/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include <stdexcept>

#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "MRCPP/trees/FunctionNode.h"
#include "MRCPP/operators/BSOperator.h"
#include "MRCPP/utils/Plotter.h"

#include "XCFunctional.h"

using mrcpp::DerivativeOperator;
using mrcpp::FunctionNode;
using mrcpp::FunctionTree;
using mrcpp::Plotter;
using mrcpp::FunctionTreeVector;
using mrcpp::Printer;
using mrcpp::Timer;

using Eigen::MatrixXi;
using Eigen::VectorXi;
using Eigen::MatrixXd;
using Eigen::VectorXd;

static int xc_iteration = 0;

namespace {
MatrixXi build_output_mask(bool is_lda, bool is_spin_sep, int order);
VectorXi build_density_mask(bool is_lda, bool is_spin_sep, int order);
void fill_output_mask(MatrixXi &mask, int start);


MatrixXi build_output_mask(bool is_lda, bool is_spin_sep, int order) {
    int start = 2;
    bool is_gga = not is_lda;
    MatrixXi mask(1,1) ;
    mask << 1;
    switch (order) {
    case 0:
        break;
    case 1:
        if (is_lda and is_spin_sep) {
            mask.resize(2,1);
            mask << 1, 2;
        } else if (is_gga and not is_spin_sep) {
            mask.resize(4,1);
            mask << 1, 2, 3, 4;
        } else if (is_gga and is_spin_sep) {
            mask.resize(8,1);
            mask << 1, 2, 3, 4, 5, 6, 7, 8;
        }
        break;
    case 2:
        if (is_lda and is_spin_sep) {
            start = 3;
            mask.resize(2,2);
        } else if (is_gga and not is_spin_sep) {
            start = 5;
            mask.resize(4,4);
        } else if (is_gga and is_spin_sep) {
            start = 9;
            mask.resize(8,8);
        }
        fill_output_mask(mask, start);
        break;
    default:
        MSG_ABORT("Not implemented");
    }
    return mask;
}

VectorXi build_density_mask(bool is_lda, bool is_spin_sep, int order) {
    bool is_gga = not is_lda;
    VectorXi mask(1) ;
    switch (order) {
    case 0:
    case 1:
        mask(0) = -1;
        break;
    case 2:
        mask(0) = 0;
        if (is_lda and is_spin_sep) {
            mask.resize(2);
            mask << 0, 1;
        } else if (is_gga and not is_spin_sep) {
            mask.resize(4);
            mask << 0, 1, 2, 3;
        } else if (is_gga and is_spin_sep) {
            mask.resize(8);
            mask << 0, 1, 2, 3, 4, 5, 6, 7;
        }
        break;
    default:
        MSG_ABORT("Not implemented");
    }
    return mask;
}

void fill_output_mask(MatrixXi &mask, int value) {
    for (int i = 0; i < mask.rows(); i++) {
        mask(i,i) = value;
        value++;
        for (int j = i+1; j < mask.cols(); j++) {
            mask(i,j) = value;
            mask(j,i) = value;
            value++;
        }
    }
}
}

namespace mrdft {

/** @brief Constructor
 *
 * Initializes the new functional
 *
 * @param[in] mra Computational domain for the density grid
 * @param[in] spin Use spin-separated functionals
 *
 * The MRA and spin parameters are const and cannot be changed after construction.
 * The constructor allocates a derivative operator (must be ABGV_00 since the grids
 * are fixed) and the density FunctionTrees (rho_a and rho_b for spin DFT, otherwise
 * only rho_t).
 */
XCFunctional::XCFunctional(mrcpp::MultiResolutionAnalysis<3> &mra, bool spin)
        : spin_separated(spin)
        , MRA(mra)
        , use_gamma(false)
        , cutoff(-1.0)
        , functional(xc_new_functional())
        , derivative(nullptr) {
    derivative = new mrcpp::BSOperator<3>(MRA, 1);
}

/** @brief Destructor */
XCFunctional::~XCFunctional() {
    xc_free_functional(functional);
    if (rho_a.size() > 0) mrcpp::clear(rho_a, true);
    if (rho_b.size() > 0) mrcpp::clear(rho_b, true);
    if (rho_t.size() > 0) mrcpp::clear(rho_t, true);
    if (derivative != nullptr) delete derivative;
}

/** @brief User-friendly setup of the xcfun calculation
 *
 *  @param[in] order Functional derivative order (1 for potential, 2 for hessian, ...)
 *
 * Prepare the XCFun object for evaluation based on the chosen parameters.
 */
void XCFunctional::evalSetup(int ord) {
    unsigned int mode = 1; //!< only partial derivative mode implemented
    order = ord; //!< update the order parameter in the object
    unsigned int func_type = isGGA(); //!< only LDA and GGA supported for now
    unsigned int dens_type = 1 + isSpinSeparated();  //!< only n (dens_type = 1) or alpha & beta (denst_type = 2) supported now.
    unsigned int laplacian = 0; //!< no laplacian
    unsigned int kinetic = 0;   //!< no kinetic energy density
    unsigned int current = 0;   //!< no current density
    unsigned int exp_derivative = not useGamma(); //!< use gamma or explicit derivatives
    if (isLDA()) exp_derivative = 0; //!< fall back to gamma-type derivatives if LDA (bad hack: no der are actually needed here!)
    xc_user_eval_setup(functional, order, func_type, dens_type, mode, laplacian, kinetic, current, exp_derivative);
}

/** @brief Set density functional
 *
 * @param[in] name The name of the chosen functional
 * @param[in] coef The amount of the chosen functional
 *
 * For each functional part in calculation a corresponding token is
 * created in xcfun. All functionals are added on top of each other
 * with the given coefficient.
 */
void XCFunctional::setFunctional(const std::string &name, double coef) {
    xc_set(functional, name.c_str(), coef);
}

/** @brief Check whether the density has been computed
 *
 * Checks if the required input densities have been computed and are valid
 * function representations. For spin separated functionals both the alpha
 * and beta densities are tested, otherwise only the total density.
 */
bool XCFunctional::hasDensity(int n_dens) const {
    bool out = true;
    if (isSpinSeparated()) {
        out = checkDensity(rho_a, n_dens);
        out = checkDensity(rho_b, n_dens);
    } else {
        out = checkDensity(rho_t, n_dens);
    }
    return out;
}

/** @brief Allocate the elements of the density vector arrays
 *
 * The density arrays which will contain the GS and the perturbed
 * densities are allocated as empty FunctinTree objects. The pointers
 * ae then stored in the corresponding FunctionTreeVector
 */
void XCFunctional::allocateDensities() {
    for (int i = 0; i < nDensities; i++) {
        if (isSpinSeparated()) {
            FunctionTree<3> *temp_a = new FunctionTree<3>(MRA);
            FunctionTree<3> *temp_b = new FunctionTree<3>(MRA);
            rho_a.push_back(std::make_tuple(1.0, temp_a));
            rho_b.push_back(std::make_tuple(1.0, temp_b));
        } else {
            FunctionTree<3> *temp_t = new FunctionTree<3>(MRA);
            rho_t.push_back(std::make_tuple(1.0, temp_t));
        }
    }
    
}

/** @brief Check whether the density vector has all required components
 *
 * Checks if the required input densities have been computed and are
 * valid function representations. The density vector shall contain a
 * number of densities up to the requested order (or higher)
 *
 * LUCA: I envisage a possible problem here. We might want to
 * initialize the functional with a higher order although that will
 * not be used in reality.
 *
 */
bool XCFunctional::checkDensity(FunctionTreeVector<3> density, int n_dens) const {
    bool out = true;
    if (density.size() < n_dens) out = false; //Not all required densities are present and we should make them
    for (int i = 0; i < n_dens; i++) {
        FunctionTree<3> &dens = mrcpp::get_func(density, i);
        if (dens.getSquareNorm() < 0.0) out = false;
    }
    return out;
}

/** @brief Return FunctionTree for the input density
 *
 * @param[in] type Which density to return (alpha, beta or total)
 *
 * Returns a reference to the internal vector of density functions so that it can be
 * computed by the host program. This needs to be done before setup().
 */
FunctionTreeVector<3> &XCFunctional::getDensityVector(DENSITY::DensityType type) {
    switch (type) {
        case DENSITY::DensityType::Total:
            return rho_t;
        case DENSITY::DensityType::Alpha:
            return rho_a;
        case DENSITY::DensityType::Beta:
            return rho_b;
        default:
            MSG_ABORT("Invalid density type");
    }
    MSG_ABORT("Total density functions not properly initialized");
}

/** @brief Return FunctionTree for the input density
 *
 * @param[in] type Which density to return (alpha, beta or total)
 *
 * Returns a reference to the internal vector of density functions so that it can be
 * computed by the host program. This needs to be done before setup().
 */
FunctionTree<3> &XCFunctional::getDensity(DENSITY::DensityType type, int index) {
    FunctionTreeVector<3> dens_vec = getDensityVector(type);
    if(index >= dens_vec.size()) MSG_ABORT("Vector out of bounds");
    return mrcpp::get_func(dens_vec, index);
}

void XCFunctional::setDensity(FunctionTree<3> &density, DENSITY::DensityType spin, int index) {
    FunctionTreeVector<3> &dens_vec = getDensityVector(spin);
    if (dens_vec.size() != index) {
        MSG_ABORT("Index mismatch");
    }
    dens_vec.push_back(std::make_tuple(1.0, &density));
}

/** @brief Return the number of nodes in the density grid (including branch nodes)*/
int XCFunctional::getNNodes() const {
    int nodes = 0;
    if (isSpinSeparated()) {
        const FunctionTree<3> &rho_gs_a = mrcpp::get_func(rho_a, 0);
        const FunctionTree<3> &rho_gs_b = mrcpp::get_func(rho_b, 0);
        nodes = rho_gs_a.getNNodes();
        std::cout << "nnodes " << rho_gs_a.getNNodes()    << std::endl;
        std::cout << "nnodes " << rho_gs_b.getNNodes()    << std::endl;        
        if (nodes != rho_gs_b.getNNodes()) MSG_ERROR("Alpha and beta grids not equal");
    } else {
        const FunctionTree<3> &rho_gs_t = mrcpp::get_func(rho_t, 0);
        nodes = rho_gs_t.getNNodes();
    }
    return nodes;
}

/** @brief Return the number of grid points in the density grid (only leaf nodes)*/
int XCFunctional::getNPoints() const {
    int nodes = 0;
    int points = 0;
    if (isSpinSeparated()) {
        const FunctionTree<3> &rho_gs_a = mrcpp::get_func(rho_a, 0);
        const FunctionTree<3> &rho_gs_b = mrcpp::get_func(rho_b, 0);
        nodes = rho_gs_a.getNEndNodes();
        points = rho_gs_a.getTDim()*rho_gs_a.getKp1_d();
        std::cout << "enodes  " << rho_gs_a.getNEndNodes() << std::endl;
        std::cout << "enodes  " << rho_gs_b.getNEndNodes() << std::endl;
        if (nodes != rho_gs_b.getNEndNodes()) MSG_ERROR("Alpha and beta grids not equal");
    } else {
        const FunctionTree<3> &rho_gs_t = mrcpp::get_func(rho_t, 0);
        nodes = rho_gs_t.getNEndNodes();
        points = rho_gs_t.getTDim() * rho_gs_t.getKp1_d();
    }
    return nodes * points;
}

/** @brief Construct initial density grid based on nuclear charge and position
 *
 * @param[in] Z Nuclear charge
 * @param[in] R Nuclear position
 *
 * This will _extend_ the density grid with extra nodes surrounding a nuclear site.
 * The refinement level depends on the nuclear charge (higher charge -> more
 * refinement). This is meant only to provide an initial guess for the grid, there
 * are no guarantees regarding the precision.
 */
void XCFunctional::buildGrid(double Z, const mrcpp::Coord<3> &R) {
    for (int i = 0; i < 5; i++) {
        mrcpp::GaussFunc<3> gauss(std::pow(Z, i), 1.0, R);
        if (isSpinSeparated()) {
            FunctionTree<3> &rho_gs_a = mrcpp::get_func(rho_a, 0);
            FunctionTree<3> &rho_gs_b = mrcpp::get_func(rho_b, 0);
            mrcpp::build_grid(rho_gs_a, gauss);
            mrcpp::build_grid(rho_gs_b, gauss);
            copyGrid(rho_a);
            copyGrid(rho_b);
        } else {
            FunctionTree<3> &rho_gs_t = mrcpp::get_func(rho_t, 0);
            mrcpp::build_grid(rho_gs_t, gauss);
            copyGrid(rho_t);
        }
    }
}

/** @brief Copy the grid from the ground state densities to the perturbed ones
 *
 * In contracted mode XCFun needs pointwise information about all
 * functions, therefore all grids must be identical to the ground
 * state grid
 */
void XCFunctional::copyGrid(FunctionTreeVector<3> densities) {
    mrcpp::FunctionTree<3> &rho0 = mrcpp::get_func(densities, 0);
    for (int i = 1; i < densities.size(); i++) {
        mrcpp::FunctionTree<3> &rho = mrcpp::get_func(densities, i);
        mrcpp::build_grid(rho, rho0);
    }
}

    /*
      void XCFunctional::refineGrid(double prec, bool abs_prec) {
      if (not hasDensity()) return;
      
      double scale = 1.0;
      if (isSpinSeparated()) {
      if (rho_a.size() == 0) MSG_ABORT("Uninitialized alpha density");
      if (rho_b.size() == 0) MSG_ABORT("Uninitialized beta density");
      if (abs_prec) scale = rho_a[0]->integrate() + rho_b[0]->integrate();
      mrcpp::refine_grid(rho_a, prec / scale);
      mrcpp::refine_grid(rho_b, prec / scale);
      
      // Extend to union grid
      int nNodes = 1;
      while (nNodes > 0) {
      int nAlpha = mrcpp::refine_grid(rho_a, rho_b);
      int nBeta = mrcpp::refine_grid(rho_b, rho_a);
      nNodes = nAlpha + nBeta;
      }
      } else {
      if (rho_t == nullptr) MSG_ABORT("Uninitialized total density");
      if (abs_prec) scale = rho_t->integrate();
      mrcpp::refine_grid(*rho_t, prec / scale);
      }
      }
    */


/** @brief Remove all grid refinement for a given density vector
 *
 * This will _remove_ all existing grid refinement and leave only root nodes
 * in the density grids.
 */
void XCFunctional::clearGrid(FunctionTreeVector<3> densities) {
    for (int i = 0; i < densities.size(); i++) {
        FunctionTree<3> &rho = mrcpp::get_func(densities, i);
        rho.clear();
    }
}

void XCFunctional::setup() {
    if (not hasDensity(order)) {
        MSG_ABORT("Not enough density functions initialized");
    }
    setupGradient();
    setupXCInput();
    setupXCOutput();
    setupXCDensityVariables();
}

void XCFunctional::setupGradient() {
    if(isLDA()) return;
    for (int i = 0; i < nDensities; i++) {
        if(isSpinSeparated()) {
            FunctionTreeVector<3> temp_a = mrcpp::gradient(*derivative, getDensity(DENSITY::DensityType::Alpha, i));
            FunctionTreeVector<3> temp_b = mrcpp::gradient(*derivative, getDensity(DENSITY::DensityType::Beta, i));
            grad_a.insert(grad_a.end(), temp_a.begin(), temp_a.end());
            grad_b.insert(grad_b.end(), temp_b.begin(), temp_b.end());
        } else {
            FunctionTreeVector<3> temp_t = mrcpp::gradient(*derivative, getDensity(DENSITY::DensityType::Total, i));
            grad_t.insert(grad_t.end(), temp_t.begin(), temp_t.end());
        }
    }
    
}

/** @brief Cleanup after xcfun evaluation
 *
 * This deallocates the memory used by both the input and output functions
 * to xcfun. The output must be collected (calcEnergy() and calcPotential())
 * before this function is called. The density grids will keep their current
 * refinement, but their MW coefs will be cleared.
 */
void XCFunctional::clear() {
    mrcpp::clear(xcInput, false);
    mrcpp::clear(xcDensity, false);
    mrcpp::clear(xcOutput, true);
    mrcpp::clear(grad_a, true);
    mrcpp::clear(grad_b, true);
    mrcpp::clear(grad_t, true);
    mrcpp::clear(gamma, true);
    clearGrid(rho_a); // We want to keep empty the density trees
    clearGrid(rho_b); // We want to keep empty the density trees
    clearGrid(rho_t); // We want to keep empty the density trees
}

/** @brief Allocate input arrays for xcfun
 *
 * Based on the xcfun setup, the requested array of FunctionTrees(s)
 * is allocared and its pointers assigned to the required input
 * functions.
 */
void XCFunctional::setupXCInput() {
    if (xcInput.size() != 0) MSG_ERROR("XC input not empty");

    int nInp = getInputLength();
    int nUsed = setupXCInputDensity();
    if (isGGA()) nUsed += setupXCInputGradient();

    if (nInp != nUsed) MSG_ERROR("Mismatch between used vs requested " << nUsed << " : " << nInp);
}

/** @brief Allocate input arrays for xcfun
 *
 * Based on the xcfun setup, the requested array of FunctionTrees(s)
 * is allocared and its pointers assigned to the required input
 * functions.
 */
void XCFunctional::setupXCDensityVariables() {
    if (xcDensity.size() != 0) MSG_ERROR("XC density vars not empty");

    int n_dens = getDensityLength();
    
    if(getDensityLength() > 0) {
        if (isSpinSeparated()) {
            xcDensity.push_back(rho_a[1]);
            xcDensity.push_back(rho_b[1]);
        } else {
            xcDensity.push_back(rho_t[1]);
        }
        if (isGGA()) {
            if (isSpinSeparated()) {
                xcDensity.insert(xcDensity.begin() + 2, grad_a.begin() + 3, grad_a.begin() + 6);
                xcDensity.insert(xcDensity.begin() + 5, grad_b.begin() + 3, grad_b.begin() + 6);
            } else {
                xcDensity.insert(xcDensity.begin() + 1, grad_t.begin() + 3, grad_t.begin() + 6);
            }
        }
    }
    //	std::cout << "Plotting densities" << std::endl;
    //	plot_function_tree_vector(xcDensity, "Dens_");

    if (n_dens != xcDensity.size()) MSG_ABORT("Mismatch between used vs requested " << n_dens << " : " << xcDensity.size());

}

    /*
void XCFunctional::plot_function_tree_vector(FunctionTreeVector<3> &functions, std::string prefix) {
    
    int nPts = 10000;                               // Number of points
    double a[3] = { 0.0,  0.0, 16.0};               // Start point of plot
	double b[3] = { 0.0, 16.0,  0.0};               // End point of plot
    double o[3] = { 0.0, -8.0, -8.0};               // Origin of plot
    mrcpp::Plotter<3> plot;                         // Plotter of 3D functions
    plot.setNPoints(nPts);                          // Set number of points
	plot.setRange(a, b, o);                         // Set plot range

    for (int i = 0; i < functions.size(); i++) {
        mrcpp::FunctionTree<3> &func = mrcpp::get_func(functions, i);
        std::string name= prefix + std::to_string(i) + "_iter_" + std::to_string(xc_iteration);
        std::cout << name << std::endl;
        std::cout << func << std::endl;
        plot.surfPlot(func, name);
    }
    
    xc_iteration++;

}
    */

/** @brief Sets xcInput pointers for the density
 *
 * Returns the nr. of pointers used for sanity checking.
 */
int XCFunctional::setupXCInputDensity() {
    int nUsed = 0;
    if (isSpinSeparated()) {
        if (rho_a.size() <= 0) MSG_ABORT("Invalid alpha density");
        if (rho_b.size() <= 0) MSG_ABORT("Invalid beta density");
        xcInput.push_back(rho_a[0]);
        xcInput.push_back(rho_b[0]);
        nUsed = 2;
    } else {
        if (rho_t.size() <= 0) MSG_ABORT("Invalid total density");
        xcInput.push_back(rho_t[0]);
        nUsed = 1;
    }
    return nUsed;
}

/** @brief Sets xcInput pointers for the gradient(s)
 *
 * Returns the nr. of pointers used for sanity checking.
 */
int XCFunctional::setupXCInputGradient() {
    int nUsed = 0;
    int nFetch = 3;
    if (useGamma()) {
        NOT_IMPLEMENTED_ABORT;
    } else {
        if (isSpinSeparated()) {
            xcInput.insert(xcInput.end(), grad_a.begin(), grad_a.begin() + nFetch);
            xcInput.insert(xcInput.end(), grad_b.begin(), grad_b.begin() + nFetch);
            nUsed += 2 * nFetch;
        } else {
            xcInput.insert(xcInput.end(), grad_t.begin(), grad_t.begin() + nFetch);
            nUsed += nFetch;
        }
    }
    return nUsed;
}

/** @brief Allocate output arrays for xcfun
 *
 * Based on the xcfun setup, the requested array of FunctionTrees(s)
 * is allocated and the function objects are created, borrowing the
 * grid from the input density.
 */
void XCFunctional::setupXCOutput() {
    if (xcOutput.size() != 0) MSG_ERROR("XC output not empty");
    if (xcInput.size() == 0) MSG_ERROR("XC input not initialized");

    // Fetch density grid to copy
    FunctionTree<3> &grid = mrcpp::get_func(xcInput, 0);

    int nOut = getContractedLength() + 1; //We always keep the functional as well
    for (int i = 0; i < nOut; i++) {
        auto *tmp = new FunctionTree<3>(MRA);
        mrcpp::copy_grid(*tmp, grid);
        xcOutput.push_back(std::make_tuple(1.0, tmp));
    }
}

/** @brief Evaluation of the functional and its derivatives
 *
 * The data contained in the xcInput is converted in matrix form and fed to
 * the functional. The output matrix is then converted back to function form.
 */
void XCFunctional::evaluate() {
    if (xcInput.size() == 0) MSG_ERROR("XC input not initialized");
    if (xcOutput.size() == 0) MSG_ERROR("XC output not initialized");

    Timer timer;
    println(2, "Evaluating");

    int nInp = getInputLength();          // Input parameters to XCFun
    int nOut = getOutputLength();         // Output parameters from XCFun
    int nCon = getContractedLength();     // Contracted parameters to XCPotential
    int nFcs = nCon + 1;                  // One extra function for the energy density
    int nPts = getNodeLength();           // Number of gridpoints in a node

#pragma omp parallel firstprivate(nInp, nOut)
    {
        int nNodes = mrcpp::get_func(xcInput, 0).getNEndNodes();
#pragma omp for schedule(guided)
        for (int n_idx = 0; n_idx < nNodes; n_idx++) {
            MatrixXd inpData, outData, conData, denData;
            inpData = MatrixXd::Zero(nPts, nInp);
            outData = MatrixXd::Zero(nPts, nOut);
            conData = MatrixXd::Zero(nPts, nFcs);
            compressNodeData(n_idx, nInp, xcInput, inpData);
            evaluateBlock(inpData, outData);
            conData.col(0) = outData.col(0); // we always keep the energy functional
            contractNodeData(n_idx, nPts, outData, conData);
            expandNodeData(n_idx, nFcs, xcOutput, conData);
        }
    }
    for (int i = 0; i < nFcs; i++) {
        FunctionTree<3> &func = mrcpp::get_func(xcOutput, i);
        func.mwTransform(mrcpp::BottomUp);
        func.calcSquareNorm();
        std::cout << "Potential norm " << i << " " << func.getSquareNorm() 
				  << " nEndNodes " << func.getNEndNodes() << std::endl;
    }

    for (int i = 0; i < xcDensity.size(); i++) {
        FunctionTree<3> &func = mrcpp::get_func(xcDensity, i);
        std::cout << "Density norm " << i << " " << func.getSquareNorm() 
				  << " nEndNodes " << func.getNEndNodes() << std::endl;
    }
    auto n = mrcpp::get_n_nodes(xcOutput);
    auto m = mrcpp::get_size_nodes(xcOutput);
    auto t = timer.elapsed();
    mrcpp::print::tree(2, "XC evaluate xcfun", n, m, t);
    printout(2, std::endl);
}

void XCFunctional::contractNodeData(int node_index, int n_points, MatrixXd &out_data, MatrixXd &con_data) {

    MatrixXi output_mask = build_output_mask(isLDA(), isSpinSeparated(), order);
    VectorXi density_mask = build_density_mask(isLDA(), isSpinSeparated(), order);
    if(output_mask.cols() != density_mask.size()) MSG_ABORT("Inconsistent lengths");

    VectorXd cont_i;
    VectorXd cont_ij;
    for (int i = 0; i < output_mask.rows(); i++) {
        cont_i = VectorXd::Zero(n_points);
        for (int j = 0; j < output_mask.cols(); j++) {
            cont_ij = VectorXd::Zero(n_points);
            int out_index = output_mask(i,j);
            int den_index = density_mask(j);
            //if(den_index >= 2) {  //WARNING: This HACK works for Open shell only
            //        continue;
            //} else if(den_index >= 0) {
            if(den_index >= 0) {
                FunctionTree<3> &dens_func = mrcpp::get_func(xcDensity, den_index);
                FunctionNode<3> &dens_node = dens_func.getEndFuncNode(node_index);
                VectorXd dens_i;
                dens_node.getValues(dens_i);
                cont_ij = out_data.col(out_index).array() * dens_i.array();
            } else {
                cont_ij = out_data.col(out_index);
            }
            cont_i += cont_ij;
        }
        con_data.col(i+1) = cont_i; //The first column contains the energy functional
    }
}


/** \brief Evaluates XC functional and derivatives
 *
 * Computes the alpha and beta exchange-correlation functionals and
 * their derivatives.  The electronic density (total/alpha/beta) and their gradients are
 * given as input. Results are then stored in the xcfun output
 * functions. Higher order derivatives can be computed changing the parameter k.
 *
 * XCFunctional output (with k=1 and explicit derivatives):
 *
 * LDA: \f$ \left(F_{xc}, \frac{\partial F_{xc}}{\partial \rho}\right) \f$
 *
 * GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho},
 *  \frac{\partial F_{xc}}{\partial \rho_x},
 *  \frac{\partial F_{xc}}{\partial \rho_y},
 *  \frac{\partial F_{xc}}{\partial \rho_z}\right) \f$
 *
 * Spin LDA: \f$ \left(F_{xc}, \frac{\partial F_{xc}}{\partial \rho^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho^\beta}\right) \f$
 *
 * Spin GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_x^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_y^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_z^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_x^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_y^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_z^\beta}
 *  \right) \f$
 *
 * XCFunctional output (with k=1 and gamma-type derivatives):
 *
 * GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho},
 *  \frac{\partial F_{xc}}{\partial \gamma} \f$
 *
 * Spin GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho^\beta },
 *  \frac{\partial F_{xc}}{\partial \gamma^{\alpha \alpha}},
 *  \frac{\partial F_{xc}}{\partial \gamma^{\alpha \beta }},
 *  \frac{\partial F_{xc}}{\partial \gamma^{\beta  \beta }}
 *  \right) \f$
 *
 * The points are passed with through a matrix of dimension nInp x nPts
 * where nInp is the number of input data required for a single evaluation
 * and nPts is the number of points requested. Similarly the output is provided
 * as a matrix nOut x nPts.
 *
 * param[in] inp Input values
 * param[out] out Output values
 *
 */
void XCFunctional::evaluateBlock(MatrixXd &inp, MatrixXd &out) const {
    if (inp.cols() != getInputLength()) MSG_ERROR("Invalid input");

    int nInp = getInputLength();
    int nOut = getOutputLength();

    int nPts = inp.rows();
    out = MatrixXd::Zero(nPts, nOut);

    double iDat[nInp];
    double oDat[nOut];

    for (int i = 0; i < nPts; i++) {
        if (inp(i, 0) > cutoff) {
            for (int j = 0; j < nInp; j++) iDat[j] = inp(i, j);
            xc_eval(functional, iDat, oDat);
            for (int j = 0; j < nOut; j++) out(i, j) = oDat[j];
        } else {
            for (int j = 0; j < nOut; j++) out(i, j) = 0.0;
        }
    }
}

/** @brief Converts data from a FunctionNode to a matrix
 *
 * The FunctionNode(s) row data is packed into a matrix whose
 * dimensions are the overall number of grid points (nCoefs) and the
 * number of functions (nFuncs).
 *
 * parma[in] n the Index of the requested node
 * param[in] nFuncs The number of functions
 * param[in] trees The array of FunctionTree(s)
 * param[in] data The matrix object
 */
void XCFunctional::compressNodeData(int n, int nFuncs, FunctionTreeVector<3> trees, MatrixXd &data) {
    if (trees.size() == 0) MSG_ERROR("Invalid input");

    FunctionTree<3> &tree = mrcpp::get_func(trees, 0);
    int nCoefs = tree.getTDim() * tree.getKp1_d();
    data = MatrixXd::Zero(nCoefs, nFuncs);

    for (int i = 0; i < nFuncs; i++) {
        FunctionNode<3> &node = mrcpp::get_func(trees, i).getEndFuncNode(n);
        VectorXd col_i;
        node.getValues(col_i);
        data.col(i) = col_i;
    }
}

/** @brief Converts data from a matrix to a FunctionNode
 *
 * The matrix containing the output from xcfun is converted back to the corresponding FunctionNode(s). The matrix
 * dimensions are the overall number of grid points (nCoefs) and the number of functions (nFuncs).
 *
 * parma[in] n the Index of the requested node
 * param[in] nFuncs The number of functions
 * param[in] trees The array of FunctionTree(s)
 * param[in] data The matrix object
 */
void XCFunctional::expandNodeData(int n, int nFuncs, FunctionTreeVector<3> trees, MatrixXd &data) {
    if (trees.size() == 0) MSG_ERROR("Invalid input");

    for (int i = 0; i < nFuncs; i++) {
        VectorXd col_i = data.col(i);
        FunctionNode<3> &node = mrcpp::get_func(trees, i).getEndFuncNode(n);
        node.setValues(col_i);
    }
}

/** @brief Computes the XC energy as the integral of the energy density */
double XCFunctional::calcEnergy() {
    if (xcOutput.size() == 0) MSG_ERROR("XC output not initialized");

    Timer timer;
    FunctionTree<3> &E_dens = mrcpp::get_func(xcOutput, 0);
    double energy = E_dens.integrate();
    mrcpp::print::tree(2, "XC energy", E_dens, timer);
    return energy;
}

/** @brief Compute the XC potential(s)
 *
 * Combines the xcfun output functions into the final XC potential functions.
 * Different calls for LDA and GGA, and for gamma-type vs explicit derivatives.
 */
FunctionTreeVector<3> XCFunctional::calcPotential() {
    FunctionTreeVector<3> xc_pot;
    if (xcOutput.size() == 0) MSG_ERROR("XC output not initialized");

    Timer timer;
    if (isLDA()) {
        calcPotentialLDA(xc_pot);
    } else if (isGGA()) {
        calcPotentialGGA(xc_pot);
    } else {
        MSG_ABORT("Invalid functional type");
    }
    auto n = mrcpp::get_n_nodes(xc_pot);
    auto m = mrcpp::get_size_nodes(xc_pot);
    auto t = timer.elapsed();
    mrcpp::print::tree(2, "XC potential", n, m, t);
    //	std::cout << "Plotting potentials" << std::endl;
    //	plot_function_tree_vector(xc_pot, "Potential_");
    return xc_pot;
}

/** @brief Potential calculation for LDA functionals
 *
 * The potential conicides with the xcfun output, which is then
 * deep copied into the corresponding potential functions.
 */
void XCFunctional::calcPotentialLDA(FunctionTreeVector<3> &potentials) {
    int nPotentials = isSpinSeparated()? 2 : 1;
    int iStart = 1;
    for (int i = 0; i < nPotentials; i++) {
        FunctionTree<3> &out_i = mrcpp::get_func(xcOutput, i + iStart);
        auto *pot = new FunctionTree<3>(MRA);
        mrcpp::copy_grid(*pot, out_i);
        mrcpp::copy_func(*pot, out_i);
        potentials.push_back(std::make_tuple(1.0, pot));
    }
}

/** @brief Potential calculation for GGA functionals
 *
 * The potential functions are assembled from the xcfun output functions.
 * The method used depends on whether the functional is spin-separated
 * and whether explicit or gamma-type derivatives have been used in xcfun.
 */
void XCFunctional::calcPotentialGGA(FunctionTreeVector<3> &potentials) {
    if (useGamma()) NOT_IMPLEMENTED_ABORT;
    FunctionTree<3> *pot;
    if (isSpinSeparated()) {
        FunctionTree<3> &df_da = mrcpp::get_func(xcOutput, 1);
        FunctionTree<3> &df_db = mrcpp::get_func(xcOutput, 2);
        FunctionTreeVector<3> df_dga(xcOutput.begin() + 3, xcOutput.begin() + 6);
        FunctionTreeVector<3> df_dgb(xcOutput.begin() + 6, xcOutput.begin() + 9);
        pot = calcPotentialGGA(df_da, df_dga);
        potentials.push_back(std::make_tuple(1.0, pot));
        pot = calcPotentialGGA(df_db, df_dgb);
        potentials.push_back(std::make_tuple(1.0, pot));
    } else {
        FunctionTree<3> &df_dt = mrcpp::get_func(xcOutput, 1);
        FunctionTreeVector<3> df_dgt(xcOutput.begin() + 2, xcOutput.begin() + 5);
        pot = calcPotentialGGA(df_dt, df_dgt);
        potentials.push_back(std::make_tuple(1.0, pot));
    }
    pot = nullptr;
}

/** @brief XC potential calculation
 *
 * @param[in] df_drho Functional derivative wrt rho
 * @param[in] df_dgr  Functional_derivative wrt grad_rho
 *
 * Computes the XC potential for explicit derivatives.
 */
FunctionTree<3> *XCFunctional::calcPotentialGGA(FunctionTree<3> &df_drho, FunctionTreeVector<3> &df_dgr) {
    FunctionTreeVector<3> funcs;
    funcs.push_back(std::make_tuple(1.0, &df_drho));

    auto *tmp = new FunctionTree<3>(MRA);
    mrcpp::divergence(*tmp, *derivative, df_dgr);
    funcs.push_back(std::make_tuple(-1.0, tmp));

    auto *V = new FunctionTree<3>(MRA);
    mrcpp::build_grid(*V, funcs);
    mrcpp::add(-1.0, *V, funcs);
    delete tmp;
    return V;
}

int XCFunctional::getNodeLength() {
    FunctionTree<3> &tree = mrcpp::get_func(xcInput, 0);
    return tree.getTDim() * tree.getKp1_d();
}

int XCFunctional::getContractedLength() const {
    int length = -1;
    if(isLDA()) length = 1; //only one function needed for LDA
    if(isGGA()) length = 4; //four functions for LDA
    if(isSpinSeparated()) length *= 2; //twice as many for spin-separated functionals;
    return length;
}

int XCFunctional::getDensityLength() const {
    int length = -1;
    if(this->order < 2) {
        length = 0; // no contractions for order 0 or 1
    } else if (order == 2) {
        length = getContractedLength(); //same length as contracted fcns for order = 2
    }
    return length;
}

} // namespace mrdft


