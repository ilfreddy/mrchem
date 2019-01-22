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

#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "MRCPP/trees/FunctionNode.h"

#include "XCFunctional.h"

using mrcpp::DerivativeOperator;
using mrcpp::FunctionNode;
using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;
using mrcpp::Printer;
using mrcpp::Timer;

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
    derivative = new mrcpp::ABGVOperator<3>(MRA, 0.0, 0.0);
    if (isSpinSeparated()) {
        rho_a = new FunctionTree<3>(MRA);
        rho_b = new FunctionTree<3>(MRA);
    } else {
        rho_t = new FunctionTree<3>(MRA);
    }
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
    if (mod != 1 and mod != 3) {
        MSG_FATAL("Only contracted or partial derivative modes are available")
    }
    mode = mod;
    order = ord; //!< update the order parameter in the object
    unsigned int func_type = isGGA(); //!< only LDA and GGA supported for now
    unsigned int dens_type = 1 + isSpinSeparated();  //!< only n (dens_type = 1) or alpha & beta (denst_type = 2) supported now.
    unsigned int laplacian = 0; //!< no laplacian
    unsigned int kinetic = 0;   //!< no kinetic energy density
    unsigned int current = 0;   //!< no current density
    unsigned int exp_derivative = not useGamma(); //!< use gamma or explicit derivatives
    if (isLDA()) exp_derivative = 0; //!< fall back to gamma-type derivatives if LDA (bad hack: no der are actually needed here!)
    xc_user_eval_setup(functional, order, func_type, dens_type, mode_type, laplacian, kinetic, current, exp_derivative);
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
bool XCFunctional::hasDensity() const {
    if (isSpinSeparated()) {
        out = checkDensity(rho_a);
        out = checkDensity(rho_b);
    } else {
        out = checkDensity(rho_t);
    }
    return out;
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
bool checkDensity(FunctionTreeVector<3> density) const {
    bool out = true;
    if (density.size() < nDensities) out = false;
    for (int i = 0; i < order; i++) {
        FunctionTree<3> dens = mrcpp::get_func(density, i);
        if (temp->getSquareNorm() < 0.0) out = false;
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
FunctionTreeVector<3> & XCFunctional::getDensity(DensityType type) {
    switch (type) {
        case DensityType::Total:
            if (checkDensity(rho_t)) return rho_t;
        case DensityType::Alpha:
            if (checkDensity(rho_a)) return rho_a;
        case DensityType::Beta:
            if (checkDensity(rho_b)) return rho_b;
        default:
            MSG_FATAL("Invalid density type");
    }
    MSG_FATAL("Total density functions not properly initialized");
}

/** @brief Return the number of nodes in the density grid (including branch nodes)*/
int XCFunctional::getNNodes() const {
    int nodes = 0;
    if (isSpinSeparated()) {
        FunctionTree<3> &rho_gs_a = mrcpp::get_func(rho_a, 0);
        FunctionTree<3> &rho_gs_b = mrcpp::get_func(rho_b, 0);
        nodes = rho_gs_a->getNNodes();
        if (nodes != rho_gs_b->getNNodes()) MSG_ERROR("Alpha and beta grids not equal");
    } else {
        FunctionTree<3> &rho_gs_t = mrcpp::get_func(rho_t, 0);
        nodes = rho_gs_t->getNNodes();
    }
    return nodes;
}

/** @brief Return the number of grid points in the density grid (only leaf nodes)*/
int XCFunctional::getNPoints() const {
    int nodes = 0;
    int points = 0;
    if (isSpinSeparated()) {
        FunctionTree<3> &rho_gs_a = mrcpp::get_func(rho_a, 0);
        FunctionTree<3> &rho_gs_b = mrcpp::get_func(rho_b, 0);
        nodes = rho_gs_a->getNEndNodes();
        points = rho_gs_a->getTDim()*rho_a->getKp1_d();
        if (nodes != rho_gs_b->getNEndNodes()) MSG_ERROR("Alpha and beta grids not equal");
    } else {
        FunctionTree<3> &rho_gs_t = mrcpp::get_func(rho_t, 0);
        nodes = rho_t->getNEndNodes();
        points = rho_t->getTDim() * rho_t->getKp1_d();
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
            distributeGrid(rho_a)
            distributeGrid(rho_b)
        } else {
            FunctionTree<3> &rho_gs_t = mrcpp::get_func(rho_t, 0);
            mrcpp::build_grid(rho_gs_t, gauss);
            distributeGrid(rho_t)
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

/** @brief Remove all grid refinement for all density vectors
 *
 * This will _remove_ all existing grid refinement and leave only root nodes
 * in the density grids.
 */
void XCFunctional::refineGrid(double prec, bool abs_prec) {
    if (not hasDensity()) return;

    double scale = 1.0;
    if (isSpinSeparated()) {
        if (rho_a == nullptr) MSG_ABORT("Uninitialized alpha density");
        if (rho_b == nullptr) MSG_ABORT("Uninitialized beta density");
        if (abs_prec) scale = rho_a->integrate() + rho_b->integrate();
        mrcpp::refine_grid(*rho_a, prec / scale);
        mrcpp::refine_grid(*rho_b, prec / scale);

        // Extend to union grid
        int nNodes = 1;
        while (nNodes > 0) {
            int nAlpha = mrcpp::refine_grid(*rho_a, *rho_b);
            int nBeta = mrcpp::refine_grid(*rho_b, *rho_a);
            nNodes = nAlpha + nBeta;
        }
    } else {
        if (rho_t == nullptr) MSG_ABORT("Uninitialized total density");
        if (abs_prec) scale = rho_t->integrate();
        mrcpp::refine_grid(*rho_t, prec / scale);
    }
}

/** @brief Remove all grid refinement for a given density vector
 *
 * This will _remove_ all existing grid refinement and leave only root nodes
 * in the density grids.
 */
void XCFunctional::clearGrid(densities) {
    for (int i = 0; i < densities.size(); i++) {
        FunctionTree<3> &rho = mrcpp::get_func(densities, i);
        rho_a.clear();
    }
}

void XCFunctional::setup() {
    switch (mode) {
    case 1:
        setup_partial();
        break;
    case 3:
        setup_contracted();
        break;
    default:
        MSG_FATAL("Not implemented: abort");
    }
}


void XCFunctional::setup_contracted() {
    MSG_FATAL("Not implemented: abort");
}
/** @brief Prepare for xcfun evaluation
 *
 * This computes the necessary input functions (gradients and gamma) and
 * constructs empty grids to hold the output functions of xcfun, and
 * collects them in the xcInput and xcOutput vectors, respectively.
 * Assumes that the density functions rho_t or rho_a/rho_b have already
 * been computed.
 */
void XCFunctional::setup_partial() {
    if (isSpinSeparated()) {
        FunctionTreeVector &rho_gs_a = mrcpp::get_function(rho_a, 0);
        FunctionTreeVector &rho_gs_b = mrcpp::get_function(rho_b, 0);
        setup_partial(rho_gs_a, rho_gs_b);
    } else {
        FunctionTreeVector &rho_gs_t = mrcpp::get_function(rho_t, 0);
        setup_partial(rho_gs_t);
    }
}

void XCFunctional::setup_partial(FunctionTree<3> &rho_a, FunctionTree<3> &rho_b) {
    if (isGGA()) {
        grad_a = mrcpp::gradient(*derivative, rho_a);
        grad_b = mrcpp::gradient(*derivative, rho_b);
    }
    if (useGamma()) {
        FunctionTree<3> *gamma_aa = new FunctionTree<3>(MRA);
        FunctionTree<3> *gamma_ab = new FunctionTree<3>(MRA);
        FunctionTree<3> *gamma_bb = new FunctionTree<3>(MRA);
        mrcpp::build_grid(*gamma_aa, rho_a);
        mrcpp::build_grid(*gamma_ab, rho_a);
        mrcpp::build_grid(*gamma_bb, rho_a);
        mrcpp::dot(-1.0, *gamma_aa, grad_a, grad_a);
        mrcpp::dot(-1.0, *gamma_ab, grad_a, grad_b);
        mrcpp::dot(-1.0, *gamma_bb, grad_b, grad_b);
        gamma.push_back(std::make_tuple(1.0, gamma_aa));
        gamma.push_back(std::make_tuple(1.0, gamma_ab));
        gamma.push_back(std::make_tuple(1.0, gamma_bb));
    }
    setupXCInput();
    setupXCOutput();
}

void XCFunctional::setup_partial(FunctionTree<3> &rho_t) {
    if (isGGA()) {
        grad_t = mrcpp::gradient(*derivative, rho_t);
    }
    if (useGamma()) {
        FunctionTree<3> *gamma_tt = new FunctionTree<3>(MRA);
        mrcpp::build_grid(*gamma_tt, rho_t);
        mrcpp::dot(-1.0, *gamma_tt, grad_t, grad_t);
        gamma.push_back(std::make_tuple(1.0, gamma_tt));
    }
    setupXCInput();
    setupXCOutput();
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
    mrcpp::clear(xcOutput, true);
    mrcpp::clear(grad_a, true);
    mrcpp::clear(grad_b, true);
    mrcpp::clear(grad_t, true);
    mrcpp::clear(gamma, true);
    clearGrid();

    // Clear MW coefs but keep the grid
    // if (rho_a != nullptr) mrcpp::clear_grid(*rho_a);
    // if (rho_b != nullptr) mrcpp::clear_grid(*rho_b);
    // if (rho_t != nullptr) mrcpp::clear_grid(*rho_t);
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

/** @brief Sets xcInput pointers for the density
 *
 * Returns the nr. of pointers used for sanity checking.
 */
int XCFunctional::setupXCInputDensity() {
    int nUsed = 0;
    if (isSpinSeparated()) {
        if (rho_a == nullptr) MSG_FATAL("Invalid alpha density");
        if (rho_b == nullptr) MSG_FATAL("Invalid beta density");
        xcInput.push_back(rho_a[0]);
        xcInput.push_back(rho_b[0]);
        nUsed = 2;
    } else {
        if (rho_t == nullptr) MSG_FATAL("Invalid total density");
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
    if (useGamma()) {
        xcInput.insert(xcInput.end(), gamma.begin(), gamma.end());
        nUsed += gamma.size();
    } else {
        if (isSpinSeparated()) {
            xcInput.insert(xcInput.end(), grad_a.begin(), grad_a.end());
            xcInput.insert(xcInput.end(), grad_b.begin(), grad_b.end());
            nUsed += grad_a.size() + grad_b.size();
        } else {
            xcInput.insert(xcInput.end(), grad_t.begin(), grad_t.end());
            nUsed += grad_t.size();
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

    int nOut = getOutputLength();
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
    int nInp = getInputLength();
    int nOut = getOutputLength();

#pragma omp parallel firstprivate(nInp, nOut)
    {
        int nNodes = mrcpp::get_func(xcInput, 0).getNEndNodes();
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MatrixXd inpData, outData;
            compressNodeData(n, nInp, xcInput, inpData);
            evaluateBlock(inpData, outData);
            expandNodeData(n, nOut, xcOutput, outData);
        }
    }
    for (int i = 0; i < nOut; i++) {
        mrcpp::get_func(xcOutput, i).mwTransform(mrcpp::BottomUp);
        mrcpp::get_func(xcOutput, i).calcSquareNorm();
    }

    auto n = mrcpp::get_n_nodes(xcOutput);
    auto m = mrcpp::get_size_nodes(xcOutput);
    auto t = timer.elapsed();
    mrcpp::print::tree(2, "XC evaluate xcfun", n, m, t);
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

    return xc_pot;
}

/** @brief Potential calculation for LDA functionals
 *
 * The potential conicides with the xcfun output, which is then
 * deep copied into the corresponding potential functions.
 */
void XCFunctional::calcPotentialLDA(FunctionTreeVector<3> &potentials) {
    int nPotentials = 1;
    int nStart = this->order; // PROBLEM: if I use a higher order than necessary this wil fail miserably!
    if (isSpinSeparated()) {
        nPotentials = this->order + 1;
        nStart = this->order * (this->order + 1) / 2;
    }
    for (int i = 0; i < nPotentials; i++) {
        FunctionTree<3> &out_i = mrcpp::get_func(xcOutput, nStart + i);
        auto *pot = new FunctionTree<3>(MRA);
        mrcpp::copy_grid(*pot, out_i);
        mrcpp::copy_func(*pot, out_i);
        potentials.push_back(std::make_tuple(1.0, pot));
    }
}

/** @brief Potential calculation for GGA functionals
 *
 */
void XCFunctional::calcPotentialGGA(FunctionTreeVector<3> &potentials) {
    switch (this->order) {
        case 1:
            calcGradientGGA(potentials);
            break;
        case 2:
            calcHessianGGA(potentials);
            break;
        default:
            NOT_IMPLEMENTED_ABORT;
    }
}

/** @brief Potential calculation for GGA functionals
 *
 * The potential functions are assembled from the xcfun output functions.
 * The method used depends on whether the functional is spin-separated
 * and whether explicit or gamma-type derivatives have been used in xcfun.
 */
void XCFunctional::calcGradientGGA(FunctionTreeVector<3> &potentials) {
    FunctionTree<3> *pot;
    if (isSpinSeparated()) {
        FunctionTree<3> &df_da = mrcpp::get_func(xcOutput, 1);
        FunctionTree<3> &df_db = mrcpp::get_func(xcOutput, 2);
        if (useGamma()) {
            FunctionTree<3> &df_dgaa = mrcpp::get_func(xcOutput, 3);
            FunctionTree<3> &df_dgab = mrcpp::get_func(xcOutput, 4);
            FunctionTree<3> &df_dgbb = mrcpp::get_func(xcOutput, 5);
            pot = calcGradientGGA(df_da, df_dgaa, df_dgab, grad_a, grad_b);
            potentials.push_back(std::make_tuple(1.0, pot));
            pot = calcGradientGGA(df_db, df_dgbb, df_dgab, grad_b, grad_a);
            potentials.push_back(std::make_tuple(1.0, pot));
        } else {
            FunctionTreeVector<3> df_dga;
            FunctionTreeVector<3> df_dgb;
            df_dga.push_back(xcOutput[3]);
            df_dga.push_back(xcOutput[4]);
            df_dga.push_back(xcOutput[5]);
            df_dgb.push_back(xcOutput[6]);
            df_dgb.push_back(xcOutput[7]);
            df_dgb.push_back(xcOutput[8]);
            pot = calcGradientGGA(df_da, df_dga);
            potentials.push_back(std::make_tuple(1.0, pot));
            pot = calcGradientGGA(df_db, df_dgb);
            potentials.push_back(std::make_tuple(1.0, pot));
        }
    } else {
        FunctionTree<3> &df_dt = mrcpp::get_func(xcOutput, 1);
        if (useGamma()) {
            FunctionTree<3> &df_dgamma = mrcpp::get_func(xcOutput, 2);
            pot = calcGradientGGA(df_dt, df_dgamma, grad_t);
            potentials.push_back(std::make_tuple(1.0, pot));
        } else {
            FunctionTreeVector<3> df_dgt;
            df_dgt.push_back(xcOutput[2]);
            df_dgt.push_back(xcOutput[3]);
            df_dgt.push_back(xcOutput[4]);
            pot = calcGradientGGA(df_dt, df_dgt);
            potentials.push_back(std::make_tuple(1.0, pot));
        }
    }
    pot = nullptr;
}

void XCFunctional::calcHessianGGA(FunctionTreeVector<3> &potentials) {
    if (isSpinSeparated()) NOT_IMPLEMENTED_ABORT;
    if (useGamma()) {
        calcHessianGGAgamma(potentials);
    } else {
        calcHessianGGAgrad(potentials);
    }
}

void XCFunctional::calcHessianGGAgamma(FunctionTreeVector<3> &potentials) {
    // 0 1    2    3      4       5
    // f dfdr dfdg df2dr2 df2drdg df2dg2

    for (int i = 2; i < 6; i++) {
        FunctionTree<3> &out_i = mrcpp::get_func(xcOutput, i);
        auto *pot = new FunctionTree<3>(MRA);
        mrcpp::copy_grid(*pot, out_i);
        mrcpp::copy_func(*pot, out_i);
        potentials.push_back(std::make_tuple(1.0, pot));
    }
}

void XCFunctional::calcHessianGGAgrad(FunctionTreeVector<3> &potentials) {
    // 0 1    2-4  5      6-8     9-14
    // f dfdr dfdg df2dr2 df2drdg df2dg2

    for (int i = 5; i < 15; i++) {
        FunctionTree<3> &out_i = mrcpp::get_func(xcOutput, i);
        auto *pot = new FunctionTree<3>(MRA);
        mrcpp::copy_grid(*pot, out_i);
        mrcpp::copy_func(*pot, out_i);
        potentials.push_back(std::make_tuple(1.0, pot));
    }
}

/** @brief XC potential calculation
 *
 * @param[in] df_drho Functional derivative wrt rho
 * @param[in] df_dgamma Functional_derivative wrt gamma
 * @param[in] grad_rho Gradient of rho
 *
 * Computes the XC potential for a non-spin separated functional and
 * gamma-type derivatives.
 */
FunctionTree<3> *XCFunctional::calcGradientGGA(FunctionTree<3> &df_drho,
                                               FunctionTree<3> &df_dgamma,
                                               FunctionTreeVector<3> grad_rho) {
    FunctionTreeVector<3> funcs;
    funcs.push_back(std::make_tuple(1.0, &df_drho));

    FunctionTree<3> *tmp = calcGradDotPotDensVec(df_dgamma, grad_rho);
    funcs.push_back(std::make_tuple(-2.0, tmp));

    auto *V = new FunctionTree<3>(MRA);
    mrcpp::build_grid(*V, funcs);
    mrcpp::add(-1.0, *V, funcs);
    delete tmp;
    return V;
}

/** @brief XC potential calculation
 *
 * @param[in] df_drhoa Functional derivative wrt rhoa
 * @param[in] df_dgaa  Functional_derivative wrt gamma_aa
 * @param[in] df_dgab  Functional_derivative wrt gamma_ab
 * @param[in] df_dgbb  Functional_derivative wrt gamma_bb
 * @param[in] grad_rhoa Gradient of rho_a
 * @param[in] grad_rhob Gradient of rho_b
 *
 * Computes the XC potential for a spin separated functional and
 * gamma-type derivatives. Can be used both for alpha and beta
 * potentials by permuting the spin parameter.
 */
FunctionTree<3> *XCFunctional::calcGradientGGA(FunctionTree<3> &df_drhoa,
                                               FunctionTree<3> &df_dgaa,
                                               FunctionTree<3> &df_dgab,
                                               FunctionTreeVector<3> grad_rhoa,
                                               FunctionTreeVector<3> grad_rhob) {
    FunctionTreeVector<3> funcs;
    funcs.push_back(std::make_tuple(1.0, &df_drhoa));

    FunctionTree<3> *tmp1 = calcGradDotPotDensVec(df_dgaa, grad_rhoa);
    funcs.push_back(std::make_tuple(-2.0, tmp1));

    FunctionTree<3> *tmp2 = calcGradDotPotDensVec(df_dgab, grad_rhob);
    funcs.push_back(std::make_tuple(-1.0, tmp2));

    auto *V = new FunctionTree<3>(MRA);
    mrcpp::build_grid(*V, funcs);
    mrcpp::add(-1.0, *V, funcs);
    delete tmp1;
    delete tmp2;
    return V;
}

/** @brief XC potential calculation
 *
 * @param[in] df_drho Functional derivative wrt rho
 * @param[in] df_dgr  Functional_derivative wrt grad_rho
 *
 * Computes the XC potential for explicit derivatives.
 */
FunctionTree<3> *XCFunctional::calcGradientGGA(FunctionTree<3> &df_drho, FunctionTreeVector<3> &df_dgr) {
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

/** @brief Helper function to compute divergence of a vector field times a function
 *
 * @param[in] V Function (derivative of the functional wrt gamma)
 * @param[in] rho Vector field (density gradient)
 *
 */
FunctionTree<3> *XCFunctional::calcGradDotPotDensVec(FunctionTree<3> &V, FunctionTreeVector<3> &rho) {
    FunctionTreeVector<3> vec;
    for (int d = 0; d < rho.size(); d++) {
        Timer timer;
        FunctionTree<3> &rho_d = mrcpp::get_func(rho, d);
        auto *Vrho = new FunctionTree<3>(MRA);
        mrcpp::copy_grid(*Vrho, rho_d);
        mrcpp::multiply(-1.0, *Vrho, 1.0, V, rho_d);
        vec.push_back(std::make_tuple(1.0, Vrho));
        mrcpp::print::tree(2, "Multiply", *Vrho, timer);
    }

    Timer timer;
    auto *result = new FunctionTree<3>(MRA);
    mrcpp::divergence(*result, *derivative, vec);
    mrcpp::clear(vec, true);
    mrcpp::print::tree(2, "Divergence", *result, timer);
    return result;
}

} // namespace mrdft
