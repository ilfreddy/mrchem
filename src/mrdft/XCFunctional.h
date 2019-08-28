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

#pragma once

#include <Eigen/Core>

#include "MRCPP/MWFunctions"
#include "XCFun/xcfun.h"
#include "mrenum.h"

/**
 *  @class XCFunctional
 *  @brief Compute XC functional with XCFun
 *
 *  Interface class for the XCFun library
 *
 * Depending on the mode chosen, xcfun needs either the gamma
 * functions or the explicit gradients. The first mode is possibly
 * more efficient (fewer functions to compute/handle), whereas the
 * other is simpler to implement. We keep both options open and
 * compute the gradient invariants if and when necessary.
 *
 * The XCFunctional keeps track of the density grid, whose size is
 * initially defined through the interface function buildGrid().
 * The grid is kept _fixed_ for all internal calculations
 * within the module.
 *
 * Typical usage within one SCF cycle:
 *
 * 1) getDensity() and compute density on given grid
 * 2) setup()
 * 3) evaluate()
 * 4) calcEnergy()
 * 5) calcPotential()
 * 7) clear()
 *
 */


namespace mrdft {

class XCFunctional final {
public:
    XCFunctional(mrcpp::MultiResolutionAnalysis<3> &mra, bool spin);
    XCFunctional(const XCFunctional &func) = delete;
    XCFunctional &operator=(const XCFunctional &func) = delete;
    ~XCFunctional();

    bool hasDensity(int n_dens_a = 1
                    ) const;
    bool checkDensity(mrcpp::FunctionTreeVector<3> density, int n_dens = 1) const;
    mrcpp::FunctionTreeVector<3> &getDensityVector(DENSITY::DensityType spin);
    mrcpp::FunctionTree<3> &getDensity(DENSITY::DensityType type, int index = 0);
    void setDensity(mrcpp::FunctionTree<3> &density, DENSITY::DensityType spin, int index = 0);

    int getNNodes() const;
    int getNPoints() const;
    void allocateDensities();

    void buildGrid(double Z, const mrcpp::Coord<3> &R);
    void copyGrid(mrcpp::FunctionTreeVector<3> densities);

    void setNDensities(int n) { this->nDensities = n; }
    int getNDensities() { return this->nDensities; }
    void setDensityCutoff(double thrs) { this->cutoff = thrs; }
    void setFunctional(const std::string &name, double coef = 1.0);

    void setUseGamma(bool use) { this->use_gamma = use; }
    bool useGamma() const { return this->use_gamma; }

    int getOrder() const { return this->order; }
    bool isLDA() const { return (not(isGGA() or isMetaGGA())); }
    bool isGGA() const { return (xc_is_gga(this->functional)); }
    bool isMetaGGA() const { return (xc_is_metagga(this->functional)); }
    bool isHybrid() const { return (std::abs(this->amountEXX()) > mrcpp::MachineZero); }
    bool isSpinSeparated() const { return this->spin_separated; }

    double amountEXX() const {
        double exx = 0.0;
        xc_get(this->functional, "exx", &exx);
        return exx;
    }

    void evalSetup(int order);
    void setup();
    void clear();

    void evaluate();
    double calcEnergy();
    mrcpp::FunctionTreeVector<3> calcPotential();

private:
    int order;
    int nDensities;
    const bool spin_separated;                  ///< Spin polarization
    const mrcpp::MultiResolutionAnalysis<3> MRA;///< Computational domain

    bool use_gamma; ///< Whether gamma-type or explicit derivatives are used
    double cutoff;  ///< Below the cutoff value, the density will be considered zero

    xc_functional functional;                 ///< The functional in the XCFun library (struct from xcfun library)
    mrcpp::DerivativeOperator<3> *derivative; ///< Derivative operator

    mrcpp::FunctionTreeVector<3> rho_a;  ///< Alpha densities
    mrcpp::FunctionTreeVector<3> rho_b;  ///< Beta densities
    mrcpp::FunctionTreeVector<3> rho_t;  ///< Total densities
    mrcpp::FunctionTreeVector<3> grad_a; ///< Gradient of the alpha densities
    mrcpp::FunctionTreeVector<3> grad_b; ///< Gradient of the beta  densities
    mrcpp::FunctionTreeVector<3> grad_t; ///< Gradient of the total densities
    mrcpp::FunctionTreeVector<3> gamma;  ///< Gamma function(s)

    mrcpp::FunctionTreeVector<3> xcInput;   ///< Bookkeeping array to feed XCFun
    mrcpp::FunctionTreeVector<3> xcOutput;  ///< Bookkeeping array returned by XCFun
    mrcpp::FunctionTreeVector<3> xcDensity; ///< Bookkeeping array with density variables

    void clearGrid();
    void clearGrid(mrcpp::FunctionTreeVector<3> densities);
    int getInputLength() const { return xc_input_length(this->functional); }
    int getOutputLength() const { return xc_output_length(this->functional); }
    int getContractedLength() const; 
    int getDensityLength() const;
    int getNodeLength();
    
    void setupGradient();
        
    void setupXCInput();
    void setupXCOutput();
    int setupXCInputDensity();
    int setupXCInputGradient();
    void setupXCDensityVariables();

    void evaluateBlock(Eigen::MatrixXd &inp, Eigen::MatrixXd &out) const;
    void compressNodeData(int n, int nFuncs, mrcpp::FunctionTreeVector<3> trees, Eigen::MatrixXd &data);
    void expandNodeData(int n, int nFuncs, mrcpp::FunctionTreeVector<3> trees, Eigen::MatrixXd &data);
    void contractNodeData(int n, int n_coefs, Eigen::MatrixXd &out_data, Eigen::MatrixXd &con_data);

    void calcPotentialLDA(mrcpp::FunctionTreeVector<3> &potentials);
    void calcPotentialGGA(mrcpp::FunctionTreeVector<3> &potentials);

    mrcpp::FunctionTree<3> *calcPotentialGGA(mrcpp::FunctionTree<3> &df_drho, mrcpp::FunctionTreeVector<3> &df_dgr);

    mrcpp::FunctionTree<3> *calcGradDotPotDensVec(mrcpp::FunctionTree<3> &V, mrcpp::FunctionTreeVector<3> &rho);
};

} // namespace mrdft
