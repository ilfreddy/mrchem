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

enum class DensityType { Total, Alpha, Beta };

class XCFunctional final {
public:
    XCFunctional(mrcpp::MultiResolutionAnalysis<3> &mra, bool spin);
    XCFunctional(const XCFunctional &func) = delete;
    XCFunctional &operator=(const XCFunctional &func) = delete;
    ~XCFunctional();

    bool hasDensity() const;
    bool checkDensity(FunctionTreeVector<3> density) const;
    mrcpp::FunctionTreeVector<3> &getDensity(DensityType type);

    int getNNodes() const;
    int getNPoints() const;

    void buildGrid(double Z, const mrcpp::Coord<3> &R);
    void copyGrid(FunctionTreeVector<3> densities);
    void clearGrid();

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

protected:
    int order;
    int nDensities;
    unsigned int mode;
    const bool spin_separated;                  ///< Spin polarization
    const mrcpp::MultiResolutionAnalysis<3> MRA;///< Computational domain

    bool use_gamma; ///< Whether gamma-type or explicit derivatives are used
    double cutoff;  ///< Below the cutoff value, the density will be considered zero

    xc_functional functional;                 ///< The functional in the XCFun library (struct from xcfun library)
    mrcpp::DerivativeOperator<3> *derivative; ///< Derivative operator

    mrcpp::FunctionTreeVector<3> rho_a;  ///< Alpha densities
    mrcpp::FunctionTreeVector<3> rho_b;  ///< Beta densities
    mrcpp::FunctionTreeVecroe<3> rho_t;  ///< Total densities
    mrcpp::FunctionTreeVector<3> grad_a; ///< Gradient of the alpha densities
    mrcpp::FunctionTreeVector<3> grad_b; ///< Gradient of the beta  densities
    mrcpp::FunctionTreeVector<3> grad_t; ///< Gradient of the total densities
    mrcpp::FunctionTreeVector<3> gamma;  ///< Gamma function(s)

    mrcpp::FunctionTreeVector<3> xcInput;  ///< Bookkeeping array to feed XCFun
    mrcpp::FunctionTreeVector<3> xcOutput; ///< Bookkeeping array returned by XCFun

    void clearGrid();
    void clearGrid(FunctionTreeVector<3> densities);
    int getInputLength() const { return xc_input_length(this->functional); }
    int getOutputLength() const { return xc_output_length(this->functional); }

    void setup_partial();
    void setup_partial(FunctionTree<3> &rho_a, FunctionTree<3> &rho_b);
    void setup_partial(FunctionTree<3> &rho_t);
    void setup_contracted();
    void setupXCInput();
    void setupXCOutput();
    int setupXCInputDensity();
    int setupXCInputGradient();

    void evaluateBlock(Eigen::MatrixXd &inp, Eigen::MatrixXd &out) const;
    void compressNodeData(int n, int nFuncs, mrcpp::FunctionTreeVector<3> trees, Eigen::MatrixXd &data);
    void expandNodeData(int n, int nFuncs, mrcpp::FunctionTreeVector<3> trees, Eigen::MatrixXd &data);

    void calcPotentialLDA(mrcpp::FunctionTreeVector<3> &potentials);
    void calcPotentialGGA(mrcpp::FunctionTreeVector<3> &potentials);
    void calcGradientGGA(mrcpp::FunctionTreeVector<3> &potentials);
    void calcHessianGGA(mrcpp::FunctionTreeVector<3> &potentials);
    void calcHessianGGAgamma(mrcpp::FunctionTreeVector<3> &potentials);
    void calcHessianGGAgrad(mrcpp::FunctionTreeVector<3> &potentials);

    mrcpp::FunctionTree<3> *calcGradientGGA(mrcpp::FunctionTree<3> &df_drho, mrcpp::FunctionTreeVector<3> &df_dgr);
    mrcpp::FunctionTree<3> *calcGradientGGA(mrcpp::FunctionTree<3> &df_drho,
                                            mrcpp::FunctionTree<3> &df_dgamma,
                                            mrcpp::FunctionTreeVector<3> grad_rho);
    mrcpp::FunctionTree<3> *calcGradientGGA(mrcpp::FunctionTree<3> &df_drhoa,
                                            mrcpp::FunctionTree<3> &df_dgaa,
                                            mrcpp::FunctionTree<3> &df_dgab,
                                            mrcpp::FunctionTreeVector<3> grad_rhoa,
                                            mrcpp::FunctionTreeVector<3> grad_rhob);

    mrcpp::FunctionTree<3> *calcGradDotPotDensVec(mrcpp::FunctionTree<3> &V, mrcpp::FunctionTreeVector<3> &rho);
};

} // namespace mrdft
