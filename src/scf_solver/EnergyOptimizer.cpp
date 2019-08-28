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

#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "EnergyOptimizer.h"
#include "HelmholtzVector.h"

#include "chemistry/Molecule.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/CoulombOperator.h"
#include "qmoperators/two_electron/ExchangeOperator.h"
#include "qmoperators/two_electron/FockOperator.h"
#include "qmoperators/two_electron/XCOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Run energy optimization
 *
 * Optimize orbitals until convergence thresholds are met. This algorithm does NOT
 * use any kinetic energy operator, but computes the Fock matrix as a potential
 * update from the previous iteration. No KAIN accelerator and no diagonalization
 * or localization within the SCF iterations. Main points of the algorithm:
 *
 * Pre SCF: diagonalize/localize orbitals
 *
 * 1) Setup Fock operator
 * 2) Compute current SCF energy
 * 3) Apply Helmholtz operator on all orbitals
 * 4) Compute orbital updates
 * 5) Compute errors and check for convergence
 * 6) Compute Fock matrix update
 * 7) Orthonormalize orbitals (Löwdin)
 *
 * Post SCF: diagonalize/localize orbitals
 */
bool EnergyOptimizer::optimize(Molecule &mol, FockOperator &F_n) {
    printParameters("Optimize energy (kinetic free)");

    SCFEnergy &E_n = mol.getSCFEnergy();
    OrbitalVector &Phi_n = mol.getOrbitals();
    ComplexMatrix &F_mat_n = mol.getFockMatrix();

    DoubleVector errors = DoubleVector::Ones(Phi_n.size());
    double err_o = errors.maxCoeff();
    double err_t = errors.norm();
    double err_p = 1.0;

    this->error.push_back(err_t);
    this->energy.push_back(E_n);
    this->property.push_back(E_n.getTotalEnergy());

    auto plevel = Printer::getPrintLevel();
    if (plevel < 1) {
        printConvergenceHeader("Total energy");
        printConvergenceRow(0);
    }

    int nIter = 0;
    bool converged = false;
    while (nIter++ < this->maxIter or this->maxIter < 0) {
        std::stringstream o_header;
        o_header << "SCF cycle " << nIter;
        mrcpp::print::header(1, o_header.str(), 0, '#');
        mrcpp::print::separator(2, ' ', 1);

        // Initialize SCF cycle
        Timer t_lap;
        double orb_prec = adjustPrecision(err_o);
        if (nIter < 2) F_n.setup(orb_prec);

        // Init Helmholtz operator
        HelmholtzVector H(orb_prec, F_mat_n.real().diagonal());

        // Setup argument
        Timer t_arg;
        mrcpp::print::header(2, "Computing Helmholtz argument");
        ComplexMatrix L_mat = H.getLambdaMatrix();
        OrbitalVector Psi = orbital::rotate(L_mat - F_mat_n, Phi_n);
        mrcpp::print::time(2, "Rotating orbitals", t_arg);
        mrcpp::print::footer(2, t_arg, 2);
        if (plevel == 1) mrcpp::print::time(1, "Computing Helmholtz argument", t_arg);

        // Apply Helmholtz operator
        OrbitalVector Phi_np1 = H.apply(F_n.potential(), Phi_n, Psi);
        Psi.clear();

        // Compute orbital updates
        OrbitalVector dPhi_n = orbital::add(1.0, Phi_np1, -1.0, Phi_n);

        // Compute errors
        errors = orbital::get_norms(dPhi_n);
        err_o = errors.maxCoeff();
        err_t = errors.norm();
        err_p = calcPropertyError();
        converged = checkConvergence(err_o, err_p);

        // Compute Fock matrix
        ComplexMatrix dS_1 = orbital::calc_overlap_matrix(dPhi_n, Phi_n);
        ComplexMatrix dF_mat_1 = dS_1 * F_mat_n;

        ComplexMatrix dS_2 = orbital::calc_overlap_matrix(Phi_np1, dPhi_n);
        ComplexMatrix L_mat_n = H.getLambdaMatrix();
        ComplexMatrix dF_mat_2 = dS_2 * L_mat_n;

        ComplexMatrix dF_mat_3 = calcFockMatrixUpdate(orb_prec, F_n, Phi_np1, Phi_n, dPhi_n);
        dPhi_n.clear();
        F_n.clear();

        // Symmetrizing
        ComplexMatrix F_mat_np1 = F_mat_n + dF_mat_1 + dF_mat_2 + dF_mat_3;
        ComplexMatrix F_mat_sym = 0.5 * (F_mat_np1 + F_mat_np1.transpose());

        // Rotate orbitals
        ComplexMatrix U_mat = orbital::calc_lowdin_matrix(Phi_np1);
        Phi_n = orbital::rotate(U_mat, Phi_np1, orb_prec);
        F_mat_n = U_mat * F_mat_sym * U_mat.adjoint();
        Phi_np1.clear();

        // Compute energy
        F_n.setup(orb_prec);
        E_n = F_n.trace(Phi_n, F_mat_n);

        if (converged) {
            if (localize) {
                orbital::localize(orb_prec, Phi_n, F_mat_n);
            } else {
                orbital::diagonalize(orb_prec, Phi_n, F_mat_n);
            }
        }

        // Collect convergence data
        this->error.push_back(err_t);
        this->energy.push_back(E_n);
        this->property.push_back(E_n.getTotalEnergy());

        if (plevel < 1) printConvergenceRow(nIter);
        printOrbitals(F_mat_n.real().diagonal(), errors, Phi_n, 0);
        printProperty();
        mrcpp::print::separator(2, ' ', 1);
        mrcpp::print::footer(1, t_lap, 2, '#');
        mrcpp::print::separator(2, ' ', 2);

        if (converged) break;
    }
    F_n.clear();

    printConvergence(converged, "Total energy");
    reset();

    return converged;
}

/** @brief Compute Fock matrix update
 *
 * The Fock matrix is computed as a direct update from the previous iteration.
 * This update is exact, provided that the orbital relation between the orbitals
 * at iterations n and n+1 is exactly the application of the Helmholtz operator.
 *
 * \Delta F = \Delta S_1 F + \Delta S_2 \Lambda + \Delta F_{pot}
 *
 * (\Delta S_1^n)_{ij} = <\Delta \phi_i^{n}   | \phi_j^n>
 * (\Delta S_2^n)_{ij} = <\Delta \phi_i^{n+1} | \Delta \phi_j^n>
 *
 * (\Delta F_^n{pot})_{ij} = <\phi_i^{n+1} | V^n | \Delta \phi_j^n>
 *                         + <\phi_i^{n+1} | \Delta V^n | \phi_j^n>
 *
 */
ComplexMatrix EnergyOptimizer::calcFockMatrixUpdate(double prec,
                                                    FockOperator &F_n,
                                                    OrbitalVector &Phi_np1,
                                                    OrbitalVector &Phi_n,
                                                    OrbitalVector &dPhi_n) {
    Timer t_tot, t_lap;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing Fock matrix update");

    auto v_n = F_n.getNuclearOperator();
    auto j_n = F_n.getCoulombOperator();
    auto k_n = F_n.getExchangeOperator();
    auto xc_n = F_n.getXCOperator();

    auto exx = 1.0;
    if (xc_n != nullptr) exx = xc_n->getFunctional()->amountEXX();

    ComplexMatrix dV_mat_n;
    { // Nuclear potential matrix is computed explicitly
        t_lap.start();
        dV_mat_n = (*v_n)(Phi_np1, dPhi_n);
        mrcpp::print::time(2, "Nuclear potential matrix", t_lap);
        mrcpp::print::separator(2, '-');
    }

    ComplexMatrix W_mat_n = ComplexMatrix::Zero(Phi_n.size(), Phi_n.size());
    { // Computing two-electron part of Fock matrix
        t_lap.start();
        if (j_n != nullptr) W_mat_n += (*j_n)(Phi_np1, Phi_n);
        if (k_n != nullptr) W_mat_n -= exx * (*k_n)(Phi_np1, Phi_n);
        if (xc_n != nullptr) W_mat_n += (*xc_n)(Phi_np1, Phi_n);
        mrcpp::print::time(2, "Fock matrix n", t_lap);
        mrcpp::print::separator(2, '-');
    }

    { // The n+1 Fock operator needs orthonormalized orbitals
        t_lap.start();
        ComplexMatrix U_mat = orbital::calc_lowdin_matrix(Phi_np1);
        Phi_np1 = orbital::rotate(U_mat, Phi_np1, prec);
        mrcpp::print::time(1, "Orthonormalize n+1", t_lap);
        mrcpp::print::separator(2, '-');
    }

    std::shared_ptr<CoulombOperator> j_np1{nullptr};
    std::shared_ptr<ExchangeOperator> k_np1{nullptr};
    std::shared_ptr<XCOperator> xc_np1{nullptr};
    { // Setup Fock operator n+1
        t_lap.start();
        auto Phi_p = std::make_shared<OrbitalVector>();
        *Phi_p = Phi_np1;

        if (j_n != nullptr) {
            j_np1 = std::make_shared<CoulombOperator>(j_n->getPoisson(), Phi_p);
            j_np1->setup(prec);
        }
        if (k_n != nullptr) {
            // Do not setup internal exchange, it must be applied on the fly anyway
            k_np1 = std::make_shared<ExchangeOperator>(k_n->getPoisson(), Phi_p);
            k_np1->setup(prec);
        }
        if (xc_n != nullptr) {
            xc_np1 = std::make_shared<XCOperator>(xc_n->getFunctional(), Phi_p);
            xc_np1->setup(prec);
        }
        mrcpp::print::time(2, "Fock operator n+1", t_lap);
        mrcpp::print::separator(2, '-');
    }

    ComplexMatrix W_mat_np1;
    { // Computing potential matrix excluding nuclear part
        t_lap.start();
        ComplexMatrix W_mat_1 = ComplexMatrix::Zero(Phi_n.size(), Phi_n.size());
        if (j_np1 != nullptr) W_mat_1 += (*j_np1)(Phi_n, Phi_n);
        if (k_np1 != nullptr) W_mat_1 -= exx * (*k_np1)(Phi_n, Phi_n);
        if (xc_np1 != nullptr) W_mat_1 += (*xc_np1)(Phi_n, Phi_n);

        ComplexMatrix W_mat_2 = ComplexMatrix::Zero(Phi_n.size(), Phi_n.size());
        if (j_np1 != nullptr) W_mat_2 += (*j_np1)(Phi_n, dPhi_n);
        if (k_np1 != nullptr) W_mat_2 -= exx * (*k_np1)(Phi_n, dPhi_n);
        if (xc_np1 != nullptr) W_mat_2 += (*xc_np1)(Phi_n, dPhi_n);

        W_mat_np1 = W_mat_1 + W_mat_2 + W_mat_2.transpose();
        mrcpp::print::time(2, "Fock matrix n+1", t_lap);
    }
    if (j_np1 != nullptr) j_np1->clear();
    if (k_np1 != nullptr) k_np1->clear();
    if (xc_np1 != nullptr) xc_np1->clear();

    // Re-computing non-orthogonal phi_np1
    Phi_np1 = orbital::add(1.0, Phi_n, 1.0, dPhi_n);

    // Adding up the pieces
    ComplexMatrix dF_mat_n = dV_mat_n + W_mat_np1 - W_mat_n;

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Calculating Fock matrix update", t_tot);

    return dF_mat_n;
}

} // namespace mrchem
