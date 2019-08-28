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

#include "HelmholtzVector.h"
#include "KAIN.h"
#include "OrbitalOptimizer.h"

#include "chemistry/Molecule.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

#include "qmoperators/two_electron/FockOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Run orbital optimization
 *
 * @param F: FockOperator defining the SCF equations
 * @param Phi_n: Orbitals to optimize
 * @param F_mat: Fock matrix
 *
 * Optimize orbitals until convergence thresholds are met. This algorithm computes
 * the Fock matrix explicitly using the kinetic energy operator, and uses a KAIN
 * accelerator to improve convergence. Diagonalization or localization may be performed
 * during the SCF iterations. Main points of the algorithm:
 *
 * Pre SCF: setup Fock operator and compute Fock matrix
 *
 *  1) Diagonalize/localize orbitals
 *  2) Compute current SCF energy
 *  3) Apply Helmholtz operator on all orbitals
 *  4) Orthonormalize orbitals (Löwdin)
 *  5) Compute orbital updates
 *  6) Compute KAIN update
 *  7) Compute errors and check for convergence
 *  8) Add orbital updates
 *  9) Orthonormalize orbitals (Löwdin)
 * 10) Setup Fock operator
 * 11) Compute Fock matrix
 *
 * Post SCF: diagonalize/localize orbitals
 */
bool OrbitalOptimizer::optimize(Molecule &mol, FockOperator &F) {
    printParameters("Optimize molecular orbitals");

    KAIN kain(this->history);
    SCFEnergy &E_n = mol.getSCFEnergy();
    OrbitalVector &Phi_n = mol.getOrbitals();
    ComplexMatrix &F_mat = mol.getFockMatrix();

    DoubleVector errors = DoubleVector::Ones(Phi_n.size());
    double err_o = errors.maxCoeff();
    double err_t = errors.norm();

    this->error.push_back(err_t);
    this->energy.push_back(E_n);
    this->property.push_back(E_n.getTotalEnergy());

    auto plevel = Printer::getPrintLevel();
    if (plevel < 1) {
        printConvergenceHeader();
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
        double helm_prec = getHelmholtzPrec();
        if (nIter < 2) F.setup(orb_prec);

        // Apply Helmholtz operator
        HelmholtzVector H(helm_prec, F_mat.real().diagonal());
        OrbitalVector Psi = H.rotate(F_mat, Phi_n);
        OrbitalVector Phi_np1 = H.apply(F.potential(), Phi_n, Psi);
        Psi.clear();
        F.clear();

        orbital::orthonormalize(orb_prec, Phi_np1, F_mat);

        // Compute orbital updates
        OrbitalVector dPhi_n = orbital::add(1.0, Phi_np1, -1.0, Phi_n);
        Phi_np1.clear();

        // Employ KAIN accelerator
        kain.accelerate(orb_prec, Phi_n, dPhi_n);

        // Compute errors
        errors = orbital::get_norms(dPhi_n);
        err_o = errors.maxCoeff();
        err_t = errors.norm();

        // Update orbitals
        Phi_n = orbital::add(1.0, Phi_n, 1.0, dPhi_n);
        dPhi_n.clear();

        orbital::orthonormalize(orb_prec, Phi_n, F_mat);

        // Compute Fock matrix and energy
        F.setup(orb_prec);
        F_mat = F(Phi_n, Phi_n);
        E_n = F.trace(Phi_n, F_mat);

        // Collect convergence data
        this->error.push_back(err_t);
        this->energy.push_back(E_n);
        this->property.push_back(E_n.getTotalEnergy());
        auto err_p = calcPropertyError();
        converged = checkConvergence(err_o, err_p);

        // Rotate orbitals
        if (needLocalization(nIter, converged)) {
            ComplexMatrix U_mat = orbital::localize(orb_prec, Phi_n, F_mat);
            F.rotate(U_mat);
            kain.clear();
        } else if (needDiagonalization(nIter, converged)) {
            ComplexMatrix U_mat = orbital::diagonalize(orb_prec, Phi_n, F_mat);
            F.rotate(U_mat);
            kain.clear();
        }

        // Finalize SCF cycle
        if (plevel < 1) printConvergenceRow(nIter);
        mrcpp::print::separator(1, '-');
        printOrbitals(F_mat.real().diagonal(), errors, Phi_n, 0);
        printProperty();
        printMemory();
        mrcpp::print::footer(1, t_lap, 2, '#');
        mrcpp::print::separator(2, ' ', 2);

        if (converged) break;
    }

    F.clear();

    printConvergence(converged);
    reset();

    return converged;
}

} // namespace mrchem
