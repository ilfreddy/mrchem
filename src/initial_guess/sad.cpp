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

#include <Eigen/Eigenvalues>

#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "core.h"
#include "gto.h"
#include "parallel.h"
#include "sad.h"
#include "utils/math_utils.h"

#include "chemistry/Molecule.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/OrbitalIterator.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"

#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/CoulombOperator.h"
#include "qmoperators/two_electron/XCOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace initial_guess {
namespace sad {

ComplexMatrix diagonalize_fock(KineticOperator &T, RankZeroTensorOperator &V, OrbitalVector &Phi, int spin);
OrbitalVector rotate_orbitals(double prec, ComplexMatrix &U, OrbitalVector &Phi, int N, int spin);
void project_atomic_densities(double prec, Density &rho_tot, const Molecule &mol);

} // namespace sad
} // namespace initial_guess

OrbitalVector initial_guess::sad::setup(double prec, const Molecule &mol, bool restricted, int zeta) {
    std::stringstream o_prec, o_zeta;
    o_prec << std::setprecision(5) << std::scientific << prec;
    o_zeta << zeta;
    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation ", "Diagonalize Hamiltonian matrix");
    print_utils::text(0, "Precision   ", o_prec.str());
    print_utils::text(0, "Restricted  ", (restricted) ? "True" : "False");
    print_utils::text(0, "Hamiltonian ", "Superposition of Atomic Densities (SAD)");
    print_utils::text(0, "Functional  ", "LDA (SVWN5)");
    print_utils::text(0, "AO basis    ", "Hydrogenic orbitals");
    print_utils::text(0, "Zeta quality", o_zeta.str());
    mrcpp::print::separator(0, '~', 2);

    Timer t_tot, t_lap;
    auto plevel = Printer::getPrintLevel();
    if (plevel == 1) mrcpp::print::header(1, "SAD Initial Guess");

    // Figure out number of occupied orbitals
    int mult = mol.getMultiplicity(); // multiplicity
    int Ne = mol.getNElectrons();     // total electrons
    int Nd = Ne - (mult - 1);         // doubly occupied electrons
    if (Nd % 2 != 0) MSG_ABORT("Invalid multiplicity");
    int Na = Nd / 2 + (mult - 1); // alpha orbitals
    int Nb = Nd / 2;              // beta orbitals

    // Make Fock operator contributions
    t_lap.start();
    auto P_p = std::make_shared<mrcpp::PoissonOperator>(*MRA, prec);
    auto D_p = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);
    auto xcfun_p = std::make_shared<mrdft::XCFunctional>(*MRA, not(restricted));
    xcfun_p->setFunctional("SLATERX");
    xcfun_p->setFunctional("VWN5C");
    xcfun_p->evalSetup(MRDFT::Gradient);
    xcfun.setNDensities(1);
    xcfun.allocateDensities();
    KineticOperator T(D_p);
    NuclearOperator V_nuc(mol.getNuclei(), prec);
    CoulombOperator J(P_p);
    XCOperator XC(xcfun_p);
    RankZeroTensorOperator V = V_nuc + J + XC;
    if (plevel == 1) mrcpp::print::time(1, "Projecting nuclear potential", t_lap);

    // Compute Coulomb density
    t_lap.start();
    Density &rho_j = J.getDensity();
    initial_guess::sad::project_atomic_densities(prec, rho_j, mol);

    // Compute XC density
    if (restricted) {
        mrcpp::FunctionTree<3> &rho_xc = XC.getDensity(DENSITY::DensityType::Total);
        mrcpp::copy_grid(rho_xc, rho_j.real());
        mrcpp::copy_func(rho_xc, rho_j.real());
    } else {
        mrcpp::FunctionTree<3> &rho_a = XC.getDensity(DENSITY::DensityType::Alpha);
        mrcpp::FunctionTree<3> &rho_b = XC.getDensity(DENSITY::DensityType::Beta);
        mrcpp::add(prec, rho_a, 1.0, rho_j.real(), -1.0 * Nb / Ne, rho_j.real());
        mrcpp::add(prec, rho_b, 1.0, rho_j.real(), -1.0 * Na / Ne, rho_j.real());

        // Extend to union grid
        int nNodes = 1;
        while (nNodes > 0) {
            int nAlpha = mrcpp::refine_grid(rho_a, rho_b);
            int nBeta = mrcpp::refine_grid(rho_b, rho_a);
            nNodes = nAlpha + nBeta;
        }
    }
    if (plevel == 1) mrcpp::print::time(1, "Projecting GTO density", t_lap);

    // Project AO basis of hydrogen functions
    t_lap.start();
    OrbitalVector Phi = initial_guess::core::project_ao(prec, mol.getNuclei(), SPIN::Paired, zeta);
    if (plevel == 1) mrcpp::print::time(1, "Projecting Hydrogen AOs", t_lap);

    mrcpp::print::header(2, "Building Fock operator");
    t_lap.start();
    T.setup(prec);
    V.setup(prec);
    mrcpp::print::footer(2, t_lap, 2);
    if (plevel == 1) mrcpp::print::time(1, "Building Fock operator", t_lap);

    // Compute Fock matrix
    mrcpp::print::header(2, "Diagonalizing Fock matrix");
    t_lap.start();
    OrbitalVector Psi;
    if (restricted) {
        if (mult != 1) MSG_ABORT("Restricted open-shell not available");
        int Np = Nd / 2; // paired orbitals
        ComplexMatrix U = initial_guess::sad::diagonalize_fock(T, V, Phi, SPIN::Paired);
        Psi = initial_guess::sad::rotate_orbitals(prec, U, Phi, Np, SPIN::Paired);
    } else {
        int Na = Nd / 2 + (mult - 1); // alpha orbitals
        int Nb = Nd / 2;              // beta orbitals
        ComplexMatrix U_a = initial_guess::sad::diagonalize_fock(T, V, Phi, SPIN::Alpha);
        OrbitalVector Psi_a = initial_guess::sad::rotate_orbitals(prec, U_a, Phi, Na, SPIN::Alpha);
        mrcpp::print::separator(2, '-');
        ComplexMatrix U_b = initial_guess::sad::diagonalize_fock(T, V, Phi, SPIN::Beta);
        OrbitalVector Psi_b = initial_guess::sad::rotate_orbitals(prec, U_b, Phi, Nb, SPIN::Beta);

        Psi = orbital::adjoin(Psi_a, Psi_b);
    }
    T.clear();
    V.clear();
    mrcpp::print::footer(2, t_lap, 2);
    if (plevel == 1) mrcpp::print::footer(1, t_tot, 2);

    return Psi;
}

ComplexMatrix initial_guess::sad::diagonalize_fock(KineticOperator &T,
                                                   RankZeroTensorOperator &V,
                                                   OrbitalVector &Phi,
                                                   int spin) {
    Timer t1;
    for (auto &i : Phi) i.setSpin(spin);
    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);
    ComplexMatrix f_tilde = T(Phi, Phi) + V(Phi, Phi);
    ComplexMatrix f = S_m12.adjoint() * f_tilde * S_m12;
    mrcpp::print::time(1, "Computing Fock matrix", t1);

    Timer t2;
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(f.cols());
    es.compute(f);
    ComplexMatrix ei_vec = es.eigenvectors();
    ComplexMatrix U = ei_vec.transpose() * S_m12;
    mrcpp::print::time(1, "Diagonalizing Fock matrix", t2);

    return U;
}

OrbitalVector initial_guess::sad::rotate_orbitals(double prec, ComplexMatrix &U, OrbitalVector &Phi, int N, int spin) {
    Timer t_tot;
    OrbitalVector Psi;
    for (int i = 0; i < N; i++) Psi.push_back(Orbital(spin));
    mpi::distribute(Psi);

    OrbitalIterator iter(Phi);
    while (iter.next()) {
        for (int i = 0; i < Psi.size(); i++) {
            if (not mpi::my_orb(Psi[i])) continue;
            QMFunctionVector func_vec;
            ComplexVector coef_vec(iter.get_size());
            for (int j = 0; j < iter.get_size(); j++) {
                int idx_j = iter.idx(j);
                Orbital &recv_j = iter.orbital(j);
                coef_vec[j] = U(i, idx_j);
                func_vec.push_back(recv_j);
            }
            Orbital tmp_i = Psi[i].paramCopy();
            qmfunction::linear_combination(tmp_i, coef_vec, func_vec, prec);
            Psi[i].add(1.0, tmp_i); // In place addition
            Psi[i].crop(prec);
        }
    }
    mrcpp::print::time(1, "Rotating orbitals", t_tot);
    return Psi;
}

void initial_guess::sad::project_atomic_densities(double prec, Density &rho_tot, const Molecule &mol) {
    auto pprec = Printer::getPrecision();
    auto w0 = Printer::getWidth() - 1;
    auto w1 = 5;
    auto w2 = 8;
    auto w3 = w0 / 3;
    auto w4 = w0 - (w1 + w2 + 2 * w3);

    std::stringstream o_head;
    o_head << std::setw(w1) << "N";
    o_head << std::setw(w2) << "Atom";
    o_head << std::setw(w4) << " ";
    o_head << std::setw(w3) << "Nuclear charge";
    o_head << std::setw(w3) << "Electron charge";

    mrcpp::print::header(2, "Projecting GTO density");
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    auto crop_prec = (mpi::numerically_exact) ? -1.0 : prec;
    std::string sad_path = SAD_BASIS_DIR;

    Timer t_tot;
    Density rho_loc(false);
    rho_loc.alloc(NUMBER::Real);
    rho_loc.real().setZero();

    auto tot_nuc = 0.0;
    auto tot_rho = 0.0;

    Timer t_loc;
    const Nuclei &nucs = mol.getNuclei();
    for (int k = 0; k < nucs.size(); k++) {
        if (mpi::orb_rank != k % mpi::orb_size) continue;

        const std::string &sym = nucs[k].getElement().getSymbol();
        std::stringstream o_bas, o_dens;
        o_bas << sad_path << sym << ".bas";
        o_dens << sad_path << sym << ".dens";

        Density rho_k = initial_guess::gto::project_density(prec, nucs[k], o_bas.str(), o_dens.str());
        rho_loc.add(1.0, rho_k);
        rho_loc.crop(crop_prec);

        auto nuc_charge = nucs[k].getCharge();
        auto rho_charge = rho_k.integrate().real();
        tot_nuc += nuc_charge;
        tot_rho += rho_charge;

        std::stringstream o_row;
        o_row << std::setw(w1) << k;
        o_row << std::setw(w2) << sym;
        o_row << std::setw(w4) << " ";
        o_row << std::setw(w3) << std::setprecision(2 * pprec) << std::fixed << nuc_charge;
        o_row << std::setw(w3) << std::setprecision(2 * pprec) << std::fixed << rho_charge;
        println(2, o_row.str());
    }
    t_loc.stop();

    Timer t_com;
    density::allreduce_density(prec, rho_tot, rho_loc);
    t_com.stop();

    std::stringstream o_row;
    o_row << " Total charge";
    o_row << std::string(w1 + w2 + w4 - 13, ' ');
    o_row << std::setw(w3) << std::setprecision(2 * pprec) << std::fixed << tot_nuc;
    o_row << std::setw(w3) << std::setprecision(2 * pprec) << std::fixed << tot_rho;

    mrcpp::print::separator(2, '-');
    println(2, o_row.str());
    mrcpp::print::separator(2, '-');
    print_utils::qmfunction(2, "Local density", rho_loc, t_loc);
    print_utils::qmfunction(2, "Allreduce density", rho_tot, t_com);
    mrcpp::print::footer(2, t_tot, 2);
}

} // namespace mrchem
