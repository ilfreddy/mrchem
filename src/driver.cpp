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

#include "driver.h"

#include "chemistry/Molecule.h"
#include "chemistry/Nucleus.h"

#include "initial_guess/core.h"
#include "initial_guess/gto.h"
#include "initial_guess/sad.h"

#include "utils/MolPlotter.h"
#include "utils/math_utils.h"

#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"

#include "qmoperators/one_electron/ElectricFieldOperator.h"
#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/CoulombOperator.h"
#include "qmoperators/two_electron/ExchangeOperator.h"
#include "qmoperators/two_electron/FockOperator.h"
#include "qmoperators/two_electron/XCOperator.h"

#include "qmoperators/one_electron/H_BB_dia.h"
#include "qmoperators/one_electron/H_BM_dia.h"
#include "qmoperators/one_electron/H_B_dip.h"
#include "qmoperators/one_electron/H_B_spin.h"
#include "qmoperators/one_electron/H_E_dip.h"
#include "qmoperators/one_electron/H_M_fc.h"
#include "qmoperators/one_electron/H_M_pso.h"
#include "qmoperators/one_electron/NuclearGradientOperator.h"

#include "scf_solver/EnergyOptimizer.h"
#include "scf_solver/KAIN.h"
#include "scf_solver/LinearResponseSolver.h"
#include "scf_solver/OrbitalOptimizer.h"

#include "mrdft/XCFunctional.h"

using mrcpp::ABGVOperator;
using mrcpp::Coord;
using mrcpp::PoissonOperator;
using mrcpp::Printer;
using mrcpp::Timer;

using mrdft::XCFunctional;

using nlohmann::json;

using DerivativeOperator = mrcpp::DerivativeOperator<3>;
using DerivativeOperator_p = std::shared_ptr<mrcpp::DerivativeOperator<3>>;

using PoissonOperator = mrcpp::PoissonOperator;
using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;

extern mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

namespace mrchem {

namespace driver {
void build_fock_operator(const json &input, Molecule &mol, FockOperator &F, int order);

void calc_scf_properties(const json &input, Molecule &mol);
void calc_rsp_properties(const json &input, Molecule &mol, int dir, double omega);

void plot_scf_quantities(const json &input, Molecule &mol);

RankOneTensorOperator<3> get_perturbation(const json &input);
DerivativeOperator_p get_derivative(const std::string &name);
} // namespace driver

/** @brief Initialize a molecule from input
 *
 * This function expects the "molecule" subsection of the input.
 */
void driver::init_molecule(const json &json_mol, Molecule &mol) {
    print_utils::headline(0, "Initializing Molecule");

    auto charge = json_mol["charge"].get<int>();
    auto multiplicity = json_mol["multiplicity"].get<int>();

    mol.setCharge(charge);
    mol.setMultiplicity(multiplicity);

    Nuclei &nuclei = mol.getNuclei();
    for (const auto &coord : json_mol["coords"].get<json>()) {
        auto atom = coord["atom"].get<std::string>();
        auto xyz = coord["xyz"].get<std::array<double, 3>>();
        nuclei.push_back(atom, xyz);
    }

    mol.printGeometry();
}

/** @brief Run initial guess calculation for the orbitals
 *
 * This function will update the ground state orbitals and the Fock
 * matrix of the molecule, based on the chosen initial guess method.
 * The orbital vector is initialized with the appropriate particle
 * number and spin. The Fock matrix is initialized to the zero matrix
 * of the appropriate size.
 *
 * This function expects the "initial_guess" subsection of the input.
 */
bool driver::run_guess(const json &json_guess, Molecule &mol) {
    print_utils::headline(0, "Computing Initial Guess Wavefunction");

    auto &Phi = mol.getOrbitals();
    auto &F_mat = mol.getFockMatrix();

    auto method = json_guess["method"].get<std::string>();
    if (method == "mw") {
        auto start_orbs = json_guess["start_orbitals"].get<std::string>();
        mrcpp::print::separator(0, '~');
        print_utils::text(0, "Calculation ", "Read orbitals from file (MW)");
        print_utils::text(0, "File name   ", start_orbs);
        mrcpp::print::separator(0, '~', 2);
        Phi = orbital::load_orbitals(start_orbs);
    } else if (method == "core") {
        auto guess_prec = json_guess["guess_prec"].get<double>();
        auto restricted = json_guess["restricted"].get<bool>();
        auto zeta = json_guess["zeta"].get<int>();
        Phi = initial_guess::core::setup(guess_prec, mol, restricted, zeta);
    } else if (method == "sad") {
        auto guess_prec = json_guess["guess_prec"].get<double>();
        auto restricted = json_guess["restricted"].get<bool>();
        auto zeta = json_guess["zeta"].get<int>();
        Phi = initial_guess::sad::setup(guess_prec, mol, restricted, zeta);
    } else if (method == "gto") {
        auto guess_prec = json_guess["guess_prec"].get<double>();
        auto restricted = json_guess["restricted"].get<bool>();
        auto file_basis = json_guess["file_basis"].get<std::string>();
        auto file_moa = json_guess["file_moa"].get<std::string>();
        auto file_mob = json_guess["file_mob"].get<std::string>();
        if (restricted) {
            Phi = initial_guess::gto::setup(guess_prec, mol, file_basis, file_moa);
        } else {
            Phi = initial_guess::gto::setup(guess_prec, mol, file_basis, file_moa, file_mob);
        }
    } else {
        MSG_ERROR("Invalid initial guess");
        return false;
    }
    orbital::print(Phi);

    auto write_orbs = json_guess["write_orbitals"].get<bool>();
    auto final_orbs = json_guess["final_orbitals"].get<std::string>();
    if (write_orbs) orbital::save_orbitals(Phi, final_orbs);

    F_mat = ComplexMatrix::Zero(Phi.size(), Phi.size());

    return true;
}

/** @brief Run ground-state SCF calculation
 *
 * This function will update the ground state orbitals and the Fock
 * matrix of the molecule based on the chosen electronic structure
 * method. The resulting orbitals will be either diagonalized or
 * localized at exit. Returns true if the calculation converges.
 *
 * After convergence the requested ground-state properties are computed.
 *
 * This function expects the "scf_calculation" subsection of the input.
 */
bool driver::run_scf(const json &json_scf, Molecule &mol) {
    print_utils::headline(0, "Computing Ground State Wavefunction");

    const auto &json_fock = json_scf["fock_operator"].get<json>();
    FockOperator F;
    driver::build_fock_operator(json_fock, mol, F, 0);

    auto success = true;
    auto &Phi = mol.getOrbitals();
    auto &F_mat = mol.getFockMatrix();

    // Calc inital energy if present in input JSON
    auto initial_energy = json_scf.find("initial_energy");
    if (initial_energy != json_scf.end()) {
        auto prec = (*initial_energy)["prec"].get<double>();
        auto localize = (*initial_energy)["localize"].get<bool>();
        auto method = (*initial_energy)["method_name"].get<std::string>();

        std::stringstream o_prec;
        o_prec << std::setprecision(5) << std::scientific << prec;
        mrcpp::print::separator(0, '~');
        print_utils::text(0, "Calculation", "Compute initial energy");
        print_utils::text(0, "Method", method);
        print_utils::text(0, "Precision", o_prec.str());
        print_utils::text(0, "Localization", (localize) ? "On" : "Off");
        mrcpp::print::separator(0, '~', 2);

        Timer timer;
        auto plevel = Printer::getPrintLevel();
        if (plevel == 1) mrcpp::print::header(1, "Calculating Molecular Energy");
        F.setup(prec);
        F_mat = F(Phi, Phi);
        mol.getSCFEnergy() = F.trace(Phi, F_mat);
        F.clear();

        if (localize) {
            orbital::localize(prec, Phi, F_mat);
        } else {
            orbital::diagonalize(prec, Phi, F_mat);
        }
        if (plevel == 1) mrcpp::print::footer(1, timer, 2);

        mol.getSCFEnergy().print();
    }

    // Run OrbitalOptimizer if present in input JSON
    auto orbital_solver = json_scf.find("orbital_solver");
    if (orbital_solver != json_scf.end()) {
        auto method = (*orbital_solver)["method_name"].get<std::string>();
        auto kain = (*orbital_solver)["kain"].get<int>();
        auto max_iter = (*orbital_solver)["max_iter"].get<int>();
        auto rotation = (*orbital_solver)["rotation"].get<int>();
        auto localize = (*orbital_solver)["localize"].get<bool>();
        auto start_prec = (*orbital_solver)["start_prec"].get<double>();
        auto final_prec = (*orbital_solver)["final_prec"].get<double>();
        auto orbital_thrs = (*orbital_solver)["orbital_thrs"].get<double>();
        auto property_thrs = (*orbital_solver)["property_thrs"].get<double>();
        auto helmholtz_prec = (*orbital_solver)["helmholtz_prec"].get<double>();

        OrbitalOptimizer solver;
        solver.setMethodName(method);
        solver.setHistory(kain);
        solver.setMaxIterations(max_iter);
        solver.setRotation(rotation);
        solver.setLocalize(localize);
        solver.setHelmholtzPrec(helmholtz_prec);
        solver.setOrbitalPrec(start_prec, final_prec);
        solver.setThreshold(orbital_thrs, property_thrs);
        success = solver.optimize(mol, F);
    }

    // Run EnergyOptimizer if present in input JSON
    auto energy_solver = json_scf.find("energy_solver");
    if (energy_solver != json_scf.end()) {
        auto method = (*orbital_solver)["method_name"].get<std::string>();
        auto max_iter = (*energy_solver)["max_iter"].get<int>();
        auto localize = (*energy_solver)["localize"].get<bool>();
        auto start_prec = (*energy_solver)["start_prec"].get<double>();
        auto final_prec = (*energy_solver)["final_prec"].get<double>();
        auto orbital_thrs = (*energy_solver)["orbital_thrs"].get<double>();
        auto property_thrs = (*energy_solver)["property_thrs"].get<double>();
        auto helmholtz_prec = (*energy_solver)["helmholtz_prec"].get<double>();

        EnergyOptimizer solver;
        solver.setMethodName(method);
        solver.setMaxIterations(max_iter);
        solver.setRotation(1);
        solver.setLocalize(localize);
        solver.setHelmholtzPrec(helmholtz_prec);
        solver.setOrbitalPrec(start_prec, final_prec);
        solver.setThreshold(orbital_thrs, property_thrs);
        success = solver.optimize(mol, F);
    }

    if (success) {
        auto final_orbs = json_scf["final_orbitals"].get<std::string>();
        auto write_orbs = json_scf["write_orbitals"].get<bool>();
        if (write_orbs) orbital::save_orbitals(Phi, final_orbs);

        auto json_prop = json_scf["properties"].get<json>();
        calc_scf_properties(json_prop, mol);

        auto json_plot = json_scf["cube_plot"].get<json>();
        plot_scf_quantities(json_plot, mol);
    }

    return success;
}

/** @brief Run linear response SCF calculation
 *
 * This function will update the perturbed orbitals of the molecule
 * based on the chosen electronic structure method and perturbation
 * operator. Each response calculation corresponds to one particular
 * perturbation operator (could be a vector operator with several
 * components). Returns true if the calculation converges.
 *
 * After convergence the requested linear response properties are computed.
 *
 * This function expects a single subsection entry in the "rsp_calculations"
 * vector of the input.
 */
bool driver::run_rsp(const json &json_rsp, Molecule &mol) {
    mrcpp::print::header(0, "Response input");
    println(0, json_rsp.dump(2));
    mrcpp::print::separator(0, '=', 2);

    auto success = true;
    auto rsp_prec = json_rsp["rsp_prec"].get<double>();
    auto dynamic = json_rsp["dynamic"].get<bool>();
    auto localize = json_rsp["localize"].get<bool>();

    mol.initPerturbedOrbitals(dynamic);

    // Setup Fock operators
    const auto &json_fock_0 = json_rsp["fock_operator_0"].get<json>();
    const auto &json_fock_1 = json_rsp["fock_operator_1"].get<json>();
    FockOperator F_0, F_1;
    driver::build_fock_operator(json_fock_0, mol, F_0, 0);
    driver::build_fock_operator(json_fock_1, mol, F_1, 1);

    F_1.getXCOperator()->setupDensity(rsp_prec);
    //    F_1.getXCOperator()->setupPotential(rsp_prec);

    auto &F_mat = mol.getFockMatrix();
    auto &Phi = mol.getOrbitals();
    auto &X = mol.getOrbitalsX();
    auto &Y = mol.getOrbitalsY();

    // Setup perturbation operator
    const auto &json_pert = json_rsp["perturbation"].get<json>();
    auto h_1 = driver::get_perturbation(json_pert);

    ///////////////////////////////////////////////////////////
    /////////////////   Running RSP Solver  ///////////////////
    ///////////////////////////////////////////////////////////

    auto rsp_solver_it = json_rsp.find("rsp_solver");
    if (rsp_solver_it != json_rsp.end()) {
        const auto &json_solver = json_rsp["rsp_solver"].get<json>();
        auto kain = json_solver["kain"].get<int>();
        auto omega = json_solver["frequency"].get<double>();
        auto max_iter = json_solver["max_iter"].get<int>();
        auto directions = json_solver["directions"].get<std::array<int, 3>>();
        auto start_prec = json_solver["start_prec"].get<double>();
        auto final_prec = json_solver["final_prec"].get<double>();
        auto orbital_thrs = json_solver["orbital_thrs"].get<double>();
        auto property_thrs = json_solver["property_thrs"].get<double>();

        LinearResponseSolver solver(dynamic, F_0, Phi, F_mat);
        solver.setHistory(kain);
        solver.setMaxIterations(max_iter);
        solver.setLocalize(localize);
        solver.setOrbitalPrec(start_prec, final_prec);
        solver.setThreshold(orbital_thrs, property_thrs);

        F_0.setup(rsp_prec);
        for (int d = 0; d < 3; d++) {
            if (directions[d]) {
                F_1.perturbation() = h_1[d];
                X = orbital::param_copy(Phi);
                Y = orbital::param_copy(Phi);
                success = solver.optimize(omega, F_1, X, Y);
                if (success) {
                    const auto &json_prop = json_rsp["properties"].get<json>();
                    driver::calc_rsp_properties(json_prop, mol, d, omega);

                    auto final_orbs = json_rsp["final_orbitals"].get<std::string>();
                    auto write_orbs = json_rsp["write_orbitals"].get<bool>();
                    if (write_orbs) {
                        orbital::save_orbitals(X, final_orbs);
                        if (dynamic) orbital::save_orbitals(Y, final_orbs);
                    }
                }
                X.clear();
                Y.clear();
            }
            if (not(success)) break;
        }
        F_0.clear();
    }

    return success;
}

/** @brief Plot ground-state quantities
 *
 * This function expects the "cube_plot" subsection of the
 * "scf_calculation" input section.
 */
void driver::plot_scf_quantities(const json &json_plot, Molecule &mol) {
    Timer t_tot, t_lap;
    print_utils::headline(1, "Plotting Ground State Quantities");
    mrcpp::print::header(1, "CubePlot");

    auto npts = json_plot["plotter"]["points"].get<std::array<int, 3>>();
    auto O = json_plot["plotter"]["O"].get<mrcpp::Coord<3>>();
    auto A = json_plot["plotter"]["A"].get<mrcpp::Coord<3>>();
    auto B = json_plot["plotter"]["B"].get<mrcpp::Coord<3>>();
    auto C = json_plot["plotter"]["C"].get<mrcpp::Coord<3>>();
    auto dens_plot = json_plot["density"].get<bool>();
    auto orb_idx = json_plot["orbital"].get<std::vector<int>>();

    auto &Phi = mol.getOrbitals();
    MolPlotter plt(mol, O);
    plt.setRange(A, B, C);

    if (dens_plot) {
        Density rho(false);

        t_lap.start();
        std::string fname = "plots/rho_t";
        density::compute(-1.0, rho, Phi, DENSITY::DensityType::Total);
        plt.cubePlot(npts, rho, fname);
        rho.free(NUMBER::Total);
        mrcpp::print::time(1, fname, t_lap);

        if (orbital::size_singly(Phi) > 0) {
            t_lap.start();
            fname = "plots/rho_s";
            density::compute(-1.0, rho, Phi, DENSITY::DensityType::Spin);
            plt.cubePlot(npts, rho, fname);
            mrcpp::print::time(1, fname, t_lap);
            rho.free(NUMBER::Total);

            t_lap.start();
            fname = "plots/rho_a";
            density::compute(-1.0, rho, Phi, DENSITY::DensityType::Alpha);
            plt.cubePlot(npts, rho, fname);
            mrcpp::print::time(1, fname, t_lap);
            rho.free(NUMBER::Total);

            t_lap.start();
            fname = "plots/rho_b";
            density::compute(-1.0, rho, Phi, DENSITY::DensityType::Beta);
            plt.cubePlot(npts, rho, fname);
            rho.free(NUMBER::Total);
            mrcpp::print::time(1, fname, t_lap);
        }
    }

    // Plotting NO orbitals
    if (orb_idx.size() == 0) return;

    if (orb_idx[0] < 0) {
        // Plotting ALL orbitals
        for (auto i = 0; i < Phi.size(); i++) {
            if (not mpi::my_orb(Phi[i])) continue;
            t_lap.start();
            std::stringstream name;
            name << "plots/phi_" << i;
            plt.cubePlot(npts, Phi[i], name.str());
            mrcpp::print::time(1, name.str(), t_lap);
        }
    } else {
        // Plotting some orbitals
        for (auto &i : orb_idx) {
            if (not mpi::my_orb(Phi[i])) continue;
            t_lap.start();
            std::stringstream name;
            name << "plots/phi_" << i;
            plt.cubePlot(npts, Phi[i], name.str());
            mrcpp::print::time(1, name.str(), t_lap);
        }
    }

    mrcpp::print::footer(1, t_tot, 2);
}

/** @brief Compute ground-state properties
 *
 * This function expects the "properties" subsection of the "scf_calculation"
 * input section, and will compute all properties which are present in this input.
 * This includes the diamagnetic contributions to the magnetic response properties.
 */
void driver::calc_scf_properties(const json &json_prop, Molecule &mol) {
    Timer t_tot, t_lap;
    auto plevel = Printer::getPrintLevel();
    print_utils::headline(1, "Computing Ground State Properties");
    if (plevel == 1) mrcpp::print::header(1, "Calculating Molecular Properties");

    auto &nuclei = mol.getNuclei();
    auto &Phi = mol.getOrbitals();

    auto json_dipole = json_prop.find("dipole_moment");
    if (json_dipole != json_prop.end()) {
        t_lap.start();
        mrcpp::print::header(2, "Dipole moment");
        auto prec = (*json_dipole)["setup_prec"].get<double>();
        auto r_O = (*json_dipole)["origin"].get<Coord<3>>();

        DipoleMoment &mu = mol.getDipoleMoment();
        mu.getOrigin() = r_O;

        H_E_dip h(r_O);
        h.setup(prec);
        mu.getNuclear() = h.trace(nuclei).real();
        mu.getElectronic() = h.trace(Phi).real();
        h.clear();
        mrcpp::print::footer(2, t_lap, 2);
        if (plevel == 1) mrcpp::print::time(1, "Dipole moment", t_lap);
    }

    auto json_quadrupole = json_prop.find("quadrupole_moment");
    if (json_quadrupole != json_prop.end()) MSG_ERROR("Quadrupole moment not implemented");

    auto json_hyperpolarizability = json_prop.find("hyperpolarizability");
    if (json_hyperpolarizability != json_prop.end()) MSG_ERROR("Hyperpolarizability not implemented");

    auto json_gradient = json_prop.find("nuclear_gradient");
    if (json_gradient != json_prop.end()) MSG_ERROR("Nuclear gradient not implemented");

    auto json_mag = json_prop.find("magnetizability");
    if (json_mag != json_prop.end()) {
        t_lap.start();
        mrcpp::print::header(2, "Magnetizability (dia)");
        auto prec = (*json_mag)["setup_prec"].get<double>();
        auto r_O = (*json_mag)["origin"].get<Coord<3>>();

        Magnetizability &khi = mol.getMagnetizability();

        H_BB_dia h(r_O);
        h.setup(prec);
        khi.getOrigin() = r_O;
        khi.getDiamagnetic() = -h.trace(Phi).real();
        h.clear();
        mrcpp::print::footer(2, t_lap, 2);
        if (plevel == 1) mrcpp::print::time(1, "Magnetizability (dia)", t_lap);
    }

    auto json_nmr = json_prop.find("nmr_shielding");
    if (json_nmr != json_prop.end()) {
        t_lap.start();
        mrcpp::print::header(2, "NMR shielding (dia)");
        auto prec = (*json_nmr)["setup_prec"].get<double>();
        auto r_O = (*json_nmr)["origin"].get<Coord<3>>();
        auto nucleus_k = (*json_nmr)["nucleus_k"].get<std::vector<int>>();

        for (int k = 0; k < nucleus_k.size(); k++) {
            if (k != 0) mrcpp::print::separator(2, '-');
            NMRShielding &sigma_k = mol.getNMRShielding(nucleus_k[k]);
            const auto &r_K = sigma_k.getNucleus().getCoord();

            H_BM_dia h(r_O, r_K);
            h.setup(prec);
            sigma_k.getOrigin() = r_O;
            sigma_k.getDiamagnetic() = h.trace(Phi).real();
            h.clear();
        }
        mrcpp::print::footer(2, t_lap, 2);
        if (plevel == 1) mrcpp::print::time(1, "NMR shielding (dia)", t_lap);
    }
    if (plevel == 1) mrcpp::print::footer(1, t_tot, 2);
}

/** @brief Compute linear response properties
 *
 * This function expects the "properties" subsection of the "rsp_calculations"
 * input section, and will compute all properties which are present in this input.
 */
void driver::calc_rsp_properties(const json &json_prop, Molecule &mol, int dir, double omega) {
    auto &Phi = mol.getOrbitals();
    auto &X = mol.getOrbitalsX();
    auto &Y = mol.getOrbitalsY();

    auto json_pol = json_prop.find("polarizability");
    if (json_pol != json_prop.end()) {
        mrcpp::print::header(1, "Calculating polarizability");
        auto prec = (*json_pol)["setup_prec"].get<double>();
        auto r_O = (*json_pol)["origin"].get<Coord<3>>();

        Polarizability &alpha = mol.getPolarizability(omega);

        Timer timer;
        H_E_dip h(r_O);
        h.setup(prec);
        alpha.getTensor().row(dir) = -h.trace(Phi, X, Y).real();
        h.clear();
        timer.stop();
        mrcpp::print::footer(1, timer, 2);
    }

    auto json_mag = json_prop.find("magnetizability");
    if (json_mag != json_prop.end()) {
        mrcpp::print::header(1, "Calculating paramagnetic magnetizability");
        auto prec = (*json_mag)["setup_prec"].get<double>();
        auto r_O = (*json_mag)["origin"].get<Coord<3>>();
        auto pert_diff = (*json_mag)["derivative"].get<std::string>();
        auto D = driver::get_derivative(pert_diff);

        Magnetizability &khi = mol.getMagnetizability();

        Timer timer;
        H_B_dip h(D, r_O);
        h.setup(prec);
        khi.getParamagnetic().row(dir) = -h.trace(Phi, X, Y).real();
        h.clear();
        timer.stop();
        mrcpp::print::footer(1, timer, 2);
    }
}

/** @brief Build Fock operator based on input parameters
 *
 * This function expects the "fock_operator" subsection of input, and will
 * construct all operator which are present in this input. Option to set
 * perturbation order of the operators.
 */
void driver::build_fock_operator(const json &json_fock, Molecule &mol, FockOperator &F, int order) {
    auto &nuclei = mol.getNuclei();
    auto Phi_p = mol.getOrbitals_p();
    auto X_p = mol.getOrbitalsX_p();
    auto Y_p = mol.getOrbitalsY_p();

    ///////////////////////////////////////////////////////////
    //////////////////   Kinetic Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_kinetic = json_fock.find("kinetic_operator");
    if (json_kinetic != json_fock.end()) {
        auto kin_diff = (*json_kinetic)["derivative"].get<std::string>();
        auto D_p = driver::get_derivative(kin_diff);
        auto T_p = std::make_shared<KineticOperator>(D_p);
        F.getKineticOperator() = T_p;
    }
    ///////////////////////////////////////////////////////////
    //////////////////   Nuclear Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_nuclear = json_fock.find("nuclear_operator");
    if (json_nuclear != json_fock.end()) {
        auto proj_prec = (*json_nuclear)["proj_prec"].get<double>();
        auto smooth_prec = (*json_nuclear)["smooth_prec"].get<double>();
        auto shared_memory = (*json_nuclear)["shared_memory"].get<bool>();
        auto V_p = std::make_shared<NuclearOperator>(nuclei, proj_prec, smooth_prec, shared_memory);
        F.getNuclearOperator() = V_p;
    }
    ///////////////////////////////////////////////////////////
    //////////////////   Coulomb Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_coulomb = json_fock.find("coulomb_operator");
    if (json_coulomb != json_fock.end()) {
        auto poisson_prec = (*json_coulomb)["poisson_prec"].get<double>();
        auto shared_memory = (*json_coulomb)["shared_memory"].get<bool>();
        auto P_p = std::make_shared<PoissonOperator>(*MRA, poisson_prec);
        if (order == 0) {
            auto J_p = std::make_shared<CoulombOperator>(P_p, Phi_p, shared_memory);
            F.getCoulombOperator() = J_p;
        } else if (order == 1) {
            auto J_p = std::make_shared<CoulombOperator>(P_p, Phi_p, X_p, Y_p, shared_memory);
            F.getCoulombOperator() = J_p;
        } else {
            MSG_ABORT("Invalid perturbation orde>r");
        }
    }
    ///////////////////////////////////////////////////////////
    ////////////////////   XC Operator   //////////////////////
    ///////////////////////////////////////////////////////////
    double exx = 1.0;
    auto json_xc = json_fock.find("xc_operator");
    if (json_xc != json_fock.end()) {
        auto grid_prec = (*json_xc)["grid_prec"].get<double>();
        auto shared_memory = (*json_xc)["shared_memory"].get<bool>();
        auto json_xcfunc = (*json_xc)["xc_functional"].get<json>();
        auto xc_spin = json_xcfunc["spin"].get<bool>();
        auto xc_gamma = json_xcfunc["gamma"].get<bool>();
        auto xc_cutoff = json_xcfunc["cutoff"].get<double>();
        auto xc_diff = json_xcfunc["derivative"].get<std::string>();
        auto xc_funcs = json_xcfunc["functionals"].get<json>();

        int xc_order = order + 1;
        auto xcfun_p = std::make_shared<XCFunctional>(*MRA, xc_spin);
        for (const auto &f : xc_funcs) {
            auto name = f["name"].get<std::string>();
            auto coef = f["coef"].get<double>();
            xcfun_p->setFunctional(name, coef);
        }
        xcfun_p->setUseGamma(xc_gamma);
        xcfun_p->setDensityCutoff(xc_cutoff);
        xcfun_p->setNDensities(xc_order); // Nr of dens is the same as xc_order
        xcfun_p->evalSetup(xc_order);
        xcfun_p->allocateDensities();
        exx = xcfun_p->amountEXX();

        if (order == 0) {
            auto XC_p = std::make_shared<XCOperator>(xcfun_p, Phi_p, shared_memory);
            F.getXCOperator() = XC_p;
        } else if (order == 1) {
            auto XC_p = std::make_shared<XCOperator>(xcfun_p, Phi_p, X_p, Y_p, shared_memory);
            F.getXCOperator() = XC_p;
        } else {
            MSG_ABORT("Invalid perturbation order");
        }
    }
    ///////////////////////////////////////////////////////////
    /////////////////   Exchange Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_exchange = json_fock.find("exchange_operator");
    if (json_exchange != json_fock.end() and exx > mrcpp::MachineZero) {
        auto poisson_prec = (*json_exchange)["poisson_prec"].get<double>();
        auto screen_prec = (*json_exchange)["screen"].get<bool>();
        auto P_p = std::make_shared<PoissonOperator>(*MRA, poisson_prec);
        if (order == 0) {
            auto K_p = std::make_shared<ExchangeOperator>(P_p, Phi_p, screen_prec);
            F.getExchangeOperator() = K_p;
        } else {
            MSG_ABORT("Invalid perturbation order");
        }
    }
    ///////////////////////////////////////////////////////////
    /////////////////   External Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_external = json_fock.find("external_operator");
    if (json_external != json_fock.end()) {
        auto field = (*json_external)["electric_field"].get<std::array<double, 3>>();
        auto r_O = (*json_external)["origin"].get<Coord<3>>();
        auto V_ext = std::make_shared<ElectricFieldOperator>(field, r_O);
        F.getExtOperator() = V_ext;
    }
    F.build(exx);
}

/** @brief Construct perturbation operator based on input keyword */
RankOneTensorOperator<3> driver::get_perturbation(const json &json_pert) {
    RankOneTensorOperator<3> h_1;
    auto pert_oper = json_pert["operator"].get<std::string>();
    if (pert_oper == "h_e_dip") {
        auto r_O = json_pert["origin"].get<mrcpp::Coord<3>>();
        h_1 = H_E_dip(r_O);
    }
    if (pert_oper == "h_b_dip") {
        auto r_O = json_pert["origin"].get<mrcpp::Coord<3>>();
        auto pert_diff = json_pert["derivative"].get<std::string>();
        auto D = driver::get_derivative(pert_diff);
        h_1 = H_B_dip(D, r_O);
    }
    return h_1;
}

/** @brief Construct derivative operator based on input keyword */
DerivativeOperator_p driver::get_derivative(const std::string &name) {
    DerivativeOperator_p D = nullptr;
    if (name == "abgv_00") {
        D = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);
    } else if (name == "abgv_55") {
        D = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5);
        //    } else if (name == "ph") {
        //        D = std::make_shared<mrcpp::PHOperator<3>>(*MRA, 1);
        //    } else if (name == "bspline") {
        //        D = std::make_shared<mrcpp::BSOperator<3>>(*MRA, 1);
    } else {
        MSG_ERROR("Invalid derivative operator");
    }
    return D;
}

void driver::print_properties(const Molecule &mol) {
    print_utils::headline(0, "Printing Molecular Properties");
    mol.printGeometry();
    mol.printProperties();
}

} // namespace mrchem
