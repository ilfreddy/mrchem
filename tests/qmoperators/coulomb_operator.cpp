/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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

#include "catch.hpp"

#include "MRCPP/MWOperators"

#include "mrchem.h"
#include "parallel.h"

#include "analyticfunctions/HydrogenFunction.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/CoulombOperator.h"

using namespace mrchem;
using namespace orbital;

namespace nuclear_potential {

TEST_CASE("CoulombOperator", "[coulomb_operator]") {
    const double prec = 1.0e-5;
    const double thrs = prec * prec;

    const int nShells = 2;
    std::vector<int> ns;
    std::vector<int> ls;
    std::vector<int> ms;

    OrbitalVector Phi;
    for (int n = 1; n <= nShells; n++) {
        int L = n;
        for (int l = 0; l < L; l++) {
            int M = 2 * l + 1;
            for (int m = 0; m < M; m++) {
                ns.push_back(n);
                ls.push_back(l);
                ms.push_back(m);
                Phi.push_back(SPIN::Paired);
            }
        }
    }
    mpi::distribute(Phi);

    for (int i = 0; i < Phi.size(); i++) {
        HydrogenFunction f(ns[i], ls[i], ms[i]);
        if (mpi::my_orb(Phi[i])) Phi[i].alloc(NUMBER::Real);
        if (mpi::my_orb(Phi[i])) mrcpp::project(prec, Phi[i].real(), f);
    }

    // reference values
    int i = 0;
    DoubleVector E_P(Phi.size());
    for (int n = 1; n <= nShells; n++) {
        int L = n;
        double E_n = 1.0 / (2.0 * n * n); //E_n = Z^2/(2*n^2)
        for (int l = 0; l < L; l++) {
            int M = 2 * l + 1;
            for (int m = 0; m < M; m++) {
                E_P(i++) = -2.0 * E_n; //virial theorem: 2<E_K> = -<E_P>
            }
        }
    }

    Nuclei nucs;
    nucs.push_back("H", {0.0, 0.0, 0.0});

    mrcpp::PoissonOperator* P = new mrcpp::PoissonOperator(*MRA, prec);
    CoulombOperator V(P, &Phi);

    V.setup(prec);
    SECTION("apply") {
        Orbital Vphi_0 = V(Phi[0]);
        ComplexDouble V_00 = orbital::dot(Phi[0], Vphi_0);
        if (mpi::my_orb(Phi[0])) {
            REQUIRE(V_00.real() == Approx(V_00.real()).epsilon(prec));
            REQUIRE(V_00.imag() < thrs);
        } else {
            REQUIRE(V_00.real() < thrs);
            REQUIRE(V_00.imag() < thrs);
        }
        Vphi_0.free();
    }
    SECTION("vector apply") {
        OrbitalVector VPhi = V(Phi);
        for (int i = 0; i < Phi.size(); i++) {
            ComplexDouble V_ii = orbital::dot(Phi[i], VPhi[i]);
            if (mpi::my_orb(Phi[i])) {
                REQUIRE(V_ii.real() == Approx(V_ii.real()).epsilon(prec));
                REQUIRE(V_ii.imag() < thrs);
            } else {
                REQUIRE(V_ii.real() < thrs);
                REQUIRE(V_ii.imag() < thrs);
            }
        }
        free(VPhi);
    }
    SECTION("expectation value") {
        ComplexDouble V_00 = V(Phi[0], Phi[0]);
        if (mpi::my_orb(Phi[0])) {
            REQUIRE(V_00.real() == Approx(V_00.real()).epsilon(prec));
            REQUIRE(V_00.imag() < thrs);
        } else {
            REQUIRE(V_00.real() < thrs);
            REQUIRE(V_00.imag() < thrs);
        }
    }
    SECTION("expectation matrix ") {
        ComplexMatrix v = V(Phi, Phi);
        for (int i = 0; i < Phi.size(); i++) {
            for (int j = 0; j <= i; j++) {
                REQUIRE(v(i, j).real() == Approx(v(i, j).real()).epsilon(prec));
                REQUIRE(v(i, j).imag() < thrs);
            }
        }
    }
    V.clear();
    free(Phi);
}

} // namespace nuclear_potential

/* reference values
 0 0  3.1265429248
 1 0  0.2473632188
 1 1  1.6931905293
 2 2  1.8996591964
 3 3  1.8996591964
 4 4  1.8996591964
*/


