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

#include "SCFSolver.h"

/** @class LinearResponseSolver
 *
 */

namespace mrchem {

class FockOperator;

class LinearResponseSolver final : public SCFSolver {
public:
    LinearResponseSolver(bool dyn, FockOperator &F_0, OrbitalVector &Phi_0, ComplexMatrix &F_mat_0);
    ~LinearResponseSolver() override = default;

    bool optimize(double omega, FockOperator &F_1, OrbitalVector &X, OrbitalVector &Y);
    void setOrthPrec(double prec) { this->orth_prec = prec; }

protected:
    bool dynamic{false};
    double orth_prec{mrcpp::MachineZero};

    FockOperator *f_oper_0{nullptr};
    ComplexMatrix *f_mat_0{nullptr};
    OrbitalVector *phi_0{nullptr};

    void printProperty() const;
    void printParameters(double omega, const std::string &oper) const;
};

} // namespace mrchem
