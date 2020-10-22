/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "mrchem.h"
#include "qmfunctions/qmfunction_fwd.h"
#include "qmoperators/qmoperator_fwd.h"

/** @class HelmholtzVector
 *
 * @brief Container of HelmholtzOperators for a corresponding OrbtialVector
 *
 * This class assigns one HelmholtzOperator to each orbital in an OrbitalVector.
 * The operators are produced on the fly based on a vector of lambda parameters.
 */

class QMPotential;

namespace mrchem {

class HelmholtzVector final {
public:
    HelmholtzVector(double pr, const DoubleVector &l);

    DoubleMatrix getLambdaMatrix() const { return this->lambda.asDiagonal(); }

    OrbitalVector apply(RankZeroTensorOperator &V, OrbitalVector &Phi, OrbitalVector &Psi) const;
    OrbitalVector apply_zora(RankZeroTensorOperator &V,
                             RankZeroTensorOperator &GlnkG,
                             RankZeroTensorOperator &zora,
                             OrbitalVector &Phi,
                             OrbitalVector &Psi) const;
    OrbitalVector operator()(OrbitalVector &Phi) const;

private:
    double prec;         ///< Precision for construction and application of Helmholtz operators
    DoubleVector lambda; ///< Helmholtz parameter, mu_i = sqrt(-2.0*lambda_i)

    Orbital apply(int i, Orbital &phi) const;
};

} // namespace mrchem
