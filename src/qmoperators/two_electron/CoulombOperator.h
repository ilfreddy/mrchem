#pragma once

#include "CoulombPotential.h"
#include "CoulombPotentialD1.h"
#include "CoulombPotentialD2.h"
#include "qmoperators/RankZeroTensorOperator.h"

/** @class CoulombOperator
 *
 * @brief Operator containing a single CoulombPotential
 *
 * This class is a simple TensorOperator realization of @class CoulombPotential.
 *
 */

namespace mrchem {

class CoulombOperator final : public RankZeroTensorOperator {
public:
    CoulombOperator(mrcpp::PoissonOperator *P)
            : potential(new CoulombPotential(P)) {
        RankZeroTensorOperator &J = (*this);
        J = *this->potential;
    }
    CoulombOperator(mrcpp::PoissonOperator *P, OrbitalVector *Phi)
            : potential(new CoulombPotentialD1(P, Phi)) {
        RankZeroTensorOperator &J = (*this);
        J = *this->potential;
    }
    CoulombOperator(mrcpp::PoissonOperator *P, OrbitalVector *Phi, OrbitalVector *X, OrbitalVector *Y)
            : potential(new CoulombPotentialD2(P, Phi, X, Y)) {
        RankZeroTensorOperator &J = (*this);
        J = *this->potential;
    }
    ~CoulombOperator() {
        if (this->potential != nullptr) delete this->potential;
    }

    Density &getDensity() {
        if (potential != nullptr) return this->potential->getDensity();
        MSG_FATAL("Coulomb operator not properly initialized");
    }

    ComplexDouble trace(OrbitalVector &Phi) { return 0.5 * RankZeroTensorOperator::trace(Phi); }

private:
    CoulombPotential *potential;
};

} //namespace mrchem
