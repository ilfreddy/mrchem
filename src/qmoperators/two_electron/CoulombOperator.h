#pragma once

#include "qmoperators/RankZeroTensorOperator.h"
#include "CoulombPotential.h"
#include "CoulombPotentialD1.h"
#include "CoulombPotentialD2.h"

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
            : potential(nullptr) {
        this->potential = new CoulombPotential(P);
        RankZeroTensorOperator &J = (*this);
        J = *this->potential;
    }
    CoulombOperator(mrcpp::PoissonOperator *P, OrbitalVector *Phi)
            : potential(nullptr) {
        this->potential = new CoulombPotentialD1(P, Phi);
        RankZeroTensorOperator &J = (*this);
        J = *this->potential;
    }
    CoulombOperator(mrcpp::PoissonOperator *P,
                    OrbitalVector *Phi,
                    OrbitalVector *X,
                    OrbitalVector *Y)
            : potential(nullptr) {
        this->potential = new CoulombPotentialD2(P, Phi, X, Y);
        RankZeroTensorOperator &J = (*this);
        J = *this->potential;
    }
    ~CoulombOperator() { if (this->potential != nullptr) delete this->potential; }

    Density &getDensity() {
        if (potential != nullptr) return this->potential->getDensity();
        MSG_FATAL("Coulomb operator not properly initialized");
    }
        
    ComplexDouble trace(OrbitalVector &Phi) { return 0.5*RankZeroTensorOperator::trace(Phi); }

protected:
    CoulombPotential *potential;
};

} //namespace mrchem
