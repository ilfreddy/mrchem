#pragma once

#include "CoulombPotential.h"

namespace mrchem {

class CoulombPotentialD1 final : public CoulombPotential {
public:
    CoulombPotentialD1(mrcpp::PoissonOperator *P, OrbitalVector *Phi);

protected:
    OrbitalVector *orbitals; ///< Unperturbed orbitals defining the ground-state electron density

    void setupDensity(double prec);
};

} //namespace mrchem
