#pragma once

#include "CoulombPotential.h"

namespace mrchem {

class CoulombPotentialD1 final : public CoulombPotential {
public:
    CoulombPotentialD1(std::shared_ptr<mrcpp::PoissonOperator> P,
                       std::shared_ptr<OrbitalVector> Phi,
                       bool mpi_share = false)
            : CoulombPotential(P, Phi, mpi_share) {}

private:
    void setupLocalDensity(double prec) override;
    void setupGlobalDensity(double prec) override;
};

} // namespace mrchem
