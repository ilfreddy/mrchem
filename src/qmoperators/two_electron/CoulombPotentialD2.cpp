#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombPotentialD2.h"
#include "qmfunctions/density_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator = mrcpp::PoissonOperator;
using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {

CoulombPotentialD2::CoulombPotentialD2(PoissonOperator_p P,
                                       OrbitalVector_p Phi,
                                       OrbitalVector_p X,
                                       OrbitalVector_p Y,
                                       bool mpi_share)
        : CoulombPotential(P, Phi, mpi_share)
        , orbitals_x(X)
        , orbitals_y(Y) {}

void CoulombPotentialD2::setupGlobalDensity(double prec) {
    if (hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_x == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_y == nullptr) MSG_ERROR("Orbitals not initialized");

    Density &rho = this->density;
    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &X = *this->orbitals_x;
    OrbitalVector &Y = *this->orbitals_y;

    Timer timer;
    density::compute(prec, rho, Phi, X, Y, DENSITY::DensityType::Total);
    print_utils::qmfunction(2, "Perturbed Coulomb density", rho, timer);
}

void CoulombPotentialD2::setupLocalDensity(double prec) {
    if (hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_x == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_y == nullptr) MSG_ERROR("Orbitals not initialized");

    Density &rho = this->density;
    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &X = *this->orbitals_x;
    OrbitalVector &Y = *this->orbitals_y;

    Timer timer;
    density::compute_local(prec, rho, Phi, X, Y, DENSITY::DensityType::Total);
    print_utils::qmfunction(2, "Perturbed Coulomb density", rho, timer);
}

} // namespace mrchem
