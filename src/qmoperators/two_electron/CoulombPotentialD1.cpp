#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombPotentialD1.h"
#include "qmfunctions/density_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

CoulombPotentialD1::CoulombPotentialD1(mrcpp::PoissonOperator *P, OrbitalVector *Phi)
        : CoulombPotential(P, Phi) {}

/** @brief compute electron density
 *
 * @param[in] prec: apply precision
 *
 * This will compute the electron density as the sum of squares of the orbitals.
 */
void CoulombPotentialD1::setupGlobalDensity(double prec) {
    if (hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");

    OrbitalVector &Phi = *this->orbitals;
    Density &rho = this->density;

    Timer timer;
    density::compute(prec, rho, Phi, DENSITY::DensityType::Total);
    timer.stop();
    double t = timer.getWallTime();
    int n = rho.getNNodes(NUMBER::Total);
    Printer::printTree(0, "Coulomb density", n, t);
}

/** @brief compute local electron density (only own MPI orbitals)
 *
 * @param[in] prec: apply precision
 *
 * This will compute the electron density as the sum of squares of the orbitals.
 */
void CoulombPotentialD1::setupLocalDensity(double prec) {
    if (hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");

    OrbitalVector &Phi = *this->orbitals;
    Density &rho = this->density;

    Timer timer;
    density::compute_local(prec, rho, Phi, DENSITY::DensityType::Total);
    timer.stop();
    double t = timer.getWallTime();
    int n = rho.getNNodes(NUMBER::Total);
    Printer::printTree(0, "Coulomb density", n, t);
}

} // namespace mrchem
