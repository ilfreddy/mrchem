#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombPotentialD1.h"
#include "qmfunctions/density_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator = mrcpp::PoissonOperator;
using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;

namespace mrchem {

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
    print_utils::qmfunction(2, "Coulomb density", rho, timer);
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
    density::compute_local(prec, rho, Phi, DENSITY::Total);
    print_utils::qmfunction(2, "Coulomb density", rho, timer);
}

} // namespace mrchem
