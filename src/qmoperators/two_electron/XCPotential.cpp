#include "XCPotential.h"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "parallel.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Compute electron density
 *
 * The density is computed on the grid provided by the MRDFT module. The grid
 * is kept as is, e.i. no additional refinement at this point, since the grid
 * size is determined inside the module.
 */
void XCPotential::setupDensity(double prec) {
    if (this->functional->hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    OrbitalVector &Phi = *this->orbitals;
    if (this->functional->isSpinSeparated()) {
        buildDensity(Phi, DENSITY::DensityType::Alpha, prec);
        buildDensity(Phi, DENSITY::DensityType::Beta, prec);
        FunctionTree<3> &func_a = this->getDensity(DENSITY::DensityType::Alpha);
        FunctionTree<3> &func_b = this->getDensity(DENSITY::DensityType::Beta);
        println(5, static_cast<mrcpp::MWTree<3> &>(func_a));
        println(5, static_cast<mrcpp::MWTree<3> &>(func_b));
        while (mrcpp::refine_grid(func_a, func_b)) {}
        while (mrcpp::refine_grid(func_b, func_a)) {}
        println(5, static_cast<mrcpp::MWTree<3> &>(func_a));
        println(5, static_cast<mrcpp::MWTree<3> &>(func_b));
    } else {
        buildDensity(Phi, DENSITY::DensityType::Total, prec);
    }
}

/** @brief Clears all data in the XCPotential object */
void XCPotential::clear() {
    this->energy = 0.0;
    mrcpp::clear(this->potentials, true);
    clearApplyPrec();
}

void XCPotential::buildDensity(OrbitalVector &Phi, DENSITY::DensityType spin, double prec) {
    Timer time;
    FunctionTree<3> &func = this->getDensity(spin);
    Density rho(false);
    rho.setReal(&func);
    density::compute(prec, rho, Phi, spin);
    rho.setReal(nullptr);
}

mrcpp::FunctionTree<3> &XCPotential::getDensity(DENSITY::DensityType spin, int index) {
    return this->functional->getDensity(spin, index);
}

/** @brief Return FunctionTree for the XC spin potential
 *
 * @param[in] type Which spin potential to return (alpha, beta or total)
 */
FunctionTree<3> &XCPotential::getPotential(int spin) {
    bool spinFunctional = this->functional->isSpinSeparated();
    int pot_idx = -1;
    if (spinFunctional and spin == SPIN::Alpha) {
        pot_idx = 0;
    } else if (spinFunctional and spin == SPIN::Beta) {
        pot_idx = 1;
    } else if (not spinFunctional) {
        pot_idx = 0;
    } else {
        NOT_IMPLEMENTED_ABORT;
    }
    return mrcpp::get_func(this->potentials, pot_idx);
}

/** @brief XCPotentialD1 application
 *
 * @param[in] phi Orbital to which the potential is applied
 *
 * The operator is applied by choosing the correct potential function
 * which is then assigned to the real function part of the operator
 * base-class before the base class function is called.
 */
Orbital XCPotential::apply(Orbital phi) {
    QMPotential &V = *this;
    if (V.hasImag()) MSG_ERROR("Imaginary part of XC potential non-zero");

    FunctionTree<3> &pot = getPotential(phi.spin());
    V.setReal(&pot);
    Orbital Vphi = QMPotential::apply(phi);
    V.setReal(nullptr);

    if (phi.spin() == SPIN::Alpha) pot.setName("v_xc alpha");
    if (phi.spin() == SPIN::Beta) pot.setName("v_xc beta");
    println(5, static_cast<mrcpp::MWTree<3> &>(pot));
    if (phi.hasReal()) println(5, static_cast<mrcpp::MWTree<3> &>(phi.real()));
    if (Vphi.hasReal()) println(5, static_cast<mrcpp::MWTree<3> &>(Vphi.real()));

    return Vphi;
}

// NOTE AFTER DISCUSSION WITH STIG: Need to move stuff that is
// iteration-independent out of the response loop, so that all required
// functions, which only depend on the GS density are computed
// once. Comment: still the grid for rho_1 is borrowed from rho_0 and
// rho_0 should still be available.

} // namespace mrchem
