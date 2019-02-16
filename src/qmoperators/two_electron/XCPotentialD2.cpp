#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "XCPotential.h"
#include "XCPotentialD2.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using mrcpp::ABGVOperator;
using mrcpp::DerivativeOperator;
using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

using XCFunctional = mrdft::XCFunctional;
using XCFunctional_p = std::shared_ptr<mrdft::XCFunctional>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {

/** @brief Constructor
 *
 * @param[in] F XCFunctional pointer
 * @param[in] Phi Vector of orbitals
 *
 * Based on the order and spin the correct nr. of potential functions is determined.
 * Then the functional is set up for subsequent calculations, fixing some internals of
 * xcfun when F.evalSetup is invoked.
 */
XCPotentialD2::XCPotentialD2(XCFunctional_p F,
                             OrbitalVector_p Phi,
                             OrbitalVector_p X,
                             OrbitalVector_p Y,
                             bool mpi_shared)
        : XCPotential(F, Phi, mpi_shared)
        , orbitals_x(X)
        , orbitals_y(Y)
        , pertDensity_t(nullptr)
        , pertDensity_a(nullptr)
        , pertDensity_b(nullptr) {}

XCPotentialD2::~XCPotentialD2() {
    mrcpp::clear(this->potentials, true);
    if (this->pertDensity_t != nullptr) MSG_ABORT("Operator not properly cleared");
    if (this->pertDensity_a != nullptr) MSG_ABORT("Operator not properly cleared");
    if (this->pertDensity_b != nullptr) MSG_ABORT("Operator not properly cleared");
}

/** @brief Prepare the operator for application
 *
 * @param[in] prec Apply precision
 *
 * Sequence of steps required to compute the XC potentials:
 *
 * 1) Compute density
 * 2) Setup xcfun input functions (gradients etc.)
 * 3) Evaluate xcfun
 * 4) Compute XC potential(s) from xcfun output
 *
 */
void XCPotentialD2::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);
    setupDensity(prec);
    setupPerturbedDensity(prec);
    setupPotential(prec);
}

/** @brief Clears all data in the XCPotentialD2 object */
void XCPotentialD2::clear() {
    this->energy = 0.0;
    if (this->pertDensity_t != nullptr) delete this->pertDensity_t;
    if (this->pertDensity_a != nullptr) delete this->pertDensity_a;
    if (this->pertDensity_b != nullptr) delete this->pertDensity_b;
    this->pertDensity_t = nullptr;
    this->pertDensity_a = nullptr;
    this->pertDensity_b = nullptr;
    clearApplyPrec();
}

// LUCA This does not work in the case of a non spin separated functional used for an open-shell system!!
void XCPotentialD2::setupPerturbedDensity(double prec) {
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_x == nullptr) MSG_ERROR("X-Orbitals not initialized");
    if (this->orbitals_y == nullptr) MSG_ERROR("Y-Orbitals not initialized");
    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &X = *this->orbitals_x;
    OrbitalVector &Y = *this->orbitals_y;

    if (this->functional->isSpinSeparated()) {
        this->pertDensity_a = buildPerturbedDensity(Phi, X, Y, DENSITY::Alpha);
        this->pertDensity_b = buildPerturbedDensity(Phi, X, Y, DENSITY::Beta);
    } else {
        this->pertDensity_t = buildPerturbedDensity(Phi, X, Y, DENSITY::Total);
    }
}

Density *XCPotentialD2::buildPerturbedDensity(OrbitalVector &Phi,
                                              OrbitalVector &X,
                                              OrbitalVector &Y,
                                              int density_spin) {
    Timer time;
    FunctionTree<3> &rho = this->getDensity(density_spin);
    pert_dens = new Density(false);
    Density &dRho = *pert_dens;
    mrcpp::build_grid(dRho.real(), rho);
    density::compute(-1.0, dRho, Phi, X, Y, density_spin);
    time.stop();
    Printer::printTree(0, "XC perturbed density", dRho_b.getNNodes(NUMBER::Total), time_b.getWallTime());
    FunctionTreeVector<3> & dv = this->functional->getDensityVector(density_spin);
    dv.push_back(std::make_tuple(1.0, dRho.real()));
}

/** @brief Compute XC potential(s)
 *
 * @param prec Precision used in refinement of density grid
 *
 * This will invoke a sequence of steps in the XCFunctional to compute the final
 * XC potential(s) that define this operator. Assuming the density has alredy been
 * computed:
 *
 * 1) Setup xcfun input functions (gradients etc.)
 * 2) Evaluate xcfun
 * 3) Compute XC energy by integrating energy density
 * 4) Compute XC potential(s) from xcfun output functions
 * 5) Remove excess grid nodes based on precision
 * 6) Add extra grid nodes based on precision
 * 7) Clear internal functions in XCFunctional (density grid is kept)
 *
 */
void XCPotentialD2::setupPotential(double prec) {
    if (this->functional == nullptr) MSG_ERROR("XCFunctional not initialized");
    if (not this->functional->hasDensity()) MSG_ERROR("XC density not initialized");
    if (this->potentials.size() != 0) MSG_ERROR("Potential not properly cleared");

    int inpNodes = this->functional->getNNodes();
    int inpPoints = this->functional->getNPoints();

    this->functional->setup();
    this->functional->evaluate();
    this->potentials = this->functional->calcPotential();
    this->functional->clear();

    int newNodes = this->functional->getNNodes() - inpNodes;
    int newPoints = this->functional->getNPoints() - inpPoints;

    println(0, " XC grid size   " << std::setw(26) << inpNodes << std::setw(17) << inpPoints);
    println(0, " XC grid change " << std::setw(26) << newNodes << std::setw(17) << newPoints);
}

} // namespace mrchem
