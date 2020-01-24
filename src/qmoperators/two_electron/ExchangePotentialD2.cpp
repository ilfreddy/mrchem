#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "ExchangePotentialD2.h"
#include "parallel.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/OrbitalIterator.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator = mrcpp::PoissonOperator;
using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {

/** @brief constructor
 *
 * @param[in] P Poisson operator (does not take ownership)
 * @param[in] Phi vector of orbitals which define the exchange operator
 */
ExchangePotentialD2::ExchangePotentialD2(PoissonOperator_p P,
                                         OrbitalVector_p Phi,
                                         OrbitalVector_p X,
                                         OrbitalVector_p Y,
                                         bool s)
        : ExchangePotential(P, Phi, s)
        , orbitals_x(X)
        , orbitals_y(Y) {
    if (&X == &Y) useOnlyX = true;
}

/** @brief Computes the exchange potential on the fly
 *
 *  \param[in] inp input orbital
 *
 * The exchange potential is computed and applied on the fly to the given orbital.
 */
Orbital ExchangePotentialD2::calcExchange(Orbital phi_p) {
    Timer timer;
    Orbital ex_p;
    if (useOnlyX) {
        ex_p = calcExchange_X(phi_p);
    } else {
        ex_p = calcExchange_XY(phi_p);
    }

    //    print_utils::qmfunction(3, "Applied exchange", ex_p, timer);
    return ex_p;
}

/** @brief Computes the diagonal part of the internal exchange potential
 *
 *  \param[in] i orbital index
 *
 * The diagonal term K_ii is computed.
 */
void ExchangePotentialD2::calcInternal(int i) {
    if (useOnlyX) {
        calcInternal_X(i);
    } else {
        calcInternal_XY(i);
    }
}

/** @brief computes the off-diagonal part of the exchange potential
 *
 *  \param[in] i first orbital index
 *  \param[in] j second orbital index
 *
 * The off-diagonal terms K_ij and K_ji are computed.
 */
void ExchangePotentialD2::calcInternal(int i, int j) {
    if (useOnlyX) {
        calcInternal_X(i);
    } else {
        calcInternal_XY(i);
    }
}

void ExchangePotentialD2::calcInternal_X(int i) {
    NOT_IMPLEMENTED_ABORT;
}
void ExchangePotentialD2::calcInternal_X(int i, int j) {
    NOT_IMPLEMENTED_ABORT;
}
void ExchangePotentialD2::calcInternal_XY(int i) {
    NOT_IMPLEMENTED_ABORT;
}
void ExchangePotentialD2::calcInternal_XY(int i, int j) {
    NOT_IMPLEMENTED_ABORT;
}
Orbital ExchangePotentialD2::calcExchange_X(Orbital phi_p) {
    NOT_IMPLEMENTED_ABORT;
}
Orbital ExchangePotentialD2::calcExchange_XY(Orbital phi_p) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
