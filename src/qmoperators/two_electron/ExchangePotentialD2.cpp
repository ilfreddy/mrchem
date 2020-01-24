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
    Timer timer;

    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &X = *this->orbitals_x;

    ComplexVector coef_vec(Phi.size());
    QMFunctionVector func_vec;

    for (int i = 0; i < Phi.size(); i++) {
        Orbital &phi_i = Phi[i];
        Orbital &x_i = X[i];
        double spin_fac = getSpinFactor(phi_i, phi_p);
        if (std::abs(spin_fac) >= mrcpp::MachineZero) {
            coef_vec(i) = spin_fac / phi_i.squaredNorm();
            Orbital phi_iip = calcExchangeComponent(phi_p, phi_i, x_i);
            func_vec.push_back(phi_iip);
        }
    }

    // compute ex_p = sum_i c_i*phi_iip
    Orbital ex_p = phi_p.paramCopy();
    qmfunction::linear_combination(ex_p, coef_vec, func_vec, -1.0);
    print_utils::qmfunction(3, "Applied exchange", ex_p, timer);

    return ex_p;
}

Orbital ExchangePotentialD2::calcExchange_XY(Orbital phi_p) {
    NOT_IMPLEMENTED_ABORT;
}

Orbital ExchangePotentialD2::calcExchangeComponent(Orbital phi_p, Orbital phi_i, Orbital x_i) {
    double prec = this->apply_prec;
    mrcpp::PoissonOperator &P = *this->poisson;
    // compute phi_ip = phi_i^dag * phi_p
    Orbital phi_ip = phi_p.paramCopy();
    qmfunction::multiply(phi_ip, phi_i.dagger(), phi_p, -1.0);

    // compute V_ip = P[phi_ip]
    Orbital V_ip = phi_p.paramCopy();
    if (phi_ip.hasReal()) {
        V_ip.alloc(NUMBER::Real);
        mrcpp::apply(prec, V_ip.real(), P, phi_ip.real());
    }
    if (phi_ip.hasImag()) {
        V_ip.alloc(NUMBER::Imag);
        mrcpp::apply(prec, V_ip.imag(), P, phi_ip.imag());
    }
    phi_ip.release();

    // compute phi_iip = x_i * V_ip
    Orbital phi_iip = phi_p.paramCopy();
    qmfunction::multiply(phi_iip, x_i, V_ip, -1.0);
    return phi_iip;
}

} // namespace mrchem
