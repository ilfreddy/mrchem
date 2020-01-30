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
    if (X == Y) useOnlyX = true;
}

/** @brief Test if a given contribution has been precomputed
 *
 * @param[in] phi_p orbital for which the check is performed
 *
 * We always return -1 for response calculations
 *
 */
int ExchangePotentialD2::testPreComputed(Orbital phi_p) const {
    int out = -1;
    return out;
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
        calcInternal_X(i, j);
    } else {
        calcInternal_XY(i, j);
    }
}

void ExchangePotentialD2::setupInternal(double prec) {
    setApplyPrec(prec);
    if (this->exchange.size() != 0) MSG_ERROR("Exchange not properly cleared");
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
    std::vector<std::complex<double>> coef_vec;
    QMFunctionVector func_vec;

    for (int i = 0; i < Phi.size(); i++) {
        Orbital &phi_i = Phi[i];
        Orbital &x_i = X[i];
        double spin_fac = getSpinFactor(phi_i, phi_p);
        if (std::abs(spin_fac) >= mrcpp::MachineZero) {
            Orbital phi_apb = calcExchangeComponent(phi_p, x_i, phi_i);
            Orbital phi_bpa = calcExchangeComponent(phi_p, phi_i, x_i);
            func_vec.push_back(phi_apb);
            func_vec.push_back(phi_bpa);
            coef_vec.push_back(spin_fac / phi_i.squaredNorm());
            coef_vec.push_back(spin_fac / phi_i.squaredNorm());
        }
    }

    // compute ex_p = sum_i c_i* (phi_apb + phi_bpa)
    Orbital ex_p = phi_p.paramCopy();
    Eigen::Map<ComplexVector> coefs(coef_vec.data(), coef_vec.size());
    qmfunction::linear_combination(ex_p, coefs, func_vec, -1.0);
    print_utils::qmfunction(3, "Applied exchange", ex_p, timer);

    return ex_p;
}

Orbital ExchangePotentialD2::calcExchange_XY(Orbital phi_p) {
    NOT_IMPLEMENTED_ABORT;
}

Orbital ExchangePotentialD2::calcExchangeComponent(Orbital phi_p, Orbital phi_a, Orbital phi_b) {

    double prec = this->apply_prec;
    mrcpp::PoissonOperator &P = *this->poisson;
    // compute phi_pb = phi_p * phi_b.dagger()
    Orbital phi_pb = phi_p.paramCopy();
    qmfunction::multiply(phi_pb, phi_p, phi_b.dagger(), -1.0);

    // compute V_pb = P[phi_pb]
    Orbital V_pb = phi_p.paramCopy();
    if (phi_pb.hasReal()) {
        V_pb.alloc(NUMBER::Real);
        mrcpp::apply(prec, V_pb.real(), P, phi_pb.real());
    }
    if (phi_pb.hasImag()) {
        V_pb.alloc(NUMBER::Imag);
        mrcpp::apply(prec, V_pb.imag(), P, phi_pb.imag());
    }
    phi_pb.release();

    // compute phi_apb = phi_a * V_pb
    Orbital phi_apb = phi_p.paramCopy();
    qmfunction::multiply(phi_apb, phi_a, V_pb, -1.0);
    return phi_apb;
}

} // namespace mrchem
