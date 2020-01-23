#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "ExchangePotentialD1.h"
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
ExchangePotentialD1::ExchangePotentialD1(PoissonOperator_p P, OrbitalVector_p Phi, bool s)
        : ExchangePotential(P, Phi, s) {}

/** @brief Computes the exchange potential on the fly
 *
 *  \param[in] inp input orbital
 *
 * The exchange potential is computed and applied on the fly to the given orbital.
 */
Orbital ExchangePotentialD1::calcExchange(Orbital phi_p) {
    Timer timer;

    double prec = this->apply_prec;
    OrbitalVector &Phi = *this->orbitals;
    mrcpp::PoissonOperator &P = *this->poisson;

    ComplexVector coef_vec(Phi.size());
    QMFunctionVector func_vec;

    OrbitalIterator iter(Phi);
    while (iter.next()) {
        for (int i = 0; i < iter.get_size(); i++) {
            Orbital &phi_i = iter.orbital(i);

            double spin_fac = getSpinFactor(phi_i, phi_p);
            if (std::abs(spin_fac) < mrcpp::MachineZero) continue;

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

            // compute phi_iip = phi_i * V_ip
            Orbital phi_iip = phi_p.paramCopy();
            qmfunction::multiply(phi_iip, phi_i, V_ip, -1.0);

            coef_vec(i) = spin_fac / phi_i.squaredNorm();
            func_vec.push_back(phi_iip);
        }
    }

    // compute ex_p = sum_i c_i*phi_iip
    Orbital ex_p = phi_p.paramCopy();
    qmfunction::linear_combination(ex_p, coef_vec, func_vec, -1.0);
    print_utils::qmfunction(3, "Applied exchange", ex_p, timer);

    return ex_p;
}

/** @brief Computes the diagonal part of the internal exchange potential
 *
 *  \param[in] i orbital index
 *
 * The diagonal term K_ii is computed.
 */
void ExchangePotentialD1::calcInternal(int i) {
    Orbital &phi_i = (*this->orbitals)[i];

    if (mpi::my_orb(phi_i)) {
        double prec = std::min(getScaledPrecision(i, i), 1.0e-1);
        mrcpp::PoissonOperator &P = *this->poisson;

        // compute phi_ii = phi_i^dag * phi_i
        Orbital phi_ii = phi_i.paramCopy();
        qmfunction::multiply(phi_ii, phi_i.dagger(), phi_i, prec);

        // compute V_ii = P[phi_ii]
        Orbital V_ii = phi_i.paramCopy();
        if (phi_ii.hasReal()) {
            V_ii.alloc(NUMBER::Real);
            mrcpp::apply(prec, V_ii.real(), P, phi_ii.real());
        }
        if (phi_ii.hasImag()) {
            V_ii.alloc(NUMBER::Imag);
            mrcpp::apply(prec, V_ii.imag(), P, phi_ii.imag());
        }
        phi_ii.release();

        // compute phi_iii = phi_i * V_ii
        Orbital phi_iii = phi_i.paramCopy();
        qmfunction::multiply(phi_iii, phi_i, V_ii, prec);
        phi_iii.rescale(1.0 / phi_i.squaredNorm());
        this->part_norms(i, i) = phi_iii.norm();
        this->exchange.push_back(phi_iii);
    } else {
        // put empty orbital to fill the exchange vector
        Orbital phi_iii = phi_i.paramCopy();
        this->exchange.push_back(phi_iii);
    }
}

/** @brief computes the off-diagonal part of the exchange potential
 *
 *  \param[in] i first orbital index
 *  \param[in] j second orbital index
 *
 * The off-diagonal terms K_ij and K_ji are computed.
 */
void ExchangePotentialD1::calcInternal(int i, int j, Orbital &phi_i, Orbital &phi_j) {
    mrcpp::PoissonOperator &P = *this->poisson;
    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &Ex = this->exchange;

    if (i == j) MSG_ABORT("Cannot handle diagonal term");
    if (Ex.size() != Phi.size()) MSG_ABORT("Size mismatch");
    if (phi_i.hasImag() or phi_j.hasImag()) MSG_ABORT("Orbitals must be real");

    double i_fac = getSpinFactor(phi_i, phi_j);
    double j_fac = getSpinFactor(phi_j, phi_i);

    double thrs = mrcpp::MachineZero;
    if (std::abs(i_fac) < thrs or std::abs(j_fac) < thrs) {
        this->part_norms(i, j) = 0.0;
        return;
    }

    // set correctly scaled precision for components ij and ji
    double prec = std::min(getScaledPrecision(i, j), getScaledPrecision(j, i));
    if (prec > 1.0e00) return;     // orbital does not contribute within the threshold
    prec = std::min(prec, 1.0e-1); // very low precision does not work properly

    // compute phi_ij = phi_i^dag * phi_j (dagger NOT used, orbitals must be real!)
    Orbital phi_ij = phi_i.paramCopy();
    qmfunction::multiply(phi_ij, phi_i, phi_j, prec);

    // compute V_ij = P[phi_ij]
    Orbital V_ij = phi_i.paramCopy();
    if (phi_ij.hasReal()) {
        V_ij.alloc(NUMBER::Real);
        mrcpp::apply(prec, V_ij.real(), P, phi_ij.real());
    }
    if (phi_ij.hasImag()) {
        MSG_ABORT("Orbitals must be real");
        V_ij.alloc(NUMBER::Imag);
        mrcpp::apply(prec, V_ij.imag(), P, phi_ij.imag());
    }
    phi_ij.release();

    // compute phi_jij = phi_j * V_ij
    Orbital phi_jij = phi_j.paramCopy();
    qmfunction::multiply(phi_jij, phi_j, V_ij, prec);
    phi_jij.rescale(1.0 / phi_j.squaredNorm());
    this->part_norms(j, i) = phi_jij.norm();

    // compute phi_iij = phi_i * V_ij
    Orbital phi_iij = phi_i.paramCopy();
    qmfunction::multiply(phi_iij, phi_i, V_ij, prec);
    phi_iij.rescale(1.0 / phi_i.squaredNorm());
    this->part_norms(i, j) = phi_iij.norm();
    V_ij.release();

    // compute x_i += phi_jij
    Ex[i].add(i_fac, phi_jij);
    phi_jij.release();

    // compute x_j += phi_iij
    Ex[j].add(j_fac, phi_iij);
    phi_iij.release();
}

} // namespace mrchem
