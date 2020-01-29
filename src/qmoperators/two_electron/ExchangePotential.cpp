#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "ExchangePotential.h"
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
ExchangePotential::ExchangePotential(PoissonOperator_p P, OrbitalVector_p Phi, bool s)
        : screen(s)
        , orbitals(Phi)
        , poisson(P) {
    int nOrbs = this->orbitals->size();
    this->tot_norms = DoubleVector::Zero(nOrbs);
    this->part_norms = DoubleMatrix::Zero(nOrbs, nOrbs);
}

/** @brief Test if a given contribution has been precomputed
 *
 * @param[in] phi_p orbital for which the check is performed
 *
 * If the given contribution has been precomputed, it is simply copied,
 * without additional recalculation.
 */
int ExchangePotential::testPreComputed(Orbital phi_p) const {
    const OrbitalVector &Phi = *this->orbitals;
    const OrbitalVector &Ex = this->exchange;

    int out = -1;
    if (Ex.size() == Phi.size()) {
        for (int i = 0; i < Phi.size(); i++) {
            if (&Phi[i].real() == &phi_p.real() and &Phi[i].imag() == &phi_p.imag()) {
                out = i;
                break;
            }
        }
    }
    return out;
}

/** @brief precomputes the exchange potential
 *
 *  @param[in] phi_p input orbital
 *
 * The exchange potential is (pre)computed among the orbitals that define the operator
 */
void ExchangePotential::setupInternal(double prec) {
    setApplyPrec(prec);

    if (this->exchange.size() != 0) MSG_ERROR("Exchange not properly cleared");

    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &Ex = this->exchange;

    Timer timer;
    // Diagonal must come first because it's NOT in-place
    for (int i = 0; i < Phi.size(); i++) calcInternal(i);

    // Off-diagonal must come last because it IS in-place
    OrbitalIterator iter(Phi, true); // symmetric iterator
    Orbital ex_rcv;
    while (iter.next(1)) { // one orbital at the time
        if (iter.get_size() > 0) {
            Orbital &phi_i = iter.orbital(0);
            int idx = iter.idx(0);
            for (int j = 0; j < Phi.size(); j++) {
                Orbital &phi_j = (*this->orbitals)[j];
                if (idx == j) continue; // skip diagonal terms
                if (mpi::my_orb(phi_i) and mpi::my_orb(phi_j)) calcInternal(idx, j);
            }
            // must send exchange_i to owner and receive exchange computed by other
            if (iter.get_step(0) and not mpi::my_orb(phi_i))
                mpi::send_function(Ex[idx], phi_i.rankID(), idx, mpi::comm_orb);

            if (iter.get_sent_size()) {
                // get exchange from where we sent orbital to
                int idx_sent = iter.get_idx_sent(0);
                int sent_rank = iter.get_rank_sent(0);
                mpi::recv_function(ex_rcv, sent_rank, idx_sent, mpi::comm_orb);
                Ex[idx_sent].add(1.0, ex_rcv);
            }

            if (not iter.get_step(0) and not mpi::my_orb(phi_i))
                mpi::send_function(Ex[idx], phi_i.rankID(), idx, mpi::comm_orb);
            if (not mpi::my_orb(Ex[idx])) Ex[idx].free(NUMBER::Total);
        } else {
            if (iter.get_sent_size()) { // must receive exchange computed by other
                // get exchange from where we sent orbital to
                int idx_sent = iter.get_idx_sent(0);
                int sent_rank = iter.get_rank_sent(0);
                mpi::recv_function(ex_rcv, sent_rank, idx_sent, mpi::comm_orb);
                Ex[idx_sent].add(1.0, ex_rcv);
            }
        }
        ex_rcv.free(NUMBER::Total);
    }

    // Collect info from the calculation
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) this->tot_norms(i) = Ex[i].norm();
    }

    mpi::allreduce_vector(this->tot_norms, mpi::comm_orb);  // to be checked
    mpi::allreduce_matrix(this->part_norms, mpi::comm_orb); // to be checked

    auto n = orbital::get_n_nodes(Ex);
    auto m = orbital::get_size_nodes(Ex);
    auto t = timer.elapsed();
    mrcpp::print::tree(2, "Hartree-Fock exchange", n, m, t);
}

/** @brief Perform a unitary transformation among the precomputed exchange contributions
 *
 * @param[in] U unitary matrix defining the rotation
 */
void ExchangePotential::rotate(const ComplexMatrix &U) {
    if (this->exchange.size() == 0) return;
    this->exchange = orbital::rotate(this->exchange, U, this->apply_prec);

    // NOTE: The following MPI point is currently NOT implemented!
    //
    // the last parameter, 1, means MPI will send only one orbital at a time
    // (because Exchange orbitals can be large for large molecules).
    // OrbitalAdder add(this->apply_prec, this->max_scale, 1);
    // add.rotate(this->exchange, U);
}

/** @brief determines the exchange factor to be used in the calculation of the exact exchange
 *
 * @param [in] phi_i orbital defining the K operator
 * @param [in] phi_j orbital to which K is applied
 *
 * The factor is computed in terms of the occupancy of the two orbitals and in terms of the spin
 * 0.5 factors are used in order to preserve occupancy of the set of doubly occupied orbitals
 * this-> is the orbital defining the operator whereas the input orbital (orb) is the one
 * the operator is applied to
 *
 * Occupancy: Single/Double
 * Spin: alpha/beta
 *
 * K (this->) | orb (input) | factor
 * alpha      | alpha       | 1.0
 * alpha      | beta        | 0.0
 * alpha      | double      | 0.5
 * -------------------------------
 * beta       | alpha       | 0.0
 * beta       | beta        | 1.0
 * beta       | double      | 0.5
 * -------------------------------
 * double     | alpha       | 1.0
 * double     | beta        | 1.0
 * double     | double      | 1.0
 *
 */
double ExchangePotential::getSpinFactor(Orbital phi_i, Orbital phi_j) const {
    double out = 0.0;
    if (phi_i.spin() == SPIN::Paired)
        out = 1.0;
    else if (phi_j.spin() == SPIN::Paired)
        out = 0.5;
    else if (phi_i.spin() == phi_j.spin())
        out = 1.0;
    return out;
}

/** @brief Prepare operator for application
 *
 * @param[in] prec reqested precision
 *
 * This will NOT precompute the internal exchange between the orbtials defining
 * the operator, which is done explicitly using setupInternal().
 */
void ExchangePotential::setup(double prec) {
    setApplyPrec(prec);

    int nOrbs = this->orbitals->size();
    if (tot_norms.size() != nOrbs) this->tot_norms = DoubleVector::Zero(nOrbs);
    if (part_norms.rows() != nOrbs) this->part_norms = DoubleMatrix::Zero(nOrbs, nOrbs);
    if (part_norms.cols() != nOrbs) this->part_norms = DoubleMatrix::Zero(nOrbs, nOrbs);
}

/** @brief Clears the Exchange Operator
 *
 *  Clears deletes the precomputed exchange contributions.
 */
void ExchangePotential::clear() {
    this->exchange.clear();
    clearApplyPrec();
}

/** @brief Applies operator potential
 *
 *  @param[in] inp input orbital
 *
 * The exchange potential is applied to the given orbital. Checks first if this
 * particular exchange contribution has been precomputed.
 */
Orbital ExchangePotential::apply(Orbital inp) {
    if (this->apply_prec < 0.0) {
        MSG_ERROR("Uninitialized operator");
        return inp.paramCopy();
    }
    int i = testPreComputed(inp);
    // if (i < 0) {
    if (true) {
        println(4, "On-the-fly exchange");
        return calcExchange(inp);
    } else {
        println(4, "Precomputed exchange");
        Orbital out = this->exchange[i].paramCopy();
        qmfunction::deep_copy(out, this->exchange[i]);
        return out;
    }
}

/** @brief Applies the adjoint of the operator
 *  \param[in] inp input orbital
 *
 * NOT IMPLEMENTED
 */
Orbital ExchangePotential::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

/** @brief scale the relative precision based on norm
 *
 * The internal norms are saved between SCF iterations so that they can
 * be used to estimate the size of the different contributions to the total
 * exchange. The relative precision of the Poisson terms is scaled to give a
 * consistent _absolute_ pecision in the final output.
 */
double ExchangePotential::getScaledPrecision(int i, int j) const {
    double scaled_prec = this->apply_prec;
    if (this->screen) {
        double tNorm = this->tot_norms(i);
        double pNorm = std::max(this->part_norms(i, j), this->part_norms(j, i));
        if (tNorm > 0.0) scaled_prec *= tNorm / pNorm;
    }
    return scaled_prec;
}

} // namespace mrchem
