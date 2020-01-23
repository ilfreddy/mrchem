#pragma once

#include "ExchangePotentialD1.h"
#include "ExchangePotentialD2.h"
#include "qmoperators/RankZeroTensorOperator.h"

/**
 * @class ExchangeOperator
 *
 * @brief Operator containing a single ExchangePotential
 *
 * The operator is defined as the Hartree-Fock exchange arising from a particular
 * set of orbitals. The OrbitalVector defining the operator is fixed throughout the
 * operators life time, but the orbitals themselves are allowed to change in between
 * each application. The internal exchange potentials (the operator applied to it's
 * own orbitals) can be precomputed and stored for fast retrieval.
 */

namespace mrchem {

class ExchangeOperator final : public RankZeroTensorOperator {
public:
    ExchangeOperator(std::shared_ptr<mrcpp::PoissonOperator> P,
                     std::shared_ptr<OrbitalVector> Phi,
                     bool screen = false) {
        exchange = std::make_shared<ExchangePotentialD1>(P, Phi, screen);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &K = (*this);
        K = exchange;
        K.name() = "K";
    }

    ExchangeOperator(std::shared_ptr<mrcpp::PoissonOperator> P,
                     std::shared_ptr<OrbitalVector> Phi,
                     std::shared_ptr<OrbitalVector> X,
                     std::shared_ptr<OrbitalVector> Y,
                     bool screen = false) {
        exchange = std::make_shared<ExchangePotentialD2>(P, Phi, X, Y, screen);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &K = (*this);
        K = exchange;
        K.name() = "K";
    }

    ~ExchangeOperator() override = default;

    auto &getPoisson() { return exchange->getPoisson(); }
    void setupInternal(double prec) { exchange->setupInternal(prec); }
    void rotate(const ComplexMatrix &U) { exchange->rotate(U); }

    ComplexDouble trace(OrbitalVector &Phi) { return 0.5 * RankZeroTensorOperator::trace(Phi); }

private:
    std::shared_ptr<ExchangePotential> exchange{nullptr};
};

} // namespace mrchem
