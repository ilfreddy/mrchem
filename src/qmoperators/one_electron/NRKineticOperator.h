#pragma once

#include "MomentumOperator.h"
#include "qmoperators/one_electron/KineticOperator.h"

/** @class NRKineticOperator
 *
 * @brief Operator for kinetic energy
 *
 * This operator is constructed as the square of the more fundamental
 * MomentumOperator. The general base class functions for calculation of
 * expectation values are overwritten, as they can be improved due to
 * symmetry.
 *
 */

namespace mrchem {

class NRKineticOperator final : public KineticOperator {
public:
    NRKineticOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D)
              : p(D) {
        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &t = (*this);
        t = 0.5 * (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        t.name() = "T";
    }

    ComplexMatrix operator()(OrbitalVector &bra, OrbitalVector &ket) override;
    ComplexMatrix dagger(OrbitalVector &bra, OrbitalVector &ket) override;

    using RankZeroTensorOperator::operator();
    using RankZeroTensorOperator::dagger;

private:
    MomentumOperator p;
};

} // namespace mrchem
