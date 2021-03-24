#pragma once

#include "MomentumOperator.h"
#include "ZoraOperator.h"
#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/ZoraOperator.h"

/** @class KineticOperator
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

class ZoraKineticOperator final : public KineticOperator {
public:
    ZoraKineticOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D, 
                        const Nuclei &nuclei, 
                        double zora_factor, 
                        double proj_prec,  
                        double smooth_prec, 
                        bool shared_memory)
                : p(D)
                , vz(nuclei, zora_factor, proj_prec, smooth_prec, shared_memory) {

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &t = (*this);
        t = 0.5 * (p[0] * vz * p[0] + p[1] * vz * p[1] + p[2] * vz * p[2]);
        t.name() = "T";
    }

    ComplexMatrix operator()(OrbitalVector &bra, OrbitalVector &ket) override;
    ComplexMatrix dagger(OrbitalVector &bra, OrbitalVector &ket) override;

    using RankZeroTensorOperator::operator();
    using RankZeroTensorOperator::dagger;

private:
    MomentumOperator p;
    ZoraOperator vz;
};

} // namespace mrchem
