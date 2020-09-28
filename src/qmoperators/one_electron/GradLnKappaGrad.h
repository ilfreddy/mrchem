#pragma once

#include "MomentumOperator.h"
#include "ZoraOperator.h"
#include "qmoperators/RankZeroTensorOperator.h"

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

class GradLnKappaGrad final : public RankZeroTensorOperator<3> {
public:
    GradLnKappaGrad(std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                    const Nuclei &nucs,
                    double proj_prec,
                    double smooth_prec = -1.0,
                    bool mpi_share = false)
            : p(D)
            , lnkappa(nucs, proj_prec, smooth_prec, mpi_share, 1) {
        // Invoke operator= to assign *this operator

        RankZeroTensorOperator &k = (*this);
        k = (p[0] * vz * p[0] + p[1] * vz * p[1] + p[2] * vz * p[2]);
        k.name() = "grad ln k grad";
    }

    ComplexMatrix operator()(OrbitalVector &bra, OrbitalVector &ket);
    ComplexMatrix dagger(OrbitalVector &bra, OrbitalVector &ket);

    using RankZeroTensorOperator::operator();
    using RankZeroTensorOperator::dagger;

private:
    MomentumOperator p;
    ZoraOperator vz;
};

} // namespace mrchem
