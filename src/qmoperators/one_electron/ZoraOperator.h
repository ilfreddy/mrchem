#pragma once

#include "analyticfunctions/NuclearFunction.h"
#include "qmoperators/RankZeroTensorOperator.h"
#include "qmoperators/one_electron/QMPotential.h"
#include "qmoperators/one_electron/ZoraPotential.h"

namespace mrchem {

class ZoraOperator final : public RankZeroTensorOperator {
public:
    ZoraOperator(const Nuclei &nucs,
                 double proj_prec,
                 double smooth_prec = -1.0,
                 bool mpi_share = false,
                 bool inverse = false) {
        r_m1 = std::make_shared<ZoraPotential>(nucs, proj_prec, smooth_prec, mpi_share, inverse);
        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &v = (*this);
        v = r_m1;
        v.name() = "V_nuc";
    }

private:
    std::shared_ptr<ZoraPotential> r_m1{nullptr};
};

} // namespace mrchem
