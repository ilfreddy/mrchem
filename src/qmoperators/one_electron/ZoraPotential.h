#pragma once

#include "analyticfunctions/NuclearFunction.h"
#include "qmoperators/RankZeroTensorOperator.h"
#include "qmoperators/one_electron/QMPotential.h"

namespace mrchem {

class ZoraPotential final : public QMPotential {
public:
    ZoraPotential(const Nuclei &nucs, double proj_prec, double smooth_prec = -1.0, bool mpi_share = false);
    ~ZoraPotential() override { free(NUMBER::Total); }

    void setup(double prec) override { setApplyPrec(prec); }
    void clear() override { clearApplyPrec(); }
    void computeZora(double prec = -1.0);
    void computeKappa(double prec = -1.0);

private:
    NuclearFunction func;
    void allreducePotential(double prec, QMFunction &V_loc);
};

} // namespace mrchem
