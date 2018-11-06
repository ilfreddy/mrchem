#pragma once

#include "RankOneTensorOperator.h"
#include "DeltaOperator.h"
#include "H_B_spin.h"

namespace mrchem {

class H_M_fc final : public RankOneTensorOperator<3> {
public:
    H_M_fc(const mrcpp::Coord<3> &o)
            : delta(o, 1.0e6) {
        const double coef = -(8.0 / 3.0) * MATHCONST::pi;
        const double alpha_2 = PHYSCONST::alpha * PHYSCONST::alpha;

        RankOneTensorOperator<3> &h = (*this);
        h[0] = (coef * alpha_2) * delta * s[0];
        h[1] = (coef * alpha_2) * delta * s[1];
        h[2] = (coef * alpha_2) * delta * s[2];
    }

protected:
    H_B_spin s;
    DeltaOperator delta;
};

} //namespace mrchem
