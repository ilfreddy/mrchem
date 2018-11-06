#pragma once

#include "QMPotential.h"
#include "qmoperators/RankZeroTensorOperator.h"

namespace mrchem {

class QMDelta final : public QMPotential {
public:
    QMDelta(const mrcpp::Coord<3> &o, double expo);

protected:
    mrcpp::GaussFunc<3> func;

    void setup(double prec);
    void clear();
};

class DeltaOperator final : public RankZeroTensorOperator {
public:
    DeltaOperator(const mrcpp::Coord<3> &o, double expo = 1.0e6)
            : delta(o, expo) {
        RankZeroTensorOperator &d = (*this);
        d = delta;
    }

protected:
    QMDelta delta;
};

} //namespace mrchem
