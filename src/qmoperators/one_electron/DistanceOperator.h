#pragma once

#include "qmoperators/one_electron/QMPotential.h"
#include "qmoperators/RankZeroTensorOperator.h"
#include "analyticfunctions/NuclearFunction.h"

namespace mrchem {

class DistancePotential final : public QMPotential {
public:
    DistancePotential(double pow, const mrcpp::Coord<3> &r_K, double S);

    void setup(double prec);
    void clear();

protected:
    const double power;
    NuclearFunction func;
};

class DistanceOperator final : public RankZeroTensorOperator {
public:
    DistanceOperator(double pow, const mrcpp::Coord<3> &R, double S = 1.0e-7)
            : r_pow(pow, R, S) {
        RankZeroTensorOperator &v = (*this);
        v = r_pow;
    }

protected:
    DistancePotential r_pow;
};


} //namespace mrchem
