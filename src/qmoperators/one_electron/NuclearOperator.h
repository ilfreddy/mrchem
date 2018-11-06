#pragma once

#include "qmoperators/one_electron/QMPotential.h"
#include "qmoperators/RankZeroTensorOperator.h"
#include "analyticfunctions/NuclearFunction.h"

namespace mrchem {

class NuclearPotential final : public QMPotential {
public:
    NuclearPotential(const Nuclei &nucs, double prec);

    void setup(double prec);
    void clear();

    Nuclei &getNuclei() { return this->func.getNuclei(); }
    const Nuclei &getNuclei() const { return this->func.getNuclei(); }
    double evalf(const mrcpp::Coord<3> &r)  { return this->func.evalf(r); }
    
protected:
    NuclearFunction func;
};

class NuclearOperator final : public RankZeroTensorOperator {
public:
    NuclearOperator(const Nuclei &nucs, double prec)
            : r_m1(nucs, prec) {
        RankZeroTensorOperator &v = (*this);
        v = r_m1;
    }

    Nuclei &getNuclei() { return this->r_m1.getNuclei(); }
    const Nuclei &getNuclei() const { return this->r_m1.getNuclei(); }

    double trace(const Nuclei &nucs); 

    using RankZeroTensorOperator::trace;
    
 protected:
    NuclearPotential r_m1;
};

} //namespace mrchem
