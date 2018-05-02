#pragma once

#include "RankZeroTensorOperator.h"
#include "XCPotential.h"

/**
 * \class XCOperator
 * \brief Exchange and Correlation operator
 *
 * Notes for myself: I need a way to "select" the correct potential to apply
 */
namespace mrchem {
    
class XCOperator final : public RankZeroTensorOperator {
 public:
 XCOperator(XCFunctional &F, OrbitalVector &Phi, int k)
     : xcPotential(0) {
        this->xcPotential = new XCPotential(F, Phi, k);

        RankZeroTensorOperator &XC = (*this);
        XC = *xcPotential;
    }
    ~XCOperator() { delete this->xcPotential; }

    double getEnergy() {return xcPotential->getEnergy();}
 protected:
    XCPotential *xcPotential;
};
    
} //namespace mrchem
