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
 XCOperator(mrcpp::XCFunctional &F, OrbitalVector &Phi, mrcpp::DerivativeOperator D, int k)
     : potential(0) {
        this->xcPotential = new XCPotential(F, Phi, D, k);

        RankZeroTensorOperator &XC = (*this);
        XC = *xcPotential;
    }
    ~XCOperator() { delete this->xcPotential; }
    
 protected:
    XCPotential *xcPotential;
};
    
} //namespace mrchem
