#pragma once

#include "QMPotential.h"

namespace mrchem {

/** 
 *  \class XCPotential
 *  \brief Exchange and Correlation potential
 *
 * Ideally, this should handle an arbitrary-order operator (k>1) for
 * response calculations Currently nly the potential is computed
 * (k=1). LDA and GGA functionals are allowed as well as two different ways
 * to compute the XC potentials: either with explicit derivatives or gamma-type derivatives
 * The calss handles the bookkeeping for the input/output of the xcfun library through arrays of
 * FunctionTree(s). 
 *
 *  \author Luca Frediani
 *  \date 2018
 *  
 */
class XCPotential final : public QMPotential {
public:
    XCPotential(XCFunctional &F, OrbitalVector &Phi, mrcpp::DerivativeOperator<3> &D, int k);
    ~XCPotential();

    void setup(double prec);
    void clear();

    double getEnergy() const { return this->energy; }
    
protected:
    int order;                           ///< Order of kernel derivative
    int nPotentials;                           ///< Number of potential energy functions
    XCFunctional *functional;                  ///< External XC functional to be used
    mrcpp::DerivativeOperator<3> *derivative;  ///< External derivative operator
    OrbitalVector *orbitals;                   ///< External set of orbitals used to build the density
    double energy;                             ///< XC energy
    mrcpp::FunctionTreeVector<3> potentials;   ///< XC Potential functions collected in a vector

    void evaluateXCFunctional();
    int getPotentialFunctionIndex(const Orbital & orb);
    Orbital apply (Orbital phi);
};


} //namespace mrchem
