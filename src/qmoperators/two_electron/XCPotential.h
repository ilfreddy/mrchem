#pragma once

#include "QMPotential.h"

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

namespace mrchem {
class XCFunctional;

class XCPotential final : public QMPotential {
public:
    XCPotential(XCFunctional &F, OrbitalVector &Phi, int k);
    ~XCPotential() { }

    void setup(double prec);
    void clear();

    double getEnergy() const { return energy; }
    
protected:
    int order;                                 ///< Order of kernel derivative
    double energy;                             ///< XC energy
    XCFunctional *functional;                  ///< External XC functional to be used
    OrbitalVector *orbitals;                   ///< External set of orbitals used to build the density
    mrcpp::FunctionTreeVector<3> potentials;   ///< XC Potential functions collected in a vector

    void evaluateXCFunctional();
    Orbital apply(Orbital phi);
};


} //namespace mrchem
