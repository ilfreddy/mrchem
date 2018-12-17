#pragma once

#include "qmoperators/two_electron/XCPotential.h"

/** @class XCPotential
 *
 * @brief Exchange-Correlation potential defined by a particular (spin) density
 *
 * The XC potential is computed by mapping of the density through a XC functional,
 * provided by the XCFun library. There are two ways of defining the density:
 *
 *  1) Use getDensity() prior to setup() and build the density as you like.
 *  2) Provide a default set of orbitals in the constructor that is used to
 *     compute the density on-the-fly in setup().
 *
 * If a set of orbitals has NOT been given in the constructor, the density
 * MUST be explicitly computed prior to setup(). The density will be computed
 * on-the-fly in setup() ONLY if it is not already available. After setup() the
 * operator will be fixed until clear(), which deletes both the density and the
 * potential.
 *
 * LDA and GGA functionals are supported as well as two different ways to compute
 * the XC potentials: either with explicit derivatives or gamma-type derivatives.
 */

namespace mrchem {

class XCPotentialD1 final : public XCPotential {
public:
    XCPotentialD1(mrdft::XCFunctional *F, OrbitalVector *Phi = nullptr);

private:
    mrcpp::FunctionTreeVector<3> potentials; ///< XC Potential functions collected in a vector

    void setup(double prec);
    void clear();

    void setupPotential(double prec);
    mrcpp::FunctionTree<3> &getPotential(int spin);

    Orbital apply(Orbital phi);
};

} //namespace mrchem
