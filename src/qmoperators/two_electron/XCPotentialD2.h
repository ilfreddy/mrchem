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

#include "MRCPP/MWFunctions"
namespace mrchem {

class XCPotentialD2 final : public XCPotential {
public:
    XCPotentialD2(std::shared_ptr<mrdft::XCFunctional> F,
                  std::shared_ptr<OrbitalVector> Phi,
                  std::shared_ptr<OrbitalVector> X,
                  std::shared_ptr<OrbitalVector> Y,
                  bool mpi_shared = false);
    ~XCPotentialD2() override;

private:
    std::shared_ptr<OrbitalVector> *orbitals_x;               ///< 1st external set of perturbed orbitals used to build the density
    std::shared_ptr<OrbitalVector> *orbitals_y;               ///< 2nd external set of perturbed orbitals used to build the density

    void setup(double prec) override;
    void setupPotential(double prec) override;
    void buildPerturbedDensity(OrbitalVector &Phi,
                               OrbitalVector &X,
                               OrbitalVector &Y,
                               DENSITY::DensityType spin);

    void setupPerturbedDensity(double prec = -1.0);
};

} // namespace mrchem
