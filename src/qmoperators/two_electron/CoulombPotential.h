#pragma once

#include "qmfunctions/Density.h"
#include "qmoperators/one_electron/QMPotential.h"

/** @class CoulombPotential
 *
 * @brief Coulomb potential defined by a particular electron density
 *
 * The Coulomb potential is computed by application of the Poisson operator
 * an electron density. There are two ways of defining the density:
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
 */

namespace mrchem {

class CoulombPotential : public QMPotential {
public:
    explicit CoulombPotential(std::shared_ptr<mrcpp::PoissonOperator> P,
                              std::shared_ptr<OrbitalVector> Phi = nullptr,
                              bool mpi_share = false);
    ~CoulombPotential() override = default;

    friend class CoulombOperator;

protected:
    Density density; ///< Ground-state electron density

    std::shared_ptr<OrbitalVector> orbitals; ///< Unperturbed orbitals defining the ground-state electron density
    std::shared_ptr<mrcpp::PoissonOperator> poisson; ///< Operator used to compute the potential

    auto &getPoisson() { return this->poisson; }
    auto &getDensity() { return this->density; }

    bool hasDensity() const { return (this->density.squaredNorm() < 0.0) ? false : true; }

    void setup(double prec) override;
    void clear() override;

    virtual void setupGlobalDensity(double prec) {}
    virtual void setupLocalDensity(double prec) {}

    void setupGlobalPotential(double prec);
    QMFunction setupLocalPotential(double prec);
    void allreducePotential(double prec, QMFunction &V_loc);
};

} // namespace mrchem
