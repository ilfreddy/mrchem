#pragma once

#include <memory>

#include "qmfunctions/qmfunction_fwd.h"
#include "qmoperators/QMOperator.h"

namespace mrchem {

/** @class ExchangePotential
 *
 *  @brief Hartree-Fock exchange potential defined by a particular set of orbitals
 *
 * The operator is defined as the Hartree-Fock exchange arising from a particular
 * set of orbitals. The OrbitalVector defining the operator is fixed throughout the
 * operators life time, but the orbitals themselves are allowed to change in between
 * each application. The internal exchange potentials (the operator applied to it's
 * own orbitals) can be precomputed and stored for fast retrieval. Option to use
 * screening based on previous calculations of the internal exchange (make sure that
 * the internal orbitals haven't been significantly changed since the last time the
 * operator was set up, e.g. through an orbital rotation).
 */

class ExchangePotential : public QMOperator {
public:
    ExchangePotential(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi, bool s = false);
    ~ExchangePotential() override = default;

    friend class ExchangeOperator;

protected:
    bool screen;             ///< Apply screening in exchange evaluation
    DoubleVector tot_norms;  ///< Total norms for use in screening
    DoubleMatrix part_norms; ///< Partial norms for use in screening
    OrbitalVector exchange;  ///< Precomputed exchange orbitals from the occupied orbital set

    std::shared_ptr<OrbitalVector> orbitals;         ///< Orbitals defining the exchange operator
    std::shared_ptr<mrcpp::PoissonOperator> poisson; ///< Poisson operator to compute orbital contributions

    auto &getPoisson() { return this->poisson; }
    double getSpinFactor(Orbital phi_i, Orbital phi_j) const;

    void rotate(const ComplexMatrix &U);
    void setup(double prec) override;
    void clear() override;

    Orbital apply(Orbital phi_p) override;
    Orbital dagger(Orbital phi_p) override;

    using QMOperator::apply;
    using QMOperator::dagger;

    double getScaledPrecision(int i, int j) const;

    virtual Orbital calcExchange(Orbital phi_p) = 0;
    virtual void calcInternal(int i) = 0;
    virtual void calcInternal(int i, int j) = 0;
    virtual int testPreComputed(Orbital phi_p) const = 0;
    virtual void setupInternal(double prec) = 0;
};

} // namespace mrchem
