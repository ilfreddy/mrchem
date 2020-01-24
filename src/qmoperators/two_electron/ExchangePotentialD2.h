#pragma once

#include <memory>

#include "ExchangePotential.h"
#include "qmfunctions/qmfunction_fwd.h"
#include "qmoperators/QMOperator.h"

namespace mrchem {

/** @class ExchangePotentialD2
 *
 *  @brief Hartree-Fock exchange potential defined by a set of unperturbed orbitals
 *
 * The operator is defined as the Hartree-Fock exchange arising from a
 * set of unperturbed orbitals. The OrbitalVector defining the
 * operator is fixed throughout the operator life time, but the
 * orbitals themselves are allowed to change in between each
 * application. The internal exchange potentials (the operator applied
 * to it's own orbitals) can be precomputed and stored for fast
 * retrieval. Option to use screening based on previous calculations
 * of the internal exchange (make sure that the internal orbitals
 * haven't been significantly changed since the last time the operator
 * was set up, e.g. through an orbital rotation).
 */

class ExchangePotentialD2 final : public ExchangePotential {
public:
    ExchangePotentialD2(std::shared_ptr<mrcpp::PoissonOperator> P,
                        std::shared_ptr<OrbitalVector> Phi,
                        std::shared_ptr<OrbitalVector> X,
                        std::shared_ptr<OrbitalVector> Y,
                        bool s = false);
    ~ExchangePotentialD2() override = default;

    friend class ExchangeOperator;

private:
    std::shared_ptr<OrbitalVector> orbitals_x; ///< first set of perturbed orbitals defining the exchange operator
    std::shared_ptr<OrbitalVector> orbitals_y; ///< second set of perturbed orbitals defining the exchange operator
    bool useOnlyX;                             ///< true if X and Y are the same set of orbitals

    Orbital calcExchange(Orbital phi_p);
    void calcInternal(int i);
    void calcInternal(int i, int j);

    Orbital calcExchange_X(Orbital phi_p);
    void calcInternal_X(int i);
    void calcInternal_X(int i, int j);

    Orbital calcExchange_XY(Orbital phi_p);
    void calcInternal_XY(int i);
    void calcInternal_XY(int i, int j);

    Orbital calcExchangeComponent(Orbital phi_p, Orbital phi_i, Orbital x_i);
};

} // namespace mrchem
