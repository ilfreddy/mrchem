#pragma once

#include "qmoperators/QMOperator.h"
#include "qmoperators/RankZeroTensorOperator.h"

namespace mrchem {

class QMIdentity final : public QMOperator {
public:
    QMIdentity()
            : QMOperator() {}

private:
    Orbital apply(Orbital inp) override;
    Orbital dagger(Orbital inp) override;
};

class IdentityOperator final : public RankZeroTensorOperator {
public:
    IdentityOperator() {
        I = std::make_shared<QMIdentity>();

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &h = (*this);
        h = I;
        h.name() = "I";
    }

    ComplexDouble operator()(Orbital bra, Orbital ket);
    ComplexDouble dagger(Orbital bra, Orbital ket);

    ComplexMatrix operator()(OrbitalVector &bra, OrbitalVector &ket);
    ComplexMatrix dagger(OrbitalVector &bra, OrbitalVector &ket);

    // Necessary in order to pick up base class definitions for overloaded functions
    using RankZeroTensorOperator::operator();
    using RankZeroTensorOperator::dagger;

private:
    std::shared_ptr<QMIdentity> I;
};

} // namespace mrchem
