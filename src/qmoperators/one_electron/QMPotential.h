#pragma once

#include "qmfunctions/QMFunction.h"
#include "qmoperators/QMOperator.h"

/** @class QMPotential
 *
 * @brief Operator defining a multiplicative potential
 *
 * Inherits the general features of a complex function from QMFunction and
 * implements the multiplication of this function with an Orbital. The actual
 * function representing the operator needs to be implemented in the derived
 * classes, where the *re and *im FunctionTree pointers should be assigned in
 * the setup() function and deallocated in the clear() function.
 *
 */

namespace mrchem {

class QMPotential : public QMFunction, public QMOperator {
public:
    explicit QMPotential(int adap, bool shared = false);
    QMPotential(const QMPotential &pot) = delete;
    QMPotential &operator=(const QMPotential &pot) = delete;
    ~QMPotential() override;

protected:
    int adap_build;

    ComplexDouble evalf(const mrcpp::Coord<3> &r) const override {
        ComplexDouble out(0.0, 0.0), i(0.0, 1.0);
        if (this->hasReal()) out += this->real().evalf(r);
        if (this->hasImag()) out += i * this->imag().evalf(r);
        return out;
    }

    Orbital apply(Orbital inp) override;
    Orbital dagger(Orbital inp) override;

    void calcRealPart(Orbital &out, Orbital &inp, bool dagger);
    void calcImagPart(Orbital &out, Orbital &inp, bool dagger);
};

} // namespace mrchem
