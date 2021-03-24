#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "NablaOperator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

#include "qmoperators/RankOneTensorOperator.h"
#include "qmoperators/one_electron/QMPotential.h"

using mrcpp::DerivativeOperator;
using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

QMNabla::QMNabla(int d, std::shared_ptr<DerivativeOperator<3>> D)
        : QMOperator()
        , apply_dir(d)
        , derivative(D) {}

Orbital QMNabla::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    if (this->derivative == nullptr) MSG_ERROR("No derivative operator");

    auto dir = this->apply_dir;
    auto &D = *this->derivative;

    Orbital out = inp.paramCopy();

    // Calc real part
    if (inp.hasReal()) {
        out.alloc(NUMBER::Real);
        mrcpp::apply(out.real(), D, inp.real(), dir);
    }
    // Calc imag part
    if (inp.hasImag()) {
        out.alloc(NUMBER::Imag);
        mrcpp::apply(out.imag(), D, inp.imag(), dir);
        if (inp.conjugate()) out.imag().rescale(-1.0);
    }

    return out;
}

/* Overload operator() to allow for an intuitive interface to
   compute the gradient of a QMPotential.  */
RankOneTensorOperator<3> NablaOperator::operator()(QMPotential &V) {
    auto dx = std::make_shared<QMPotential>(1);
    auto dy = std::make_shared<QMPotential>(1);
    auto dz = std::make_shared<QMPotential>(1);

    if (V.hasReal()) {
        mrcpp::FunctionTreeVector<3> dV = mrcpp::gradient(*(this->derivativeOperator), V.real());
        dx->setReal(&mrcpp::get_func(dV, 0));
        dy->setReal(&mrcpp::get_func(dV, 1));
        dz->setReal(&mrcpp::get_func(dV, 2));
    }

    if (V.hasImag()) {
        mrcpp::FunctionTreeVector<3> dV = mrcpp::gradient(*(this->derivativeOperator), V.imag());
        dx->setImag(&mrcpp::get_func(dV, 0));
        dy->setImag(&mrcpp::get_func(dV, 1));
        dz->setImag(&mrcpp::get_func(dV, 2));
    }

    RankOneTensorOperator<3> out;
    out[0] = dx;
    out[1] = dy;
    out[2] = dz;

    return out;
}

Orbital QMNabla::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
