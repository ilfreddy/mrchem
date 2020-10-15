#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "GradZoraOperator.h"
#include "ZoraPotential.h"

#include "utils/print_utils.h"

using mrcpp::DerivativeOperator;
using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

void GradZoraOperator::computeGradComponent(QMFunction component,
                                            int dir,
                                            std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                                            ZoraPotential Z) {

    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    if (this->derivative == nullptr) MSG_ERROR("No derivative operator");
    if (zora.hasReal()) {}
    // Calc real part
    if (zora.hasReal()) {
        component.alloc(NUMBER::Real);
        mrcpp::apply(component.real(), D, Z.real(), dir);
    }
    // Calc imag part
    if (zora.hasImag()) {
        component.alloc(NUMBER::Imag);
        mrcpp::apply(component.imag(), D, Z.imag(), dir);
        if (zora.conjugate()) component.imag().rescale(-1.0);
    }
} // namespace mrchem
