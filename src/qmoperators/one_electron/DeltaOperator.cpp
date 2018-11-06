#include "MRCPP/Gaussians"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "DeltaOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

QMDelta::QMDelta(const mrcpp::Coord<3> &o, double expo)
        : QMPotential(1) {
    // Delta function should integrate to one
    double coef = pow(expo/MATHCONST::pi, 3.0/2.0);

    this->func.setCoef(coef);
    this->func.setExp(expo);
    this->func.setPos(o.data());
}

void QMDelta::setup(double prec) {
    if (this->isSetup(prec)) return;
    this->setApplyPrec(prec);

    if (this->hasReal()) MSG_ERROR("Potential not properly cleared");
    if (this->hasImag()) MSG_ERROR("Potential not properly cleared");

    Timer timer;
    this->alloc(NUMBER::Real);
    mrcpp::build_grid(this->real(), this->func);
    mrcpp::project(this->apply_prec, this->real(), this->func);
    timer.stop();

    int n = this->getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(1, "Position operator", n, t);
}

void QMDelta::clear() {
    this->free();
    this->clearApplyPrec();
}

} //namespace mrchem
