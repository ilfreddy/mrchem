#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombPotentialD2.h"
#include "qmfunctions/density_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

CoulombPotentialD2::CoulombPotentialD2(mrcpp::PoissonOperator *P,
                                       OrbitalVector *Phi,
                                       OrbitalVector *X,
                                       OrbitalVector *Y)
        : CoulombPotential(P),
          orbitals(Phi),
          orbitals_x(X),
          orbitals_y(Y) {
}

void CoulombPotentialD2::setupDensity(double prec) {
    if (hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_x == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_y == nullptr) MSG_ERROR("Orbitals not initialized");

    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &X = *this->orbitals_x;
    OrbitalVector &Y = *this->orbitals_y;
    Density &rho = this->density;

    Timer timer;
    density::compute(prec, rho, Phi, X, Y, DENSITY::Total);
    timer.stop();
    double t = timer.getWallTime();
    int n = rho.getNNodes();
    Printer::printTree(0, "Coulomb density", n, t);
}

} //namespace mrchem
