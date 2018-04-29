#include "mrchem.h"
#include "constants.h"
#include "Orbital.h"
#include "qmfunctions.h"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "XCFunctional.h"
#include "XCPotential.h"

using namespace mrcpp;
using namespace std;
using namespace Eigen;

extern MultiResolutionAnalysis<3> *MRA;

namespace mrchem {

/** @brief constructor
 *
 * @param[in] k order of the operator
 * @param[in] F XCFunctional pointer
 * @param[in] Phi vector of orbitals
 *
 * Based on the order and spin the correct nr. of potential functions is determined
 * Then the functional is set up for subsequent calculations, fixing some internals of
 * xcfun when F.evalSetup is invoked.
 *
 */

XCPotential::XCPotential(XCFunctional &F, OrbitalVector &Phi, int k) 
    : QMPotential(1), //HACK is the correct value of adap here?
      orbitals(&Phi),
      functional(&F),
      energy(0.0),
      order(k) {};
    
/** @brief destructor
 *
 */
XCPotential::~XCPotential() {
    this->functional = 0;
    this->orbitals = 0;
}

/** @brief setup the XCPotential
 * 
 * @param[in] prec precision 
 *
 * Sequence of steps required to compute the XC potentials The moat
 * important steps are evaluateXCFunctional and calcPotential where
 * the functional derivatives are computed and the potential assembled
 * respectively
 *
 */
void XCPotential::setup(double prec) {
    setApplyPrec(prec);
    evaluateXCFunctional();
}

/** @brief XCPotential application
 *
 * The operator is applied by choosing the correct potential function
 * which is then assigned to the real function part of the operator
 * base-class before the base class function is called.
 *
 * @param[in] phi orbital to which the potential is applied.
 *
 */
Orbital XCPotential::apply(Orbital phi) {
    int potentialIndex =  this->functional->getPotentialFunctionIndex(phi);
    FunctionTree<3> * potential = this->functional->getPotentialFunction(potentialIndex);
    this->setReal(potential);
    this->setImag(0);
    Orbital Vphi = QMPotential::apply(phi); 
    this->setReal(0);
    return Vphi;
}

/** @brief clears all data in the XCPotential object
 *
 */
void XCPotential::clear() {
    energy = 0.0;
    nPotentials = -1;
    order = -1;
    clearApplyPrec();
}

/** @brief evaluation of the functional and its derivatives
 *
 */
void XCPotential::evaluateXCFunctional() {
    this->functional->setup(this->order);
}

} //namespace mrchem
