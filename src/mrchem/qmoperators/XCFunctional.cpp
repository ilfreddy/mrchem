#include "XCFunctional.h"
#include "TelePrompter.h"

using namespace std;
using namespace Eigen;

XCFunctional::XCFunctional(bool s) : spin(s) {
    this->functional = xc_new_functional();
    if (this->spin) {
        xc_set_mode(this->functional, XC_VARS_AB);
    } else {
        xc_set_mode(this->functional, XC_VARS_N);
    }
}

XCFunctional::~XCFunctional() {
    xc_free_functional(this->functional);
}

void XCFunctional::setFunctional(const string &name, double coef) {
    int param = getParamFromName(name);
    xc_set_param(this->functional, param, coef);
}

int XCFunctional::getParamFromName(const string &name) {
    int param = -1;
    if (name == "SLATERX") {
        MSG_WARN("Functional not tested");
        param = XC_SLATERX;
    } else if (name == "VWN5C") {
        MSG_WARN("Functional not tested");
        param = XC_VWN5C;
    } else if (name == "BECKEX") {
        MSG_WARN("Functional not tested");
        param = XC_BECKEX;
    } else if (name == "BECKECORRX") {
        MSG_WARN("Functional not tested");
        param = XC_BECKECORRX;
    } else if (name == "BECKESRX") {
        MSG_WARN("Functional not tested");
        param = XC_BECKESRX;
    } else if (name == "OPTX") {
        MSG_WARN("Functional not tested");
        param = XC_OPTX;
    } else if (name == "LYPC") {
        MSG_WARN("Functional not tested");
        param = XC_LYPC;
    } else if (name == "PBEX") {
        MSG_WARN("Functional not tested");
        param = XC_PBEX;
    } else if (name == "REVPBEX") {
        MSG_WARN("Functional not tested");
        param = XC_REVPBEX;
    } else if (name == "RPBEX") {
        MSG_WARN("Functional not tested");
        param = XC_RPBEX;
    } else if (name == "PBEC") {
        MSG_WARN("Functional not tested");
        param = XC_PBEC;
    } else if (name == "SPBEC") {
        MSG_WARN("Functional not tested");
        param = XC_SPBEC;
    } else if (name == "VWN_PBEC") {
        MSG_WARN("Functional not tested");
        param = XC_VWN_PBEC;
    } else if (name == "LDAERFX") {
        MSG_WARN("Functional not tested");
        param = XC_LDAERFX;
    } else if (name == "LDAERFC") {
        MSG_WARN("Functional not tested");
        param = XC_LDAERFC;
    } else if (name == "LDAERFC_JT") {
        MSG_WARN("Functional not tested");
        param = XC_LDAERFC_JT;
    } else if (name == "KTX") {
        MSG_WARN("Functional not tested");
        param = XC_KTX;
    } else if (name == "TFK") {
        MSG_WARN("Functional not tested");
        param = XC_TFK;
    } else if (name == "PW91X") {
        MSG_WARN("Functional not tested");
        param = XC_PW91X;
    } else if (name == "PW91K") {
        MSG_WARN("Functional not tested");
        param = XC_PW91K;
    } else if (name == "PW92C") {
        MSG_WARN("Functional not tested");
        param = XC_PW92C;
    } else if (name == "MO5X") {
        MSG_WARN("Functional not tested");
        param = XC_M05X;
    } else if (name == "MO5X2X") {
        MSG_WARN("Functional not tested");
        param = XC_M05X2X;
    } else if (name == "MO6X") {
        MSG_WARN("Functional not tested");
        param = XC_M06X;
    } else if (name == "MO6X2X") {
        MSG_WARN("Functional not tested");
        param = XC_M06X2X;
    } else if (name == "MO6LX") {
        MSG_WARN("Functional not tested");
        param = XC_M06LX;
    } else if (name == "MO6HFX") {
        MSG_WARN("Functional not tested");
        param = XC_M06HFX;
    } else if (name == "BRX") {
        MSG_WARN("Functional not tested");
        param = XC_BRX;
    } else if (name == "MO5X2C") {
        MSG_WARN("Functional not tested");
        param = XC_M05X2C;
    } else if (name == "MO5C") {
        MSG_WARN("Functional not tested");
        param = XC_M05C;
    } else if (name == "MO6C") {
        MSG_WARN("Functional not tested");
        param = XC_M06C;
    } else if (name == "MO6HFC") {
        MSG_WARN("Functional not tested");
        param = XC_M06HFC;
    } else if (name == "MO6LC") {
        MSG_WARN("Functional not tested");
        param = XC_M06LC;
    } else if (name == "MO6X2C") {
        MSG_WARN("Functional not tested");
        param = XC_M06X2C;
    } else if (name == "TPSSC") {
        MSG_WARN("Functional not tested");
        param = XC_TPSSC;
    } else if (name == "TPSSX") {
        MSG_WARN("Functional not tested");
        param = XC_TPSSX;
    } else if (name == "REVTPSSC") {
        MSG_WARN("Functional not tested");
        param = XC_REVTPSSC;
    } else if (name == "REVTPSSX") {
        MSG_WARN("Functional not tested");
        param = XC_REVTPSSX;
    } else {
        MSG_ERROR("Invalid functional");
    }
    return param;
}

/** Computes the alpha and beta exchange-correlation potentials
 * from the xcfun output functions. For LDA's these are the second
 * and third output functions, respectively. For GGA's the potentials
 * must be computed through
 * \f$ v_{xc}^\sigma = \frac{\partial F_{xc}}{\partial \rho^\sigma(r)}
 *  - \nabla\cdot\frac{\partial F_{xc}}{\partial(\nabla\rho^\sigma)} \f$
 *
 * XCFunctional output:
 *
 * LDA: \f$ \left(F_{xc}, \frac{\partial F_{xc}}{\partial \rho}\right) \f$
 *
 * GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho},
 *  \frac{\partial F_{xc}}{\partial \rho_x},
 *  \frac{\partial F_{xc}}{\partial \rho_y},
 *  \frac{\partial F_{xc}}{\partial \rho_z}\right) \f$
 *
 * Spin LDA: \f$ \left(F_{xc}, \frac{\partial F_{xc}}{\partial \rho^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho^\beta}\right) \f$
 *
 * Spin GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_x^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_y^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_z^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_x^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_y^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_z^\beta}
 *  \right) \f$
 */

void XCFunctional::evaluate(int k, MatrixXd &inp, MatrixXd &out) const {
    if (inp.cols() != getInputLength()) MSG_ERROR("Invalid input");

    int nPts = inp.rows();
    out = MatrixXd::Zero(nPts, getOutputLength(k));

    for (int i = 0; i < nPts; i++) {
        xc_eval(this->functional, k, inp.row(i).data(), out.row(i).data());
    }
}
