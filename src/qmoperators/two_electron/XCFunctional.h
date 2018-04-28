#pragma once

#include <Eigen/Core>

#include "MRCPP/MWFunctions"
#include "MRCPP/Printer"
#include "xcfun.h"
#include "Density.h"

/** 
 *  \class XCFunctional
 *  \brief Compute XC functional with XCFun
 *
 *  Interface class for the XCFun library
 *
 *  \author Stig Rune Jensen
 *  \date 2015
 *  
 */

        
namespace mrchem {

class XCFunctional {
public:
    XCFunctional(bool s, bool e, double thrs, OrbitalVector &phi, mrcpp::DerivativeOperator<3> *D);
    virtual ~XCFunctional();

    void setDensityCutoff(double thrs) { this->cutoff = thrs; }
    void setFunctional(const std::string &name, double coef = 1.0);

    int getInputLength() const { return xc_input_length(this->functional); }
    int getOutputLength() const { return xc_output_length(this->functional); }
    double getEnergy() const { return energy; }

    bool isLDA() const { return (!(this->isGGA() || this->isMetaGGA())); }
    bool isGGA() const { return (xc_is_gga(this->functional)); }
    bool isMetaGGA() const { return (xc_is_metagga(this->functional)); }
    
    bool isSpinSeparated() const { return this->spin; }
    bool needsGamma() const { return (expDerivatives == 0);};

    void evaluate(int k, DoubleMatrix &inp, DoubleMatrix &out) const;
    void setup(const int order);
    void evalSetup(const int order);

    int getPotentialFunctionIndex(const Orbital & orb);
    mrcpp::FunctionTree<3>* getPotentialFunction(int index) {return potentialFunction[index];};

 protected:
    Density density;                                ///< Unperturbed density
    mrcpp::FunctionTree<3> **xcInput;               ///< Bookkeeping array to feed XCFun
    mrcpp::FunctionTree<3> **xcOutput;              ///< Bookkeeping array returned by XCFun
    mrcpp::FunctionTreeVector<3> potentialFunction; ///< Storage of the computed potential functions
    mrcpp::FunctionTreeVector<3> grad_a;            ///< Gradient of the alpha density        
    mrcpp::FunctionTreeVector<3> grad_b;            ///< Gradient of the beta  density        
    mrcpp::FunctionTreeVector<3> grad_t;            ///< Gradient of the total density        
    mrcpp::FunctionTreeVector<3> gamma;             ///< Gamma function(s) (three fcns for spin separated calculations)       
    mrcpp::DerivativeOperator<3> *derivative;       ///< External derivative operator

    mrcpp::FunctionTree<3> * calcPotentialGGA(mrcpp::FunctionTree<3> & df_drho,
                                              mrcpp::FunctionTree<3> & df_dgamma,
                                              mrcpp::FunctionTreeVector<3> grad_rho,
                                              mrcpp::DerivativeOperator<3> *derivative);
    mrcpp::FunctionTree<3> * calcPotentialGGA(mrcpp::FunctionTree<3> & df_drhoa,
                                              mrcpp::FunctionTree<3> & df_dgaa,
                                              mrcpp::FunctionTree<3> & df_dgab,
                                              mrcpp::FunctionTreeVector<3> grad_rhoa,
                                              mrcpp::FunctionTreeVector<3> grad_rhob,
                                              mrcpp::DerivativeOperator<3> *derivative);
    mrcpp::FunctionTree<3> * calcPotentialGGA(mrcpp::FunctionTree<3> & df_drho,
                                              mrcpp::FunctionTreeVector<3> & df_dgr,
                                              mrcpp::DerivativeOperator<3> *derivative);
    mrcpp::FunctionTree<3> * addPotentialContributions(mrcpp::FunctionTreeVector<3> & contributions);
    mrcpp::FunctionTree<3> * calcDivergence(mrcpp::FunctionTreeVector<3> &inp,
                                            mrcpp::DerivativeOperator<3> *derivative);
    mrcpp::FunctionTree<3> * calcGradDotPotDensVec(mrcpp::FunctionTree<3> &V,
                                                   mrcpp::FunctionTreeVector<3> &rho,
                                                   mrcpp::DerivativeOperator<3> *derivative);
    
    void evaluate(OrbitalVector * orbitals);
    void setupXCInput();
    void setupXCOutput();
    int setupXCInputDensity(int nUsed);
    int setupXCInputGradient(int nUsed);
    void clearXCInput();
    void clearXCOutput();
    int calcDensityGradient();
    void calcGamma();
    mrcpp::FunctionTreeVector<3> calcPotential();
    bool cropPotential(double prec);
    void calcPotentialLDA();
    void calcPotentialGGA();
    void calcEnergy();
    void evaluateFunctional();

    void compressNodeData(int n, int nFuncs, mrcpp::FunctionTree<3> **trees, Eigen::MatrixXd &data);
    void expandNodeData(int n, int nFuncs, mrcpp::FunctionTree<3> **trees, Eigen::MatrixXd &data);

    mrcpp::FunctionTreeVector<3> calcGradient(mrcpp::FunctionTree<3> &inp);
    mrcpp::FunctionTree<3>* calcDotProduct(mrcpp::FunctionTreeVector<3> &vec_a, mrcpp::FunctionTreeVector<3> &vec_b);


    template<class T>
    int sumNodes(T **trees, int nTrees) {
        int nNodes = 0;
        for (int i = 0; i < nTrees; i++) {
            if (trees[i] != 0) {
                nNodes += trees[i]->getNNodes();
            }
        }
        return nNodes;
    }

    template<class T>
    T** allocPtrArray(int n_funcs) {
        T **ptr = new T*[n_funcs];
        for (int i = 0; i < n_funcs; i++) {
            ptr[i] = 0;
        }
        return ptr;
    }

    template<class T>
    void clearPtrArray(int n_funcs, T **ptr) {
        if (ptr == 0) MSG_ERROR("Clearing NULL pointer");
        for (int i = 0; i < n_funcs; i++) {
            if (ptr[i] != 0) {
                delete ptr[i];
            }
            ptr[i] = 0;
        }
    }

    template<class T>
    T** deletePtrArray(int n_funcs, T ***ptr) {
        if (*ptr != 0) {
            for (int i = 0; i < n_funcs; i++) {
                if ((*ptr)[i] != 0) {
                    delete (*ptr)[i];
                }
                (*ptr)[i] = 0;
            }
            delete[] *ptr;
        }
        return 0;
    }

private:
    bool spin;                  ///< Spin polarization
    unsigned int expDerivatives;///< whether gamma-type or explicit derivatives are used
    double cutoff;              ///< Below the cutoff value, the density will be considered zero
    xc_functional functional;   ///< The functional in the XCFun library (struct from xcfun library)
    OrbitalVector * orbitals;   ///< Set of orbitals used to compute the density defining the functional
    double energy;              ///< XC energy
};
 
} //namespace mrchem
