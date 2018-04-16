#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "FunctionTree.h"
#include "FunctionTreeVector.h"
#include "xcfun.h"

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
class XCFunctional {
public:
    XCFunctional(bool s, bool e, double thrs = 0.0);
    virtual ~XCFunctional();

    void setDensityCutoff(double thrs) { this->cutoff = thrs; }
    void setFunctional(const std::string &name, double coef = 1.0);

    int getInputLength() const { return xc_input_length(this->functional); }
    int getOutputLength() const { return xc_output_length(this->functional); }

    bool isLDA() const { return (!(this->isGGA() || this->isMetaGGA())); }
    bool isGGA() const { return (xc_is_gga(this->functional)); }
    bool isMetaGGA() const { return (xc_is_metagga(this->functional)); }
    
    bool isSpinSeparated() const { return this->spin; }
    bool needsGamma() const { return (expDerivatives == 0);};

    void evaluate(int k, Eigen::MatrixXd &inp, Eigen::MatrixXd &out) const;
    void evalSetup(const int order);

    FunctionTree<3> * calcPotentialGGA(FunctionTree<3> & df_drho, FunctionTree<3> & df_dgamma,
                                       FunctionTreeVector<3> grad_rho, DerivativeOperator<3> *derivative,
                                       int maxScale);
    FunctionTree<3> * calcPotentialGGA(FunctionTree<3> & df_drhoa, FunctionTree<3> & df_dgaa,
                                       FunctionTree<3> & df_dgab, FunctionTreeVector<3> grad_rhoa,
                                       FunctionTreeVector<3> grad_rhob, DerivativeOperator<3> *derivative,
                                       int maxScale);
    FunctionTree<3> * calcPotentialGGA(FunctionTree<3> & df_drho, FunctionTreeVector<3> & df_dgr,
                                       DerivativeOperator<3> *derivative, int maxScale);
 protected:
    FunctionTree<3> * addPotentialContributions(FunctionTreeVector<3> & contributions,
                                                int maxScale);
    FunctionTree<3> * calcDivergence(FunctionTreeVector<3> &inp,
                                    DerivativeOperator<3> *derivative,
                                    int maxScale);
    FunctionTree<3> * calcGradDotPotDensVec(FunctionTree<3> &V,
                                           FunctionTreeVector<3> &rho,
                                           DerivativeOperator<3> *derivative,
                                           int maxScale);
    DerivativeOperator<3> *derivative;  ///< External derivative operator
    double energy;                      ///< XC energy
    Density density;                    ///< Unperturbed density

    FunctionTree<3> **xcInput;          ///< Bookkeeping array to feed XCFun
    FunctionTree<3> **xcOutput;         ///< Bookkeeping array returned by XCFun

    FunctionTreeVector<3> potentialFunction; ///< Storage of the computed potential functions
    FunctionTreeVector<3> grad_a; ///< Gradient of the alpha density        
    FunctionTreeVector<3> grad_b; ///< Gradient of the beta  density        
    FunctionTreeVector<3> grad_t; ///< Gradient of the total density        
    FunctionTreeVector<3> gamma;  ///< Gamma function(s) (three fcns for spin separated calculations)       

    void evaluate(OrbitalVector * orbitals);
    void setupXCInput();
    void setupXCOutput();
    int setupXCInputDensity(int nUsed);
    int setupXCInputGradient(int nUsed);
    void clearXCInput();
    void clearXCOutput();
    void calcDensity(OrbitalVector * orbitals);
    int calcDensityGradient();
    void calcGamma();

    FunctionTreeVector<3> calcPotential();
    bool cropPotential(double prec);

    void calcPotentialLDA();
    
    void calcPotentialGGA();

    void calcEnergy();
    void evaluateXCFunctional();

    void compressTreeData(int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);
    void expandTreeData(int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);

    void compressNodeData(int n, int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);
    void expandNodeData(int n, int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);

    FunctionTreeVector<3> calcGradient(FunctionTree<3> &inp);
    FunctionTree<3>* calcDotProduct(FunctionTreeVector<3> &vec_a, FunctionTreeVector<3> &vec_b);

    int getPotentialFunctionIndex(const Orbital & orb);
    
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
        if (ptr == 0) MSG_FATAL("Clearing NULL pointer");
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

};

