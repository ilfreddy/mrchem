#ifndef FUNCTIONTREEVECTOR_H
#define FUNCTIONTREEVECTOR_H

#include <vector>

#include "FunctionTree.h"

template<int D>
class FunctionTreeVector {
public:
    FunctionTreeVector() { }
    virtual ~FunctionTreeVector() { }

    FunctionTreeVector(const FunctionTreeVector<D> &vec) {
        this->clear();
        for (int i = 0; i < vec.size(); i++) {
            this->coefs.push_back(vec.coefs[i]);
            this->funcs.push_back(vec.funcs[i]);
        }
    }
    FunctionTreeVector& operator=(const FunctionTreeVector<D> &vec) {
        this->clear();
        for (int i = 0; i < vec.size(); i++) {
            this->coefs.push_back(vec.coefs[i]);
            this->funcs.push_back(vec.funcs[i]);
        }
    }

    int size() const {
        return this->funcs.size();
    }
    void push_back(double c, FunctionTree<D> &f) {
        this->coefs.push_back(c);
        this->funcs.push_back(&f);
    }
    void push_back(FunctionTree<D> &f) {
        this->coefs.push_back(1.0);
        this->funcs.push_back(&f);
    }
    void clear() {
        this->coefs.clear();
        this->funcs.clear();
    }
    double getCoef(int i) const {
        if (i < 0 or i >= this->coefs.size()) MSG_ERROR("Out of bounds");
        return this->coefs[i];
    }
    FunctionTree<D> &getFunc(int i) {
        if (i < 0 or i >= this->funcs.size()) MSG_ERROR("Out of bounds");
        return *this->funcs[i];
    }
    const FunctionTree<D> &getFunc(int i) const {
        if (i < 0 or i >= this->funcs.size()) MSG_ERROR("Out of bounds");
        return *this->funcs[i];
    }
protected:
    std::vector<double> coefs;
    std::vector<FunctionTree<D> *> funcs;
};

#endif // FUNCTIONTREEVECTOR_H
