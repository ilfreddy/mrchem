#pragma once

#include "qmoperators/QMOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class QMMomentum final : public QMOperator {
public:
    QMMomentum(int d, std::shared_ptr<mrcpp::DerivativeOperator<3>> D);

private:
    const int apply_dir;
    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative;

    Orbital apply(Orbital inp) override;
    Orbital dagger(Orbital inp) override;
};

class MomentumOperator final : public RankOneTensorOperator<3> {
public:
    MomentumOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D) {
        p_x = std::make_shared<QMMomentum>(0, D);
        p_y = std::make_shared<QMMomentum>(1, D);
        p_z = std::make_shared<QMMomentum>(2, D);

        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &p = (*this);
        p[0] = p_x;
        p[1] = p_y;
        p[2] = p_z;
        p[0].name() = "p[x]";
        p[1].name() = "p[y]";
        p[2].name() = "p[z]";
    }

private:
    std::shared_ptr<QMMomentum> p_x{nullptr};
    std::shared_ptr<QMMomentum> p_y{nullptr};
    std::shared_ptr<QMMomentum> p_z{nullptr};
};

} // namespace mrchem
