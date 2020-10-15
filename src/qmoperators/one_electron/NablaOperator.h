#pragma once

#include "qmoperators/QMOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class QMNabla final : public QMOperator {
public:
    QMNabla(int d, std::shared_ptr<mrcpp::DerivativeOperator<3>> D);

private:
    const int apply_dir;
    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative;

    Orbital apply(Orbital inp) override;
    Orbital dagger(Orbital inp) override;
};

class NablaOperator final : public RankOneTensorOperator<3> {
public:
    NablaOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D) {
        d_x = std::make_shared<QMNabla>(0, D);
        d_y = std::make_shared<QMNabla>(1, D);
        d_z = std::make_shared<QMNabla>(2, D);

        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &d = (*this);
        d[0] = d_x;
        d[1] = d_y;
        d[2] = d_z;
        d[0].name() = "del[x]";
        d[1].name() = "del[y]";
        d[2].name() = "del[z]";
    }

private:
    std::shared_ptr<QMNabla> d_x{nullptr};
    std::shared_ptr<QMNabla> d_y{nullptr};
    std::shared_ptr<QMNabla> d_z{nullptr};
};

} // namespace mrchem
