#pragma once

#include "qmoperators/QMOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class GradZoraOperator final : public RankOneTensorOperator<3> {
public:
    GradZoraOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D, std::shared_ptr<ZoraPotential> Z) {

        computeGradComponent(dkdx, 0, D, Z);
        computeGradComponent(dkdy, 1, D, Z);
        computeGradComponent(dkdz, 2, D, Z);

        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &d = (*this);
        //        p_x = std::make_shared<QMMomentum>(0, D);
        //        p_y = std::make_shared<QMMomentum>(1, D);
        //        p_z = std::make_shared<QMMomentum>(2, D);
        d[0] = std::make_shared<QMFunction>(&dkdx);
        d[0] = std::make_shared<QMFunction>(&dkdy);
        d[0] = std::make_shared<QMFunction>(&dkdz);
        d[0].name() = "grad_x_kappa";
        d[1].name() = "grad_y_kappa";
        d[2].name() = "grad_z_kappa";
    }

private:
    std::shared_ptr<ZoraPotential> zora;
    QMFunction dkdx;
    QMFunction dkdy;
    QMFunction dkdz;

    void computeGradComponent(QMFunction component,
                              int dir,
                              std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                              std::shared_ptr<ZoraPotential> Z);
};

} // namespace mrchem
