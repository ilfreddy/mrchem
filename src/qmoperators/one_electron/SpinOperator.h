#pragma once

#include "qmoperators/QMOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class QMSpin final : public QMOperator {
public:
    QMSpin(int d)
            : D(d) {}

private:
    const int D;

    Orbital apply(Orbital inp) override;
    Orbital dagger(Orbital inp) override;
};

class SpinOperator final : public RankOneTensorOperator<3> {
public:
    SpinOperator() {
        s_x = std::make_shared<QMSpin>(0);
        s_y = std::make_shared<QMSpin>(1);
        s_z = std::make_shared<QMSpin>(2);

        // Invoke operator= to assign *this operator
        RankOneTensorOperator &s = *this;
        s[0] = s_x;
        s[1] = s_y;
        s[2] = s_z;
        s[0].name() = "s[x]";
        s[1].name() = "s[y]";
        s[2].name() = "s[z]";
    }

private:
    std::shared_ptr<QMSpin> s_x{nullptr};
    std::shared_ptr<QMSpin> s_y{nullptr};
    std::shared_ptr<QMSpin> s_z{nullptr};
};

} // namespace mrchem
