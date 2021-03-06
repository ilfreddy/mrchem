#pragma once

#include "H_E_dip.h"
#include "qmoperators/RankZeroTensorOperator.h"

/** @class ElectricFieldOperator
 *
 * @brief External electric field operator
 *
 * An external electric field interacts with the molecular dipole
 * moment. The operator is simply implemented as a scalar product of
 * the dipole moment operator with the electric field vector.
 *
 * It implements also explicit trace functions for the nuclear contributions
 *
 */

namespace mrchem {

/** @brief constructor
 *
 * @param[in] f the external electric field
 *
 */
class ElectricFieldOperator final : public RankZeroTensorOperator {
public:
    ElectricFieldOperator(const std::array<double, 3> &f, const mrcpp::Coord<3> &o)
            : field(f[0], f[1], f[2])
            , dipole(o) {
        RankZeroTensorOperator &d_x = this->dipole[0];
        RankZeroTensorOperator &d_y = this->dipole[1];
        RankZeroTensorOperator &d_z = this->dipole[2];

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &HEF = (*this);
        HEF = -f[0] * d_x - f[1] * d_y - f[2] * d_z;
        HEF.name() = "E . mu_E";
    }
    ElectricFieldOperator(const Eigen::Vector3d &f, const mrcpp::Coord<3> &o)
            : field(f)
            , dipole(o) {
        RankZeroTensorOperator &d_x = this->dipole[0];
        RankZeroTensorOperator &d_y = this->dipole[1];
        RankZeroTensorOperator &d_z = this->dipole[2];

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &HEF = (*this);
        HEF = -f[0] * d_x - f[1] * d_y - f[2] * d_z;
        HEF.name() = "E . mu_E";
    }

private:
    Eigen::Vector3d field; ///< the external field vector
    H_E_dip dipole;        ///< the dipole moment operator
};

} // namespace mrchem
