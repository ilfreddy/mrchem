#include "qmoperators/RankZeroTensorOperator.h"
#include "qmoperators/RankOneTensorOperator.h"
#include "qmoperators/one_electron/QMPotential.h"
#include "qmoperators/qmoperator_utils.h"

namespace mrchem {

RankZeroTensorOperator qmoperator::dot(RankOneTensorOperator<3> &A, RankOneTensorOperator<3> &B) {
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

} // namespace mrchem
