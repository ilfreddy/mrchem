#include "qmoperators/one_electron/QMPotential.h"
#include "qmoperators/RankZeroTensorOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {
namespace qmoperator {

RankZeroTensorOperator dot(RankOneTensorOperator<3> &A, RankOneTensorOperator<3> &B);

} // namespace qmoperator
} // namespace mrchem
