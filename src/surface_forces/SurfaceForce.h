#ifndef SURFACEFORCE_H
#define SURFACEFORCE_H

#include <vector>
#include "qmfunctions/Orbital.h"
#include "chemistry/Molecule.h"
#include <nlohmann/json.hpp>
#include <Eigen/Dense>

namespace surface_force {

// Function declaration
Eigen::MatrixXd surface_forces(mrchem::Molecule &mol, mrchem::OrbitalVector &Phi, double prec, const json &json_fock);

} // namespace surface_force

#endif // SURFACEFORCE_H

