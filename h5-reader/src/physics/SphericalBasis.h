// SphericalBasis — Reader's canonical irreducible decomposition of a Cartesian
// rank-2 tensor. T2 uses the project order [xy, yz, zz, xz, xx-yy] and an
// isometric normalization, so its Euclidean norm equals the Frobenius norm of
// the symmetric-traceless Cartesian tensor.

#pragma once

#include "../model/Types.h"  // Vec3, Mat3, SphericalTensor

namespace h5reader::physics {

// Decompose a generally asymmetric 3x3 tensor into T0, T1, and T2.
model::SphericalTensor DecomposeLibrary(const model::Mat3& sigma);

// Reconstruct the symmetric matrix whose library-basis T0/T2 components are
// exactly `t0` and `t2`. This is the inverse of the T2 convention documented
// above; T1 is antisymmetric and is not represented in an EFG T2 tensor.
model::Mat3 ReconstructLibraryT2Matrix(double t0, const std::array<double, 5>& t2);
model::Mat3 ReconstructLibraryT2Matrix(const std::array<double, 5>& t2);

}  // namespace h5reader::physics
