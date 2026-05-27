// SphericalDecomposition — irreducible decomposition of a 3x3 shielding tensor
// into a SphericalTensor (T0/T1/T2). Pure tensor math (Eigen only, no Qt/VTK).

#pragma once

#include "Types.h"  // Mat3, SphericalTensor

namespace h5reader::model {

// Decompose a (generally ASYMMETRIC) 3x3 shielding tensor sigma:
//   T0 = isotropic shielding = trace/3 (ppm) — the conventional iso and our
//        fail-loud anchor (must equal ORCA's printed iso=).
//   T1 = the antisymmetric part as the Cartesian pseudovector (A_yz, A_zx, A_xy),
//        A_ij = (sigma_ij - sigma_ji)/2 — real and non-trivial for shielding,
//        so we keep it rather than symmetrising it away.
//   T2 = the symmetric-traceless part in a real basis normalised so that
//        |T2| (SphericalTensor::T2Magnitude) equals the Frobenius norm of that
//        part — i.e. the anisotropy magnitude, basis-independent.
//
// NOTE: T0 and |T2| are exact and meaningful. The per-component T2 ordering is a
// valid rank-2 representation but is NOT claimed to match the library's exact
// e3nn component order — moot here, since kernel shielding is not compared
// against DFT for the trajectory. If that comparison returns, reconcile the T2
// basis against the library's Decompose then.
SphericalTensor DecomposeShielding(const Mat3& sigma);

}  // namespace h5reader::model
