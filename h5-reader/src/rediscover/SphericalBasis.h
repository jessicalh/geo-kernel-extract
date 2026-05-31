// SphericalBasis — irreducible decomposition in the LIBRARY's T2 component
// order, so reader-side DFT tensors land in the SAME basis as the H5 kernel
// tensors (which the producer decomposed this way).
//
// The reader's existing `model::SphericalDecomposition` uses a DIFFERENT T2
// component order ([xx-yy, zz, xy, xz, yz]); its own header notes it is not
// claimed to match the library and is fine for the viewer's display, which
// never cross-compares. The rediscover substrate DOES cross-compare DFT
// tensors against the producer's kernels, so it needs the library order.
//
// This is a faithful port of nmr::SphericalTensor::Decompose
// (src/Types.cpp:25-60): T0 = trace/3; T1 = Cartesian antisymmetric
// pseudovector (v_x=A_yz, v_y=A_zx, v_z=A_xy); T2 = isometric real-SH
// basis [xy, yz, zz, xz, xx-yy] (m = -2..+2), Frobenius-norm preserving.
// The reader does not link the library, so the convention is re-implemented,
// pinned by a fixture against analytic values.

#pragma once

#include "../model/Types.h"  // Vec3, Mat3, SphericalTensor

namespace h5reader::rediscover {

// Decompose a (generally asymmetric) 3x3 tensor into T0/T1/T2 in the library
// component order + isometric normalization. Matches
// nmr::SphericalTensor::Decompose component-for-component.
model::SphericalTensor DecomposeLibrary(const model::Mat3& sigma);

}  // namespace h5reader::rediscover
