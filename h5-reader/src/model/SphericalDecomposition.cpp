#include "SphericalDecomposition.h"

namespace h5reader::model {

SphericalTensor DecomposeShielding(const Mat3& s) {
    SphericalTensor t;

    // T0: isotropic shielding (ppm).
    const double iso = (s(0, 0) + s(1, 1) + s(2, 2)) / 3.0;
    t.T0 = iso;

    // T1: antisymmetric part as the Cartesian pseudovector (A_yz, A_zx, A_xy).
    t.T1[0] = 0.5 * (s(1, 2) - s(2, 1));  // A_yz
    t.T1[1] = 0.5 * (s(2, 0) - s(0, 2));  // A_zx
    t.T1[2] = 0.5 * (s(0, 1) - s(1, 0));  // A_xy

    // T2: symmetric-traceless part. Component normalisation chosen so that
    // sum(T2^2) == Frobenius norm^2 of the symmetric-traceless tensor, i.e.
    // |T2| is the anisotropy magnitude (verified: the diagonal pair below
    // contributes Sxx^2+Syy^2+Szz^2 and the off-diagonals 2*(Sxy^2+Sxz^2+Syz^2)).
    const double Sxx = s(0, 0) - iso;
    const double Syy = s(1, 1) - iso;
    const double Szz = s(2, 2) - iso;
    const double Sxy = 0.5 * (s(0, 1) + s(1, 0));
    const double Sxz = 0.5 * (s(0, 2) + s(2, 0));
    const double Syz = 0.5 * (s(1, 2) + s(2, 1));

    constexpr double kSqrt2  = 1.41421356237309515;  // sqrt(2)
    constexpr double kSqrt32 = 1.22474487139158905;  // sqrt(3/2)
    t.T2[0] = (Sxx - Syy) / kSqrt2;
    t.T2[1] = Szz * kSqrt32;
    t.T2[2] = Sxy * kSqrt2;
    t.T2[3] = Sxz * kSqrt2;
    t.T2[4] = Syz * kSqrt2;

    return t;
}

}  // namespace h5reader::model
