#include "SphericalBasis.h"

#include <cmath>

namespace h5reader::rediscover {

model::SphericalTensor DecomposeLibrary(const model::Mat3& sigma) {
    model::SphericalTensor st;

    // T0: isotropic = trace / 3.
    st.T0 = sigma.trace() / 3.0;

    // T1: antisymmetric pseudovector (Cartesian Levi-Civita dual).
    //   A_ij = (s_ij - s_ji)/2; v_x = A_yz, v_y = A_zx, v_z = A_xy.
    st.T1[0] = 0.5 * (sigma(1, 2) - sigma(2, 1));
    st.T1[1] = 0.5 * (sigma(2, 0) - sigma(0, 2));
    st.T1[2] = 0.5 * (sigma(0, 1) - sigma(1, 0));

    // Traceless symmetric part S_ij = (s_ij + s_ji)/2 - (trace/3) delta_ij.
    const double Sxx = sigma(0, 0) - st.T0;
    const double Syy = sigma(1, 1) - st.T0;
    const double Szz = sigma(2, 2) - st.T0;
    const double Sxy = 0.5 * (sigma(0, 1) + sigma(1, 0));
    const double Sxz = 0.5 * (sigma(0, 2) + sigma(2, 0));
    const double Syz = 0.5 * (sigma(1, 2) + sigma(2, 1));

    // T2: isometric real spherical-harmonic basis, m = -2..+2,
    // component order [xy, yz, zz, xz, xx-yy]; sum|T2_m|^2 = sum S_ij^2.
    const double kSqrt2          = std::sqrt(2.0);
    const double kSqrtThreeHalves = std::sqrt(3.0 / 2.0);
    st.T2[0] = kSqrt2 * Sxy;                  // m = -2
    st.T2[1] = kSqrt2 * Syz;                  // m = -1
    st.T2[2] = kSqrtThreeHalves * Szz;        // m =  0
    st.T2[3] = kSqrt2 * Sxz;                  // m = +1
    st.T2[4] = (1.0 / kSqrt2) * (Sxx - Syy);  // m = +2

    return st;
}

}  // namespace h5reader::rediscover
