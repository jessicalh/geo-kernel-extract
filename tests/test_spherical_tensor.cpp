#include <gtest/gtest.h>
#include "Types.h"
#include <cmath>

using namespace nmr;

// ============================================================================
// Decompose a known Mat3, reconstruct, verify roundtrip.
// ============================================================================

TEST(SphericalTensor, RoundtripSymmetric) {
    // Symmetric matrix (T1 should be zero)
    Mat3 sigma;
    sigma << 10.0, 1.0, 2.0,
              1.0, 20.0, 3.0,
              2.0, 3.0, 30.0;

    auto st = SphericalTensor::Decompose(sigma);
    Mat3 reconstructed = st.Reconstruct();

    // Roundtrip should be exact (floating point)
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(reconstructed(i,j), sigma(i,j), 1e-12)
                << "Mismatch at (" << i << "," << j << ")";
}


TEST(SphericalTensor, RoundtripAsymmetric) {
    // General asymmetric matrix
    Mat3 sigma;
    sigma << 1.0, 0.5, 0.3,
             0.2, 2.0, 0.7,
             0.1, 0.4, 3.0;

    auto st = SphericalTensor::Decompose(sigma);
    Mat3 reconstructed = st.Reconstruct();

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(reconstructed(i,j), sigma(i,j), 1e-12)
                << "Mismatch at (" << i << "," << j << ")";
}


TEST(SphericalTensor, RoundtripIdentity) {
    Mat3 sigma = Mat3::Identity();
    auto st = SphericalTensor::Decompose(sigma);
    Mat3 reconstructed = st.Reconstruct();

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(reconstructed(i,j), sigma(i,j), 1e-12);
}


// ============================================================================
// Verify T0 against hand calculation: T0 = trace / 3
// ============================================================================

TEST(SphericalTensor, T0IsTrace) {
    Mat3 sigma;
    sigma << 10.0, 1.0, 2.0,
              1.0, 20.0, 3.0,
              2.0, 3.0, 30.0;

    auto st = SphericalTensor::Decompose(sigma);
    double expected_T0 = (10.0 + 20.0 + 30.0) / 3.0;
    EXPECT_NEAR(st.T0, expected_T0, 1e-12);
}


// ============================================================================
// Verify T1 against hand calculation.
// T1 = antisymmetric pseudovector: v_x = (s_yz - s_zy)/2, etc.
// For a symmetric matrix, T1 should be zero.
// ============================================================================

TEST(SphericalTensor, T1ZeroForSymmetric) {
    Mat3 sigma;
    sigma << 10.0, 1.0, 2.0,
              1.0, 20.0, 3.0,
              2.0, 3.0, 30.0;

    auto st = SphericalTensor::Decompose(sigma);
    for (int i = 0; i < 3; ++i)
        EXPECT_NEAR(st.T1[i], 0.0, 1e-12);
}


TEST(SphericalTensor, T1HandCalculation) {
    Mat3 sigma;
    sigma << 1.0, 0.5, 0.3,
             0.2, 2.0, 0.7,
             0.1, 0.4, 3.0;

    auto st = SphericalTensor::Decompose(sigma);

    // v_x = (s_yz - s_zy)/2 = (0.7 - 0.4)/2 = 0.15
    EXPECT_NEAR(st.T1[0], 0.15, 1e-12);
    // v_y = (s_zx - s_xz)/2 = (0.1 - 0.3)/2 = -0.1
    EXPECT_NEAR(st.T1[1], -0.1, 1e-12);
    // v_z = (s_xy - s_yx)/2 = (0.5 - 0.2)/2 = 0.15
    EXPECT_NEAR(st.T1[2], 0.15, 1e-12);
}


// ============================================================================
// Verify T2 against hand calculation for a diagonal matrix.
// For sigma = diag(a, b, c): trace = a+b+c, T0 = trace/3.
// Traceless symmetric: Sxx = a-T0, Syy = b-T0, Szz = c-T0.
// Off-diagonal S = 0.
// T2[0] = sqrt(2)*Sxy = 0
// T2[1] = sqrt(2)*Syz = 0
// T2[2] = sqrt(3/2)*Szz
// T2[3] = sqrt(2)*Sxz = 0
// T2[4] = (Sxx-Syy)/sqrt(2)
// ============================================================================

TEST(SphericalTensor, T2DiagonalMatrix) {
    double a = 10.0, b = 20.0, c = 30.0;
    Mat3 sigma = Mat3::Zero();
    sigma(0,0) = a; sigma(1,1) = b; sigma(2,2) = c;

    auto st = SphericalTensor::Decompose(sigma);

    double T0 = (a + b + c) / 3.0;
    double Szz = c - T0;
    double Sxx = a - T0;
    double Syy = b - T0;

    EXPECT_NEAR(st.T2[0], 0.0, 1e-12);   // m=-2: sqrt(2)*Sxy=0
    EXPECT_NEAR(st.T2[1], 0.0, 1e-12);   // m=-1: sqrt(2)*Syz=0
    EXPECT_NEAR(st.T2[2], std::sqrt(3.0/2.0) * Szz, 1e-12);  // m=0
    EXPECT_NEAR(st.T2[3], 0.0, 1e-12);   // m=+1: sqrt(2)*Sxz=0
    EXPECT_NEAR(st.T2[4], (Sxx - Syy) / std::sqrt(2.0), 1e-12);  // m=+2
}


// ============================================================================
// Isometric norm preservation: sum|T2_m|^2 == sum S_ij^2
// ============================================================================

TEST(SphericalTensor, IsometricNormPreservation) {
    Mat3 sigma;
    sigma << 1.0, 0.5, 0.3,
             0.2, 2.0, 0.7,
             0.1, 0.4, 3.0;

    auto st = SphericalTensor::Decompose(sigma);

    // Compute Frobenius norm of traceless symmetric part
    double T0 = sigma.trace() / 3.0;
    Mat3 symm = 0.5 * (sigma + sigma.transpose());
    Mat3 S = symm - T0 * Mat3::Identity();
    double S_frobenius_sq = S.squaredNorm();

    double T2_norm_sq = 0.0;
    for (double v : st.T2) T2_norm_sq += v * v;

    EXPECT_NEAR(T2_norm_sq, S_frobenius_sq, 1e-10);
}


// ============================================================================
// Zero matrix decomposes to all zeros
// ============================================================================

TEST(SphericalTensor, ZeroMatrix) {
    auto st = SphericalTensor::Decompose(Mat3::Zero());
    EXPECT_NEAR(st.T0, 0.0, 1e-15);
    for (int i = 0; i < 3; ++i) EXPECT_NEAR(st.T1[i], 0.0, 1e-15);
    for (int i = 0; i < 5; ++i) EXPECT_NEAR(st.T2[i], 0.0, 1e-15);
}


// ============================================================================
// T2Magnitude
// ============================================================================

TEST(SphericalTensor, T2Magnitude) {
    Mat3 sigma;
    sigma << 10.0, 0.0, 0.0,
              0.0, 20.0, 0.0,
              0.0, 0.0, 30.0;
    auto st = SphericalTensor::Decompose(sigma);
    double mag = st.T2Magnitude();
    EXPECT_GT(mag, 0.0);

    // For identity, T2 should be zero
    auto st_id = SphericalTensor::Decompose(Mat3::Identity());
    EXPECT_NEAR(st_id.T2Magnitude(), 0.0, 1e-12);
}


// ============================================================================
// Cross-check: our T2 basis is the real spherical harmonic basis Y^2_m
// up to a single fixed normalization constant.
//
// This pins the basis against a known reference (hand-coded Y^2_m
// closed forms — equivalent to sphericart's SphericalHarmonics<double>
// at l_max=2 up to the same fixed constant, which we cross-checked by
// hand) so that any future edit to the SphericalTensor basis has to
// break this test before it corrupts downstream consumers that depend
// on the real-spherical-harmonic convention.
//
// For a unit vector r̂ = (x, y, z), decompose the outer product minus
// isotropic part:
//   M = r̂ r̂ᵀ - I/3
// M is symmetric and traceless, so T0 = 0, T1 = 0, and all information
// lives in T2. The five T2 components are proportional to Y^2_m(r̂):
//
//   T2[m] = k · Y^2_m(r̂)    for all m = -2..+2
//
// with k = 2·√(2π/15) ≈ 1.2947. This test verifies the ratio is
// constant across all 5 m's for many random unit vectors.
// ============================================================================

namespace {

// Real spherical harmonics Y^2_m evaluated at a unit vector (x, y, z).
// Signs + prefactors follow the common chemistry convention (no Condon-
// Shortley phase twist for real SH).
//
//   Y^2_{-2} = (1/2)√(15/π) · x·y
//   Y^2_{-1} = (1/2)√(15/π) · y·z
//   Y^2_{ 0} = (1/4)√(5/π)  · (3z² - 1)
//   Y^2_{+1} = (1/2)√(15/π) · x·z
//   Y^2_{+2} = (1/4)√(15/π) · (x² - y²)
std::array<double, 5> Y2m(double x, double y, double z) {
    const double PI = 3.14159265358979323846;
    const double c15 = 0.5  * std::sqrt(15.0 / PI);
    const double c5  = 0.25 * std::sqrt( 5.0 / PI);
    const double c15_4 = 0.25 * std::sqrt(15.0 / PI);
    return {
        c15   * x * y,          // m = -2
        c15   * y * z,          // m = -1
        c5    * (3.0*z*z - 1.0),// m =  0
        c15   * x * z,          // m = +1
        c15_4 * (x*x - y*y),    // m = +2
    };
}

}  // namespace


TEST(SphericalTensor, T2BasisMatchesY2m) {
    // Expected constant ratio T2[m] / Y^2_m = 2·√(2π/15).
    const double PI = 3.14159265358979323846;
    const double k_expected = 2.0 * std::sqrt(2.0 * PI / 15.0);

    // A spread of unit vectors covering all octants + axes + generic
    // orientations. Each is normalised to unit length.
    const std::vector<Vec3> unit_vectors = {
        Vec3(1.0, 0.0, 0.0).normalized(),
        Vec3(0.0, 1.0, 0.0).normalized(),
        Vec3(0.0, 0.0, 1.0).normalized(),
        Vec3(1.0, 1.0, 0.0).normalized(),
        Vec3(1.0, 0.0, 1.0).normalized(),
        Vec3(0.0, 1.0, 1.0).normalized(),
        Vec3(1.0, 1.0, 1.0).normalized(),
        Vec3(-1.0, 2.0, -3.0).normalized(),
        Vec3(0.3, -0.7, 0.5).normalized(),
        Vec3(0.9, 0.1, -0.4).normalized(),
    };

    for (const Vec3& r : unit_vectors) {
        // M = r̂ r̂ᵀ - I/3 is symmetric traceless.
        Mat3 M = r * r.transpose() - Mat3::Identity() / 3.0;
        const auto st = SphericalTensor::Decompose(M);

        // T0 and T1 must vanish on a symmetric traceless matrix.
        EXPECT_NEAR(st.T0, 0.0, 1e-13)
            << "T0 should vanish for symmetric traceless input";
        for (int i = 0; i < 3; ++i)
            EXPECT_NEAR(st.T1[i], 0.0, 1e-13)
                << "T1[" << i << "] should vanish for symmetric input";

        // T2[m] / Y^2_m should equal k_expected for every m.
        const auto y = Y2m(r.x(), r.y(), r.z());
        for (int m = 0; m < 5; ++m) {
            // Skip m's where Y^2_m is near-zero for this r̂ (e.g.
            // Y^2_{-2} ∝ xy vanishes when x or y is zero); the ratio
            // would be 0/0. The other m's still constrain the basis.
            if (std::abs(y[m]) < 1e-12) continue;

            const double ratio = st.T2[m] / y[m];
            EXPECT_NEAR(ratio, k_expected, 1e-10)
                << "T2[" << m << "]/Y^2_{" << (m - 2)
                << "} ratio breaks the basis convention"
                   " for r̂=(" << r.x() << "," << r.y() << "," << r.z() << ")";
        }
    }
}


// ============================================================================
// Orthogonality of the T2 basis: five symmetric-traceless basis matrices
// each decompose to a T2 vector supported on exactly one m with the
// expected normalization, and zero on all other m's.
// ============================================================================

TEST(SphericalTensor, T2BasisOrthogonality) {
    // Build each basis matrix B_m such that our Decompose(B_m) returns
    // T2 with a single unit component at index m. These are the
    // real-spherical-harmonic basis matrices in our isometric
    // normalization — inverting the coefficients in Decompose.
    const double SQRT2   = std::sqrt(2.0);
    const double SQRT3_2 = std::sqrt(3.0 / 2.0);
    const double INV_SQRT2 = 1.0 / SQRT2;

    auto make = [](double Sxx, double Syy, double Szz,
                   double Sxy, double Sxz, double Syz) {
        Mat3 m;
        m << Sxx, Sxy, Sxz,
             Sxy, Syy, Syz,
             Sxz, Syz, Szz;
        return m;
    };

    // Inverse of Decompose's T2 coefficients: each basis matrix has
    // one nonzero S component sized to give unit T2[m] on decomposition.
    //
    //   T2[0] = √2 · Sxy                → Sxy = 1/√2
    //   T2[1] = √2 · Syz                → Syz = 1/√2
    //   T2[2] = √(3/2) · Szz            → Szz = √(2/3); matrix must be
    //                                     traceless, so Sxx = Syy = -Szz/2
    //                                     = -1/√6 (also gives Sxx-Syy=0).
    //   T2[3] = √2 · Sxz                → Sxz = 1/√2
    //   T2[4] = (Sxx - Syy)/√2          → Sxx = √2/2, Syy = -√2/2, Szz=0.
    const double SQRT2_3 = std::sqrt(2.0 / 3.0);   // = 1/SQRT3_2
    const double NEG_INV_SQRT6 = -1.0 / std::sqrt(6.0);  // = -√(2/3)/2
    const std::array<Mat3, 5> basis = {
        make(0, 0, 0,   INV_SQRT2,         0,         0),   // m=-2
        make(0, 0, 0,           0,         0, INV_SQRT2),   // m=-1
        make(NEG_INV_SQRT6, NEG_INV_SQRT6, SQRT2_3,
              0, 0, 0),                                      // m= 0
        make(0, 0, 0,           0, INV_SQRT2,         0),   // m=+1
        make( SQRT2/2.0, -SQRT2/2.0, 0,
              0, 0, 0),                                      // m=+2
    };

    for (int m = 0; m < 5; ++m) {
        const auto st = SphericalTensor::Decompose(basis[m]);
        EXPECT_NEAR(st.T0, 0.0, 1e-13)
            << "basis[" << m << "] T0 should vanish";
        for (int i = 0; i < 3; ++i)
            EXPECT_NEAR(st.T1[i], 0.0, 1e-13)
                << "basis[" << m << "] T1[" << i << "] should vanish";
        for (int k = 0; k < 5; ++k) {
            const double expected = (k == m) ? 1.0 : 0.0;
            EXPECT_NEAR(st.T2[k], expected, 1e-12)
                << "basis[" << m << "] T2[" << k
                << "] should be " << expected;
        }
    }
}
