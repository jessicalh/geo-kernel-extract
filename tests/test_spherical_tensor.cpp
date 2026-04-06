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
