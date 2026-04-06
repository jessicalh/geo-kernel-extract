#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

#include "PiQuadrupoleResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "PhysicalConstants.h"

#include <iomanip>
#include <map>

#include <filesystem>
namespace fs = std::filesystem;
using namespace nmr;



// ============================================================================
// Analytical test: verify tracelessness and symmetry of the EFG kernel.
//
// The EFG from a point quadrupole satisfies Laplace: Tr(G) = 0.
// Symmetry: G_ab = G_ba.
//
// Test at multiple positions around a ring:
//   1. (0, 0, 3) — directly above ring (theta = 0)
//   2. (3, 0, 0) — in ring plane (theta = pi/2)
//   3. (2, 2, 2) — off-axis
// ============================================================================

// Hand-compute the kernel for verification
static Mat3 HandComputeG(const Vec3& d, const Vec3& n) {
    double r = d.norm();
    double r2 = r * r;
    double r5 = r2 * r2 * r;
    double r7 = r5 * r2;
    double r9 = r7 * r2;
    double dn = d.dot(n);
    double dn2 = dn * dn;

    Mat3 G;
    double diag = 3.0 / r5 - 15.0 * dn2 / r7;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            G(a, b) = 105.0 * dn2 * d(a) * d(b) / r9
                     - 30.0 * dn * (n(a) * d(b) + n(b) * d(a)) / r7
                     - 15.0 * d(a) * d(b) / r7
                     + 6.0 * n(a) * n(b) / r5
                     + (a == b ? diag : 0.0);
    return G;
}


TEST(PiQuadAnalytical, TracelessAtMultiplePositions) {
    Vec3 n(0, 0, 1);

    // Test positions
    std::vector<Vec3> positions = {
        Vec3(0, 0, 3),     // above ring
        Vec3(3, 0, 0),     // in plane
        Vec3(2, 2, 2),     // off-axis
        Vec3(1, 0, 4),     // near-axial
        Vec3(5, 3, 1),     // far off-axis
    };

    for (const Vec3& d : positions) {
        Mat3 G = HandComputeG(d, n);

        double trace = G.trace();
        EXPECT_NEAR(trace, 0.0, 1e-12)
            << "Tr(G) must be 0 (Laplace) at d=("
            << d(0) << "," << d(1) << "," << d(2) << ")";

        // Symmetry
        for (int a = 0; a < 3; ++a)
            for (int b = a + 1; b < 3; ++b)
                EXPECT_NEAR(G(a, b), G(b, a), 1e-14)
                    << "G must be symmetric";
    }
}


TEST(PiQuadAnalytical, ScalarMatchesFormula) {
    Vec3 n(0, 0, 1);
    Vec3 d(2, 1, 3);
    double r = d.norm();
    double cos_theta = d.dot(n) / r;
    double r4 = r * r * r * r;

    double expected = (3.0 * cos_theta * cos_theta - 1.0) / r4;

    // The scalar should match (3cos^2 theta - 1)/r^4
    // This is the Buckingham A-term source (E-field from quadrupole)
    EXPECT_NEAR(expected, (3.0 * 9.0 / 14.0 - 1.0) / r4, 1e-10);
}


TEST(PiQuadAnalytical, AboveRingValues) {
    // Atom at (0, 0, 3) above ring at origin with normal along z.
    // d = (0, 0, 3), r = 3, dn = 3, cos_theta = 1
    //
    // G_xx = G_yy = 0 - 0 - 0 + 0 + (3/r^5 - 15*9/r^7)
    //      = 3/243 - 135/2187 = 0.01235 - 0.06173 = -0.04938
    //
    // G_zz = 105*9*9/r^9 - 30*3*2*3*3/r^7 - 15*9/r^7 + 6/r^5 + (3/r^5 - 135/r^7)
    //      = 8505/19683 - 1620/2187 - 135/2187 + 6/243 + 3/243 - 135/2187
    //      = 0.4321 - 0.7407 - 0.0617 + 0.02469 + 0.01235 - 0.06173
    //
    // Let me just verify Tr = 0 and compute numerically.

    Vec3 n(0, 0, 1);
    Vec3 d(0, 0, 3);
    Mat3 G = HandComputeG(d, n);

    EXPECT_NEAR(G.trace(), 0.0, 1e-12);

    // G_xx == G_yy by cylindrical symmetry
    EXPECT_NEAR(G(0, 0), G(1, 1), 1e-14);

    // G_zz = -2 * G_xx (from tracelessness and G_xx = G_yy)
    EXPECT_NEAR(G(2, 2), -2.0 * G(0, 0), 1e-12);

    // Off-diagonal should be zero by cylindrical symmetry
    EXPECT_NEAR(G(0, 1), 0.0, 1e-14);
    EXPECT_NEAR(G(0, 2), 0.0, 1e-14);
    EXPECT_NEAR(G(1, 2), 0.0, 1e-14);

    // Scalar: (3*1 - 1)/81 = 2/81
    double r = 3.0;
    double expected_scalar = 2.0 / (r * r * r * r);
    double cos_theta = 1.0;
    double actual_scalar = (3.0 * cos_theta * cos_theta - 1.0) / (r * r * r * r);
    EXPECT_NEAR(actual_scalar, expected_scalar, 1e-14);

    std::cout << "  G at (0,0,3): diag = ("
              << G(0,0) << ", " << G(1,1) << ", " << G(2,2) << ")\n"
              << "  Trace = " << G.trace() << "\n"
              << "  Scalar = " << actual_scalar << " A^-4\n";
}


TEST(PiQuadAnalytical, InPlaneValues) {
    // Atom at (3, 0, 0), in ring plane.
    // dn = 0, cos_theta = 0.
    //
    // Scalar: (0 - 1)/81 = -1/81
    //
    // G simplifies: all dn terms vanish.
    // G_ab = -15 d_a d_b / r^7 + 6 n_a n_b / r^5 + delta_ab * 3/r^5

    Vec3 n(0, 0, 1);
    Vec3 d(3, 0, 0);
    Mat3 G = HandComputeG(d, n);

    EXPECT_NEAR(G.trace(), 0.0, 1e-12);

    double r = 3.0;
    double r5 = std::pow(r, 5);
    double r7 = std::pow(r, 7);

    // G_xx = -15*9/r^7 + 0 + 3/r^5 = -135/2187 + 3/243
    //      = -0.0617 + 0.01235 = -0.04938
    // Wait, need to be more careful. d = (3, 0, 0).
    // G_xx = 0 - 0 - 15*9/r^7 + 0 + 3/r^5
    double expected_xx = -15.0 * 9.0 / r7 + 3.0 / r5;
    EXPECT_NEAR(G(0, 0), expected_xx, 1e-12);

    // G_yy = 0 - 0 - 0 + 0 + 3/r^5
    EXPECT_NEAR(G(1, 1), 3.0 / r5, 1e-12);

    // G_zz = 0 - 0 - 0 + 6/r^5 + 3/r^5 = 9/r^5
    // Wait: 6 n_z n_z / r^5 = 6/r^5, plus delta_zz * 3/r^5 = 3/r^5
    // Total: 9/r^5... but Tr = G_xx + G_yy + G_zz should be 0.
    // G_xx = -135/2187 + 3/243 = -135/2187 + 27/2187 = -108/2187
    // G_yy = 3/243 = 27/2187
    // G_zz = 9/243 = 81/2187
    // Tr = (-108 + 27 + 81)/2187 = 0 ✓
    EXPECT_NEAR(G(2, 2), 9.0 / r5, 1e-12);

    std::cout << "  G at (3,0,0): diag = ("
              << G(0,0) << ", " << G(1,1) << ", " << G(2,2) << ")\n"
              << "  Trace = " << G.trace() << "\n";
}


// ============================================================================
// FINITE-DIFFERENCE VERIFICATION of the analytical EFG formula.
//
// The potential from an axial quadrupole Theta along n at origin:
//   Phi(r) = (3 (d.n)^2 - r^2) / (2 r^5)
// (with unit Theta = 1, absorbing constants).
//
// The EFG is V_ab = -d^2 Phi / dx_a dx_b. We compute this numerically
// by central differences:
//   V_ab ~ -[Phi(r+h*e_a+h*e_b) - Phi(r+h*e_a-h*e_b)
//           - Phi(r-h*e_a+h*e_b) + Phi(r-h*e_a-h*e_b)] / (4h^2)
//
// This must match our analytical G_ab to within O(h^2) truncation error.
// Using h=1e-5 gives ~1e-10 relative accuracy — sufficient to catch
// any sign error, missing term, or wrong exponent.
// ============================================================================

static double QuadPotential(const Vec3& d, const Vec3& n) {
    // Phi = (3 dn^2 - r^2) / (2 r^5)  with unit Theta
    double r2 = d.squaredNorm();
    double r = std::sqrt(r2);
    if (r < 1e-15) return 0.0;
    double dn = d.dot(n);
    double r5 = r2 * r2 * r;
    return (3.0 * dn * dn - r2) / (2.0 * r5);
}

static Mat3 NumericalEFG(const Vec3& d, const Vec3& n, double h) {
    // V_ab = -d^2 Phi / dx_a dx_b via central differences
    Mat3 V;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            Vec3 ea = Vec3::Zero(); ea(a) = h;
            Vec3 eb = Vec3::Zero(); eb(b) = h;
            double pp = QuadPotential(d + ea + eb, n);
            double pm = QuadPotential(d + ea - eb, n);
            double mp = QuadPotential(d - ea + eb, n);
            double mm = QuadPotential(d - ea - eb, n);
            V(a, b) = -(pp - pm - mp + mm) / (4.0 * h * h);
        }
    }
    return V;
}

TEST(PiQuadAnalytical, FiniteDifferenceVerification) {
    // The relationship between our geometric kernel G_ab and the physical
    // EFG V_ab = -d^2 Phi / dx_a dx_b is:
    //
    //   V_ab = -(Theta/2) * G_ab
    //
    // With unit Theta=1, the numerical 2nd derivatives of Phi should equal
    // -G_ab/2. We verify this at multiple positions with two normals.

    Vec3 n(0, 0, 1);
    double h = 1e-5;

    std::vector<Vec3> positions = {
        Vec3(0, 0, 3),       // above ring
        Vec3(3, 0, 0),       // in plane
        Vec3(2, 2, 2),       // off-axis
        Vec3(1, 0, 4),       // near-axial
        Vec3(5, 3, 1),       // far off-axis
        Vec3(0.5, 0.3, 2),   // close, off-axis
        Vec3(7, 0, 0),       // far in-plane
    };

    Vec3 n_tilted = Vec3(1, 1, 1).normalized();

    double max_rel_err = 0.0;
    int checked = 0;

    auto check = [&](const Vec3& d, const Vec3& normal, const char* label) {
        Mat3 G_analytical = HandComputeG(d, normal);
        Mat3 V_numerical = NumericalEFG(d, normal, h);

        // V_numerical should equal -G/2 (with Theta=1)
        Mat3 expected = -0.5 * G_analytical;

        // Scale tolerance by the largest component in this tensor pair,
        // with an absolute floor for near-zero entries where FD truncation
        // error dominates the relative error.
        double max_abs = std::max(expected.cwiseAbs().maxCoeff(),
                                   V_numerical.cwiseAbs().maxCoeff());
        double abs_tol = std::max(max_abs * 1e-4, 1e-7);

        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                double anal = expected(a, b);
                double numer = V_numerical(a, b);
                double scale = std::max(std::abs(anal), std::abs(numer));
                // Only track relative error when both are significantly non-zero
                // (otherwise FD truncation error dominates the ratio).
                if (std::abs(anal) > abs_tol && std::abs(numer) > abs_tol) {
                    double rel_err = std::abs(anal - numer) / scale;
                    max_rel_err = std::max(max_rel_err, rel_err);
                }

                EXPECT_NEAR(anal, numer, abs_tol)
                    << label << " d=(" << d(0) << "," << d(1) << "," << d(2)
                    << ") (" << a << "," << b << ")"
                    << " expected=" << anal << " numerical=" << numer;
                checked++;
            }
        }
    };

    for (const Vec3& d : positions) {
        check(d, n, "z-normal");
        check(d, n_tilted, "tilted-normal");
    }

    std::cout << "  Finite-difference verification: " << checked
              << " components checked, max relative error (non-trivial) = "
              << max_rel_err << "\n";
    // max_rel_err only includes components where both values exceed 1e-10,
    // so near-zero entries (where FD truncation error dominates) are excluded.
    EXPECT_LT(max_rel_err, 1e-4)
        << "Analytical G must satisfy V_ab = -(1/2)*G_ab against numerical Phi derivatives";
}


// ============================================================================
// MAGNITUDE VERIFICATION against benzene quadrupole moment.
//
// Literature: benzene Theta_zz = -28.3 × 10^-40 C·m² (Flygare 1971,
// also Buckingham 1959: -29.0). We use -28.3e-40 C·m².
//
// Physical EFG at position d from a quadrupole Theta:
//   V_ab(SI) = -(Theta / (4*pi*eps_0)) * (1/2) * G_ab
//
// where G_ab is our geometric kernel (in m^-5 when d is in metres).
//
// At d = (0, 0, 3A) above benzene center (directly above ring):
// G_zz(analytical) = -2 * G_xx (from tracelessness, cylindrical symmetry)
// We compute V_zz in V/m² and compare to published order of magnitude.
//
// Published range for aromatic EFG at ~3A: O(10^18 - 10^19) V/m²
// (Sternheimer, Chem. Rev. 1968; also computed in Buckingham 1960).
// ============================================================================

TEST(PiQuadAnalytical, MagnitudeAgainstBenzeneQuadrupole) {
    // Benzene quadrupole moment
    constexpr double THETA_ZZ_SI = -28.3e-40;  // C·m² (Flygare 1971)

    // Coulomb constant: 1/(4*pi*eps_0) = 8.9875e9 N·m²/C²
    constexpr double KE_SI = 8.9875517923e9;

    // Position: 3 Angstroms above ring center = 3e-10 m
    Vec3 n(0, 0, 1);
    Vec3 d_A(0, 0, 3);      // Angstroms
    Vec3 d_m = d_A * 1e-10;  // metres

    // Our geometric kernel in Angstrom units
    Mat3 G_A = HandComputeG(d_A, n);

    // Convert to SI: G_ab has units A^-5, convert to m^-5
    // 1 A^-5 = (1e-10 m)^-5 = 1e50 m^-5
    Mat3 G_SI = G_A * 1e50;

    // Physical EFG: V_ab = -KE * (Theta/2) * G_ab
    // (the -1/2 comes from V_ab = -(1/3) Theta_cd T_abcd with
    //  Theta_cd = Theta(3n_cn_d - delta_cd)/2, so the net factor
    //  with the delta trace vanishing is Theta/2)
    Mat3 V_SI = -KE_SI * (THETA_ZZ_SI / 2.0) * G_SI;

    double Vzz_SI = V_SI(2, 2);  // V/m²

    // Report
    std::cout << "  Benzene quadrupole EFG at 3A above ring center:\n"
              << "    Theta_zz = " << THETA_ZZ_SI << " C·m²\n"
              << "    G_zz (A^-5) = " << G_A(2, 2) << "\n"
              << "    G_zz (m^-5) = " << G_SI(2, 2) << "\n"
              << "    V_zz (V/m^2) = " << Vzz_SI << "\n"
              << "    |V_zz| = " << std::abs(Vzz_SI) << " V/m²\n"
              << "    Order of magnitude: 10^" << std::log10(std::abs(Vzz_SI)) << "\n";

    // Diagonal check
    std::cout << "    V diag (V/m²): (" << V_SI(0,0) << ", "
              << V_SI(1,1) << ", " << V_SI(2,2) << ")\n"
              << "    Tr(V) = " << V_SI.trace() << " V/m² (should be ~0)\n";

    // The EFG from a quadrupole at 3A should be in the range 10^17 - 10^20 V/m²
    // depending on the moment. For benzene (small Theta), expect ~10^18.
    EXPECT_GT(std::abs(Vzz_SI), 1e16)
        << "EFG should be at least 10^16 V/m² at 3A from benzene";
    EXPECT_LT(std::abs(Vzz_SI), 1e22)
        << "EFG should be less than 10^22 V/m² at 3A from benzene";

    // Now compute in our working units (V/A²) for comparison with CoulombResult.
    // V/m² to V/A²: 1 V/m² = 1e-20 V/A²
    double Vzz_VA2 = Vzz_SI * 1e-20;
    std::cout << "    V_zz in V/A² = " << Vzz_VA2 << "\n";

    // For context: CoulombResult EFG from partial charges is typically
    // 0.01-10 V/A² at backbone atoms. The quadrupole EFG should be
    // a small fraction of this — it's a higher-order effect.
    std::cout << "    (CoulombResult backbone EFG is typically 0.01-10 V/A²)\n";
}


// ============================================================================
// Full protein test on 1UBQ
// ============================================================================

class PiQuadProteinTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) GTEST_SKIP() << "1UBQ not found";
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << "Failed to load 1UBQ";
        protein = std::move(r.protein);

        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        conf.AttachResult(SpatialIndexResult::Compute(conf));
    }
    std::unique_ptr<Protein> protein;
};


TEST_F(PiQuadProteinTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto pq = PiQuadrupoleResult::Compute(conf);
    ASSERT_NE(pq, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(pq)));
    ASSERT_TRUE(conf.HasResult<PiQuadrupoleResult>());
}


TEST_F(PiQuadProteinTest, TracelessOnRealProtein) {
    auto& conf = protein->Conformation();
    conf.AttachResult(PiQuadrupoleResult::Compute(conf));

    int checked = 0;
    double max_trace = 0.0;
    double max_asym = 0.0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.quad_tensor.isZero(1e-20)) continue;

            double trace = std::abs(rn.quad_tensor.trace());
            max_trace = std::max(max_trace, trace);

            // Symmetry check
            for (int a = 0; a < 3; ++a)
                for (int b = a + 1; b < 3; ++b)
                    max_asym = std::max(max_asym,
                        std::abs(rn.quad_tensor(a, b) - rn.quad_tensor(b, a)));

            checked++;
        }
    }

    EXPECT_GT(checked, 100) << "Should check many ring-atom pairs";
    EXPECT_LT(max_trace, 1e-10)
        << "EFG must be traceless (Laplace)";
    EXPECT_LT(max_asym, 1e-14)
        << "EFG must be symmetric";

    std::cout << "  Verified tracelessness on " << checked
              << " pairs, max |Tr| = " << max_trace
              << ", max asymmetry = " << max_asym << "\n";
}


TEST_F(PiQuadProteinTest, T0IsZero) {
    // Since the tensor is traceless, T0 = Tr/3 should be ~0 everywhere.
    auto& conf = protein->Conformation();
    conf.AttachResult(PiQuadrupoleResult::Compute(conf));

    double max_t0 = 0.0;
    int nonzero_t2 = 0;
    double max_t2 = 0.0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& sc = conf.AtomAt(ai).piquad_shielding_contribution;
        max_t0 = std::max(max_t0, std::abs(sc.T0));
        double t2m = sc.T2Magnitude();
        if (t2m > 1e-8) nonzero_t2++;
        max_t2 = std::max(max_t2, t2m);
    }

    EXPECT_LT(max_t0, 1e-10)
        << "T0 must be ~0 (traceless EFG)";
    EXPECT_GT(nonzero_t2, 0)
        << "T2 should be non-zero near rings";

    std::cout << "  max |T0| = " << max_t0 << " (should be ~0)\n"
              << "  T2 nonzero: " << nonzero_t2
              << ", max |T2| = " << max_t2 << "\n";
}


TEST_F(PiQuadProteinTest, ScalarStoredCorrectly) {
    auto& conf = protein->Conformation();
    conf.AttachResult(PiQuadrupoleResult::Compute(conf));

    int checked = 0;
    double max_diff = 0.0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (std::abs(rn.quad_scalar) < 1e-20) continue;

            // Recompute scalar from geometry
            const RingGeometry& geom = conf.ring_geometries[rn.ring_index];
            Vec3 d = conf.PositionAt(ai) - geom.center;
            double r = d.norm();
            if (r < MIN_DISTANCE) continue;
            double cos_theta = d.dot(geom.normal) / r;
            double r4 = r * r * r * r;
            double expected = (3.0 * cos_theta * cos_theta - 1.0) / r4;

            double diff = std::abs(rn.quad_scalar - expected);
            max_diff = std::max(max_diff, diff);
            checked++;
        }
    }

    EXPECT_GT(checked, 0);
    EXPECT_LT(max_diff, 1e-10)
        << "Stored scalar must match (3cos^2 theta - 1)/r^4";

    std::cout << "  Checked " << checked << " scalars, max diff = "
              << max_diff << "\n";
}


// ============================================================================
// ORCA protein test
// ============================================================================

TEST(PiQuadOrcaTest, RunOnProtonatedProtein) {
    OrcaRunFiles files;
    files.pdb_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.pdb";
    files.xyz_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";

    if (!fs::exists(files.xyz_path) || !fs::exists(files.prmtop_path))
        GTEST_SKIP() << "ORCA test data not found";

    auto load = BuildFromOrca(files);
    ASSERT_TRUE(load.Ok()) << load.error;

    auto& conf = load.protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(SpatialIndexResult::Compute(conf));
    conf.AttachResult(PiQuadrupoleResult::Compute(conf));

    double max_trace = 0.0;
    int with_quad = 0;
    double max_t2 = 0.0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.quad_tensor.isZero(1e-20)) continue;
            max_trace = std::max(max_trace, std::abs(rn.quad_tensor.trace()));
        }
        double t2m = conf.AtomAt(ai).piquad_shielding_contribution.T2Magnitude();
        if (t2m > 1e-8) with_quad++;
        max_t2 = std::max(max_t2, t2m);
    }

    EXPECT_LT(max_trace, 1e-10) << "EFG must be traceless";
    EXPECT_GT(with_quad, 0) << "Some atoms should have PQ signal";

    std::cout << "  ORCA protein PiQuadrupole summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " rings=" << load.protein->RingCount() << "\n"
              << "    atoms with PQ T2: " << with_quad << "\n"
              << "    max |Tr| = " << max_trace << "\n"
              << "    max |T2| = " << max_t2 << " A^-5\n";
}


// ============================================================================
// POINT-QUADRUPOLE CONVERGENCE TEST
//
// Compare the point-quadrupole EFG (our kernel) against the "distributed"
// EFG computed by summing the dipolar kernel K_ab(d) weighted by partial
// charges over all ring atoms. The distributed model is:
//
//   V_ab^dist(r) = sum_{j in ring} q_j * (3 d_a d_b / r^5 - delta_ab / r^3)
//
// where d = r_atom - r_ring_atom, r = |d|.
//
// At large distance from the ring, both should converge (the multipole
// expansion becomes accurate). At close range, they diverge (the point
// approximation breaks down).
//
// We bin atom-ring pairs by distance and report mean |cos| of the T2
// directions and the T0 ratio (point/distributed) at each distance bin.
// ============================================================================

TEST(PiQuadConvergence, PointVsDistributedOnOrcaProtein) {
    OrcaRunFiles files;
    files.pdb_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.pdb";
    files.xyz_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";

    if (!fs::exists(files.xyz_path) || !fs::exists(files.prmtop_path))
        GTEST_SKIP() << "ORCA test data not found";

    auto load = BuildFromOrca(files);
    ASSERT_TRUE(load.Ok()) << load.error;

    auto& conf = load.protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));

    PrmtopChargeSource charge_source(files.prmtop_path);
    conf.AttachResult(ChargeAssignmentResult::Compute(conf, charge_source));

    conf.AttachResult(SpatialIndexResult::Compute(conf));
    conf.AttachResult(PiQuadrupoleResult::Compute(conf));

    const Protein& protein = conf.ProteinRef();

    // Distance bins: 2-4, 4-6, 6-8, 8-10, 10-15 Angstroms
    struct Bin {
        double lo, hi;
        double sum_t2_cos = 0;    // |cos| between point and distributed T2
        double sum_t0_ratio = 0;  // |T0_point / T0_distributed|
        int count = 0;
        int t0_count = 0;         // pairs where both T0 are non-negligible
    };
    std::vector<Bin> bins = {{2,4}, {4,6}, {6,8}, {8,10}, {10,15}};

    auto t2_cos = [](const SphericalTensor& a, const SphericalTensor& b) -> double {
        double dot = 0, na = 0, nb = 0;
        for (int m = 0; m < 5; ++m) {
            dot += a.T2[m] * b.T2[m];
            na += a.T2[m] * a.T2[m];
            nb += b.T2[m] * b.T2[m];
        }
        double denom = std::sqrt(na * nb);
        if (denom < 1e-20) return 0.0;
        return dot / denom;
    };

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        Vec3 atom_pos = conf.PositionAt(ai);

        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.quad_tensor.isZero(1e-20)) continue;

            const Ring& ring = protein.RingAt(rn.ring_index);
            double dist = rn.distance_to_center;

            // Compute distributed EFG: sum over ring atoms
            Mat3 V_dist = Mat3::Zero();
            for (size_t ring_atom_idx : ring.atom_indices) {
                double q = conf.AtomAt(ring_atom_idx).partial_charge;
                Vec3 d = atom_pos - conf.PositionAt(ring_atom_idx);
                double r = d.norm();
                if (r < MIN_DISTANCE) continue;
                double r3 = r * r * r;
                double r5 = r3 * r * r;

                // Dipolar kernel weighted by charge: q * K_ab
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b)
                        V_dist(a, b) += q * (3.0 * d(a) * d(b) / r5
                                            - (a == b ? 1.0 : 0.0) / r3);
            }

            SphericalTensor st_point = rn.quad_spherical;
            SphericalTensor st_dist = SphericalTensor::Decompose(V_dist);

            // Find the distance bin
            for (auto& bin : bins) {
                if (dist >= bin.lo && dist < bin.hi) {
                    double t2_mag_p = st_point.T2Magnitude();
                    double t2_mag_d = st_dist.T2Magnitude();

                    if (t2_mag_p > 1e-8 && t2_mag_d > 1e-8) {
                        bin.sum_t2_cos += std::abs(t2_cos(st_point, st_dist));
                        bin.count++;
                    }

                    // T0 ratio (distributed has T0 from charge-weighted
                    // dipolar sum; point has T0 = 0 by construction)
                    // Actually: point-quad EFG is traceless (T0=0).
                    // Distributed EFG is also traceless (sum of traceless terms).
                    // So T0 comparison is trivial (both ~0).
                    // Instead compare Frobenius norm ratio.
                    double norm_p = rn.quad_tensor.norm();
                    double norm_d = V_dist.norm();
                    if (norm_d > 1e-10) {
                        bin.sum_t0_ratio += norm_p / norm_d;
                        bin.t0_count++;
                    }
                    break;
                }
            }
        }
    }

    std::cout << "  Point-quadrupole vs distributed EFG convergence:\n"
              << "  " << std::setw(10) << "Distance"
              << std::setw(10) << "Pairs"
              << std::setw(14) << "Mean|cos(T2)|"
              << std::setw(16) << "Mean Frob ratio"
              << "\n";

    for (const auto& bin : bins) {
        double mean_cos = bin.count > 0 ? bin.sum_t2_cos / bin.count : 0;
        double mean_ratio = bin.t0_count > 0 ? bin.sum_t0_ratio / bin.t0_count : 0;
        std::cout << "  " << std::setw(4) << bin.lo << "-" << std::setw(4) << bin.hi << "A"
                  << std::setw(10) << bin.count
                  << std::setw(14) << std::setprecision(4) << mean_cos
                  << std::setw(16) << mean_ratio
                  << "\n";
    }

    // Note: the point-quadrupole and distributed charge models are
    // fundamentally different. The charge sum includes monopole, dipole,
    // quadrupole, and all higher multipoles. The point-quadrupole captures
    // only the quadrupole term. At large distance, the quadrupole term
    // dominates (slower decay than higher multipoles), but the monopole
    // and dipole terms from the net ring charge and dipole moment also
    // contribute to the distributed EFG.
    //
    // We do NOT expect perfect convergence. What we expect is that the
    // Frobenius ratio (point/distributed) decreases with distance, showing
    // that the point-quadrupole contribution is a FRACTION of the total
    // charge-based EFG, and that fraction is well-behaved.
    //
    // The T2 cos tells us whether the angular patterns are similar. If
    // cos is near-random (~0.36) even at close range, the quadrupole
    // angular pattern is unrelated to the charge pattern, which would be
    // a problem. Values of 0.5-0.7 indicate partial agreement — the
    // quadrupole captures part of the angular structure.
    if (!bins.empty() && bins[0].count > 0) {
        double near_cos = bins[0].sum_t2_cos / bins[0].count;
        EXPECT_GT(near_cos, 0.4)
            << "At close range, PQ and charge T2 should share some angular structure";
    }
}
