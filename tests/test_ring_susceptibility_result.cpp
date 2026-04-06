#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

#include "RingSusceptibilityResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "PhysicalConstants.h"

#include <filesystem>
namespace fs = std::filesystem;
using namespace nmr;



// ============================================================================
// Analytical test: atom at known position relative to ring center,
// hand-calculable result.
//
// Same formula as McConnell analytical test (lesson 19) with b_hat → n_hat.
// Ring center at origin, normal along z. Atom at (3, 0, 0).
//
// d = (3, 0, 0), r = 3, d_hat = (1, 0, 0)
// n_hat = (0, 0, 1)
// cos_theta = d_hat . n_hat = 0
//
// Scalar f = (3*0 - 1) / 27 = -1/27
//
// Full tensor M / r³:
//   Term 1: 9*0 * d_hat ⊗ n_hat / r³ = 0 (cos_theta = 0)
//   Term 2: -3 * n_hat ⊗ n_hat / r³ → only (2,2) = -3/27 = -1/9
//   Term 3: -(3 d_hat ⊗ d_hat - I) / r³ = -K
//
//   M_11 / r³ = 0 + 0 + (-2/27) = -2/27
//   M_22 / r³ = 0 + 0 + (1/27)  = 1/27
//   M_33 / r³ = 0 + (-3/27) + (1/27) = -2/27
//
//   Trace = -2/27 + 1/27 + (-2/27) = -3/27 = -1/9
//   T0 = Trace/3 = -1/27 = f  ✓
// ============================================================================

TEST(RingChiAnalytical, T0EqualsFOnRealProtein) {
    // Verify T0 = f identity on a real protein. This is the same identity
    // as McConnell (Trace(M)/3 = (3cos²θ-1)/r³) applied to the ring
    // susceptibility tensor with n_hat instead of b_hat.
    //
    // We verify by recomputing M from scratch for each ring-atom pair
    // and checking Trace(M/r³)/3 == f at machine precision.

    if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) GTEST_SKIP() << "1UBQ not found";
    auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    if (!r.Ok()) GTEST_SKIP() << r.error;

    auto& conf = r.protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(SpatialIndexResult::Compute(conf));
    conf.AttachResult(RingSusceptibilityResult::Compute(conf));

    ASSERT_GT(r.protein->RingCount(), 0u) << "1UBQ should have rings";

    int checked = 0;
    double max_diff = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (std::abs(rn.chi_scalar) < 1e-15) continue;

            // Reconstruct M/r³ from geometry and verify trace
            const RingGeometry& geom = conf.ring_geometries[rn.ring_index];
            Vec3 d = conf.PositionAt(ai) - geom.center;
            double rm = d.norm();
            if (rm < MIN_DISTANCE) continue;

            Vec3 d_hat = d / rm;
            double cos_theta = d_hat.dot(geom.normal);
            double r3 = rm * rm * rm;

            Mat3 M;
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    M(a, b) = (9.0 * cos_theta * d_hat(a) * geom.normal(b)
                               - 3.0 * geom.normal(a) * geom.normal(b)
                               - (3.0 * d_hat(a) * d_hat(b) - (a==b ? 1.0 : 0.0)))
                              / r3;

            double t0 = M.trace() / 3.0;
            double f = (3.0 * cos_theta * cos_theta - 1.0) / r3;
            double diff = std::abs(t0 - f);
            max_diff = std::max(max_diff, diff);

            // Also verify stored values match
            EXPECT_NEAR(rn.chi_scalar, f, 1e-10);
            EXPECT_NEAR(rn.chi_spherical.T0, f, 1e-10);
            checked++;
        }
    }

    EXPECT_GT(checked, 100) << "Should verify many ring-atom pairs";
    EXPECT_LT(max_diff, 1e-10) << "T0 must equal f at machine precision";

    std::cout << "  Verified T0 = f on " << checked
              << " ring-atom pairs, max |T0-f| = " << max_diff << "\n";
}


// ============================================================================
// Full protein test
// ============================================================================

class RingChiProteinTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ.pdb not found";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << "Failed to load 1UBQ";
        protein = std::move(r.protein);

        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        conf.AttachResult(SpatialIndexResult::Compute(conf));
    }

    std::unique_ptr<Protein> protein;
};


TEST_F(RingChiProteinTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto rchi = RingSusceptibilityResult::Compute(conf);
    ASSERT_NE(rchi, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(rchi)));
    ASSERT_TRUE(conf.HasResult<RingSusceptibilityResult>());
}


TEST_F(RingChiProteinTest, T0EqualsScalarF) {
    // For each atom-ring pair, T0 of the full tensor M/r³ should equal f.
    auto& conf = protein->Conformation();
    conf.AttachResult(RingSusceptibilityResult::Compute(conf));

    int checked = 0;
    double max_diff = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (std::abs(rn.chi_scalar) < 1e-15) continue;
            double diff = std::abs(rn.chi_spherical.T0 - rn.chi_scalar);
            max_diff = std::max(max_diff, diff);
            checked++;
        }
    }

    EXPECT_GT(checked, 0) << "Should have checked some ring-atom pairs";
    EXPECT_LT(max_diff, 1e-10)
        << "T0 of full tensor must equal scalar f";

    std::cout << "  Checked " << checked << " pairs, max |T0 - f| = "
              << max_diff << "\n";
}


TEST_F(RingChiProteinTest, ShieldingContributionHasT0AndT2) {
    auto& conf = protein->Conformation();
    conf.AttachResult(RingSusceptibilityResult::Compute(conf));

    int nonzero_t0 = 0, nonzero_t2 = 0;
    double max_t0 = 0, max_t2 = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& sc = conf.AtomAt(ai).ringchi_shielding_contribution;
        if (std::abs(sc.T0) > 1e-8) nonzero_t0++;
        double t2 = sc.T2Magnitude();
        if (t2 > 1e-8) nonzero_t2++;
        max_t0 = std::max(max_t0, std::abs(sc.T0));
        max_t2 = std::max(max_t2, t2);
    }

    EXPECT_GT(nonzero_t0, 0) << "Full tensor must have non-zero T0";
    EXPECT_GT(nonzero_t2, 0) << "Full tensor must have non-zero T2";

    std::cout << "  T0 nonzero: " << nonzero_t0 << ", max |T0| = " << max_t0 << "\n";
    std::cout << "  T2 nonzero: " << nonzero_t2 << ", max |T2| = " << max_t2 << "\n";
}


// ============================================================================
// ORCA protein test
// ============================================================================

TEST(RingChiOrcaTest, RunOnProtonatedProtein) {
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

    auto rchi = RingSusceptibilityResult::Compute(conf);
    ASSERT_NE(rchi, nullptr);
    conf.AttachResult(std::move(rchi));

    // Summary
    double max_t0 = 0, max_t2 = 0;
    int with_rings = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& sc = conf.AtomAt(ai).ringchi_shielding_contribution;
        max_t0 = std::max(max_t0, std::abs(sc.T0));
        max_t2 = std::max(max_t2, sc.T2Magnitude());
        if (!conf.AtomAt(ai).ring_neighbours.empty()) with_rings++;
    }

    std::cout << "  ORCA protein RingSusceptibility summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " rings=" << load.protein->RingCount() << "\n"
              << "    atoms with ring neighbours: " << with_rings << "\n"
              << "    max |T0| = " << max_t0 << " A^-3\n"
              << "    max |T2| = " << max_t2 << " A^-3\n";

    EXPECT_GT(with_rings, 0) << "Some atoms should see rings";
    EXPECT_GT(max_t0, 0.001) << "T0 should be non-zero";
    EXPECT_GT(max_t2, 0.001) << "T2 should be non-zero";
}
