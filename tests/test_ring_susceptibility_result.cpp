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
// Analytical test: verify the kept scalar rescue.
// ============================================================================

TEST(RingChiAnalytical, ScalarMatchesRingContributionColumn7Formula) {
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

            const RingGeometry& geom = conf.ring_geometries[rn.ring_index];
            Vec3 d = conf.PositionAt(ai) - geom.center;
            double rm = d.norm();
            if (rm < MIN_DISTANCE) continue;

            Vec3 d_hat = d / rm;
            double cos_theta = d_hat.dot(geom.normal);
            double r3 = rm * rm * rm;

            double f = (3.0 * cos_theta * cos_theta - 1.0) / r3;
            double col7_formula = 0.0;
            if (rn.distance_to_center > 1e-12) {
                double cos_th = rn.z / rn.distance_to_center;
                double rr3 = rn.distance_to_center * rn.distance_to_center
                           * rn.distance_to_center;
                if (rr3 > 1e-30)
                    col7_formula = (3.0 * cos_th * cos_th - 1.0) / rr3;
            }
            double diff = std::abs(col7_formula - f);
            max_diff = std::max(max_diff, diff);

            EXPECT_NEAR(rn.chi_scalar, f, 1e-10);
            EXPECT_NEAR(col7_formula, f, 1e-10);
            checked++;
        }
    }

    EXPECT_GT(checked, 100) << "Should verify many ring-atom pairs";
    EXPECT_LT(max_diff, 1e-10) << "Column-7 formula must equal f at machine precision";

    std::cout << "  Verified ring-chi scalar rescue on " << checked
              << " ring-atom pairs, max |col7-f| = " << max_diff << "\n";
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


TEST_F(RingChiProteinTest, ScalarMatchesColumn7Formula) {
    auto& conf = protein->Conformation();
    conf.AttachResult(RingSusceptibilityResult::Compute(conf));

    int checked = 0;
    double max_diff = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (std::abs(rn.chi_scalar) < 1e-15) continue;
            double col7_formula = 0.0;
            if (rn.distance_to_center > 1e-12) {
                double cos_th = rn.z / rn.distance_to_center;
                double r3 = rn.distance_to_center * rn.distance_to_center
                          * rn.distance_to_center;
                if (r3 > 1e-30)
                    col7_formula = (3.0 * cos_th * cos_th - 1.0) / r3;
            }
            double diff = std::abs(col7_formula - rn.chi_scalar);
            max_diff = std::max(max_diff, diff);
            checked++;
        }
    }

    EXPECT_GT(checked, 0) << "Should have checked some ring-atom pairs";
    EXPECT_LT(max_diff, 1e-10)
        << "Ring contribution column-7 formula must equal scalar f";

    std::cout << "  Checked " << checked << " pairs, max |col7 - f| = "
              << max_diff << "\n";
}


TEST_F(RingChiProteinTest, ScalarIsPopulated) {
    auto& conf = protein->Conformation();
    conf.AttachResult(RingSusceptibilityResult::Compute(conf));

    int nonzero_scalar = 0;
    double max_scalar = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            const double s = std::abs(rn.chi_scalar);
            if (s > 1e-8) nonzero_scalar++;
            max_scalar = std::max(max_scalar, s);
        }
    }

    EXPECT_GT(nonzero_scalar, 0) << "Ring chi scalar should be non-zero";

    std::cout << "  Scalar nonzero: " << nonzero_scalar
              << ", max |scalar| = " << max_scalar << "\n";
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
    double max_scalar = 0;
    int with_rings = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& atom = conf.AtomAt(ai);
        if (!atom.ring_neighbours.empty()) with_rings++;
        for (const auto& rn : atom.ring_neighbours)
            max_scalar = std::max(max_scalar, std::abs(rn.chi_scalar));
    }

    std::cout << "  ORCA protein RingSusceptibility summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " rings=" << load.protein->RingCount() << "\n"
              << "    atoms with ring neighbours: " << with_rings << "\n"
              << "    max |scalar| = " << max_scalar << " A^-3\n";

    EXPECT_GT(with_rings, 0) << "Some atoms should see rings";
    EXPECT_GT(max_scalar, 0.001) << "Scalar should be non-zero";
}
