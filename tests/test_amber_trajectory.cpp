// PHASE 0 trajectory loader tests — AMBER-ff GROMACS straight-relaxation
// MD trajectories. Verifies that the existing libgromacs C++ path
// (FullSystemReader::ReadTopology + GromacsFrameHandler) handles AMBER
// TPRs cleanly when the force-field tag is set explicitly.
//
// GROMACS interaction is via the linked C++ API (read_tpx_state /
// xdrfile) — no subprocess `gmx dump` or `gmx trjconv` calls. See the
// "GROMACS doctrine" note in
// spec/plan/amber-implementation-plan-2026-04-29.md.
//
// Tests SKIP when fixtures are absent (the AMBER trajectory tree is
// gitignored at tests/data/fleet_amber/ and generated via the fleet
// preparation pipeline).
//

#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "BuildResult.h"
#include "ChargeSource.h"
#include "ForceFieldChargeTable.h"
#include "FullSystemReader.h"
#include "Protein.h"
#include "ProteinBuildContext.h"
#include "ProteinConformation.h"

#include <cmath>
#include <filesystem>

using namespace nmr;

class AmberTrajectoryFixtureTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (nmr::test::TestEnvironment::FleetAmberData().empty()) {
            GTEST_SKIP() << "fleet_amber not configured in testpaths.toml";
        }
    }

    nmr::test::AmberTrajectoryFixture FixtureFor(
            const std::string& protein_id) {
        auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(protein_id);
        if (fix.tpr_path.empty()) {
            ADD_FAILURE() << "no fleet_amber subpath for " << protein_id;
        }
        return fix;
    }
};


TEST_F(AmberTrajectoryFixtureTest, FixturesPathsResolveToRealFiles_1P9J) {
    auto fix = FixtureFor("1P9J_5801");
    if (fix.tpr_path.empty()) GTEST_SKIP();
    EXPECT_TRUE(std::filesystem::exists(fix.tpr_path)) << fix.tpr_path;
    EXPECT_TRUE(std::filesystem::exists(fix.xtc_path)) << fix.xtc_path;
    EXPECT_TRUE(std::filesystem::exists(fix.edr_path)) << fix.edr_path;
}

TEST_F(AmberTrajectoryFixtureTest, FixturesPathsResolveToRealFiles_1Z9B) {
    auto fix = FixtureFor("1Z9B_6577");
    if (fix.tpr_path.empty()) GTEST_SKIP();
    EXPECT_TRUE(std::filesystem::exists(fix.tpr_path)) << fix.tpr_path;
    EXPECT_TRUE(std::filesystem::exists(fix.xtc_path)) << fix.xtc_path;
    EXPECT_TRUE(std::filesystem::exists(fix.edr_path)) << fix.edr_path;
}


// First real test: libgromacs's read_tpx_state must accept the AMBER
// TPR end-to-end through FullSystemReader. The TPR format itself is
// force-field-agnostic; this confirms there's no AMBER-specific quirk
// in our parsing path.
TEST_F(AmberTrajectoryFixtureTest, FullSystemReaderReadsAmberTpr_1P9J) {
    auto fix = FixtureFor("1P9J_5801");
    if (fix.tpr_path.empty()) GTEST_SKIP();
    if (!std::filesystem::exists(fix.tpr_path)) {
        GTEST_SKIP() << "tpr not on disk: " << fix.tpr_path;
    }

    FullSystemReader reader;
    ASSERT_TRUE(reader.ReadTopology(fix.tpr_path)) << reader.error();

    const auto& topo = reader.Topology();
    EXPECT_GT(topo.protein_count, 0u);
    EXPECT_GT(topo.water_count, 0u);  // solvated production
}


TEST_F(AmberTrajectoryFixtureTest, FullSystemReaderBuildsAmberProtein_1P9J) {
    auto fix = FixtureFor("1P9J_5801");
    if (fix.tpr_path.empty()) GTEST_SKIP();
    if (!std::filesystem::exists(fix.tpr_path)) GTEST_SKIP();

    FullSystemReader reader;
    ASSERT_TRUE(reader.ReadTopology(fix.tpr_path)) << reader.error();

    auto build = reader.BuildProtein(fix.protein_id, ForceField::Amber_ff14SB);
    ASSERT_TRUE(build.Ok()) << build.error;

    EXPECT_GT(build.protein->AtomCount(), 100u);
    EXPECT_GT(build.protein->ResidueCount(), 10u);

    // Charges came from the TPR via the libgromacs C++ path; the typed
    // table must record AMBER ff14SB as the source force field.
    ASSERT_TRUE(build.protein->HasForceFieldCharges());
    const auto& table = build.protein->ForceFieldCharges();
    EXPECT_EQ(table.SourceForceField(), ForceField::Amber_ff14SB);
    EXPECT_EQ(table.AtomCount(), build.protein->AtomCount());

    // Net charge should be near-integer (whole-system MD with neutralised box).
    double total = table.TotalCharge();
    double frac = std::abs(total - std::round(total));
    EXPECT_LT(frac, 0.05) << "Total charge " << total
        << " not near-integer (NetCharge integrity check)";
}


TEST_F(AmberTrajectoryFixtureTest, FullSystemReaderBuildsAmberProtein_1Z9B) {
    auto fix = FixtureFor("1Z9B_6577");
    if (fix.tpr_path.empty()) GTEST_SKIP();
    if (!std::filesystem::exists(fix.tpr_path)) GTEST_SKIP();

    FullSystemReader reader;
    ASSERT_TRUE(reader.ReadTopology(fix.tpr_path)) << reader.error();

    auto build = reader.BuildProtein(fix.protein_id, ForceField::Amber_ff14SB);
    ASSERT_TRUE(build.Ok()) << build.error;

    EXPECT_GT(build.protein->AtomCount(), 100u);
    EXPECT_GT(build.protein->ResidueCount(), 10u);
    ASSERT_TRUE(build.protein->HasForceFieldCharges());
    EXPECT_EQ(build.protein->ForceFieldCharges().SourceForceField(),
              ForceField::Amber_ff14SB);

    double total = build.protein->ForceFieldCharges().TotalCharge();
    double frac = std::abs(total - std::round(total));
    EXPECT_LT(frac, 0.05) << "Total charge " << total
        << " not near-integer (NetCharge integrity check)";
}


// 1P9J and 1Z9B are distinct proteins. A test that happens to match
// one fixture's numerics by coincidence cannot also match the other.
TEST_F(AmberTrajectoryFixtureTest, ProteinsDifferBetweenFixtures) {
    auto fix1 = FixtureFor("1P9J_5801");
    auto fix2 = FixtureFor("1Z9B_6577");
    if (fix1.tpr_path.empty() || fix2.tpr_path.empty()) GTEST_SKIP();
    if (!std::filesystem::exists(fix1.tpr_path) ||
            !std::filesystem::exists(fix2.tpr_path)) {
        GTEST_SKIP() << "fixtures not on disk";
    }

    FullSystemReader r1, r2;
    ASSERT_TRUE(r1.ReadTopology(fix1.tpr_path)) << r1.error();
    ASSERT_TRUE(r2.ReadTopology(fix2.tpr_path)) << r2.error();

    auto b1 = r1.BuildProtein(fix1.protein_id, ForceField::Amber_ff14SB);
    auto b2 = r2.BuildProtein(fix2.protein_id, ForceField::Amber_ff14SB);
    ASSERT_TRUE(b1.Ok());
    ASSERT_TRUE(b2.Ok());

    EXPECT_NE(b1.protein->AtomCount(), b2.protein->AtomCount());
    EXPECT_NE(b1.protein->ResidueCount(), b2.protein->ResidueCount());
}
