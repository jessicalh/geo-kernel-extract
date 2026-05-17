//
// test_hbond_count_welford: discipline + integration for
// HBondCountWelfordTrajectoryResult. AV-pattern scalar Welford on
// per-atom H-bond pair count.
//

#include "HBondResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "HBondCountWelfordTrajectoryResult.h"
#include "TestConfig.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryAtom.h"
#include "TrajectoryProtein.h"
#include "Types.h"

#include <gtest/gtest.h>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

namespace {

constexpr const char* kFixtureProtein = "1P9J_5801";

std::string TrrPathFor(const std::string& p) { return fs::path(p).replace_extension(".trr").string(); }
std::string ProductionDirFor(const std::string& p) { return fs::path(p).parent_path().string(); }
bool FixtureAvailable(const nmr::test::AmberTrajectoryFixture& fix) {
    return !fix.tpr_path.empty() && fs::exists(fix.tpr_path)
        && fs::exists(TrrPathFor(fix.tpr_path)) && fs::exists(fix.edr_path);
}
void LoadCalculatorConfig() {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(std::string(NMR_TEST_DATA_DIR) + "/../data/calculator_params.toml");
}

// Per-file bridge to the shared TestConfig profile system. HBondResult
// requires DsspResult; KernelWithDssp profile is the right pick. See
// spec/plan/test-suite-realignment-deferred-2026-05-17.md.
void BuildBaseConfig(nmr::RunConfiguration& config, const std::string& name, std::size_t stride) {
    config = nmr::test::BuildTestConfig(nmr::test::TestProfile::KernelWithDssp, name, stride);
    config.RequireConformationResult(typeid(nmr::HBondResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::HBondCountWelfordTrajectoryResult::Create(tp_in);
        });
}

}  // namespace


TEST(HBondCountWelford, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    BuildBaseConfig(config, "HBondCountWelfordFrame0Semantics", 99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    ASSERT_EQ(traj.FrameCount(), 1u);

    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& ta = tp.AtomAt(i);
        EXPECT_EQ(ta.hbond_count_welford.n_frames, 1u);
        EXPECT_DOUBLE_EQ(ta.hbond_count_welford.count.std, 0.0);
        EXPECT_EQ(ta.hbond_count_welford.delta_n, 0u);
    }
}


TEST(HBondCountWelford, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    BuildBaseConfig(config, "HBondCountWelfordFinalizeIdempotency", 99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::HBondCountWelfordTrajectoryResult>();
    const size_t probe = tp.AtomCount() / 2;
    const double mean_first = tp.AtomAt(probe).hbond_count_welford.count.mean;
    const double std_first  = tp.AtomAt(probe).hbond_count_welford.count.std;

    tr.Finalize(tp, traj);

    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).hbond_count_welford.count.mean, mean_first);
    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).hbond_count_welford.count.std,  std_first);
}


TEST(HBondCountWelford, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    BuildBaseConfig(config, "HBondCountWelfordH5RoundTrip", 99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::HBondCountWelfordTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("hbond_count_welford_h5_roundtrip_" + std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr.WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/hbond_count_welford"));
    auto grp = reopen.getGroup("/trajectory/hbond_count_welford");

    std::string units;
    grp.getAttribute("units").read(units);
    EXPECT_EQ(units, "pairs");

    double radius;
    grp.getAttribute("source_radius_A").read(radius);
    EXPECT_DOUBLE_EQ(radius, 3.5);

    const auto dims = grp.getDataSet("mean").getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 1u);
    EXPECT_EQ(dims[0], tp.AtomCount());

    fs::remove(h5_path);
}


TEST(HBondCountWelford, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    BuildBaseConfig(config, "HBondCountWelfordIntegration", 300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_GE(traj.FrameCount(), 2u);

    size_t populated = 0;
    double max_mean = 0.0;
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& ta = tp.AtomAt(i);
        EXPECT_TRUE(std::isfinite(ta.hbond_count_welford.count.mean));
        EXPECT_TRUE(std::isfinite(ta.hbond_count_welford.count.std));
        EXPECT_GE(ta.hbond_count_welford.count.mean, 0.0);  // count is non-negative
        EXPECT_EQ(ta.hbond_count_welford.n_frames, traj.FrameCount());
        if (ta.hbond_count_welford.count.mean > 1e-12) ++populated;
        max_mean = std::max(max_mean, ta.hbond_count_welford.count.mean);
    }
    // Catches silent-zero wiring bugs (e.g. DSSP starved → HBondResult
    // ran with empty hbond set → every atom got 0 counts every frame).
    EXPECT_GT(populated, 0u)
        << "HBondCount all-zero — HBondResult or DsspResult not firing";
    std::cout << "HBondCountWelford integration: populated=" << populated
              << "/" << tp.AtomCount()
              << " max mean=" << max_mean << " pairs\n";
}
