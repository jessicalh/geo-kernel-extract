//
// test_mc_welford: discipline + integration for McConnellWelfordTrajectoryResult.
// AV-pattern Welford on the McConnell shielding kernel. Mirrors the BS
// Welford test shape — no synthetic test; the shared MomentsUpdate math
// is exercised by every AV TR. Discipline trio + Integration1P9J.
//

#include "McConnellResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "McConnellWelfordTrajectoryResult.h"
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


std::string TrrPathFor(const std::string& tpr_path) {
    return fs::path(tpr_path).replace_extension(".trr").string();
}

std::string ProductionDirFor(const std::string& tpr_path) {
    return fs::path(tpr_path).parent_path().string();
}

bool FixtureAvailable(const nmr::test::AmberTrajectoryFixture& fix) {
    return !fix.tpr_path.empty() &&
           fs::exists(fix.tpr_path) &&
           fs::exists(TrrPathFor(fix.tpr_path)) &&
           fs::exists(fix.edr_path);
}

void LoadCalculatorConfig() {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../data/calculator_params.toml");
}

}  // namespace


// ============================================================================
// DISCIPLINE: Frame-0 semantics.
// ============================================================================

TEST(McConnellWelford, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "McConnellWelfordFrame0Semantics", 99999);
    config.RequireConformationResult(typeid(nmr::McConnellResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::McConnellWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    ASSERT_EQ(traj.FrameCount(), 1u);

    ASSERT_TRUE(tp.HasResult<nmr::McConnellWelfordTrajectoryResult>());
    const auto& tr = tp.Result<nmr::McConnellWelfordTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);

    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& ta = tp.AtomAt(i);
        EXPECT_EQ(ta.mc_welford.n_frames, 1u);
        EXPECT_DOUBLE_EQ(ta.mc_welford.t0.std, 0.0);
        EXPECT_EQ(ta.mc_welford.delta_n, 0u);
    }
}


// ============================================================================
// DISCIPLINE: Finalize idempotency.
// ============================================================================

TEST(McConnellWelford, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "McConnellWelfordFinalizeIdempotency", 99999);
    config.RequireConformationResult(typeid(nmr::McConnellResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::McConnellWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::McConnellWelfordTrajectoryResult>();
    const size_t probe = tp.AtomCount() / 2;
    const double mean_first = tp.AtomAt(probe).mc_welford.t0.mean;
    const double std_first  = tp.AtomAt(probe).mc_welford.t0.std;

    tr.Finalize(tp, traj);

    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).mc_welford.t0.mean, mean_first);
    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).mc_welford.t0.std,  std_first);
}


// ============================================================================
// DISCIPLINE: H5 round-trip.
// ============================================================================

TEST(McConnellWelford, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "McConnellWelfordH5RoundTrip", 99999);
    config.RequireConformationResult(typeid(nmr::McConnellResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::McConnellWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::McConnellWelfordTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("mc_welford_h5_roundtrip_" + std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr.WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/mc_welford"));
    auto grp = reopen.getGroup("/trajectory/mc_welford");

    std::string units;
    grp.getAttribute("units").read(units);
    EXPECT_EQ(units, "Angstrom^-3");

    const auto dims = grp.getDataSet("t0_mean").getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 1u);
    EXPECT_EQ(dims[0], tp.AtomCount());

    fs::remove(h5_path);
}


// ============================================================================
// INTEGRATION: real McConnell kernel through Trajectory::Run.
// ============================================================================

TEST(McConnellWelford, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "McConnellWelfordIntegration", 300);
    config.RequireConformationResult(typeid(nmr::McConnellResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::McConnellWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_GE(traj.FrameCount(), 2u);

    size_t populated = 0;
    double max_abs_t0 = 0.0;
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& ta = tp.AtomAt(i);
        EXPECT_TRUE(std::isfinite(ta.mc_welford.t0.mean));
        EXPECT_TRUE(std::isfinite(ta.mc_welford.t0.std));
        EXPECT_TRUE(std::isfinite(ta.mc_welford.t2magnitude.mean));
        EXPECT_EQ(ta.mc_welford.n_frames, traj.FrameCount());
        if (std::abs(ta.mc_welford.t0.mean) > 1e-12) ++populated;
        max_abs_t0 = std::max(max_abs_t0, std::abs(ta.mc_welford.t0.mean));
    }
    EXPECT_GT(populated, 0u)
        << "McConnell Welford all-zero — calculator regression or attach miss";
    std::cout << "McConnellWelford integration: populated=" << populated
              << "/" << tp.AtomCount()
              << " max|t0_mean|=" << max_abs_t0 << " A^-3\n";
}
