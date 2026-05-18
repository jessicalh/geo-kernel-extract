//
// test_bs_welford: discipline + integration for BsWelfordTrajectoryResult.
// AV-pattern Welford on the Biot-Savart shielding kernel. BsWelford is
// the exemplar TR for the Welford family (Phase 2b refactor 2026-05-17):
// every other Welford TR clones this shape. It is also the upstream
// writer for BsAnomalousAtomMarker's cross-result read of
// `ta.bs_welford.{t0.mean, t0.m2, n_frames}` (PATTERNS §17).
//
// Discipline trio (Frame0Semantics + FinalizeIdempotency + H5RoundTrip)
// + Integration1P9J. The integration test asserts T0 + per-component
// T1/T2 channel coverage to protect against silent regression of the
// per-component channels (T2 Completeness; analogous to the McConnell
// T1 omission caught on 2026-05-17).
//

#include "BiotSavartResult.h"

#include "BsWelfordTrajectoryResult.h"
#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "TestConfig.h"
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
// DISCIPLINE: Frame-0 semantics — single-frame fixture exercises Compute once.
// ============================================================================

TEST(BsWelford, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "BsWelfordFrame0Semantics", 99999);
    config.RequireConformationResult(typeid(nmr::BiotSavartResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::BsWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    ASSERT_EQ(traj.FrameCount(), 1u);

    ASSERT_TRUE(tp.HasResult<nmr::BsWelfordTrajectoryResult>());
    const auto& tr = tp.Result<nmr::BsWelfordTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);

    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& ta = tp.AtomAt(i);
        EXPECT_EQ(ta.bs_welford.n_frames, 1u);
        EXPECT_DOUBLE_EQ(ta.bs_welford.t0.std, 0.0);
        EXPECT_EQ(ta.bs_welford.delta_n, 0u);
    }
}


// ============================================================================
// DISCIPLINE: Finalize idempotency — second Finalize must not corrupt state.
// ============================================================================

TEST(BsWelford, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "BsWelfordFinalizeIdempotency", 99999);
    config.RequireConformationResult(typeid(nmr::BiotSavartResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::BsWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::BsWelfordTrajectoryResult>();
    const size_t probe = tp.AtomCount() / 2;
    const double mean_first = tp.AtomAt(probe).bs_welford.t0.mean;
    const double std_first  = tp.AtomAt(probe).bs_welford.t0.std;

    tr.Finalize(tp, traj);

    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).bs_welford.t0.mean, mean_first);
    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).bs_welford.t0.std,  std_first);
}


// ============================================================================
// DISCIPLINE: H5 round-trip — schema attributes + dataset shape survive.
// ============================================================================

TEST(BsWelford, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "BsWelfordH5RoundTrip", 99999);
    config.RequireConformationResult(typeid(nmr::BiotSavartResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::BsWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::BsWelfordTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("bs_welford_h5_roundtrip_" + std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr.WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/bs_welford"));
    auto grp = reopen.getGroup("/trajectory/bs_welford");

    std::string units;
    grp.getAttribute("units").read(units);
    EXPECT_EQ(units, "ppm_T_per_nA");

    // T1 storage convention: Cartesian Levi-Civita dual, NOT real-Y_1m.
    std::string irrep_t1;
    grp.getAttribute("irrep_layout_t1").read(irrep_t1);
    EXPECT_EQ(irrep_t1, "v_x,v_y,v_z");

    const auto dims_t0 = grp.getDataSet("t0_mean").getSpace().getDimensions();
    ASSERT_EQ(dims_t0.size(), 1u);
    EXPECT_EQ(dims_t0[0], tp.AtomCount());

    // Per-component T1/T2 shapes: 2D (N, K).
    const auto dims_t1 = grp.getDataSet("t1_mean").getSpace().getDimensions();
    ASSERT_EQ(dims_t1.size(), 2u);
    EXPECT_EQ(dims_t1[0], tp.AtomCount());
    EXPECT_EQ(dims_t1[1], 3u);

    const auto dims_t2 = grp.getDataSet("t2_mean").getSpace().getDimensions();
    ASSERT_EQ(dims_t2.size(), 2u);
    EXPECT_EQ(dims_t2[0], tp.AtomCount());
    EXPECT_EQ(dims_t2[1], 5u);

    fs::remove(h5_path);
}


// ============================================================================
// INTEGRATION: real BS kernel through Trajectory::Run, multi-frame.
// ============================================================================

TEST(BsWelford, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "BsWelfordIntegration", 300);
    config.RequireConformationResult(typeid(nmr::BiotSavartResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::BsWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_GE(traj.FrameCount(), 2u);

    // T0 + per-component T1/T2 coverage. The T1/T2 floors protect against
    // silent regression of the per-component channels (T2 Completeness;
    // analogous protection to the McConnell T1 omission caught 2026-05-17).
    size_t populated  = 0;
    size_t t1_populated = 0, t2_populated = 0;
    double max_abs_t0 = 0.0, max_abs_t1 = 0.0, max_abs_t2 = 0.0;
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& ta = tp.AtomAt(i);
        EXPECT_TRUE(std::isfinite(ta.bs_welford.t0.mean));
        EXPECT_TRUE(std::isfinite(ta.bs_welford.t0.std));
        EXPECT_TRUE(std::isfinite(ta.bs_welford.t2magnitude.mean));
        EXPECT_EQ(ta.bs_welford.n_frames, traj.FrameCount());
        if (std::abs(ta.bs_welford.t0.mean) > 1e-12) ++populated;
        max_abs_t0 = std::max(max_abs_t0, std::abs(ta.bs_welford.t0.mean));

        bool atom_t1 = false;
        for (size_t k = 0; k < 3; ++k) {
            EXPECT_TRUE(std::isfinite(ta.bs_welford.t1[k].mean));
            const double v = std::abs(ta.bs_welford.t1[k].mean);
            max_abs_t1 = std::max(max_abs_t1, v);
            if (v > 1e-12) atom_t1 = true;
        }
        if (atom_t1) ++t1_populated;

        bool atom_t2 = false;
        for (size_t k = 0; k < 5; ++k) {
            EXPECT_TRUE(std::isfinite(ta.bs_welford.t2[k].mean));
            const double v = std::abs(ta.bs_welford.t2[k].mean);
            max_abs_t2 = std::max(max_abs_t2, v);
            if (v > 1e-12) atom_t2 = true;
        }
        if (atom_t2) ++t2_populated;
    }
    EXPECT_GT(populated, 0u)
        << "BS Welford T0 all-zero — calculator regression or attach miss";
    EXPECT_GT(t1_populated, 0u)
        << "BS Welford T1 all-zero — per-component regression";
    EXPECT_GT(t2_populated, 0u)
        << "BS Welford T2 per-component all-zero — T2 Completeness regression";
    std::cout << "BsWelford integration: populated=" << populated
              << "/" << tp.AtomCount()
              << " max|t0|=" << max_abs_t0
              << " t1_pop=" << t1_populated << " max|t1|=" << max_abs_t1
              << " t2_pop=" << t2_populated << " max|t2|=" << max_abs_t2
              << " (ppm_T_per_nA)\n";
}
