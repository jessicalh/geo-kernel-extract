//
// test_hm_welford: discipline + integration for HmWelfordTrajectoryResult.
// AV-pattern Welford on the HaighMallion shielding kernel. Mirrors the
// BS Welford test shape — no synthetic test (the shared MomentsUpdate
// math is exercised by every AV TR's integration path; pinning it again
// per TR is redundant). Discipline trio + Integration1P9J only.
//

#include "HaighMallionResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "HmWelfordTrajectoryResult.h"
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
// DISCIPLINE: Frame-0 semantics — stride ≥ fixture length → one Compute call.
// ============================================================================

TEST(HmWelford, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "HmWelfordFrame0Semantics", 99999);
    config.RequireConformationResult(typeid(nmr::HaighMallionResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::HmWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    ASSERT_EQ(traj.FrameCount(), 1u);

    ASSERT_TRUE(tp.HasResult<nmr::HmWelfordTrajectoryResult>());
    const auto& tr = tp.Result<nmr::HmWelfordTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);

    // Single-frame sanity: n_frames == 1 for every atom, std == 0 (no spread),
    // delta_n == 0 (no prior frame).
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& ta = tp.AtomAt(i);
        EXPECT_EQ(ta.hm_welford.n_frames, 1u);
        EXPECT_DOUBLE_EQ(ta.hm_welford.t0.std, 0.0);
        EXPECT_EQ(ta.hm_welford.delta_n, 0u);
    }
}


// ============================================================================
// DISCIPLINE: Finalize idempotency.
// ============================================================================

TEST(HmWelford, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "HmWelfordFinalizeIdempotency", 99999);
    config.RequireConformationResult(typeid(nmr::HaighMallionResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::HmWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::HmWelfordTrajectoryResult>();
    const size_t probe = tp.AtomCount() / 2;
    const double mean_first = tp.AtomAt(probe).hm_welford.t0.mean;
    const double std_first  = tp.AtomAt(probe).hm_welford.t0.std;

    tr.Finalize(tp, traj);

    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).hm_welford.t0.mean, mean_first);
    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).hm_welford.t0.std,  std_first);
}


// ============================================================================
// DISCIPLINE: H5 round-trip.
// ============================================================================

TEST(HmWelford, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "HmWelfordH5RoundTrip", 99999);
    config.RequireConformationResult(typeid(nmr::HaighMallionResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::HmWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::HmWelfordTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("hm_welford_h5_roundtrip_" + std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr.WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/hm_welford"));
    auto grp = reopen.getGroup("/trajectory/hm_welford");

    std::string units;
    grp.getAttribute("units").read(units);
    EXPECT_EQ(units, "Angstrom^-1");

    const auto dims = grp.getDataSet("t0_mean").getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 1u);
    EXPECT_EQ(dims[0], tp.AtomCount());

    // Keep file at fixed location for manual inspection
    std::cout << "INSPECT_DEBUG: copying " << h5_path << " to /tmp/hm_welford_inspect.h5\n";
    try {
        fs::copy_file(h5_path, "/tmp/hm_welford_inspect.h5",
                      fs::copy_options::overwrite_existing);
        std::cout << "INSPECT_DEBUG: copy OK\n";
    } catch (const std::exception& e) {
        std::cout << "INSPECT_DEBUG: copy failed: " << e.what() << "\n";
    }
    fs::remove(h5_path);
}


// ============================================================================
// INTEGRATION: real HM kernel through Trajectory::Run.
// ============================================================================

TEST(HmWelford, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "HmWelfordIntegration", 300);
    config.RequireConformationResult(typeid(nmr::HaighMallionResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::HmWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_GE(traj.FrameCount(), 2u);

    // T0 / T1 / T2 channel coverage — protects against silent regression
    // of the per-component T1 / T2 channels (T2 Completeness, PATTERNS).
    // A future refactor that quiets all three per-component T1 channels
    // would not trip the t0 assertion alone.
    size_t populated  = 0;
    size_t t1_populated = 0, t2_populated = 0;
    double max_abs_t0 = 0.0, max_abs_t1 = 0.0, max_abs_t2 = 0.0;
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& ta = tp.AtomAt(i);
        EXPECT_TRUE(std::isfinite(ta.hm_welford.t0.mean));
        EXPECT_TRUE(std::isfinite(ta.hm_welford.t0.std));
        EXPECT_TRUE(std::isfinite(ta.hm_welford.t2magnitude.mean));
        EXPECT_EQ(ta.hm_welford.n_frames, traj.FrameCount());
        if (std::abs(ta.hm_welford.t0.mean) > 1e-12) ++populated;
        max_abs_t0 = std::max(max_abs_t0, std::abs(ta.hm_welford.t0.mean));

        bool atom_t1 = false;
        for (size_t k = 0; k < 3; ++k) {
            EXPECT_TRUE(std::isfinite(ta.hm_welford.t1[k].mean));
            const double v = std::abs(ta.hm_welford.t1[k].mean);
            max_abs_t1 = std::max(max_abs_t1, v);
            if (v > 1e-12) atom_t1 = true;
        }
        if (atom_t1) ++t1_populated;

        bool atom_t2 = false;
        for (size_t k = 0; k < 5; ++k) {
            EXPECT_TRUE(std::isfinite(ta.hm_welford.t2[k].mean));
            const double v = std::abs(ta.hm_welford.t2[k].mean);
            max_abs_t2 = std::max(max_abs_t2, v);
            if (v > 1e-12) atom_t2 = true;
        }
        if (atom_t2) ++t2_populated;
    }
    EXPECT_GT(populated, 0u)
        << "HM Welford T0 all-zero — calculator regression or attach miss";
    EXPECT_GT(t1_populated, 0u)
        << "HM Welford T1 all-zero — per-component regression";
    EXPECT_GT(t2_populated, 0u)
        << "HM Welford T2 per-component all-zero — T2 Completeness regression";

    // Codex 2026-05-18: dxdt_n must equal delta_n on a well-formed
    // trajectory (no duplicated-timestamp frames). Regression guard
    // for the separate-counter fix that skips zero-dt frames rather
    // than zero-filling the rate accumulator.
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& w = tp.AtomAt(i).hm_welford;
        EXPECT_EQ(w.dxdt_n, w.delta_n)
            << "dxdt_n != delta_n on atom " << i
            << " — investigate timestamp duplication in fixture";
    }

    std::cout << "HmWelford integration: populated=" << populated
              << "/" << tp.AtomCount()
              << " max|t0|=" << max_abs_t0
              << " t1_pop=" << t1_populated << " max|t1|=" << max_abs_t1
              << " t2_pop=" << t2_populated << " max|t2|=" << max_abs_t2
              << " (A^-1)\n";
}
