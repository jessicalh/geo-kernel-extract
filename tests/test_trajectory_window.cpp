//
// test_trajectory_window: the restored window-local extraction bound on
// Trajectory::Run's single dispatch loop. The window is a [start,
// start+len) bound over raw TRR frame indices; stride still applies WITHIN
// the window; window_len == 0 means "no window" = full trajectory.
//
// Discipline (PATTERNS.md §15): these tests drive Trajectory::Run with a
// narrow RunConfiguration, never a shadow-reimplementation of the
// eight-phase loop. Setup (fixture path, minimal per-frame stack,
// RunConfiguration usage) matches test_rmsd_tracking.cpp verbatim;
// RmsdTrackingTrajectoryResult is the attached TR because it reads
// positions only (Dependencies() == {}) — no AIMNet2/CUDA, minimal noise —
// and it exposes the reference-vs-first-dispatched-frame semantics the
// window most directly affects.
//
// Asserts:
//   (a) every recorded frame_index ∈ [start, start+len);
//   (b) captured frame count == the strided window count;
//   (c) conf0 / RMSD reference == frame `start` (rmsd[0] == 0 and the
//       first recorded frame index == start);
//   (d) a no-window run is unchanged (starts at frame 0, dispatches the
//       whole strided trajectory — strictly more frames than the window).
//

#include "CalculatorConfig.h"
#include "OperationLog.h"
#include "RmsdTrackingTrajectoryResult.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "Types.h"

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

namespace {

// 1P9J's TRR fixture is 1501 frames (10 ps/frame), so a window in the low
// hundreds is comfortably in-bounds.
constexpr const char* kFixtureProtein = "1P9J_5801";

std::string TrrPathFor(const std::string& p) {
    return fs::path(p).replace_extension(".trr").string();
}
std::string ProductionDirFor(const std::string& p) {
    return fs::path(p).parent_path().string();
}
bool FixtureAvailable(const nmr::test::AmberTrajectoryFixture& fix) {
    return !fix.tpr_path.empty() && fs::exists(fix.tpr_path)
        && fs::exists(TrrPathFor(fix.tpr_path)) && fs::exists(fix.edr_path);
}

// Minimal per-frame stack (positions-only TR), optionally windowed. Mirrors
// test_rmsd_tracking.cpp::BuildConfig; adds SetWindow. window_len == 0 means
// no window (the default, byte-for-byte the prior full-trajectory run).
nmr::RunConfiguration BuildConfig(std::size_t stride,
                                  std::size_t window_start,
                                  std::size_t window_len) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    // RmsdTracking reads positions only.
    opts.skip_mopac = true;
    opts.skip_coulomb = true;
    opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::RmsdTrackingTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    config.SetWindow(window_start, window_len);
    return config;
}

// Strided window count: frames dispatched at start, start+stride, ...,
// while < start+len. That is ceil(len / stride).
std::size_t StridedWindowCount(std::size_t len, std::size_t stride) {
    return (len + stride - 1) / stride;
}

}  // namespace


// ── Window bounds the dispatch loop; stride applies within ───────────

TEST(TrajectoryWindow, WindowLocalDispatch) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    // Window = TRR frames [100, 110); stride 2 within → {100,102,104,106,108}.
    constexpr std::size_t kStart = 100;
    constexpr std::size_t kLen = 10;
    constexpr std::size_t kStride = 2;

    auto config = BuildConfig(kStride, kStart, kLen);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& frame_indices = traj.FrameIndices();

    // (b) captured frame count == strided window count.
    EXPECT_EQ(traj.FrameCount(), StridedWindowCount(kLen, kStride));
    EXPECT_EQ(frame_indices.size(), StridedWindowCount(kLen, kStride));

    // (a) every recorded frame index ∈ [start, start+len).
    for (std::size_t fi : frame_indices) {
        EXPECT_GE(fi, kStart) << "dispatched frame below window start";
        EXPECT_LT(fi, kStart + kLen) << "dispatched frame at/past window end";
    }

    // (c) conf0 / RMSD reference == frame `start`.
    ASSERT_FALSE(frame_indices.empty());
    EXPECT_EQ(frame_indices.front(), kStart)
        << "first dispatched frame must be the window start (the reference)";

    const auto& tr = tp.Result<nmr::RmsdTrackingTrajectoryResult>();
    ASSERT_GE(tr.NumFrames(), 1u);
    // The window start is its own Kabsch reference → RMSD exactly 0.
    EXPECT_DOUBLE_EQ(tr.RmsdAtSampleIndex(0), 0.0);

    // The strided window frames are exactly {100,102,104,106,108}.
    const std::vector<std::size_t> expected = {100, 102, 104, 106, 108};
    ASSERT_EQ(frame_indices.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        EXPECT_EQ(frame_indices[i], expected[i]) << "frame slot " << i;
    }
}


// ── No-window run is unchanged (full trajectory) ─────────────────────

TEST(TrajectoryWindow, NoWindowIsFullTrajectory) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    constexpr std::size_t kStride = 2;

    // window_len == 0 → no window. Same stride as the windowed run so the
    // counts are directly comparable.
    auto config = BuildConfig(kStride, /*window_start=*/0, /*window_len=*/0);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& frame_indices = traj.FrameIndices();
    ASSERT_FALSE(frame_indices.empty());

    // (d) Full-trajectory behaviour: starts at frame 0, dispatches far more
    // than a 10-frame window's worth of strided frames, and reaches frames
    // beyond the window's upper bound (110). The reference is frame 0.
    EXPECT_EQ(frame_indices.front(), 0u)
        << "no-window run must start at trajectory frame 0";
    EXPECT_GT(traj.FrameCount(), StridedWindowCount(10, kStride))
        << "no-window run must dispatch more frames than a 10-frame window";
    EXPECT_GT(frame_indices.back(), 110u)
        << "no-window run must extend past the windowed run's upper bound";

    const auto& tr = tp.Result<nmr::RmsdTrackingTrajectoryResult>();
    EXPECT_DOUBLE_EQ(tr.RmsdAtSampleIndex(0), 0.0);
}


// ── window_start past EOF fails loud ─────────────────────────────────

TEST(TrajectoryWindow, WindowStartBeyondTrajectoryFailsLoud) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    // 1P9J is 1501 frames; a window starting past the end cannot seed.
    auto config = BuildConfig(/*stride=*/1, /*window_start=*/10'000'000,
                              /*window_len=*/5);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    // Loud, not silent: Phase 2 returns kFrameReadFailed when window_start
    // runs off the end of the trajectory.
    EXPECT_EQ(traj.Run(tp, config, session), nmr::kFrameReadFailed);
}


// ── window END past EOF fails loud (graceful hard error, not a short window) ──

TEST(TrajectoryWindow, WindowEndBeyondTrajectoryFailsLoud) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    // 1P9J is 1501 frames. A window that STARTS in-bounds but whose end runs
    // past EOF cannot be fully honoured: captured (1495..1500 = 6) < requested
    // (20). The run must fail loud rather than silently emit the short prefix
    // it managed to capture.
    auto config = BuildConfig(/*stride=*/1, /*window_start=*/1495,
                              /*window_len=*/20);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    EXPECT_EQ(traj.Run(tp, config, session), nmr::kFrameReadFailed);
}
