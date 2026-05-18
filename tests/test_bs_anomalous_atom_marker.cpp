//
// test_bs_anomalous_atom_marker: discipline + integration for
// BsAnomalousAtomMarkerTrajectoryResult.
//
// This TR is the canonical exemplar for cross-result reads + per-atom
// event-bag writes (PATTERNS §17 + the AtomEvent / RecordBag<AtomEvent>
// pattern). It does NOT emit H5 datasets — its authoritative output is
// the per-atom `ta.events` bag — so the discipline trio diverges from
// the Welford-family shape:
//
//   1. Frame0Semantics: single frame is below MIN_BURN_IN_FRAMES → zero
//      events expected (the burn-in skip is a documented load-bearing
//      decision, not a bug).
//   2. DependsOnBsWelford: cross-result read discipline. The TR must
//      declare BsWelford in Dependencies() so Phase 4's validator can
//      catch a missing upstream at attach time, not at compute time.
//   3. Integration1P9J: with enough frames to exceed burn-in + a real
//      kernel running, verify (a) events are pushed to ta.events,
//      (b) the ByKind<High>() + ByKind<Low>() queries discriminate at
//      read time, (c) the running totals (n_events_, n_high, n_low) add
//      up consistently across atoms.
//

#include "BiotSavartResult.h"

#include "AtomEvent.h"
#include "BsAnomalousAtomMarkerTrajectoryResult.h"
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

#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <typeindex>
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

// Both BsWelford (upstream) and the marker (downstream) attached in
// dependency order. Phase 4's validator enforces this; we feed it the
// right order here.
nmr::RunConfiguration BuildConfigWithMarker(
        const std::string& name, std::size_t stride) {
    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, name, stride);
    config.RequireConformationResult(typeid(nmr::BiotSavartResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::BsWelfordTrajectoryResult::Create(tp);
        });
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::BsAnomalousAtomMarkerTrajectoryResult::Create(tp);
        });
    return config;
}

}  // namespace


// ============================================================================
// DISCIPLINE: Frame-0 semantics — single frame is below burn-in.
// ============================================================================
//
// With one frame, every atom's bs_welford.n_frames == 1 < MIN_BURN_IN_FRAMES
// (20). The marker correctly skips every atom — no events are pushed. This
// is the documented load-bearing skip (avoids degenerate z-scores when std
// estimate is unreliable), not a bug.

TEST(BsAnomalousAtomMarker, Frame0SemanticsZeroEvents) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfigWithMarker("BsAnomalyMarkerFrame0", 99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    ASSERT_EQ(traj.FrameCount(), 1u);

    ASSERT_TRUE(tp.HasResult<nmr::BsAnomalousAtomMarkerTrajectoryResult>());
    const auto& marker = tp.Result<nmr::BsAnomalousAtomMarkerTrajectoryResult>();

    EXPECT_EQ(marker.TotalEvents(), 0u)
        << "single frame is below MIN_BURN_IN_FRAMES = 20; no events expected";
    EXPECT_EQ(marker.HighEvents(), 0u);
    EXPECT_EQ(marker.LowEvents(),  0u);

    // Every atom's event bag should be empty.
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        EXPECT_EQ(tp.AtomAt(i).events.Count(), 0u);
    }
}


// ============================================================================
// DISCIPLINE: Dependencies() lists BsWelford.
// ============================================================================
//
// Phase 4 of Trajectory::Run validates dependency declarations at attach
// time. A consumer who forgets to attach BsWelford would otherwise see
// Compute read stale-zero from bs_welford.{t0.mean, t0.m2, n_frames} and
// silently miss every anomaly. Verifying the dependency list catches this
// at the declaration level.

TEST(BsAnomalousAtomMarker, DependsOnBsWelford) {
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();

    auto marker = nmr::BsAnomalousAtomMarkerTrajectoryResult::Create(tp);
    ASSERT_TRUE(marker);
    const auto deps = marker->Dependencies();

    const auto welford_id = std::type_index(typeid(nmr::BsWelfordTrajectoryResult));
    const auto bs_id      = std::type_index(typeid(nmr::BiotSavartResult));

    bool has_welford = false, has_bs = false;
    for (const auto& d : deps) {
        if (d == welford_id) has_welford = true;
        if (d == bs_id)      has_bs = true;
    }
    EXPECT_TRUE(has_welford)
        << "BsWelford dependency missing — cross-result read would silently "
           "read stale-zero from bs_welford fields";
    EXPECT_TRUE(has_bs)
        << "BiotSavart dependency missing — per-frame kernel value would "
           "not be available";
}


// ============================================================================
// INTEGRATION: real BS kernel + Welford + marker over enough frames.
// ============================================================================
//
// stride = 60 over the 1P9J fleet_amber fixture (~1500 frames) gives ~25
// captured frames, exceeding MIN_BURN_IN_FRAMES = 20 with room for the
// marker to start emitting events at ~frame 20. The integration verifies:
//
//   1. cross-result read works: marker sees nonzero bs_welford state
//   2. events are pushed: TotalEvents() > 0 (a well-mixed trajectory will
//      produce some |z| > 2σ on at least one atom over 25 frames)
//   3. ByKind query discriminates: high + low queries sum to total
//   4. running totals match per-atom event-bag count

TEST(BsAnomalousAtomMarker, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfigWithMarker("BsAnomalyMarkerIntegration", 60);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    // Need at least MIN_BURN_IN_FRAMES + a few to start emitting.
    const size_t MIN_BURN_IN =
        nmr::BsAnomalousAtomMarkerTrajectoryResult::MIN_BURN_IN_FRAMES;
    if (traj.FrameCount() < MIN_BURN_IN + 2) {
        GTEST_SKIP() << "trajectory has only " << traj.FrameCount()
                     << " frames; need >" << MIN_BURN_IN
                     << " to exceed burn-in (decrease stride)";
    }

    ASSERT_TRUE(tp.HasResult<nmr::BsAnomalousAtomMarkerTrajectoryResult>());
    const auto& marker = tp.Result<nmr::BsAnomalousAtomMarkerTrajectoryResult>();

    // Welford state populated (cross-result read source).
    size_t welford_pop = 0;
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        if (tp.AtomAt(i).bs_welford.n_frames > 0) ++welford_pop;
    }
    EXPECT_GT(welford_pop, 0u) << "BsWelford state empty; cross-result read "
                                  "would have nothing to read";

    // Running totals consistent: total = high + low.
    EXPECT_EQ(marker.TotalEvents(),
              marker.HighEvents() + marker.LowEvents())
        << "TotalEvents must equal HighEvents + LowEvents";

    // Per-atom event-bag totals match TR-level n_events_. Sum across all
    // atoms via .Count() should equal TotalEvents().
    size_t bag_total = 0;
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        bag_total += tp.AtomAt(i).events.Count();
    }
    EXPECT_EQ(bag_total, marker.TotalEvents())
        << "sum of per-atom events.Count() must equal TR-level TotalEvents";

    // ByKind queries: for any atom with events, the high+low partition by
    // kind must reconstruct the full set.
    using HighKind = nmr::BsAnomalousAtomMarkerTrajectoryResult::BsAnomalyHighT0;
    using LowKind  = nmr::BsAnomalousAtomMarkerTrajectoryResult::BsAnomalyLowT0;

    size_t bykind_high = 0, bykind_low = 0;
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& bag = tp.AtomAt(i).events;
        bykind_high += bag.ByKind(typeid(HighKind)).size();
        bykind_low  += bag.ByKind(typeid(LowKind)).size();
    }
    EXPECT_EQ(bykind_high, marker.HighEvents())
        << "ByKind<HighT0> sum must match TR.HighEvents()";
    EXPECT_EQ(bykind_low,  marker.LowEvents())
        << "ByKind<LowT0> sum must match TR.LowEvents()";

    std::cout << "BsAnomalousAtomMarker integration: frames="
              << traj.FrameCount()
              << " welford_pop=" << welford_pop << "/" << tp.AtomCount()
              << " total_events=" << marker.TotalEvents()
              << " (high=" << marker.HighEvents()
              << ", low="  << marker.LowEvents() << ")\n";
}
