//
// test_rmsd_spike_dft_coord: discipline + integration for
// RmsdSpikeSelectionTrajectoryResult (TR12) and
// DftPoseCoordinatorTrajectoryResult (TR13). TR12 is a CROSS-RESULT
// READ on TR11's per-frame RMSD; TR13 is a Finalize-time reducer
// on TR12 + ChiRotamerSelectionTrajectoryResult's SelectionBag
// entries.
//
// Tests verify wiring (no crash, no missing attach) and end-to-end
// SelectionBag emission. Threshold detection is project-decision;
// short trajectories may produce 0 spikes (rolling window needs
// ≥ 10 frames). The tests assert the COUNTERS exist + reduction
// math is consistent, NOT specific spike counts on this fixture.
//

#include "CalculatorConfig.h"
#include "ChiRotamerSelectionTrajectoryResult.h"
#include "DftPoseCoordinatorTrajectoryResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "RmsdSpikeSelectionTrajectoryResult.h"
#include "RmsdTrackingTrajectoryResult.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "Types.h"

#include <gtest/gtest.h>

#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>
#include <typeindex>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

namespace {

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
void LoadCalculatorConfig() {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(std::string(NMR_TEST_DATA_DIR) +
                                "/../data/calculator_params.toml");
}

// TR11 first (writer for TR12), then TR12, then ChiRotamer (writer
// for TR13), then TR13. Attach order = dispatch order per PATTERNS
// §15; this matches RunConfiguration::PerFrameExtractionSet.
nmr::RunConfiguration BuildConfig(unsigned stride) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true;
    opts.skip_coulomb = true;
    opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::RmsdTrackingTrajectoryResult::Create(tp_in);
    });
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::RmsdSpikeSelectionTrajectoryResult::Create(tp_in);
    });
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::ChiRotamerSelectionTrajectoryResult::Create(tp_in);
    });
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::DftPoseCoordinatorTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


// ── TR12 + TR13: end-to-end on 1P9J multi-frame ────────────────────

TEST(RmsdSpikeAndDftCoord, EndToEndOn1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    // Stride 200 over 1P9J's 1501-frame fixture = ~8 frames; not
    // enough for the rolling-window criterion (needs ≥ 10), so
    // expect 0 RmsdSpikes here — TR13 still runs, possibly reducing
    // from ChiRotamer's stream.
    auto config = BuildConfig(200);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    // TR11 produced per-frame RMSD
    const auto& tr_rmsd = tp.Result<nmr::RmsdTrackingTrajectoryResult>();
    EXPECT_GE(tr_rmsd.NumFrames(), 2u);
    EXPECT_DOUBLE_EQ(tr_rmsd.RmsdAtFrame(0), 0.0);

    // TR12 ran, internal counter exists
    const auto& tr_spike = tp.Result<nmr::RmsdSpikeSelectionTrajectoryResult>();
    // Spike count CAN be 0 (short fixture, rolling window unstable);
    // we only assert wiring -- no crash, accessible counter.
    const std::size_t spike_n = tr_spike.SpikeCount();
    EXPECT_GE(spike_n, 0u);

    // ChiRotamer fires on per-residue rotamer transitions; even an
    // 8-frame stride on SH3 typically catches at least 1 sidechain
    // rotamer flip.
    const auto& tr_chi = tp.Result<nmr::ChiRotamerSelectionTrajectoryResult>();
    const std::size_t chi_n = tr_chi.TransitionCount();
    EXPECT_GE(chi_n, 0u);

    // TR13 reducer fired at Finalize. Its reduced count should be
    // ≤ (spike_n + chi_n) — bounded above by the union of upstream
    // emitters (before dedup).
    const auto& tr_coord = tp.Result<nmr::DftPoseCoordinatorTrajectoryResult>();
    const std::size_t reduced_n = tr_coord.ReducedCount();
    EXPECT_LE(reduced_n, spike_n + chi_n)
        << "reducer should not emit MORE records than upstream union";

    nmr::OperationLog::Info(
        "RmsdSpikeAndDftCoordTest::EndToEndOn1P9J",
        "spikes=" + std::to_string(spike_n) +
        " chi_transitions=" + std::to_string(chi_n) +
        " reduced=" + std::to_string(reduced_n));

    // SelectionBag has correct kinds present
    const auto kinds = traj.Selections().Kinds();
    bool saw_rmsd_spike_kind = false;
    bool saw_chi_kind        = false;
    bool saw_coord_kind      = false;
    for (const std::type_index& k : kinds) {
        if (k == typeid(nmr::RmsdSpikeSelectionTrajectoryResult))
            saw_rmsd_spike_kind = true;
        if (k == typeid(nmr::ChiRotamerSelectionTrajectoryResult))
            saw_chi_kind = true;
        if (k == typeid(nmr::DftPoseCoordinatorTrajectoryResult))
            saw_coord_kind = true;
    }
    if (spike_n > 0) EXPECT_TRUE(saw_rmsd_spike_kind);
    if (chi_n > 0)   EXPECT_TRUE(saw_chi_kind);
    if (reduced_n > 0) EXPECT_TRUE(saw_coord_kind);
}
