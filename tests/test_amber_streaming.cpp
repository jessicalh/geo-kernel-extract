//
// test_amber_streaming: AMBER trajectory streaming via TRR (libgromacs
// gmx_trr_* path) and the trajectory-scope object model.
//
// Replaces the retired test_gromacs_streaming.cpp (CHARMM/XTC,
// fleet_test_fullsys/1ZR7_6721) which became structurally stale after
// the 2026-05-02 CHARMM-retired / TRR-replaces-XTC decision. This file
// covers the same integration-pattern surface against the AMBER+TRR
// fixtures (1P9J_5801, 1Z9B_6577): TrajectoryProtein::BuildFromTrajectory,
// GromacsFrameHandler reading TRR frames, Trajectory::Run driving the
// per-frame loop, and the BondLengthStats discipline trio.
//
// Per-TR full coverage (BsWelford, PositionsTimeSeries, BsAnomalyMarker,
// BsT0Autocorrelation, ChiRotamerSelection) lives in their own dedicated
// test files; this file is the integration-pattern anchor for the
// trajectory-scope architecture and the AMBER+TRR path. BondLengthStats
// is the chosen TR for discipline tests because Dependencies() == {} —
// minimal noise from the per-frame stack.
//

// AIMNet2Result.h MUST come before GROMACS headers — GROMACS defines
// DIM as a macro which poisons PyTorch template parameters.
#include "AIMNet2Result.h"

// Nanoflann-using headers (SpatialIndexResult, via BiotSavartResult +
// McConnell) MUST come before GROMACS headers too — gromacs vectypes.h
// #define DIM 3 collides with nanoflann's DIM template parameter.
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "BiotSavartResult.h"

#include "TrajectoryProtein.h"
#include "Trajectory.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "BondLengthStatsTrajectoryResult.h"
#include "CalculatorConfig.h"
#include "OperationRunner.h"
#include "OperationLog.h"
#include "TestEnvironment.h"

#include <gtest/gtest.h>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

namespace {

// fixture.xtc_path is a misnomer post-TRR-replaces-XTC; the AMBER
// fixtures carry production.{tpr,trr,edr}. Until TestEnvironment grows
// a typed trr_path field, derive .trr from .tpr by extension swap (same
// prefix, both files coexist in the fixture).
std::string TrrPathFor(const std::string& tpr_path) {
    return fs::path(tpr_path).replace_extension(".trr").string();
}

// Trajectory directory = parent of production.tpr (TrajectoryProtein::
// BuildFromTrajectory consumes a directory, then composes
// <dir>/production.tpr and <dir>/../topol.top per the documented layout).
std::string ProductionDirFor(const std::string& tpr_path) {
    return fs::path(tpr_path).parent_path().string();
}

bool FixtureAvailable(const nmr::test::AmberTrajectoryFixture& fix) {
    return !fix.tpr_path.empty() &&
           fs::exists(fix.tpr_path) &&
           fs::exists(TrrPathFor(fix.tpr_path)) &&
           fs::exists(fix.edr_path);
}

void ConfigureBondLengthStats(nmr::RunConfiguration& config,
                              const std::string& test_name) {
    config.SetName(test_name);
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::BondLengthStatsTrajectoryResult::Create(tp);
        });
}

void LoadCalculatorConfig() {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../data/calculator_params.toml");
}

const std::string kFixtureProtein = "1P9J_5801";

}  // namespace


// ============================================================================
// TrajectoryProtein builds + Trajectory::Run drives the per-frame loop.
// ============================================================================

TEST(AmberStreaming, TrajectoryBuildAndRun) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fleet_amber " << kFixtureProtein
                                              << " fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();

    EXPECT_GT(tp.ProteinRef().AtomCount(), 100u);
    EXPECT_GT(tp.ProteinRef().ResidueCount(), 10u);

    // Minimal per-frame run: GeometryResult + SpatialIndexResult only.
    nmr::RunConfiguration config;
    config.SetName("AmberStreamingBuildAndRun");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));

    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path,
                         fix.edr_path);
    nmr::Session session;

    const nmr::Status s = traj.Run(tp, config, session);
    ASSERT_EQ(s, nmr::kOk);
    EXPECT_TRUE(traj.IsComplete());
    EXPECT_TRUE(tp.IsFinalized());
    EXPECT_GT(traj.FrameCount(), 1u);
    EXPECT_GT(tp.ProteinRef().BondCount(), 0u);
}


// ============================================================================
// BondLengthStats end-to-end. Validates trajectory-scope AV pattern:
// per-frame Compute updates Welford accumulators, Finalize lands the
// stats. Physical plausibility: bond lengths are 0.5–3.0 Å.
// ============================================================================

TEST(AmberStreaming, BondLengthStatsEndToEnd) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    ConfigureBondLengthStats(config, "AmberStreamingBondLengthStatsEndToEnd");

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;

    const nmr::Status s = traj.Run(tp, config, session);
    ASSERT_EQ(s, nmr::kOk);
    EXPECT_TRUE(traj.IsComplete());
    EXPECT_TRUE(tp.IsFinalized());

    ASSERT_GT(tp.ProteinRef().BondCount(), 0u);
    const size_t n_frames = traj.FrameCount();
    ASSERT_GT(n_frames, 1u);

    ASSERT_TRUE(tp.HasResult<nmr::BondLengthStatsTrajectoryResult>());
    const auto& bonds =
        tp.Result<nmr::BondLengthStatsTrajectoryResult>().PerBond();
    ASSERT_EQ(bonds.size(), tp.ProteinRef().BondCount());

    bool any_moved = false;
    for (const auto& pb : bonds) {
        EXPECT_EQ(pb.n_frames, n_frames);
        EXPECT_GT(pb.length_mean, 0.5);
        EXPECT_LT(pb.length_mean, 3.0);
        EXPECT_GE(pb.length_std, 0.0);
        EXPECT_LE(pb.length_min, pb.length_mean);
        EXPECT_GE(pb.length_max, pb.length_mean);
        if (pb.length_std > 1e-5) any_moved = true;
    }
    EXPECT_TRUE(any_moved)
        << "expected some bonds to fluctuate in length across frames";

    std::cout << "BondLengthStats: " << bonds.size()
              << " bonds × " << n_frames << " frames; sample bond 0 mean="
              << bonds[0].length_mean << "Å std=" << bonds[0].length_std
              << "Å\n";
}


// ============================================================================
// BondLengthStats discipline trio. Pattern from spec/INDEX.md:
//   - frame-0 semantics (stride ≥ fixture length so only frame 0 dispatches)
//   - Finalize idempotency
//   - H5 round-trip via a temp file
// ============================================================================

TEST(AmberStreaming, BondLengthStatsFrame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    ConfigureBondLengthStats(config, "AmberStreamingFrame0Semantics");
    config.SetStride(99999);  // > fixture length → only frame 0 dispatches

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    const nmr::Status s = traj.Run(tp, config, session);
    ASSERT_EQ(s, nmr::kOk);
    ASSERT_EQ(traj.FrameCount(), 1u)
        << "stride > fixture length should leave only frame 0 dispatched";

    const auto& bonds =
        tp.Result<nmr::BondLengthStatsTrajectoryResult>().PerBond();
    ASSERT_EQ(bonds.size(), tp.ProteinRef().BondCount());

    // AV-pattern frame-0 semantics: one sample per bond, no deltas
    // (delta needs a prior frame), min == mean == max, std == 0.
    for (const auto& pb : bonds) {
        EXPECT_EQ(pb.n_frames, 1u);
        EXPECT_EQ(pb.delta_n, 0u)
            << "delta tracker needs a prior frame — must be 0 at frame 0";
        EXPECT_DOUBLE_EQ(pb.length_min, pb.length_mean);
        EXPECT_DOUBLE_EQ(pb.length_mean, pb.length_max);
        EXPECT_DOUBLE_EQ(pb.length_std, 0.0)
            << "n=1 has no variance";
        EXPECT_DOUBLE_EQ(pb.delta_std, 0.0);
    }
}


TEST(AmberStreaming, BondLengthStatsFinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    ConfigureBondLengthStats(config, "AmberStreamingFinalizeIdempotency");

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    const nmr::Status s = traj.Run(tp, config, session);
    ASSERT_EQ(s, nmr::kOk);

    auto& blst = tp.Result<nmr::BondLengthStatsTrajectoryResult>();
    const auto snapshot = blst.PerBond();   // copy post-first-Finalize

    // Re-run Finalize; state must be unchanged.
    blst.Finalize(tp, traj);

    const auto& after = blst.PerBond();
    ASSERT_EQ(after.size(), snapshot.size());
    for (size_t i = 0; i < after.size(); ++i) {
        EXPECT_DOUBLE_EQ(after[i].length_mean,  snapshot[i].length_mean)  << "bond " << i;
        EXPECT_DOUBLE_EQ(after[i].length_m2,    snapshot[i].length_m2)    << "bond " << i;
        EXPECT_DOUBLE_EQ(after[i].length_std,   snapshot[i].length_std)   << "bond " << i;
        EXPECT_DOUBLE_EQ(after[i].length_min,   snapshot[i].length_min)   << "bond " << i;
        EXPECT_DOUBLE_EQ(after[i].length_max,   snapshot[i].length_max)   << "bond " << i;
        EXPECT_EQ(after[i].n_frames,            snapshot[i].n_frames)     << "bond " << i;
        EXPECT_DOUBLE_EQ(after[i].delta_mean,   snapshot[i].delta_mean)   << "bond " << i;
        EXPECT_DOUBLE_EQ(after[i].delta_m2,     snapshot[i].delta_m2)     << "bond " << i;
        EXPECT_DOUBLE_EQ(after[i].delta_std,    snapshot[i].delta_std)    << "bond " << i;
        EXPECT_EQ(after[i].delta_n,             snapshot[i].delta_n)      << "bond " << i;
    }
}


TEST(AmberStreaming, BondLengthStatsH5RoundTrip) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    ConfigureBondLengthStats(config, "AmberStreamingH5RoundTrip");

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    const nmr::Status s = traj.Run(tp, config, session);
    ASSERT_EQ(s, nmr::kOk);

    const auto& blst = tp.Result<nmr::BondLengthStatsTrajectoryResult>();
    const auto pre = blst.PerBond();

    const std::string h5_path = (fs::temp_directory_path() /
        ("amber_streaming_blst_h5_roundtrip_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        blst.WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    // Read it back via HighFive and compare a representative sample.
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/bond_length_stats"));

    auto length_mean_ds = reopen.getDataSet(
        "/trajectory/bond_length_stats/length_mean");
    std::vector<double> length_mean_read;
    length_mean_ds.read(length_mean_read);
    ASSERT_EQ(length_mean_read.size(), pre.size());
    for (size_t i = 0; i < pre.size(); ++i) {
        EXPECT_DOUBLE_EQ(length_mean_read[i], pre[i].length_mean)
            << "bond " << i;
    }

    fs::remove(h5_path);
}
