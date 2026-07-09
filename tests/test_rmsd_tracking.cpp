//
// test_rmsd_tracking: discipline + integration for
// RmsdTrackingTrajectoryResult (TR11). Scalar AV (T,) double — first
// scalar-AV TR shape in the codebase. Backbone heavy atoms (N/CA/C/O)
// vs the first dispatched frame, Kabsch-aligned.
//

#include "CalculatorConfig.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RmsdTrackingTrajectoryResult.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "Types.h"

#include <gtest/gtest.h>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>
#include <unistd.h>

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
nmr::RunConfiguration BuildConfig(unsigned stride) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    // TR11 reads positions only.
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
    return config;
}

}  // namespace


// ── Frame 0 semantics ───────────────────────────────────────────────

TEST(RmsdTracking, Frame0Semantics) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_EQ(traj.FrameCount(), 1u);

    const auto& tr = tp.Result<nmr::RmsdTrackingTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);
    EXPECT_GT(tr.NumAlignmentAtoms(), 0u);
    // Frame 0 is its own reference -- RMSD must be exactly 0.
    EXPECT_DOUBLE_EQ(tr.RmsdAtSampleIndex(0), 0.0);
}


// ── Finalize idempotency ────────────────────────────────────────────

TEST(RmsdTracking, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::RmsdTrackingTrajectoryResult>();
    const std::size_t T = tr.NumFrames();
    tr.Finalize(tp, traj);  // double-Finalize must be safe
    EXPECT_EQ(tr.NumFrames(), T);
    EXPECT_DOUBLE_EQ(tr.RmsdAtSampleIndex(0), 0.0);
}


// ── H5 round-trip + dataset shapes ─────────────────────────────────

TEST(RmsdTracking, H5RoundTrip) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::RmsdTrackingTrajectoryResult>();
    const std::size_t M = tr.NumAlignmentAtoms();
    const std::size_t T = tr.NumFrames();

    const std::string h5_path = (fs::temp_directory_path() /
        ("rmsd_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/rmsd_tracking"));
    auto grp = reopen.getGroup("/trajectory/rmsd_tracking");

    // Convention attrs
    std::string alignment, sel, units, policy;
    grp.getAttribute("alignment_method").read(alignment);
    grp.getAttribute("atom_selection").read(sel);
    grp.getAttribute("units").read(units);
    grp.getAttribute("source_attached_policy").read(policy);
    EXPECT_EQ(alignment, "kabsch_svd");
    EXPECT_EQ(sel, "backbone_heavy_atoms_NCACO");
    EXPECT_EQ(units, "Angstrom");
    EXPECT_NE(policy.find("always_attached"), std::string::npos);

    // Dataset shapes
    {
        ASSERT_TRUE(grp.exist("rmsd"));
        const auto dims = grp.getDataSet("rmsd").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 1u);
        EXPECT_EQ(dims[0], T);
    }
    {
        ASSERT_TRUE(grp.exist("atom_indices"));
        const auto dims = grp.getDataSet("atom_indices").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 1u);
        EXPECT_EQ(dims[0], M);
    }
    EXPECT_TRUE(grp.exist("frame_indices"));
    EXPECT_TRUE(grp.exist("frame_times"));
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));

    // Content read-back
    std::vector<double> rmsd;
    grp.getDataSet("rmsd").read(rmsd);
    ASSERT_EQ(rmsd.size(), T);
    EXPECT_DOUBLE_EQ(rmsd[0], 0.0);
    for (std::size_t t = 0; t < T; ++t) {
        EXPECT_TRUE(std::isfinite(rmsd[t])) << "non-finite RMSD at frame " << t;
        EXPECT_GE(rmsd[t], 0.0) << "negative RMSD at frame " << t;
    }
}


// ── Insufficient alignment set sentinel ─────────────────────────────

TEST(RmsdTracking, InsufficientAlignmentAtomsSentinelAndH5) {
    nmr::Protein protein;
    nmr::ProteinConformation conf0(
        &protein,
        {nmr::Vec3(0.0, 0.0, 0.0), nmr::Vec3(1.0, 0.0, 0.0)},
        "rmsd-insufficient-0");
    nmr::ProteinConformation conf1(
        &protein,
        {nmr::Vec3(0.0, 0.0, 0.0), nmr::Vec3(2.0, 0.0, 0.0)},
        "rmsd-insufficient-1");
    nmr::TrajectoryProtein tp;
    nmr::Trajectory traj("", "", "");

    auto tr = nmr::RmsdTrackingTrajectoryResult::CreateForTesting({0, 1});
    tr->Compute(conf0, tp, traj, 7, 1.25);
    tr->Compute(conf1, tp, traj, 11, 2.50);

    ASSERT_EQ(tr->NumFrames(), 2u);
    EXPECT_TRUE(tr->InsufficientAlignmentAtoms());
    EXPECT_TRUE(std::isnan(tr->RmsdAtSampleIndex(0)));
    EXPECT_TRUE(std::isnan(tr->RmsdAtSampleIndex(1)));
    EXPECT_TRUE(std::isnan(tr->LatestRmsd()));

    const std::string h5_path = (fs::temp_directory_path() /
        ("rmsd_insufficient_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/rmsd_tracking");

    bool insufficient = false;
    std::size_t min_alignment_atoms = 0;
    std::size_t reference_frame = 0;
    std::string convention;
    grp.getAttribute("insufficient_alignment_atoms").read(insufficient);
    grp.getAttribute("min_alignment_atoms").read(min_alignment_atoms);
    grp.getAttribute("reference_frame_trr_index").read(reference_frame);
    grp.getAttribute("rmsd_frame_0_convention").read(convention);
    EXPECT_TRUE(insufficient);
    EXPECT_EQ(min_alignment_atoms, 3u);
    EXPECT_EQ(reference_frame, 7u);
    EXPECT_NE(convention.find("NaN when n_atoms < 3"), std::string::npos);

    std::vector<std::uint8_t> mask;
    grp.getDataSet("insufficient_alignment_atoms_mask").read(mask);
    ASSERT_EQ(mask.size(), 1u);
    EXPECT_EQ(mask[0], 1u);

    std::vector<double> rmsd;
    grp.getDataSet("rmsd").read(rmsd);
    ASSERT_EQ(rmsd.size(), 2u);
    EXPECT_TRUE(std::isnan(rmsd[0]));
    EXPECT_TRUE(std::isnan(rmsd[1]));

    fs::remove(h5_path);
}


TEST(RmsdTracking, SufficientAlignmentAtomsKeepFrame0ZeroAndMaskClear) {
    nmr::Protein protein;
    nmr::ProteinConformation conf0(
        &protein,
        {nmr::Vec3(0.0, 0.0, 0.0),
         nmr::Vec3(1.0, 0.0, 0.0),
         nmr::Vec3(0.0, 1.0, 0.0)},
        "rmsd-sufficient-0");
    nmr::TrajectoryProtein tp;
    nmr::Trajectory traj("", "", "");

    auto tr = nmr::RmsdTrackingTrajectoryResult::CreateForTesting({0, 1, 2});
    tr->Compute(conf0, tp, traj, 3, 0.75);

    ASSERT_EQ(tr->NumFrames(), 1u);
    EXPECT_FALSE(tr->InsufficientAlignmentAtoms());
    EXPECT_DOUBLE_EQ(tr->RmsdAtSampleIndex(0), 0.0);

    const std::string h5_path = (fs::temp_directory_path() /
        ("rmsd_sufficient_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/rmsd_tracking");

    bool insufficient = true;
    std::size_t reference_frame = 0;
    grp.getAttribute("insufficient_alignment_atoms").read(insufficient);
    grp.getAttribute("reference_frame_trr_index").read(reference_frame);
    EXPECT_FALSE(insufficient);
    EXPECT_EQ(reference_frame, 3u);

    std::vector<std::uint8_t> mask;
    grp.getDataSet("insufficient_alignment_atoms_mask").read(mask);
    ASSERT_EQ(mask.size(), 1u);
    EXPECT_EQ(mask[0], 0u);

    std::vector<double> rmsd;
    grp.getDataSet("rmsd").read(rmsd);
    ASSERT_EQ(rmsd.size(), 1u);
    EXPECT_DOUBLE_EQ(rmsd[0], 0.0);

    fs::remove(h5_path);
}


// ── Multi-frame physical-range sanity ──────────────────────────────

TEST(RmsdTracking, RmsdRange1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    // Stride 300 over 1P9J's 1501-frame fixture = ~5 frames spanning
    // the trajectory.
    auto config = BuildConfig(300);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::RmsdTrackingTrajectoryResult>();
    const auto& rmsd = tr.Rmsd();
    ASSERT_GE(rmsd.size(), 2u);

    EXPECT_DOUBLE_EQ(rmsd[0], 0.0);
    // SH3 in 15 ns MD: backbone RMSD typically grows monotonically(-ish)
    // over the first few ns, plateauing at 1-3 Å. Loose bounds:
    //   - non-zero by frame 1 (geometry moves)
    //   - sub-10 Å through the whole run (otherwise the protein has
    //     unfolded; not the case for 1P9J at 300K)
    double max_rmsd = 0.0;
    for (std::size_t t = 1; t < rmsd.size(); ++t) {
        EXPECT_GT(rmsd[t], 0.0) << "RMSD zero past frame 0 — Kabsch identity?";
        EXPECT_LT(rmsd[t], 10.0) << "RMSD > 10A at frame " << t
                                  << " — unfolding? Inspect.";
        if (rmsd[t] > max_rmsd) max_rmsd = rmsd[t];
    }
    nmr::OperationLog::Info(
        "RmsdTrackingTest::RmsdRange1P9J",
        "n_frames=" + std::to_string(rmsd.size()) +
        " max_rmsd=" + std::to_string(max_rmsd) + " A");
}
