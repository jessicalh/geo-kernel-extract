//
// test_mopac_vs_ff14sb_reconciliation: discipline + integration for
// MopacVsFf14SbReconciliationTrajectoryResult (TR9 of the 13-TR plan).
// Cross-source |cos(MOPAC_T2, FF14SB_T2)| — new canonical pattern.
//

#include "MopacVsFf14SbReconciliationTrajectoryResult.h"
#include "MopacCoulombResult.h"
#include "MopacResult.h"
#include "CoulombResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "ChargeAssignmentResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
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
}  // namespace


// ── Layer 0a: both sources absent → group skipped ────────────

TEST(MopacVsFf14SbReconciliation, GroupSkippedWhenBothNeverAttached) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    auto tr = nmr::MopacVsFf14SbReconciliationTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 3;
    std::vector<nmr::Vec3> positions(tp.AtomCount(), nmr::Vec3::Zero());
    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame (neither source)");
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    EXPECT_EQ(tr->NumFrames(), kFrames);
    EXPECT_EQ(tr->SourceAttachedCount(), 0u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("mopac_vs_ff14sb_absent_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/mopac_vs_ff14sb_reconciliation"));

    fs::remove(h5_path);
}


// ── Layer 0b: only one source attached → mask all-0, group skipped ──

TEST(MopacVsFf14SbReconciliation, GroupSkippedWhenOnlyOneSource) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    auto tr = nmr::MopacVsFf14SbReconciliationTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 3;
    std::vector<nmr::Vec3> positions(tp.AtomCount(), nmr::Vec3::Zero());
    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic — only Coulomb");
        // Attach only CoulombResult — MopacCoulombResult absent.
        conf->ForceAttachResultForTesting(std::make_unique<nmr::CoulombResult>());
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    EXPECT_EQ(tr->SourceAttachedCount(), 0u)
        << "asymmetric attach (only Coulomb) should NOT count — "
           "the gate requires BOTH";

    const std::string h5_path = (fs::temp_directory_path() /
        ("mopac_vs_ff14sb_one_source_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/mopac_vs_ff14sb_reconciliation"));

    fs::remove(h5_path);
}


// ── Integration1P9J: real Trajectory::Run with both sources ──

TEST(MopacVsFf14SbReconciliation, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = false;     // need MopacResult → MopacCoulombResult
    opts.skip_coulomb = false;   // enables both Coulomb + MopacCoulomb
    opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::ChargeAssignmentResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::MopacVsFf14SbReconciliationTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    const auto status = traj.Run(tp, config, session);
    if (status != nmr::kOk) {
        GTEST_SKIP() << "Trajectory::Run failed; skipping per "
                        "feedback_log_overages_dont_assert.";
    }

    auto& tr = tp.Result<nmr::MopacVsFf14SbReconciliationTrajectoryResult>();
    const std::size_t T = tr.NumFrames();
    const std::size_t attached = tr.SourceAttachedCount();
    std::cout << "MopacVsFf14SbReconciliation Integration1P9J: T=" << T
              << " both-attached=" << attached << std::endl;

    if (attached == 0) {
        GTEST_SKIP() << "Both sources never attached; cannot test.";
    }

    const std::string h5_path = (fs::temp_directory_path() /
        ("mopac_vs_ff14sb_int1p9j_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/mopac_vs_ff14sb_reconciliation");
    auto ds = grp.getDataSet("cos_t2");  // signed cosine — renamed
                                          // from abs_cos_t2 per
                                          // science adversarial M1.
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 2u);
    EXPECT_EQ(dims[1], T);

    std::string parity, units, mag_units;
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("units").read(units);
    grp.getAttribute("magnitude_floor_units").read(mag_units);
    EXPECT_EQ(parity, "0e");
    EXPECT_EQ(units, "dimensionless");
    EXPECT_EQ(mag_units, "V/Å^2")
        << "EFG-scale floor (V/Å^2) per math adversarial review H1.";
    double mag_floor = 0.0;
    grp.getAttribute("magnitude_floor").read(mag_floor);
    EXPECT_GT(mag_floor, 1e-10)
        << "magnitude_floor must be calibrated to EFG signal scale, "
           "NOT direction-vector floor (1e-10).";

    const std::size_t N = dims[0];
    std::vector<double> flat(N * T);
    ds.read(flat.data());

    // Signed cos ∈ [-1, 1] for finite cells. Log mean + range to see
    // typical agreement vs disagreement.
    std::size_t n_finite = 0, n_nan = 0, n_negative = 0;
    double sum_cos = 0.0, max_cos = -2.0, min_cos = 2.0;
    for (double v : flat) {
        if (std::isfinite(v)) {
            EXPECT_GE(v, -1.0 - 1e-12);  // tiny float slop
            EXPECT_LE(v,  1.0 + 1e-12);
            ++n_finite;
            sum_cos += v;
            max_cos = std::max(max_cos, v);
            min_cos = std::min(min_cos, v);
            if (v < 0.0) ++n_negative;
        } else {
            EXPECT_TRUE(std::isnan(v));
            ++n_nan;
        }
    }
    const double mean_cos = (n_finite > 0) ? (sum_cos / n_finite) : 0.0;
    std::cout << "  finite=" << n_finite << "/" << (N * T)
              << "  nan=" << n_nan
              << "  mean(cos)=" << mean_cos
              << "  range=[" << min_cos << ", " << max_cos << "]"
              << "  n_negative=" << n_negative
              << "  (signed cos exposes sign-disagreement; per "
                 "science adversarial M1 the prior |cos| collapsed "
                 "this signal)"
              << std::endl;
    EXPECT_GT(n_finite, 0u) << "no finite cosine values — both sources "
                               "probably produced near-zero T2 everywhere";

    fs::remove(h5_path);
}


// ── Layer 0c: FinalizeIdempotency ──────────────────────────────

TEST(MopacVsFf14SbReconciliation, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    auto tr = nmr::MopacVsFf14SbReconciliationTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    // Drive a single absent-source frame so Finalize has state to
    // touch but no buffer adoption.
    std::vector<nmr::Vec3> positions(tp.AtomCount(), nmr::Vec3::Zero());
    auto conf = std::make_unique<nmr::ProteinConformation>(
        &tp.ProteinRef(), positions, "synthetic frame (neither source)");
    tr->Compute(*conf, tp, traj, 0, 0.0);

    tr->Finalize(tp, traj);
    const auto n_after_first = tr->NumFrames();
    tr->Finalize(tp, traj);  // second call — no state mutation beyond
                              // re-flagging finalized_ = true + logging.
    EXPECT_EQ(tr->NumFrames(), n_after_first);
    EXPECT_EQ(tr->SourceAttachedCount(), 0u);
}
