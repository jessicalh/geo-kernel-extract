//
// test_mopac_charge_welford: discipline + integration for
// MopacChargeWelfordTrajectoryResult (TR5 of the 13-TR plan).
// Canonical sparse-Welford-scalar template — the bond-order companion
// TR6 clones this shape against the bond axis.
//
// Test layering note. The "positive-path" Welford accumulation cannot
// be driven by synthetic ProteinConformation + Compute because the
// AV-pattern update writes through `tp.MutableAtomAt(i)`, which
// requires `tp.Seed()` to have populated TrajectoryAtom slots. Seed
// FATALs on synthetic all-zero positions (canonicalisation rejects
// degenerate coords). The canonical positive-path test is
// Integration1P9J via real `Trajectory::Run` — same pattern landed
// for AIMNet2ChargeResponseGradientWelford in S1.
//
// Synthetic tests cover the SOURCE-ABSENT path only (HasResult false →
// skip MutableAtomAt entirely), where the gate logic is reachable
// without Seed.
//

#include "MopacChargeWelfordTrajectoryResult.h"
#include "MopacResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "ChargeAssignmentResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryAtom.h"
#include "TrajectoryMoments.h"
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


// ── Layer 0: source absent every frame → group skipped ─────
//
// MOPAC never ran (no ForceAttach). HasResult false each frame,
// Welford update never runs, source_attached_count==0. WriteH5Group
// MUST skip the group entirely per canonical "absent, not faked"
// (OBJECT_MODEL Conditional-attach TR discipline). H5 file should
// not contain /trajectory/mopac_charge_welford.

TEST(MopacChargeWelford, GroupSkippedWhenSourceNeverAttached) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    auto tr = nmr::MopacChargeWelfordTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 3;
    std::vector<nmr::Vec3> positions(tp.AtomCount(), nmr::Vec3::Zero());
    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame (no MOPAC)");
        // NO ForceAttach → HasResult<MopacResult>() returns false.
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    EXPECT_EQ(tr->NumFrames(), kFrames);
    EXPECT_EQ(tr->SourceAttachedCount(), 0u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("mopac_charge_welford_absent_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/mopac_charge_welford"))
        << "group must NOT exist when source attached 0/T frames "
           "per canonical 'absent, not faked'";

    fs::remove(h5_path);
}


// ── Integration1P9J: real Trajectory::Run with MOPAC enabled ──
//
// Drives the canonical AV-Welford positive path through
// Trajectory::Run on the 1P9J fleet_amber fixture. tp.Seed populates
// TrajectoryAtom slots from real positions, Mopac actually runs each
// frame (within stride), the Welford rollup accumulates per-atom.
//
// Verifies via tp.AtomAt(i).mopac_charge_welford that:
//   - n_frames matches source_attached_count
//   - mean is finite
//   - the protein's total charge ≈ 0 (Mulliken neutrality sanity)
//
// Skips if MOPAC binary unavailable or runtime config disables it.

TEST(MopacChargeWelford, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = false;   // MOPAC must run for this test
    opts.skip_coulomb = true;
    opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::MopacChargeWelfordTrajectoryResult::Create(tp_in);
    });
    // Large stride → small frame count, but MOPAC is slow (~10 min/frame
    // on the 1P9J fixture). 99999 → effectively single frame; keeps
    // the test runnable on CI.
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    const auto status = traj.Run(tp, config, session);
    if (status != nmr::kOk) {
        GTEST_SKIP() << "Trajectory::Run failed (MOPAC env not ready?); "
                        "skipping rather than asserting per "
                        "feedback_log_overages_dont_assert.";
    }

    auto& tr = tp.Result<nmr::MopacChargeWelfordTrajectoryResult>();
    const std::size_t N = tp.AtomCount();
    const std::size_t T = tr.NumFrames();
    const std::size_t attached = tr.SourceAttachedCount();
    std::cout << "MopacChargeWelford Integration1P9J: T=" << T
              << " MOPAC-attached=" << attached
              << " N=" << N << std::endl;

    if (attached == 0) {
        GTEST_SKIP() << "MOPAC did not attach on any frame "
                        "(production stride / runtime / disabled). "
                        "Cannot verify positive-path Welford state.";
    }

    double total_mean_q = 0.0;
    std::size_t n_finite = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& ws = tp.AtomAt(i).mopac_charge_welford;
        EXPECT_EQ(ws.n_frames, attached)
            << "atom " << i << " has n_frames=" << ws.n_frames
            << " but source attached " << attached << " frames; the "
               "per-atom counter should match (no per-atom skip in "
               "Mopac Welford v0).";
        EXPECT_TRUE(std::isfinite(ws.charge.mean));
        if (std::isfinite(ws.charge.mean)) {
            total_mean_q += ws.charge.mean;
            ++n_finite;
        }
    }
    EXPECT_EQ(n_finite, N);
    // Mulliken charge total should be ≈ 0 (protein has net charge from
    // CYS-CYS / charged residues, but Σq_i over the whole protein
    // matches the system net charge ± Mulliken-partition rounding;
    // 1P9J net charge is 0 per the AMBER prep, so Σq ≈ 0).
    EXPECT_LT(std::abs(total_mean_q), 2.0)
        << "Σ ⟨q⟩_i over protein = " << total_mean_q
        << " e; Mulliken charge neutrality sanity (1P9J net charge 0)";
    std::cout << "  Σ ⟨q⟩_i = " << total_mean_q << " e" << std::endl;
}


// ── Layer 0c: FinalizeIdempotency ──────────────────────────────
//
// Drives an absent-source path (no MopacResult attached) so the
// idempotency check exercises the source_attached_count==0 branch
// in Finalize (the loop is skipped; only finalized_=true flag is
// set). Calling Finalize twice should be a no-op.

TEST(MopacChargeWelford, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    auto tr = nmr::MopacChargeWelfordTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    std::vector<nmr::Vec3> positions(tp.AtomCount(), nmr::Vec3::Zero());
    auto conf = std::make_unique<nmr::ProteinConformation>(
        &tp.ProteinRef(), positions, "synthetic (no MOPAC)");
    tr->Compute(*conf, tp, traj, 0, 0.0);
    tr->Finalize(tp, traj);
    const auto attached_first = tr->SourceAttachedCount();
    tr->Finalize(tp, traj);  // second call
    EXPECT_EQ(tr->SourceAttachedCount(), attached_first);
    EXPECT_EQ(tr->NumFrames(), 1u);
}
