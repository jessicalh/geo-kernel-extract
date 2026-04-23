//
// test_gromacs_streaming: tests for the GROMACS trajectory path
// using the new TrajectoryProtein / TrajectoryResult model.
//
// Byte-parity with the pre-refactor output is NOT asserted in this
// session — only BsWelfordTrajectoryResult exists as a concrete
// TrajectoryResult, so most old Welford fields are no longer
// accumulated. Follow-up sessions migrate the remaining fields per
// spec/WIP_OBJECT_MODEL.md Appendix F; byte-parity diffs against the
// pre-refactor fixture restore at that point.
//
// The tests here verify:
//   - TrajectoryProtein builds from a real TPR
//   - GromacsFrameHandler reads frames and finalizes the Protein
//   - BsWelfordTrajectoryResult attaches, dispatches per frame,
//     finalizes, and produces non-zero per-atom fields
//   - Trajectory::Run drives the 5-phase loop end-to-end
//

// AIMNet2Result.h MUST come before GROMACS headers — GROMACS defines
// DIM as a macro which poisons PyTorch template parameters.
#include "AIMNet2Result.h"

// Nanoflann-using headers (SpatialIndexResult, via BiotSavartResult +
// McConnell + others) MUST come before GROMACS headers too — GROMACS
// vectypes.h #define DIM 3 collides with nanoflann's DIM template
// parameter.
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "BiotSavartResult.h"

#include "TrajectoryProtein.h"
#include "Trajectory.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "GromacsFrameHandler.h"
#include "BsWelfordTrajectoryResult.h"
#include "BsShieldingTimeSeriesTrajectoryResult.h"
#include "BsAnomalousAtomMarkerTrajectoryResult.h"
#include "BsT0AutocorrelationTrajectoryResult.h"
#include "BondLengthStatsTrajectoryResult.h"
#include "PositionsTimeSeriesTrajectoryResult.h"
#include "ChiRotamerSelectionTrajectoryResult.h"
#include "DenseBuffer.h"
#include "CalculatorConfig.h"
#include "OperationRunner.h"
#include "OperationLog.h"

#include <gtest/gtest.h>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <set>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

// walker_0/ has md.tpr, md.xtc, md.edr — production layout.
static const std::string TRAJ_DIR =
    std::string(NMR_TEST_DATA_DIR) + "/fleet_test_fullsys/1ZR7_6721/walker_0";


// ============================================================================
// TrajectoryProtein build + handler Open finalizes Protein
// ============================================================================

TEST(GromacsStreaming, TrajectoryBuildAndScan) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);

    if (!fs::exists(TRAJ_DIR + "/md.tpr") ||
        !fs::exists(TRAJ_DIR + "/md.xtc")) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    // --- 1. Build TrajectoryProtein from directory ---
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(TRAJ_DIR)) << tp.Error();

    // Protein not finalized yet — no bonds, no rings
    EXPECT_EQ(tp.ProteinRef().AtomCount(), 479u);
    EXPECT_EQ(tp.ProteinRef().BondCount(), 0u);

    // --- 2. Open handler (pure reader — no frame read yet) ---
    nmr::GromacsFrameHandler handler(tp);
    const std::string tpr = TRAJ_DIR + "/md.tpr";
    const std::string xtc = TRAJ_DIR + "/md.xtc";
    ASSERT_TRUE(handler.Open(xtc, tpr)) << handler.error();

    // Still not finalized — Open does not read.
    EXPECT_EQ(tp.ProteinRef().BondCount(), 0u);

    // --- 3. Read frame 0 + Seed (the canonical-conformation step) ---
    ASSERT_TRUE(handler.ReadNextFrame()) << handler.error();
    tp.Seed(handler.ProteinPositions(), handler.Time());

    // Now Protein is finalized: bonds, rings, conf0, TrajectoryAtoms.
    EXPECT_GT(tp.ProteinRef().BondCount(), 0u);
    EXPECT_GT(tp.ProteinRef().RingCount(), 0u);
    EXPECT_EQ(tp.ProteinRef().ConformationCount(), 1u);
    EXPECT_EQ(tp.AtomCount(), 479u);
    EXPECT_EQ(tp.Atoms().size(), 479u);

    std::cout << "Trajectory: " << tp.ProteinRef().AtomCount() << " atoms, "
              << tp.ProteinRef().BondCount() << " bonds, "
              << tp.ProteinRef().RingCount() << " rings\n";

    // --- 4. Scan 10 more frames (pure reader, no calculators) ---
    size_t scanned = 0;
    for (int i = 0; i < 10; ++i) {
        if (!handler.ReadNextFrame()) break;
        ++scanned;
    }
    ASSERT_GT(scanned, 0u);
    std::cout << "Scanned " << scanned << " frames after frame 0\n";

    // With no TrajectoryResults attached and no OperationRunner calls,
    // no Finalize has run. This test only exercises frame reading and
    // Protein seeding.
    EXPECT_FALSE(tp.IsFinalized());
}


// ============================================================================
// BsWelfordTrajectoryResult end-to-end: attach, dispatch, finalize
// ============================================================================

TEST(GromacsStreaming, BsWelfordAttachAndFinalize) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../../data/calculator_params.toml");

    const std::string tpr = TRAJ_DIR + "/md.tpr";
    const std::string xtc = TRAJ_DIR + "/md.xtc";

    if (!fs::exists(tpr) || !fs::exists(xtc)) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    // Build TrajectoryProtein
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(TRAJ_DIR)) << tp.Error();

    // Per-frame opts — BS only, nothing requiring charges/DSSP/APBS/etc.
    nmr::RunOptions opts;
    opts.charge_source = tp.Charges();
    opts.skip_mopac    = true;
    opts.skip_coulomb  = true;
    opts.skip_apbs     = true;
    opts.skip_dssp     = true;

    // Manual orchestration matching Trajectory::Run's phase order:
    //   1) Open handler (pure reader).
    //   2) Read frame 0, tp.Seed (finalize Protein, create conf0).
    //   3) Attach TrajectoryResults (factories see finalized Protein).
    //   4) Run OperationRunner on conf0, dispatch TRs for frame 0.
    //   5) Per-frame loop: read, tick, run calcs, dispatch.
    //   6) FinalizeAllResults.

    nmr::GromacsFrameHandler handler(tp);
    ASSERT_TRUE(handler.Open(xtc, tpr)) << handler.error();

    ASSERT_TRUE(handler.ReadNextFrame()) << handler.error();
    tp.Seed(handler.ProteinPositions(), handler.Time());

    ASSERT_TRUE(tp.AttachResult(
        nmr::BsWelfordTrajectoryResult::Create(tp)));
    ASSERT_TRUE(tp.HasResult<nmr::BsWelfordTrajectoryResult>());

    // Stand up a Trajectory so the dispatch API has a run-scope
    // context to pass through. BsWelford doesn't use it; future
    // Results that do (selection emitters, env readers) need it.
    nmr::Trajectory traj(xtc, tpr, TRAJ_DIR + "/md.edr");

    // Frame 0 — run per-frame calculators on the canonical
    // conformation, then dispatch TRs.
    {
        auto& conf0 = tp.CanonicalConformation();
        auto rr = nmr::OperationRunner::Run(conf0, opts);
        ASSERT_TRUE(rr.Ok()) << rr.error;
        tp.DispatchCompute(conf0, traj, 0, handler.Time());
    }

    // Run 10 more frames — read, tick, calc, dispatch.
    size_t total = 1;
    while (total < 11 && handler.ReadNextFrame()) {
        auto conf = tp.TickConformation(handler.ProteinPositions());
        auto rr = nmr::OperationRunner::Run(*conf, opts);
        ASSERT_TRUE(rr.Ok()) << rr.error;
        tp.DispatchCompute(*conf, traj, handler.Index(), handler.Time());
        ++total;
    }

    // Finalize: converts Welford M2 → std.
    tp.FinalizeAllResults(traj);
    EXPECT_TRUE(tp.IsFinalized());

    // ── Singleton check on the attach discipline ──
    // A second BsWelfordTrajectoryResult attach must be rejected.
    EXPECT_FALSE(tp.AttachResult(
        nmr::BsWelfordTrajectoryResult::Create(tp)))
        << "second attach should be rejected (singleton per type)";

    // ── Bridge shape: frame counts are consistent per atom ──
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& a = tp.AtomAt(i);
        EXPECT_EQ(a.bs_n_frames, total) << "atom " << i;
        EXPECT_EQ(a.bs_t0_delta_n, total - 1) << "atom " << i;
    }

    // ── Numerical invariants across every atom ──
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& a = tp.AtomAt(i);
        // Finite: no NaN / Inf leaking through the Welford update.
        EXPECT_TRUE(std::isfinite(a.bs_t0_mean))  << "atom " << i;
        EXPECT_TRUE(std::isfinite(a.bs_t0_std))   << "atom " << i;
        EXPECT_TRUE(std::isfinite(a.bs_t2mag_mean)) << "atom " << i;
        EXPECT_TRUE(std::isfinite(a.bs_t2mag_std))  << "atom " << i;
        // Non-negative magnitudes.
        EXPECT_GE(a.bs_t0_std, 0.0)    << "atom " << i;
        EXPECT_GE(a.bs_t2mag_std, 0.0) << "atom " << i;
        EXPECT_GE(a.bs_t2mag_mean, 0.0) << "atom " << i
            << " — |T2| cannot be negative";
        EXPECT_GE(a.bs_t2mag_min, 0.0) << "atom " << i;
        // Order: min ≤ mean ≤ max.
        EXPECT_LE(a.bs_t0_min, a.bs_t0_mean) << "atom " << i;
        EXPECT_LE(a.bs_t0_mean, a.bs_t0_max) << "atom " << i;
        EXPECT_LE(a.bs_t2mag_min, a.bs_t2mag_mean) << "atom " << i;
        EXPECT_LE(a.bs_t2mag_mean, a.bs_t2mag_max) << "atom " << i;
        // min/max frame indices point at an actual frame.
        EXPECT_LT(a.bs_t0_min_frame, total) << "atom " << i;
        EXPECT_LT(a.bs_t0_max_frame, total) << "atom " << i;
    }

    // ── Evidence the bridge carried ConformationAtom → TrajectoryAtom ──
    // At least one atom must see measurable BS contribution AND show
    // non-zero variance across frames (otherwise the trajectory moved
    // but we silently accumulated constants). |T2| magnitude is the
    // cleanest probe — always ≥ 0 and ring-proximity-dependent.
    bool found_active_atom = false;
    bool found_varying_atom = false;
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& a = tp.AtomAt(i);
        if (a.bs_t2mag_mean > 0.01) found_active_atom = true;
        if (a.bs_t2mag_std  > 1e-6) found_varying_atom = true;
        if (std::abs(a.bs_t0_mean) > 0.01 && a.bs_t0_std > 1e-8) {
            std::cout << "Atom " << i << " ("
                      << tp.ProteinRef().AtomAt(i).pdb_atom_name << "): "
                      << "bs_t0_mean=" << a.bs_t0_mean
                      << "  bs_t0_std=" << a.bs_t0_std
                      << "  bs_t2mag_mean=" << a.bs_t2mag_mean
                      << "  bs_t0_delta_mean=" << a.bs_t0_delta_mean << "\n";
        }
    }
    EXPECT_TRUE(found_active_atom)
        << "expected at least one atom with measurable BS |T2|";
    EXPECT_TRUE(found_varying_atom)
        << "expected at least one atom with non-trivial |T2| variance "
           "across frames (otherwise bridge read stale values)";
}


// ============================================================================
// Stride cadence: RunContext::SetStride controls the frame loop gap.
// Tests only Trajectory's Skip/Next interleaving — no AIMNet2, no APBS.
// Attaches BsWelfordTrajectoryResult explicitly so we can count frames
// dispatched without pulling in the full PerFrameExtractionSet.
// ============================================================================

TEST(GromacsStreaming, TrajectoryStrideSkipsFramesCorrectly) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../../data/calculator_params.toml");

    const std::string tpr = TRAJ_DIR + "/md.tpr";
    const std::string xtc = TRAJ_DIR + "/md.xtc";
    const std::string edr = TRAJ_DIR + "/md.edr";

    if (!fs::exists(tpr) || !fs::exists(xtc)) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    // Ad-hoc RunConfiguration for this test: BiotSavart only, no
    // AIMNet2. Mimics ScanForDftPointSet's shape without using
    // production factories (those add more ConformationResults than
    // needed to test the stride).
    nmr::RunConfiguration config;
    config.SetName("StrideTest");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_apbs    = true;
    opts.skip_coulomb = true;
    opts.skip_dssp    = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::BiotSavartResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::BsWelfordTrajectoryResult::Create(tp);
        });

    // Run with stride=3 — dispatch every 3rd frame (frames 0, 3, 6, …).
    // Stride lives on the RunConfiguration (the named run shape).
    config.SetStride(3);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(TRAJ_DIR)) << tp.Error();

    nmr::Trajectory traj(xtc, tpr, edr);
    nmr::Session session;  // no aimnet2 needed here; StrideTest config doesn't require it
    const nmr::Status s = traj.Run(tp, config, session, /*extras=*/{}, /*output_dir=*/"/tmp");
    ASSERT_EQ(s, nmr::kOk);

    EXPECT_TRUE(traj.IsComplete());
    EXPECT_GT(traj.FrameCount(), 1u) << "stride=3 should still produce "
                                        "multiple dispatched frames";

    // Frame indices should be at stride multiples: 0, 3, 6, 9, …
    const auto& idx = traj.FrameIndices();
    ASSERT_EQ(idx.size(), traj.FrameCount());
    EXPECT_EQ(idx[0], 0u) << "frame 0 always dispatched";
    for (size_t k = 1; k < idx.size(); ++k) {
        EXPECT_EQ(idx[k] - idx[k - 1], 3u)
            << "between dispatched frames " << (k - 1) << " and " << k;
    }

    // BsWelford should have been dispatched once per frame in idx.
    EXPECT_EQ(tp.AtomAt(0).bs_n_frames, traj.FrameCount());

    std::cout << "Stride 3: " << traj.FrameCount()
              << " frames dispatched, first indices = [";
    for (size_t k = 0; k < std::min<size_t>(5, idx.size()); ++k) {
        std::cout << (k ? ", " : "") << idx[k];
    }
    std::cout << (idx.size() > 5 ? ", …" : "") << "]\n";
}


// ============================================================================
// Trajectory::Run 5-phase drive (end-to-end via PerFrameExtractionSet)
// ============================================================================

// ============================================================================
// PositionsTimeSeriesTrajectoryResult — dense-buffer worked example
//
// The Finalize-only time-series pattern. After a short run,
// TrajectoryProtein holds a DenseBuffer<Vec3> owned by this Result's
// type_index. The buffer's atom slice for each atom equals the
// positions seen in each frame, atom-major.
// ============================================================================

TEST(GromacsStreaming, PositionsTimeSeriesEndToEnd) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../../data/calculator_params.toml");

    const std::string tpr = TRAJ_DIR + "/md.tpr";
    const std::string xtc = TRAJ_DIR + "/md.xtc";
    if (!fs::exists(tpr) || !fs::exists(xtc)) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(TRAJ_DIR)) << tp.Error();

    nmr::RunOptions opts;
    opts.charge_source = tp.Charges();
    opts.skip_mopac    = true;
    opts.skip_coulomb  = true;
    opts.skip_apbs     = true;
    opts.skip_dssp     = true;

    nmr::GromacsFrameHandler handler(tp);
    ASSERT_TRUE(handler.Open(xtc, tpr)) << handler.error();
    ASSERT_TRUE(handler.ReadNextFrame()) << handler.error();
    tp.Seed(handler.ProteinPositions(), handler.Time());

    // Record frame 0's positions BEFORE attach (for later verification).
    const std::vector<nmr::Vec3> positions_frame0 = handler.ProteinPositions();

    ASSERT_TRUE(tp.AttachResult(
        nmr::PositionsTimeSeriesTrajectoryResult::Create(tp)));

    nmr::Trajectory traj(xtc, tpr, TRAJ_DIR + "/md.edr");

    // Frame 0.
    {
        auto& conf0 = tp.CanonicalConformation();
        auto rr = nmr::OperationRunner::Run(conf0, opts);
        ASSERT_TRUE(rr.Ok()) << rr.error;
        tp.DispatchCompute(conf0, traj, 0, handler.Time());
    }

    // 9 more frames.
    size_t total = 1;
    while (total < 10 && handler.ReadNextFrame()) {
        auto conf = tp.TickConformation(handler.ProteinPositions());
        auto rr = nmr::OperationRunner::Run(*conf, opts);
        ASSERT_TRUE(rr.Ok()) << rr.error;
        tp.DispatchCompute(*conf, traj, handler.Index(), handler.Time());
        ++total;
    }

    tp.FinalizeAllResults(traj);
    EXPECT_TRUE(tp.IsFinalized());

    // ── Post-Finalize: dense buffer is owned by tp ──
    auto* buffer = tp.GetDenseBuffer<nmr::Vec3>(
        std::type_index(
            typeid(nmr::PositionsTimeSeriesTrajectoryResult)));
    ASSERT_NE(buffer, nullptr)
        << "PositionsTimeSeriesTrajectoryResult did not transfer its buffer";

    EXPECT_EQ(buffer->AtomCount(),     tp.AtomCount());
    EXPECT_EQ(buffer->StridePerAtom(), total)
        << "stride should equal number of frames dispatched";

    // Frame 0's positions survive round-trip.
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const nmr::Vec3& got = buffer->At(i, 0);
        const nmr::Vec3& exp = positions_frame0[i];
        EXPECT_DOUBLE_EQ(got.x(), exp.x()) << "atom " << i << " frame 0 x";
        EXPECT_DOUBLE_EQ(got.y(), exp.y()) << "atom " << i << " frame 0 y";
        EXPECT_DOUBLE_EQ(got.z(), exp.z()) << "atom " << i << " frame 0 z";
    }

    // Atoms actually move across frames — at least one atom has some
    // non-trivial displacement between frames 0 and N-1.
    double max_frame_displacement = 0.0;
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const nmr::Vec3 d = buffer->At(i, total - 1) - buffer->At(i, 0);
        const double dist = d.norm();
        if (dist > max_frame_displacement) max_frame_displacement = dist;
    }
    EXPECT_GT(max_frame_displacement, 1e-6)
        << "expected some atomic motion across the sampled frames";

    std::cout << "PositionsTimeSeries: " << total << " frames x "
              << tp.AtomCount() << " atoms stored in DenseBuffer<Vec3>; "
              << "max atom displacement over window: "
              << max_frame_displacement << " Å\n";
}


// ============================================================================
// BsShieldingTimeSeriesTrajectoryResult — SphericalTensor time-series
// worked example. Same Finalize-only dense-buffer pattern as
// PositionsTimeSeries, but the payload is a 9-component
// SphericalTensor per atom per frame. Emission carries the e3nn
// metadata (irrep_layout / normalization / parity) so Python
// consumers can construct Irreps("0e+1o+2e") directly.
// ============================================================================

TEST(GromacsStreaming, BsShieldingTimeSeriesEndToEnd) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../../data/calculator_params.toml");

    const std::string tpr = TRAJ_DIR + "/md.tpr";
    const std::string xtc = TRAJ_DIR + "/md.xtc";
    if (!fs::exists(tpr) || !fs::exists(xtc)) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(TRAJ_DIR)) << tp.Error();

    // BiotSavart must run per frame so that
    // ca.bs_shielding_contribution is populated. It doesn't require
    // DSSP, charges, APBS, or MOPAC — skip everything except what BS
    // needs (Geometry, SpatialIndex, Enrichment).
    nmr::RunOptions opts;
    opts.charge_source = tp.Charges();
    opts.skip_mopac    = true;
    opts.skip_coulomb  = true;
    opts.skip_apbs     = true;
    opts.skip_dssp     = true;

    nmr::GromacsFrameHandler handler(tp);
    ASSERT_TRUE(handler.Open(xtc, tpr)) << handler.error();
    ASSERT_TRUE(handler.ReadNextFrame()) << handler.error();
    tp.Seed(handler.ProteinPositions(), handler.Time());

    ASSERT_TRUE(tp.AttachResult(
        nmr::BsShieldingTimeSeriesTrajectoryResult::Create(tp)));

    nmr::Trajectory traj(xtc, tpr, TRAJ_DIR + "/md.edr");

    // Frame 0.
    nmr::SphericalTensor expected_atom_frame0;
    {
        auto& conf0 = tp.CanonicalConformation();
        auto rr = nmr::OperationRunner::Run(conf0, opts);
        ASSERT_TRUE(rr.Ok()) << rr.error;
        // Snapshot frame-0 bs_shielding for one atom BEFORE dispatch
        // so we can verify round-trip into the dense buffer.
        expected_atom_frame0 = conf0.AtomAt(0).bs_shielding_contribution;
        tp.DispatchCompute(conf0, traj, 0, handler.Time());
    }

    // 5 more frames — enough to exercise the growing buffer + final
    // transfer without running long.
    size_t total = 1;
    while (total < 6 && handler.ReadNextFrame()) {
        auto conf = tp.TickConformation(handler.ProteinPositions());
        auto rr = nmr::OperationRunner::Run(*conf, opts);
        ASSERT_TRUE(rr.Ok()) << rr.error;
        tp.DispatchCompute(*conf, traj, handler.Index(), handler.Time());
        ++total;
    }

    tp.FinalizeAllResults(traj);
    EXPECT_TRUE(tp.IsFinalized());

    // ── Dense buffer owned by tp, indexed by SphericalTensor ──
    auto* buffer = tp.GetDenseBuffer<nmr::SphericalTensor>(
        std::type_index(
            typeid(nmr::BsShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buffer, nullptr)
        << "BsShieldingTimeSeriesTrajectoryResult did not transfer";

    EXPECT_EQ(buffer->AtomCount(),     tp.AtomCount());
    EXPECT_EQ(buffer->StridePerAtom(), total);

    // Frame 0's SphericalTensor survives round-trip bit-identical.
    const nmr::SphericalTensor& got = buffer->At(0, 0);
    EXPECT_DOUBLE_EQ(got.T0, expected_atom_frame0.T0);
    for (int k = 0; k < 3; ++k)
        EXPECT_DOUBLE_EQ(got.T1[k], expected_atom_frame0.T1[k])
            << "atom 0 frame 0 T1[" << k << "]";
    for (int k = 0; k < 5; ++k)
        EXPECT_DOUBLE_EQ(got.T2[k], expected_atom_frame0.T2[k])
            << "atom 0 frame 0 T2[" << k << "]";

    // At least one atom shows a non-trivial |T2| across the window
    // (the protein moves; ring-proximity-dependent BS shielding varies).
    bool found_varying = false;
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& st0 = buffer->At(i, 0);
        const auto& stN = buffer->At(i, total - 1);
        double delta_t2 = 0.0;
        for (int k = 0; k < 5; ++k) {
            const double d = stN.T2[k] - st0.T2[k];
            delta_t2 += d * d;
        }
        if (std::sqrt(delta_t2) > 1e-8) { found_varying = true; break; }
    }
    EXPECT_TRUE(found_varying)
        << "expected non-trivial T2 change across frames "
           "(ring-proximity-dependent shielding varies with motion)";

    std::cout << "BsShieldingTimeSeries: " << total
              << " frames x " << tp.AtomCount()
              << " atoms of SphericalTensor stored in "
                 "DenseBuffer<SphericalTensor>\n";
}


// ============================================================================
// ChiRotamerSelectionTrajectoryResult — scan-mode emitter worked example
//
// Pushes SelectionRecords to traj.MutableSelections() when any
// residue's chi angle crosses a rotamer bin boundary. Over a short
// window we expect at least a few transitions and they must show up
// in the run-scope bag under this Result's kind.
// ============================================================================

TEST(GromacsStreaming, ChiRotamerSelectionEmitsToBag) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../../data/calculator_params.toml");

    const std::string tpr = TRAJ_DIR + "/md.tpr";
    const std::string xtc = TRAJ_DIR + "/md.xtc";
    if (!fs::exists(tpr) || !fs::exists(xtc)) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(TRAJ_DIR)) << tp.Error();

    nmr::RunOptions opts;
    opts.charge_source = tp.Charges();
    opts.skip_mopac    = true;
    opts.skip_coulomb  = true;
    opts.skip_apbs     = true;
    opts.skip_dssp     = true;

    nmr::GromacsFrameHandler handler(tp);
    ASSERT_TRUE(handler.Open(xtc, tpr)) << handler.error();
    ASSERT_TRUE(handler.ReadNextFrame()) << handler.error();
    tp.Seed(handler.ProteinPositions(), handler.Time());

    ASSERT_TRUE(tp.AttachResult(
        nmr::ChiRotamerSelectionTrajectoryResult::Create(tp)));

    nmr::Trajectory traj(xtc, tpr, TRAJ_DIR + "/md.edr");

    // Bag is empty before any frame runs.
    EXPECT_EQ(traj.Selections().Count(), 0u);

    // Frame 0: establishes prior bins, no transitions emitted.
    {
        auto& conf0 = tp.CanonicalConformation();
        auto rr = nmr::OperationRunner::Run(conf0, opts);
        ASSERT_TRUE(rr.Ok()) << rr.error;
        tp.DispatchCompute(conf0, traj, 0, handler.Time());
    }
    EXPECT_EQ(traj.Selections().Count(), 0u)
        << "frame 0 should not emit — no prior bin to compare against";

    // Run enough subsequent frames to observe some rotamer transitions.
    size_t total = 1;
    while (total < 60 && handler.ReadNextFrame()) {
        auto conf = tp.TickConformation(handler.ProteinPositions());
        auto rr = nmr::OperationRunner::Run(*conf, opts);
        ASSERT_TRUE(rr.Ok()) << rr.error;
        tp.DispatchCompute(*conf, traj, handler.Index(), handler.Time());
        ++total;
    }

    tp.FinalizeAllResults(traj);

    // The bag has events of exactly one kind — ChiRotamerSelection's.
    const auto kinds = traj.Selections().Kinds();
    ASSERT_EQ(kinds.size(), 1u) << "only one emitter is attached";
    EXPECT_EQ(kinds[0],
        std::type_index(typeid(nmr::ChiRotamerSelectionTrajectoryResult)));

    // There should be at least ONE transition across ~60 frames with
    // hundreds of chi angles on the protein.
    const auto records =
        traj.Selections().ByKind<nmr::ChiRotamerSelectionTrajectoryResult>();
    EXPECT_GT(records.size(), 0u)
        << "expected at least one rotamer transition across "
        << total << " frames";

    // Each record carries frame + time as explicit top-level fields,
    // and emitter-specific metadata.
    for (const auto* rec : records) {
        EXPECT_LT(rec->frame_idx, total);
        EXPECT_GT(rec->time_ps, 0.0);
        EXPECT_FALSE(rec->reason.empty())
            << "reason should describe the transition";
        EXPECT_TRUE(rec->metadata.count("residue_index"))
            << "metadata should carry residue_index";
        EXPECT_TRUE(rec->metadata.count("chi_index"));
        EXPECT_TRUE(rec->metadata.count("bin_before"));
        EXPECT_TRUE(rec->metadata.count("bin_after"));
        EXPECT_NE(rec->metadata.at("bin_before"),
                  rec->metadata.at("bin_after"))
            << "a transition must change the bin";
    }

    // Sample affordance: ByKindSinceFrame filter.
    const auto late =
        traj.Selections().ByKindSinceFrame<
            nmr::ChiRotamerSelectionTrajectoryResult>(total / 2);
    for (const auto* rec : late) {
        EXPECT_GE(rec->frame_idx, total / 2);
    }
    EXPECT_LE(late.size(), records.size());

    std::cout << "ChiRotamerSelection: " << records.size()
              << " transitions emitted across " << total
              << " frames (first-half: " << (records.size() - late.size())
              << ", second-half: " << late.size() << ")\n";
}


// ============================================================================
// BondLengthStatsTrajectoryResult — per-bond scope worked example.
//
// Per-bond Welford state lives INTERNAL to the Result (option b per
// WIP §3), not on a first-class TrajectoryBond store. Assert on
// physical plausibility: bond lengths are ~1-2 Å, std >= 0, atoms
// actually move across frames so std > 0 for some bonds.
// ============================================================================

TEST(GromacsStreaming, BondLengthStatsEndToEnd) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../../data/calculator_params.toml");

    const std::string tpr = TRAJ_DIR + "/md.tpr";
    const std::string xtc = TRAJ_DIR + "/md.xtc";
    if (!fs::exists(tpr) || !fs::exists(xtc)) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(TRAJ_DIR)) << tp.Error();

    nmr::RunOptions opts;
    opts.charge_source = tp.Charges();
    opts.skip_mopac    = true;
    opts.skip_coulomb  = true;
    opts.skip_apbs     = true;
    opts.skip_dssp     = true;

    nmr::GromacsFrameHandler handler(tp);
    ASSERT_TRUE(handler.Open(xtc, tpr)) << handler.error();
    ASSERT_TRUE(handler.ReadNextFrame()) << handler.error();
    tp.Seed(handler.ProteinPositions(), handler.Time());

    // After Seed, Protein is finalized — BondCount is valid here.
    ASSERT_GT(tp.ProteinRef().BondCount(), 0u);

    ASSERT_TRUE(tp.AttachResult(
        nmr::BondLengthStatsTrajectoryResult::Create(tp)));

    nmr::Trajectory traj(xtc, tpr, TRAJ_DIR + "/md.edr");

    // Frame 0 + 9 more.
    auto& conf0 = tp.CanonicalConformation();
    auto rr0 = nmr::OperationRunner::Run(conf0, opts);
    ASSERT_TRUE(rr0.Ok()) << rr0.error;
    tp.DispatchCompute(conf0, traj, 0, handler.Time());
    size_t total = 1;
    while (total < 10 && handler.ReadNextFrame()) {
        auto conf = tp.TickConformation(handler.ProteinPositions());
        auto rr = nmr::OperationRunner::Run(*conf, opts);
        ASSERT_TRUE(rr.Ok()) << rr.error;
        tp.DispatchCompute(*conf, traj, handler.Index(), handler.Time());
        ++total;
    }
    tp.FinalizeAllResults(traj);

    ASSERT_TRUE(tp.HasResult<nmr::BondLengthStatsTrajectoryResult>());
    const auto& bonds =
        tp.Result<nmr::BondLengthStatsTrajectoryResult>().PerBond();
    ASSERT_EQ(bonds.size(), tp.ProteinRef().BondCount());

    // Physical plausibility. Every covalent bond length in a protein
    // is 0.8 < L < 2.5 Å (C-H to S-S), and every bond was observed
    // every frame.
    bool any_moved = false;
    for (const auto& pb : bonds) {
        EXPECT_EQ(pb.n_frames, total);
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
              << " bonds x " << total << " frames; sample bond 0 mean="
              << bonds[0].length_mean << "Å std=" << bonds[0].length_std
              << "Å\n";
}


// ============================================================================
// BsAnomalousAtomMarkerTrajectoryResult — cross-Result + atom events.
//
// Exercises: attach-order dependency on another TrajectoryResult,
// reading that Result's stashed TrajectoryAtom fields during Compute,
// pushing to the per-atom events bag with two kinds (HighT0, LowT0)
// from one emitter.
// ============================================================================

TEST(GromacsStreaming, BsAnomalousAtomMarkerPushesEvents) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../../data/calculator_params.toml");

    const std::string tpr = TRAJ_DIR + "/md.tpr";
    const std::string xtc = TRAJ_DIR + "/md.xtc";
    if (!fs::exists(tpr) || !fs::exists(xtc)) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(TRAJ_DIR)) << tp.Error();

    nmr::RunOptions opts;
    opts.charge_source = tp.Charges();
    opts.skip_mopac    = true;
    opts.skip_coulomb  = true;
    opts.skip_apbs     = true;
    opts.skip_dssp     = true;

    nmr::GromacsFrameHandler handler(tp);
    ASSERT_TRUE(handler.Open(xtc, tpr)) << handler.error();
    ASSERT_TRUE(handler.ReadNextFrame()) << handler.error();
    tp.Seed(handler.ProteinPositions(), handler.Time());

    // Attach order matters: BsWelford first so its fields on
    // TrajectoryAtom exist when BsAnomalous reads them.
    ASSERT_TRUE(tp.AttachResult(
        nmr::BsWelfordTrajectoryResult::Create(tp)));
    ASSERT_TRUE(tp.AttachResult(
        nmr::BsAnomalousAtomMarkerTrajectoryResult::Create(tp)));

    nmr::Trajectory traj(xtc, tpr, TRAJ_DIR + "/md.edr");

    // Frame 0 + 60 more. Burn-in of 20 frames per atom means the
    // first 20 frames accumulate but produce no events.
    auto& conf0 = tp.CanonicalConformation();
    auto rr0 = nmr::OperationRunner::Run(conf0, opts);
    ASSERT_TRUE(rr0.Ok()) << rr0.error;
    tp.DispatchCompute(conf0, traj, 0, handler.Time());
    size_t total = 1;
    while (total < 60 && handler.ReadNextFrame()) {
        auto conf = tp.TickConformation(handler.ProteinPositions());
        auto rr = nmr::OperationRunner::Run(*conf, opts);
        ASSERT_TRUE(rr.Ok()) << rr.error;
        tp.DispatchCompute(*conf, traj, handler.Index(), handler.Time());
        ++total;
    }
    tp.FinalizeAllResults(traj);

    const auto& marker =
        tp.Result<nmr::BsAnomalousAtomMarkerTrajectoryResult>();
    const size_t n_events = marker.TotalEvents();

    // We ran 60 frames, MIN_BURN_IN_FRAMES is 20, so frames 20..59
    // contribute. Across ~479 atoms × 40 frames × ~4.5% two-tail of
    // a 2σ threshold, expect roughly hundreds of events. Exact count
    // depends on motion; assert non-zero and that events fell AFTER
    // the burn-in window.
    EXPECT_GT(n_events, 0u)
        << "expected some anomalies across " << total
        << " frames with burn-in=20";
    EXPECT_EQ(n_events, marker.HighEvents() + marker.LowEvents())
        << "totals should split between High and Low kinds";

    // Inspect atom events bag on one atom that got flagged.
    size_t inspected = 0;
    for (size_t i = 0; i < tp.AtomCount() && inspected < 3; ++i) {
        const auto& bag = tp.AtomAt(i).events;
        if (bag.Count() == 0) continue;
        ++inspected;

        // Every event on this atom came from BsAnomalousAtomMarker.
        for (const auto& rec : bag.All()) {
            EXPECT_EQ(rec.emitter, std::type_index(
                typeid(nmr::BsAnomalousAtomMarkerTrajectoryResult)));
            EXPECT_TRUE(
                rec.kind == std::type_index(typeid(
                    nmr::BsAnomalousAtomMarkerTrajectoryResult::
                        BsAnomalyHighT0)) ||
                rec.kind == std::type_index(typeid(
                    nmr::BsAnomalousAtomMarkerTrajectoryResult::
                        BsAnomalyLowT0)));
            EXPECT_GE(rec.frame_idx, 20u) << "event before burn-in end";
            EXPECT_LT(rec.frame_idx, total);
            EXPECT_TRUE(rec.metadata.count("z_sigma"));
            EXPECT_TRUE(rec.metadata.count("current_t0"));
            EXPECT_TRUE(rec.metadata.count("running_mean"));
            EXPECT_TRUE(rec.metadata.count("running_n"));
        }

        // The per-atom bag's own affordance queries. Count by kind
        // should partition the bag.
        const size_t high = bag.CountByKind<
            nmr::BsAnomalousAtomMarkerTrajectoryResult::BsAnomalyHighT0>();
        const size_t low = bag.CountByKind<
            nmr::BsAnomalousAtomMarkerTrajectoryResult::BsAnomalyLowT0>();
        EXPECT_EQ(high + low, bag.Count());
    }
    EXPECT_GT(inspected, 0u) << "no atom had any events — unexpected";

    std::cout << "BsAnomalousAtomMarker: " << n_events
              << " total events (" << marker.HighEvents()
              << " high, " << marker.LowEvents()
              << " low) across " << total << " frames\n";
}


// ============================================================================
// BsT0AutocorrelationTrajectoryResult — rolling-window autocorrelation.
//
// Exercises: non-trivial circular-buffer internal state, Finalize-
// heavy normalization, DenseBuffer<double> emission with
// per-atom × lag shape.
// ============================================================================

TEST(GromacsStreaming, BsT0AutocorrelationEndToEnd) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../../data/calculator_params.toml");

    const std::string tpr = TRAJ_DIR + "/md.tpr";
    const std::string xtc = TRAJ_DIR + "/md.xtc";
    if (!fs::exists(tpr) || !fs::exists(xtc)) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(TRAJ_DIR)) << tp.Error();

    nmr::RunOptions opts;
    opts.charge_source = tp.Charges();
    opts.skip_mopac    = true;
    opts.skip_coulomb  = true;
    opts.skip_apbs     = true;
    opts.skip_dssp     = true;

    nmr::GromacsFrameHandler handler(tp);
    ASSERT_TRUE(handler.Open(xtc, tpr)) << handler.error();
    ASSERT_TRUE(handler.ReadNextFrame()) << handler.error();
    tp.Seed(handler.ProteinPositions(), handler.Time());

    ASSERT_TRUE(tp.AttachResult(
        nmr::BsT0AutocorrelationTrajectoryResult::Create(tp)));

    nmr::Trajectory traj(xtc, tpr, TRAJ_DIR + "/md.edr");

    // Need > a handful of frames to get any lag > 0 populated; run
    // 30 which gives lag 0..29 (N_LAGS=120 but only the first 30 get
    // counts).
    auto& conf0 = tp.CanonicalConformation();
    auto rr0 = nmr::OperationRunner::Run(conf0, opts);
    ASSERT_TRUE(rr0.Ok()) << rr0.error;
    tp.DispatchCompute(conf0, traj, 0, handler.Time());
    size_t total = 1;
    while (total < 30 && handler.ReadNextFrame()) {
        auto conf = tp.TickConformation(handler.ProteinPositions());
        auto rr = nmr::OperationRunner::Run(*conf, opts);
        ASSERT_TRUE(rr.Ok()) << rr.error;
        tp.DispatchCompute(*conf, traj, handler.Index(), handler.Time());
        ++total;
    }
    tp.FinalizeAllResults(traj);

    auto* rho = tp.GetDenseBuffer<double>(std::type_index(
        typeid(nmr::BsT0AutocorrelationTrajectoryResult)));
    ASSERT_NE(rho, nullptr);
    EXPECT_EQ(rho->AtomCount(), tp.AtomCount());
    EXPECT_EQ(rho->StridePerAtom(),
        nmr::BsT0AutocorrelationTrajectoryResult::N_LAGS);

    // ρ(0) must be exactly 1 for every atom that has non-zero
    // variance (the lag-0 autocorrelation is self-correlation). For
    // atoms with constant T0 (e.g. buried atoms far from rings), the
    // variance is 0 and ρ(0) falls through to 0 — our Finalize
    // special-cases that.
    size_t n_rho0_unit = 0;
    size_t n_rho0_zero = 0;
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const double r0 = rho->At(i, 0);
        if (std::fabs(r0 - 1.0) < 1e-9) ++n_rho0_unit;
        else if (std::fabs(r0) < 1e-9)  ++n_rho0_zero;
        else ADD_FAILURE() << "atom " << i << " rho(0)=" << r0
                           << " — expected 1.0 or 0.0";
    }
    EXPECT_GT(n_rho0_unit, 0u)
        << "expected some atoms to show self-correlation ρ(0)=1";

    // ρ(k) ∈ [−1, +1] at every k and every N — guaranteed by the
    // biased estimator via Cauchy-Schwarz. This is why we picked
    // biased over unbiased: the physical-meaning bound holds
    // exactly at finite sample size. Check every lag, every atom.
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        for (size_t k = 0; k <
             nmr::BsT0AutocorrelationTrajectoryResult::N_LAGS; ++k) {
            const double v = rho->At(i, k);
            EXPECT_GE(v, -1.0 - 1e-9)
                << "atom " << i << " lag " << k << " out of bounds";
            EXPECT_LE(v,  1.0 + 1e-9)
                << "atom " << i << " lag " << k << " out of bounds";
        }
    }

    // Lags beyond (total - 1) have zero pairs contributed — ρ(k)
    // there is 0 by construction.
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        for (size_t k = total; k <
             nmr::BsT0AutocorrelationTrajectoryResult::N_LAGS; ++k) {
            EXPECT_NEAR(rho->At(i, k), 0.0, 1e-12)
                << "atom " << i << " lag " << k
                << " should be 0 (beyond trajectory length)";
        }
    }

    std::cout << "BsT0Autocorrelation: " << tp.AtomCount()
              << " atoms × " << rho->StridePerAtom()
              << " lags; " << n_rho0_unit
              << " atoms with ρ(0)=1, " << n_rho0_zero
              << " with constant T0 (ρ(0)=0)\n";
}


TEST(GromacsStreaming, TrajectoryRunDrivesLoop) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../../data/calculator_params.toml");

    const std::string tpr = TRAJ_DIR + "/md.tpr";
    const std::string xtc = TRAJ_DIR + "/md.xtc";
    const std::string edr = TRAJ_DIR + "/md.edr";
    const std::string jpt = std::string(NMR_TEST_DATA_DIR) +
                            "/../../data/models/aimnet2_wb97m_0.jpt";

    if (!fs::exists(tpr) || !fs::exists(xtc)) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    // PerFrameExtractionSet makes AIMNet2 mandatory — the
    // configuration's RequiresAimnet2() flag is honoured by
    // Trajectory::Run Phase 2. Load the model into Session once; the
    // run reads session.Aimnet2Model() there.
    ASSERT_TRUE(fs::exists(jpt))
        << "AIMNet2 model not found at " << jpt
        << " — PerFrameExtractionSet requires it";

    nmr::Session session;
    ASSERT_EQ(session.LoadAimnet2Model(jpt), nmr::kOk)
        << session.LastError();

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(TRAJ_DIR)) << tp.Error();

    nmr::Trajectory traj(xtc, tpr, edr);
    nmr::RunConfiguration config = nmr::RunConfiguration::PerFrameExtractionSet();

    // Full-trajectory run: every per-frame ConformationResult (APBS,
    // AIMNet2, classical kernels, SASA, water, hydration, EEQ,
    // bonded energy, DSSP, BS, HM, McConnell, …) plus the attached
    // BsWelfordTrajectoryResult. Stride 1 (every frame) — production
    // analysis mode uses STRIDE=2, not supported by Trajectory::Run
    // in this session's landing. This is long-running by design.
    const nmr::Status s = traj.Run(tp, config, session, /*extras=*/{}, /*output_dir=*/"/tmp");
    ASSERT_EQ(s, nmr::kOk);

    EXPECT_TRUE(traj.IsComplete());
    EXPECT_GT(traj.FrameCount(), 0u);
    EXPECT_TRUE(tp.IsFinalized());
    EXPECT_TRUE(tp.HasResult<nmr::BsWelfordTrajectoryResult>());
    EXPECT_EQ(tp.AtomAt(0).bs_n_frames, traj.FrameCount());

    std::cout << "Trajectory::Run drove " << traj.FrameCount()
              << " frames, " << tp.AtomCount() << " atoms, "
              << traj.Selections().Count() << " selections\n";
}
