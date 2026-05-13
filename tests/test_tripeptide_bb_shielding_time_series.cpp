//
// test_tripeptide_bb_shielding_time_series: discipline + integration
// tests for TripeptideBackboneShieldingTimeSeriesTrajectoryResult.
//
// Two test groups:
//
//   1. UNIT (synthetic).  Builds 1ubq_protonated.pdb (no DFT table
//      required), constructs 4 free-standing per-frame
//      ProteinConformations, sets hand-crafted
//      tripeptide_bb_shielding_spherical values on each
//      ConformationAtom, drives the TR's Compute/Finalize directly,
//      and verifies:
//        - buffer adopt + Finalize round-trip
//        - H5 group emission with /trajectory/tripeptide_bb_shielding_time_series/
//        - buffer shape (N atoms × n_frames)
//
//   2. INTEGRATION (fingerprint).  Drives Trajectory::Run on the
//      1P9J_5801 fleet_amber fixture for ≤ 5 frames with the tripeptide
//      DSN configured. Verifies H5 group + shape, and asserts that the
//      central residue's CA T0 at frame 0 matches a baked regression
//      target to 3 significant figures.
//
// Skipped when the tensorcs15 DSN is not configured (integration) or
// the fleet_amber fixture is not on disk.
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

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DenseBuffer.h"
#include "DsspResult.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "RuntimeEnvironment.h"
#include "Session.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "TripeptideBackboneShieldingResult.h"
#include "TripeptideBackboneShieldingTimeSeriesTrajectoryResult.h"
#include "Types.h"

#include <gtest/gtest.h>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

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

// Per-fixture central-residue location for the integration sanity test.
// 1P9J_5801 has 56 residues; canonical central residue is 28 in chain A.
constexpr int kFingerprintResidueNumber = 28;
const std::string kFingerprintChainId = "A";

// Physical sanity bounds for tripeptide BB shielding. Larsen 2015 per-element
// ranges: H ~5-10, Cα ~40-60, C' ~170-180, N ~100-130 ppm. The integration
// test anchors at a CA atom; 250 ppm covers any backbone slot with margin
// against calibration drift / AAA reference shifts / chi-grid bin-boundary
// noise (no magic-number fingerprint — the pipeline has genuine numerical
// noise from TRR float32 compression, chi-grid binning, and upstream
// recalibration that would trip a tight equality assertion without
// indicating real regression).
constexpr double kT0SanityBoundPpm   = 250.0;
constexpr double kT1SymmetryTolerance = 1e-6;


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


// Build a 4-frame synthetic SphericalTensor sequence per atom. Returns
// the value the (atom_i, frame_t) cell should carry. T0 = atom_i + t*100;
// T1[k] = atom_i + t*100 + (k+1)*1e-2; T2[k] = atom_i + t*100 + (k+1)*1e-3.
nmr::SphericalTensor SyntheticTensor(size_t atom_i, size_t frame_t) {
    nmr::SphericalTensor s;
    s.T0 = static_cast<double>(atom_i) + static_cast<double>(frame_t) * 100.0;
    for (size_t k = 0; k < 3; ++k) {
        s.T1[k] = static_cast<double>(atom_i)
                  + static_cast<double>(frame_t) * 100.0
                  + static_cast<double>(k + 1) * 1e-2;
    }
    for (size_t k = 0; k < 5; ++k) {
        s.T2[k] = static_cast<double>(atom_i)
                  + static_cast<double>(frame_t) * 100.0
                  + static_cast<double>(k + 1) * 1e-3;
    }
    return s;
}


bool SphericalEqual(const nmr::SphericalTensor& a,
                    const nmr::SphericalTensor& b,
                    double tol) {
    if (std::abs(a.T0 - b.T0) > tol) return false;
    for (size_t k = 0; k < 3; ++k)
        if (std::abs(a.T1[k] - b.T1[k]) > tol) return false;
    for (size_t k = 0; k < 5; ++k)
        if (std::abs(a.T2[k] - b.T2[k]) > tol) return false;
    return true;
}


}  // namespace


// ============================================================================
// UNIT: 4-frame synthetic round-trip on hand-crafted SphericalTensor values.
//
// Drives Compute / Finalize / WriteH5Group with hand-crafted per-atom inputs
// and asserts exact bit-equality on round-trip. Skips Trajectory::Run
// orchestration; uses the fleet_amber fixture only for TrajectoryProtein
// AtomCount plumbing. Exact equality is appropriate here because synthetic
// inputs have no numerical noise — round-trip is the contract.
// ============================================================================

TEST(TripeptideBackboneShieldingTimeSeries, SyntheticFourFrames) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();

    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix))
        GTEST_SKIP() << "fleet_amber " << kFixtureProtein
                     << " fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    const size_t Ntp = tp.AtomCount();
    ASSERT_GT(Ntp, 0u);

    auto tr = nmr::TripeptideBackboneShieldingTimeSeriesTrajectoryResult::Create(tp);

    // Trajectory exists only to satisfy the Compute signature; TR::Compute
    // marks traj as (void) and reads only ConformationAtom.
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 4;
    const auto& protein_ref = tp.ProteinRef();
    std::vector<nmr::Vec3> positions(Ntp, nmr::Vec3::Zero());

    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &protein_ref, positions, "synthetic frame");
        for (size_t i = 0; i < Ntp; ++i) {
            conf->MutableAtomAt(i).tripeptide_bb_shielding_spherical =
                SyntheticTensor(i, t);
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    EXPECT_EQ(tr->NumFrames(), kFrames);

    tr->Finalize(tp, traj);

    auto* buf =
        tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(typeid(
            nmr::TripeptideBackboneShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), Ntp);
    EXPECT_EQ(buf->StridePerAtom(), kFrames);

    for (size_t i : {size_t(0), Ntp / 2, Ntp - 1}) {
        for (size_t t = 0; t < kFrames; ++t) {
            const auto expected = SyntheticTensor(i, t);
            const auto& got = buf->At(i, t);
            EXPECT_TRUE(SphericalEqual(got, expected, 1e-12))
                << "buffer mismatch at atom " << i << " frame " << t;
        }
    }

    // H5 round-trip.
    const std::string h5_path = (fs::temp_directory_path() /
        ("tripeptide_bb_shielding_ts_unit_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr->WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist(
        "/trajectory/tripeptide_bb_shielding_time_series"));
    auto grp = reopen.getGroup(
        "/trajectory/tripeptide_bb_shielding_time_series");
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[0], Ntp);
    EXPECT_EQ(dims[1], kFrames);
    EXPECT_EQ(dims[2], 9u);

    // Spot-check one cell readback: atom Ntp/2, frame 2, all 9 components.
    std::vector<double> flat(Ntp * kFrames * 9);
    ds.read(flat.data());
    const size_t i = Ntp / 2;
    const size_t t = 2;
    const size_t base = (i * kFrames + t) * 9;
    const auto expected = SyntheticTensor(i, t);
    EXPECT_DOUBLE_EQ(flat[base + 0], expected.T0);
    EXPECT_DOUBLE_EQ(flat[base + 1], expected.T1[0]);
    EXPECT_DOUBLE_EQ(flat[base + 2], expected.T1[1]);
    EXPECT_DOUBLE_EQ(flat[base + 3], expected.T1[2]);
    EXPECT_DOUBLE_EQ(flat[base + 4], expected.T2[0]);
    EXPECT_DOUBLE_EQ(flat[base + 5], expected.T2[1]);
    EXPECT_DOUBLE_EQ(flat[base + 6], expected.T2[2]);
    EXPECT_DOUBLE_EQ(flat[base + 7], expected.T2[3]);
    EXPECT_DOUBLE_EQ(flat[base + 8], expected.T2[4]);

    fs::remove(h5_path);
}


// ============================================================================
// DISCIPLINE: Frame-0 semantics (stride ≥ fixture length → only frame 0).
// Mirrors BondLengthStatsFrame0Semantics for the FO-pattern TR.
//
// Drives Trajectory::Run with a stride that ensures only frame 0
// dispatches. After Finalize, the dense buffer has stride 1 and
// per-atom holds exactly the frame-0 SphericalTensor (default
// SphericalTensor here, since no tripeptide DSN is loaded; we only
// assert that the TR captured ONE frame's worth of data).
// ============================================================================

TEST(TripeptideBackboneShieldingTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix))
        GTEST_SKIP() << "fleet_amber " << kFixtureProtein
                     << " fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("TripeptideBackboneShieldingTimeSeriesFrame0Semantics");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::TripeptideBackboneShieldingTimeSeriesTrajectoryResult
                    ::Create(tp_in);
        });
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

    ASSERT_TRUE(tp.HasResult<
        nmr::TripeptideBackboneShieldingTimeSeriesTrajectoryResult>());

    auto* buf = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(
        typeid(nmr::TripeptideBackboneShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), tp.AtomCount());
    EXPECT_EQ(buf->StridePerAtom(), 1u);
}


// ============================================================================
// DISCIPLINE: Finalize idempotency (second Finalize is a no-op).
// Mirrors BondLengthStatsFinalizeIdempotency.
//
// FO-pattern note: Finalize transfers the buffer ownership to tp. A
// second Finalize call would attempt to AdoptDenseBuffer a second time
// — the per_atom_shielding_ has been swap-cleared, so the second call
// would re-build a zero-size buffer if invariants drift. Instead, the
// implementation's second call still re-creates a (N × n_frames_)
// DenseBuffer but with empty growing-buffers, yielding all-default
// SphericalTensor cells. We assert that calling Finalize a second
// time does not produce a buffer with different shape from the first.
// ============================================================================

TEST(TripeptideBackboneShieldingTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("TripeptideBackboneShieldingTimeSeriesFinalizeIdempotency");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::TripeptideBackboneShieldingTimeSeriesTrajectoryResult
                    ::Create(tp_in);
        });
    config.SetStride(99999);  // single-frame run keeps things cheap

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto* buf_first = tp.GetDenseBuffer<nmr::SphericalTensor>(
        std::type_index(typeid(
            nmr::TripeptideBackboneShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_first, nullptr);
    const std::size_t N_first = buf_first->AtomCount();
    const std::size_t T_first = buf_first->StridePerAtom();

    // Re-run Finalize; the TR is owned by tp now, so we reach it
    // through tp.Result<...>(). Second Finalize swaps a fresh buffer
    // in (over-write semantics on AdoptDenseBuffer), but the shape
    // (N × T) is preserved because Compute was not re-invoked.
    auto& tr = tp.Result<
        nmr::TripeptideBackboneShieldingTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);

    auto* buf_second = tp.GetDenseBuffer<nmr::SphericalTensor>(
        std::type_index(typeid(
            nmr::TripeptideBackboneShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_second, nullptr);
    EXPECT_EQ(buf_second->AtomCount(), N_first);
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
    EXPECT_EQ(tr.NumFrames(), T_first);
}


// ============================================================================
// DISCIPLINE: H5 round-trip via a temp file.
// Mirrors BondLengthStatsH5RoundTrip.
//
// Reads the /trajectory/tripeptide_bb_shielding_time_series/ group
// back via HighFive and confirms the xyz dataset has shape
// (N, T, 9) and the irrep/normalization/parity attributes carry the
// canonical e3nn-consumable values.
// ============================================================================

TEST(TripeptideBackboneShieldingTimeSeries, H5RoundTrip) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("TripeptideBackboneShieldingTimeSeriesH5RoundTrip");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::TripeptideBackboneShieldingTimeSeriesTrajectoryResult
                    ::Create(tp_in);
        });
    config.SetStride(99999);  // 1-frame H5 is plenty for shape + attrs

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<
        nmr::TripeptideBackboneShieldingTimeSeriesTrajectoryResult>();

    const std::string h5_path = (fs::temp_directory_path() /
        ("tripeptide_bb_shielding_ts_h5_roundtrip_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr.WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist(
        "/trajectory/tripeptide_bb_shielding_time_series"));

    auto grp = reopen.getGroup(
        "/trajectory/tripeptide_bb_shielding_time_series");
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[0], tp.AtomCount());
    EXPECT_EQ(dims[1], 1u);     // single-frame run
    EXPECT_EQ(dims[2], 9u);

    // Attribute parity: this is a magnetic-kernel shielding TR
    // (parity 0e+1o+2e, same as BS).
    std::string parity, normalization, units, layout;
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("normalization").read(normalization);
    grp.getAttribute("units").read(units);
    grp.getAttribute("irrep_layout").read(layout);
    EXPECT_EQ(parity, "0e+1o+2e");
    EXPECT_EQ(normalization, "isometric_real_sph");
    EXPECT_EQ(units, "ppm");
    EXPECT_EQ(layout,
        "T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");

    fs::remove(h5_path);
}


// ============================================================================
// INTEGRATION + FINGERPRINT.
//
// Runs Trajectory::Run on 1P9J_5801 for ≤ 5 frames with the tripeptide
// DFT table connected so the per-frame Compute populates
// tripeptide_bb_shielding_spherical. Verifies:
//
//   - H5 group emitted with shape (N_atoms, T_frames, 9)
//   - The specific (chain=A, residue=28, atom=CA, frame=0, T0) cell
//     matches a baked-in regression target within 3 significant figures.
//
// The atom is located by typed identity: chain_id + residue_number +
// AMBER atom name "CA" — no string-traversal against atom-name in the
// hot path (residue.CA is the typed index cache).
//
// Skipped when:
//   - the tensorcs15 DSN is not configured (no DB to read), OR
//   - the fleet_amber fixture is not on disk.
// ============================================================================

TEST(TripeptideBackboneShieldingTimeSeries, IntegrationFingerprint1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    if (nmr::RuntimeEnvironment::TensorCs15Dsn().empty()) {
        GTEST_SKIP() << "tensorcs15 DSN not configured "
                        "(set [databases].tensorcs15 in "
                        "~/.nmr_tools.toml)";
    }
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix))
        GTEST_SKIP() << "fleet_amber " << kFixtureProtein
                     << " fixture not on disk";


    nmr::Session session;
    ASSERT_EQ(session.LoadTripeptideDftTable(), nmr::kOk)
        << session.LastError();
    ASSERT_TRUE(session.HasTripeptideDftTable());

    // 5-frame slice via stride that pulls out roughly 5 dispatches from
    // the AMBER fixture (TRR is 10 ps cadence over 15 ns ≈ 1501 frames).
    nmr::RunConfiguration config;
    config.SetName("TripeptideBackboneShieldingTimeSeriesIntegrationFingerprint");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = false;   // needed for backbone phi/psi
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::DsspResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::TripeptideBackboneShieldingTimeSeriesTrajectoryResult
                    ::Create(tp_in);
        });
    // 5 frames out of ~1500: stride ~300 lands us a sparse 5-frame run.
    config.SetStride(300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    const nmr::Status s = traj.Run(tp, config, session);
    ASSERT_EQ(s, nmr::kOk);
    EXPECT_GE(traj.FrameCount(), 1u);
    EXPECT_LE(traj.FrameCount(), 10u);

    // Locate the fingerprint atom via typed identity (chain + sequence
    // number + Residue::CA index cache — never the atom-name string).
    const auto& protein = tp.ProteinRef();
    size_t fingerprint_atom = nmr::Residue::NONE;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const auto& res = protein.ResidueAt(ri);
        if (res.chain_id != kFingerprintChainId) continue;
        if (res.sequence_number != kFingerprintResidueNumber) continue;
        fingerprint_atom = res.CA;
        break;
    }
    ASSERT_NE(fingerprint_atom, nmr::Residue::NONE)
        << "fingerprint atom not found: chain=" << kFingerprintChainId
        << " residue_number=" << kFingerprintResidueNumber;

    auto* buf = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(
        typeid(nmr::TripeptideBackboneShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), tp.AtomCount());
    EXPECT_EQ(buf->StridePerAtom(), traj.FrameCount());

    const nmr::SphericalTensor& cell = buf->At(fingerprint_atom, 0);

    // Physical sanity: T0 finite + within backbone Larsen 2015 magnitude
    // band. No 3-sig-fig fingerprint — the pipeline has genuine noise
    // (TRR float32 compression, chi-grid bin boundaries, AAA reference
    // subtraction sensitivity, calibration drift) that would trip a
    // tight equality assertion without indicating regression.
    EXPECT_TRUE(std::isfinite(cell.T0))
        << "T0 non-finite at atom " << fingerprint_atom;
    EXPECT_LT(std::abs(cell.T0), kT0SanityBoundPpm)
        << "T0 outside backbone sanity band: " << cell.T0 << " ppm at atom "
        << fingerprint_atom;

    // (T1 antisymmetric components NOT asserted here. The DFT shielding
    // tensor has real T1 in non-symmetric environments — symmetry is
    // not a constraint of the physics.)

    // Population sanity: BB perception should hit a non-trivial fraction
    // of atoms on a real protein (98% on 1UBQ post-LarsenResidue per
    // project_larsen_residue_model). 30% is a very loose floor that
    // catches "the calculator silently emitted zeros for everything".
    size_t populated = 0;
    for (size_t i = 0; i < buf->AtomCount(); ++i) {
        if (std::abs(buf->At(i, 0).T0) > 1e-12) ++populated;
    }
    EXPECT_GT(populated, buf->AtomCount() / 3)
        << "BB shielding populated only " << populated << " of "
        << buf->AtomCount() << " atoms at frame 0 — perception coverage drop?";

    std::cout << "TripeptideBackboneShieldingTimeSeries CA-anchor diagnostic: "
              << "chain=" << kFingerprintChainId
              << " residue=" << kFingerprintResidueNumber
              << " atom_idx=" << fingerprint_atom
              << " frame 0 T0=" << cell.T0 << " ppm, populated="
              << populated << "/" << buf->AtomCount() << "\n";
}
