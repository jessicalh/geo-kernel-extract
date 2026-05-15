//
// test_tripeptide_bb_residual_vec_time_series: discipline + integration
// tests for TripeptideBackboneResidualVecTimeSeriesTrajectoryResult.
//
// Four-test pattern mirrors test_tripeptide_bb_shielding_time_series.cpp
// against the Vec3 residual buffer:
//
//   1. SyntheticFourFrames — hand-crafted Vec3 values, exact-equality
//      round-trip through Compute / Finalize / WriteH5Group. No DFT
//      table needed.
//   2. Frame0Semantics — stride > fixture length asserts only frame 0
//      dispatches and the dense buffer arrives with stride 1.
//   3. FinalizeIdempotency — second Finalize is a bounds-check no-op,
//      buffer shape preserved.
//   4. H5RoundTrip — single-frame Trajectory::Run + HighFive read-back
//      confirms (N, T, 3) shape and the cartesian/1o/angstrom attrs.
//   5. IntegrationLogOverages1P9J — real fleet_amber fixture with DSN
//      connected. Per feedback_log_overages_dont_assert: assert
//      isfinite + population floor; LOG magnitude overages without
//      asserting (chi-grid coarseness routinely produces 1-4 Å Vec3
//      norms on deep sidechains — legitimate, not regression).
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
#include "TripeptideBackboneResidualVecTimeSeriesTrajectoryResult.h"
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

// Match the BB shielding test's central-residue anchor for diagnostics.
constexpr int kFingerprintResidueNumber = 28;
const std::string kFingerprintChainId = "A";

// Physical sanity threshold for LOGGING (NOT asserted). Backbone atoms
// typically ≤ 0.5 Å residual; deep sidechain rotamer mismatches reach
// 1-4 Å routinely. > 5 Å is unusual but legitimate (chi-fallback,
// perception edges, AAA reference subtraction sensitivity). Per
// feedback_log_overages_dont_assert: the pipeline produces imperfect
// output for legitimate reasons and calibration absorbs; log the
// overage, don't fail the test.
constexpr double kResidualMagWarnAngstrom = 5.0;


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


// Build a 4-frame synthetic Vec3 sequence per atom. Returns the value
// the (atom_i, frame_t) cell should carry. x = atom_i + t*100;
// y = atom_i + t*100 + 0.1; z = atom_i + t*100 + 0.2. Distinct per
// component so component-mapping errors surface as wrong values, not
// equal-looking values.
nmr::Vec3 SyntheticVec3(size_t atom_i, size_t frame_t) {
    const double base = static_cast<double>(atom_i)
                      + static_cast<double>(frame_t) * 100.0;
    return nmr::Vec3(base, base + 0.1, base + 0.2);
}


bool Vec3Equal(const nmr::Vec3& a, const nmr::Vec3& b, double tol) {
    return std::abs(a.x() - b.x()) <= tol
        && std::abs(a.y() - b.y()) <= tol
        && std::abs(a.z() - b.z()) <= tol;
}


}  // namespace


// ============================================================================
// UNIT: 4-frame synthetic round-trip on hand-crafted Vec3 values.
//
// Drives Compute / Finalize / WriteH5Group with hand-crafted per-atom
// inputs and asserts exact bit-equality on round-trip. Skips
// Trajectory::Run orchestration; uses the fleet_amber fixture only
// for TrajectoryProtein AtomCount plumbing. Exact equality is
// appropriate here because synthetic inputs have no numerical noise.
// ============================================================================

TEST(TripeptideBackboneResidualVecTimeSeries, SyntheticFourFrames) {
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

    auto tr =
        nmr::TripeptideBackboneResidualVecTimeSeriesTrajectoryResult::Create(tp);

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
            conf->MutableAtomAt(i).tripeptide_bb_residual_vec =
                SyntheticVec3(i, t);
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    EXPECT_EQ(tr->NumFrames(), kFrames);

    tr->Finalize(tp, traj);

    auto* buf = tp.GetDenseBuffer<nmr::Vec3>(std::type_index(typeid(
        nmr::TripeptideBackboneResidualVecTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), Ntp);
    EXPECT_EQ(buf->StridePerAtom(), kFrames);

    for (size_t i : {size_t(0), Ntp / 2, Ntp - 1}) {
        for (size_t t = 0; t < kFrames; ++t) {
            const auto expected = SyntheticVec3(i, t);
            const auto& got = buf->At(i, t);
            EXPECT_TRUE(Vec3Equal(got, expected, 1e-12))
                << "buffer mismatch at atom " << i << " frame " << t;
        }
    }

    // H5 round-trip.
    const std::string h5_path = (fs::temp_directory_path() /
        ("tripeptide_bb_residual_vec_ts_unit_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr->WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist(
        "/trajectory/tripeptide_bb_residual_vec_time_series"));
    auto grp = reopen.getGroup(
        "/trajectory/tripeptide_bb_residual_vec_time_series");
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[0], Ntp);
    EXPECT_EQ(dims[1], kFrames);
    EXPECT_EQ(dims[2], 3u);

    // Spot-check one cell readback: atom Ntp/2, frame 2, all 3 components.
    std::vector<double> flat(Ntp * kFrames * 3);
    ds.read(flat.data());
    const size_t i = Ntp / 2;
    const size_t t = 2;
    const size_t base = (i * kFrames + t) * 3;
    const auto expected = SyntheticVec3(i, t);
    EXPECT_DOUBLE_EQ(flat[base + 0], expected.x());
    EXPECT_DOUBLE_EQ(flat[base + 1], expected.y());
    EXPECT_DOUBLE_EQ(flat[base + 2], expected.z());

    fs::remove(h5_path);
}


// ============================================================================
// DISCIPLINE: Frame-0 semantics (stride ≥ fixture length → only frame 0).
// ============================================================================

TEST(TripeptideBackboneResidualVecTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix))
        GTEST_SKIP() << "fleet_amber " << kFixtureProtein
                     << " fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("TripeptideBackboneResidualVecTimeSeriesFrame0Semantics");
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
            return nmr::TripeptideBackboneResidualVecTimeSeriesTrajectoryResult
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
        nmr::TripeptideBackboneResidualVecTimeSeriesTrajectoryResult>());

    auto* buf = tp.GetDenseBuffer<nmr::Vec3>(std::type_index(
        typeid(nmr::TripeptideBackboneResidualVecTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), tp.AtomCount());
    EXPECT_EQ(buf->StridePerAtom(), 1u);
}


// ============================================================================
// DISCIPLINE: Finalize idempotency (second Finalize is a no-op).
//
// The bounds-check idempotency pattern (feedback_bounds_check_over_state_flag)
// makes the second call short-circuit when per_atom_residual_ is post-swap
// empty. Shape (N × T) is preserved because Compute was not re-invoked.
// ============================================================================

TEST(TripeptideBackboneResidualVecTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("TripeptideBackboneResidualVecTimeSeriesFinalizeIdempotency");
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
            return nmr::TripeptideBackboneResidualVecTimeSeriesTrajectoryResult
                    ::Create(tp_in);
        });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto* buf_first = tp.GetDenseBuffer<nmr::Vec3>(
        std::type_index(typeid(
            nmr::TripeptideBackboneResidualVecTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_first, nullptr);
    const std::size_t N_first = buf_first->AtomCount();
    const std::size_t T_first = buf_first->StridePerAtom();

    auto& tr = tp.Result<
        nmr::TripeptideBackboneResidualVecTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);

    auto* buf_second = tp.GetDenseBuffer<nmr::Vec3>(
        std::type_index(typeid(
            nmr::TripeptideBackboneResidualVecTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_second, nullptr);
    EXPECT_EQ(buf_second->AtomCount(), N_first);
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
    EXPECT_EQ(tr.NumFrames(), T_first);
}


// ============================================================================
// DISCIPLINE: H5 round-trip via a temp file.
// Verifies the cartesian/1o/angstrom attribute triple.
// ============================================================================

TEST(TripeptideBackboneResidualVecTimeSeries, H5RoundTrip) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("TripeptideBackboneResidualVecTimeSeriesH5RoundTrip");
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
            return nmr::TripeptideBackboneResidualVecTimeSeriesTrajectoryResult
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
        nmr::TripeptideBackboneResidualVecTimeSeriesTrajectoryResult>();

    const std::string h5_path = (fs::temp_directory_path() /
        ("tripeptide_bb_residual_vec_ts_h5_roundtrip_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr.WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist(
        "/trajectory/tripeptide_bb_residual_vec_time_series"));

    auto grp = reopen.getGroup(
        "/trajectory/tripeptide_bb_residual_vec_time_series");
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[0], tp.AtomCount());
    EXPECT_EQ(dims[1], 1u);
    EXPECT_EQ(dims[2], 3u);

    // Attribute triple: cartesian xyz layout, polar vector (1o), in Å.
    std::string parity, normalization, units, layout;
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("normalization").read(normalization);
    grp.getAttribute("units").read(units);
    grp.getAttribute("irrep_layout").read(layout);
    EXPECT_EQ(parity, "1o");
    EXPECT_EQ(normalization, "cartesian");
    EXPECT_EQ(units, "angstrom");
    EXPECT_EQ(layout, "x,y,z");

    fs::remove(h5_path);
}


// ============================================================================
// INTEGRATION + LOG-OVERAGES.
//
// Runs Trajectory::Run on 1P9J_5801 with tripeptide DFT table connected so
// the per-frame Compute populates tripeptide_bb_residual_vec. Verifies:
//
//   - H5 buffer shape (N_atoms, T_frames, 3) consistent
//   - Each component isfinite
//   - LOG warnings for |residual_vec| > kResidualMagWarnAngstrom — chi-grid
//     coarseness on deep sidechains routinely produces 1-4 Å Vec3 norms,
//     and >5 Å is legitimate (chi-fallback, perception edges); calibration
//     absorbs (feedback_log_overages_dont_assert)
//   - Population floor (33% of atoms touched at frame 0; matches BB
//     Shielding TR integration sanity)
//
// Skipped when DSN not configured or fixture missing.
// ============================================================================

TEST(TripeptideBackboneResidualVecTimeSeries, IntegrationLogOverages1P9J) {
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

    // 5-frame slice via stride 300 over the ~1500-frame AMBER fixture.
    nmr::RunConfiguration config;
    config.SetName("TripeptideBackboneResidualVecTimeSeriesIntegration");
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
            return nmr::TripeptideBackboneResidualVecTimeSeriesTrajectoryResult
                    ::Create(tp_in);
        });
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

    auto* buf = tp.GetDenseBuffer<nmr::Vec3>(std::type_index(
        typeid(nmr::TripeptideBackboneResidualVecTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), tp.AtomCount());
    EXPECT_EQ(buf->StridePerAtom(), traj.FrameCount());

    // Walk every cell: assert finiteness; tally overages and frame-0
    // populated count for diagnostic emission.
    std::size_t populated_at_frame0 = 0;
    std::size_t overage_count = 0;
    double max_seen_magnitude = 0.0;
    for (std::size_t i = 0; i < buf->AtomCount(); ++i) {
        for (std::size_t t = 0; t < buf->StridePerAtom(); ++t) {
            const nmr::Vec3& v = buf->At(i, t);
            EXPECT_TRUE(std::isfinite(v.x())
                     && std::isfinite(v.y())
                     && std::isfinite(v.z()))
                << "non-finite residual at atom " << i << " frame " << t;
            const double mag = v.norm();
            if (mag > max_seen_magnitude) max_seen_magnitude = mag;
            if (mag > kResidualMagWarnAngstrom) ++overage_count;
            if (t == 0 && mag > 1e-12) ++populated_at_frame0;
        }
    }

    // Overage diagnostic — chi-grid coarseness on deep sidechains
    // routinely pushes Vec3 norm to 1-4 Å; > 5 Å is unusual but
    // legitimate (deep pathology, AAA reference subtraction sensitivity).
    // Logged, never asserted.
    if (overage_count > 0) {
        std::cout << "[warn] TripeptideBackboneResidualVec: " << overage_count
                  << " cells exceed " << kResidualMagWarnAngstrom
                  << " Å (max norm " << max_seen_magnitude
                  << " Å) — chi-grid / perception edges, not asserted\n";
    }

    // Population floor — same 33% threshold as the BB shielding test.
    // Catches the regression "calculator silently emitted zeros for
    // everything" without committing to a specific fingerprint.
    EXPECT_GT(populated_at_frame0, buf->AtomCount() / 3)
        << "residual_vec populated only " << populated_at_frame0
        << " of " << buf->AtomCount() << " atoms at frame 0 — "
        << "perception coverage drop?";

    // Locate the central-residue CA atom for a diagnostic readout
    // (symmetric with the BB shielding test's fingerprint anchor; we
    // print but do NOT assert specific values).
    const auto& protein = tp.ProteinRef();
    size_t anchor_atom = nmr::Residue::NONE;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const auto& res = protein.ResidueAt(ri);
        if (res.chain_id != kFingerprintChainId) continue;
        if (res.sequence_number != kFingerprintResidueNumber) continue;
        anchor_atom = res.CA;
        break;
    }
    if (anchor_atom != nmr::Residue::NONE) {
        const nmr::Vec3& v0 = buf->At(anchor_atom, 0);
        std::cout << "TripeptideBackboneResidualVecTimeSeries diagnostic: "
                  << "chain=" << kFingerprintChainId
                  << " residue=" << kFingerprintResidueNumber
                  << " atom_idx=" << anchor_atom
                  << " frame 0 residual=(" << v0.x() << ", " << v0.y() << ", "
                  << v0.z() << ") |v|=" << v0.norm() << " Å, populated="
                  << populated_at_frame0 << "/" << buf->AtomCount()
                  << " overages=" << overage_count
                  << " max=" << max_seen_magnitude << " Å\n";
    }
}
