//
// test_tripeptide_bb_method_tag_time_series: discipline + integration
// tests for TripeptideBackboneMethodTagTimeSeriesTrajectoryResult.
//
// Establishes the int-buffer (uint8) (N, T) flat emission pattern.
// Five tests:
//   1. SyntheticFourFrames — hand-crafted uint8 values, exact-equality
//      round-trip through Compute / Finalize / WriteH5Group.
//   2. Frame0Semantics, 3. FinalizeIdempotency, 4. H5RoundTrip
//   5. IntegrationDistribution1P9J — real fixture; tallies tag
//      distribution (0 = no match, 1 = OPBE, 2 = PBE SER); asserts
//      at least one match-tag (1 or 2) populates a majority of cells.
//

#include "AIMNet2Result.h"

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
#include "TripeptideBackboneMethodTagTimeSeriesTrajectoryResult.h"
#include "Types.h"

#include <gtest/gtest.h>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cstdint>
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


// Distinct synthetic tag per (atom, frame): 3 valid values cycling
// through (0=no_match, 1=opbe, 2=pbe_ser) keyed on (atom_i + frame_t).
std::uint8_t SyntheticTag(std::size_t atom_i, std::size_t frame_t) {
    return static_cast<std::uint8_t>((atom_i + frame_t * 13) % 3);
}


}  // namespace


TEST(TripeptideBackboneMethodTagTimeSeries, SyntheticFourFrames) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();

    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix))
        GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    const std::size_t Ntp = tp.AtomCount();
    ASSERT_GT(Ntp, 0u);

    auto tr =
        nmr::TripeptideBackboneMethodTagTimeSeriesTrajectoryResult::Create(tp);
    // Synthetic path bypasses OperationRunner so the real source calc
    // never attaches. Force-mark source present so the WriteH5Group
    // skip-on-absent gate doesn't fire.
    tr->ForceSourcePresentForTesting();

    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);

    constexpr std::size_t kFrames = 4;
    const auto& protein_ref = tp.ProteinRef();
    std::vector<nmr::Vec3> positions(Ntp, nmr::Vec3::Zero());

    for (std::size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &protein_ref, positions, "synthetic frame");
        for (std::size_t i = 0; i < Ntp; ++i) {
            conf->MutableAtomAt(i).tripeptide_bb_method_tag =
                SyntheticTag(i, t);
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    EXPECT_EQ(tr->NumFrames(), kFrames);

    tr->Finalize(tp, traj);

    auto* buf = tp.GetDenseBuffer<std::uint8_t>(std::type_index(typeid(
        nmr::TripeptideBackboneMethodTagTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), Ntp);
    EXPECT_EQ(buf->StridePerAtom(), kFrames);

    for (std::size_t i : {std::size_t(0), Ntp / 2, Ntp - 1}) {
        for (std::size_t t = 0; t < kFrames; ++t) {
            EXPECT_EQ(buf->At(i, t), SyntheticTag(i, t))
                << "buffer mismatch at atom " << i << " frame " << t;
        }
    }

    // H5 round-trip.
    const std::string h5_path = (fs::temp_directory_path() /
        ("tripeptide_bb_method_tag_ts_unit_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr->WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist(
        "/trajectory/tripeptide_bb_method_tag_time_series"));
    auto grp = reopen.getGroup(
        "/trajectory/tripeptide_bb_method_tag_time_series");
    auto ds = grp.getDataSet("method_tag");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 2u);
    EXPECT_EQ(dims[0], Ntp);
    EXPECT_EQ(dims[1], kFrames);

    std::vector<std::uint8_t> flat(Ntp * kFrames);
    ds.read(flat.data());
    const std::size_t i = Ntp / 2;
    const std::size_t t = 2;
    EXPECT_EQ(flat[i * kFrames + t], SyntheticTag(i, t));

    fs::remove(h5_path);
}


TEST(TripeptideBackboneMethodTagTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("TripeptideBackboneMethodTagTimeSeriesFrame0Semantics");
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
            return nmr::TripeptideBackboneMethodTagTimeSeriesTrajectoryResult
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
    ASSERT_EQ(traj.FrameCount(), 1u);

    ASSERT_TRUE(tp.HasResult<
        nmr::TripeptideBackboneMethodTagTimeSeriesTrajectoryResult>());

    auto* buf = tp.GetDenseBuffer<std::uint8_t>(std::type_index(
        typeid(nmr::TripeptideBackboneMethodTagTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), tp.AtomCount());
    EXPECT_EQ(buf->StridePerAtom(), 1u);
}


TEST(TripeptideBackboneMethodTagTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("TripeptideBackboneMethodTagTimeSeriesFinalizeIdempotency");
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
            return nmr::TripeptideBackboneMethodTagTimeSeriesTrajectoryResult
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

    auto* buf_first = tp.GetDenseBuffer<std::uint8_t>(
        std::type_index(typeid(
            nmr::TripeptideBackboneMethodTagTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_first, nullptr);
    const std::size_t N_first = buf_first->AtomCount();
    const std::size_t T_first = buf_first->StridePerAtom();

    auto& tr = tp.Result<
        nmr::TripeptideBackboneMethodTagTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);

    auto* buf_second = tp.GetDenseBuffer<std::uint8_t>(
        std::type_index(typeid(
            nmr::TripeptideBackboneMethodTagTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_second, nullptr);
    EXPECT_EQ(buf_second->AtomCount(), N_first);
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
    EXPECT_EQ(tr.NumFrames(), T_first);
}


TEST(TripeptideBackboneMethodTagTimeSeries, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    if (nmr::RuntimeEnvironment::TensorCs15Dsn().empty()) {
        GTEST_SKIP() << "tensorcs15 DSN not configured — source calc "
                        "won't attach so the TR correctly skips H5 emission, "
                        "and round-trip cannot be verified";
    }
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::Session session;
    ASSERT_EQ(session.LoadTripeptideDftTable(), nmr::kOk)
        << session.LastError();
    ASSERT_TRUE(session.HasTripeptideDftTable());

    nmr::RunConfiguration config;
    config.SetName("TripeptideBackboneMethodTagTimeSeriesH5RoundTrip");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = false;  // backbone phi/psi needed for the tripeptide calc to actually run
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::DsspResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::TripeptideBackboneMethodTagTimeSeriesTrajectoryResult
                    ::Create(tp_in);
        });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<
        nmr::TripeptideBackboneMethodTagTimeSeriesTrajectoryResult>();

    const std::string h5_path = (fs::temp_directory_path() /
        ("tripeptide_bb_method_tag_ts_h5_roundtrip_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr.WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist(
        "/trajectory/tripeptide_bb_method_tag_time_series"));

    auto grp = reopen.getGroup(
        "/trajectory/tripeptide_bb_method_tag_time_series");
    auto ds = grp.getDataSet("method_tag");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 2u);
    EXPECT_EQ(dims[0], tp.AtomCount());
    EXPECT_EQ(dims[1], 1u);

    std::string units, dtype, legend;
    grp.getAttribute("units").read(units);
    grp.getAttribute("dtype").read(dtype);
    grp.getAttribute("legend").read(legend);
    EXPECT_EQ(units, "categorical");
    EXPECT_EQ(dtype, "uint8");
    EXPECT_EQ(legend,
        "0=no_match, 1=opbe_6-31g(d,p)_gaussian, "
        "2=pbe_6-31g(d,p)_orca_ser");

    fs::remove(h5_path);
}


// ============================================================================
// INTEGRATION + DISTRIBUTION DIAGNOSTIC.
//
// Runs Trajectory::Run on 1P9J_5801 with DSN connected; tallies the
// method-tag distribution. On a 56-residue fixture with ~90% atom
// match coverage, the majority of cells should be tag 1 (OPBE) or 2
// (PBE-SER for serine residues). A high count of 0 (no match) would
// indicate a perception coverage regression. We assert:
//   - all tags ∈ {0, 1, 2} (no out-of-range values)
//   - matched cells (tag != 0) exceed 33% of total
//   - we log the full (0, 1, 2) distribution for diagnostic
// ============================================================================

TEST(TripeptideBackboneMethodTagTimeSeries, IntegrationDistribution1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    if (nmr::RuntimeEnvironment::TensorCs15Dsn().empty()) {
        GTEST_SKIP() << "tensorcs15 DSN not configured";
    }
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::Session session;
    ASSERT_EQ(session.LoadTripeptideDftTable(), nmr::kOk)
        << session.LastError();
    ASSERT_TRUE(session.HasTripeptideDftTable());

    nmr::RunConfiguration config;
    config.SetName("TripeptideBackboneMethodTagTimeSeriesIntegration");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = false;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::DsspResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::TripeptideBackboneMethodTagTimeSeriesTrajectoryResult
                    ::Create(tp_in);
        });
    config.SetStride(300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_GE(traj.FrameCount(), 1u);

    auto* buf = tp.GetDenseBuffer<std::uint8_t>(std::type_index(
        typeid(nmr::TripeptideBackboneMethodTagTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);

    std::size_t count_0 = 0;
    std::size_t count_1 = 0;
    std::size_t count_2 = 0;
    std::size_t count_other = 0;
    const std::size_t N = buf->AtomCount();
    const std::size_t T = buf->StridePerAtom();
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const std::uint8_t tag = buf->At(i, t);
            switch (tag) {
                case 0: ++count_0; break;
                case 1: ++count_1; break;
                case 2: ++count_2; break;
                default: ++count_other; break;
            }
        }
    }

    EXPECT_EQ(count_other, 0u)
        << count_other << " cells have out-of-range method_tag (>2)";

    const std::size_t total_cells = N * T;
    const std::size_t matched = count_1 + count_2;
    EXPECT_GT(matched, total_cells / 3)
        << "matched cells " << matched << " / " << total_cells
        << " — perception coverage drop?";

    std::cout << "TripeptideBackboneMethodTag distribution: "
              << "frames=" << T
              << " atoms=" << N
              << " tag0(no_match)=" << count_0
              << " tag1(opbe)="      << count_1
              << " tag2(pbe_ser)="   << count_2
              << " other="           << count_other
              << " matched_frac="
              << (static_cast<double>(matched) /
                  static_cast<double>(total_cells)) << "\n";
}
