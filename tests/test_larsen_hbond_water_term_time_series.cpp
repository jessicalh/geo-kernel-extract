//
// test_larsen_hbond_water_term_time_series: discipline + integration
// tests for LarsenHBondWaterTermTimeSeriesTrajectoryResult.
//
// Establishes the scalar-double (N, T) 2D H5 emission pattern for the
// Larsen bundle. Five tests: SyntheticFourFrames + Frame0Semantics +
// FinalizeIdempotency + H5RoundTrip + IntegrationDistribution1P9J.
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
#include "LarsenHBondShieldingResult.h"
#include "LarsenHBondWaterTermTimeSeriesTrajectoryResult.h"
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

// Larsen Δσ_w canonical value (2.07 ppm); per-cell value is either
// 0.0 (non-HN or HN with pairs) or exactly 2.07 (HN with no pairs).
constexpr double kLarsenWaterTermPpm = 2.07;


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


// Synthetic scalar value per (atom, frame). Distinct so component-
// mapping errors surface as wrong values.
double SyntheticWaterTerm(std::size_t atom_i, std::size_t frame_t) {
    return static_cast<double>(atom_i)
         + static_cast<double>(frame_t) * 100.0
         + 0.5;
}


}  // namespace


TEST(LarsenHBondWaterTermTimeSeries, SyntheticFourFrames) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();

    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    const std::size_t Ntp = tp.AtomCount();
    ASSERT_GT(Ntp, 0u);

    auto tr = nmr::LarsenHBondWaterTermTimeSeriesTrajectoryResult::Create(tp);
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
            conf->MutableAtomAt(i).larsen_hbond_water_term =
                SyntheticWaterTerm(i, t);
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    EXPECT_EQ(tr->NumFrames(), kFrames);

    tr->Finalize(tp, traj);

    auto* buf = tp.GetDenseBuffer<double>(std::type_index(typeid(
        nmr::LarsenHBondWaterTermTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), Ntp);
    EXPECT_EQ(buf->StridePerAtom(), kFrames);

    for (std::size_t i : {std::size_t(0), Ntp / 2, Ntp - 1}) {
        for (std::size_t t = 0; t < kFrames; ++t) {
            EXPECT_DOUBLE_EQ(buf->At(i, t), SyntheticWaterTerm(i, t))
                << "buffer mismatch at atom " << i << " frame " << t;
        }
    }

    const std::string h5_path = (fs::temp_directory_path() /
        ("larsen_hbond_water_term_ts_unit_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr->WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist(
        "/trajectory/larsen_hbond_water_term_time_series"));
    auto grp = reopen.getGroup(
        "/trajectory/larsen_hbond_water_term_time_series");
    auto ds = grp.getDataSet("water_term");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 2u);
    EXPECT_EQ(dims[0], Ntp);
    EXPECT_EQ(dims[1], kFrames);

    std::vector<double> flat(Ntp * kFrames);
    ds.read(flat.data());
    const std::size_t i = Ntp / 2;
    const std::size_t t = 2;
    EXPECT_DOUBLE_EQ(flat[i * kFrames + t], SyntheticWaterTerm(i, t));

    fs::remove(h5_path);
}


TEST(LarsenHBondWaterTermTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("LarsenHBondWaterTermTimeSeriesFrame0Semantics");
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
            return nmr::LarsenHBondWaterTermTimeSeriesTrajectoryResult
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
        nmr::LarsenHBondWaterTermTimeSeriesTrajectoryResult>());

    auto* buf = tp.GetDenseBuffer<double>(std::type_index(
        typeid(nmr::LarsenHBondWaterTermTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), tp.AtomCount());
    EXPECT_EQ(buf->StridePerAtom(), 1u);
}


TEST(LarsenHBondWaterTermTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("LarsenHBondWaterTermTimeSeriesFinalizeIdempotency");
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
            return nmr::LarsenHBondWaterTermTimeSeriesTrajectoryResult
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

    auto* buf_first = tp.GetDenseBuffer<double>(
        std::type_index(typeid(
            nmr::LarsenHBondWaterTermTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_first, nullptr);
    const std::size_t N_first = buf_first->AtomCount();
    const std::size_t T_first = buf_first->StridePerAtom();

    auto& tr = tp.Result<
        nmr::LarsenHBondWaterTermTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);

    auto* buf_second = tp.GetDenseBuffer<double>(
        std::type_index(typeid(
            nmr::LarsenHBondWaterTermTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_second, nullptr);
    EXPECT_EQ(buf_second->AtomCount(), N_first);
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
    EXPECT_EQ(tr.NumFrames(), T_first);
}


TEST(LarsenHBondWaterTermTimeSeries, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    nmr::Session session;
    if (session.LoadLarsenHBondGrid() != nmr::kOk
        || !session.HasLarsenHBondGrid()) {
        GTEST_SKIP() << "larsen_hbond_grids not configured — source calc "
                        "won't attach so the TR correctly skips H5 emission, "
                        "and round-trip cannot be verified";
    }
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("LarsenHBondWaterTermTimeSeriesH5RoundTrip");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = false;  // backbone phi/psi needed for the larsen H-bond calc to actually run
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::DsspResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::LarsenHBondWaterTermTimeSeriesTrajectoryResult
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
        nmr::LarsenHBondWaterTermTimeSeriesTrajectoryResult>();

    const std::string h5_path = (fs::temp_directory_path() /
        ("larsen_hbond_water_term_ts_h5_roundtrip_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr.WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist(
        "/trajectory/larsen_hbond_water_term_time_series"));

    auto grp = reopen.getGroup(
        "/trajectory/larsen_hbond_water_term_time_series");
    auto ds = grp.getDataSet("water_term");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 2u);
    EXPECT_EQ(dims[0], tp.AtomCount());
    EXPECT_EQ(dims[1], 1u);

    std::string irrep_layout, parity, units;
    grp.getAttribute("irrep_layout").read(irrep_layout);
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("units").read(units);
    EXPECT_EQ(irrep_layout, "T0");
    EXPECT_EQ(parity, "0e");
    EXPECT_EQ(units, "ppm");

    fs::remove(h5_path);
}


// ============================================================================
// INTEGRATION + DISTRIBUTION DIAGNOSTIC.
//
// Runs Trajectory::Run on 1P9J_5801 with Larsen H-bond grids loaded.
// The water_term should be EXACTLY 0.0 or EXACTLY 2.07 ppm on every
// cell (Larsen Δσ_w is a fixed NMA-water complex constant; the
// calculator applies it as 0 or 2.07, no interpolation). We:
//   - assert no NaN, no Inf
//   - assert every cell is in {0.0, 2.07} (no off-quantum values)
//   - tally count(0.0) and count(2.07); the 2.07 fraction should be
//     small but nonzero (amide H atoms with no H-bond pair this frame)
//   - log the distribution
// ============================================================================

TEST(LarsenHBondWaterTermTimeSeries, IntegrationDistribution1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();

    nmr::Session session;
    if (session.LoadLarsenHBondGrid() != nmr::kOk
        || !session.HasLarsenHBondGrid()) {
        GTEST_SKIP() << "larsen_hbond_grids not configured in "
                        "~/.nmr_tools.toml — Larsen calc cannot run, "
                        "integration assertion would be vacuous";
    }
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("LarsenHBondWaterTermTimeSeriesIntegration");
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
            return nmr::LarsenHBondWaterTermTimeSeriesTrajectoryResult
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

    auto* buf = tp.GetDenseBuffer<double>(std::type_index(
        typeid(nmr::LarsenHBondWaterTermTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);

    const std::size_t N = buf->AtomCount();
    const std::size_t T = buf->StridePerAtom();
    std::size_t zero_count = 0;
    std::size_t larsen_count = 0;
    std::size_t other_count = 0;
    double max_other = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const double v = buf->At(i, t);
            ASSERT_TRUE(std::isfinite(v))
                << "non-finite water_term at atom " << i << " frame " << t;
            if (v == 0.0) {
                ++zero_count;
            } else if (std::abs(v - kLarsenWaterTermPpm) < 1e-9) {
                ++larsen_count;
            } else {
                ++other_count;
                if (std::abs(v) > std::abs(max_other)) max_other = v;
            }
        }
    }

    EXPECT_EQ(other_count, 0u)
        << "found " << other_count << " water_term cells outside "
        << "{0, " << kLarsenWaterTermPpm << "} (max=" << max_other << ")";

    // Real-data coverage floor: the grid is loaded, source calc must
    // have evaluated something. zero_count + larsen_count covers all
    // {0.0, 2.07} cells (other_count is asserted == 0 just above). With
    // N atoms × T frames cells, the calc should have written into at
    // least N×T/10 of them — anything less suggests the source calc
    // returned all-zero (grid loaded but never queried), which is the
    // exact regression the source-attached gate was added to catch.
    const std::size_t covered = zero_count + larsen_count;
    const std::size_t coverage_floor = (N * T) / 10;
    EXPECT_GT(covered, coverage_floor)
        << "Larsen water term covered only " << covered << " of " << (N * T)
        << " cells (floor=" << coverage_floor << ") — calc likely never ran "
        << "despite grid being loaded; possible source-attached regression";

    // 2.07-population floor: 1P9J has ~76 amide H atoms; in any given
    // frame some fraction lack H-bond pairs and pick up Δσ_w. Across T
    // frames we expect at least a handful. This catches the
    // "grid loaded but every cell stayed at default 0.0" regression
    // (vs. "calc ran and legitimately produced 0.0 for paired HN").
    EXPECT_GT(larsen_count, 0u)
        << "no Larsen Δσ_w (2.07 ppm) cells across " << T
        << " frames — calc likely returned all zeros";

    std::cout << "LarsenHBondWaterTerm distribution: "
              << "frames=" << T
              << " atoms=" << N
              << " zero=" << zero_count
              << " larsen(2.07)=" << larsen_count
              << " other=" << other_count << "\n";
}
