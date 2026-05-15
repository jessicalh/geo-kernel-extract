//
// test_larsen_hbond_count_time_series: discipline + integration tests
// for LarsenHBondCountTimeSeriesTrajectoryResult.
//
// Five tests: SyntheticFourFrames + Frame0Semantics + FinalizeIdempotency
// + H5RoundTrip + IntegrationDistribution1P9J.
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
#include "LarsenHBondCountTimeSeriesTrajectoryResult.h"
#include "Types.h"

#include <gtest/gtest.h>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

namespace {

constexpr const char* kFixtureProtein = "1P9J_5801";

// Upper bound for physical pair counts — log-only, not asserted.
// Larsen Table 2 contribution classes mean an atom can in principle
// see contributions from multiple flanking residues; we set a generous
// 16 ceiling for the "log overage" warning.
constexpr int kHBondCountSanityCap = 16;


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


int SyntheticCount(std::size_t atom_i, std::size_t frame_t) {
    return static_cast<int>((atom_i + frame_t * 7) % 5);
}


}  // namespace


TEST(LarsenHBondCountTimeSeries, SyntheticFourFrames) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    const std::size_t Ntp = tp.AtomCount();
    ASSERT_GT(Ntp, 0u);

    auto tr = nmr::LarsenHBondCountTimeSeriesTrajectoryResult::Create(tp);
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
            conf->MutableAtomAt(i).larsen_hbond_n_pairs =
                SyntheticCount(i, t);
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    EXPECT_EQ(tr->NumFrames(), kFrames);

    tr->Finalize(tp, traj);

    auto* buf = tp.GetDenseBuffer<int>(std::type_index(typeid(
        nmr::LarsenHBondCountTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), Ntp);
    EXPECT_EQ(buf->StridePerAtom(), kFrames);

    for (std::size_t i : {std::size_t(0), Ntp / 2, Ntp - 1}) {
        for (std::size_t t = 0; t < kFrames; ++t) {
            EXPECT_EQ(buf->At(i, t), SyntheticCount(i, t));
        }
    }

    const std::string h5_path = (fs::temp_directory_path() /
        ("larsen_hbond_count_ts_unit_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr->WriteH5Group(tp, file);
    }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/larsen_hbond_count_time_series");
    auto ds = grp.getDataSet("count");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 2u);
    EXPECT_EQ(dims[0], Ntp);
    EXPECT_EQ(dims[1], kFrames);

    std::vector<int> flat(Ntp * kFrames);
    ds.read(flat.data());
    const std::size_t i = Ntp / 2;
    const std::size_t t = 2;
    EXPECT_EQ(flat[i * kFrames + t], SyntheticCount(i, t));

    fs::remove(h5_path);
}


TEST(LarsenHBondCountTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("LarsenHBondCountTimeSeriesFrame0Semantics");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true;
    opts.skip_apbs = true;  opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::LarsenHBondCountTimeSeriesTrajectoryResult
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

    auto* buf = tp.GetDenseBuffer<int>(std::type_index(
        typeid(nmr::LarsenHBondCountTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), tp.AtomCount());
    EXPECT_EQ(buf->StridePerAtom(), 1u);
}


TEST(LarsenHBondCountTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("LarsenHBondCountTimeSeriesFinalizeIdempotency");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true;
    opts.skip_apbs = true;  opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::LarsenHBondCountTimeSeriesTrajectoryResult
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

    auto* buf_first = tp.GetDenseBuffer<int>(
        std::type_index(typeid(
            nmr::LarsenHBondCountTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_first, nullptr);
    const std::size_t N_first = buf_first->AtomCount();
    const std::size_t T_first = buf_first->StridePerAtom();

    auto& tr = tp.Result<
        nmr::LarsenHBondCountTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);

    auto* buf_second = tp.GetDenseBuffer<int>(
        std::type_index(typeid(
            nmr::LarsenHBondCountTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_second, nullptr);
    EXPECT_EQ(buf_second->AtomCount(), N_first);
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
    EXPECT_EQ(tr.NumFrames(), T_first);
}


TEST(LarsenHBondCountTimeSeries, H5RoundTrip) {
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
    config.SetName("LarsenHBondCountTimeSeriesH5RoundTrip");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true;
    opts.skip_apbs = true;  opts.skip_dssp = false;  // backbone phi/psi needed for the larsen H-bond calc to actually run
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::DsspResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::LarsenHBondCountTimeSeriesTrajectoryResult
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
        nmr::LarsenHBondCountTimeSeriesTrajectoryResult>();

    const std::string h5_path = (fs::temp_directory_path() /
        ("larsen_hbond_count_ts_h5_roundtrip_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr.WriteH5Group(tp, file);
    }

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/larsen_hbond_count_time_series");
    auto ds = grp.getDataSet("count");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 2u);
    EXPECT_EQ(dims[0], tp.AtomCount());
    EXPECT_EQ(dims[1], 1u);

    std::string units, dtype;
    grp.getAttribute("units").read(units);
    grp.getAttribute("dtype").read(dtype);
    EXPECT_EQ(units, "pairs");
    EXPECT_EQ(dtype, "int32");

    fs::remove(h5_path);
}


// Integration: tally count distribution; log if any cell exceeds the
// physical sanity cap.
TEST(LarsenHBondCountTimeSeries, IntegrationDistribution1P9J) {
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
    config.SetName("LarsenHBondCountTimeSeriesIntegration");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true;
    opts.skip_apbs = true;  opts.skip_dssp = false;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::DsspResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::LarsenHBondCountTimeSeriesTrajectoryResult
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

    auto* buf = tp.GetDenseBuffer<int>(std::type_index(
        typeid(nmr::LarsenHBondCountTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);

    const std::size_t N = buf->AtomCount();
    const std::size_t T = buf->StridePerAtom();
    int min_count = std::numeric_limits<int>::max();
    int max_count = std::numeric_limits<int>::min();
    long long sum_count = 0;
    std::size_t overage_count = 0;
    std::size_t negative_count = 0;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const int c = buf->At(i, t);
            if (c < min_count) min_count = c;
            if (c > max_count) max_count = c;
            sum_count += c;
            if (c > kHBondCountSanityCap) ++overage_count;
            if (c < 0) ++negative_count;
        }
    }

    EXPECT_EQ(negative_count, 0u)
        << "found " << negative_count << " cells with negative count "
        << "(min=" << min_count << ")";

    // Real-data coverage floor: 1P9J has many H-bond donors and over
    // T frames at least some H-bond pairs are formed. max_count==0
    // means the calc didn't run (or ran and produced all zeros), which
    // is the exact regression the source-attached gate was added to
    // catch.
    EXPECT_GT(max_count, 0)
        << "no H-bonds found across " << T << " frames on 1P9J — "
        << "calc likely never ran despite grid being loaded";

    if (overage_count > 0) {
        std::cout << "[warn] LarsenHBondCount: " << overage_count
                  << " cells exceed " << kHBondCountSanityCap
                  << " (max=" << max_count << ") — chains of strong "
                  << "H-bonds, not asserted\n";
    }

    const double mean_count = static_cast<double>(sum_count)
                            / static_cast<double>(N * T);
    std::cout << "LarsenHBondCount diagnostic: frames=" << T
              << " atoms=" << N
              << " min=" << min_count << " max=" << max_count
              << " mean=" << mean_count
              << " overages=" << overage_count << "\n";
}
