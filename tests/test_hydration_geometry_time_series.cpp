//
// test_hydration_geometry_time_series: discipline + integration for
// HydrationGeometryTimeSeriesTrajectoryResult. Per-atom 6-channel TR
// (Vec3 dipole + Vec3 normal + uint32 count + 3 scalars).
//

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "GeometryResult.h"
#include "HydrationGeometryResult.h"
#include "HydrationGeometryTimeSeriesTrajectoryResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "SasaResult.h"
#include "Session.h"
#include "SpatialIndexResult.h"
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
void LoadCalculatorConfig() {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(std::string(NMR_TEST_DATA_DIR) + "/../data/calculator_params.toml");
}

nmr::RunConfiguration BuildConfig(unsigned stride) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::SasaResult));
    config.RequireConformationResult(typeid(nmr::HydrationGeometryResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::HydrationGeometryTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


TEST(HydrationGeometryTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_EQ(traj.FrameCount(), 1u);
    const auto& tr = tp.Result<nmr::HydrationGeometryTimeSeriesTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);
}


TEST(HydrationGeometryTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    auto& tr = tp.Result<nmr::HydrationGeometryTimeSeriesTrajectoryResult>();
    const std::size_t T = tr.NumFrames();
    tr.Finalize(tp, traj);
    EXPECT_EQ(tr.NumFrames(), T);
}


TEST(HydrationGeometryTimeSeries, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::HydrationGeometryTimeSeriesTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("hyd_geo_ts_h5_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/hydration_geometry_time_series"));
    auto grp = reopen.getGroup("/trajectory/hydration_geometry_time_series");

    std::string layout, polar_chans;
    grp.getAttribute("dipole_vector_layout").read(layout);
    grp.getAttribute("polarisation_signal_channels").read(polar_chans);
    EXPECT_EQ(layout, "x,y,z");
    EXPECT_EQ(polar_chans, "dipole_alignment,dipole_coherence,half_shell_asymmetry");

    const std::size_t N = tp.AtomCount();
    const std::size_t T = tr.NumFrames();
    const auto dipole_dims = grp.getDataSet("dipole_vector").getSpace().getDimensions();
    ASSERT_EQ(dipole_dims.size(), 3u);
    EXPECT_EQ(dipole_dims[0], N); EXPECT_EQ(dipole_dims[1], T); EXPECT_EQ(dipole_dims[2], 3u);

    EXPECT_TRUE(grp.exist("dipole_alignment"));
    EXPECT_TRUE(grp.exist("dipole_coherence"));
    EXPECT_TRUE(grp.exist("half_shell_asymmetry"));
    EXPECT_TRUE(grp.exist("first_shell_count"));
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));

    fs::remove(h5_path);
}


TEST(HydrationGeometryTimeSeries, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(300);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::HydrationGeometryTimeSeriesTrajectoryResult>();
    EXPECT_GE(tr.NumFrames(), 2u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("hyd_geo_ts_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/hydration_geometry_time_series");

    const std::size_t N = tp.AtomCount();
    const std::size_t T = tr.NumFrames();
    std::vector<std::vector<double>> alignment;
    grp.getDataSet("dipole_alignment").read(alignment);
    ASSERT_EQ(alignment.size(), N);

    std::size_t populated = 0;
    double max_abs_align = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            EXPECT_TRUE(std::isfinite(alignment[i][t]));
            // alignment is cos angle ∈ [-1, 1]
            EXPECT_GE(alignment[i][t], -1.001);
            EXPECT_LE(alignment[i][t],  1.001);
            max_abs_align = std::max(max_abs_align, std::abs(alignment[i][t]));
        }
        if (std::abs(alignment[i][0]) > 1e-12) ++populated;
    }
    EXPECT_GT(populated, 0u) << "dipole_alignment all zero on a solvated protein";
    std::cout << "HydrationGeometryTimeSeries: " << T << " frames; "
              << "populated=" << populated << "/" << N
              << " max|alignment|=" << max_abs_align << "\n";

    fs::remove(h5_path);
}


TEST(HydrationGeometryTimeSeries, SyntheticAllAbsentSkipsGroup) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::HydrationGeometryTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < 3; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    const std::string h5_path = (fs::temp_directory_path() /
        ("hyd_geo_ts_absent_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/hydration_geometry_time_series"))
        << "All-absent run should skip group emission.";

    fs::remove(h5_path);
}
