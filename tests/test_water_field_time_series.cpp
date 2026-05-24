//
// test_water_field_time_series: discipline + integration for
// WaterFieldTimeSeriesTrajectoryResult. Per-atom 6-channel TR (Vec3 + Vec3
// + SphericalTensor + SphericalTensor + int + int) over WaterFieldResult.
//

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "GeometryResult.h"
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
#include "WaterFieldResult.h"
#include "WaterFieldTimeSeriesTrajectoryResult.h"

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

nmr::RunConfiguration BuildWaterFieldConfig(unsigned stride) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::SasaResult));
    config.RequireConformationResult(typeid(nmr::WaterFieldResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::WaterFieldTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


// ============================================================================
// SYNTHETIC: source-absent path — verifies all-absent → group skipped.
// R3 codex F5 2026-05-18.
// ============================================================================

TEST(WaterFieldTimeSeries, SyntheticAllAbsentSkipsGroup) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::WaterFieldTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < 3; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        // No WaterFieldResult attached — source absent every frame.
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    const std::string h5_path = (fs::temp_directory_path() /
        ("water_ts_allabsent_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/water_field_time_series"))
        << "All-absent run should skip group emission entirely.";

    fs::remove(h5_path);
}


TEST(WaterFieldTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildWaterFieldConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_EQ(traj.FrameCount(), 1u);

    const auto& tr = tp.Result<nmr::WaterFieldTimeSeriesTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);
}


TEST(WaterFieldTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildWaterFieldConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::WaterFieldTimeSeriesTrajectoryResult>();
    const std::size_t T = tr.NumFrames();
    tr.Finalize(tp, traj);
    EXPECT_EQ(tr.NumFrames(), T);
}


TEST(WaterFieldTimeSeries, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildWaterFieldConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::WaterFieldTimeSeriesTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("water_field_ts_h5_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/water_field_time_series"));
    auto grp = reopen.getGroup("/trajectory/water_field_time_series");

    std::string efield_units, efg_layout, count_units, efg_parity;
    grp.getAttribute("efield_units").read(efield_units);
    grp.getAttribute("efg_irrep_layout").read(efg_layout);
    grp.getAttribute("efg_parity").read(efg_parity);
    grp.getAttribute("count_units").read(count_units);
    EXPECT_EQ(efield_units, "V/Angstrom");
    // Water EFG: T0 (trace) and T1 (antisym pseudovector) are both
    // structural zeros — only T2 (symmetric-traceless, 5 components)
    // is emitted. Parity "2e" reflects T2 only.
    EXPECT_EQ(efg_layout, "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");
    EXPECT_EQ(efg_parity, "2e");
    EXPECT_EQ(count_units, "dimensionless");

    bool t0_zero = false, t1_zero = false;
    grp.getAttribute("efg_t0_structural_zero").read(t0_zero);
    grp.getAttribute("efg_t1_structural_zero").read(t1_zero);
    EXPECT_TRUE(t0_zero);
    EXPECT_TRUE(t1_zero);

    // Shapes
    const std::size_t N = tp.AtomCount();
    const std::size_t T = tr.NumFrames();
    const auto efield_dims = grp.getDataSet("efield").getSpace().getDimensions();
    ASSERT_EQ(efield_dims.size(), 3u);
    EXPECT_EQ(efield_dims[0], N); EXPECT_EQ(efield_dims[1], T); EXPECT_EQ(efield_dims[2], 3u);
    const auto efg_dims = grp.getDataSet("efg").getSpace().getDimensions();
    ASSERT_EQ(efg_dims.size(), 3u);
    EXPECT_EQ(efg_dims[0], N); EXPECT_EQ(efg_dims[1], T); EXPECT_EQ(efg_dims[2], 5u);
    const auto n_first_dims = grp.getDataSet("n_first").getSpace().getDimensions();
    ASSERT_EQ(n_first_dims.size(), 2u);
    EXPECT_EQ(n_first_dims[0], N); EXPECT_EQ(n_first_dims[1], T);

    fs::remove(h5_path);
}


TEST(WaterFieldTimeSeries, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildWaterFieldConfig(300);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::WaterFieldTimeSeriesTrajectoryResult>();
    EXPECT_GE(tr.NumFrames(), 2u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("water_field_ts_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/water_field_time_series");

    const std::size_t N = tp.AtomCount();
    const std::size_t T = tr.NumFrames();

    // HighFive supports nested-vector read for multi-dim datasets.
    std::vector<std::vector<std::vector<double>>> efield;
    grp.getDataSet("efield").read(efield);
    std::vector<std::vector<std::uint32_t>> n_first;
    grp.getDataSet("n_first").read(n_first);
    ASSERT_EQ(efield.size(),  N);
    ASSERT_EQ(efield[0].size(), T);
    ASSERT_EQ(efield[0][0].size(), 3u);
    ASSERT_EQ(n_first.size(),  N);
    ASSERT_EQ(n_first[0].size(), T);

    // Populated check: a real solvated protein should see nonzero E-field
    // on every atom (the cutoff is 15 Å and water surrounds it).
    std::size_t populated = 0;
    double max_efield_mag = 0.0;
    std::uint32_t max_n_first = 0;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const double Ex = efield[i][t][0];
            const double Ey = efield[i][t][1];
            const double Ez = efield[i][t][2];
            EXPECT_TRUE(std::isfinite(Ex));
            EXPECT_TRUE(std::isfinite(Ey));
            EXPECT_TRUE(std::isfinite(Ez));
            const double mag = std::sqrt(Ex*Ex + Ey*Ey + Ez*Ez);
            max_efield_mag = std::max(max_efield_mag, mag);
            max_n_first = std::max(max_n_first, n_first[i][t]);
        }
        if (std::abs(efield[i][0][0]) > 1e-12 ||
            std::abs(efield[i][0][1]) > 1e-12 ||
            std::abs(efield[i][0][2]) > 1e-12) {
            ++populated;
        }
    }
    EXPECT_GT(populated, N / 2)
        << "Water E-field populated on < 50% of atoms — solvent not loaded?";
    std::cout << "WaterFieldTimeSeries: " << T << " frames; "
              << "max |E|=" << max_efield_mag << " V/A; "
              << "max n_first=" << max_n_first << "; "
              << "populated=" << populated << "/" << N << "\n";

    fs::remove(h5_path);
}
