//
// test_gromacs_energy_time_series: discipline + integration for
// GromacsEnergyTimeSeriesTrajectoryResult. System-scalar per-frame TR;
// breaks the per-atom DenseBuffer pattern that the BS/HM/Sasa TS TRs
// follow. Storage is a flat std::vector<GromacsEnergy> indexed by frame.
//

#include "CalculatorConfig.h"
#include "GeometryResult.h"
#include "GromacsEnergyResult.h"
#include "GromacsEnergyTimeSeriesTrajectoryResult.h"
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
void LoadCalculatorConfig() {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(std::string(NMR_TEST_DATA_DIR) + "/../data/calculator_params.toml");
}

}  // namespace


// ============================================================================
// SYNTHETIC: attach a fake GromacsEnergyResult on each of three frames and
// confirm the TR rolls up the snapshots in order with the correct values.
// System-scalar TR cannot use the per-atom-field shortcut other TRs do;
// the source is GromacsEnergyResult itself, so we attach it.
// ============================================================================

TEST(GromacsEnergyTimeSeries, SyntheticThreeFrames) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    constexpr std::size_t kFrames = 3;
    for (std::size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");

        // Synthesize a GromacsEnergy snapshot with a recognizable pattern.
        nmr::GromacsEnergy ge;
        ge.time_ps      = static_cast<double>(t);
        ge.coulomb_sr   = -1000.0 - t;
        ge.coulomb_recip = -500.0 - t;
        ge.potential    = -5000.0 - 10.0 * t;
        ge.kinetic      = 2000.0 + t;
        ge.total_energy = ge.potential + ge.kinetic;
        ge.temperature  = 300.0 + 0.1 * t;
        ge.pressure     = 1.0 + 0.01 * t;
        ge.volume       = 100.0 + t;
        ge.density      = 1000.0 - t;
        ge.box_x        = ge.box_y = ge.box_z = std::cbrt(ge.volume);
        for (int k = 0; k < 9; ++k) ge.vir[k]  = static_cast<double>(t * 100 + k);
        for (int k = 0; k < 9; ++k) ge.pres[k] = static_cast<double>(t * 10  + k);

        conf->AttachResult(nmr::GromacsEnergyResult::Compute(*conf, ge));
        tr->Compute(*conf, tp, traj, t, ge.time_ps);
    }
    tr->Finalize(tp, traj);
    EXPECT_EQ(tr->NumFrames(), kFrames);

    // H5 round-trip the snapshots and read back the channels.
    const std::string h5_path = (fs::temp_directory_path() /
        ("gromacs_energy_ts_unit_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/gromacs_energy_time_series"));
    auto grp = reopen.getGroup("/trajectory/gromacs_energy_time_series");

    std::vector<double> potential;
    grp.getDataSet("potential").read(potential);
    ASSERT_EQ(potential.size(), kFrames);
    for (std::size_t t = 0; t < kFrames; ++t)
        EXPECT_DOUBLE_EQ(potential[t], -5000.0 - 10.0 * t);

    std::vector<std::vector<double>> vir;
    grp.getDataSet("virial").read(vir);
    ASSERT_EQ(vir.size(), kFrames);
    ASSERT_EQ(vir[0].size(), 9u);
    EXPECT_DOUBLE_EQ(vir[2][3], 2 * 100 + 3);  // (frame=2, k=3)

    std::string units, tensor_layout;
    grp.getAttribute("units").read(units);
    grp.getAttribute("tensor_layout").read(tensor_layout);
    EXPECT_EQ(units, "kJ/mol");
    EXPECT_EQ(tensor_layout, "XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ");

    fs::remove(h5_path);
}


// ============================================================================
// SYNTHETIC: source-absent path. Mix 2 attached + 2 absent frames, verify
// the source-attached gate behavior — H5 group lists source_attached_count=2
// and absent frames carry NaN on the energy datasets. R2 review 2026-05-18:
// the gate logic was uncovered until this test landed.
// ============================================================================

TEST(GromacsEnergyTimeSeries, SyntheticSourceAbsentFrames) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr std::size_t kFrames = 4;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        // Attach GromacsEnergyResult only on EVEN frames (0, 2). Odd
        // frames (1, 3) are source-absent.
        if (t % 2 == 0) {
            nmr::GromacsEnergy ge;
            ge.time_ps     = static_cast<double>(t);
            ge.potential   = -1000.0 - 10.0 * t;
            ge.temperature = 300.0;
            conf->AttachResult(nmr::GromacsEnergyResult::Compute(*conf, ge));
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    // H5 round-trip + verify NaN on absent frames.
    const std::string h5_path = (fs::temp_directory_path() /
        ("gromacs_energy_ts_absent_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/gromacs_energy_time_series");

    std::size_t source_attached_count = 0;
    grp.getAttribute("source_attached_count").read(source_attached_count);
    EXPECT_EQ(source_attached_count, 2u)
        << "Expected 2/4 frames source-attached (even frames only)";

    std::vector<std::uint8_t> mask;
    grp.getDataSet("source_attached_per_frame").read(mask);
    ASSERT_EQ(mask.size(), kFrames);
    EXPECT_EQ(mask[0], 1u); EXPECT_EQ(mask[1], 0u);
    EXPECT_EQ(mask[2], 1u); EXPECT_EQ(mask[3], 0u);

    std::vector<double> potential;
    grp.getDataSet("potential").read(potential);
    ASSERT_EQ(potential.size(), kFrames);
    // Frame 0, 2: finite values (attached). Frame 1, 3: NaN.
    EXPECT_TRUE(std::isfinite(potential[0]));
    EXPECT_TRUE(std::isnan(potential[1]));
    EXPECT_TRUE(std::isfinite(potential[2]));
    EXPECT_TRUE(std::isnan(potential[3]));
    EXPECT_DOUBLE_EQ(potential[0], -1000.0);
    EXPECT_DOUBLE_EQ(potential[2], -1020.0);

    // energy_frame_times_ps NaN-filled too.
    std::vector<double> et;
    grp.getDataSet("energy_frame_times_ps").read(et);
    EXPECT_TRUE(std::isfinite(et[0]));
    EXPECT_TRUE(std::isnan(et[1]));

    fs::remove(h5_path);
}


// ============================================================================
// SYNTHETIC: ALL frames source-absent → WriteH5Group skips emission entirely.
// ============================================================================

TEST(GromacsEnergyTimeSeries, SyntheticAllAbsentSkipsGroup) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < 3; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        // No AttachResult — source is absent for every frame.
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);
    EXPECT_EQ(tr->NumFrames(), 3u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("gromacs_energy_ts_allabsent_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/gromacs_energy_time_series"))
        << "All-absent run should skip group emission entirely.";

    fs::remove(h5_path);
}


// ============================================================================
// DISCIPLINE: Frame-0 semantics — single-frame fixture exercises Compute once.
// ============================================================================

TEST(GromacsEnergyTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::GromacsEnergyResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_EQ(traj.FrameCount(), 1u);

    const auto& tr = tp.Result<nmr::GromacsEnergyTimeSeriesTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);
}


// ============================================================================
// DISCIPLINE: Finalize idempotency — second Finalize must not corrupt state.
// ============================================================================

TEST(GromacsEnergyTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::GromacsEnergyResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::GromacsEnergyTimeSeriesTrajectoryResult>();
    const std::size_t T_first = tr.NumFrames();
    tr.Finalize(tp, traj);
    EXPECT_EQ(tr.NumFrames(), T_first);
}


// ============================================================================
// DISCIPLINE: H5 round-trip — group + all expected datasets land.
// ============================================================================

TEST(GromacsEnergyTimeSeries, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::GromacsEnergyResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::GromacsEnergyTimeSeriesTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("gromacs_energy_ts_h5_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/gromacs_energy_time_series"));
    auto grp = reopen.getGroup("/trajectory/gromacs_energy_time_series");

    std::string units;
    grp.getAttribute("units").read(units);
    EXPECT_EQ(units, "kJ/mol");

    // Expected datasets — sample 4 representative ones across categories.
    EXPECT_TRUE(grp.exist("potential"));
    EXPECT_TRUE(grp.exist("temperature"));
    EXPECT_TRUE(grp.exist("virial"));
    EXPECT_TRUE(grp.exist("frame_times"));

    // Virial / pressure_tensor shape: (T, 9).
    const auto vir_dims = grp.getDataSet("virial").getSpace().getDimensions();
    ASSERT_EQ(vir_dims.size(), 2u);
    EXPECT_EQ(vir_dims[1], 9u);

    fs::remove(h5_path);
}


// ============================================================================
// INTEGRATION: real .edr energies through Trajectory::Run, multi-frame.
// ============================================================================

TEST(GromacsEnergyTimeSeries, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::GromacsEnergyResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::GromacsEnergyTimeSeriesTrajectoryResult>();
    EXPECT_GE(tr.NumFrames(), 2u);

    // H5 round-trip to inspect populated channels — confirm finite values
    // on a sample of the physically-meaningful channels.
    const std::string h5_path = (fs::temp_directory_path() /
        ("gromacs_energy_ts_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/gromacs_energy_time_series");

    std::vector<double> total, temperature, pressure;
    grp.getDataSet("total_energy").read(total);
    grp.getDataSet("temperature").read(temperature);
    grp.getDataSet("pressure").read(pressure);
    ASSERT_EQ(total.size(), tr.NumFrames());

    for (double v : total)       EXPECT_TRUE(std::isfinite(v));
    for (double T : temperature) EXPECT_TRUE(std::isfinite(T));
    for (double P : pressure)    EXPECT_TRUE(std::isfinite(P));

    // Loose physical-sanity: thermostatted MD around 300 K within ±50 K.
    // Per feedback_log_overages_dont_assert this is logged, not asserted.
    for (double T : temperature) {
        if (T < 250.0 || T > 400.0) {
            std::cerr << "WARN: temperature out of [250,400] K: " << T << "\n";
        }
    }
    std::cout << "GromacsEnergyTimeSeries: " << tr.NumFrames() << " frames; "
              << "total_energy[0]=" << total[0]
              << " T[0]=" << temperature[0]
              << " P[0]=" << pressure[0] << "\n";

    fs::remove(h5_path);
}
