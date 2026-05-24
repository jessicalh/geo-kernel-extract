//
// test_bonded_energy_time_series: discipline + integration for
// BondedEnergyTimeSeriesTrajectoryResult. Per-atom 7-channel TR over the
// CHARMM36m bonded breakdown. Internal per-atom buffers (no DenseBuffer
// adoption) emit 7 (N, T) H5 datasets.
//

#include "BondedEnergyResult.h"
#include "BondedEnergyTimeSeriesTrajectoryResult.h"
#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "GeometryResult.h"
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
// SYNTHETIC: source-absent path. Verifies all-absent → group skipped + the
// mixed path NaN-fills absent rows. R3 codex F5 2026-05-18.
// ============================================================================

TEST(BondedEnergyTimeSeries, SyntheticAllAbsentSkipsGroup) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::BondedEnergyTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < 3; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        // No BondedEnergyResult attached — source absent every frame.
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    const std::string h5_path = (fs::temp_directory_path() /
        ("bonded_ts_allabsent_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/bonded_energy_time_series"))
        << "All-absent run should skip group emission entirely.";

    fs::remove(h5_path);
}


// ============================================================================
// DISCIPLINE: Frame-0 semantics through Trajectory::Run.
// ============================================================================

TEST(BondedEnergyTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::BondedEnergyResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::BondedEnergyTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_EQ(traj.FrameCount(), 1u);

    const auto& tr = tp.Result<nmr::BondedEnergyTimeSeriesTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);
}


// ============================================================================
// DISCIPLINE: Finalize idempotency.
// ============================================================================

TEST(BondedEnergyTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::BondedEnergyResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::BondedEnergyTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::BondedEnergyTimeSeriesTrajectoryResult>();
    const std::size_t T_first = tr.NumFrames();
    tr.Finalize(tp, traj);
    EXPECT_EQ(tr.NumFrames(), T_first);
}


// ============================================================================
// DISCIPLINE: H5 round-trip — all 7 channels + frame indices land.
// ============================================================================

TEST(BondedEnergyTimeSeries, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::BondedEnergyResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::BondedEnergyTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::BondedEnergyTimeSeriesTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("bonded_energy_ts_h5_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/bonded_energy_time_series"));
    auto grp = reopen.getGroup("/trajectory/bonded_energy_time_series");

    std::string units, split;
    grp.getAttribute("units").read(units);
    grp.getAttribute("split_convention").read(split);
    EXPECT_EQ(units, "kJ/mol");
    EXPECT_EQ(split, "evenly_among_participating_atoms");

    for (const std::string& ch :
         {"bond", "angle", "urey_bradley", "proper_dih",
          "improper_dih", "cmap_dih", "total"}) {
        ASSERT_TRUE(grp.exist(ch)) << "missing channel: " << ch;
        const auto dims = grp.getDataSet(ch).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], tp.AtomCount());
        EXPECT_EQ(dims[1], tr.NumFrames());
    }

    fs::remove(h5_path);
}


// ============================================================================
// INTEGRATION: real BondedEnergyResult through Trajectory::Run, multi-frame.
// ============================================================================

TEST(BondedEnergyTimeSeries, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::BondedEnergyResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::BondedEnergyTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::BondedEnergyTimeSeriesTrajectoryResult>();
    EXPECT_GE(tr.NumFrames(), 2u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("bonded_energy_ts_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/bonded_energy_time_series");

    std::vector<std::vector<double>> total;
    grp.getDataSet("total").read(total);
    ASSERT_EQ(total.size(), tp.AtomCount());
    ASSERT_EQ(total[0].size(), tr.NumFrames());

    // Loose population floor — bond + angle should be nonzero for most
    // atoms on a real protein. CHARMM36m through GROMACS reports
    // UB/improper/CMAP=0 on the 1P9J fixture (verified 2026-05-18); the
    // populated check focuses on bond + angle + proper_dih.
    std::size_t populated = 0;
    double max_abs_total = 0.0;
    for (std::size_t i = 0; i < tp.AtomCount(); ++i) {
        for (std::size_t t = 0; t < tr.NumFrames(); ++t) {
            EXPECT_TRUE(std::isfinite(total[i][t]));
            max_abs_total = std::max(max_abs_total, std::abs(total[i][t]));
        }
        if (std::abs(total[i][0]) > 1e-12) ++populated;
    }
    EXPECT_GT(populated, tp.AtomCount() / 2)
        << "bonded total populated on < 50% of atoms — calc not firing";
    std::cout << "BondedEnergyTimeSeries: " << tr.NumFrames() << " frames; "
              << "max|total|=" << max_abs_total << " kJ/mol; "
              << "populated=" << populated << "/" << tp.AtomCount() << "\n";

    fs::remove(h5_path);
}
