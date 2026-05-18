//
// test_hydration_shell_welford: discipline + integration for
// HydrationShellWelfordTrajectoryResult. Per-atom Welford rollup of
// the COM-based hydration shell features (4 scalar channels).
//
// COVERAGE GAP (per c55c825 + test_hydration_geometry_welford.cpp): the
// synthetic source-absent test that the sibling TimeSeries TR carries is
// NOT replicated here — same Seed-with-real-positions limitation as the
// other Welford TRs. Source-attached gate discipline IS exercised on the
// TimeSeries side; the Welford prev_valid_ invalidation across a gap
// is a documented coverage gap that requires production-path
// `Trajectory::Run` infrastructure beyond the synthetic test framework.
//

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "GeometryResult.h"
#include "HydrationShellResult.h"
#include "HydrationShellWelfordTrajectoryResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "SpatialIndexResult.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryAtom.h"
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
    config.RequireConformationResult(typeid(nmr::HydrationShellResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::HydrationShellWelfordTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


TEST(HydrationShellWelford, Frame0Semantics) {
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

    for (std::size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& w = tp.AtomAt(i).hydration_shell_welford;
        EXPECT_EQ(w.n_frames, 1u);
        EXPECT_EQ(w.delta_n,  0u);
    }
}


TEST(HydrationShellWelford, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const std::size_t probe = tp.AtomCount() / 2;
    const double mean_first = tp.AtomAt(probe).hydration_shell_welford.half_shell_asymmetry.mean;
    auto& tr = tp.Result<nmr::HydrationShellWelfordTrajectoryResult>();
    tr.Finalize(tp, traj);
    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).hydration_shell_welford.half_shell_asymmetry.mean, mean_first);
}


TEST(HydrationShellWelford, H5RoundTrip) {
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

    const auto& tr = tp.Result<nmr::HydrationShellWelfordTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("hyd_shell_welford_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/hydration_shell_welford");

    std::string ref;
    grp.getAttribute("reference_frame").read(ref);
    EXPECT_EQ(ref, "COM");

    for (const std::string& base : {"half_shell_asymmetry", "mean_water_dipole_cos",
                                     "nearest_ion_distance", "nearest_ion_charge"}) {
        EXPECT_TRUE(grp.exist(base + "_mean"))            << base;
        EXPECT_TRUE(grp.exist(base + "_delta_mean"))      << base;
        EXPECT_TRUE(grp.exist(base + "_dxdt_mean"))       << base;
        EXPECT_TRUE(grp.exist(base + "_rms_delta"))       << base;
    }
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));
    EXPECT_TRUE(grp.exist("dxdt_n_per_atom"));

    fs::remove(h5_path);
}


TEST(HydrationShellWelford, Integration1P9J) {
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
    EXPECT_GE(traj.FrameCount(), 2u);

    // half_shell_asymmetry and mean_water_dipole_cos are always finite (no
    // inf sentinel). nearest_ion_* may be inf/NaN on atoms without an ion
    // in cutoff — log but don't assert.
    std::size_t populated_asymmetry = 0, populated_ion = 0;
    for (std::size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& w = tp.AtomAt(i).hydration_shell_welford;
        EXPECT_TRUE(std::isfinite(w.half_shell_asymmetry.mean));
        EXPECT_TRUE(std::isfinite(w.mean_water_dipole_cos.mean));
        EXPECT_EQ(w.n_frames, traj.FrameCount());
        if (std::abs(w.half_shell_asymmetry.mean) > 1e-12) ++populated_asymmetry;
        if (std::isfinite(w.nearest_ion_distance.mean)) ++populated_ion;
    }
    EXPECT_GT(populated_asymmetry, 0u);
    std::cout << "HydrationShellWelford: asymmetry populated=" << populated_asymmetry
              << "/" << tp.AtomCount()
              << ", ion_distance finite=" << populated_ion << "/" << tp.AtomCount() << "\n";
}
