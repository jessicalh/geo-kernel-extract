//
// test_hydration_geometry_welford: discipline + integration for
// HydrationGeometryWelfordTrajectoryResult. Per-atom Welford rollup of the
// SASA-normal water polarisation features.
//
// COVERAGE GAP (per c55c825 docstring on test_water_field_welford.cpp): the
// synthetic source-absent test that the sibling TimeSeries TR carries is
// NOT replicated here. Reason: the Welford TR's Compute writes to
// `tp.MutableAtomAt(i).hydration_geometry_welford`, which requires
// TrajectoryAtoms allocated via `tp.Seed()`. The synthetic Seed path with
// all-zero positions FATALs at the protein canonicalisation step
// (substrate gap on missing-after-canonical atom names). Production code
// goes through `Trajectory::Run` which seeds with real frame-0 positions;
// the synthetic-positions path can't reach that state. Source-attached
// gate discipline IS exercised on the TimeSeries side; the Welford-side
// prev_valid_ invalidation across a gap remains a documented coverage gap.
//

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "GeometryResult.h"
#include "HydrationGeometryResult.h"
#include "HydrationGeometryWelfordTrajectoryResult.h"
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
    config.RequireConformationResult(typeid(nmr::SasaResult));
    config.RequireConformationResult(typeid(nmr::HydrationGeometryResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::HydrationGeometryWelfordTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


TEST(HydrationGeometryWelford, Frame0Semantics) {
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
        const auto& w = tp.AtomAt(i).hydration_geometry_welford;
        EXPECT_EQ(w.n_frames, 1u);
        EXPECT_EQ(w.delta_n,  0u);
        EXPECT_DOUBLE_EQ(w.dipole_alignment.std, 0.0);
    }
}


TEST(HydrationGeometryWelford, FinalizeIdempotency) {
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
    const double mean_first = tp.AtomAt(probe).hydration_geometry_welford.dipole_alignment.mean;
    auto& tr = tp.Result<nmr::HydrationGeometryWelfordTrajectoryResult>();
    tr.Finalize(tp, traj);
    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).hydration_geometry_welford.dipole_alignment.mean, mean_first);
}


TEST(HydrationGeometryWelford, H5RoundTrip) {
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

    const auto& tr = tp.Result<nmr::HydrationGeometryWelfordTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("hyd_geo_welford_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/hydration_geometry_welford");

    for (const std::string& base : {"half_shell_asymmetry", "dipole_alignment",
                                     "dipole_coherence", "shell_count"}) {
        EXPECT_TRUE(grp.exist(base + "_delta_mean"))     << base;
        EXPECT_TRUE(grp.exist(base + "_abs_delta_mean")) << base;
        EXPECT_TRUE(grp.exist(base + "_dxdt_mean"))      << base;
        EXPECT_TRUE(grp.exist(base + "_rms_delta"))      << base;
    }
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));
    EXPECT_TRUE(grp.exist("dxdt_n_per_atom"));

    fs::remove(h5_path);
}


TEST(HydrationGeometryWelford, Integration1P9J) {
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

    std::size_t populated = 0;
    double max_abs_alignment = 0.0;
    for (std::size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& w = tp.AtomAt(i).hydration_geometry_welford;
        EXPECT_TRUE(std::isfinite(w.dipole_alignment.mean));
        EXPECT_TRUE(std::isfinite(w.dipole_coherence.mean));
        EXPECT_TRUE(std::isfinite(w.shell_count.mean));
        EXPECT_EQ(w.n_frames, traj.FrameCount());
        EXPECT_EQ(w.dxdt_n, w.delta_n)
            << "atom " << i << ": dxdt_n != delta_n";
        // alignment ∈ [-1, 1], coherence ∈ [0, 1] when computed
        EXPECT_GE(w.dipole_alignment.mean, -1.001);
        EXPECT_LE(w.dipole_alignment.mean,  1.001);
        if (std::abs(w.dipole_alignment.mean) > 1e-12) ++populated;
        max_abs_alignment = std::max(max_abs_alignment, std::abs(w.dipole_alignment.mean));
    }
    EXPECT_GT(populated, 0u) << "dipole_alignment all zero on solvated protein";
    std::cout << "HydrationGeometryWelford: populated=" << populated << "/"
              << tp.AtomCount() << " max|alignment|=" << max_abs_alignment << "\n";
}
