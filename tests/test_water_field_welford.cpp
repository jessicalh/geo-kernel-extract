//
// test_water_field_welford: discipline + integration for
// WaterFieldWelfordTrajectoryResult. Per-atom Welford rollup of the
// explicit-water E-field + EFG kernel. Discipline trio + Integration1P9J
// with population floor on the primary scalars + dxdt_n alignment check.
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
#include "TrajectoryAtom.h"
#include "TrajectoryProtein.h"
#include "Types.h"
#include "WaterFieldResult.h"
#include "WaterFieldWelfordTrajectoryResult.h"

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

nmr::RunConfiguration BuildWaterFieldWelfordConfig(unsigned stride) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::SasaResult));
    config.RequireConformationResult(typeid(nmr::WaterFieldResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::WaterFieldWelfordTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


// Source-absent path is exercised on GromacsEnergyTimeSeries via its
// synthetic test pair (SyntheticSourceAbsentFrames + SyntheticAllAbsentSkipsGroup).
// The WaterFieldWelford-specific path that invalidates prev_valid_ across
// a gap requires `Seed()` with real first-frame positions to allocate
// TrajectoryAtoms — the synthetic-zero-positions path FATALs at
// canonicalization. The discipline is identical to GromacsEnergyTS; the
// WaterFieldWelford-specific cache invalidation is a coverage gap noted
// in the R2 review followup commit. Integration1P9J below verifies the
// happy path (all frames source-attached) including dxdt_n == delta_n.


TEST(WaterFieldWelford, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildWaterFieldWelfordConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_EQ(traj.FrameCount(), 1u);

    for (std::size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& w = tp.AtomAt(i).water_field_welford;
        EXPECT_EQ(w.n_frames, 1u);
        EXPECT_EQ(w.delta_n,  0u);
        EXPECT_DOUBLE_EQ(w.efg_t2magnitude.std, 0.0);
    }
}


TEST(WaterFieldWelford, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildWaterFieldWelfordConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const std::size_t probe = tp.AtomCount() / 2;
    const double mean_first = tp.AtomAt(probe).water_field_welford.efg_t2magnitude.mean;
    const double std_first  = tp.AtomAt(probe).water_field_welford.efg_t2magnitude.std;
    auto& tr = tp.Result<nmr::WaterFieldWelfordTrajectoryResult>();
    tr.Finalize(tp, traj);
    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).water_field_welford.efg_t2magnitude.mean, mean_first);
    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).water_field_welford.efg_t2magnitude.std,  std_first);
}


TEST(WaterFieldWelford, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildWaterFieldWelfordConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::WaterFieldWelfordTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("water_field_welford_h5_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/water_field_welford");

    std::string units, ef, et2;
    grp.getAttribute("units").read(units);
    grp.getAttribute("irrep_layout_efield").read(ef);
    grp.getAttribute("irrep_layout_efg_t2").read(et2);
    EXPECT_EQ(units, "V/Angstrom");
    EXPECT_EQ(ef, "v_x,v_y,v_z");
    EXPECT_EQ(et2, "m-2,m-1,m0,m+1,m+2");

    // efg_t1_* datasets must NOT be emitted (structurally zero).
    EXPECT_FALSE(grp.exist("efg_t1_mean"));
    EXPECT_FALSE(grp.exist("efg_first_t1_mean"));
    // efg_t2 per-component shape is (N, 5).
    const auto t2_dims = grp.getDataSet("efg_t2_mean").getSpace().getDimensions();
    ASSERT_EQ(t2_dims.size(), 2u); EXPECT_EQ(t2_dims[1], 5u);
    bool t1_zero = false;
    grp.getAttribute("efg_t1_structural_zero").read(t1_zero);
    EXPECT_TRUE(t1_zero);

    // Delta variants on the 3 primary scalars (efg_t0 deltas absent —
    // T0 structurally zero, channel removed per 2026-05-18 review).
    for (const std::string& base : {"efield_magnitude", "n_first", "n_second"}) {
        EXPECT_TRUE(grp.exist(base + "_delta_mean")) << base;
        EXPECT_TRUE(grp.exist(base + "_abs_delta_mean")) << base;
        EXPECT_TRUE(grp.exist(base + "_dxdt_mean")) << base;
        EXPECT_TRUE(grp.exist(base + "_rms_delta")) << base;
    }
    // efg_t0 datasets MUST NOT be emitted (structural zero).
    EXPECT_FALSE(grp.exist("efg_t0_mean"));
    EXPECT_FALSE(grp.exist("efg_t0_delta_mean"));
    // Provenance flag.
    bool t0_zero = false;
    grp.getAttribute("efg_t0_structural_zero").read(t0_zero);
    EXPECT_TRUE(t0_zero);
    EXPECT_TRUE(grp.exist("dxdt_n_per_atom"));

    fs::remove(h5_path);
}


TEST(WaterFieldWelford, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildWaterFieldWelfordConfig(300);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_GE(traj.FrameCount(), 2u);

    std::size_t populated = 0, t2_populated = 0;
    double max_abs_efield_mag = 0.0, max_abs_t2 = 0.0;
    for (std::size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& w = tp.AtomAt(i).water_field_welford;
        EXPECT_TRUE(std::isfinite(w.efield_magnitude.mean));
        EXPECT_TRUE(std::isfinite(w.efg_t2magnitude.mean));
        EXPECT_EQ(w.n_frames, traj.FrameCount());
        if (std::abs(w.efield_magnitude.mean) > 1e-12) ++populated;
        max_abs_efield_mag = std::max(max_abs_efield_mag, std::abs(w.efield_magnitude.mean));
        for (std::size_t k = 0; k < 5; ++k) {
            EXPECT_TRUE(std::isfinite(w.efg_t2[k].mean));
            const double v = std::abs(w.efg_t2[k].mean);
            max_abs_t2 = std::max(max_abs_t2, v);
            if (v > 1e-12) { ++t2_populated; break; }
        }
        // dxdt_n must equal delta_n on a well-formed trajectory.
        EXPECT_EQ(w.dxdt_n, w.delta_n);
    }
    EXPECT_GT(populated, tp.AtomCount() / 2)
        << "Water E-field magnitude populated on < 50% of atoms";
    EXPECT_GT(t2_populated, 0u)
        << "Water EFG T2 all-zero — per-component regression";
    std::cout << "WaterFieldWelford: populated=" << populated << "/" << tp.AtomCount()
              << " max|E_mag|=" << max_abs_efield_mag
              << " max|t2|=" << max_abs_t2 << "\n";
}
