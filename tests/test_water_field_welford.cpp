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
#include "PdbFileReader.h"
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
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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


TEST(WaterFieldWelford, H5DirectionalMetadataZeroCountSynthetic) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto build = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto tp_owner = nmr::TrajectoryProtein::CreateForTesting(
        std::move(build.protein));
    ASSERT_NE(tp_owner, nullptr);
    auto& tp = *tp_owner;
    ASSERT_GT(tp.AtomCount(), 1u);

    auto tr = nmr::WaterFieldWelfordTrajectoryResult::Create(tp);
    tr->ForceSourcePresentForTesting();
    nmr::Trajectory traj({}, {}, {});
    auto conf = tp.TickConformation(tp.CanonicalConformation().Positions());
    tr->Compute(*conf, tp, traj, 7, 2.5);
    tp.MutableAtomAt(0).water_field_welford = {};
    tr->Finalize(tp, traj);

    const std::string h5_path = (fs::temp_directory_path() /
        ("water_field_welford_directional_metadata_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr->WriteH5Group(tp, file);
    }
    {
        HighFive::File file(h5_path, HighFive::File::ReadOnly);
        auto grp = file.getGroup("/trajectory/water_field_welford");

        std::string irrep_scope, mean_law, component_law;
        std::string zero_count_validity;
        grp.getAttribute("irrep_metadata_scope").read(irrep_scope);
        grp.getAttribute("directional_mean_transformation").read(mean_law);
        grp.getAttribute("componentwise_statistic_transformation")
            .read(component_law);
        grp.getAttribute("zero_count_sentinel_validity")
            .read(zero_count_validity);
        EXPECT_EQ(irrep_scope,
                  "only assembled component means carry directional irrep metadata");
        EXPECT_EQ(mean_law,
                  "assembled efield_{x,y,z}_mean and "
                  "efield_first_{x,y,z}_mean are polar: v'=R v; assembled "
                  "efg_t2_mean and efg_first_t2_mean are even rank-2: "
                  "T'=R T R^T");
        EXPECT_EQ(component_law,
                  "componentwise m2,std,min,max,min_frame,max_frame have no "
                  "closed irrep transformation law");
        EXPECT_EQ(zero_count_validity,
                  "when n_frames_per_atom=0, mean,m2,std are NaN and min=+inf,"
                  "max=-inf,min_frame=0,max_frame=0 are invalid sentinels; "
                  "n_frames_per_atom gates validity");

        std::vector<std::size_t> counts;
        grp.getDataSet("n_frames_per_atom").read(counts);
        ASSERT_EQ(counts.size(), tp.AtomCount());
        EXPECT_EQ(counts[0], 0u);
        for (std::size_t i = 1; i < counts.size(); ++i) {
            EXPECT_EQ(counts[i], 1u);
        }

        auto expect_count_gated_zero = [&](const std::string& prefix,
                                            std::size_t width) {
            const std::size_t size = tp.AtomCount() * width;
            std::vector<double> mean(size), m2(size), stddev(size), min(size),
                                max(size);
            std::vector<std::size_t> min_frame(size), max_frame(size);
            grp.getDataSet(prefix + "_mean").read(mean.data());
            grp.getDataSet(prefix + "_m2").read(m2.data());
            grp.getDataSet(prefix + "_std").read(stddev.data());
            grp.getDataSet(prefix + "_min").read(min.data());
            grp.getDataSet(prefix + "_max").read(max.data());
            grp.getDataSet(prefix + "_min_frame").read(min_frame.data());
            grp.getDataSet(prefix + "_max_frame").read(max_frame.data());
            ASSERT_GE(mean.size(), 2 * width);
            for (std::size_t k = 0; k < width; ++k) {
                EXPECT_TRUE(std::isnan(mean[k]));
                EXPECT_TRUE(std::isnan(m2[k]));
                EXPECT_TRUE(std::isnan(stddev[k]));
                EXPECT_TRUE(std::isinf(min[k]) && min[k] > 0.0);
                EXPECT_TRUE(std::isinf(max[k]) && max[k] < 0.0);
                EXPECT_EQ(min_frame[k], 0u);
                EXPECT_EQ(max_frame[k], 0u);

                const std::size_t valid = width + k;
                EXPECT_DOUBLE_EQ(mean[valid], 0.0);
                EXPECT_DOUBLE_EQ(m2[valid], 0.0);
                EXPECT_DOUBLE_EQ(stddev[valid], 0.0);
                EXPECT_DOUBLE_EQ(min[valid], 0.0);
                EXPECT_DOUBLE_EQ(max[valid], 0.0);
                EXPECT_EQ(min_frame[valid], 7u);
                EXPECT_EQ(max_frame[valid], 7u);
            }
        };
        for (const std::string& prefix : {
                 "efield_x", "efield_y", "efield_z",
                 "efield_first_x", "efield_first_y", "efield_first_z"}) {
            expect_count_gated_zero(prefix, 1);
        }
        expect_count_gated_zero("efg_t2", 5);
        expect_count_gated_zero("efg_first_t2", 5);
    }
    fs::remove(h5_path);
}


TEST(WaterFieldWelford, H5RoundTrip) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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

    std::string units, ef, et2_basis, et2_order, et2;
    bool legacy_deprecated = false;
    grp.getAttribute("units").read(units);
    grp.getAttribute("irrep_layout_efield").read(ef);
    grp.getAttribute("efg_t2_basis").read(et2_basis);
    grp.getAttribute("efg_t2_component_order").read(et2_order);
    grp.getAttribute("irrep_layout_efg_t2").read(et2);
    grp.getAttribute("legacy_irrep_attrs_deprecated").read(legacy_deprecated);
    EXPECT_EQ(units, "V/Angstrom");
    EXPECT_EQ(ef, "v_x,v_y,v_z");
    EXPECT_EQ(et2_basis, "project_native_t2_isometric_real_tesseral_v1");
    EXPECT_EQ(et2_order, "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");
    EXPECT_EQ(et2, "m-2,m-1,m0,m+1,m+2");
    EXPECT_TRUE(legacy_deprecated);

    std::string irrep_scope, mean_law, component_law, zero_count_validity;
    grp.getAttribute("irrep_metadata_scope").read(irrep_scope);
    grp.getAttribute("directional_mean_transformation").read(mean_law);
    grp.getAttribute("componentwise_statistic_transformation").read(component_law);
    grp.getAttribute("zero_count_sentinel_validity").read(zero_count_validity);
    EXPECT_EQ(irrep_scope,
              "only assembled component means carry directional irrep metadata");
    EXPECT_EQ(mean_law,
              "assembled efield_{x,y,z}_mean and efield_first_{x,y,z}_mean "
              "are polar: v'=R v; assembled efg_t2_mean and "
              "efg_first_t2_mean are even rank-2: T'=R T R^T");
    EXPECT_EQ(component_law,
              "componentwise m2,std,min,max,min_frame,max_frame have no closed "
              "irrep transformation law");
    EXPECT_EQ(zero_count_validity,
              "when n_frames_per_atom=0, mean,m2,std are NaN and min=+inf,"
              "max=-inf,min_frame=0,max_frame=0 are invalid sentinels; "
              "n_frames_per_atom gates validity");

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
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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
