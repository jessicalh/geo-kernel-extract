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

#include <cstring>

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "GeometryResult.h"
#include "HydrationGeometryResult.h"
#include "HydrationGeometryWelfordTrajectoryResult.h"
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
constexpr const char* kSurfaceNormalTransformation =
    "outward polar vector: v'=R v in the continuum limit; live upstream "
    "finite lab-fixed Fibonacci estimator has no exact O(3) law and is "
    "only approximately rotation-covariant within the recorded test envelope";
constexpr const char* kSourceComponentSemantics =
    "dipole_vector and surface_normal are independently accumulated "
    "Cartesian x,y,z source components; they are not SphericalTensor T1 "
    "components";
constexpr const char* kIrrepMetadataScope =
    "only assembled dipole component means carry exact irrep metadata; "
    "assembled surface-normal means carry only the declared "
    "continuum/finite-grid approximate directional law";
constexpr const char* kDirectionalMeanTransformation =
    "assembled dipole_vector_{x,y,z}_mean: polar vector v'=R v; assembled "
    "surface_normal_{x,y,z}_mean: outward polar vector v'=R v in the "
    "continuum limit, but the live upstream finite lab-fixed Fibonacci "
    "estimator has no exact O(3) law and is only approximately "
    "rotation-covariant within the recorded test envelope";
constexpr const char* kScalarStatisticTransformation =
    "for B in {half_shell_asymmetry,dipole_alignment}: "
    "B_{mean,m2,std,min,max,min_frame,max_frame}, "
    "B_{delta,abs_delta,delta_squared,dxdt}_"
    "{mean,m2,std,min,max,min_frame,max_frame}, and B_rms_delta are "
    "continuum rotation invariants that inherit the live upstream finite "
    "lab-fixed Fibonacci estimator: no exact O(3) law and only approximate "
    "rotation covariance within the recorded test envelope";
constexpr const char* kDirectionalMetadataScope =
    "114 directional-or-finite-grid-qualified paths: "
    "42 vector-component paths "
    "dipole_vector_{x,y,z}_{mean,m2,std,min,max,min_frame,max_frame} and "
    "surface_normal_{x,y,z}_{mean,m2,std,min,max,min_frame,max_frame}; "
    "plus 72 scalar paths for B in "
    "{half_shell_asymmetry,dipole_alignment}: "
    "B_{mean,m2,std,min,max,min_frame,max_frame}, "
    "B_{delta,abs_delta,delta_squared,dxdt}_"
    "{mean,m2,std,min,max,min_frame,max_frame}, and B_rms_delta; "
    "excludes exact rotation-invariant dipole_magnitude_"
    "{mean,m2,std,min,max,min_frame,max_frame} and, for E in "
    "{dipole_coherence,mean_net_dipole_eA,dipole_order_parameter,"
    "shell_count}, E_{mean,m2,std,min,max,min_frame,max_frame}, "
    "E_{delta,abs_delta,delta_squared,dxdt}_"
    "{mean,m2,std,min,max,min_frame,max_frame}, and E_rms_delta; also "
    "excludes "
    "n_frames_per_atom,delta_n_per_atom,dxdt_n_per_atom,"
    "source_attached_per_frame and group provenance";

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
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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


TEST(HydrationGeometryWelford, H5DirectionalMetadataZeroCountSynthetic) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto build = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto tp_owner = nmr::TrajectoryProtein::CreateForTesting(
        std::move(build.protein));
    ASSERT_NE(tp_owner, nullptr);
    auto& tp = *tp_owner;
    ASSERT_GT(tp.AtomCount(), 1u);

    auto tr = nmr::HydrationGeometryWelfordTrajectoryResult::Create(tp);
    tr->ForceSourcePresentForTesting();
    nmr::Trajectory traj({}, {}, {});
    auto conf = tp.TickConformation(tp.CanonicalConformation().Positions());
    tr->Compute(*conf, tp, traj, 7, 2.5);
    tp.MutableAtomAt(0).hydration_geometry_welford = {};
    tr->Finalize(tp, traj);

    const std::string h5_path = (fs::temp_directory_path() /
        ("hydration_geometry_welford_directional_metadata_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr->WriteH5Group(tp, file);
    }
    {
        HighFive::File file(h5_path, HighFive::File::ReadOnly);
        auto grp = file.getGroup(
            "/trajectory/hydration_geometry_welford");

        std::string dipole_frame, dipole_parity, dipole_law;
        std::string normal_frame, normal_parity, normal_law;
        std::string component_semantics, scalar_law, directional_scope;
        std::string irrep_scope, mean_law, component_law;
        std::string zero_count_validity;
        grp.getAttribute("dipole_vector_coordinate_frame")
            .read(dipole_frame);
        grp.getAttribute("dipole_vector_parity").read(dipole_parity);
        grp.getAttribute("dipole_vector_transformation").read(dipole_law);
        grp.getAttribute("surface_normal_coordinate_frame")
            .read(normal_frame);
        grp.getAttribute("surface_normal_parity").read(normal_parity);
        grp.getAttribute("surface_normal_transformation").read(normal_law);
        grp.getAttribute("source_component_semantics")
            .read(component_semantics);
        grp.getAttribute("scalar_statistic_transformation")
            .read(scalar_law);
        grp.getAttribute("directional_metadata_scope")
            .read(directional_scope);
        grp.getAttribute("irrep_metadata_scope").read(irrep_scope);
        grp.getAttribute("directional_mean_transformation").read(mean_law);
        grp.getAttribute("componentwise_statistic_transformation")
            .read(component_law);
        grp.getAttribute("zero_count_sentinel_validity")
            .read(zero_count_validity);
        EXPECT_EQ(dipole_frame, "conformation_cartesian_xyz");
        EXPECT_EQ(dipole_parity, "1o");
        EXPECT_EQ(dipole_law, "polar vector: v'=R v");
        EXPECT_EQ(normal_frame, "conformation_cartesian_xyz");
        EXPECT_EQ(normal_parity, "mixed");
        EXPECT_EQ(normal_law, kSurfaceNormalTransformation);
        EXPECT_EQ(component_semantics, kSourceComponentSemantics);
        EXPECT_EQ(scalar_law, kScalarStatisticTransformation);
        EXPECT_EQ(directional_scope, kDirectionalMetadataScope);
        EXPECT_EQ(irrep_scope, kIrrepMetadataScope);
        EXPECT_EQ(mean_law, kDirectionalMeanTransformation);
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

        auto expect_count_gated_zero = [&](const std::string& prefix) {
            std::vector<double> mean, m2, stddev, min, max;
            std::vector<std::size_t> min_frame, max_frame;
            grp.getDataSet(prefix + "_mean").read(mean);
            grp.getDataSet(prefix + "_m2").read(m2);
            grp.getDataSet(prefix + "_std").read(stddev);
            grp.getDataSet(prefix + "_min").read(min);
            grp.getDataSet(prefix + "_max").read(max);
            grp.getDataSet(prefix + "_min_frame").read(min_frame);
            grp.getDataSet(prefix + "_max_frame").read(max_frame);
            ASSERT_GE(mean.size(), 2u);
            EXPECT_TRUE(std::isnan(mean[0]));
            EXPECT_TRUE(std::isnan(m2[0]));
            EXPECT_TRUE(std::isnan(stddev[0]));
            EXPECT_TRUE(std::isinf(min[0]) && min[0] > 0.0);
            EXPECT_TRUE(std::isinf(max[0]) && max[0] < 0.0);
            EXPECT_EQ(min_frame[0], 0u);
            EXPECT_EQ(max_frame[0], 0u);

            EXPECT_DOUBLE_EQ(mean[1], 0.0);
            EXPECT_DOUBLE_EQ(m2[1], 0.0);
            EXPECT_DOUBLE_EQ(stddev[1], 0.0);
            EXPECT_DOUBLE_EQ(min[1], 0.0);
            EXPECT_DOUBLE_EQ(max[1], 0.0);
            EXPECT_EQ(min_frame[1], 7u);
            EXPECT_EQ(max_frame[1], 7u);
        };
        for (const std::string& prefix : {
                 "dipole_vector_x", "dipole_vector_y", "dipole_vector_z",
                 "surface_normal_x", "surface_normal_y", "surface_normal_z"}) {
            expect_count_gated_zero(prefix);
        }
    }
    fs::remove(h5_path);
}


TEST(HydrationGeometryWelford, H5RoundTrip) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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
                                     "dipole_coherence", "mean_net_dipole_eA",
                                     "dipole_order_parameter", "shell_count"}) {
        EXPECT_TRUE(grp.exist(base + "_delta_mean"))     << base;
        EXPECT_TRUE(grp.exist(base + "_abs_delta_mean")) << base;
        EXPECT_TRUE(grp.exist(base + "_dxdt_mean"))      << base;
        EXPECT_TRUE(grp.exist(base + "_rms_delta"))      << base;
    }
    EXPECT_TRUE(grp.exist("mean_net_dipole_eA_mean"));
    EXPECT_TRUE(grp.exist("dipole_order_parameter_mean"));
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));
    EXPECT_TRUE(grp.exist("dxdt_n_per_atom"));

    std::string polar_chans, alias_name, formula;
    grp.getAttribute("polarisation_signal_channels").read(polar_chans);
    grp.getAttribute("dipole_coherence_alias").read(alias_name);
    grp.getAttribute("dipole_order_parameter_formula").read(formula);
    EXPECT_EQ(polar_chans,
              "dipole_alignment,dipole_coherence,mean_net_dipole_eA,dipole_order_parameter,half_shell_asymmetry");
    EXPECT_EQ(alias_name, "mean_net_dipole_eA");
    EXPECT_EQ(formula, "|sum d_i| / sum |d_i|");

    std::string irrep_scope, mean_law, component_law, zero_count_validity;
    grp.getAttribute("irrep_metadata_scope").read(irrep_scope);
    grp.getAttribute("directional_mean_transformation").read(mean_law);
    grp.getAttribute("componentwise_statistic_transformation").read(component_law);
    grp.getAttribute("zero_count_sentinel_validity").read(zero_count_validity);
    EXPECT_EQ(irrep_scope, kIrrepMetadataScope);
    EXPECT_EQ(mean_law, kDirectionalMeanTransformation);
    EXPECT_EQ(component_law,
              "componentwise m2,std,min,max,min_frame,max_frame have no closed "
              "irrep transformation law");
    EXPECT_EQ(zero_count_validity,
              "when n_frames_per_atom=0, mean,m2,std are NaN and min=+inf,"
              "max=-inf,min_frame=0,max_frame=0 are invalid sentinels; "
              "n_frames_per_atom gates validity");

    auto read_units = [&](const std::string& dataset) {
        std::string units;
        grp.getDataSet(dataset).getAttribute("units").read(units);
        return units;
    };
    auto expect_double_alias = [&](const std::string& alias,
                                   const std::string& source) {
        ASSERT_TRUE(grp.exist(alias)) << alias;
        ASSERT_TRUE(grp.exist(source)) << source;
        std::vector<double> alias_values;
        std::vector<double> source_values;
        grp.getDataSet(alias).read(alias_values);
        grp.getDataSet(source).read(source_values);
        // Byte-identity is the spec's truth for an alias (the payload is a copy
        // of the source) and is NaN-safe, unlike EXPECT_EQ's element-wise ==,
        // which is false for the all-NaN dxdt/rms_delta channels even though the
        // bits are identical. Compare raw payload bytes.
        ASSERT_EQ(alias_values.size(), source_values.size())
            << alias << " vs " << source;
        EXPECT_EQ(0, std::memcmp(alias_values.data(), source_values.data(),
                                 alias_values.size() * sizeof(double)))
            << alias << " vs " << source << " (byte-identical payloads)";
        std::string alias_of;
        grp.getDataSet(alias).getAttribute("alias_of").read(alias_of);
        EXPECT_EQ(alias_of, source);
        EXPECT_EQ(read_units(alias), read_units(source));
    };
    auto expect_size_alias = [&](const std::string& alias,
                                 const std::string& source) {
        ASSERT_TRUE(grp.exist(alias)) << alias;
        ASSERT_TRUE(grp.exist(source)) << source;
        std::vector<std::size_t> alias_values;
        std::vector<std::size_t> source_values;
        grp.getDataSet(alias).read(alias_values);
        grp.getDataSet(source).read(source_values);
        EXPECT_EQ(alias_values, source_values) << alias << " vs " << source;
        std::string alias_of;
        grp.getDataSet(alias).getAttribute("alias_of").read(alias_of);
        EXPECT_EQ(alias_of, source);
        EXPECT_EQ(read_units(alias), read_units(source));
    };

    for (const std::string& prefix : {
             "mean_net_dipole_eA",
             "mean_net_dipole_eA_delta",
             "mean_net_dipole_eA_abs_delta",
             "mean_net_dipole_eA_delta_squared",
             "mean_net_dipole_eA_dxdt"}) {
        const std::string source = (prefix == "mean_net_dipole_eA")
            ? "dipole_coherence"
            : "dipole_coherence" + prefix.substr(std::string("mean_net_dipole_eA").size());
        for (const std::string& suffix : {"_mean", "_m2", "_std", "_min", "_max"}) {
            expect_double_alias(prefix + suffix, source + suffix);
        }
        expect_size_alias(prefix + "_min_frame", source + "_min_frame");
        expect_size_alias(prefix + "_max_frame", source + "_max_frame");
    }
    expect_double_alias("mean_net_dipole_eA_rms_delta",
                        "dipole_coherence_rms_delta");

    EXPECT_EQ(read_units("dipole_order_parameter_mean"), "dimensionless");
    EXPECT_EQ(read_units("dipole_order_parameter_m2"), "dimensionless^2");
    EXPECT_EQ(read_units("dipole_order_parameter_delta_mean"), "dimensionless");
    EXPECT_EQ(read_units("dipole_order_parameter_delta_m2"), "dimensionless^2");
    EXPECT_EQ(read_units("dipole_order_parameter_abs_delta_mean"), "dimensionless");
    EXPECT_EQ(read_units("dipole_order_parameter_abs_delta_m2"), "dimensionless^2");
    EXPECT_EQ(read_units("dipole_order_parameter_delta_squared_mean"), "dimensionless^2");
    EXPECT_EQ(read_units("dipole_order_parameter_delta_squared_m2"), "dimensionless^4");
    EXPECT_EQ(read_units("dipole_order_parameter_dxdt_mean"), "dimensionless/ps");
    EXPECT_EQ(read_units("dipole_order_parameter_dxdt_m2"), "dimensionless^2/ps^2");
    EXPECT_EQ(read_units("dipole_order_parameter_rms_delta"), "dimensionless");

    fs::remove(h5_path);
}


TEST(HydrationGeometryWelford, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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
        EXPECT_TRUE(std::isfinite(w.dipole_order_parameter.mean));
        EXPECT_TRUE(std::isfinite(w.shell_count.mean));
        EXPECT_EQ(w.n_frames, traj.FrameCount());
        EXPECT_EQ(w.dxdt_n, w.delta_n)
            << "atom " << i << ": dxdt_n != delta_n";
        // alignment is a cos angle; order parameter is dimensionless [0, 1].
        EXPECT_GE(w.dipole_alignment.mean, -1.001);
        EXPECT_LE(w.dipole_alignment.mean,  1.001);
        EXPECT_GE(w.dipole_order_parameter.mean, -0.001);
        EXPECT_LE(w.dipole_order_parameter.mean,  1.001);
        if (std::abs(w.dipole_alignment.mean) > 1e-12) ++populated;
        max_abs_alignment = std::max(max_abs_alignment, std::abs(w.dipole_alignment.mean));
    }
    EXPECT_GT(populated, 0u) << "dipole_alignment all zero on solvated protein";
    std::cout << "HydrationGeometryWelford: populated=" << populated << "/"
              << tp.AtomCount() << " max|alignment|=" << max_abs_alignment << "\n";
}
