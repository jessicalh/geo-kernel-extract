//
// test_apbs_efield_time_series: discipline + integration for
// ApbsEfieldTimeSeriesTrajectoryResult. Vec3 DenseBuffer pattern.
//

#include "ApbsFieldResult.h"
#include "ChargeAssignmentResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DenseBuffer.h"
#include "ApbsEfieldTimeSeriesTrajectoryResult.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
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
}  // namespace


TEST(ApbsEfieldTimeSeries, SyntheticFourFrames) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto build = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto tp_owner = nmr::TrajectoryProtein::CreateForTesting(
        std::move(build.protein));
    ASSERT_NE(tp_owner, nullptr);
    auto& tp = *tp_owner;
    const size_t N = tp.AtomCount();
    auto tr = nmr::ApbsEfieldTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj({}, {}, {});

    constexpr size_t kFrames = 4;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        conf->ForceAttachResultForTesting(
            std::make_unique<nmr::ApbsFieldResult>());
        for (size_t i = 0; i < N; ++i) {
            conf->MutableAtomAt(i).apbs_efield =
                nmr::Vec3(static_cast<double>(i) + t * 1.0,
                          static_cast<double>(i) + t * 1.0 + 0.1,
                          static_cast<double>(i) + t * 1.0 + 0.2);
            conf->MutableAtomAt(i).apbs_efield_clamp_mask =
                (i == 0 && t == 1) ? 1u : 0u;
            conf->MutableAtomAt(i).apbs_efield_clamp_scale =
                (i == 0 && t == 1) ? 0.5 : 1.0;
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    auto* buf = tp.GetDenseBuffer<nmr::Vec3>(std::type_index(typeid(
        nmr::ApbsEfieldTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->StridePerAtom(), kFrames);
    for (size_t i : {size_t(0), N / 2, N - 1})
        for (size_t t = 0; t < kFrames; ++t) {
            const nmr::Vec3& v = buf->At(i, t);
            EXPECT_DOUBLE_EQ(v.x(), static_cast<double>(i) + t * 1.0);
            EXPECT_DOUBLE_EQ(v.y(), static_cast<double>(i) + t * 1.0 + 0.1);
            EXPECT_DOUBLE_EQ(v.z(), static_cast<double>(i) + t * 1.0 + 0.2);
        }

    const std::string h5_path = (fs::temp_directory_path() /
        ("apbs_efield_ts_unit_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/apbs_efield_time_series");
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[0], N);
    EXPECT_EQ(dims[1], kFrames);
    EXPECT_EQ(dims[2], 3u);
    std::string parity, layout, directional_scope;
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("irrep_layout").read(layout);
    grp.getAttribute("directional_metadata_scope").read(directional_scope);
    EXPECT_EQ(parity, "1o");
    EXPECT_EQ(layout, "x,y,z");
    EXPECT_EQ(directional_scope,
        "irrep_layout,normalization,parity,units describe xyz only; "
        "clamp_* and apbs_grid_*_per_frame datasets carry per-dataset "
        "scalar/lab-axis contracts");
    std::string xyz_frame, xyz_transformation, xyz_parity;
    ds.getAttribute("coordinate_frame").read(xyz_frame);
    ds.getAttribute("transformation").read(xyz_transformation);
    ds.getAttribute("parity").read(xyz_parity);
    EXPECT_EQ(xyz_frame, "conformation_cartesian_xyz");
    EXPECT_EQ(xyz_transformation,
              "continuum polar_vector: v'=R v; translation invariant. The "
              "live axis-aligned finite-difference APBS solve has no exact "
              "O(3) law; transformed production reruns use the recorded "
              "1.8e-2 V/Angstrom absolute + 5e-2 relative finite-grid envelope");
    EXPECT_EQ(xyz_parity, "1o");

    auto source_ds = grp.getDataSet("source_attached_per_frame");
    std::vector<std::uint8_t> source_mask(kFrames);
    source_ds.read(source_mask.data());
    for (const auto value : source_mask) EXPECT_EQ(value, 1u);

    auto clamp_mask_ds = grp.getDataSet("clamp_mask");
    auto clamp_scale_ds = grp.getDataSet("clamp_scale");
    EXPECT_EQ(clamp_mask_ds.getSpace().getDimensions(),
              (std::vector<std::size_t>{N, kFrames}));
    EXPECT_EQ(clamp_scale_ds.getSpace().getDimensions(),
              (std::vector<std::size_t>{N, kFrames}));
    std::vector<std::uint8_t> clamp_mask(N * kFrames);
    std::vector<double> clamp_scale(N * kFrames);
    clamp_mask_ds.read(clamp_mask.data());
    clamp_scale_ds.read(clamp_scale.data());
    EXPECT_EQ(clamp_mask[1], 1u);
    EXPECT_DOUBLE_EQ(clamp_scale[1], 0.5);
    EXPECT_EQ(clamp_mask[0], 0u);
    EXPECT_DOUBLE_EQ(clamp_scale[0], 1.0);

    for (const std::string& name : {
            "apbs_grid_dims_per_frame", "apbs_grid_lengths_A_per_frame",
            "apbs_grid_origin_A_per_frame", "apbs_grid_spacing_A_per_frame"}) {
        EXPECT_EQ(grp.getDataSet(name).getSpace().getDimensions(),
                  (std::vector<std::size_t>{kFrames, 3u}));
    }
    const auto expect_grid_contract = [&](const std::string& name,
            const std::string& expected_units,
            const std::string& expected_parity,
            const std::string& expected_transformation) {
        auto grid = grp.getDataSet(name);
        std::string order, grid_units, frame, grid_parity, transformation;
        grid.getAttribute("component_order").read(order);
        grid.getAttribute("units").read(grid_units);
        grid.getAttribute("coordinate_frame").read(frame);
        grid.getAttribute("parity").read(grid_parity);
        grid.getAttribute("transformation").read(transformation);
        EXPECT_EQ(order, "x,y,z");
        EXPECT_EQ(grid_units, expected_units);
        EXPECT_EQ(frame, "apbs_lab_axis_aligned_grid_xyz");
        EXPECT_EQ(grid_parity, expected_parity);
        EXPECT_EQ(transformation, expected_transformation);
    };
    expect_grid_contract("apbs_grid_dims_per_frame", "grid_points", "even",
        "configured cubic grid point counts (n,n,n); invariant under "
        "proper/improper orthogonal transforms and translations");
    expect_grid_contract("apbs_grid_lengths_A_per_frame", "Angstrom", "mixed",
        "lab-axis grid extents; translation invariant; no closed O(3) "
        "transformation law under arbitrary orthogonal transforms");
    expect_grid_contract("apbs_grid_origin_A_per_frame", "Angstrom", "mixed",
        "lab-axis grid origin; o'=o+t under pure translations; no "
        "affine-position/O(3) law under arbitrary orthogonal transforms "
        "because grid axes remain lab-fixed");
    expect_grid_contract("apbs_grid_spacing_A_per_frame", "Angstrom", "mixed",
        "lab-axis grid step lengths; translation invariant; no closed O(3) "
        "transformation law under arbitrary orthogonal transforms");
    std::string field_quantity, grid_mode, reference_subtracts,
                clamp_units, rank3_policy;
    grp.getAttribute("field_quantity").read(field_quantity);
    grp.getAttribute("apbs_grid_mode").read(grid_mode);
    grp.getAttribute("reference_subtracts").read(reference_subtracts);
    grp.getAttribute("efield_clamp_units").read(clamp_units);
    grp.getAttribute("rank3_policy").read(rank3_policy);
    EXPECT_EQ(field_quantity,
              "reaction_field_total_minus_homogeneous_vacuum_reference");
    EXPECT_EQ(grid_mode, "single_manual");
    EXPECT_EQ(reference_subtracts, "all_solute_charges");
    EXPECT_EQ(clamp_units, "V/Angstrom");
    EXPECT_EQ(rank3_policy, "not_emitted_no_local_frame");
    int max_rank = 0;
    bool higher_present = true;
    grp.getAttribute("max_potential_derivative_rank").read(max_rank);
    grp.getAttribute("higher_derivatives_present").read(higher_present);
    EXPECT_EQ(max_rank, 2);
    EXPECT_FALSE(higher_present);
    fs::remove(h5_path);
}


TEST(ApbsEfieldTimeSeries, SourceAbsentNanFillMaskAndMetadata) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    const size_t N = tp.AtomCount();
    auto tr = nmr::ApbsEfieldTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 2;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic absent APBS frame");
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    auto* buffer = tp.GetDenseBuffer<nmr::Vec3>(std::type_index(typeid(
        nmr::ApbsEfieldTimeSeriesTrajectoryResult)));
    ASSERT_NE(buffer, nullptr);
    for (size_t i : {size_t(0), N / 2, N - 1})
        for (size_t t = 0; t < kFrames; ++t)
            EXPECT_FALSE(buffer->At(i, t).allFinite());

    const std::string h5_path = (fs::temp_directory_path() /
        ("apbs_efield_ts_absent_" + std::to_string(::getpid()) + ".h5"))
        .string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/apbs_efield_time_series");

    std::vector<std::uint8_t> source(kFrames);
    grp.getDataSet("source_attached_per_frame").read(source.data());
    for (const auto value : source) EXPECT_EQ(value, 0u);

    std::vector<std::uint8_t> clamp_mask(N * kFrames);
    std::vector<double> clamp_scale(N * kFrames);
    grp.getDataSet("clamp_mask").read(clamp_mask.data());
    grp.getDataSet("clamp_scale").read(clamp_scale.data());
    for (const auto value : clamp_mask) EXPECT_EQ(value, 0u);
    for (const auto value : clamp_scale) EXPECT_TRUE(std::isnan(value));

    std::vector<std::uint64_t> grid_dims(kFrames * 3);
    std::vector<double> grid_lengths(kFrames * 3);
    grp.getDataSet("apbs_grid_dims_per_frame").read(grid_dims.data());
    grp.getDataSet("apbs_grid_lengths_A_per_frame").read(
        grid_lengths.data());
    for (const auto value : grid_dims) EXPECT_EQ(value, 0u);
    for (const auto value : grid_lengths) EXPECT_TRUE(std::isnan(value));
    fs::remove(h5_path);
}


TEST(ApbsEfieldTimeSeries, Frame0Semantics) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = false;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::ChargeAssignmentResult));
    config.RequireConformationResult(typeid(nmr::ApbsFieldResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::ApbsEfieldTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_EQ(traj.FrameCount(), 1u);
}


TEST(ApbsEfieldTimeSeries, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = false;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::ChargeAssignmentResult));
    config.RequireConformationResult(typeid(nmr::ApbsFieldResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::ApbsEfieldTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);  // single-frame keeps the test cheap

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto* buf_first = tp.GetDenseBuffer<nmr::Vec3>(std::type_index(typeid(
        nmr::ApbsEfieldTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_first, nullptr);
    const std::size_t N_first = buf_first->AtomCount();
    const std::size_t T_first = buf_first->StridePerAtom();

    auto& tr = tp.Result<nmr::ApbsEfieldTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);
    auto* buf_second = tp.GetDenseBuffer<nmr::Vec3>(std::type_index(typeid(
        nmr::ApbsEfieldTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_second, nullptr);
    EXPECT_EQ(buf_second->AtomCount(), N_first);
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
}


TEST(ApbsEfieldTimeSeries, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = false;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::ChargeAssignmentResult));
    config.RequireConformationResult(typeid(nmr::ApbsFieldResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::ApbsEfieldTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto* buf = tp.GetDenseBuffer<nmr::Vec3>(std::type_index(typeid(
        nmr::ApbsEfieldTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);

    double max_mag = 0.0;
    const double clamp_threshold = nmr::CalculatorConfig::Get(
        "efield_magnitude_sanity_clamp");
    for (size_t i = 0; i < buf->AtomCount(); ++i) {
        for (size_t t = 0; t < buf->StridePerAtom(); ++t) {
            const nmr::Vec3& v = buf->At(i, t);
            EXPECT_TRUE(std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z()));
            // Canonical reaction E is clamped in V/A by calculator config.
            EXPECT_LE(v.norm(), clamp_threshold + 1e-12);
            max_mag = std::max(max_mag, v.norm());
        }
    }
    EXPECT_GT(max_mag, 0.01) << "APBS E-field all near zero — calc not firing";
    std::cout << "ApbsEfieldTimeSeries max|E|=" << max_mag << " V/A\n";

    auto& tr_result = tp.Result<
        nmr::ApbsEfieldTimeSeriesTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("apbs_efield_ts_integration_" + std::to_string(::getpid()) +
         ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr_result.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/apbs_efield_time_series");
    double temperature = 0.0;
    double thermal_voltage = 0.0;
    double padding = 0.0;
    double min_dim = 0.0;
    double emitted_clamp = 0.0;
    grp.getAttribute("apbs_temperature_K").read(temperature);
    grp.getAttribute("apbs_thermal_voltage_V").read(thermal_voltage);
    grp.getAttribute("apbs_manual_grid_padding_A").read(padding);
    grp.getAttribute("apbs_manual_grid_min_dim_A").read(min_dim);
    grp.getAttribute("efield_clamp_threshold").read(emitted_clamp);
    EXPECT_DOUBLE_EQ(temperature,
        nmr::CalculatorConfig::Get("apbs_temperature_K"));
    EXPECT_DOUBLE_EQ(thermal_voltage,
        nmr::ApbsThermalVoltage(temperature));
    EXPECT_DOUBLE_EQ(padding,
        nmr::CalculatorConfig::Get("apbs_manual_grid_padding_A"));
    EXPECT_DOUBLE_EQ(min_dim,
        nmr::CalculatorConfig::Get("apbs_manual_grid_min_dim_A"));
    EXPECT_DOUBLE_EQ(emitted_clamp, clamp_threshold);
    EXPECT_EQ(grp.getDataSet("apbs_grid_dims_per_frame")
                  .getSpace().getDimensions(),
              (std::vector<std::size_t>{buf->StridePerAtom(), 3u}));
    fs::remove(h5_path);
}
