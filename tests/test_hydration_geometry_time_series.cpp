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
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "SasaResult.h"
#include "Session.h"
#include "SolventEnvironment.h"
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
#include <unistd.h>
#include <utility>
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
constexpr const char* kScalarTransformation =
    "half_shell_asymmetry,dipole_alignment: continuum rotation invariants "
    "inheriting the finite lab-fixed Fibonacci surface-normal approximation; "
    "dipole_coherence,mean_net_dipole_eA,dipole_order_parameter,"
    "first_shell_count: exact rotation invariants";
constexpr const char* kDirectionalMetadataScope =
    "dipole_vector,surface_normal,half_shell_asymmetry,dipole_alignment,"
    "dipole_coherence,mean_net_dipole_eA,dipole_order_parameter,"
    "first_shell_count; excludes frame_indices,frame_times,"
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
        return nmr::HydrationGeometryTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


TEST(HydrationGeometryTimeSeries,
     H5DirectionalMetadataIsExactWithoutFleetFixture) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();

    auto build = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto tp_owner = nmr::TrajectoryProtein::CreateForTesting(
        std::move(build.protein));
    ASSERT_NE(tp_owner, nullptr);
    auto& tp = *tp_owner;
    ASSERT_GT(tp.AtomCount(), 0u);

    auto tr = nmr::HydrationGeometryTimeSeriesTrajectoryResult::Create(tp);
    tr->ForceSourcePresentForTesting();
    nmr::Trajectory traj({}, {}, {});
    std::vector<nmr::Vec3> positions(tp.AtomCount(), nmr::Vec3::Zero());
    auto conf = std::make_unique<nmr::ProteinConformation>(
        &tp.ProteinRef(), positions,
        "synthetic hydration geometry metadata frame");
    auto& atom = conf->MutableAtomAt(0);
    atom.water_dipole_vector = nmr::Vec3(1.25, -2.5, 3.75);
    atom.water_surface_normal = nmr::Vec3(-0.25, 0.5, 0.75);
    atom.sasa_half_shell_asymmetry = 0.6;
    atom.sasa_dipole_alignment = -0.4;
    atom.sasa_dipole_coherence = 0.8;
    atom.sasa_dipole_order_parameter = 0.3;
    atom.sasa_first_shell_count = 7;

    tr->Compute(*conf, tp, traj, 0, 0.0);
    tr->Finalize(tp, traj);

    const std::string h5_path = (fs::temp_directory_path() /
        ("hydration_geometry_ts_directional_metadata_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr->WriteH5Group(tp, file);
    }

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup(
        "/trajectory/hydration_geometry_time_series");
    std::string dipole_layout;
    std::string dipole_frame;
    std::string dipole_parity;
    std::string dipole_transformation;
    std::string normal_layout;
    std::string normal_frame;
    std::string normal_parity;
    std::string normal_transformation;
    std::string scalar_transformation;
    std::string metadata_scope;
    grp.getAttribute("dipole_vector_layout").read(dipole_layout);
    grp.getAttribute("dipole_vector_coordinate_frame").read(dipole_frame);
    grp.getAttribute("dipole_vector_parity").read(dipole_parity);
    grp.getAttribute("dipole_vector_transformation").read(
        dipole_transformation);
    grp.getAttribute("surface_normal_layout").read(normal_layout);
    grp.getAttribute("surface_normal_coordinate_frame").read(normal_frame);
    grp.getAttribute("surface_normal_parity").read(normal_parity);
    grp.getAttribute("surface_normal_transformation").read(
        normal_transformation);
    grp.getAttribute("scalar_transformation").read(scalar_transformation);
    grp.getAttribute("directional_metadata_scope").read(metadata_scope);

    EXPECT_EQ(dipole_layout, "x,y,z");
    EXPECT_EQ(dipole_frame, "conformation_cartesian_xyz");
    EXPECT_EQ(dipole_parity, "1o");
    EXPECT_EQ(dipole_transformation, "polar vector: v'=R v");
    EXPECT_EQ(normal_layout, "x,y,z");
    EXPECT_EQ(normal_frame, "conformation_cartesian_xyz");
    EXPECT_EQ(normal_parity, "mixed");
    EXPECT_EQ(normal_transformation, kSurfaceNormalTransformation);
    EXPECT_EQ(scalar_transformation, kScalarTransformation);
    EXPECT_EQ(metadata_scope, kDirectionalMetadataScope);

    std::vector<double> dipole(tp.AtomCount() * 3u);
    std::vector<double> normal(tp.AtomCount() * 3u);
    grp.getDataSet("dipole_vector").read(dipole.data());
    grp.getDataSet("surface_normal").read(normal.data());
    EXPECT_DOUBLE_EQ(dipole[0], 1.25);
    EXPECT_DOUBLE_EQ(dipole[1], -2.5);
    EXPECT_DOUBLE_EQ(dipole[2], 3.75);
    EXPECT_DOUBLE_EQ(normal[0], -0.25);
    EXPECT_DOUBLE_EQ(normal[1], 0.5);
    EXPECT_DOUBLE_EQ(normal[2], 0.75);

    fs::remove(h5_path);
}

TEST(HydrationGeometryResult, ConstructedAlignedWaterDipoles) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();

    nmr::Protein protein;
    nmr::Residue residue;
    residue.type = nmr::AminoAcid::ALA;
    residue.sequence_number = 1;
    const size_t ri = protein.AddResidue(std::move(residue));
    auto atom = nmr::Atom::Create(nmr::Element::C);
    atom->residue_index = ri;
    const size_t ai = protein.AddAtom(std::move(atom));
    protein.MutableResidueAt(ri).atom_indices.push_back(ai);
    auto& conf = protein.AddConformation({nmr::Vec3(0.0, 0.0, 0.0)},
                                         "aligned water dipoles");
    conf.MutableAtomAt(0).sasa_normal = nmr::Vec3(1.0, 0.0, 0.0);

    nmr::SolventEnvironment solvent;
    for (const nmr::Vec3& o : {nmr::Vec3(1.0, 0.0, 0.0),
                               nmr::Vec3(1.0, 1.0, 0.0)}) {
        nmr::WaterMolecule water;
        water.O_pos = o;
        water.H1_pos = o + nmr::Vec3(1.0, 0.0, 0.0);
        water.H2_pos = o + nmr::Vec3(1.0, 0.0, 0.0);
        water.O_charge = -0.834;
        water.H_charge = 0.417;
        solvent.waters.push_back(water);
        solvent.water_O_positions.push_back(water.O_pos);
    }

    auto result = nmr::HydrationGeometryResult::Compute(conf, solvent);
    ASSERT_NE(result, nullptr);
    const auto& out = conf.AtomAt(0);
    EXPECT_EQ(out.sasa_first_shell_count, 2);
    EXPECT_NEAR(out.sasa_dipole_coherence, 0.834, 1e-12);
    EXPECT_NEAR(out.sasa_dipole_order_parameter, 1.0, 1e-12);
}

TEST(HydrationGeometryResult, ConstructedCancellingWaterDipoles) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();

    nmr::Protein protein;
    nmr::Residue residue;
    residue.type = nmr::AminoAcid::ALA;
    residue.sequence_number = 1;
    const size_t ri = protein.AddResidue(std::move(residue));
    auto atom = nmr::Atom::Create(nmr::Element::C);
    atom->residue_index = ri;
    const size_t ai = protein.AddAtom(std::move(atom));
    protein.MutableResidueAt(ri).atom_indices.push_back(ai);
    auto& conf = protein.AddConformation({nmr::Vec3(0.0, 0.0, 0.0)},
                                         "cancelling water dipoles");
    conf.MutableAtomAt(0).sasa_normal = nmr::Vec3(1.0, 0.0, 0.0);

    nmr::SolventEnvironment solvent;
    const std::vector<std::pair<nmr::Vec3, nmr::Vec3>> waters = {
        {nmr::Vec3(1.0, 0.0, 0.0), nmr::Vec3(1.0, 0.0, 0.0)},
        {nmr::Vec3(1.0, 1.0, 0.0), nmr::Vec3(-1.0, 0.0, 0.0)},
    };
    for (const auto& entry : waters) {
        nmr::WaterMolecule water;
        water.O_pos = entry.first;
        water.H1_pos = entry.first + entry.second;
        water.H2_pos = entry.first + entry.second;
        water.O_charge = -0.834;
        water.H_charge = 0.417;
        solvent.waters.push_back(water);
        solvent.water_O_positions.push_back(water.O_pos);
    }

    auto result = nmr::HydrationGeometryResult::Compute(conf, solvent);
    ASSERT_NE(result, nullptr);
    const auto& out = conf.AtomAt(0);
    EXPECT_EQ(out.sasa_first_shell_count, 2);
    EXPECT_NEAR(out.sasa_dipole_coherence, 0.0, 1e-12);
    EXPECT_NEAR(out.sasa_dipole_order_parameter, 0.0, 1e-12);
}


TEST(HydrationGeometryTimeSeries, Frame0Semantics) {
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
    const auto& tr = tp.Result<nmr::HydrationGeometryTimeSeriesTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);
}


TEST(HydrationGeometryTimeSeries, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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

    const auto& tr = tp.Result<nmr::HydrationGeometryTimeSeriesTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("hyd_geo_ts_h5_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/hydration_geometry_time_series"));
    auto grp = reopen.getGroup("/trajectory/hydration_geometry_time_series");

    std::string layout, polar_chans, alias, formula;
    grp.getAttribute("dipole_vector_layout").read(layout);
    grp.getAttribute("polarisation_signal_channels").read(polar_chans);
    grp.getAttribute("dipole_coherence_alias").read(alias);
    grp.getAttribute("dipole_order_parameter_formula").read(formula);
    EXPECT_EQ(layout, "x,y,z");
    EXPECT_EQ(polar_chans, "dipole_alignment,dipole_coherence,mean_net_dipole_eA,dipole_order_parameter,half_shell_asymmetry");
    EXPECT_EQ(alias, "mean_net_dipole_eA");
    EXPECT_EQ(formula, "|sum d_i| / sum |d_i|");

    const std::size_t N = tp.AtomCount();
    const std::size_t T = tr.NumFrames();
    const auto dipole_dims = grp.getDataSet("dipole_vector").getSpace().getDimensions();
    ASSERT_EQ(dipole_dims.size(), 3u);
    EXPECT_EQ(dipole_dims[0], N); EXPECT_EQ(dipole_dims[1], T); EXPECT_EQ(dipole_dims[2], 3u);

    EXPECT_TRUE(grp.exist("dipole_alignment"));
    EXPECT_TRUE(grp.exist("dipole_coherence"));
    EXPECT_TRUE(grp.exist("mean_net_dipole_eA"));
    EXPECT_TRUE(grp.exist("dipole_order_parameter"));
    EXPECT_TRUE(grp.exist("half_shell_asymmetry"));
    EXPECT_TRUE(grp.exist("first_shell_count"));
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));

    std::vector<std::vector<double>> coherence;
    std::vector<std::vector<double>> mean_net;
    grp.getDataSet("dipole_coherence").read(coherence);
    grp.getDataSet("mean_net_dipole_eA").read(mean_net);
    ASSERT_EQ(coherence.size(), mean_net.size());
    for (std::size_t i = 0; i < coherence.size(); ++i) {
        ASSERT_EQ(coherence[i].size(), mean_net[i].size());
        for (std::size_t t = 0; t < coherence[i].size(); ++t) {
            EXPECT_DOUBLE_EQ(coherence[i][t], mean_net[i][t]);
        }
    }

    fs::remove(h5_path);
}


TEST(HydrationGeometryTimeSeries, Integration1P9J) {
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
    std::vector<std::vector<double>> order;
    grp.getDataSet("dipole_alignment").read(alignment);
    grp.getDataSet("dipole_order_parameter").read(order);
    ASSERT_EQ(alignment.size(), N);
    ASSERT_EQ(order.size(), N);

    std::size_t populated = 0;
    double max_abs_align = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            EXPECT_TRUE(std::isfinite(alignment[i][t]));
            // alignment is cos angle ∈ [-1, 1]
            EXPECT_GE(alignment[i][t], -1.001);
            EXPECT_LE(alignment[i][t],  1.001);
            EXPECT_TRUE(std::isfinite(order[i][t]));
            EXPECT_GE(order[i][t], -0.001);
            EXPECT_LE(order[i][t],  1.001);
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
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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
