//
// test_water_field_time_series: discipline + integration for
// WaterFieldTimeSeriesTrajectoryResult. Per-atom 6-channel TR (Vec3 + Vec3
// + SphericalTensor + SphericalTensor + int + int) over WaterFieldResult.
//

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "GeometryResult.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"
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
#include <cstdint>
#include <filesystem>
#include <fstream>
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

nmr::ProteinConformation& BuildOneAtomWaterFixture(nmr::Protein& protein) {
    nmr::Residue residue;
    residue.type = nmr::AminoAcid::ALA;
    residue.sequence_number = 1;
    const size_t ri = protein.AddResidue(std::move(residue));
    auto atom = nmr::Atom::Create(nmr::Element::C);
    atom->residue_index = ri;
    const size_t ai = protein.AddAtom(std::move(atom));
    protein.MutableResidueAt(ri).atom_indices.push_back(ai);
    return protein.AddConformation({nmr::Vec3(0.0, 0.0, 0.0)},
                                   "one atom water field");
}

void RunWaterClampIndependenceInChild() {
    const fs::path toml = fs::temp_directory_path() /
        ("water_clamp_" + std::to_string(::getpid()) + ".toml");
    {
        std::ofstream out(toml);
        out << "efield_magnitude_sanity_clamp = 1.0\n";
    }
    nmr::CalculatorConfig::Load(toml.string());
    nmr::Protein protein;
    auto& conf = BuildOneAtomWaterFixture(protein);

    nmr::SolventEnvironment solvent;
    nmr::WaterMolecule first;
    first.O_pos = nmr::Vec3(2.0, 0.0, 0.0);
    first.H1_pos = first.O_pos;
    first.H2_pos = first.O_pos;
    first.O_charge = 0.001;
    first.H_charge = 0.0;
    nmr::WaterMolecule outer;
    outer.O_pos = nmr::Vec3(4.0, 0.0, 0.0);
    outer.H1_pos = outer.O_pos;
    outer.H2_pos = outer.O_pos;
    outer.O_charge = 2.0;
    outer.H_charge = 0.0;
    solvent.waters = {first, outer};
    solvent.water_O_positions = {first.O_pos, outer.O_pos};

    auto result = nmr::WaterFieldResult::Compute(conf, solvent);
    const auto& atom = conf.AtomAt(0);
    const nmr::Vec3 expected_first =
        (nmr::COULOMB_KE * 0.001 / 8.0) * nmr::Vec3(-2.0, 0.0, 0.0);
    const bool ok = result != nullptr
        && atom.water_efield_clamp_mask == 1
        && atom.water_efield_clamp_scale < 1.0
        && std::abs(atom.water_efield.norm() - 1.0) < 1e-9
        && atom.water_efield_first_clamp_mask == 0
        && atom.water_efield_first_clamp_scale == 1.0
        && (atom.water_efield_first - expected_first).norm() < 1e-12;
    fs::remove(toml);
    _exit(ok ? 0 : 1);
}

}  // namespace


// ============================================================================
// SYNTHETIC: source-absent path — verifies all-absent → group skipped.
// R3 codex F5 2026-05-18.
// ============================================================================

TEST(WaterFieldTimeSeries, SyntheticAllAbsentSkipsGroup) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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

TEST(WaterFieldTimeSeries, ClampIndependence) {
    EXPECT_EXIT({ RunWaterClampIndependenceInChild(); },
                ::testing::ExitedWithCode(0), "");
}

TEST(WaterFieldResult, UsesTargetMinusSourceCoulombConvention) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::Protein protein;
    auto& conf = BuildOneAtomWaterFixture(protein);

    nmr::SolventEnvironment solvent;
    nmr::WaterMolecule water;
    water.O_pos = nmr::Vec3(2.0, 0.0, 0.0);
    water.H1_pos = water.O_pos;
    water.H2_pos = water.O_pos;
    water.O_charge = 0.5;
    water.H_charge = 0.0;
    solvent.waters = {water};
    solvent.water_O_positions = {water.O_pos};

    ASSERT_NE(nmr::WaterFieldResult::Compute(conf, solvent), nullptr);
    const auto& atom = conf.AtomAt(0);
    const nmr::Vec3 r = conf.PositionAt(0) - water.O_pos;
    const double rmag = r.norm();
    const double r3 = rmag * rmag * rmag;
    const double r5 = r3 * rmag * rmag;
    const nmr::Vec3 expected_E =
        nmr::COULOMB_KE * water.O_charge * r / r3;
    const nmr::Mat3 expected_EFG = nmr::COULOMB_KE * water.O_charge *
        (3.0 * r * r.transpose() / r5 - nmr::Mat3::Identity() / r3);

    EXPECT_LT((atom.water_efield - expected_E).norm(), 1e-12);
    EXPECT_LT((atom.water_efg - expected_EFG).norm(), 1e-12);
    EXPECT_LT(std::abs(atom.water_efg.trace()), 1e-12);
}


TEST(WaterFieldTimeSeries, ClampH5RoundTripMixedSource) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::WaterFieldTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    auto present = std::make_unique<nmr::ProteinConformation>(
        &tp.ProteinRef(), positions, "synthetic present");
    present->ForceAttachResultForTesting(std::make_unique<nmr::WaterFieldResult>());
    tr->Compute(*present, tp, traj, 0, 0.0);

    auto absent = std::make_unique<nmr::ProteinConformation>(
        &tp.ProteinRef(), positions, "synthetic absent");
    tr->Compute(*absent, tp, traj, 1, 1.0);
    tr->Finalize(tp, traj);

    const std::string h5_path = (fs::temp_directory_path() /
        ("water_ts_clamp_mixed_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/water_field_time_series"));
    auto grp = reopen.getGroup("/trajectory/water_field_time_series");

    for (const std::string& name : {"efield_clamp_mask", "efield_clamp_scale",
                                    "efield_first_clamp_mask",
                                    "efield_first_clamp_scale"}) {
        ASSERT_TRUE(grp.exist(name)) << name;
        const auto dims = grp.getDataSet(name).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], N);
        EXPECT_EQ(dims[1], 2u);
    }
    std::uint8_t absent_sentinel = 0;
    grp.getDataSet("efield_clamp_mask")
        .getAttribute("absent_sentinel").read(absent_sentinel);
    EXPECT_EQ(absent_sentinel, 255u);

    std::vector<std::vector<std::uint8_t>> total_mask;
    std::vector<std::vector<std::uint8_t>> first_mask;
    std::vector<std::vector<double>> total_scale;
    std::vector<std::vector<double>> first_scale;
    grp.getDataSet("efield_clamp_mask").read(total_mask);
    grp.getDataSet("efield_first_clamp_mask").read(first_mask);
    grp.getDataSet("efield_clamp_scale").read(total_scale);
    grp.getDataSet("efield_first_clamp_scale").read(first_scale);
    ASSERT_EQ(total_mask.size(), N);
    for (std::size_t i = 0; i < N; ++i) {
        EXPECT_EQ(total_mask[i][0], 0u);
        EXPECT_EQ(first_mask[i][0], 0u);
        EXPECT_DOUBLE_EQ(total_scale[i][0], 1.0);
        EXPECT_DOUBLE_EQ(first_scale[i][0], 1.0);
        EXPECT_EQ(total_mask[i][1], 255u);
        EXPECT_EQ(first_mask[i][1], 255u);
        EXPECT_TRUE(std::isnan(total_scale[i][1]));
        EXPECT_TRUE(std::isnan(first_scale[i][1]));
    }

    fs::remove(h5_path);
}


TEST(WaterFieldTimeSeries, Frame0Semantics) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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

    std::string efield_units, efg_basis, efg_order, efg_parity;
    std::string efg_export, count_units, legacy_layout;
    bool legacy_deprecated = false;
    grp.getAttribute("efield_units").read(efield_units);
    grp.getAttribute("efg_t2_basis").read(efg_basis);
    grp.getAttribute("efg_t2_component_order").read(efg_order);
    grp.getAttribute("efg_t2_parity").read(efg_parity);
    grp.getAttribute("efg_e3nn_export").read(efg_export);
    grp.getAttribute("legacy_irrep_attrs_deprecated").read(legacy_deprecated);
    grp.getAttribute("efg_irrep_layout").read(legacy_layout);
    grp.getAttribute("count_units").read(count_units);
    EXPECT_EQ(efield_units, "V/Angstrom");
    // Water EFG: T0 (trace) and T1 (antisym pseudovector) are both
    // structural zeros — only T2 (symmetric-traceless, 5 components)
    // is emitted.
    EXPECT_EQ(efg_basis, "project_native_t2_isometric_real_tesseral_v1");
    EXPECT_EQ(efg_order, "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");
    EXPECT_EQ(efg_parity, "even");
    EXPECT_NE(efg_export.find("project_t2_to_e3nn"), std::string::npos);
    EXPECT_TRUE(legacy_deprecated);
    EXPECT_EQ(legacy_layout, "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");
    EXPECT_EQ(count_units, "dimensionless");

    bool t0_zero = false, t1_zero = false;
    grp.getAttribute("efg_t0_structural_zero").read(t0_zero);
    grp.getAttribute("efg_t1_structural_zero").read(t1_zero);
    EXPECT_TRUE(t0_zero);
    EXPECT_TRUE(t1_zero);

    int max_rank = 0;
    bool higher = true;
    std::string rank3_policy;
    grp.getAttribute("max_potential_derivative_rank").read(max_rank);
    grp.getAttribute("higher_derivatives_present").read(higher);
    grp.getAttribute("rank3_policy").read(rank3_policy);
    EXPECT_EQ(max_rank, 2);
    EXPECT_FALSE(higher);
    EXPECT_EQ(rank3_policy, "not_emitted_no_local_frame");

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
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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
