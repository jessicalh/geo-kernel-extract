//
// test_hydration_shell_time_series: discipline + integration for
// HydrationShellTimeSeriesTrajectoryResult. Per-atom 4-scalar TR
// (COM-based water shell features + nearest ion).
//

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "GeometryResult.h"
#include "HydrationGeometryResult.h"
#include "HydrationShellResult.h"
#include "HydrationShellTimeSeriesTrajectoryResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
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
nmr::RunConfiguration BuildConfig(unsigned stride) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::HydrationShellResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::HydrationShellTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


TEST(HydrationShellTimeSeries, Frame0Semantics) {
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
    const auto& tr = tp.Result<nmr::HydrationShellTimeSeriesTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);
}

TEST(HydrationShellTimeSeries, hydration_coincident_water_guard) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();

    nmr::Protein protein;
    nmr::Residue residue;
    residue.type = nmr::AminoAcid::ALA;
    residue.sequence_number = 1;
    const size_t ri = protein.AddResidue(std::move(residue));
    for (nmr::Element element : {nmr::Element::C, nmr::Element::N}) {
        auto atom = nmr::Atom::Create(element);
        atom->residue_index = ri;
        const size_t ai = protein.AddAtom(std::move(atom));
        protein.MutableResidueAt(ri).atom_indices.push_back(ai);
    }

    auto& conf = protein.AddConformation(
        {nmr::Vec3(0.0, 0.0, 0.0), nmr::Vec3(10.0, 0.0, 0.0)},
        "hydration coincident water");
    conf.MutableAtomAt(0).sasa_normal = nmr::Vec3(1.0, 0.0, 0.0);

    nmr::SolventEnvironment solvent;
    nmr::WaterMolecule water;
    water.O_pos = nmr::Vec3(0.0, 0.0, 0.0);
    water.H1_pos = nmr::Vec3(0.1, 0.0, 0.0);
    water.H2_pos = nmr::Vec3(0.0, 0.1, 0.0);
    water.O_charge = -0.834;
    water.H_charge = 0.417;
    solvent.waters.push_back(water);
    solvent.water_O_positions.push_back(water.O_pos);

    auto shell = nmr::HydrationShellResult::Compute(conf, solvent);
    ASSERT_NE(shell, nullptr);
    auto geometry = nmr::HydrationGeometryResult::Compute(conf, solvent);
    ASSERT_NE(geometry, nullptr);

    const auto& atom0 = conf.AtomAt(0);
    EXPECT_TRUE(std::isfinite(atom0.half_shell_asymmetry));
    EXPECT_TRUE(std::isfinite(atom0.mean_water_dipole_cos));
    EXPECT_TRUE(std::isfinite(atom0.sasa_half_shell_asymmetry));
    EXPECT_TRUE(std::isfinite(atom0.sasa_dipole_alignment));
    EXPECT_TRUE(std::isfinite(atom0.sasa_dipole_coherence));
    EXPECT_EQ(atom0.sasa_first_shell_count, 0);

    bool saw_shell_guard = false;
    bool saw_geometry_guard = false;
    auto has_rejected_one = [](const nmr::GeometryChoice& choice) {
        for (const auto& number : choice.Numbers()) {
            if (number.name == "waters_rejected" && number.value == 1.0)
                return true;
        }
        return false;
    };
    for (const auto& choice : conf.geometry_choices) {
        if (choice.Label() == "hydration shell singularity guard")
            saw_shell_guard = has_rejected_one(choice);
        if (choice.Label() == "hydration geometry singularity guard")
            saw_geometry_guard = has_rejected_one(choice);
    }
    EXPECT_TRUE(saw_shell_guard);
    EXPECT_TRUE(saw_geometry_guard);
}


TEST(HydrationShellTimeSeries, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    auto& tr = tp.Result<nmr::HydrationShellTimeSeriesTrajectoryResult>();
    const std::size_t T = tr.NumFrames();
    tr.Finalize(tp, traj);
    EXPECT_EQ(tr.NumFrames(), T);
}


TEST(HydrationShellTimeSeries, H5RoundTrip) {
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

    const auto& tr = tp.Result<nmr::HydrationShellTimeSeriesTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("hyd_shell_ts_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/hydration_shell_time_series"));
    auto grp = reopen.getGroup("/trajectory/hydration_shell_time_series");

    std::string ref_frame;
    grp.getAttribute("reference_frame").read(ref_frame);
    EXPECT_EQ(ref_frame, "COM");

    for (const std::string& ch : {"half_shell_asymmetry", "mean_water_dipole_cos",
                                   "nearest_ion_distance", "nearest_ion_charge"}) {
        EXPECT_TRUE(grp.exist(ch)) << ch;
        const auto dims = grp.getDataSet(ch).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], tp.AtomCount());
        EXPECT_EQ(dims[1], tr.NumFrames());
    }
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));

    fs::remove(h5_path);
}


TEST(HydrationShellTimeSeries, Integration1P9J) {
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

    const auto& tr = tp.Result<nmr::HydrationShellTimeSeriesTrajectoryResult>();
    EXPECT_GE(tr.NumFrames(), 2u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("hyd_shell_ts_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/hydration_shell_time_series");

    std::vector<std::vector<double>> ion_dist;
    grp.getDataSet("nearest_ion_distance").read(ion_dist);
    ASSERT_EQ(ion_dist.size(), tp.AtomCount());

    std::size_t with_ion = 0, without_ion = 0;
    for (std::size_t i = 0; i < tp.AtomCount(); ++i) {
        if (std::isfinite(ion_dist[i][0])) ++with_ion; else ++without_ion;
    }
    std::cout << "HydrationShellTimeSeries: " << tr.NumFrames() << " frames; "
              << "atoms with ion in cutoff at frame 0: " << with_ion
              << "/" << tp.AtomCount() << " (without: " << without_ion << ")\n";

    fs::remove(h5_path);
}


TEST(HydrationShellTimeSeries, SyntheticAllAbsentSkipsGroup) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::HydrationShellTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < 3; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    const std::string h5_path = (fs::temp_directory_path() /
        ("hyd_shell_ts_absent_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/hydration_shell_time_series"))
        << "All-absent run should skip group emission.";

    fs::remove(h5_path);
}
