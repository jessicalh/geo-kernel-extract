//
// test_bs_shielding_time_series: focused H5 metadata test for
// BsShieldingTimeSeriesTrajectoryResult.
//

#include "BsShieldingTimeSeriesTrajectoryResult.h"
#include "PdbFileReader.h"
#include "ProteinConformation.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "Types.h"

#include <gtest/gtest.h>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <filesystem>
#include <memory>
#include <string>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

namespace {

constexpr const char* kFixtureProtein = "1P9J_5801";

std::string TrrPathFor(const std::string& tpr_path) {
    return fs::path(tpr_path).replace_extension(".trr").string();
}

std::string ProductionDirFor(const std::string& tpr_path) {
    return fs::path(tpr_path).parent_path().string();
}

bool FixtureAvailable(const nmr::test::AmberTrajectoryFixture& fix) {
    return !fix.tpr_path.empty() &&
           fs::exists(fix.tpr_path) &&
           fs::exists(TrrPathFor(fix.tpr_path)) &&
           fs::exists(fix.edr_path);
}

}  // namespace

TEST(BsShieldingTimeSeries, H5LayoutUsesCartesianT1Labels) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto build = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto tp_owner = nmr::TrajectoryProtein::CreateForTesting(
        std::move(build.protein));
    ASSERT_NE(tp_owner, nullptr);
    auto& tp = *tp_owner;

    auto tr = nmr::BsShieldingTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj({}, {}, {});
    std::vector<nmr::Vec3> positions(tp.AtomCount(), nmr::Vec3::Zero());
    auto conf = std::make_unique<nmr::ProteinConformation>(
        &tp.ProteinRef(), positions, "synthetic bs metadata frame");

    tr->Compute(*conf, tp, traj, 0, 0.0);
    tr->Finalize(tp, traj);

    const std::string h5_path = (fs::temp_directory_path() /
        ("bs_shielding_ts_h5_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/bs_shielding_time_series");
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[0], tp.AtomCount());
    EXPECT_EQ(dims[1], 1u);
    EXPECT_EQ(dims[2], 9u);

    std::string layout, parity, basis, order, frame, tensor_parity;
    std::string transformation, t1_semantics, structural_zeros;
    std::string e3nn_export, normalization_scope;
    bool t1_structural_zero = true;
    grp.getAttribute("irrep_layout").read(layout);
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("tensor_basis").read(basis);
    grp.getAttribute("tensor_component_order").read(order);
    grp.getAttribute("tensor_frame").read(frame);
    grp.getAttribute("tensor_parity").read(tensor_parity);
    grp.getAttribute("tensor_transformation").read(transformation);
    grp.getAttribute("tensor_t1_semantics").read(t1_semantics);
    grp.getAttribute("tensor_t1_structural_zero").read(t1_structural_zero);
    grp.getAttribute("tensor_structural_zero_components").read(structural_zeros);
    grp.getAttribute("e3nn_export").read(e3nn_export);
    grp.getAttribute("normalization_scope").read(normalization_scope);
    EXPECT_EQ(layout,
        "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");
    EXPECT_EQ(parity, "0e+1e+2e");
    EXPECT_EQ(basis, "project_native_full9_spherical_tensor_v1");
    EXPECT_EQ(order,
        "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");
    EXPECT_EQ(frame, "conformation_cartesian_xyz");
    EXPECT_EQ(tensor_parity, "even");
    EXPECT_EQ(transformation, "even_rank2: T'=R T R^T");
    EXPECT_EQ(t1_semantics,
        "Cartesian Levi-Civita dual x,y,z (not real-Y1m); axial "
        "a'=det(R) R a; generically nonzero");
    EXPECT_FALSE(t1_structural_zero);
    EXPECT_EQ(structural_zeros, "none");
    EXPECT_EQ(e3nn_export,
        "explicit project-basis to e3nn conversion required before use");
    EXPECT_EQ(normalization_scope,
        "T2 uses isometric real-tesseral normalization; T1 is Cartesian "
        "Levi-Civita dual");

    fs::remove(h5_path);
}
