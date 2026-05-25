//
// test_dssp8_transition: discipline + integration for
// Dssp8TransitionTrajectoryResult. AV companion to Dssp8TimeSeries.
//

#include "CalculatorConfig.h"
#include "DsspResult.h"
#include "Dssp8TransitionTrajectoryResult.h"
#include "GeometryResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "RunConfiguration.h"
#include "Session.h"
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
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true;
    opts.skip_dssp = false;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::DsspResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::Dssp8TransitionTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


TEST(Dssp8Transition, H5RoundTripSingleFrame) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::Dssp8TransitionTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    EXPECT_EQ(tr.NumFrames(), 1u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("dssp8_trans_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/dssp8_transition"));
    auto grp = reopen.getGroup("/trajectory/dssp8_transition");

    for (const std::string& name : {"ss8_transition_count",
                                     "ss8_dominant",
                                     "n_frames_observed"}) {
        ASSERT_TRUE(grp.exist(name)) << name;
        const auto dims = grp.getDataSet(name).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 1u);
        EXPECT_EQ(dims[0], R);
    }
    {
        const auto dims =
            grp.getDataSet("ss8_occupancy").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], 8u);
    }
    {
        const auto dims =
            grp.getDataSet("ss8_transition_matrix").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 3u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], 8u);
        EXPECT_EQ(dims[2], 8u);
    }

    // Single-frame run: zero transitions, occupancy sums == n_obs.
    std::vector<std::uint32_t> trans_count;
    grp.getDataSet("ss8_transition_count").read(trans_count);
    for (auto c : trans_count) EXPECT_EQ(c, 0u);

    fs::remove(h5_path);
}


TEST(Dssp8Transition, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(300);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::Dssp8TransitionTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    const std::size_t T = tr.NumFrames();
    EXPECT_GE(T, 2u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("dssp8_trans_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/dssp8_transition");

    // Invariants:
    // (1) Sum of ss8_occupancy across all 8 states == n_frames_observed.
    // (2) Transition matrix diagonal is zero (self-transitions excluded).
    // (3) Sum of off-diagonal transition matrix == ss8_transition_count.
    std::vector<std::uint32_t> trans_count, n_obs;
    grp.getDataSet("ss8_transition_count").read(trans_count);
    grp.getDataSet("n_frames_observed").read(n_obs);
    std::vector<std::vector<std::uint32_t>> occ;
    grp.getDataSet("ss8_occupancy").read(occ);
    std::vector<std::vector<std::vector<std::uint32_t>>> mat;
    grp.getDataSet("ss8_transition_matrix").read(mat);

    std::size_t total_trans = 0, total_obs = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        std::uint32_t occ_sum = 0;
        for (std::size_t s = 0; s < 8; ++s) occ_sum += occ[ri][s];
        EXPECT_EQ(occ_sum, n_obs[ri])
            << "occupancy sum != n_frames_observed at ri=" << ri;

        std::uint32_t mat_off_diag_sum = 0;
        for (std::size_t p = 0; p < 8; ++p) {
            for (std::size_t c = 0; c < 8; ++c) {
                if (p == c) {
                    EXPECT_EQ(mat[ri][p][c], 0u)
                        << "diagonal nonzero at ri=" << ri
                        << " state=" << p;
                } else {
                    mat_off_diag_sum += mat[ri][p][c];
                }
            }
        }
        EXPECT_EQ(mat_off_diag_sum, trans_count[ri])
            << "matrix off-diagonal sum != transition_count at ri=" << ri;

        total_trans += trans_count[ri];
        total_obs   += n_obs[ri];
    }
    std::cout << "Dssp8Transition: " << R << " residues, " << T
              << " frames; total_obs=" << total_obs
              << " total_transitions=" << total_trans << "\n";

    fs::remove(h5_path);
}
