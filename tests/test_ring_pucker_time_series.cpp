//
// test_ring_pucker_time_series: discipline + integration for
// RingPuckerTimeSeriesTrajectoryResult.
//

#include "CalculatorConfig.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "LegacyAmberTopology.h"
#include "OperationLog.h"
#include "PlanarGeometryResult.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Ring.h"
#include "RingPuckerTimeSeriesTrajectoryResult.h"
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
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::EnrichmentResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::RingPuckerTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


TEST(RingPuckerTimeSeries, H5RoundTripSingleFrame) {
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

    const auto& tr = tp.Result<nmr::RingPuckerTimeSeriesTrajectoryResult>();
    const auto& topo = tp.ProteinRef().LegacyAmber();
    const std::size_t S = topo.SaturatedRingCount();
    const std::size_t A = topo.AromaticRingCount();
    EXPECT_EQ(tr.NumFrames(), 1u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("ring_pucker_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/ring_pucker_time_series"));
    auto grp = reopen.getGroup("/trajectory/ring_pucker_time_series");

    // Per-ring lookups always present.
    EXPECT_TRUE(grp.exist("saturated_parent_residue_index"));
    EXPECT_TRUE(grp.exist("aromatic_parent_residue_index"));

    if (S > 0) {
        EXPECT_TRUE(grp.exist("pucker_Q"));
        EXPECT_TRUE(grp.exist("pucker_theta"));
        const auto dims = grp.getDataSet("pucker_Q").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], S);
        EXPECT_EQ(dims[1], tr.NumFrames());
    }
    if (A > 0) {
        EXPECT_TRUE(grp.exist("aromatic_chi2"));
        const auto dims =
            grp.getDataSet("aromatic_chi2").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], A);
        EXPECT_EQ(dims[1], tr.NumFrames());
    }

    EXPECT_TRUE(grp.exist("source_attached_per_frame"));

    // Convention attrs
    std::string conv, Q_units, theta_units;
    grp.getAttribute("pucker_convention").read(conv);
    grp.getAttribute("pucker_Q_units").read(Q_units);
    grp.getAttribute("pucker_theta_units").read(theta_units);
    EXPECT_NE(conv.find("Cremer-Pople"), std::string::npos);
    EXPECT_EQ(Q_units, "Angstrom");
    EXPECT_EQ(theta_units, "degrees");

    std::cout << "RingPuckerTimeSeries: S=" << S << " A=" << A << " T="
              << tr.NumFrames() << "\n";

    fs::remove(h5_path);
}


TEST(RingPuckerTimeSeries, Integration1P9J) {
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

    const auto& tr = tp.Result<nmr::RingPuckerTimeSeriesTrajectoryResult>();
    const auto& topo = tp.ProteinRef().LegacyAmber();
    const std::size_t S = topo.SaturatedRingCount();
    const std::size_t A = topo.AromaticRingCount();
    const std::size_t T = tr.NumFrames();
    EXPECT_GE(T, 2u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("ring_pucker_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/ring_pucker_time_series");

    if (A > 0) {
        std::vector<std::vector<double>> chi2;
        grp.getDataSet("aromatic_chi2").read(chi2);
        std::size_t finite_chi2 = 0;
        for (std::size_t ai = 0; ai < A; ++ai) {
            for (std::size_t f = 0; f < T; ++f) {
                if (std::isfinite(chi2[ai][f])) {
                    ++finite_chi2;
                    EXPECT_GE(chi2[ai][f], -M_PI);
                    EXPECT_LE(chi2[ai][f],  M_PI);
                }
            }
        }
        EXPECT_GT(finite_chi2, 0u);
        std::cout << "  aromatic_chi2 finite: " << finite_chi2
                  << " / " << (A*T) << "\n";
    }

    if (S > 0) {
        std::vector<std::vector<double>> Q;
        grp.getDataSet("pucker_Q").read(Q);
        std::size_t finite_Q = 0;
        for (std::size_t si = 0; si < S; ++si) {
            for (std::size_t f = 0; f < T; ++f) {
                if (std::isfinite(Q[si][f])) {
                    ++finite_Q;
                    EXPECT_GE(Q[si][f], 0.0)
                        << "pucker_Q (amplitude) must be non-negative";
                    EXPECT_LE(Q[si][f], 1.0)
                        << "pucker_Q implausibly large for a 5-ring";
                }
            }
        }
        std::cout << "  pucker_Q finite: " << finite_Q
                  << " / " << (S*T) << "\n";
    }

    fs::remove(h5_path);
}
