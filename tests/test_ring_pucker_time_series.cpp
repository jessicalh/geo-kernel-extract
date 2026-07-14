//
// test_ring_pucker_time_series: discipline + integration for
// RingPuckerTimeSeriesTrajectoryResult.
//

#include "CalculatorConfig.h"
#include "DirectionalTestHelpers.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "LegacyAmberTopology.h"
#include "OperationLog.h"
#include "PlanarGeometryResult.h"
#include "PdbFileReader.h"
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
#include <cstdint>
#include <filesystem>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
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

template <typename T>
std::vector<T> ReadH5Flat(const std::string& path,
                          const std::string& dataset,
                          std::vector<std::size_t>* dimensions = nullptr) {
    HighFive::File file(path, HighFive::File::ReadOnly);
    auto data_set = file.getDataSet(dataset);
    const auto dims = data_set.getSpace().getDimensions();
    if (dimensions) *dimensions = dims;
    const std::size_t count = std::accumulate(
        dims.begin(), dims.end(), std::size_t{1},
        std::multiplies<std::size_t>());
    std::vector<T> values(count);
    if (!values.empty()) data_set.read(values.data());
    return values;
}

void ExpectDirectionalPuckerMetadata(const std::string& path) {
    HighFive::File file(path, HighFive::File::ReadOnly);
    const std::string group = "/trajectory/ring_pucker_time_series/";
    auto expect_strings = [&](const char* dataset_name,
                              const std::vector<std::pair<
                                  std::string, std::string>>& expected) {
        const std::string dataset_path = group + dataset_name;
        auto dataset = file.getDataSet(dataset_path);
        for (const auto& [name, value] : expected) {
            std::string actual;
            dataset.getAttribute(name).read(actual);
            EXPECT_EQ(actual, value) << dataset_path << " @" << name;
        }
    };
    expect_strings("pucker_Q", {
        {"coordinate_frame", "intrinsic_cremer_pople_5ring"},
        {"irrep_layout", "0e"},
        {"parity", "even"},
        {"transformation",
         "exact O(3)- and translation-invariant pucker amplitude: Q'=Q"},
        {"validity",
         "source_attached_per_frame gates source absence; NaN marks an "
         "incomplete ring or degenerate mean-plane construction"},
    });
    expect_strings("pucker_theta", {
        {"coordinate_frame", "intrinsic_oriented_cremer_pople_5ring_phase"},
        {"irrep_layout", "none_periodic_phase_with_improper_offset"},
        {"parity", "mixed"},
        {"transformation",
         "translation/proper-rotation invariant phase; under an improper "
         "orthogonal transform theta'=(theta+180 degrees) mod 360"},
        {"validity",
         "source_attached_per_frame gates source absence; NaN marks an "
         "incomplete/degenerate ring or pucker_Q < 1e-6 Angstrom; compare "
         "finite values with circular differences"},
    });
    expect_strings("aromatic_chi2", {
        {"coordinate_frame", "intrinsic_signed_dihedral"},
        {"irrep_layout", "0o"},
        {"parity", "odd"},
        {"transformation",
         "wrapped signed-dihedral pseudoscalar: chi2'=det(R) chi2 modulo "
         "2pi; translation invariant"},
        {"validity",
         "source_attached_per_frame gates source absence; NaN marks an "
         "unavailable or geometrically degenerate parent-residue chi2; "
         "compare finite values with circular differences"},
    });
}

nmr::RunConfiguration BuildConfig(unsigned stride) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::EnrichmentResult));
    config.RequireConformationResult(typeid(nmr::PlanarGeometryResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::RingPuckerTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


TEST(RingPuckerTimeSeries,
     DirectionalPhasesRerunPlanarOwnerO3SerializedH5) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    const std::string pdb = nmr::test::TestEnvironment::UbqProtonated();
    if (pdb.empty() || !fs::exists(pdb))
        GTEST_SKIP() << "1UBQ protonated fixture unavailable";
    auto build = nmr::BuildFromProtonatedPdb(pdb);
    ASSERT_TRUE(build.Ok()) << build.error;
    auto protein = std::move(build.protein);
    auto& original = protein->Conformation();
    auto attach_planar = [](nmr::ProteinConformation& conf) {
        auto geometry = nmr::GeometryResult::Compute(conf);
        ASSERT_NE(geometry, nullptr);
        ASSERT_TRUE(conf.AttachResult(std::move(geometry)));
        auto enrichment = nmr::EnrichmentResult::Compute(conf);
        ASSERT_NE(enrichment, nullptr);
        ASSERT_TRUE(conf.AttachResult(std::move(enrichment)));
        auto result = nmr::PlanarGeometryResult::Compute(conf);
        ASSERT_NE(result, nullptr);
        ASSERT_TRUE(conf.AttachResult(std::move(result)));
    };
    attach_planar(original);

    const auto proper = nmr::test::directional::SeededTransform(
        0x5055434B45525453ULL, false);
    const auto improper = nmr::test::directional::SeededTransform(
        0x5055434B45525453ULL, true);
    auto& moved_proper = protein->AddConformation(
        nmr::test::directional::Positions(proper, original.Positions()),
        "ring pucker proper owner rerun");
    attach_planar(moved_proper);
    auto& moved_improper = protein->AddConformation(
        nmr::test::directional::Positions(improper, original.Positions()),
        "ring pucker improper owner rerun");
    attach_planar(moved_improper);

    auto tp = nmr::TrajectoryProtein::CreateForTesting(std::move(protein));
    ASSERT_NE(tp, nullptr);
    auto result = nmr::RingPuckerTimeSeriesTrajectoryResult::Create(*tp);
    ASSERT_NE(result, nullptr);
    nmr::Trajectory dummy("", "", "");
    for (std::size_t frame = 0; frame < 3; ++frame) {
        result->Compute(tp->ProteinRef().ConformationAt(frame), *tp, dummy,
                        51u + frame, 0.5 * frame);
    }
    result->Finalize(*tp, dummy);

    const std::string h5_path =
        nmr::test::TestEnvironment::TempPath(
            "ring_pucker_directional_covariance.h5");
    fs::remove(h5_path);
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        result->WriteH5Group(*tp, file);
    }
    // The one H5 payload carries source, proper, and improper production
    // reruns as separate frame columns; assert every directional dataset's
    // exact contract on that serialized output.
    ExpectDirectionalPuckerMetadata(h5_path);
    const std::string group = "/trajectory/ring_pucker_time_series/";
    const auto& topology = tp->ProteinRef().LegacyAmber();
    const std::size_t saturated =
        topology.SaturatedRingCount();
    const std::size_t aromatic =
        topology.AromaticRingCount();
    ASSERT_GT(saturated, 0u);
    ASSERT_GT(aromatic, 0u);
    constexpr double kScalarAbsTolerance = 2.0e-10;

    std::vector<std::size_t> dims;
    const auto q = ReadH5Flat<double>(
        h5_path, group + "pucker_Q", &dims);
    ASSERT_EQ(dims,
              (std::vector<std::size_t>{saturated, 3u}));
    const auto theta = ReadH5Flat<double>(
        h5_path, group + "pucker_theta", &dims);
    ASSERT_EQ(dims,
              (std::vector<std::size_t>{saturated, 3u}));
    std::size_t finite_puckers = 0;
    for (std::size_t ring = 0; ring < saturated; ++ring) {
        const double source_q = q[ring * 3u];
        const double source_theta = theta[ring * 3u];
        if (std::isnan(source_q) || std::isnan(source_theta)) {
            EXPECT_TRUE(std::isnan(q[ring * 3u + 1u]));
            EXPECT_TRUE(std::isnan(q[ring * 3u + 2u]));
            EXPECT_TRUE(std::isnan(theta[ring * 3u + 1u]));
            EXPECT_TRUE(std::isnan(theta[ring * 3u + 2u]));
            continue;
        }
        EXPECT_NEAR(q[ring * 3u + 1u], source_q,
                    kScalarAbsTolerance);
        EXPECT_NEAR(q[ring * 3u + 2u], source_q,
                    kScalarAbsTolerance);
        EXPECT_NEAR(theta[ring * 3u + 1u], source_theta,
                    kScalarAbsTolerance);
        const double expected_improper =
            std::fmod(source_theta + 180.0, 360.0);
        const double theta_error = std::remainder(
            theta[ring * 3u + 2u] - expected_improper, 360.0);
        EXPECT_NEAR(theta_error, 0.0, kScalarAbsTolerance)
            << group << "pucker_theta ring=" << ring;
        ++finite_puckers;
    }
    EXPECT_GT(finite_puckers, 0u);

    const auto chi2 = ReadH5Flat<double>(
        h5_path, group + "aromatic_chi2", &dims);
    ASSERT_EQ(dims,
              (std::vector<std::size_t>{aromatic, 3u}));
    std::size_t finite_chi2 = 0;
    for (std::size_t ring = 0; ring < aromatic; ++ring) {
        const double source = chi2[ring * 3u];
        if (std::isnan(source)) {
            EXPECT_TRUE(std::isnan(chi2[ring * 3u + 1u]));
            EXPECT_TRUE(std::isnan(chi2[ring * 3u + 2u]));
            continue;
        }
        EXPECT_NEAR(nmr::test::directional::CircularDifference(
                        chi2[ring * 3u + 1u], source),
                    0.0, kScalarAbsTolerance);
        EXPECT_NEAR(nmr::test::directional::CircularDifference(
                        chi2[ring * 3u + 2u], -source),
                    0.0, kScalarAbsTolerance);
        ++finite_chi2;
    }
    EXPECT_GT(finite_chi2, 0u);
    EXPECT_EQ(ReadH5Flat<std::uint8_t>(
                  h5_path, group + "source_attached_per_frame"),
              (std::vector<std::uint8_t>{1u, 1u, 1u}));
    EXPECT_EQ(ReadH5Flat<std::size_t>(h5_path, group + "frame_indices"),
              (std::vector<std::size_t>{51u, 52u, 53u}));
    EXPECT_EQ(ReadH5Flat<double>(h5_path, group + "frame_times"),
              (std::vector<double>{0.0, 0.5, 1.0}));
    std::vector<std::int32_t> expected_saturated_parent(saturated);
    for (std::size_t ring = 0; ring < saturated; ++ring) {
        expected_saturated_parent[ring] = static_cast<std::int32_t>(
            topology.SaturatedRingAt(ring).parent_residue_index);
    }
    std::vector<std::int32_t> expected_aromatic_parent(aromatic);
    for (std::size_t ring = 0; ring < aromatic; ++ring) {
        expected_aromatic_parent[ring] = static_cast<std::int32_t>(
            topology.AromaticRingAt(ring).parent_residue_index);
    }
    EXPECT_EQ(ReadH5Flat<std::int32_t>(
                  h5_path, group + "saturated_parent_residue_index"),
              expected_saturated_parent);
    EXPECT_EQ(ReadH5Flat<std::int32_t>(
                  h5_path, group + "aromatic_parent_residue_index"),
              expected_aromatic_parent);
    EXPECT_TRUE(fs::remove(h5_path));
}


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
