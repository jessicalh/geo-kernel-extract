//
// test_tripeptide_neighbor_residual_vec_next_time_series: discipline +
// integration tests for TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult.
//
// Mirrors the prev-direction test; identical NaN-tolerant integration
// pattern against the i+1 source field.
//

#include "AIMNet2Result.h"

#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "BiotSavartResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DenseBuffer.h"
#include "DsspResult.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "RuntimeEnvironment.h"
#include "Session.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "TripeptideNeighborShieldingResult.h"
#include "TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult.h"
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
constexpr double kResidualMagWarnAngstrom = 5.0;

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

void LoadCalculatorConfig() {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../data/calculator_params.toml");
}

nmr::Vec3 SyntheticVec3(size_t atom_i, size_t frame_t) {
    const double base = static_cast<double>(atom_i)
                      + static_cast<double>(frame_t) * 100.0;
    return nmr::Vec3(base, base + 0.1, base + 0.2);
}

bool Vec3Equal(const nmr::Vec3& a, const nmr::Vec3& b, double tol) {
    return std::abs(a.x() - b.x()) <= tol
        && std::abs(a.y() - b.y()) <= tol
        && std::abs(a.z() - b.z()) <= tol;
}

}  // namespace


TEST(TripeptideNeighborResidualVecNextTimeSeries, SyntheticFourFrames) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix))
        GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    const size_t Ntp = tp.AtomCount();
    ASSERT_GT(Ntp, 0u);

    auto tr =
        nmr::TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult::Create(tp);
    // Synthetic path bypasses OperationRunner so the real source calc
    // never attaches. Force-mark source present so the WriteH5Group
    // skip-on-absent gate doesn't fire.
    tr->ForceSourcePresentForTesting();

    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 4;
    const auto& protein_ref = tp.ProteinRef();
    std::vector<nmr::Vec3> positions(Ntp, nmr::Vec3::Zero());

    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &protein_ref, positions, "synthetic frame");
        for (size_t i = 0; i < Ntp; ++i) {
            conf->MutableAtomAt(i).tripeptide_neighbor_residual_vec_next =
                SyntheticVec3(i, t);
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    EXPECT_EQ(tr->NumFrames(), kFrames);

    tr->Finalize(tp, traj);

    auto* buf = tp.GetDenseBuffer<nmr::Vec3>(std::type_index(typeid(
        nmr::TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), Ntp);
    EXPECT_EQ(buf->StridePerAtom(), kFrames);

    for (size_t i : {size_t(0), Ntp / 2, Ntp - 1}) {
        for (size_t t = 0; t < kFrames; ++t) {
            const auto expected = SyntheticVec3(i, t);
            const auto& got = buf->At(i, t);
            EXPECT_TRUE(Vec3Equal(got, expected, 1e-12))
                << "buffer mismatch at atom " << i << " frame " << t;
        }
    }

    const std::string h5_path = (fs::temp_directory_path() /
        ("tripeptide_neighbor_residual_vec_next_ts_unit_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr->WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist(
        "/trajectory/tripeptide_neighbor_residual_vec_next_time_series"));
    auto grp = reopen.getGroup(
        "/trajectory/tripeptide_neighbor_residual_vec_next_time_series");
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[0], Ntp);
    EXPECT_EQ(dims[1], kFrames);
    EXPECT_EQ(dims[2], 3u);

    std::vector<double> flat(Ntp * kFrames * 3);
    ds.read(flat.data());
    const size_t i = Ntp / 2;
    const size_t t = 2;
    const size_t base = (i * kFrames + t) * 3;
    const auto expected = SyntheticVec3(i, t);
    EXPECT_DOUBLE_EQ(flat[base + 0], expected.x());
    EXPECT_DOUBLE_EQ(flat[base + 1], expected.y());
    EXPECT_DOUBLE_EQ(flat[base + 2], expected.z());

    fs::remove(h5_path);
}


TEST(TripeptideNeighborResidualVecNextTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("TripeptideNeighborResidualVecNextTimeSeriesFrame0Semantics");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult
                    ::Create(tp_in);
        });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    ASSERT_EQ(traj.FrameCount(), 1u);

    ASSERT_TRUE(tp.HasResult<
        nmr::TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult>());

    auto* buf = tp.GetDenseBuffer<nmr::Vec3>(std::type_index(
        typeid(nmr::TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), tp.AtomCount());
    EXPECT_EQ(buf->StridePerAtom(), 1u);
}


TEST(TripeptideNeighborResidualVecNextTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("TripeptideNeighborResidualVecNextTimeSeriesFinalizeIdempotency");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult
                    ::Create(tp_in);
        });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto* buf_first = tp.GetDenseBuffer<nmr::Vec3>(
        std::type_index(typeid(
            nmr::TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_first, nullptr);
    const std::size_t N_first = buf_first->AtomCount();
    const std::size_t T_first = buf_first->StridePerAtom();

    auto& tr = tp.Result<
        nmr::TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);

    auto* buf_second = tp.GetDenseBuffer<nmr::Vec3>(
        std::type_index(typeid(
            nmr::TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_second, nullptr);
    EXPECT_EQ(buf_second->AtomCount(), N_first);
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
    EXPECT_EQ(tr.NumFrames(), T_first);
}


TEST(TripeptideNeighborResidualVecNextTimeSeries, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    if (nmr::RuntimeEnvironment::TensorCs15Dsn().empty()) {
        GTEST_SKIP() << "tensorcs15 DSN not configured — source calc "
                        "won't attach so the TR correctly skips H5 emission, "
                        "and round-trip cannot be verified";
    }
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::Session session;
    ASSERT_EQ(session.LoadTripeptideDftTable(), nmr::kOk)
        << session.LastError();
    ASSERT_TRUE(session.HasTripeptideDftTable());

    nmr::RunConfiguration config;
    config.SetName("TripeptideNeighborResidualVecNextTimeSeriesH5RoundTrip");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = false;  // backbone phi/psi needed for the tripeptide calc to actually run
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::DsspResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult
                    ::Create(tp_in);
        });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<
        nmr::TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult>();

    const std::string h5_path = (fs::temp_directory_path() /
        ("tripeptide_neighbor_residual_vec_next_ts_h5_roundtrip_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr.WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist(
        "/trajectory/tripeptide_neighbor_residual_vec_next_time_series"));

    auto grp = reopen.getGroup(
        "/trajectory/tripeptide_neighbor_residual_vec_next_time_series");
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[0], tp.AtomCount());
    EXPECT_EQ(dims[1], 1u);
    EXPECT_EQ(dims[2], 3u);

    std::string parity, normalization, units, layout;
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("normalization").read(normalization);
    grp.getAttribute("units").read(units);
    grp.getAttribute("irrep_layout").read(layout);
    EXPECT_EQ(parity, "1o");
    EXPECT_EQ(normalization, "cartesian");
    EXPECT_EQ(units, "angstrom");
    EXPECT_EQ(layout, "x,y,z");

    fs::remove(h5_path);
}


TEST(TripeptideNeighborResidualVecNextTimeSeries, IntegrationLogOverages1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    if (nmr::RuntimeEnvironment::TensorCs15Dsn().empty()) {
        GTEST_SKIP() << "tensorcs15 DSN not configured";
    }
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::Session session;
    ASSERT_EQ(session.LoadTripeptideDftTable(), nmr::kOk)
        << session.LastError();
    ASSERT_TRUE(session.HasTripeptideDftTable());

    nmr::RunConfiguration config;
    config.SetName("TripeptideNeighborResidualVecNextTimeSeriesIntegration");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = false;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::DsspResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult
                    ::Create(tp_in);
        });
    config.SetStride(300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_GE(traj.FrameCount(), 1u);

    auto* buf = tp.GetDenseBuffer<nmr::Vec3>(std::type_index(
        typeid(nmr::TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);

    std::size_t finite_cells = 0;
    std::size_t nan_cells = 0;
    std::size_t finite_nonzero_at_frame0 = 0;
    std::size_t overage_count = 0;
    double max_seen_magnitude = 0.0;
    for (std::size_t i = 0; i < buf->AtomCount(); ++i) {
        for (std::size_t t = 0; t < buf->StridePerAtom(); ++t) {
            const nmr::Vec3& v = buf->At(i, t);
            const bool all_finite = std::isfinite(v.x())
                                 && std::isfinite(v.y())
                                 && std::isfinite(v.z());
            const bool all_nan = std::isnan(v.x())
                              && std::isnan(v.y())
                              && std::isnan(v.z());
            EXPECT_TRUE(all_finite || all_nan)
                << "atom " << i << " frame " << t << " has mixed "
                << "finite/NaN/Inf components";
            if (all_finite) {
                ++finite_cells;
                const double mag = v.norm();
                if (mag > max_seen_magnitude) max_seen_magnitude = mag;
                if (mag > kResidualMagWarnAngstrom) ++overage_count;
                if (t == 0 && mag > 1e-12) ++finite_nonzero_at_frame0;
            } else {
                ++nan_cells;
            }
        }
    }

    if (overage_count > 0) {
        std::cout << "[warn] TripeptideNeighborResidualVecNext: "
                  << overage_count << " finite cells exceed "
                  << kResidualMagWarnAngstrom << " Å (max norm "
                  << max_seen_magnitude << " Å) — chi-grid / perception "
                  << "edges, not asserted\n";
    }

    // Population floor on finite, nonzero cells. C-terminal residue has
    // no i+1 neighbour; still expect a majority finite.
    EXPECT_GT(finite_nonzero_at_frame0, buf->AtomCount() / 3)
        << "next residual_vec finite-nonzero at frame 0: "
        << finite_nonzero_at_frame0 << " / " << buf->AtomCount()
        << " — perception or i+1 coverage drop?";

    std::cout << "TripeptideNeighborResidualVecNext diagnostic: "
              << "frames=" << traj.FrameCount()
              << " finite_cells=" << finite_cells
              << " nan_cells=" << nan_cells
              << " finite_nonzero_frame0=" << finite_nonzero_at_frame0
              << "/" << buf->AtomCount()
              << " max=" << max_seen_magnitude << " Å"
              << " overages=" << overage_count << "\n";
}
