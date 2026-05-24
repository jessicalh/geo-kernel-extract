//
// test_larsen_hbond_2pHB_shielding_time_series: clone of 1pHB tests
// against larsen_hbond_2pHB_spherical.
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
#include "LarsenHBondShieldingResult.h"
#include "LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult.h"
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
constexpr double kPerClassT0SanityPpm = 50.0;

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

nmr::SphericalTensor SyntheticTensor(size_t atom_i, size_t frame_t) {
    nmr::SphericalTensor s;
    s.T0 = static_cast<double>(atom_i) + static_cast<double>(frame_t) * 100.0;
    for (size_t k = 0; k < 3; ++k)
        s.T1[k] = static_cast<double>(atom_i)
                  + static_cast<double>(frame_t) * 100.0
                  + static_cast<double>(k + 1) * 1e-2;
    for (size_t k = 0; k < 5; ++k)
        s.T2[k] = static_cast<double>(atom_i)
                  + static_cast<double>(frame_t) * 100.0
                  + static_cast<double>(k + 1) * 1e-3;
    return s;
}
bool SphericalEqual(const nmr::SphericalTensor& a,
                    const nmr::SphericalTensor& b, double tol) {
    if (std::abs(a.T0 - b.T0) > tol) return false;
    for (size_t k = 0; k < 3; ++k)
        if (std::abs(a.T1[k] - b.T1[k]) > tol) return false;
    for (size_t k = 0; k < 5; ++k)
        if (std::abs(a.T2[k] - b.T2[k]) > tol) return false;
    return true;
}

}  // namespace


TEST(LarsenHBond2pHBShieldingTimeSeries, SyntheticFourFrames) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    const size_t Ntp = tp.AtomCount();
    ASSERT_GT(Ntp, 0u);

    auto tr = nmr::LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult::Create(tp);
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
        for (size_t i = 0; i < Ntp; ++i)
            conf->MutableAtomAt(i).larsen_hbond_2pHB_spherical =
                SyntheticTensor(i, t);
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }

    tr->Finalize(tp, traj);

    auto* buf = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(typeid(
        nmr::LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), Ntp);
    EXPECT_EQ(buf->StridePerAtom(), kFrames);

    for (size_t i : {size_t(0), Ntp / 2, Ntp - 1}) {
        for (size_t t = 0; t < kFrames; ++t) {
            EXPECT_TRUE(SphericalEqual(buf->At(i, t), SyntheticTensor(i, t),
                                       1e-12));
        }
    }

    const std::string h5_path = (fs::temp_directory_path() /
        ("larsen_hbond_2pHB_ts_unit_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr->WriteH5Group(tp, file);
    }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup(
        "/trajectory/larsen_hbond_2pHB_shielding_time_series");
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[0], Ntp);
    EXPECT_EQ(dims[1], kFrames);
    EXPECT_EQ(dims[2], 9u);

    fs::remove(h5_path);
}


TEST(LarsenHBond2pHBShieldingTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("LarsenHBond2pHBShieldingTimeSeriesFrame0Semantics");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true;
    opts.skip_apbs = true;  opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult
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

    auto* buf = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(
        typeid(nmr::LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), tp.AtomCount());
    EXPECT_EQ(buf->StridePerAtom(), 1u);
}


TEST(LarsenHBond2pHBShieldingTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("LarsenHBond2pHBShieldingTimeSeriesFinalizeIdempotency");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true;
    opts.skip_apbs = true;  opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult
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

    auto* buf_first = tp.GetDenseBuffer<nmr::SphericalTensor>(
        std::type_index(typeid(
            nmr::LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_first, nullptr);
    const std::size_t N_first = buf_first->AtomCount();
    const std::size_t T_first = buf_first->StridePerAtom();

    auto& tr = tp.Result<
        nmr::LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);

    auto* buf_second = tp.GetDenseBuffer<nmr::SphericalTensor>(
        std::type_index(typeid(
            nmr::LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_second, nullptr);
    EXPECT_EQ(buf_second->AtomCount(), N_first);
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
    EXPECT_EQ(tr.NumFrames(), T_first);
}


TEST(LarsenHBond2pHBShieldingTimeSeries, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    nmr::Session session;
    if (session.LoadLarsenHBondGrid() != nmr::kOk
        || !session.HasLarsenHBondGrid()) {
        GTEST_SKIP() << "larsen_hbond_grids not configured — source calc "
                        "won't attach so the TR correctly skips H5 emission, "
                        "and round-trip cannot be verified";
    }
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("LarsenHBond2pHBShieldingTimeSeriesH5RoundTrip");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true;
    opts.skip_apbs = true;  opts.skip_dssp = false;  // backbone phi/psi needed for the larsen H-bond calc to actually run
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::DsspResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult
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
        nmr::LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult>();

    const std::string h5_path = (fs::temp_directory_path() /
        ("larsen_hbond_2pHB_ts_h5_roundtrip_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr.WriteH5Group(tp, file);
    }

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup(
        "/trajectory/larsen_hbond_2pHB_shielding_time_series");
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[0], tp.AtomCount());
    EXPECT_EQ(dims[1], 1u);
    EXPECT_EQ(dims[2], 9u);

    std::string parity, normalization, units;
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("normalization").read(normalization);
    grp.getAttribute("units").read(units);
    EXPECT_EQ(parity, "0e+1o+2e");
    EXPECT_EQ(normalization, "isometric_real_sph");
    EXPECT_EQ(units, "ppm");

    fs::remove(h5_path);
}


TEST(LarsenHBond2pHBShieldingTimeSeries, IntegrationLogOverages1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();

    nmr::Session session;
    if (session.LoadLarsenHBondGrid() != nmr::kOk
        || !session.HasLarsenHBondGrid()) {
        GTEST_SKIP() << "larsen_hbond_grids not configured in "
                        "~/.nmr_tools.toml — Larsen calc cannot run, "
                        "integration assertion would be vacuous";
    }
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("LarsenHBond2pHBShieldingTimeSeriesIntegration");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true;
    opts.skip_apbs = true;  opts.skip_dssp = false;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::DsspResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult
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

    auto* buf = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(
        typeid(nmr::LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);

    std::size_t overage_count = 0;
    double max_t0 = 0.0;
    std::size_t nonzero_at_frame0 = 0;
    for (std::size_t i = 0; i < buf->AtomCount(); ++i) {
        for (std::size_t t = 0; t < buf->StridePerAtom(); ++t) {
            const nmr::SphericalTensor& st = buf->At(i, t);
            EXPECT_TRUE(std::isfinite(st.T0)
                     && std::isfinite(st.T1[0]) && std::isfinite(st.T1[1])
                     && std::isfinite(st.T1[2])
                     && std::isfinite(st.T2[0]) && std::isfinite(st.T2[1])
                     && std::isfinite(st.T2[2]) && std::isfinite(st.T2[3])
                     && std::isfinite(st.T2[4]))
                << "non-finite component at atom " << i << " frame " << t;
            const double abs_t0 = std::abs(st.T0);
            if (abs_t0 > max_t0) max_t0 = abs_t0;
            if (abs_t0 > kPerClassT0SanityPpm) ++overage_count;
            if (t == 0 && abs_t0 > 1e-12) ++nonzero_at_frame0;
        }
    }

    // Real-data coverage floor at frame 0: with the grid loaded, the
    // per-class calc should produce nonzero T0 on at least a small
    // fraction of atoms. Catches "every cell stayed at default 0.0"
    // regression that the source-attached gate guards against.
    const std::size_t coverage_floor = buf->AtomCount() / 100;
    EXPECT_GT(nonzero_at_frame0, coverage_floor)
        << "Larsen 2pHB produced nonzero T0 on only " << nonzero_at_frame0
        << " of " << buf->AtomCount() << " atoms at frame 0 (floor="
        << coverage_floor << ") — calc likely never ran despite grid loaded";

    if (overage_count > 0) {
        std::cout << "[warn] LarsenHBond2pHB: " << overage_count
                  << " cells exceed " << kPerClassT0SanityPpm << " ppm "
                  << "(max=" << max_t0 << ")\n";
    }
    std::cout << "LarsenHBond2pHB diagnostic: max_T0=" << max_t0
              << " nonzero_frame0=" << nonzero_at_frame0
              << " overages=" << overage_count << "\n";
}
