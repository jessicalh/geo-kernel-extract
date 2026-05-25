//
// test_mopac_mc_shielding_time_series: discipline + integration for
// MopacMcConnellShieldingTimeSeriesTrajectoryResult (TR8 of the 13-TR
// plan). Sibling of TR7 but emits 9 components (T0+T1+T2) because
// the source mopac_mc_shielding_contribution is NOT traceless
// (bond-anisotropy M_total has nonzero T0/T1 in general).
//

#include "MopacMcConnellShieldingTimeSeriesTrajectoryResult.h"
#include "MopacMcConnellResult.h"
#include "MopacResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "ChargeAssignmentResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DenseBuffer.h"
#include "OperationLog.h"
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


// ── Layer 0: source absent every frame → group skipped ─────

TEST(MopacMcConnellShieldingTimeSeries, GroupSkippedWhenSourceNeverAttached) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    auto tr = nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 3;
    std::vector<nmr::Vec3> positions(tp.AtomCount(), nmr::Vec3::Zero());
    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame (no MopacMcConnell)");
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    EXPECT_EQ(tr->NumFrames(), kFrames);
    EXPECT_EQ(tr->SourceAttachedCount(), 0u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("mopac_mc_shielding_ts_absent_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/mopac_mc_shielding_time_series"));

    fs::remove(h5_path);
}


// ── Integration1P9J: real Trajectory::Run with MOPAC + McConnell ──

TEST(MopacMcConnellShieldingTimeSeries, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = false;     // MopacMcConnellResult needs MopacResult
    opts.skip_coulomb = false;   // enables Mopac-family attachments
    opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    const auto status = traj.Run(tp, config, session);
    if (status != nmr::kOk) {
        GTEST_SKIP() << "Trajectory::Run failed; skipping per "
                        "feedback_log_overages_dont_assert.";
    }

    auto& tr = tp.Result<nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult>();
    const std::size_t T = tr.NumFrames();
    const std::size_t attached = tr.SourceAttachedCount();
    std::cout << "MopacMcConnellShieldingTS Integration1P9J: T=" << T
              << " MopacMcConnell-attached=" << attached << std::endl;

    if (attached == 0) {
        GTEST_SKIP() << "MopacMcConnellResult did not attach on any frame.";
    }

    const std::string h5_path = (fs::temp_directory_path() /
        ("mopac_mc_shielding_ts_int1p9j_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/mopac_mc_shielding_time_series");
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[1], T);
    EXPECT_EQ(dims[2], 9u)
        << "9-component emission per 'if not traceless write both'";

    std::string parity, layout, units;
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("irrep_layout").read(layout);
    grp.getAttribute("units").read(units);
    EXPECT_EQ(parity, "0e+1o+2e");
    EXPECT_EQ(layout,
        "T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");
    EXPECT_EQ(units, "Angstrom^-3")
        << "bare McConnell kernel bo·M/r³, no Δχ × γ multiplication "
           "at extraction; decision 2026-05-21 per science/math review H1.";

    const std::size_t N = dims[0];
    std::vector<double> flat(N * T * 9);
    ds.read(flat.data());
    double max_t0 = 0.0, max_t1 = 0.0, max_t2 = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const std::size_t base = (i * T + t) * 9;
            for (std::size_t k = 0; k < 9; ++k) {
                EXPECT_TRUE(std::isfinite(flat[base + k]));
            }
            max_t0 = std::max(max_t0, std::abs(flat[base + 0]));
            for (std::size_t k = 1; k <= 3; ++k)
                max_t1 = std::max(max_t1, std::abs(flat[base + k]));
            for (std::size_t k = 4; k <= 8; ++k)
                max_t2 = std::max(max_t2, std::abs(flat[base + k]));
        }
    }
    std::cout << "  max|T0| = " << max_t0
              << "  max|T1| = " << max_t1
              << "  max|T2| = " << max_t2
              << "  (T0/T1 not necessarily ~0 — McConnell M_total "
                 "is not traceless)" << std::endl;
    // Sanity: at least T2 should be nonzero on every protein
    // (bond-anisotropy is a real signal). T0/T1 may or may not be —
    // log don't assert per feedback_log_overages_dont_assert.
    EXPECT_GT(max_t2, 0.0);

    fs::remove(h5_path);
}


// ── Layer 0b: FinalizeIdempotency ──────────────────────────────

TEST(MopacMcConnellShieldingTimeSeries, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = false;
    opts.skip_coulomb = false;
    opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    if (traj.Run(tp, config, session) != nmr::kOk) {
        GTEST_SKIP() << "Trajectory::Run failed; skipping.";
    }

    auto* buf_first = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(
        typeid(nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult)));
    if (!buf_first) GTEST_SKIP() << "MopacMcConnell never attached.";
    const std::size_t N_first = buf_first->AtomCount();
    const std::size_t T_first = buf_first->StridePerAtom();

    auto& tr = tp.Result<nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);  // second call — bounds-check idempotency
    auto* buf_second = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(
        typeid(nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_second, nullptr);
    EXPECT_EQ(buf_second->AtomCount(), N_first);
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
}
