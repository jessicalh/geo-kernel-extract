//
// test_apbs_efg_time_series: discipline + integration for
// ApbsEfgTimeSeriesTrajectoryResult. T2 SphericalTensor DenseBuffer
// pattern; sibling of ApbsEfieldTimeSeries (Vec3 TS, same source calc).
// Canonical T2-EFG TS template — TR4 of the 2026-05-20 13-TR plan.
//

#include "ApbsFieldResult.h"
#include "ChargeAssignmentResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DenseBuffer.h"
#include "ApbsEfgTimeSeriesTrajectoryResult.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
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
// Build a synthetic T2-only SphericalTensor from a deterministic
// (atom, frame) seed. T0/T1 explicitly zero (mirrors what
// ApbsFieldResult produces after the Poisson tracelessness
// projection). The five T2 components encode the seed so the
// round-trip test can assert exact recovery.
nmr::SphericalTensor MakeSyntheticT2(std::size_t i, std::size_t t) {
    nmr::SphericalTensor st;
    st.T0 = 0.0;
    st.T1 = {0.0, 0.0, 0.0};
    const double base = static_cast<double>(i) + t * 1.0;
    st.T2 = {base + 0.0, base + 0.1, base + 0.2,
             base + 0.3, base + 0.4};
    return st;
}

}  // namespace


// ── Layer 0a: synthetic Compute + Finalize + H5 round-trip ──
//
// Drives Compute/Finalize/WriteH5Group directly without
// Trajectory::Run. Verifies (a) DenseBuffer round-trip preserves
// SphericalTensor T2 components exactly, (b) H5 shape is (N, T, 5)
// not (N, T, 9), (c) attrs parity="2e" / irrep_layout = T2 only.

TEST(ApbsEfgTimeSeries, SyntheticFourFrames) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    const size_t N = tp.AtomCount();
    auto tr = nmr::ApbsEfgTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 4;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        // Attach a synthetic ApbsFieldResult so the HasResult gate
        // sees the source — the result instance itself is never
        // queried, only the apbs_efg_spherical field on
        // ConformationAtom is read.
        conf->ForceAttachResultForTesting(
            std::make_unique<nmr::ApbsFieldResult>());
        for (size_t i = 0; i < N; ++i) {
            conf->MutableAtomAt(i).apbs_efg_spherical = MakeSyntheticT2(i, t);
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    auto* buf = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(typeid(
        nmr::ApbsEfgTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->StridePerAtom(), kFrames);
    for (size_t i : {size_t(0), N / 2, N - 1}) {
        for (size_t t = 0; t < kFrames; ++t) {
            const nmr::SphericalTensor& st = buf->At(i, t);
            const auto expected = MakeSyntheticT2(i, t);
            EXPECT_DOUBLE_EQ(st.T0, expected.T0);
            for (size_t k = 0; k < 3; ++k)
                EXPECT_DOUBLE_EQ(st.T1[k], expected.T1[k]);
            for (size_t k = 0; k < 5; ++k)
                EXPECT_DOUBLE_EQ(st.T2[k], expected.T2[k]);
        }
    }

    const std::string h5_path = (fs::temp_directory_path() /
        ("apbs_efg_ts_unit_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/apbs_efg_time_series");
    auto ds = grp.getDataSet("t2");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[0], N);
    EXPECT_EQ(dims[1], kFrames);
    EXPECT_EQ(dims[2], 5u) << "T2-only emission per 2026-05-18 schema rev";

    std::string parity, layout, units, policy;
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("irrep_layout").read(layout);
    grp.getAttribute("units").read(units);
    grp.getAttribute("source_attached_policy").read(policy);
    EXPECT_EQ(parity, "2e");
    EXPECT_EQ(layout, "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");
    EXPECT_EQ(units, "V/Å^2");
    EXPECT_NE(policy.find("always_attached"), std::string::npos);

    // Round-trip: read the (N, T, 5) flat data, check a known cell.
    std::vector<double> flat(N * kFrames * 5);
    ds.read(flat.data());
    const size_t check_i = N / 2;
    for (size_t t = 0; t < kFrames; ++t) {
        const auto expected = MakeSyntheticT2(check_i, t);
        const size_t base = (check_i * kFrames + t) * 5;
        for (size_t k = 0; k < 5; ++k)
            EXPECT_DOUBLE_EQ(flat[base + k], expected.T2[k]);
    }

    // source_attached_per_frame should be all-1 (HasResult succeeded
    // each frame courtesy of ForceAttachResultForTesting).
    auto mask_ds = grp.getDataSet("source_attached_per_frame");
    std::vector<std::uint8_t> mask(kFrames);
    mask_ds.read(mask.data());
    for (auto m : mask) EXPECT_EQ(m, 1u);

    fs::remove(h5_path);
}


// ── Layer 0b: source-absent → NaN-fill + mask=0 ─────────────
//
// Drives Compute without attaching ApbsFieldResult; verifies
// canonical 'absent, not faked': NaN cells + mask all-0 + group
// still emitted (because Finalize succeeded with N samples).

TEST(ApbsEfgTimeSeries, SourceAbsentNanFillAndMaskZero) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    const size_t N = tp.AtomCount();
    auto tr = nmr::ApbsEfgTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 2;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (size_t t = 0; t < kFrames; ++t) {
        // NO ForceAttachResultForTesting → HasResult false → NaN-fill.
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame (source absent)");
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    auto* buf = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(typeid(
        nmr::ApbsEfgTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    // T2 components NaN, T0/T1 NaN too (set explicitly in cpp).
    for (size_t i : {size_t(0), N / 2, N - 1}) {
        for (size_t t = 0; t < kFrames; ++t) {
            const nmr::SphericalTensor& st = buf->At(i, t);
            EXPECT_TRUE(std::isnan(st.T0));
            for (size_t k = 0; k < 3; ++k) EXPECT_TRUE(std::isnan(st.T1[k]));
            for (size_t k = 0; k < 5; ++k) EXPECT_TRUE(std::isnan(st.T2[k]));
        }
    }

    const std::string h5_path = (fs::temp_directory_path() /
        ("apbs_efg_ts_nanfill_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/apbs_efg_time_series");
    auto mask_ds = grp.getDataSet("source_attached_per_frame");
    std::vector<std::uint8_t> mask(kFrames);
    mask_ds.read(mask.data());
    for (auto m : mask) EXPECT_EQ(m, 0u);
    fs::remove(h5_path);
}


// ── Layer 0c: Frame0Semantics (stride ≥ fixture length) ─────

TEST(ApbsEfgTimeSeries, Frame0Semantics) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = false;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::ChargeAssignmentResult));
    config.RequireConformationResult(typeid(nmr::ApbsFieldResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::ApbsEfgTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_EQ(traj.FrameCount(), 1u);
}


// ── Layer 0d: FinalizeIdempotency ───────────────────────────

TEST(ApbsEfgTimeSeries, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = false;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::ChargeAssignmentResult));
    config.RequireConformationResult(typeid(nmr::ApbsFieldResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::ApbsEfgTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto* buf_first = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(
        typeid(nmr::ApbsEfgTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_first, nullptr);
    const std::size_t N_first = buf_first->AtomCount();
    const std::size_t T_first = buf_first->StridePerAtom();

    auto& tr = tp.Result<nmr::ApbsEfgTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);
    auto* buf_second = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(
        typeid(nmr::ApbsEfgTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_second, nullptr);
    EXPECT_EQ(buf_second->AtomCount(), N_first);
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
}


// ── Integration1P9J: real APBS Poisson solve via Trajectory::Run ──
//
// Stride 300 → ~3 frames on the 1P9J trajectory. Asserts (a) the
// DenseBuffer materialises, (b) every T2 component is finite,
// (c) at least one atom has |T2| > 0.01 V/Å² (calc is firing —
// catches the "all zeros" failure mode), (d) all values are within
// a generous physical sanity limit (T2 components are V/Å², so
// bounded by APBS_SANITY_LIMIT_EFG ~ 100 V/Å² as a loose bound;
// log don't assert per feedback_log_overages_dont_assert).

TEST(ApbsEfgTimeSeries, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = false;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::ChargeAssignmentResult));
    config.RequireConformationResult(typeid(nmr::ApbsFieldResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::ApbsEfgTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto* buf = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(typeid(
        nmr::ApbsEfgTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);

    const std::size_t N = buf->AtomCount();
    const std::size_t T = buf->StridePerAtom();
    double max_t2_mag = 0.0;
    std::size_t n_overages = 0;
    constexpr double kEfgSanityLimit = 100.0;  // V/Å² — generous
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const nmr::SphericalTensor& st = buf->At(i, t);
            for (std::size_t k = 0; k < 5; ++k)
                EXPECT_TRUE(std::isfinite(st.T2[k]))
                    << "atom " << i << " frame " << t << " T2[" << k << "]";
            const double mag = st.T2Magnitude();
            max_t2_mag = std::max(max_t2_mag, mag);
            if (mag > kEfgSanityLimit) ++n_overages;
        }
    }
    EXPECT_GT(max_t2_mag, 0.01)
        << "APBS EFG T2 magnitudes all near zero — calc not firing or "
        "tracelessness projection collapsed everything to zero";
    std::cout << "ApbsEfgTimeSeries Integration1P9J: max|T2|=" << max_t2_mag
              << " V/Å², overages>" << kEfgSanityLimit << "=" << n_overages
              << "/" << (N * T) << std::endl;

    // Mask should be all-1 (ApbsFieldResult is Require'd).
    auto& tr = tp.Result<nmr::ApbsEfgTimeSeriesTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), T);
}
