//
// Smoke tests for the AIMNet2 fleet TR pair:
//   - AIMNet2EmbeddingTimeSeriesTrajectoryResult (256-dim float32)
//   - AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult (Vec3 + scalar)
//
// Synthetic-driven (no model load required). Verifies shape, attrs,
// round-trip, and per-atom values via the public WriteH5Group surface.
// Integration-on-1P9J path is gated by the AIMNet2 model and skipped
// when unavailable (same pattern as test_aimnet2_charge_time_series).
//

#include "AIMNet2EmbeddingTimeSeriesTrajectoryResult.h"
#include "AIMNet2ChargeResponseGradientResult.h"
#include "AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult.h"
#include "AIMNet2ChargeResponseGradientWelfordTrajectoryResult.h"
#include "AIMNet2Result.h"
#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "GeometryResult.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "SpatialIndexResult.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryAtom.h"
#include "TrajectoryProtein.h"
#include "TrajectoryResult.h"
#include "Types.h"

#include <gtest/gtest.h>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cmath>
#include <filesystem>
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


// ============================================================================
// Embedding TS (256-dim float32)
// ============================================================================

TEST(AIMNet2EmbeddingTimeSeries, SyntheticThreeFramesH5RoundTrip) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::AIMNet2EmbeddingTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr std::size_t kFrames = 3;
    constexpr std::size_t kDim    = nmr::AIMNET2_AIM_DIMS;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        for (std::size_t i = 0; i < N; ++i) {
            auto& vec = conf->MutableAtomAt(i).aimnet2_aim;
            for (std::size_t d = 0; d < kDim; ++d) {
                vec[d] = static_cast<float>(0.0001 * static_cast<double>(i)
                                            + 0.01 * static_cast<double>(t)
                                            + 0.001 * static_cast<double>(d));
            }
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);
    EXPECT_EQ(tr->NumFrames(), kFrames);

    const std::string h5_path = (fs::temp_directory_path() /
        ("aimnet2_embedding_ts_unit_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/aimnet2_embedding_time_series"));
    auto grp = reopen.getGroup("/trajectory/aimnet2_embedding_time_series");

    auto ds = grp.getDataSet("embedding");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[0], N);
    EXPECT_EQ(dims[1], kFrames);
    EXPECT_EQ(dims[2], kDim);

    // Spot-check a single (atom, frame, dim) cell. In synthetic mode the
    // ConformationResult is NOT attached, so the TR's HasResult gate
    // correctly emits the NaN placeholder ("absent, not faked",
    // codex review 2026-05-20); synthetic aimnet2_aim values on
    // MutableAtom are ignored.
    std::vector<float> buf(N * kFrames * kDim);
    ds.read(buf.data());
    const std::size_t i = N / 2;
    const std::size_t t = 1;
    const std::size_t d = 5;
    EXPECT_TRUE(std::isnan(buf[(i * kFrames + t) * kDim + d]));

    // Attr checks
    std::string source, policy, irrep, parity;
    bool optional_large = false;
    std::size_t embedding_dim = 0;
    grp.getAttribute("source").read(source);
    grp.getAttribute("source_attached_policy").read(policy);
    grp.getAttribute("optional_large").read(optional_large);
    grp.getAttribute("embedding_dim").read(embedding_dim);
    grp.getAttribute("irrep_layout").read(irrep);
    grp.getAttribute("parity").read(parity);
    EXPECT_NE(source.find("AIMNet2Result.aimnet2_aim"), std::string::npos);
    EXPECT_EQ(policy, "always_attached");
    EXPECT_TRUE(optional_large);
    EXPECT_EQ(embedding_dim, kDim);
    EXPECT_EQ(irrep, "feature_vector");
    EXPECT_EQ(parity, "0e");

    // source_attached_per_frame canonical mask (SDK contract). In this
    // synthetic-driven test the conformation has NO ConformationResult
    // attached, so the TR's HasResult<...>() gate correctly records
    // mask=0 for every frame. Production runs (where AIMNet2Result /
    // AIMNet2ChargeResponseGradientResult ARE attached) land mask=1. This
    // asserts the gate logic works as advertised.
    std::vector<std::uint8_t> mask;
    grp.getDataSet("source_attached_per_frame").read(mask);
    ASSERT_EQ(mask.size(), kFrames);
    for (auto m : mask) EXPECT_EQ(m, 0u);

    fs::remove(h5_path);
}


TEST(AIMNet2EmbeddingTimeSeries, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::AIMNet2EmbeddingTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr std::size_t kFrames = 3;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic");
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);
    const std::size_t T_first = tr->NumFrames();
    tr->Finalize(tp, traj);  // second call must not corrupt state
    EXPECT_EQ(tr->NumFrames(), T_first);
}


// ============================================================================
// Polarisability TS (Vec3 + scalar)
// ============================================================================

TEST(AIMNet2ChargeResponseGradientTimeSeries, SyntheticThreeFramesH5RoundTrip) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr std::size_t kFrames = 3;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        for (std::size_t i = 0; i < N; ++i) {
            auto& ca = conf->MutableAtomAt(i);
            ca.aimnet2_charge_response_gradient_vector =
                nmr::Vec3(0.01 * static_cast<double>(i),
                          0.02 * static_cast<double>(t),
                          0.03 * static_cast<double>(i + t));
            ca.aimnet2_charge_response_gradient_scalar =
                ca.aimnet2_charge_response_gradient_vector.norm();
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);
    EXPECT_EQ(tr->NumFrames(), kFrames);

    const std::string h5_path = (fs::temp_directory_path() /
        ("aimnet2_charge_response_gradient_ts_unit_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/aimnet2_charge_response_gradient_time_series"));
    auto grp = reopen.getGroup("/trajectory/aimnet2_charge_response_gradient_time_series");

    auto ds_vec = grp.getDataSet("charge_response_gradient_vector");
    auto ds_scl = grp.getDataSet("charge_response_gradient_scalar");
    const auto vdims = ds_vec.getSpace().getDimensions();
    const auto sdims = ds_scl.getSpace().getDimensions();
    ASSERT_EQ(vdims.size(), 3u);
    EXPECT_EQ(vdims[0], N);
    EXPECT_EQ(vdims[1], kFrames);
    EXPECT_EQ(vdims[2], std::size_t(3));
    ASSERT_EQ(sdims.size(), 2u);
    EXPECT_EQ(sdims[0], N);
    EXPECT_EQ(sdims[1], kFrames);

    // Synthetic mode: no ConformationResult attached, so the gate
    // emits NaN placeholders for both vector and scalar ("absent,
    // not faked", codex review 2026-05-20).
    std::vector<double> vbuf(N * kFrames * 3);
    std::vector<double> sbuf(N * kFrames);
    ds_vec.read(vbuf.data());
    ds_scl.read(sbuf.data());
    const std::size_t i = N / 2;
    const std::size_t t = 1;
    EXPECT_TRUE(std::isnan(vbuf[(i * kFrames + t) * 3 + 0]));
    EXPECT_TRUE(std::isnan(vbuf[(i * kFrames + t) * 3 + 1]));
    EXPECT_TRUE(std::isnan(vbuf[(i * kFrames + t) * 3 + 2]));
    EXPECT_TRUE(std::isnan(sbuf[i * kFrames + t]));

    // Attr checks. Vec3 metadata follows existing TR convention:
    // layout + normalization + parity emitted as separate attrs.
    std::string source, policy, uv, us, ilv, norm_v, plv, ils, pls;
    grp.getAttribute("source").read(source);
    grp.getAttribute("source_attached_policy").read(policy);
    grp.getAttribute("units_vector").read(uv);
    grp.getAttribute("units_scalar").read(us);
    grp.getAttribute("irrep_layout_vector").read(ilv);
    grp.getAttribute("normalization_vector").read(norm_v);
    grp.getAttribute("parity_vector").read(plv);
    grp.getAttribute("irrep_layout_scalar").read(ils);
    grp.getAttribute("parity_scalar").read(pls);
    EXPECT_NE(source.find("AIMNet2ChargeResponseGradientResult"), std::string::npos);
    EXPECT_EQ(policy, "always_attached");
    EXPECT_EQ(uv, "e^2/Å");
    EXPECT_EQ(us, "e^2/Å");
    EXPECT_EQ(ilv, "x,y,z");
    EXPECT_EQ(norm_v, "cartesian");
    EXPECT_EQ(plv, "1o");
    EXPECT_EQ(ils, "T0");
    EXPECT_EQ(pls, "0e");

    // source_attached_per_frame canonical mask (SDK contract). In this
    // synthetic-driven test the conformation has NO ConformationResult
    // attached, so the TR's HasResult<...>() gate correctly records
    // mask=0 for every frame. Production runs (where AIMNet2Result /
    // AIMNet2ChargeResponseGradientResult ARE attached) land mask=1. This
    // asserts the gate logic works as advertised.
    std::vector<std::uint8_t> mask;
    grp.getDataSet("source_attached_per_frame").read(mask);
    ASSERT_EQ(mask.size(), kFrames);
    for (auto m : mask) EXPECT_EQ(m, 0u);

    fs::remove(h5_path);
}


TEST(AIMNet2ChargeResponseGradientTimeSeries, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr std::size_t kFrames = 3;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic");
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);
    const std::size_t T_first = tr->NumFrames();
    tr->Finalize(tp, traj);
    EXPECT_EQ(tr->NumFrames(), T_first);
}


// ============================================================================
// Polarisability Welford (TR #3 — AV companion)
// ============================================================================

TEST(AIMNet2ChargeResponseGradientWelford, SyntheticThreeFramesSkipsGroupOnAbsence) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::AIMNet2ChargeResponseGradientWelfordTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr std::size_t kFrames = 3;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic");
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);
    EXPECT_EQ(tr->NumFrames(), kFrames);
    // Synthetic mode = no ConformationResult attached → 0 source-attached
    // frames → WriteH5Group skips the group entirely (Welford accumulator
    // contract: no frames means no honest mean/M2).
    EXPECT_EQ(tr->SourceAttachedCount(), 0u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("aimnet2_charge_response_gradient_welford_skipped_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/aimnet2_charge_response_gradient_welford"));
    fs::remove(h5_path);
}


// ============================================================================
// INTEGRATION: real AIMNet2 CRG kernel through Trajectory::Run, multi-frame.
// Codex F4 2026-05-20 — exercises the actual Welford accumulation path
// against per-frame AIMNet2 backward-pass gradients. Skips when CUDA /
// model unavailable. The synthetic-positions accumulation path cannot
// be exercised here because TrajectoryProtein::Seed FATALs on all-zero
// positions (substrate canonicalisation); see sibling
// test_hydration_geometry_welford coverage-gap note for the same
// pattern.
// ============================================================================

TEST(AIMNet2ChargeResponseGradientWelford, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    const std::string& model_path = nmr::test::TestEnvironment::Aimnet2Model();
    if (model_path.empty() || !fs::exists(model_path)) {
        GTEST_SKIP() << "AIMNet2 model not available";
    }

    nmr::Session session;
    ASSERT_EQ(session.LoadAimnet2Model(model_path), nmr::kOk)
        << session.LastError();

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true;
    opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::AIMNet2Result));
    config.RequireConformationResult(typeid(nmr::AIMNet2ChargeResponseGradientResult));
    config.SetRequiresAimnet2(true);
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
                return nmr::AIMNet2ChargeResponseGradientWelfordTrajectoryResult::Create(tp_in);
            });
    config.SetStride(500);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    const std::size_t F = traj.FrameCount();
    EXPECT_GE(F, 2u);

    // Per-atom Welford state populated by the live accumulation path.
    // Verifies the F1 non-finite guard didn't reject all frames AND the
    // canonical mean/m2/std are emitted at the canonical shape with
    // physically-sensible magnitudes.
    std::size_t populated = 0;
    double max_abs_mean = 0.0;
    for (std::size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& ws = tp.AtomAt(i).aimnet2_charge_response_gradient_welford;
        EXPECT_EQ(ws.n_frames, F);
        EXPECT_TRUE(std::isfinite(ws.charge_response_gradient_scalar.mean));
        EXPECT_TRUE(std::isfinite(ws.charge_response_gradient_scalar.std));
        for (std::size_t c = 0; c < 3; ++c) {
            EXPECT_TRUE(std::isfinite(ws.charge_response_gradient_vector[c].mean));
            EXPECT_TRUE(std::isfinite(ws.charge_response_gradient_vector[c].std));
            max_abs_mean = std::max(max_abs_mean,
                std::abs(ws.charge_response_gradient_vector[c].mean));
        }
        if (ws.charge_response_gradient_scalar.mean > 1e-12) ++populated;
    }
    EXPECT_GT(populated, tp.AtomCount() / 2)
        << "Most atoms should have nonzero charge-response gradient L2 norm";
    EXPECT_GT(max_abs_mean, 1e-6) << "Welford accumulation at noise floor";

    // H5 write-back coverage is provided by the synthetic-skip test
    // (covers attr/dataset layout when source is absent) + the per-atom
    // state check above (covers the live accumulation path on real data
    // through Trajectory::Run). Direct H5 readback from this test would
    // require capturing the run-attached TR by type_index from
    // TrajectoryProtein, which the existing API doesn't expose for this
    // result family — that surface is added in a later commit if a
    // downstream integration test needs to round-trip H5 from the
    // production accumulation path.
    std::cout << "AIMNet2ChargeResponseGradientWelford max|grad_x|="
              << max_abs_mean << " e²/Å, populated=" << populated
              << "/" << tp.AtomCount() << "\n";
}
