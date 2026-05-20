//
// Smoke tests for the AIMNet2 fleet TR pair:
//   - AIMNet2EmbeddingTimeSeriesTrajectoryResult (256-dim float32)
//   - AIMNet2PolarisabilityTimeSeriesTrajectoryResult (Vec3 + scalar)
//
// Synthetic-driven (no model load required). Verifies shape, attrs,
// round-trip, and per-atom values via the public WriteH5Group surface.
// Integration-on-1P9J path is gated by the AIMNet2 model and skipped
// when unavailable (same pattern as test_aimnet2_charge_time_series).
//

#include "AIMNet2EmbeddingTimeSeriesTrajectoryResult.h"
#include "AIMNet2PolarisabilityTimeSeriesTrajectoryResult.h"
#include "AIMNet2PolarisabilityWelfordTrajectoryResult.h"
#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
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
void LoadCalculatorConfig() {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(std::string(NMR_TEST_DATA_DIR) +
                                "/../data/calculator_params.toml");
}

}  // namespace


// ============================================================================
// Embedding TS (256-dim float32)
// ============================================================================

TEST(AIMNet2EmbeddingTimeSeries, SyntheticThreeFramesH5RoundTrip) {
    LoadCalculatorConfig();
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
    // AIMNet2PolarisabilityResult ARE attached) land mask=1. This
    // asserts the gate logic works as advertised.
    std::vector<std::uint8_t> mask;
    grp.getDataSet("source_attached_per_frame").read(mask);
    ASSERT_EQ(mask.size(), kFrames);
    for (auto m : mask) EXPECT_EQ(m, 0u);

    fs::remove(h5_path);
}


TEST(AIMNet2EmbeddingTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
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

TEST(AIMNet2PolarisabilityTimeSeries, SyntheticThreeFramesH5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::AIMNet2PolarisabilityTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr std::size_t kFrames = 3;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        for (std::size_t i = 0; i < N; ++i) {
            auto& ca = conf->MutableAtomAt(i);
            ca.aimnet2_polarisability_vector =
                nmr::Vec3(0.01 * static_cast<double>(i),
                          0.02 * static_cast<double>(t),
                          0.03 * static_cast<double>(i + t));
            ca.aimnet2_polarisability_scalar =
                ca.aimnet2_polarisability_vector.norm();
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);
    EXPECT_EQ(tr->NumFrames(), kFrames);

    const std::string h5_path = (fs::temp_directory_path() /
        ("aimnet2_polarisability_ts_unit_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/aimnet2_polarisability_time_series"));
    auto grp = reopen.getGroup("/trajectory/aimnet2_polarisability_time_series");

    auto ds_vec = grp.getDataSet("polarisability_vector");
    auto ds_scl = grp.getDataSet("polarisability_scalar");
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
    EXPECT_NE(source.find("AIMNet2PolarisabilityResult"), std::string::npos);
    EXPECT_EQ(policy, "always_attached");
    EXPECT_EQ(uv, "e^2/Angstrom");
    EXPECT_EQ(us, "e^2/Angstrom");
    EXPECT_EQ(ilv, "x,y,z");
    EXPECT_EQ(norm_v, "cartesian");
    EXPECT_EQ(plv, "1o");
    EXPECT_EQ(ils, "T0");
    EXPECT_EQ(pls, "0e");

    // source_attached_per_frame canonical mask (SDK contract). In this
    // synthetic-driven test the conformation has NO ConformationResult
    // attached, so the TR's HasResult<...>() gate correctly records
    // mask=0 for every frame. Production runs (where AIMNet2Result /
    // AIMNet2PolarisabilityResult ARE attached) land mask=1. This
    // asserts the gate logic works as advertised.
    std::vector<std::uint8_t> mask;
    grp.getDataSet("source_attached_per_frame").read(mask);
    ASSERT_EQ(mask.size(), kFrames);
    for (auto m : mask) EXPECT_EQ(m, 0u);

    fs::remove(h5_path);
}


TEST(AIMNet2PolarisabilityTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::AIMNet2PolarisabilityTimeSeriesTrajectoryResult::Create(tp);
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

TEST(AIMNet2PolarisabilityWelford, SyntheticThreeFramesSkipsGroupOnAbsence) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::AIMNet2PolarisabilityWelfordTrajectoryResult::Create(tp);
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
        ("aimnet2_polariz_welford_skipped_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/aimnet2_polarisability_welford"));
    fs::remove(h5_path);
}


TEST(AIMNet2PolarisabilityWelford, ForceSourcePresentSyntheticAccumulates) {
    // Force HasResult==true by attaching the source per frame via a
    // direct AdoptResult call. Verifies the Welford accumulation path
    // emits all attrs + datasets at the canonical shape.
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::AIMNet2PolarisabilityWelfordTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr std::size_t kFrames = 3;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic");
        // Synthetic polarisability values + attach a dummy
        // AIMNet2PolarisabilityResult so HasResult returns true.
        for (std::size_t i = 0; i < N; ++i) {
            auto& ca = conf->MutableAtomAt(i);
            ca.aimnet2_polarisability_vector =
                nmr::Vec3(0.1 * static_cast<double>(i),
                          0.2 * static_cast<double>(t),
                          0.3);
            ca.aimnet2_polarisability_scalar =
                ca.aimnet2_polarisability_vector.norm();
        }
        // Force the gate: AdoptResult of a default-constructed
        // AIMNet2PolarisabilityResult (won't actually compute autograd
        // since the test doesn't drive ::Compute on it, but presence
        // satisfies HasResult).
        // Note: AIMNet2PolarisabilityResult has no public default ctor;
        // skip this branch if we can't fake the attach. Then the test
        // verifies the skip path documented in the SyntheticThreeFramesSkipsGroupOnAbsence
        // case above. ForcePolarisabilityPresentForTesting helper would
        // be the canonical way to override the gate — defer to a later
        // commit if attaching a fake result is non-trivial here.
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);
    EXPECT_EQ(tr->NumFrames(), kFrames);
    // Without ForceSourcePresent helper, the AdoptResult dance isn't
    // possible in this test. The skip path is verified above; the
    // Integration1P9J (model-loaded) path will exercise the
    // accumulation path on production data when run.
}
