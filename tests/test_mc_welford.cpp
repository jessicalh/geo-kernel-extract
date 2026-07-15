//
// test_mc_welford: discipline + integration for McConnellWelfordTrajectoryResult.
// AV-pattern Welford on the McConnell shielding kernel. Mirrors the BS
// Welford test shape — no synthetic test; the shared MomentsUpdate math
// is exercised by every AV TR. Discipline trio + Integration1P9J.
//

#include "McConnellResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "McConnellWelfordTrajectoryResult.h"
#include "TestConfig.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryAtom.h"
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


// ============================================================================
// DISCIPLINE: Frame-0 semantics.
// ============================================================================

TEST(McConnellWelford, Frame0Semantics) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "McConnellWelfordFrame0Semantics", 99999);
    config.RequireConformationResult(typeid(nmr::McConnellResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::McConnellWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    ASSERT_EQ(traj.FrameCount(), 1u);

    ASSERT_TRUE(tp.HasResult<nmr::McConnellWelfordTrajectoryResult>());
    const auto& tr = tp.Result<nmr::McConnellWelfordTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);

    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& ta = tp.AtomAt(i);
        if (ta.mc_welford.n_frames == 0u) {
            EXPECT_TRUE(std::isnan(ta.mc_welford.t0.mean));
            EXPECT_TRUE(std::isnan(ta.mc_welford.t0.std));
        } else {
            EXPECT_EQ(ta.mc_welford.n_frames, 1u);
            EXPECT_DOUBLE_EQ(ta.mc_welford.t0.std, 0.0);
        }
        EXPECT_EQ(ta.mc_welford.delta_n, 0u);
    }
}


// ============================================================================
// DISCIPLINE: Finalize idempotency.
// ============================================================================

TEST(McConnellWelford, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "McConnellWelfordFinalizeIdempotency", 99999);
    config.RequireConformationResult(typeid(nmr::McConnellResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::McConnellWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::McConnellWelfordTrajectoryResult>();
    const size_t probe = tp.AtomCount() / 2;
    const double mean_first = tp.AtomAt(probe).mc_welford.t0.mean;
    const double std_first  = tp.AtomAt(probe).mc_welford.t0.std;

    tr.Finalize(tp, traj);

    if (std::isnan(mean_first)) {
        EXPECT_TRUE(std::isnan(tp.AtomAt(probe).mc_welford.t0.mean));
        EXPECT_TRUE(std::isnan(tp.AtomAt(probe).mc_welford.t0.std));
    } else {
        EXPECT_DOUBLE_EQ(tp.AtomAt(probe).mc_welford.t0.mean, mean_first);
        EXPECT_DOUBLE_EQ(tp.AtomAt(probe).mc_welford.t0.std,  std_first);
    }
}


TEST(McConnellWelford,
     UnevaluableFrameIsOmittedAndDoesNotBridgeDeltas) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto build = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto tp_owner = nmr::TrajectoryProtein::CreateForTesting(
        std::move(build.protein));
    ASSERT_NE(tp_owner, nullptr);
    auto& tp = *tp_owner;
    ASSERT_GT(tp.AtomCount(), 1u);

    auto tr = nmr::McConnellWelfordTrajectoryResult::Create(tp);
    nmr::Trajectory traj({}, {}, {});
    std::vector<nmr::Vec3> positions(tp.AtomCount(), nmr::Vec3::Zero());
    auto conf = std::make_unique<nmr::ProteinConformation>(
        &tp.ProteinRef(), positions, "conditional McConnell Welford");

    nmr::SphericalTensor first;
    first.T0 = 1.0;
    first.T1 = {0.1, 0.2, 0.3};
    first.T2 = {0.4, 0.5, 0.6, 0.7, 0.8};
    conf->MutableAtomAt(0).mc_shielding_contribution = first;
    tr->Compute(*conf, tp, traj, 10, 0.0);

    nmr::SphericalTensor unavailable;
    unavailable.T0 = std::nan("");
    unavailable.T1.fill(std::nan(""));
    unavailable.T2.fill(std::nan(""));
    conf->MutableAtomAt(0).mc_shielding_contribution = unavailable;
    tr->Compute(*conf, tp, traj, 11, 1.0);

    nmr::SphericalTensor third = first;
    third.T0 = 3.0;
    third.T1 = {0.3, 0.4, 0.5};
    third.T2 = {0.6, 0.7, 0.8, 0.9, 1.0};
    conf->MutableAtomAt(0).mc_shielding_contribution = third;
    tr->Compute(*conf, tp, traj, 12, 2.0);
    tr->Finalize(tp, traj);

    const auto& w = tp.AtomAt(0).mc_welford;
    EXPECT_EQ(tr->NumFrames(), 3u);
    EXPECT_EQ(w.n_frames, 2u);
    EXPECT_DOUBLE_EQ(w.t0.mean, 2.0);
    EXPECT_NEAR(w.t0.std, std::sqrt(2.0), 1e-15);
    EXPECT_DOUBLE_EQ(w.t0.min, 1.0);
    EXPECT_DOUBLE_EQ(w.t0.max, 3.0);
    EXPECT_EQ(w.t0.min_frame, 10u);
    EXPECT_EQ(w.t0.max_frame, 12u);
    EXPECT_EQ(w.delta_n, 0u);
    EXPECT_EQ(w.dxdt_n, 0u);
    EXPECT_TRUE(std::isnan(w.t0_delta.mean));

    for (std::size_t i = 1; i < tp.AtomCount(); ++i) {
        EXPECT_EQ(tp.AtomAt(i).mc_welford.n_frames, 3u);
        EXPECT_DOUBLE_EQ(tp.AtomAt(i).mc_welford.t0.mean, 0.0);
    }
}


TEST(McConnellWelford, H5DirectionalMetadataZeroCountSynthetic) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto build = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto tp_owner = nmr::TrajectoryProtein::CreateForTesting(
        std::move(build.protein));
    ASSERT_NE(tp_owner, nullptr);
    auto& tp = *tp_owner;
    ASSERT_GT(tp.AtomCount(), 1u);

    auto tr = nmr::McConnellWelfordTrajectoryResult::Create(tp);
    nmr::Trajectory traj({}, {}, {});
    tr->Compute(tp.CanonicalConformation(), tp, traj, 7, 2.5);
    tp.MutableAtomAt(0).mc_welford = {};
    tr->Finalize(tp, traj);

    const std::string h5_path = (fs::temp_directory_path() /
        ("mc_welford_directional_metadata_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr->WriteH5Group(tp, file);
    }
    {
    HighFive::File file(h5_path, HighFive::File::ReadOnly);
    auto grp = file.getGroup("/trajectory/mc_welford");

    std::string irrep_scope, mean_law, component_law, validity,
                zero_count_validity;
    grp.getAttribute("irrep_metadata_scope").read(irrep_scope);
    grp.getAttribute("directional_mean_transformation").read(mean_law);
    grp.getAttribute("componentwise_statistic_transformation").read(component_law);
    grp.getAttribute("validity").read(validity);
    grp.getAttribute("zero_count_sentinel_validity").read(zero_count_validity);
    EXPECT_EQ(irrep_scope,
              "only assembled component means carry directional irrep metadata");
    EXPECT_EQ(mean_law,
              "t0_mean is invariant; assembled t1_mean is axial: "
              "a'=det(R) R a; assembled t2_mean is even rank-2: T'=R T R^T");
    EXPECT_EQ(component_law,
              "componentwise m2,std,min,max,min_frame,max_frame have no closed "
              "irrep transformation law");
    EXPECT_EQ(validity,
              "moments condition on complete finite McConnell tensors; "
              "n_frames_per_atom is the evaluable-sample count, and physical "
              "zero means an evaluable empty or cancelled fixed-channel sum");
    EXPECT_EQ(zero_count_validity,
              "when n_frames_per_atom=0, mean,m2,std are NaN and min=+inf,"
              "max=-inf,min_frame=0,max_frame=0 are invalid sentinels; "
              "n_frames_per_atom gates validity");

    std::vector<std::size_t> counts;
    grp.getDataSet("n_frames_per_atom").read(counts);
    ASSERT_EQ(counts.size(), tp.AtomCount());
    EXPECT_EQ(counts[0], 0u);
    for (std::size_t i = 1; i < counts.size(); ++i) EXPECT_EQ(counts[i], 1u);

    auto expect_zero_count_components = [&](const std::string& prefix,
                                             std::size_t width) {
        const std::size_t size = tp.AtomCount() * width;
        std::vector<double> mean(size), m2(size), stddev(size), min(size),
                            max(size);
        std::vector<std::size_t> min_frame(size), max_frame(size);
        grp.getDataSet(prefix + "_mean").read(mean.data());
        grp.getDataSet(prefix + "_m2").read(m2.data());
        grp.getDataSet(prefix + "_std").read(stddev.data());
        grp.getDataSet(prefix + "_min").read(min.data());
        grp.getDataSet(prefix + "_max").read(max.data());
        grp.getDataSet(prefix + "_min_frame").read(min_frame.data());
        grp.getDataSet(prefix + "_max_frame").read(max_frame.data());
        ASSERT_GE(mean.size(), 2 * width);
        for (std::size_t k = 0; k < width; ++k) {
            EXPECT_TRUE(std::isnan(mean[k]));
            EXPECT_TRUE(std::isnan(m2[k]));
            EXPECT_TRUE(std::isnan(stddev[k]));
            EXPECT_TRUE(std::isinf(min[k]) && min[k] > 0.0);
            EXPECT_TRUE(std::isinf(max[k]) && max[k] < 0.0);
            EXPECT_EQ(min_frame[k], 0u);
            EXPECT_EQ(max_frame[k], 0u);

            const std::size_t valid = width + k;
            EXPECT_DOUBLE_EQ(mean[valid], 0.0);
            EXPECT_DOUBLE_EQ(m2[valid], 0.0);
            EXPECT_DOUBLE_EQ(stddev[valid], 0.0);
            EXPECT_DOUBLE_EQ(min[valid], 0.0);
            EXPECT_DOUBLE_EQ(max[valid], 0.0);
            EXPECT_EQ(min_frame[valid], 7u);
            EXPECT_EQ(max_frame[valid], 7u);
        }
    };
    expect_zero_count_components("t1", 3);
    expect_zero_count_components("t2", 5);

    }
    fs::remove(h5_path);
}


// ============================================================================
// DISCIPLINE: H5 round-trip.
// ============================================================================

TEST(McConnellWelford, H5RoundTrip) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "McConnellWelfordH5RoundTrip", 99999);
    config.RequireConformationResult(typeid(nmr::McConnellResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::McConnellWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::McConnellWelfordTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("mc_welford_h5_roundtrip_" + std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr.WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/mc_welford"));
    auto grp = reopen.getGroup("/trajectory/mc_welford");

    std::string units;
    grp.getAttribute("units").read(units);
    EXPECT_EQ(units, "Angstrom^-3");

    std::string irrep_scope, mean_law, component_law, zero_count_validity;
    grp.getAttribute("irrep_metadata_scope").read(irrep_scope);
    grp.getAttribute("directional_mean_transformation").read(mean_law);
    grp.getAttribute("componentwise_statistic_transformation").read(component_law);
    grp.getAttribute("zero_count_sentinel_validity").read(zero_count_validity);
    EXPECT_EQ(irrep_scope,
              "only assembled component means carry directional irrep metadata");
    EXPECT_EQ(mean_law,
              "t0_mean is invariant; assembled t1_mean is axial: "
              "a'=det(R) R a; assembled t2_mean is even rank-2: T'=R T R^T");
    EXPECT_EQ(component_law,
              "componentwise m2,std,min,max,min_frame,max_frame have no closed "
              "irrep transformation law");
    EXPECT_EQ(zero_count_validity,
              "when n_frames_per_atom=0, mean,m2,std are NaN and min=+inf,"
              "max=-inf,min_frame=0,max_frame=0 are invalid sentinels; "
              "n_frames_per_atom gates validity");

    const auto dims = grp.getDataSet("t0_mean").getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 1u);
    EXPECT_EQ(dims[0], tp.AtomCount());

    fs::remove(h5_path);
}


// ============================================================================
// INTEGRATION: real McConnell kernel through Trajectory::Run.
// ============================================================================

TEST(McConnellWelford, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = nmr::test::BuildTestConfig(
        nmr::test::TestProfile::KernelOnly, "McConnellWelfordIntegration", 300);
    config.RequireConformationResult(typeid(nmr::McConnellResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::McConnellWelfordTrajectoryResult::Create(tp_in);
        });

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_GE(traj.FrameCount(), 2u);

    // T0 / T1 / T2 channel coverage — protects against silent regression of
    // the 2026-05-17 PM McConnell T1 omission (PATTERNS Lesson 19: McConnell's
    // asymmetric three-term form has nonzero antisymmetric part, T1 ≠ 0). A
    // future refactor that quiets all three per-component T1 channels
    // would not trip the t0 assertion alone.
    size_t populated  = 0;
    size_t t1_populated = 0, t2_populated = 0;
    size_t unavailable = 0;
    double max_abs_t0 = 0.0, max_abs_t1 = 0.0, max_abs_t2 = 0.0;
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& ta = tp.AtomAt(i);
        EXPECT_LE(ta.mc_welford.n_frames, traj.FrameCount());
        if (ta.mc_welford.n_frames == 0u) {
            ++unavailable;
            EXPECT_TRUE(std::isnan(ta.mc_welford.t0.mean));
            EXPECT_TRUE(std::isnan(ta.mc_welford.t0.std));
            EXPECT_TRUE(std::isnan(ta.mc_welford.t2magnitude.mean));
            for (size_t k = 0; k < 3; ++k)
                EXPECT_TRUE(std::isnan(ta.mc_welford.t1[k].mean));
            for (size_t k = 0; k < 5; ++k)
                EXPECT_TRUE(std::isnan(ta.mc_welford.t2[k].mean));
            continue;
        }
        EXPECT_TRUE(std::isfinite(ta.mc_welford.t0.mean));
        EXPECT_TRUE(std::isfinite(ta.mc_welford.t0.std));
        EXPECT_TRUE(std::isfinite(ta.mc_welford.t2magnitude.mean));
        if (std::abs(ta.mc_welford.t0.mean) > 1e-12) ++populated;
        max_abs_t0 = std::max(max_abs_t0, std::abs(ta.mc_welford.t0.mean));

        bool atom_t1 = false;
        for (size_t k = 0; k < 3; ++k) {
            EXPECT_TRUE(std::isfinite(ta.mc_welford.t1[k].mean));
            const double v = std::abs(ta.mc_welford.t1[k].mean);
            max_abs_t1 = std::max(max_abs_t1, v);
            if (v > 1e-12) atom_t1 = true;
        }
        if (atom_t1) ++t1_populated;

        bool atom_t2 = false;
        for (size_t k = 0; k < 5; ++k) {
            EXPECT_TRUE(std::isfinite(ta.mc_welford.t2[k].mean));
            const double v = std::abs(ta.mc_welford.t2[k].mean);
            max_abs_t2 = std::max(max_abs_t2, v);
            if (v > 1e-12) atom_t2 = true;
        }
        if (atom_t2) ++t2_populated;
    }
    EXPECT_GT(populated, 0u)
        << "McConnell Welford T0 all-zero — calculator regression or attach miss";
    EXPECT_GT(t1_populated, 0u)
        << "McConnell Welford T1 all-zero — Lesson 19 regression (was the "
           "2026-05-17 critical bug; McConnell's three-term form has "
           "nonzero antisymmetric part)";
    EXPECT_GT(t2_populated, 0u)
        << "McConnell Welford T2 per-component all-zero — T2 Completeness regression";
    std::cout << "McConnell Welford zero-evaluable-sample atoms=" << unavailable
              << " / " << tp.AtomCount() << "\n";

    // Codex 2026-05-18: dxdt_n must equal delta_n on a well-formed
    // trajectory (no duplicated-timestamp frames). Regression guard
    // for the separate-counter fix that skips zero-dt frames rather
    // than zero-filling the rate accumulator.
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& w = tp.AtomAt(i).mc_welford;
        EXPECT_EQ(w.dxdt_n, w.delta_n)
            << "dxdt_n != delta_n on atom " << i
            << " — investigate timestamp duplication in fixture";
    }

    std::cout << "McConnellWelford integration: populated=" << populated
              << "/" << tp.AtomCount()
              << " max|t0|=" << max_abs_t0
              << " t1_pop=" << t1_populated << " max|t1|=" << max_abs_t1
              << " t2_pop=" << t2_populated << " max|t2|=" << max_abs_t2
              << " (A^-3)\n";
}
