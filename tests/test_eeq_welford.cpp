//
// test_eeq_welford: discipline + integration for EeqWelfordTrajectoryResult.
// AV-pattern scalar Welford on Eeq atomic charge. Discipline trio +
// Integration1P9J — no synthetic test (shared Welford math, see
// test_hm_welford.cpp header).
//

#include "EeqResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "EeqWelfordTrajectoryResult.h"
#include "TestConfig.h"
#include "OperationLog.h"
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

std::string TrrPathFor(const std::string& p) { return fs::path(p).replace_extension(".trr").string(); }
std::string ProductionDirFor(const std::string& p) { return fs::path(p).parent_path().string(); }
bool FixtureAvailable(const nmr::test::AmberTrajectoryFixture& fix) {
    return !fix.tpr_path.empty() && fs::exists(fix.tpr_path)
        && fs::exists(TrrPathFor(fix.tpr_path)) && fs::exists(fix.edr_path);
}
void LoadCalculatorConfig() {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(std::string(NMR_TEST_DATA_DIR) + "/../data/calculator_params.toml");
}

// Per-file bridge to the shared TestConfig profile system. See
// spec/plan/test-suite-realignment-deferred-2026-05-17.md.
void BuildBaseConfig(nmr::RunConfiguration& config, const std::string& name, std::size_t stride) {
    config = nmr::test::BuildTestConfig(nmr::test::TestProfile::KernelOnly, name, stride);
    config.RequireConformationResult(typeid(nmr::EeqResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::EeqWelfordTrajectoryResult::Create(tp_in);
        });
}

}  // namespace


TEST(EeqWelford, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    BuildBaseConfig(config, "EeqWelfordFrame0Semantics", 99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    ASSERT_EQ(traj.FrameCount(), 1u);

    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& ta = tp.AtomAt(i);
        EXPECT_EQ(ta.eeq_welford.n_frames, 1u);
        EXPECT_DOUBLE_EQ(ta.eeq_welford.charge.std, 0.0);
        EXPECT_EQ(ta.eeq_welford.delta_n, 0u);
    }
}


TEST(EeqWelford, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    BuildBaseConfig(config, "EeqWelfordFinalizeIdempotency", 99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::EeqWelfordTrajectoryResult>();
    const size_t probe = tp.AtomCount() / 2;
    const double mean_first = tp.AtomAt(probe).eeq_welford.charge.mean;
    const double std_first  = tp.AtomAt(probe).eeq_welford.charge.std;

    tr.Finalize(tp, traj);

    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).eeq_welford.charge.mean, mean_first);
    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).eeq_welford.charge.std,  std_first);
}


TEST(EeqWelford, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    BuildBaseConfig(config, "EeqWelfordH5RoundTrip", 99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::EeqWelfordTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("eeq_welford_h5_roundtrip_" + std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr.WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/eeq_welford"));
    auto grp = reopen.getGroup("/trajectory/eeq_welford");

    std::string units;
    grp.getAttribute("units").read(units);
    EXPECT_EQ(units, "elementary_charge");

    const auto dims = grp.getDataSet("mean").getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 1u);
    EXPECT_EQ(dims[0], tp.AtomCount());

    fs::remove(h5_path);
}


TEST(EeqWelford, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    BuildBaseConfig(config, "EeqWelfordIntegration", 300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_GE(traj.FrameCount(), 2u);

    size_t populated = 0;
    double max_abs = 0.0;
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& ta = tp.AtomAt(i);
        EXPECT_TRUE(std::isfinite(ta.eeq_welford.charge.mean));
        EXPECT_TRUE(std::isfinite(ta.eeq_welford.charge.std));
        EXPECT_EQ(ta.eeq_welford.n_frames, traj.FrameCount());
        if (std::abs(ta.eeq_welford.charge.mean) > 1e-12) ++populated;
        max_abs = std::max(max_abs, std::abs(ta.eeq_welford.charge.mean));
    }
    EXPECT_GT(populated, 0u);

    // Codex 2026-05-18: dxdt_n must equal delta_n on a well-formed
    // trajectory (no duplicated-timestamp frames). Regression guard
    // for the separate-counter fix that skips zero-dt frames rather
    // than zero-filling the rate accumulator.
    for (size_t i = 0; i < tp.AtomCount(); ++i) {
        const auto& w = tp.AtomAt(i).eeq_welford;
        EXPECT_EQ(w.dxdt_n, w.delta_n)
            << "dxdt_n != delta_n on atom " << i
            << " — investigate timestamp duplication in fixture";
    }

    std::cout << "EeqWelford integration: populated=" << populated
              << "/" << tp.AtomCount()
              << " max|mean|=" << max_abs << " e\n";
}
