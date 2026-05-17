//
// test_eeq_welford: discipline + integration for EeqWelfordTrajectoryResult.
// AV-pattern scalar Welford on Eeq atomic charge. Discipline trio +
// Integration1P9J — no synthetic test (shared Welford math, see
// test_hm_welford.cpp header).
//

#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "EeqResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "EeqWelfordTrajectoryResult.h"
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

void BuildBaseConfig(nmr::RunConfiguration& config, const std::string& name, std::size_t stride) {
    config.SetName(name);
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::EeqResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::EeqWelfordTrajectoryResult::Create(tp_in);
        });
    config.SetStride(stride);
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
        EXPECT_EQ(ta.eeq_charge_n_frames, 1u);
        EXPECT_DOUBLE_EQ(ta.eeq_charge_std, 0.0);
        EXPECT_EQ(ta.eeq_charge_delta_n, 0u);
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
    const double mean_first = tp.AtomAt(probe).eeq_charge_mean;
    const double std_first  = tp.AtomAt(probe).eeq_charge_std;

    tr.Finalize(tp, traj);

    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).eeq_charge_mean, mean_first);
    EXPECT_DOUBLE_EQ(tp.AtomAt(probe).eeq_charge_std,  std_first);
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
        EXPECT_TRUE(std::isfinite(ta.eeq_charge_mean));
        EXPECT_TRUE(std::isfinite(ta.eeq_charge_std));
        EXPECT_EQ(ta.eeq_charge_n_frames, traj.FrameCount());
        if (std::abs(ta.eeq_charge_mean) > 1e-12) ++populated;
        max_abs = std::max(max_abs, std::abs(ta.eeq_charge_mean));
    }
    EXPECT_GT(populated, 0u);
    std::cout << "EeqWelford integration: populated=" << populated
              << "/" << tp.AtomCount()
              << " max|mean|=" << max_abs << " e\n";
}
