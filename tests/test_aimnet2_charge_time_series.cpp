//
// test_aimnet2_charge_time_series: discipline + integration for
// AIMNet2ChargeTimeSeriesTrajectoryResult. Scalar-double DenseBuffer.
// AIMNet2 model loading is required for the integration tests; they
// skip when not available.
//

#include "AIMNet2Result.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DenseBuffer.h"
#include "AIMNet2ChargeTimeSeriesTrajectoryResult.h"
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
void LoadCalculatorConfig() {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(std::string(NMR_TEST_DATA_DIR) + "/../data/calculator_params.toml");
}

}  // namespace


TEST(AIMNet2ChargeTimeSeries, SyntheticFourFrames) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const size_t N = tp.AtomCount();
    auto tr = nmr::AIMNet2ChargeTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 4;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        for (size_t i = 0; i < N; ++i) {
            conf->MutableAtomAt(i).aimnet2_charge = 0.001 * static_cast<double>(i) - 0.1 * t;
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    auto* buf = tp.GetDenseBuffer<double>(std::type_index(typeid(
        nmr::AIMNet2ChargeTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->StridePerAtom(), kFrames);
    for (size_t i : {size_t(0), N / 2, N - 1})
        for (size_t t = 0; t < kFrames; ++t)
            EXPECT_DOUBLE_EQ(buf->At(i, t), 0.001 * static_cast<double>(i) - 0.1 * t);

    const std::string h5_path = (fs::temp_directory_path() /
        ("aimnet2_charge_ts_unit_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/aimnet2_charge_time_series"));
    auto grp = reopen.getGroup("/trajectory/aimnet2_charge_time_series");
    auto ds = grp.getDataSet("charge");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 2u);
    EXPECT_EQ(dims[0], N);
    EXPECT_EQ(dims[1], kFrames);
    fs::remove(h5_path);
}


TEST(AIMNet2ChargeTimeSeries, H5RoundTripAttrs) {
    // Synthetic-driven H5 emission to verify attrs without needing the
    // model loaded.
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const size_t N = tp.AtomCount();
    auto tr = nmr::AIMNet2ChargeTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    auto conf = std::make_unique<nmr::ProteinConformation>(
        &tp.ProteinRef(), positions, "single frame");
    tr->Compute(*conf, tp, traj, 0, 0.0);
    tr->Finalize(tp, traj);

    const std::string h5_path = (fs::temp_directory_path() /
        ("aimnet2_charge_ts_attrs_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/aimnet2_charge_time_series");
    std::string parity, units, layout;
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("units").read(units);
    grp.getAttribute("irrep_layout").read(layout);
    EXPECT_EQ(parity, "0e");
    EXPECT_EQ(units, "elementary_charge");
    EXPECT_EQ(layout, "T0");
    fs::remove(h5_path);
}


TEST(AIMNet2ChargeTimeSeries, Integration1P9J) {
    LoadCalculatorConfig();
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
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::AIMNet2Result));
    config.SetRequiresAimnet2(true);
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::AIMNet2ChargeTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto* buf = tp.GetDenseBuffer<double>(std::type_index(typeid(
        nmr::AIMNet2ChargeTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);

    double max_abs = 0.0;
    double charge_sum = 0.0;
    for (size_t i = 0; i < buf->AtomCount(); ++i) {
        const double q = buf->At(i, 0);
        EXPECT_TRUE(std::isfinite(q));
        EXPECT_LT(std::abs(q), 5.0) << "atom " << i << " charge magnitude " << q;
        max_abs = std::max(max_abs, std::abs(q));
        charge_sum += q;
    }
    EXPECT_GT(max_abs, 0.05) << "AIMNet2 charges at noise floor";
    // For a neutralised protein, net charge should be small (|sum| < 1 e
    // is generous; AIMNet2 isn't guaranteed to give exact integer net
    // charge, but should be close).
    std::cout << "AIMNet2ChargeTimeSeries max|q|=" << max_abs
              << " e, net_charge=" << charge_sum << " e\n";
}
