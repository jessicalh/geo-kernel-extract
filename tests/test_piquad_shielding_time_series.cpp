//
// test_piquad_shielding_time_series: discipline + integration for
// PiQuadrupoleShieldingTimeSeriesTrajectoryResult.
//

#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "PiQuadrupoleResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DenseBuffer.h"
#include "PiQuadrupoleShieldingTimeSeriesTrajectoryResult.h"
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
constexpr int kFingerprintResidueNumber = 28;
const std::string kFingerprintChainId = "A";
constexpr double kT0SanityBoundPpm = 50.0;

std::string TrrPathFor(const std::string& tpr_path) {
    return fs::path(tpr_path).replace_extension(".trr").string();
}
std::string ProductionDirFor(const std::string& tpr_path) {
    return fs::path(tpr_path).parent_path().string();
}
bool FixtureAvailable(const nmr::test::AmberTrajectoryFixture& fix) {
    return !fix.tpr_path.empty() && fs::exists(fix.tpr_path)
        && fs::exists(TrrPathFor(fix.tpr_path)) && fs::exists(fix.edr_path);
}
void LoadCalculatorConfig() {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(std::string(NMR_TEST_DATA_DIR) + "/../data/calculator_params.toml");
}
nmr::SphericalTensor SyntheticTensor(size_t i, size_t t) {
    nmr::SphericalTensor s;
    s.T0 = static_cast<double>(i) + static_cast<double>(t) * 100.0;
    for (size_t k = 0; k < 3; ++k)
        s.T1[k] = static_cast<double>(i) + t * 100.0 + (k + 1) * 1e-2;
    for (size_t k = 0; k < 5; ++k)
        s.T2[k] = static_cast<double>(i) + t * 100.0 + (k + 1) * 1e-3;
    return s;
}
bool SphericalEqual(const nmr::SphericalTensor& a, const nmr::SphericalTensor& b, double tol) {
    if (std::abs(a.T0 - b.T0) > tol) return false;
    for (size_t k = 0; k < 3; ++k) if (std::abs(a.T1[k] - b.T1[k]) > tol) return false;
    for (size_t k = 0; k < 5; ++k) if (std::abs(a.T2[k] - b.T2[k]) > tol) return false;
    return true;
}

}  // namespace


TEST(PiQuadrupoleShieldingTimeSeries, SyntheticFourFrames) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const size_t Ntp = tp.AtomCount();
    auto tr = nmr::PiQuadrupoleShieldingTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 4;
    std::vector<nmr::Vec3> positions(Ntp, nmr::Vec3::Zero());
    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        for (size_t i = 0; i < Ntp; ++i)
            conf->MutableAtomAt(i).piquad_shielding_contribution = SyntheticTensor(i, t);
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    auto* buf = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(typeid(
        nmr::PiQuadrupoleShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), Ntp);
    EXPECT_EQ(buf->StridePerAtom(), kFrames);
    for (size_t i : {size_t(0), Ntp / 2, Ntp - 1})
        for (size_t t = 0; t < kFrames; ++t)
            EXPECT_TRUE(SphericalEqual(buf->At(i, t), SyntheticTensor(i, t), 1e-12));

    const std::string h5_path = (fs::temp_directory_path() /
        ("piquad_shielding_ts_unit_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/piquad_shielding_time_series"));
    fs::remove(h5_path);
}


TEST(PiQuadrupoleShieldingTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("PiQuadrupoleShieldingTimeSeriesFrame0Semantics");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::PiQuadrupoleResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::PiQuadrupoleShieldingTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    ASSERT_EQ(traj.FrameCount(), 1u);
}


TEST(PiQuadrupoleShieldingTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("PiQuadrupoleShieldingTimeSeriesFinalizeIdempotency");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::PiQuadrupoleResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::PiQuadrupoleShieldingTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto* buf_first = tp.GetDenseBuffer<nmr::SphericalTensor>(
        std::type_index(typeid(nmr::PiQuadrupoleShieldingTimeSeriesTrajectoryResult)));
    const std::size_t T_first = buf_first->StridePerAtom();
    auto& tr = tp.Result<nmr::PiQuadrupoleShieldingTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);
    auto* buf_second = tp.GetDenseBuffer<nmr::SphericalTensor>(
        std::type_index(typeid(nmr::PiQuadrupoleShieldingTimeSeriesTrajectoryResult)));
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
}


TEST(PiQuadrupoleShieldingTimeSeries, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("PiQuadrupoleShieldingTimeSeriesH5RoundTrip");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::PiQuadrupoleResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::PiQuadrupoleShieldingTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::PiQuadrupoleShieldingTimeSeriesTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("piquad_shielding_ts_h5_roundtrip_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/piquad_shielding_time_series");
    std::string parity, units;
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("units").read(units);
    EXPECT_EQ(parity, "0e+1o+2e");
    EXPECT_EQ(units, "Angstrom^-4");
    fs::remove(h5_path);
}


TEST(PiQuadrupoleShieldingTimeSeries, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("PiQuadrupoleShieldingTimeSeriesIntegration");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::PiQuadrupoleResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::PiQuadrupoleShieldingTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& protein = tp.ProteinRef();
    size_t fingerprint_atom = nmr::Residue::NONE;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const auto& res = protein.ResidueAt(ri);
        if (res.chain_id == kFingerprintChainId && res.sequence_number == kFingerprintResidueNumber) {
            fingerprint_atom = res.CA; break;
        }
    }
    ASSERT_NE(fingerprint_atom, nmr::Residue::NONE);

    auto* buf = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(
        typeid(nmr::PiQuadrupoleShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    const nmr::SphericalTensor& cell = buf->At(fingerprint_atom, 0);
    EXPECT_TRUE(std::isfinite(cell.T0));
    EXPECT_LT(std::abs(cell.T0), kT0SanityBoundPpm);

    // PiQuad kernel is pure-T2 by Laplace's equation (analytically
    // traceless EFG; see PiQuadrupoleResult.cpp:40-44 and the asserted
    // tests/test_pi_quadrupole_result.cpp:475-476 EXPECT_LT(max_t0, 1e-10)).
    // T0 is structurally zero; the signal lives in T2. Check max |T2[k]|.
    size_t populated_t2 = 0;
    double max_abs_t2 = 0.0;
    for (size_t i = 0; i < buf->AtomCount(); ++i) {
        double atom_max_t2 = 0.0;
        for (size_t t = 0; t < buf->StridePerAtom(); ++t) {
            const nmr::SphericalTensor& st = buf->At(i, t);
            EXPECT_TRUE(std::isfinite(st.T0)
                     && std::isfinite(st.T1[0]) && std::isfinite(st.T1[1]) && std::isfinite(st.T1[2])
                     && std::isfinite(st.T2[0]) && std::isfinite(st.T2[1])
                     && std::isfinite(st.T2[2]) && std::isfinite(st.T2[3]) && std::isfinite(st.T2[4]));
            for (int k = 0; k < 5; ++k) {
                const double a = std::abs(st.T2[k]);
                max_abs_t2 = std::max(max_abs_t2, a);
                atom_max_t2 = std::max(atom_max_t2, a);
            }
        }
        if (atom_max_t2 > 1e-10) ++populated_t2;
    }
    // Assert kernel fires somewhere above noise on T2 channels.
    EXPECT_GT(max_abs_t2, 1e-6) << "PiQuad T2 channels at noise floor — kernel may be broken";
    std::cout << "PiQuadrupoleShieldingTimeSeries CA-anchor T0=" << cell.T0
              << " A^-4 (structural zero), max|T2|=" << max_abs_t2 << " A^-4, populated_t2="
              << populated_t2 << "/" << buf->AtomCount()
              << " (pure-T2 traceless EFG; 1/r^4 falloff)\n";
}
