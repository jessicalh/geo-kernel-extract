//
// test_j_coupling_time_series: discipline + integration for
// JCouplingTimeSeriesTrajectoryResult. Per-residue per-frame Karplus
// 3J observables; thin transform on positions + Residue cache; no
// source ConformationResult dependency.
//

#include "AminoAcidType.h"
#include "CalculatorConfig.h"
#include "GeometryResult.h"
#include "JCouplingTimeSeriesTrajectoryResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "SpatialIndexResult.h"
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
    nmr::CalculatorConfig::Load(std::string(NMR_TEST_DATA_DIR) +
                                "/../data/calculator_params.toml");
}

nmr::RunConfiguration BuildConfig(unsigned stride) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::JCouplingTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


TEST(JCouplingTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::JCouplingTimeSeriesTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);
}


TEST(JCouplingTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::JCouplingTimeSeriesTrajectoryResult>();
    const std::size_t T = tr.NumFrames();
    tr.Finalize(tp, traj);
    EXPECT_EQ(tr.NumFrames(), T);
}


TEST(JCouplingTimeSeries, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::JCouplingTimeSeriesTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();

    const std::string h5_path = (fs::temp_directory_path() /
        ("j_coupling_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/j_coupling_time_series"));
    auto grp = reopen.getGroup("/trajectory/j_coupling_time_series");

    for (const std::string& name : {"J_HN_Halpha", "J_N_Cgamma",
                                     "J_Cprime_Cgamma"}) {
        ASSERT_TRUE(grp.exist(name)) << name;
        const auto dims = grp.getDataSet(name).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], tr.NumFrames());
    }

    EXPECT_TRUE(grp.exist("J_HN_Halpha_exists"));
    EXPECT_TRUE(grp.exist("J_chi1_exists"));
    EXPECT_TRUE(grp.exist("residue_index_per_atom"));
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));

    // Convention attrs
    std::string karplus, hn_ha, n_cg, units;
    grp.getAttribute("karplus_form").read(karplus);
    grp.getAttribute("J_HN_Halpha_coefficients").read(hn_ha);
    grp.getAttribute("J_N_Cgamma_coefficients").read(n_cg);
    grp.getAttribute("units").read(units);
    EXPECT_NE(karplus.find("cos^2"), std::string::npos);
    EXPECT_NE(hn_ha.find("Wang"), std::string::npos);
    EXPECT_NE(n_cg.find("Perez"), std::string::npos);
    EXPECT_EQ(units, "Hz");

    fs::remove(h5_path);
}


TEST(JCouplingTimeSeries, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(300);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::JCouplingTimeSeriesTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    const std::size_t T = tr.NumFrames();
    EXPECT_GE(T, 2u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("j_coupling_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/j_coupling_time_series");

    std::vector<std::vector<double>> j_hn, j_n_cg, j_cp_cg;
    grp.getDataSet("J_HN_Halpha").read(j_hn);
    grp.getDataSet("J_N_Cgamma").read(j_n_cg);
    grp.getDataSet("J_Cprime_Cgamma").read(j_cp_cg);
    ASSERT_EQ(j_hn.size(), R);

    std::vector<std::uint8_t> hn_exists, chi1_exists;
    grp.getDataSet("J_HN_Halpha_exists").read(hn_exists);
    grp.getDataSet("J_chi1_exists").read(chi1_exists);
    ASSERT_EQ(hn_exists.size(), R);
    ASSERT_EQ(chi1_exists.size(), R);

    // Per-residue invariants:
    // (1) PRO has no HN, so J_HN_Halpha exists mask == 0 and the
    //     corresponding row is all-NaN.
    // (2) GLY/ALA have no chi1, so J_chi1_exists == 0 and both
    //     chi1-dependent J channels are all-NaN.
    // (3) Karplus arithmetic bounds. 3J(θ) = A·cos²(θ) + B·cos(θ) + C
    //     is a closed-form quadratic in cos(θ); over cos(θ) ∈ [-1, 1]
    //     and A > 0, B < 0 (for all three published parametrizations)
    //     the global max is f(-1) = A - B + C = A + |B| + C and the
    //     global min is at the vertex u* = -B/(2A) ∈ (0, 1), evaluating
    //     to f(u*) = C - B²/(4A). Numerical values:
    //       J(HN,Hα)  : [1.48, 9.87] Hz (Wang & Bax 1996)
    //       J(N,Cγ)   : [0.32, 2.15] Hz (Pérez 2001)
    //       J(C',Cγ)  : [0.20, 2.56] Hz (Pérez 2001)
    //     A violation indicates wrong coefficients or wrong dihedral
    //     definition — these aren't "physical range overages" (no
    //     TRR-float32 / chi-fallback path enters here). Bounds below
    //     have small slack for float epsilon and to leave room for
    //     refit Karplus parametrizations.
    std::size_t pro_count = 0, gly_ala_count = 0;
    std::size_t hn_finite_obs = 0, n_cg_finite_obs = 0, cp_cg_finite_obs = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        const auto& res = tp.ProteinRef().ResidueAt(ri);
        const bool is_pro = (res.type == nmr::AminoAcid::PRO);
        const bool is_gly_ala = (res.type == nmr::AminoAcid::GLY ||
                                   res.type == nmr::AminoAcid::ALA);
        if (is_pro) {
            ++pro_count;
            EXPECT_EQ(hn_exists[ri], 0u);
            for (std::size_t t = 0; t < T; ++t) {
                EXPECT_TRUE(std::isnan(j_hn[ri][t]))
                    << "PRO ri=" << ri << " J_HN_Halpha should be NaN";
            }
        }
        if (is_gly_ala) {
            ++gly_ala_count;
            EXPECT_EQ(chi1_exists[ri], 0u);
            for (std::size_t t = 0; t < T; ++t) {
                EXPECT_TRUE(std::isnan(j_n_cg[ri][t]));
                EXPECT_TRUE(std::isnan(j_cp_cg[ri][t]));
            }
        }
        for (std::size_t t = 0; t < T; ++t) {
            if (std::isfinite(j_hn[ri][t])) {
                ++hn_finite_obs;
                EXPECT_GE(j_hn[ri][t], 1.0)  << "J(HN,Hα) below Karplus min ~1.48";
                EXPECT_LE(j_hn[ri][t], 10.0) << "J(HN,Hα) above Karplus max ~9.87";
            }
            if (std::isfinite(j_n_cg[ri][t])) {
                ++n_cg_finite_obs;
                EXPECT_GE(j_n_cg[ri][t], 0.0) << "J(N,Cγ) below Karplus min ~0.32";
                EXPECT_LE(j_n_cg[ri][t], 2.5) << "J(N,Cγ) above Karplus max ~2.15";
            }
            if (std::isfinite(j_cp_cg[ri][t])) {
                ++cp_cg_finite_obs;
                EXPECT_GE(j_cp_cg[ri][t], 0.0) << "J(C',Cγ) below Karplus min ~0.20";
                EXPECT_LE(j_cp_cg[ri][t], 3.0) << "J(C',Cγ) above Karplus max ~2.56";
            }
        }
    }
    std::cout << "JCouplingTimeSeries: R=" << R << " T=" << T
              << " PRO=" << pro_count
              << " GLY+ALA=" << gly_ala_count
              << " | J(HN,Hα) finite obs=" << hn_finite_obs
              << " | J(N,Cγ) finite obs=" << n_cg_finite_obs
              << " | J(C',Cγ) finite obs=" << cp_cg_finite_obs << "\n";
    EXPECT_GT(hn_finite_obs, 0u);
    EXPECT_GT(n_cg_finite_obs, 0u);
    EXPECT_GT(cp_cg_finite_obs, 0u);

    fs::remove(h5_path);
}
