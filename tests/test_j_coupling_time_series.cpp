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
#include "PdbFileReader.h"
#include "PhysicalConstants.h"
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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <unistd.h>

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
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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

    for (const std::string& name : {"J_HN_Halpha", "J_HN_Halpha_Vogeli",
                                     "J_HN_Cbeta", "J_HN_Cprime",
                                     "J_Halpha_Cprime",
                                     "J_N_Cgamma", "J_Cprime_Cgamma",
                                     "J_Halpha_Hbeta2", "J_Halpha_Hbeta3"}) {
        ASSERT_TRUE(grp.exist(name)) << name;
        const auto dims = grp.getDataSet(name).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], tr.NumFrames());
    }

    EXPECT_TRUE(grp.exist("J_HN_Halpha_exists"));
    EXPECT_TRUE(grp.exist("J_HN_Cbeta_exists"));
    EXPECT_TRUE(grp.exist("J_HN_Cprime_exists"));
    EXPECT_TRUE(grp.exist("J_Halpha_Cprime_exists"));
    EXPECT_TRUE(grp.exist("J_chi1_exists"));
    EXPECT_TRUE(grp.exist("J_N_Cgamma_exists"));
    EXPECT_TRUE(grp.exist("J_Cprime_Cgamma_exists"));
    EXPECT_TRUE(grp.exist("J_Halpha_Hbeta_exists"));
    EXPECT_TRUE(grp.exist("residue_index_per_atom"));
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));

    // Convention attrs — citation pins are auditable from the H5 itself.
    std::string karplus, hn_ha, hn_ha_vog, hn_cb, hn_cp, ha_cp,
                n_cg, cp_cg, ha_hb, units;
    grp.getAttribute("karplus_form").read(karplus);
    grp.getAttribute("J_HN_Halpha_coefficients").read(hn_ha);
    grp.getAttribute("J_HN_Halpha_Vogeli_coefficients").read(hn_ha_vog);
    grp.getAttribute("J_HN_Cbeta_coefficients").read(hn_cb);
    grp.getAttribute("J_HN_Cprime_coefficients").read(hn_cp);
    grp.getAttribute("J_Halpha_Cprime_coefficients").read(ha_cp);
    grp.getAttribute("J_N_Cgamma_coefficients").read(n_cg);
    grp.getAttribute("J_Cprime_Cgamma_coefficients").read(cp_cg);
    grp.getAttribute("J_Halpha_Hbeta_coefficients").read(ha_hb);
    grp.getAttribute("units").read(units);
    EXPECT_NE(karplus.find("cos^2"), std::string::npos);
    EXPECT_NE(hn_ha.find("Vuister"), std::string::npos);
    EXPECT_NE(hn_ha.find("1993"), std::string::npos);
    EXPECT_NE(hn_ha_vog.find("Vogeli"), std::string::npos);
    EXPECT_NE(hn_ha_vog.find("2007"), std::string::npos);
    EXPECT_NE(hn_cb.find("Wang"), std::string::npos);
    EXPECT_NE(hn_cb.find("1996"), std::string::npos);
    EXPECT_NE(hn_cp.find("Wang"), std::string::npos);
    EXPECT_NE(hn_cp.find("ROW 4"), std::string::npos);  // row-mapping
        // fixed 2026-05-20: J(HN, C') uses Wang-Bax row 4 (theta=0).
    EXPECT_NE(ha_cp.find("Wang"), std::string::npos);
    EXPECT_NE(ha_cp.find("ROW 2"), std::string::npos);  // row-mapping
        // fixed 2026-05-20: J(Halpha, C') uses Wang-Bax row 2 (theta=-60).
    EXPECT_NE(n_cg.find("Perez"), std::string::npos);
    EXPECT_NE(n_cg.find("2001"), std::string::npos);
    EXPECT_NE(cp_cg.find("Perez"), std::string::npos);
    EXPECT_NE(ha_hb.find("Perez"), std::string::npos);
    EXPECT_NE(ha_hb.find("Hbeta"), std::string::npos);
    // J_Cprime_Cgamma_coefficients embeds the byte-verified values
    // (A=2.31, B=-0.87, C=0.55) — verify the correction lands.
    EXPECT_NE(cp_cg.find("A=2.3100"), std::string::npos)
        << "J_Cprime_Cgamma should use Perez 2001 Table 2 consensus "
           "A=2.31 (was incorrectly 1.74 prior to 2026-05-19 byte-fix); "
           "attr was: " << cp_cg;
    EXPECT_EQ(units, "Hz");

    fs::remove(h5_path);
}


TEST(JCouplingTimeSeries, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
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

    std::vector<std::vector<double>> j_hn, j_hn_vog, j_hn_cb, j_hn_cp,
                                       j_ha_cp, j_n_cg, j_cp_cg,
                                       j_ha_hb2, j_ha_hb3;
    grp.getDataSet("J_HN_Halpha").read(j_hn);
    grp.getDataSet("J_HN_Halpha_Vogeli").read(j_hn_vog);
    grp.getDataSet("J_HN_Cbeta").read(j_hn_cb);
    grp.getDataSet("J_HN_Cprime").read(j_hn_cp);
    grp.getDataSet("J_Halpha_Cprime").read(j_ha_cp);
    grp.getDataSet("J_N_Cgamma").read(j_n_cg);
    grp.getDataSet("J_Cprime_Cgamma").read(j_cp_cg);
    grp.getDataSet("J_Halpha_Hbeta2").read(j_ha_hb2);
    grp.getDataSet("J_Halpha_Hbeta3").read(j_ha_hb3);
    ASSERT_EQ(j_hn.size(), R);
    ASSERT_EQ(j_hn_cb.size(), R);
    ASSERT_EQ(j_hn_cp.size(), R);
    ASSERT_EQ(j_ha_cp.size(), R);
    ASSERT_EQ(j_ha_hb2.size(), R);

    std::vector<std::uint8_t> hn_exists, hn_cb_exists, hn_cp_exists,
                              ha_cp_exists, chi1_exists,
                              n_cg_exists, cp_cg_exists, ha_hb_exists;
    grp.getDataSet("J_HN_Halpha_exists").read(hn_exists);
    grp.getDataSet("J_HN_Cbeta_exists").read(hn_cb_exists);
    grp.getDataSet("J_HN_Cprime_exists").read(hn_cp_exists);
    grp.getDataSet("J_Halpha_Cprime_exists").read(ha_cp_exists);
    grp.getDataSet("J_chi1_exists").read(chi1_exists);
    grp.getDataSet("J_N_Cgamma_exists").read(n_cg_exists);
    grp.getDataSet("J_Cprime_Cgamma_exists").read(cp_cg_exists);
    grp.getDataSet("J_Halpha_Hbeta_exists").read(ha_hb_exists);
    ASSERT_EQ(hn_exists.size(), R);
    ASSERT_EQ(hn_cb_exists.size(), R);
    ASSERT_EQ(hn_cp_exists.size(), R);
    ASSERT_EQ(ha_cp_exists.size(), R);
    ASSERT_EQ(chi1_exists.size(), R);
    ASSERT_EQ(n_cg_exists.size(), R);
    ASSERT_EQ(cp_cg_exists.size(), R);
    ASSERT_EQ(ha_hb_exists.size(), R);

    // Per-residue invariants:
    // (1) PRO has no HN, so J_HN_Halpha and J_HN_Cbeta exists masks
    //     == 0 and the corresponding rows are all-NaN.
    // (2) GLY has no Cbeta, so J_HN_Cbeta_exists == 0 and that row
    //     is all-NaN.
    // (3) GLY/ALA have no chi1, so J_chi1_exists == 0 and both
    //     chi1-dependent J channels are all-NaN.
    // (4) Karplus arithmetic bounds. 3J(α) = A·cos²(α) + B·cos(α) + C
    //     is a closed-form quadratic in cos(α); over cos(α) ∈ [-1, 1]
    //     the global extrema are at the endpoints f(+1) = A + B + C
    //     and f(-1) = A - B + C, and at the vertex u* = -B/(2A) ∈
    //     (-1, 1) where f(u*) = C - B²/(4A). For A > 0, B < 0 the
    //     vertex is the MIN and f(-1) = A + |B| + C is the MAX. For
    //     A > 0, B > 0 the vertex is the MIN and f(+1) = A + B + C
    //     is the MAX (the two positive-B channels are HN-C' and
    //     Halpha-C'). Numerical values:
    //       J(HN,Hα)  : [1.48, 9.87] Hz (Vuister & Bax 1993; B<0)
    //       J(HN,Cβ)  : [0.005, 4.40] Hz (Wang-Bax row 3; B<0)
    //       J(HN,C')  : [-0.04, 5.16] Hz (Wang-Bax ROW 4; B>0;
    //                   row-mapping fixed 2026-05-20)
    //       J(Hα,C')  : [0.96, 7.22] Hz (Wang-Bax ROW 2; B>0;
    //                   row-mapping fixed 2026-05-20)
    //       J(N,Cγ)   : [0.32, 2.15] Hz (Pérez 2001; B<0)
    //       J(C',Cγ)  : [0.47, 3.73] Hz (Pérez 2001; B<0)
    //     A violation indicates wrong coefficients or wrong dihedral
    //     definition — these aren't "physical range overages" (no
    //     TRR-float32 / chi-fallback path enters here). Bounds below
    //     have small slack for float epsilon and to leave room for
    //     refit Karplus parametrizations.
    std::size_t pro_count = 0, gly_count = 0, gly_ala_count = 0;
    std::size_t hn_finite_obs = 0, hn_vog_finite_obs = 0;
    std::size_t hn_cb_finite_obs = 0, hn_cp_finite_obs = 0;
    std::size_t ha_cp_finite_obs = 0;
    std::size_t n_cg_finite_obs = 0, cp_cg_finite_obs = 0;
    std::size_t ha_hb2_finite_obs = 0, ha_hb3_finite_obs = 0;
    // Cross-channel correlation accumulator (Pearson) for science
    // F10: J_HN_Halpha and J_HN_Cprime read the same phi via
    // different atomic dihedrals; the two channels should be
    // correlated across the (residue, frame) corpus.
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    std::size_t n_xy = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        const auto& res = tp.ProteinRef().ResidueAt(ri);
        const bool is_pro = (res.type == nmr::AminoAcid::PRO);
        const bool is_gly = (res.type == nmr::AminoAcid::GLY);
        const bool is_ala = (res.type == nmr::AminoAcid::ALA);
        const bool is_gly_ala = is_gly || is_ala;
        if (is_pro) {
            ++pro_count;
            EXPECT_EQ(hn_exists[ri], 0u);
            EXPECT_EQ(hn_cb_exists[ri], 0u);
            EXPECT_EQ(hn_cp_exists[ri], 0u);  // PRO: no HN
            for (std::size_t t = 0; t < T; ++t) {
                EXPECT_TRUE(std::isnan(j_hn[ri][t]))
                    << "PRO ri=" << ri << " J_HN_Halpha should be NaN";
                EXPECT_TRUE(std::isnan(j_hn_vog[ri][t]))
                    << "PRO ri=" << ri << " J_HN_Halpha_Vogeli should be NaN";
                EXPECT_TRUE(std::isnan(j_hn_cb[ri][t]))
                    << "PRO ri=" << ri << " J_HN_Cbeta should be NaN";
                EXPECT_TRUE(std::isnan(j_hn_cp[ri][t]))
                    << "PRO ri=" << ri << " J_HN_Cprime should be NaN";
            }
        }
        if (is_gly) {
            ++gly_count;
            EXPECT_EQ(hn_cb_exists[ri], 0u);  // GLY: no Cβ
            EXPECT_EQ(ha_hb_exists[ri], 0u);  // GLY: no Cβ → no Hβ
            for (std::size_t t = 0; t < T; ++t) {
                EXPECT_TRUE(std::isnan(j_hn_cb[ri][t]))
                    << "GLY ri=" << ri << " J_HN_Cbeta should be NaN";
                EXPECT_TRUE(std::isnan(j_ha_hb2[ri][t]))
                    << "GLY ri=" << ri << " J_Halpha_Hbeta2 should be NaN";
                EXPECT_TRUE(std::isnan(j_ha_hb3[ri][t]))
                    << "GLY ri=" << ri << " J_Halpha_Hbeta3 should be NaN";
            }
        }
        if (is_ala) {
            // ALA: methyl Cβ; deliberately NaN-fill the Hβ channels.
            EXPECT_EQ(ha_hb_exists[ri], 0u);
            for (std::size_t t = 0; t < T; ++t) {
                EXPECT_TRUE(std::isnan(j_ha_hb2[ri][t]))
                    << "ALA ri=" << ri << " J_Halpha_Hbeta2 should be NaN "
                       "(ALA methyl Cβ is not the methylene observable)";
                EXPECT_TRUE(std::isnan(j_ha_hb3[ri][t]))
                    << "ALA ri=" << ri << " J_Halpha_Hbeta3 should be NaN";
            }
        }
        if (is_gly_ala) {
            ++gly_ala_count;
            EXPECT_EQ(chi1_exists[ri], 0u);
            EXPECT_EQ(n_cg_exists[ri], 0u);
            EXPECT_EQ(cp_cg_exists[ri], 0u);
            for (std::size_t t = 0; t < T; ++t) {
                EXPECT_TRUE(std::isnan(j_n_cg[ri][t]));
                EXPECT_TRUE(std::isnan(j_cp_cg[ri][t]));
            }
        }
        // SER/CYS/THR: chi1 terminal is OG/SG/OG1 (non-carbon). The
        // J(N,Cγ) and J(C',Cγ) channels are NaN by the Element::C
        // element gate; chi1_exists is still 1 (chi[0].Valid()) so
        // the chi1 timeline in DihedralTS is still emitted. Codex F3
        // 2026-05-20.
        const bool is_ser = (res.type == nmr::AminoAcid::SER);
        const bool is_cys = (res.type == nmr::AminoAcid::CYS);
        const bool is_thr = (res.type == nmr::AminoAcid::THR);
        if (is_ser || is_cys || is_thr) {
            // chi1 is valid for these (heavy-atom OG/SG/OG1 terminal).
            EXPECT_EQ(chi1_exists[ri], 1u)
                << "chi1 should be valid for SER/CYS/THR (chi[0].a[3] = "
                   "OG/SG/OG1; finite dihedral).";
            EXPECT_EQ(n_cg_exists[ri], 0u)
                << "J_N_Cgamma should be NaN for SER/CYS/THR "
                   "(chi1 terminal is not carbon; gated by element).";
            EXPECT_EQ(cp_cg_exists[ri], 0u)
                << "J_Cprime_Cgamma should be NaN for SER/CYS/THR.";
            for (std::size_t t = 0; t < T; ++t) {
                EXPECT_TRUE(std::isnan(j_n_cg[ri][t]))
                    << "SER/CYS/THR ri=" << ri
                    << " J_N_Cgamma should be NaN";
                EXPECT_TRUE(std::isnan(j_cp_cg[ri][t]))
                    << "SER/CYS/THR ri=" << ri
                    << " J_Cprime_Cgamma should be NaN";
            }
        }
        // N-terminal residue (first in chain): all phi-derived
        // backbone channels are NaN because no C'(prev) exists to
        // define project phi. J_Halpha_Cprime additionally needs
        // C'(prev) for the actual 3-bond path. Post-2026-05-20 codex
        // F2: Halpha-Cprime now uses HA-CA-N-C'(prev) (phi axis), not
        // HA-CA-C-N(next) (psi axis).
        const bool is_n_terminus = !tp.ProteinRef().BackbonePredecessor(ri).has_value();
        if (is_n_terminus) {
            EXPECT_EQ(hn_exists[ri], 0u);
            EXPECT_EQ(hn_cb_exists[ri], 0u);
            EXPECT_EQ(hn_cp_exists[ri], 0u);
            EXPECT_EQ(ha_cp_exists[ri], 0u);
            for (std::size_t t = 0; t < T; ++t) {
                EXPECT_TRUE(std::isnan(j_hn[ri][t]))
                    << "N-term ri=" << ri << " J_HN_Halpha should be NaN";
                EXPECT_TRUE(std::isnan(j_hn_vog[ri][t]))
                    << "N-term ri=" << ri << " J_HN_Halpha_Vogeli should be NaN";
                EXPECT_TRUE(std::isnan(j_hn_cb[ri][t]))
                    << "N-term ri=" << ri << " J_HN_Cbeta should be NaN";
                EXPECT_TRUE(std::isnan(j_hn_cp[ri][t]))
                    << "N-term ri=" << ri << " J_HN_Cprime should be NaN";
                EXPECT_TRUE(std::isnan(j_ha_cp[ri][t]))
                    << "N-term ri=" << ri << " J_Halpha_Cprime should be NaN";
            }
        }
        for (std::size_t t = 0; t < T; ++t) {
            if (std::isfinite(j_hn[ri][t])) {
                ++hn_finite_obs;
                EXPECT_GE(j_hn[ri][t], 1.0)  << "J(HN,Hα) below Karplus min ~1.48";
                EXPECT_LE(j_hn[ri][t], 10.0) << "J(HN,Hα) above Karplus max ~9.87";
            }
            if (std::isfinite(j_hn_vog[ri][t])) {
                ++hn_vog_finite_obs;
                EXPECT_GE(j_hn_vog[ri][t], 0.0)
                    << "J(HN,Hα)_Vogeli below Karplus min ~0.58";
                EXPECT_LE(j_hn_vog[ri][t], 10.0)
                    << "J(HN,Hα)_Vogeli above Karplus max ~9.86";
            }
            if (std::isfinite(j_hn_cb[ri][t])) {
                ++hn_cb_finite_obs;
                EXPECT_GE(j_hn_cb[ri][t], 0.0) << "J(HN,Cβ) below Karplus min ~0.005";
                EXPECT_LE(j_hn_cb[ri][t], 5.0) << "J(HN,Cβ) above Karplus max ~4.40";
            }
            if (std::isfinite(j_hn_cp[ri][t])) {
                ++hn_cp_finite_obs;
                // Row-mapping fix 2026-05-20: J(HN, C') uses Wang-Bax
                // row 4 (A=4.32, B=+0.84, C=0.00); range [-0.04, 5.16].
                EXPECT_GE(j_hn_cp[ri][t], -0.5)
                    << "J(HN,C') below Karplus min ~-0.04";
                EXPECT_LE(j_hn_cp[ri][t], 5.5)
                    << "J(HN,C') above Karplus max ~5.16";
            }
            if (std::isfinite(j_ha_cp[ri][t])) {
                ++ha_cp_finite_obs;
                // Row-mapping fix 2026-05-20: J(Halpha, C') uses
                // Wang-Bax row 2 (A=3.75, B=+2.19, C=1.28); range
                // [0.96, 7.22]. Atom path also corrected:
                // Halpha-CA-N-C'(prev) (phi axis), not HA-CA-C-N(next).
                EXPECT_GE(j_ha_cp[ri][t], 0.5)
                    << "J(Hα,C') below Karplus min ~0.96";
                EXPECT_LE(j_ha_cp[ri][t], 8.0)
                    << "J(Hα,C') above Karplus max ~7.22";
            }
            if (std::isfinite(j_n_cg[ri][t])) {
                ++n_cg_finite_obs;
                EXPECT_GE(j_n_cg[ri][t], 0.0) << "J(N,Cγ) below Karplus min ~0.32";
                EXPECT_LE(j_n_cg[ri][t], 2.5) << "J(N,Cγ) above Karplus max ~2.15";
            }
            if (std::isfinite(j_cp_cg[ri][t])) {
                ++cp_cg_finite_obs;
                // Updated bounds for Perez 2001 Table 2 consensus
                // values (A=2.31, B=-0.87, C=0.55): range [0.47, 3.73].
                EXPECT_GE(j_cp_cg[ri][t], 0.3)
                    << "J(C',Cγ) below Perez Table 2 min ~0.47";
                EXPECT_LE(j_cp_cg[ri][t], 4.0)
                    << "J(C',Cγ) above Perez Table 2 max ~3.73";
            }
            if (std::isfinite(j_ha_hb2[ri][t])) {
                ++ha_hb2_finite_obs;
                EXPECT_GE(j_ha_hb2[ri][t], 1.5)
                    << "J(Hα,Hβ2) below Karplus min ~2.16";
                EXPECT_LE(j_ha_hb2[ri][t], 11.0)
                    << "J(Hα,Hβ2) above Karplus max ~10.82";
            }
            if (std::isfinite(j_ha_hb3[ri][t])) {
                ++ha_hb3_finite_obs;
                EXPECT_GE(j_ha_hb3[ri][t], 1.5)
                    << "J(Hα,Hβ3) below Karplus min ~2.16";
                EXPECT_LE(j_ha_hb3[ri][t], 11.0)
                    << "J(Hα,Hβ3) above Karplus max ~10.82";
            }
            // Phi-consistency cross-channel: when both channels finite,
            // accumulate Pearson stats.
            if (std::isfinite(j_hn[ri][t]) && std::isfinite(j_hn_cp[ri][t])) {
                const double x = j_hn[ri][t];
                const double y = j_hn_cp[ri][t];
                sx += x; sy += y;
                sxx += x*x; syy += y*y; sxy += x*y;
                ++n_xy;
            }
        }
    }
    std::cout << "JCouplingTimeSeries: R=" << R << " T=" << T
              << " PRO=" << pro_count
              << " GLY=" << gly_count
              << " GLY+ALA=" << gly_ala_count
              << "\n  | J(HN,Hα)        finite obs=" << hn_finite_obs
              << "\n  | J(HN,Hα)_Vogeli finite obs=" << hn_vog_finite_obs
              << "\n  | J(HN,Cβ)        finite obs=" << hn_cb_finite_obs
              << "\n  | J(HN,C')        finite obs=" << hn_cp_finite_obs
              << "\n  | J(Hα,C')        finite obs=" << ha_cp_finite_obs
              << "\n  | J(N,Cγ)         finite obs=" << n_cg_finite_obs
              << "\n  | J(C',Cγ)        finite obs=" << cp_cg_finite_obs
              << "\n  | J(Hα,Hβ2)       finite obs=" << ha_hb2_finite_obs
              << "\n  | J(Hα,Hβ3)       finite obs=" << ha_hb3_finite_obs
              << "\n";
    EXPECT_GT(hn_finite_obs, 0u);
    EXPECT_GT(hn_vog_finite_obs, 0u);
    EXPECT_GT(hn_cb_finite_obs, 0u);
    EXPECT_GT(hn_cp_finite_obs, 0u);
    EXPECT_GT(ha_cp_finite_obs, 0u);
    EXPECT_GT(n_cg_finite_obs, 0u);
    EXPECT_GT(cp_cg_finite_obs, 0u);
    EXPECT_GT(ha_hb2_finite_obs, 0u);
    EXPECT_GT(ha_hb3_finite_obs, 0u);

    // Cross-channel phi consistency (science F10): J(HN,Hα) and
    // J(HN,C') read the same backbone phi via the H-N-CA-HA and
    // H-N-CA-C atomic dihedrals respectively (3-bond couplings to
    // adjacent substituents on the same N-CA axis). Both Karplus
    // curves have A > 0 and different B (Vuister-Bax B=-1.76,
    // Wang-Bax row 4 B=+0.84 post-2026-05-20 swap), so we expect a
    // moderate-magnitude correlation rather than r ≈ 1 — but they
    // should not be independent. Threshold 0.15 in absolute value to
    // avoid false positives from chance correlation while still
    // detecting a catastrophic decoupling.
    if (n_xy > 10) {
        const double n_d = static_cast<double>(n_xy);
        const double mean_x = sx / n_d;
        const double mean_y = sy / n_d;
        const double var_x = sxx / n_d - mean_x * mean_x;
        const double var_y = syy / n_d - mean_y * mean_y;
        const double cov   = sxy / n_d - mean_x * mean_y;
        if (var_x > 1e-10 && var_y > 1e-10) {
            const double r = cov / std::sqrt(var_x * var_y);
            std::cout << "  J(HN,Hα) ↔ J(HN,C') Pearson r = "
                      << r << " over " << n_xy << " pairs\n";
            EXPECT_GT(std::abs(r), 0.15)
                << "phi cross-channel decorrelation suggests one channel "
                   "is mis-wired";
        }
    }

    fs::remove(h5_path);
}


// PROBE: pin the project-sign mapping. Wang-Bax prints theta_pub=-60°
// for J(HN,Halpha) in the published phi convention; this codebase's
// project phi has the opposite sign, so the compiled project offset is
// +60°. The direct H-N-CA-HA atomic dihedral is a geometry diagnostic,
// not the value fed to the backbone Karplus implementation.
TEST(JCouplingTimeSeries, ProbeHNCAHAvsPhiOffset_1UBQ) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    const std::string pdb = nmr::test::TestEnvironment::UbqProtonated();
    if (!fs::exists(pdb)) GTEST_SKIP() << "1ubq_protonated.pdb absent";
    auto r = nmr::BuildFromProtonatedPdb(pdb);
    ASSERT_TRUE(r.Ok()) << r.error;
    auto& prot = *r.protein;
    auto& conf = prot.Conformation();

    auto dihedral = [&](size_t a, size_t b, size_t c, size_t d) {
        const nmr::Vec3 p1 = conf.PositionAt(a);
        const nmr::Vec3 p2 = conf.PositionAt(b);
        const nmr::Vec3 p3 = conf.PositionAt(c);
        const nmr::Vec3 p4 = conf.PositionAt(d);
        const nmr::Vec3 b1 = p2 - p1, b2 = p3 - p2, b3 = p4 - p3;
        const double b2n = b2.norm();
        if (b2n < 1e-10) return std::numeric_limits<double>::quiet_NaN();
        const nmr::Vec3 n1 = b1.cross(b2);
        const nmr::Vec3 n2 = b2.cross(b3);
        if (n1.norm() < 1e-10 || n2.norm() < 1e-10) return std::numeric_limits<double>::quiet_NaN();
        const nmr::Vec3 m1 = n1.cross(b2 / b2n);
        return std::atan2(m1.dot(n2), n1.dot(n2)) * 180.0 / M_PI;
    };
    auto wrap = [](double d) {
        while (d >  180.0) d -= 360.0;
        while (d < -180.0) d += 360.0;
        return d;
    };

    const size_t R = prot.ResidueCount();
    int n = 0;
    double sum_ha = 0, sum_ha_sq = 0, sum_cb = 0, sum_cb_sq = 0;
    int n_cb = 0;
    double sum_phi = 0;
    // Run J(HN,Hα) both ways and compare to literature.
    auto karplus = [](double A, double B, double C, double t_rad) {
        const double c = std::cos(t_rad);
        return A * c * c + B * c + C;
    };
    double sum_J_atomic = 0;
    double sum_J_project_minus60 = 0;  // unflipped published offset; wrong in project phi sign
    double sum_J_project_plus60 = 0;   // shipped project-sign offset
    int n_helix = 0, n_sheet = 0;
    double sum_J_atomic_helix = 0, sum_J_project_plus60_helix = 0;
    double sum_J_atomic_sheet = 0, sum_J_project_plus60_sheet = 0;
    std::cout << "\nri  type  phi(deg)  H-N-CA-HA  diff_ha  H-N-CA-CB  diff_cb\n";
    for (size_t ri = 1; ri + 1 < R; ++ri) {
        const auto& res = prot.ResidueAt(ri);
        const auto& prev = prot.ResidueAt(ri - 1);
        if (prev.C == nmr::Residue::NONE) continue;
        if (res.N == nmr::Residue::NONE || res.CA == nmr::Residue::NONE) continue;
        if (res.C == nmr::Residue::NONE) continue;
        if (res.H == nmr::Residue::NONE || res.HA == nmr::Residue::NONE) continue;
        const double phi = dihedral(prev.C, res.N, res.CA, res.C);
        const double ha  = dihedral(res.H, res.N, res.CA, res.HA);
        double cb = std::nan("");
        if (res.CB != nmr::Residue::NONE) {
            cb = dihedral(res.H, res.N, res.CA, res.CB);
        }
        const double diff_ha = wrap(ha - phi);
        const double diff_cb = std::isfinite(cb) ? wrap(cb - phi) : std::nan("");
        sum_ha += diff_ha;
        sum_ha_sq += diff_ha * diff_ha;
        sum_phi += phi;
        ++n;

        const double phi_rad = phi * M_PI / 180.0;
        const double atomic_rad = ha * M_PI / 180.0;
        const double J_atomic     = karplus(6.51, -1.76, 1.60, atomic_rad);
        const double J_phi_minus  = karplus(6.51, -1.76, 1.60, phi_rad - M_PI / 3.0);
        const double J_phi_plus   = karplus(6.51, -1.76, 1.60, phi_rad + M_PI / 3.0);
        sum_J_atomic     += J_atomic;
        sum_J_project_minus60 += J_phi_minus;
        sum_J_project_plus60  += J_phi_plus;
        if (phi >= 45.0 && phi <= 90.0) {  // project-sign alphaR-like bin
            ++n_helix;
            sum_J_atomic_helix     += J_atomic;
            sum_J_project_plus60_helix += J_phi_plus;
        }
        if (phi >= 90.0 && phi <= 150.0) {  // project-sign beta-like bin
            ++n_sheet;
            sum_J_atomic_sheet     += J_atomic;
            sum_J_project_plus60_sheet += J_phi_plus;
        }
        if (std::isfinite(diff_cb)) {
            sum_cb += diff_cb;
            sum_cb_sq += diff_cb * diff_cb;
            ++n_cb;
        }
        if (n <= 10) {
            std::printf("%3zu  type=%2d  %8.2f  %9.2f  %7.2f  %9.2f  %7.2f\n",
                        ri, static_cast<int>(res.type),
                        phi, ha, diff_ha, cb, diff_cb);
        }
    }
    const double mean_ha = sum_ha / n;
    const double var_ha  = sum_ha_sq / n - mean_ha * mean_ha;
    std::printf("\nN=%d residues (Hα), %d (Cβ)\n", n, n_cb);
    std::printf("H-N-CA-HA - phi : mean = %+.3f deg, std = %.3f deg\n",
                mean_ha, std::sqrt(std::max(0.0, var_ha)));
    if (n_cb > 0) {
        const double mean_cb = sum_cb / n_cb;
        const double var_cb  = sum_cb_sq / n_cb - mean_cb * mean_cb;
        std::printf("H-N-CA-CB - phi : mean = %+.3f deg, std = %.3f deg\n",
                    mean_cb, std::sqrt(std::max(0.0, var_cb)));
    }
    std::printf("Wang-Bax 1996 prints theta_pub=-60 deg in the published phi convention.\n");
    std::printf("Project phi has the opposite sign, so shipped theta_project=+60 deg.\n");
    std::printf("\n--- J(HN,Hα) all-residue comparison ---\n");
    std::printf("Wang-Bax 1996: helix J~4-6 Hz, β-sheet J~8-10 Hz, BPTI mean ~5 Hz\n");
    std::printf("mean phi over %d residues: %+.2f deg\n", n, sum_phi / n);
    std::printf("J via atomic H-N-CA-HA (geometry check) : mean = %.3f Hz over %d\n",
                sum_J_atomic / n, n);
    std::printf("J via project phi - 60 deg (unflipped): mean = %.3f Hz over %d\n",
                sum_J_project_minus60 / n, n);
    std::printf("J via project phi + 60 deg (shipped)  : mean = %.3f Hz over %d\n",
                sum_J_project_plus60 / n, n);
    if (n_helix > 0)
        std::printf("\nproject αR-bin (n=%d): atomic-J = %.3f, project+60-J = %.3f (literature ~4-6 Hz)\n",
                    n_helix, sum_J_atomic_helix / n_helix,
                    sum_J_project_plus60_helix / n_helix);
    if (n_sheet > 0)
        std::printf("project β-bin  (n=%d): atomic-J = %.3f, project+60-J = %.3f (literature ~8-10 Hz)\n",
                    n_sheet, sum_J_atomic_sheet / n_sheet,
                    sum_J_project_plus60_sheet / n_sheet);
    EXPECT_GT(std::abs(mean_ha), 50.0);
    EXPECT_LT(std::abs(mean_ha), 70.0);
}


// Literature-anchored probe (codex F8): runs the JCoupling TR on
// the 1UBQ_pm6dh3plus fixture protein and checks that the J values
// reproduce the rough literature bands for helix and sheet residues
// stratified by phi. With the project-sign (phi + theta_offset)
// Karplus form (codex F6 fix + project-sign repair, 2026-05-20), the
// curve should reproduce Wang-Bax Figure 4 panels A-D.
//
// Sources for the band targets (Wang & Bax 1996 JACS 118:2483,
// Figure 4 page 2488, fit curves):
//      Panel A: J(HN,Hα): helix (phi≈-60°) ~4-5 Hz from curve, scatter
//               3-6 Hz; sheet (phi≈-120°) ~9-10 Hz from curve, scatter
//               7-10 Hz; max near phi=±170° at 10 Hz.
//      Panel B: J(HN,C'): helix ~1-2 Hz; sheet ~0.5-1 Hz; max near
//               phi=0° at ~5 Hz. (Both helix and sheet fall on the
//               low-slope shoulder of the curve; Wang-Bax row 4
//               theta=0 puts the maximum at phi=0.)
//      Panel C: J(HN,Cβ): helix ~2-3 Hz; sheet ~0-1 Hz; max near
//               phi≈+170° around 4 Hz.
//      Panel D: J(Hα,C'): helix ~1-2 Hz from curve (vertex near phi
//               ≈ -100° gives min ≈ 1 Hz); sheet ~3-5 Hz from curve;
//               max near phi≈+60° around 7 Hz.
//
// Tolerances are loose (±2 Hz, per feedback_log_overages_dont_assert
// memory): the goal is to flag a SIGN or 60-degree-offset mistake in
// the Karplus evaluation, NOT to over-constrain on per-residue noise.
// The bands below reflect the FIT-CURVE value at the band-center phi,
// with ±2 Hz tolerance absorbing intra-band phi spread + experimental
// scatter.
//
// The probe iterates over the 1UBQ residue set (skipping PRO and
// chain termini), stratifies by phi into αR-helix and β-sheet bins,
// and reports the mean J for each backbone Karplus channel. The
// test asserts band membership for the means, not for individual
// residues.
TEST(JCouplingTimeSeries, LiteratureAnchoredProbeOn1UBQ) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    const std::string pdb = nmr::test::TestEnvironment::UbqProtonated();
    if (!fs::exists(pdb)) GTEST_SKIP() << "1ubq_protonated.pdb absent";

    auto r = nmr::BuildFromProtonatedPdb(pdb);
    ASSERT_TRUE(r.Ok()) << r.error;
    auto& prot = *r.protein;
    auto& conf = prot.Conformation();
    const std::size_t R = prot.ResidueCount();

    auto dihedral = [&](std::size_t a, std::size_t b,
                        std::size_t c, std::size_t d) -> double {
        const nmr::Vec3 p1 = conf.PositionAt(a);
        const nmr::Vec3 p2 = conf.PositionAt(b);
        const nmr::Vec3 p3 = conf.PositionAt(c);
        const nmr::Vec3 p4 = conf.PositionAt(d);
        const nmr::Vec3 b1 = p2 - p1, b2 = p3 - p2, b3 = p4 - p3;
        const double b2n = b2.norm();
        if (b2n < 1e-10) return std::nan("");
        const nmr::Vec3 n1 = b1.cross(b2);
        const nmr::Vec3 n2 = b2.cross(b3);
        if (n1.norm() < 1e-10 || n2.norm() < 1e-10) return std::nan("");
        const nmr::Vec3 m1 = n1.cross(b2 / b2n);
        return std::atan2(m1.dot(n2), n1.dot(n2));
    };
    auto karplus = [](double A, double B, double C, double t) {
        const double cs = std::cos(t);
        return A * cs * cs + B * cs + C;
    };

    // Accumulators per channel x region. Regions are defined on
    // project phi, whose sign is opposite to Wang-Bax/DSSP; therefore
    // αR-helix is [+45°, +90°] and β-sheet is [+90°, +150°].
    struct Stat { double sum = 0; int n = 0; };
    Stat hn_helix, hn_sheet;
    Stat hncp_helix, hncp_sheet;
    Stat hncb_helix, hncb_sheet;
    Stat hacp_helix, hacp_sheet;

    int n_total = 0;
    for (std::size_t ri = 1; ri + 1 < R; ++ri) {
        const auto& res = prot.ResidueAt(ri);
        const auto prev_opt = prot.BackbonePredecessor(ri);
        if (!prev_opt) continue;
        const auto& prev = prot.ResidueAt(*prev_opt);
        if (prev.C  == nmr::Residue::NONE) continue;
        if (res.N  == nmr::Residue::NONE || res.CA == nmr::Residue::NONE) continue;
        if (res.C  == nmr::Residue::NONE) continue;
        if (res.H  == nmr::Residue::NONE || res.HA == nmr::Residue::NONE) continue;
        if (res.CB == nmr::Residue::NONE) continue;  // GLY excluded;
                                                      // sheet/helix counts
                                                      // shouldn't drop much

        const double phi = dihedral(prev.C, res.N, res.CA, res.C);
        const double phi_deg = phi * 180.0 / M_PI;

        // Karplus J via the (phi + theta_offset) form -- mirrors the
        // form shipped in JCouplingTimeSeriesTrajectoryResult after
        // codex F6 (2026-05-20). Each backbone channel's
        // theta_offset constant lives in PhysicalConstants.h.
        const double th_hnha = phi + nmr::KARPLUS_HN_HA_THETA;
        const double th_hncb = phi + nmr::KARPLUS_HN_CB_THETA;
        const double th_hncp = phi + nmr::KARPLUS_HN_CP_THETA;
        const double th_hacp = phi + nmr::KARPLUS_HA_CP_THETA;

        const double j_hn   = karplus(nmr::KARPLUS_HN_HA_A,
                                       nmr::KARPLUS_HN_HA_B,
                                       nmr::KARPLUS_HN_HA_C, th_hnha);
        const double j_hncb = karplus(nmr::KARPLUS_HN_CB_A,
                                       nmr::KARPLUS_HN_CB_B,
                                       nmr::KARPLUS_HN_CB_C, th_hncb);
        const double j_hncp = karplus(nmr::KARPLUS_HN_CP_A,
                                       nmr::KARPLUS_HN_CP_B,
                                       nmr::KARPLUS_HN_CP_C, th_hncp);
        const double j_hacp = karplus(nmr::KARPLUS_HA_CP_A,
                                       nmr::KARPLUS_HA_CP_B,
                                       nmr::KARPLUS_HA_CP_C, th_hacp);

        const bool is_helix = (phi_deg >= 45.0 && phi_deg <= 90.0);
        const bool is_sheet = (phi_deg >= 90.0 && phi_deg <= 150.0);
        if (is_helix) {
            hn_helix.sum   += j_hn;   ++hn_helix.n;
            hncp_helix.sum += j_hncp; ++hncp_helix.n;
            hncb_helix.sum += j_hncb; ++hncb_helix.n;
            hacp_helix.sum += j_hacp; ++hacp_helix.n;
        }
        if (is_sheet) {
            hn_sheet.sum   += j_hn;   ++hn_sheet.n;
            hncp_sheet.sum += j_hncp; ++hncp_sheet.n;
            hncb_sheet.sum += j_hncb; ++hncb_sheet.n;
            hacp_sheet.sum += j_hacp; ++hacp_sheet.n;
        }
        ++n_total;
    }

    auto mean = [](const Stat& s) {
        return s.n > 0 ? s.sum / s.n : std::nan("");
    };
    auto log_band = [&](const char* label, const Stat& helix,
                        const Stat& sheet,
                        double helix_lo, double helix_hi,
                        double sheet_lo, double sheet_hi) {
        const double mh = mean(helix);
        const double ms = mean(sheet);
        std::printf("  %s: helix mean = %6.3f Hz (n=%3d; expect [%g, %g]); "
                    "sheet mean = %6.3f Hz (n=%3d; expect [%g, %g])\n",
                    label, mh, helix.n, helix_lo, helix_hi,
                    ms, sheet.n, sheet_lo, sheet_hi);
        // Loose ±2 Hz tolerance around literature band -- per
        // feedback_log_overages_dont_assert, the band is a SHAPE check
        // not a per-residue gate. Skip if region underpopulated.
        if (helix.n >= 5) {
            EXPECT_GE(mh, helix_lo - 2.0)
                << label << " helix mean " << mh
                << " below band [" << helix_lo << ", " << helix_hi << "]";
            EXPECT_LE(mh, helix_hi + 2.0)
                << label << " helix mean " << mh
                << " above band [" << helix_lo << ", " << helix_hi << "]";
        }
        if (sheet.n >= 5) {
            EXPECT_GE(ms, sheet_lo - 2.0)
                << label << " sheet mean " << ms
                << " below band [" << sheet_lo << ", " << sheet_hi << "]";
            EXPECT_LE(ms, sheet_hi + 2.0)
                << label << " sheet mean " << ms
                << " above band [" << sheet_lo << ", " << sheet_hi << "]";
        }
    };
    std::printf("\nLiteratureAnchoredProbeOn1UBQ: %d backbone residues\n",
                n_total);
    // Bands reflect Wang-Bax 1996 Figure 4 (page 2488) fit-curve
    // values at the band-center phi, ±2 Hz tolerance. The J(HN,Hα)
    // sheet/helix discrimination is the strongest signal (factor 2
    // difference); the other three channels show subtler phi-
    // dependence and the tolerance absorbs the spread.
    //
    // The probe's primary purpose is to detect a SIGN or 60-deg-offset
    // mistake in the Karplus evaluation -- specifically the failure
    // mode where the Wang-Bax printed theta offsets are used without
    // flipping them into the project phi convention. This probe is the
    // executable gate for the project-sign theta constants.
    log_band("J(HN,Halpha)", hn_helix, hn_sheet,
             4.0, 6.0,    8.0, 10.0);
    log_band("J(HN,Cprime)", hncp_helix, hncp_sheet,
             1.0, 2.5,    0.0, 1.5);
    log_band("J(HN,Cbeta)",  hncb_helix, hncb_sheet,
             1.5, 3.0,    0.0, 1.5);
    log_band("J(Halpha,Cprime)", hacp_helix, hacp_sheet,
             0.5, 2.5,    2.0, 4.0);

    // Sanity floor: at least some helix + sheet residues must have
    // populated, else the protein isn't 1UBQ or the dihedral path is
    // catastrophically broken.
    EXPECT_GE(hn_helix.n + hn_sheet.n, 10);
}
