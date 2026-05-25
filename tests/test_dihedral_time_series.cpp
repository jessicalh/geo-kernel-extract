//
// test_dihedral_time_series: discipline + integration for
// DihedralTimeSeriesTrajectoryResult. Per-residue per-frame
// phi/psi/omega/chi timelines + Ramachandran-region classification +
// static per-residue masks. Self-contained TR (positions + topology
// only — no source ConformationResult dependency).
//

#include "AminoAcidType.h"
#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DihedralTimeSeriesTrajectoryResult.h"
#include "DsspResult.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "OperationLog.h"
#include "PlanarGeometryResult.h"
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
#include <limits>
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
nmr::RunConfiguration BuildConfig(unsigned stride) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    // DihedralTS reads positions only — strip everything heavy.
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::DihedralTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

// Heavier config: enables DSSP + Enrichment so PlanarGeometryResult
// attaches. Used by the cross-result-consistency test which validates
// our independently-computed phi/psi/omega against DsspResult and
// PlanarGeometryResult at 1e-6 tolerance on internal residues.
nmr::RunConfiguration BuildConfigWithCrossResults(unsigned stride) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true;
    opts.skip_dssp = false;  // need DSSP for cross-check
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::EnrichmentResult));
    config.RequireConformationResult(typeid(nmr::DsspResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::DihedralTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


// ── Frame 0 smoke ────────────────────────────────────────────────────

TEST(DihedralTimeSeries, Frame0Semantics) {
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
    EXPECT_EQ(traj.FrameCount(), 1u);

    const auto& tr = tp.Result<nmr::DihedralTimeSeriesTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);
}


// ── Finalize idempotency (data-flow short-circuit) ───────────────────

TEST(DihedralTimeSeries, FinalizeIdempotency) {
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

    auto& tr = tp.Result<nmr::DihedralTimeSeriesTrajectoryResult>();
    const std::size_t T = tr.NumFrames();
    tr.Finalize(tp, traj);
    EXPECT_EQ(tr.NumFrames(), T);
}


// ── H5 round-trip: schema + attrs + dataset shapes ───────────────────

TEST(DihedralTimeSeries, H5RoundTrip) {
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

    const auto& tr = tp.Result<nmr::DihedralTimeSeriesTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    const std::size_t T = tr.NumFrames();

    const std::string h5_path = (fs::temp_directory_path() /
        ("dihedral_ts_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/dihedral_time_series"));
    auto grp = reopen.getGroup("/trajectory/dihedral_time_series");

    // Convention attrs.
    std::string units, periodicity, convention;
    grp.getAttribute("angle_units").read(units);
    grp.getAttribute("periodicity").read(periodicity);
    grp.getAttribute("angle_convention").read(convention);
    EXPECT_EQ(units, "radians");
    EXPECT_EQ(periodicity, "2pi");
    EXPECT_NE(convention.find("IUPAC"), std::string::npos);

    // Per-frame datasets.
    for (const std::string& name : {"phi", "psi", "omega", "omega_deviation"}) {
        ASSERT_TRUE(grp.exist(name)) << name;
        const auto dims = grp.getDataSet(name).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], T);
    }
    {
        ASSERT_TRUE(grp.exist("chi"));
        const auto dims = grp.getDataSet("chi").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 3u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], T);
        EXPECT_EQ(dims[2], 4u);
    }
    {
        ASSERT_TRUE(grp.exist("rama_region"));
        const auto dims = grp.getDataSet("rama_region").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], T);
    }

    // Static masks.
    EXPECT_TRUE(grp.exist("chi_exists"));
    EXPECT_TRUE(grp.exist("omega_is_xpro"));
    EXPECT_TRUE(grp.exist("is_glycine"));
    EXPECT_TRUE(grp.exist("is_proline"));
    EXPECT_TRUE(grp.exist("is_pre_proline"));
    EXPECT_TRUE(grp.exist("residue_terminal_state"));
    EXPECT_TRUE(grp.exist("chain_id_per_residue"));
    EXPECT_TRUE(grp.exist("residue_index_per_atom"));
    EXPECT_TRUE(grp.exist("frame_indices"));
    EXPECT_TRUE(grp.exist("frame_times"));
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));

    // n_atoms attr validation: must match dataset shape.
    std::size_t n_atoms_attr = 0;
    grp.getAttribute("n_atoms").read(n_atoms_attr);
    EXPECT_EQ(n_atoms_attr, tp.AtomCount());
    std::vector<std::int32_t> ria;
    grp.getDataSet("residue_index_per_atom").read(ria);
    EXPECT_EQ(ria.size(), n_atoms_attr);

    // chain_id_per_residue read-back must round-trip variable-length strings.
    std::vector<std::string> chain_ids;
    grp.getDataSet("chain_id_per_residue").read(chain_ids);
    ASSERT_EQ(chain_ids.size(), R);
    for (std::size_t ri = 0; ri < R; ++ri) {
        EXPECT_EQ(chain_ids[ri], tp.ProteinRef().ResidueAt(ri).chain_id)
            << "chain_id mismatch at ri=" << ri;
    }

    // Convention pin attrs — both new ones added in cleanup pass.
    std::string value_range, chunking_policy, source_policy;
    grp.getAttribute("value_range").read(value_range);
    grp.getAttribute("chunking_policy").read(chunking_policy);
    grp.getAttribute("source_attached_policy").read(source_policy);
    EXPECT_NE(value_range.find("[-pi, pi]"), std::string::npos);
    EXPECT_NE(chunking_policy.find("{R, min(T, 64)}"), std::string::npos);
    EXPECT_NE(source_policy.find("always_attached"), std::string::npos);

    fs::remove(h5_path);
}


// ── Integration on 1P9J: real data, real boundaries ──────────────────

TEST(DihedralTimeSeries, Integration1P9J) {
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

    const auto& tr = tp.Result<nmr::DihedralTimeSeriesTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    const std::size_t T = tr.NumFrames();
    EXPECT_GE(T, 2u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("dihedral_ts_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/dihedral_time_series");

    // Boundary-NaN discipline: N-terminus residue must have NaN phi.
    // C-terminus must have NaN psi and omega.
    std::vector<std::vector<double>> phi_data, psi_data, omega_data;
    grp.getDataSet("phi").read(phi_data);
    grp.getDataSet("psi").read(psi_data);
    grp.getDataSet("omega").read(omega_data);
    ASSERT_EQ(phi_data.size(), R);
    ASSERT_EQ(phi_data[0].size(), T);

    std::vector<std::uint8_t> terminal_state;
    grp.getDataSet("residue_terminal_state").read(terminal_state);
    ASSERT_EQ(terminal_state.size(), R);

    std::size_t n_term_seen = 0, c_term_seen = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        const std::uint8_t ts = terminal_state[ri];
        // 1 = NTerminus, 2 = CTerminus, 3 = NAndCTerminus
        if (ts == 1u || ts == 3u) {
            EXPECT_TRUE(std::isnan(phi_data[ri][0]))
                << "N-terminus residue " << ri << " has finite phi";
            ++n_term_seen;
        }
        if (ts == 2u || ts == 3u) {
            EXPECT_TRUE(std::isnan(psi_data[ri][0]))
                << "C-terminus residue " << ri << " has finite psi";
            EXPECT_TRUE(std::isnan(omega_data[ri][0]))
                << "C-terminus residue " << ri << " has finite omega";
            ++c_term_seen;
        }
    }
    std::cout << "DihedralTimeSeries: " << R << " residues, " << T << " frames; "
              << "termini observed: N=" << n_term_seen
              << " C=" << c_term_seen << "\n";

    // chi_exists static: GLY rows should be all-zero, PHE/TYR rows
    // should be 1 in slots 0+1.
    std::vector<std::vector<std::uint8_t>> chi_exists_data;
    grp.getDataSet("chi_exists").read(chi_exists_data);
    ASSERT_EQ(chi_exists_data.size(), R);
    std::size_t gly_count = 0, two_chi_count = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        const auto& row = chi_exists_data[ri];
        ASSERT_EQ(row.size(), 4u);
        const nmr::Residue& res = tp.ProteinRef().ResidueAt(ri);
        if (res.type == nmr::AminoAcid::GLY) {
            EXPECT_EQ(row[0] + row[1] + row[2] + row[3], 0)
                << "GLY residue " << ri << " has chi_exists != all-zero";
            ++gly_count;
        }
        if (res.type == nmr::AminoAcid::PHE ||
            res.type == nmr::AminoAcid::TYR ||
            res.type == nmr::AminoAcid::HIS ||
            res.type == nmr::AminoAcid::TRP) {
            EXPECT_EQ(row[0], 1u);
            EXPECT_EQ(row[1], 1u);
            ++two_chi_count;
        }
    }
    std::cout << "  GLY residues (all-zero chi_exists): " << gly_count
              << "; aromatic residues (>=2 chi): " << two_chi_count << "\n";

    // omega_is_xpro mask cross-check: for each PRO residue at ri > 0,
    // assert omega_is_xpro[ri-1] == 1 AND is_pre_proline[ri-1] == 1 AND
    // is_proline[ri] == 1 (verifies the static-mask wiring done in
    // Create() picked up Pro correctly).
    std::vector<std::uint8_t> omega_xpro, pre_pro, is_pro;
    grp.getDataSet("omega_is_xpro").read(omega_xpro);
    grp.getDataSet("is_pre_proline").read(pre_pro);
    grp.getDataSet("is_proline").read(is_pro);
    ASSERT_EQ(omega_xpro.size(), R);
    ASSERT_EQ(pre_pro.size(),    R);
    ASSERT_EQ(is_pro.size(),     R);
    std::size_t pro_count = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        const nmr::Residue& res = tp.ProteinRef().ResidueAt(ri);
        if (res.type == nmr::AminoAcid::PRO) {
            EXPECT_EQ(is_pro[ri], 1u);
            if (ri > 0) {
                // Pre-Pro flag should be set on i-1 IF backbone-connected.
                // 1P9J is single-chain with no numbering gaps so we expect
                // every Pro to have a connected predecessor.
                EXPECT_EQ(omega_xpro[ri - 1], 1u)
                    << "omega_is_xpro[" << (ri-1) << "] should be 1 because "
                    << "residue " << ri << " is PRO";
                EXPECT_EQ(pre_pro[ri - 1], 1u);
            }
            ++pro_count;
        }
    }
    std::cout << "  PRO residues: " << pro_count
              << "; omega_is_xpro+is_pre_proline cross-check passed\n";

    // omega_deviation now matches PlanarGeometryResult impl: emitted as
    // actual WrapPi(omega-pi) for EVERY well-defined peptide bond
    // (including X->Pro). Verify by sampling: omega_is_xpro==1 rows
    // should have finite omega_deviation when omega is finite.
    std::vector<std::vector<double>> omega_dev_data;
    grp.getDataSet("omega_deviation").read(omega_dev_data);
    std::size_t xpro_dev_finite = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        if (omega_xpro[ri] != 1u) continue;
        for (std::size_t t = 0; t < T; ++t) {
            if (std::isfinite(omega_data[ri][t]) &&
                std::isfinite(omega_dev_data[ri][t])) {
                ++xpro_dev_finite;
                EXPECT_GE(omega_dev_data[ri][t], -M_PI);
                EXPECT_LE(omega_dev_data[ri][t],  M_PI);
            }
        }
    }
    std::cout << "  X-Pro rows with finite omega_deviation: " << xpro_dev_finite
              << " (post-cleanup: actual values, not NaN-filled)\n";

    // Range sanity on emitted dihedrals. atan2 returns [-pi, pi]
    // (closed both ends). Use the inclusive bound.
    std::size_t phi_finite = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        for (std::size_t t = 0; t < T; ++t) {
            const double v = phi_data[ri][t];
            if (std::isfinite(v)) {
                ++phi_finite;
                EXPECT_GE(v, -M_PI) << "phi out of [-pi, pi] at (" << ri << "," << t << ")";
                EXPECT_LE(v,  M_PI) << "phi out of [-pi, pi] at (" << ri << "," << t << ")";
            }
        }
    }
    EXPECT_GE(phi_finite, R - 5u)
        << "Expected most rows to have at least one finite phi sample";

    fs::remove(h5_path);
}


// ── Cross-result consistency: DSSP phi/psi + PG omega ────────────────
//
// Validates that our independently-computed dihedrals agree with
// DsspResult and PlanarGeometryResult on internal residues.
//
// Two distinct convention layers checked:
//   1. PlanarGeometryResult uses IUPAC signed dihedral via the same
//      atan2(y, x) formula — bit-identical agreement expected for
//      omega (and for phi/psi if PG had them).
//   2. DsspResult forwards values from libdssp/libcifpp, which uses
//      the NEGATED IUPAC sign convention for backbone dihedrals
//      (well-known DSSP quirk; phi_DSSP = -phi_IUPAC, psi_DSSP =
//      -psi_IUPAC). We compare to -DSSP.Phi(ri) / -DSSP.Psi(ri) and
//      document the convention divergence in code + test comment so
//      downstream consumers (and reviewers) see the trap. Drift OTHER
//      than the sign flip would mean an algorithmic bug.

TEST(DihedralTimeSeries, CrossResultConsistencyDsspPlanarGeometry) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfigWithCrossResults(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& conf = tp.CanonicalConformation();
    ASSERT_TRUE(conf.HasResult<nmr::DsspResult>());
    const auto& dssp = conf.Result<nmr::DsspResult>();
    const bool has_pg = conf.HasResult<nmr::PlanarGeometryResult>();

    const auto& dts = tp.Result<nmr::DihedralTimeSeriesTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();

    const std::string h5_path = (fs::temp_directory_path() /
        ("dihedral_ts_cross_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      dts.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/dihedral_time_series");

    std::vector<std::vector<double>> phi, psi, omega;
    grp.getDataSet("phi").read(phi);
    grp.getDataSet("psi").read(psi);
    grp.getDataSet("omega").read(omega);

    // PG uses the same atan2 formula as DTS — bit-identical agreement
    // expected (1e-12 tolerance).
    //
    // DSSP has TWO convention layers vs IUPAC:
    //   (a) negated sign: phi_DSSP = -phi_IUPAC (compare to -DSSP value)
    //   (b) "undefined" sentinel is ±2π (degrees=360) — NOT NaN. Filter
    //       residues where |DSSP value| > 2π - 0.1 as undefined.
    // Tolerance ~0.01 rad (~0.6°) absorbs libdssp's algorithmic drift
    // (different cross-product accumulation order can produce O(1e-3)
    // numerical differences even though both use atan2-of-cross-product).
    // Tight enough to catch a sign flip (would be O(π) diff); loose
    // enough to ignore the implementation drift.
    constexpr double kTolPg   = 1e-12;
    constexpr double kTolDssp = 1e-2;
    constexpr double kDsspUndefSentinel = 2.0 * M_PI - 0.1;

    auto dssp_finite = [](double d) {
        return std::isfinite(d) && std::abs(d) < kDsspUndefSentinel;
    };

    std::size_t phi_match = 0, phi_skip = 0, phi_mismatch = 0;
    std::size_t psi_match = 0, psi_skip = 0, psi_mismatch = 0;
    std::size_t omega_match = 0, omega_skip = 0, omega_mismatch = 0;

    for (std::size_t ri = 0; ri < R; ++ri) {
        const double dts_phi   = phi[ri][0];
        const double dts_psi   = psi[ri][0];
        const double dts_omega = omega[ri][0];

        if (std::isfinite(dts_phi)) {
            // DSSP convention: negated IUPAC (well-known libdssp/libcifpp
            // quirk). Compare to -DSSP.Phi. Skip rows where DSSP returns
            // its "undefined" sentinel (±2π, from 360-degree internal).
            const double d = dssp.Phi(ri);
            if (dssp_finite(d) && std::abs(d) > 1e-12) {
                const double d_neg = -d;
                if (std::abs(dts_phi - d_neg) < kTolDssp) ++phi_match;
                else { ++phi_mismatch;
                    EXPECT_NEAR(dts_phi, d_neg, kTolDssp)
                        << "phi diverged at ri=" << ri
                        << " (vs -DSSP — convention is negated IUPAC; "
                        << "tolerance absorbs libdssp algorithmic drift)"; }
            } else { ++phi_skip; }
        }
        if (std::isfinite(dts_psi)) {
            const double d = dssp.Psi(ri);
            if (dssp_finite(d) && std::abs(d) > 1e-12) {
                const double d_neg = -d;
                if (std::abs(dts_psi - d_neg) < kTolDssp) ++psi_match;
                else { ++psi_mismatch;
                    EXPECT_NEAR(dts_psi, d_neg, kTolDssp)
                        << "psi diverged at ri=" << ri
                        << " (vs -DSSP — convention is negated IUPAC; "
                        << "tolerance absorbs libdssp algorithmic drift)"; }
            } else { ++psi_skip; }
        }

        if (has_pg && std::isfinite(dts_omega)) {
            const auto& pg_omega =
                conf.Result<nmr::PlanarGeometryResult>().OmegaActual();
            if (ri < pg_omega.size() && std::isfinite(pg_omega[ri])) {
                // PG uses the same atan2 formula — bit-identical match
                // expected.
                if (std::abs(dts_omega - pg_omega[ri]) < kTolPg) ++omega_match;
                else { ++omega_mismatch;
                    EXPECT_NEAR(dts_omega, pg_omega[ri], kTolPg)
                        << "omega diverged from PlanarGeometryResult at ri=" << ri; }
            } else { ++omega_skip; }
        }
    }

    std::cout << "Cross-result consistency on " << R << " residues:\n"
              << "  phi:   match=" << phi_match
              << " mismatch=" << phi_mismatch
              << " skip=" << phi_skip << " (DSSP)\n"
              << "  psi:   match=" << psi_match
              << " mismatch=" << psi_mismatch
              << " skip=" << psi_skip << " (DSSP)\n"
              << "  omega: match=" << omega_match
              << " mismatch=" << omega_mismatch
              << " skip=" << omega_skip
              << " (PlanarGeometryResult"
              << (has_pg ? "" : " — NOT attached")
              << ")\n";

    // Loose floor — at least most residues should match. Per
    // feedback_log_overages_dont_assert we don't assert exact counts;
    // the EXPECT_NEAR per-row above is what enforces correctness.
    EXPECT_EQ(phi_mismatch,   0u) << "phi divergence vs DSSP";
    EXPECT_EQ(psi_mismatch,   0u) << "psi divergence vs DSSP";
    EXPECT_EQ(omega_mismatch, 0u) << "omega divergence vs PlanarGeometryResult";
    EXPECT_GT(phi_match + psi_match, R)
        << "DSSP cross-check produced too few matches (DSSP attached but "
        << "returning NaN/zero everywhere?)";

    fs::remove(h5_path);
}


// Documented coverage gap: Protein::BackbonePredecessor / BackboneSuccessor
// walk the cifpp bond graph and return the actual covalent neighbour
// (handles cyclic-peptide wrap, ACE/NME caps, antibody insertion codes,
// and engineered chimeras with non-monotonic numbering uniformly --
// none of those branches require special-case logic at the call site).
// 1P9J is single-chain linear with monotonic numbering and no caps;
// the bond walk and the linear ri±1 path return the same answer. A
// multi-chain / cyclic / insertion-coded fixture would distinguish
// them numerically — see B's MED 4 in the dual adversarial review and
// the cyclic-peptide caveat in Protein.h.
