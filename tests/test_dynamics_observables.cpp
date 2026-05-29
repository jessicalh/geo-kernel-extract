//
// test_dynamics_observables: discipline + integration for the three
// positions-only model-free TrajectoryResults --
// ReorientationalDynamics, IRedOrderParameter, DihedralAutocorrelation.
// All read positions + the Residue backbone/chi cache only (no source
// ConformationResult), so the per-frame stack is stripped to nothing and
// the test runs without APBS / AIMNet2 / GROMACS-energy machinery.
//
// KernelDynamics + KernelCoherence depend on the classical kernel
// ConformationResults incl. APBS, so they are exercised on the production
// config verbatim -- RunConfiguration::PerFrameExtractionSet() (APBS on,
// AIMNet2 on, MOPAC the only thing skipped), which already Produces both --
// in the KernelInstrumentPerFrameExtractionSet test below. That test skips
// only if the AIMNet2 model is unavailable, like the other AIMNet2 tests.
//
// Engine numerics (biased ACF, Parzen spectrum, Legendre/circular TCF)
// are unit-verified separately; this file verifies the TRs assemble
// correct output end-to-end on the 1P9J fixture.
//

#include "DihedralAutocorrelationTrajectoryResult.h"
#include "IRedOrderParameterTrajectoryResult.h"
#include "KernelCoherenceTrajectoryResult.h"
#include "KernelDynamicsTrajectoryResult.h"
#include "OperationLog.h"
#include "ReorientationalDynamicsTrajectoryResult.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"

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

// Attaches all three positions-only model-free TRs; strips the per-frame
// stack (they need no ConformationResult).
nmr::RunConfiguration BuildConfig(unsigned stride) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true;
    opts.skip_coulomb = true;
    opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::ReorientationalDynamicsTrajectoryResult::Create(tp);
        });
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::IRedOrderParameterTrajectoryResult::Create(tp);
        });
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp) -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::DihedralAutocorrelationTrajectoryResult::Create(tp);
        });
    config.SetStride(stride);
    return config;
}

// Runs the fixture once and returns the TrajectoryProtein with the three
// TRs finalized. Caller already checked FixtureAvailable.
bool RunFixture(nmr::TrajectoryProtein& tp, unsigned stride) {
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) return false;
    auto config = BuildConfig(stride);
    if (!tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) return false;
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    return traj.Run(tp, config, session) == nmr::kOk;
}

}  // namespace


// ── Frame 0 smoke: the three TRs attach, run, and populate ──────────

TEST(DynamicsObservables, Frame0Semantics) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(RunFixture(tp, 99999u));  // huge stride -> one frame

    const auto& reo = tp.Result<nmr::ReorientationalDynamicsTrajectoryResult>();
    const auto& ired = tp.Result<nmr::IRedOrderParameterTrajectoryResult>();
    const auto& dih = tp.Result<nmr::DihedralAutocorrelationTrajectoryResult>();
    // 1P9J is a 36-residue protein: backbone N-H/CA-HA/C=O vectors exist,
    // amide N-H exist, and phi/psi/chi are defined for interior residues.
    EXPECT_GT(reo.NumVectors(), 0u);
    EXPECT_GT(ired.NumVectors(), 0u);
    EXPECT_GT(dih.NumFrames(), 0u);
}


// ── Finalize idempotency (the codex 2026-05-29 fix) ─────────────────

TEST(DynamicsObservables, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(RunFixture(tp, 300u));
    nmr::Trajectory dummy("", "", "");  // Finalize ignores traj here

    auto& reo = tp.Result<nmr::ReorientationalDynamicsTrajectoryResult>();
    auto& ired = tp.Result<nmr::IRedOrderParameterTrajectoryResult>();
    auto& dih = tp.Result<nmr::DihedralAutocorrelationTrajectoryResult>();
    const std::size_t reo_v = reo.NumVectors();
    const std::size_t ired_v = ired.NumVectors();
    const std::size_t dih_t = dih.NumFrames();
    // Second Finalize must be a safe no-op (accumulators were released).
    reo.Finalize(tp, dummy);
    ired.Finalize(tp, dummy);
    dih.Finalize(tp, dummy);
    EXPECT_EQ(reo.NumVectors(), reo_v);
    EXPECT_EQ(ired.NumVectors(), ired_v);
    EXPECT_EQ(dih.NumFrames(), dih_t);
}


// ── H5 round-trip + physics content ─────────────────────────────────

TEST(DynamicsObservables, H5RoundTripAndContent) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(RunFixture(tp, 300u));  // ~5 frames over the 1501-frame fixture

    const auto& reo = tp.Result<nmr::ReorientationalDynamicsTrajectoryResult>();
    const auto& ired = tp.Result<nmr::IRedOrderParameterTrajectoryResult>();
    const auto& dih = tp.Result<nmr::DihedralAutocorrelationTrajectoryResult>();

    const std::string h5_path = (fs::temp_directory_path() /
        ("dynobs_" + std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        reo.WriteH5Group(tp, file);
        ired.WriteH5Group(tp, file);
        dih.WriteH5Group(tp, file);
    }
    HighFive::File f(h5_path, HighFive::File::ReadOnly);

    // ── Reorientational ──────────────────────────────────────────────
    ASSERT_TRUE(f.exist("/trajectory/reorientational_dynamics"));
    auto rg = f.getGroup("/trajectory/reorientational_dynamics");
    {
        std::vector<double> s2;
        rg.getDataSet("order_parameter_S2").read(s2);
        EXPECT_GT(s2.size(), 0u);
        for (double v : s2) {
            // S^2 is a generalized order parameter: physically [0,1]; allow
            // a small finite-sampling excursion. Must be finite.
            ASSERT_TRUE(std::isfinite(v)) << "non-finite S^2";
            EXPECT_GE(v, -0.2);
            EXPECT_LE(v, 1.05);
        }
        // Internal TCF C_I(0) == 1 by construction for every vector.
        std::vector<std::vector<double>> acf;
        rg.getDataSet("bond_vector_autocorrelation").read(acf);
        for (const auto& row : acf)
            EXPECT_NEAR(row.at(0), 1.0, 1e-9) << "C_I(0) must be 1";
        // Orientation tensor shape (V, 3, 3).
        const auto od = rg.getDataSet("bond_orientation_tensor").getSpace().getDimensions();
        ASSERT_EQ(od.size(), 3u);
        EXPECT_EQ(od[1], 3u);
        EXPECT_EQ(od[2], 3u);
        // The 1P9J fixture is ~15 ns -> a few tau_m only -> tau_m NOT
        // converged. The honesty flag must say so.
        bool converged = true;
        rg.getAttribute("tau_m_converged").read(converged);
        EXPECT_FALSE(converged)
            << "15 ns fixture should report tau_m unconverged";

        // ── 15N relaxation layer (NH rows) ───────────────────────────
        // Compute-and-flag: rates ARE emitted on the unconverged fixture,
        // gated for reliability by tau_m_converged (asserted false above).
        // Verify the code path: NH rows finite + physical sign, J(w) >= 0,
        // non-NH rows NaN, reporting field round-trips.
        std::vector<std::uint8_t> kind;
        rg.getDataSet("vector_kind").read(kind);
        std::vector<double> R1, R2, NOE;
        rg.getDataSet("relaxation_R1").read(R1);
        rg.getDataSet("relaxation_R2").read(R2);
        rg.getDataSet("relaxation_NOE").read(NOE);
        std::vector<std::vector<double>> jw;
        rg.getDataSet("spectral_density_j").read(jw);
        ASSERT_EQ(R1.size(), kind.size());
        ASSERT_EQ(jw.size(), kind.size());
        double field = 0.0;
        rg.getAttribute("relaxation_field_tesla").read(field);
        EXPECT_NEAR(field, 14.1, 1e-6);
        std::size_t n_nh = 0;
        for (std::size_t v = 0; v < kind.size(); ++v) {
            if (kind[v] == 1u) {  // NH
                ++n_nh;
                ASSERT_TRUE(std::isfinite(R1[v])) << "NH R1 must be finite";
                ASSERT_TRUE(std::isfinite(R2[v])) << "NH R2 must be finite";
                ASSERT_TRUE(std::isfinite(NOE[v])) << "NH NOE must be finite";
                EXPECT_GT(R1[v], 0.0);
                EXPECT_GT(R2[v], 0.0);
                EXPECT_GE(R2[v], R1[v]);  // J(0) dominates R2's extra terms
                ASSERT_EQ(jw[v].size(), 5u);
                for (double jj : jw[v]) {
                    ASSERT_TRUE(std::isfinite(jj)) << "NH J(w) must be finite";
                    EXPECT_GE(jj, 0.0) << "spectral density is non-negative";
                }
            } else {  // CaHa / CO: no 15N params -> NaN
                EXPECT_TRUE(std::isnan(R1[v])) << "non-NH R1 must be NaN";
                EXPECT_TRUE(std::isnan(NOE[v])) << "non-NH NOE must be NaN";
            }
        }
        EXPECT_GT(n_nh, 0u) << "1P9J should have N-H vectors";
    }

    // ── iRED ─────────────────────────────────────────────────────────
    ASSERT_TRUE(f.exist("/trajectory/ired_order_parameters"));
    auto ig = f.getGroup("/trajectory/ired_order_parameters");
    {
        std::vector<double> s2i, evals;
        ig.getDataSet("s2_ired").read(s2i);
        ig.getDataSet("eigenvalues").read(evals);
        EXPECT_GT(s2i.size(), 0u);
        for (double v : s2i) {
            ASSERT_TRUE(std::isfinite(v)) << "non-finite iRED S^2";
            EXPECT_GE(v, -0.2);
            EXPECT_LE(v, 1.05);
        }
        // Eigenvalues are emitted descending.
        for (std::size_t k = 1; k < evals.size(); ++k)
            EXPECT_GE(evals[k - 1] + 1e-9, evals[k]) << "eigenvalues not descending";
    }

    // ── Dihedral autocorrelation ─────────────────────────────────────
    ASSERT_TRUE(f.exist("/trajectory/dihedral_autocorrelation"));
    auto dg = f.getGroup("/trajectory/dihedral_autocorrelation");
    {
        std::vector<std::vector<double>> phi_acf;
        dg.getDataSet("phi_acf").read(phi_acf);
        std::vector<std::uint8_t> phi_defined;
        dg.getDataSet("phi_defined").read(phi_defined);
        ASSERT_EQ(phi_acf.size(), phi_defined.size());
        // Circular ACF C(0) == 1 for every defined residue; NaN where
        // structurally undefined (terminus phi).
        std::size_t n_defined = 0;
        for (std::size_t r = 0; r < phi_acf.size(); ++r) {
            if (phi_defined[r]) {
                ++n_defined;
                EXPECT_NEAR(phi_acf[r].at(0), 1.0, 1e-9)
                    << "phi circular ACF C(0) must be 1 for defined residue " << r;
            }
        }
        EXPECT_GT(n_defined, 0u) << "1P9J should have phi-defined residues";
    }

    nmr::OperationLog::Info("DynamicsObservablesTest::H5RoundTripAndContent",
        "verified reorientational + iRED + dihedral groups on 1P9J fixture");
    std::error_code ec;
    fs::remove(h5_path, ec);
}


// ── KernelDynamics + KernelCoherence on the real PerFrameExtractionSet ──
//
// The instrument TRs read the classical kernel ConformationResults incl.
// ApbsField, so they are exercised on the production config verbatim --
// APBS on, AIMNet2 on, MOPAC skipped -- which already Produces both. This
// is the first in-harness production-shape run; it skips only when the
// AIMNet2 model is unavailable (same gate as the other AIMNet2 tests).

TEST(DynamicsObservables, KernelInstrumentPerFrameExtractionSet) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    const std::string& model_path = nmr::test::TestEnvironment::Aimnet2Model();
    if (model_path.empty() || !fs::exists(model_path))
        GTEST_SKIP() << "AIMNet2 model not available";

    nmr::Session session;
    ASSERT_EQ(session.LoadAimnet2Model(model_path), nmr::kOk) << session.LastError();

    // Production config verbatim (APBS on, AIMNet2 on, MOPAC off). It already
    // attaches KernelDynamics + KernelCoherence. Stride 300 -> ~6 frames.
    auto config = nmr::RunConfiguration::PerFrameExtractionSet();
    config.SetStride(300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& kd = tp.Result<nmr::KernelDynamicsTrajectoryResult>();
    const auto& kc = tp.Result<nmr::KernelCoherenceTrajectoryResult>();

    const std::string h5_path = (fs::temp_directory_path() /
        ("dynobs_kernel_" + std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        kd.WriteH5Group(tp, file);
        kc.WriteH5Group(tp, file);
    }
    HighFive::File f(h5_path, HighFive::File::ReadOnly);

    // ── KernelDynamics: rho(0)=1 for varying channels, PSD >= 0 ──────
    ASSERT_TRUE(f.exist("/trajectory/kernel_dynamics"));
    auto kg = f.getGroup("/trajectory/kernel_dynamics");
    {
        std::vector<std::vector<std::vector<double>>> acf;     // (N, C, L)
        kg.getDataSet("acf").read(acf);
        std::vector<std::vector<std::vector<double>>> spec;    // (N, C, F)
        kg.getDataSet("power_spectrum").read(spec);
        ASSERT_GT(acf.size(), 0u);
        ASSERT_EQ(acf[0].size(), 13u) << "13 kernel channels";
        std::size_t n_active = 0;
        for (const auto& atom : acf)
            for (const auto& ch : atom) {
                const double r0 = ch.at(0);
                if (std::abs(r0) > 1e-9) {  // non-constant channel
                    EXPECT_NEAR(r0, 1.0, 1e-9) << "rho(0) must be 1";
                    ++n_active;
                }
            }
        EXPECT_GT(n_active, 0u) << "expected varying kernel channels on 6 frames";
        // Parzen one-sided PSD is non-negative everywhere.
        for (const auto& atom : spec)
            for (const auto& ch : atom)
                for (double s : ch)
                    EXPECT_GE(s, -1e-12) << "power spectrum must be >= 0";
    }

    // ── KernelCoherence: zero-lag Pearson, diagonal 1 (or NaN constant) ──
    ASSERT_TRUE(f.exist("/trajectory/kernel_coherence"));
    auto cg = f.getGroup("/trajectory/kernel_coherence");
    {
        std::vector<std::vector<std::vector<double>>> cm;      // (N, C, C)
        cg.getDataSet("correlation_matrix").read(cm);
        ASSERT_GT(cm.size(), 0u);
        const std::size_t C = cm[0].size();
        ASSERT_EQ(C, 13u);
        std::size_t n_diag = 0;
        for (const auto& atom : cm)
            for (std::size_t a = 0; a < C; ++a) {
                const double diag = atom[a][a];
                if (std::isfinite(diag)) { EXPECT_NEAR(diag, 1.0, 1e-9); ++n_diag; }
                for (std::size_t b = 0; b < C; ++b) {
                    const double off = atom[a][b];
                    if (std::isfinite(off)) {
                        EXPECT_GE(off, -1.0 - 1e-9);
                        EXPECT_LE(off,  1.0 + 1e-9);
                    }
                }
            }
        EXPECT_GT(n_diag, 0u) << "expected finite (varying) channel diagonals";
    }

    std::error_code ec;
    fs::remove(h5_path, ec);
    nmr::OperationLog::Info("DynamicsObservablesTest::KernelInstrument",
        "verified kernel_dynamics + kernel_coherence on PerFrameExtractionSet");
}
