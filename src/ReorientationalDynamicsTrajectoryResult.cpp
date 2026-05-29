#include "ReorientationalDynamicsTrajectoryResult.h"

#include "CalculatorConfig.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <limits>

namespace nmr {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Optimal rotation aligning `current` onto `reference` (both 3xM point
// clouds, one point per entry), reflection-corrected. Cloned from
// RmsdTrackingTrajectoryResult::KabschRmsd (PATTERNS.md 17), returning the
// rotation R instead of the RMSD: R rotates the centred current frame onto
// the reference, so a bond vector in the current frame becomes R * v in the
// reference (body) frame. M < 3 is rotationally underdetermined -> identity.
Mat3 KabschRotation(const std::vector<Vec3>& current,
                    const std::vector<Vec3>& reference) {
    const std::size_t M = current.size();
    if (M < 3 || reference.size() != M) return Mat3::Identity();

    Vec3 p_cen = Vec3::Zero();
    Vec3 q_cen = Vec3::Zero();
    for (std::size_t i = 0; i < M; ++i) { p_cen += current[i]; q_cen += reference[i]; }
    p_cen /= static_cast<double>(M);
    q_cen /= static_cast<double>(M);

    Eigen::Matrix<double, 3, Eigen::Dynamic> P(3, M), Q(3, M);
    for (std::size_t i = 0; i < M; ++i) {
        P.col(i) = current[i]   - p_cen;
        Q.col(i) = reference[i] - q_cen;
    }
    const Mat3 H = P * Q.transpose();
    Eigen::JacobiSVD<Mat3> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Mat3& U = svd.matrixU();
    const Mat3& V = svd.matrixV();
    // sign(det(V U^T)) only -- passing the raw det would inject O(eps)
    // anisotropic scaling and break orthogonality (RmsdTracking note).
    const double det = (V * U.transpose()).determinant();
    Eigen::DiagonalMatrix<double, 3> D(1.0, 1.0, (det < 0.0) ? -1.0 : 1.0);
    return V * D * U.transpose();
}

}  // namespace


std::unique_ptr<ReorientationalDynamicsTrajectoryResult>
ReorientationalDynamicsTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<ReorientationalDynamicsTrajectoryResult>();
    const Protein& protein = tp.ProteinRef();
    const std::size_t R = protein.ResidueCount();

    auto add = [&](Kind k, std::size_t tail, std::size_t head,
                   std::size_t owning, std::size_t ri) {
        r->kind_.push_back(static_cast<std::uint8_t>(k));
        r->tail_atom_.push_back(tail);
        r->head_atom_.push_back(head);
        r->owning_atom_.push_back(owning);
        r->residue_index_.push_back(static_cast<std::int32_t>(ri));
    };

    for (std::size_t ri = 0; ri < R; ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.N != Residue::NONE && res.H != Residue::NONE)
            add(Kind::NH, res.N, res.H, res.H, ri);          // N->H, key H
        if (res.CA != Residue::NONE && res.HA != Residue::NONE)
            add(Kind::CaHa, res.CA, res.HA, res.HA, ri);     // CA->HA, key HA
        if (res.C != Residue::NONE && res.O != Residue::NONE)
            add(Kind::CO, res.C, res.O, res.O, ri);          // C->O, key O
        // Backbone-heavy alignment set (N/CA/C/O), per RmsdTracking.
        if (res.N  != Residue::NONE) r->align_atoms_.push_back(res.N);
        if (res.CA != Residue::NONE) r->align_atoms_.push_back(res.CA);
        if (res.C  != Residue::NONE) r->align_atoms_.push_back(res.C);
        if (res.O  != Residue::NONE) r->align_atoms_.push_back(res.O);
    }

    // Clamp to >= 1: a misconfigured dynamics_n_lags <= 0 must not size the
    // TCF accumulators to zero (Finalize would then read an empty vector).
    // codex review 2026-05-29.
    const double n_lags_raw = CalculatorConfig::Get("dynamics_n_lags");
    const std::size_t n_lags =
        (n_lags_raw >= 1.0) ? static_cast<std::size_t>(n_lags_raw) : 1;
    r->n_lags_ = n_lags;
    r->n_vectors_ = r->kind_.size();
    r->body_acc_.assign(r->n_vectors_, LegendreTcfAccumulator(n_lags));
    r->lab_acc_.assign(r->n_vectors_, LegendreTcfAccumulator(n_lags));
    r->order_sum_.assign(r->n_vectors_ * 6, 0.0);

    OperationLog::Info("ReorientationalDynamicsTrajectoryResult::Create",
        "backbone vectors: " + std::to_string(r->n_vectors_) + " (N-H/CA-HA/C=O), "
        "alignment set: " + std::to_string(r->align_atoms_.size()) + " heavy atoms");
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Per frame: the Kabsch rotation onto frame 0 (computed once from the
// alignment set), applied to each bond unit vector to reach the body
// frame. Push body and lab vectors to their TCF accumulators and add the
// body-frame order-tensor products.

void ReorientationalDynamicsTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj; (void)frame_idx;
    const std::size_t M = align_atoms_.size();

    std::vector<Vec3> current;
    current.reserve(M);
    for (std::size_t ai : align_atoms_) current.push_back(conf.PositionAt(ai));

    Mat3 R = Mat3::Identity();
    if (!reference_captured_) {
        reference_positions_ = current;        // frame 0 is its own reference
        reference_captured_ = true;
    } else {
        R = KabschRotation(current, reference_positions_);
    }

    for (std::size_t v = 0; v < n_vectors_; ++v) {
        const Vec3 d = conf.PositionAt(head_atom_[v]) - conf.PositionAt(tail_atom_[v]);
        const double nrm = d.norm();
        // Backbone covalent bonds are never zero-length; this guard is
        // defensive. Skip a degenerate frame rather than push a zero
        // "unit" vector that would poison the TCF and the order tensor.
        if (nrm < 1e-9) continue;
        const Vec3 lab = d / nrm;
        const Vec3 body = R * lab;

        lab_acc_[v].Push(lab.x(), lab.y(), lab.z());
        body_acc_[v].Push(body.x(), body.y(), body.z());

        double* os = &order_sum_[v * 6];
        os[0] += body.x() * body.x();
        os[1] += body.y() * body.y();
        os[2] += body.z() * body.z();
        os[3] += body.x() * body.y();
        os[4] += body.x() * body.z();
        os[5] += body.y() * body.z();
    }

    frame_times_.push_back(time_ps);
    ++n_frames_;
}


// ── Finalize ──────────────────────────────────────────────────────

void ReorientationalDynamicsTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                       Trajectory& traj) {
    (void)tp; (void)traj;
    if (finalized_) return;
    const std::size_t V = n_vectors_;
    const std::size_t L = n_lags_;

    sample_interval_ps_ = 0.0;
    if (frame_times_.size() >= 2) {
        std::vector<double> dd;
        dd.reserve(frame_times_.size() - 1);
        for (std::size_t t = 1; t < frame_times_.size(); ++t)
            dd.push_back(frame_times_[t] - frame_times_[t - 1]);
        std::sort(dd.begin(), dd.end());
        sample_interval_ps_ = dd[dd.size() / 2];
    }
    const double dt = sample_interval_ps_;

    acf_internal_.assign(V * L, kNaN);
    acf_lab_.assign(V * L, kNaN);
    s2_.assign(V, kNaN);
    tau_e_.assign(V, kNaN);
    orient_tensor_.assign(V * 9, kNaN);

    const double T = static_cast<double>(n_frames_);
    for (std::size_t v = 0; v < V; ++v) {
        const std::vector<double> ci = body_acc_[v].Finalize();
        const std::vector<double> cl = lab_acc_[v].Finalize();
        for (std::size_t k = 0; k < L; ++k) {
            acf_internal_[v * L + k] = ci[k];
            acf_lab_[v * L + k] = cl[k];
        }
        if (n_frames_ < 2) continue;

        const double* os = &order_sum_[v * 6];
        const double mxx = os[0] / T, myy = os[1] / T, mzz = os[2] / T;
        const double mxy = os[3] / T, mxz = os[4] / T, myz = os[5] / T;
        // S^2 = (3/2)[ <x^2>^2 + <y^2>^2 + <z^2>^2
        //            + 2(<xy>^2 + <xz>^2 + <yz>^2) ] - 1/2  (Henry-Szabo 1985)
        const double s2 = 1.5 * (mxx * mxx + myy * myy + mzz * mzz
                                 + 2.0 * (mxy * mxy + mxz * mxz + myz * myz))
                          - 0.5;
        s2_[v] = s2;

        double* ot = &orient_tensor_[v * 9];
        ot[0] = mxx; ot[1] = mxy; ot[2] = mxz;
        ot[3] = mxy; ot[4] = myy; ot[5] = myz;
        ot[6] = mxz; ot[7] = myz; ot[8] = mzz;

        // tau_e area method: integral [C_I(t) - S^2] / (1 - S^2) dt, summed
        // from lag 0 to the first lag where C_I has decayed to S^2.
        if (dt > 0.0 && std::isfinite(s2)) {
            if ((1.0 - s2) > 1e-6) {
                double area = 0.0;
                for (std::size_t k = 0; k < L; ++k) {
                    if (!std::isfinite(ci[k]) || ci[k] <= s2) break;
                    area += (ci[k] - s2) / (1.0 - s2);
                }
                tau_e_[v] = area * dt;
            } else {
                // Effectively rigid (S^2 ~ 1): no internal decay, tau_e -> 0.
                // The (1-S^2) term in J vanishes regardless, but tau_e must be
                // finite (not the skipped NaN) so the relaxation layer runs
                // for rigid vectors -- exactly where relaxation is best
                // defined (codex review 2026-05-29).
                tau_e_[v] = 0.0;
            }
        }
    }

    // Global tau_m: the trajectory-averaged lab-frame N-H second-rank TCF
    // C_O_est(k), integrated to its first zero. This is the mean N-H
    // reorientational correlation time (overall tumbling plus any residual
    // internal motion) -- honest on long runs, unreliable on short ones, so
    // it ships with the convergence ratio and flag, never trusted silently.
    tau_m_ps_ = kNaN;
    traj_len_over_tau_m_ = kNaN;
    tau_m_converged_ = false;
    if (dt > 0.0 && n_frames_ >= 2) {
        std::vector<double> co(L, 0.0);
        std::vector<std::size_t> cnt(L, 0);
        for (std::size_t v = 0; v < V; ++v) {
            if (kind_[v] != static_cast<std::uint8_t>(Kind::NH)) continue;
            for (std::size_t k = 0; k < L; ++k) {
                const double c = acf_lab_[v * L + k];
                if (std::isfinite(c)) { co[k] += c; ++cnt[k]; }
            }
        }
        if (cnt[0] > 0) {
            double area = 0.0;
            for (std::size_t k = 0; k < L; ++k) {
                if (cnt[k] == 0) break;
                const double mean = co[k] / static_cast<double>(cnt[k]);
                if (mean <= 0.0) break;
                area += mean;
            }
            tau_m_ps_ = area * dt;
            const double traj_len = (T - 1.0) * dt;  // T samples span T-1 intervals
            if (tau_m_ps_ > 0.0) {
                traj_len_over_tau_m_ = traj_len / tau_m_ps_;
                tau_m_converged_ = traj_len_over_tau_m_ >= 5.0;
            }
        }
    }

    // ── 15N relaxation layer (N-H vectors) ───────────────────────────
    // Lipari-Szabo spectral density + the 15N dipolar + CSA rate equations.
    // Rides on the SAME global tau_m and per-vector S^2/tau_e computed above,
    // so it inherits tau_m_converged: on a run too short for tau_m (1P9J at
    // 15 ns) the rates are computed and emitted but flagged unreliable.
    // Constants + citations: PhysicalConstants.h; field: CalculatorConfig.
    relax_freqs_.assign(5, kNaN);
    spectral_density_j_.assign(V * 5, kNaN);
    r1_.assign(V, kNaN);
    r2_.assign(V, kNaN);
    noe_.assign(V, kNaN);

    relaxation_field_tesla_ = CalculatorConfig::Get("relaxation_field_tesla");
    const double B0 = relaxation_field_tesla_;
    // Larmor angular frequencies. J(omega) is even, so only |omega| enters J
    // -- BUT the 15N combination frequencies must use the SIGNED Larmor
    // frequencies (Kay-Torchia-Bax 1989). gamma_N < 0, so the KTB term
    // "J(omega_H - omega_N)" is the HIGH combination |omega_H| + |omega_N|,
    // and "J(omega_H + omega_N)" is the LOW one |omega_H| - |omega_N|. Using
    // magnitudes for the combinations silently swaps the 1x/6x spectral
    // weights and inflates the steady-state NOE toward 1 (~0.89 instead of
    // the correct ~0.82-0.84). Caught by codex review 2026-05-29.
    const double wH_signed = GAMMA_H   * B0;          // > 0
    const double wN_signed = GAMMA_N15 * B0;          // < 0 for 15N
    const double wH    = std::abs(wH_signed);         // |omega_H|
    const double wN    = std::abs(wN_signed);         // |omega_N|
    const double w_zq  = std::abs(wH_signed - wN_signed);  // |omega_H - omega_N| (high)
    const double w_dq  = std::abs(wH_signed + wN_signed);  // |omega_H + omega_N| (low)
    proton_larmor_mhz_ = wH / (2.0 * PI) / 1.0e6;
    // relax_freqs_ slot order matches spectral_density_j: [0, omega_N,
    // |omega_H-omega_N|, omega_H, |omega_H+omega_N|]. For 15N the third slot
    // (high) exceeds the fifth (low) -- the signed-combination ordering.
    const double freqs5[5] = {0.0, wN, w_zq, wH, w_dq};
    for (int i = 0; i < 5; ++i) relax_freqs_[i] = freqs5[i];

    // Dipolar coupling constant d = mu0 hbar gamma_H gamma_N / (4 pi r^3)
    // [rad/s]; d enters the rates squared, so gamma_N's sign is immaterial
    // here. Units check: d^2 [s^-2] * J [s] -> rate [s^-1].
    const double r_nh_m = NH_DIPOLAR_BOND_LENGTH_A * ANGSTROMS_TO_METRES;
    const double d = (VACUUM_PERMEABILITY * REDUCED_PLANCK * GAMMA_H * GAMMA_N15)
                     / (4.0 * PI * r_nh_m * r_nh_m * r_nh_m);
    const double d2 = d * d;
    // CSA constant c = wN * Dsigma / sqrt(3) [rad/s]; squared in the rates.
    const double dsigma = N15_CSA_PPM * 1.0e-6;   // ppm -> dimensionless
    const double c = wN * dsigma / std::sqrt(3.0);
    const double c2 = c * c;

    const double tau_m_s = (std::isfinite(tau_m_ps_) && tau_m_ps_ > 0.0)
                           ? tau_m_ps_ * 1.0e-12 : kNaN;
    if (std::isfinite(tau_m_s)) {
        for (std::size_t v = 0; v < V; ++v) {
            if (kind_[v] != static_cast<std::uint8_t>(Kind::NH)) continue;
            const double s2 = s2_[v];
            const double tau_e_s = tau_e_[v] * 1.0e-12;   // ps -> s
            if (!std::isfinite(s2) || !std::isfinite(tau_e_s)) continue;
            // Effective internal correlation time: 1/tau = 1/tau_m + 1/tau_e.
            // tau_e = 0 (rigid internal) -> tau = 0; the (1-S^2) term vanishes.
            const double tau_s = (tau_e_s > 0.0)
                ? (tau_m_s * tau_e_s) / (tau_m_s + tau_e_s) : 0.0;
            // J(w) = (2/5)[ S^2 tau_m/(1+(w tau_m)^2)
            //             + (1-S^2) tau/(1+(w tau)^2) ]   (Lipari-Szabo 1982)
            auto J = [&](double w) {
                const double a = s2 * tau_m_s
                                 / (1.0 + (w * tau_m_s) * (w * tau_m_s));
                const double b = (1.0 - s2) * tau_s
                                 / (1.0 + (w * tau_s) * (w * tau_s));
                return 0.4 * (a + b);
            };
            // KTB terms in signed notation: jHmN = J(omega_H - omega_N) is
            // the HIGH combination (w_zq) and jHpN = J(omega_H + omega_N) is
            // the LOW one (w_dq) for 15N (gamma_N < 0). J is even, so we feed
            // the magnitudes computed above.
            const double j0   = J(0.0);
            const double jN   = J(wN);
            const double jHmN = J(w_zq);   // J(omega_H - omega_N), high
            const double jH   = J(wH);
            const double jHpN = J(w_dq);   // J(omega_H + omega_N), low
            double* jr = &spectral_density_j_[v * 5];
            jr[0] = j0; jr[1] = jN; jr[2] = jHmN; jr[3] = jH; jr[4] = jHpN;

            // 15N R1/R2 and steady-state {1H}-NOE (Kay, Torchia & Bax 1989).
            const double R1 = (d2 / 4.0) * (jHmN + 3.0 * jN + 6.0 * jHpN)
                              + c2 * jN;
            const double R2 = (d2 / 8.0) * (4.0 * j0 + jHmN + 3.0 * jN
                                            + 6.0 * jH + 6.0 * jHpN)
                              + (c2 / 6.0) * (4.0 * j0 + 3.0 * jN);
            // Cross-relaxation sigma; NOE = 1 + (gamma_H/gamma_N)(sigma/R1).
            // gamma_H/gamma_N < 0 (gamma_N < 0) is exactly what drives the
            // 15N NOE below 1 (negative for fast internal motion).
            const double sigma = (d2 / 4.0) * (6.0 * jHpN - jHmN);
            r1_[v] = R1;
            r2_[v] = R2;
            noe_[v] = (R1 > 0.0)
                ? 1.0 + (GAMMA_H / GAMMA_N15) * (sigma / R1) : kNaN;
        }
    }

    std::vector<LegendreTcfAccumulator>().swap(body_acc_);
    std::vector<LegendreTcfAccumulator>().swap(lab_acc_);
    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "ReorientationalDynamicsTrajectoryResult::Finalize",
        std::to_string(V) + " backbone vectors, " + std::to_string(n_frames_) +
        " frames; tau_m ~ " + std::to_string(tau_m_ps_) + " ps, length/tau_m = " +
        std::to_string(traj_len_over_tau_m_) +
        (tau_m_converged_ ? " (converged)" : " (NOT converged)"));
}


// ── WriteH5Group ──────────────────────────────────────────────────

void ReorientationalDynamicsTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    (void)tp;
    if (!finalized_) {
        OperationLog::Warn(
            "ReorientationalDynamicsTrajectoryResult::WriteH5Group",
            "Finalize not called; nothing to write");
        return;
    }
    if (n_vectors_ == 0) {
        OperationLog::Warn(
            "ReorientationalDynamicsTrajectoryResult::WriteH5Group",
            "no backbone vectors; group omitted");
        return;
    }
    const std::size_t V = n_vectors_;
    const std::size_t L = n_lags_;

    auto grp = file.createGroup("/trajectory/reorientational_dynamics");
    grp.createAttribute("result_name",        Name());
    grp.createAttribute("n_vectors",          V);
    grp.createAttribute("n_lags",             L);
    grp.createAttribute("n_frames",           n_frames_);
    grp.createAttribute("finalized",          finalized_);
    grp.createAttribute("sample_interval_ps", sample_interval_ps_);
    grp.createAttribute("vector_kind_legend", std::string("1=NH, 2=CaHa, 3=CO"));
    grp.createAttribute("superposition",      std::string("kabsch_svd backbone_NCACO"));
    grp.createAttribute("reference",          std::string("trajectory_frame_0"));
    grp.createAttribute("tcf_estimator",      std::string("legendre_P2_unbiased"));
    grp.createAttribute("s2_estimator",       std::string("henry_szabo_order_tensor_body_frame"));
    grp.createAttribute("tau_e_method",       std::string("area_to_first_S2_crossing"));
    grp.createAttribute("tau_m_ps",           tau_m_ps_);
    grp.createAttribute("tau_m_provenance",   std::string(
        "area of the trajectory-averaged lab-frame N-H second-rank TCF; "
        "overall tumbling plus residual internal motion, not a pure "
        "tumbling separation"));
    grp.createAttribute("trajectory_length_over_tau_m", traj_len_over_tau_m_);
    grp.createAttribute("tau_m_converged",    tau_m_converged_);
    grp.createAttribute("tau_m_converged_criterion",
                        std::string("trajectory_length >= 5 * tau_m"));
    grp.createAttribute("deferred", std::string(
        "sidechain X-H, methyl symmetry axis (with the 1/9 fast-rotation "
        "factor), aromatic C-H; and dipole-CSA cross-correlated rates "
        "(eta_xy/eta_z)."));

    // Relaxation-layer provenance (15N-1H, NH rows only).
    grp.createAttribute("relaxation_field_tesla",        relaxation_field_tesla_);
    grp.createAttribute("relaxation_proton_larmor_MHz",  proton_larmor_mhz_);
    grp.createAttribute("relaxation_nh_bond_length_A",
                        NH_DIPOLAR_BOND_LENGTH_A);
    grp.createAttribute("relaxation_n15_csa_ppm",        N15_CSA_PPM);
    grp.createAttribute("relaxation_spectral_density_model",
                        std::string("lipari_szabo_1982"));
    grp.createAttribute("relaxation_equations", std::string(
        "Kay-Torchia-Bax 1989 15N dipolar + CSA: "
        "R1=(d^2/4)[J(wH-wN)+3J(wN)+6J(wH+wN)]+c^2 J(wN); "
        "R2=(d^2/8)[4J0+J(wH-wN)+3J(wN)+6J(wH)+6J(wH+wN)]+(c^2/6)[4J0+3J(wN)]; "
        "NOE=1+(gamma_H/gamma_N)(d^2/4)[6J(wH+wN)-J(wH-wN)]/R1. "
        "d=mu0 hbar gamma_H gamma_N/(4 pi r_NH^3); c=wN Dsigma/sqrt(3). "
        "wH-wN and wH+wN are SIGNED Larmor combinations; for 15N (gamma_N<0) "
        "|wH-wN|=wH+|wN| (high) and |wH+wN|=wH-|wN| (low)."));
    grp.createAttribute("relaxation_reliability", std::string(
        "R1/R2/NOE inherit tau_m_converged: when false (run shorter than "
        "~5 tau_m) the rates are computed but NOT reliable. The S(f)/J(omega) "
        "and S^2 quantities are field-independent and do not carry this caveat."));

    auto emit_v_l = [&](const std::string& name,
                        const std::vector<double>& flat) -> HighFive::DataSet {
        std::vector<std::size_t> dims = {V, L};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<double>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("layout", std::string("(vector, lag)"));
        return ds;
    };
    {
        auto ds = emit_v_l("bond_vector_autocorrelation", acf_internal_);
        ds.createAttribute("frame", std::string("body (tumbling removed)"));
        ds.createAttribute("note", std::string(
            "internal second-rank TCF C_I(k)=<P2(u(t).u(t+k))>; C_I(0)=1, "
            "plateau ~ S^2; mean_subtracted=false"));
    }
    {
        auto ds = emit_v_l("bond_vector_autocorrelation_lab", acf_lab_);
        ds.createAttribute("frame", std::string("lab (carries tumbling)"));
    }

    grp.createDataSet("order_parameter_S2", s2_)
       .createAttribute("units", std::string("dimensionless"));
    grp.createDataSet("lipari_szabo_tau_e", tau_e_)
       .createAttribute("units", std::string("ps"));
    {
        std::vector<std::size_t> dims = {V, std::size_t(3), std::size_t(3)};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<double>("bond_orientation_tensor", space);
        ds.write_raw(orient_tensor_.data());
        ds.createAttribute("layout", std::string("(vector, 3, 3) body-frame <u (x) u>"));
    }

    // 15N relaxation observables (NH rows finite; CaHa/CO rows NaN).
    {
        std::vector<std::size_t> dims = {V, std::size_t(5)};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<double>("spectral_density_j", space);
        ds.write_raw(spectral_density_j_.data());
        ds.createAttribute("units", std::string("s"));
        ds.createAttribute("layout", std::string(
            "(vector, 5) J at the KTB combination frequencies "
            "[0, omega_N, omega_H-omega_N, omega_H, omega_H+omega_N] (SIGNED; "
            "for 15N |omega_H-omega_N| > |omega_H+omega_N|). The magnitude of "
            "each slot is in relaxation_larmor_freqs_rad_per_s."));
    }
    grp.createDataSet("relaxation_R1", r1_)
       .createAttribute("units", std::string("s^-1"));
    grp.createDataSet("relaxation_R2", r2_)
       .createAttribute("units", std::string("s^-1"));
    grp.createDataSet("relaxation_NOE", noe_)
       .createAttribute("units", std::string("dimensionless"));
    grp.createDataSet("relaxation_larmor_freqs_rad_per_s", relax_freqs_)
       .createAttribute("layout", std::string(
           "|angular frequency| rad/s for each spectral_density_j slot: "
           "[0, omega_N, |omega_H-omega_N|, omega_H, |omega_H+omega_N|]"));

    grp.createDataSet("vector_kind", kind_)
       .createAttribute("legend", std::string("1=NH, 2=CaHa, 3=CO"));
    auto to_i32 = [](const std::vector<std::size_t>& v) {
        std::vector<std::int32_t> out(v.size());
        for (std::size_t i = 0; i < v.size(); ++i) out[i] = static_cast<std::int32_t>(v[i]);
        return out;
    };
    grp.createDataSet("owning_atom", to_i32(owning_atom_));
    grp.createDataSet("tail_atom",   to_i32(tail_atom_));
    grp.createDataSet("head_atom",   to_i32(head_atom_));
    grp.createDataSet("residue_index", residue_index_);

    std::vector<std::uint64_t> lag_frames(L);
    std::vector<double> lag_times(L);
    for (std::size_t k = 0; k < L; ++k) {
        lag_frames[k] = k;
        lag_times[k] = static_cast<double>(k) * sample_interval_ps_;
    }
    grp.createDataSet("lag_frames", lag_frames);
    grp.createDataSet("lag_times_ps", lag_times);
}

}  // namespace nmr
