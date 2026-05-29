#include "DihedralAutocorrelationTrajectoryResult.h"

#include "CalculatorConfig.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "TrajectoryProtein.h"
#include "Types.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace nmr {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Signed dihedral, atan2 form, NaN at degeneracy. Cloned from
// DihedralTimeSeriesTrajectoryResult (PATTERNS.md 17, 10):
//   D = atan2( (n1 x b2hat) . n2, n1 . n2 ),
//   b1=p2-p1, b2=p3-p2, b3=p4-p3, n1=b1xb2, n2=b2xb3.
double Dihedral(const Vec3& p1, const Vec3& p2,
                const Vec3& p3, const Vec3& p4) {
    const Vec3 b1 = p2 - p1;
    const Vec3 b2 = p3 - p2;
    const Vec3 b3 = p4 - p3;
    const double b2n = b2.norm();
    if (b2n < 1e-10) return kNaN;
    const Vec3 n1 = b1.cross(b2);
    const Vec3 n2 = b2.cross(b3);
    if (n1.norm() < 1e-10 || n2.norm() < 1e-10) return kNaN;
    const Vec3 m1 = n1.cross(b2 / b2n);
    return std::atan2(m1.dot(n2), n1.dot(n2));
}

// Push one frame's angle into a circular-ACF accumulator. A finite value
// advances the hold; a non-finite value (rare degenerate geometry) repeats
// the last finite value so the lag-to-time mapping is preserved; an angle
// never yet finite (structurally undefined) is skipped entirely.
void PushAngle(CircularAcfAccumulator& acc, double& last,
               std::uint8_t& has_last, double val) {
    if (std::isfinite(val)) {
        acc.Push(val);
        last = val;
        has_last = 1;
    } else if (has_last) {
        acc.Push(last);
    }
}

// 1/e decorrelation time: dt * the (interpolated) lag where C(k) first
// falls to e^-1. NaN if the curve is undefined (fewer than two samples);
// the full window (a lower bound) if C never decays that far. Units ps.
double OneOverETime(const std::vector<double>& c, double dt) {
    const std::size_t L = c.size();
    if (!(dt > 0.0) || L < 2) return kNaN;
    if (!std::isfinite(c[0]) || !std::isfinite(c[1])) return kNaN;
    const double thresh = std::exp(-1.0);
    std::size_t last_finite_k = 1;
    for (std::size_t k = 1; k < L; ++k) {
        if (!std::isfinite(c[k])) break;
        last_finite_k = k;
        if (c[k] <= thresh) {
            const double c0 = c[k - 1], c1 = c[k];
            const double frac = (c0 > c1) ? (c0 - thresh) / (c0 - c1) : 0.0;
            return dt * (static_cast<double>(k - 1) + frac);
        }
    }
    // No crossing within the finite window: lower bound = the last finite
    // lag (NOT L-1, which over-reports when the curve is only defined to
    // T-1 because the run is shorter than the lag count).
    return dt * static_cast<double>(last_finite_k);
}

}  // namespace


std::unique_ptr<DihedralAutocorrelationTrajectoryResult>
DihedralAutocorrelationTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<DihedralAutocorrelationTrajectoryResult>();
    const Protein& protein = tp.ProteinRef();
    const std::size_t R = protein.ResidueCount();
    const std::size_t N = tp.AtomCount();
    r->n_residues_ = R;
    r->n_atoms_ = N;

    // Clamp to >= 1: a misconfigured dynamics_n_lags <= 0 must not size the
    // accumulators to zero (Finalize would then read an empty vector).
    // codex review 2026-05-29.
    const double n_lags_raw = CalculatorConfig::Get("dynamics_n_lags");
    const std::size_t n_lags =
        (n_lags_raw >= 1.0) ? static_cast<std::size_t>(n_lags_raw) : 1;
    r->n_lags_ = n_lags;
    r->phi_acc_.assign(R, CircularAcfAccumulator(n_lags));
    r->psi_acc_.assign(R, CircularAcfAccumulator(n_lags));
    r->chi_acc_.assign(R * 4, CircularAcfAccumulator(n_lags));

    r->phi_defined_.assign(R, 0);
    r->psi_defined_.assign(R, 0);
    r->chi_defined_.assign(R * 4, 0);
    r->phi_last_.assign(R, 0.0);
    r->psi_last_.assign(R, 0.0);
    r->chi_last_.assign(R * 4, 0.0);
    r->phi_has_last_.assign(R, 0);
    r->psi_has_last_.assign(R, 0);
    r->chi_has_last_.assign(R * 4, 0);

    for (std::size_t ri = 0; ri < R; ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        const auto prev = protein.BackbonePredecessor(ri);
        const auto next = protein.BackboneSuccessor(ri);
        if (prev && res.N != Residue::NONE && res.CA != Residue::NONE &&
            res.C != Residue::NONE) {
            r->phi_defined_[ri] = 1;
        }
        if (next && res.N != Residue::NONE && res.CA != Residue::NONE &&
            res.C != Residue::NONE) {
            r->psi_defined_[ri] = 1;
        }
        for (int k = 0; k < 4; ++k) {
            if (res.chi[k].Valid()) r->chi_defined_[ri * 4 + k] = 1;
        }
    }

    r->residue_index_per_atom_.assign(N, -1);
    for (std::size_t ai = 0; ai < N; ++ai) {
        r->residue_index_per_atom_[ai] =
            static_cast<std::int32_t>(protein.AtomAt(ai).residue_index);
    }
    return r;
}


// ── Compute ──────────────────────────────────────────────────────

void DihedralAutocorrelationTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj; (void)frame_idx;
    const Protein& protein = conf.ProteinRef();
    const std::size_t R = protein.ResidueCount();

    for (std::size_t ri = 0; ri < R; ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        const auto prev = protein.BackbonePredecessor(ri);
        const auto next = protein.BackboneSuccessor(ri);

        double phi = kNaN;
        if (prev && res.CA != Residue::NONE && res.C != Residue::NONE) {
            const Residue& rp = protein.ResidueAt(*prev);
            phi = Dihedral(conf.PositionAt(rp.C), conf.PositionAt(res.N),
                           conf.PositionAt(res.CA), conf.PositionAt(res.C));
        }
        PushAngle(phi_acc_[ri], phi_last_[ri], phi_has_last_[ri], phi);

        double psi = kNaN;
        if (next && res.N != Residue::NONE && res.CA != Residue::NONE &&
            res.C != Residue::NONE) {
            const Residue& rn = protein.ResidueAt(*next);
            psi = Dihedral(conf.PositionAt(res.N), conf.PositionAt(res.CA),
                           conf.PositionAt(res.C), conf.PositionAt(rn.N));
        }
        PushAngle(psi_acc_[ri], psi_last_[ri], psi_has_last_[ri], psi);

        for (int k = 0; k < 4; ++k) {
            double chi = kNaN;
            if (res.chi[k].Valid()) {
                chi = Dihedral(conf.PositionAt(res.chi[k].a[0]),
                               conf.PositionAt(res.chi[k].a[1]),
                               conf.PositionAt(res.chi[k].a[2]),
                               conf.PositionAt(res.chi[k].a[3]));
            }
            const std::size_t c = ri * 4 + k;
            PushAngle(chi_acc_[c], chi_last_[c], chi_has_last_[c], chi);
        }
    }

    frame_times_.push_back(time_ps);
    ++n_frames_;
}


// ── Finalize ──────────────────────────────────────────────────────

void DihedralAutocorrelationTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                       Trajectory& traj) {
    (void)tp; (void)traj;
    if (finalized_) return;
    const std::size_t R = n_residues_;
    const std::size_t L = n_lags_;

    sample_interval_ps_ = 0.0;
    if (frame_times_.size() >= 2) {
        std::vector<double> d;
        d.reserve(frame_times_.size() - 1);
        for (std::size_t t = 1; t < frame_times_.size(); ++t) {
            d.push_back(frame_times_[t] - frame_times_[t - 1]);
        }
        std::sort(d.begin(), d.end());
        sample_interval_ps_ = d[d.size() / 2];
    }

    phi_acf_.assign(R * L, kNaN);
    psi_acf_.assign(R * L, kNaN);
    chi_acf_.assign(R * 4 * L, kNaN);
    phi_corr_time_.assign(R, kNaN);
    psi_corr_time_.assign(R, kNaN);
    chi_corr_time_.assign(R * 4, kNaN);

    auto finalize_one = [&](CircularAcfAccumulator& acc, double* acf_row,
                            double& corr_time) {
        const std::vector<double> c = acc.Finalize();
        for (std::size_t k = 0; k < L; ++k) acf_row[k] = c[k];
        corr_time = OneOverETime(c, sample_interval_ps_);
    };

    for (std::size_t ri = 0; ri < R; ++ri) {
        finalize_one(phi_acc_[ri], &phi_acf_[ri * L], phi_corr_time_[ri]);
        finalize_one(psi_acc_[ri], &psi_acf_[ri * L], psi_corr_time_[ri]);
        for (int k = 0; k < 4; ++k) {
            const std::size_t c = ri * 4 + k;
            finalize_one(chi_acc_[c], &chi_acf_[c * L], chi_corr_time_[c]);
        }
    }

    std::vector<CircularAcfAccumulator>().swap(phi_acc_);
    std::vector<CircularAcfAccumulator>().swap(psi_acc_);
    std::vector<CircularAcfAccumulator>().swap(chi_acc_);
    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "DihedralAutocorrelationTrajectoryResult::Finalize",
        "circular ACF (phi/psi/chi) across " + std::to_string(R) +
        " residues, " + std::to_string(n_frames_) + " frames; sample "
        "interval ~ " + std::to_string(sample_interval_ps_) + " ps");
}


// ── WriteH5Group ──────────────────────────────────────────────────

void DihedralAutocorrelationTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    (void)tp;
    if (!finalized_) {
        OperationLog::Warn(
            "DihedralAutocorrelationTrajectoryResult::WriteH5Group",
            "Finalize not called; nothing to write");
        return;
    }
    const std::size_t R = n_residues_;
    const std::size_t L = n_lags_;

    auto grp = file.createGroup("/trajectory/dihedral_autocorrelation");
    grp.createAttribute("result_name",        Name());
    grp.createAttribute("n_residues",         R);
    grp.createAttribute("n_atoms",            n_atoms_);
    grp.createAttribute("n_lags",             L);
    grp.createAttribute("n_frames",           n_frames_);
    grp.createAttribute("finalized",          finalized_);
    grp.createAttribute("sample_interval_ps", sample_interval_ps_);
    grp.createAttribute("estimator",          std::string("circular_unbiased"));
    grp.createAttribute("acf_definition", std::string(
        "C(k) = <cos(theta(t+k) - theta(t))>, C(0)=1, plateau = "
        "<cos>^2 + <sin>^2 (circular order). Unbiased (T-k pairs per lag)."));
    grp.createAttribute("correlation_time_definition", std::string(
        "1/e decorrelation time: dt * interpolated lag where C(k) first "
        "<= exp(-1); the full window (a lower bound) when C never decays "
        "that far; NaN where the angle is structurally undefined. Units ps."));
    grp.createAttribute("angle_convention", std::string(
        "IUPAC signed dihedral atan2(y,x). phi = C(i-1)-N-CA-C; "
        "psi = N-CA-C-N(i+1); chi_k from Residue.chi[k]. Backbone "
        "adjacency via Protein::BackbonePredecessor/Successor (bond graph)."));
    grp.createAttribute("deferred", std::string(
        "torsional power spectrum and rotamer-state survival ACF (jump "
        "time) are planned follow-ups; rotamer transition counts are in "
        "/trajectory/dihedral_bin_transition/."));

    auto emit_2d = [&](const std::string& name,
                       const std::vector<double>& flat) {
        std::vector<std::size_t> dims = {R, L};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<double>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("dimensionless"));
        ds.createAttribute("layout", std::string("(residue, lag)"));
    };
    emit_2d("phi_acf", phi_acf_);
    emit_2d("psi_acf", psi_acf_);
    {
        std::vector<std::size_t> dims = {R, std::size_t(4), L};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<double>("chi_acf", space);
        ds.write_raw(chi_acf_.data());
        ds.createAttribute("units", std::string("dimensionless"));
        ds.createAttribute("layout", std::string("(residue, chi_index, lag)"));
    }

    grp.createDataSet("phi_corr_time_ps", phi_corr_time_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("psi_corr_time_ps", psi_corr_time_)
       .createAttribute("units", std::string("ps"));
    {
        std::vector<std::size_t> dims = {R, std::size_t(4)};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<double>("chi_corr_time_ps", space);
        ds.write_raw(chi_corr_time_.data());
        ds.createAttribute("units", std::string("ps"));
        ds.createAttribute("layout", std::string("(residue, chi_index)"));
    }

    grp.createDataSet("phi_defined", phi_defined_);
    grp.createDataSet("psi_defined", psi_defined_);
    {
        std::vector<std::size_t> dims = {R, std::size_t(4)};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<std::uint8_t>("chi_defined", space);
        ds.write_raw(chi_defined_.data());
        ds.createAttribute("layout", std::string("(residue, chi_index)"));
    }
    grp.createDataSet("residue_index_per_atom", residue_index_per_atom_)
       .createAttribute("units", std::string("residue_index"));

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
