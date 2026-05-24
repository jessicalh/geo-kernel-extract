#include "BsT0AutocorrelationTrajectoryResult.h"
#include "BiotSavartResult.h"
#include "DenseBuffer.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <algorithm>
#include <cmath>
#include <typeinfo>

namespace nmr {


std::vector<std::type_index>
BsT0AutocorrelationTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(BiotSavartResult)) };
}


std::unique_ptr<BsT0AutocorrelationTrajectoryResult>
BsT0AutocorrelationTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<BsT0AutocorrelationTrajectoryResult>();
    r->per_atom_history_.assign(tp.AtomCount(), std::vector<double>{});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Append the current frame's bs_t0 to each atom's history buffer.
// Record the simulation time in frame_times_ so Finalize can derive
// the sample interval for the lag_times_ps emission without reading
// traj.FrameTimes() (which is only populated by Trajectory::Run —
// manually-orchestrated callers don't fill it).

void BsT0AutocorrelationTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj; (void)frame_idx;
    const std::size_t N = tp.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_history_[i].push_back(
            conf.AtomAt(i).bs_shielding_contribution.T0);
    }
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
//
// Canonical biased autocorrelation with full-range mean.
//
//   μ    = (1/N) Σ x_t                      (over all N frames)
//   C(0) = Var(x) = (1/N) Σ (x_t − μ)²
//   C(k) = (1/N) Σ_{t=0..N-k-1} (x_t − μ)(x_{t+k} − μ)
//   ρ(k) = C(k) / C(0)                      (∈ [−1, +1] by
//                                            Cauchy-Schwarz on the
//                                            biased covariance)
//
// Why biased with full-range μ:
//   1. |ρ(k)| ≤ 1 at all finite N (property we want).
//   2. FT[C(·)] is a non-negative power spectrum (Wiener-Khinchin),
//      which is the right input for J(ω) / relaxation-time analysis
//      in NMR.
//   3. Matches the convention used by MDAnalysis, numpy.correlate
//      (1/N normalisation), and the MD-NMR literature (Lipari-Szabo
//      order parameters, spectral density mapping).
//
// Constant-signal atoms (C(0) < 1e-15): ρ is undefined and we emit
// zero at all lags. Tests that check ρ(0) = 1 partition atoms into
// variance-bearing and constant.
//
// Lag-range: k ∈ [0, N_LAGS). Entries where k ≥ N (lag exceeds
// trajectory length) have no pairs to contribute; we emit zero.
// Entries where k < N but N−k is small (say N−k < 10) are noisy but
// mathematically valid — the |ρ| ≤ 1 bound still holds.

void BsT0AutocorrelationTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                   Trajectory& traj) {
    (void)traj;
    const std::size_t N_atoms = tp.AtomCount();

    // Sample interval: median of consecutive Δt from our own
    // frame_times_ record. Falls back to 0 if fewer than 2 frames.
    sample_interval_ps_ = 0.0;
    if (frame_times_.size() >= 2) {
        std::vector<double> deltas;
        deltas.reserve(frame_times_.size() - 1);
        for (std::size_t t = 1; t < frame_times_.size(); ++t) {
            deltas.push_back(frame_times_[t] - frame_times_[t - 1]);
        }
        std::sort(deltas.begin(), deltas.end());
        sample_interval_ps_ = deltas[deltas.size() / 2];
    }

    auto buffer = std::make_unique<DenseBuffer<double>>(N_atoms, N_LAGS);

    for (std::size_t i = 0; i < N_atoms; ++i) {
        const auto& x = per_atom_history_[i];
        const std::size_t T = x.size();

        if (T < 2) {
            for (std::size_t k = 0; k < N_LAGS; ++k) buffer->At(i, k) = 0.0;
            continue;
        }

        // μ = full-range mean.
        double sum = 0.0;
        for (double v : x) sum += v;
        const double mean = sum / static_cast<double>(T);

        // C(0) = (1/T) Σ (x_t − μ)². Also caches (x_t − μ) for reuse.
        std::vector<double> dev(T);
        double c0 = 0.0;
        for (std::size_t t = 0; t < T; ++t) {
            dev[t] = x[t] - mean;
            c0 += dev[t] * dev[t];
        }
        c0 /= static_cast<double>(T);

        if (c0 < 1e-15) {
            // Constant (or numerically near-constant) signal.
            for (std::size_t k = 0; k < N_LAGS; ++k) buffer->At(i, k) = 0.0;
            continue;
        }

        // C(k) = (1/T) Σ_{t=0..T-k-1} dev[t] · dev[t+k].  ρ(k) = C(k) / C(0).
        // Lags beyond the trajectory length are left at 0.
        for (std::size_t k = 0; k < N_LAGS; ++k) {
            if (k >= T) { buffer->At(i, k) = 0.0; continue; }
            double ck = 0.0;
            const std::size_t stop = T - k;
            for (std::size_t t = 0; t < stop; ++t) {
                ck += dev[t] * dev[t + k];
            }
            ck /= static_cast<double>(T);
            buffer->At(i, k) = ck / c0;
        }

        // Release the per-atom history now that we've consumed it.
        std::vector<double>().swap(per_atom_history_[i]);
    }

    tp.AdoptDenseBuffer<double>(
        std::move(buffer),
        std::type_index(typeid(BsT0AutocorrelationTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "BsT0AutocorrelationTrajectoryResult::Finalize",
        "biased ACF across " + std::to_string(N_atoms) + " atoms, " +
        std::to_string(n_frames_) + " frames, " +
        std::to_string(N_LAGS) + " lags; sample interval ≈ " +
        std::to_string(sample_interval_ps_) + " ps");
}


// ── WriteH5Group ─────────────────────────────────────────────────

void BsT0AutocorrelationTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<double>(std::type_index(
            typeid(BsT0AutocorrelationTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "BsT0AutocorrelationTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t L = buffer->StridePerAtom();

    auto grp = file.createGroup("/trajectory/bs_t0_autocorrelation");
    grp.createAttribute("result_name",         Name());
    grp.createAttribute("n_atoms",             N);
    grp.createAttribute("n_lags",              L);
    grp.createAttribute("n_frames",            n_frames_);
    grp.createAttribute("finalized",           finalized_);
    grp.createAttribute("sample_interval_ps",  sample_interval_ps_);
    grp.createAttribute("units",               std::string("dimensionless"));
    grp.createAttribute("estimator",           std::string("biased"));
    grp.createAttribute("mean_convention",     std::string("full_range"));

    // Flat (N, L) double via explicit indexed access (no reinterpret).
    std::vector<double> flat(N * L);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t k = 0; k < L; ++k) {
            flat[i * L + k] = buffer->At(i, k);
        }
    }
    std::vector<std::size_t> dims = {N, L};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("rho", space);
    ds.write_raw(flat.data());

    std::vector<std::uint64_t> lag_frames(L);
    std::vector<double>        lag_times(L);
    for (std::size_t k = 0; k < L; ++k) {
        lag_frames[k] = k;
        lag_times[k]  = static_cast<double>(k) * sample_interval_ps_;
    }
    grp.createDataSet("lag_frames",   lag_frames);
    grp.createDataSet("lag_times_ps", lag_times);
}

}  // namespace nmr
