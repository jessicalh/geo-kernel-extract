#include "EeqWelfordTrajectoryResult.h"
#include "EeqResult.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "TrajectoryAtom.h"
#include "TrajectoryMoments.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "Types.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <highfive/H5DataSpace.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cmath>
#include <functional>
#include <limits>
#include <typeinfo>

namespace nmr {


std::vector<std::type_index> EeqWelfordTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(EeqResult)) };
}


std::unique_ptr<EeqWelfordTrajectoryResult>
EeqWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto result = std::make_unique<EeqWelfordTrajectoryResult>();
    const size_t N = tp.AtomCount();
    result->prev_value_.assign(N, 0.0);
    result->prev_valid_.assign(N, false);
    return result;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Phase 2b expansion (2026-05-17): scalar source — no T1/T2 channels.
// Per PATTERNS Lesson 25, the three frame-to-frame delta variants
// distinguish drift (signed Δ telescopes to (x_N - x_0)/(N-1) — endpoint
// drift) from dynamics (|Δ| and Δ² don't telescope; rms_delta = sqrt(<Δ²>)
// is the standard fluctuation measure).

void EeqWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj; (void)time_ps;
    const size_t N = tp.AtomCount();

    for (size_t i = 0; i < N; ++i) {
        const double q = conf.AtomAt(i).eeq_charge;

        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        EeqWelfordState& w = ta.eeq_welford;
        const size_t n_new = w.n_frames + 1;

        WelfordUpdate(w.charge, q, n_new, frame_idx);
        w.n_frames = n_new;

        // Frame-to-frame delta variants (skip first frame per atom)
        if (prev_valid_[i]) {
            const double delta_x   = q - prev_value_[i];
            const double abs_delta = std::abs(delta_x);
            const double delta_sq  = delta_x * delta_x;
            const size_t dn_new    = w.delta_n + 1;

            WelfordUpdate(w.charge_delta,         delta_x,   dn_new, frame_idx);
            WelfordUpdate(w.charge_abs_delta,     abs_delta, dn_new, frame_idx);
            WelfordUpdate(w.charge_delta_squared, delta_sq,  dn_new, frame_idx);
            w.delta_n = dn_new;
        }
        prev_value_[i] = q;
        prev_valid_[i] = true;
    }

    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
// Welford m2 → unbiased std for every channel. Derive charge_rms_delta
// from charge_delta_squared.mean. Derive mean_dt_ps from frame_times.

void EeqWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                         Trajectory& traj) {
    const size_t N = tp.AtomCount();
    for (size_t i = 0; i < N; ++i) {
        EeqWelfordState& w = tp.MutableAtomAt(i).eeq_welford;

        WelfordFinalize(w.charge,               w.n_frames);
        WelfordFinalize(w.charge_delta,         w.delta_n);
        WelfordFinalize(w.charge_abs_delta,     w.delta_n);
        WelfordFinalize(w.charge_delta_squared, w.delta_n);

        // RMS fluctuation = sqrt(<Δ²>). The squared-delta accumulator's
        // mean is the per-atom <Δ²>; sqrt at Finalize gives RMS.
        // NaN when uncomputable (no delta samples); sqrt(0) = 0 when
        // atom is truly static. Downstream isfinite() distinguishes.
        w.charge_rms_delta = (w.delta_n == 0)
                           ? std::nan("")
                           : std::sqrt(w.charge_delta_squared.mean);
    }

    // Cadence metadata: mean Δt across captured frames. Lets downstream
    // convert delta-per-stride to dx/dt in physical units.
    const auto& times = traj.FrameTimes();
    if (times.size() >= 2) {
        mean_dt_ps_ = (times.back() - times.front())
                    / static_cast<double>(times.size() - 1);
    }

    finalized_ = true;

    OperationLog::Info(LogCalcOther, "EeqWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) + " frames, " +
        std::to_string(N) + " atoms, mean_dt_ps=" +
        std::to_string(mean_dt_ps_));
}


// ── WriteFeatures (NPY) ──────────────────────────────────────────
// eeq_welford.npy (N, 7) kept unchanged from pre-Phase-2b for backward
// compatibility with the existing SDK ArraySpec. The expanded H5 schema
// is the canonical interface for new consumers.

int EeqWelfordTrajectoryResult::WriteFeatures(
        const TrajectoryProtein& tp,
        const std::string& output_dir) const {
    const size_t N = tp.AtomCount();
    constexpr size_t K = 7;

    std::vector<double> data(N * K);
    for (size_t i = 0; i < N; ++i) {
        const EeqWelfordState& w = tp.AtomAt(i).eeq_welford;
        data[i * K + 0] = w.charge.mean;
        data[i * K + 1] = w.charge.std;
        data[i * K + 2] = w.charge.min;
        data[i * K + 3] = w.charge.max;
        data[i * K + 4] = w.charge_delta.mean;
        data[i * K + 5] = w.charge_delta.std;
        data[i * K + 6] = static_cast<double>(w.n_frames);
    }

    NpyWriter::WriteFloat64(output_dir + "/eeq_welford.npy",
                            data.data(), N, K);
    return 1;
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// /trajectory/eeq_welford/ — expanded schema (Phase 2b). Datasets:
//
//   Scalar channel (1D shape (N,)):
//     mean, std, min, max         — legacy (kept for backward compat)
//     m2, min_frame, max_frame    — provenance
//
//   Frame-to-frame delta variants (1D):
//     delta_mean, delta_std       — legacy signed Δ (kept)
//     abs_delta_{mean,std,min,max}      — |Δ| (fluctuation amplitude)
//     delta_squared_{mean,std,min,max}  — Δ² (Finalize → rms_delta)
//     rms_delta                          — Finalize-derived sqrt(<Δ²>)
//
//   Provenance (1D):
//     n_frames_per_atom, delta_n_per_atom
//
//   H5 group attributes:
//     result_name, n_frames, finalized
//     ddof = 1                          (unbiased std convention)
//     mean_dt_ps                        (cadence)
//     units = "elementary_charge"
//
// Old dataset names (mean, std, min, max, delta_mean, delta_std,
// n_frames_per_atom) stay emitted for backward compatibility.

void EeqWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const size_t N = tp.AtomCount();
    auto grp = file.createGroup("/trajectory/eeq_welford");

    // ── Group attributes ────────────────────────────────────────
    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_frames",    n_frames_);
    grp.createAttribute("finalized",   finalized_);
    grp.createAttribute("ddof",        static_cast<int>(1));
    grp.createAttribute("mean_dt_ps",  mean_dt_ps_);
    grp.createAttribute("units",       std::string("elementary_charge"));

    // ── Helper: emit one WelfordMoments channel as 7 1D datasets ──
    auto emit_1d = [&](const std::string& prefix,
                       std::function<const WelfordMoments&(size_t)> get) {
        std::vector<double> mean(N), m2(N), std_(N), min_(N), max_(N);
        std::vector<size_t> min_frame(N), max_frame(N);
        for (size_t i = 0; i < N; ++i) {
            const WelfordMoments& w = get(i);
            mean[i]      = w.mean;
            m2[i]        = w.m2;
            std_[i]      = w.std;
            min_[i]      = w.min;
            max_[i]      = w.max;
            min_frame[i] = w.min_frame;
            max_frame[i] = w.max_frame;
        }
        grp.createDataSet(prefix + "_mean",      mean);
        grp.createDataSet(prefix + "_m2",        m2);
        grp.createDataSet(prefix + "_std",       std_);
        grp.createDataSet(prefix + "_min",       min_);
        grp.createDataSet(prefix + "_max",       max_);
        grp.createDataSet(prefix + "_min_frame", min_frame);
        grp.createDataSet(prefix + "_max_frame", max_frame);
    };

    // ── Legacy datasets (unchanged from pre-Phase-2b) ───────────
    // These read straight off the charge channel; emit them by name
    // for backward compatibility before the canonical "charge_*" block.
    std::vector<double> mean(N), std_(N), min_(N), max_(N);
    std::vector<double> delta_mean(N), delta_std(N);
    std::vector<size_t> n_frames(N), delta_n(N);
    for (size_t i = 0; i < N; ++i) {
        const EeqWelfordState& w = tp.AtomAt(i).eeq_welford;
        mean[i]       = w.charge.mean;
        std_[i]       = w.charge.std;
        min_[i]       = w.charge.min;
        max_[i]       = w.charge.max;
        delta_mean[i] = w.charge_delta.mean;
        delta_std[i]  = w.charge_delta.std;
        n_frames[i]   = w.n_frames;
        delta_n[i]    = w.delta_n;
    }
    grp.createDataSet("mean",              mean);
    grp.createDataSet("std",               std_);
    grp.createDataSet("min",               min_);
    grp.createDataSet("max",               max_);
    grp.createDataSet("delta_mean",        delta_mean);
    grp.createDataSet("delta_std",         delta_std);
    grp.createDataSet("n_frames_per_atom", n_frames);
    grp.createDataSet("delta_n_per_atom",  delta_n);

    // ── New canonical channels via emit_1d helper ────────────────
    // Note we don't re-emit "charge_*" prefixed datasets to avoid
    // doubling storage on the primary channel (legacy names above
    // already cover charge mean/std/min/max). What's new and exposed:
    // - provenance m2/min_frame/max_frame on charge
    // - |Δ| and Δ² channels with full 7-stat blocks
    // - rms_delta scalar
    emit_1d("abs_delta",     [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).eeq_welford.charge_abs_delta; });
    emit_1d("delta_squared", [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).eeq_welford.charge_delta_squared; });

    // Provenance for the primary charge channel (m2 + min/max frame indices).
    std::vector<double> charge_m2(N);
    std::vector<size_t> charge_min_frame(N), charge_max_frame(N);
    for (size_t i = 0; i < N; ++i) {
        const WelfordMoments& w = tp.AtomAt(i).eeq_welford.charge;
        charge_m2[i]        = w.m2;
        charge_min_frame[i] = w.min_frame;
        charge_max_frame[i] = w.max_frame;
    }
    grp.createDataSet("m2",        charge_m2);
    grp.createDataSet("min_frame", charge_min_frame);
    grp.createDataSet("max_frame", charge_max_frame);

    // RMS-Δ per atom (Finalize-derived).
    std::vector<double> rms_delta(N);
    for (size_t i = 0; i < N; ++i) {
        rms_delta[i] = tp.AtomAt(i).eeq_welford.charge_rms_delta;
    }
    grp.createDataSet("rms_delta", rms_delta);
}

}  // namespace nmr
