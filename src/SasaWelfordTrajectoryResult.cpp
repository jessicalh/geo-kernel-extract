#include "SasaWelfordTrajectoryResult.h"
#include "SasaResult.h"
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


std::vector<std::type_index> SasaWelfordTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(SasaResult)) };
}


std::unique_ptr<SasaWelfordTrajectoryResult>
SasaWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto result = std::make_unique<SasaWelfordTrajectoryResult>();
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

void SasaWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj; (void)time_ps;
    const size_t N = tp.AtomCount();

    for (size_t i = 0; i < N; ++i) {
        const double s = conf.AtomAt(i).atom_sasa;

        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        SasaWelfordState& w = ta.sasa_welford;
        const size_t n_new = w.n_frames + 1;

        WelfordUpdate(w.sasa, s, n_new, frame_idx);
        w.n_frames = n_new;

        // Frame-to-frame delta variants (skip first frame per atom)
        if (prev_valid_[i]) {
            const double delta_x   = s - prev_value_[i];
            const double abs_delta = std::abs(delta_x);
            const double delta_sq  = delta_x * delta_x;
            const size_t dn_new    = w.delta_n + 1;

            WelfordUpdate(w.sasa_delta,         delta_x,   dn_new, frame_idx);
            WelfordUpdate(w.sasa_abs_delta,     abs_delta, dn_new, frame_idx);
            WelfordUpdate(w.sasa_delta_squared, delta_sq,  dn_new, frame_idx);
            w.delta_n = dn_new;
        }
        prev_value_[i] = s;
        prev_valid_[i] = true;
    }

    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
// Welford m2 → unbiased std for every channel. Derive sasa_rms_delta
// from sasa_delta_squared.mean. Derive mean_dt_ps from frame_times.

void SasaWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                          Trajectory& traj) {
    const size_t N = tp.AtomCount();
    for (size_t i = 0; i < N; ++i) {
        SasaWelfordState& w = tp.MutableAtomAt(i).sasa_welford;

        WelfordFinalize(w.sasa,               w.n_frames);
        WelfordFinalize(w.sasa_delta,         w.delta_n);
        WelfordFinalize(w.sasa_abs_delta,     w.delta_n);
        WelfordFinalize(w.sasa_delta_squared, w.delta_n);

        // RMS fluctuation = sqrt(<Δ²>). The squared-delta accumulator's
        // mean is the per-atom <Δ²>; sqrt at Finalize gives RMS.
        w.sasa_rms_delta = (w.sasa_delta_squared.mean > 0.0)
                         ? std::sqrt(w.sasa_delta_squared.mean)
                         : 0.0;
    }

    // Cadence metadata: mean Δt across captured frames. Lets downstream
    // convert delta-per-stride to dx/dt in physical units.
    const auto& times = traj.FrameTimes();
    if (times.size() >= 2) {
        mean_dt_ps_ = (times.back() - times.front())
                    / static_cast<double>(times.size() - 1);
    }

    finalized_ = true;

    OperationLog::Info(LogCalcOther, "SasaWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) + " frames, " +
        std::to_string(N) + " atoms, mean_dt_ps=" +
        std::to_string(mean_dt_ps_));
}


// ── WriteFeatures (NPY) ──────────────────────────────────────────
// sasa_welford.npy (N, 7) kept unchanged from pre-Phase-2b for backward
// compatibility with the existing SDK ArraySpec.

int SasaWelfordTrajectoryResult::WriteFeatures(
        const TrajectoryProtein& tp,
        const std::string& output_dir) const {
    const size_t N = tp.AtomCount();
    constexpr size_t K = 7;

    std::vector<double> data(N * K);
    for (size_t i = 0; i < N; ++i) {
        const SasaWelfordState& w = tp.AtomAt(i).sasa_welford;
        data[i * K + 0] = w.sasa.mean;
        data[i * K + 1] = w.sasa.std;
        data[i * K + 2] = w.sasa.min;
        data[i * K + 3] = w.sasa.max;
        data[i * K + 4] = w.sasa_delta.mean;
        data[i * K + 5] = w.sasa_delta.std;
        data[i * K + 6] = static_cast<double>(w.n_frames);
    }

    NpyWriter::WriteFloat64(output_dir + "/sasa_welford.npy",
                            data.data(), N, K);
    return 1;
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// /trajectory/sasa_welford/ — expanded schema (Phase 2b). Datasets:
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
//     units = "Angstrom^2"

void SasaWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const size_t N = tp.AtomCount();
    auto grp = file.createGroup("/trajectory/sasa_welford");

    // ── Group attributes ────────────────────────────────────────
    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_frames",    n_frames_);
    grp.createAttribute("finalized",   finalized_);
    grp.createAttribute("ddof",        static_cast<int>(1));
    grp.createAttribute("mean_dt_ps",  mean_dt_ps_);
    grp.createAttribute("units",       std::string("Angstrom^2"));

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
    std::vector<double> mean(N), std_(N), min_(N), max_(N);
    std::vector<double> delta_mean(N), delta_std(N);
    std::vector<size_t> n_frames(N), delta_n(N);
    for (size_t i = 0; i < N; ++i) {
        const SasaWelfordState& w = tp.AtomAt(i).sasa_welford;
        mean[i]       = w.sasa.mean;
        std_[i]       = w.sasa.std;
        min_[i]       = w.sasa.min;
        max_[i]       = w.sasa.max;
        delta_mean[i] = w.sasa_delta.mean;
        delta_std[i]  = w.sasa_delta.std;
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
    emit_1d("abs_delta",     [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).sasa_welford.sasa_abs_delta; });
    emit_1d("delta_squared", [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).sasa_welford.sasa_delta_squared; });

    // Provenance for the primary sasa channel (m2 + min/max frame indices).
    std::vector<double> sasa_m2(N);
    std::vector<size_t> sasa_min_frame(N), sasa_max_frame(N);
    for (size_t i = 0; i < N; ++i) {
        const WelfordMoments& w = tp.AtomAt(i).sasa_welford.sasa;
        sasa_m2[i]        = w.m2;
        sasa_min_frame[i] = w.min_frame;
        sasa_max_frame[i] = w.max_frame;
    }
    grp.createDataSet("m2",        sasa_m2);
    grp.createDataSet("min_frame", sasa_min_frame);
    grp.createDataSet("max_frame", sasa_max_frame);

    // RMS-Δ per atom (Finalize-derived).
    std::vector<double> rms_delta(N);
    for (size_t i = 0; i < N; ++i) {
        rms_delta[i] = tp.AtomAt(i).sasa_welford.sasa_rms_delta;
    }
    grp.createDataSet("rms_delta", rms_delta);
}

}  // namespace nmr
