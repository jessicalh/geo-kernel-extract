#include "HBondCountWelfordTrajectoryResult.h"
#include "HBondResult.h"
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


std::vector<std::type_index> HBondCountWelfordTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(HBondResult)) };
}


std::unique_ptr<HBondCountWelfordTrajectoryResult>
HBondCountWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto result = std::make_unique<HBondCountWelfordTrajectoryResult>();
    const size_t N = tp.AtomCount();
    result->prev_value_.assign(N, 0.0);
    result->prev_valid_.assign(N, false);
    return result;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Phase 2b expansion (2026-05-17): scalar source + occupancy_fraction
// companion channel. The count channel accumulates the integer H-bond
// count promoted to double; the occupancy_fraction channel accumulates
// the binary indicator (count > 0 ? 1.0 : 0.0). The latter's mean lives
// in [0, 1] and captures "what fraction of frames is this atom engaged
// in any H-bond?" — a distinct physics question from ⟨N⟩ expected count.
//
// Per PATTERNS Lesson 25, the three frame-to-frame delta variants on
// count distinguish drift (signed Δ telescopes) from dynamics (|Δ|, Δ²).

void HBondCountWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj; (void)time_ps;
    const size_t N = tp.AtomCount();

    for (size_t i = 0; i < N; ++i) {
        const double c = static_cast<double>(conf.AtomAt(i).hbond_count_within_3_5A);

        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        HBondCountWelfordState& w = ta.hbond_count_welford;
        const size_t n_new = w.n_frames + 1;

        WelfordUpdate(w.count, c, n_new, frame_idx);

        // Occupancy indicator: 1.0 if any H-bond this frame, else 0.0.
        // Mean → fraction of frames this atom is engaged at all.
        const double indicator = (c > 0.0) ? 1.0 : 0.0;
        WelfordUpdate(w.occupancy_fraction, indicator, n_new, frame_idx);

        w.n_frames = n_new;

        // Frame-to-frame count delta variants (skip first frame per atom)
        if (prev_valid_[i]) {
            const double delta_x   = c - prev_value_[i];
            const double abs_delta = std::abs(delta_x);
            const double delta_sq  = delta_x * delta_x;
            const size_t dn_new    = w.delta_n + 1;

            WelfordUpdate(w.count_delta,         delta_x,   dn_new, frame_idx);
            WelfordUpdate(w.count_abs_delta,     abs_delta, dn_new, frame_idx);
            WelfordUpdate(w.count_delta_squared, delta_sq,  dn_new, frame_idx);
            w.delta_n = dn_new;
        }
        prev_value_[i] = c;
        prev_valid_[i] = true;
    }

    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
// Welford m2 → unbiased std for every channel including occupancy_fraction.
// Derive count_rms_delta from count_delta_squared.mean.
// Derive mean_dt_ps from frame_times.

void HBondCountWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                 Trajectory& traj) {
    const size_t N = tp.AtomCount();
    for (size_t i = 0; i < N; ++i) {
        HBondCountWelfordState& w = tp.MutableAtomAt(i).hbond_count_welford;

        WelfordFinalize(w.count,               w.n_frames);
        WelfordFinalize(w.occupancy_fraction,  w.n_frames);
        WelfordFinalize(w.count_delta,         w.delta_n);
        WelfordFinalize(w.count_abs_delta,     w.delta_n);
        WelfordFinalize(w.count_delta_squared, w.delta_n);

        // NaN when uncomputable (no delta samples); sqrt(0) = 0 when
        // atom is truly static. Downstream isfinite() distinguishes.
        w.count_rms_delta = (w.delta_n == 0)
                          ? std::nan("")
                          : std::sqrt(w.count_delta_squared.mean);
    }

    // Cadence metadata.
    const auto& times = traj.FrameTimes();
    if (times.size() >= 2) {
        mean_dt_ps_ = (times.back() - times.front())
                    / static_cast<double>(times.size() - 1);
    }

    finalized_ = true;

    OperationLog::Info(LogCalcOther, "HBondCountWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) + " frames, " +
        std::to_string(N) + " atoms, mean_dt_ps=" +
        std::to_string(mean_dt_ps_));
}


// ── WriteFeatures (NPY) ──────────────────────────────────────────
// hbond_count_welford.npy (N, 7) kept unchanged from pre-Phase-2b.

int HBondCountWelfordTrajectoryResult::WriteFeatures(
        const TrajectoryProtein& tp,
        const std::string& output_dir) const {
    const size_t N = tp.AtomCount();
    constexpr size_t K = 7;

    std::vector<double> data(N * K);
    for (size_t i = 0; i < N; ++i) {
        const HBondCountWelfordState& w = tp.AtomAt(i).hbond_count_welford;
        data[i * K + 0] = w.count.mean;
        data[i * K + 1] = w.count.std;
        data[i * K + 2] = w.count.min;
        data[i * K + 3] = w.count.max;
        data[i * K + 4] = w.count_delta.mean;
        data[i * K + 5] = w.count_delta.std;
        data[i * K + 6] = static_cast<double>(w.n_frames);
    }

    NpyWriter::WriteFloat64(output_dir + "/hbond_count_welford.npy",
                            data.data(), N, K);
    return 1;
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// /trajectory/hbond_count_welford/ — expanded schema (Phase 2b).
//
//   Primary count channel (1D shape (N,)):
//     mean, std, min, max         — legacy (kept for backward compat)
//     m2, min_frame, max_frame    — provenance
//
//   Occupancy indicator companion channel (1D, mean ∈ [0, 1]):
//     occupancy_fraction_{mean,m2,std,min,max,min_frame,max_frame}
//
//   Frame-to-frame count delta variants (1D):
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
//     units = "pairs"
//     source_radius_A = 3.5

void HBondCountWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const size_t N = tp.AtomCount();
    auto grp = file.createGroup("/trajectory/hbond_count_welford");

    // ── Group attributes ────────────────────────────────────────
    grp.createAttribute("result_name",     Name());
    grp.createAttribute("n_frames",        n_frames_);
    grp.createAttribute("finalized",       finalized_);
    grp.createAttribute("ddof",            static_cast<int>(1));
    grp.createAttribute("mean_dt_ps",      mean_dt_ps_);
    grp.createAttribute("units",           std::string("pairs"));
    grp.createAttribute("source_radius_A", 3.5);

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
        const HBondCountWelfordState& w = tp.AtomAt(i).hbond_count_welford;
        mean[i]       = w.count.mean;
        std_[i]       = w.count.std;
        min_[i]       = w.count.min;
        max_[i]       = w.count.max;
        delta_mean[i] = w.count_delta.mean;
        delta_std[i]  = w.count_delta.std;
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
    emit_1d("occupancy_fraction", [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).hbond_count_welford.occupancy_fraction; });
    emit_1d("abs_delta",          [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).hbond_count_welford.count_abs_delta; });
    emit_1d("delta_squared",      [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).hbond_count_welford.count_delta_squared; });

    // Provenance for the primary count channel (m2 + min/max frame indices).
    std::vector<double> count_m2(N);
    std::vector<size_t> count_min_frame(N), count_max_frame(N);
    for (size_t i = 0; i < N; ++i) {
        const WelfordMoments& w = tp.AtomAt(i).hbond_count_welford.count;
        count_m2[i]        = w.m2;
        count_min_frame[i] = w.min_frame;
        count_max_frame[i] = w.max_frame;
    }
    grp.createDataSet("m2",        count_m2);
    grp.createDataSet("min_frame", count_min_frame);
    grp.createDataSet("max_frame", count_max_frame);

    // RMS-Δ per atom (Finalize-derived).
    std::vector<double> rms_delta(N);
    for (size_t i = 0; i < N; ++i) {
        rms_delta[i] = tp.AtomAt(i).hbond_count_welford.count_rms_delta;
    }
    grp.createDataSet("rms_delta", rms_delta);
}

}  // namespace nmr
