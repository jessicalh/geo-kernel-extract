#include "McConnellWelfordTrajectoryResult.h"
#include "McConnellResult.h"
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


std::vector<std::type_index> McConnellWelfordTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(McConnellResult)) };
}


std::unique_ptr<McConnellWelfordTrajectoryResult>
McConnellWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto result = std::make_unique<McConnellWelfordTrajectoryResult>();
    const size_t N = tp.AtomCount();
    result->prev_t0_.assign(N, 0.0);
    result->prev_valid_.assign(N, false);
    return result;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Phase 2b expansion (2026-05-17): clones the BS pattern WITHOUT T1.
// McConnell-form has T1 ≡ 0 by construction (PATTERNS Lesson 19);
// the three-term asymmetric form lacks the antisymmetric rank-1
// channel that BS / HM rank-1 outer products carry. T0 + per-component
// T2[5] + |T2| amplitude + three delta variants on T0 follow BS
// exactly.

void McConnellWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj; (void)time_ps;
    const size_t N = tp.AtomCount();

    for (size_t i = 0; i < N; ++i) {
        const SphericalTensor& st = conf.AtomAt(i).mc_shielding_contribution;
        const double t0    = st.T0;
        const double t2mag = st.T2Magnitude();

        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        McConnellWelfordState& w = ta.mc_welford;
        const size_t n_new = w.n_frames + 1;

        // Scalar channels
        WelfordUpdate(w.t0,          t0,    n_new, frame_idx);
        WelfordUpdate(w.t2magnitude, t2mag, n_new, frame_idx);

        // NO per-component T1 — McConnell-form has T1 ≡ 0 (PATTERNS Lesson 19).

        // Per-component T2 (5 components)
        for (size_t k = 0; k < 5; ++k)
            WelfordUpdate(w.t2[k], st.T2[k], n_new, frame_idx);

        w.n_frames = n_new;

        // Frame-to-frame T0 delta variants (skip first frame per atom)
        if (prev_valid_[i]) {
            const double delta_x   = t0 - prev_t0_[i];
            const double abs_delta = std::abs(delta_x);
            const double delta_sq  = delta_x * delta_x;
            const size_t dn_new    = w.delta_n + 1;

            WelfordUpdate(w.t0_delta,         delta_x,   dn_new, frame_idx);
            WelfordUpdate(w.t0_abs_delta,     abs_delta, dn_new, frame_idx);
            WelfordUpdate(w.t0_delta_squared, delta_sq,  dn_new, frame_idx);
            w.delta_n = dn_new;
        }
        prev_t0_[i] = t0;
        prev_valid_[i] = true;
    }

    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
// Welford m2 → unbiased std for every channel. Derive t0_rms_delta
// from t0_delta_squared.mean. Derive mean_dt_ps from frame_times.
// No T1 finalize — McConnell-form has no T1 channel.

void McConnellWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                Trajectory& traj) {
    const size_t N = tp.AtomCount();
    for (size_t i = 0; i < N; ++i) {
        McConnellWelfordState& w = tp.MutableAtomAt(i).mc_welford;

        WelfordFinalize(w.t0,          w.n_frames);
        WelfordFinalize(w.t2magnitude, w.n_frames);
        for (size_t k = 0; k < 5; ++k) WelfordFinalize(w.t2[k], w.n_frames);

        WelfordFinalize(w.t0_delta,         w.delta_n);
        WelfordFinalize(w.t0_abs_delta,     w.delta_n);
        WelfordFinalize(w.t0_delta_squared, w.delta_n);

        // RMS fluctuation = sqrt(<Δ²>). The squared-delta accumulator's
        // mean is the per-atom <Δ²>; sqrt at Finalize gives RMS.
        w.t0_rms_delta = (w.t0_delta_squared.mean > 0.0)
                       ? std::sqrt(w.t0_delta_squared.mean)
                       : 0.0;
    }

    // Cadence metadata: mean Δt across captured frames. Lets
    // downstream convert delta-per-stride to dx/dt in physical units.
    const auto& times = traj.FrameTimes();
    if (times.size() >= 2) {
        mean_dt_ps_ = (times.back() - times.front())
                    / static_cast<double>(times.size() - 1);
    }

    finalized_ = true;

    OperationLog::Info(LogCalcOther, "McConnellWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) + " frames, " +
        std::to_string(N) + " atoms, mean_dt_ps=" +
        std::to_string(mean_dt_ps_));
}


// ── WriteFeatures (NPY) ──────────────────────────────────────────
// mc_welford.npy (N, 11) kept unchanged from pre-Phase-2b for
// backward compatibility with the existing SDK ArraySpec. The
// expanded H5 schema is the canonical interface for new consumers.

int McConnellWelfordTrajectoryResult::WriteFeatures(
        const TrajectoryProtein& tp,
        const std::string& output_dir) const {
    const size_t N = tp.AtomCount();
    constexpr size_t K = 11;

    std::vector<double> data(N * K);
    for (size_t i = 0; i < N; ++i) {
        const McConnellWelfordState& w = tp.AtomAt(i).mc_welford;
        data[i * K + 0]  = w.t0.mean;
        data[i * K + 1]  = w.t0.std;
        data[i * K + 2]  = w.t0.min;
        data[i * K + 3]  = w.t0.max;
        data[i * K + 4]  = w.t2magnitude.mean;
        data[i * K + 5]  = w.t2magnitude.std;
        data[i * K + 6]  = w.t2magnitude.min;
        data[i * K + 7]  = w.t2magnitude.max;
        data[i * K + 8]  = w.t0_delta.mean;
        data[i * K + 9]  = w.t0_delta.std;
        data[i * K + 10] = static_cast<double>(w.n_frames);
    }

    NpyWriter::WriteFloat64(output_dir + "/mc_welford.npy",
                            data.data(), N, K);
    return 1;
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// /trajectory/mc_welford/ — expanded schema (Phase 2b). Same shape
// as BS minus T1 (McConnell-form has T1 ≡ 0 by construction; PATTERNS
// Lesson 19). Datasets:
//
//   Scalar channels (1D shape (N,)):
//     t0_{mean,m2,std,min,max,min_frame,max_frame}
//     t2magnitude_{mean,m2,std,min,max,min_frame,max_frame}
//
//   Per-component channel (2D shape (N, 5)):
//     t2_{mean,m2,std,min,max,min_frame,max_frame}    — irrep layout m = -2, -1, 0, +1, +2
//
//   Frame-to-frame delta variants (1D, on T0):
//     t0_delta_{mean,m2,std,min,max}      — signed Δ (telescopes → drift)
//     t0_abs_delta_{mean,m2,std,min,max}  — |Δ| (fluctuation amplitude)
//     t0_delta_squared_{mean,m2,std,min,max}  — Δ² (Finalize → rms_delta)
//     t0_rms_delta                         — Finalize-derived sqrt(<Δ²>)
//
//   Provenance (1D):
//     n_frames_per_atom, delta_n_per_atom
//
//   H5 group attributes:
//     result_name, n_frames, finalized
//     ddof = 1                                    (unbiased std convention)
//     mean_dt_ps                                  (cadence)
//     irrep_layout_t2                             (component ordering; no t1)
//     units = "Angstrom^-3"                        (McConnell kernel unit)
//
// Old dataset names (t0_mean, t2mag_mean, t0_delta_mean, t0_delta_std,
// n_frames_per_atom) stay emitted for backward compatibility — old
// consumers continue to work; new consumers use the canonical names.

void McConnellWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const size_t N = tp.AtomCount();
    auto grp = file.createGroup("/trajectory/mc_welford");

    // ── Group attributes ────────────────────────────────────────
    // No irrep_layout_t1 — McConnell-form has no T1 channel.
    grp.createAttribute("result_name",     Name());
    grp.createAttribute("n_frames",        n_frames_);
    grp.createAttribute("finalized",       finalized_);
    grp.createAttribute("ddof",            static_cast<int>(1));
    grp.createAttribute("mean_dt_ps",      mean_dt_ps_);
    grp.createAttribute("irrep_layout_t2", std::string("m-2,m-1,m0,m+1,m+2"));
    grp.createAttribute("units",           std::string("Angstrom^-3"));

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

    // ── Helper: emit one per-component channel as 7 2D datasets ──
    auto emit_2d = [&](const std::string& prefix, size_t K,
                       std::function<const WelfordMoments&(size_t, size_t)> get) {
        std::vector<double> mean(N * K), m2(N * K), std_(N * K), min_(N * K), max_(N * K);
        std::vector<size_t> min_frame(N * K), max_frame(N * K);
        for (size_t i = 0; i < N; ++i) {
            for (size_t k = 0; k < K; ++k) {
                const WelfordMoments& w = get(i, k);
                mean[i * K + k]      = w.mean;
                m2[i * K + k]        = w.m2;
                std_[i * K + k]      = w.std;
                min_[i * K + k]      = w.min;
                max_[i * K + k]      = w.max;
                min_frame[i * K + k] = w.min_frame;
                max_frame[i * K + k] = w.max_frame;
            }
        }
        HighFive::DataSpace space({N, K});
        grp.createDataSet<double>(prefix + "_mean",       space).write_raw(mean.data());
        grp.createDataSet<double>(prefix + "_m2",         space).write_raw(m2.data());
        grp.createDataSet<double>(prefix + "_std",        space).write_raw(std_.data());
        grp.createDataSet<double>(prefix + "_min",        space).write_raw(min_.data());
        grp.createDataSet<double>(prefix + "_max",        space).write_raw(max_.data());
        grp.createDataSet<size_t>(prefix + "_min_frame",  space).write_raw(min_frame.data());
        grp.createDataSet<size_t>(prefix + "_max_frame",  space).write_raw(max_frame.data());
    };

    // ── Scalar channels (T0, |T2|) ──────────────────────────────
    emit_1d("t0",          [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).mc_welford.t0; });
    emit_1d("t2magnitude", [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).mc_welford.t2magnitude; });

    // ── Per-component channel (T2 only; no T1 for McConnell) ────
    emit_2d("t2", 5, [&](size_t i, size_t k) -> const WelfordMoments& { return tp.AtomAt(i).mc_welford.t2[k]; });

    // ── Delta variants (T0) ──────────────────────────────────────
    emit_1d("t0_delta",         [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).mc_welford.t0_delta; });
    emit_1d("t0_abs_delta",     [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).mc_welford.t0_abs_delta; });
    emit_1d("t0_delta_squared", [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).mc_welford.t0_delta_squared; });

    // ── Single scalar: RMS-Δ per atom (Finalize-derived) ─────────
    std::vector<double> t0_rms_delta(N);
    std::vector<size_t> n_frames(N), delta_n(N);
    for (size_t i = 0; i < N; ++i) {
        const McConnellWelfordState& w = tp.AtomAt(i).mc_welford;
        t0_rms_delta[i] = w.t0_rms_delta;
        n_frames[i]     = w.n_frames;
        delta_n[i]      = w.delta_n;
    }
    grp.createDataSet("t0_rms_delta",      t0_rms_delta);
    grp.createDataSet("n_frames_per_atom", n_frames);
    grp.createDataSet("delta_n_per_atom",  delta_n);

    // ── Backward-compatibility aliases ───────────────────────────
    // Pre-Phase-2b consumers read t2mag_*; emit those as aliases to
    // t2magnitude_* (same data, deprecated naming).
    std::vector<double> t2mag_mean(N), t2mag_std(N), t2mag_min(N), t2mag_max(N);
    for (size_t i = 0; i < N; ++i) {
        const WelfordMoments& w = tp.AtomAt(i).mc_welford.t2magnitude;
        t2mag_mean[i] = w.mean;
        t2mag_std[i]  = w.std;
        t2mag_min[i]  = w.min;
        t2mag_max[i]  = w.max;
    }
    grp.createDataSet("t2mag_mean", t2mag_mean);
    grp.createDataSet("t2mag_std",  t2mag_std);
    grp.createDataSet("t2mag_min",  t2mag_min);
    grp.createDataSet("t2mag_max",  t2mag_max);
}

}  // namespace nmr
