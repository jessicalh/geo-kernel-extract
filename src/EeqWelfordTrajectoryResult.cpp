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
    result->prev_time_.assign(N, 0.0);
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
    (void)traj;
    const size_t N = tp.AtomCount();

    // dxdt uses its OWN counter (eeq_welford.dxdt_n) to skip zero-dt
    // frames rather than zero-fill them — codex 2026-05-18: zero-fill
    // biases the rate mean toward zero and inflates variance.
    // 1e-12 ps is well below any physical sampling rate.
    constexpr double MIN_DT_PS = 1e-12;

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

            // Cadence-normalized rate (codex HIGH finding 2026-05-17):
            // signed Δ alone is sample-rate dependent. dxdt is physically
            // meaningful across runs of different stride. Uses its own
            // dxdt_n counter — zero-dt frames are SKIPPED, not zero-filled
            // (codex 2026-05-18: zero-fill biases mean and inflates variance).
            const double dt = time_ps - prev_time_[i];
            if (std::abs(dt) > MIN_DT_PS) {
                const size_t dxdt_n_new = w.dxdt_n + 1;
                WelfordUpdate(w.charge_dxdt, delta_x / dt, dxdt_n_new, frame_idx);
                w.dxdt_n = dxdt_n_new;
            }
        }
        prev_value_[i] = q;
        prev_time_[i]  = time_ps;
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
        // charge_dxdt is finalised against its OWN counter — only frames
        // with dt > MIN_DT_PS contributed (codex 2026-05-18). NaN-fills
        // mean/m2/std when n==0, surfacing the all-zero-dt edge case.
        WelfordFinalize(w.charge_dxdt,          w.dxdt_n);

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

    // Frame-index span this rollup covered. Emitted as H5 group attribute.
    const auto& fidx = traj.FrameIndices();
    if (!fidx.empty()) {
        frame_index_range_ = {fidx.front(), fidx.back()};
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
// /trajectory/eeq_welford/ — expanded schema (Phase 2b + Commit C).
//
// Codex review (2026-05-17) surfaced four MEDIUM findings closed
// here:
//
//   (1) Per-dataset `units` attributes are authoritative — the group-
//       level `units` attribute is kept for backward-compat but is
//       wrong for squared-units (`charge_delta_squared_*` = e²) and
//       rate (`charge_dxdt_*` = e/ps). Per-dataset is honest.
//
//   (2) `frame_index_range` group attribute added (span of trajectory
//       covered).
//
//   (3) The legacy signed-delta channel previously emitted ONLY
//       `delta_mean` + `delta_std`, while `abs_delta` and
//       `delta_squared` emit via `emit_1d` (full 7-stat block).
//       Now routed uniformly via `emit_1d` for shape consistency.
//
//   (4) `charge_dxdt_*` channel added — codex HIGH finding: signed Δ
//       is sample-rate dependent (stride=2 vs stride=300 → 150× change
//       on identical physics); dxdt = Δx/Δt is physically meaningful.
//
// Canonical channels (prefixed):
//   charge_*, charge_delta_*, charge_abs_delta_*,
//   charge_delta_squared_*, charge_dxdt_*
// Each block emits {mean, m2, std, min, max, min_frame, max_frame}
// (7 datasets) with a per-dataset `units` attribute.
//
// Legacy datasets (unprefixed):
//   mean, std, min, max, delta_mean, delta_std, n_frames_per_atom
// Still emitted for pre-Phase-2b SDK ArraySpec readback. Each carries
// a `units` attribute and a `deprecated_use` attribute pointing to
// the canonical prefixed name.
//
// Group attributes:
//   result_name, n_frames, finalized, ddof, mean_dt_ps,
//   frame_index_range, units (group-level, "elementary_charge").

void EeqWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const size_t N = tp.AtomCount();
    auto grp = file.createGroup("/trajectory/eeq_welford");

    // ── Group attributes ────────────────────────────────────────
    grp.createAttribute("result_name",       Name());
    grp.createAttribute("n_frames",          n_frames_);
    grp.createAttribute("finalized",         finalized_);
    grp.createAttribute("ddof",              static_cast<int>(1));
    grp.createAttribute("mean_dt_ps",        mean_dt_ps_);
    grp.createAttribute("frame_index_range", frame_index_range_);
    // Group-level `units` kept for backward-compat with consumers
    // reading the group attribute. Per-dataset `units` attributes
    // below are authoritative (they distinguish e vs e² vs e/ps).
    grp.createAttribute("units",             std::string("elementary_charge"));

    // ── Helper: emit one WelfordMoments channel as 7 1D datasets ──
    // base_units is the unit of the underlying samples (mean / std / min / max);
    // m2_units is the unit of the Welford M2 accumulator (sum of (sample-mean)²).
    // Caller passes both explicitly because `base² ≠ "base^2"` when base already
    // carries an exponent (e.g. for the *_delta_squared channel the sample is
    // already squared so m2 needs base^4) or when base is a compound rate unit
    // (e.g. "elementary_charge_per_ps" → m2 must be "elementary_charge^2_per_ps^2"
    // with the ^2 distributed over each token). Per codex MEDIUM finding
    // 2026-05-17.
    auto emit_1d = [&](const std::string& prefix,
                       const std::string& base_units,
                       const std::string& m2_units,
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
        auto with_units = [&](HighFive::DataSet ds, const std::string& u) {
            ds.createAttribute("units", u);
        };
        with_units(grp.createDataSet(prefix + "_mean",      mean),      base_units);
        with_units(grp.createDataSet(prefix + "_m2",        m2),        m2_units);
        with_units(grp.createDataSet(prefix + "_std",       std_),      base_units);
        with_units(grp.createDataSet(prefix + "_min",       min_),      base_units);
        with_units(grp.createDataSet(prefix + "_max",       max_),      base_units);
        with_units(grp.createDataSet(prefix + "_min_frame", min_frame), std::string("frame_index"));
        with_units(grp.createDataSet(prefix + "_max_frame", max_frame), std::string("frame_index"));
    };

    // ── Canonical prefixed channels ──────────────────────────────
    // Each block: 7 datasets with explicit (base, m2) units.
    //   charge_*:         e,      e²
    //   charge_delta_*:   e,      e²       (signed Δ; same dimension as charge)
    //   charge_abs_delta: e,      e²
    //   charge_delta_squared: e², e⁴      (Welford on Δ² → m2 has e⁴)
    //   charge_dxdt:      e/ps,   e²/ps²  (per-token exponents — unambiguous)
    emit_1d("charge",               "elementary_charge",          "elementary_charge^2",          [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).eeq_welford.charge; });
    emit_1d("charge_delta",         "elementary_charge",          "elementary_charge^2",          [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).eeq_welford.charge_delta; });
    emit_1d("charge_abs_delta",     "elementary_charge",          "elementary_charge^2",          [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).eeq_welford.charge_abs_delta; });
    emit_1d("charge_delta_squared", "elementary_charge^2",        "elementary_charge^4",          [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).eeq_welford.charge_delta_squared; });
    emit_1d("charge_dxdt",          "elementary_charge_per_ps",   "elementary_charge^2_per_ps^2", [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).eeq_welford.charge_dxdt; });

    // ── Single scalars and provenance ────────────────────────────
    std::vector<double> rms_delta(N);
    std::vector<size_t> n_frames(N), delta_n(N), dxdt_n(N);
    for (size_t i = 0; i < N; ++i) {
        const EeqWelfordState& w = tp.AtomAt(i).eeq_welford;
        rms_delta[i] = w.charge_rms_delta;
        n_frames[i]  = w.n_frames;
        delta_n[i]   = w.delta_n;
        dxdt_n[i]    = w.dxdt_n;
    }
    auto rms_ds = grp.createDataSet("rms_delta", rms_delta);
    rms_ds.createAttribute("units", std::string("elementary_charge"));

    auto nf_ds = grp.createDataSet("n_frames_per_atom", n_frames);
    nf_ds.createAttribute("units", std::string("frame_count"));

    auto dn_ds = grp.createDataSet("delta_n_per_atom", delta_n);
    dn_ds.createAttribute("units", std::string("frame_count"));

    // dxdt_n_per_atom: per-atom count of VALID-dt charge_dxdt samples.
    // Codex 2026-05-18 — equals delta_n on well-formed trajectories,
    // diverges when zero-dt rows sneak in. Provenance for rate stats.
    auto dxdtn_ds = grp.createDataSet("dxdt_n_per_atom", dxdt_n);
    dxdtn_ds.createAttribute("units", std::string("frame_count"));

    // ── Legacy unprefixed dataset aliases ────────────────────────
    // Pre-Phase-2b SDK ArraySpec entries may already read these.
    // Each one duplicates the canonical prefixed dataset above and
    // carries `units` + `deprecated_use` pointer to its replacement.
    std::vector<double> legacy_mean(N), legacy_std(N), legacy_min(N), legacy_max(N);
    std::vector<double> legacy_delta_mean(N), legacy_delta_std(N);
    for (size_t i = 0; i < N; ++i) {
        const EeqWelfordState& w = tp.AtomAt(i).eeq_welford;
        legacy_mean[i]       = w.charge.mean;
        legacy_std[i]        = w.charge.std;
        legacy_min[i]        = w.charge.min;
        legacy_max[i]        = w.charge.max;
        legacy_delta_mean[i] = w.charge_delta.mean;
        legacy_delta_std[i]  = w.charge_delta.std;
    }
    auto attach_legacy = [&](HighFive::DataSet ds,
                             const std::string& u,
                             const std::string& canonical) {
        ds.createAttribute("units", u);
        ds.createAttribute("deprecated_use",
            std::string("canonical name is ") + canonical +
            "; this is a pre-Phase-2b alias");
    };
    attach_legacy(grp.createDataSet("mean",       legacy_mean),       std::string("elementary_charge"), "charge_mean");
    attach_legacy(grp.createDataSet("std",        legacy_std),        std::string("elementary_charge"), "charge_std");
    attach_legacy(grp.createDataSet("min",        legacy_min),        std::string("elementary_charge"), "charge_min");
    attach_legacy(grp.createDataSet("max",        legacy_max),        std::string("elementary_charge"), "charge_max");
    attach_legacy(grp.createDataSet("delta_mean", legacy_delta_mean), std::string("elementary_charge"), "charge_delta_mean");
    attach_legacy(grp.createDataSet("delta_std",  legacy_delta_std),  std::string("elementary_charge"), "charge_delta_std");
}

}  // namespace nmr
