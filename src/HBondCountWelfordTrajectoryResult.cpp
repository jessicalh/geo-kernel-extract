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
    result->prev_time_.assign(N, 0.0);
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
    (void)traj;
    const size_t N = tp.AtomCount();

    // Guard against zero dt (frame duplication or stride misconfig);
    // MIN_DT_PS prevents division blowup.
    constexpr double MIN_DT_PS = 1e-12;

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

            // Cadence-normalized rate (codex HIGH finding 2026-05-17):
            // signed Δ alone is sample-rate dependent. dxdt is physically
            // meaningful across runs of different stride.
            const double dt   = time_ps - prev_time_[i];
            const double dxdt = (std::abs(dt) > MIN_DT_PS)
                              ? (delta_x / dt)
                              : 0.0;
            WelfordUpdate(w.count_dxdt, dxdt, dn_new, frame_idx);

            w.delta_n = dn_new;
        }
        prev_value_[i] = c;
        prev_time_[i]  = time_ps;
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
        WelfordFinalize(w.count_dxdt,          w.delta_n);

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

    // Frame-index span this rollup covered. Emitted as H5 group attribute.
    const auto& fidx = traj.FrameIndices();
    if (!fidx.empty()) {
        frame_index_range_ = {fidx.front(), fidx.back()};
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
// /trajectory/hbond_count_welford/ — expanded schema (Phase 2b + Commit C).
//
// Codex review (2026-05-17) surfaced four MEDIUM findings closed
// here:
//
//   (1) Per-dataset `units` attributes — occupancy_fraction_* is
//       "dimensionless" (NOT "pairs"); squared (`*_m2` = pairs²) and
//       rate (`count_dxdt_*` = pairs/ps) get honest unit strings.
//       Group-level `units = "pairs"` kept for backward-compat.
//
//   (2) `frame_index_range` group attribute added.
//
//   (3) Legacy signed-delta channel previously emitted ONLY
//       `delta_mean` + `delta_std`. Now routed via `emit_1d` (full
//       7-stat block as `count_delta_*`) for shape consistency.
//
//   (4) `count_dxdt_*` channel added — codex HIGH finding: signed Δ
//       is sample-rate dependent.
//
//   (5) Verified .h docstring says "expected count ⟨N⟩" — was fixed
//       in commit 3463f5f; no further change needed.
//
// Canonical channels (prefixed):
//   count_*, count_delta_*, count_abs_delta_*, count_delta_squared_*,
//   count_dxdt_*, occupancy_fraction_*
//
// Legacy datasets (unprefixed): mean, std, min, max, delta_mean,
// delta_std — kept for pre-Phase-2b SDK ArraySpec readback, each
// with `units` + `deprecated_use` attributes.
//
// Group attributes: result_name, n_frames, finalized, ddof,
// mean_dt_ps, frame_index_range, units (group-level "pairs"),
// source_radius_A = 3.5.

void HBondCountWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const size_t N = tp.AtomCount();
    auto grp = file.createGroup("/trajectory/hbond_count_welford");

    // ── Group attributes ────────────────────────────────────────
    grp.createAttribute("result_name",       Name());
    grp.createAttribute("n_frames",          n_frames_);
    grp.createAttribute("finalized",         finalized_);
    grp.createAttribute("ddof",              static_cast<int>(1));
    grp.createAttribute("mean_dt_ps",        mean_dt_ps_);
    grp.createAttribute("frame_index_range", frame_index_range_);
    // Group-level `units = "pairs"` reflects the primary count
    // channel. Per-dataset `units` below are authoritative — the
    // occupancy_fraction channel is dimensionless, not "pairs".
    grp.createAttribute("units",             std::string("pairs"));
    grp.createAttribute("source_radius_A",   3.5);

    // ── Helper: emit one WelfordMoments channel as 7 1D datasets ──
    // base_units is the unit of the underlying samples; m2_units is the unit
    // of the Welford M2 accumulator (sum of (sample-mean)²). Caller passes both
    // explicitly because m2's exponent must distribute over every token of a
    // compound unit (e.g. base "pairs_per_ps" → m2 "pairs^2_per_ps^2"). Per
    // codex MEDIUM finding 2026-05-17.
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
    // HBondCount source is integer pairs (promoted to double).
    //   count_*:                pairs,      pairs²
    //   count_delta_*:          pairs,      pairs²
    //   count_abs_delta_*:      pairs,      pairs²
    //   count_delta_squared_*:  pairs²,     pairs⁴   (Welford on Δ²)
    //   count_dxdt_*:           pairs/ps,   pairs²/ps²
    //   occupancy_fraction_*:   dimensionless (indicator in {0, 1})
    // The group-level `units = "pairs"` attribute describes count_*; the
    // per-dataset `units` attributes are authoritative everywhere else
    // (especially occupancy_fraction, which is dimensionless not pairs).
    emit_1d("count",               "pairs",          "pairs^2",          [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).hbond_count_welford.count; });
    emit_1d("count_delta",         "pairs",          "pairs^2",          [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).hbond_count_welford.count_delta; });
    emit_1d("count_abs_delta",     "pairs",          "pairs^2",          [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).hbond_count_welford.count_abs_delta; });
    emit_1d("count_delta_squared", "pairs^2",        "pairs^4",          [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).hbond_count_welford.count_delta_squared; });
    emit_1d("count_dxdt",          "pairs_per_ps",   "pairs^2_per_ps^2", [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).hbond_count_welford.count_dxdt; });
    emit_1d("occupancy_fraction",  "dimensionless",  "dimensionless",    [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).hbond_count_welford.occupancy_fraction; });

    // ── Single scalars and provenance ────────────────────────────
    std::vector<double> rms_delta(N);
    std::vector<size_t> n_frames(N), delta_n(N);
    for (size_t i = 0; i < N; ++i) {
        const HBondCountWelfordState& w = tp.AtomAt(i).hbond_count_welford;
        rms_delta[i] = w.count_rms_delta;
        n_frames[i]  = w.n_frames;
        delta_n[i]   = w.delta_n;
    }
    auto rms_ds = grp.createDataSet("rms_delta", rms_delta);
    rms_ds.createAttribute("units", std::string("pairs"));

    auto nf_ds = grp.createDataSet("n_frames_per_atom", n_frames);
    nf_ds.createAttribute("units", std::string("frame_count"));

    auto dn_ds = grp.createDataSet("delta_n_per_atom", delta_n);
    dn_ds.createAttribute("units", std::string("frame_count"));

    // ── Legacy unprefixed dataset aliases ────────────────────────
    std::vector<double> legacy_mean(N), legacy_std(N), legacy_min(N), legacy_max(N);
    std::vector<double> legacy_delta_mean(N), legacy_delta_std(N);
    for (size_t i = 0; i < N; ++i) {
        const HBondCountWelfordState& w = tp.AtomAt(i).hbond_count_welford;
        legacy_mean[i]       = w.count.mean;
        legacy_std[i]        = w.count.std;
        legacy_min[i]        = w.count.min;
        legacy_max[i]        = w.count.max;
        legacy_delta_mean[i] = w.count_delta.mean;
        legacy_delta_std[i]  = w.count_delta.std;
    }
    auto attach_legacy = [&](HighFive::DataSet ds,
                             const std::string& u,
                             const std::string& canonical) {
        ds.createAttribute("units", u);
        ds.createAttribute("deprecated_use",
            std::string("canonical name is ") + canonical +
            "; this is a pre-Phase-2b alias");
    };
    attach_legacy(grp.createDataSet("mean",       legacy_mean),       std::string("pairs"), "count_mean");
    attach_legacy(grp.createDataSet("std",        legacy_std),        std::string("pairs"), "count_std");
    attach_legacy(grp.createDataSet("min",        legacy_min),        std::string("pairs"), "count_min");
    attach_legacy(grp.createDataSet("max",        legacy_max),        std::string("pairs"), "count_max");
    attach_legacy(grp.createDataSet("delta_mean", legacy_delta_mean), std::string("pairs"), "count_delta_mean");
    attach_legacy(grp.createDataSet("delta_std",  legacy_delta_std),  std::string("pairs"), "count_delta_std");
}

}  // namespace nmr
