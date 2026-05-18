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

void SasaWelfordTrajectoryResult::Compute(
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

            // Cadence-normalized rate (codex HIGH finding 2026-05-17):
            // signed Δ alone is sample-rate dependent. dxdt is physically
            // meaningful across runs of different stride.
            const double dt   = time_ps - prev_time_[i];
            const double dxdt = (std::abs(dt) > MIN_DT_PS)
                              ? (delta_x / dt)
                              : 0.0;
            WelfordUpdate(w.sasa_dxdt, dxdt, dn_new, frame_idx);

            w.delta_n = dn_new;
        }
        prev_value_[i] = s;
        prev_time_[i]  = time_ps;
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
        WelfordFinalize(w.sasa_dxdt,          w.delta_n);

        // NaN when uncomputable (no delta samples); sqrt(0) = 0 when
        // atom is truly static. Downstream isfinite() distinguishes.
        w.sasa_rms_delta = (w.delta_n == 0)
                         ? std::nan("")
                         : std::sqrt(w.sasa_delta_squared.mean);
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
// /trajectory/sasa_welford/ — expanded schema (Phase 2b + Commit C).
//
// Codex review (2026-05-17) surfaced four MEDIUM findings closed
// here:
//
//   (1) Per-dataset `units` attributes — squared-units (`*_m2` = Å⁴)
//       and rate (`sasa_dxdt_*` = Å²/ps) get honest unit strings;
//       group-level `units` kept for backward-compat.
//
//   (2) `frame_index_range` group attribute added.
//
//   (3) Legacy signed-delta channel previously emitted ONLY
//       `delta_mean` + `delta_std`. Now routed via `emit_1d` (full
//       7-stat block as `sasa_delta_*`) for shape consistency.
//
//   (4) `sasa_dxdt_*` channel added — codex HIGH finding: signed Δ
//       is sample-rate dependent.
//
// Canonical channels (prefixed):
//   sasa_*, sasa_delta_*, sasa_abs_delta_*, sasa_delta_squared_*,
//   sasa_dxdt_*
//
// Legacy datasets (unprefixed): mean, std, min, max, delta_mean,
// delta_std — kept for pre-Phase-2b SDK ArraySpec readback, each
// with `units` + `deprecated_use` attributes pointing to canonical.
//
// Group attributes: result_name, n_frames, finalized, ddof,
// mean_dt_ps, frame_index_range, units (group-level "Angstrom^2").

void SasaWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const size_t N = tp.AtomCount();
    auto grp = file.createGroup("/trajectory/sasa_welford");

    // ── Group attributes ────────────────────────────────────────
    grp.createAttribute("result_name",       Name());
    grp.createAttribute("n_frames",          n_frames_);
    grp.createAttribute("finalized",         finalized_);
    grp.createAttribute("ddof",              static_cast<int>(1));
    grp.createAttribute("mean_dt_ps",        mean_dt_ps_);
    grp.createAttribute("frame_index_range", frame_index_range_);
    grp.createAttribute("units",             std::string("Angstrom^2"));

    // ── Helper: emit one WelfordMoments channel as 7 1D datasets ──
    // base_units is the unit of the underlying samples; m2_units is the unit
    // of the Welford M2 accumulator (sum of (sample-mean)²). Caller passes both
    // explicitly because m2's exponent must distribute over every token of a
    // compound unit (e.g. base "Angstrom^2_per_ps" → m2 "Angstrom^4_per_ps^2").
    // Per codex MEDIUM finding 2026-05-17.
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
    // SASA sample is in Å²; rate sample in Å²/ps; squared-Δ sample in Å⁴.
    //   sasa_*:                Å²,     Å⁴
    //   sasa_delta_*:          Å²,     Å⁴
    //   sasa_abs_delta_*:      Å²,     Å⁴
    //   sasa_delta_squared_*:  Å⁴,     Å⁸    (Welford on Δ²; m2 = (Δ²)² = Å⁸)
    //   sasa_dxdt_*:           Å²/ps,  Å⁴/ps²
    emit_1d("sasa",               "Angstrom^2",         "Angstrom^4",          [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).sasa_welford.sasa; });
    emit_1d("sasa_delta",         "Angstrom^2",         "Angstrom^4",          [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).sasa_welford.sasa_delta; });
    emit_1d("sasa_abs_delta",     "Angstrom^2",         "Angstrom^4",          [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).sasa_welford.sasa_abs_delta; });
    emit_1d("sasa_delta_squared", "Angstrom^4",         "Angstrom^8",          [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).sasa_welford.sasa_delta_squared; });
    emit_1d("sasa_dxdt",          "Angstrom^2_per_ps",  "Angstrom^4_per_ps^2", [&](size_t i) -> const WelfordMoments& { return tp.AtomAt(i).sasa_welford.sasa_dxdt; });

    // ── Single scalars and provenance ────────────────────────────
    std::vector<double> rms_delta(N);
    std::vector<size_t> n_frames(N), delta_n(N);
    for (size_t i = 0; i < N; ++i) {
        const SasaWelfordState& w = tp.AtomAt(i).sasa_welford;
        rms_delta[i] = w.sasa_rms_delta;
        n_frames[i]  = w.n_frames;
        delta_n[i]   = w.delta_n;
    }
    auto rms_ds = grp.createDataSet("rms_delta", rms_delta);
    rms_ds.createAttribute("units", std::string("Angstrom^2"));

    auto nf_ds = grp.createDataSet("n_frames_per_atom", n_frames);
    nf_ds.createAttribute("units", std::string("frame_count"));

    auto dn_ds = grp.createDataSet("delta_n_per_atom", delta_n);
    dn_ds.createAttribute("units", std::string("frame_count"));

    // ── Legacy unprefixed dataset aliases ────────────────────────
    std::vector<double> legacy_mean(N), legacy_std(N), legacy_min(N), legacy_max(N);
    std::vector<double> legacy_delta_mean(N), legacy_delta_std(N);
    for (size_t i = 0; i < N; ++i) {
        const SasaWelfordState& w = tp.AtomAt(i).sasa_welford;
        legacy_mean[i]       = w.sasa.mean;
        legacy_std[i]        = w.sasa.std;
        legacy_min[i]        = w.sasa.min;
        legacy_max[i]        = w.sasa.max;
        legacy_delta_mean[i] = w.sasa_delta.mean;
        legacy_delta_std[i]  = w.sasa_delta.std;
    }
    auto attach_legacy = [&](HighFive::DataSet ds,
                             const std::string& u,
                             const std::string& canonical) {
        ds.createAttribute("units", u);
        ds.createAttribute("deprecated_use",
            std::string("canonical name is ") + canonical +
            "; this is a pre-Phase-2b alias");
    };
    attach_legacy(grp.createDataSet("mean",       legacy_mean),       std::string("Angstrom^2"), "sasa_mean");
    attach_legacy(grp.createDataSet("std",        legacy_std),        std::string("Angstrom^2"), "sasa_std");
    attach_legacy(grp.createDataSet("min",        legacy_min),        std::string("Angstrom^2"), "sasa_min");
    attach_legacy(grp.createDataSet("max",        legacy_max),        std::string("Angstrom^2"), "sasa_max");
    attach_legacy(grp.createDataSet("delta_mean", legacy_delta_mean), std::string("Angstrom^2"), "sasa_delta_mean");
    attach_legacy(grp.createDataSet("delta_std",  legacy_delta_std),  std::string("Angstrom^2"), "sasa_delta_std");
}

}  // namespace nmr
