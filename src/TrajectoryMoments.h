#pragma once
//
// TrajectoryMoments: inline scalar Welford and min/max maintenance.
//
// Free functions, not a class. The ConformationAtom → TrajectoryAtom
// bridge in each TrajectoryResult::Compute does per-atom online
// scalar accumulation (mean, M2, min, max). These three primitives
// centralise the numerical discipline so the 40-plus TrajectoryResult
// subclasses that follow BsWelfordTrajectoryResult share one
// implementation of the Welford update rather than copy-pasting the
// four-line inline block.
//
// Not a helper class — the project's PATTERNS.md rules out Adapter /
// Wrapper / Proxy / Helper / Bridge naming. These are three one-line
// functions that operate on references into caller-owned
// TrajectoryAtom fields. The Welford state lives on TrajectoryAtom
// in named `WelfordMoments` substructs (one per channel) per
// PATTERNS.md Lesson 25 corollary — see
// spec/plan/welford-data-shape-design-2026-05-17.md.
//
// Numerical notes: the scalar Welford formula is from Welford 1962
// (Technometrics 4(3):419), unchanged since. We compute unbiased
// variance (n-1 denominator) at Finalize, not here.
//

#include <cmath>
#include <cstddef>
#include <limits>

namespace nmr {

// WelfordMoments: per-channel running statistics. One instance per
// scalar channel of a Welford-pattern TR rollup. Used by named
// substruct on TrajectoryAtom (e.g. `bs_welford.t0`, `bs_welford.t1[k]`)
// per the differentiated-structure principle in PATTERNS.md Lesson 25.
//
// Layout: 5 doubles + 2 size_t = 56 bytes per channel. Storage is
// per-atom × per-channel, so 6 Welford TRs × ~10 channels each ×
// 56 bytes × 1500 atoms = ~5 MB per protein for the rollup state.
//
// Field discipline:
//   - mean, m2: running Welford accumulators, updated every frame
//   - std: Finalize-only — derived from m2 / (n-1) where n is the
//     enclosing Welford state's frame count
//   - min/max: running extrema; min_frame/max_frame are the frame
//     indices where the extrema occurred (lets a downstream reader
//     go back to the trajectory)
//
// The shared frame-count denominator lives on the enclosing per-Welford
// state struct (e.g. `bs_welford.n_frames`), not inside the
// WelfordMoments — channels within one Welford share one n.
struct WelfordMoments {
    double mean = 0.0;
    double m2   = 0.0;
    double std  = 0.0;
    double min  =  std::numeric_limits<double>::infinity();
    double max  = -std::numeric_limits<double>::infinity();
    std::size_t min_frame = 0;
    std::size_t max_frame = 0;
};

// Online update of one WelfordMoments channel. Reads x and updates
// mean / m2 (Welford 1962) plus min / max with frame-index
// provenance. Caller provides the frame count after this sample
// (old_n + 1) so the formula is consistent across channels sharing
// one n.
inline void WelfordUpdate(WelfordMoments& w,
                          double x,
                          std::size_t n_new,
                          std::size_t frame_idx) {
    const double delta    = x - w.mean;
    const double new_mean = w.mean + delta / static_cast<double>(n_new);
    w.m2  += delta * (x - new_mean);
    w.mean = new_mean;
    if (x < w.min) { w.min = x; w.min_frame = frame_idx; }
    if (x > w.max) { w.max = x; w.max_frame = frame_idx; }
}

// Finalize a WelfordMoments channel: derive std from m2 / (n-1).
// Idempotent — calling Finalize twice produces the same result
// because std is computed from m2 each time, not accumulated.
inline void WelfordFinalize(WelfordMoments& w, std::size_t n) {
    w.std = (n <= 1) ? 0.0 : std::sqrt(w.m2 / static_cast<double>(n - 1));
}


// Online scalar Welford update. n_new is the frame count AFTER this
// sample (caller must compute it as old_n + 1). Writes to mean/m2
// in place. No branches; safe to use in the hot per-atom loop.
inline void MomentsUpdate(double& mean, double& m2,
                          double x, std::size_t n_new) {
    const double delta     = x - mean;
    const double new_mean  = mean + delta / static_cast<double>(n_new);
    m2   += delta * (x - new_mean);
    mean  = new_mean;
}

// Online min/max with source-frame tracking. Updates both the value
// and the frame index on which the extremum was observed.
inline void MinMaxUpdate(double& min_val, double& max_val,
                         std::size_t& min_frame, std::size_t& max_frame,
                         double x, std::size_t frame_idx) {
    if (x < min_val) { min_val = x; min_frame = frame_idx; }
    if (x > max_val) { max_val = x; max_frame = frame_idx; }
}

// Lighter overload: min/max values only, no frame indices. For
// per-atom fields that don't carry provenance (e.g. delta-trackers
// where the "min delta frame" isn't meaningful).
inline void MinMaxUpdate(double& min_val, double& max_val, double x) {
    if (x < min_val) min_val = x;
    if (x > max_val) max_val = x;
}

// Unbiased standard deviation from accumulated Welford M2 and
// sample count. Returns 0 for n <= 1.
inline double MomentsStd(double m2, std::size_t n) {
    if (n <= 1) return 0.0;
    return std::sqrt(m2 / static_cast<double>(n - 1));
}

}  // namespace nmr
