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
// TrajectoryAtom fields. The accumulator state lives where it did
// before: on TrajectoryAtom (finalised doubles only, never Welford
// structs).
//
// Numerical notes: the scalar Welford formula is from Welford 1962
// (Technometrics 4(3):419), unchanged since. We compute unbiased
// variance (n-1 denominator) at Finalize, not here.
//

#include <cmath>
#include <cstddef>

namespace nmr {

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
