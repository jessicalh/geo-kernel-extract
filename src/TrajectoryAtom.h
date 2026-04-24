#pragma once
//
// TrajectoryAtom: per-atom trajectory-scope data store. See
// OBJECT_MODEL.md (trajectory-scope) and PATTERNS.md §13. Private
// constructor; only TrajectoryProtein constructs via friend.
//
// Invariant worth stating: accumulator implementation state
// (Welford, DeltaTracker, TransitionCounter, rolling window) lives
// inside the owning TR, not on this struct. This struct holds
// finalized output fields with one writer per field, plus the
// per-atom event bag.
//

#include "AtomEvent.h"
#include "RecordBag.h"

#include <cstddef>
#include <limits>

namespace nmr {

class TrajectoryProtein;

class TrajectoryAtom {
    friend class TrajectoryProtein;
public:
    // =================================================================
    // Written by BsWelfordTrajectoryResult.
    // Always-valid mid-stream: Compute() updates these in-place each
    // frame. bs_t0_std is Finalize-only — undefined before Finalize.
    // =================================================================

    // BiotSavart isotropic shielding contribution (ppm).
    double bs_t0_mean = 0.0;
    double bs_t0_m2 = 0.0;     // Welford sum-of-squared-deviations
    double bs_t0_std = 0.0;    // Finalize-only: sqrt(m2 / (n-1))
    double bs_t0_min = std::numeric_limits<double>::infinity();
    double bs_t0_max = -std::numeric_limits<double>::infinity();
    size_t bs_t0_min_frame = 0;
    size_t bs_t0_max_frame = 0;
    size_t bs_n_frames = 0;

    // BiotSavart |T2| magnitude (ppm).
    double bs_t2mag_mean = 0.0;
    double bs_t2mag_m2 = 0.0;
    double bs_t2mag_std = 0.0;
    double bs_t2mag_min = std::numeric_limits<double>::infinity();
    double bs_t2mag_max = -std::numeric_limits<double>::infinity();
    size_t bs_t2mag_min_frame = 0;
    size_t bs_t2mag_max_frame = 0;

    // Frame-to-frame BiotSavart T0 delta (ppm/frame).
    double bs_t0_delta_mean = 0.0;
    double bs_t0_delta_m2 = 0.0;
    double bs_t0_delta_std = 0.0;
    double bs_t0_delta_min = std::numeric_limits<double>::infinity();
    double bs_t0_delta_max = -std::numeric_limits<double>::infinity();
    size_t bs_t0_delta_n = 0;  // count of delta samples (= n_frames - 1)

    // =================================================================
    // Pattern C — per-atom event bag.
    // Push via events.Push({emitter, kind, frame, time, metadata}) from
    // any TrajectoryResult::Compute or ::Finalize that emits per-atom
    // events at this atom. Queries are the standard RecordBag verbs.
    // =================================================================
    RecordBag<AtomEvent> events;

private:
    explicit TrajectoryAtom() = default;
};

}  // namespace nmr
