#pragma once
//
// TrajectoryAtom: per-atom trajectory-scope DATA STORE.
//
// ⚠ ARCHITECTURAL NOTE — READ BEFORE MODIFYING ⚠
//
// This class replaces GromacsProteinAtom (moved to learn/bones/). The
// old struct held Welford / DeltaTracker / TransitionCounter instances
// as fields — accumulator IMPLEMENTATION STATE baked into a per-atom
// data struct. That pattern was wrong and must not be reintroduced
// under new names.
//
// TrajectoryAtom holds finalized OUTPUT fields only — typed summary
// values (double, int, SphericalTensor, small arrays).
//
// DO NOT add as fields on this class:
//   - Welford, DeltaTracker, TransitionCounter, or any other
//     accumulator-state object.
//   - std::vector<std::unique_ptr<Accumulator>> or similar
//     "polymorphic accumulator list" — that is the anti-pattern
//     under virtual dispatch.
//   - Any field whose writer is not clearly identifiable as a
//     specific TrajectoryResult subclass (singleton-per-type).
//
// DO add as fields on this class:
//   - Finalized means / stds / deltas / counts (doubles, ints).
//   - Finalized tensor summaries (SphericalTensor, Mat3, Vec3).
//   - Rich per-source structured vectors populated at Finalize.
//
// Each field has EXACTLY ONE WRITER: the TrajectoryResult subclass
// that owns the physics producing it. Enforced by attach discipline
// at TrajectoryProtein::AttachResult.
//
// For the full rationale, see spec/WIP_OBJECT_MODEL.md §3
// anti-patterns subsection.
//
// Identity (element, residue, bonds, NmrAtomIdentity) goes through
// the wrapped Protein via TrajectoryProtein::Protein().AtomAt(i),
// NOT duplicated here.
//
// Three coexisting shapes hold per-atom trajectory data:
//
//   Pattern A — typed struct vectors for known-shape per-source data.
//               The trajectory-scope parallel to ConformationAtom's
//               ring_neighbours / bond_neighbours / spatial_neighbours.
//               Lands as classes like RingNeighbourhoodTrajectoryStats
//               as the catalog fills in.
//
//   Pattern B — typed accumulator fields, one writer per field
//               (bs_t0_mean, bs_t0_std, bs_t2mag_mean, …). What this
//               file currently holds. Matches the ConformationAtom
//               accumulator discipline lifted across frames.
//
//   Pattern C — the per-atom event bag. Open-shape events emitted by
//               scan-mode detectors and lifetime/transition
//               accumulators (rotamer transitions, hbond form/break,
//               ring flips). Queried with the same grammar as the
//               run-scope SelectionBag on Trajectory — ByKind,
//               ByKindSinceFrame, ByKindSinceTime, MostRecent,
//               CountByKind — at the atom axis.
//
// The three patterns do not compete. A single TrajectoryResult may
// write into any combination of them, depending on the shape of its
// output.
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
