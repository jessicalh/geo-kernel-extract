// TrajectoryAtom: per-atom trajectory-scope DATA STORE.
//
// ⚠ ARCHITECTURAL NOTE — READ BEFORE MODIFYING ⚠
//
// This class replaces GromacsProteinAtom. The old struct held
// Welford / DeltaTracker / TransitionCounter instances as fields
// (accumulator IMPLEMENTATION STATE baked into a per-atom data
// struct). That pattern was wrong and must not be reintroduced under
// new names.
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
//   - Rich per-source structured vectors (paralleling
//     ConformationAtom::ring_neighbours) populated at Finalize.
//
// Each field has EXACTLY ONE WRITER: the TrajectoryResult subclass
// that owns the physics producing it. Enforced by attach discipline
// at TrajectoryProtein::AttachResult.
//
// For the full rationale, see
// spec/WIP_OBJECT_MODEL.md §3 anti-patterns subsection.

#ifndef TRAJECTORY_FRAMEWORK_SKETCH_TRAJECTORY_ATOM_H
#define TRAJECTORY_FRAMEWORK_SKETCH_TRAJECTORY_ATOM_H

#include <cstddef>
#include <limits>

class TrajectoryProtein;

// Per WIP §3: private ctor, owned exclusively by TrajectoryProtein's
// atoms_ vector. Public typed fields. No Position() member (frame
// positions are conformation-scoped).
//
// Field inventory below is the narrative example from WIP §3 under
// "`TrajectoryAtom`: per-atom running buffer". The spec's literal
// declared field set is "defined by the accumulator-kin fields in
// GromacsProteinAtom today, redistributed as the TrajectoryResult
// subclasses land." Since this is a framework sketch — not a production
// catalogue — only the fields written by BsWelfordTrajectoryResult (the
// one worked-example Result required by the task) are declared here.
// Adding a new TrajectoryResult means adding its fields here, per the
// "declared upfront" rule.
class TrajectoryAtom {
    friend class TrajectoryProtein;

public:
    // Written by BsWelfordTrajectoryResult. Always-valid mid-stream:
    // Compute() updates these fields directly each frame.
    double bs_t0_mean = 0.0;
    double bs_t0_m2 = 0.0;   // running sum of squared deviations
    double bs_t0_min = std::numeric_limits<double>::infinity();
    double bs_t0_max = -std::numeric_limits<double>::infinity();
    std::size_t bs_n_frames = 0;

    // Finalize-only. Valid after Trajectory::Run's Finalize phase
    // has invoked BsWelfordTrajectoryResult::Finalize().
    double bs_t0_std = 0.0;

private:
    TrajectoryAtom() = default;
};

#endif // TRAJECTORY_FRAMEWORK_SKETCH_TRAJECTORY_ATOM_H
