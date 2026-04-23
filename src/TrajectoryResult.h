#pragma once
//
// TrajectoryResult: base class for per-trajectory modular calculators.
// Parallel to ConformationResult at conformation scope.
//
// ⚠ ARCHITECTURAL NOTE — READ BEFORE MODIFYING ⚠
//
// TrajectoryResult subclasses OWN their per-frame accumulator state.
// Welford, DeltaTracker, TransitionCounter, rolling windows, FFT
// buffers — all internal state of the subclass, populated during
// Compute, dissolved or transferred at Finalize.
//
// They do NOT write accumulator state onto TrajectoryAtom. They
// write finalized OUTPUT fields onto TrajectoryAtom, one writer per
// field, enforced by singleton-per-type discipline at attach time.
//
// If you find yourself writing
//   struct MyTrajectoryResult { std::vector<Welford> per_atom_; };
// that is CORRECT (accumulator state owned by the result).
//
// If you find yourself writing
//   struct TrajectoryAtom { Welford bs_T0; /* ... */ };
// that is WRONG (accumulator state on the per-atom struct). Stop;
// re-read spec/WIP_OBJECT_MODEL.md §3 anti-patterns.
//
// For the full pattern, see spec/WIP_OBJECT_MODEL.md §4 and the
// worked example BsWelfordTrajectoryResult.
//

#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace HighFive { class File; }

namespace nmr {

class ProteinConformation;
class Trajectory;
class TrajectoryProtein;

class TrajectoryResult {
public:
    virtual ~TrajectoryResult() = default;

    // Human-readable name, for logging and AttachResult diagnostics.
    virtual std::string Name() const = 0;

    // Declared dependencies. Type indices of:
    //   - other TrajectoryResult types that must be attached first, AND/OR
    //   - ConformationResult types that must be run per frame (i.e., the
    //     RunConfiguration's per-frame calculator set must include them).
    // Attach-time check fires if any declared dependency is missing.
    virtual std::vector<std::type_index> Dependencies() const = 0;

    // THE work method. Called once per frame during streaming.
    //
    // Arguments:
    //   conf       — this frame's ProteinConformation, with its
    //                ConformationResults already attached and its
    //                ConformationAtoms populated.
    //   tp         — the TrajectoryProtein running buffer. Reads from
    //                and writes to tp.AtomAt(i) fields this result
    //                owns, and to ta.events on any atom for per-atom
    //                event emission. Reads from tp.AtomAt(i) fields
    //                owned by earlier-attached TrajectoryResults.
    //   traj       — the Trajectory driving this run. Provides access
    //                to run-scope state: env (current frame's solvent,
    //                energy pointer), the selection bag
    //                (traj.MutableSelections().Push(...) for scan-mode
    //                event emission), EDR energy lookups, frame
    //                record. Only valid for the duration of Run.
    //   frame_idx  — zero-based frame index in the current traversal.
    //   time_ps    — simulation time of this frame.
    //
    // No return value. Errors are logged via OperationLog, not thrown;
    // a failing Compute either skips its update for this frame or
    // marks its own output invalid via a state flag.
    virtual void Compute(const ProteinConformation& conf,
                         TrajectoryProtein& tp,
                         Trajectory& traj,
                         std::size_t frame_idx,
                         double time_ps) = 0;

    // End-of-stream synthesis. Called once after the last frame.
    // Default no-op for results whose TrajectoryAtom fields are
    // always-valid mid-stream (pure Welfords, transition counters).
    // Override for results that need end-of-stream work:
    //   - FFT over a buffered window
    //   - division by (n-1) for unbiased variance
    //   - transferring dense-buffer ownership to TrajectoryProtein
    //   - computing derived fields from accumulated state
    //   - reading traj.Selections() and pushing a reduced set back
    //     (the DftPoseCoordinator pattern)
    virtual void Finalize(TrajectoryProtein& tp, Trajectory& traj) {
        (void)tp; (void)traj;
    }

    // Self-serialisation to NPY arrays. Same discipline as
    // ConformationResult — each result knows its own fields and
    // writes them. Default returns 0 (no arrays written).
    virtual int WriteFeatures(const TrajectoryProtein& tp,
                              const std::string& output_dir) const {
        (void)tp; (void)output_dir;
        return 0;
    }

    // H5 group emission — parallel to WriteFeatures but targeting the
    // analysis H5 schema. Each result writes its own group.
    virtual void WriteH5Group(const TrajectoryProtein& tp,
                              HighFive::File& file) const {
        (void)tp; (void)file;
    }
};

}  // namespace nmr
