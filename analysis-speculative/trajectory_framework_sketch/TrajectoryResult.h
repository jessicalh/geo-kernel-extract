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
// re-read §3 of the WIP object model doc.
//
// For the full pattern, see spec/WIP_OBJECT_MODEL.md §4 and the
// worked example BsWelfordTrajectoryResult.

#ifndef TRAJECTORY_FRAMEWORK_SKETCH_TRAJECTORY_RESULT_H
#define TRAJECTORY_FRAMEWORK_SKETCH_TRAJECTORY_RESULT_H

#include <cstddef>
#include <string>
#include <typeindex>
#include <vector>

#include "Stubs.h"

class TrajectoryProtein;

// Per WIP §4. Four virtuals the spec commits to on the base:
// Name, Dependencies, Compute, Finalize. Plus two optional
// serialisation hooks (WriteFeatures, WriteH5Group) that default
// to no-ops so subclasses opt in.
class TrajectoryResult {
public:
    virtual ~TrajectoryResult() = default;

    // Human-readable name, for logging and AttachResult diagnostics.
    virtual std::string Name() const = 0;

    // Declared dependencies. Type indices may refer to:
    //   - other TrajectoryResult types that must be attached first, AND/OR
    //   - ConformationResult types that must run per frame (RunConfiguration's
    //     declared ConformationResult set must include them).
    // Trajectory::Run Phase 2 validates that every dep is satisfied by
    // either already-attached TrajectoryResults or the RunConfiguration's
    // required ConformationResult set. Per WIP §12 item 4 RESOLVED.
    virtual std::vector<std::type_index> Dependencies() const = 0;

    // Per-frame work. Called once per frame during Trajectory::Run
    // phase 4 iteration, via polymorphic dispatch through this base.
    //
    // conf      : this frame's ProteinConformation, with its
    //             ConformationResults already attached.
    // tp        : the TrajectoryProtein running buffer. This result
    //             reads and writes tp.AtomAt(i) fields it owns;
    //             reads tp.AtomAt(i) fields owned by earlier-attached
    //             TrajectoryResults; reads and writes its own internal state.
    // frame_idx : zero-based frame index in the current traversal.
    // time_ps   : simulation time of this frame.
    //
    // No return value. Errors are logged, not thrown.
    virtual void Compute(const ProteinConformation& conf,
                         TrajectoryProtein& tp,
                         std::size_t frame_idx,
                         double time_ps) = 0;

    // End-of-stream synthesis. Called once after the last frame by
    // Trajectory::Run phase 5.
    // Default no-op for results whose TrajectoryAtom fields are
    // always-valid mid-stream (pure Welfords, transition counters).
    virtual void Finalize(TrajectoryProtein& /*tp*/) {}

    // Self-serialisation. Each result knows its own fields and writes
    // them. The top-level writer iterates attached results and calls
    // these. Default returns 0 / no-op so subclasses opt in.
    virtual int WriteFeatures(const TrajectoryProtein& /*tp*/,
                              const std::string& /*output_dir*/) const {
        return 0;
    }

    virtual void WriteH5Group(const TrajectoryProtein& /*tp*/,
                              HighFive::File& /*file*/) const {}
};

#endif // TRAJECTORY_FRAMEWORK_SKETCH_TRAJECTORY_RESULT_H
