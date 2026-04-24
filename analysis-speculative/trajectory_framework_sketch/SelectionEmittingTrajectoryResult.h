// SelectionEmittingTrajectoryResult.h
//
// Opt-in mixin interface per WIP_OBJECT_MODEL.md Appendix D.
//
// TrajectoryResults that detect events during per-frame Compute
// (rotamer transitions, RMSD spikes, chi1 bin crossings, ring-flip
// hints) and need to record FrameSelectionRecords for downstream DFT
// submission inherit from this interface in addition to TrajectoryResult.
//
// Trajectory::Run's Phase 5 Finalize sweep iterates attached results,
// dynamic_casts each to this interface, and merges returned records
// into Trajectory::selections_. See Trajectory.h.
//
// Keeping this separate from TrajectoryResult keeps the base class
// free of process-scope concerns (the reason given in Appendix D).

#ifndef TRAJECTORY_FRAMEWORK_SKETCH_SELECTION_EMITTING_H
#define TRAJECTORY_FRAMEWORK_SKETCH_SELECTION_EMITTING_H

#include <vector>

#include "FrameSelectionRecord.h"

// Per Appendix D, the spec shows the return type as std::span<const
// FrameSelectionRecord>. C++17 has no std::span; we return
// const std::vector<FrameSelectionRecord>& for equivalent semantics.
// See README's "Choices the spec did not fully specify".
class SelectionEmittingTrajectoryResult {
public:
    virtual ~SelectionEmittingTrajectoryResult() = default;
    virtual const std::vector<FrameSelectionRecord>& SelectionRecords() const = 0;
};

#endif // TRAJECTORY_FRAMEWORK_SKETCH_SELECTION_EMITTING_H
