// FrameSelectionRecord.h
//
// Per WIP_OBJECT_MODEL.md §5 and Appendix D: the record shape that
// SelectionEmittingTrajectoryResult subclasses emit during scan-mode
// runs (ScanForDftPointSet). Trajectory::Run collects these at end of
// Finalize via dynamic_cast through SelectionEmittingTrajectoryResult.

#ifndef TRAJECTORY_FRAMEWORK_SKETCH_FRAME_SELECTION_RECORD_H
#define TRAJECTORY_FRAMEWORK_SKETCH_FRAME_SELECTION_RECORD_H

#include <cstddef>
#include <map>
#include <string>

struct FrameSelectionRecord {
    std::size_t frame_idx = 0;
    double time_ps = 0.0;
    std::string reason;             // e.g. "rotamer_transition_chi1_LYS_42"
    std::string selector;           // which TrajectoryResult flagged it
    // Free-form metadata for downstream consumers (Appendix D).
    std::map<std::string, std::string> metadata;
};

#endif // TRAJECTORY_FRAMEWORK_SKETCH_FRAME_SELECTION_RECORD_H
