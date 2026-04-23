#pragma once
//
// SelectionEmittingTrajectoryResult: opt-in mixin for TrajectoryResults
// that record per-frame selection events during Compute (rotamer
// transitions, RMSD spikes, ring-flip candidates, χ₁ bin crossings,
// etc.). Trajectory collects these at end of Run via dynamic_cast.
//
// Pattern per spec/WIP_OBJECT_MODEL.md Appendix D. TrajectoryResults
// that emit selections inherit from BOTH TrajectoryResult AND
// SelectionEmittingTrajectoryResult; those that do not inherit only
// from TrajectoryResult.
//
// Why a separate interface instead of a parameter on Compute:
//   - Most TrajectoryResults do not emit selections; adding a
//     parameter they must accept but never use is signature pollution.
//   - The base TrajectoryResult class should not depend on
//     FrameSelectionRecord, Trajectory, or any process-scope type.
//   - dynamic_cast in one collection loop is the explicit "who knows
//     about selections" check. Easy to audit, easy to grep, not
//     scattered across every TrajectoryResult subclass.
//

#include <map>
#include <string>
#include <vector>

namespace nmr {

// Selection record: what was flagged, when, why, by whom.
struct FrameSelectionRecord {
    size_t frame_idx = 0;
    double time_ps = 0.0;
    std::string reason;           // e.g., "rotamer_transition_chi1_LYS_42"
    std::string selector;         // which TrajectoryResult flagged it
    std::map<std::string, std::string> metadata;  // free-form for downstream
};

class SelectionEmittingTrajectoryResult {
public:
    virtual ~SelectionEmittingTrajectoryResult() = default;

    // Post-Run accessor. Trajectory::Run iterates attached
    // TrajectoryResults, dynamic_casts to this interface, and merges
    // the returned records into Trajectory::selections_.
    virtual const std::vector<FrameSelectionRecord>& SelectionRecords() const = 0;
};

}  // namespace nmr
