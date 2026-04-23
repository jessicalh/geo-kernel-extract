#pragma once
//
// SelectionRecord: one entry in the run-scope SelectionBag on
// Trajectory. Emitted by TrajectoryResults (typically during scan-mode
// Compute) to flag that this frame is interesting — rotamer
// transition, RMSD spike, DFT pose candidate, ring flip, SS8 change,
// etc.
//
// Reducer TrajectoryResults (DftPoseCoordinator and future kin) read
// the bag at Finalize via traj.Selections().ByKind<...>(), deduplicate
// and budget, and push their own reduced set back to the bag under
// their own kind. The final H5 writer walks traj.Selections().Kinds()
// and emits one group per kind.
//
// Fields:
//   kind       — typeid of the emitting TrajectoryResult.
//   frame_idx  — original XTC frame index (matches
//                Trajectory::FrameIndices).
//   time_ps    — simulation time of the frame.
//   reason     — human-readable cause, e.g.
//                "chi1_transition_LYS42_A_to_B".
//   metadata   — free-form k/v for emitter-specific extras (atom_idx,
//                bin_before, bin_after, ...). frame_idx and time_ps
//                are NOT stashed here — they are top-level so the bag
//                can do windowed queries without string lookups.
//

#include <cstddef>
#include <map>
#include <string>
#include <typeindex>
#include <utility>

namespace nmr {

struct SelectionRecord {
    std::type_index kind;
    std::size_t frame_idx = 0;
    double time_ps = 0.0;
    std::string reason;
    std::map<std::string, std::string> metadata;

    SelectionRecord(std::type_index k,
                    std::size_t fi,
                    double t,
                    std::string r = {},
                    std::map<std::string, std::string> m = {})
        : kind(k), frame_idx(fi), time_ps(t),
          reason(std::move(r)), metadata(std::move(m)) {}
};

}  // namespace nmr
