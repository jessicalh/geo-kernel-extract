// Trajectory.h
//
// The process entity representing one traversal run. Per
// WIP_OBJECT_MODEL.md §5. Holds source file paths, frame metadata and
// selection records during and after the run, and drives the traversal
// via Run().
//
// Trajectory owns no per-atom state — that's TrajectoryProtein's job.
// Trajectory owns what's about the run itself.

#ifndef TRAJECTORY_FRAMEWORK_SKETCH_TRAJECTORY_H
#define TRAJECTORY_FRAMEWORK_SKETCH_TRAJECTORY_H

#include <cstddef>
#include <filesystem>
#include <memory>
#include <vector>

#include "FrameSelectionRecord.h"
#include "Stubs.h"

class TrajectoryProtein;
class RunContext;
class GromacsFrameHandler;

class Trajectory {
public:
    // Constructor takes source file paths. Records them; does not
    // open them until Run() is called. (§5 constructor contract.)
    Trajectory(std::filesystem::path xtc_path,
               std::filesystem::path tpr_path,
               std::filesystem::path edr_path);
    ~Trajectory();

    // The main method. Drives the five phases per §5:
    //   1. Attach TrajectoryResults per ctx.Configuration().
    //   2. Validate deps against the per-frame ConformationResult set.
    //   3. Open source, read frame 0, finalize protein, tick frame 0.
    //   4. Per-frame loop: each attached TrajectoryResult's Compute
    //      is called once per frame.
    //   5. Finalize: each attached TrajectoryResult's Finalize is
    //      called once. Collect selections from SelectionEmitting
    //      interfaces. Mark tp finalized.
    void Run(TrajectoryProtein& tp, RunContext& ctx);

    // Record fields (valid post-Run).
    std::size_t FrameCount() const { return frame_count_; }
    const std::vector<double>& FrameTimes() const { return frame_times_; }
    const std::vector<std::size_t>& FrameIndices() const { return frame_indices_; }
    const std::vector<FrameSelectionRecord>& Selections() const { return selections_; }
    const std::filesystem::path& XtcPath() const { return xtc_path_; }
    const std::filesystem::path& TprPath() const { return tpr_path_; }
    const std::filesystem::path& EdrPath() const { return edr_path_; }

    bool IsComplete() const { return state_ == State::Complete; }

    // Trajectory-scope H5 metadata. Per §5 and §7: source paths, frame
    // times, frame indices, selections. Does NOT emit per-atom data —
    // that is TrajectoryProtein's business, composed by a top-level
    // writer in production code.
    void WriteH5(HighFive::File& file) const;

private:
    enum class State { Constructed, Running, Complete };
    State state_ = State::Constructed;

    // Sources.
    std::filesystem::path xtc_path_;
    std::filesystem::path tpr_path_;
    std::filesystem::path edr_path_;

    // Record fields.
    std::size_t frame_count_ = 0;
    std::vector<double> frame_times_;
    std::vector<std::size_t> frame_indices_;
    std::vector<FrameSelectionRecord> selections_;

    // Process/active state, valid during Run() only.
    std::unique_ptr<GromacsFrameHandler> handler_;
};

#endif // TRAJECTORY_FRAMEWORK_SKETCH_TRAJECTORY_H
