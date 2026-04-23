#pragma once
//
// RunContext: caller-constructed description of a specific trajectory
// run instance. Class (not struct): has a constructor that takes the
// pieces the run needs, plus methods for optional extras.
//
// Replaces the role GromacsRunContext played as a "bag of run-time
// state" — though note those are different scopes. GromacsRunContext
// (dissolved) held bonded params, preloaded EDR, and a frame cursor;
// RunContext here holds a RunConfiguration, an output directory, any
// ad-hoc TrajectoryResult extras, and optional runtime resources
// (AIMNet2 model). Preloaded EDR migrates to Trajectory; bonded
// params migrate to TrajectoryProtein; the cursor is internal to the
// frame handler.
//
// Reading Construction as intent rather than assembly:
//
//   RunContext ctx(RunConfiguration::PerFrameExtractionSet(),
//                  "output/1ubq");
//   ctx.SetAimnet2Model(aimnet2_model);
//   ctx.AttachExtraResult(std::make_unique<CustomTrajectoryResult>());
//
//   Trajectory traj(xtc_path, tpr_path, edr_path);
//   TrajectoryProtein tp;
//   tp.BuildFromTrajectory(dir);
//   traj.Run(tp, ctx);
//

#include "RunConfiguration.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace nmr {

class TrajectoryResult;
struct AIMNet2Model;

class RunContext {
public:
    // Minimal constructor: caller supplies a RunConfiguration and an
    // output directory. Every other parameter is an optional setter.
    RunContext(RunConfiguration config,
               std::filesystem::path output_dir);

    // Attach an ad-hoc TrajectoryResult beyond the configuration's
    // default factories. Common for tests, ad-hoc analyses, and
    // μs-harvester custom selectors.
    void AttachExtraResult(std::unique_ptr<TrajectoryResult> extra);

    // Optional AIMNet2 model (torch-loaded outside this class,
    // shared across conformations).
    void SetAimnet2Model(AIMNet2Model* model) { aimnet2_model_ = model; }

    // Stride: process every N-th frame, skipping (stride-1) between
    // each dispatched frame. Production analysis mode uses stride=2
    // on ~1200-frame trajectories to harvest ~600 frames. Default 1
    // = every frame.
    //
    // Honoured by Trajectory::Run Phase 4: Skip() is called
    // (stride-1) times between each Next(), so the underlying XTC
    // is still read linearly (xdrfile doesn't support seeks) but the
    // per-frame calculator pipeline only fires on sampled frames.
    // Trajectory::FrameIndices() records the original XTC index so
    // downstream consumers see the gap.
    void SetStride(size_t s);
    size_t Stride() const { return stride_; }

    // Accessors.
    const RunConfiguration& Configuration() const { return config_; }
    const std::filesystem::path& OutputDir() const { return output_dir_; }

    // Per-frame RunOptions composed on demand: configuration's
    // per_frame_opts_ plus aimnet2 patched in. Returned by value so
    // no mutable cache state is needed.
    RunOptions PerFrameRunOptions() const;

    // Transfer the extras to Trajectory::Run at the start of Phase 1.
    std::vector<std::unique_ptr<TrajectoryResult>> MoveOutExtraResults();

private:
    RunConfiguration config_;
    std::filesystem::path output_dir_;
    std::vector<std::unique_ptr<TrajectoryResult>> extras_;
    AIMNet2Model* aimnet2_model_ = nullptr;
    size_t stride_ = 1;
};

}  // namespace nmr
