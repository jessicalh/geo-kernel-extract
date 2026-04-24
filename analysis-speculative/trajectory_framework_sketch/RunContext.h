// RunContext.h
//
// Per WIP_OBJECT_MODEL.md §6: "Class not struct. Takes its required
// arguments via constructor; lets caller attach ad-hoc extras via
// methods."
//
// Explicitly NOT a struct. Holds a RunConfiguration and an output
// destination, plus optional extras that get moved out by
// Trajectory::Run's Phase 1 after the config's factories fire.
//
// No mutable cache (§6 private section: "No mutable cache.
// PerFrameRunOptions() returns by value, composed fresh at call time").

#ifndef TRAJECTORY_FRAMEWORK_SKETCH_RUN_CONTEXT_H
#define TRAJECTORY_FRAMEWORK_SKETCH_RUN_CONTEXT_H

#include <filesystem>
#include <memory>
#include <utility>
#include <vector>

#include "RunConfiguration.h"
#include "Stubs.h"
#include "TrajectoryResult.h"

class RunContext {
public:
    // Caller supplies a configuration and an output destination.
    // Optional extras can be added after construction.
    RunContext(RunConfiguration config,
               std::filesystem::path output_dir)
        : config_(std::move(config)),
          output_dir_(std::move(output_dir)) {}

    // Add an ad-hoc TrajectoryResult beyond the config's defaults.
    void AttachExtraResult(std::unique_ptr<TrajectoryResult> extra) {
        extras_.push_back(std::move(extra));
    }

    // Set AIMNet2 model (optional).
    void SetAimnet2Model(AIMNet2Model* model) { aimnet2_model_ = model; }

    // Accessors.
    const RunConfiguration& Configuration() const { return config_; }
    const std::filesystem::path& OutputDir() const { return output_dir_; }

    // Per §6: returns by value, composed fresh at call time from
    // config_ + aimnet2_model_. Trivial copy cost.
    RunOptions PerFrameRunOptions() const {
        RunOptions opts = config_.PerFrameRunOptions();
        opts.aimnet2_model = aimnet2_model_;
        return opts;
    }

    // Moved out by Trajectory::Run when it attaches extras.
    std::vector<std::unique_ptr<TrajectoryResult>> MoveOutExtraResults() {
        return std::move(extras_);
    }

private:
    RunConfiguration config_;
    std::filesystem::path output_dir_;
    std::vector<std::unique_ptr<TrajectoryResult>> extras_;
    AIMNet2Model* aimnet2_model_ = nullptr;
};

#endif // TRAJECTORY_FRAMEWORK_SKETCH_RUN_CONTEXT_H
