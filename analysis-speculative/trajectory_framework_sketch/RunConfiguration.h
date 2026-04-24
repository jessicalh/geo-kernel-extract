// RunConfiguration.h
//
// First-class typed description of a specific run shape. Per
// WIP_OBJECT_MODEL.md §6. Three named configurations produced by static
// factory methods; each encodes:
//   - per-frame RunOptions (skip flags etc. for the ConformationResult
//     pipeline that runs each frame via OperationRunner),
//   - a list of TrajectoryResult factories (the TrajectoryResults this
//     configuration wants attached), and
//   - the ConformationResult type_indices the configuration requires
//     per frame (for Trajectory::Run Phase 2 dep validation).

#ifndef TRAJECTORY_FRAMEWORK_SKETCH_RUN_CONFIGURATION_H
#define TRAJECTORY_FRAMEWORK_SKETCH_RUN_CONFIGURATION_H

#include <functional>
#include <memory>
#include <string>
#include <typeindex>
#include <unordered_set>
#include <vector>

#include "Stubs.h"
#include "TrajectoryResult.h"

class TrajectoryProtein;

class RunConfiguration {
public:
    using TrajectoryResultFactory =
        std::function<std::unique_ptr<TrajectoryResult>(const TrajectoryProtein&)>;

    // The three named shapes, per §6.
    static RunConfiguration ScanForDftPointSet();
    static RunConfiguration PerFrameExtractionSet();
    static RunConfiguration FullFatFrameExtraction();

    // Accessors used by Trajectory::Run and validators.
    const std::vector<TrajectoryResultFactory>& TrajectoryResultFactories() const {
        return traj_factories_;
    }
    const RunOptions& PerFrameRunOptions() const { return per_frame_opts_; }
    bool RequiresConformationResult(std::type_index tid) const {
        return required_conf_result_types_.count(tid) > 0;
    }
    const std::string& Name() const { return name_; }

private:
    std::string name_;
    RunOptions per_frame_opts_;
    std::vector<TrajectoryResultFactory> traj_factories_;
    std::unordered_set<std::type_index> required_conf_result_types_;
};

#endif // TRAJECTORY_FRAMEWORK_SKETCH_RUN_CONFIGURATION_H
