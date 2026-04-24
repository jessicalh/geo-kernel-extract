#pragma once
//
// RunConfiguration: typed description of a trajectory-run shape.
// See OBJECT_MODEL.md (trajectory-scope). Three named static
// factories: ScanForDftPointSet, PerFrameExtractionSet,
// FullFatFrameExtraction. Each encodes per-frame RunOptions, the
// TrajectoryResult factory list (attach order = dispatch order),
// the required ConformationResult types (validated at Trajectory::Run
// Phase 4), stride, and mandatory session resources.
//

#include "OperationRunner.h"   // RunOptions

#include <functional>
#include <memory>
#include <string>
#include <typeindex>
#include <unordered_set>
#include <vector>

namespace nmr {

class TrajectoryResult;
class TrajectoryProtein;

class RunConfiguration {
public:
    using TrajectoryResultFactory =
        std::function<std::unique_ptr<TrajectoryResult>(const TrajectoryProtein&)>;

    // Named shapes (static factories).
    static RunConfiguration ScanForDftPointSet();
    static RunConfiguration PerFrameExtractionSet();
    static RunConfiguration FullFatFrameExtraction();

    // Accessors used by Trajectory::Run.
    const std::string& Name() const { return name_; }
    const std::vector<TrajectoryResultFactory>& TrajectoryResultFactories() const {
        return traj_factories_;
    }
    const RunOptions& PerFrameRunOptions() const { return per_frame_opts_; }

    bool RequiresConformationResult(std::type_index tid) const {
        return required_conf_result_types_.count(tid) > 0;
    }
    const std::unordered_set<std::type_index>& RequiredConformationResultTypes() const {
        return required_conf_result_types_;
    }

    // If true, Trajectory::Run Phase 4 returns kConfigRequiresAimnet2
    // when the Session has no AIMNet2 model loaded. Rationale:
    // silently skipping a mandatory per-frame driver produces a
    // different extraction, which is worse than stopping.
    bool RequiresAimnet2() const { return requires_aimnet2_; }
    void SetRequiresAimnet2(bool b) { requires_aimnet2_ = b; }

    // Process every N-th frame; Skip() the (N-1) in between. Default 1.
    void SetStride(std::size_t s) { stride_ = (s == 0) ? 1 : s; }
    std::size_t Stride() const { return stride_; }

    // Mutable access for the factories themselves to populate.
    RunOptions& MutablePerFrameRunOptions() { return per_frame_opts_; }
    void AddTrajectoryResultFactory(TrajectoryResultFactory f) {
        traj_factories_.push_back(std::move(f));
    }
    void RequireConformationResult(std::type_index tid) {
        required_conf_result_types_.insert(tid);
    }
    void SetName(std::string name) { name_ = std::move(name); }

private:
    std::string name_;
    RunOptions per_frame_opts_;
    std::vector<TrajectoryResultFactory> traj_factories_;
    std::unordered_set<std::type_index> required_conf_result_types_;
    bool requires_aimnet2_ = false;
    std::size_t stride_ = 1;
};

}  // namespace nmr
