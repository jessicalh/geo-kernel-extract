#pragma once
//
// RunConfiguration: typed description of a trajectory-run shape.
// See OBJECT_MODEL.md (trajectory-scope). Two named static shapes:
// PerFrameExtractionSet and FullFatFrameExtraction (the latter is the
// former with MOPAC-family trajectory results and vacuum Coulomb
// enabled). Each encodes per-frame RunOptions, the list of results to
// build (attach order = dispatch order), the required ConformationResult
// types (validated at Trajectory::Run Phase 4), stride, and mandatory
// session resources.
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
    // How to build one TrajectoryResult, deferred. Deferred because a
    // result sizes its buffers from the TrajectoryProtein, which isn't
    // seeded until Trajectory::Run Phase 2 — so a shape lists how to
    // build each result, not the results themselves.
    using DeferredResult =
        std::function<std::unique_ptr<TrajectoryResult>(const TrajectoryProtein&)>;

    // Named shapes (static).
    static RunConfiguration PerFrameExtractionSet();
    static RunConfiguration FullFatFrameExtraction();

    // Accessors used by Trajectory::Run.
    const std::string& Name() const { return name_; }
    const std::vector<DeferredResult>& ResultsToBuild() const {
        return results_to_build_;
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

    // The one frame-selection knob. Process every N-th frame; Skip() the
    // (N-1) in between. Default 1 = every frame. This is the SINGLE source
    // of truth for trajectory cadence: MOPAC (when enabled), per-frame NPY
    // and per-frame PDB emission all run on exactly the dispatched frames.
    // There is no separate mopac/npy/pdb stride — that granularity was
    // removed 2026-05-31 (it let dispatch=2 hide behind emit=1).
    void SetStride(std::size_t s) { stride_ = (s == 0) ? 1 : s; }
    std::size_t Stride() const { return stride_; }

    // Mutable access for the shapes themselves to populate.
    RunOptions& MutablePerFrameRunOptions() { return per_frame_opts_; }
    void AddTrajectoryResultFactory(DeferredResult f) {
        results_to_build_.push_back(std::move(f));
    }
    void RequireConformationResult(std::type_index tid) {
        required_conf_result_types_.insert(tid);
    }
    void SetName(std::string name) { name_ = std::move(name); }

private:
    std::string name_;
    RunOptions per_frame_opts_;
    std::vector<DeferredResult> results_to_build_;
    std::unordered_set<std::type_index> required_conf_result_types_;
    bool requires_aimnet2_ = false;
    std::size_t stride_ = 1;
};

}  // namespace nmr
