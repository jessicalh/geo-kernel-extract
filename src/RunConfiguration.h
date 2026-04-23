#pragma once
//
// RunConfiguration: typed description of a trajectory-run shape.
//
// Three named configurations as static factories:
//   ScanForDftPointSet      — cheap per-frame, selection-emitting
//                              TrajectoryResults for DFT pose choice
//   PerFrameExtractionSet   — full classical stack every frame, time
//                              series + Welford rollups for analysis H5
//   FullFatFrameExtraction  — PerFrameExtractionSet + MOPAC family on
//                              selected frames (use case D's pass 2)
//
// Each configuration encodes:
//   - per-frame RunOptions for OperationRunner::Run
//   - TrajectoryResult factory list (attached to TrajectoryProtein
//     before the loop starts)
//   - required ConformationResult types (validated against the
//     RunOptions by Trajectory::Run before the loop starts — a
//     TrajectoryResult that declares a ConformationResult type in
//     its Dependencies() must find that type in this set)
//
// "What does this run do?" is answered by reading one factory
// function. Changing a configuration means changing one function.
// Adding a configuration means writing a new static factory.
//
// This session populates only the BsWelfordTrajectoryResult factory
// in PerFrameExtractionSet (and in ScanForDftPointSet, since BS runs
// per-frame there too). The other TrajectoryResult factories remain
// empty placeholders until follow-up sessions migrate their fields.
// See spec/WIP_OBJECT_MODEL.md §6 + Appendix F.
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

    // Caller-supplied runtime resource requirements. True means
    // Trajectory::Run Phase 2 will throw if the corresponding
    // pointer in the resolved RunOptions is null. Rationale: these
    // are MANDATORY per-frame drivers, not optional accelerators —
    // silently proceeding without them produces a different
    // extraction, which is worse than stopping. See user feedback
    // 2026-04-23: AIMNet2 had been silently switched off by a
    // predecessor when misconfigured; the correct behaviour is to
    // force the caller to set it up.
    bool RequiresAimnet2() const { return requires_aimnet2_; }
    void SetRequiresAimnet2(bool b) { requires_aimnet2_ = b; }

    // Stride: process every N-th frame, Skip() the (N-1) in between.
    // Stride is a property of the run shape (how you sample), not of
    // runtime resources or per-frame plumbing — it lives here, beside
    // the named shape that says "this is what a run of this kind does."
    // Default 1 = every frame.
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
