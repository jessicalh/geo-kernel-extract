// Trajectory.cpp
//
// Five-phase Run() per WIP_OBJECT_MODEL.md §5.
//
// The spec's §5 concrete Run() body uses `throw std::runtime_error` in
// Phase 1 / Phase 2 dependency-validation failures. PATTERNS.md §9 says
// "Return codes, not exceptions" for calculator/result code, and §"Exception
// hierarchies" explicitly forbids exception hierarchies in calculator,
// result, or pipeline code — but §9 exempts "external library boundaries"
// and the WIP spec §5 itself uses throws. We follow the WIP spec's §5
// shape literally; see README's "Spec contradictions noted".
//
// GromacsFrameHandler is stubbed — its Open/Next/LastConformation/
// LastTimePs/LastXtcIndex surface is forward-declared in a minimal
// shape the sandbox can reference. Real body is process-scope
// infrastructure elsewhere per WIP §9.

#include "Trajectory.h"

#include <stdexcept>
#include <utility>

#include "RunConfiguration.h"
#include "RunContext.h"
#include "SelectionEmittingTrajectoryResult.h"
#include "TrajectoryProtein.h"
#include "TrajectoryResult.h"

// Minimal GromacsFrameHandler interface for the sketch. Real body is
// in src/GromacsFrameHandler.{h,cpp}. The spec §5 calls these methods
// by name; we declare only those the Run() body invokes.
class GromacsFrameHandler {
public:
    GromacsFrameHandler(TrajectoryProtein& /*tp*/,
                        const std::filesystem::path& /*xtc*/,
                        const std::filesystem::path& /*tpr*/) {}

    // Reads frame 0 internally, builds conf0, adds it to
    // tp.Protein().conformations_, runs this frame's ConformationResults
    // per opts, exposes it via LastConformation(). Per WIP §5 Phase 3.
    void Open(const RunOptions& /*opts*/) {}

    // Advances to the next frame; returns false at end of stream.
    // Per WIP §5 Phase 4 `while (handler_->Next(...))`.
    bool Next(const RunOptions& /*opts*/) { return false; }

    const ProteinConformation& LastConformation() const {
        static ProteinConformation sentinel;
        return sentinel;
    }
    double LastTimePs() const { return 0.0; }
    std::size_t LastXtcIndex() const { return 0; }
};

Trajectory::Trajectory(std::filesystem::path xtc_path,
                       std::filesystem::path tpr_path,
                       std::filesystem::path edr_path)
    : xtc_path_(std::move(xtc_path)),
      tpr_path_(std::move(tpr_path)),
      edr_path_(std::move(edr_path)) {}

Trajectory::~Trajectory() = default;

void Trajectory::Run(TrajectoryProtein& tp, RunContext& ctx) {
    if (state_ == State::Complete) {
        throw std::logic_error("Trajectory::Run called on already-completed run");
    }
    state_ = State::Running;

    // === Phase 1: attach TrajectoryResults per the RunConfiguration ===
    for (const auto& factory : ctx.Configuration().TrajectoryResultFactories()) {
        auto result = factory(tp);
        if (!tp.AttachResult(std::move(result))) {
            throw std::runtime_error("Failed to attach TrajectoryResult");
        }
    }
    // Any ad-hoc extras from RunContext (after config factories; §5
    // Phase 1 second part).
    for (auto& extra : ctx.MoveOutExtraResults()) {
        if (!tp.AttachResult(std::move(extra))) {
            throw std::runtime_error("Failed to attach extra TrajectoryResult");
        }
    }

    // === Phase 2: validate ConformationResult deps in the run config ===
    // Per §5 Phase 2: for each attached TrajectoryResult, every declared
    // dependency must be satisfied by either (a) the config's required
    // ConformationResult set, or (b) another already-attached
    // TrajectoryResult.
    for (auto* result : tp.ResultsInAttachOrder()) {
        for (const std::type_index& dep : result->Dependencies()) {
            if (ctx.Configuration().RequiresConformationResult(dep)) {
                continue;  // declared by config
            }
            if (tp.AllResults().count(dep) > 0) {
                continue;  // attached
            }
            throw std::runtime_error(
                std::string("TrajectoryResult ") + result->Name() +
                " declares unmet dependency");
        }
    }

    // === Phase 3: open source, read frame 0, tick frame 0 ===
    handler_ = std::make_unique<GromacsFrameHandler>(tp, xtc_path_, tpr_path_);
    handler_->Open(ctx.PerFrameRunOptions());
    {
        const ProteinConformation& conf0 = handler_->LastConformation();
        for (auto* result : tp.ResultsInAttachOrder()) {
            result->Compute(conf0, tp, 0, handler_->LastTimePs());
        }
        frame_times_.push_back(handler_->LastTimePs());
        frame_indices_.push_back(0);
        frame_count_ = 1;
    }

    // === Phase 4: per-frame loop ===
    while (handler_->Next(ctx.PerFrameRunOptions())) {
        const ProteinConformation& conf = handler_->LastConformation();
        const std::size_t idx = frame_count_;
        const double time_ps = handler_->LastTimePs();

        // Each TrajectoryResult acts on this frame. Iteration is
        // attach-order, which respects declared dependencies (§5
        // "Ordering guarantees").
        for (auto* result : tp.ResultsInAttachOrder()) {
            result->Compute(conf, tp, idx, time_ps);
        }

        frame_times_.push_back(time_ps);
        frame_indices_.push_back(handler_->LastXtcIndex());
        ++frame_count_;
    }

    // === Phase 5: Finalize ===
    for (auto* result : tp.ResultsInAttachOrder()) {
        result->Finalize(tp);
    }

    // Collect selections from SelectionEmittingTrajectoryResult
    // interfaces. Per Appendix D "Collection in Trajectory::Run".
    for (auto* result : tp.ResultsInAttachOrder()) {
        if (auto* emitter =
                dynamic_cast<SelectionEmittingTrajectoryResult*>(result)) {
            const auto& records = emitter->SelectionRecords();
            selections_.insert(selections_.end(), records.begin(), records.end());
        }
    }

    tp.SetFinalized(true);
    state_ = State::Complete;
    handler_.reset();  // Active state cleared; record fields remain.
}

void Trajectory::WriteH5(HighFive::File& /*file*/) const {
    // Per WIP §5 / §7: emits source paths as attributes, frame_times_
    // and frame_indices_ as datasets, selections_ as a record group.
    // Real body depends on HighFive; the shape of the responsibility is
    // what the spec fixes ("writes source metadata and frame records").
    // Framework sketch intentionally leaves the body empty — what goes
    // where is called out in §7, but the HDF5 emission is not framework
    // design; it is writer concern.
}
