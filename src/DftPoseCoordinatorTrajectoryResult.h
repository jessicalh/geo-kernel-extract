#pragma once
//
// DftPoseCoordinatorTrajectoryResult: cross-TR SelectionBag aggregator.
// TR13 of the 13-TR plan; the canonical SelectionBag REDUCER worked
// example.
//
// At Finalize, walks the run-scope SelectionBag for two upstream
// emitter kinds — RmsdSpikeSelectionTrajectoryResult (TR12) and
// ChiRotamerSelectionTrajectoryResult — and pushes a deduplicated
// REDUCED set back into the same bag under THIS class's own kind.
//
// Dedup key: `(residue_index, floor(time_ps / 1000 ps))`. residue_index
// comes from a record's `metadata["residue_index"]` if present
// (ChiRotamer records carry it; RmsdSpike doesn't, so we substitute `-1` —
// distinguishing whole-protein events from per-residue events).
//
// Within each (residue_index, ns_bucket) cell, keep the FIRST record
// encountered (push-order iteration of the SelectionBag — stable).
// This spreads DFT pose candidates across both time and residue
// dimensions without prejudice.
//
// Lifecycle: FO (Finalize-only). Compute is a no-op; the reducer
// fires at Finalize after all emitters have completed.
//
// CROSS-RESULT READ (reader side):
//   Reads `traj.Selections().ByKind<RmsdSpikeSelectionTrajectoryResult>()`
//   and `traj.Selections().ByKind<ChiRotamerSelectionTrajectoryResult>()`
//   at Finalize. Dependencies() lists both writer TRs so Phase 4 validates
//   that they are attached; RunConfiguration attaches them before this
//   reducer. Per PATTERNS.md §17, the cross-result read is warranted:
//   the whole purpose of this TR is to reduce
//   the OTHER TRs' output streams; duplicating their accumulation
//   would be wasteful AND the semantic coupling is explicit.
//
// Emission: pushes SelectionRecord entries onto
// `traj.MutableSelections()` with kind=typeid(this). The Trajectory's
// own WriteH5 walks selections_.Kinds() and emits the reduced set
// under `/trajectory/selections/<kind>/`. No private WriteH5Group.
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class DftPoseCoordinatorTrajectoryResult : public TrajectoryResult {
public:
    // Project-decision parameter (per 13-TR plan, 2026-05-21):
    //   1 ns dedup bucket = 1000 ps absolute trajectory time. Keyed
    //   on `record.time_ps`, not the raw TRR frame index — the
    //   previous frame-index form was wrong at any TRR cadence other
    //   than 20 ps/frame (codex round 2 2026-05-21 HIGH finding).
    static constexpr double kNsBucketPs = 1000.0;

    std::string Name() const override {
        return "DftPoseCoordinatorTrajectoryResult";
    }

    // CROSS-RESULT READ: dedupes records from both upstream emitters.
    // Phase 4 validates both are attached.
    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<DftPoseCoordinatorTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp, Trajectory& traj) override;

    // Diagnostic accessor for tests: total reduced records pushed.
    std::size_t ReducedCount() const { return n_reduced_; }

private:
    std::size_t n_reduced_ = 0;
    // Finalize idempotency guard (codex round 2 2026-05-21 MEDIUM):
    // unlike most TRs, the data this reducer emits lives in the
    // SHARED SelectionBag rather than on TR-private state, so the
    // canonical data-flow short-circuit (have_we_emitted check)
    // requires inspecting the bag for our kind. A simple state flag
    // is the right shape here per `feedback_bounds_check_over_state_flag`'s
    // shared-bag escape hatch — the data flow IS shared-bag-shaped.
    bool finalized_ = false;
};

}  // namespace nmr
