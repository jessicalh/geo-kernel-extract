#pragma once
//
// BsAnomalousAtomMarkerTrajectoryResult: per-atom per-frame outlier
// marker against the BiotSavart T0 distribution. Worked example of:
//
//   1. TrajectoryResult depending on another TrajectoryResult (not
//      just on per-frame ConformationResults). BsWelfordTrajectoryResult
//      must be attached first; Phase 4's validator enforces this via
//      Dependencies().
//
//   2. Cross-Result read: this Result reads TrajectoryAtom fields
//      owned by BsWelfordTrajectoryResult. See the CROSS-RESULT READ
//      block below.
//
//   3. Pushing to the per-atom event bag (`TrajectoryAtom::events`,
//      Pattern C). The run-scope SelectionBag on Trajectory gets per-
//      frame-level selections (chi transitions, pose candidates, etc);
//      the per-atom bag gets events attributed to a specific atom.
//      Both bags share the same RecordBag<Record> template and the
//      same query grammar.
//
//   4. One emitter producing multiple kinds. Z > +threshold → a
//      BsAnomalyHighT0 kind; Z < -threshold → a BsAnomalyLowT0 kind.
//      Same class emits both; the bag's ByKind<T>() query discriminates
//      at read time.
//
// ── CROSS-RESULT READ ──
//
// This Result reads TrajectoryAtom fields owned by another Result:
//
//   tp.AtomAt(i).bs_t0_mean      (owned by BsWelfordTrajectoryResult)
//   tp.AtomAt(i).bs_t0_m2
//   tp.AtomAt(i).bs_n_frames
//
// Cross-Result reads are allowed by the open-buffer discipline but
// are not the default shape: every cross-read is a human-maintained
// dependency between two Results. The writer-side markers live in
// BsWelfordTrajectoryResult.h; inline read-point markers live in
// this Result's Compute in the .cpp. grep "CROSS-RESULT READ" to
// enumerate every cross-Result dependency in the codebase.
//
// Reason for the cross-read here: the anomaly z-score needs the
// running mean and variance of bs_t0 across the frames seen so far.
// BsWelford already maintains those as its canonical output;
// recomputing them locally would duplicate state and create two
// sources of truth for the same statistic. The cross-read is the
// correct call in this case; if a later design lets anomalies be
// detected against a reference distribution (e.g. from calibration
// data) rather than the running distribution, this Result's
// Dependencies() loses BsWelfordTrajectoryResult and the cross-read
// lifts out cleanly.
//
// ── Attach-order and the running-mean bias ──
//
// Attach order is dispatch order, so BsWelford's Compute runs before
// this Result's Compute on each frame. The running mean includes the
// current frame's contribution when we compare against it; at `n`
// frames the z-score magnitude is biased by ~1/n toward zero. With
// MIN_BURN_IN_FRAMES = 20 the bias is below the threshold noise
// floor. Documented rather than worked around; the alternatives
// (attach BsWelford after, or two-pass compute) complicate the
// exemplar without changing the result.
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class BsAnomalousAtomMarkerTrajectoryResult : public TrajectoryResult {
public:
    // Burn-in: skip frames until BsWelford has accumulated enough
    // samples for a meaningful std estimate.
    static constexpr std::size_t MIN_BURN_IN_FRAMES = 20;

    // z-score threshold for flagging. 2.0σ ≈ two-tail 4.5% tail weight.
    static constexpr double Z_THRESHOLD = 2.0;

    // Event-kind tag structs. These exist only to be `typeid`'d;
    // no instances are ever constructed. One emitter, two kinds.
    struct BsAnomalyHighT0 {};
    struct BsAnomalyLowT0  {};

    std::string Name() const override {
        return "BsAnomalousAtomMarkerTrajectoryResult";
    }

    // Declares deps on TWO TrajectoryResult-or-ConformationResult
    // types. Phase 4 of Trajectory::Run checks:
    //   - typeid(BsWelfordTrajectoryResult): satisfied iff BsWelford
    //     is in tp.AllResults() (i.e. attached earlier).
    //   - typeid(BiotSavartResult): satisfied iff the per-frame
    //     RunConfiguration.RequiredConformationResultTypes() contains
    //     it (i.e. BS runs per frame).
    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<BsAnomalousAtomMarkerTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    // Running totals for diagnostics / tests. Not emitted to H5 —
    // the authoritative record is the per-atom events bag.
    std::size_t TotalEvents() const { return n_events_; }
    std::size_t HighEvents() const { return n_high_events_; }
    std::size_t LowEvents()  const { return n_low_events_;  }

private:
    std::size_t n_events_      = 0;
    std::size_t n_high_events_ = 0;
    std::size_t n_low_events_  = 0;
};

}  // namespace nmr
