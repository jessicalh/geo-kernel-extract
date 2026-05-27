#pragma once
//
// RmsdSpikeSelectionTrajectoryResult: per-frame RMSD-spike detector.
// TR12 of the 13-TR plan; near-mechanical clone of
// ChiRotamerSelectionTrajectoryResult (scan-mode SelectionBag emitter)
// with the rotamer-bin transition criterion replaced by a dual-
// threshold RMSD spike criterion.
//
// Dual-threshold spike criterion (OR), with cooldown:
//   - ABSOLUTE: current frame's backbone RMSD vs frame 0 > 1.5 Å.
//     The 1.5 Å threshold is calibrated for small (<100-residue)
//     globular domains at 300 K — SH3 / β-barrel folds typically
//     equilibrate around 1.0-1.5 Å backbone-RMSD; values above 1.5
//     reflect a real conformational excursion rather than equilibrium
//     drift. (Pre-2026-05-21 science review: 2.0 Å was too aggressive
//     for our 60-residue 1P9J SH3 fixture; equilibrated values
//     routinely sat in [1.0, 1.5] Å.) For larger proteins consider
//     gating on per-residue maximum displacement rather than
//     whole-protein RMSD.
//   - LOCAL: |rmsd - rolling_100_frame_mean| > 0.5 Å.
//     Catches sudden excursions from the local plateau, which is
//     what "spike" connotes in trajectory analysis.
//
// Either trigger fires a SelectionRecord push to
// `traj.MutableSelections()`. After firing, the detector enters a
// 100-frame COOLDOWN during which no further spikes are emitted —
// prevents oversampling a single conformational event.
//
// Rolling mean window: 100 frames. The absolute trigger is suppressed
// until at least 10 frames have been observed; the local-delta trigger
// waits for the full 100-frame window.
//
// Per the 13-TR plan: thresholds + cooldown are PROJECT-DECISION
// values (not literature constants); they balance pose diversity
// against trajectory length for DFT pose coordination. The cooldown
// is roughly the conformational-event lifetime at our stride (at
// 20 ps stride: 100 frames ≈ 2 ns).
//
// CROSS-RESULT READ (reader side):
//   Reads `RmsdTrackingTrajectoryResult::LatestRmsd()` during each
//   Compute call. RmsdTracking writes its per-frame rmsd_ vector
//   in-place (AV); attach order = dispatch order so TR11 has already
//   pushed the current frame's RMSD by the time TR12.Compute fires.
//   Alternative considered: this TR computes its own RMSD locally —
//   rejected because (a) it would duplicate the Kabsch SVD per frame,
//   which IS expensive at fleet scale, and (b) the semantic coupling
//   is explicit (this TR's whole purpose depends on the same RMSD
//   distribution that TR11 produces). See PATTERNS.md §17
//   cross-result-read marker discipline.
//
//   Note: an earlier draft used TR11's `RmsdAtFrame(frame_idx)` which
//   silently returned NaN at any stride > 1 (TR11 stores DENSELY by
//   sample order; `frame_idx` is the original TRR frame index, not a
//   positional sample index). Codex round 1 2026-05-21 CRITICAL
//   finding; fixed via the new `LatestRmsd()` accessor.
//
//   Dependencies() returns typeid(RmsdTrackingTrajectoryResult) per
//   Phase 4 validation; attach order = dispatch order so TR11 runs
//   first per frame.
//
// No ConformationResult dependency at this layer (positions used by
// TR11 already).
//
// Emission: pushes SelectionRecord entries to
// `traj.MutableSelections()` with kind=typeid(this); the Trajectory's
// own WriteH5 walks selections_.Kinds() and emits the records under
// `/trajectory/selections/<kind>/`. No private WriteH5Group on this
// TR (the SelectionBag is the H5 surface).
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <deque>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class RmsdSpikeSelectionTrajectoryResult : public TrajectoryResult {
public:
    // Project-decision parameters (per 13-TR plan, 2026-05-21;
    // absolute threshold revised 2.0 -> 1.5 per science adversarial
    // review on 2026-05-21).
    static constexpr double      kAbsoluteThresholdA = 1.5;
    static constexpr double      kLocalDeltaThresholdA = 0.5;
    static constexpr std::size_t kRollingWindowFrames = 100;
    static constexpr std::size_t kCooldownFrames = 100;
    static constexpr std::size_t kMinFramesForRollingMean = 10;

    std::string Name() const override {
        return "RmsdSpikeSelectionTrajectoryResult";
    }

    // CROSS-RESULT READ: requires RmsdTrackingTrajectoryResult.
    // RunConfiguration attaches TR11 before TR12 so it dispatches first.
    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<RmsdSpikeSelectionTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    // Diagnostic accessor — total spikes pushed across the run.
    std::size_t SpikeCount() const { return n_spikes_; }

private:
    // Rolling window of RMSD values. push_back at each Compute;
    // pop_front when len > kRollingWindowFrames.
    std::deque<double> rolling_rmsd_;

    std::size_t frames_until_cooldown_clear_ = 0;
    std::size_t n_spikes_ = 0;
};

}  // namespace nmr
