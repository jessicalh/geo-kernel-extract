#include "RmsdSpikeSelectionTrajectoryResult.h"

#include "OperationLog.h"
#include "RmsdTrackingTrajectoryResult.h"
#include "SelectionRecord.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"

#include <cmath>
#include <typeinfo>

namespace nmr {

std::vector<std::type_index>
RmsdSpikeSelectionTrajectoryResult::Dependencies() const {
    return {typeid(RmsdTrackingTrajectoryResult)};
}

std::unique_ptr<RmsdSpikeSelectionTrajectoryResult>
RmsdSpikeSelectionTrajectoryResult::Create(const TrajectoryProtein& tp) {
    (void)tp;
    return std::make_unique<RmsdSpikeSelectionTrajectoryResult>();
}

// ── Compute ──────────────────────────────────────────────────────
//
// Per frame:
//   1. CROSS-RESULT READ: get current RMSD from RmsdTrackingTrajectoryResult.
//      RunConfiguration attaches TR11 before TR12 so it dispatches first.
//   2. Decrement cooldown counter if active.
//   3. Push to rolling window.
//   4. If not in cooldown AND rolling window has >= 10 frames:
//      check absolute and local thresholds; on either trigger, push
//      a SelectionRecord and enter cooldown.

void RmsdSpikeSelectionTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)conf;

    // CROSS-RESULT READ: per-frame RMSD scalar from TR11.
    // `LatestRmsd()` returns the value TR11.Compute just pushed for
    // this same frame (attach order = dispatch order, Phase 6/7).
    // Earlier code used `RmsdAtFrame(frame_idx)` which silently
    // returned NaN at any stride > 1 because TR11 stores DENSELY by
    // sample order while `frame_idx` is the original TRR frame index
    // (codex round 1 2026-05-21 CRITICAL finding).
    const auto& rmsd_tracker =
        tp.Result<RmsdTrackingTrajectoryResult>();
    const double rmsd = rmsd_tracker.LatestRmsd();
    if (!std::isfinite(rmsd)) {
        // TR11's `LatestRmsd` / `RmsdAtSampleIndex` return NaN when
        // the alignment set has fewer than 3 atoms (Kabsch rotation
        // underdetermined) OR the per-frame rmsd vector is empty
        // (Compute not yet fired). Skip the spike check for that
        // frame — there is no valid RMSD to compare against thresholds.
        return;
    }

    // Decrement cooldown.
    if (frames_until_cooldown_clear_ > 0) {
        --frames_until_cooldown_clear_;
    }

    // Update rolling window.
    rolling_rmsd_.push_back(rmsd);
    if (rolling_rmsd_.size() > kRollingWindowFrames) {
        rolling_rmsd_.pop_front();
    }

    if (frames_until_cooldown_clear_ > 0) {
        // Cooling down; no spike detection this frame.
        return;
    }

    // Spike-detection criteria. The absolute trigger can fire as soon
    // as the rolling window has at least kMinFramesForRollingMean
    // entries (= 10 frames), but the LOCAL-Δ trigger requires a FULL
    // 100-frame window — during early equilibration the rolling mean
    // is non-stationary (RMSD typically rises monotonically), and a
    // 10-30-frame mean would lag the true plateau enough to fire
    // spurious local-Δ spikes. Per math adversarial review 2026-05-21.
    if (rolling_rmsd_.size() < kMinFramesForRollingMean) {
        // Below the absolute-trigger gate too.
        return;
    }

    // Compute rolling mean (used by both triggers' logging; local
    // trigger only consults it when the window is fully populated).
    double sum = 0.0;
    for (double v : rolling_rmsd_) sum += v;
    const double mean = sum / static_cast<double>(rolling_rmsd_.size());
    const double local_delta = std::fabs(rmsd - mean);

    const bool abs_trigger   = (rmsd        > kAbsoluteThresholdA);
    const bool local_trigger =
        (rolling_rmsd_.size() >= kRollingWindowFrames) &&
        (local_delta > kLocalDeltaThresholdA);
    if (!abs_trigger && !local_trigger) return;

    // Spike. Emit SelectionRecord, enter cooldown.
    const char* triggers =
        (abs_trigger && local_trigger) ? "abs+local"
        : abs_trigger                  ? "abs"
                                       : "local";

    std::string reason =
        std::string("rmsd_spike_frame_") + std::to_string(frame_idx) +
        "_triggers_" + triggers +
        "_rmsd_" + std::to_string(rmsd) + "_A";

    traj.MutableSelections().Push(SelectionRecord(
        std::type_index(typeid(RmsdSpikeSelectionTrajectoryResult)),
        frame_idx, time_ps,
        std::move(reason),
        {
            {"rmsd_A",            std::to_string(rmsd)},
            {"rolling_mean_A",    std::to_string(mean)},
            {"local_delta_A",     std::to_string(local_delta)},
            {"window_frames",     std::to_string(rolling_rmsd_.size())},
            {"abs_threshold_A",   std::to_string(kAbsoluteThresholdA)},
            {"local_threshold_A", std::to_string(kLocalDeltaThresholdA)},
            {"trigger",           triggers},
        }));
    ++n_spikes_;

    // Enter cooldown.
    frames_until_cooldown_clear_ = kCooldownFrames;
}

}  // namespace nmr
