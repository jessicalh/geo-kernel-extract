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
//      Phase 4 guarantees TR11 is attached first per frame.
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
    const auto& rmsd_tracker =
        tp.Result<RmsdTrackingTrajectoryResult>();
    const double rmsd = rmsd_tracker.RmsdAtFrame(frame_idx);
    if (!std::isfinite(rmsd)) {
        // RmsdTracking returned NaN (e.g. fewer than 3 alignment atoms;
        // shouldn't happen on real proteins but defensively skip).
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
    if (rolling_rmsd_.size() < kMinFramesForRollingMean) {
        // Rolling mean unstable; wait.
        return;
    }

    // Compute rolling mean.
    double sum = 0.0;
    for (double v : rolling_rmsd_) sum += v;
    const double mean = sum / static_cast<double>(rolling_rmsd_.size());
    const double local_delta = std::fabs(rmsd - mean);

    const bool abs_trigger   = (rmsd        > kAbsoluteThresholdA);
    const bool local_trigger = (local_delta > kLocalDeltaThresholdA);
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
