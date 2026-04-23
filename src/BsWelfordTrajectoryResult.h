#pragma once
//
// BsWelfordTrajectoryResult: running mean / variance / min / max of
// BiotSavart shielding contribution per atom, accumulated across all
// frames of a trajectory.
//
// First concrete TrajectoryResult per spec/WIP_OBJECT_MODEL.md §4
// worked example. Illustrates the always-valid-mid-stream pattern:
// Compute() updates TrajectoryAtom fields in place each frame; the
// mean/m2/min/max are valid at any point after Compute; std is
// finalised at end-of-stream (division by n-1 for unbiased variance).
//
// Source: each frame, after BiotSavartResult::Compute has run on the
// frame's ProteinConformation, this result reads:
//   conf.AtomAt(i).bs_shielding_contribution.T0       (ppm, isotropic)
//   conf.AtomAt(i).bs_shielding_contribution.T2Magnitude()  (|T2|)
// and updates three Welford accumulators on TrajectoryAtom:
//   bs_t0_*         — isotropic ppm
//   bs_t2mag_*      — |T2| ppm
//   bs_t0_delta_*   — frame-to-frame T0 change
//
// Internal state: per-atom previous-frame T0 cache for the delta
// tracker. Small (one double per atom), stays on the result.
//

#include "TrajectoryResult.h"

#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class BsWelfordTrajectoryResult : public TrajectoryResult {
public:
    // Human-readable name.
    std::string Name() const override {
        return "BsWelfordTrajectoryResult";
    }

    // Per-frame dependency: BiotSavartResult must run each frame so
    // that ConformationAtom::bs_shielding_contribution is populated
    // before this result reads it.
    std::vector<std::type_index> Dependencies() const override;

    // Factory. Constructs with zeroed running state sized to the
    // TrajectoryProtein's atom count.
    static std::unique_ptr<BsWelfordTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    // Per-frame compute. Updates TrajectoryAtom fields in place.
    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    // End-of-stream: finalise std from m2 / (n-1) across all atoms.
    void Finalize(TrajectoryProtein& tp, Trajectory& traj) override;

    size_t NumFrames() const { return n_frames_; }

    // NPY output: bs_welford.npy (N, K) with K columns
    // (t0_mean, t0_std, t0_min, t0_max,
    //  t2mag_mean, t2mag_std, t2mag_min, t2mag_max,
    //  t0_delta_mean, t0_delta_std,
    //  n_frames).
    int WriteFeatures(const TrajectoryProtein& tp,
                      const std::string& output_dir) const override;

    // H5 group: /trajectory/bs_welford/ with per-column datasets.
    void WriteH5Group(const TrajectoryProtein& tp,
                      HighFive::File& file) const override;

private:
    // Per-atom previous-frame T0 value, used to compute delta each
    // frame (delta = T0_now - T0_prev). Valid for atom i only after
    // the first Compute call for that atom (prev_valid_[i] == true).
    std::vector<double> prev_t0_;
    std::vector<bool> prev_valid_;

    size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
