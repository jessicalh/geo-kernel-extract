#pragma once
//
// BsWelfordTrajectoryResult: running mean / variance / min / max of
// BiotSavart shielding contribution per atom, accumulated across all
// frames of a trajectory.
//
// CROSS-RESULT READ (writer side). These TrajectoryAtom fields,
// written by this Result, are also read by other TrajectoryResult
// classes during their own Compute. Renames or semantics changes to
// the fields below require the reader side to update in lockstep.
//
//   bs_welford.t0.mean, bs_welford.t0.m2, bs_welford.n_frames
//       read by BsAnomalousAtomMarkerTrajectoryResult during Compute
//       (anomaly z-score against the running distribution).
//
// AV-pattern exemplar (see PATTERNS.md §14 + Lesson 25): Compute()
// updates TrajectoryAtom fields in place each frame; mean/m2/min/max
// are valid after any Compute; std is finalised at end-of-stream
// (division by n-1 for unbiased variance).
//
// Source: each frame, after BiotSavartResult::Compute has run on the
// frame's ProteinConformation, this result reads the full
// SphericalTensor at:
//   conf.AtomAt(i).bs_shielding_contribution
// and updates the following Welford channels on
// TrajectoryAtom::bs_welford:
//   t0                 — isotropic scalar
//   t1[3]              — antisymmetric rank-1 (m = -1, 0, +1)
//   t2[5]              — symmetric traceless (m = -2..+2)
//   t2magnitude        — |T2| Frobenius amplitude
//   t0_delta           — signed Δ (telescopes → drift)
//   t0_abs_delta       — |Δ| (fluctuation amplitude)
//   t0_delta_squared   — Δ² (sqrt(mean) → RMS fluctuation, Finalize)
//
// Internal state: per-atom previous-frame T0 cache for the delta
// trackers. Small (one double per atom), stays on the result.
//

#include "TrajectoryResult.h"

#include <array>
#include <cstddef>
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

    // Per-atom previous-frame time (ps), used to compute cadence-
    // normalized dxdt = (x_now - x_prev) / (time_now - time_prev).
    // Valid for atom i only when prev_valid_[i] is true.
    std::vector<double> prev_time_;

    size_t n_frames_ = 0;
    bool finalized_ = false;

    // Cadence metadata — mean Δt between captured frames (ps). Derived
    // at Finalize from traj.FrameTimes(). Emitted as H5 attribute so
    // downstream can convert delta-per-stride to dx/dt in physical
    // units. 0.0 before Finalize.
    double mean_dt_ps_ = 0.0;

    // Frame-index span this rollup covered: [first, last] from
    // traj.FrameIndices(). Lets downstream verify cadence and time
    // alignment against the trajectory. Captured at Finalize, emitted
    // as H5 group attribute. {0, 0} before Finalize.
    std::array<std::size_t, 2> frame_index_range_ = {0, 0};
};

}  // namespace nmr
