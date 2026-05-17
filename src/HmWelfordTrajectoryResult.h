#pragma once
//
// HmWelfordTrajectoryResult: running mean / variance / min / max of
// HaighMallion shielding kernel per atom, accumulated across all frames
// of a trajectory. AV-pattern clone of BsWelfordTrajectoryResult; the
// HaighMallion kernel coexists with the Biot-Savart kernel and the
// rollups stay side-by-side per `feedback_methods_accumulate`.
//
// Source: `ConformationAtom::hm_shielding_contribution` (units = Å⁻¹,
// rank-1 same as BS but no PPM_FACTOR multiplier — see OBJECT_MODEL.md
// contract drift table). HaighMallionResult is unconditionally
// attached in `PerFrameExtractionSet`, so the dep is enforced via
// `Dependencies()`.
//
// Phase 2b expansion (2026-05-17): identical channel shape to BS —
// per-component T1[3] + T2[5] preserve tensor orientation that |T2|
// amplitude rollup discards; three frame-to-frame delta variants on
// T0 distinguish drift / |Δ| / Δ² → RMS fluctuation. Per PATTERNS
// Lesson 25 (Export Everything Upstream). See
// spec/plan/welford-data-shape-design-2026-05-17.md.
//

#include "TrajectoryResult.h"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class HmWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "HmWelfordTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<HmWelfordTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp, Trajectory& traj) override;

    size_t NumFrames() const { return n_frames_; }

    int WriteFeatures(const TrajectoryProtein& tp,
                      const std::string& output_dir) const override;

    void WriteH5Group(const TrajectoryProtein& tp,
                      HighFive::File& file) const override;

private:
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
