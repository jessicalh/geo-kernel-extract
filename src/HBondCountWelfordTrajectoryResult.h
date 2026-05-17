#pragma once
//
// HBondCountWelfordTrajectoryResult: running mean / variance / min / max
// of per-atom H-bond count (`hbond_count_within_3_5A`, integer source,
// pairs within 3.5 Å), accumulated across all frames.
//
// The running mean is fractional — it represents the trajectory-averaged
// H-bond occupancy frequency at each atom (e.g., 1.20 means this atom
// averages 1.2 H-bonds across the frames). AV-pattern scalar Welford.
//
// HBondResult is unconditionally attached in `PerFrameExtractionSet`.
//

#include "TrajectoryResult.h"

#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class HBondCountWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "HBondCountWelfordTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<HBondCountWelfordTrajectoryResult> Create(
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
    std::vector<double> prev_value_;
    std::vector<bool> prev_valid_;
    size_t n_frames_ = 0;
    bool finalized_ = false;

    // Cadence metadata — mean Δt between captured frames (ps). Derived
    // at Finalize from traj.FrameTimes(). Emitted as H5 attribute so
    // downstream can convert delta-per-stride to dx/dt in physical
    // units. 0.0 before Finalize. See PATTERNS Lesson 25.
    double mean_dt_ps_ = 0.0;
};

}  // namespace nmr
