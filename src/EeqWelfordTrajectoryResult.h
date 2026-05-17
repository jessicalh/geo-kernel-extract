#pragma once
//
// EeqWelfordTrajectoryResult: running mean / variance / min / max of
// the Eeq geometry-dependent atomic charge per atom, accumulated
// across all frames. AV-pattern scalar Welford — single-channel source
// (`ConformationAtom::eeq_charge`, units = e), no T0/T2 distinction.
//
// EeqResult is unconditionally attached in `PerFrameExtractionSet`, so
// the dep is enforced via `Dependencies()`.
//

#include "TrajectoryResult.h"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class EeqWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "EeqWelfordTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<EeqWelfordTrajectoryResult> Create(
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

    // Per-atom previous-frame time (ps), used to compute cadence-
    // normalized dxdt = (x_now - x_prev) / (time_now - time_prev).
    // Valid for atom i only when prev_valid_[i] is true.
    std::vector<double> prev_time_;

    size_t n_frames_ = 0;
    bool finalized_ = false;

    // Cadence metadata — mean Δt between captured frames (ps). Derived
    // at Finalize from traj.FrameTimes(). Emitted as H5 attribute so
    // downstream can convert delta-per-stride to dx/dt in physical
    // units. 0.0 before Finalize. See PATTERNS Lesson 25.
    double mean_dt_ps_ = 0.0;

    // Frame-index span this rollup covered. [first, last] from
    // traj.FrameIndices() at Finalize. Emitted as H5 group attribute.
    std::array<std::size_t, 2> frame_index_range_ = {0, 0};
};

}  // namespace nmr
