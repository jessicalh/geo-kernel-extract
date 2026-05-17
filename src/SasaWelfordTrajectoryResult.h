#pragma once
//
// SasaWelfordTrajectoryResult: running mean / variance / min / max of
// Shrake-Rupley solvent-accessible surface area per atom (`atom_sasa`,
// units = Å²), accumulated across all frames. AV-pattern scalar Welford.
//
// Pairs with the landed FO-pattern `SasaTimeSeriesTrajectoryResult`:
// the time-series captures every frame; this Welford captures running
// stats. Both consume the same source field.
//
// SasaResult is unconditionally attached in `PerFrameExtractionSet`.
//

#include "TrajectoryResult.h"

#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class SasaWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "SasaWelfordTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<SasaWelfordTrajectoryResult> Create(
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
};

}  // namespace nmr
