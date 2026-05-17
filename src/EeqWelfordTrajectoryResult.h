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
    size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
