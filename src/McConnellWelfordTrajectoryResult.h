#pragma once
//
// McConnellWelfordTrajectoryResult: running mean / variance / min / max
// of McConnell shielding tensor per atom, accumulated across all frames
// of a trajectory. AV-pattern clone of BsWelfordTrajectoryResult.
//
// Source: `ConformationAtom::mc_shielding_contribution` (units = Å⁻³,
// full McConnell asymmetric non-traceless three-term form — see
// PATTERNS.md Lesson 19; T0 = (3cos²θ-1)/r³ is meaningful). McConnellResult
// is unconditionally attached in `PerFrameExtractionSet`, so the dep
// is enforced via `Dependencies()`.
//

#include "TrajectoryResult.h"

#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class McConnellWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "McConnellWelfordTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<McConnellWelfordTrajectoryResult> Create(
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
    size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
