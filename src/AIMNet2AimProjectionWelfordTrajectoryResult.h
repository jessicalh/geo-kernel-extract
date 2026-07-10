#pragma once
//
// H5-only per-atom Welford rollup of the stored 32-component AIMNet2 aim
// projection. AIMNet2Result::Compute owns projection arithmetic and writes
// ConformationAtom::aimnet2_aim_projection exactly once. This trajectory
// result is deliberately a pure reader of that field.

#include "TrajectoryResult.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class AIMNet2AimProjectionWelfordTrajectoryResult final
    : public TrajectoryResult {
public:
    std::string Name() const override {
        return "AIMNet2AimProjectionWelfordTrajectoryResult";
    }

    // The source-presence gate is intentional: defensive custom trajectory
    // configurations can attach this reducer without AIMNet2Result, in which
    // case the frame is masked and no Welford state is updated.
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<AIMNet2AimProjectionWelfordTrajectoryResult>
    Create(const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp, Trajectory& traj) override;

    // H5-only: inherit TrajectoryResult::WriteFeatures (returns 0).
    void WriteH5Group(const TrajectoryProtein& tp,
                      HighFive::File& file) const override;

    std::size_t NumFrames() const { return n_frames_; }
    std::size_t SourceAttachedCount() const {
        return source_attached_count_;
    }

private:
    std::vector<std::uint64_t> frame_indices_;
    std::vector<double> frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    std::size_t n_frames_ = 0;
    std::size_t source_attached_count_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
