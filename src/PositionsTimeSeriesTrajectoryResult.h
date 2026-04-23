#pragma once
//
// PositionsTimeSeriesTrajectoryResult: per-atom per-frame position
// time series. The Finalize-only dense-buffer worked example,
// parallel to BsWelfordTrajectoryResult which is the always-valid-
// mid-stream rollup worked example.
//
// Every future time-series TrajectoryResult follows this shape:
//
//   Compute   — append this frame's per-atom data to internal
//               per-atom growing buffers. No TrajectoryAtom field
//               writes; the buffers are owned internally until
//               Finalize.
//   Finalize  — flatten internal buffers into a DenseBuffer<T>,
//               transfer ownership to TrajectoryProtein via
//               AdoptDenseBuffer<T>, record the frame_indices /
//               frame_times so WriteH5Group has its provenance.
//   WriteH5Group — retrieve the buffer via tp.GetDenseBuffer<T>,
//               emit /trajectory/positions/{xyz,frame_indices,
//               frame_times} with shape attributes.
//
// Reads positions directly from ProteinConformation per frame; no
// per-frame ConformationResult dependency.
//

#include "DenseBuffer.h"
#include "TrajectoryResult.h"
#include "Types.h"

#include <cstddef>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class PositionsTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "PositionsTimeSeriesTrajectoryResult";
    }

    // No ConformationResult dependency — positions are always present
    // on every ProteinConformation.
    std::vector<std::type_index> Dependencies() const override {
        return {};
    }

    // Factory. Allocates the per-atom growing buffers sized to the
    // TrajectoryProtein's atom count.
    static std::unique_ptr<PositionsTimeSeriesTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    // Transfer internal per-atom buffers into a contiguous
    // DenseBuffer<Vec3> owned by TrajectoryProtein.
    void Finalize(TrajectoryProtein& tp, Trajectory& traj) override;

    void WriteH5Group(const TrajectoryProtein& tp,
                      HighFive::File& file) const override;

    std::size_t NumFrames() const { return n_frames_; }

private:
    // Per-atom growing buffers: one vector per atom, Vec3 per frame.
    // At Finalize, flattened into the atom-major DenseBuffer layout.
    std::vector<std::vector<Vec3>> per_atom_positions_;

    // Frame-axis metadata — populated per Compute so WriteH5Group can
    // emit the frame list alongside the xyz dataset.
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;

    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
