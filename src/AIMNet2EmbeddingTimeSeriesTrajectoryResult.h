#pragma once
//
// AIMNet2EmbeddingTimeSeriesTrajectoryResult: per-atom per-frame time
// series of the 256-dim AIMNet2 'aim' embedding tensor
// (ConformationAtom::aimnet2_aim, std::array<float, AIMNET2_AIM_DIMS>).
// FO dense-buffer pattern; clone of AIMNet2ChargeTimeSeriesTrajectoryResult
// with the embedding axis added.
//
// AIMNet2Result is in PerFrameExtractionSet (RunConfiguration.cpp); when
// the AIMNet2 model is not Session-loaded, OperationRunner aborts the
// run before any TR Compute fires — no per-frame conditional gate at
// the TR layer (always_attached policy).
//
// Emission:
//
//   /trajectory/aimnet2_embedding_time_series/
//     embedding      (N, T, 256) float32 — AIMNet2 'aim' tensor
//     frame_indices  (T,)        uint64
//     frame_times    (T,)        float64 — ps
//     attrs:
//       result_name             = "AIMNet2EmbeddingTimeSeriesTrajectoryResult"
//       n_atoms, n_frames, finalized
//       embedding_dim           = 256
//       units                   = "dimensionless"
//       source                  = "AIMNet2Result.aimnet2_aim (AIMNET2_AIM_DIMS=256)"
//       source_attached_policy  = "always_attached"
//       optional_large          = true
//
// Storage discipline: float32 native (per `feedback_embedding_float32`).
// At fleet scale the dataset is ~3.8 GB/protein uncompressed; chunked
// (N, frame_chunk=64, 256) for movie-target seek. Compression is the
// HighFive default for now; revisit with Blosc2 if storage budget
// requires.
//

#include "TrajectoryResult.h"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class AIMNet2EmbeddingTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "AIMNet2EmbeddingTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<AIMNet2EmbeddingTimeSeriesTrajectoryResult>
    Create(const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp, Trajectory& traj) override;

    void WriteH5Group(const TrajectoryProtein& tp,
                      HighFive::File& file) const override;

    std::size_t NumFrames() const { return n_frames_; }

private:
    // (N, T, 256) — outer per-atom, inner per-frame array slot.
    std::vector<std::vector<std::array<float, 256>>> per_atom_embedding_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double>      frame_times_;
    std::size_t n_frames_  = 0;
    bool        finalized_ = false;
};

}  // namespace nmr
