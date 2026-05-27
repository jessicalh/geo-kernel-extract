#pragma once
//
// AIMNet2EmbeddingTimeSeriesTrajectoryResult: per-atom per-frame time
// series of the 256-dim AIMNet2 'aim' embedding tensor
// (ConformationAtom::aimnet2_aim, std::array<float, AIMNET2_AIM_DIMS>).
// Per-atom buffered writer with the embedding axis added.
//
// AIMNet2Result is in PerFrameExtractionSet (RunConfiguration.cpp); when
// the AIMNet2 model is not Session-loaded, OperationRunner aborts the
// run before any TR Compute fires. Compute still gates on
// HasResult<AIMNet2Result>() and emits NaN-fill + mask=0 if a custom
// configuration omits the required source.
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
//       source                  describes AIMNet2Result.aimnet2_aim and AIMNET2_AIM_DIMS
//       source_attached_policy  = "always_attached" — but Compute's
//                                  HasResult<AIMNet2Result>() gate
//                                  emits NaN-fill + source_attached
//                                  _per_frame=0 on absent frames
//                                  (codex review 2026-05-20; "absent,
//                                  not faked").
//       optional_large          = true
//
// Storage discipline: float32 native (per `feedback_embedding_float32`).
// At fleet scale the dataset is ~3.8 GB/protein uncompressed; chunked
// (1, frame_chunk=64, 256) for per-atom writes and movie-target seek.
// No explicit compression is configured here.
//

#include "TrajectoryResult.h"

#include <array>
#include <cstddef>
#include <cstdint>
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
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    // Per-frame source-attached mask. Always-attached policy means the
    // mask is normally all-1, but a custom config that wires the TR
    // without RequireConformationResult'ing AIMNet2Result would land
    // mask=0 for those frames (codex review 2026-05-20).
    std::vector<std::uint8_t> source_attached_per_frame_;
    std::size_t n_frames_  = 0;
    bool        finalized_ = false;
};

}  // namespace nmr
