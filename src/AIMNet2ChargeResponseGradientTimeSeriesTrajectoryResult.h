#pragma once
//
// AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult: per-atom per-frame
// time series of the AIMNet2 charge-response gradient. Two
// emissions per atom per frame:
//
//   - vector (Vec3): aimnet2_charge_response_gradient_vector,
//                    ConformationAtom field.
//                    Gradient of L = Σ_j q_j² (AIMNet2 charges, in units
//                    e²) with respect to atomic coordinates (autograd
//                    backward through AIMNet2 charge head). Units e²/Å.
//   - scalar (double): aimnet2_charge_response_gradient_scalar, L2 norm
//                    of the vector — emit both rather than recompute
//                    downstream. Units e²/Å.
//
// Per-atom buffered writer; pairs with AIMNet2EmbeddingTS. AIMNet2ChargeResponseGradientResult is registered as
// a required ConformationResult in the trajectory PerFrameExtractionSet
// (RunConfiguration.cpp); OperationRunner aborts the run if the model
// cannot Compute it on any frame. Compute still gates on
// HasResult<AIMNet2ChargeResponseGradientResult>() and emits NaN-fill +
// mask=0 if a custom configuration omits the required source.
//
// Emission:
//
//   /trajectory/aimnet2_charge_response_gradient_time_series/
//     charge_response_gradient_vector       (N, T, 3) float64 — e²/Å
//     charge_response_gradient_scalar       (N, T)    float64 — e²/Å (L2 norm)
//     frame_indices               (T,)      uint64
//     frame_times                 (T,)      float64 — ps
//     source_attached_per_frame   (T,)      uint8   — per-frame source mask
//     attrs:
//       result_name             = "AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult"
//       n_atoms, n_frames, finalized
//       units_vector            = "e^2/Å"
//       units_scalar            = "e^2/Å"
//       irrep_layout_vector     = "x,y,z"      (Cartesian component order)
//       normalization_vector    = "cartesian"
//       parity_vector           = "1o"         (odd parity, rank-1)
//       irrep_layout_scalar     = "T0"         (rank-0 invariant)
//       parity_scalar           = "0e"
//       source                  describes the AIMNet2ChargeResponseGradientResult vector/scalar fields
//       source_attached_policy  = "always_attached" — but Compute's
//                                  HasResult gate emits NaN-fill +
//                                  source_attached_per_frame=0 on
//                                  absent frames (codex review
//                                  2026-05-20; "absent, not faked").
//
// Per `feedback_methods_accumulate`, emit both vector AND scalar even
// though scalar is derivable from vector — both are downstream-useful
// and the cost is trivial.
//

#include "TrajectoryResult.h"
#include "Types.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult>
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
    std::vector<std::vector<Vec3>>   per_atom_vector_;  // (N, T) Vec3
    std::vector<std::vector<double>> per_atom_scalar_;  // (N, T)
    std::vector<std::size_t>         frame_indices_;
    std::vector<double>              frame_times_;
    // Per-frame source-attached mask (canonical gate; codex review
    // 2026-05-20). Normally all-1 under the always-attached policy.
    std::vector<std::uint8_t>        source_attached_per_frame_;
    std::size_t n_frames_  = 0;
    bool        finalized_ = false;
};

}  // namespace nmr
