#pragma once
//
// AIMNet2PolarisabilityTimeSeriesTrajectoryResult: per-atom per-frame
// time series of the AIMNet2 charge-polarisation gradient. Two
// emissions per atom per frame:
//
//   - vector (Vec3): aimnet2_polarisability_vector, ConformationAtom field.
//                    Gradient of L = Σ_j q_j² (AIMNet2 charges, in units
//                    e²) with respect to atomic coordinates (autograd
//                    backward through AIMNet2 charge head). Units e²/Å.
//   - scalar (double): aimnet2_polarisability_scalar, L2 norm of the
//                    vector — emit both rather than recompute downstream.
//                    Units e²/Å.
//
// FO dense-buffer pattern; pairs with AIMNet2EmbeddingTS (TR #1 of the
// AIMNet2 fleet trio). AIMNet2PolarisabilityResult is now registered as
// a required ConformationResult in the trajectory PerFrameExtractionSet
// (RunConfiguration.cpp); OperationRunner aborts the run if the model
// cannot Compute it on any frame — no per-frame conditional gate at the
// TR layer (always_attached policy).
//
// Emission:
//
//   /trajectory/aimnet2_polarisability_time_series/
//     polarisability_vector       (N, T, 3) float64 — e²/Å
//     polarisability_scalar       (N, T)    float64 — e²/Å (L2 norm)
//     frame_indices               (T,)      uint64
//     frame_times                 (T,)      float64 — ps
//     source_attached_per_frame   (T,)      uint8   — always 1 (canonical
//                                                     SDK contract for the
//                                                     "Conditional-attach TR
//                                                     discipline" subsection
//                                                     of OBJECT_MODEL.md)
//     attrs:
//       result_name             = "AIMNet2PolarisabilityTimeSeriesTrajectoryResult"
//       n_atoms, n_frames, finalized
//       units_vector            = "e^2/Angstrom"
//       units_scalar            = "e^2/Angstrom"
//       irrep_layout_vector     = "1o"  (Vec3 odd parity, rank-1)
//       parity_vector           = "1o"
//       irrep_layout_scalar     = "T0"  (rank-0 invariant)
//       parity_scalar           = "0e"
//       source                  = "AIMNet2PolarisabilityResult.{vector,scalar}"
//       source_attached_policy  = "always_attached"
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

class AIMNet2PolarisabilityTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "AIMNet2PolarisabilityTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<AIMNet2PolarisabilityTimeSeriesTrajectoryResult>
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
