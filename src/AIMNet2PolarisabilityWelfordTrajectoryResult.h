#pragma once
//
// AIMNet2PolarisabilityWelfordTrajectoryResult: AV companion to
// AIMNet2PolarisabilityTimeSeriesTrajectoryResult (TR #3 of the
// AIMNet2 fleet trio). Per-atom Welford rollup of the four
// channels — three Cartesian Vec3 components + one scalar L2 norm
// of the polarisability gradient ∂L/∂r_i where L = Σ_j q_j².
//
// AV (Always-valid) pattern: each Compute updates the per-atom
// WelfordMoments on TrajectoryAtom::aimnet2_polarisability_welford;
// Finalize writes the means/variances to H5. Source-attached gate
// matches the TS pair: Compute checks
// HasResult<AIMNet2PolarisabilityResult>() and skips the update on
// absent frames (records mask=0 for SDK provenance). Production
// trajectory pipelines RequireConformationResult'd it at line 167
// of RunConfiguration.cpp so all frames should attach.
//
// Emission (/trajectory/aimnet2_polarisability_welford/):
//   vector_mean       (N, 3) float64 — per-component mean (e²/Å)
//   vector_m2         (N, 3) float64 — per-component sum-of-squared
//                                       deviations (Welford M2)
//   vector_n          (N, 3) uint64  — per-component sample count
//   scalar_mean       (N,)   float64 — L2-norm mean (e²/Å)
//   scalar_m2         (N,)   float64 — L2-norm M2
//   scalar_n          (N,)   uint64  — L2-norm sample count
//   source_attached_per_frame (T,) uint8 — per-frame attach mask
//   frame_indices            (T,) uint64
//   frame_times              (T,) float64
//
// Attrs:
//   result_name             = "AIMNet2PolarisabilityWelfordTrajectoryResult"
//   n_atoms, n_frames, source_attached_count, finalized
//   units_vector            = "e^2/Angstrom"
//   units_scalar            = "e^2/Angstrom"
//   irrep_layout_vector     = "x,y,z"
//   normalization_vector    = "cartesian"
//   parity_vector           = "1o"
//   irrep_layout_scalar     = "T0"
//   parity_scalar           = "0e"
//   source                  = "AIMNet2PolarisabilityResult.{vector,scalar}"
//   source_attached_policy  = "always_attached" with HasResult gate
//
// Minimum-viable design (no delta variants in v0). Delta-variant
// pattern from HydrationGeometryWelfordTrajectoryResult is available
// if a calibration finding requests dx/dt or rms_delta later.
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class AIMNet2PolarisabilityWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "AIMNet2PolarisabilityWelfordTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<AIMNet2PolarisabilityWelfordTrajectoryResult>
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
    std::size_t SourceAttachedCount() const { return source_attached_count_; }

private:
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    std::size_t               n_frames_              = 0;
    std::size_t               source_attached_count_ = 0;
    bool                      finalized_             = false;
};

}  // namespace nmr
