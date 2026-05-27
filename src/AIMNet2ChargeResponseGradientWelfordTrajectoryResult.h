#pragma once
//
// AIMNet2ChargeResponseGradientWelfordTrajectoryResult: AV companion to
// AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult. Per-atom
// Welford rollup of the four
// channels — three Cartesian Vec3 components + one scalar L2 norm
// of the charge-response gradient ∂L/∂r_i where L = Σ_j q_j².
//
// AV (Always-valid) pattern: each Compute updates the per-atom
// WelfordMoments on
// TrajectoryAtom::aimnet2_charge_response_gradient_welford; Finalize
// derives std and WriteH5Group writes the means/variances to H5. Source-attached gate
// matches the TS pair: Compute checks
// HasResult<AIMNet2ChargeResponseGradientResult>() and skips the update on
// absent frames (records mask=0 for SDK provenance). Production
// trajectory pipelines require this source so all frames should attach.
//
// Emission (/trajectory/aimnet2_charge_response_gradient_welford/) — full
// canonical Welford row per sibling TR convention (HydrationGeometry /
// BsWelford / HmWelford); WelfordFinalize derives std + NaN-fills n=0.
//   vector_mean         (N, 3) float64 — per-component mean (e²/Å)
//   vector_std          (N, 3) float64 — per-component std
//                                         (= sqrt(m2/(n-1)) for n≥2;
//                                          0 for n=1; NaN for n=0)
//   vector_m2           (N, 3) float64 — per-component Welford M2
//   vector_min          (N, 3) float64 — per-component running min
//   vector_max          (N, 3) float64 — per-component running max
//   vector_min_frame    (N, 3) uint64  — frame index at which each
//                                         component min was attained
//   vector_max_frame    (N, 3) uint64  — frame index at which each
//                                         component max was attained
//   scalar_mean         (N,)   float64 — L2-norm mean (e²/Å)
//   scalar_std          (N,)   float64 — L2-norm std
//   scalar_m2           (N,)   float64 — L2-norm M2
//   scalar_min          (N,)   float64 — L2-norm running min
//   scalar_max          (N,)   float64 — L2-norm running max
//   scalar_min_frame    (N,)   uint64  — frame at which scalar_min was hit
//   scalar_max_frame    (N,)   uint64  — frame at which scalar_max was hit
//   n_per_atom          (N,)   uint64  — shared sample count across
//                                         the 4 channels
//   source_attached_per_frame (T,) uint8 — per-frame attach mask
//   frame_indices            (T,) uint64
//   frame_times              (T,) float64
//
// Attrs:
//   result_name             = "AIMNet2ChargeResponseGradientWelfordTrajectoryResult"
//   n_atoms, n_frames, source_attached_count, finalized
//   units_vector            = "e^2/Å"
//   units_scalar            = "e^2/Å"
//   irrep_layout_vector     = "x,y,z"
//   normalization_vector    = "cartesian"
//   parity_vector           = "1o"
//   irrep_layout_scalar     = "T0"
//   parity_scalar           = "0e"
//   source                  describes the AIMNet2ChargeResponseGradientResult vector/scalar fields
//   source_attached_policy  = "always_attached" with HasResult gate
//
// Minimum-viable design (no delta variants in v0; mean/std/m2/min/max
// canonical row only). Delta-variant pattern from
// HydrationGeometryWelfordTrajectoryResult is available if a
// calibration finding requests dx/dt or rms_delta later.
//
// PHYSICS NOTE: The emitted quantity is ∂(Σ_j q_j²)/∂r_i (gradient
// of a sum-of-squared-AIMNet2-charges scalar with respect to atomic
// coordinates), NOT a Buckingham α tensor (α_ab = ∂μ_a/∂E_b, the
// dipole-response-to-field quantity that conventionally appears in
// NMR shielding theory). The L = Σ q² objective is a
// computationally-cheap, autograd-friendly proxy for per-atom charge
// sensitivity; physical-observable connection to NMR shielding is
// exploratory and calibration-ridge decides if signal carries beyond
// the AIMNet2 charge channel.
//
// NUMERICAL CAVEAT: the per-atom Welford std includes variation from
// the AIMNet2 autograd backward pass as well as trajectory motion.
// Treat very small std values as numerical-plus-physical variation,
// especially for tightly constrained atoms.
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class AIMNet2ChargeResponseGradientWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "AIMNet2ChargeResponseGradientWelfordTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<AIMNet2ChargeResponseGradientWelfordTrajectoryResult>
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
