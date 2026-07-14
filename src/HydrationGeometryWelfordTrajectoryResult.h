#pragma once
//
// HydrationGeometryWelfordTrajectoryResult: per-atom Welford rollup of the
// SASA-normal water polarisation features from HydrationGeometryResult.
// Clones the WaterFieldWelford shape: per-component Welford on Vec3 channels,
// scalar Welford + signed/abs/squared/dxdt delta variants on the
// polarisation scalars.  Dipole magnitude/coherence/order and shell count
// are exact rotation invariants; alignment and asymmetry are continuum
// invariants that inherit the live finite SASA-stencil approximation.
//
// Channels (mirroring HydrationGeometryTimeSeries source fields):
//
//   dipole_vector         Vec3 (e·Å, raw sum)          per-component[3] + magnitude
//   surface_normal        Vec3 (unit vector)           per-component[3]
//   half_shell_asymmetry  scalar (fraction)            full Welford + delta variants
//   dipole_alignment      scalar (cos angle)           full Welford + delta variants
//   dipole_coherence      scalar (e·Å, legacy `|Σd|/n`) full Welford + delta variants
//   dipole_order_parameter scalar ([0,1], `|Σd|/Σ|d|`) full Welford + delta variants
//   shell_count           int (dimensionless)          full Welford + delta variants
//
// R6 codex 2026-05-18: dipole_coherence is NOT a dimensionless order
// parameter despite the name — source formula `|Σd_i| / n_shell` has
// e·Å units. A true coherence (dimensionless [0,1]) would divide by
// `Σ |d_i|` instead. Consumers can post-process if needed.
//
// Delta variants on: half_shell_asymmetry, dipole_alignment,
// dipole_coherence, dipole_order_parameter, shell_count.
//
// Emission:
//
//   /trajectory/hydration_geometry_welford/
//     dipole_vector_{x,y,z}_{mean,m2,std,min,max,min_frame,max_frame}  (N,)
//     dipole_magnitude_{...}                                            (N,)
//     surface_normal_{x,y,z}_{...}                                      (N,)
//     half_shell_asymmetry_{...}, dipole_alignment_{...},
//     dipole_coherence_{...}, shell_count_{...}                         (N,)
//     <scalar>_{delta,abs_delta,delta_squared,dxdt}_{...}              for the 4 scalars
//     <scalar>_rms_delta                                                (N,) Finalize-derived
//     n_frames_per_atom, delta_n_per_atom, dxdt_n_per_atom              (N,)
//     source_attached_per_frame                                         (T,) uint8
//     attrs: result_name, n_frames, source_attached_count, finalized,
//            ddof, mean_dt_ps (attached-subset), frame_index_range
//            (attached-subset), irrep_layout_dipole, irrep_layout_normal
//            dipole_vector_coordinate_frame = "conformation_cartesian_xyz"
//            dipole_vector_parity = "1o"
//            dipole_vector_transformation = "polar vector: v'=R v"
//            surface_normal_coordinate_frame = "conformation_cartesian_xyz"
//            surface_normal_parity = "mixed"
//            surface_normal_transformation = continuum-polar with the live
//                finite lab-fixed Fibonacci/no-exact-O(3) qualifier
//            scalar_statistic_transformation = all 72 normal-derived scalar
//                statistic paths inherit that finite-grid qualifier
//            directional_metadata_scope = exactly 114 paths: 42 vector
//                component statistics + 72 normal-derived scalar statistics
//
// The three per-source channels are independently accumulated Cartesian
// x/y/z components, not SphericalTensor T1 components.  Reassembling the
// three dipole component means yields an exact polar vector.  Reassembling
// the three surface-normal means yields a continuum-polar vector that still
// inherits SasaResult's finite lab-fixed Fibonacci approximation.  No
// componentwise m2/std/min/max/extremum dataset is a covariant vector.  Every
// Welford/delta descendant of half_shell_asymmetry and dipole_alignment also
// inherits the surface normal's finite-grid approximation; only the dipole
// magnitude/coherence/net-dipole/order and count/provenance families are
// exact rotation invariants outside the 114-path directional scope.
//
// Dependencies: HydrationGeometryResult — REQUIRED by PerFrameExtractionSet
// but conditionally attached by OperationRunner when solvent is loaded.
// Follows "absent, not faked" — frames failing the source-present check skip
// Welford updates, invalidate prev_* caches for clean delta restart on the
// next attached frame, and WriteH5Group skips the group when
// source_attached_count == 0.
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

class HydrationGeometryWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "HydrationGeometryWelfordTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<HydrationGeometryWelfordTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp, Trajectory& traj) override;

    void WriteH5Group(const TrajectoryProtein& tp,
                      HighFive::File& file) const override;

    std::size_t NumFrames() const { return n_frames_; }

    void ForceSourcePresentForTesting() {
        force_source_present_for_testing_ = true;
    }

private:
    // Per-atom previous-frame caches for the delta trackers
    std::vector<double> prev_half_shell_asymmetry_;
    std::vector<double> prev_dipole_alignment_;
    std::vector<double> prev_dipole_coherence_;
    std::vector<double> prev_dipole_order_parameter_;
    std::vector<double> prev_shell_count_;
    std::vector<bool>   prev_valid_;
    std::vector<double> prev_time_;

    std::vector<std::uint8_t> source_attached_per_frame_;
    bool                      force_source_present_for_testing_ = false;

    std::size_t              n_frames_  = 0;
    bool                     finalized_ = false;
    double                   mean_dt_ps_ = 0.0;
    std::array<std::size_t, 2> frame_index_range_ = {0, 0};
};

}  // namespace nmr
