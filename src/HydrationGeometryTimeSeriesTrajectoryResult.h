#pragma once
//
// HydrationGeometryTimeSeriesTrajectoryResult: per-atom per-frame timeline
// of the SASA-normal water polarisation features from HydrationGeometryResult.
// Six channels per atom per frame:
//
//   water_dipole_vector       Vec3      net first-shell water dipole (sum of d_i)
//   water_surface_normal      Vec3      SASA outward normal (unit vector)
//   sasa_first_shell_count    uint32    first-shell water O count (statistical weight)
//   sasa_half_shell_asymmetry scalar    fraction of waters on solvent-exposed side
//   sasa_dipole_alignment     scalar    cos(net dipole, surface normal)
//   sasa_dipole_coherence     scalar    |Σ d_i| / n
//
// The alignment / coherence / asymmetry trio IS the polarisation signal —
// calibration weight on these channels is what we want to learn.
//
// Dependencies: HydrationGeometryResult — REQUIRED by PerFrameExtractionSet
// but conditionally attached by OperationRunner: if `opts.solvent` is null
// or `opts.solvent->Empty()` (no solvent environment loaded), the source
// ConformationResult is silently skipped. Follows "absent, not faked":
//   - Per-frame `conf.HasResult<HydrationGeometryResult>()` check in Compute
//   - `source_attached_per_frame` mask emitted as H5 provenance
//   - NaN-fill atom-axis data for source-absent frames
//   - WriteH5Group skips entire group emission when source attached zero times
//
// Emission:
//
//   /trajectory/hydration_geometry_time_series/
//     dipole_vector             (N, T, 3)  float64  Debye-like (unnormalised)
//     surface_normal            (N, T, 3)  float64  unit vector
//     first_shell_count         (N, T)     uint32   shell-occupancy count
//     half_shell_asymmetry      (N, T)     float64  fraction
//     dipole_alignment          (N, T)     float64  cos angle
//     dipole_coherence          (N, T)     float64  order parameter
//     frame_indices             (T,)       uint64
//     frame_times               (T,)       float64  ps
//     source_attached_per_frame (T,)       uint8    1=attached, 0=absent
//     attrs:
//       result_name, n_atoms, n_frames, source_attached_count, finalized
//       dipole_vector_layout       = "x,y,z"
//       dipole_vector_parity       = "1o"   (polar vector)
//       surface_normal_layout      = "x,y,z"
//       surface_normal_parity      = "1o"
//       polarisation_signal_channels = "alignment,coherence,asymmetry"
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

class HydrationGeometryTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "HydrationGeometryTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<HydrationGeometryTimeSeriesTrajectoryResult> Create(
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

    // Test-only bypass for the source-attached gate.
    void ForceSourcePresentForTesting() {
        force_source_present_for_testing_ = true;
    }

private:
    std::vector<std::vector<Vec3>>          dipole_vector_;
    std::vector<std::vector<Vec3>>          surface_normal_;
    std::vector<std::vector<std::uint32_t>> first_shell_count_;
    std::vector<std::vector<double>>        half_shell_asymmetry_;
    std::vector<std::vector<double>>        dipole_alignment_;
    std::vector<std::vector<double>>        dipole_coherence_;
    std::vector<std::size_t>                frame_indices_;
    std::vector<double>                     frame_times_;
    std::vector<std::uint8_t>               source_attached_per_frame_;
    bool                                    force_source_present_for_testing_ = false;
    std::size_t                             n_frames_  = 0;
    bool                                    finalized_ = false;
};

}  // namespace nmr
