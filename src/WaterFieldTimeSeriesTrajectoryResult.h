#pragma once
//
// WaterFieldTimeSeriesTrajectoryResult: per-atom per-frame timeline of
// the explicit-water E-field + EFG kernel from WaterFieldResult. Six
// channels per atom per frame:
//
//   water_efield           Vec3 (V/Å)            — total field (cutoff sphere)
//   water_efield_first     Vec3 (V/Å)            — first-shell-only (< 3.5 Å)
//   water_efg_spherical    SphericalTensor       — total EFG (T0+T1+T2)
//   water_efg_first_spherical SphericalTensor    — first-shell-only EFG
//   water_n_first          int                   — water O count, 0..3.5 Å
//   water_n_second         int                   — water O count, 3.5..5.5 Å
//
// Storage: separate per-atom growing buffers, six in parallel. H5 emits
// (N, T, ...) datasets per channel. No DenseBuffer adoption — water-field
// time-series is a leaf consumer for downstream calibration / ML; no other
// TR reads per-frame water E-field channels.
//
// Use case: the explicit-water side of the per-atom electrostatic
// environment time series. Pairs with `ApbsEfieldTimeSeriesTrajectoryResult`
// (continuum dielectric) so the calibration step can read both and
// learn which approximation best matches reference shifts.
//
// Export-everything-upstream (PATTERNS Lesson 25): emit both total and
// first-shell-only variants for each E-field / EFG channel; emit shell
// occupancy counts separately. Removed columns are forever.
//
// Emission:
//
//   /trajectory/water_field_time_series/
//     efield                    (N, T, 3)  float64  V/Angstrom
//     efield_first              (N, T, 3)  float64  V/Angstrom
//     efg                       (N, T, 5)  float64  V/Angstrom^2  (T2_m-2..T2_m+2)
//     efg_first                 (N, T, 5)  float64  V/Angstrom^2
//     n_first                   (N, T)     uint32   shell-occupancy count
//     n_second                  (N, T)     uint32   shell-occupancy count
//     frame_indices             (T,)       uint64
//     frame_times               (T,)       float64  ps
//     source_attached_per_frame (T,)       uint8    provenance mask
//     attrs:
//       result_name, n_atoms, n_frames, source_attached_count, finalized
//       efield_layout       = "x,y,z" / efield_parity = "1o"
//       efg_irrep_layout    = "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
//       efg_parity          = "2e" (water EFG: T0=0 trace + T1=0 antisym both structural)
//       efg_t0_structural_zero, efg_t1_structural_zero  (true)
//       efg_units           = "V/Angstrom^2"
//       count_units         = "dimensionless"
//
// Dependencies: WaterFieldResult — REQUIRED by PerFrameExtractionSet but
// conditionally attached by OperationRunner: if `opts.solvent` is null
// or `opts.solvent->Empty()` (no solvent environment loaded), the source
// ConformationResult is silently skipped. Follows "absent, not faked":
//   - Per-frame `conf.HasResult<WaterFieldResult>()` check in Compute
//   - `source_attached_per_frame` mask emitted as H5 provenance
//   - NaN-fill atom-axis data for source-absent frames
//   - WriteH5Group skips entire group emission when source attached zero times
//   - Protein-only extractions (no water in TPR) ⇒ group absent, not
//     a hidden all-zero H5 group that mimics solvated runs.
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

class WaterFieldTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "WaterFieldTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<WaterFieldTimeSeriesTrajectoryResult> Create(
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

    // Test-only bypass for the source-attached gate. Production code
    // never calls this — only synthetic tests that don't run through
    // OperationRunner.
    void ForceSourcePresentForTesting() {
        force_source_present_for_testing_ = true;
    }

private:
    std::vector<std::vector<Vec3>>             efield_;
    std::vector<std::vector<Vec3>>             efield_first_;
    std::vector<std::vector<SphericalTensor>>  efg_;
    std::vector<std::vector<SphericalTensor>>  efg_first_;
    std::vector<std::vector<std::uint32_t>>    n_first_;
    std::vector<std::vector<std::uint32_t>>    n_second_;
    std::vector<std::size_t>                   frame_indices_;
    std::vector<double>                        frame_times_;
    std::vector<std::uint8_t>                  source_attached_per_frame_;
    bool                                       force_source_present_for_testing_ = false;
    std::size_t                                n_frames_ = 0;
    bool                                       finalized_ = false;
};

}  // namespace nmr
