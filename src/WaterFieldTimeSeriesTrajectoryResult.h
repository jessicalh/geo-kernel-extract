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
//     efg                       (N, T, 9)  float64  V/Angstrom^2  (T0..T2_+2)
//     efg_first                 (N, T, 9)  float64  V/Angstrom^2
//     n_first                   (N, T)     uint32   shell-occupancy count
//     n_second                  (N, T)     uint32   shell-occupancy count
//     frame_indices             (T,)       uint64
//     frame_times               (T,)       float64  ps
//     attrs:
//       result_name, n_atoms, n_frames, finalized
//       efield_layout       = "x,y,z"
//       efield_units        = "V/Angstrom"
//       efg_irrep_layout    = "T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
//       efg_units           = "V/Angstrom^2"
//       count_units         = "dimensionless"
//
// Dependencies: WaterFieldResult (unconditionally attached in
// PerFrameExtractionSet — no source-attached gate).
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

private:
    std::vector<std::vector<Vec3>>             efield_;
    std::vector<std::vector<Vec3>>             efield_first_;
    std::vector<std::vector<SphericalTensor>>  efg_;
    std::vector<std::vector<SphericalTensor>>  efg_first_;
    std::vector<std::vector<std::uint32_t>>    n_first_;
    std::vector<std::vector<std::uint32_t>>    n_second_;
    std::vector<std::size_t>                   frame_indices_;
    std::vector<double>                        frame_times_;
    std::size_t                                n_frames_ = 0;
    bool                                       finalized_ = false;
};

}  // namespace nmr
