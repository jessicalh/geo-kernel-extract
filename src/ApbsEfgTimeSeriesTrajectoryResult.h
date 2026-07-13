#pragma once
//
// ApbsEfgTimeSeriesTrajectoryResult: per-atom per-frame time series of
// the canonical APBS reaction electric field gradient
// (ConformationAtom::apbs_efg_spherical, T2 SphericalTensor, V/Å²).
// Tensor TS sibling of ApbsEfieldTimeSeriesTrajectoryResult (Vec3 TS
// for the E-field) — same source ConformationResult (ApbsFieldResult),
// different read field.
//
// FO dense-buffer pattern (DenseBuffer<SphericalTensor> via
// AdoptDenseBuffer at Finalize). ApbsFieldResult is unconditionally
// attached in PerFrameExtractionSet, so in production the
// HasResult<ApbsFieldResult>() gate is defensive and
// source_attached_per_frame is all-1. The gate + per-frame mask +
// NaN-fill on absence is still emitted for SDK uniformity per
// OBJECT_MODEL.md "Conditional-attach TR discipline" subsection
// (2026-05-19; all-1 + source_attached_policy="always_attached").
//
// Emission (T2-only — APBS reaction EFG is the symmetrized, trace-projected
// gradient of E; T0 and T1 are structurally zero after the source
// symmetrization and trace projection):
//
//   /trajectory/apbs_efg_time_series/
//     t2             (N, T, 5)  float64 — T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64 — ps
//     source_attached_per_frame (T,) uint8 — canonical SDK contract;
//                                            all-1 in production
//     attrs:
//       result_name             = "ApbsEfgTimeSeriesTrajectoryResult"
//       n_atoms, n_frames, finalized
//       t2_basis                = "project_native_t2_isometric_real_tesseral_v1"
//       t2_component_order      = "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
//       t2_frame                = "cartesian_xyz_emitted_frame"
//       t2_parity               = "even"
//       efg_t0_structural_zero  = true
//       efg_t1_structural_zero  = true
//       e3nn_export             = explicit conversion-note string
//       irrep_layout/normalization/parity = deprecated legacy attrs
//       legacy_irrep_attrs_deprecated = true
//       units                   = "V/Å^2"
//       source                  = "ApbsFieldResult.apbs_efg_spherical
//                                   (SphericalTensor, T2 components 0..4)"
//       source_attached_policy  = "always_attached" — HasResult gate
//                                  emits NaN-fill + mask=0 on absent
//                                  frames per canonical 'absent, not
//                                  faked'.
//
// t2_parity "even": EFG is a rank-2 even-parity tensor derived from the
// gradient of the polar E-field. T0 and T1 are structurally zero from
// source symmetrization and trace projection in ApbsFieldResult.cpp;
// only 5 T2 components are emitted.
// Group tensor/irrep attributes apply only to t2; lab-axis APBS grid
// diagnostics carry dataset-local mixed/non-covariant contracts.
//

#include "DenseBuffer.h"
#include "TrajectoryResult.h"
#include "Types.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class ApbsEfgTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "ApbsEfgTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<ApbsEfgTimeSeriesTrajectoryResult> Create(
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
    // Per-atom growing buffer of SphericalTensor. Transferred into
    // DenseBuffer<SphericalTensor> on Finalize.
    std::vector<std::vector<SphericalTensor>> per_atom_efg_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;
    // Per-frame source-attached mask. All-1 in production because
    // ApbsFieldResult is RequireConformationResult'd; the gate is
    // defensive against non-PerFrameExtractionSet configurations.
    std::vector<std::uint8_t> source_attached_per_frame_;
    std::vector<std::array<std::uint64_t, 3>> grid_dims_per_frame_;
    std::vector<std::array<double, 3>> grid_lengths_A_per_frame_;
    std::vector<std::array<double, 3>> grid_origin_A_per_frame_;
    std::vector<std::array<double, 3>> grid_spacing_A_per_frame_;
    double apbs_manual_grid_padding_A_ =
        std::numeric_limits<double>::quiet_NaN();
    double apbs_manual_grid_min_dim_A_ =
        std::numeric_limits<double>::quiet_NaN();
    double apbs_temperature_K_ =
        std::numeric_limits<double>::quiet_NaN();
    double apbs_thermal_voltage_V_ =
        std::numeric_limits<double>::quiet_NaN();
    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
