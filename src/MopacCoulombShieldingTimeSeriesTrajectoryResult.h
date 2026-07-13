#pragma once
//
// MopacCoulombShieldingTimeSeriesTrajectoryResult: per-atom per-frame
// time series of the MOPAC-Coulomb T2 EFG kernel
// (ConformationAtom.mopac_coulomb_shielding_contribution, T2
// SphericalTensor, units V/Å²). NOTE: despite the historical
// "shielding_contribution" field name, the stored value is the raw
// EFG kernel — no γ multiplication is applied at the source. The
// MopacCoulombResult source comment "gamma converts this to
// shielding" is forward-looking: shielding is recovered downstream
// by multiplying by per-element γ at calibration time.
//
// TR7 of the 13-TR plan. Combines:
//   - TR4 pattern: T2-only (N, T, 5) emission via DenseBuffer<SphericalTensor>
//   - TR5 gate:    conditional-source HasResult<MopacCoulombResult>() skip,
//                  NaN-fill on absent frames, source_attached_per_frame
//                  mask, group-absent when source never attached.
//
// SOURCE STRUCTURE: MopacCoulombResult stores
// `SphericalTensor::Decompose(EFG_total)`, where EFG_total is the
// trace-projected dipolar EFG from the legacy F15.6 projection of MOPAC
// Coulson charges. T0 and T1
// are structurally zero, so (N, T, 5) emission is information-preserving.
//
// CONDITIONAL SOURCE: MopacCoulombResult attaches via TimedAttach, not
// RequireConformationResult. There is no MOPAC-specific cadence; when MOPAC
// is enabled, attached samples are on the same frames dispatched by the
// trajectory's single --stride and any dispatch window. WriteH5Group skips
// both the canonical /trajectory/mopac_coulomb_efg_time_series/ group and
// the legacy alias group when no frame attached the source.
//
// Emission:
//   /trajectory/mopac_coulomb_efg_time_series/
//     t2             (N, T, 5)  float64 — T2_m-2..T2_m+2 (V/Å²)
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64 — ps
//     source_attached_per_frame (T,) uint8
//     attrs:
//       result_name             = "MopacCoulombEfgTimeSeriesTrajectoryResult"
//       n_atoms, n_frames, source_attached_count, finalized
//       quantity                = "raw_efg_t2"
//       historical_field_name   = "mopac_coulomb_shielding_contribution"
//       gamma_applied           = false
//       t2_basis                = "project_native_t2_isometric_real_tesseral_v1"
//       t2_component_order      = "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
//       t2_frame                = "cartesian_xyz_emitted_frame"
//       t2_parity               = "even"
//       legacy_irrep_attrs_deprecated = true
//       units                   = "V/Å^2"  (EFG kernel, pre-γ)
//       source_result/source_field/source_operation/source_tensor pin the
//                                  symbol-level provenance
//       source_attached_policy  = "conditional -- MopacCoulombResult attaches
//                                  only on dispatched frames where MOPAC
//                                  is enabled and Compute succeeds ..."
//   /trajectory/mopac_coulomb_shielding_time_series/
//     Legacy duplicate of the canonical group, with alias_of and
//     legacy_name_deprecated attrs.
//

#include "DenseBuffer.h"
#include "TrajectoryResult.h"
#include "Types.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class MopacCoulombShieldingTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "MopacCoulombEfgTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<MopacCoulombShieldingTimeSeriesTrajectoryResult>
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
    // Per-atom growing buffer of SphericalTensor. Transferred to
    // DenseBuffer<SphericalTensor> at Finalize.
    std::vector<std::vector<SphericalTensor>> per_atom_;
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    std::size_t n_frames_              = 0;
    std::size_t source_attached_count_ = 0;
    bool        finalized_             = false;
};

}  // namespace nmr
