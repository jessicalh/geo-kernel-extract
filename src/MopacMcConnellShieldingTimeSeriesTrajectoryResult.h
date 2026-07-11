#pragma once
//
// MopacMcConnellShieldingTimeSeriesTrajectoryResult: per-atom per-frame
// time series of the MOPAC-weighted McConnell bond-anisotropy kernel
// (ConformationAtom.mopac_mc_shielding_contribution, SphericalTensor,
// units Å⁻³ -- bond-order-weighted D(r)Qhat response).
// NOTE: despite the historical "shielding_contribution" field name,
// the stored value is the bare kernel — no Δχ × γ multiplication
// at extraction. Sibling FF14SB `mc_shielding_contribution` carries
// the same convention; `McConnellShieldingTimeSeriesTrajectoryResult`
// emits `units = "Angstrom^-3"`.
//
// TR8 of the 13-TR plan. Similar shape to TR7 (conditional-source T2 TS)
// but emits ALL 9 components (T0+T1+T2) because the source field is
// generally not symmetric-traceless; T0/T1/T2 are all meaningful.
//
// SOURCE STRUCTURE: McConnellResult sets
//   ca.mopac_mc_shielding_contribution =
//       SphericalTensor::Decompose(sum bo_s * D(r_is) * Qhat_s);
// where missing/below-floor bond order zeros only the BO channel.
//
// CONDITIONAL SOURCE: MopacMcConnellResult attaches via TimedAttach, not
// RequireConformationResult. Same gate as TR5/TR6/TR7; attached samples
// ride the trajectory's dispatched frames rather than a MOPAC-specific
// cadence.
//
// Emission:
//   /trajectory/mopac_mc_shielding_time_series/
//     xyz            (N, T, 9)  float64 — T0, T1_x, T1_y, T1_z,
//                                          T2_m-2, T2_m-1, T2_m0, T2_m+1,
//                                          T2_m+2  (Å⁻³)
//     T0 channel: trace(DQhat)/3, the PCS scalar branch.
//     T1 channel: even antisymmetric pseudovector from non-commuting D
//       and Qhat.
//     T2 channel: symmetric traceless branch of the source response.
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64 — ps
//     source_attached_per_frame (T,) uint8
//     attrs:
//       result_name             = "MopacMcConnellShieldingTimeSeriesTrajectoryResult"
//       n_atoms, n_frames, source_attached_count, finalized
//       tensor_basis            = "project_native_full9_spherical_tensor_v1"
//       tensor_component_order  = "T0,T1_x,T1_y,T1_z,T2_m-2,...,T2_m+2"
//       tensor_frame            = "conformation_cartesian_xyz"
//       tensor_parity           = "even"
//       e3nn_export             = "explicit project-basis ... conversion required"
//       normalization           = "isometric_real_sph"
//       units                   = "Angstrom^-3"  (bare kernel, pre-Δχ)
//       source                  = "MopacMcConnellResult.mopac_mc_shielding_contribution
//                                  (unscaled BO channel from D(r)Qhat)"
//       source_attached_policy  = "conditional -- TimedAttach ..."
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

class MopacMcConnellShieldingTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "MopacMcConnellShieldingTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<MopacMcConnellShieldingTimeSeriesTrajectoryResult>
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
    std::vector<std::vector<SphericalTensor>> per_atom_;
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    std::size_t n_frames_              = 0;
    std::size_t source_attached_count_ = 0;
    bool        finalized_             = false;
};

}  // namespace nmr
