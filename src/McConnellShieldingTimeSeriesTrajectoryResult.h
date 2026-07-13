#pragma once
//
// McConnellShieldingTimeSeriesTrajectoryResult: per-atom per-frame
// time series of the McConnell bond-anisotropy shielding contribution
// as a SphericalTensor (geometric kernel, Å⁻³). FO dense-buffer
// pattern, clones BsShieldingTimeSeriesTrajectoryResult against
// ConformationAtom::mc_shielding_contribution.
//
// Source is McConnellResult, unconditionally attached in
// PerFrameExtractionSet. No source-attached gate.
//
// Units: Angstrom^-3 (NOT ppm). McConnell's field on ConformationAtom
// stores the unparameterised D(r)Qhat source-shape response
// (T0+T1+T2 all physical). Susceptibility scale, sign, and unit
// conversion are learned downstream in the calibration pipeline.
//
// Emission:
//
//   /trajectory/mc_shielding_time_series/
//     xyz            (N, T, 9)  float64
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name             = "McConnellShieldingTimeSeriesTrajectoryResult"
//       tensor_basis            = "project_native_full9_spherical_tensor_v1"
//       tensor_component_order  = "T0,T1_x,T1_y,T1_z,T2_m-2,...,T2_m+2"
//       tensor_frame            = "conformation_cartesian_xyz"
//       coordinate_frame        = "conformation_cartesian_xyz"
//       parity                 = "0e+1e+2e"
//       tensor_parity           = "even"
//       tensor_transformation   = "even_rank2: T'=R T R^T"
//       tensor_t1_semantics     = "Cartesian Levi-Civita dual ..."
//       tensor_t1_structural_zero = false
//       tensor_structural_zero_components = "none"
//       e3nn_export             = "explicit project-basis ... conversion required"
//       normalization           = "isometric_real_sph"
//       normalization_scope     = "xyz tensor payload: T2 uses ..."
//       units                   = "Angstrom^-3"
//
// McConnellResult forms D(r) Qhat.  Under every orthogonal transform R,
// D' = R D R^T and Qhat' = R Qhat R^T, hence the (generally nonsymmetric)
// product transforms as T' = R T R^T.  Its antisymmetric T1 channel is the
// Cartesian Levi-Civita dual from SphericalTensor::Decompose (an axial 1e
// vector), not a real-Y1m vector, and is not structurally zero.
//

#include "DenseBuffer.h"
#include "TrajectoryResult.h"
#include "Types.h"

#include <cstddef>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class McConnellShieldingTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "McConnellShieldingTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<McConnellShieldingTimeSeriesTrajectoryResult> Create(
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
    std::vector<std::vector<SphericalTensor>> per_atom_shielding_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;
    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
