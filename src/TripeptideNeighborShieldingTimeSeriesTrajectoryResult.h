#pragma once
//
// TripeptideNeighborShieldingTimeSeriesTrajectoryResult: per-atom
// per-frame time series of the TripeptideNeighborShieldingResult
// Δσ_BB^{i±1} contribution as a SphericalTensor. Finalize-only
// dense-buffer pattern, mirrors BsShieldingTimeSeriesTrajectoryResult
// against the tripeptide_neighbor_shielding_spherical source field
// on ConformationAtom.
//
// Per Larsen 2015 Eq 3 the per-residue contribution
// Δσ_BB^{i-1}(i) + Δσ_BB^{i+1}(i) is read at the flanking ALA cap
// atoms of the (i±1)-centered tripeptides, with the AAA reference
// at standard angles (φ_std=-120°, ψ_std=140°) subtracted; the per-
// frame ConformationResult writes the SUM onto each central atom's
// tripeptide_neighbor_shielding_spherical field. This TR captures
// that field across the trajectory.
//
// Emission uses the SphericalTensor::PackFull9 payload order:
//
//   /trajectory/tripeptide_neighbor_shielding_time_series/
//     xyz            (N, T, 9)  float64
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name    = "TripeptideNeighborShieldingTimeSeriesTrajectoryResult"
//       irrep_layout   = "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
//       normalization  = "isometric_real_sph"
//       parity         = "mixed"
//       coordinate_frame = "conformation_cartesian_xyz"
//       tensor_basis   = "project_native_full9_spherical_tensor_v1"
//       tensor_component_order = "T0,T1_x,T1_y,T1_z,T2_m-2,...,T2_m+2"
//       tensor_frame   = "conformation_cartesian_xyz"
//       tensor_t1_semantics = "Cartesian Levi-Civita dual ..."
//       tensor_t1_structural_zero = false
//       tensor_structural_zero_components = "none"
//       e3nn_export    = "explicit project-basis ... conversion required"
//       normalization_scope = "xyz tensor payload: T2 uses ..."
//       transformation = "even_rank2 under proper rotations: ..."
//       units          = "ppm"
//       n_atoms, n_frames, finalized
//
// The packed tensor obeys T'=R T R^T under proper rotations.  Its typed
// lookup/proper-Kabsch alignment is conditioned on an unchanged chiral
// L-amino-acid DFT source, so no improper-transform parity is claimed.
// The T1 entries are the exact Cartesian Levi-Civita dual
// ((T_yz-T_zy)/2, (T_zx-T_xz)/2, (T_xy-T_yx)/2), not real-Y1m, and are
// generically nonzero. T2 alone uses the isometric real-tesseral basis.
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

class TripeptideNeighborShieldingTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "TripeptideNeighborShieldingTimeSeriesTrajectoryResult";
    }

    // No trajectory-scope dependencies; the underlying ConformationResult
    // (TripeptideNeighborShieldingResult) is conditionally attached when
    // the [databases].tensorcs15 table is configured. This TR captures
    // the per-atom tensor and `tripeptide_neighbor_has_match`
    // applicability state each frame. H5 rows are NaN when the source did
    // not attach or the atom had no neighboring tripeptide contribution.
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<TripeptideNeighborShieldingTimeSeriesTrajectoryResult> Create(
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

    // Test-only: bypass the per-frame `conf.HasResult<SourceCalc>()`
    // check inside Compute. Call BEFORE the first Compute. Use ONLY
    // in synthetic unit tests that don't go through Trajectory::Run /
    // OperationRunner (the production path that attaches the source
    // ConformationResult). Production code never calls this.
    void ForceSourcePresentForTesting() {
        force_source_present_for_testing_ = true;
    }

private:
    // Per-atom growing buffers of SphericalTensor. Flattened into an
    // atom-major DenseBuffer<SphericalTensor> at Finalize.
    std::vector<std::vector<SphericalTensor>> per_atom_shielding_;
    // Internal atom-by-frame applicability, captured from
    // ConformationAtom::tripeptide_neighbor_has_match. It is intentionally
    // not emitted as a new H5 dataset; it only controls NaN filling of xyz.
    std::vector<std::vector<std::uint8_t>> per_atom_has_match_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;

    // Per-frame source-present mask. 1 if the source ConformationResult
    // (TripeptideNeighborShieldingResult) was attached this frame, or if
    // the test-only override forces it present; 0 otherwise. Emitted as
    // the `source_attached_per_frame` H5 dataset for downstream provenance.
    // When all-zero (calc never ran), WriteH5Group skips emission
    // entirely per the "absent, not faked" discipline.
    std::vector<std::uint8_t> source_present_per_frame_;
    bool force_source_present_for_testing_ = false;

    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
