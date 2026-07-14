#pragma once
//
// LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult: per-atom per-frame
// time series of the Larsen 2°HαB shielding contribution
// (ConformationAtom::larsen_hbond_2pHaB_spherical, SphericalTensor in
// ppm, protein lab frame). One of four per-class Larsen TRs cloned
// from TripeptideBackboneShieldingTimeSeriesTrajectoryResult — see
// LarsenContribDispatch in LarsenHBondShieldingResult.h.
// H5 contract: source_attached_count (uint64),
// source_attached_policy="conditional_larsen_grid_source",
// atom_axis="protein_atom_index", frame_axis="trajectory_frame_row", and
// irrep_layout="PackFull9: [T0, T1_cartesian_xyz,
// T2_real_tesseral_m-2..m+2]", parity="mixed",
// coordinate_frame/tensor_frame="conformation_cartesian_xyz",
// tensor_basis="project_native_full9_spherical_tensor_v1",
// tensor_component_order="T0,T1_x,T1_y,T1_z,T2_m-2,...,T2_m+2",
// tensor_t1_structural_zero=false,
// tensor_structural_zero_components="none", explicit e3nn conversion and
// normalization-scope attrs, and an explicit proper-rotation-only
// transformation law. The T1 convention is the exact axial Cartesian
// Levi-Civita dual a=((T_yz-T_zy)/2,(T_zx-T_xz)/2,(T_xy-T_yx)/2), not
// real-Y1m, and is generically nonzero. The signed-rho DFT-grid lookup is
// chirality-conditioned, so no improper-transform parity is claimed.
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

class LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult
    : public TrajectoryResult {
public:
    std::string Name() const override {
        return "LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult>
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

    // Test-only: override the source-present state for subsequent frames.
    // May be toggled between synthetic Compute calls to exercise mixed
    // absent/present provenance without an external trajectory or grid.
    // Production code never calls this.
    void ForceSourcePresentForTesting(bool present = true) {
        force_source_present_for_testing_ = present;
    }

private:
    std::vector<std::vector<SphericalTensor>> per_atom_shielding_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;

    // Per-frame source-present mask. 1 if the source ConformationResult
    // (LarsenHBondShieldingResult) was attached this frame, or if the
    // test-only override forces it present; 0 otherwise. Emitted as the
    // `source_attached_per_frame` H5 dataset for downstream provenance.
    // When all-zero (calc never ran), WriteH5Group skips emission
    // entirely per the "absent, not faked" discipline.
    std::vector<std::uint8_t> source_present_per_frame_;
    bool force_source_present_for_testing_ = false;

    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
