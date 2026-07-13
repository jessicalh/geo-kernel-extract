#pragma once
//
// LarsenHBond1pHBShieldingTimeSeriesTrajectoryResult: per-atom per-frame
// time series of the Larsen 1°HB shielding contribution
// (ConformationAtom::larsen_hbond_1pHB_spherical, SphericalTensor in
// ppm, protein lab frame). Finalize-only dense-buffer pattern;
// mirrors TripeptideBackboneShieldingTimeSeriesTrajectoryResult.
//
// The 1°HB contribution is one of four per-class fields stored
// separately for ML feature stratification — see LarsenContribDispatch
// in LarsenHBondShieldingResult.h and Larsen Table 2.
//
// H5 attrs on /trajectory/larsen_hbond_1pHB_shielding_time_series:
//   source_attached_count  = uint64 count of present-source frames
//   source_attached_policy = "conditional_larsen_grid_source"
//   atom_axis              = "protein_atom_index"
//   frame_axis             = "trajectory_frame_row"
//   irrep_layout  = "PackFull9: [T0, T1_cartesian_xyz, T2_real_tesseral_m-2..m+2]"
//   normalization = "isometric_real_sph"
//   parity        = "mixed"
//   coordinate_frame = "conformation_cartesian_xyz"
//   transformation = "even_rank2 under proper rotations: ..."
//   units         = "ppm"
//
// The signed-rho DFT-grid lookup is chirality-conditioned.  Its packed
// shielding tensor obeys T'=R T R^T under proper rotations, but no
// improper-transform parity is claimed.
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

class LarsenHBond1pHBShieldingTimeSeriesTrajectoryResult
    : public TrajectoryResult {
public:
    std::string Name() const override {
        return "LarsenHBond1pHBShieldingTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<LarsenHBond1pHBShieldingTimeSeriesTrajectoryResult>
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
