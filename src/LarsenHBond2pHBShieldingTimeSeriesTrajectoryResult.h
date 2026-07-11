#pragma once
//
// LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult: per-atom per-frame
// time series of the Larsen 2°HB shielding contribution
// (ConformationAtom::larsen_hbond_2pHB_spherical, SphericalTensor in
// ppm, protein lab frame). One of four per-class Larsen TRs cloned
// from TripeptideBackboneShieldingTimeSeriesTrajectoryResult — see
// LarsenContribDispatch in LarsenHBondShieldingResult.h.
// H5 contract: source_attached_count (uint64),
// source_attached_policy="conditional_larsen_grid_source",
// atom_axis="protein_atom_index", frame_axis="trajectory_frame_row", and
// irrep_layout="PackFull9: [T0, T1_cartesian_xyz,
// T2_real_tesseral_m-2..m+2]".
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

class LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult
    : public TrajectoryResult {
public:
    std::string Name() const override {
        return "LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult>
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
