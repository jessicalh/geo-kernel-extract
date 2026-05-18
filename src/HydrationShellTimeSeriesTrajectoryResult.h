#pragma once
//
// HydrationShellTimeSeriesTrajectoryResult: per-atom per-frame timeline
// of HydrationShellResult's 4 scalar channels. The COM-based older sibling
// of HydrationGeometryTimeSeries — separate physics call, kept side-by-side
// per the methods-accumulate discipline (feedback_methods_accumulate).
//
// Channels:
//   half_shell_asymmetry   scalar (fraction)  COM-reference frame
//   mean_water_dipole_cos  scalar             water orientation order parameter
//   nearest_ion_distance   scalar (Å)         distance to nearest ion in cutoff;
//                                              +infinity sentinel when no ion
//                                              within ion_cutoff (default 20 Å)
//   nearest_ion_charge     scalar (e)         charge of nearest ion; 0 when no ion
//
// "+infinity sentinel" on nearest_ion_distance is meaningful — distinct from
// "0 distance" and from "missing measurement." Buried atoms with no ion
// within cutoff carry +inf forever; consumers should `np.isfinite()` filter
// the channel before naive aggregation. The Welford rollup propagates inf
// → NaN once a second inf sample arrives (IEEE 754 inf-inf = NaN), which
// is the right signal for "this atom's distribution isn't well-defined."
//
// Dependencies: HydrationShellResult — REQUIRED by PerFrameExtractionSet
// but conditionally attached by OperationRunner: if `opts.solvent` is
// null/empty, the source ConformationResult is silently skipped. Follows
// "absent, not faked":
//   - Per-frame `conf.HasResult<HydrationShellResult>()` check in Compute
//   - `source_attached_per_frame` uint8 mask emitted as H5 provenance
//   - NaN-fill float channels on source-absent frames
//   - WriteH5Group skips entire group emission when source attached zero times
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class HydrationShellTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "HydrationShellTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<HydrationShellTimeSeriesTrajectoryResult> Create(
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

    void ForceSourcePresentForTesting() {
        force_source_present_for_testing_ = true;
    }

private:
    std::vector<std::vector<double>> half_shell_asymmetry_;
    std::vector<std::vector<double>> mean_water_dipole_cos_;
    std::vector<std::vector<double>> nearest_ion_distance_;
    std::vector<std::vector<double>> nearest_ion_charge_;
    std::vector<std::size_t>         frame_indices_;
    std::vector<double>              frame_times_;
    std::vector<std::uint8_t>        source_attached_per_frame_;
    bool                             force_source_present_for_testing_ = false;
    std::size_t                      n_frames_ = 0;
    bool                             finalized_ = false;
};

}  // namespace nmr
