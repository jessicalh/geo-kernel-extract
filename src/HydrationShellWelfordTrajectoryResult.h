#pragma once
//
// HydrationShellWelfordTrajectoryResult: per-atom Welford rollup of the
// COM-based water shell features. Four scalar channels, each with full
// signed/abs/squared/dxdt delta variants. Sibling of
// HydrationGeometryWelford (SASA-normal frame) — both methods accumulate
// per feedback_methods_accumulate; calibration weights them separately.
//
// Channels (mirroring HydrationShellTimeSeries source fields):
//   half_shell_asymmetry   scalar (fraction)
//   mean_water_dipole_cos  scalar (cos angle)
//   nearest_ion_distance   scalar (Å)         +inf sentinel on no-ion-in-cutoff
//   nearest_ion_charge     scalar (e)
//
// All four scalars carry per-channel-distinct dynamics questions; full
// delta variants on each.
//
// Inf-sentinel behavior: when an atom never has an ion within cutoff,
// every frame's nearest_ion_distance is +inf. The Welford running mean
// over inf-only data degrades to NaN at the second sample (IEEE 754:
// inf - inf = NaN); downstream `np.isfinite(welford.mean)` filters these
// atoms cleanly. Mixed atoms (some frames in cutoff, some out) also
// produce NaN, signaling "the distribution isn't well-defined."
//
// Dependencies: HydrationShellResult — conditionally attached gated on
// `opts.solvent`. Source-absent frames skip Welford updates entirely.
//

#include "TrajectoryResult.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class HydrationShellWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "HydrationShellWelfordTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<HydrationShellWelfordTrajectoryResult> Create(
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
    std::vector<double> prev_half_shell_asymmetry_;
    std::vector<double> prev_mean_water_dipole_cos_;
    std::vector<double> prev_nearest_ion_distance_;
    std::vector<double> prev_nearest_ion_charge_;
    std::vector<bool>   prev_valid_;
    std::vector<double> prev_time_;

    std::vector<std::uint8_t> source_attached_per_frame_;
    bool                      force_source_present_for_testing_ = false;

    std::size_t              n_frames_  = 0;
    bool                     finalized_ = false;
    double                   mean_dt_ps_ = 0.0;
    std::array<std::size_t, 2> frame_index_range_ = {0, 0};
};

}  // namespace nmr
