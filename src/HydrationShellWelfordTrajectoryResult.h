#pragma once
//
// HydrationShellWelfordTrajectoryResult: per-atom Welford rollup of the
// COM-based water shell features. Four scalar channels, each with full
// signed/abs/squared/dxdt delta variants. Sibling of
// HydrationGeometryWelford (SASA-normal frame) — both methods accumulate
// per feedback_methods_accumulate; calibration weights them separately.
//
// Channels (mirroring HydrationShellTimeSeries source fields):
//   half_shell_asymmetry   scalar (fraction)         full Welford + delta variants
//   mean_water_dipole_cos  scalar (cos angle)        full Welford + delta variants
//   nearest_ion_distance   scalar (Å)                conditional Welford (R6)
//   nearest_ion_charge     scalar (e)                full Welford + delta variants
//   ion_present_fraction   scalar [0,1] (R6)         Welford only — no deltas
//
// half_shell_asymmetry / mean_water_dipole_cos / nearest_ion_charge carry
// per-channel-distinct dynamics questions; full delta variants on each.
//
// Nearest-ion conditional Welford (R6 codex 2026-05-18): naively
// accumulating `nearest_ion_distance = +inf` into a Welford NaN-poisoned
// any atom that was MIXED contact/no-contact across frames — and at
// the 20 Å cutoff most protein atoms ARE mixed. Schema rev:
//   - `nearest_ion_distance` Welford accumulates ONLY when the sample
//     is finite. Per-atom counter `n_ion_present` becomes the
//     Finalize divisor. Atoms that never see an ion in cutoff finalize
//     to NaN via WelfordFinalize(_, 0) — consistent with the
//     "uncomputable" sentinel discipline.
//   - Delta variants on `nearest_ion_distance` accumulate only when
//     BOTH the current and previous frames had a finite ion distance
//     (per-atom counter `n_ion_delta`). dxdt likewise (`n_ion_dxdt`).
//   - New channel `ion_present_fraction`: Welford on the 1.0/0.0
//     indicator that nearest_ion_distance is finite. Its mean ∈ [0, 1]
//     is interpretable as Pr(ion in cutoff this frame).
// Conditional Welford + presence indicator give the actually-meaningful
// summary: "P(ion present)" + "E(distance | ion present)" + standard
// deviation of the conditional distribution.
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
    // R6: separately track whether the previous attached frame had a
    // FINITE nearest_ion_distance. Distinct from prev_valid_ which
    // tracks source-attached gate; this tracks the inf-vs-finite branch
    // ONLY for the nearest_ion_distance channel's delta computation.
    std::vector<bool>   prev_ion_finite_;
    std::vector<double> prev_time_;

    std::vector<std::uint8_t> source_attached_per_frame_;
    bool                      force_source_present_for_testing_ = false;

    std::size_t              n_frames_  = 0;
    bool                     finalized_ = false;
    double                   mean_dt_ps_ = 0.0;
    std::array<std::size_t, 2> frame_index_range_ = {0, 0};
};

}  // namespace nmr
