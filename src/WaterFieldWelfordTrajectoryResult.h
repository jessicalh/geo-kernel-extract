#pragma once
//
// WaterFieldWelfordTrajectoryResult: per-atom Welford rollup of the
// explicit-water E-field + EFG kernel from WaterFieldResult. Clones the
// BS Welford 5-part shape: per-component Welford on vector channels,
// T2-only on the EFG (T0 and T1 are structural zeros — see schema note
// below), plus signed/abs/squared/dxdt delta variants on the primary
// rotationally-invariant scalar channels (E-field magnitude, shell-
// occupancy counts).
//
// Channels (mirroring WaterFieldTimeSeries source fields):
//
//   efield                Vec3 V/Å              per-component[3] + magnitude
//   efield_first          Vec3 V/Å              per-component[3] + magnitude
//   efg                   T2-only V/Å²          T2[5] + |T2|  (5+1 channels)
//   efg_first             T2-only V/Å²          T2[5] + |T2|  (5+1 channels)
//   n_first               int (dimensionless)   full scalar Welford
//   n_second              int (dimensionless)   full scalar Welford
//
// Delta variants on: efield_magnitude, n_first, n_second. Three primary
// rotationally-invariant scalars carry the per-channel-distinct dynamics
// question. efg_t0 deltas are NOT emitted — T0 is structurally zero.
//
// EFG T0+T1 schema rev (2026-05-18 codex F4): water EFG is built from
// symmetric r⊗r outer products (WaterFieldResult.cpp:130) and explicitly
// traceless-projected (WaterFieldResult.cpp:147-150). After projection
// T0 = trace = 0; T1 = antisymmetric pseudovector = 0 because the
// decomposition reads `0.5*(s_ij - s_ji)` (Types.cpp:31) which is bit-
// exact zero on symmetric input. Only T2 (symmetric-traceless,
// 5 real-spherical-tesseral components m=-2..+2) carries signal.
//
// Emission:
//
//   /trajectory/water_field_welford/
//     efield_{x,y,z}_{mean,m2,std,min,max,min_frame,max_frame}     (N,)
//     efield_magnitude_{mean,...}                                  (N,)
//     efield_first_{x,y,z}_{...}                                   (N,)
//     efield_first_magnitude_{...}                                 (N,)
//     efg_t2_{...}                                                 (N, 5)  m=-2..+2
//     efg_t2magnitude_{...}                                        (N,)
//     efg_first_t2_{...}                                           (N, 5)
//     efg_first_t2magnitude_{...}                                  (N,)
//     n_first_{...}, n_second_{...}                                (N,)
//     <scalar>_delta_{...}, <scalar>_abs_delta_{...},
//     <scalar>_delta_squared_{...}, <scalar>_dxdt_{...}           (N,) for the 3 primary scalars
//     n_frames_per_atom, delta_n_per_atom, dxdt_n_per_atom        (N,)
//     source_attached_per_frame                                    (T,) uint8 provenance
//     attrs: result_name, n_frames, source_attached_count, finalized,
//            ddof, mean_dt_ps (attached-subset), frame_index_range
//            (attached-subset), irrep_layout_efield, irrep_layout_efg_t2,
//            efg_t0_structural_zero=true, efg_t1_structural_zero=true
//
// Dependencies: WaterFieldResult — REQUIRED by PerFrameExtractionSet but
// conditionally attached by OperationRunner: if `opts.solvent` is null
// or `opts.solvent->Empty()` (no solvent environment loaded), the source
// ConformationResult is silently skipped. Follows "absent, not faked":
//   - Per-frame `conf.HasResult<WaterFieldResult>()` check in Compute
//   - Source-absent frames SKIP Welford updates entirely (no biased zero
//     accumulation) and SKIP prev_* cache updates (next-frame delta is
//     not computed across a gap)
//   - `source_attached_per_frame` mask emitted as H5 provenance
//   - WriteH5Group skips group emission when source attached zero times
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

class WaterFieldWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "WaterFieldWelfordTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<WaterFieldWelfordTrajectoryResult> Create(
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

    // Test-only bypass for the source-attached gate.
    void ForceSourcePresentForTesting() {
        force_source_present_for_testing_ = true;
    }

private:
    // Per-atom previous-frame caches for the delta trackers. No
    // prev_efg_t0_ cache — `efg_t0` channel removed per 2026-05-18
    // adversarial review (structurally zero after the traceless
    // projection in `WaterFieldResult.cpp:147-150`).
    std::vector<double>      prev_efield_mag_;
    std::vector<double>      prev_n_first_;
    std::vector<double>      prev_n_second_;
    std::vector<bool>        prev_valid_;
    std::vector<double>      prev_time_;

    std::vector<std::uint8_t>  source_attached_per_frame_;
    bool                       force_source_present_for_testing_ = false;

    std::size_t              n_frames_  = 0;
    bool                     finalized_ = false;
    double                   mean_dt_ps_ = 0.0;
    std::array<std::size_t, 2> frame_index_range_ = {0, 0};
};

}  // namespace nmr
