#pragma once
//
// MopacChargeWelfordTrajectoryResult: AV (always-valid) per-atom
// Welford rollup of MOPAC Mulliken charges (ConformationAtom.mopac_charge,
// units e — elementary charges). TR5 of the 13-TR plan; canonical
// sparse-Welford-scalar template (the bond-order companion TR6
// clones this shape against the bond axis).
//
// SPARSE CADENCE — "absent, not faked": MopacResult is attached
// conditionally by OperationRunner (gated by `!opts.skip_mopac` and
// a non-null Compute return; NOT RequireConformationResult). MOPAC runs
// every ~20 ps in
// production (CLI-driven; the TR is cadence-agnostic). On
// frames where MopacResult is absent, the per-frame `HasResult` gate
// skips the Welford update and does not increment source_attached_count_;
// the per-frame mask records 0. When source_attached_count == 0 across
// the whole trajectory (e.g. MOPAC disabled), WriteH5Group skips the
// /trajectory/mopac_charge_welford/ group entirely (canonical
// "absent, not faked" — see OBJECT_MODEL.md "Conditional-attach TR
// discipline").
//
// MINIMUM-VIABLE v0 (no delta variants): single channel, full
// canonical Welford row (mean / std / m2 / min / max / min_frame /
// max_frame) + n_per_atom. Mirrors the AIMNet2ChargeResponseGradient
// Welford v0 precedent — MOPAC charges are continuous (not discrete
// count transitions like HBondCount), so per-atom std already
// captures fluctuation. Delta variants (charge_delta_*, charge_dxdt_*)
// can be added later if a calibration finding requests them.
//
// Emission at /trajectory/mopac_charge_welford/:
//
//   charge_mean          (N,) float64 — e
//   charge_std           (N,) float64 — e  (= sqrt(m2/(n-1)) for n>=2;
//                                            0 for n=1; NaN for n=0)
//   charge_m2            (N,) float64 — e²
//   charge_min           (N,) float64 — e
//   charge_max           (N,) float64 — e
//   charge_min_frame     (N,) uint64  — frame_index at which min was hit
//   charge_max_frame     (N,) uint64  — frame_index at which max was hit
//   n_per_atom           (N,) uint64  — per-atom samples accumulated
//   frame_indices        (T,) uint64
//   frame_times          (T,) float64 — ps
//   source_attached_per_frame (T,) uint8 — canonical SDK contract
//
// Attrs:
//   result_name             = "MopacChargeWelfordTrajectoryResult"
//   n_atoms, n_frames, source_attached_count, finalized
//   units                   = "e" (elementary charge)
//   source                  describes MopacResult.mopac_charge (Mulliken)
//   source_attached_policy  = "conditional -- MopacResult attaches
//                              sparsely per Mopac cadence
//                              (OperationRunner TimedAttach not
//                              Require). Compute's HasResult gate
//                              ..."
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class MopacChargeWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "MopacChargeWelfordTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<MopacChargeWelfordTrajectoryResult> Create(
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
    std::size_t SourceAttachedCount() const { return source_attached_count_; }

private:
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    std::size_t               n_frames_              = 0;
    std::size_t               source_attached_count_ = 0;
    bool                      finalized_             = false;
};

}  // namespace nmr
