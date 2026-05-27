#pragma once
//
// MopacBondOrderWelfordTrajectoryResult: AV (always-valid) per-bond
// Welford rollup of MOPAC Wiberg bond orders. TR6 of the 13-TR plan;
// clone of MopacChargeWelfordTrajectoryResult with the bond axis
// substituted for the atom axis.
//
// Bond-scope TR. Per-bond Welford state lives INSIDE this Result as
// `std::vector<WelfordMoments> per_bond_` parallel to
// `tp.ProteinRef().Bonds()` (same axis as `bonds.npy` emitted by
// TopologySidecar). Sized at Create from `protein.BondCount()` —
// valid because Trajectory::Run Phase 2 (Seed) precedes Phase 3
// (factory invocation), so the Protein is finalised.
//
// SPARSE CADENCE — same "absent, not faked" gate as TR5: MopacResult
// is attached when OperationRunner allows MOPAC and Compute returns
// non-null; it is NOT a RequireConformationResult. Frames without MOPAC skip the Welford
// update and record mask=0. When source_attached_count == 0 across
// the whole trajectory, WriteH5Group skips
// /trajectory/mopac_bond_order_welford/ entirely.
//
// Source: MopacResult.TopologyBondOrders() — std::vector<double>
// parallel to protein.Bonds(). MopacResult sets bond order to 0.0
// (exact) for bonds MOPAC didn't report (NOT NaN). NOTE: the
// MopacResult parser itself filters at `bo > 0.01`, so any Wiberg
// order in (0.0, 0.01] is
// dropped at the parser and arrives here as the 0.0 sentinel —
// indistinguishable from "MOPAC didn't print this bond at all".
// For typical MD this fuses two cases:
//   (1) MOZYME-merged interior bond (MOPAC genuinely didn't report),
//   (2) Transient bond at extension/breaking with Wiberg order
//       below the 0.01 parser threshold (MOPAC printed a small
//       value, the parser dropped it).
// Both contribute to `bo == 0.0` here. The Welford treats the
// fused event as "no observation"; for production-stable bonds
// (case 1 only) this is what the sentinel-aware design wants,
// and for case 2 the bond is by definition not a meaningful
// covalent observation for the calibration target.
//
// SENTINEL-AWARE WELFORD (per `feedback_conditional_welford_for_sentinels`,
// codex R6 2026-05-18): naive accumulation of "no observation"
// sentinels biases the running mean toward 0. Instead, we accumulate
// the order Welford ONLY on frames where the bond was reported
// (`bo != 0.0`) AND emit a companion `order_present_fraction`
// indicator-Welford on the binary "MOPAC reported this bond" event
// (1.0 if `bo != 0.0`, else 0.0). Mirrors the HydrationShellWelford
// `ion_present_fraction` pattern landed for nearest_ion_distance.
//
// MINIMUM-VIABLE v0 (no delta variants on order or present_fraction):
// single channel per Welford, full canonical row.
//
// Emission at /trajectory/mopac_bond_order_welford/:
//   order_mean                     (B,) float64 — dimensionless (Wiberg bond order)
//   order_std                      (B,) float64 — dimensionless
//   order_m2                       (B,) float64 — dimensionless (squared)
//   order_min                      (B,) float64
//   order_max                      (B,) float64
//   order_min_frame                (B,) uint64
//   order_max_frame                (B,) uint64
//   n_per_bond                     (B,) uint64  — order Welford divisor
//                                                 (frames where MOPAC
//                                                 reported the bond)
//   order_present_fraction_mean    (B,) float64 — dimensionless in [0, 1]
//                                                 (Pr(MOPAC reports bond))
//   order_present_fraction_std     (B,) float64
//   order_present_fraction_m2      (B,) float64
//   order_present_fraction_min     (B,) float64
//   order_present_fraction_max     (B,) float64
//   order_present_fraction_min_frame (B,) uint64
//   order_present_fraction_max_frame (B,) uint64
//   n_total_per_bond               (B,) uint64  — present_fraction
//                                                 Welford divisor =
//                                                 source_attached_count
//                                                 (MOPAC-attached frames)
//   frame_indices                  (T,) uint64
//   frame_times                    (T,) float64 — ps
//   source_attached_per_frame      (T,) uint8
//
// Attrs:
//   result_name             = "MopacBondOrderWelfordTrajectoryResult"
//   n_bonds, n_frames, source_attached_count, finalized
//   units                   = "dimensionless"
//   bond_axis               = "bonds.npy" (canonical sidecar axis)
//   source                  = "MopacResult.TopologyBondOrders() ..."
//   source_attached_policy  = "conditional -- MopacResult ..."
//

#include "TrajectoryMoments.h"  // WelfordMoments
#include "TrajectoryResult.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class MopacBondOrderWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "MopacBondOrderWelfordTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<MopacBondOrderWelfordTrajectoryResult> Create(
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
    std::size_t BondCount() const { return per_bond_.size(); }

private:
    // Per-bond Welford state (axis parallel to protein.Bonds() /
    // bonds.npy). Sized at Create from protein.BondCount().
    //
    // per_bond_              : Welford on bo (only updated when bo != 0)
    // per_bond_n_present_    : count of frames where bo != 0
    //                          (= divisor for per_bond_ Welford)
    // per_bond_present_      : indicator Welford on (bo != 0 ? 1.0 : 0.0)
    //                          (mean ∈ [0, 1] = Pr(MOPAC reports bond))
    //                          divisor = source_attached_count_ (every
    //                          MOPAC-attached frame contributes an
    //                          indicator sample)
    std::vector<WelfordMoments> per_bond_;
    std::vector<std::size_t>    per_bond_n_present_;
    std::vector<WelfordMoments> per_bond_present_;
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    std::size_t n_frames_              = 0;
    std::size_t source_attached_count_ = 0;
    bool        finalized_             = false;
};

}  // namespace nmr
