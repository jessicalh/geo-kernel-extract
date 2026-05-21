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
// is TimedAttach'd (OperationRunner.cpp:142), not Required. Frames
// without MOPAC skip the Welford update and record mask=0. When
// source_attached_count == 0 across the whole trajectory,
// WriteH5Group skips /trajectory/mopac_bond_order_welford/ entirely.
//
// Source: MopacResult.TopologyBondOrders() — std::vector<double>
// parallel to protein.Bonds(). MopacResult sets bond order to 0
// for bonds MOPAC didn't report (NOT NaN), so we don't gate per
// bond — every bond gets the same number of Welford updates as
// the protein-wide source_attached_count.
//
// MINIMUM-VIABLE v0 (no delta variants): single channel, full
// canonical Welford row. Mirrors TR5 + AIMNet2 CRG Welford v0
// precedent.
//
// Emission at /trajectory/mopac_bond_order_welford/:
//   order_mean          (B,) float64 — dimensionless (Wiberg bond order)
//   order_std           (B,) float64 — dimensionless
//   order_m2            (B,) float64 — dimensionless (squared)
//   order_min           (B,) float64
//   order_max           (B,) float64
//   order_min_frame     (B,) uint64
//   order_max_frame     (B,) uint64
//   n_per_bond          (B,) uint64
//   frame_indices       (T,) uint64
//   frame_times         (T,) float64 — ps
//   source_attached_per_frame (T,) uint8
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
    std::vector<WelfordMoments> per_bond_;
    std::vector<std::size_t>    per_bond_n_;
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    std::size_t n_frames_              = 0;
    std::size_t source_attached_count_ = 0;
    bool        finalized_             = false;
};

}  // namespace nmr
