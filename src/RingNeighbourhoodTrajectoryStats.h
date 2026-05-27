#pragma once
//
// RingNeighbourhoodTrajectoryStats -- FO per-(atom, aromatic-ring)
// geometric residual across the trajectory. TR10 of the 13-TR plan;
// S3 locked-scope A3.2 (the geometric-residual-only TR; see also
// A3.0 gaussian_density delete + A3.1 ring_membership_per_atom
// emitted in this TR's H5 group).
//
// Per-frame, per (atom, ring-in-static-snapshot): 4-channel geometric
// residual:
//   0 distance       (A)      atom-to-ring-center
//   1 rho            (A)      in-plane radial distance
//   2 z              (A)      out-of-plane projection (signed; z > 0
//                             means same side as ring normal vector)
//   3 in_plane_angle (rad)    azimuth in ring plane from vertex 0
//                             toward atom, in [0, 2*pi); NaN when the
//                             atom is on the ring axis (rho <
//                             MIN_DISTANCE = 0.1 Å -- the project-wide
//                             singularity-guard threshold).
//
// Static (atom, ring) cutoff snapshot frozen at first Compute call
// (frame 0): for each atom, the list of aromatic-ring indices within
// `ring_current_spatial_cutoff` (15 A) at frame 0 geometry. The pair
// set is fixed FOR THE TRAJECTORY; subsequent frames recompute the
// geometric channels for those same pairs (a ring that drifts past
// 15 A during the run still has its geometry emitted -- consumer
// applies their own analysis-time cutoff via the distance channel).
//
// Aromatic only -- TR filters SpatialIndexResult's full-ring result on
// `ri < protein.RingCount()`; Protein::RingCount is the aromatic-ring count.
// ProPyrrolidine excluded (no
// ring-current physics; emit via `RingPuckerTimeSeries` separately).
//
// Per-frame geometry computed FRESH from `conf.ring_geometries[ri]`
// + `conf.PositionAt(atom_idx)`. Does NOT read
// `ConformationAtom::ring_neighbours` fat-union (the 5 ring calcs
// write that struct with their own calc-specific fields; TR10 stays
// independent for cleanup-phase auditability). Does NOT touch the
// 5 ring calculators.
//
// Lifecycle: FO (Finalize-only). DenseBuffer<double>, atom-major,
// stride per atom = T * R_per_atom_max * 4. Per-atom flat offset:
//   offset(frame, r_slot, channel) = frame * R_per_atom_max * 4
//                                  + r_slot * 4
//                                  + channel
// r_slot >= per_atom_ring_count[atom] is NaN-padded.
//
// Source policy: "always_attached" -- positions present at `tp.Seed`,
// SpatialIndexResult + GeometryResult attached every frame in
// PerFrameExtractionSet. `source_attached_per_frame` emitted as all-1
// + `source_attached_policy="always_attached"` group attribute for
// SDK uniformity (per "Conditional-attach TR discipline" in
// OBJECT_MODEL).
//
// Emission at /trajectory/ring_neighbourhood_trajectory_stats/:
//   geometry                  (N, T, R_per_atom_max, 4) float64
//                             dynamic per-frame channels
//   ring_membership_per_atom  (N, R_per_atom_max)        int32
//                             static substrate-snapshot membership;
//                             -1 sentinel for unfilled slots
//   frame_indices             (T,)                       uint64
//   frame_times               (T,)                       float64 (ps)
//   source_attached_per_frame (T,)                       uint8 (all 1)
//
// Group attrs document channel_layout, units, cutoff, NaN semantics,
// aromatic-only convention, static-snapshot origin.
//
// Group is SKIPPED entirely when `r_per_atom_max_ == 0` (no
// aromatic-ring/atom pairs in the frame-0 cutoff set). Reader contract:
// group absence means "no aromatic ring neighbourhoods to emit,"
// analogous to the conditional-source skip pattern.
//
// Locked-scope provenance: project_ring_neighbourhood_debt_2026-05-21.
// Fat-union + parallel-writer pattern on `ca.ring_neighbours` stays
// as cleanup-phase target -- per `feedback_calculator_inclusion_two_use_cases`
// (cross-result-read markers across 5 calcs would have been wasteful)
// and `feedback_extractor_untouchable` (trajectory-TR landing window
// is not when to refactor 5 ring calculators).
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

class RingNeighbourhoodTrajectoryStats : public TrajectoryResult {
public:
    // Channel ordering documented above. Use this constant when
    // packing / unpacking rather than hardcoded 4.
    static constexpr std::size_t kChannelCount = 4;

    std::string Name() const override {
        return "RingNeighbourhoodTrajectoryStats";
    }

    // Declares the ConformationResult dependencies explicitly so
    // Phase 4 validates them against the active RunConfiguration's
    // `required_conf_result_types_` set. Compute reads
    // `conf.Result<SpatialIndexResult>()` (init snapshot) +
    // `conf.ring_geometries[ri]` (populated by GeometryResult).
    // PerFrameExtractionSet requires both; declaring deps here means
    // non-canonical configs fail loud at Phase 4 instead of crashing
    // at the first frame's Compute (codex round 1 2026-05-21 MED).
    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<RingNeighbourhoodTrajectoryStats> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp, Trajectory& traj) override;

    void WriteH5Group(const TrajectoryProtein& tp,
                      HighFive::File& file) const override;

    // Accessors for tests + integration probes.
    std::size_t NumAtoms()    const { return n_atoms_; }
    std::size_t NumFrames()   const { return n_frames_; }
    std::size_t RPerAtomMax() const { return r_per_atom_max_; }
    const std::vector<std::size_t>& RingListForAtom(std::size_t ai) const {
        return per_atom_ring_list_[ai];
    }

private:
    void InitStaticSnapshot_(const ProteinConformation& conf,
                              const TrajectoryProtein& tp);

    // Static per-atom (atom, aromatic-ring) snapshot from frame 0.
    // per_atom_ring_list_[atom_idx] is the sorted-ascending list of
    // aromatic ring indices within ring_current_spatial_cutoff at
    // conf0. Sized at first Compute call.
    std::vector<std::vector<std::size_t>> per_atom_ring_list_;
    std::size_t r_per_atom_max_ = 0;

    // Per-atom growing buffer. data_[atom_idx] grows by R_max * 4
    // doubles per frame; unused slots NaN-padded.
    std::vector<std::vector<double>> per_atom_data_;

    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;

    std::size_t n_atoms_  = 0;
    std::size_t n_frames_ = 0;
    bool        finalized_ = false;
};

}  // namespace nmr
