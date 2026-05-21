#pragma once
//
// RmsdTrackingTrajectoryResult: AV per-frame backbone-heavy-atom RMSD
// vs trajectory frame 0. TR11 of the 13-TR plan; first scalar AV (T,)
// double TR in the codebase.
//
// Reference frame: frame 0 of the trajectory (captured at first
// Compute call). Reference geometry is the position of the backbone
// heavy atoms (N, CA, C, O) for every residue where each is present.
// ACE/NME caps contribute their own backbone atoms when present.
//
// RMSD definition: Kabsch-aligned RMSD over the selected backbone atom
// set. Per frame:
//   1. Collect current positions of atom_indices_
//   2. Compute optimal rotation R via SVD of cross-covariance
//      against frame-0 reference (translation removed by centroid
//      subtraction).
//   3. Apply rotation to current positions, compute RMSD vs reference.
//
// Lifecycle: AV (always valid mid-stream). The per-frame rmsd_ vector
// is appended in Compute; consumers can read rmsd_[k] at any time
// after frame k's Compute. Finalize is trivial.
//
// Source policy: "always_attached" -- positions present at tp.Seed.
// `source_attached_per_frame` emitted as all-1 + group attr for SDK
// uniformity (per Conditional-attach TR discipline).
//
// Emission at /trajectory/rmsd_tracking/:
//   rmsd                        (T,)   float64 Angstrom
//   atom_indices                (M,)   int32   protein atom indices
//                                              used for alignment
//   frame_indices               (T,)   uint64
//   frame_times                 (T,)   float64 ps
//   source_attached_per_frame   (T,)   uint8 (all 1)
//
// Attrs:
//   result_name             = "RmsdTrackingTrajectoryResult"
//   n_atoms                 = M (alignment set size)
//   n_frames                = T
//   finalized               = bool
//   alignment_method        = "kabsch_svd"
//   atom_selection          = "backbone_heavy_atoms_NCACO"
//   reference_frame_origin  = "trajectory_frame_0"
//   units                   = "Angstrom"
//   source_attached_policy  = "always_attached"
//
// Pairs with:
//   - TR12 RmsdSpikeSelectionTrajectoryResult, which CROSS-RESULT READs
//     this TR's per-frame rmsd to detect spikes (thresholds: 2.0A vs
//     frame 0 + 0.5A vs rolling 100-frame mean; cooldown 100 frames).
//
// CROSS-RESULT READ (writer side):
//   Fields read by other TRs during their Compute:
//     - `RmsdAtFrame(frame_idx)` const accessor returns the AV
//       per-frame RMSD scalar (Angstrom). Read by
//       RmsdSpikeSelectionTrajectoryResult.

#include "TrajectoryResult.h"
#include "Types.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class RmsdTrackingTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "RmsdTrackingTrajectoryResult";
    }

    // No TR dependency. Reads positions only; AtomIndices captured at
    // Seed time from typed Residue cache (Residue.N / CA / C / O).
    std::vector<std::type_index> Dependencies() const override {
        return {};
    }

    static std::unique_ptr<RmsdTrackingTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp, Trajectory& traj) override;

    void WriteH5Group(const TrajectoryProtein& tp,
                      HighFive::File& file) const override;

    // CROSS-RESULT READ (writer side): read by
    // RmsdSpikeSelectionTrajectoryResult during its own Compute.
    // Returns AV (immediately-valid mid-stream) per-frame RMSD scalar
    // in Angstroms. Returns NaN if frame_idx >= NumFrames() (frame
    // not yet computed in the current per-frame dispatch order).
    double RmsdAtFrame(std::size_t frame_idx) const;
    std::size_t NumFrames() const { return n_frames_; }
    std::size_t NumAlignmentAtoms() const { return atom_indices_.size(); }
    const std::vector<double>& Rmsd() const { return rmsd_; }

private:
    // Per-residue backbone atom indices captured at Create time from
    // Residue.N / CA / C / O typed slots. Atoms with NONE-marked slots
    // are skipped (ACE: no N/CA; NME: no C/O; etc.). Backbone-heavy-
    // only (no hydrogens).
    std::vector<std::size_t> atom_indices_;

    // Reference positions for the alignment set, captured at first
    // Compute call from frame-0 geometry.
    std::vector<Vec3> reference_positions_;
    bool reference_captured_ = false;

    // Per-frame Kabsch-aligned RMSD (Angstrom). Appended in Compute.
    std::vector<double>       rmsd_;
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;

    std::size_t n_frames_  = 0;
    bool        finalized_ = false;
};

}  // namespace nmr
