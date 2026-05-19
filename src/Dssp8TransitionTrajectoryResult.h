#pragma once
//
// Dssp8TransitionTrajectoryResult: per-residue DSSP 8-state transition
// statistics. AV companion to `Dssp8TimeSeriesTrajectoryResult`.
//
// Tracks how the secondary structure of each residue moves around the
// 8-state alphabet (H/G/I/E/B/T/S/C) over the trajectory: full
// transition-count matrix (R, 8, 8), occupancy fractions (R, 8),
// dominant state (R,), and the count of frames where SS was observed.
//
// Same conditional-source pattern as Dssp8TimeSeries: depends on
// DsspResult attaching per frame; source-attached gate + skip-emission
// when source never attached.
//
// Pure AV: Compute updates running counters in place, no DenseBuffer.
// Per-residue state internal to the TR (same as
// DihedralBinTransitionTrajectoryResult).
//
// Emission: /trajectory/dssp8_transition/
//   Per-residue per-frame (none -- this is pure AV; per-frame SS
//   labels live in Dssp8TimeSeriesTrajectoryResult).
//
//   Per-residue stats (R,):
//     ss8_transition_count   uint32  total transitions (any-to-any,
//                                    excluding self-edges and including
//                                    only consecutive observed pairs)
//     ss8_dominant           uint8   argmax of ss8_occupancy
//                                    (255 if no observation)
//     n_frames_observed      uint32  (frames where source attached AND
//                                    DSSP reported a code != unassigned)
//
//   Per-residue per state (R, 8):
//     ss8_occupancy          uint32  frame count per state
//
//   Per-residue transition matrix (R, 8, 8):
//     ss8_transition_matrix  uint32  count of (prev=row, curr=col)
//                                    transitions, prev != curr
//
//   Per-frame metadata:
//     frame_indices, frame_times, source_attached_per_frame
//
// Transition gate: BOTH prev and curr observed (non-unassigned) for
// the consecutive pair to count. Diagonals are EXCLUDED from
// ss8_transition_count and ss8_transition_matrix (no-transition
// frames bump occupancy only).
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

class Dssp8TransitionTrajectoryResult : public TrajectoryResult {
public:
    static constexpr std::size_t kSSCount = 8;  // H/G/I/E/B/T/S/C
    static constexpr std::uint8_t kSSUnassigned = 255;

    std::string Name() const override {
        return "Dssp8TransitionTrajectoryResult";
    }

    // Same as Dssp8TS: DsspResult is conditionally attached, no declared
    // dependency.
    std::vector<std::type_index> Dependencies() const override {
        return {};
    }

    static std::unique_ptr<Dssp8TransitionTrajectoryResult>
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

    // Test-only: bypass the per-frame `conf.HasResult<DsspResult>()`
    // check. SAFETY (codex review 2026-05-19): When this flag is true,
    // Compute calls `conf.Result<DsspResult>()` which throws if the
    // Result is not actually attached. Test caller MUST attach
    // DsspResult to every ProteinConformation before setting this flag.
    void ForceSourcePresentForTesting() {
        force_source_present_for_testing_ = true;
    }

private:
    std::vector<std::uint8_t> prev_ss_;  // kSSUnassigned = no prev obs

    std::vector<std::uint32_t> ss8_transition_count_;
    std::vector<std::array<std::uint32_t, kSSCount>> ss8_occupancy_;
    std::vector<
        std::array<std::array<std::uint32_t, kSSCount>, kSSCount>>
        ss8_transition_matrix_;
    std::vector<std::uint32_t> n_frames_observed_;

    // Finalize-derived.
    std::vector<std::uint8_t> ss8_dominant_;

    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    bool force_source_present_for_testing_ = false;

    std::size_t n_frames_  = 0;
    bool        finalized_ = false;
};

}  // namespace nmr
