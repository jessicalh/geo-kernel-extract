// BsWelfordTrajectoryResult.h
//
// Running mean, variance, min, max of the BiotSavart T0 shielding
// contribution per atom, accumulated across all frames of a trajectory.
//
// The worked example from WIP_OBJECT_MODEL.md §4. Present here so the
// framework sketch has at least one concrete TrajectoryResult showing
// the pattern.
//
// Source: each frame, after BiotSavartResult::Compute has run on the
// frame's ProteinConformation, this result reads
// conf.atoms[i].bs_shielding_contribution.T0 and ticks a Welford update
// per atom. Running mean / M2 / min / max live as fields on
// TrajectoryAtom; they are always valid mid-stream.
//
// Finalize divides M2 by (n-1) and writes the final std to
// TrajectoryAtom.bs_t0_std.

#ifndef TRAJECTORY_FRAMEWORK_SKETCH_BS_WELFORD_TRAJECTORY_RESULT_H
#define TRAJECTORY_FRAMEWORK_SKETCH_BS_WELFORD_TRAJECTORY_RESULT_H

#include <cstddef>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

#include "Stubs.h"
#include "TrajectoryResult.h"

class TrajectoryProtein;

class BsWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override { return "BsWelfordTrajectoryResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: create empty with running Welford state zeroed.
    static std::unique_ptr<BsWelfordTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 std::size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp) override;

    // Post-finalize convenience query methods (§4 example).
    double MeanAtAtom(std::size_t atom_idx) const;
    double StdAtAtom(std::size_t atom_idx) const;  // valid only after Finalize
    std::size_t NumFrames() const { return n_frames_; }

private:
    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

#endif // TRAJECTORY_FRAMEWORK_SKETCH_BS_WELFORD_TRAJECTORY_RESULT_H
