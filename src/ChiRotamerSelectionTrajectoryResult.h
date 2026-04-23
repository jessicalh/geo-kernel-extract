#pragma once
//
// ChiRotamerSelectionTrajectoryResult: per-residue chi-angle rotamer
// transition detector. Scan-mode selection emitter worked example.
//
// On each frame, for every residue whose chi[k] angle atom indices
// are valid (set by Protein::FinalizeConstruction at load time),
// compute the dihedral chi_k from positions, classify into a rotamer
// bin (g+, t, g-), and compare to the bin observed on the prior
// frame for that residue + chi index. When the bin changes, push a
// SelectionRecord onto the run-scope bag via traj.MutableSelections().
//
// Every future scan-mode selection emitter follows this shape:
//
//   - Internal state: per-residue (or per-atom, per-bond) previous
//     bin/value + valid flag. Owned by the Result; never stashed on
//     TrajectoryAtom.
//   - Compute reads conf, compares to prior state, pushes records
//     directly to traj.MutableSelections() on transitions. The push
//     carries {kind=typeid(this), frame_idx, time_ps, reason,
//     metadata} — frame and time explicit; no selector string.
//   - Finalize optional (none here — transitions are emitted per
//     frame; no end-of-stream synthesis).
//   - WriteH5Group optional — the SelectionBag's own writer emits
//     /trajectory/selections/<kind>/ with this Result's records.
//
// No ConformationResult dependency: dihedral computation uses only
// positions + Residue.chi[k] atom indices (set once at FinalizeConstruction).
//

#include "TrajectoryResult.h"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class ChiRotamerSelectionTrajectoryResult : public TrajectoryResult {
public:
    // Three-bin rotamer classification. Gauche+, trans, gauche-.
    enum class RotamerBin : int { Gplus = 0, Trans = 1, Gminus = 2 };

    std::string Name() const override {
        return "ChiRotamerSelectionTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override {
        return {};
    }

    static std::unique_ptr<ChiRotamerSelectionTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    // Running count of transitions pushed — exposed for tests and
    // diagnostics. The authoritative record is the SelectionBag.
    std::size_t TransitionCount() const { return n_transitions_; }

private:
    // Per-residue previous bin per chi index. bin_valid_[r][k] == false
    // until the first observation (first frame skips — a transition
    // needs a prior to compare against).
    std::vector<std::array<RotamerBin, 4>> prev_bin_;
    std::vector<std::array<bool, 4>> bin_valid_;

    std::size_t n_transitions_ = 0;
};

}  // namespace nmr
