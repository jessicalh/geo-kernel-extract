#pragma once
//
// BondLengthStatsTrajectoryResult: per-bond length mean/std/min/max
// + frame-to-frame length delta stats, across all frames of a run.
//
// Bond-scope TrajectoryResult. Per-bond accumulator state is
// INTERNAL TO THIS RESULT as `std::vector<PerBondWelford>`, not a
// first-class TrajectoryBond peer-of-TrajectoryAtom store.
//
// Bonds are identified by index into `tp.ProteinRef().Bonds()`. The
// Result sizes `per_bond_` from `protein.BondCount()` at Create —
// valid because Trajectory::Run Phase 2 (Seed) precedes Phase 3
// (factory invocation), so the Protein is finalised.
//
// Emission: /trajectory/bond_length_stats/
//   length_mean     (B,) float64
//   length_std      (B,)
//   length_min      (B,)
//   length_max      (B,)
//   length_delta_mean (B,)
//   length_delta_std  (B,)
//   atom_a          (B,) uint64  — topology passthrough
//   atom_b          (B,) uint64
//   order           (B,) int8    — BondOrder enum value
//   category        (B,) int8    — BondCategory enum value
//   attrs: result_name, n_bonds, n_frames, finalized, units="Å"
//
// Topology metadata (atom_a, atom_b, order, category) is frame-
// invariant and technically belongs in a future /topology/bonds/
// group emitted once by TrajectoryProtein::WriteH5. Emitting it here
// keeps this group self-contained for downstream consumers; when the
// canonical /topology/ group lands, bond-scope Results can drop the
// duplication.
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class BondLengthStatsTrajectoryResult : public TrajectoryResult {
public:
    struct PerBondWelford {
        double length_mean = 0.0;
        double length_m2 = 0.0;
        double length_std = 0.0;   // Finalize-only
        double length_min = std::numeric_limits<double>::infinity();
        double length_max = -std::numeric_limits<double>::infinity();
        std::size_t n_frames = 0;

        double delta_mean = 0.0;
        double delta_m2 = 0.0;
        double delta_std = 0.0;    // Finalize-only
        double prev_length = 0.0;
        bool has_prev = false;
        std::size_t delta_n = 0;
    };

    std::string Name() const override {
        return "BondLengthStatsTrajectoryResult";
    }

    // No ConformationResult dependency — reads positions directly.
    std::vector<std::type_index> Dependencies() const override {
        return {};
    }

    static std::unique_ptr<BondLengthStatsTrajectoryResult> Create(
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
    const std::vector<PerBondWelford>& PerBond() const { return per_bond_; }

private:
    std::vector<PerBondWelford> per_bond_;
    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
