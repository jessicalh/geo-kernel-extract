#pragma once
//
// McConnellWelfordTrajectoryResult: running mean / variance / min / max
// of McConnell shielding tensor per atom, accumulated across all frames
// of a trajectory. AV-pattern clone of BsWelfordTrajectoryResult.
//
// Source: `ConformationAtom::mc_shielding_contribution` (units = Å⁻³,
// full McConnell asymmetric non-traceless three-term form — see
// PATTERNS.md Lesson 19; T0 = (3cos²θ-1)/r³ is meaningful). McConnellResult
// is unconditionally attached in `PerFrameExtractionSet`, so the dep
// is enforced via `Dependencies()`.
//
// Phase 2b expansion (2026-05-17): clones BS shape with one critical
// exception — McConnell-form has T1 ≡ 0 by construction (PATTERNS
// Lesson 19), so T1 channels are intentionally absent (no t1 array
// in McConnellWelfordState, no per-component T1 in Compute/Finalize,
// no t1_* H5 datasets, no irrep_layout_t1 attribute). T0 / T2[5] /
// |T2| / delta variants follow BS exactly. Per PATTERNS Lesson 25
// (Export Everything Upstream). See
// spec/plan/welford-data-shape-design-2026-05-17.md.
//

#include "TrajectoryResult.h"

#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class McConnellWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "McConnellWelfordTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<McConnellWelfordTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp, Trajectory& traj) override;

    size_t NumFrames() const { return n_frames_; }

    int WriteFeatures(const TrajectoryProtein& tp,
                      const std::string& output_dir) const override;

    void WriteH5Group(const TrajectoryProtein& tp,
                      HighFive::File& file) const override;

private:
    std::vector<double> prev_t0_;
    std::vector<bool> prev_valid_;
    size_t n_frames_ = 0;
    bool finalized_ = false;

    // Cadence metadata — mean Δt between captured frames (ps). Derived
    // at Finalize from traj.FrameTimes(). Emitted as H5 attribute so
    // downstream can convert delta-per-stride to dx/dt in physical
    // units. 0.0 before Finalize.
    double mean_dt_ps_ = 0.0;
};

}  // namespace nmr
