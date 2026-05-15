#pragma once
//
// LarsenHBondCountTimeSeriesTrajectoryResult: per-atom per-frame time
// series of ConformationAtom::larsen_hbond_n_pairs (int) — the number
// of H-bond pair contributions this atom received this frame, summed
// across all four Larsen classes (1°HB + 2°HB + 1°HαB + 2°HαB).
// Finalize-only dense-buffer pattern; (N, T) int32 emission.
//
// Emission shape:
//
//   /trajectory/larsen_hbond_count_time_series/
//     count          (N, T)     int32
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name = "LarsenHBondCountTimeSeriesTrajectoryResult"
//       units       = "pairs"
//       dtype       = "int32"
//       description = "Per-atom H-bond pair count (sum over 1°HB, 2°HB, 1°HαB, 2°HαB)"
//       n_atoms, n_frames, finalized
//
// No parity/irrep_layout attrs: this is a discrete count, not a
// tensor irrep. Downstream ML consumers may pair it with the per-class
// shielding TRs to feature-engineer "pair density per donor class."
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

class LarsenHBondCountTimeSeriesTrajectoryResult
    : public TrajectoryResult {
public:
    std::string Name() const override {
        return "LarsenHBondCountTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<LarsenHBondCountTimeSeriesTrajectoryResult>
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

    // Test-only: bypass the per-frame `conf.HasResult<SourceCalc>()`
    // check inside Compute. Call BEFORE the first Compute. Use ONLY
    // in synthetic unit tests that don't go through Trajectory::Run /
    // OperationRunner (the production path that attaches the source
    // ConformationResult). Production code never calls this.
    void ForceSourcePresentForTesting() {
        force_source_present_for_testing_ = true;
    }

private:
    std::vector<std::vector<int>> per_atom_count_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;

    // Per-frame source-attached mask. 1 if the source ConformationResult
    // (LarsenHBondShieldingResult) was attached this frame; 0 if
    // not. Emitted as the `source_attached_per_frame` H5 dataset for
    // downstream provenance. When all-zero (calc never ran), WriteH5Group
    // skips emission entirely per the "absent, not faked" discipline.
    std::vector<std::uint8_t> source_present_per_frame_;
    bool force_source_present_for_testing_ = false;

    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
