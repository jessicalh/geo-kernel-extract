#pragma once
//
// AIMNet2ChargeTimeSeriesTrajectoryResult: per-atom per-frame time
// series of the AIMNet2 neural-network Hirshfeld charge
// (ConformationAtom::aimnet2_charge, elementary charge). FO
// dense-buffer pattern. AIMNet2Result is in PerFrameExtractionSet
// (RunConfiguration.cpp:136) and requires a Session-loaded model;
// when the model is not loaded, OperationRunner aborts the run
// before any TR Compute is called — no per-frame conditional gate
// at the TR layer.
//
// Emission:
//
//   /trajectory/aimnet2_charge_time_series/
//     charge         (N, T)     float64  — Hirshfeld charge (e)
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name    = "AIMNet2ChargeTimeSeriesTrajectoryResult"
//       irrep_layout   = "T0"
//       parity         = "0e"
//       units          = "elementary_charge"
//       n_atoms, n_frames, finalized
//

#include "DenseBuffer.h"
#include "TrajectoryResult.h"
#include "Types.h"

#include <cstddef>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class AIMNet2ChargeTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "AIMNet2ChargeTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<AIMNet2ChargeTimeSeriesTrajectoryResult> Create(
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

private:
    std::vector<std::vector<double>> per_atom_charge_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;
    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
