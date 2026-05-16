#pragma once
//
// SasaTimeSeriesTrajectoryResult: per-atom per-frame time series of
// the Shrake-Rupley solvent-accessible surface area
// (ConformationAtom::atom_sasa, A^2). FO dense-buffer pattern, clones
// the scalar-double shape of LarsenHBondWaterTermTimeSeries against
// an unconditional source — SasaResult is in PerFrameExtractionSet
// (RunConfiguration.cpp:134); no source-attached gate needed.
//
// Emission:
//
//   /trajectory/sasa_time_series/
//     sasa           (N, T)     float64  — Shrake-Rupley SASA (A^2)
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name    = "SasaTimeSeriesTrajectoryResult"
//       irrep_layout   = "T0"
//       parity         = "0e"
//       units          = "Angstrom^2"
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

class SasaTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "SasaTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<SasaTimeSeriesTrajectoryResult> Create(
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
    std::vector<std::vector<double>> per_atom_sasa_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;
    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
