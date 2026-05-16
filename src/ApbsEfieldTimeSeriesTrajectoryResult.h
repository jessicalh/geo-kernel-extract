#pragma once
//
// ApbsEfieldTimeSeriesTrajectoryResult: per-atom per-frame time
// series of the APBS solvated E-field
// (ConformationAtom::apbs_efield, V/Angstrom). FO dense-buffer
// pattern with Vec3 payload (clones PositionsTimeSeriesTrajectoryResult
// shape against the apbs_efield field). ApbsFieldResult is
// unconditionally attached in PerFrameExtractionSet
// (RunConfiguration.cpp:126); no source-attached gate.
//
// Emission:
//
//   /trajectory/apbs_efield_time_series/
//     xyz            (N, T, 3)  float64  — E_x, E_y, E_z (V/A)
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name    = "ApbsEfieldTimeSeriesTrajectoryResult"
//       irrep_layout   = "x,y,z"
//       normalization  = "cartesian"
//       parity         = "1o"
//       units          = "V/Angstrom"
//       n_atoms, n_frames, finalized
//
// Parity "1o": E-field is a polar (true) vector, parity-odd under
// inversion. Same parity convention as PositionsTimeSeries.
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

class ApbsEfieldTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "ApbsEfieldTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<ApbsEfieldTimeSeriesTrajectoryResult> Create(
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
    std::vector<std::vector<Vec3>> per_atom_efield_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;
    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
