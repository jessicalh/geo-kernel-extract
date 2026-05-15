#pragma once
//
// TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult: per-atom
// per-frame time series of the i+1-direction post-Kabsch positional
// residual produced by TripeptideNeighborShieldingResult
// (Vec3, units Å, protein lab frame). Finalize-only dense-buffer
// pattern; mirrors the prev-direction TR against the next field.
//
// **NaN-fill contract:** TripeptideNeighborShieldingResult.cpp:162-169
// NaN-initialises tripeptide_neighbor_residual_vec_next at the top of
// each per-frame Compute. NaN means "no i+1 neighbour contribution
// this frame"; finite means "i+1 tripeptide aligned and stamped a
// residual." NaN is signal; this TR captures the field as-is.
//
// Emission shape:
//
//   /trajectory/tripeptide_neighbor_residual_vec_next_time_series/
//     xyz            (N, T, 3)  float64
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name    = "TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult"
//       irrep_layout   = "x,y,z"
//       normalization  = "cartesian"
//       parity         = "1o"
//       units          = "angstrom"
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

class TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult
    : public TrajectoryResult {
public:
    std::string Name() const override {
        return "TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult>
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

private:
    std::vector<std::vector<Vec3>> per_atom_residual_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;
    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
