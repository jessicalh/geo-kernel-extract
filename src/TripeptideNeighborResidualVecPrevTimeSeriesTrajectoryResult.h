#pragma once
//
// TripeptideNeighborResidualVecPrevTimeSeriesTrajectoryResult: per-atom
// per-frame time series of the i-1-direction post-Kabsch positional
// residual produced by TripeptideNeighborShieldingResult
// (Vec3, units Å, protein lab frame). Finalize-only dense-buffer
// pattern; mirrors TripeptideBackboneResidualVecTimeSeriesTrajectoryResult
// against the neighbor prev field.
//
// **NaN-fill contract:** The producing calculator
// (TripeptideNeighborShieldingResult.cpp:162-169) explicitly
// NaN-initialises tripeptide_neighbor_residual_vec_prev at the top of
// each per-frame Compute. A NaN cell means "no i-1 neighbour
// contribution this frame for this atom" — finite cells mean "i-1
// tripeptide aligned and stamped a residual here." NaN is signal,
// not corruption; downstream consumers distinguish via isnan
// rather than against a sentinel value. This TR captures the field
// as-is; NaN propagates through the DenseBuffer<Vec3> and the H5
// dataset unchanged.
//
// Emission shape:
//
//   /trajectory/tripeptide_neighbor_residual_vec_prev_time_series/
//     xyz            (N, T, 3)  float64
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name    = "TripeptideNeighborResidualVecPrevTimeSeriesTrajectoryResult"
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

class TripeptideNeighborResidualVecPrevTimeSeriesTrajectoryResult
    : public TrajectoryResult {
public:
    std::string Name() const override {
        return "TripeptideNeighborResidualVecPrevTimeSeriesTrajectoryResult";
    }

    // No trajectory-scope dependencies; the underlying
    // TripeptideNeighborShieldingResult is conditionally attached when
    // the [databases].tensorcs15 DSN is configured. Captures whatever
    // is in tripeptide_neighbor_residual_vec_prev each frame — NaN
    // means no i-1 contribution that frame.
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<TripeptideNeighborResidualVecPrevTimeSeriesTrajectoryResult>
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
