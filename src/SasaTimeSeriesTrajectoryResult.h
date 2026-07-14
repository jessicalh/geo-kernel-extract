#pragma once
//
// SasaTimeSeriesTrajectoryResult: per-atom per-frame time series of
// the Shrake-Rupley solvent-accessible surface area
// (ConformationAtom::atom_sasa, A^2). FO dense-buffer pattern, clones
// the scalar-double shape of LarsenHBondWaterTermTimeSeries against
// an unconditional source — SasaResult is in PerFrameExtractionSet;
// no source-attached gate needed.
//
// Emission:
//
//   /trajectory/sasa_time_series/
//     sasa           (N, T)     float64  — Shrake-Rupley SASA (A^2)
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name    = "SasaTimeSeriesTrajectoryResult"
//       irrep_layout   = "raw_scalar_no_exact_o3_irrep"
//       parity         = "mixed"
//       units          = "Angstrom^2"
//       coordinate_frame = "lab_fixed_fibonacci_sampling_grid"
//       transformation = continuum scalar law qualified by the live finite,
//                        lab-fixed Fibonacci estimator's recorded envelope
//       directional_metadata_scope = "sasa dataset only"
//       n_atoms, n_frames, finalized
//
// The physical continuum SASA is rotation-invariant.  SasaResult's finite
// sampling directions are fixed in lab axes, so this captured estimator has
// no exact O(3) irrep at finite point count.
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
