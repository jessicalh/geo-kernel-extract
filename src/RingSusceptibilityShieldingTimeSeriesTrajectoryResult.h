#pragma once
//
// RingSusceptibilityShieldingTimeSeriesTrajectoryResult: per-atom
// per-frame time series of the ring-susceptibility geometric kernel
// as a SphericalTensor (Angstrom^-3). FO dense-buffer pattern,
// clones BsShieldingTimeSeriesTrajectoryResult against
// ConformationAtom::ringchi_shielding_contribution. Despite the field
// name, the source calc stores the decomposed GEOMETRIC KERNEL, not
// parameterised ppm shielding.
//
// Source is RingSusceptibilityResult, unconditionally attached in
// PerFrameExtractionSet. No source-attached gate.
//
// Emission:
//
//   /trajectory/ringchi_shielding_time_series/
//     xyz            (N, T, 9)  float64
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name    = "RingSusceptibilityShieldingTimeSeriesTrajectoryResult"
//       irrep_layout   = "T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
//       normalization  = "isometric_real_sph"
//       parity         = "0e+1o+2e"
//       units          = "Angstrom^-3"
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

class RingSusceptibilityShieldingTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "RingSusceptibilityShieldingTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<RingSusceptibilityShieldingTimeSeriesTrajectoryResult> Create(
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
    std::vector<std::vector<SphericalTensor>> per_atom_shielding_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;
    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
