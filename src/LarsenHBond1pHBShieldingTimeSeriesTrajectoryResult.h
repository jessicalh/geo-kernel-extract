#pragma once
//
// LarsenHBond1pHBShieldingTimeSeriesTrajectoryResult: per-atom per-frame
// time series of the Larsen 1°HB shielding contribution
// (ConformationAtom::larsen_hbond_1pHB_spherical, SphericalTensor in
// ppm, protein lab frame). Finalize-only dense-buffer pattern;
// mirrors TripeptideBackboneShieldingTimeSeriesTrajectoryResult.
//
// The 1°HB contribution is one of four per-class fields stored
// separately for ML feature stratification — see LarsenContribDispatch
// in LarsenHBondShieldingResult.h and Larsen Table 2.
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

class LarsenHBond1pHBShieldingTimeSeriesTrajectoryResult
    : public TrajectoryResult {
public:
    std::string Name() const override {
        return "LarsenHBond1pHBShieldingTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<LarsenHBond1pHBShieldingTimeSeriesTrajectoryResult>
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
    std::vector<std::vector<SphericalTensor>> per_atom_shielding_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;
    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
