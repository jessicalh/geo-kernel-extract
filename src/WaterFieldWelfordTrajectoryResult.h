#pragma once
//
// WaterFieldWelfordTrajectoryResult: per-atom Welford rollup of the
// explicit-water E-field + EFG kernel from WaterFieldResult. Clones the
// BS Welford 5-part shape: per-component Welford on vector / tensor
// channels, plus signed/abs/squared/dxdt delta variants on the primary
// scalar channels (E-field magnitude, EFG T0, shell-occupancy counts).
//
// Channels (mirroring WaterFieldTimeSeries source fields):
//
//   efield                Vec3 V/Å              per-component[3] + magnitude
//   efield_first          Vec3 V/Å              per-component[3] + magnitude
//   efg                   SphericalTensor V/Å²  T0 + T1[3] + T2[5] + |T2|
//   efg_first             SphericalTensor V/Å²  T0 + T1[3] + T2[5] + |T2|
//   n_first               int (dimensionless)   full scalar Welford
//   n_second              int (dimensionless)   full scalar Welford
//
// Delta variants on: efield_magnitude, efg_t0, n_first, n_second.
// These four scalars carry the per-channel-distinct dynamics question.
//
// Emission:
//
//   /trajectory/water_field_welford/
//     efield_{x,y,z}_{mean,m2,std,min,max,min_frame,max_frame}     (N,)
//     efield_magnitude_{mean,...}                                  (N,)
//     efield_first_{x,y,z}_{...}                                   (N,)
//     efield_first_magnitude_{...}                                 (N,)
//     efg_t0_{...}                                                 (N,)
//     efg_t1_{...}                                                 (N, 3)
//     efg_t2_{...}                                                 (N, 5)
//     efg_t2magnitude_{...}                                        (N,)
//     efg_first_{t0,t1,t2,t2magnitude}_{...}                       same shapes
//     n_first_{...}, n_second_{...}                                (N,)
//     <scalar>_delta_{...}, <scalar>_abs_delta_{...},
//     <scalar>_delta_squared_{...}, <scalar>_dxdt_{...}           (N,) for the 4 scalars
//     n_frames_per_atom, delta_n_per_atom, dxdt_n_per_atom        (N,)
//     attrs: result_name, n_frames, finalized, ddof, mean_dt_ps,
//            frame_index_range, irrep_layout_t1, irrep_layout_t2
//
// Dependencies: WaterFieldResult (unconditionally attached in
// PerFrameExtractionSet — no source-attached gate).
//

#include "TrajectoryResult.h"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class WaterFieldWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "WaterFieldWelfordTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<WaterFieldWelfordTrajectoryResult> Create(
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
    // Per-atom previous-frame caches for the delta trackers. No
    // prev_efg_t0_ cache — `efg_t0` channel removed per 2026-05-18
    // adversarial review (structurally zero after the traceless
    // projection in `WaterFieldResult.cpp:147-150`).
    std::vector<double>      prev_efield_mag_;
    std::vector<double>      prev_n_first_;
    std::vector<double>      prev_n_second_;
    std::vector<bool>        prev_valid_;
    std::vector<double>      prev_time_;

    std::size_t              n_frames_  = 0;
    bool                     finalized_ = false;
    double                   mean_dt_ps_ = 0.0;
    std::array<std::size_t, 2> frame_index_range_ = {0, 0};
};

}  // namespace nmr
