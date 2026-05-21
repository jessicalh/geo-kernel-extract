#pragma once
//
// MopacMcConnellShieldingTimeSeriesTrajectoryResult: per-atom per-frame
// time series of the MOPAC-weighted McConnell bond-anisotropy
// shielding contribution (ConformationAtom.mopac_mc_shielding_contribution,
// SphericalTensor, units ppm — chemical shielding from bond-order-weighted
// dipolar kernels).
//
// TR8 of the 13-TR plan. Similar shape to TR7 (sparse-cadence T2 TS)
// but emits ALL 9 components (T0+T1+T2) because the source field is
// NOT traceless. Per user direction 2026-05-21: "if not traceless
// write both" — and the methods-accumulate principle.
//
// SOURCE STRUCTURE: MopacMcConnellResult.cpp:265 sets
//   ca.mopac_mc_shielding_contribution = SphericalTensor::Decompose(M_total);
// where M_total accumulates the bond-anisotropy kernel
//   M_ab = 9*cos_theta*d_hat_a*b_hat_b - 3*b_hat_a*b_hat_b - (3*d_hat_a*d_hat_b - delta_ab)
// weighted by Wiberg bond order. M_total is NOT symmetric-traceless —
// the antisymmetric part (T1) and trace (T0) are nonzero in general.
// The per-category T2 totals (T2_backbone, T2_sidechain, T2_aromatic)
// ARE explicitly symmetrized at lines 252-262, but the overall
// shielding contribution is not. So emit all 9 components and let
// downstream readers separate T0/T1/T2 channels as needed.
//
// SPARSE CADENCE: MopacMcConnellResult attaches via TimedAttach in
// OperationRunner.cpp:185 (NOT Require). Same gate as TR5/TR6/TR7.
//
// Emission:
//   /trajectory/mopac_mc_shielding_time_series/
//     xyz            (N, T, 9)  float64 — T0, T1_m-1, T1_m0, T1_m+1,
//                                          T2_m-2, T2_m-1, T2_m0, T2_m+1,
//                                          T2_m+2  (ppm)
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64 — ps
//     source_attached_per_frame (T,) uint8
//     attrs:
//       result_name             = "MopacMcConnellShieldingTimeSeriesTrajectoryResult"
//       n_atoms, n_frames, source_attached_count, finalized
//       irrep_layout            = "T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
//       normalization           = "isometric_real_sph"
//       parity                  = "0e+1o+2e" (matches FF14SB
//                                  McConnellShieldingTimeSeries
//                                  convention — bond-anisotropy
//                                  kernel sources are non-T2 even
//                                  for symmetric polar bond directions)
//       units                   = "ppm"
//       source                  = "MopacMcConnellResult.mopac_mc_shielding_contribution
//                                  (NOT traceless; emit all 9 components per user
//                                  'if not traceless write both' 2026-05-21)"
//       source_attached_policy  = "conditional -- TimedAttach ..."
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

class MopacMcConnellShieldingTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "MopacMcConnellShieldingTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<MopacMcConnellShieldingTimeSeriesTrajectoryResult>
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
    std::size_t SourceAttachedCount() const { return source_attached_count_; }

private:
    std::vector<std::vector<SphericalTensor>> per_atom_;
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    std::size_t n_frames_              = 0;
    std::size_t source_attached_count_ = 0;
    bool        finalized_             = false;
};

}  // namespace nmr
