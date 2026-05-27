#pragma once
//
// MopacCoulombShieldingTimeSeriesTrajectoryResult: per-atom per-frame
// time series of the MOPAC-Coulomb T2 EFG kernel
// (ConformationAtom.mopac_coulomb_shielding_contribution, T2
// SphericalTensor, units V/Å²). NOTE: despite the historical
// "shielding_contribution" field name, the stored value is the raw
// EFG kernel — no γ multiplication is applied at the source. The
// MopacCoulombResult source comment "gamma converts this to
// shielding" is forward-looking: shielding is recovered downstream
// by multiplying by per-element γ at calibration time.
//
// TR7 of the 13-TR plan. Combines:
//   - TR4 pattern: T2-only (N, T, 5) emission via DenseBuffer<SphericalTensor>
//   - TR5 gate:    sparse-cadence HasResult<MopacCoulombResult>() skip,
//                  NaN-fill on absent frames, source_attached_per_frame
//                  mask, group-absent when source never attached.
//
// SOURCE STRUCTURE: MopacCoulombResult stores
// `SphericalTensor::Decompose(EFG_total)`, where EFG_total is the
// trace-projected dipolar EFG from MOPAC Mulliken charges. T0 and T1
// are structurally zero, so (N, T, 5) emission is information-preserving.
//
// SPARSE CADENCE: MopacCoulombResult attaches via TimedAttach, not
// RequireConformationResult. Same CLI-driven Mopac cadence as TR5/TR6.
// WriteH5Group skips the entire
// /trajectory/mopac_coulomb_shielding_time_series/ group when no
// frame attached the source.
//
// Emission:
//   /trajectory/mopac_coulomb_shielding_time_series/
//     t2             (N, T, 5)  float64 — T2_m-2..T2_m+2 (V/Å²)
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64 — ps
//     source_attached_per_frame (T,) uint8
//     attrs:
//       result_name             = "MopacCoulombShieldingTimeSeriesTrajectoryResult"
//       n_atoms, n_frames, source_attached_count, finalized
//       irrep_layout            = "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
//       normalization           = "isometric_real_sph"
//       parity                  = "2e"
//       units                   = "V/Å^2"  (EFG kernel, pre-γ)
//       source                  = "MopacCoulombResult.mopac_coulomb_shielding_contribution
//                                  (T2-only per source comment 'Pure T2
//                                  (EFG is traceless)')"
//       source_attached_policy  = "conditional -- MopacCoulombResult attaches
//                                  sparsely per the Mopac cadence ..."
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

class MopacCoulombShieldingTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "MopacCoulombShieldingTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<MopacCoulombShieldingTimeSeriesTrajectoryResult>
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
    // Per-atom growing buffer of SphericalTensor. Transferred to
    // DenseBuffer<SphericalTensor> at Finalize.
    std::vector<std::vector<SphericalTensor>> per_atom_;
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    std::size_t n_frames_              = 0;
    std::size_t source_attached_count_ = 0;
    bool        finalized_             = false;
};

}  // namespace nmr
