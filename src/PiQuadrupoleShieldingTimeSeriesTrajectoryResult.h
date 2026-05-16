#pragma once
//
// PiQuadrupoleShieldingTimeSeriesTrajectoryResult: per-atom per-frame
// time series of the pi-electron quadrupole geometric kernel as a
// SphericalTensor (Angstrom^-5). FO dense-buffer pattern, clones
// BsShieldingTimeSeriesTrajectoryResult against
// ConformationAtom::piquad_shielding_contribution. Note: despite the
// field name, the source calc (PiQuadrupoleResult) stores the
// decomposed GEOMETRIC KERNEL, not the parameterised ppm shielding.
// Pre-existing naming drift across most classical calcs — see
// OBJECT_MODEL.md "Calculator Shielding Contribution Contract drift".
// Units: G tensor is Å⁻⁵ per PiQuadrupoleResult.cpp:43 (the related
// Buckingham scalar A-term is Å⁻⁴; do not confuse).
//
// **T0 is structurally zero.** The EFG kernel is analytically traceless
// by Laplace's equation (Stone, *Theory of Intermolecular Forces*,
// OUP 2013, Ch. 3; see PiQuadrupoleResult.cpp:40-44, asserted at
// tests/test_pi_quadrupole_result.cpp:475-476). All physical signal
// lives in T2. Downstream consumers reading this TR's H5 output
// should compute statistics on the 5 T2 components, not T0.
//
// Source is PiQuadrupoleResult, unconditionally attached in
// PerFrameExtractionSet. No source-attached gate.
//
// Emission:
//
//   /trajectory/piquad_shielding_time_series/
//     xyz            (N, T, 9)  float64
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name    = "PiQuadrupoleShieldingTimeSeriesTrajectoryResult"
//       irrep_layout   = "T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
//       normalization  = "isometric_real_sph"
//       parity         = "0e+1o+2e"
//       units          = "Angstrom^-5"
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

class PiQuadrupoleShieldingTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "PiQuadrupoleShieldingTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<PiQuadrupoleShieldingTimeSeriesTrajectoryResult> Create(
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
