#pragma once
//
// BsShieldingTimeSeriesTrajectoryResult: per-atom per-frame time
// series of the BiotSavart shielding contribution as a
// SphericalTensor. Finalize-only dense-buffer pattern, same as
// PositionsTimeSeriesTrajectoryResult but with a tensor payload.
//
// Canonical worked example for SphericalTensor time-series emission.
// Every *ShieldingTimeSeriesTrajectoryResult that follows (McConnell,
// HaighMallion, RingSusceptibility, PiQuadrupole, Dispersion, HBond,
// Coulomb, APBS, AIMNet2) clones this shape against its own source
// field on ConformationAtom and its own calculator-specific parity.
//
// Emission pins the e3nn-compatible convention:
//
//   /trajectory/bs_shielding_time_series/
//     xyz            (N, T, 9)  float64
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name    = "BsShieldingTimeSeriesTrajectoryResult"
//       irrep_layout   = "T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
//       normalization  = "isometric_real_sph"
//       parity         = "0e+1o+2e"
//       units          = "ppm"
//       n_atoms, n_frames, finalized
//
// Why parity "0e+1o+2e" for BiotSavart: the shielding kernel
// G_ab = -n_b · B_a · PPM_FACTOR is the outer product of a ring
// normal with a magnetic field. B is an axial (pseudo) vector, n is
// an axial vector, but the rank-2 shielding tensor's antisymmetric
// part is a parity-odd pseudovector in the e3nn Irreps convention
// (1o). Every future magnetic-kernel shielding TR has the same
// parity; Coulomb / APBS EFG (which are symmetric traceless, T1=0)
// use parity "0e+1e+2e".
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

class BsShieldingTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "BsShieldingTimeSeriesTrajectoryResult";
    }

    // Per-frame dependency: BiotSavartResult must run so that
    // ConformationAtom::bs_shielding_contribution is populated before
    // this result reads it.
    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<BsShieldingTimeSeriesTrajectoryResult> Create(
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
    // Per-atom growing buffers of SphericalTensor. Flattened into an
    // atom-major DenseBuffer<SphericalTensor> at Finalize.
    std::vector<std::vector<SphericalTensor>> per_atom_shielding_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;
    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
