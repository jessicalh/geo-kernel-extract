#pragma once
//
// TripeptideNeighborShieldingTimeSeriesTrajectoryResult: per-atom
// per-frame time series of the TripeptideNeighborShieldingResult
// Δσ_BB^{i±1} contribution as a SphericalTensor. Finalize-only
// dense-buffer pattern, mirrors BsShieldingTimeSeriesTrajectoryResult
// against the tripeptide_neighbor_shielding_spherical source field
// on ConformationAtom.
//
// Per Larsen 2015 Eq 3 the per-residue contribution
// Δσ_BB^{i-1}(i) + Δσ_BB^{i+1}(i) is read at the flanking ALA cap
// atoms of the (i±1)-centered tripeptides, with the AAA reference
// at standard angles (φ_std=-120°, ψ_std=140°) subtracted; the per-
// frame ConformationResult writes the SUM onto each central atom's
// tripeptide_neighbor_shielding_spherical field. This TR captures
// that field across the trajectory.
//
// Emission pins the e3nn-compatible convention:
//
//   /trajectory/tripeptide_neighbor_shielding_time_series/
//     xyz            (N, T, 9)  float64
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name    = "TripeptideNeighborShieldingTimeSeriesTrajectoryResult"
//       irrep_layout   = "T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
//       normalization  = "isometric_real_sph"
//       parity         = "0e+1o+2e"
//       units          = "ppm"
//       n_atoms, n_frames, finalized
//
// Parity "0e+1o+2e": the tripeptide neighbour Δσ tensor is a DFT-
// derived shielding tensor restamped into the central-atom frame
// by the cap-side Kabsch alignment of the (i±1)-centered tripeptide
// onto residue i's N/CA/C backbone. The underlying object is a
// magnetic shielding tensor (its T1 part is a parity-odd
// pseudovector in the e3nn Irreps convention), so the parity matches
// the BS / HM / McConnell magnetic-kernel TimeSeries TRs. Coulomb /
// APBS EFG (symmetric traceless, T1=0) use parity "0e+1e+2e".
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

class TripeptideNeighborShieldingTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "TripeptideNeighborShieldingTimeSeriesTrajectoryResult";
    }

    // No trajectory-scope dependencies; the underlying ConformationResult
    // (TripeptideNeighborShieldingResult) is always attached when DSN is
    // configured (project precondition). This TR captures whatever is in
    // tripeptide_neighbor_shielding_spherical each frame — zero-default
    // SphericalTensor if the calc did not attach.
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<TripeptideNeighborShieldingTimeSeriesTrajectoryResult> Create(
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

    // Test-only: bypass the per-frame `conf.HasResult<SourceCalc>()`
    // check inside Compute. Call BEFORE the first Compute. Use ONLY
    // in synthetic unit tests that don't go through Trajectory::Run /
    // OperationRunner (the production path that attaches the source
    // ConformationResult). Production code never calls this.
    void ForceSourcePresentForTesting() {
        force_source_present_for_testing_ = true;
    }

private:
    // Per-atom growing buffers of SphericalTensor. Flattened into an
    // atom-major DenseBuffer<SphericalTensor> at Finalize.
    std::vector<std::vector<SphericalTensor>> per_atom_shielding_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;

    // Per-frame source-attached mask. 1 if the source ConformationResult
    // (TripeptideNeighborShieldingResult) was attached this frame; 0 if
    // not. Emitted as the `source_attached_per_frame` H5 dataset for
    // downstream provenance. When all-zero (calc never ran), WriteH5Group
    // skips emission entirely per the "absent, not faked" discipline.
    std::vector<std::uint8_t> source_present_per_frame_;
    bool force_source_present_for_testing_ = false;

    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
