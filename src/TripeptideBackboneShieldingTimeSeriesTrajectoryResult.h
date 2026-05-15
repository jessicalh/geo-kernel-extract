#pragma once
//
// TripeptideBackboneShieldingTimeSeriesTrajectoryResult: per-atom
// per-frame time series of the tripeptide-backbone DFT shielding
// (σ_BB^i per Larsen 2015) as a SphericalTensor. Finalize-only
// dense-buffer pattern, mirrors BsShieldingTimeSeriesTrajectoryResult
// against ConformationAtom::tripeptide_bb_shielding_spherical.
//
// Per the Trajectory.cpp:160-200 comment, the per-frame
// TripeptideBackboneShieldingResult already populates the source field
// on each ConformationAtom when the [databases].tensorcs15 DSN is
// configured; this TR is the H5/NPY emission surface for the per-atom
// tensor time series.
//
// Emission pins the e3nn-compatible convention (identical layout to
// BsShieldingTimeSeriesTrajectoryResult):
//
//   /trajectory/tripeptide_bb_shielding_time_series/
//     xyz            (N, T, 9)  float64
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name    = "TripeptideBackboneShieldingTimeSeriesTrajectoryResult"
//       irrep_layout   = "T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
//       normalization  = "isometric_real_sph"
//       parity         = "0e+1o+2e"
//       units          = "ppm"
//       n_atoms, n_frames, finalized
//
// Why parity "0e+1o+2e": the tripeptide DFT tensor σ_BB^i is a
// rank-2 magnetic-shielding tensor (same shape and parity as the
// Biot-Savart kernel). Both magnetic-kernel shielding TRs share this
// parity.
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

class TripeptideBackboneShieldingTimeSeriesTrajectoryResult
    : public TrajectoryResult {
public:
    std::string Name() const override {
        return "TripeptideBackboneShieldingTimeSeriesTrajectoryResult";
    }

    // No trajectory-scope dependencies; the underlying ConformationResult
    // (TripeptideBackboneShieldingResult) is always attached when DSN is
    // configured, which is a project precondition. This TR captures
    // whatever is in tripeptide_bb_shielding_spherical each frame —
    // zero-default SphericalTensor if the calc did not attach.
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<TripeptideBackboneShieldingTimeSeriesTrajectoryResult>
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
    // (TripeptideBackboneShieldingResult) was attached this frame; 0 if
    // not. Emitted as the `source_attached_per_frame` H5 dataset for
    // downstream provenance. When all-zero (calc never ran), WriteH5Group
    // skips emission entirely per the "absent, not faked" discipline.
    std::vector<std::uint8_t> source_present_per_frame_;
    bool force_source_present_for_testing_ = false;

    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
