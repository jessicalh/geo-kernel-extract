#pragma once
//
// LarsenHBondWaterTermTimeSeriesTrajectoryResult: per-atom per-frame
// time series of ConformationAtom::larsen_hbond_water_term (double,
// ppm). This is Larsen's isotropic Δσ_w contribution (2.07 ppm,
// NMA-water complex value) applied on amide H atoms that received
// ZERO H-bond pair contributions in this frame; zero elsewhere.
//
// Finalize-only dense-buffer pattern. Establishes the scalar-double
// (N, T) 2D H5 emission shape for the bundle's scalar TRs. Same
// (N, T) flat layout as BsT0AutocorrelationTrajectoryResult but
// with frame_indices/frame_times as the second axis instead of lag
// indices.
//
// Emission shape:
//
//   /trajectory/larsen_hbond_water_term_time_series/
//     water_term     (N, T)     float64
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name    = "LarsenHBondWaterTermTimeSeriesTrajectoryResult"
//       irrep_layout   = "T0"
//       parity         = "0e"
//       units          = "ppm"
//       description    = (Larsen water-term semantics)
//       n_atoms, n_frames, finalized
//
// Why parity "0e": this is an isotropic L0 shielding contribution
// (Larsen NMA-water complex value). A scalar in ppm, l=0 even
// parity — directly consumable as an e3nn L0 input.
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

class LarsenHBondWaterTermTimeSeriesTrajectoryResult
    : public TrajectoryResult {
public:
    std::string Name() const override {
        return "LarsenHBondWaterTermTimeSeriesTrajectoryResult";
    }

    // No trajectory-scope dependencies. The underlying
    // LarsenHBondShieldingResult is conditionally attached when the
    // Larsen H-bond grids are configured (project precondition).
    // Capture-as-is: default 0.0 means either non-HN atom or HN with
    // ≥1 H-bond pair this frame; positive 2.07 ppm means amide H with
    // no pairs.
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<LarsenHBondWaterTermTimeSeriesTrajectoryResult>
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
    // Per-atom growing buffers of double water-term scalars. Flattened
    // into an atom-major DenseBuffer<double> at Finalize.
    std::vector<std::vector<double>> per_atom_water_term_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;

    // Per-frame source-attached mask. 1 if the source ConformationResult
    // (LarsenHBondShieldingResult) was attached this frame; 0 if
    // not. Emitted as the `source_attached_per_frame` H5 dataset for
    // downstream provenance. When all-zero (calc never ran), WriteH5Group
    // skips emission entirely per the "absent, not faked" discipline.
    std::vector<std::uint8_t> source_present_per_frame_;
    bool force_source_present_for_testing_ = false;

    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
