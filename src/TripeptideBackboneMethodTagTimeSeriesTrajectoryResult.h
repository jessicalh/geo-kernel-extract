#pragma once
//
// TripeptideBackboneMethodTagTimeSeriesTrajectoryResult: per-atom
// per-frame time series of the tripeptide-backbone method-tag carried
// by ConformationAtom::tripeptide_bb_method_tag (uint8_t). The tag is
// a categorical discriminator carrying the per-frame DFT match state
// and source orientation:
//
//   0 = no match (tripeptide_bb_has_match == false)
//   1 = gaussian_standard_orientation (OPBE/6-31G(d,p), 19 residues)
//   2 = orca_input_orientation        (PBE/6-31G(d,p), SER regen)
//
// Finalize-only dense-buffer pattern. Establishes the int-buffer +
// (N, T) 2D H5 emission shape for the bundle's other discrete-categorical
// TRs (e.g. LarsenHBondCountTimeSeries).
//
// Emission shape:
//
//   /trajectory/tripeptide_bb_method_tag_time_series/
//     method_tag     (N, T)     uint8
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name = "TripeptideBackboneMethodTagTimeSeriesTrajectoryResult"
//       units       = "categorical"
//       dtype       = "uint8"
//       legend      = "0=no_match, 1=opbe_gaussian, 2=pbe_orca_ser"
//       n_atoms, n_frames, finalized
//
// No parity/irrep_layout attributes: this is a discrete categorical,
// not a tensor irrep. Downstream consumers read the legend attribute
// to interpret the integer values.
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

class TripeptideBackboneMethodTagTimeSeriesTrajectoryResult
    : public TrajectoryResult {
public:
    std::string Name() const override {
        return "TripeptideBackboneMethodTagTimeSeriesTrajectoryResult";
    }

    // No trajectory-scope dependencies; capture-as-is project
    // precondition. uint8_t field defaults to 0 ("no match") naturally
    // — no separate sentinel needed.
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<TripeptideBackboneMethodTagTimeSeriesTrajectoryResult>
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
    // Per-atom growing buffers of uint8 method tags. Flattened into an
    // atom-major DenseBuffer<uint8_t> at Finalize.
    std::vector<std::vector<std::uint8_t>> per_atom_tag_;
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
