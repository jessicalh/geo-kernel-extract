#pragma once
//
// TripeptideBackboneResidualVecTimeSeriesTrajectoryResult: per-atom
// per-frame time series of the post-Kabsch positional residual
// produced by TripeptideBackboneShieldingResult (Vec3, units Å, in
// the protein's lab frame). Finalize-only dense-buffer pattern,
// mirrors PositionsTimeSeriesTrajectoryResult for the Vec3 layout
// and TripeptideBackboneShieldingTimeSeriesTrajectoryResult for the
// no-deps capture-as-is source-field shape.
//
// The residual is defined as `aligned_DFT_position - protein_position`
// in TripeptidePoseAssembler; sign matters. It encodes chi-grid
// coarseness and DFT geometry mismatch as a per-atom ML feature alongside
// the σ_BB^i tensor. Deep sidechain atoms can routinely carry 2-4 Å
// residuals without being rejected on the central path.
//
// No-match contract:
//   The ConformationAtom field defaults to Vec3::Zero() and is only
//   written for atoms with `tripeptide_bb_has_match == true`. This TR
//   captures the match bit with the vector and NaN-fills unmatched H5
//   xyz rows. A matched residual that is physically zero remains zero.
//
// Emission shape:
//
//   /trajectory/tripeptide_bb_residual_vec_time_series/
//     xyz            (N, T, 3)  float64
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     attrs:
//       result_name    = "TripeptideBackboneResidualVecTimeSeriesTrajectoryResult"
//       irrep_layout   = "x,y,z"
//       normalization  = "cartesian"
//       parity         = "1o"
//       units          = "angstrom"
//       n_atoms, n_frames, finalized
//
// Why parity "1o": the residual is a polar 3-vector (spatial
// displacement in the protein's lab frame), behaves under rotation
// and inversion like positions and the electric field. Cartesian
// layout: the three doubles are literally x, y, z in the protein's
// lab frame, not the e3nn real-spherical (m=-1, m=0, m=+1) ordering;
// the attribute pair (irrep_layout="x,y,z", normalization="cartesian")
// pins this for downstream consumers.
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

class TripeptideBackboneResidualVecTimeSeriesTrajectoryResult
    : public TrajectoryResult {
public:
    std::string Name() const override {
        return "TripeptideBackboneResidualVecTimeSeriesTrajectoryResult";
    }

    // No trajectory-scope dependencies; the underlying ConformationResult
    // (TripeptideBackboneShieldingResult) is conditionally attached when
    // the [databases].tensorcs15 DSN is configured. This TR captures
    // the per-atom residual and `tripeptide_bb_has_match` applicability
    // state each frame. H5 rows are NaN when the source did not attach
    // or the atom had no DFT match.
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<TripeptideBackboneResidualVecTimeSeriesTrajectoryResult>
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
    // Per-atom growing buffers of Vec3. Flattened into an atom-major
    // DenseBuffer<Vec3> at Finalize.
    std::vector<std::vector<Vec3>> per_atom_residual_;
    // Internal atom-by-frame applicability, captured from
    // ConformationAtom::tripeptide_bb_has_match. It is intentionally not
    // emitted as a new H5 dataset; it only controls NaN filling of xyz.
    std::vector<std::vector<std::uint8_t>> per_atom_has_match_;
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;

    // Per-frame source-present mask. 1 if the source ConformationResult
    // (TripeptideBackboneShieldingResult) was attached this frame, or if
    // the test-only override forces it present; 0 otherwise. Emitted as
    // the `source_attached_per_frame` H5 dataset for downstream provenance.
    // When all-zero (calc never ran), WriteH5Group skips emission
    // entirely per the "absent, not faked" discipline.
    std::vector<std::uint8_t> source_present_per_frame_;
    bool force_source_present_for_testing_ = false;

    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
