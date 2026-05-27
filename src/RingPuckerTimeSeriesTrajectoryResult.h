#pragma once
//
// RingPuckerTimeSeriesTrajectoryResult: per-ring per-frame timelines
// of Cremer-Pople saturated-ring pucker (Q, θ) and aromatic-ring χ₂.
//
// Sources are per-frame `PlanarGeometryResult` arrays:
//   - PuckerQ()      length = SaturatedRingCount() (Q in Å)
//   - PuckerTheta()  length = SaturatedRingCount() (θ in degrees [0,360))
//   - AromaticChi2() length = AromaticRingCount()  (radians)
//
// CROSS-RESULT READ (reader side): PlanarGeometryResult exposes the
// per-frame pucker arrays + aromatic-ring χ₂ array. This TR copies
// them into per-ring growing buffers. OperationRunner calls PG each frame;
// PG attaches only if PlanarGeometryResult::Compute succeeds, which requires
// the LegacyAmber AtomSemanticTable substrate. When PG never attached, the
// H5 group is skipped per the 2026-05-15 conditional-attach discipline.
//
// Two axes (saturated and aromatic rings) -- these are DIFFERENT
// counts indexed differently inside LegacyAmberTopology. Per-ring
// metadata (parent_residue_index for each ring) emitted alongside
// for the SDK/viewer broadcast.
//
// Emission: /trajectory/ring_pucker_time_series/
//   Per-saturated-ring per-frame (S, T):
//     pucker_Q       float64  amplitude (Å)
//     pucker_theta   float64  phase (degrees, [0, 360))
//
//   Per-aromatic-ring per-frame (A, T):
//     aromatic_chi2  float64  radians (IUPAC sign convention)
//
//   Per-ring static lookups:
//     saturated_parent_residue_index   (S,) int32
//     aromatic_parent_residue_index    (A,) int32
//
//   Per-frame metadata:
//     frame_indices, frame_times, source_attached_per_frame
//
// Convention pins (group attrs):
//   pucker_convention      = "Cremer-Pople 1975 J. Am. Chem. Soc. 97, 1354"
//   pucker_Q_units         = "Angstrom"
//   pucker_theta_units     = "degrees"
//   pucker_theta_range     = "[0, 360)"
//   pucker_5ring_endvtwist = "theta mod 72 gives envelope (E) vs twist (T)
//                              endo/exo classification"
//   aromatic_chi2_units    = "radians"
//   aromatic_chi2_convention = "IUPAC signed dihedral (Ca-Cb-Cg-Cd1
//                                for PHE/TYR; Ca-Cb-Cg-Nd1 for HIS;
//                                Ca-Cb-Cg-Cd1 for TRP); matches
//                                DihedralTimeSeries chi[1] for the
//                                parent residue."
//   source                 = "PlanarGeometryResult"
//   source_attached_policy = "conditional -- PG attaches when
//                              PlanarGeometryResult::Compute succeeds."
//
// NaN-fill: PG emits NaN for rings where the geometry was degenerate
// (5-ring CP requires a well-defined ring plane). We pass NaN through
// unchanged.
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class RingPuckerTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "RingPuckerTimeSeriesTrajectoryResult";
    }

    // No declared dependency: PlanarGeometryResult attaches when the
    // LegacyAmber substrate is populated. We capture per-frame whether
    // it was attached via the source_attached gate.
    std::vector<std::type_index> Dependencies() const override {
        return {};
    }

    static std::unique_ptr<RingPuckerTimeSeriesTrajectoryResult>
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

    // Test-only: bypass the per-frame `conf.HasResult<PlanarGeometryResult>()`
    // check. SAFETY (codex review 2026-05-19): When this flag is true,
    // Compute calls `conf.Result<PlanarGeometryResult>()` which throws
    // if the Result is not actually attached. Test caller MUST attach
    // PlanarGeometryResult to every ProteinConformation before setting
    // this flag.
    void ForceSourcePresentForTesting() {
        force_source_present_for_testing_ = true;
    }

private:
    // Per-ring per-frame growing buffers.
    std::vector<std::vector<double>> pucker_Q_;        // (S outer, T inner)
    std::vector<std::vector<double>> pucker_theta_;    // (S outer, T inner)
    std::vector<std::vector<double>> aromatic_chi2_;   // (A outer, T inner)

    // Per-ring static metadata (parent residue index for each ring).
    std::vector<std::int32_t> saturated_parent_residue_index_;
    std::vector<std::int32_t> aromatic_parent_residue_index_;

    // Per-frame metadata.
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    bool force_source_present_for_testing_ = false;

    std::size_t n_frames_  = 0;
    bool        finalized_ = false;
};

}  // namespace nmr
