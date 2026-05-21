#pragma once
//
// MopacVsFf14SbReconciliationTrajectoryResult: per-atom per-frame
// scalar |cos(MOPAC Coulomb T2, FF14SB Coulomb T2)| — the absolute
// cosine similarity between the MOPAC-charge-derived Coulomb shielding
// tensor and the FF14SB-charge-derived Coulomb shielding tensor, in
// the T2 (symmetric traceless) subspace.
//
// TR9 of the 13-TR plan. New cross-source pattern; its own canonical
// (only TR of this shape).
//
// PHYSICS RATIONALE: both source calcs produce a T2 shielding
// contribution from a Coulomb EFG, but with charges derived from
// different methods — MOPAC PM7+MOZYME Mulliken vs FF14SB
// parameterised partial charges. |cos(T2_MOPAC, T2_FF14SB)| ∈ [0, 1]
// measures the orientational agreement between the two methods at
// each atom each frame. Aligned (≈1): both methods agree on the
// tensor's principal axes. Perpendicular (≈0): methods disagree.
// Per-atom-per-frame so calibration can stratify by atom type and
// observe whether agreement degrades dynamically.
//
// CROSS-SOURCE GATE: HasResult<MopacCoulombResult>() AND
// HasResult<CoulombResult>() must both be true for that frame to
// contribute. Either-absent → emit NaN for all atoms that frame,
// source_attached_per_frame=0. When NO frame had both, WriteH5Group
// skips the entire /trajectory/mopac_vs_ff14sb_reconciliation/
// group per canonical "absent, not faked".
//
// ZERO-MAGNITUDE HANDLING: cosine is undefined when either
// |T2_MOPAC| or |T2_FF14SB| is near zero (atom has no Coulomb
// shielding contribution from one or both methods — typical for
// remote-from-charged-groups atoms). Per-atom NaN under that
// condition; threshold from CalculatorConfig::Get("near_zero_vector_norm_threshold").
// SDK readers MUST use isfinite() to distinguish "real measurement"
// from "below noise floor / both-absent / either-absent".
//
// Emission:
//   /trajectory/mopac_vs_ff14sb_reconciliation/
//     abs_cos_t2     (N, T) float64 — |cos(T_MOPAC, T_FF14SB)| ∈ [0, 1]
//     frame_indices  (T,)   uint64
//     frame_times    (T,)   float64 — ps
//     source_attached_per_frame (T,) uint8 — 1 iff BOTH sources attached
//     attrs:
//       result_name             = "MopacVsFf14SbReconciliationTrajectoryResult"
//       n_atoms, n_frames, source_attached_count, finalized
//       parity                  = "0e"  (rotation-invariant scalar)
//       units                   = "dimensionless"  ([0, 1] cosine)
//       sources                 = "MopacCoulombResult.mopac_coulomb_shielding_contribution
//                                   + CoulombResult.coulomb_shielding_contribution
//                                   (both T2 SphericalTensor; |cos|
//                                   in the T2 5-vector subspace)"
//       source_attached_policy  = "conditional -- requires BOTH ..."
//       zero_magnitude_threshold = (value of near_zero_vector_norm_threshold)
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class MopacVsFf14SbReconciliationTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "MopacVsFf14SbReconciliationTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<MopacVsFf14SbReconciliationTrajectoryResult>
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
    // Per-atom growing buffer of |cos|. NaN cells = either-absent
    // source or either-side zero-magnitude T2.
    std::vector<std::vector<double>> per_atom_abs_cos_;
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    double zero_magnitude_threshold_ = 0.0;  // cached at Create from CalculatorConfig
    std::size_t n_frames_              = 0;
    std::size_t source_attached_count_ = 0;
    bool        finalized_             = false;
};

}  // namespace nmr
