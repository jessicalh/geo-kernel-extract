#pragma once
//
// MopacVsFf14SbReconciliationTrajectoryResult: per-atom per-frame
// scalar cos(MOPAC Coulomb T2, FF14SB Coulomb T2) — the SIGNED cosine
// similarity between the MOPAC-charge-derived Coulomb T2 EFG kernel
// and the FF14SB-charge-derived Coulomb T2 EFG kernel, in the T2
// (symmetric traceless) 5-vector subspace.
//
// TR9 of the 13-TR plan. New cross-source pattern; its own canonical
// (only TR of this shape).
//
// SIGNED COSINE in [-1, 1]: the T2 5-vector representation
// (isometric_real_sph) is sign-deterministic — a sign-flipped T2
// means a physically different (opposite-polarisation) tensor, NOT
// an eigenvector-convention ambiguity (the latter applies only when
// comparing principal axes). cos ≈ +1: methods agree on tensor
// orientation. cos ≈ -1: methods produce opposite-polarisation
// tensors (e.g., a charge-sign disagreement at a chemically
// distinctive group like SER OG or ARG NH2 can flip the EFG sign).
// cos ≈ 0: methods produce orthogonal tensor orientations. All
// three signals are diagnostically meaningful and the calibration
// ridge MUST see the signed value, not |cos|, to expose
// chemistry-driven disagreement. (Decision 2026-05-21 per science
// adversarial review M1.)
//
// PHYSICS RATIONALE: both source calcs produce a T2 EFG kernel from
// per-atom Coulomb-summed Hessian-of-φ contributions, with charges
// derived from different methods — MOPAC PM7+MOZYME Mulliken vs
// FF14SB parameterised partial charges. cos(T2_MOPAC, T2_FF14SB)
// measures the orientational agreement between the two methods at
// each atom each frame. Per-atom-per-frame so calibration can
// stratify by atom type and observe whether agreement degrades
// dynamically.
//
// CROSS-SOURCE GATE: HasResult<MopacCoulombResult>() AND
// HasResult<CoulombResult>() must both be true for that frame to
// contribute. Either-absent → emit NaN for all atoms that frame,
// source_attached_per_frame=0. When NO frame had both, WriteH5Group
// skips the entire /trajectory/mopac_vs_ff14sb_reconciliation/
// group per canonical "absent, not faked".
//
// MAGNITUDE FLOOR: cosine is undefined when either |T2| <
// `coulomb_efg_t2_magnitude_floor` (CalculatorConfig, V/Å² —
// calibrated to the EFG signal scale, NOT the project-wide
// direction-vector floor 1e-10 which is seven orders below this
// EFG-scale floor and would let FP-noise-dominated atoms through.
// Decision 2026-05-21 per math adversarial review H1.) Per-atom NaN
// under that condition; SDK readers MUST use isfinite() to gate.
//
// Emission:
//   /trajectory/mopac_vs_ff14sb_reconciliation/
//     cos_t2         (N, T) float64 — cos(T_MOPAC, T_FF14SB) ∈ [-1, 1]
//     frame_indices  (T,)   uint64
//     frame_times    (T,)   float64 — ps
//     source_attached_per_frame (T,) uint8 — 1 iff BOTH sources attached
//     attrs:
//       result_name             = "MopacVsFf14SbReconciliationTrajectoryResult"
//       n_atoms, n_frames, source_attached_count, finalized
//       parity                  = "0e"  (rotation-invariant scalar)
//       units                   = "dimensionless"  ([-1, 1] cosine)
//       sources                 describes the MopacCoulombResult and CoulombResult T2 EFG fields
//       source_attached_policy  = "conditional -- requires BOTH ..."
//       magnitude_floor         = (value of coulomb_efg_t2_magnitude_floor)
//       magnitude_floor_units   = "V/Å^2"
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
    // Per-atom growing buffer of signed cos. NaN cells = either-absent
    // source or either-side |T2| below magnitude floor.
    std::vector<std::vector<double>> per_atom_cos_;
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    double magnitude_floor_ = 0.0;  // cached at Create from CalculatorConfig
    std::size_t n_frames_              = 0;
    std::size_t source_attached_count_ = 0;
    bool        finalized_             = false;
};

}  // namespace nmr
