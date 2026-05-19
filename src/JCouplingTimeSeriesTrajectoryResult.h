#pragma once
//
// JCouplingTimeSeriesTrajectoryResult: per-residue per-frame Karplus
// 3J-coupling observables. Thin downstream transform on positions +
// the residue backbone-cache + AminoAcidType.chi_angles — no source
// ConformationResult dependency (positions and the typed substrate
// are always present at tp.Seed time, same as DihedralTimeSeries).
//
// Why a J-coupling TR? 3J couplings are direct experimental NMR
// observables for backbone phi (³J(HN, Hα)) and sidechain chi1
// rotamer populations (³J(N, Cγ), ³J(C', Cγ)). These close the
// loop between the dihedral-family TR slice's persisted phi/psi/chi
// timelines and the experiment-side measurements an NMR advisor
// works with. TALOS-N basis + Vögeli / Wang-Bax / Pérez rotamer
// disambiguation start here.
//
// Channels (R, T) float64, Hz:
//
//   J_HN_Halpha     phi observable via the H-N-Cα-Hα dihedral.
//                   Wang & Bax 1996 JACS 118:2483 parametrization
//                   for the Karplus equation:
//                     ³J(HN, Hα) = 6.51·cos²(θ) - 1.76·cos(θ) + 1.60
//                   θ = signed dihedral H-N-CA-HA (IUPAC). PRO: NaN
//                   (no amide H). GLY: uses Residue.HA which is HA2
//                   by Residue.h convention; HA3 is not separately
//                   measured here. The slight averaging error vs
//                   GLY's two Hα is documented; consumers needing
//                   strict pro-R/pro-S resolution should compute
//                   directly.
//
//   J_N_Cgamma      chi1 rotamer observable via the N-Cα-Cβ-Cγ
//                   dihedral (= chi1 directly). Pérez et al. 2001
//                   JACS 123:7081 parametrization:
//                     ³J(N, Cγ) = 1.29·cos²(θ) - 0.49·cos(θ) + 0.37
//                   GLY / ALA: NaN (no chi1).
//
//   J_Cprime_Cgamma chi1 rotamer observable via the C'-Cα-Cβ-Cγ
//                   dihedral (chi1 with C' as the leading atom
//                   instead of N — 120° offset on the same Cβ-Cα
//                   axis). Pérez et al. 2001:
//                     ³J(C', Cγ) = 1.74·cos²(θ) - 0.57·cos(θ) + 0.25
//                   GLY / ALA: NaN.
//
// Static per-residue masks (R,):
//
//   j_hn_halpha_exists  uint8 — 1 if H + N + CA + HA all cached
//                                (i.e. residue can have ³J(HN, Hα)).
//   j_chi1_exists       uint8 — 1 if chi[0].Valid() (residue has
//                                chi1 defined).
//
// Per-atom lookup:
//
//   residue_index_per_atom  (N,) int32 — atom_i → residue_i broadcast
//                                         for the SDK / viewer.
//
// Per-frame metadata:
//
//   frame_indices, frame_times, source_attached_per_frame (trivially
//   all-1; positions always present at tp.Seed time —
//   source_attached_policy="always_attached"). Per OBJECT_MODEL.md
//   "Conditional-attach TR discipline (2026-05-15)" canonical-statement
//   for SDK uniformity.
//
// Convention pins (group attrs):
//   karplus_form = "3J = A·cos²(θ) + B·cos(θ) + C; θ in radians via
//                    IUPAC signed dihedral atan2"
//   J_HN_Halpha_coefficients = "A=6.51, B=-1.76, C=1.60 (Wang & Bax
//                                1996 JACS 118:2483)"
//   J_N_Cgamma_coefficients  = "A=1.29, B=-0.49, C=0.37 (Pérez et al.
//                                2001 JACS 123:7081)"
//   J_Cprime_Cgamma_coefficients = "A=1.74, B=-0.57, C=0.25 (Pérez
//                                    et al. 2001)"
//   dihedral_convention = "H-N-CA-HA for J(HN,Hα); N-CA-CB-CG for
//                           J(N,Cγ); C-CA-CB-CG for J(C',Cγ).
//                           IUPAC sign convention via atan2."
//   GLY_caveat = "GLY HA uses Residue.HA (HA2 by Residue.h cache
//                  convention); HA3 not separately measured. Consumers
//                  needing pro-R/pro-S resolution should compute
//                  directly via residue's two Hα atom indices."
//
// Movie-target note: the J-coupling values are scalars in Hz and make
// natural color-by overlays for the viewer's per-residue glyph
// rendering. The atom-axis broadcast via residue_index_per_atom is
// the SDK consumer pattern.
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class JCouplingTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "JCouplingTimeSeriesTrajectoryResult";
    }

    // No declared dependencies; reads positions + Residue backbone-
    // cache + AminoAcidType.chi_angles (via Residue.chi[k]). Always-on
    // source (positions present from tp.Seed).
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<JCouplingTimeSeriesTrajectoryResult>
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

private:
    // Per-residue per-frame growing buffers (R outer, T inner).
    std::vector<std::vector<double>> j_hn_halpha_;
    std::vector<std::vector<double>> j_n_cgamma_;
    std::vector<std::vector<double>> j_cprime_cgamma_;

    // Static per-residue masks (computed once at Create).
    std::vector<std::uint8_t> j_hn_halpha_exists_;
    std::vector<std::uint8_t> j_chi1_exists_;

    // Static atom→residue lookup.
    std::vector<std::int32_t> residue_index_per_atom_;

    // Per-frame metadata.
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;

    std::size_t n_frames_  = 0;
    bool        finalized_ = false;
};

}  // namespace nmr
