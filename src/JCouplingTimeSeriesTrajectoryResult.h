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
//                   Vuister & Bax 1993 JACS 115:7772 parametrization
//                   (DOI 10.1021/ja00070a024):
//                     ³J(HN, Hα) = 6.51·cos²(θ) - 1.76·cos(θ) + 1.60
//                   θ = signed dihedral H-N-CA-HA (IUPAC); equals
//                   (φ - 60°) under ideal tetrahedral L-amino-acid
//                   backbone geometry. PRO: NaN (no amide H). GLY:
//                   uses Residue.HA which is HA2 by Residue.h
//                   convention; HA3 is not separately measured here.
//                   The slight averaging error vs GLY's two Hα is
//                   absorbed in the Vuister-Bax parametrization (the
//                   fit included glycine 3J values); consumers needing
//                   strict pro-R / pro-S resolution should compute
//                   directly from the two Hα indices.
//
//   J_HN_Cbeta      phi observable via the H-N-Cα-Cβ dihedral.
//                   Wang & Bax 1996 JACS 118:2483 NMR/X-ray refined
//                   parametrization (DOI 10.1021/ja9535524, Table 1
//                   row 3, page 2487):
//                     ³J(HN, Cβ) = 3.39·cos²(θ) - 0.94·cos(θ) + 0.07
//                   θ = signed dihedral H-N-CA-CB (IUPAC); equals
//                   (φ + 60°) under ideal tetrahedral L-amino-acid
//                   geometry. PRO: NaN (no amide H). GLY: NaN (no
//                   Cβ). Orthogonal phi probe to ³J(HN, Hα); the
//                   pair together overdetermines phi.
//
//   J_N_Cgamma      chi1 rotamer observable via the N-Cα-Cβ-Cγ
//                   dihedral (= chi1 directly). Pérez, Löhr,
//                   Rüterjans & Schmidt 2001 JACS 123:7081
//                   parametrization (DOI 10.1021/ja003724j):
//                     ³J(N, Cγ) = 1.29·cos²(θ) - 0.49·cos(θ) + 0.37
//                   GLY / ALA: NaN (no chi1).
//
//   J_Cprime_Cgamma chi1 rotamer observable via the C'-Cα-Cβ-Cγ
//                   dihedral (chi1 with C' as the leading atom
//                   instead of N — 120° offset on the same Cβ-Cα
//                   axis). Pérez et al. 2001 (same paper):
//                     ³J(C', Cγ) = 1.74·cos²(θ) - 0.57·cos(θ) + 0.25
//                   GLY / ALA: NaN.
//
// All Karplus coefficients live in `PhysicalConstants.h` with full
// literature citations + reference PDFs in `references/`. Caveat:
// the Pérez 2001 (N,Cγ) and (C',Cγ) coefficients are widely-circulated
// but the 2026-05-19 reference audit could not byte-verify them
// against the paywalled original paper -- circulated unchanged in
// downstream NMR software (TALOS-N, NMRViewJ) but flagged for
// institutional-access verification.
//
// Static per-residue masks (R,):
//
//   J_HN_Halpha_exists  uint8 — 1 if H + N + CA + HA all cached
//                                (i.e. residue can have ³J(HN, Hα)).
//   J_HN_Cbeta_exists   uint8 — 1 if H + N + CA + CB all cached
//                                (PRO=0 via H; GLY=0 via CB).
//   J_chi1_exists       uint8 — 1 if chi[0].Valid() (residue has
//                                chi1 defined; GLY/ALA → 0).
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
// Convention pins (group attrs) — see emitter for verbatim strings;
// coefficient pins always include the full citation + DOI; the
// emitter reads the numeric values directly from `PhysicalConstants.h`
// so the H5 attrs and the compiled binary stay in lockstep.
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
    std::vector<std::vector<double>> j_hn_cbeta_;
    std::vector<std::vector<double>> j_n_cgamma_;
    std::vector<std::vector<double>> j_cprime_cgamma_;

    // Static per-residue masks (computed once at Create).
    std::vector<std::uint8_t> j_hn_halpha_exists_;
    std::vector<std::uint8_t> j_hn_cbeta_exists_;
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
