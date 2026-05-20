#pragma once
//
// JCouplingTimeSeriesTrajectoryResult: per-residue per-frame Karplus
// 3J-coupling observables. Thin downstream transform on positions +
// the residue backbone-cache + AminoAcidType.chi_angles — no source
// ConformationResult dependency (positions and the typed substrate
// are always present at tp.Seed time, same as DihedralTimeSeries).
//
// Why a J-coupling TR? 3J couplings are direct experimental NMR
// observables for backbone phi (³J(HN, Hα), ³J(HN, Cβ), ³J(HN, C'),
// ³J(Hα, C')) and sidechain chi1 rotamer populations (³J(N, Cγ),
// ³J(C', Cγ), ³J(Hα, Hβ)). These close the loop between the
// dihedral-family TR slice's persisted phi/psi/chi timelines and the
// experiment-side measurements an NMR advisor works with. TALOS-N
// basis + Vögeli / Wang-Bax / Pérez rotamer disambiguation start here.
//
// Convention: backbone channels are evaluated as
//   J = A*cos^2(phi_project + theta_project)
//       + B*cos(phi_project + theta_project) + C
// where phi_project is the same IUPAC-signed C(prev)-N-CA-C angle
// emitted by DihedralTimeSeries. The project phi sign is opposite to
// the Wang-Bax / Vogeli plotting convention, so the theta constants
// in PhysicalConstants.h are stored in the project convention
// (HN-Halpha=+pi/3, HN-Cbeta=-pi/3, HN-C'=0, Halpha-C'=+pi/3).
// Chi1 channels use the actual atom-by-atom 4-atom dihedral directly;
// Perez 2001 Table 2 folds the substituent offset into the per-channel
// coefficients.
//
// Channels (R, T) float64, Hz:
//
//   J_HN_Halpha     phi observable via the H-N-Cα-Hα dihedral.
//                   Vuister & Bax 1993 JACS 115:7772 parametrization
//                   (DOI 10.1021/ja00070a024):
//                     ³J(HN, Hα) = 6.51·cos²(θ) - 1.76·cos(θ) + 1.60
//                   PRO: NaN (no amide H). GLY: uses Residue.HA which
//                   is HA2 by Residue.h convention; HA3 is not
//                   separately measured here.
//
//   J_HN_Halpha_Vogeli
//                   phi observable via the H-N-Cα-Hα dihedral — same
//                   atomic dihedral as J_HN_Halpha, alternate Karplus
//                   parametrization. Vögeli, Ying, Grishaev & Bax 2007
//                   JACS 129:9377 (DOI 10.1021/ja070324o, Table 1
//                   "rigid" row, page 9383):
//                     ³J(HN, Hα) = 7.97·cos²(θ) - 1.26·cos(θ) + 0.63
//                   Methods-accumulate alternate to Vuister-Bax 1993
//                   (feedback_methods_accumulate: both stay; the
//                   difference is reportable).
//
//   J_HN_Cbeta      phi observable via the H-N-Cα-Cβ dihedral.
//                   Wang & Bax 1996 JACS 118:2483 NMR/X-ray refined
//                   parametrization (DOI 10.1021/ja9535524, Table 1
//                   row 3, page 2487):
//                     ³J(HN, Cβ) = 3.39·cos²(θ) - 0.94·cos(θ) + 0.07
//                   PRO: NaN (no amide H). GLY: NaN (no Cβ).
//                   Orthogonal phi probe to ³J(HN, Hα).
//
//   J_HN_Cprime     phi observable via the H-N-Cα-C dihedral.
//                   Wang & Bax 1996 JACS 118:2483 NMR/X-ray refined
//                   parametrization (DOI 10.1021/ja9535524, Table 1
//                   row 4, theta=0 deg, page 2487):
//                     ³J(HN, C') = 4.32·cos²(θ) + 0.84·cos(θ) + 0.00
//                   Note: B is POSITIVE; J can be slightly negative
//                   at the vertex (~-0.04 Hz). PRO: NaN (no amide H).
//                   (Row-mapping fixed 2026-05-20 per codex F1: prior
//                   bundle attributed row 2 values to this channel.)
//
//   J_Halpha_Cprime phi observable via the Hα-Cα-N-C'(prev) dihedral
//                   -- 3-bond path crosses the peptide bond at the
//                   previous residue's C', rotation around N-Cα
//                   (phi axis), per Vuister teaching lecture sect 6.1
//                   + Vogeli 2007 page 9384 (3J(C'(i-1), Hα) listed
//                   among the six phi-related couplings). Wang & Bax
//                   1996 JACS 118:2483 NMR/X-ray refined fit row 2
//                   (DOI 10.1021/ja9535524, theta=-60 deg, page 2487):
//                     ³J(Hα, C') = 3.75·cos²(θ) + 2.19·cos(θ) + 1.28
//                   Note: B is POSITIVE. N-terminal residue: NaN (no
//                   C'(prev)). GLY: uses HA2. (Atom path corrected
//                   2026-05-20 per codex F2: prior bundle used
//                   HA-Cα-C-N(next), which is rotation around CA-C
//                   = psi axis, NOT phi.)
//
//   J_N_Cgamma      chi1 rotamer observable via the N-Cα-Cβ-Cγ
//                   dihedral (= chi1 directly). Pérez, Löhr,
//                   Rüterjans & Schmidt 2001 JACS 123:7081 (DOI
//                   10.1021/ja003724j, Table 2 page 7086, consensus
//                   row):
//                     ³J(N, Cγ) = 1.29·cos²(θ) - 0.49·cos(θ) + 0.37
//                   GLY / ALA: NaN (no chi1). SER / CYS / THR: NaN
//                   (chi1 terminal is OG/SG/OG1, not Cγ -- gated by
//                   element of chi[0].a[3]; per-residue heteroatom
//                   terminal coefficients in Pérez Table 2 are a
//                   different chemical observable, not shipped here).
//
//   J_Cprime_Cgamma chi1 rotamer observable via the C-Cα-Cβ-Cγ
//                   dihedral. Pérez 2001 (Table 2 consensus row,
//                   byte-verified 2026-05-19):
//                     ³J(C', Cγ) = 2.31·cos²(θ) - 0.87·cos(θ) + 0.55
//                   (Prior commit carried (1.74, -0.57, 0.25); fixed
//                   at byte-verification 2026-05-19.)
//                   GLY / ALA: NaN. SER / CYS / THR: NaN (same
//                   element gate as J_N_Cgamma).
//
//   J_Halpha_Hbeta2 chi1 rotamer observable via the HA-Cα-Cβ-HB2
//                   dihedral. Pérez 2001 Table 2 consensus row:
//                     ³J(Hα, Hβ) = 7.23·cos²(θ) - 1.37·cos(θ) + 2.22
//                   Most residues carry a prochiral Hβ pair (HB2,
//                   HB3 in IUPAC ordering); we emit them as two
//                   separate channels.
//                   GLY: NaN (no Cβ). ALA: NaN (Cβ has a methyl, not
//                   methylene — emitting HB1/HB2/HB3 averaged is not
//                   the same chemical observable; consumers can
//                   compute directly if needed).
//                   Ile/Val/Thr: have a SINGLE Hβ (methine); we emit
//                   that Hβ in BOTH J_Halpha_Hbeta2 AND
//                   J_Halpha_Hbeta3 (the channels carry the same
//                   value and the mask is uniform).
//
//   J_Halpha_Hbeta3 As J_Halpha_Hbeta2 but with HB3.
//
// All Karplus coefficients live in `PhysicalConstants.h` with full
// literature citations + reference PDFs in `references/`. All channel
// families are byte-verified against the source PDFs as of 2026-05-19.
//
// Static per-residue masks (R,):
//
//   J_HN_Halpha_exists  uint8 — 1 if C(prev) + H + N + CA + C + HA all cached
//                                (i.e. residue can have ³J(HN, Hα));
//                                also gates J_HN_Halpha_Vogeli (same
//                                atomic dihedral).
//   J_HN_Cbeta_exists   uint8 — 1 if C(prev) + H + N + CA + C + CB all cached
//                                (PRO=0 via H; GLY=0 via CB).
//   J_HN_Cprime_exists  uint8 — 1 if C(prev) + H + N + CA + C all cached
//                                (PRO=0 via H).
//   J_Halpha_Cprime_exists  uint8 — 1 if C(prev) + N + CA + C + HA all
//                                    cached. The physical atom path is
//                                    HA-CA-N-C'(prev), but the emitted
//                                    value is phi-derived and therefore
//                                    also needs current-residue C.
//                                    N-terminus=0 via C(prev), looked
//                                    up via canonical
//                                    Protein::BackbonePredecessor.
//   J_chi1_exists       uint8 — 1 if chi[0].Valid() (residue has
//                                chi1 defined; GLY/ALA → 0). Necessary
//                                but NOT sufficient for J_N_Cgamma /
//                                J_Cprime_Cgamma (those further
//                                require Element::C at chi[0].a[3]).
//   J_N_Cgamma_exists   uint8 — 1 if chi[0].Valid() AND chi[0].a[3]
//                                element is Carbon. SER/CYS/THR → 0.
//   J_Cprime_Cgamma_exists  uint8 — 1 if J_N_Cgamma_exists AND C/CA
//                                    are cached. SER/CYS/THR → 0.
//   J_Halpha_Hbeta_exists  uint8 — 1 if HA + CA + CB + at-least-one-HB
//                                   present (GLY=0; ALA may carry
//                                   methyl HBs but we deliberately
//                                   NaN — see channel docstring).
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
    std::vector<std::vector<double>> j_hn_halpha_vogeli_;
    std::vector<std::vector<double>> j_hn_cbeta_;
    std::vector<std::vector<double>> j_hn_cprime_;
    std::vector<std::vector<double>> j_halpha_cprime_;
    std::vector<std::vector<double>> j_n_cgamma_;
    std::vector<std::vector<double>> j_cprime_cgamma_;
    std::vector<std::vector<double>> j_halpha_hbeta2_;
    std::vector<std::vector<double>> j_halpha_hbeta3_;

    // Static per-residue masks (computed once at Create).
    std::vector<std::uint8_t> j_hn_halpha_exists_;
    std::vector<std::uint8_t> j_hn_cbeta_exists_;
    std::vector<std::uint8_t> j_hn_cprime_exists_;
    std::vector<std::uint8_t> j_halpha_cprime_exists_;
    std::vector<std::uint8_t> j_chi1_exists_;
    // STRICTER chi1 masks: also require chi[0].a[3] = Element::C
    // (NaN for SER/CYS/THR whose chi1 terminal is OG/SG/OG1).
    std::vector<std::uint8_t> j_n_cgamma_exists_;
    std::vector<std::uint8_t> j_cprime_cgamma_exists_;
    std::vector<std::uint8_t> j_halpha_hbeta_exists_;

    // Per-residue Hbeta atom-index caches (resolved at Create time
    // by walking residue atoms and matching pdb_atom_name). NONE for
    // residues without that Hbeta.
    std::vector<std::size_t> hb2_index_;
    std::vector<std::size_t> hb3_index_;

    // Per-residue C(prev) cache (Halpha-Cprime needs the previous
    // residue's C' atom; NONE for N-terminus / chain start). Resolved
    // via Protein::BackbonePredecessor at Create time (typed bond-
    // graph query; correct on ACE/NME caps, insertion codes, and
    // cyclic peptides).
    std::vector<std::size_t> c_prev_index_;

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
