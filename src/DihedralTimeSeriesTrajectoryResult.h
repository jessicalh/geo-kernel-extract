#pragma once
//
// DihedralTimeSeriesTrajectoryResult: per-residue, per-frame timelines of
// every backbone and sidechain dihedral angle, plus a Ramachandran-region
// classification and the static residue-property masks consumers need to
// interpret per-row semantics.
//
// Self-contained: reads positions + Residue.chi[k] cached atom indices +
// Protein chain structure. No declared `Dependencies()` — positions and
// topology are always available once `tp.Seed` has run. The TR runs in
// every PerFrameExtractionSet frame and emits a complete record.
//
// Per-frame channels (per residue, R rows × T columns):
//
//   phi             N-terminal residue NaN; chain-break-in NaN
//   psi             C-terminal residue NaN; chain-break-out NaN
//   omega           C-terminal residue NaN; chain-break-out NaN
//   omega_deviation omega − π wrapped to (−π, π]; NaN at X→Pro (mask)
//   chi             (R, T, 4) — NaN where chi[k] not defined
//   rama_region     uint8 R × T (0=unassigned, 1=αR, 2=β, 3=αL,
//                                4=PPII, 5=other)
//
// Per-residue static masks (R rows):
//
//   chi_exists                 (R, 4) uint8 — chi[k] valid for this AA
//   omega_is_xpro              (R,)   uint8 — bond i→(i+1) is X→Pro
//   is_glycine                 (R,)   uint8 — Rama special: full allowed map
//   is_proline                 (R,)   uint8 — Rama special: phi constrained
//   is_pre_proline             (R,)   uint8 — psi constrained because i+1=Pro
//   residue_terminal_state     (R,)   uint8 — see ResidueTerminalState legend
//   chain_id_per_residue       (R,)   string — variable-length chain IDs
//
// Per-atom lookup (N rows):
//
//   residue_index_per_atom     (N,) int32 — atom_i → residue_i for the
//                              SDK / viewer's atom-axis broadcast pattern.
//                              Same index space as the topology sidecar's
//                              residues.npy row index.
//
// Per-frame:
//
//   frame_indices              (T,) uint64
//   frame_times                (T,) float64 (ps)
//   source_attached_per_frame  (T,) uint8 — emitted as all-1 for SDK
//                              uniformity; positions are always present.
//
// Convention pins (recorded as group attrs):
//   - angles in radians
//   - dihedral sign: IUPAC, atan2(y,x) with the standard 4-atom signed
//     formulation in PlanarGeometryResult.cpp:35-46 / ChiRotamerSelection
//     TrajectoryResult.cpp:22-40 — single implementation reused here.
//   - phi   = C(i−1)-N(i)-Cα(i)-C(i)
//     psi   = N(i)-Cα(i)-C(i)-N(i+1)
//     omega = Cα(i)-C(i)-N(i+1)-Cα(i+1)
//     chi_k from AminoAcidType.chi_angles (Residue.chi[k] pre-cached atom
//                                          indices, IUPAC sidechain order)
//   - omega_deviation is the wrapped deviation from π; emitted for every
//     well-defined peptide bond INCLUDING X→Pro (cis-trans isomerism is
//     real signal, not a deviation). The `omega_is_xpro` static mask flags
//     X→Pro rows so consumers interpret deviation appropriately.
//   - Ramachandran region binning is a simple 5-region literature grid
//     (boundaries written into the `rama_region_boundaries` attr); proper
//     re-binning against Lovell-Richardson maps is left to downstream.
//
// Chi symmetry caveats (consumer responsibility; group attr documents):
//   - PHE/TYR χ₂ has 180° ring-flip symmetry (CD1↔CD2 swap). Raw χ₂ here
//     is the IUPAC-signed value; rotamer counters that need mod-π must
//     apply it themselves.
//   - ASP χ₂ / GLU χ₃ have carboxylate-end symmetry; same caveat.
//
// Movie/playback note: this TR is the dihedral-state source for the
// viewer's color-by / glyph rendering. Per-frame resolution; no
// downsampling. Atom-axis broadcasts via `residue_index_per_atom` at
// render time (option C from the design discussion).
//

#include "TrajectoryResult.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class DihedralTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "DihedralTimeSeriesTrajectoryResult";
    }

    // No declared dependencies — reads positions and Residue.chi[k] which
    // are always available after tp.Seed (PATTERNS §15).
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<DihedralTimeSeriesTrajectoryResult>
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
    // Per-residue growing buffers, R outer × T inner. Finalize flattens
    // these into the (R, T) / (R, T, 4) H5 datasets atom-major (residue-
    // major here): residue-i-frame-0, residue-i-frame-1, ..., residue-i-
    // frame-T-1, residue-(i+1)-frame-0, ...
    std::vector<std::vector<double>>  phi_;
    std::vector<std::vector<double>>  psi_;
    std::vector<std::vector<double>>  omega_;
    std::vector<std::vector<double>>  omega_deviation_;
    std::vector<std::vector<std::array<double, 4>>> chi_;
    std::vector<std::vector<std::uint8_t>> rama_region_;

    // Static per-residue masks (computed once at Create, stable across
    // frames because Protein is invariant after tp.Seed).
    std::vector<std::uint8_t>  chi_exists_;            // R * 4 row-major
    std::vector<std::uint8_t>  omega_is_xpro_;
    std::vector<std::uint8_t>  is_glycine_;
    std::vector<std::uint8_t>  is_proline_;
    std::vector<std::uint8_t>  is_pre_proline_;
    std::vector<std::uint8_t>  residue_terminal_state_;
    std::vector<std::string>   chain_id_per_residue_;

    // Static atom→residue lookup for the viewer/SDK atom-axis broadcast.
    std::vector<std::int32_t>  residue_index_per_atom_;

    // Per-frame metadata.
    std::vector<std::size_t>   frame_indices_;
    std::vector<double>        frame_times_;
    std::vector<std::uint8_t>  source_attached_per_frame_;

    std::size_t n_frames_  = 0;
    bool        finalized_ = false;
};

}  // namespace nmr
