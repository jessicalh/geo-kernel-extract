#pragma once
//
// TripeptidePoseAssembler: the validation layer between
// TripeptideDftRecord (the DB row) and the per-atom shielding output.
//
// Why a separate "layer" (in the software-engineering sense): any
// single mismapping (wrong DFT atom → wrong protein atom; wrong
// Kabsch input) silently puts a tensor on the wrong atom. The
// kernels feed downstream ML; one bad atom mapping invalidates the
// protein. This layer is the place atom identity is checked end-to-end.
//
// Discipline: typed-identity matching across both substrates.
//
//   DFT side — TripeptideDftRecord carries a perceived LarsenTripeptide
//   (built in TripeptideDftTable::QueryNearest by
//   PerceiveLarsenTripeptide). Each atom holds an AtomMechanicalIdentity
//   derived from a bond-graph isomorphism against canonical chemistry
//   from AminoAcidType — NOT from positional ordering in the DFT
//   record. The DFT atom_idx is preserved for log correlation but is
//   never used as the matching key.
//
//   Protein side — LegacyAmberTopology::SemanticAt(ai) returns an
//   AtomSemanticTable with the same AtomMechanicalIdentity tuple
//   (Element + Locant + BranchAddress + DiastereotopicIndex +
//   BackboneRole) that perception emits on the DFT side. Both sides
//   speak the same typed vocabulary.
//
//   Matching — identity equality (strict, then DiastereotopicIndex/
//   BranchAddress-relaxed for prochiral methylene + aromatic-symmetric
//   equivalence classes), with nearest-spatial tiebreak within an
//   equivalence class. NO element-walks, NO element+distance heuristics
//   in the central-residue assembly path. The cap path reads BB+CB
//   atoms from typed LarsenResidue slots directly. Perception-or-nothing:
//   absent rec.larsen → decline the residue with structured
//   OperationLog::Warn.
//
// Per-atom residual is emitted as a Vec3 (the displacement
// aligned_position − protein_atom_position). It is NOT a rejection
// gate — large residuals stay in the output for the ML model to
// consume as a feature alongside the rotated shielding tensor (per
// feedback_residual_as_ml_feature). The n_above_threshold counter is
// a diagnostic stat, not a filter.
//
// History: the previous element-pattern + 5 Å radius + threshold-
// rejection design was retired 2026-05-11 alongside the introduction
// of LarsenResidue perception. See
// spec/plan/larsen-residue-design-2026-05-11.md for the durable design
// and the retired-design pointer.
//

#include "Types.h"
#include "TripeptideDftTable.h"

#include <cstdint>
#include <string>
#include <vector>

namespace nmr {

class Protein;
class ProteinConformation;


// Which side of the tripeptide the data is read from.
//
//   Central: σ_BB^i — read at the central residue's atoms. The central
//            residue's chemistry matches the protein residue's; full
//            sidechain mapping applies.
//   NTerm:   Δσ_BB^{i+1}(i) — read at the N-terminal ALA cap atoms of
//            the (i+1)-centered tripeptide. The cap represents the
//            residue BEFORE central in tripeptide, i.e. protein
//            residue i. Cap is always ALA; only BB+CB+HBs.
//   CTerm:   Δσ_BB^{i-1}(i) — read at the C-terminal ALA cap atoms of
//            the (i-1)-centered tripeptide. Cap represents the
//            residue AFTER central, i.e. protein residue i. Cap is
//            always ALA.
enum class TripeptidePoseSide : std::uint8_t {
    Central = 0,
    NTerm   = 1,
    CTerm   = 2,
};


// One DFT atom mapped to a protein atom after pose alignment.
struct AlignedDftAtom {
    int         dft_atom_idx     = -1;          // into TripeptideDftRecord::atoms
    std::size_t protein_atom_idx = 0;            // into Protein::atoms
    Element     element          = Element::Unknown;

    Vec3 aligned_position = Vec3::Zero();         // R*(p_dft - cs) + ct
    Vec3 residual_vec     = Vec3::Zero();         // aligned - protein_pos
    double residual_distance = 0.0;                // residual_vec.norm()

    Mat3            shielding_tensor_aligned    = Mat3::Zero();    // R σ R^T, ppm
    SphericalTensor shielding_spherical_aligned;                    // Decompose

    // Substrate cross-check outcome. True iff the protein atom's
    // typed substrate role (BackboneRole / Locant) matches the role
    // implied by the canonical DFT cap position. When false, the
    // atom is excluded from the output and the caller is warned.
    bool substrate_role_agrees = false;
};


// Output of one Assemble call.
struct AssembledTripeptide {
    bool ok = false;                          // false on hard structural failure
    int  calc_id = 0;
    std::string frame_type;
    TripeptidePoseSide side = TripeptidePoseSide::Central;
    int  n_chi_axes_used = 0;                 // caller-set diagnostic

    Mat3   rotation        = Mat3::Identity();
    Vec3   source_centroid = Vec3::Zero();
    Vec3   target_centroid = Vec3::Zero();
    double backbone_kabsch_rmsd = 0.0;         // 3-point Kabsch fit, Å

    std::vector<AlignedDftAtom> aligned_atoms;  // atoms that passed both
                                                // path 1 + path 2 checks

    // Aggregate diagnostics (across atoms that passed validation).
    int    n_substrate_disagreements = 0;     // path 2 failures (rejected)
    int    n_above_threshold = 0;             // residual > threshold (rejected)
    double max_residual_A = 0.0;
    double mean_residual_A = 0.0;
};


// Assemble a tripeptide DFT record's data onto a protein residue's
// pose. See top-of-file for the two-path validation discipline.
//
// validation_threshold_A: the threshold's role differs between the
// cap path and the central path.
//
//   * Cap path (`AssembleAlaCap` → `EmitAlignedAtom`): atoms whose
//     residual_distance exceeds this ARE excluded from
//     aligned_atoms. The cap path's seven slots are tightly
//     constrained by typed identity so a large residual on a
//     BB/Cβ atom indicates a Kabsch failure or a substrate
//     disagreement worth excising.
//   * Central path (`AssembleCentralTyped`): atoms whose
//     residual_distance exceeds this are kept but counted in
//     `out.n_above_threshold` as a diagnostic stat. The residual
//     itself is the load-bearing ML feature
//     (per `feedback_residual_as_ml_feature`); rejection on
//     residual magnitude would discard signal the calibration
//     expects to consume. Chi-grid coarseness puts deep-sidechain
//     atoms (Arg/Lys terminal N region) at 2-4 Å residual
//     routinely, which is normal not a Kabsch failure.
//
// Default 3.0 Å — chosen for the cap-path gate; on the central path
// it is the threshold above which n_above_threshold ticks.
//
// substrate_check_strict: when true, atoms whose typed substrate role
// doesn't match the canonical cap-slot role are excluded from
// aligned_atoms (and the disagreement is logged). When false they
// are kept with substrate_role_agrees=false. Default true: the
// substrate cross-check is the load-bearing independent validation.
// Applies to the cap path's per-slot dispatch in `EmitAlignedAtom`;
// the central path uses typed-identity-equality at the candidate
// step and substrate_role_agrees is set to true by construction
// (a non-empty candidate set means the chemistry agreed).
//
// Returns AssembledTripeptide with ok=false on hard structural error:
//   - Protein residue's backbone N/CA/C cache incomplete
//   - DFT record's cap atoms can't be identified by element pattern
//   - Kabsch input is degenerate
AssembledTripeptide AssembleTripeptide(
    const Protein&             protein,
    const ProteinConformation& conf,
    std::size_t                protein_residue_idx,
    const TripeptideDftRecord& rec,
    TripeptidePoseSide         side,
    double                     validation_threshold_A = 3.0,
    bool                       substrate_check_strict = true);


}  // namespace nmr
