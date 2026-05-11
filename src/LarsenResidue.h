#pragma once
//
// LarsenResidue: typed model of one residue or cap inside a Larsen
// ProCS15 tripeptide DFT calculation.
//
// Identity is derived from bond-graph perception against canonical
// AminoAcidType chemistry at TripeptideDftRecord load time. The DFT
// row position (`dft_atom_idx`) is preserved per atom for log
// correlation but is NEVER consulted by runtime calculator code:
// chemistry questions resolve through the typed AtomMechanicalIdentity.
//
// The five pieces inside one tripeptide DFT calculation are:
//
//   ACE (acetyl cap)   -- 6 atoms, no peptide-amide H, no N
//   NCapAla            -- 10 atoms, canonical ALA (N-terminal flanking)
//   Central            -- variable (any of 20 amino acids; HIS picks
//                         HID/HIE/HIP variant by atom count)
//   CCapAla            -- 10 atoms, canonical ALA (C-terminal flanking)
//   NME (methylamide)  -- 6 atoms, no carbonyl-O, has amide H
//
// Perception algorithm (in detail in
// spec/plan/larsen-residue-design-2026-05-11.md). The Python script
// at scripts/perceive_larsen_tripeptide.py is the original algorithm
// prototype that validated the bond-graph + WL approach against
// every (tripeptide, frame_type) DB combination and the AAA Gaussian
// log — it is NOT a normative spec for the C++ object model. Names
// die at canonical construction; runtime identity is the typed
// AtomMechanicalIdentity stamped by the generated topology table.
// The Python script still returns canonical atom names from its
// `match_piece` because it predates the typed-substrate work; treat
// it as a diagnostic tool, not as the contract this code must match.
//
//   1. Build bond graph from (position, element) with element-pair-
//      specific distance cutoffs.
//   2. Identify peptide-amide bonds: sp2 C with C=O double-bond and
//      an N neighbour that itself has another heavy neighbour.
//   3. Cut amides → 5 connected components.
//   4. Order components along the directed amide chain (ACE → NME).
//   5. Per piece, run graph isomorphism against canonical chemistry
//      (AminoAcidType for the 19 non-HIS residues; HIS_VARIANTS for
//      HID/HIE/HIP; ACE/NME identities hand-coded here).
//   6. Emit typed AtomMechanicalIdentity per atom.
//
// On failure (any of: amide count != 4, segmentation != 5 pieces,
// atom-count mismatch with canonical, graph isomorphism fail), the
// perception returns nullopt and the caller (typically
// TripeptideDftTable::QueryNearest) leaves the record's `larsen`
// field empty. Downstream calculators decline records without a
// perceived LarsenTripeptide.
//
// See spec/plan/larsen-residue-design-2026-05-11.md for the full
// design rationale (supersedes the retired typed-tripeptide-topology
// design from 2026-05-10).
//

#include "SemanticEnums.h"
#include "Types.h"

#include <array>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace nmr {

struct TripeptideDftRecord;  // forward; defined in TripeptideDftTable.h


// One residue's worth of perceived atoms in a Larsen tripeptide.
class LarsenResidue {
public:
    enum class Kind : std::uint8_t {
        AceCap   = 0,  ///< Acetyl cap atoms: -C(=O)-CH3.
        NCapAla  = 1,  ///< N-terminal flanking ALA residue.
        Central  = 2,  ///< Central residue (residue field carries which).
        CCapAla  = 3,  ///< C-terminal flanking ALA residue.
        NmeCap   = 4,  ///< N-methylamide cap: -NH-CH3.
    };

    // One perceived atom. The `identity` is the typed chemistry
    // assignment from graph isomorphism; `dft_atom_idx` is the
    // 1-based row index in the Gaussian / ORCA Standard orientation
    // output (or equivalent), preserved for log correlation only.
    struct PerAtom {
        AtomMechanicalIdentity identity = {};
        int                    dft_atom_idx = 0;
        Element                element  = Element::Unknown;
        Vec3                   position = Vec3::Zero();

        // Pre-decomposed shielding tensor (irreps split at ingest
        // and carried on the DB row; preserved here for direct use
        // by calculators).
        Mat3                  shielding_tensor = Mat3::Zero();
        double                isotropic        = 0.0;
        double                anisotropy       = 0.0;
        std::array<double, 5> t2_components    = {};

        // True iff the canonical name was chosen within a K=3
        // Weisfeiler-Lehman signature class of size ≥ 2 — the bond
        // graph cannot distinguish this atom from a graph-automorphic
        // sibling (PHE/TYR CD1↔CD2, ARG NH1↔NH2, ASN/GLN HD21↔HD22,
        // methyl-Hs whose DiastereotopicIndex collapses). The
        // chemical identity tuple (Element, Locant, BackboneRole) is
        // sound; only the within-class label (BranchAddress and
        // DiastereotopicIndex) is interchangeable. Matchers should
        // drop those two fields and resolve by nearest-spatial when
        // this flag is true. For singleton WL classes (the common
        // case after K=3 — chemistry-distinct branches like ILE
        // CG1/CG2 land here), strict identity match is correct.
        bool canonical_assignment_ambiguous = false;
    };

    // One bond in the perceived covalent graph. Indices are into
    // `atoms` (this LarsenResidue's local atoms list). External
    // amide bonds are NOT recorded here — they belong to the
    // LarsenTripeptide-level chain.
    struct Bond {
        int          a     = -1;
        int          b     = -1;
        std::uint8_t order = 1;  ///< 1=single, 2=double, etc. (heuristic)
    };

    Kind                 kind     = Kind::Central;
    AminoAcid            residue  = AminoAcid::Unknown;
    std::vector<PerAtom> atoms;
    std::vector<Bond>    bonds;

    // Role-pinned slot cache. -1 if absent for this piece's kind.
    int N_idx  = -1;
    int H_idx  = -1;
    int CA_idx = -1;
    int HA_idx = -1;
    int CB_idx = -1;
    int C_idx  = -1;
    int O_idx  = -1;

    // Find the local index of an atom with the given typed identity.
    // Returns -1 if not found. For chemically-equivalent atom sets
    // (e.g., methyl Hs that collapse to one identity), returns the
    // first match.
    int LookupByIdentity(const AtomMechanicalIdentity& id) const;

    // True iff this piece has all required role-pinned slots filled
    // for its Kind. ACE: C, O, CH3-C only (no BB slots in the same
    // sense). NME: N + amide-H only.
    bool HasAllRequiredSlots() const;
};


// The full 5-piece Larsen tripeptide.
struct LarsenTripeptide {
    LarsenResidue ace;
    LarsenResidue n_cap;
    LarsenResidue central;
    LarsenResidue c_cap;
    LarsenResidue nme;

    // For HIS central residues, the variant index that perception
    // matched: 0=HID, 1=HIE, 2=HIP. -1 for non-HIS centrals or when
    // no variant has been determined. Downstream calculators can
    // cross-check this against the protein's `Residue::protonation_variant_index`
    // and warn on mismatch.
    int central_variant_index = -1;

    int TotalAtoms() const;

    // Find a (piece*, local_idx) by the global dft_atom_idx (1-based
    // in the source Gaussian / ORCA output). Returns (nullptr, -1)
    // if not found.
    std::pair<const LarsenResidue*, int> FindByDftIdx(int dft_atom_idx) const;
};


// Perceive a full LarsenTripeptide from a TripeptideDftRecord.
//
// expected_central is the amino-acid expected at the central
// position (typically `AminoAcidFromOneLetterCode(rec.tripeptide[1])`
// from the calling site). The perception verifies that the central
// piece's atom inventory matches the expected residue's canonical
// chemistry; mismatches fail the perception.
//
// his_variant_hint applies only when expected_central == HIS. Values
// are the canonical protonation variant indices: 0=HID, 1=HIE, 2=HIP.
// -1 means "no hint" (default). When a hint is provided, perception
// REQUIRES that variant to match — non-matching DB rows fail loud.
// Without a hint, perception tries HID/HIE/HIP in order and accepts
// the first that fits, which can silently misassign Hε/Hδ atoms if
// the protein's actual variant differs from the perceived one. Pass
// the protein's `Residue::protonation_variant_index` to keep the
// match honest.
//
// On failure, returns nullopt. The reason is logged via OperationLog.
std::optional<LarsenTripeptide> PerceiveLarsenTripeptide(
    const TripeptideDftRecord& rec,
    AminoAcid                  expected_central,
    int                        his_variant_hint = -1);


}  // namespace nmr
