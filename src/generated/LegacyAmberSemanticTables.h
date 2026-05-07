#pragma once
//
// LegacyAmberSemanticTables.h -- runtime API for the typed
// chemistry-substrate tables emitted by tools/topology/build_semantic_tables.
//
// String barrier: this header contains typed-enum declarations, two
// inline composition helpers, and an inline atom-name parser
// (`ParseAtomName` + `ComputeAtomMechanicalIdentity`). No chemistry-
// string libraries (cifpp, RDKit, gemmi) are referenced. Strings are
// CONSUMED by the parser at the runtime PDB-loading boundary and DIE
// at function return; the parser produces typed-enum outputs only.
// Safe to include from any runtime translation unit that links
// libnmr_shielding.
//
// The header is HAND-WRITTEN, not generated -- its contents are fixed
// (declarations + inline helpers), independent of the chemistry data
// emitted into the companion .cpp. Companion:
// src/generated/LegacyAmberSemanticTables.cpp (generator output).
//
// Per spec/plan/topology-encoding-dependencies-2026-05-05.md §H.
//

#include <cctype>
#include <cstdint>
#include <string>

#include "../SemanticEnums.h"

namespace nmr::topology_generated {

// ============================================================================
// Variant-index sentinel
// ============================================================================
//
// LookupBy selects per-residue tables by `variant_idx`. The chain
// (non-variant) table is selected by the sentinel `kBaseVariantIdx`;
// numeric values 0/1/2/... select the residue's variants in
// declaration order. Any other value returns nullptr (fail-fast on
// unknown indices, per dependencies §H.5).
//
inline constexpr std::uint8_t kBaseVariantIdx = 255;


// ============================================================================
// LookupBy -- structural lookup over the per-residue / per-variant tables
// ============================================================================
//
// Returns the AtomSemanticTable entry whose mechanical identity matches
// `identity` within the table selected by (residue, variant_idx).
// Returns nullptr if:
//   * no atom in the selected table has matching identity, or
//   * (residue, variant_idx) is not a valid combination (fail-fast).
//
// Mechanical-identity lookup is the canonical runtime mechanism;
// atom_local_idx integration is unsound and is forbidden -- see
// dependencies §H.
//
const ::nmr::AtomSemanticTable*
LookupBy(::nmr::AminoAcid residue,
         std::uint8_t variant_idx,
         const ::nmr::AtomMechanicalIdentity& identity);


// ============================================================================
// LookupCap -- structural lookup over the four terminal-state cap tables
// ============================================================================
//
// Returns the cap-table entry whose mechanical identity matches
// `identity` within the table selected by `state`. Returns nullptr
// if no atom matches.
//
// Cap-table entries are SEMANTIC DELTAS, not full rows. Use
// ApplyCapDelta() to compose with the chain entry; do not assign the
// cap entry whole-row over a chain entry.
//
const ::nmr::AtomSemanticTable*
LookupCap(::nmr::TerminalState state,
          const ::nmr::AtomMechanicalIdentity& identity);


// ============================================================================
// ApplyCapDelta -- canonical cap-on-chain composition primitive
// ============================================================================
//
// Compose a terminal-state cap delta on top of a chain-derived
// AtomSemanticTable. The runtime composition rule from dependencies
// §H.5 specifies that cap entries control ONLY the chemistry-level
// fields enumerated below; identity fields (Element, Locant,
// BranchAddress, DiastereotopicIndex, BackboneRole) and RDKit-derived
// fields (aromatic, equivalence_class) stay from chain, as do
// chemistry fields that do not vary at termini (prochiral,
// planar_stereo).
//
// The cap entries' non-delta fields carry placeholder values (typically
// the type's default-constructed state). Those values are NEVER read
// by ApplyCapDelta. This is by design: cap tables are space-efficient
// AtomSemanticTable rows where only the delta fields carry meaning.
//
// For cap-only atoms (H1/H2/H3 in NTERM, OXT/HXT in CTERM -- no chain
// entry), the runtime initialises `result` from the cap entry directly
// (the cap entry's identity fields are valid; RDKit-derived fields
// default to 0/false because cap atoms are not RDKit-perceived). Use
// ApplyCapDelta for the override case (chain entry exists AND cap
// entry exists for the same atom -- the terminal N at NTERM, or the
// terminal C/O at CTERM).
//
inline void ApplyCapDelta(::nmr::AtomSemanticTable& chain_then_result,
                          const ::nmr::AtomSemanticTable& cap_delta) {
    chain_then_result.planar_group  = cap_delta.planar_group;
    chain_then_result.polar_h       = cap_delta.polar_h;
    chain_then_result.formal_charge = cap_delta.formal_charge;
    chain_then_result.pseudoatom    = cap_delta.pseudoatom;
    chain_then_result.is_exchangeable =
        (chain_then_result.polar_h != ::nmr::PolarHKind::NotPolar);
    // NOT overridden by cap (preserved from chain):
    //   element, locant, branch, di_index, backbone_role (identity);
    //   aromatic, equivalence_class (RDKit perception);
    //   prochiral, planar_stereo (chemistry not affected by termini
    //   per Section 4 of the residue reference doc);
    //   ring_position. The chain row's ring membership is the truth.
    //   Cap atoms are NOT in any ring (cap entries carry
    //   RingSystemKind::None / RingPositionLabel::NotInRing as
    //   placeholder defaults). For PRO at NTERM, the chain N row
    //   correctly carries Pyrrolidine_Pro/ProRingNitrogen (the Pro-
    //   specific position label introduced in Slice A; replaces the
    //   generic Heteroatom_NoH that earlier versions used); the
    //   backbone cap-delta path must preserve that ring membership
    //   rather than clobber it with the cap entry's NotInRing
    //   default. (For non-Pro residues, chain N's ring_position is
    //   already NotInRing, so the bug was silent before — the PRO
    //   case surfaced it.)
}


// ============================================================================
// Atom-name parser -- mechanical fields from a PDB / AMBER atom name
// ============================================================================
//
// Lifted from tools/topology/build_semantic_tables.cpp::ParseAtomName so
// the runtime composition path and the generator share a single source
// of truth for the (Locant, BranchAddress, DiastereotopicIndex,
// BackboneRole) mapping. The function consumes std::string atom names
// at the load boundary and returns a typed `ParsedAtomName`; strings
// die at function return. No chemistry-string library is involved.
//
// Mirrors the documented locant convention from Markley, Bax, Arata,
// Hilbers, Kaptein, Sykes, Wright, Wuethrich, J. Biomol. NMR 12 (1998)
// 1-23 §2.1.1 (Greek-letter side-chain locants); BackboneRole assignment
// follows the peptide-amide chemistry of Pauling, Corey, Branson, PNAS
// 37 (1951) 205-211.
//
// The parser does NOT take Element as input; it derives Element from
// the first character of the atom name in the same way the generator
// does. Runtime callers that already have an authoritative Element on
// the runtime atom (PDB column-77 element field) should use
// `ComputeAtomMechanicalIdentity` below, which prefers that authority
// over name-derived element.
//
struct ParsedAtomName {
    ::nmr::Element             element       = ::nmr::Element::Unknown;
    bool                       is_backbone   = false;  ///< N/CA/C/O/H/HA (and HN alias).
    bool                       is_n_terminus = false;  ///< H1/H2/H3 (extra ammonium Hs).
    bool                       is_c_terminus = false;  ///< OXT/HXT.
    /// Cap-only-NTERM flag: this atom has NO entry in the residue's chain
    /// table and resolves only via LookupCap(NtermCharged|NtermNeutral).
    /// After 2026-05-06 (CharmmLegacy cleanup), this is 1:1 with
    /// `is_n_terminus`. The intent-explicit name lets ComposeAtomSemantic
    /// dispatch on `parsed.is_cap_only_n || parsed.is_cap_only_c` instead
    /// of re-reading the atom-name string via IsCapOnlyAtomName(string)
    /// (codex-review Finding 3 — "parse once, then no string work in
    /// composition").
    bool                       is_cap_only_n = false;
    /// Cap-only-CTERM flag: same intent, for OXT/HXT.
    bool                       is_cap_only_c = false;
    ::nmr::Locant              locant        = ::nmr::Locant::None;
    ::nmr::BranchAddress       branch        = {};
    ::nmr::DiastereotopicIndex di_index      = ::nmr::DiastereotopicIndex::None;
    ::nmr::BackboneRole        backbone_role = ::nmr::BackboneRole::None;
};

inline ::nmr::Element ElementFromAtomName(const std::string& name) {
    if (name.empty()) return ::nmr::Element::Unknown;
    switch (name[0]) {
        case 'H': return ::nmr::Element::H;
        case 'C': return ::nmr::Element::C;
        case 'N': return ::nmr::Element::N;
        case 'O': return ::nmr::Element::O;
        case 'S': return ::nmr::Element::S;
        default:  return ::nmr::Element::Unknown;
    }
}

inline ::nmr::Locant LocantLetterToEnum(char c) {
    // Note: 'H' as locant suffix (Arg/Tyr eta), NOT element.
    switch (c) {
        case 'A': return ::nmr::Locant::Alpha;
        case 'B': return ::nmr::Locant::Beta;
        case 'G': return ::nmr::Locant::Gamma;
        case 'D': return ::nmr::Locant::Delta;
        case 'E': return ::nmr::Locant::Epsilon;
        case 'Z': return ::nmr::Locant::Zeta;
        case 'H': return ::nmr::Locant::Eta;
        default:  return ::nmr::Locant::None;
    }
}

// Parse an atom name plus its parent name (the first heavy-atom
// neighbour, used to disambiguate H atoms with branch indices like
// HG21 on Thr CG2). For non-H atoms `parent_name` is empty.
//
// This implementation is a 1:1 transcription of
// tools/topology/build_semantic_tables.cpp::ParseAtomName circa
// commit ee1f1b4 ("structural-matching lookup + cap-table separation")
// and remains in sync with that file. Lifted into the runtime header
// so the StringBarrier discipline holds: there is no private copy in
// the generator anymore -- both the runtime and the generator call
// THIS function. If the parser changes, edit it here and
// build_semantic_tables.cpp picks up the change automatically (it
// includes this header).
inline ParsedAtomName ParseAtomName(const std::string& name,
                                    const std::string& parent_name) {
    ParsedAtomName p;
    p.element = ElementFromAtomName(name);

    // Pure backbone names. BackboneRole identifies which slot in the
    // peptide-amide unit this atom occupies (the (Element, Locant,
    // Branch, DiIndex) tuple all share Locant::None for these atoms,
    // so it cannot disambiguate them on its own). HA carries
    // BackboneRole::AlphaHydrogen AND Locant::Alpha; only Gly's
    // HA2/HA3 stay BackboneRole::None and rely on Locant::Alpha.
    if (name == "N" || name == "CA" || name == "C" || name == "O" ||
        name == "H" || name == "HN" || name == "HA") {
        p.is_backbone = true;
        if      (name == "N")  p.backbone_role = ::nmr::BackboneRole::Nitrogen;
        else if (name == "CA") p.backbone_role = ::nmr::BackboneRole::AlphaCarbon;
        else if (name == "C")  p.backbone_role = ::nmr::BackboneRole::CarbonylCarbon;
        else if (name == "O")  p.backbone_role = ::nmr::BackboneRole::CarbonylOxygen;
        else if (name == "H" || name == "HN") p.backbone_role = ::nmr::BackboneRole::AmideHydrogen;
        else if (name == "HA") {
            p.backbone_role = ::nmr::BackboneRole::AlphaHydrogen;
            p.locant = ::nmr::Locant::Alpha;  // HA is a Cα-H.
        }
        return p;
    }
    // N-terminal extra ammonium Hs (cap-only at NTERM_CHARGED /
    // NTERM_NEUTRAL — no chain-table entry).
    if (name == "H1" || name == "H2" || name == "H3") {
        p.is_n_terminus = true;
        p.is_cap_only_n = true;
        return p;
    }
    // C-terminal carboxyl atoms (cap-only at CTERM_DEPROTONATED /
    // CTERM_PROTONATED — no chain-table entry).
    if (name == "OXT" || name == "HXT") {
        p.is_c_terminus = true;
        p.is_cap_only_c = true;
        return p;
    }
    // NOTE: H2N (NTERM amine literal) and OT1/OT2 (CHARMM-style C-term
    // oxygens) were removed 2026-05-06 with the CharmmLegacy cleanup
    // (codex-review Finding 2). No active load path emits them; if a
    // future CHARMM-input path arises, they belong with a concrete
    // source tag in the canonicalisation rule layer (NamingApplicator),
    // not as silent terminus-flag aliases here.
    // Glycine alpha-Hs (only diastereotopic CH2 in the standard 20).
    if (name == "HA2" || name == "HA3") {
        p.is_backbone = true;
        p.locant = ::nmr::Locant::Alpha;
        p.di_index = (name.back() == '2') ? ::nmr::DiastereotopicIndex::Position2
                                          : ::nmr::DiastereotopicIndex::Position3;
        return p;
    }

    // Sidechain: <element-letter><locant-letter><digits...>.
    if (name.size() < 2) return p;
    p.locant = LocantLetterToEnum(name[1]);
    if (p.locant == ::nmr::Locant::None) return p;

    // Trailing digits form the numeric suffix.
    const std::string suffix = name.substr(2);
    if (suffix.empty()) return p;  // single side-chain atom (CB, SG, OH).

    // Heavy atoms with a digit suffix: digit is the BranchAddress::outer.
    if (p.element != ::nmr::Element::H) {
        if (suffix.size() == 1 && std::isdigit(static_cast<unsigned char>(suffix[0]))) {
            p.branch.outer = static_cast<uint8_t>(suffix[0] - '0');
            return p;
        }
        // Multi-digit on heavy atom is unusual; record as outer.
        p.branch.outer = static_cast<uint8_t>(suffix[0] - '0');
        return p;
    }

    // H atoms: use the heavy-atom parent map to disambiguate.
    bool parent_has_digit = false;
    char parent_digit = 0;
    if (!parent_name.empty()) {
        for (char c : parent_name) {
            if (std::isdigit(static_cast<unsigned char>(c))) {
                parent_has_digit = true;
                parent_digit = c;
                break;
            }
        }
    }

    if (!parent_has_digit) {
        // Diastereotopic methylene-style: HB2 vs HB3, HG2 vs HG3, etc.
        // Methyl-on-unbranched (Ala HB1/HB2/HB3) lands here too;
        // PseudoatomMembership downstream carries methyl-pair semantics.
        if (suffix.size() == 1 && std::isdigit(static_cast<unsigned char>(suffix[0]))) {
            const char d = suffix[0];
            if (d == '2') p.di_index = ::nmr::DiastereotopicIndex::Position2;
            else if (d == '3') p.di_index = ::nmr::DiastereotopicIndex::Position3;
            // d == '1' (Ala HB1) leaves di_index = None.
        }
        return p;
    }

    // Parent has branch digit: outer = parent's branch, inner = remaining
    // digit on the H name (if any).
    p.branch.outer = static_cast<uint8_t>(parent_digit - '0');
    for (char c : suffix) {
        if (!std::isdigit(static_cast<unsigned char>(c))) continue;
        if (c == parent_digit && p.branch.inner == 0 && suffix.size() == 1) {
            // Single matching digit (HG1 on CG1 = "the H on CG1").
            return p;
        }
        if (c != parent_digit) {
            p.branch.inner = static_cast<uint8_t>(c - '0');
            break;
        }
    }
    if (p.branch.inner == 0) {
        p.branch.inner = static_cast<uint8_t>(suffix.back() - '0');
        if (p.branch.inner == p.branch.outer && suffix.size() == 1) {
            p.branch.inner = 0;
        }
    }
    return p;
}


// ============================================================================
// ComputeAtomMechanicalIdentity -- mechanical-identity tuple from an atom
// ============================================================================
//
// Wraps `ParseAtomName` and prefers the runtime atom's authoritative
// Element over the name-derived element. The runtime carries Element as
// a typed enum from the PDB load boundary; using it here keeps Element
// authority on the typed runtime atom, not on the parser.
//
// Returns the typed `AtomMechanicalIdentity` consumed by `LookupBy` and
// `LookupCap`.
//
inline ::nmr::AtomMechanicalIdentity
ComputeAtomMechanicalIdentity(::nmr::Element element,
                              const std::string& atom_name,
                              const std::string& parent_atom_name) {
    const ParsedAtomName p = ParseAtomName(atom_name, parent_atom_name);
    ::nmr::AtomMechanicalIdentity ident;
    ident.element       = element;  // Authority: typed runtime Element.
    ident.locant        = p.locant;
    ident.branch        = p.branch;
    ident.di_index      = p.di_index;
    ident.backbone_role = p.backbone_role;
    return ident;
}


// ============================================================================
// IsCapOnlyAtomName -- partition standard chain atoms from cap-only atoms
// ============================================================================
//
// Cap-only atoms have NO chain-table entry and resolve only via
// `LookupCap`. The set is closed: NTERM_CHARGED contributes H1/H2/H3,
// NTERM_NEUTRAL contributes H1/H2, CTERM_DEPROTONATED contributes OXT,
// CTERM_PROTONATED contributes OXT/HXT.
//
// History note: H2N / OT1 / OT2 were previously listed here as
// CHARMM-port alternates; they were removed 2026-05-06 (codex-review
// Finding 2) — no active load path emits them, the corresponding
// CharmmLegacy rules had no live emitter, and CHARMM-the-force-field
// is retired (memory `project_charmm_retired_amber_only_2026-05-02`).
//
// Implementation note (codex-review Finding 3): this function is a
// THIN WRAPPER around `ParseAtomName` and consumes the typed
// `is_cap_only_n` / `is_cap_only_c` flags that ParseAtomName already
// populates in a single pass. Composition consumers (ComposeAtomSemantic)
// should `ParseAtomName` once per atom and read the typed flags
// directly, NOT call this string-taking wrapper. The wrapper exists
// for diagnostic / spot-check use cases where a parsed object isn't
// already in scope.
//
inline bool IsCapOnlyAtomName(const std::string& atom_name) {
    const ParsedAtomName p = ParseAtomName(atom_name, /*parent_name=*/"");
    return p.is_cap_only_n || p.is_cap_only_c;
}


// ============================================================================
// IsBackboneCapOverlayAtom -- decides which backbone atoms get cap-deltas
// ============================================================================
//
// Cap entries OVERRIDE chain-derived chemistry on a small set of
// backbone atoms at the terminus:
//   * NTERM (charged or neutral): backbone N gets cap-delta
//     (formal_charge, polar_h on the new H1/H2/H3, planar_group).
//   * CTERM (deprotonated or protonated): backbone C and O get cap-delta
//     (planar_group flips PeptideAmide -> Carboxylate; on OXT the
//     formal_charge flips per state).
//
// Returns true if `parsed.backbone_role` is one of those slots and the
// terminus is the corresponding state. Used by FinalizeConstruction to
// decide whether to call ApplyCapDelta after LookupBy.
//
// Glycine HA2/HA3 keep BackboneRole::None and are intentionally NOT
// cap-overridden (their chemistry is not affected by the terminus).
//
inline bool IsBackboneCapOverlayAtom(::nmr::BackboneRole role,
                                     ::nmr::TerminalState n_state,
                                     ::nmr::TerminalState c_state) {
    const bool n_cap = (n_state == ::nmr::TerminalState::NtermCharged ||
                        n_state == ::nmr::TerminalState::NtermNeutral);
    const bool c_cap = (c_state == ::nmr::TerminalState::CtermDeprotonated ||
                        c_state == ::nmr::TerminalState::CtermProtonated);
    if (n_cap && role == ::nmr::BackboneRole::Nitrogen) return true;
    if (c_cap && (role == ::nmr::BackboneRole::CarbonylCarbon ||
                  role == ::nmr::BackboneRole::CarbonylOxygen)) return true;
    return false;
}


}  // namespace nmr::topology_generated
