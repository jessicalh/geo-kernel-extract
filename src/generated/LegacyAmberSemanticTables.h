#pragma once
//
// LegacyAmberSemanticTables.h -- runtime API for the typed
// chemistry-substrate tables emitted by tools/topology/build_semantic_tables.
//
// String barrier: this header contains only typed-enum declarations
// and an inline composition helper. No chemistry-string libraries
// (cifpp, RDKit, gemmi) are referenced. Safe to include from any
// runtime translation unit that links libnmr_shielding.
//
// The header is HAND-WRITTEN, not generated -- its contents are fixed
// (declarations + one inline helper), independent of the chemistry
// data emitted into the companion .cpp. Companion:
// src/generated/LegacyAmberSemanticTables.cpp (generator output).
//
// Per spec/plan/topology-encoding-dependencies-2026-05-05.md §H.
//

#include <cstdint>

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
    chain_then_result.ring_position = cap_delta.ring_position;
    chain_then_result.is_exchangeable =
        (chain_then_result.polar_h != ::nmr::PolarHKind::NotPolar);
    // NOT overridden by cap (preserved from chain):
    //   element, locant, branch, di_index, backbone_role (identity);
    //   aromatic, equivalence_class (RDKit perception);
    //   prochiral, planar_stereo (chemistry not affected by termini
    //   per Section 4 of the residue reference doc).
}

}  // namespace nmr::topology_generated
