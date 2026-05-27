// QtNpyRecords — packed POD structs that mirror the byte layout of the
// 5 sidecar NPY files exactly. Each struct's sizeof() matches the
// writer-side kXRecordSize constant (src/TopologySidecar.cpp +
// src/CategoryInfoProjection.cpp) so a raw memcpy from NPY data into a
// std::vector<QtNpyXRow> reconstitutes the typed rows.
//
// The fixed-width string columns (chain_id S2, atom_name S8, etc.) are
// raw char[N] arrays here — they are zero-padded by the writer. The
// loader's job is to char[N] → QString conversion at the point where
// the Qt model is populated (not here; this is the raw-byte boundary).
//
// IMPORTANT: every struct uses `#pragma pack(push, 1)` so the compiler
// does NOT inject padding. The static_asserts at the bottom of this
// header are the build-time guard: if the byte layout drifts from the
// writer (e.g. someone adds a field to one side only), the build fails
// here rather than mis-decoding silently.

#pragma once

#include <cstddef>
#include <cstdint>

namespace h5reader::io {

#pragma pack(push, 1)

// ============================================================================
// residues.npy — 42 bytes per row, 20 fields.
// Writer: src/TopologySidecar.cpp:104..127 (kResiduesDtypeDescr +
//         kResidueRecordSize).
// ============================================================================

struct QtNpyResidueRow {
    int32_t residue_index;          // 4
    char chain_id[2];               // S2
    int32_t residue_number;         // 4
    char insertion_code[1];         // S1
    int8_t residue_type;            // AminoAcid ordinal
    char amber_residue_3letter[4];  // S4 (e.g. "ALA", "HIE", "CYX")
    char iupac_residue_3letter[4];  // S4
    char one_letter[1];             // S1
    int8_t protonation_variant_index;
    int8_t terminal_state;       // TerminalState ordinal
    int32_t prev_residue_index;  // -1 if none
    int32_t next_residue_index;  // -1 if none
    int8_t prev_residue_type;
    int8_t next_residue_type;
    int32_t atom_count;
    int8_t is_proline;
    int8_t is_aromatic;
    int8_t is_titratable;
    int8_t has_amide_h;
    int8_t is_xpro_context;
};


// ============================================================================
// bonds.npy — 18 bytes per row, 9 fields.
// Writer: src/TopologySidecar.cpp:251..262.
// ============================================================================

struct QtNpyBondRow {
    int32_t bond_index;    // 4
    int32_t atom_index_a;  // 4
    int32_t atom_index_b;  // 4
    int8_t bond_order;     // BondOrder ordinal
    int8_t bond_category;  // BondCategory ordinal
    int8_t is_rotatable;
    int8_t is_aromatic;
    int8_t is_peptide;
    int8_t is_backbone;
};


// ============================================================================
// rings.npy — 24 bytes per row, 9 fields (one is a _pad0 alignment byte).
// Writer: src/TopologySidecar.cpp:303..314.
// ============================================================================

struct QtNpyRingRow {
    int32_t ring_id;         // absolute row index (4)
    int8_t ring_kind;        // 0 = aromatic, 1 = saturated
    int8_t ring_type_index;  // RingTypeIndex ordinal
    int8_t atom_count;
    int8_t _pad0;               // alignment to int32
    int32_t native_axis_index;  // index within aromatic-only or saturated-only axis
    int32_t parent_residue_index;
    int32_t parent_residue_number;
    int32_t fused_partner_ring_id;  // -1 if not fused
};


// ============================================================================
// ring_membership.npy — 12 bytes per row, 6 fields.
// Writer: src/TopologySidecar.cpp:392..400.
// ============================================================================

struct QtNpyRingMembershipRow {
    int32_t ring_id;
    int32_t atom_index;
    int8_t ring_atom_order;
    int8_t is_vertex;
    int8_t is_substituent;
    int8_t _pad0;
};


// ============================================================================
// atoms_category_info.npy — 83 bytes per row, 37 fields.
// Writer: src/CategoryInfoProjection.cpp WriteStructuredNpy + the typed
// field layout it composes per atom. (No single descriptor literal there;
// the field order is the literal sequence of memcpys.)
//
// Field order verified against the manifest dump on the 1P9J fixture
// 2026-05-23 — see the conversation transcript for the h5py dtype dump.
// ============================================================================

struct QtNpyAtomCategoryRow {
    int32_t atom_index;             // 4
    int32_t residue_index;          // 4
    int8_t element;                 // atomic number (1/6/7/8/16), NOT enum ordinal
    char amber_atom_name[8];        // S8
    char iupac_atom_name[8];        // S8
    char bmrb_atom_name[8];         // S8
    char amber_residue_3letter[4];  // S4
    char iupac_residue_3letter[4];  // S4
    char bmrb_residue_3letter[4];   // S4
    char residue_1letter[1];        // S1
    int8_t residue_type;            // AminoAcid ordinal
    int8_t residue_variant_index;
    int8_t terminal_state;  // TerminalState ordinal
    int8_t locant;          // Locant ordinal
    int8_t branch_outer;
    int8_t branch_inner;
    int8_t di_index;                 // DiastereotopicIndex ordinal
    int8_t backbone_role;            // BackboneRole ordinal
    int8_t prochiral;                // ProchiralStereo ordinal
    int8_t planar_group;             // PlanarGroupKind ordinal
    int8_t planar_stereo;            // PlanarStereo ordinal
    int8_t polar_h_kind;             // PolarHKind ordinal
    int8_t ring_position_primary;    // RingPositionLabel ordinal
    int8_t ring_position_secondary;  // RingPositionLabel ordinal
    int8_t pseudoatom_kind;          // PseudoatomKind ordinal
    int8_t in_super_group;
    int8_t aromatic;
    int8_t formal_charge;
    int8_t is_exchangeable;
    int8_t iupac_naming_provenance;  // NamingProvenance ordinal
    int8_t bmrb_naming_provenance;   // NamingProvenance ordinal
    char chain_id[2];                // S2
    int32_t residue_number;          // 4
    char insertion_code[1];          // S1
    int32_t parent_atom_index;       // 4 (-1 if not H)
    char ff_atom_type_string[4];     // S4
    int8_t equivalence_class;
};

#pragma pack(pop)

// ── Build-time guards: byte layouts must match the writer constants ──
// Sizes verified 2026-05-23 against the 1P9J baseline fixture via
// numpy.dtype.itemsize. Writer-side kXRecordSize constants match.
static_assert(sizeof(QtNpyResidueRow) == 42, "QtNpyResidueRow size diverged from src/TopologySidecar.cpp kResidueRecordSize");
static_assert(sizeof(QtNpyBondRow) == 18, "QtNpyBondRow size diverged from src/TopologySidecar.cpp kBondRecordSize");
static_assert(sizeof(QtNpyRingRow) == 24, "QtNpyRingRow size diverged from src/TopologySidecar.cpp kRingRecordSize");
static_assert(sizeof(QtNpyRingMembershipRow) == 12,
              "QtNpyRingMembershipRow size diverged from src/TopologySidecar.cpp kRingMembershipRecordSize");
static_assert(sizeof(QtNpyAtomCategoryRow) == 83,
              "QtNpyAtomCategoryRow size diverged from src/CategoryInfoProjection.cpp atom-row layout");

}  // namespace h5reader::io
