// QtSemanticEnums.h — chemistry-substrate typed vocabulary.
//
// Mirrors `nmr/src/SemanticEnums.h` field-for-field, ordinal-
// compatible. The library splits core typed primitives (`Types.h`)
// from chemistry substrate (`SemanticEnums.h`) for module-boundary
// reasons; we keep the same split because the reader echoes the
// library's organisation.
//
// Sidecar storage for these enums is `atoms_category_info.npy`'s int8
// columns. Manifest `enum_vocab` declares `terminal_state` ordinals;
// the other substrate enums are not in the manifest vocab block
// (no library-side enum_vocab key for them), so the loader trusts the
// hardcoded ordinals here matched against the library's
// SemanticEnums.h header.
//
// No-strings discipline: zero string fields in this file. Every value
// is a typed enum or POD struct with typed-enum members. The
// `Marker`/`citation` style of the library file (long docstrings
// citing chemistry literature) is preserved as comments — these are
// design intent, not runtime overhead.

#pragma once

#include "Types.h"  // Element

#include <array>
#include <cstdint>

namespace h5reader::model {


// ============================================================================
// TerminalState — chain-end residue classification.
//
// Library mirror: nmr::ResidueTerminalState (src/Residue.h:17) — the
// PROJECTION variant used by residues.npy + atoms_category_info.npy.
// Manifest enum_vocab confirms: {0:Internal, 1:NTerminus, 2:CTerminus,
// 3:NAndCTerminus, 4:Unknown}.
//
// Note: nmr/src/SemanticEnums.h has a finer-grained nmr::TerminalState
// (NtermCharged, NtermNeutral, CtermDeprotonated, CtermProtonated)
// used internally for the four cap-atom tables. The h5-reader
// consumes the projection enum from the manifest, not the internal
// variant; both share the {Internal, ...Term, Unknown} positions but
// differ in the ordinals 1-4. Use the projection ordinals here.
// ============================================================================

enum class TerminalState : int8_t {
    Internal = 0,
    NTerminus = 1,
    CTerminus = 2,
    NAndCTerminus = 3,
    Unknown = 4,
};


// ============================================================================
// BackboneRole — which canonical backbone slot an atom occupies.
// Library mirror: nmr::BackboneRole (src/SemanticEnums.h:89).
// Sidecar storage: atoms_category_info.backbone_role (int8).
//
// Pauling, Corey, Branson, PNAS 37 (1951) 205-211 (planar peptide
// bond) + Ramachandran & Sasisekharan, Adv. Protein Chem. 23 (1968).
// Glycine HA2/HA3 carry `Locant::Alpha` + `BackboneRole::None` —
// disambiguated by (Locant, DiastereotopicIndex), not BackboneRole.
// ============================================================================

enum class BackboneRole : int8_t {
    None = 0,  // sidechain atom or cap
    Nitrogen = 1,
    AlphaCarbon = 2,
    CarbonylCarbon = 3,
    CarbonylOxygen = 4,
    AmideHydrogen = 5,
    AlphaHydrogen = 6,  // non-GLY HA only; GLY HA2/HA3 stay None
};


// ============================================================================
// Locant — Greek-letter / IUPAC sidechain position label.
// Library mirror: nmr::Locant (src/SemanticEnums.h:119).
// Sidecar storage: atoms_category_info.locant (int8).
//
// Markley, Bax, Arata, Hilbers, Kaptein, Sykes, Wright, Wuethrich,
// J. Biomol. NMR 12 (1998) 1-23 §2.1.1. The Greek letters walk
// outward from C-alpha along the side-chain heavy-atom chain.
// Backbone atoms (N, CA, C, O, H, HA) carry Locant::None.
// ============================================================================

enum class Locant : int8_t {
    None = 0,
    Alpha = 1,
    Beta = 2,
    Gamma = 3,
    Delta = 4,
    Epsilon = 5,
    Zeta = 6,
    Eta = 7,
};


// ============================================================================
// BranchAddress — two-level disambiguator for atoms sharing a locant.
// Library mirror: nmr::BranchAddress (src/SemanticEnums.h:152).
// Sidecar storage: atoms_category_info.branch_outer / branch_inner
// (int8 each). For 99% of atoms only `outer` is non-zero; arginine
// side-chain H-eta atoms use both indices.
//
// Markley 1998 Figure 1 caption: CIP-clockwise rule for sight-down-
// highest-priority-substituent disambiguation.
// ============================================================================

struct BranchAddress {
    int8_t outer = 0;  // heavy-atom branch (1, 2). 0 = not applicable.
    int8_t inner = 0;  // H-within-group branch (1, 2). 0 = not applicable.

    constexpr bool IsBranched() const { return outer != 0; }
    constexpr bool HasInner() const { return inner != 0; }
};
constexpr bool operator==(BranchAddress a, BranchAddress b) {
    return a.outer == b.outer && a.inner == b.inner;
}


// ============================================================================
// DiastereotopicIndex — IUPAC 2/3 label on prochiral methylene Hs.
// Library mirror: nmr::DiastereotopicIndex (src/SemanticEnums.h:178).
// Sidecar storage: atoms_category_info.di_index (int8).
// Markley 1998 Figure 1 caption.
// ============================================================================

enum class DiastereotopicIndex : int8_t {
    None = 0,
    Position2 = 2,
    Position3 = 3,
};


// ============================================================================
// ProchiralStereo — CIP R/S designation for prochiral substituents.
// Library mirror: nmr::ProchiralStereo (src/SemanticEnums.h:206).
// Sidecar storage: atoms_category_info.prochiral (int8).
// Cahn, Ingold, Prelog Angew. Chem. Int. Ed. 5 (1966) 385; RDKit
// CIPLabeler is the algorithmic source.
// ============================================================================

enum class ProchiralStereo : int8_t {
    NotProchiral = 0,
    ProR = 1,
    ProS = 2,
    Unassigned = 3,
};


// ============================================================================
// PlanarGroupKind — which planar/sp2 functional group an atom is in.
// Library mirror: nmr::PlanarGroupKind (src/SemanticEnums.h:231).
// Sidecar storage: atoms_category_info.planar_group (int8).
// ============================================================================

enum class PlanarGroupKind : int8_t {
    None = 0,
    PeptideAmide = 1,      // backbone peptide bond plane
    SidechainAmide = 2,    // ASN, GLN
    Guanidinium = 3,       // ARG
    Imidazole = 4,         // HIS (variant-dependent)
    Aromatic6Ring = 5,     // PHE, TYR, TRP benzene
    Aromatic5Ring = 6,     // TRP pyrrole, HIS variant 5-ring
    Carboxylate = 7,       // ASP, GLU, C-term
    AromaticHydroxyl = 8,  // TYR -OH
    AromaticOxide = 9,     // TYM phenolate
};


// ============================================================================
// PlanarStereo — canonical IUPAC E/Z label at planar centres.
// Library mirror: nmr::PlanarStereo (src/SemanticEnums.h:319).
// Sidecar storage: atoms_category_info.planar_stereo (int8).
// ============================================================================

enum class PlanarStereo : int8_t {
    NotApplicable = 0,
    E = 1,
    Z = 2,
    Unspecified = 3,
};


// ============================================================================
// PseudoatomKind — Markley 1998 Table 1 pseudoatom letter.
// Library mirror: nmr::PseudoatomKind (src/SemanticEnums.h:354).
// Sidecar storage: atoms_category_info.pseudoatom_kind (int8).
// Markley 1998 Table 1 IUPAC pseudoatom taxonomy: M / Q / R only.
// ============================================================================

enum class PseudoatomKind : int8_t {
    None = 0,
    M = 1,  // methyl group member
    Q = 2,  // equivalent-H group member (non-methyl)
    R = 3,  // ring-atom group member
};


// ============================================================================
// PseudoatomMembership — atom's pseudoatom-group membership record.
// Library mirror: nmr::PseudoatomMembership (src/SemanticEnums.h:376).
//
// Reader-side projection on QtAtom uses two fields:
// `pseudoatomKind` (the kind) and `inSuperGroup` (whether the atom is
// also a member of a higher-order Q aggregator like Val QG, Leu QD,
// Arg QH). The locant + branch are not duplicated here — they live on
// QtAtom directly.
// ============================================================================


// ============================================================================
// PolarHKind — functional-group classification of polar hydrogens.
// Library mirror: nmr::PolarHKind (src/SemanticEnums.h:406).
// Sidecar storage: atoms_category_info.polar_h_kind (int8).
// Wuethrich 1986; Englander 2008 for H/D exchange framework.
// SHIFTX2 (Han et al. 2011) trains separate models per environment.
// ============================================================================

enum class PolarHKind : int8_t {
    NotPolar = 0,
    BackboneAmide = 1,
    SidechainPrimaryAmide = 2,
    IndoleNH = 3,
    AmmoniumNH = 4,
    GuanidiniumNH = 5,
    ImidazoleNH = 6,
    CarboxylOH = 7,
    HydroxylOH_Aliphatic = 8,
    HydroxylOH_Aromatic = 9,
    ThiolSH = 10,
    AmineNH = 11,
    OtherPolarH = 12,
};


// ============================================================================
// RingSystemKind — which named ring system an atom belongs to.
// Library mirror: nmr::RingSystemKind (src/SemanticEnums.h:496).
// Sidecar projection: indirectly via atoms_category_info's
// ring_position_primary / _secondary fields linking atoms to ring
// rows in rings.npy. Not directly a column.
// ============================================================================

enum class RingSystemKind : int8_t {
    NotInRing = 0,
    Benzene_Phe = 1,
    Benzene_Tyr = 2,
    Imidazole_His = 3,
    Indole_Trp_5 = 4,
    Indole_Trp_6 = 5,
    Pyrrolidine_Pro = 6,
    Indole_Trp_9 = 7,  // 9-atom indole perimeter (Case 1995)
};


// ============================================================================
// RingPositionLabel — position of an atom within its ring system.
// Library mirror: nmr::RingPositionLabel (src/SemanticEnums.h:533).
// Sidecar storage: atoms_category_info.ring_position_primary,
// ring_position_secondary (int8 each).
// Vollhardt & Schore (benzene); Joule & Mills (heterocycles).
// ============================================================================

enum class RingPositionLabel : int8_t {
    NotInRing = 0,
    Ipso = 1,
    Ortho1 = 2,
    Ortho2 = 3,
    Meta1 = 4,
    Meta2 = 5,
    Para = 6,
    PyrroleAlpha = 7,
    PyrroleBeta = 8,
    BridgeFusion = 9,
    Heteroatom_NH = 10,
    Heteroatom_NoH = 11,
    Heteroatom_OH = 12,
    Saturated = 13,
    ProRingNitrogen = 14,
    ProRingAlphaCarbon = 15,
    ProRingBeta = 16,
    ProRingPuckerPivot = 17,
    ProRingDelta = 18,
    PerimeterMember = 19,
};


// ============================================================================
// RingMembership — one ring's worth of context for an atom.
// Library mirror: nmr::RingMembership (src/SemanticEnums.h:643).
//
// The reader's `QtRingMembership` (in QtRingMembership.h) is the
// per-(ring, vertex-atom) record from ring_membership.npy — a
// different shape from the library's `RingMembership` (which is an
// atom-side view). We don't expose this library shape on the reader
// side because the atoms_category_info projection already collapses
// it into `ring_position_primary/secondary` labels on QtAtom.
// ============================================================================


// ============================================================================
// AtomMechanicalIdentity — the 5-tuple lookup key.
// Library mirror: nmr::AtomMechanicalIdentity (src/SemanticEnums.h:985).
//
// Reader doesn't need to LOOK UP semantics by mechanical identity
// (loader already populated them), but the tuple is useful for
// equality checks across atoms (e.g., "do these two atoms share a
// mechanical identity, modulo branch?").
// ============================================================================

struct AtomMechanicalIdentity {
    Element element = Element::Unknown;
    Locant locant = Locant::None;
    BranchAddress branch = {};
    DiastereotopicIndex di_index = DiastereotopicIndex::None;
    BackboneRole backbone_role = BackboneRole::None;
};

constexpr bool operator==(const AtomMechanicalIdentity& a, const AtomMechanicalIdentity& b) {
    return a.element == b.element && a.locant == b.locant && a.branch == b.branch && a.di_index == b.di_index
           && a.backbone_role == b.backbone_role;
}


}  // namespace h5reader::model
