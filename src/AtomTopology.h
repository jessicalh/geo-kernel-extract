#pragma once
//
// AtomTopology: per-atom symbolic topology stamped from the IUPAC-IUBMB-IUPAB
// nomenclature recommendations (Markley et al. 1998, J. Biomol. NMR 12:1-23).
//
// Populated at protein load time by Protein::StampAtomTopology() against the
// AminoAcidType table. Invariant for the protein's lifetime; describes what
// the atom IS in IUPAC structural terms, independent of geometry.
//
// Two fields encode prochirality because they are NOT the same concept:
//   diastereotopic_index — IUPAC 2/3 numbering. PURELY STRUCTURAL.
//                          Determined by main-chain-precedence + clockwise
//                          rule (Figure 1 caption). Does NOT alternate with
//                          chain depth.
//   prochiral_stereo     — Cahn-Ingold-Prelog pro-R / pro-S descriptor.
//                          ALTERNATES with chain depth: Hβ3 is pro-R for most
//                          residues but Hγ3 is pro-S, Hδ3 is pro-R again, etc.
//
// PDB / BMRB use the IUPAC 2/3 names ("HB2", "HB3"). The CIP descriptor is a
// separate axis. Conflating them produces silently wrong values at half the
// chain depths.
//
// Reference PDF: references/markley-1998-iupac-nmr-nomenclature-recommendations.pdf
//

#include <array>
#include <cstdint>

namespace nmr {

// Greek-letter atom locant (Markley 1998 Sec 2.1.1, Figure 1).
// Roman fallback: A=Alpha, B=Beta, G=Gamma, D=Delta, E=Epsilon, Z=Zeta, H=Eta.
// PDB strings use the Roman codes (HB2, CG, ...).
enum class Locant : uint8_t {
    Backbone,   // N, HN/H, C, O — peptide main chain (no Greek position)
    Alpha,      // Cα, Hα — start of the Greek-letter chain
    Beta,
    Gamma,
    Delta,
    Epsilon,
    Zeta,
    Eta,
    Terminal,   // H1/H2/H3 (N-term amine), O'/O''/H'' (C-term carboxyl)
    Special     // unusual / not covered by IUPAC tabulation
};

// Branch index for atoms with branched parent: Cγ1/Cγ2 (Val, Ile), Cδ1/Cδ2
// (Leu, Phe, Tyr, Trp, His), Nη1/Nη2 (Arg), Cε3/Cζ3 (Trp), etc. None for
// atoms without a sibling at the same Greek position. Three appears in Trp
// (Cε3, Cζ3) where the indole numbering produces a third branch.
enum class BranchIndex : uint8_t { None, One, Two, Three };

// IUPAC 2/3 numbering for diastereotopic substituents (Hβ2/Hβ3, Hγ2/Hγ3, etc.).
// Determined by main-chain-precedence + clockwise rule (Figure 1 caption).
// Purely structural; does not alternate with chain depth. PDB strings carry
// these numbers verbatim ("HB2", "HB3").
enum class DiastereotopicIndex : uint8_t {
    None,    // not at a diastereotopic position
    Two,     // IUPAC index 2 (HB2, HG2, ...)
    Three    // IUPAC index 3 (HB3, HG3, ...)
};

// Cahn-Ingold-Prelog pro-R / pro-S descriptor.
// CRITICAL: alternates with chain depth. The mapping from IUPAC 2/3 to CIP
// pro-R/pro-S is NOT constant across the chain. Hand-curated per the
// geometric rule applied at each prochiral center.
enum class ProchiralStereo : uint8_t {
    NotProchiral,
    Equivalent,    // chemically equivalent under fast averaging (methyl Hs, etc.)
    ProR,
    ProS
};

// Cis/trans (E/Z) descriptor for amide hydrogens on Asn / Gln / Arg side-chain
// NH2 groups. Markley 1998 Section 2.1.1: "distinguished by numbers (1 or 2)
// on the basis of their relationship (cis or trans, respectively) to the
// heavy atom three bonds closer to the main chain."
//   Asn Hδ21: cis to Cβ; Hδ22: trans
//   Gln Hε21: cis to Cγ; Hε22: trans
//   Arg Hη11/Hη12 (on Nη1, cis to Cδ via Nε); Hη21/Hη22 (on Nη2, trans).
//   Arg Nη1 itself is cis (Z) to Cδ; Nη2 is trans (E).
enum class PlanarStereo : uint8_t {
    None,
    Cis,    // index 1 — cis (Z) to reference atom
    Trans   // index 2 — trans (E) to reference atom
};

// IUPAC pseudoatom class (Markley 1998 Sec 3.1, Table 1).
// Identifies symmetry-equivalent groups under fast averaging.
//   M = methyl (3 chemically equivalent Hs on a single carbon)
//   Q = other (methylene pair, NH2 pair, ring super-class, etc.)
enum class PseudoatomClass : uint8_t {
    None,
    // Methyls (all three Hs on one carbon are equivalent)
    MB,    // β-methyl (Ala HB1/HB2/HB3)
    MG,    // γ-methyl (Thr HG21/22/23, Ile HG21/22/23)
    MG1,   // γ1-methyl (Val HG11/12/13)
    MG2,   // γ2-methyl (Val HG21/22/23)
    MD,    // δ-methyl (Ile HD11/12/13)
    MD1,   // δ1-methyl (Leu HD11/12/13)
    MD2,   // δ2-methyl (Leu HD21/22/23)
    ME,    // ε-methyl (Met HE1/2/3)
    // Methylenes / amino groups
    QA,    // α-methylene (Gly Hα2/Hα3)
    QB,    // β-methylene (generic Hβ2/Hβ3 super-class)
    QG,    // γ-methylene or all-γ super-class
    QD,    // δ-methylene or full δ super-class (Leu QD = all six δ-methyl Hs)
    QE,    // ε-methylene or NH2 (Gln Hε21/Hε22)
    QZ,    // ζ-amino (Lys NH3+ Hζ1/2/3)
    // Arg guanidinium NH2 groups
    QH1,   // Arg Hη11/Hη12 (one nitrogen's NH2 pair)
    QH2,   // Arg Hη21/Hη22 (other nitrogen's pair)
    QH,    // Arg all four Hη super-class
    // Ring super-class
    QR     // Phe/Tyr all-ring super-class
};

// Role of an atom within an aromatic ring.
enum class RingPosition : uint8_t {
    NotInRing,
    Substituent,  // ring atom bonded to a non-ring sidechain atom
                  //   (Phe Cγ, Tyr Cγ, His Cγ, Trp Cγ)
    Member,       // generic ring vertex
    Junction      // atom at a ring fusion (Trp Cδ2 / Cε2 are shared between
                  //   TRP5 / TRP6 / TRP9)
};

// Kind of polar / labile / heteroatom-bonded hydrogen.
//
// Set on H atoms whose chemistry is materially distinct from generic C-H:
// ring N-H of His variants and Trp indole, hydroxyl O-H of Ser/Thr/Tyr,
// carboxyl acid O-H of Ash/Glh, sulfanyl S-H of Cys. Default None covers
// all C-H atoms and the backbone/side-chain amide H atoms whose role is
// already captured by EnrichmentResult flags (`is_amide_H` etc.).
//
// Variant-conditional Hs (HD1/HE2 of HIS variants, HD2 of ASP-ASH, HE2 of
// GLU-GLH) carry the chemistry they exhibit when present — the topology
// is invariant given existence; only existence depends on the variant.
enum class PolarHKind : uint8_t {
    None,            // generic C-H or backbone/sidechain amide N-H (use atom_flags)
    RingNH,          // His HD1 (HID/HIP) / HE2 (HIE/HIP); Trp HE1 (indole NH)
    Hydroxyl,        // Ser HG; Thr HG1; Tyr HH (alcohol / phenol O-H)
    CarboxylAcid,    // Asp HD2 (ASH); Glu HE2 (GLH) — protonated carboxylate
    Sulfanyl,        // Cys HG (thiol S-H; absent in CYX/CYM variants)
};

// Per-atom symbolic topology. ~16 bytes. One value per Atom, populated at
// protein load via name lookup; const after stamping.
struct AtomTopology {
    Locant locant = Locant::Special;
    BranchIndex branch_index = BranchIndex::None;
    DiastereotopicIndex diastereotopic_index = DiastereotopicIndex::None;
    ProchiralStereo prochiral_stereo = ProchiralStereo::NotProchiral;
    PlanarStereo planar_stereo = PlanarStereo::None;
    PseudoatomClass pseudoatom_class = PseudoatomClass::None;
    RingPosition ring_position = RingPosition::NotInRing;
    PolarHKind polar_h_kind = PolarHKind::None;

    // Chi participation: chi_position[i] is this atom's 0..3 position in
    // chi(i+1), or -1 if not in that chi. An atom can participate in
    // multiple chi angles (e.g., Lys CB is in chi1 at position 2 and
    // chi2 at position 1).
    std::array<int8_t, 4> chi_position = {-1, -1, -1, -1};

    // True after a successful StampAtomTopology lookup. False = the
    // (residue_type, atom_name) pair was not found in the IUPAC table —
    // a loud diagnostic was emitted at load and this atom needs attention.
    bool stamped = false;
};

// ---------------------------------------------------------------------------
// constexpr builder helpers — terse table construction.
// ---------------------------------------------------------------------------

constexpr AtomTopology TopologyBackbone() {
    AtomTopology t;
    t.locant = Locant::Backbone;
    t.stamped = true;
    return t;
}

constexpr AtomTopology TopologyAlpha() {
    AtomTopology t;
    t.locant = Locant::Alpha;
    t.stamped = true;
    return t;
}

constexpr AtomTopology TopologyAlphaProchiral(DiastereotopicIndex idx,
                                               ProchiralStereo cip) {
    AtomTopology t = TopologyAlpha();
    t.diastereotopic_index = idx;
    t.prochiral_stereo = cip;
    t.pseudoatom_class = PseudoatomClass::QA;
    return t;
}

constexpr AtomTopology TopologySidechain(Locant loc,
                                          BranchIndex br = BranchIndex::None) {
    AtomTopology t;
    t.locant = loc;
    t.branch_index = br;
    t.stamped = true;
    return t;
}

constexpr AtomTopology TopologyMethyl(Locant loc, BranchIndex br,
                                       PseudoatomClass pseudo) {
    AtomTopology t;
    t.locant = loc;
    t.branch_index = br;
    t.prochiral_stereo = ProchiralStereo::Equivalent;
    t.pseudoatom_class = pseudo;
    t.stamped = true;
    return t;
}

constexpr AtomTopology TopologyDiastereo(Locant loc, DiastereotopicIndex idx,
                                          ProchiralStereo cip,
                                          PseudoatomClass pseudo = PseudoatomClass::None,
                                          BranchIndex br = BranchIndex::None) {
    AtomTopology t;
    t.locant = loc;
    t.branch_index = br;
    t.diastereotopic_index = idx;
    t.prochiral_stereo = cip;
    t.pseudoatom_class = pseudo;
    t.stamped = true;
    return t;
}

constexpr AtomTopology TopologyRing(Locant loc, BranchIndex br,
                                     RingPosition pos,
                                     PseudoatomClass pseudo = PseudoatomClass::None) {
    AtomTopology t;
    t.locant = loc;
    t.branch_index = br;
    t.ring_position = pos;
    t.pseudoatom_class = pseudo;
    t.stamped = true;
    return t;
}

constexpr AtomTopology TopologyPlanar(Locant loc, BranchIndex br,
                                       PlanarStereo stereo,
                                       PseudoatomClass pseudo = PseudoatomClass::None) {
    AtomTopology t;
    t.locant = loc;
    t.branch_index = br;
    t.planar_stereo = stereo;
    t.pseudoatom_class = pseudo;
    t.stamped = true;
    return t;
}

constexpr AtomTopology TopologyTerminalH() {
    AtomTopology t;
    t.locant = Locant::Terminal;
    t.prochiral_stereo = ProchiralStereo::Equivalent;
    t.stamped = true;
    return t;
}

// Polar / heteroatom-bonded H — heteroatom H atoms with chemistry distinct
// from generic sidechain C-H. The kind argument is the value picked up by
// post-hoc analysis when classifying NMR resonances; locant + branch carry
// the structural identity. Use TopologySidechain for ordinary C-H atoms.
constexpr AtomTopology TopologyPolarH(Locant loc, BranchIndex br,
                                       PolarHKind kind) {
    AtomTopology t;
    t.locant = loc;
    t.branch_index = br;
    t.polar_h_kind = kind;
    t.stamped = true;
    return t;
}

constexpr AtomTopology TopologyTerminalO() {
    AtomTopology t;
    t.locant = Locant::Terminal;
    t.stamped = true;
    return t;
}

}  // namespace nmr
