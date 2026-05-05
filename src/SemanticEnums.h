#pragma once
//
// SemanticEnums.h -- typed atom-level chemistry vocabulary for
// LegacyAmberTopology.
//
// This header defines 14 typed-enum / typed-struct fields that
// together specify the chemistry-substrate classification of every
// atom in a standard amino-acid residue. Values are populated at
// `Protein::FinalizeConstruction` from a generated static table
// (`src/generated/LegacyAmberSemanticTables.cpp`); after that point
// every chemistry question on an atom resolves to a typed lookup,
// not a string comparison.
//
// String barrier:
// - This header includes only <array>, <cstdint>. No cifpp, no RDKit,
//   no gemmi. No std::string fields.
// - The generator binary (tools/topology/build_semantic_tables) reads
//   chemistry strings from CCD via cifpp and runs RDKit to perceive
//   stereo/aromatic/hybridisation, then emits typed-enum literals
//   into src/generated/LegacyAmberSemanticTables.cpp.
// - The runtime library libnmr_shielding.a does not link RDKit at
//   all. The string barrier is enforced by the linker, not by code
//   review. See spec/plan/topology-substrate-implementation-plan-2026-05-05.md.
//
// Citation pattern:
// - Every typed-enum value has a docstring with a chemistry citation,
//   following the precedent set by Ring::Intensity() with its Case
//   1995 J. Biomol. NMR citation in src/Ring.cpp.
// - Synthesised fields (PlanarGroupKind, PolarHKind, RingPositionLabel)
//   cite the literature that grounds the taxonomy. These are first-
//   class chemistry contributions, not weaker provenance.
//
// Provenance:
// - Each atom-field carries a typed SemanticProvenance witness
//   recording the source(s) and confidence.
// - Witness stash holds up to four source-value pairs so post-hoc
//   stats can ask "for atoms where source A and source B disagreed,
//   what is the shielding distribution?"
//

#include <array>
#include <cstdint>

namespace nmr {


// ============================================================================
// Locant -- Greek-letter / IUPAC position label on the side chain
// ============================================================================
//
// Markley, Bax, Arata, Hilbers, Kaptein, Sykes, Wright, Wuethrich,
// J. Biomol. NMR 12 (1998) 1-23, "Recommendations for the presentation
// of NMR structures of proteins and nucleic acids", §2.1.1.
// Original IUPAC-IUB tentative rules of 1969 (cited as ref. 20 in
// Markley 1998).
//
// The Greek letters (alpha, beta, gamma, ...) walk outward from
// C-alpha along the side-chain heavy-atom chain. Backbone atoms (N,
// CA, C, O, H, HA) carry no locant and use Locant::None. The Roman
// counterparts (A, B, G, D, E, Z, H) are ASCII renderings; storage
// is the integer value, rendering picks Greek (UTF-8) or Roman per
// the output target.
//
enum class Locant : uint8_t {
    None    = 0,   ///< Backbone or terminal atom; no Greek-letter locant.
    Alpha   = 1,   ///< C-alpha (Greek alpha). Glycine HA2/HA3 are alpha-Hs.
    Beta    = 2,   ///< C-beta -- first side-chain heavy atom.
    Gamma   = 3,   ///< C-gamma.
    Delta   = 4,   ///< C-delta.
    Epsilon = 5,   ///< C-epsilon.
    Zeta    = 6,   ///< C-zeta.
    Eta     = 7,   ///< C-eta (Arg side-chain N-eta, Tyr O-eta proton).
};


// ============================================================================
// BranchAddress -- two-level disambiguator when atoms share a locant
// ============================================================================
//
// Markley 1998, Figure 1 caption (text-2:357-372): when atoms
// share a Greek-letter locant, the branch index resolves them per
// the CIP-clockwise rule (sight from highest-priority substituent
// toward the chiral centre; remaining diastereotopic atoms numbered
// 2/3 such that Y, Z2, Z3 follow a clockwise orientation).
//
// Two-level address handles arginine side-chain hydrogens (text-2:
// 367-371): "The hydrogen atoms of the side-chain -NH2 groups of
// Asn, Gln, and Arg are distinguished by numbers (1 or 2) on the
// basis of their relationship (cis or trans, respectively) to the
// heavy atom three bonds closer to the main chain". Arg H-eta11 has
// outer=1 (which N-eta) and inner=1 (which H within that N-eta).
//
// For 99% of atoms only `outer` is non-zero; for Arg side-chain Hs
// both indices apply. Both 0 means the atom needs no branch
// disambiguator.
//
struct BranchAddress {
    uint8_t outer = 0;   ///< Heavy-atom branch (1, 2). 0 = NotApplicable.
    uint8_t inner = 0;   ///< H-within-group branch (1, 2). 0 = NotApplicable.

    constexpr bool IsBranched() const { return outer != 0; }
    constexpr bool HasInner() const { return inner != 0; }
};

constexpr bool operator==(BranchAddress a, BranchAddress b) {
    return a.outer == b.outer && a.inner == b.inner;
}


// ============================================================================
// DiastereotopicIndex -- the IUPAC 2/3 label on prochiral methylene Hs
// ============================================================================
//
// Markley 1998, Figure 1 caption: prochiral methylene hydrogens
// (HB2/HB3, HG2/HG3, HD2/HD3, HE2/HE3 plus glycine HA2/HA3) carry
// numeric labels 2 and 3 by the CIP-derived clockwise rule.
//
// This is the IUPAC numeric label only; the corresponding pro-R /
// pro-S CIP designation lives in ProchiralStereo. Most residues
// alternate (3=R, 2=S) but Gly inverts (HA2=R, HA3=S) -- see Markley
// Fig 1 marks. Per-residue tables required.
//
enum class DiastereotopicIndex : uint8_t {
    None      = 0,   ///< Atom is not part of a prochiral methylene pair.
    Position2 = 2,   ///< IUPAC numeric label "2".
    Position3 = 3,   ///< IUPAC numeric label "3".
};


// ============================================================================
// ProchiralStereo -- Cahn-Ingold-Prelog R/S designation
// ============================================================================
//
// Cahn, Ingold, Prelog, Angew. Chem. Int. Ed. 5 (1966) 385.
// The CIP rules formalise prochiral assignment via priority
// substitution. Markley 1998 §28 (text-7:148-165) restates them:
// for tetrahedral X with substituents A > B > C = C', sight down
// A-X with C' replacing C as the heavier isotope; if B, C', C
// (away from viewer) are clockwise, C' = pro-R, C = pro-S.
//
// The thesis methodology cites "RDKit version 2023.09.6 CIPLabeler"
// as the algorithmic source for this field. RDKit is the primary
// reconciliation source per the project precedence table; Markley
// Figure 1 marks are the cross-check. Disagreements are typed-logged
// in the SemanticProvenance witness stash.
//
// CCD `_chem_comp_atom.pdbx_stereo_config` allowed values are
// {R, S, N} -- this enum maps directly.
// https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_atom.pdbx_stereo_config.html
//
enum class ProchiralStereo : uint8_t {
    NotProchiral = 0,   ///< Atom is not a prochiral centre/substituent.
    ProR         = 1,   ///< CIP rectus (right-handed) per RDKit CIPLabeler.
    ProS         = 2,   ///< CIP sinister (left-handed) per RDKit CIPLabeler.
    Unassigned   = 3,   ///< Prochiral but no CIP determination available;
                        ///< should be empty after the RDKit pass on the
                        ///< standard-20 substrate.
};


// ============================================================================
// PlanarGroupKind -- which chemical planar/sp2 group an atom is in
// ============================================================================
//
// Synthesised classification grounded in standard biochemistry:
// every sp2/planar functional group present in the standard 20
// amino acids gets an enum value. Citation per value below.
//
// The substrate stores the topological membership (which atoms are
// in the group + canonical IUPAC E/Z label per Markley convention).
// Per-frame angles (peptide omega, planarity deviation, ring
// puckering) live on PlanarGeometryResult, the conformation-side
// companion ConformationResult. Substrate-vs-conformation split
// per `feedback_object_model_scope_discipline`.
//
enum class PlanarGroupKind : uint8_t {
    None              = 0,   ///< Atom is not in a planar/sp2 group.

    /// Backbone peptide bond plane. Atoms: C-alpha(i-1), C(i-1),
    /// O(i-1), N(i), H(i), C-alpha(i). Pauling, Corey, Branson,
    /// PNAS 37 (1951) 205-211; Ramachandran & Sasisekharan, Adv.
    /// Protein Chem. 23 (1968) 283-437. SHIFTX2 (Han et al., J.
    /// Biomol. NMR 50 (2011) 43-57) gives omega 5-14% of backbone
    /// shielding signal as a continuous feature.
    PeptideAmide      = 1,

    /// Side-chain primary amide (Asn C-gamma + O-delta1 + N-delta2;
    /// Gln C-delta + O-epsilon1 + N-epsilon2). Resonance-stabilised
    /// planar amide; Asn/Gln side-chain flip is a known crystal
    /// ambiguity per Word, Lovell, Richardson, Richardson, J. Mol.
    /// Biol. 285 (1999) 1735-1747.
    SidechainAmide    = 2,

    /// Arginine guanidinium group. C-zeta + N-epsilon + N-eta1 +
    /// N-eta2 + their hydrogens. Six-electron three-N delocalised
    /// pi system; near three-fold symmetric. Markley 1998 §2.1.4
    /// distinguishes N-eta1/N-eta2 by cis/trans relationship to
    /// C-delta.
    Guanidinium       = 3,

    /// Histidine imidazole ring. Five-membered aromatic with two
    /// nitrogens. Variant-dependent protonation: HID has H on
    /// N-delta1, HIE on N-epsilon2, HIP on both. Hueckel-aromatic
    /// (6 pi electrons including lone pair) per Joule & Mills,
    /// "Heterocyclic Chemistry" 5e (2010) ch. 7.
    Imidazole         = 4,

    /// Aromatic six-membered ring: Phe sidechain, Tyr sidechain,
    /// Trp benzene component (the six-ring of the fused indole).
    /// Standard benzene-like aromaticity. Ring-current intensity
    /// from Christensen, Sauer, Jensen, J. Chem. Theory Comput. 7
    /// (2011) 2078-2084 calibrated against benzene baseline.
    Aromatic6Ring     = 5,

    /// Aromatic five-membered ring: Trp pyrrole component (the
    /// five-ring fused with the benzene), and (variant-dependent)
    /// His imidazole 5-ring. Lower ring-current intensity than
    /// benzene per Christensen et al. (~60% for pyrrole).
    Aromatic5Ring     = 6,

    /// Carboxylate group: Asp C-gamma + O-delta1 + O-delta2;
    /// Glu C-delta + O-epsilon1 + O-epsilon2; C-terminal C + O' + O''.
    /// Two-oxygen-equivalent delocalised; CCD value_order = "delo".
    Carboxylate       = 7,

    /// Tyrosine para-OH. The hydroxyl proton is conventionally
    /// described as in-plane with the ring (Bovey, "Nuclear Magnetic
    /// Resonance Spectroscopy" 2e (1988) ch. 4). C-zeta-O-H rotation
    /// is a real conformation degree of freedom captured by
    /// PlanarGeometryResult, not by this substrate field.
    AromaticHydroxyl  = 8,

    /// Aryloxide. The deprotonated phenolate -O- on TYM (Tyr
    /// deprotonated form). The phenolate conjugates more strongly
    /// with the aromatic ring pi system than neutral phenol-OH does,
    /// redistributing ring current and shifting the para-Czeta
    /// environment by 5-10 ppm relative to neutral Tyr.
    /// Distinguished from `AromaticHydroxyl` because the H is absent
    /// and the conjugation is qualitatively different. Vollhardt &
    /// Schore, "Organic Chemistry" 8e (2018) ch. 22 for phenolate
    /// resonance; Bovey 1988 ch. 4 for shift consequences.
    AromaticOxide     = 9,
};


// ============================================================================
// PlanarStereo -- canonical IUPAC E/Z label
// ============================================================================
//
// CCD `_chem_comp_bond.pdbx_stereo_config` allowed values {E, Z, N}.
// https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_bond.pdbx_stereo_config.html
//
// Encodes the canonical (per-IUPAC-convention) configuration around
// a planar centre. Markley 1998 §2.1.4 (text-2:367-371) gives the
// E/Z conventions for Asn/Gln side-chain primary amides and Arg
// guanidinium hydrogens (cis/trans relative to a heavy atom three
// bonds closer to the main chain).
//
// The actual configuration in a given structure (Asn/Gln may be
// flipped per Word et al. 1999) is captured at load time and may
// disagree with the canonical label. The substrate carries the
// canonical; the conformation side carries the actual.
//
enum class PlanarStereo : uint8_t {
    NotApplicable = 0,   ///< Atom is not at a defined E/Z centre.
    E             = 1,   ///< CIP entgegen (opposite). German convention.
    Z             = 2,   ///< CIP zusammen (together). German convention.
    Unspecified   = 3,   ///< Planar centre exists but configuration unset.
};


// ============================================================================
// PseudoatomKind -- Markley 1998 Table 1 IUPAC pseudoatom letter
// ============================================================================
//
// Markley 1998 Table 1 (text-3:181-242) defines the IUPAC pseudoatom
// taxonomy for proteins. The full set is M, Q, R only:
//
//   M = methyl group (3 equivalent H atoms collapse to one centre)
//   Q = group of equivalent Hs that aren't methyl (typically 2 Hs:
//       methylene; or larger groups -- QH for guanidinium 4Hs)
//   R = ring (used as locant-suffix for ring-atom groupings, e.g.
//       Phe QR = Q across all ring Hs)
//
// XPLOR / CARA / DIANA introduce additional letters (X, Y, S, T)
// but those are NOT IUPAC and are not adopted here. See research
// dossier `spec/plan/topology-fields-research-2026-05-05.md`
// "Field 6" for evidence.
//
// Terminus convention (synthesised, not Markley-canonical):
// N-terminal H1/H2/H3 (NTERM_CHARGED) and H1/H2 (NTERM_NEUTRAL)
// are encoded as `Q` with `Locant::None`. Markley Table 1 does not
// list a pseudoatom for terminus-specific Hs; `Q` is the closest
// fit because they form an equivalent-H group on a single nitrogen
// (rapid exchange with water on the NMR timescale). We do NOT
// introduce a new letter ("T" or otherwise) -- staying within the
// IUPAC M/Q/R taxonomy is a discipline call.
//
enum class PseudoatomKind : uint8_t {
    None = 0,   ///< Atom is not a member of any pseudoatom group.
    M    = 1,   ///< Methyl-group member (3 H pseudoatom).
    Q    = 2,   ///< General equivalent-H-group member (methylene, NH2, NH3).
    R    = 3,   ///< Ring-atom group member (used in the QR aggregator).
};


// ============================================================================
// PseudoatomMembership -- per-atom pseudoatom membership record
// ============================================================================
//
// An atom can belong to a primary pseudoatom (e.g. Leu HD11/HD12/HD13
// belong to MD1, the delta1-methyl) and optionally to a higher-order
// super-group (e.g. all six Leu delta-methyl Hs belong to QD).
//
// Markley 1998 Table 1 (residue table); locant letter is the Greek
// position; branch is 1/2 for diastereotopic methyls (Leu MD1 vs MD2,
// Val MG1 vs MG2). `in_super_group` is true when the atom is also a
// member of the higher-order Q super-group (Val QG over 6 H, Leu QD
// over 6 H, Arg QH over 4 H).
//
struct PseudoatomMembership {
    PseudoatomKind kind            = PseudoatomKind::None;
    uint8_t        locant          = 0;          ///< Maps to the Locant enum integer.
    uint8_t        branch          = 0;          ///< 0/1/2 for branched methyls.
    bool           in_super_group  = false;      ///< True if also in a QG/QD/QH/QR aggregator.

    constexpr bool IsMember() const { return kind != PseudoatomKind::None; }
};


// ============================================================================
// PolarHKind -- functional-group classification of polar hydrogens
// ============================================================================
//
// Synthesised taxonomy grounded in protein physical chemistry.
// Each value corresponds to a well-defined exchangeable-hydrogen
// environment with distinct NMR signature and exchange behaviour.
// Wuethrich, "NMR of Proteins and Nucleic Acids" (1986), ch. 2;
// Englander, Annu. Rev. Biophys. Biomol. Struct. 39 (2008) 289-307
// for the H/D exchange framework.
//
// SHIFTX2 (Han et al. 2011) Table 2 trains separate models for
// HD21/HD22/HE21/HE22/HH11/HH12/HH21/HH22, confirming the chemistry
// distinction matters at the shift level.
//
// Coupled to protonation variant: Cys H-gamma exists only in
// reduced (not CYX); His H-delta1 only in HID/HIP; ASH H-delta2 only
// when protonated. The generator resolves variant interactions per
// (residue, variant) at table-build time.
//
enum class PolarHKind : uint8_t {
    NotPolar             = 0,    ///< C-H or aliphatic H; not exchangeable.

    /// Backbone amide N-H (HN/H). The canonical hydrogen-bond donor
    /// of the polypeptide chain. Chemical shift ~6-10 ppm, dominated
    /// by H-bond geometry and secondary structure. Yi, McDermott,
    /// J. Phys. Chem. B 128 (2024) -- the H-bond angle theta carries
    /// most of the shift sensitivity. Absent in proline (secondary
    /// amine).
    BackboneAmide        = 1,

    /// Side-chain primary amide N-H (Asn HD21, HD22; Gln HE21, HE22).
    /// Resonance-stabilised planar primary amide; shifts cluster
    /// around 6.8-7.8 ppm. Crystal structures may have flipped
    /// O/N assignment per Word, Lovell, Richardson, Richardson, J.
    /// Mol. Biol. 285 (1999) 1735-1747.
    SidechainPrimaryAmide = 2,

    /// Indole N-H (Trp HE1). Aromatic secondary amine; characteristic
    /// shift ~10 ppm, well-separated from the amide region. Slow
    /// exchange under typical solution NMR conditions.
    IndoleNH             = 3,

    /// Ammonium N-H (Lys HZ1/HZ2/HZ3 in protonated state; N-terminal
    /// H1/H2/H3). Quaternary-nitrogen-like protons; rapid exchange
    /// with water at neutral pH; rarely observable in solution NMR.
    AmmoniumNH           = 4,

    /// Guanidinium N-H (Arg HE; HH11, HH12, HH21, HH22). Delocalised
    /// pi-system; intermediate exchange. Markley 1998 §2.1.4 E/Z
    /// distinguishes the four eta-Hs by relationship to C-delta.
    GuanidiniumNH        = 5,

    /// Imidazole ring N-H. Variant-specific: HID has it on
    /// N-delta1, HIE on N-epsilon2, HIP on both. Tautomer state
    /// shifts the chemistry significantly.
    ImidazoleNH          = 6,

    /// Carboxyl O-H. Asp HD2 (in ASH variant only), Glu HE2 (in GLH
    /// variant only), C-terminal H'' (when carboxyl is protonated).
    /// Very acidic; typically not observed in solution NMR at
    /// neutral pH.
    CarboxylOH           = 7,

    /// Aliphatic hydroxyl O-H (Ser HG, Thr HG1). Slow exchange in
    /// folded state; observable under controlled exchange conditions.
    HydroxylOH_Aliphatic = 8,

    /// Aromatic hydroxyl O-H (Tyr HH). The Tyr para-OH; in-plane
    /// with the ring per the AromaticHydroxyl PlanarGroupKind.
    /// Distinguished from the aliphatic OH because the para-pi
    /// delocalisation perturbs the chemistry significantly.
    HydroxylOH_Aromatic  = 9,

    /// Thiol S-H (Cys HG, in reduced CYS only -- absent in CYX
    /// disulfide). Slow exchange; characteristic of reduced
    /// cysteines.
    ThiolSH              = 10,

    /// Neutral primary-amine N-H. Lys neutral form (LYN: HZ2, HZ3)
    /// and the NTERM_NEUTRAL state (H1, H2) on the backbone N.
    /// Distinguished from `AmmoniumNH` because the chemistry differs:
    /// pKa ~9-10 for primary amine vs always-deprotonated water-
    /// exchange for charged ammonium; SHIFTX2 (Han et al., J. Biomol.
    /// NMR 50 (2011) 43-57) trains separate models for these
    /// environments. Wuethrich, "NMR of Proteins and Nucleic Acids"
    /// (1986) ch. 2; Englander, Annu. Rev. Biophys. Biomol. Struct.
    /// 39 (2008) 289-307 for the H/D exchange framework.
    AmineNH              = 11,

    /// Catch-all for non-standard residues (modified amino acids,
    /// caps with novel polar Hs not in the standard 20).
    OtherPolarH          = 12,
};


// ============================================================================
// RingSystemKind -- which ring system an atom belongs to
// ============================================================================
//
// Classifies by the combined (residue, ring-component) since Trp
// has two coexisting rings on its sidechain (pyrrole 5-ring + benzene
// 6-ring fused). Pro is included because its sidechain forms a
// non-aromatic saturated ring that affects shielding and rotamer
// statistics.
//
// Joule & Mills, "Heterocyclic Chemistry" 5e (2010) for heterocycle
// chemistry; Vollhardt & Schore, "Organic Chemistry" 8e (2018) for
// benzene-positioning conventions.
//
enum class RingSystemKind : uint8_t {
    NotInRing       = 0,
    Benzene_Phe     = 1,    ///< Phe sidechain six-membered aromatic.
    Benzene_Tyr     = 2,    ///< Tyr sidechain six-membered aromatic.
    Imidazole_His   = 3,    ///< His imidazole; aromaticity depends on variant.
    Indole_Trp_5    = 4,    ///< Trp pyrrole five-membered ring.
    Indole_Trp_6    = 5,    ///< Trp benzene six-membered ring fused with pyrrole.
    Pyrrolidine_Pro = 6,    ///< Pro saturated five-membered ring (non-aromatic).
};


// ============================================================================
// RingPositionLabel -- position of an atom within its ring system
// ============================================================================
//
// Carries the ipso/ortho/meta/para labels for benzene-like rings
// (Vollhardt & Schore, "Organic Chemistry" 8e (2018) ch. 15) plus
// the pyrrole-position labels (Joule & Mills, "Heterocyclic Chemistry"
// 5e (2010) ch. 13) plus heteroatom-presence subtypes.
//
// For fused-ring systems (Trp), atoms at the bridge belong to both
// rings; they get RingPositionLabel::BridgeFusion in their primary
// ring and a separate RingMembership entry for the secondary ring.
//
// The chemistry-shift statistics motivating granularity: Phe Cdelta
// ~131 ppm (ortho), Cepsilon ~129 ppm (meta), Czeta ~128 ppm (para);
// Tyr Cdelta ~133, Cepsilon ~118, Czeta ~157 (the OH effect at
// para). Trp HE1 (indole NH) ~10 ppm vs amide ~8 ppm. His Hdelta2
// vs Hepsilon1 distinguishable; protonation state shifts them by
// several ppm. These distinctions earn their keep at the labelling
// level.
//
enum class RingPositionLabel : uint8_t {
    NotInRing       = 0,

    /// Connected to the rest of the side chain (typically C-gamma).
    /// Vollhardt & Schore.
    Ipso            = 1,

    /// Adjacent to ipso, branch 1 (per |chi2|-priority rule, Markley
    /// text-2:368-371). Extension for Trp 6-ring perimeter
    /// (synthesised; Markley does not publish ipso/ortho/meta/para
    /// for the indole 6-ring): `Ortho1` = C-epsilon3 (perimeter atom
    /// adjacent to the C-delta2 bridgehead).
    Ortho1          = 2,

    /// Adjacent to ipso, branch 2. For Trp 6-ring (synthesised):
    /// `Ortho2` = C-zeta2 (perimeter atom adjacent to the C-epsilon2
    /// bridgehead).
    Ortho2          = 3,

    /// Two bonds from ipso, branch 1. For Trp 6-ring (synthesised):
    /// `Meta1` = C-zeta3.
    Meta1           = 4,

    /// Two bonds from ipso, branch 2. For Trp 6-ring (synthesised):
    /// `Meta2` = C-eta2.
    Meta2           = 5,

    /// Para to ipso. For Tyr this is C-zeta with the OH; for Phe the
    /// non-substituted aromatic carbon.
    Para            = 6,

    /// Five-ring alpha position (between two heteroatoms or between
    /// a heteroatom and a bridgehead). His C-epsilon1 (between
    /// N-delta1 and N-epsilon2); pyrrole C-2 in standard numbering.
    /// Joule & Mills ch. 13.
    PyrroleAlpha    = 7,

    /// Five-ring beta position (next to a heteroatom but not at
    /// alpha). His C-delta2; Trp C-delta1 (the only non-bridgehead
    /// pyrrole carbon in indole).
    PyrroleBeta     = 8,

    /// Shared atom in a fused-ring system. Trp C-delta2 and
    /// C-epsilon2 belong to both the pyrrole 5-ring and the benzene
    /// 6-ring of the indole.
    BridgeFusion    = 9,

    /// Ring nitrogen carrying an H (Trp N-epsilon1, His N-delta1
    /// in HID/HIP, His N-epsilon2 in HIE/HIP).
    Heteroatom_NH   = 10,

    /// Ring nitrogen without H (His N-epsilon2 in HID, N-delta1 in
    /// HIE).
    Heteroatom_NoH  = 11,

    /// Ring oxygen with H (none in standard 20; reserved for
    /// non-standard residues).
    Heteroatom_OH   = 12,

    /// Saturated ring carbon (Pro C-beta, C-gamma, C-delta).
    Saturated       = 13,
};


// ============================================================================
// RingMembership -- one ring's worth of context for an atom
// ============================================================================
//
// An atom's primary ring is the ring it belongs to (the smaller ring
// if the atom is at a fused-ring bridge). For non-bridge atoms, the
// secondary RingMembership has ring == NotInRing.
//
struct RingMembership {
    RingSystemKind    ring          = RingSystemKind::NotInRing;
    RingPositionLabel position      = RingPositionLabel::NotInRing;
    uint8_t           ring_size     = 0;       ///< 5 or 6, 0 if NotInRing.
    bool              aromatic      = false;   ///< True for aromatic rings; false for Pro.
    bool              planar        = false;   ///< Aromatic + His-protonation-state-dependent.
    uint8_t           n_heteroatoms = 0;       ///< 0 (Phe), 1 (Trp pyrrole), 2 (His).
};


// ============================================================================
// RingPosition -- combined primary + secondary ring memberships
// ============================================================================
//
// For atoms in a single ring, `secondary.ring == NotInRing`. For Trp
// bridgehead atoms (C-delta2, C-epsilon2), both `primary` (the
// 5-ring per the smaller-ring convention) and `secondary` (the
// 6-ring) are populated.
//
// Locant and RingPosition are ORTHOGONAL: `Locant` records the
// backbone-vs-side-chain Greek-letter position; `RingPosition`
// records ring membership. An atom can be in a ring AND have
// `Locant::None`. Pro C-alpha is the canonical example: it is in
// the pyrrolidine ring (`Pyrrolidine_Pro/Saturated/5/f/1`) AND has
// `Locant::None` because Locant is reserved for side-chain atoms.
//
struct RingPosition {
    RingMembership primary;
    RingMembership secondary;

    constexpr bool InAnyRing() const { return primary.ring != RingSystemKind::NotInRing; }
    constexpr bool InTwoRings() const { return secondary.ring != RingSystemKind::NotInRing; }
};


// ============================================================================
// BondOrderToNeighbour -- per-bond order to a specific neighbour
// ============================================================================
//
// CCD `_chem_comp_bond.value_order` allowed values: sing, doub, trip,
// quad, arom, delo, pi, poly. Maps to this enum:
// https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_bond.value_order.html
//
// `Delocalised` is critical for guanidinium (Arg C-zeta to three Ns)
// and carboxylate (Asp/Glu C to two Os), where the formal Lewis
// double bond is symmetrised across multiple bonds.
//
enum class BondOrderToNeighbour : uint8_t {
    None         = 0,   ///< No bond / unset slot.
    Single       = 1,   ///< CCD "sing".
    Double       = 2,   ///< CCD "doub".
    Triple       = 3,   ///< CCD "trip".
    Aromatic     = 4,   ///< CCD "arom" -- classical aromatic bond.
    Delocalised  = 5,   ///< CCD "delo" -- guanidinium / carboxylate.
    Quadruple    = 6,   ///< CCD "quad" -- not present in standard 20.
    Pi           = 7,   ///< CCD "pi" -- reserved.
    Polymer      = 8,   ///< CCD "poly" -- reserved.
};


// ============================================================================
// BondOrderMask -- packed per-neighbour bond orders for an atom
// ============================================================================
//
// Up to 4 neighbours (4-bit each = 16 bits total) suffices for
// every standard amino-acid atom; quaternary carbons have 4 single
// bonds, sp2 atoms have 3 neighbours, etc. Slots are filled in the
// order that CovalentTopology::BondIndicesFor emits.
//
struct BondOrderMask {
    std::array<BondOrderToNeighbour, 4> orders = {
        BondOrderToNeighbour::None,
        BondOrderToNeighbour::None,
        BondOrderToNeighbour::None,
        BondOrderToNeighbour::None,
    };
};


// ============================================================================
// SemanticSource -- typed enum identifying which source produced a value
// ============================================================================
//
// Lists every chemistry-data source contributing to the substrate
// table. The thesis methodology cites algorithmic sources (RDKit,
// cifpp/CCD) primarily because they are version-pinned and
// reproducible.
//
enum class SemanticSource : uint8_t {
    None                      = 0,
    cifpp_CCD_pdbx_aromatic   = 1,    ///< CCD pdbx_aromatic_flag (Y/N).
    cifpp_CCD_pdbx_stereo     = 2,    ///< CCD pdbx_stereo_config (R/S/N or E/Z/N).
    cifpp_CCD_value_order     = 3,    ///< CCD bond value_order.
    cifpp_CCD_leaving_atom    = 4,    ///< CCD pdbx_leaving_atom_flag.
    cifpp_CCD_ordinal         = 5,    ///< CCD pdbx_ordinal.
    RDKit_CIPLabeler          = 6,    ///< RDKit CIP labeler (R/S).
    RDKit_CanonicalRank       = 7,    ///< RDKit Chem::CanonicalRankAtoms.
    RDKit_Hybridisation       = 8,    ///< RDKit Atom::getHybridization.
    RDKit_AromaticPerception  = 9,    ///< RDKit aromatic perception.
    RDKit_BondType            = 10,   ///< RDKit bond perception.
    Markley1998_Table1        = 11,   ///< Markley 1998 Table 1 (pseudoatoms).
    Markley1998_Figure1       = 12,   ///< Markley 1998 Fig 1 R/S marks.
    IUPAC_1969                = 13,   ///< IUPAC-IUB tentative rules 1969.
    BMRB_NomenclaturePage     = 14,   ///< https://bmrb.io/referenc/nomenclature
    BMRB_PseudoatomTable      = 15,   ///< https://bmrb.io/ref_info/pseudoatom_nom.txt
    AmberVariantTable         = 16,   ///< Project AminoAcidType variant table.
    SynthesizedFromChemistry  = 17,   ///< Project synthesis grounded in cited literature.
    AlternationDefault        = 18,   ///< Pre-CIP fallback; should be empty after RDKit pass.
    ManualOverride            = 19,   ///< Hand-edited override; see comment at edit site.
};


// ============================================================================
// SemanticConfidence -- typed enum classifying the confidence of a value
// ============================================================================
//
enum class SemanticConfidence : uint8_t {
    None                    = 0,
    AlgorithmAuthoritative  = 1,   ///< Single algorithmic source, no cross-check needed.
    CIPVerified             = 2,   ///< RDKit CIP and Markley table agree.
    SourcesAgree            = 3,   ///< Multiple sources, all agreed.
    DisagreementLogged      = 4,   ///< Sources disagreed; primary chosen per precedence.
    ManuallyOverridden      = 5,
};


// ============================================================================
// SemanticSourceWitness -- one source's value for a given field
// ============================================================================
//
// Per-(atom, field) provenance can hold up to four witnesses. Empty
// slots have source == None. The value is packed as the underlying
// enum width (uint8_t for all our enums) so this struct is fixed-
// size at 2 bytes.
//
struct SemanticSourceWitness {
    SemanticSource source       = SemanticSource::None;
    uint8_t        value_packed = 0;
};


// ============================================================================
// SemanticProvenance -- per-(atom, field) provenance record
// ============================================================================
//
// Stores up to four source witnesses, the chosen-primary index, the
// confidence classification, and a CIP-verification flag. Strings
// only appear when this is rendered for human consumption at the
// H5 emission boundary; runtime carries enums.
//
struct SemanticProvenance {
    std::array<SemanticSourceWitness, 4> witnesses = {};
    uint8_t            primary_witness_idx = 0;       ///< Index into witnesses (0..3).
    SemanticConfidence confidence          = SemanticConfidence::None;
    bool               sources_agree       = true;    ///< False if any disagreement.
    bool               cip_verified        = false;   ///< RDKit CIP confirmed value.
};


// ============================================================================
// AtomSemanticTable -- the runtime per-atom record consumed by the
// LegacyAmberTopology populator
// ============================================================================
//
// One record per (residue, variant, atom_local_idx). The generated
// table at src/generated/LegacyAmberSemanticTables.cpp emits a
// constexpr std::array of these per residue; the populator looks up
// entries at construction time and writes into the typed fields on
// LegacyAmberTopology.
//
// Provenance is NOT carried in the runtime record (it lives in the
// generation log, src/generated/LegacyAmberSemanticTables.log.txt,
// and is committed alongside as the audit trail). Runtime code never
// needs to inspect provenance; downstream stats analysis joins via
// (residue, atom_id) keys against the log.
//
struct AtomSemanticTable {
    Locant                locant            = Locant::None;
    BranchAddress         branch            = {};
    DiastereotopicIndex   di_index          = DiastereotopicIndex::None;
    ProchiralStereo       prochiral         = ProchiralStereo::NotProchiral;
    PlanarGroupKind       planar_group      = PlanarGroupKind::None;
    PlanarStereo          planar_stereo     = PlanarStereo::NotApplicable;
    PseudoatomMembership  pseudoatom        = {};
    PolarHKind            polar_h           = PolarHKind::NotPolar;
    RingPosition          ring_position     = {};
    bool                  aromatic          = false;
    int8_t                formal_charge     = 0;
    bool                  is_exchangeable   = false;
    uint8_t               equivalence_class = 0;
    // Note: Hybridisation (from Types.h) and BondOrderMask are not
    // included in the runtime record at the moment; they are
    // available via the generation log if needed for analysis.
    // The runtime is intentionally minimal until a calculator
    // requires a specific field on the substrate.
};


}  // namespace nmr
