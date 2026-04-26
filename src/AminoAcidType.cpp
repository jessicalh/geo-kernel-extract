// AminoAcidType: the complete amino acid chemistry table.
// Single source of truth for all amino acid properties.

#include "AminoAcidType.h"

namespace nmr {

using A = AminoAcidAtom;
using E = Element;

// Standard backbone atoms (most amino acids)
//
// Backbone nitrogen N, peptide H, carbonyl C, carbonyl O carry Locant::Backbone.
// Cα and Hα carry Locant::Alpha (start of the IUPAC Greek-letter chain).
// Glycine substitutes prochiral Hα2/Hα3 for the single HA (pro-R/pro-S per
// Markley 1998 Figure 1; Hα2 is pro-R for Gly).
//
// Terminal atoms (H1/H2/H3 on protonated N-term, OXT on C-term carboxylate)
// are listed in every BB macro because reduce / pdb2gmx can produce them on
// the terminal residue regardless of amino-acid type. Atoms not present in a
// given residue are simply not in its atom_indices and not stamped.
//
// Backbone topology is identical across all 20 amino acids, so the macros
// pre-populate it for every residue automatically.
#define BB  A{"N",E::N,true,TopologyBackbone()}, \
            A{"CA",E::C,true,TopologyAlpha()}, \
            A{"C",E::C,true,TopologyBackbone()}, \
            A{"O",E::O,true,TopologyBackbone()}, \
            A{"H",E::H,true,TopologyBackbone()}, \
            A{"HA",E::H,true,TopologyAlpha()}, \
            A{"H1",E::H,true,TopologyTerminalH()}, \
            A{"H2",E::H,true,TopologyTerminalH()}, \
            A{"H3",E::H,true,TopologyTerminalH()}, \
            A{"OXT",E::O,true,TopologyTerminalO()}
#define BB_PRO  A{"N",E::N,true,TopologyBackbone()}, \
                A{"CA",E::C,true,TopologyAlpha()}, \
                A{"C",E::C,true,TopologyBackbone()}, \
                A{"O",E::O,true,TopologyBackbone()}, \
                A{"HA",E::H,true,TopologyAlpha()}, \
                A{"H2",E::H,true,TopologyTerminalH()}, \
                A{"H3",E::H,true,TopologyTerminalH()}, \
                A{"OXT",E::O,true,TopologyTerminalO()}
#define BB_GLY  A{"N",E::N,true,TopologyBackbone()}, \
                A{"CA",E::C,true,TopologyAlpha()}, \
                A{"C",E::C,true,TopologyBackbone()}, \
                A{"O",E::O,true,TopologyBackbone()}, \
                A{"H",E::H,true,TopologyBackbone()}, \
                A{"HA2",E::H,true,TopologyAlphaProchiral(DiastereotopicIndex::Two, ProchiralStereo::ProR)}, \
                A{"HA3",E::H,true,TopologyAlphaProchiral(DiastereotopicIndex::Three, ProchiralStereo::ProS)}, \
                A{"H1",E::H,true,TopologyTerminalH()}, \
                A{"H2",E::H,true,TopologyTerminalH()}, \
                A{"H3",E::H,true,TopologyTerminalH()}, \
                A{"OXT",E::O,true,TopologyTerminalO()}

static const std::vector<AminoAcidType> AMINO_ACID_TYPES = {

// ALA — β-methyl. HB1/HB2/HB3 are equivalent under fast rotation (pseudoatom MB).
{ AminoAcid::ALA, "ALA", 'A', false, false, true, 0, 0,
  {BB,
   {"CB", E::C, false, TopologySidechain(Locant::Beta)},
   {"HB1",E::H, false, TopologyMethyl(Locant::Beta, BranchIndex::None, PseudoatomClass::MB)},
   {"HB2",E::H, false, TopologyMethyl(Locant::Beta, BranchIndex::None, PseudoatomClass::MB)},
   {"HB3",E::H, false, TopologyMethyl(Locant::Beta, BranchIndex::None, PseudoatomClass::MB)}},
  {}, {}, {} },

// ARG — Markley 1998 Sec 2.1.1: Nη1 is cis (Z) to Cδ, Nη2 is trans (E).
// HHx1/HHx2 distinguished cis/trans relative to Nε (heavy atom three bonds
// closer to main chain than the H). Pseudoatoms QH1 (Hη11/Hη12), QH2 (Hη21/Hη22).
{ AminoAcid::ARG, "ARG", 'R', false, true, true, 4, +1,
  {BB,
   {"CB",  E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2", E::H, false, TopologyDiastereo(Locant::Beta,    DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3", E::H, false, TopologyDiastereo(Locant::Beta,    DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"CG",  E::C, false, TopologySidechain(Locant::Gamma)},
   {"HG2", E::H, false, TopologyDiastereo(Locant::Gamma,   DiastereotopicIndex::Two,   ProchiralStereo::ProR, PseudoatomClass::QG)},
   {"HG3", E::H, false, TopologyDiastereo(Locant::Gamma,   DiastereotopicIndex::Three, ProchiralStereo::ProS, PseudoatomClass::QG)},
   {"CD",  E::C, false, TopologySidechain(Locant::Delta)},
   {"HD2", E::H, false, TopologyDiastereo(Locant::Delta,   DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QD)},
   {"HD3", E::H, false, TopologyDiastereo(Locant::Delta,   DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QD)},
   {"NE",  E::N, false, TopologySidechain(Locant::Epsilon)},
   {"HE",  E::H, false, TopologySidechain(Locant::Epsilon)},
   {"CZ",  E::C, false, TopologySidechain(Locant::Zeta)},
   {"NH1", E::N, false, TopologyPlanar(Locant::Eta, BranchIndex::One, PlanarStereo::Cis)},
   {"HH11",E::H, false, TopologyPlanar(Locant::Eta, BranchIndex::One, PlanarStereo::Cis,   PseudoatomClass::QH1)},
   {"HH12",E::H, false, TopologyPlanar(Locant::Eta, BranchIndex::One, PlanarStereo::Trans, PseudoatomClass::QH1)},
   {"NH2", E::N, false, TopologyPlanar(Locant::Eta, BranchIndex::Two, PlanarStereo::Trans)},
   {"HH21",E::H, false, TopologyPlanar(Locant::Eta, BranchIndex::Two, PlanarStereo::Cis,   PseudoatomClass::QH2)},
   {"HH22",E::H, false, TopologyPlanar(Locant::Eta, BranchIndex::Two, PlanarStereo::Trans, PseudoatomClass::QH2)}},
  {},
  {{"N","CA","CB","CG"}, {"CA","CB","CG","CD"}, {"CB","CG","CD","NE"}, {"CG","CD","NE","CZ"}},
  {{"ARN", "deprotonated arginine", 0, "deprotonated"}} },

// ASN — δ2-amido NH2 (Hδ21 cis to Cβ, Hδ22 trans). Pseudoatom QD = HD21+HD22.
{ AminoAcid::ASN, "ASN", 'N', false, false, true, 2, 0,
  {BB,
   {"CB",  E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2", E::H, false, TopologyDiastereo(Locant::Beta,  DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3", E::H, false, TopologyDiastereo(Locant::Beta,  DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"CG",  E::C, false, TopologySidechain(Locant::Gamma)},
   {"OD1", E::O, false, TopologySidechain(Locant::Delta, BranchIndex::One)},
   {"ND2", E::N, false, TopologySidechain(Locant::Delta, BranchIndex::Two)},
   {"HD21",E::H, false, TopologyPlanar(Locant::Delta, BranchIndex::Two, PlanarStereo::Cis,   PseudoatomClass::QD)},
   {"HD22",E::H, false, TopologyPlanar(Locant::Delta, BranchIndex::Two, PlanarStereo::Trans, PseudoatomClass::QD)}},
  {}, {{"N","CA","CB","CG"}, {"CA","CB","CG","OD1"}}, {} },

// ASP — Cδ1/Cδ2 distinguished by smaller |χ2| (carboxylate Os). HD2 atom
// included for the ASH protonated variant; absent in default deprotonated
// state and not stamped if missing. When present (ASH), HD2 is the
// carboxyl-acid O-H proton — chemistry-tagged via PolarHKind::CarboxylAcid
// so post-hoc analysis distinguishes it from a generic sidechain H.
{ AminoAcid::ASP, "ASP", 'D', false, true, true, 2, -1,
  {BB,
   {"CB",  E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2", E::H, false, TopologyDiastereo(Locant::Beta, DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3", E::H, false, TopologyDiastereo(Locant::Beta, DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"CG",  E::C, false, TopologySidechain(Locant::Gamma)},
   {"OD1", E::O, false, TopologySidechain(Locant::Delta, BranchIndex::One)},
   {"OD2", E::O, false, TopologySidechain(Locant::Delta, BranchIndex::Two)},
   {"HD2", E::H, false, TopologyPolarH(Locant::Delta, BranchIndex::Two, PolarHKind::CarboxylAcid)}},  // ASH variant only
  {}, {{"N","CA","CB","CG"}, {"CA","CB","CG","OD1"}},
  {{"ASH", "protonated aspartate", 0, "protonated"}} },

// CYS — γ-thiol HG is a single H (not prochiral). CYX (disulfide) and CYM
// (deprotonated thiolate) variants modify the SG/HG state. HG when present
// is the thiol S-H — PolarHKind::Sulfanyl marks it as labile.
{ AminoAcid::CYS, "CYS", 'C', false, true, true, 1, -1,
  {BB,
   {"CB",  E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2", E::H, false, TopologyDiastereo(Locant::Beta, DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3", E::H, false, TopologyDiastereo(Locant::Beta, DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"SG",  E::S, false, TopologySidechain(Locant::Gamma)},
   {"HG",  E::H, false, TopologyPolarH(Locant::Gamma, BranchIndex::None, PolarHKind::Sulfanyl)}},
  {}, {{"N","CA","CB","SG"}},
  {{"CYX", "disulfide bonded", 0, "disulfide"}, {"CYM", "deprotonated thiolate", -1, "deprotonated"}} },

// GLN — ε2-amido NH2 (Hε21 cis to Cγ, Hε22 trans). Pseudoatom QE = HE21+HE22.
{ AminoAcid::GLN, "GLN", 'Q', false, false, true, 3, 0,
  {BB,
   {"CB",  E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2", E::H, false, TopologyDiastereo(Locant::Beta,    DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3", E::H, false, TopologyDiastereo(Locant::Beta,    DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"CG",  E::C, false, TopologySidechain(Locant::Gamma)},
   {"HG2", E::H, false, TopologyDiastereo(Locant::Gamma,   DiastereotopicIndex::Two,   ProchiralStereo::ProR, PseudoatomClass::QG)},
   {"HG3", E::H, false, TopologyDiastereo(Locant::Gamma,   DiastereotopicIndex::Three, ProchiralStereo::ProS, PseudoatomClass::QG)},
   {"CD",  E::C, false, TopologySidechain(Locant::Delta)},
   {"OE1", E::O, false, TopologySidechain(Locant::Epsilon, BranchIndex::One)},
   {"NE2", E::N, false, TopologySidechain(Locant::Epsilon, BranchIndex::Two)},
   {"HE21",E::H, false, TopologyPlanar(Locant::Epsilon, BranchIndex::Two, PlanarStereo::Cis,   PseudoatomClass::QE)},
   {"HE22",E::H, false, TopologyPlanar(Locant::Epsilon, BranchIndex::Two, PlanarStereo::Trans, PseudoatomClass::QE)}},
  {}, {{"N","CA","CB","CG"}, {"CA","CB","CG","CD"}, {"CB","CG","CD","OE1"}}, {} },

// GLU — Cε1/Cε2 distinguished by smaller |χ2| (carboxylate Os). HE2 atom
// included for the GLH protonated variant; absent in default deprotonated
// state and not stamped if missing. When present (GLH), HE2 is the
// carboxyl-acid O-H proton — chemistry-tagged via PolarHKind::CarboxylAcid.
{ AminoAcid::GLU, "GLU", 'E', false, true, true, 3, -1,
  {BB,
   {"CB",  E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2", E::H, false, TopologyDiastereo(Locant::Beta,    DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3", E::H, false, TopologyDiastereo(Locant::Beta,    DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"CG",  E::C, false, TopologySidechain(Locant::Gamma)},
   {"HG2", E::H, false, TopologyDiastereo(Locant::Gamma,   DiastereotopicIndex::Two,   ProchiralStereo::ProR, PseudoatomClass::QG)},
   {"HG3", E::H, false, TopologyDiastereo(Locant::Gamma,   DiastereotopicIndex::Three, ProchiralStereo::ProS, PseudoatomClass::QG)},
   {"CD",  E::C, false, TopologySidechain(Locant::Delta)},
   {"OE1", E::O, false, TopologySidechain(Locant::Epsilon, BranchIndex::One)},
   {"OE2", E::O, false, TopologySidechain(Locant::Epsilon, BranchIndex::Two)},
   {"HE2", E::H, false, TopologyPolarH(Locant::Epsilon, BranchIndex::Two, PolarHKind::CarboxylAcid)}},  // GLH variant only
  {}, {{"N","CA","CB","CG"}, {"CA","CB","CG","CD"}, {"CB","CG","CD","OE1"}},
  {{"GLH", "protonated glutamate", 0, "protonated"}} },

// GLY
{ AminoAcid::GLY, "GLY", 'G', false, false, true, 0, 0,
  {BB_GLY}, {}, {}, {} },

// HIS — imidazole ring (5 vertices). HD1 present in HID/HIP; HE2 present in
// HIE/HIP. Both atoms listed in template; absence per variant handled at load.
// The default base entry covers the HisImidazole ring; HID/HIE/HIP variants
// modify the ring TYPE (HidImidazole/HieImidazole) but not the atom names.
//
// HD1 (when present in HID/HIP) and HE2 (when present in HIE/HIP) are
// imidazole ring N-H protons — chemistry-tagged via PolarHKind::RingNH so
// post-hoc analysis distinguishes them from generic ring C-H or sidechain H.
// HD2 (always present, on Cδ2) and HE1 (always present, on Cε1) are
// aromatic C-H — these stay as generic sidechain H.
{ AminoAcid::HIS, "HIS", 'H', true, true, true, 2, +1,
  {BB,
   {"CB", E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2",E::H, false, TopologyDiastereo(Locant::Beta, DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3",E::H, false, TopologyDiastereo(Locant::Beta, DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"CG", E::C, false, TopologyRing(Locant::Gamma,   BranchIndex::None, RingPosition::Substituent)},
   {"ND1",E::N, false, TopologyRing(Locant::Delta,   BranchIndex::One,  RingPosition::Member)},
   {"HD1",E::H, false, TopologyPolarH(Locant::Delta, BranchIndex::One, PolarHKind::RingNH)},   // HID/HIP only
   {"CD2",E::C, false, TopologyRing(Locant::Delta,   BranchIndex::Two,  RingPosition::Member)},
   {"HD2",E::H, false, TopologySidechain(Locant::Delta,   BranchIndex::Two)},
   {"CE1",E::C, false, TopologyRing(Locant::Epsilon, BranchIndex::One,  RingPosition::Member)},
   {"HE1",E::H, false, TopologySidechain(Locant::Epsilon, BranchIndex::One)},
   {"NE2",E::N, false, TopologyRing(Locant::Epsilon, BranchIndex::Two,  RingPosition::Member)},
   {"HE2",E::H, false, TopologyPolarH(Locant::Epsilon, BranchIndex::Two, PolarHKind::RingNH)}},  // HIE/HIP only
  {{RingTypeIndex::HisImidazole, {"CG","ND1","CE1","NE2","CD2"}}},
  {{"N","CA","CB","CG"}, {"CA","CB","CG","ND1"}},
  {{"HID", "Nd-protonated (delta)", 0, "delta"},
   {"HIE", "Ne-protonated (epsilon)", 0, "epsilon"},
   {"HIP", "doubly protonated", +1, "doubly"}} },

// ILE — Cβ is real chiral (no Hβ2/Hβ3, single Hβ). CG1 is a methylene branch
// (prochiral Hγ12/Hγ13). CG2 and CD1 are methyl carbons (pseudoatoms MG, MD).
{ AminoAcid::ILE, "ILE", 'I', false, false, true, 2, 0,
  {BB,
   {"CB",  E::C, false, TopologySidechain(Locant::Beta)},
   {"HB",  E::H, false, TopologySidechain(Locant::Beta)},
   {"CG1", E::C, false, TopologySidechain(Locant::Gamma, BranchIndex::One)},
   {"HG12",E::H, false, TopologyDiastereo(Locant::Gamma, DiastereotopicIndex::Two,   ProchiralStereo::ProR, PseudoatomClass::QG, BranchIndex::One)},
   {"HG13",E::H, false, TopologyDiastereo(Locant::Gamma, DiastereotopicIndex::Three, ProchiralStereo::ProS, PseudoatomClass::QG, BranchIndex::One)},
   {"CG2", E::C, false, TopologySidechain(Locant::Gamma, BranchIndex::Two)},
   {"HG21",E::H, false, TopologyMethyl(Locant::Gamma, BranchIndex::Two, PseudoatomClass::MG)},
   {"HG22",E::H, false, TopologyMethyl(Locant::Gamma, BranchIndex::Two, PseudoatomClass::MG)},
   {"HG23",E::H, false, TopologyMethyl(Locant::Gamma, BranchIndex::Two, PseudoatomClass::MG)},
   {"CD1", E::C, false, TopologySidechain(Locant::Delta, BranchIndex::One)},
   {"HD11",E::H, false, TopologyMethyl(Locant::Delta, BranchIndex::One, PseudoatomClass::MD)},
   {"HD12",E::H, false, TopologyMethyl(Locant::Delta, BranchIndex::One, PseudoatomClass::MD)},
   {"HD13",E::H, false, TopologyMethyl(Locant::Delta, BranchIndex::One, PseudoatomClass::MD)}},
  {}, {{"N","CA","CB","CG1"}, {"CA","CB","CG1","CD1"}}, {} },

// LEU — Cγ is real chiral (single Hγ). Cδ1 and Cδ2 are diastereotopic methyl
// carbons; pseudoatom QD is the union of all 6 δ-methyl Hs.
{ AminoAcid::LEU, "LEU", 'L', false, false, true, 2, 0,
  {BB,
   {"CB",  E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2", E::H, false, TopologyDiastereo(Locant::Beta,  DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3", E::H, false, TopologyDiastereo(Locant::Beta,  DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"CG",  E::C, false, TopologySidechain(Locant::Gamma)},
   {"HG",  E::H, false, TopologySidechain(Locant::Gamma)},
   {"CD1", E::C, false, TopologySidechain(Locant::Delta, BranchIndex::One)},
   {"HD11",E::H, false, TopologyMethyl(Locant::Delta, BranchIndex::One, PseudoatomClass::MD1)},
   {"HD12",E::H, false, TopologyMethyl(Locant::Delta, BranchIndex::One, PseudoatomClass::MD1)},
   {"HD13",E::H, false, TopologyMethyl(Locant::Delta, BranchIndex::One, PseudoatomClass::MD1)},
   {"CD2", E::C, false, TopologySidechain(Locant::Delta, BranchIndex::Two)},
   {"HD21",E::H, false, TopologyMethyl(Locant::Delta, BranchIndex::Two, PseudoatomClass::MD2)},
   {"HD22",E::H, false, TopologyMethyl(Locant::Delta, BranchIndex::Two, PseudoatomClass::MD2)},
   {"HD23",E::H, false, TopologyMethyl(Locant::Delta, BranchIndex::Two, PseudoatomClass::MD2)}},
  {}, {{"N","CA","CB","CG"}, {"CA","CB","CG","CD1"}}, {} },

// LYS — four prochiral methylenes (β/γ/δ/ε). Nζ-NH3+ in protonated state has
// equivalent Hζ1/Hζ2/Hζ3 (pseudoatom QZ); LYN variant drops Hζ3.
{ AminoAcid::LYS, "LYS", 'K', false, true, true, 4, +1,
  {BB,
   {"CB", E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2",E::H, false, TopologyDiastereo(Locant::Beta,    DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3",E::H, false, TopologyDiastereo(Locant::Beta,    DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"CG", E::C, false, TopologySidechain(Locant::Gamma)},
   {"HG2",E::H, false, TopologyDiastereo(Locant::Gamma,   DiastereotopicIndex::Two,   ProchiralStereo::ProR, PseudoatomClass::QG)},
   {"HG3",E::H, false, TopologyDiastereo(Locant::Gamma,   DiastereotopicIndex::Three, ProchiralStereo::ProS, PseudoatomClass::QG)},
   {"CD", E::C, false, TopologySidechain(Locant::Delta)},
   {"HD2",E::H, false, TopologyDiastereo(Locant::Delta,   DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QD)},
   {"HD3",E::H, false, TopologyDiastereo(Locant::Delta,   DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QD)},
   {"CE", E::C, false, TopologySidechain(Locant::Epsilon)},
   {"HE2",E::H, false, TopologyDiastereo(Locant::Epsilon, DiastereotopicIndex::Two,   ProchiralStereo::ProR, PseudoatomClass::QE)},
   {"HE3",E::H, false, TopologyDiastereo(Locant::Epsilon, DiastereotopicIndex::Three, ProchiralStereo::ProS, PseudoatomClass::QE)},
   {"NZ", E::N, false, TopologySidechain(Locant::Zeta)},
   {"HZ1",E::H, false, TopologyMethyl(Locant::Zeta, BranchIndex::None, PseudoatomClass::QZ)},
   {"HZ2",E::H, false, TopologyMethyl(Locant::Zeta, BranchIndex::None, PseudoatomClass::QZ)},
   {"HZ3",E::H, false, TopologyMethyl(Locant::Zeta, BranchIndex::None, PseudoatomClass::QZ)}},  // LYN drops HZ3
  {},
  {{"N","CA","CB","CG"}, {"CA","CB","CG","CD"}, {"CB","CG","CD","CE"}, {"CG","CD","CE","NZ"}},
  {{"LYN", "deprotonated lysine", 0, "deprotonated"}} },

// MET — Sδ thioether sulfur. Cε is the methyl carbon; Hε1/Hε2/Hε3 are the
// equivalent thiomethyl protons (pseudoatom ME).
{ AminoAcid::MET, "MET", 'M', false, false, true, 3, 0,
  {BB,
   {"CB", E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2",E::H, false, TopologyDiastereo(Locant::Beta,    DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3",E::H, false, TopologyDiastereo(Locant::Beta,    DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"CG", E::C, false, TopologySidechain(Locant::Gamma)},
   {"HG2",E::H, false, TopologyDiastereo(Locant::Gamma,   DiastereotopicIndex::Two,   ProchiralStereo::ProR, PseudoatomClass::QG)},
   {"HG3",E::H, false, TopologyDiastereo(Locant::Gamma,   DiastereotopicIndex::Three, ProchiralStereo::ProS, PseudoatomClass::QG)},
   {"SD", E::S, false, TopologySidechain(Locant::Delta)},
   {"CE", E::C, false, TopologySidechain(Locant::Epsilon)},
   {"HE1",E::H, false, TopologyMethyl(Locant::Epsilon, BranchIndex::None, PseudoatomClass::ME)},
   {"HE2",E::H, false, TopologyMethyl(Locant::Epsilon, BranchIndex::None, PseudoatomClass::ME)},
   {"HE3",E::H, false, TopologyMethyl(Locant::Epsilon, BranchIndex::None, PseudoatomClass::ME)}},
  {}, {{"N","CA","CB","CG"}, {"CA","CB","CG","SD"}, {"CB","CG","SD","CE"}}, {} },

// PHE — benzene ring. HB2/HB3 prochiral β-methylene (pro-S/pro-R). Ring protons
// HD1/HD2 form pseudoatom QD; HE1/HE2 form QE; HZ is unique. Cδ1/Cδ2 (and
// Cε1/Cε2) distinguished by smaller |χ2| convention (Markley 1998 Fig 1 caption).
{ AminoAcid::PHE, "PHE", 'F', true, false, true, 2, 0,
  {BB,
   {"CB", E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2",E::H, false, TopologyDiastereo(Locant::Beta, DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3",E::H, false, TopologyDiastereo(Locant::Beta, DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"CG", E::C, false, TopologyRing(Locant::Gamma,   BranchIndex::None, RingPosition::Substituent)},
   {"CD1",E::C, false, TopologyRing(Locant::Delta,   BranchIndex::One,  RingPosition::Member)},
   {"HD1",E::H, false, TopologySidechain(Locant::Delta,   BranchIndex::One)},
   {"CD2",E::C, false, TopologyRing(Locant::Delta,   BranchIndex::Two,  RingPosition::Member)},
   {"HD2",E::H, false, TopologySidechain(Locant::Delta,   BranchIndex::Two)},
   {"CE1",E::C, false, TopologyRing(Locant::Epsilon, BranchIndex::One,  RingPosition::Member)},
   {"HE1",E::H, false, TopologySidechain(Locant::Epsilon, BranchIndex::One)},
   {"CE2",E::C, false, TopologyRing(Locant::Epsilon, BranchIndex::Two,  RingPosition::Member)},
   {"HE2",E::H, false, TopologySidechain(Locant::Epsilon, BranchIndex::Two)},
   {"CZ", E::C, false, TopologyRing(Locant::Zeta,    BranchIndex::None, RingPosition::Member)},
   {"HZ", E::H, false, TopologySidechain(Locant::Zeta,    BranchIndex::None)}},
  {{RingTypeIndex::PheBenzene, {"CG","CD1","CE1","CZ","CE2","CD2"}}},
  {{"N","CA","CB","CG"}, {"CA","CB","CG","CD1"}}, {} },

// PRO — pyrrolidine ring (N–Cα–Cβ–Cγ–Cδ–N). No backbone H (BB_PRO). Three
// prochiral methylenes (β/γ/δ); pucker conformations defined by χ1..χ4.
{ AminoAcid::PRO, "PRO", 'P', false, false, false, 2, 0,
  {BB_PRO,
   {"CB", E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2",E::H, false, TopologyDiastereo(Locant::Beta,  DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3",E::H, false, TopologyDiastereo(Locant::Beta,  DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"CG", E::C, false, TopologySidechain(Locant::Gamma)},
   {"HG2",E::H, false, TopologyDiastereo(Locant::Gamma, DiastereotopicIndex::Two,   ProchiralStereo::ProR, PseudoatomClass::QG)},
   {"HG3",E::H, false, TopologyDiastereo(Locant::Gamma, DiastereotopicIndex::Three, ProchiralStereo::ProS, PseudoatomClass::QG)},
   {"CD", E::C, false, TopologySidechain(Locant::Delta)},
   {"HD2",E::H, false, TopologyDiastereo(Locant::Delta, DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QD)},
   {"HD3",E::H, false, TopologyDiastereo(Locant::Delta, DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QD)}},
  {}, {{"N","CA","CB","CG"}, {"CA","CB","CG","CD"}}, {} },

// SER — γ-hydroxyl. Hγ is a single hydroxyl O-H (not prochiral).
{ AminoAcid::SER, "SER", 'S', false, false, true, 1, 0,
  {BB,
   {"CB", E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2",E::H, false, TopologyDiastereo(Locant::Beta, DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3",E::H, false, TopologyDiastereo(Locant::Beta, DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"OG", E::O, false, TopologySidechain(Locant::Gamma)},
   {"HG", E::H, false, TopologyPolarH(Locant::Gamma, BranchIndex::None, PolarHKind::Hydroxyl)}},
  {}, {{"N","CA","CB","OG"}}, {} },

// THR — Cβ is real chiral (single Hβ). Oγ1 = hydroxyl O; Cγ2 = methyl C.
// Hγ1 is the hydroxyl O-H — chemistry-tagged via PolarHKind::Hydroxyl.
{ AminoAcid::THR, "THR", 'T', false, false, true, 1, 0,
  {BB,
   {"CB",  E::C, false, TopologySidechain(Locant::Beta)},
   {"HB",  E::H, false, TopologySidechain(Locant::Beta)},
   {"OG1", E::O, false, TopologySidechain(Locant::Gamma, BranchIndex::One)},
   {"HG1", E::H, false, TopologyPolarH(Locant::Gamma, BranchIndex::One, PolarHKind::Hydroxyl)},
   {"CG2", E::C, false, TopologySidechain(Locant::Gamma, BranchIndex::Two)},
   {"HG21",E::H, false, TopologyMethyl(Locant::Gamma, BranchIndex::Two, PseudoatomClass::MG)},
   {"HG22",E::H, false, TopologyMethyl(Locant::Gamma, BranchIndex::Two, PseudoatomClass::MG)},
   {"HG23",E::H, false, TopologyMethyl(Locant::Gamma, BranchIndex::Two, PseudoatomClass::MG)}},
  {}, {{"N","CA","CB","OG1"}}, {} },

// TRP — indole: pyrrole (5-ring) + benzene (6-ring) fused, with the 9-atom
// perimeter as the third (TrpPerimeter) ring. Cδ2 / Cε2 are ring junctions
// (members of all three rings). Nε1 is the indole NH; Hε1 carries
// PolarHKind::RingNH so post-hoc analysis sees the indole-proton chemistry.
{ AminoAcid::TRP, "TRP", 'W', true, false, true, 2, 0,
  {BB,
   {"CB", E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2",E::H, false, TopologyDiastereo(Locant::Beta, DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3",E::H, false, TopologyDiastereo(Locant::Beta, DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"CG", E::C, false, TopologyRing(Locant::Gamma,   BranchIndex::None,  RingPosition::Substituent)},
   {"CD1",E::C, false, TopologyRing(Locant::Delta,   BranchIndex::One,   RingPosition::Member)},
   {"HD1",E::H, false, TopologySidechain(Locant::Delta, BranchIndex::One)},
   {"CD2",E::C, false, TopologyRing(Locant::Delta,   BranchIndex::Two,   RingPosition::Junction)},
   {"NE1",E::N, false, TopologyRing(Locant::Epsilon, BranchIndex::One,   RingPosition::Member)},
   {"HE1",E::H, false, TopologyPolarH(Locant::Epsilon, BranchIndex::One, PolarHKind::RingNH)},  // indole NH
   {"CE2",E::C, false, TopologyRing(Locant::Epsilon, BranchIndex::Two,   RingPosition::Junction)},
   {"CE3",E::C, false, TopologyRing(Locant::Epsilon, BranchIndex::Three, RingPosition::Member)},
   {"HE3",E::H, false, TopologySidechain(Locant::Epsilon, BranchIndex::Three)},
   {"CZ2",E::C, false, TopologyRing(Locant::Zeta,    BranchIndex::Two,   RingPosition::Member)},
   {"HZ2",E::H, false, TopologySidechain(Locant::Zeta,    BranchIndex::Two)},
   {"CZ3",E::C, false, TopologyRing(Locant::Zeta,    BranchIndex::Three, RingPosition::Member)},
   {"HZ3",E::H, false, TopologySidechain(Locant::Zeta,    BranchIndex::Three)},
   {"CH2",E::C, false, TopologyRing(Locant::Eta,     BranchIndex::Two,   RingPosition::Member)},
   {"HH2",E::H, false, TopologySidechain(Locant::Eta,     BranchIndex::Two)}},
  {{RingTypeIndex::TrpBenzene,   {"CD2","CE2","CZ2","CH2","CZ3","CE3"}},
   {RingTypeIndex::TrpPyrrole,   {"CG","CD1","NE1","CE2","CD2"}},
   {RingTypeIndex::TrpPerimeter, {"CG","CD1","NE1","CE2","CZ2","CH2","CZ3","CE3","CD2"}}},
  {{"N","CA","CB","CG"}, {"CA","CB","CG","CD1"}}, {} },

// TYR — phenol ring with hydroxyl at Cζ. Same Cδ1/Cδ2 + Cε1/Cε2 ordering as
// PHE (smaller |χ2| → 1). HH absent in TYM (deprotonated tyrosinate) variant.
{ AminoAcid::TYR, "TYR", 'Y', true, true, true, 2, -1,
  {BB,
   {"CB", E::C, false, TopologySidechain(Locant::Beta)},
   {"HB2",E::H, false, TopologyDiastereo(Locant::Beta, DiastereotopicIndex::Two,   ProchiralStereo::ProS, PseudoatomClass::QB)},
   {"HB3",E::H, false, TopologyDiastereo(Locant::Beta, DiastereotopicIndex::Three, ProchiralStereo::ProR, PseudoatomClass::QB)},
   {"CG", E::C, false, TopologyRing(Locant::Gamma,   BranchIndex::None, RingPosition::Substituent)},
   {"CD1",E::C, false, TopologyRing(Locant::Delta,   BranchIndex::One,  RingPosition::Member)},
   {"HD1",E::H, false, TopologySidechain(Locant::Delta,   BranchIndex::One)},
   {"CD2",E::C, false, TopologyRing(Locant::Delta,   BranchIndex::Two,  RingPosition::Member)},
   {"HD2",E::H, false, TopologySidechain(Locant::Delta,   BranchIndex::Two)},
   {"CE1",E::C, false, TopologyRing(Locant::Epsilon, BranchIndex::One,  RingPosition::Member)},
   {"HE1",E::H, false, TopologySidechain(Locant::Epsilon, BranchIndex::One)},
   {"CE2",E::C, false, TopologyRing(Locant::Epsilon, BranchIndex::Two,  RingPosition::Member)},
   {"HE2",E::H, false, TopologySidechain(Locant::Epsilon, BranchIndex::Two)},
   {"CZ", E::C, false, TopologyRing(Locant::Zeta,    BranchIndex::None, RingPosition::Member)},
   {"OH", E::O, false, TopologySidechain(Locant::Eta)},
   {"HH", E::H, false, TopologyPolarH(Locant::Eta, BranchIndex::None, PolarHKind::Hydroxyl)}},  // phenol O-H; absent in TYM
  {{RingTypeIndex::TyrPhenol, {"CG","CD1","CE1","CZ","CE2","CD2"}}},
  {{"N","CA","CB","CG"}, {"CA","CB","CG","CD1"}},
  {{"TYM", "deprotonated tyrosinate", -1, "deprotonated"}} },

// VAL — Cβ is real chiral (single Hβ). Cγ1 and Cγ2 are diastereotopic methyl
// carbons; pseudoatoms MG1 (Hγ11/12/13), MG2 (Hγ21/22/23), QG (all six γ-methyl).
{ AminoAcid::VAL, "VAL", 'V', false, false, true, 1, 0,
  {BB,
   {"CB",  E::C, false, TopologySidechain(Locant::Beta)},
   {"HB",  E::H, false, TopologySidechain(Locant::Beta)},
   {"CG1", E::C, false, TopologySidechain(Locant::Gamma, BranchIndex::One)},
   {"HG11",E::H, false, TopologyMethyl(Locant::Gamma, BranchIndex::One, PseudoatomClass::MG1)},
   {"HG12",E::H, false, TopologyMethyl(Locant::Gamma, BranchIndex::One, PseudoatomClass::MG1)},
   {"HG13",E::H, false, TopologyMethyl(Locant::Gamma, BranchIndex::One, PseudoatomClass::MG1)},
   {"CG2", E::C, false, TopologySidechain(Locant::Gamma, BranchIndex::Two)},
   {"HG21",E::H, false, TopologyMethyl(Locant::Gamma, BranchIndex::Two, PseudoatomClass::MG2)},
   {"HG22",E::H, false, TopologyMethyl(Locant::Gamma, BranchIndex::Two, PseudoatomClass::MG2)},
   {"HG23",E::H, false, TopologyMethyl(Locant::Gamma, BranchIndex::Two, PseudoatomClass::MG2)}},
  {}, {{"N","CA","CB","CG1"}}, {} },

// Unknown
{ AminoAcid::Unknown, "UNK", 'X', false, false, false, 0, 0,
  {}, {}, {}, {} },

};

#undef BB
#undef BB_PRO
#undef BB_GLY

const AminoAcidType& GetAminoAcidType(AminoAcid aa) {
    int i = static_cast<int>(aa);
    if (i < 0 || i >= static_cast<int>(AMINO_ACID_TYPES.size()))
        return AMINO_ACID_TYPES.back();
    return AMINO_ACID_TYPES[static_cast<size_t>(i)];
}

const std::vector<AminoAcidType>& AllAminoAcidTypes() {
    return AMINO_ACID_TYPES;
}

const AminoAcidType& AminoAcidTypeFromCode(const std::string& code) {
    for (const auto& type : AMINO_ACID_TYPES)
        if (code == type.three_letter_code) return type;
    return AMINO_ACID_TYPES.back();
}


void ValidateVariantIndices() {
    // Verify the variant ordering contract. A swap here silently breaks
    // every protonation assignment in the system.
    auto check = [](AminoAcid aa, int idx, const char* expected_name) {
        const auto& aat = GetAminoAcidType(aa);
        if (idx >= static_cast<int>(aat.variants.size()) ||
            std::string(aat.variants[idx].name) != expected_name) {
            std::fprintf(stderr,
                "FATAL: AminoAcidType variant index contract violated: "
                "%s variants[%d] expected \"%s\", got \"%s\"\n",
                aat.three_letter_code, idx, expected_name,
                (idx < static_cast<int>(aat.variants.size()))
                    ? aat.variants[idx].name : "(out of range)");
            std::abort();
        }
    };

    check(AminoAcid::HIS, 0, "HID");
    check(AminoAcid::HIS, 1, "HIE");
    check(AminoAcid::HIS, 2, "HIP");
    check(AminoAcid::ASP, 0, "ASH");
    check(AminoAcid::GLU, 0, "GLH");
    check(AminoAcid::CYS, 0, "CYX");
    check(AminoAcid::CYS, 1, "CYM");
    check(AminoAcid::LYS, 0, "LYN");
    check(AminoAcid::ARG, 0, "ARN");
    check(AminoAcid::TYR, 0, "TYM");
}

}  // namespace nmr
