// Types.h — core typed primitives for the h5-reader object model.
//
// Every enum here is ordinal-compatible with a counterpart in the
// nmr-shielding library's src/Types.h. The H5 writer emits library enum
// ordinals directly (int32 or int8); the loader casts back to these.
// If you add a value here or there without matching the other, every
// atom, bond, or ring using that enum will mis-decode silently.
//
// The round-trip test in tests/EnumOrdinalTest.cpp fails on drift —
// run it when touching anything here.
//
// Vec3 / Mat3 are Eigen types, same as the library. SphericalTensor is
// the irreducible-representation decomposition (T0 isotropic, T1
// antisymmetric, T2 traceless symmetric). Decompose/Reconstruct helpers
// will land in a separate .cpp when a renderer needs them; the struct
// is declared here so slab accessors can return one by value.

#pragma once

#include <Eigen/Dense>
#include <array>
#include <cstdint>
#include <string>

namespace h5reader::model {

// ============================================================================
// Linear algebra — same Eigen types the library uses.
// ============================================================================

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;


// ============================================================================
// Element — atomic species for protein NMR.
//
// H5 stores atomic number directly (atoms/element dataset, int32):
//   1 → H, 6 → C, 7 → N, 8 → O, 16 → S.
//
// Library mirror: nmr::Element. Decoder: AtomicNumberToElement() in
// io/QtProteinLoader.cpp.
// ============================================================================

enum class Element { H, C, N, O, S, Unknown };

inline int AtomicNumberForElement(Element e) {
    switch (e) {
        case Element::H: return 1;
        case Element::C: return 6;
        case Element::N: return 7;
        case Element::O: return 8;
        case Element::S: return 16;
        default:         return 0;
    }
}

inline const char* SymbolForElement(Element e) {
    switch (e) {
        case Element::H: return "H";
        case Element::C: return "C";
        case Element::N: return "N";
        case Element::O: return "O";
        case Element::S: return "S";
        default:         return "?";
    }
}

// Covalent radius in Angstroms. Bondi 1964 / Cordero 2008.
inline double CovalentRadiusForElement(Element e) {
    switch (e) {
        case Element::H: return 0.31;
        case Element::C: return 0.76;
        case Element::N: return 0.71;
        case Element::O: return 0.66;
        case Element::S: return 1.05;
        default:         return 0.0;
    }
}

// Pauling electronegativity.
inline double ElectronegativityForElement(Element e) {
    switch (e) {
        case Element::H: return 2.20;
        case Element::C: return 2.55;
        case Element::N: return 3.04;
        case Element::O: return 3.44;
        case Element::S: return 2.58;
        default:         return 0.0;
    }
}


// ============================================================================
// Hybridisation — OpenBabel-assigned in the library's EnrichmentResult.
//
// Ordinal-compatible with nmr::Hybridisation (src/Types.h:91).
// ============================================================================

enum class Hybridisation {
    sp        = 0,
    sp2       = 1,
    sp3       = 2,
    Unassigned = 3
};

inline const char* NameForHybridisation(Hybridisation h) {
    switch (h) {
        case Hybridisation::sp:         return "sp";
        case Hybridisation::sp2:        return "sp2";
        case Hybridisation::sp3:        return "sp3";
        case Hybridisation::Unassigned: return "unassigned";
    }
    return "?";
}


// ============================================================================
// AtomRole — NMR-relevant classification of each atom.
//
// Ordinal-compatible with nmr::AtomRole (src/Types.h:98..122).
// H5: atoms/atom_role (int32). The library's Unknown sits at the end
// (value 16); we keep that contract.
// ============================================================================

enum class AtomRole : int32_t {
    BackboneN   = 0,
    BackboneCA  = 1,
    BackboneC   = 2,
    BackboneO   = 3,
    SidechainC  = 4,
    SidechainN  = 5,
    SidechainO  = 6,
    SidechainS  = 7,
    AromaticC   = 8,
    AromaticN   = 9,
    AmideH      = 10,
    AlphaH      = 11,
    MethylH     = 12,
    AromaticH   = 13,
    HydroxylH   = 14,
    OtherH      = 15,
    Unknown     = 16,
};

inline const char* NameForAtomRole(AtomRole r) {
    switch (r) {
        case AtomRole::BackboneN:  return "BackboneN";
        case AtomRole::BackboneCA: return "BackboneCA";
        case AtomRole::BackboneC:  return "BackboneC";
        case AtomRole::BackboneO:  return "BackboneO";
        case AtomRole::SidechainC: return "SidechainC";
        case AtomRole::SidechainN: return "SidechainN";
        case AtomRole::SidechainO: return "SidechainO";
        case AtomRole::SidechainS: return "SidechainS";
        case AtomRole::AromaticC:  return "AromaticC";
        case AtomRole::AromaticN:  return "AromaticN";
        case AtomRole::AmideH:     return "AmideH";
        case AtomRole::AlphaH:     return "AlphaH";
        case AtomRole::MethylH:    return "MethylH";
        case AtomRole::AromaticH:  return "AromaticH";
        case AtomRole::HydroxylH:  return "HydroxylH";
        case AtomRole::OtherH:     return "OtherH";
        case AtomRole::Unknown:    return "Unknown";
    }
    return "?";
}

inline bool IsBackboneRole(AtomRole r) {
    return r == AtomRole::BackboneN || r == AtomRole::BackboneCA
        || r == AtomRole::BackboneC || r == AtomRole::BackboneO;
}

inline bool IsHydrogenRole(AtomRole r) {
    return r == AtomRole::AmideH || r == AtomRole::AlphaH
        || r == AtomRole::MethylH || r == AtomRole::AromaticH
        || r == AtomRole::HydroxylH || r == AtomRole::OtherH;
}


// ============================================================================
// BondOrder — covalent bond multiplicity.
//
// Ordinal-compatible with nmr::BondOrder (src/Types.h:129..136).
// ============================================================================

enum class BondOrder : int32_t {
    Single   = 0,
    Double   = 1,
    Triple   = 2,
    Aromatic = 3,
    Peptide  = 4,
    Unknown  = 5,
};


// ============================================================================
// BondCategory — finer classification than BondOrder, used by McConnell
// and for visual bond styling.
//
// Ordinal-compatible with nmr::BondCategory (src/Types.h:143..152).
// ============================================================================

enum class BondCategory : int32_t {
    PeptideCO      = 0,
    PeptideCN      = 1,
    BackboneOther  = 2,
    SidechainCO    = 3,
    Aromatic       = 4,
    Disulfide      = 5,
    SidechainOther = 6,
    Unknown        = 7,
};

inline const char* NameForBondCategory(BondCategory c) {
    switch (c) {
        case BondCategory::PeptideCO:      return "PeptideCO";
        case BondCategory::PeptideCN:      return "PeptideCN";
        case BondCategory::BackboneOther:  return "BackboneOther";
        case BondCategory::SidechainCO:    return "SidechainCO";
        case BondCategory::Aromatic:       return "Aromatic";
        case BondCategory::Disulfide:      return "Disulfide";
        case BondCategory::SidechainOther: return "SidechainOther";
        case BondCategory::Unknown:        return "Unknown";
    }
    return "?";
}


// ============================================================================
// RingAromaticity — PiElectronCount-derived classification.
//
// Ordinal-compatible with nmr::RingAromaticity (src/Types.h:168).
// ============================================================================

enum class RingAromaticity { Full, Reduced, Weak };


// ============================================================================
// RingTypeIndex — 8 concrete ring types.
//
// Ordinal-compatible with nmr::RingTypeIndex (src/Types.h:175..185).
// H5: topology/ring_type (int32) and per_ring/ring_type (int8).
// The factory CreateQtRing(RingTypeIndex) in QtRing.cpp instantiates
// the matching concrete subclass.
// ============================================================================

enum class RingTypeIndex : int32_t {
    PheBenzene    = 0,
    TyrPhenol     = 1,
    TrpBenzene    = 2,
    TrpPyrrole    = 3,
    TrpPerimeter  = 4,
    HisImidazole  = 5,
    HidImidazole  = 6,
    HieImidazole  = 7,
};

constexpr int RingTypeCount = 8;


// ============================================================================
// AminoAcid — 20 standard residues.
//
// Ordinal-compatible with nmr::AminoAcid (src/Types.h:206..210).
// H5: residues/residue_name (variable-length strings); decoded by
// io/QtNamingRegistry.cpp's ThreeLetterCodeToAminoAcid().
// ============================================================================

enum class AminoAcid : int32_t {
    ALA = 0, ARG, ASN, ASP, CYS, GLN, GLU, GLY,
    HIS,     ILE, LEU, LYS, MET, PHE, PRO, SER,
    THR,     TRP, TYR, VAL,
    Unknown
};

constexpr int StandardAminoAcidCount = 20;


// ============================================================================
// ProtonationVariant — per-residue protonation state.
//
// Inferred at load time from ring subclass identity (HIS tautomers) and
// atom presence (ASP/GLU/LYS/TYR/CYS). Never from residue name string
// comparison. See io/QtProteinLoader.cpp::InferProtonationVariant.
// ============================================================================

enum class ProtonationVariant {
    Default = 0,                 // standard protonation for this residue
    HIS_delta,                   // HID — ND1 protonated
    HIS_epsilon,                 // HIE — NE2 protonated
    HIS_doubly,                  // HIP — both protonated
    ASP_protonated,              // ASH — neutral ASP
    GLU_protonated,              // GLH — neutral GLU
    LYS_deprotonated,            // LYN — neutral LYS
    ARG_deprotonated,            // ARN — neutral ARG (rare)
    TYR_deprotonated,            // TYM — anionic TYR
    CYS_disulfide,               // CYX — bonded to another CYS
    CYS_deprotonated,            // CYM — thiolate
};


// ============================================================================
// DsspCode — 8-class secondary structure.
//
// DSSP classic: H α-helix, E β-strand, G 3₁₀-helix, I π-helix, T turn,
// S bend, B β-bridge, ' ' (or '-') coil.
// H5: dssp/ss8 (int8). The library writes its own 0..7 ordinals; this
// enum mirrors that ordering.
// ============================================================================

enum class DsspCode : int8_t {
    Coil      = 0,   // ' '
    AlphaHelix= 1,   // 'H'
    BetaStrand= 2,   // 'E'
    Helix310  = 3,   // 'G'
    PiHelix   = 4,   // 'I'
    Turn      = 5,   // 'T'
    Bend      = 6,   // 'S'
    BetaBridge= 7,   // 'B'
    Unknown   = 8,
};

inline char OneLetterForDssp(DsspCode c) {
    switch (c) {
        case DsspCode::Coil:       return '-';
        case DsspCode::AlphaHelix: return 'H';
        case DsspCode::BetaStrand: return 'E';
        case DsspCode::Helix310:   return 'G';
        case DsspCode::PiHelix:    return 'I';
        case DsspCode::Turn:       return 'T';
        case DsspCode::Bend:       return 'S';
        case DsspCode::BetaBridge: return 'B';
        case DsspCode::Unknown:    return '?';
    }
    return '?';
}


// ============================================================================
// SphericalTensor — irreducible decomposition of a 3x3 tensor.
//
// T0 (scalar): isotropic component (trace/3).
// T1 (3 components): antisymmetric pseudovector (sigma_ij - sigma_ji).
// T2 (5 components, m = -2..+2): traceless symmetric, isometric
//     normalisation matching sphericart / the library convention.
//
// H5 layout per dataset attribute: "T0,T1[3],T2[5]" → 9 doubles flat.
// Slab accessor (in QtFrame) returns one of these by value per atom.
// ============================================================================

struct SphericalTensor {
    double                 T0 = 0.0;
    std::array<double, 3>  T1 = {};
    std::array<double, 5>  T2 = {};
};


}  // namespace h5reader::model
