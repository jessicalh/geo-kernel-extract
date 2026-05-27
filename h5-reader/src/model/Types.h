// Types.h — core typed primitives for the h5-reader object model.
//
// Mirrors `nmr/src/Types.h` field-for-field, ordinal-compatible. The
// chemistry-substrate vocabulary (Locant, BackboneRole, PlanarGroupKind,
// etc.) is in the sibling `QtSemanticEnums.h` — same split as the
// library (`src/Types.h` vs `src/SemanticEnums.h`).
//
// String policy: ZERO string-dispatch on chemistry identity. The
// NameForXxx() helpers in this file produce `const char*` literals for
// display surfaces (inspector dock, tooltips) — they MUST NEVER appear
// in a comparison or branch that decides physics. The design rule
// driving this discipline is in `notes/H5_READER_REWRITE_DESIGN_2026-05-23.md`
// §2 ("The No-Strings Discipline").
//
// Compile-time ordinal compatibility with the library is the runtime
// contract: the H5 + sidecar NPY format stores int8/int32 enum
// ordinals directly, and the loader casts them back to these enums.
// Drift between an enum here and its library counterpart silently
// mis-decodes every atom/bond/ring using it. The
// extraction_manifest.json's `enum_vocab` block is checked at load
// time by `io/QtEnumVocab` to catch drift early.

#pragma once

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstdint>

namespace h5reader::model {

// ============================================================================
// Linear algebra — same Eigen types the library uses.
// ============================================================================

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;


// ============================================================================
// Element — atomic species for protein NMR.
// Library mirror: nmr::Element (src/Types.h:30).
// Storage: atoms_category_info.element AND H5 /atoms/element both store
// the int8 ATOMIC NUMBER (1, 6, 7, 8, 16) — NOT the enum ordinal —
// decoded via ElementFromAtomicNumber() at the loader boundary
// (CategoryInfoProjection.cpp:417 emits AtomicNumberForElement;
// QtTopologySidecar.cpp:188 decodes). A naive static_cast<Element>(row)
// would mis-decode every atom (e.g. C atomic-number 6 -> Element ord 6 =
// Unknown), so keep the ElementFromAtomicNumber() decode.
// ============================================================================

enum class Element : int8_t { H = 0, C = 1, N = 2, O = 3, S = 4, Unknown = 5 };

inline int AtomicNumberForElement(Element e) {
    switch (e) {
    case Element::H:
        return 1;
    case Element::C:
        return 6;
    case Element::N:
        return 7;
    case Element::O:
        return 8;
    case Element::S:
        return 16;
    default:
        return 0;
    }
}
inline Element ElementFromAtomicNumber(int z) {
    switch (z) {
    case 1:
        return Element::H;
    case 6:
        return Element::C;
    case 7:
        return Element::N;
    case 8:
        return Element::O;
    case 16:
        return Element::S;
    default:
        return Element::Unknown;
    }
}
inline const char* SymbolForElement(Element e) {
    switch (e) {
    case Element::H:
        return "H";
    case Element::C:
        return "C";
    case Element::N:
        return "N";
    case Element::O:
        return "O";
    case Element::S:
        return "S";
    default:
        return "?";
    }
}
inline double CovalentRadiusForElement(Element e) {
    // Bondi 1964 / Cordero 2008.
    switch (e) {
    case Element::H:
        return 0.31;
    case Element::C:
        return 0.76;
    case Element::N:
        return 0.71;
    case Element::O:
        return 0.66;
    case Element::S:
        return 1.05;
    default:
        return 0.0;
    }
}
inline double ElectronegativityForElement(Element e) {
    switch (e) {
    case Element::H:
        return 2.20;
    case Element::C:
        return 2.55;
    case Element::N:
        return 3.04;
    case Element::O:
        return 3.44;
    case Element::S:
        return 2.58;
    default:
        return 0.0;
    }
}


// ============================================================================
// Hybridisation — library mirror nmr::Hybridisation (src/Types.h:91).
// Not in the sidecar projection; surfaced here for completeness.
// ============================================================================

enum class Hybridisation : int8_t {
    sp = 0,
    sp2 = 1,
    sp3 = 2,
    Unassigned = 3,
};


// ============================================================================
// AtomRole — NMR-relevant atom classification.
// Library mirror: nmr::AtomRole (src/Types.h:98).
// Not in atoms_category_info projection (the typed substrate
// BackboneRole + Locant + Element supersedes it). Kept for legacy
// consumers; new code should prefer the substrate fields on QtAtom.
// ============================================================================

enum class AtomRole : int8_t {
    BackboneN = 0,
    BackboneCA = 1,
    BackboneC = 2,
    BackboneO = 3,
    SidechainC = 4,
    SidechainN = 5,
    SidechainO = 6,
    SidechainS = 7,
    AromaticC = 8,
    AromaticN = 9,
    AmideH = 10,
    AlphaH = 11,
    MethylH = 12,
    AromaticH = 13,
    HydroxylH = 14,
    OtherH = 15,
    Unknown = 16,
};

inline const char* NameForAtomRole(AtomRole r) {
    switch (r) {
    case AtomRole::BackboneN:
        return "BackboneN";
    case AtomRole::BackboneCA:
        return "BackboneCA";
    case AtomRole::BackboneC:
        return "BackboneC";
    case AtomRole::BackboneO:
        return "BackboneO";
    case AtomRole::SidechainC:
        return "SidechainC";
    case AtomRole::SidechainN:
        return "SidechainN";
    case AtomRole::SidechainO:
        return "SidechainO";
    case AtomRole::SidechainS:
        return "SidechainS";
    case AtomRole::AromaticC:
        return "AromaticC";
    case AtomRole::AromaticN:
        return "AromaticN";
    case AtomRole::AmideH:
        return "AmideH";
    case AtomRole::AlphaH:
        return "AlphaH";
    case AtomRole::MethylH:
        return "MethylH";
    case AtomRole::AromaticH:
        return "AromaticH";
    case AtomRole::HydroxylH:
        return "HydroxylH";
    case AtomRole::OtherH:
        return "OtherH";
    case AtomRole::Unknown:
        return "Unknown";
    }
    return "?";
}

inline const char* NameForHybridisation(Hybridisation h) {
    switch (h) {
    case Hybridisation::sp:
        return "sp";
    case Hybridisation::sp2:
        return "sp2";
    case Hybridisation::sp3:
        return "sp3";
    case Hybridisation::Unassigned:
        return "unassigned";
    }
    return "?";
}


// ============================================================================
// BondOrder — library mirror nmr::BondOrder (src/Types.h:129).
// Sidecar storage: bonds.bond_order (int8). Manifest enum_vocab:
// {0:Single, 1:Double, 2:Triple, 3:Aromatic, 4:Peptide, 5:Unknown}.
// ============================================================================

enum class BondOrder : int8_t {
    Single = 0,
    Double = 1,
    Triple = 2,
    Aromatic = 3,
    Peptide = 4,
    Unknown = 5,
};


// ============================================================================
// BondCategory — library mirror nmr::BondCategory (src/Types.h:143).
// Sidecar storage: bonds.bond_category (int8). Manifest enum_vocab:
// {0:PeptideCO, 1:PeptideCN, 2:BackboneOther, 3:SidechainCO,
//  4:Aromatic, 5:Disulfide, 6:SidechainOther, 7:Unknown}.
// ============================================================================

enum class BondCategory : int8_t {
    PeptideCO = 0,
    PeptideCN = 1,
    BackboneOther = 2,
    SidechainCO = 3,
    Aromatic = 4,
    Disulfide = 5,
    SidechainOther = 6,
    Unknown = 7,
};


// ============================================================================
// RingAromaticity — library mirror nmr::RingAromaticity (src/Types.h:176).
// ============================================================================

enum class RingAromaticity : int8_t { Full = 0, Reduced = 1, Weak = 2, None = 3 };


// ============================================================================
// RingTypeIndex — 9 ring chemistries (8 aromatic + 1 saturated).
// Library mirror: nmr::RingTypeIndex (src/Types.h:189).
// Sidecar storage: rings.ring_type_index (int8). Manifest enum_vocab:
// {0:PheBenzene, 1:TyrPhenol, 2:TrpBenzene, 3:TrpPyrrole,
//  4:TrpPerimeter, 5:HisImidazole, 6:HidImidazole, 7:HieImidazole,
//  8:ProPyrrolidine}.
// ============================================================================

enum class RingTypeIndex : int8_t {
    PheBenzene = 0,
    TyrPhenol = 1,
    TrpBenzene = 2,
    TrpPyrrole = 3,
    TrpPerimeter = 4,
    HisImidazole = 5,
    HidImidazole = 6,
    HieImidazole = 7,
    ProPyrrolidine = 8,
};
constexpr int kRingTypeCount = 9;
constexpr int kAromaticRingTypeCount = 8;


// ============================================================================
// RingKind — aromatic vs saturated discriminator.
// Sidecar storage: rings.ring_kind (int8). Manifest enum_vocab:
// {0:aromatic, 1:saturated}.
// ============================================================================

enum class RingKind : int8_t { Aromatic = 0, Saturated = 1 };


// ============================================================================
// AminoAcid — 20 standard amino acids + Unknown.
// Library mirror: nmr::AminoAcid (src/Types.h:263).
// Sidecar storage: residues.residue_type / atoms_category_info.residue_type
// (int8). Manifest enum_vocab confirms the ordering.
// ============================================================================

enum class AminoAcid : int8_t {
    ALA = 0,
    ARG = 1,
    ASN = 2,
    ASP = 3,
    CYS = 4,
    GLN = 5,
    GLU = 6,
    GLY = 7,
    HIS = 8,
    ILE = 9,
    LEU = 10,
    LYS = 11,
    MET = 12,
    PHE = 13,
    PRO = 14,
    SER = 15,
    THR = 16,
    TRP = 17,
    TYR = 18,
    VAL = 19,
    Unknown = 20,
};
constexpr int kStandardAminoAcidCount = 20;


// ============================================================================
// DsspCode — DSSP 8-class secondary structure.
// Source: per-residue ss8 stored as uint8 in
// /trajectory/dssp8_time_series/ss8_code. The group's `ss8_legend` attr
// in the 1P9J fixture spells out the ordinals: "H=0 (alpha helix),
// G=1 (3_10 helix), I=2 (pi helix), E=3 (extended strand),
// B=4 (beta bridge), T=5 (turn), S=6 (bend), C=7 (coil)".
// ss8_unassigned_sentinel = 255 → Unknown.
// ============================================================================

enum class DsspCode : uint8_t {
    AlphaHelix = 0,
    Helix310 = 1,
    PiHelix = 2,
    ExtendedStrand = 3,
    BetaBridge = 4,
    Turn = 5,
    Bend = 6,
    Coil = 7,
    Unknown = 255,
};

inline char OneLetterForDssp(DsspCode c) {
    switch (c) {
    case DsspCode::AlphaHelix:
        return 'H';
    case DsspCode::Helix310:
        return 'G';
    case DsspCode::PiHelix:
        return 'I';
    case DsspCode::ExtendedStrand:
        return 'E';
    case DsspCode::BetaBridge:
        return 'B';
    case DsspCode::Turn:
        return 'T';
    case DsspCode::Bend:
        return 'S';
    case DsspCode::Coil:
        return 'C';
    case DsspCode::Unknown:
        return '?';
    }
    return '?';
}


// ============================================================================
// NamingProvenance — atom-name projection accuracy.
// Library mirror: nmr::NamingProvenance (src/CategoryInfoProjection.cpp:41).
// Sidecar storage: atoms_category_info.iupac_naming_provenance,
// bmrb_naming_provenance (int8 each).
// ============================================================================

enum class NamingProvenance : int8_t {
    Match = 0,       // canonical name found in atom_nom.tbl
    MissLogged = 1,  // no row; AMBER name returned as fallback
};


// ============================================================================
// NamingConvention / NamingSource — selectors for the name PROJECTIONS.
// Labels are projections, never identity: the reader is NOT label-driven,
// the typed AminoAcid/substrate is the single source of truth, and these
// pick only which display / ML-join label to surface. See QtAtomNames.h +
// QtResidueNames.h.
// ============================================================================

enum class NamingConvention : int8_t { Amber = 0, Iupac = 1, Bmrb = 2 };

// Verbatim = the string the library's CategoryInfoProjection wrote into the
//            sidecar (atom_nom.tbl-driven), read as-deposited.
// Derived  = recomputed reader-side from the typed AminoAcid (+ variant) via
//            the QtResidueNames helpers. Atom names have no Derived source
//            (the reader does not re-derive atom names; the library did).
enum class NamingSource : int8_t { Verbatim = 0, Derived = 1 };


// ============================================================================
// QtFfAtomType — typed AMBER ff14SB atom-type vocabulary.
//
// Replaces the bare std::string `ff_atom_type_string` column with a
// typed enum. The ff14SB vocabulary is finite and stable within the
// force field; new types from ff19SB or other AMBER variants would
// extend this enum + add Unknown-handling.
//
// Source: AMBER ff14SB atom types (Maier et al. 2015,
// J. Chem. Theory Comput. 11(8) 3696–3713 — Table 1 and supplementary).
// Sidecar storage: atoms_category_info.ff_atom_type_string (S4 ASCII).
// Loader maps the S4 string to this enum; mismatches fall to Unknown
// with an ErrorBus warning.
//
// Discipline note: rationale for typed enum over QString is in
// notes/H5_READER_REWRITE_DESIGN_2026-05-23.md §11.B — closes the
// last string-dispatch vulnerability.
// ============================================================================

enum class QtFfAtomType : int16_t {
    Unknown = 0,

    // Carbons
    C = 1,   // carbonyl C
    CA = 2,  // alpha C
    CB = 3,
    CC = 4,  // HIS sidechain
    CD = 5,  // ARG, LYS, PRO sidechain
    CK = 6,
    CM = 7,
    CN = 8,  // TRP indole
    CO = 9,
    CP = 10,
    CQ = 11,
    CR = 12,  // HIS imidazole
    CT = 13,  // sp3 aliphatic C
    CV = 14,  // HIE / HIP sidechain
    CW = 15,  // HID / HIP sidechain, TRP
    CX = 16,
    CY = 17,
    CZ = 18,
    Cstar = 19,  // C* — TRP CG

    // Nitrogens
    N = 20,   // amide N
    N2 = 21,  // ARG, sp2
    N3 = 22,  // sp3, N-term ammonium
    NA = 23,  // HIE, TRP NE1
    NB = 24,  // HID
    NC = 25,
    NP = 26,
    NT = 27,
    NY = 28,
    Nstar = 29,

    // Hydrogens
    H = 30,   // amide H, on N
    H1 = 31,  // on CT with one EWG
    H2 = 32,  // on CT with two EWGs
    H3 = 33,
    H4 = 34,  // on aromatic ring with one EWG
    H5 = 35,
    HA = 36,  // on sp2 aromatic ring
    HC = 37,  // on sp3 C without EWG
    HO = 38,  // on OH
    HP = 39,
    HS = 40,  // on SH
    HW = 41,  // water H
    HZ = 42,

    // Oxygens
    O = 43,   // carbonyl O
    O2 = 44,  // carboxylate
    OH = 45,
    OS = 46,
    OW = 47,  // water O
    OP = 48,
    OD = 49,

    // Sulfurs
    S = 50,   // CYS, MET sulfide
    SH = 51,  // CYS thiol
};


// ============================================================================
// SphericalTensor — irreducible decomposition of a 3x3 tensor.
// Library mirror: nmr::SphericalTensor (src/Types.h:286).
//
// 9-component layout [T0, T1[3], T2[5]]. T2 is real-spherical-tesseral
// (m = -2..+2), isometric-normalised (matches e3nn). T0 and |T2| are
// basis-independent.
//
// T1 basis — TWO paths, do not conflate:
//  * Per-frame NPY path (the result-group views over QtConformationSnapshot):
//    columns are the library's PackFull9 straight off Decompose, so T1[0..2]
//    is the CARTESIAN antisymmetric pseudovector (v_x=A_yz, v_y=A_zx,
//    v_z=A_xy; src/Types.cpp). The views pass it through raw — a T1 from a
//    group view is Cartesian.
//  * H5 per-TR `xyz` path (dense time-series): the `irrep_layout` attr labels
//    T1 "T1_m-1,T1_m0,T1_m+1" (m-basis). Whether that path re-orders to
//    m-basis or also carries Cartesian is a library convention to CONFIRM
//    against QtTrajectoryH5 — flagged, not assumed.
//
// Parity: T1 (antisymmetric pseudovector) is axial -> e3nn 1e (even); the
// whole tensor is 0e+1e+2e (all-even), per _catalog.py. The SDK _tensors.py
// declares 1o for T1 — a library-side inconsistency (a definite-parity tensor
// cannot mix 1o with 0e/2e); the catalog's 1e is correct.
// ============================================================================

struct SphericalTensor {
    double T0 = 0.0;
    std::array<double, 3> T1 = {};
    std::array<double, 5> T2 = {};

    double T2Magnitude() const {
        double s = 0.0;
        for (double v : T2)
            s += v * v;
        return std::sqrt(s);
    }
};


}  // namespace h5reader::model
