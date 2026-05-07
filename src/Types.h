#pragma once
//
// Core types for NMR shielding tensor calculations.
//
// Vec3 and Mat3 are Eigen types -- no wrapper, no indirection.
// SphericalTensor decomposes a 3x3 tensor into irreducible representations.
// Element enumerates the nuclei present in protein NMR.
// Enums for bond classification, atom roles, hybridisation, ring identity.
//

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <string>

namespace nmr {

// ============================================================================
// Linear algebra from Eigen, used everywhere.
// ============================================================================

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;


// ============================================================================
// Element (enum) -- the 5 elements in protein NMR
// ============================================================================

enum class Element { H, C, N, O, S, Unknown };

inline Element ElementFromSymbol(const std::string& sym) {
    if (sym == "H") return Element::H;
    if (sym == "C") return Element::C;
    if (sym == "N") return Element::N;
    if (sym == "O") return Element::O;
    if (sym == "S") return Element::S;
    return Element::Unknown;
}

inline std::string SymbolForElement(Element e) {
    switch (e) {
        case Element::H: return "H";
        case Element::C: return "C";
        case Element::N: return "N";
        case Element::O: return "O";
        case Element::S: return "S";
        default: return "?";
    }
}

inline int AtomicNumberForElement(Element e) {
    switch (e) {
        case Element::H: return 1;
        case Element::C: return 6;
        case Element::N: return 7;
        case Element::O: return 8;
        case Element::S: return 16;
        default: return 0;
    }
}

// Compile-time element properties (OBJECT_MODEL.md)
inline double CovalentRadiusForElement(Element e) {
    switch (e) {
        case Element::H: return 0.31;
        case Element::C: return 0.76;
        case Element::N: return 0.71;
        case Element::O: return 0.66;
        case Element::S: return 1.05;
        default: return 0.0;
    }
}

inline double ElectronegativityForElement(Element e) {
    switch (e) {
        case Element::H: return 2.20;
        case Element::C: return 2.55;
        case Element::N: return 3.04;
        case Element::O: return 3.44;
        case Element::S: return 2.58;
        default: return 0.0;
    }
}


// ============================================================================
// Hybridisation (enum)
// ============================================================================

enum class Hybridisation { sp, sp2, sp3, Unassigned };


// ============================================================================
// AtomRole (enum) -- NMR-relevant classification of each atom
// ============================================================================

enum class AtomRole {
    // Heavy backbone
    BackboneN,
    BackboneCA,
    BackboneC,
    BackboneO,

    // Heavy sidechain
    SidechainC,
    SidechainN,
    SidechainO,
    SidechainS,
    AromaticC,
    AromaticN,

    // Hydrogen (by NMR-relevant environment)
    AmideH,
    AlphaH,
    MethylH,
    AromaticH,
    HydroxylH,
    OtherH,

    Unknown
};


// ============================================================================
// BondOrder (enum)
// ============================================================================

enum class BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
    Peptide,
    Unknown
};


// ============================================================================
// BondCategory (enum) -- finer than Constitution minimum, needed for McConnell
// ============================================================================

enum class BondCategory {
    PeptideCO,
    PeptideCN,
    BackboneOther,
    SidechainCO,
    Aromatic,
    Disulfide,
    SidechainOther,
    Unknown
};


// ============================================================================
// HeuristicTier (enum)
// ============================================================================

enum class HeuristicTier { REPORT, PASS, SILENT };


// ============================================================================
// RingAromaticity and RingSize
// ============================================================================

enum class RingSize { FiveMembered = 5, SixMembered = 6 };

// Aromaticity classification for ring chemistries.
// `None` covers saturated rings (Pro pyrrolidine); `Full` covers
// six-membered aromatic carbocycles and the fully-conjugated indole
// perimeter; `Reduced` covers five-membered aromatic heterocycles
// like the indole pyrrole subring; `Weak` covers imidazoles
// where ring-current strength is variant-sensitive (His HID/HIE/HIP).
// Joule & Mills, "Heterocyclic Chemistry" 5e (2010), chapters 7
// (saturated heterocycles), 13 (pyrroles, indoles) and 17 (imidazoles).
enum class RingAromaticity { Full, Reduced, Weak, None };


// ============================================================================
// RingTypeIndex (enum) -- 9 ring types
// ============================================================================
//
// Indices 0..kAromaticRingTypeCount-1 are aromatic ring chemistries;
// indices kAromaticRingTypeCount..Count-1 are saturated. Calculator
// per-aromatic-type accumulator arrays are sized
// `kAromaticRingTypeCount` and gate Pro out via the conventional
// `if (ti < kAromaticRingTypeCount)` guard.

enum class RingTypeIndex {
    PheBenzene     = 0,
    TyrPhenol      = 1,
    TrpBenzene     = 2,
    TrpPyrrole     = 3,
    TrpPerimeter   = 4,
    HisImidazole   = 5,
    HidImidazole   = 6,
    HieImidazole   = 7,
    ProPyrrolidine = 8,   ///< Saturated 5-ring (Pro pyrrolidine).
                          ///< Aromaticity = None; Intensity = 0
                          ///< (Joule & Mills 2010 ch. 7). Falls outside
                          ///< the `< kAromaticRingTypeCount` boundary
                          ///< for per-aromatic-type calculator
                          ///< accumulation.
    Count          = 9
};

/// Boundary between aromatic ring types (indices 0..7) and saturated
/// ring types (indices 8..Count-1). Saturated rings are excluded from
/// per-aromatic-type calculator aggregation by the conventional
/// `if (ti < kAromaticRingTypeCount)` guard.
///
/// **Adoption in calculator code is deferred to the per-calculator
/// update slices.** Bundle C / Slice A is substrate-side only: the
/// constant is declared here and documented as the design seam, but
/// `BiotSavartResult.cpp`, `HaighMallionResult.cpp`,
/// `PiQuadrupoleResult.cpp`, `DispersionResult.cpp`, and
/// `ConformationAtom.h`'s `std::array<double, 8>` per-aromatic-type
/// accumulators continue to use the literal `8` until each
/// calculator's own update slice picks them up. The literal `8` in
/// `ConformationAtom.h` is also the stable NPY ABI shape (eight
/// per-aromatic-type entries per atom in the trajectory output);
/// changing it requires a coordinated NPY schema migration.
///
/// The static_asserts below pin the constant to the cross-enum
/// invariant: `kAromaticRingTypeCount` IS the index where saturated
/// ring types begin, and `Count` must remain strictly greater. If a
/// future change adds a new aromatic ring at index 8 (shifting Pro
/// to 9), the first assert fails until `kAromaticRingTypeCount` is
/// updated to 9 too — surfacing the boundary shift at the build line
/// instead of leaving the calculator-side guards silently
/// off-by-one.
inline constexpr int kAromaticRingTypeCount = 8;
static_assert(static_cast<int>(RingTypeIndex::ProPyrrolidine) == kAromaticRingTypeCount,
              "kAromaticRingTypeCount IS the index of the first saturated "
              "ring type; if you re-order RingTypeIndex (e.g. inserting a "
              "new aromatic ring before ProPyrrolidine), update the "
              "constant in lockstep so the calculator-side aromatic "
              "boundary guards remain correct.");
static_assert(kAromaticRingTypeCount < static_cast<int>(RingTypeIndex::Count),
              "kAromaticRingTypeCount must be a valid ring type index, "
              "not the Count sentinel.");

inline const char* RingTypeName(RingTypeIndex t) {
    switch (t) {
        case RingTypeIndex::PheBenzene:     return "PHE";
        case RingTypeIndex::TyrPhenol:      return "TYR";
        case RingTypeIndex::TrpBenzene:     return "TRP6";
        case RingTypeIndex::TrpPyrrole:     return "TRP5";
        case RingTypeIndex::TrpPerimeter:   return "TRP9";
        case RingTypeIndex::HisImidazole:   return "HIS";
        case RingTypeIndex::HidImidazole:   return "HID";
        case RingTypeIndex::HieImidazole:   return "HIE";
        case RingTypeIndex::ProPyrrolidine: return "PRO";
        default: return "?";
    }
}


// ============================================================================
// AminoAcid (enum) -- the 20 standard amino acids
// ============================================================================

enum class AminoAcid {
    ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY,
    HIS, ILE, LEU, LYS, MET, PHE, PRO, SER,
    THR, TRP, TYR, VAL, Unknown
};

AminoAcid AminoAcidFromThreeLetterCode(const std::string& code);
std::string ThreeLetterCodeForAminoAcid(AminoAcid aa);
bool IsAromaticAminoAcid(AminoAcid aa);


// ============================================================================
// SphericalTensor -- irreducible decomposition of a 3x3 tensor
//
// Any 3x3 tensor sigma decomposes uniquely into three parts:
//   T0: scalar (trace/3) -- the isotropic component
//   T1: pseudovector (antisymmetric part) -- 3 components
//   T2: traceless symmetric tensor -- 5 components (m = -2..+2)
//
// T2 uses isometric normalization (real spherical harmonics).
// This preserves the L2 norm: sum(|T2_m|^2) == sum(S_ij^2).
// ============================================================================

struct SphericalTensor {
    double T0 = 0.0;
    std::array<double, 3> T1 = {};
    std::array<double, 5> T2 = {};

    double Isotropic() const { return T0; }
    const std::array<double, 3>& Antisymmetric() const { return T1; }
    const std::array<double, 5>& TracelessSymmetric() const { return T2; }

    // Decompose a 3x3 tensor into T0 + T1 + T2.
    static SphericalTensor Decompose(const Mat3& tensor);

    // Reconstruct the 3x3 tensor from T0, T1, T2.
    Mat3 Reconstruct() const;

    // L2 norm of the T2 components.
    double T2Magnitude() const;
};


// ============================================================================
// FieldValue -- a calculator result attributed to a specific source
// ============================================================================

enum class CalculatorId {
    BiotSavart, HaighMallion, McConnell, Coulomb,
    PiQuadrupole, RingSusceptibility, Dispersion, HBond,
    APBS, Orca, Mopac, MopacCoulomb, MopacMcConnell, AIMNet2, SASA,
    WaterField, HydrationShell, HydrationGeometry, EEQ
};

struct FieldValue {
    Mat3 tensor = Mat3::Zero();
    SphericalTensor spherical;
    CalculatorId source_calculator = CalculatorId::BiotSavart;
    size_t source_index = 0;
};


// ============================================================================
// ProtonationTool (enum)
// ============================================================================

enum class ProtonationTool { PROPKA, KaML, TLeap, Manual };


}  // namespace nmr
