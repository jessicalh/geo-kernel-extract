// QtAtom — per-atom typed substrate, mirroring the projected fields of
// atoms_category_info.npy.
//
// This struct IS the typed identity for an atom. Every chemistry-typed
// field (Element, Locant, BranchAddress, DiastereotopicIndex,
// BackboneRole, PlanarGroupKind, PolarHKind, ProchiralStereo,
// PlanarStereo, ring positions, PseudoatomKind, aromatic,
// formal_charge, is_exchangeable, equivalence_class) is a typed enum
// or POD. The atom IS its (Element, Locant, branch, di_index,
// backbone_role) tuple plus the chemistry substrate around it.
//
// NO NAME STRINGS LIVE HERE. The three atom-name strings (amber,
// iupac, bmrb) live on `QtAtomNames` (in QtAtomNames.h), accessed via
// `QtProtein::atomNames(i)` — an explicit projection boundary. Code
// that asks chemistry questions stays on QtAtom; code that displays a
// label asks for QtAtomNames. The two surfaces never mix.
//
// ff_atom_type uses the typed `QtFfAtomType` enum (per design
// §11.B), not the raw S4 string. Loader maps the AMBER vocabulary
// string to the enum at the boundary.
//
// See notes/H5_READER_REWRITE_DESIGN_2026-05-23.md §4.2 for the design.

#pragma once

#include "QtSemanticEnums.h"
#include "Types.h"

#include <cstdint>

namespace h5reader::model {

struct QtAtom {
    // ----- Identity indices -----
    int32_t atomIndex = -1;        // row index; equals position in QtProtein.atoms()
    int32_t residueIndex = -1;     // row into QtProtein.residues()
    int32_t parentAtomIndex = -1;  // for H atoms, heavy-atom parent; -1 otherwise

    // ----- Mechanical-identity tuple (mirrors AtomMechanicalIdentity) -----
    Element element = Element::Unknown;
    Locant locant = Locant::None;
    BranchAddress branch = {};
    DiastereotopicIndex diIndex = DiastereotopicIndex::None;
    BackboneRole backboneRole = BackboneRole::None;

    // ----- Stereo / planarity -----
    ProchiralStereo prochiral = ProchiralStereo::NotProchiral;
    PlanarGroupKind planarGroup = PlanarGroupKind::None;
    PlanarStereo planarStereo = PlanarStereo::NotApplicable;
    PolarHKind polarH = PolarHKind::NotPolar;

    // ----- Ring membership labels -----
    RingPositionLabel ringPositionPrimary = RingPositionLabel::NotInRing;
    RingPositionLabel ringPositionSecondary = RingPositionLabel::NotInRing;

    // ----- Pseudoatom (Markley 1998 Table 1) -----
    PseudoatomKind pseudoatomKind = PseudoatomKind::None;
    bool inSuperGroup = false;

    // ----- Other typed flags -----
    bool aromatic = false;
    int8_t formalCharge = 0;
    bool isExchangeable = false;
    int8_t equivalenceClass = 0;
    bool hasPartialCharge = false;
    double partialCharge = 0.0;  // FF14SB/topol.top charge, elementary-charge units.

    // ----- Force-field surface (typed enum, NOT string) -----
    QtFfAtomType ffAtomType = QtFfAtomType::Unknown;

    // ----- Element-property convenience (typed; computed) -----
    double CovalentRadius() const { return CovalentRadiusForElement(element); }
    double Electronegativity() const { return ElectronegativityForElement(element); }
    int AtomicNumber() const { return AtomicNumberForElement(element); }
    bool IsHBondDonorElement() const { return element == Element::N || element == Element::O; }
    bool IsHBondAcceptorElement() const { return element == Element::N || element == Element::O; }

    // ----- Substrate predicates (mirror AtomSemanticTable's; src/SemanticEnums.h:902) -----
    constexpr bool IsBackbone() const { return backboneRole != BackboneRole::None; }
    constexpr bool IsBackboneNitrogen() const { return backboneRole == BackboneRole::Nitrogen; }
    constexpr bool IsBackboneAlphaCarbon() const { return backboneRole == BackboneRole::AlphaCarbon; }
    constexpr bool IsBackboneCarbonylCarbon() const { return backboneRole == BackboneRole::CarbonylCarbon; }
    constexpr bool IsBackboneCarbonylOxygen() const { return backboneRole == BackboneRole::CarbonylOxygen; }
    constexpr bool IsBackboneAmideHydrogen() const { return backboneRole == BackboneRole::AmideHydrogen; }
    constexpr bool IsBackboneAlphaHydrogen() const { return backboneRole == BackboneRole::AlphaHydrogen; }
    constexpr bool IsAnyAlphaHydrogen() const {
        // Non-GLY HA via BackboneRole::AlphaHydrogen; GLY HA2/HA3 via
        // Locant::Alpha + BackboneRole::None per Markley.
        return backboneRole == BackboneRole::AlphaHydrogen
               || (element == Element::H && locant == Locant::Alpha && backboneRole == BackboneRole::None);
    }
    constexpr bool IsSidechainCarboxylateOxygen() const {
        return element == Element::O && planarGroup == PlanarGroupKind::Carboxylate;
    }
    constexpr bool IsSidechainAmideOxygen() const {
        return element == Element::O && planarGroup == PlanarGroupKind::SidechainAmide;
    }
    constexpr bool IsPolarH() const { return polarH != PolarHKind::NotPolar; }

    // Aromatic ring-facing hydrogen — the rediscover ring-current stratum.
    // Typed (no string dispatch): an H atom whose AMBER ff14SB atom type is
    // one of the aromatic-ring CH types — HA (H on sp2 aromatic ring), H4
    // (aromatic ring with one EWG, e.g. HIS), or H5 (aromatic ring with two
    // EWGs). These are exactly the protons that sit in the ring-current
    // shielding cone (PHE/TYR/TRP HD/HE/HZ/HH, HIS HD2/HE1). Backbone amide
    // HN (ff type H, NotPolar=false) and aliphatic H (HC/H1/H2/H3) are
    // excluded by construction. Mirrors the HN stratum's existing
    // IsBackboneAmideHydrogen() predicate in spirit. See DESIGN.md.
    constexpr bool IsAromaticRingHydrogen() const {
        return element == Element::H
               && (ffAtomType == QtFfAtomType::HA || ffAtomType == QtFfAtomType::H4
                   || ffAtomType == QtFfAtomType::H5);
    }

    constexpr bool IsInAnyRing() const {
        return ringPositionPrimary != RingPositionLabel::NotInRing || ringPositionSecondary != RingPositionLabel::NotInRing;
    }
};

}  // namespace h5reader::model
