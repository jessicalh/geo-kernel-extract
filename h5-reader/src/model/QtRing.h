// QtRing — class hierarchy for the 9 ring types (8 aromatic + 1 saturated).
//
// Mirrors the library's nmr::Ring hierarchy: ring types ARE CLASSES,
// not enum values looked up in a table. Each subclass carries its
// physics (intensity, JB lobe offset, nitrogen count, aromaticity
// class, ring size) as virtual methods. Calculator and renderer code
// is ring-type-agnostic: ring->Intensity(), ring->NitrogenCount().
// No switch statements on ring_type outside the CreateQtRing() factory.
//
// Library mirror: src/Ring.h. Ordinal-compatible with nmr::Ring's
// type_index; the 9 types match nmr::RingTypeIndex 1:1, including
// the saturated ProPyrrolidine (RingTypeIndex 8 / RingKind::Saturated).
//
// String policy: the only string surface is `TypeName()` returning a
// `const char*` literal ("PHE", "TRP9", etc.) — display only. No
// physics or rendering code branches on the string. The discipline
// matches the library's Ring::TypeName().

#pragma once

#include "Types.h"

#include <cstdint>
#include <memory>
#include <vector>

namespace h5reader::model {

// ============================================================================
// RingGeometry — per-frame geometric state, computed at frame access
// from ring vertex positions. The current trajectory.h5 format does
// not emit per-frame ring geometry; fit it on demand.
// ============================================================================

struct RingGeometry {
    Vec3 center = Vec3::Zero();
    Vec3 normal = Vec3::Zero();
    double radius = 0.0;
};

// An orthonormal basis aligned to a ring normal. Used by overlays
// that sample in ring-local coordinates (B-field streamline grid,
// ring polygon vertices, butterfly grid orientation).
struct RingOrthoBasis {
    Vec3 u;
    Vec3 v;
    Vec3 n;
};

inline RingOrthoBasis OrthoBasisFromNormal(const Vec3& normalLike) {
    RingOrthoBasis b;
    b.n = normalLike.normalized();
    const Vec3 arbitrary = std::abs(b.n.x()) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
    b.u = b.n.cross(arbitrary).normalized();
    b.v = b.n.cross(b.u);
    return b;
}


// ============================================================================
// Abstract base — mirrors nmr::Ring (src/Ring.h:50).
// ============================================================================

class QtRing {
public:
    virtual ~QtRing() = default;

    // ----- Structural identity (decoded at load, const thereafter) -----
    int32_t ringId = -1;           // absolute row in rings.npy
    int32_t nativeAxisIndex = -1;  // aromatic-only OR saturated-only index
    RingKind ringKind = RingKind::Aromatic;
    int32_t parentResidueIndex = -1;
    int32_t parentResidueNumber = 0;   // PDB numbering, display
    int32_t fusedPartnerRingId = -1;   // -1 if not fused
    std::vector<int32_t> atomIndices;  // canonical-walk order

    // ----- Virtual physics properties (per ring type) -----
    virtual RingTypeIndex TypeIndex() const = 0;
    virtual double LiteratureIntensity() const = 0;  // nA/T, Giessner-Prettre 1969
    virtual double JohnsonBoveyLobeOffset() const = 0;         // Angstroms
    virtual int NitrogenCount() const = 0;
    virtual RingAromaticity Aromaticity() const = 0;
    virtual int RingSizeValue() const = 0;
    virtual const char* TypeName() const = 0;  // "PHE", "HIE", "TRP9", ...

    // ----- Non-virtual queries -----
    bool IsFused() const { return fusedPartnerRingId >= 0; }
    bool IsAromatic() const { return ringKind == RingKind::Aromatic; }
    int TypeIndexAsInt() const { return static_cast<int>(TypeIndex()); }
};


// ============================================================================
// Ring size categories
// ============================================================================

class QtSixMemberedRing : public QtRing {
public:
    int RingSizeValue() const override { return 6; }
};

class QtFiveMemberedRing : public QtRing {
public:
    int RingSizeValue() const override { return 5; }
};

class QtFusedRing : public QtRing {};


// ============================================================================
// Six-membered aromatic rings
// ============================================================================

class QtPheBenzeneRing final : public QtSixMemberedRing {
public:
    RingTypeIndex TypeIndex() const override { return RingTypeIndex::PheBenzene; }
    double LiteratureIntensity() const override { return -12.0; }
    double JohnsonBoveyLobeOffset() const override { return 0.64; }
    int NitrogenCount() const override { return 0; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    const char* TypeName() const override { return "PHE"; }
};

class QtTyrPhenolRing final : public QtSixMemberedRing {
public:
    RingTypeIndex TypeIndex() const override { return RingTypeIndex::TyrPhenol; }
    double LiteratureIntensity() const override { return -11.28; }
    double JohnsonBoveyLobeOffset() const override { return 0.64; }
    int NitrogenCount() const override { return 0; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    const char* TypeName() const override { return "TYR"; }
};

class QtTrpBenzeneRing final : public QtSixMemberedRing {
public:
    RingTypeIndex TypeIndex() const override { return RingTypeIndex::TrpBenzene; }
    double LiteratureIntensity() const override { return -12.48; }
    double JohnsonBoveyLobeOffset() const override { return 0.64; }
    int NitrogenCount() const override { return 0; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    const char* TypeName() const override { return "TRP6"; }
};


// ============================================================================
// Five-membered aromatic rings
// ============================================================================

class QtTrpPyrroleRing final : public QtFiveMemberedRing {
public:
    RingTypeIndex TypeIndex() const override { return RingTypeIndex::TrpPyrrole; }
    double LiteratureIntensity() const override { return -6.72; }
    double JohnsonBoveyLobeOffset() const override { return 0.52; }
    int NitrogenCount() const override { return 1; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Reduced; }
    const char* TypeName() const override { return "TRP5"; }
};

// HIS — unspecified tautomer (used when protonation state is unknown).
class QtHisImidazoleRing final : public QtFiveMemberedRing {
public:
    RingTypeIndex TypeIndex() const override { return RingTypeIndex::HisImidazole; }
    double LiteratureIntensity() const override { return -5.16; }
    double JohnsonBoveyLobeOffset() const override { return 0.50; }
    int NitrogenCount() const override { return 2; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Weak; }
    const char* TypeName() const override { return "HIS"; }
};

// HID — N-delta protonated tautomer.
class QtHidImidazoleRing final : public QtFiveMemberedRing {
public:
    RingTypeIndex TypeIndex() const override { return RingTypeIndex::HidImidazole; }
    double LiteratureIntensity() const override { return -5.16; }
    double JohnsonBoveyLobeOffset() const override { return 0.50; }
    int NitrogenCount() const override { return 2; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Weak; }
    const char* TypeName() const override { return "HID"; }
};

// HIE — N-epsilon protonated tautomer.
class QtHieImidazoleRing final : public QtFiveMemberedRing {
public:
    RingTypeIndex TypeIndex() const override { return RingTypeIndex::HieImidazole; }
    double LiteratureIntensity() const override { return -5.16; }
    double JohnsonBoveyLobeOffset() const override { return 0.50; }
    int NitrogenCount() const override { return 2; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Weak; }
    const char* TypeName() const override { return "HIE"; }
};


// ============================================================================
// Fused ring — TRP 9-atom indole perimeter
// ============================================================================

class QtIndolePerimeterRing final : public QtFusedRing {
public:
    RingTypeIndex TypeIndex() const override { return RingTypeIndex::TrpPerimeter; }
    double LiteratureIntensity() const override { return -19.2; }
    double JohnsonBoveyLobeOffset() const override { return 0.60; }
    int NitrogenCount() const override { return 1; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    int RingSizeValue() const override { return 9; }
    const char* TypeName() const override { return "TRP9"; }
};


// ============================================================================
// Saturated five-membered ring — Proline pyrrolidine
//
// Mirrors nmr::ProPyrrolidineRing (src/Ring.h:186). Aromaticity is
// None — the ring is not pi-conjugated, so ring-current intensity
// is identically zero. Joule & Mills 2010 ch. 7 (saturated
// heterocycles) is the chemistry citation. The literal 0.0 values
// for Intensity / LiteratureIntensity / JohnsonBoveyLobeOffset are therefore
// physics, not calibration parameters, and are not surfaced through
// any CalculatorConfig-equivalent.
// ============================================================================

class QtProPyrrolidineRing final : public QtFiveMemberedRing {
public:
    QtProPyrrolidineRing() { ringKind = RingKind::Saturated; }

    RingTypeIndex TypeIndex() const override { return RingTypeIndex::ProPyrrolidine; }
    double LiteratureIntensity() const override { return 0.0; }
    double JohnsonBoveyLobeOffset() const override { return 0.0; }
    int NitrogenCount() const override { return 1; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::None; }
    const char* TypeName() const override { return "PRO"; }
};


// ============================================================================
// Factory — the ONE place RingTypeIndex → concrete class is decided.
//
// Returns nullptr for unknown ordinal; caller logs via ErrorBus.
// See io/QtTopologySidecar.cpp.
// ============================================================================

std::unique_ptr<QtRing> CreateQtRing(RingTypeIndex idx);

}  // namespace h5reader::model
