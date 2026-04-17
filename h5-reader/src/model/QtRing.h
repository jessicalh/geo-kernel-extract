// QtRing — class hierarchy for the 8 aromatic ring types.
//
// Ring types are CLASSES, not enum values looked up in a table. Each
// subclass carries its physics (intensity, JB lobe offset, nitrogen
// count, aromaticity class, ring size) as virtual methods. Calculator
// and renderer code is ring-type-agnostic: ring->Intensity(),
// ring->NitrogenCount(). No switch statements on ring_type outside the
// CreateQtRing() factory.
//
// This mirrors nmr::Ring in src/Ring.h and is the canonical "objects
// answer questions about themselves" pattern for the reader.
//
// Instances are created by io::CreateQtRing(RingTypeIndex) at H5 load
// time, keyed on topology/ring_type. The factory is the ONE place the
// ordinal-to-class mapping lives.

#pragma once

#include "Types.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace h5reader::model {

// RingGeometry — per-frame geometric state, populated from
// ring_geometry/data (T, n_rings, 7) at the frame slab.
struct RingGeometry {
    Vec3   center  = Vec3::Zero();
    Vec3   normal  = Vec3::Zero();
    double radius  = 0.0;
};

// An orthonormal basis aligned to a ring normal. {u, v} span the ring
// plane; n is the normal itself. Used by overlays that sample in
// ring-local coordinates (B-field streamline grid, ring polygon
// vertices, butterfly grid orientation). Constructed via
// OrthoBasisFromNormal().
struct RingOrthoBasis {
    Vec3 u;
    Vec3 v;
    Vec3 n;
};

// Derive an orthonormal basis from an arbitrary-looking vector by
// picking a non-parallel reference axis, taking cross products. Same
// rule used throughout the library viewer and the ring-polygon /
// butterfly-field code here — centralised to avoid drift between the
// three call sites.
inline RingOrthoBasis OrthoBasisFromNormal(const Vec3& normalLike) {
    RingOrthoBasis b;
    b.n = normalLike.normalized();
    const Vec3 arbitrary =
        std::abs(b.n.x()) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
    b.u = b.n.cross(arbitrary).normalized();
    b.v = b.n.cross(b.u);
    return b;
}


// ---------------------------------------------------------------------------
// Abstract base
// ---------------------------------------------------------------------------

class QtRing {
public:
    virtual ~QtRing() = default;

    // Structural identity (decoded at load time, const thereafter).
    std::vector<size_t> atomIndices;         // indices into QtProtein atoms
    int                 parentResidueIndex = -1;
    int                 parentResidueNumber = 0;   // PDB numbering, display
    int                 fusedPartnerIndex = -1;    // into QtProtein.rings(), -1 if not fused

    // Virtual physics properties — each subclass answers for its type.
    virtual RingTypeIndex    TypeIndex()        const = 0;
    virtual double           LiteratureIntensity() const = 0;  // nA/T, Giessner-Prettre 1969
    virtual double           JBLobeOffset()     const = 0;     // Angstroms
    virtual int              NitrogenCount()    const = 0;
    virtual RingAromaticity  Aromaticity()      const = 0;
    virtual int              RingSizeValue()    const = 0;
    virtual const char*      TypeName()         const = 0;     // "PHE", "HIE", "TRP9"…

    // Non-virtual queries
    bool IsFused() const { return fusedPartnerIndex >= 0; }
    int  TypeIndexAsInt() const { return static_cast<int>(TypeIndex()); }
};


// ---------------------------------------------------------------------------
// Ring size categories
// ---------------------------------------------------------------------------

class QtSixMemberedRing : public QtRing {
public:
    int RingSizeValue() const override { return 6; }
};

class QtFiveMemberedRing : public QtRing {
public:
    int RingSizeValue() const override { return 5; }
};

class QtFusedRing : public QtRing {};


// ---------------------------------------------------------------------------
// Concrete ring types — 8 subclasses
//
// Intensities are Giessner-Prettre & Pullman 1969 literature values.
// The extractor's TOML-calibrated intensities are NOT used by the
// reader: ring_current kernel values per atom are already in the H5.
// These virtual methods exist for (a) volumetric BS/HM re-evaluation
// when the reader draws butterflies (where we sample the closed-form
// kernel at open-space grid points with the literature coefficient),
// and (b) atom-inspector display.
// ---------------------------------------------------------------------------

class QtPheBenzeneRing final : public QtSixMemberedRing {
public:
    RingTypeIndex    TypeIndex()           const override { return RingTypeIndex::PheBenzene; }
    double           LiteratureIntensity() const override { return -12.0; }
    double           JBLobeOffset()        const override { return 0.64; }
    int              NitrogenCount()       const override { return 0; }
    RingAromaticity  Aromaticity()         const override { return RingAromaticity::Full; }
    const char*      TypeName()            const override { return "PHE"; }
};

class QtTyrPhenolRing final : public QtSixMemberedRing {
public:
    RingTypeIndex    TypeIndex()           const override { return RingTypeIndex::TyrPhenol; }
    double           LiteratureIntensity() const override { return -11.28; }
    double           JBLobeOffset()        const override { return 0.64; }
    int              NitrogenCount()       const override { return 0; }
    RingAromaticity  Aromaticity()         const override { return RingAromaticity::Full; }
    const char*      TypeName()            const override { return "TYR"; }
};

class QtTrpBenzeneRing final : public QtSixMemberedRing {
public:
    RingTypeIndex    TypeIndex()           const override { return RingTypeIndex::TrpBenzene; }
    double           LiteratureIntensity() const override { return -12.48; }
    double           JBLobeOffset()        const override { return 0.64; }
    int              NitrogenCount()       const override { return 0; }
    RingAromaticity  Aromaticity()         const override { return RingAromaticity::Full; }
    const char*      TypeName()            const override { return "TRP6"; }
};

class QtTrpPyrroleRing final : public QtFiveMemberedRing {
public:
    RingTypeIndex    TypeIndex()           const override { return RingTypeIndex::TrpPyrrole; }
    double           LiteratureIntensity() const override { return -6.72; }
    double           JBLobeOffset()        const override { return 0.52; }
    int              NitrogenCount()       const override { return 1; }
    RingAromaticity  Aromaticity()         const override { return RingAromaticity::Reduced; }
    const char*      TypeName()            const override { return "TRP5"; }
};

// HIS — unspecified tautomer. Used when protonation state is unknown.
class QtHisImidazoleRing final : public QtFiveMemberedRing {
public:
    RingTypeIndex    TypeIndex()           const override { return RingTypeIndex::HisImidazole; }
    double           LiteratureIntensity() const override { return -5.16; }
    double           JBLobeOffset()        const override { return 0.50; }
    int              NitrogenCount()       const override { return 2; }
    RingAromaticity  Aromaticity()         const override { return RingAromaticity::Weak; }
    const char*      TypeName()            const override { return "HIS"; }
};

// HID — N-delta protonated tautomer.
class QtHidImidazoleRing final : public QtFiveMemberedRing {
public:
    RingTypeIndex    TypeIndex()           const override { return RingTypeIndex::HidImidazole; }
    double           LiteratureIntensity() const override { return -5.16; }
    double           JBLobeOffset()        const override { return 0.50; }
    int              NitrogenCount()       const override { return 2; }
    RingAromaticity  Aromaticity()         const override { return RingAromaticity::Weak; }
    const char*      TypeName()            const override { return "HID"; }
};

// HIE — N-epsilon protonated tautomer.
class QtHieImidazoleRing final : public QtFiveMemberedRing {
public:
    RingTypeIndex    TypeIndex()           const override { return RingTypeIndex::HieImidazole; }
    double           LiteratureIntensity() const override { return -5.16; }
    double           JBLobeOffset()        const override { return 0.50; }
    int              NitrogenCount()       const override { return 2; }
    RingAromaticity  Aromaticity()         const override { return RingAromaticity::Weak; }
    const char*      TypeName()            const override { return "HIE"; }
};

// TRP 9-atom indole perimeter — physical ring current path, sum of
// TRP5 and TRP6 intensities. Carries the actual current in MD.
class QtIndolePerimeterRing final : public QtFusedRing {
public:
    RingTypeIndex    TypeIndex()           const override { return RingTypeIndex::TrpPerimeter; }
    double           LiteratureIntensity() const override { return -19.2; }
    double           JBLobeOffset()        const override { return 0.60; }
    int              NitrogenCount()       const override { return 1; }
    RingAromaticity  Aromaticity()         const override { return RingAromaticity::Full; }
    int              RingSizeValue()       const override { return 9; }
    const char*      TypeName()            const override { return "TRP9"; }
};


// ---------------------------------------------------------------------------
// Factory — the ONE place RingTypeIndex → concrete class is decided.
//
// Returns nullptr for unknown ordinal; caller should ErrorBus::Report.
// See io/QtProteinLoader.cpp.
// ---------------------------------------------------------------------------

std::unique_ptr<QtRing> CreateQtRing(RingTypeIndex idx);

}  // namespace h5reader::model
