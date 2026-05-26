#pragma once
//
// Ring type class hierarchy: each subclass exposes its type-specific
// physical parameters (ring-current intensity, Johnson-Bovey lobe
// offset, nitrogen count, ...) through virtual methods.
//

#include "Types.h"
#include "CalculatorConfig.h"
#include <Eigen/Dense>
#include <vector>
#include <memory>

namespace nmr {

struct RingGeometry {
    Vec3              center = Vec3::Zero();
    Vec3              normal = Vec3::Zero();
    double            radius = 0.0;
    std::vector<Vec3> vertices;
};

class Ring {
public:
    // Structural identity / topology.
    std::vector<size_t> atom_indices;
    RingTypeIndex       type_index = RingTypeIndex::PheBenzene;
    size_t              parent_residue_index = 0;
    int                 parent_residue_number = 0;
    size_t              fused_partner_index = SIZE_MAX;

    virtual ~Ring() = default;

    // Ring-type physical properties (one override per type class):
    //   Intensity / LiteratureIntensity : ring-current strength, nA/T
    //                                      (negative = diamagnetic)
    //   JohnsonBoveyLobeOffset           : JB two-loop z-offset, Å
    //   NitrogenCount                    : ring N atoms
    virtual double Intensity() const = 0;
    virtual double LiteratureIntensity() const = 0;
    virtual double JohnsonBoveyLobeOffset() const = 0;
    virtual int NitrogenCount() const = 0;
    virtual RingAromaticity Aromaticity() const = 0;
    virtual int RingAtomCount() const = 0;
    virtual const char* TypeName() const = 0;

    bool IsFused() const { return fused_partner_index != SIZE_MAX; }
    int TypeIndexAsInt() const { return static_cast<int>(type_index); }

    RingGeometry ComputeGeometry(const std::vector<Vec3>& positions) const;
};


class SixMemberedRing : public Ring {
public:
    int RingAtomCount() const override { return 6; }
};

class PheBenzeneRing : public SixMemberedRing {
public:
    PheBenzeneRing() { type_index = RingTypeIndex::PheBenzene; }
    double Intensity() const override { return CalculatorConfig::Get("phe_benzene_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -12.0; }
    double JohnsonBoveyLobeOffset() const override { return CalculatorConfig::Get("phe_benzene_jb_lobe_offset"); }
    int NitrogenCount() const override { return 0; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    const char* TypeName() const override { return "PHE"; }
};

class TyrPhenolRing : public SixMemberedRing {
public:
    TyrPhenolRing() { type_index = RingTypeIndex::TyrPhenol; }
    double Intensity() const override { return CalculatorConfig::Get("tyr_phenol_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -11.28; }
    double JohnsonBoveyLobeOffset() const override { return CalculatorConfig::Get("tyr_phenol_jb_lobe_offset"); }
    int NitrogenCount() const override { return 0; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    const char* TypeName() const override { return "TYR"; }
};

class TrpBenzeneRing : public SixMemberedRing {
public:
    TrpBenzeneRing() { type_index = RingTypeIndex::TrpBenzene; }
    double Intensity() const override { return CalculatorConfig::Get("trp_benzene_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -12.48; }
    double JohnsonBoveyLobeOffset() const override { return CalculatorConfig::Get("trp_benzene_jb_lobe_offset"); }
    int NitrogenCount() const override { return 0; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    const char* TypeName() const override { return "TRP6"; }
};


class FiveMemberedRing : public Ring {
public:
    int RingAtomCount() const override { return 5; }
};

class TrpPyrroleRing : public FiveMemberedRing {
public:
    TrpPyrroleRing() { type_index = RingTypeIndex::TrpPyrrole; }
    double Intensity() const override { return CalculatorConfig::Get("trp_pyrrole_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -6.72; }
    double JohnsonBoveyLobeOffset() const override { return CalculatorConfig::Get("trp_pyrrole_jb_lobe_offset"); }
    int NitrogenCount() const override { return 1; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Reduced; }
    const char* TypeName() const override { return "TRP5"; }
};

class HisImidazoleRing : public FiveMemberedRing {
public:
    HisImidazoleRing() { type_index = RingTypeIndex::HisImidazole; }
    double Intensity() const override { return CalculatorConfig::Get("his_imidazole_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -5.16; }
    double JohnsonBoveyLobeOffset() const override { return CalculatorConfig::Get("his_imidazole_jb_lobe_offset"); }
    int NitrogenCount() const override { return 2; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Weak; }
    const char* TypeName() const override { return "HIS"; }
};

class HidImidazoleRing : public FiveMemberedRing {
public:
    HidImidazoleRing() { type_index = RingTypeIndex::HidImidazole; }
    double Intensity() const override { return CalculatorConfig::Get("hid_imidazole_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -5.16; }
    double JohnsonBoveyLobeOffset() const override { return CalculatorConfig::Get("hid_imidazole_jb_lobe_offset"); }
    int NitrogenCount() const override { return 2; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Weak; }
    const char* TypeName() const override { return "HID"; }
};

class HieImidazoleRing : public FiveMemberedRing {
public:
    HieImidazoleRing() { type_index = RingTypeIndex::HieImidazole; }
    double Intensity() const override { return CalculatorConfig::Get("hie_imidazole_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -5.16; }
    double JohnsonBoveyLobeOffset() const override { return CalculatorConfig::Get("hie_imidazole_jb_lobe_offset"); }
    int NitrogenCount() const override { return 2; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Weak; }
    const char* TypeName() const override { return "HIE"; }
};

// Pro pyrrolidine: saturated 5-ring (N + Cα + Cβ + Cγ + Cδ).
// Aromaticity is None — the ring is not π-conjugated, so the
// ring-current intensity is identically zero. Joule & Mills 2010
// ch. 7 (saturated heterocycles) is the chemistry citation: ring
// current is a property of cyclic π conjugation, which pyrrolidine
// lacks. The literal 0.0 values for Intensity / LiteratureIntensity
// / JohnsonBoveyLobeOffset are therefore physics, not calibration parameters,
// and are not surfaced through CalculatorConfig.
class ProPyrrolidineRing : public FiveMemberedRing {
public:
    ProPyrrolidineRing() { type_index = RingTypeIndex::ProPyrrolidine; }
    double Intensity() const override { return 0.0; }
    double LiteratureIntensity() const override { return 0.0; }
    double JohnsonBoveyLobeOffset() const override { return 0.0; }
    int NitrogenCount() const override { return 1; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::None; }
    const char* TypeName() const override { return "PRO"; }
};


class FusedRing : public Ring {};

class IndolePerimeterRing : public FusedRing {
public:
    IndolePerimeterRing() { type_index = RingTypeIndex::TrpPerimeter; }
    double Intensity() const override { return CalculatorConfig::Get("trp_indole_perimeter_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -19.2; }
    double JohnsonBoveyLobeOffset() const override { return CalculatorConfig::Get("trp_indole_perimeter_jb_lobe_offset"); }
    int NitrogenCount() const override { return 1; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    int RingAtomCount() const override { return 9; }
    const char* TypeName() const override { return "TRP9"; }
};


std::unique_ptr<Ring> CreateRing(RingTypeIndex type);

}  // namespace nmr
