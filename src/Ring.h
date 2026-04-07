#pragma once
//
// Ring type class hierarchy.
//
// Ring types ARE classes with physics properties baked in.
// Each type provides const properties derived from its identity.
// Calculator code is ring-type-agnostic: ring.Intensity(), ring.JBLobeOffset().
//
// 8 types in 3 size categories:
//   SixMemberedRing: PheBenzeneRing, TyrPhenolRing, TrpBenzeneRing
//   FiveMemberedRing: TrpPyrroleRing, HisImidazoleRing, HidImidazoleRing, HieImidazoleRing
//   FusedRing: IndolePerimeterRing (TRP 9-atom)
//

#include "Types.h"
#include "CalculatorConfig.h"
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <memory>

namespace nmr {

// ============================================================================
// Ring::Geometry -- conformation-dependent, computed by GeometryResult
// ============================================================================

struct RingGeometry {
    Vec3              center = Vec3::Zero();
    Vec3              normal = Vec3::Zero();
    double            radius = 0.0;
    std::vector<Vec3> vertices;
};

// ============================================================================
// Accumulated ring properties (set by ConformationResult post-pass updates)
// ============================================================================

struct RingAccumulated {
    Vec3 total_B_at_center = Vec3::Zero();
    double intensity_used = 0.0;
    double total_G_T0_diagnostic = 0.0;
    std::map<size_t, Vec3> mutual_B_from;
};

// ============================================================================
// Ring (base class)
// ============================================================================

class Ring {
public:
    // Structural identity (topology, set at construction)
    std::vector<size_t> atom_indices;
    RingTypeIndex       type_index = RingTypeIndex::PheBenzene;
    size_t              parent_residue_index = 0;
    int                 parent_residue_number = 0;
    size_t              fused_partner_index = SIZE_MAX;

    virtual ~Ring() = default;

    // Virtual const properties (overridden by each type class)
    virtual double Intensity() const = 0;
    virtual double LiteratureIntensity() const = 0;
    virtual double JBLobeOffset() const = 0;
    virtual int NitrogenCount() const = 0;
    virtual RingAromaticity Aromaticity() const = 0;
    virtual int RingSizeValue() const = 0;
    virtual const char* TypeName() const = 0;

    // Non-virtual queries
    bool IsFused() const { return fused_partner_index != SIZE_MAX; }
    int TypeIndexAsInt() const { return static_cast<int>(type_index); }

    // Compute geometry from positions (SVD normal)
    RingGeometry ComputeGeometry(const std::vector<Vec3>& positions) const;

    // Accumulated properties (set during extraction passes)
    RingAccumulated accumulated;
};


// ============================================================================
// Six-membered rings
// ============================================================================

class SixMemberedRing : public Ring {
public:
    int RingSizeValue() const override { return 6; }
};

class PheBenzeneRing : public SixMemberedRing {
public:
    PheBenzeneRing() { type_index = RingTypeIndex::PheBenzene; }
    double Intensity() const override { return CalculatorConfig::Get("phe_benzene_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -12.0; }
    double JBLobeOffset() const override { return CalculatorConfig::Get("phe_benzene_jb_lobe_offset"); }
    int NitrogenCount() const override { return 0; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    const char* TypeName() const override { return "PHE"; }
};

class TyrPhenolRing : public SixMemberedRing {
public:
    TyrPhenolRing() { type_index = RingTypeIndex::TyrPhenol; }
    double Intensity() const override { return CalculatorConfig::Get("tyr_phenol_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -11.28; }
    double JBLobeOffset() const override { return CalculatorConfig::Get("tyr_phenol_jb_lobe_offset"); }
    int NitrogenCount() const override { return 0; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    const char* TypeName() const override { return "TYR"; }
};

class TrpBenzeneRing : public SixMemberedRing {
public:
    TrpBenzeneRing() { type_index = RingTypeIndex::TrpBenzene; }
    double Intensity() const override { return CalculatorConfig::Get("trp_benzene_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -12.48; }
    double JBLobeOffset() const override { return CalculatorConfig::Get("trp_benzene_jb_lobe_offset"); }
    int NitrogenCount() const override { return 0; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    const char* TypeName() const override { return "TRP6"; }
};


// ============================================================================
// Five-membered rings
// ============================================================================

class FiveMemberedRing : public Ring {
public:
    int RingSizeValue() const override { return 5; }
};

class TrpPyrroleRing : public FiveMemberedRing {
public:
    TrpPyrroleRing() { type_index = RingTypeIndex::TrpPyrrole; }
    double Intensity() const override { return CalculatorConfig::Get("trp_pyrrole_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -6.72; }
    double JBLobeOffset() const override { return CalculatorConfig::Get("trp_pyrrole_jb_lobe_offset"); }
    int NitrogenCount() const override { return 1; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Reduced; }
    const char* TypeName() const override { return "TRP5"; }
};

class HisImidazoleRing : public FiveMemberedRing {
public:
    HisImidazoleRing() { type_index = RingTypeIndex::HisImidazole; }
    double Intensity() const override { return CalculatorConfig::Get("his_imidazole_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -5.16; }
    double JBLobeOffset() const override { return CalculatorConfig::Get("his_imidazole_jb_lobe_offset"); }
    int NitrogenCount() const override { return 2; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Weak; }
    const char* TypeName() const override { return "HIS"; }
};

class HidImidazoleRing : public FiveMemberedRing {
public:
    HidImidazoleRing() { type_index = RingTypeIndex::HidImidazole; }
    double Intensity() const override { return CalculatorConfig::Get("hid_imidazole_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -5.16; }
    double JBLobeOffset() const override { return CalculatorConfig::Get("hid_imidazole_jb_lobe_offset"); }
    int NitrogenCount() const override { return 2; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Weak; }
    const char* TypeName() const override { return "HID"; }
};

class HieImidazoleRing : public FiveMemberedRing {
public:
    HieImidazoleRing() { type_index = RingTypeIndex::HieImidazole; }
    double Intensity() const override { return CalculatorConfig::Get("hie_imidazole_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -5.16; }
    double JBLobeOffset() const override { return CalculatorConfig::Get("hie_imidazole_jb_lobe_offset"); }
    int NitrogenCount() const override { return 2; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Weak; }
    const char* TypeName() const override { return "HIE"; }
};


// ============================================================================
// Fused ring (TRP 9-atom indole perimeter)
// ============================================================================

class FusedRing : public Ring {};

class IndolePerimeterRing : public FusedRing {
public:
    IndolePerimeterRing() { type_index = RingTypeIndex::TrpPerimeter; }
    double Intensity() const override { return CalculatorConfig::Get("trp_indole_perimeter_ring_current_intensity"); }
    double LiteratureIntensity() const override { return -19.2; }
    double JBLobeOffset() const override { return CalculatorConfig::Get("trp_indole_perimeter_jb_lobe_offset"); }
    int NitrogenCount() const override { return 1; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    int RingSizeValue() const override { return 9; }
    const char* TypeName() const override { return "TRP9"; }
};


// ============================================================================
// Factory: create a Ring subclass from a RingTypeIndex
// ============================================================================

std::unique_ptr<Ring> CreateRing(RingTypeIndex type);

}  // namespace nmr
