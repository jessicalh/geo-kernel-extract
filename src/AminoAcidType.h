#pragma once
//
// AminoAcidType: amino-acid template properties and variant helpers.
//
// Each table entry stores canonical atoms, chi angle definitions,
// titratable variants, and summary flags used by loaders/calculators.
//
// Usage:
//   const AminoAcidType& type = GetAminoAcidType(AminoAcid::PHE);
//   type.is_aromatic;           // true
//   type.chi_angle_count;       // 2
//

#include "Types.h"
#include <string>
#include <vector>

namespace nmr {

// An atom in the canonical amino acid template (PDB naming).
struct AminoAcidAtom {
    const char* name;
    Element     element;
    bool        is_backbone;
};

// A protonation variant (e.g., HID for histidine, ASH for aspartate).
//
// VARIANT INDEX CONTRACT: the position of each variant in
// AminoAcidType::variants is load-bearing. ProtonationDetectionResult,
// ProtonationState, and the AMBER charge/preparation paths all use
// variant_index to identify the protonation state. Reordering silently
// breaks protonation assignment consumers.
//
// Canonical indices (checked by ValidateVariantIndices):
//   HIS: 0=HID (delta), 1=HIE (epsilon), 2=HIP (doubly)
//   ASP: 0=ASH (protonated)
//   GLU: 0=GLH (protonated)
//   CYS: 0=CYX (disulfide), 1=CYM (deprotonated)
//   LYS: 0=LYN (deprotonated)
//   ARG: 0=ARN (deprotonated, pKa ~12.5, very rare)
//   TYR: 0=TYM (deprotonated)
//
struct ProtonationVariant {
    const char* name;           // AMBER-standard name: "HID", "ASH", "CYX"
    const char* description;    // Human-readable: "Nd-protonated (delta)"
    int         formal_charge;  // Formal charge in this state
    const char* registry_key;   // NamingRegistry semantic key: "delta", "protonated"
                                // Used by ResolveForTool(canonical, context, key)
};

// A chi angle definition: four atom names defining the dihedral.
struct ChiAngleDef {
    const char* atoms[4];
};

// Ring membership is supplied by RingTopology::ConstructFromSubstrate,
// which reads typed RingPosition slots from each atom's AtomSemanticTable.
// AminoAcidType keeps canonical atom names, chi definitions, and
// protonation variants.

class AminoAcidType {
public:
    AminoAcid       index;
    const char*     three_letter_code;
    char            one_letter_code;

    bool            is_aromatic;
    bool            is_titratable;
    bool            has_amide_H;
    int             chi_angle_count;
    int             charged_formal_charge;

    std::vector<AminoAcidAtom> atoms;
    std::vector<ChiAngleDef>   chi_angles;
    std::vector<ProtonationVariant> variants;

    bool HasAtom(const char* name) const {
        for (const auto& a : atoms)
            if (std::string(a.name) == name) return true;
        return false;
    }

    const AminoAcidAtom* FindAtom(const char* name) const {
        for (const auto& a : atoms)
            if (std::string(a.name) == name) return &a;
        return nullptr;
    }
};

const AminoAcidType& GetAminoAcidType(AminoAcid aa);
const std::vector<AminoAcidType>& AllAminoAcidTypes();
const AminoAcidType& AminoAcidTypeFromCode(const std::string& code);

// Validates that variant ordering matches the documented contract.
// Aborts on mismatch.
void ValidateVariantIndices();

// Resolve a variant index from a force-field residue name.
//
// Single shared helper: callers that have a force-field residue label
// (CHARMM HSD/HSE/HSP, AMBER HID/HIE/HIP/CYX/CYM/ASH/GLH/LYN/ARN/TYM,
// CHARMM-port alternates ASPP/GLUP/CYS2) map it to the canonical
// variant index for the given amino acid type. Returns -1 when the
// label names the canonical-charged-state form (no variant) or doesn't
// match any known variant.
//
// Indices match the AminoAcidType.h canonical contract checked by
// ValidateVariantIndices():
//
//   HIS: HID/HSD = 0, HIE/HSE = 1, HIP/HSP = 2
//   ASP: ASH/ASPP = 0
//   GLU: GLH/GLUP = 0
//   CYS: CYX/CYS2 = 0, CYM = 1
//   LYS: LYN = 0
//   ARG: ARN = 0
//   TYR: TYM = 0
//
// FF-port labels that GROMACS pdb2gmx writes back (HISH / HISE / HISD
// for amber14sb, etc.) are NOT handled here — those are resolved
// upstream by reading the topol.top rtp comment line. The
// GromacsToAmberReadbackBlock is a compiler-trace shape (see memory
// feedback_readback_block_is_a_compiler_trace). Callers should
// pass canonical AMBER/CHARMM names, not GROMACS FF-port labels.
int VariantIndexFromForceFieldName(AminoAcid type, const std::string& ff_name);

// Strip an N- or C- terminal prefix from a GROMACS rtp name, returning the
// base FF-port name. The prefix is recognised only if the remaining four
// characters are not themselves a known self-canonical FF-port name (CYS2,
// ASPP, GLUP, LYSN). Pass-through for non-prefixed names.
//
//   "VAL"   → "VAL"
//   "NVAL"  → "VAL"
//   "CALA"  → "ALA"
//   "NHIP"  → "HIP"
//   "CHIE"  → "HIE"
//   "NCYX"  → "CYX"
//   "CYS2"  → "CYS2"   (self-canonical 4-letter; not stripped)
//   "ASPP"  → "ASPP"   (self-canonical 4-letter; not stripped)
std::string BaseFfPortNameFromGromacsRtp(const std::string& rtp);

// Resolve a GROMACS rtp name (possibly N-/C-terminal-prefixed FF-port form)
// to the canonical AMBER/PDB 3-letter residue code.
//
//   "VAL"   → "VAL"
//   "NVAL"  → "VAL"
//   "HIP"   → "HIS"
//   "NHIP"  → "HIS"
//   "CHIE"  → "HIS"
//   "CYX"   → "CYS"
//   "CYS2"  → "CYS"
//   "ASH"   → "ASP"
//
// Returns empty string if the name can't be resolved to a known amino acid.
std::string CanonicalThreeLetterFromGromacsRtp(const std::string& rtp);

}  // namespace nmr
