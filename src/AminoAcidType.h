#pragma once
// Amino-acid template table and force-field residue-name helpers.

#include "Types.h"
#include <string>
#include <vector>

namespace nmr {

struct AminoAcidAtom {
    const char* name;
    Element     element;
    bool        is_backbone;
};

// A protonation variant (e.g., HID for histidine, ASH for aspartate).
//
// Variant order is load-bearing: stored variant_index values refer to
// positions in AminoAcidType::variants.
//
// Canonical indices (checked by ValidateVariantIndices):
//   HIS: 0=HID (delta), 1=HIE (epsilon), 2=HIP (doubly)
//   ASP: 0=ASH (protonated)
//   GLU: 0=GLH (protonated)
//   CYS: 0=CYX (disulfide), 1=CYM (deprotonated)
//   LYS: 0=LYN (deprotonated)
//   ARG: 0=ARN (deprotonated)
//   TYR: 0=TYM (deprotonated)
//
struct ProtonationVariant {
    const char* name;           // AMBER-standard name: "HID", "ASH", "CYX"
    const char* description;    // Human-readable: "Nd-protonated (delta)"
    int         formal_charge;  // Formal charge in this state
    const char* registry_key;   // NamingRegistry semantic key: "delta", "protonated"
                                // Used by ResolveForTool(canonical, context, key)
};

struct ChiAngleDef {
    const char* atoms[4];
};

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

// Aborts if variant ordering no longer matches the contract above.
void ValidateVariantIndices();

// Resolve a variant index from a force-field residue name.
//
// Returns -1 when the label names the canonical charged form or no
// known variant for this amino-acid type. Indices match the variant
// ordering contract:
//
//   HIS: HID/HSD = 0, HIE/HSE = 1, HIP/HSP = 2
//   ASP: ASH/ASPP = 0
//   GLU: GLH/GLUP = 0
//   CYS: CYX/CYS2 = 0, CYM = 1
//   LYS: LYN = 0
//   ARG: ARN = 0
//   TYR: TYM = 0
//
// GROMACS FF-port names such as HISH/HISE/HISD are intentionally not
// handled here; callers pass the canonical AMBER/CHARMM residue name.
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
