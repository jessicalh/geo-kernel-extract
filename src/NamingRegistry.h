#pragma once
//
// NamingRegistry: the single authority for residue and atom name translation.
//
// Every tool in the pipeline has its own naming universe:
//   PDB/Standard: HIS, CYS, H (backbone amide)
//   AMBER:        HID/HIE/HIP, CYX/CYM, H
//   CHARMM36m:    HSD/HSE/HSP, CYS2/CYM, HN
//   ORCA:         standard PDB names (numbered atoms)
//   OpenMM:       AMBER-like (HIE, HID, CYX)
//
// This registry translates between them. It is the GATEKEEPER:
// unknown names are rejected. No guessing. No silent pass-through.
//
// After translation, typed objects carry their properties. The registry
// exists so that the rest of the system NEVER needs strings for identity.
//

#include "Types.h"
#include <string>
#include <map>
#include <set>
#include <vector>

namespace nmr {

// The tools we talk to. Each has its own naming conventions.
enum class ToolContext {
    Standard,    // PDB / IUPAC canonical names
    Amber,       // AMBER force fields (ff14SB, ff19SB): HID/HIE/HIP, CYX, ASH, GLH
    Charmm,      // CHARMM36m: HSD/HSE/HSP, CYS2, ASPP, GLUP, HN for backbone amide
    Orca,        // ORCA DFT (standard PDB names, numbered atoms)
    OpenMM,      // OpenMM (AMBER-like: HIE, HID, CYX)
};

inline const char* ToolContextName(ToolContext ctx) {
    switch (ctx) {
        case ToolContext::Standard: return "Standard";
        case ToolContext::Amber:    return "Amber";
        case ToolContext::Charmm:   return "CHARMM";
        case ToolContext::Orca:     return "ORCA";
        case ToolContext::OpenMM:   return "OpenMM";
    }
    return "Unknown";
}


class NamingRegistry {
public:
    NamingRegistry();

    // ================================================================
    // Residue name translation
    // ================================================================

    // Is this a known residue name (in any context)?
    bool IsKnownResidueName(const std::string& name) const;

    // Map any variant to canonical (standard 20 three-letter codes).
    // HID -> HIS, ASH -> ASP, CYX -> CYS, MSE -> MET, etc.
    // Returns empty string if unknown.
    std::string ToCanonical(const std::string& name) const;

    // Map canonical to tool-specific name.
    // For titratable residues, caller specifies the variant:
    //   ResolveForTool("HIS", ToolContext::Charmm, "delta")   -> "HSD"
    //   ResolveForTool("HIS", ToolContext::Charmm, "epsilon") -> "HSE"
    //   ResolveForTool("HIS", ToolContext::Charmm, "doubly")  -> "HSP"
    // For non-titratable residues, variant is ignored:
    //   ResolveForTool("ALA", ToolContext::Charmm) -> "ALA"
    // Returns canonical name unchanged if no mapping exists.
    std::string ResolveForTool(const std::string& canonical,
                                ToolContext context,
                                const std::string& variant = "") const;

    // ================================================================
    // Atom name translation
    // ================================================================

    // Translate an atom name between tools.
    // Backbone amide: H (PDB/AMBER) <-> HN (CHARMM)
    std::string TranslateAtomName(const std::string& atom_name,
                                   const std::string& residue_name,
                                   ToolContext from_context,
                                   ToolContext to_context) const;

private:
    // The 20 standard amino acids
    std::set<std::string> standard_residues_;

    // Any variant name -> canonical name
    std::map<std::string, std::string> to_canonical_;

    // (canonical, context, variant) -> tool-specific name
    struct ContextKey {
        std::string canonical;
        ToolContext context;
        std::string variant;
        bool operator<(const ContextKey& o) const;
    };
    std::map<ContextKey, std::string> context_map_;

    // Atom name translation: (atom_name, residue_name, from, to) -> translated
    struct AtomNameKey {
        std::string atom_name;
        std::string residue_name;  // "*" matches all
        ToolContext from_context;
        ToolContext to_context;
        bool operator<(const AtomNameKey& o) const;
    };
    std::map<AtomNameKey, std::string> atom_name_map_;

    void InitialiseStandardResidues();
    void InitialiseAmberContext();
    void InitialiseCharmmContext();
    void InitialiseAtomNameRules();

    void AddAtomNameRule(const std::string& from_name,
                          const std::string& to_name,
                          const std::string& residue_name,
                          ToolContext from_context,
                          ToolContext to_context);
};

// Global singleton (C++11 Meyers singleton, thread-safe).
NamingRegistry& GlobalNamingRegistry();

}  // namespace nmr
