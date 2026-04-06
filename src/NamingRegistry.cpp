#include "NamingRegistry.h"
#include <algorithm>
#include <cctype>

namespace nmr {

// ============================================================================
// Helpers
// ============================================================================

static std::string ToUpper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    return s;
}

static std::string Trim(std::string s) {
    s.erase(0, s.find_first_not_of(' '));
    s.erase(s.find_last_not_of(' ') + 1);
    return s;
}

bool NamingRegistry::ContextKey::operator<(const ContextKey& o) const {
    if (canonical != o.canonical) return canonical < o.canonical;
    if (context != o.context) return context < o.context;
    return variant < o.variant;
}

bool NamingRegistry::AtomNameKey::operator<(const AtomNameKey& o) const {
    if (atom_name != o.atom_name) return atom_name < o.atom_name;
    if (residue_name != o.residue_name) return residue_name < o.residue_name;
    if (from_context != o.from_context) return from_context < o.from_context;
    return to_context < o.to_context;
}


// ============================================================================
// Construction: populate all naming tables
// ============================================================================

NamingRegistry::NamingRegistry() {
    InitialiseStandardResidues();
    InitialiseAmberContext();
    InitialiseCharmmContext();
    InitialiseAtomNameRules();
}


void NamingRegistry::InitialiseStandardResidues() {
    standard_residues_ = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    };

    // Every standard name maps to itself
    for (const auto& name : standard_residues_)
        to_canonical_[name] = name;

    // Variant -> canonical mappings (accumulated from all tool universes)
    to_canonical_["HID"] = "HIS";  to_canonical_["HIE"] = "HIS";
    to_canonical_["HIP"] = "HIS";
    to_canonical_["HSD"] = "HIS";  to_canonical_["HSE"] = "HIS";
    to_canonical_["HSP"] = "HIS";

    to_canonical_["ASH"] = "ASP";
    to_canonical_["ASPP"] = "ASP";

    to_canonical_["GLH"] = "GLU";
    to_canonical_["GLUP"] = "GLU";

    to_canonical_["CYX"] = "CYS";  to_canonical_["CYM"] = "CYS";
    to_canonical_["CYS2"] = "CYS";

    to_canonical_["LYN"] = "LYS";
    to_canonical_["ARN"] = "ARG";
    to_canonical_["TYM"] = "TYR";
    to_canonical_["MSE"] = "MET";
}


void NamingRegistry::InitialiseAmberContext() {
    auto add = [&](const std::string& canonical, const std::string& variant,
                    const std::string& amber_name) {
        context_map_[{canonical, ToolContext::Amber, variant}] = amber_name;
        // OpenMM uses the same naming as AMBER
        context_map_[{canonical, ToolContext::OpenMM, variant}] = amber_name;
    };

    // HIS tautomers
    add("HIS", "delta",   "HID");
    add("HIS", "epsilon", "HIE");
    add("HIS", "doubly",  "HIP");
    add("HIS", "",        "HIS");

    // Protonated acids
    add("ASP", "protonated", "ASH");
    add("GLU", "protonated", "GLH");

    // Cysteine states
    add("CYS", "disulfide",    "CYX");
    add("CYS", "deprotonated", "CYM");

    // Neutral lysine
    add("LYS", "deprotonated", "LYN");

    // Deprotonated arginine (extremely rare, pKa ~12.5)
    add("ARG", "deprotonated", "ARN");

    // Deprotonated tyrosine
    add("TYR", "deprotonated", "TYM");
}


void NamingRegistry::InitialiseCharmmContext() {
    auto add = [&](const std::string& canonical, const std::string& variant,
                    const std::string& charmm_name) {
        context_map_[{canonical, ToolContext::Charmm, variant}] = charmm_name;
    };

    // HIS tautomers (CHARMM uses HSD/HSE/HSP)
    add("HIS", "delta",   "HSD");
    add("HIS", "epsilon", "HSE");
    add("HIS", "doubly",  "HSP");
    add("HIS", "",        "HIS");

    // Protonated acids (CHARMM uses ASPP/GLUP)
    add("ASP", "protonated", "ASPP");
    add("GLU", "protonated", "GLUP");

    // Cysteine disulfide (CHARMM uses CYS2)
    add("CYS", "disulfide",    "CYS2");
    add("CYS", "deprotonated", "CYM");
}


void NamingRegistry::InitialiseAtomNameRules() {
    // Backbone amide hydrogen: H (PDB/AMBER) <-> HN (CHARMM)
    AddAtomNameRule("H",  "HN", "*", ToolContext::Standard, ToolContext::Charmm);
    AddAtomNameRule("HN", "H",  "*", ToolContext::Charmm, ToolContext::Standard);

    // Beta methylenes (all residues: HB2,HB3 <-> HB1,HB2)
    AddAtomNameRule("HB2", "HB1", "*", ToolContext::Standard, ToolContext::Charmm);
    AddAtomNameRule("HB3", "HB2", "*", ToolContext::Standard, ToolContext::Charmm);
    AddAtomNameRule("HB1", "HB2", "*", ToolContext::Charmm, ToolContext::Standard);
    AddAtomNameRule("HB2", "HB3", "*", ToolContext::Charmm, ToolContext::Standard);

    // Gamma methylenes (all residues: HG2,HG3 <-> HG1,HG2)
    AddAtomNameRule("HG2", "HG1", "*", ToolContext::Standard, ToolContext::Charmm);
    AddAtomNameRule("HG3", "HG2", "*", ToolContext::Standard, ToolContext::Charmm);
    AddAtomNameRule("HG1", "HG2", "*", ToolContext::Charmm, ToolContext::Standard);
    AddAtomNameRule("HG2", "HG3", "*", ToolContext::Charmm, ToolContext::Standard);
}


void NamingRegistry::AddAtomNameRule(const std::string& from_name,
                                      const std::string& to_name,
                                      const std::string& residue_name,
                                      ToolContext from_context,
                                      ToolContext to_context) {
    atom_name_map_[{from_name, residue_name, from_context, to_context}] = to_name;
}


// ============================================================================
// Residue name operations
// ============================================================================

bool NamingRegistry::IsKnownResidueName(const std::string& name) const {
    std::string upper = ToUpper(Trim(name));
    return to_canonical_.count(upper) > 0;
}

std::string NamingRegistry::ToCanonical(const std::string& name) const {
    std::string upper = ToUpper(Trim(name));
    auto it = to_canonical_.find(upper);
    if (it != to_canonical_.end())
        return it->second;
    return "";  // unknown -- return empty, caller checks
}

std::string NamingRegistry::ResolveForTool(const std::string& canonical,
                                            ToolContext context,
                                            const std::string& variant) const {
    // Exact (canonical, context, variant) match
    auto it = context_map_.find({canonical, context, variant});
    if (it != context_map_.end())
        return it->second;

    // Fall back to default variant
    if (!variant.empty()) {
        it = context_map_.find({canonical, context, ""});
        if (it != context_map_.end())
            return it->second;
    }

    // No context-specific mapping: return canonical name unchanged
    return canonical;
}


// ============================================================================
// Atom name operations
// ============================================================================

std::string NamingRegistry::TranslateAtomName(const std::string& atom_name,
                                               const std::string& residue_name,
                                               ToolContext from_context,
                                               ToolContext to_context) const {
    if (from_context == to_context) return atom_name;

    // Try residue-specific rule first
    auto it = atom_name_map_.find({atom_name, residue_name, from_context, to_context});
    if (it != atom_name_map_.end())
        return it->second;

    // Try wildcard residue
    it = atom_name_map_.find({atom_name, "*", from_context, to_context});
    if (it != atom_name_map_.end())
        return it->second;

    // No rule: return unchanged
    return atom_name;
}


// ============================================================================
// Global singleton
// ============================================================================

NamingRegistry& GlobalNamingRegistry() {
    static NamingRegistry instance;
    return instance;
}

}  // namespace nmr
