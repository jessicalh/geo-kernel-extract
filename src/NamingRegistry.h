#pragma once
// Residue-name translation and atom-name canonicalisation. Unknown
// residue names return an empty canonical name; atom-name failures abort.

#include "Types.h"
#include "SemanticEnums.h"
#include <cstdint>
#include <functional>
#include <memory>
#include <set>
#include <string>
#include <string_view>
#include <map>
#include <vector>

namespace nmr {

// External naming surfaces for residue-name translation. Atom-name
// rules use NamingSource instead.
enum class ToolContext {
    Standard = 0,    // PDB / IUPAC canonical names
    Amber    = 1,    // AMBER force fields (ff14SB, ff19SB): HID/HIE/HIP, CYX, ASH, GLH
    // Slot 2 is intentionally unused; reusing it could rebind any
    // persisted records carrying the old value.
    Orca     = 3,    // ORCA DFT (standard PDB names, numbered atoms)
    OpenMM   = 4,    // OpenMM (AMBER-like: HIE, HID, CYX)
};

inline const char* ToolContextName(ToolContext ctx) {
    switch (ctx) {
        case ToolContext::Standard: return "Standard";
        case ToolContext::Amber:    return "Amber";
        case ToolContext::Orca:     return "ORCA";
        case ToolContext::OpenMM:   return "OpenMM";
    }
    return "Unknown";
}


// Source tag for atom-name rules. It names the input convention a rule
// represents, not the output ToolContext.
enum class NamingSource : uint8_t {
    Unknown                  = 0,

    /// Project-internal canonical AMBER ff14SB target vocabulary.
    AmberFf14SBCanonical     = 1,

    /// AMBER pdb2gmx writes side-chain methylenes and methyl-bearing
    /// carbons in older AMBER RTP convention; the resulting atom names
    /// deviate from canonical AMBER ff14SB / IUPAC. Verified inside
    /// fleet_amber/{1P9J,1Z9B}/topol.top mol_X rtp blocks.
    Pdb2gmxAmberRtpDeviation = 2,

    /// cifpp's PDB parsing yields mostly canonical IUPAC/AMBER
    /// post-Markley names.
    CifppPdbInput            = 3,

    /// ORCA NMR output echoes the input atom names.
    OrcaEcho                 = 4,

    /// Markley et al. 1998 J. Biomol. NMR 12:1-23 nomenclature
    /// recommendations (the IUPAC-IUB 1969 update for proteins). The
    /// Greek-letter-locant + diastereotopic-numbering canon lives
    /// here when expressed as a runtime rule.
    Markley1998              = 6,

    /// BMRB nomenclature table at https://bmrb.io/ref_info.
    BmrbAtomNomTbl           = 7,

    /// IUPAC-IUB 1969 tentative rules.
    IupacIub1969             = 8,

    /// Project-internal synthesised conventions.
    ProjectSynthesis         = 9,
};

inline const char* NamingSourceName(NamingSource src) {
    switch (src) {
        case NamingSource::Unknown:                  return "Unknown";
        case NamingSource::AmberFf14SBCanonical:     return "AmberFf14SBCanonical";
        case NamingSource::Pdb2gmxAmberRtpDeviation: return "Pdb2gmxAmberRtpDeviation";
        case NamingSource::CifppPdbInput:            return "CifppPdbInput";
        case NamingSource::OrcaEcho:                 return "OrcaEcho";
        case NamingSource::Markley1998:              return "Markley1998";
        case NamingSource::BmrbAtomNomTbl:           return "BmrbAtomNomTbl";
        case NamingSource::IupacIub1969:             return "IupacIub1969";
        case NamingSource::ProjectSynthesis:         return "ProjectSynthesis";
    }
    return "Unknown";
}


struct NamingContext {
    NamingSource  source         = NamingSource::Unknown;
    std::string   input_name;
    AminoAcid     residue_type   = AminoAcid::Unknown;

    /// Protonation variant index into AminoAcidType::variants. -1 is
    /// the unresolved-load state; resolved variants use >= 0.
    int           variant_index  = -1;

    TerminalState terminal_state = TerminalState::Internal;

    /// Sibling names in input form, snapped before any atom in the
    /// residue is rewritten. Shift-pair rules inspect this original set.
    std::set<std::string> sibling_input_names;

    /// Parent heavy atom's input name (for hydrogen disambiguation).
    /// Empty for non-hydrogen atoms or when parent is unknown.
    std::string parent_input_name;

    /// Diagnostics-only fields (not consumed by rule predicates;
    /// included in fail-loud messages for triage).
    int         residue_sequence_number = 0;
    std::string chain_id;
};


struct NamingRule {
    NamingSource     source;
    std::string_view name;       ///< Stable identifier for diagnostics + tests.
    std::string_view rationale;  ///< One-line citation of the source decision.

    /// Predicate: does this rule apply to this context?
    std::function<bool(const NamingContext&)> applies;

    /// Output: given the context, what canonical form does this rule
    /// propose?
    std::function<std::string(const NamingContext&)> output;
};


struct NamingApplication {
    const NamingRule* rule;             ///< Which rule fired.
    std::string       proposed_output;  ///< What that rule proposes.
};


namespace test {
// Forward-declared test access helper. Lives in the tests/ tree;
// production code never instantiates it.
class NamingApplicatorTestAccess;
}  // namespace test

class NamingApplicator {
public:
    NamingApplicator();

    /// Single-atom canonicalisation. Builds the per-atom application
    /// map, calls the resolver, returns the canonical output. Aborts
    /// loudly on unresolvable cases (no rule + non-canonical input,
    /// or multi-rule combination unhandled by Resolve()).
    std::string Apply(const NamingContext& ctx) const;

    /// Whole-residue canonicalisation: snapshots sibling names ONCE
    /// for the residue, then iterates atoms. Necessary for shift-pair
    /// rules to read the original sibling set, not the partially-
    /// rewritten set.
    ///
    /// `input_names` parallel to `parent_input_names`; both parallel
    /// to the resulting output. Caller is responsible for assigning
    /// outputs onto Atom::pdb_atom_name.
    std::vector<std::string> ApplyResidue(
        const std::vector<std::string>& input_names,
        const std::vector<std::string>& parent_input_names,
        AminoAcid residue_type,
        int variant_index,
        TerminalState terminal_state,
        NamingSource source,
        int residue_sequence_number,
        std::string_view chain_id) const;

    /// Diagnostics-only accessor: how many rules are loaded?
    size_t RuleCount() const { return rules_.size(); }

    /// Diagnostics-only accessor: name of the i-th rule (for tests).
    std::string_view RuleNameAt(size_t i) const {
        return (i < rules_.size()) ? rules_[i].name : std::string_view{};
    }

    /// Optional debug-logging mode (off by default).
    void SetDebugLogging(bool enabled) { debug_logging_ = enabled; }

    /// Canonicality oracle: is `ctx.input_name` already a valid
    /// canonical AMBER ff14SB atom name for `ctx.residue_type` plus
    /// `ctx.variant_index`? Reads from the AminoAcidType chain table
    /// plus per-variant atom-presence rules. Used by Resolve() to
    /// detect idempotent (already-canonical) input. Public so
    /// property tests can exercise it directly.
    bool IsCanonical(const NamingContext& ctx) const;

private:
    // Test-only constructor with custom rules. Production code reaches
    // only the default rule set through GlobalNamingApplicator().
    friend class test::NamingApplicatorTestAccess;
    struct CustomRules { std::vector<NamingRule> rules; };
    explicit NamingApplicator(CustomRules custom);

    /// Step 1: iterate rules, return all applications.
    std::vector<NamingApplication> Collect(const NamingContext& ctx) const;

    /// Step 2: resolve the per-atom map to a single output.
    std::string Resolve(const std::vector<NamingApplication>& applications,
                        const NamingContext& ctx) const;

    /// Diagnostic emit for fail-loud paths. Aborts via the project's
    /// fprintf(stderr,"FATAL: ...")+std::abort() pattern.
    [[noreturn]] void FailUnresolved(
        const NamingContext& ctx,
        const std::vector<NamingApplication>& applications,
        std::string_view reason) const;

    /// Fires when Resolve() chooses output that the canonicality oracle
    /// rejects for the resolved chemistry context.
    [[noreturn]] void FailValidator(
        const NamingContext& ctx,
        const std::vector<NamingApplication>& applications,
        const std::string& chosen_output) const;

    void InstallRules();

    std::vector<NamingRule> rules_;
    bool debug_logging_ = false;
};

/// Process-wide singleton accessor for the applicator.
NamingApplicator& GlobalNamingApplicator();


class NamingRegistry {
public:
    NamingRegistry();

    bool IsKnownResidueName(const std::string& name) const;

    // Map any variant to canonical (standard 20 three-letter codes).
    // HID -> HIS, ASH -> ASP, CYX -> CYS, MSE -> MET, etc.
    // Returns empty string if unknown.
    std::string ToCanonical(const std::string& name) const;

    // Returns canonical name unchanged if no context-specific mapping exists.
    std::string ResolveForTool(const std::string& canonical,
                                ToolContext context,
                                const std::string& variant = "") const;

private:
    std::set<std::string> standard_residues_;

    std::map<std::string, std::string> to_canonical_;

    // (canonical, context, variant) -> tool-specific name
    struct ContextKey {
        std::string canonical;
        ToolContext context;
        std::string variant;
        bool operator<(const ContextKey& o) const;
    };
    std::map<ContextKey, std::string> context_map_;

    void InitialiseStandardResidues();
    void InitialiseAmberContext();
};

// Global singleton (C++11 Meyers singleton, thread-safe).
NamingRegistry& GlobalNamingRegistry();

}  // namespace nmr
