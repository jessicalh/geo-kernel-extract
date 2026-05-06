#pragma once
//
// NamingRegistry: the single authority for residue and atom name translation.
//
// Internal architecture (2026-05-06): rule-application object model. Every
// atom-name canonicalisation request is processed by a NamingApplicator
// that iterates a flat list of NamingRule entries, accumulates a transient
// per-atom map of every rule that fired, and resolves the map to a single
// canonical output via an explicit Resolve() method. Rule sets are
// preserved as published (each tagged with a NamingSource enum value
// identifying its scientific or project source); the choice between
// sources happens explicitly in Resolve(), not by editing rules to play
// nice. See spec/plan/naming-applicator-architecture-sketch-2026-05-06.md.
//
// Every tool in the pipeline has its own naming universe:
//   PDB/Standard: HIS, CYS, H (backbone amide)
//   AMBER:        HID/HIE/HIP, CYX/CYM, H
//   ORCA:         standard PDB names (numbered atoms)
//   OpenMM:       AMBER-like (HIE, HID, CYX)
//
// CHARMM36m (formerly HSD/HSE/HSP, CYS2, ASPP, GLUP, HN) was retired
// 2026-05-02; the `Charmm` ToolContext slot was removed 2026-05-06
// (codex D1). NamingRegistry::ToCanonical still recognises HSD/HSE/HSP
// and CYS2/ASPP/GLUP at INPUT (a CHARMM-prepared file landing in the
// pipeline still maps onto the standard 3-letter codes); there is no
// canonical -> CHARMM emitter on the output side.
//
// This registry translates between them. It is the GATEKEEPER:
// unknown names are rejected. No guessing. No silent pass-through.
//
// After translation, typed objects carry their properties. The registry
// exists so that the rest of the system NEVER needs strings for identity.
//

#include "Types.h"
#include "SemanticEnums.h"   // for AminoAcid (via Types.h chain) and TerminalState
#include <cstdint>
#include <functional>
#include <memory>
#include <set>
#include <string>
#include <string_view>
#include <map>
#include <vector>

namespace nmr {

// The tools we talk to. Each has its own naming conventions.
//
// ToolContext is the EXTERNAL surface for residue-name translation
// (HIS <-> HID/HIE/HIP etc.). It is NOT the internal source-tag for
// atom-name rule routing — that role is played by NamingSource, see
// below.
//
// CHARMM force-field support was retired 2026-05-02 (memory entry
// `project_charmm_retired_amber_only_2026-05-02`); the previous
// `Charmm` enum value (formerly slot 2: HSD/HSE/HSP, CYS2, ASPP,
// GLUP, HN) was removed 2026-05-06 (codex Finding D1) after
// investigation found zero active production consumers. The numeric
// slot it occupied is deliberately RETIRED — do not re-use for a
// different concept; a future legitimate CHARMM-input path can be
// added with a fresh enum value at that time.
enum class ToolContext {
    Standard = 0,    // PDB / IUPAC canonical names
    Amber    = 1,    // AMBER force fields (ff14SB, ff19SB): HID/HIE/HIP, CYX, ASH, GLH
    // Charmm = 2  — RETIRED slot, removed 2026-05-06 (codex D1).
    // DO NOT reuse this numeric value for a different concept; if
    // any persisted records contain it they would silently re-bind
    // to an unrelated meaning.
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


// ============================================================================
// NamingSource — typed enum identifying which scientific or project source
// authored an atom-name rule. Symmetric with SemanticSource on the substrate
// side (SemanticEnums.h). New sources = new enum values. Each value names
// the source the rule body actually represents — historic ToolContext
// labels did not (e.g. the now-retired ToolContext::Charmm slot tagged
// AMBER pdb2gmx-RTP rules for fleet_amber TPRs, NOT CHARMM force-field
// naming; ToolContext::Charmm itself was removed 2026-05-06 codex D1).
// ============================================================================

enum class NamingSource : uint8_t {
    Unknown                  = 0,

    /// Project-internal canonical AMBER ff14SB. The TARGET of all
    /// canonicalisation: every rule's Output() produces a name in this
    /// vocabulary. Rules tagged AmberFf14SBCanonical are project decisions
    /// about the canonical mapping (e.g. LYN HZ1+HZ2 -> HZ2+HZ3 is an
    /// AMBER-vs-pre-Markley LYN naming decision; GLY HA -> HA2 is the
    /// Markley 1998 §2.1.2 collapsed-methylene canonical mapping).
    AmberFf14SBCanonical     = 1,

    /// AMBER pdb2gmx writes side-chain methylenes and methyl-bearing
    /// carbons in older AMBER RTP convention; the resulting atom names
    /// deviate from canonical AMBER ff14SB / IUPAC. Verified inside
    /// fleet_amber/{1P9J,1Z9B}/topol.top mol_X rtp blocks.
    ///
    /// THIS ENUM VALUE REPLACES the historic `ToolContext::Charmm`
    /// atom-name source tag that was overloaded to mean "pdb2gmx-AMBER-
    /// RTP deviation in fleet_amber TPRs". CHARMM-the-force-field is
    /// retired (2026-05-02 quarantined-legacy); the only live consumer
    /// of the rules formerly tagged Charmm is the AMBER trajectory
    /// loader reading pdb2gmx-AMBER topologies via libgromacs-direct.
    Pdb2gmxAmberRtpDeviation = 2,

    /// cifpp's PDB parsing yields IUPAC convention names (mostly
    /// canonical AMBER post-Markley). The registered cifpp surface is
    /// PdbFileReader.cpp:131. Inputs from this source rarely need
    /// rewriting; rules tagged here are typically idempotent
    /// pass-through over a set the source already gets right.
    CifppPdbInput            = 3,

    /// ORCA NMR-output echoes the input atom names (prmtop ATOM_NAME
    /// surface in OrcaRunLoader.cpp:206). prmtop-loaded names are
    /// already AMBER ff14SB canonical; rules tagged here handle
    /// any residual cases.
    OrcaEcho                 = 4,

    /// Markley et al. 1998 J. Biomol. NMR 12:1-23 nomenclature
    /// recommendations (the IUPAC-IUB 1969 update for proteins). The
    /// Greek-letter-locant + diastereotopic-numbering canon lives
    /// here when expressed as a runtime rule.
    Markley1998              = 6,

    /// BMRB nomenclature table at https://bmrb.io/ref_info. Currently
    /// reserved; no runtime rules tagged here yet.
    BmrbAtomNomTbl           = 7,

    /// IUPAC-IUB 1969 tentative rules. Currently reserved; no
    /// runtime rules tagged here yet.
    IupacIub1969             = 8,

    /// Project-internal synthesised conventions (e.g. Trp 6-ring
    /// perimeter labels per topology-encoding-dependencies §C.3).
    /// Currently reserved; no runtime rules tagged here yet.
    ProjectSynthesis         = 9,

    // NOTE: `CharmmLegacy = 5` was removed 2026-05-06 (codex-review
    // Finding 2). It was never tagged by any active load path; the
    // 3 rules carrying it (HN->H, OT1->O, OT2->OXT) had no live
    // emitter (verified: PdbFileReader / FullSystemReader /
    // OrcaRunLoader all produce canonical AMBER ff14SB; fleet_amber
    // 1Z9B and 1P9J fixtures have neither HN nor OT1/OT2 nor H2N).
    // CHARMM-the-force-field is retired (2026-05-02 memory entry
    // `project_charmm_retired_amber_only_2026-05-02`); a future
    // legitimate CHARMM-input path can be added with a non-CHARMM-
    // confused source tag at that time. The numeric value 5 is
    // retired; do NOT reuse for a different concept (would alias old
    // serialised bytes if any persisted records exist).
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


// ============================================================================
// NamingContext — input record passed to the applicator. Carries everything
// a rule's predicate might examine. New context fields land here as new
// rule types demand them.
// ============================================================================

struct NamingContext {
    NamingSource  source         = NamingSource::Unknown;
    std::string   input_name;
    AminoAcid     residue_type   = AminoAcid::Unknown;

    /// Protonation variant index into AminoAcidType::variants. -1 means
    /// "unresolved at load time" (Pass 1 over the applicator); >= 0
    /// means resolved (Pass 2, post-ResolveProtonationStates). Variant-
    /// aware rules require >= 0 in their predicates.
    int           variant_index  = -1;

    /// Terminal state from SemanticEnums.h (TerminalState; the typed
    /// enum that distinguishes NtermCharged vs NtermNeutral and
    /// CtermDeprotonated vs CtermProtonated). Defaults to Internal.
    TerminalState terminal_state = TerminalState::Internal;

    /// Sibling atom names IN THE INPUT FORM. Snapshot of the residue's
    /// atom names at the start of the current pass through the
    /// applicator, before any rules rewrite them. The snapshot is
    /// critical for shift-pair rules (e.g. Pdb2gmx ILE HD1/HD2/HD3 ↔
    /// canonical HD11/HD12/HD13): the rule's predicate examines the
    /// snapshot to determine non-canonical state, NOT the partially-
    /// renamed state. Stored as std::string (transient lifetime owned
    /// by the caller's residue snapshot vector) for simplicity over
    /// string_view aliasing.
    std::set<std::string> sibling_input_names;

    /// Parent heavy atom's input name (for hydrogen disambiguation).
    /// Empty for non-hydrogen atoms or when parent is unknown.
    std::string parent_input_name;

    /// Diagnostics-only fields (not consumed by rule predicates;
    /// included in fail-loud messages for triage).
    int         residue_sequence_number = 0;
    std::string chain_id;
};


// ============================================================================
// NamingRule — single rule. Has a source tag, an Applies predicate, an
// Output transform, plus name + rationale for diagnostics and self-docs.
// ============================================================================

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


// ============================================================================
// NamingApplication — one rule firing on one atom. Per-atom records of
// these are accumulated in the transient map.
// ============================================================================

struct NamingApplication {
    const NamingRule* rule;             ///< Which rule fired.
    std::string       proposed_output;  ///< What that rule proposes.
};


// ============================================================================
// NamingApplicator — the containing object. Owns rules, the resolver, the
// canonicality oracle. Lives in NamingRegistry.{h,cpp}; replaces the
// historic flat-map architecture (atom_name_map_).
// ============================================================================

namespace test {
// Forward-declared test access helper. Lives in the tests/ tree;
// production code never instantiates it. Its sole purpose is to keep
// the test-only NamingApplicator(custom_rules) constructor private —
// production callers cannot stumble across it on the public API
// (codex Finding D3, 2026-05-06). See
// tests/test_naming_applicator_access.h.
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
    // Test-only constructor: build an applicator with custom rules.
    // The supplied rule list REPLACES the production rule set; the
    // canonicality oracle and resolution body are unchanged. Used by
    // `tests/test_naming_registry.cpp` to exercise resolution branches
    // (notably the "multiple rules agree" branch) that the production
    // rule set does not currently reach. NOT for production use;
    // production code calls `GlobalNamingApplicator()`.
    //
    // Codex Finding D3 (2026-05-06): this constructor is private, only
    // reachable through the friend declaration below; the test code
    // gates entry through nmr::test::NamingApplicatorTestAccess. The
    // production public API does not expose a way to install a
    // hand-rolled rule list.
    friend class test::NamingApplicatorTestAccess;
    struct CustomRules { std::vector<NamingRule> rules; };
    explicit NamingApplicator(CustomRules custom);

    /// Step 1: iterate rules, return all applications.
    std::vector<NamingApplication> Collect(const NamingContext& ctx) const;

    /// Step 2: resolve the per-atom map to a single output. Body is
    /// the explicit per-case decision logic; each branch documented
    /// with a citation to the project decision authorising the choice.
    std::string Resolve(const std::vector<NamingApplication>& applications,
                        const NamingContext& ctx) const;

    /// Diagnostic emit for fail-loud paths. Aborts via the project's
    /// fprintf(stderr,"FATAL: ...")+std::abort() pattern.
    [[noreturn]] void FailUnresolved(
        const NamingContext& ctx,
        const std::vector<NamingApplication>& applications,
        std::string_view reason) const;

    /// Post-resolution validator's fail-loud emit. Fires when Resolve()
    /// returned a chosen output that the canonicality oracle rejects
    /// for the resolved chemistry context — a rule misbehaved by
    /// proposing non-canonical output. Distinct from FailUnresolved
    /// (which fires when no rule fires or a multi-rule conflict has no
    /// documented branch). The diagnostic names every rule that fired,
    /// the original input, and the bad output. Codex round-2,
    /// 2026-05-06.
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
    //   ResolveForTool("HIS", ToolContext::Amber, "delta")   -> "HID"
    //   ResolveForTool("HIS", ToolContext::Amber, "epsilon") -> "HIE"
    //   ResolveForTool("HIS", ToolContext::Amber, "doubly")  -> "HIP"
    // For non-titratable residues, variant is ignored:
    //   ResolveForTool("ALA", ToolContext::Amber) -> "ALA"
    // Returns canonical name unchanged if no mapping exists.
    std::string ResolveForTool(const std::string& canonical,
                                ToolContext context,
                                const std::string& variant = "") const;

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

    void InitialiseStandardResidues();
    void InitialiseAmberContext();
    // InitialiseCharmmContext() removed 2026-05-06 (codex D1):
    // CHARMM force-field support retired; no live consumer needed
    // canonical -> CHARMM residue-name emission.
};

// Global singleton (C++11 Meyers singleton, thread-safe).
NamingRegistry& GlobalNamingRegistry();

}  // namespace nmr
