#include "NamingRegistry.h"
#include "Residue.h"
#include "Atom.h"
#include "AminoAcidType.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <sstream>

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


// ============================================================================
// NamingRegistry: residue-name translation only. Atom-name canonicalisation
// is delegated to the NamingApplicator below.
// ============================================================================

NamingRegistry::NamingRegistry() {
    InitialiseStandardResidues();
    InitialiseAmberContext();
    InitialiseCharmmContext();
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
// Global singleton — residue-name registry
// ============================================================================

NamingRegistry& GlobalNamingRegistry() {
    static NamingRegistry instance;
    return instance;
}


// ============================================================================
// NamingApplicator — atom-name canonicalisation engine
//
// Algorithm:
//   1. Collect: iterate every rule, accumulate every NamingApplication
//      whose Applies() predicate evaluates true for this context.
//   2. Resolve: pick a single canonical output from the per-atom map
//      via the explicit per-case decision body in Resolve(). Each
//      branch documents the project decision authorising the choice
//      and is exercised by a property test.
//
// The map is transient (built per Apply() call, consumed by Resolve(),
// discarded). The chosen output is what persists, written by the caller
// onto Atom::pdb_atom_name.
//
// See spec/plan/naming-applicator-architecture-sketch-2026-05-06.md.
// ============================================================================

NamingApplicator::NamingApplicator() {
    InstallRules();
}


// ----------------------------------------------------------------------------
// Sibling-set predicates used by the rule InstallRules() function.
//
// Each predicate takes a NamingContext.sibling_input_names snapshot and
// returns true iff the residue is in the non-canonical pre-Markley state
// the rule corrects. Naming follows the convention:
//   `IsLynPreMarkley` -- siblings have HZ1+HZ2 only (LYN chemistry but
//                        non-canonical).
//   `IsCanonicalLysAmmonium` -- siblings have HZ1+HZ2+HZ3 (charged LYS,
//                                already canonical).
// ----------------------------------------------------------------------------

namespace {

bool ContainsAll(const std::set<std::string>& siblings,
                 std::initializer_list<const char*> names) {
    for (const char* n : names) {
        if (siblings.count(n) == 0) return false;
    }
    return true;
}

bool ContainsNone(const std::set<std::string>& siblings,
                  std::initializer_list<const char*> names) {
    for (const char* n : names) {
        if (siblings.count(n) != 0) return false;
    }
    return true;
}

// LYN pre-Markley shape: only HZ1+HZ2 are siblings (no HZ3). Means the
// residue carries LYN chemistry under non-canonical naming; the canonical
// AMBER ff14SB names are HZ2+HZ3.
bool IsLynPreMarkley(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HZ1", "HZ2"})
        && ContainsNone(siblings, {"HZ3"});
}

// LYN canonical shape: siblings have HZ2+HZ3 (no HZ1). Already canonical.
bool IsLynCanonical(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HZ2", "HZ3"})
        && ContainsNone(siblings, {"HZ1"});
}

// LYS canonical-charged shape: siblings have HZ1+HZ2+HZ3. NH3+ chain form.
bool IsLysAmmoniumCanonical(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HZ1", "HZ2", "HZ3"});
}

// GLY pre-Markley shape: collapsed methylene "HA" (no HA2/HA3).
bool IsGlyHaPreMarkley(const std::set<std::string>& siblings) {
    return siblings.count("HA") != 0
        && ContainsNone(siblings, {"HA2", "HA3"});
}

// GLY second pre-Markley shape: HA1+HA2 (replacing HA2+HA3).
bool IsGlyHa1Ha2PreMarkley(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HA1", "HA2"})
        && ContainsNone(siblings, {"HA3"});
}

// PRO Pdb2gmx-RTP shift: HD1+HD2 instead of canonical HD2+HD3.
bool IsProHdPdb2gmxShape(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HD1", "HD2"})
        && ContainsNone(siblings, {"HD3"});
}

// LYS Pdb2gmx-RTP δ-methylene shift: HD1+HD2 instead of HD2+HD3.
bool IsLysHdPdb2gmxShape(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HD1", "HD2"})
        && ContainsNone(siblings, {"HD3"});
}

// LYS Pdb2gmx-RTP ε-methylene shift: HE1+HE2 instead of HE2+HE3.
bool IsLysHePdb2gmxShape(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HE1", "HE2"})
        && ContainsNone(siblings, {"HE3"});
}

// ARG Pdb2gmx-RTP δ-methylene shift: HD1+HD2 instead of HD2+HD3.
bool IsArgHdPdb2gmxShape(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HD1", "HD2"})
        && ContainsNone(siblings, {"HD3"});
}

// ILE Pdb2gmx-RTP δ-methyl shape: HD1+HD2+HD3 + CD instead of canonical
// HD11+HD12+HD13 + CD1.
bool IsIleHdMethylPdb2gmxShape(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HD1", "HD2", "HD3"})
        && ContainsNone(siblings, {"HD11", "HD12", "HD13"});
}

// ILE Pdb2gmx-RTP γ-carbon shape: CD instead of CD1.
bool IsIleCdPdb2gmxShape(const std::set<std::string>& siblings) {
    return siblings.count("CD") != 0
        && ContainsNone(siblings, {"CD1"});
}

// ILE Pdb2gmx-RTP γ1-methylene shape: HG11+HG12 instead of HG12+HG13.
bool IsIleHg1Pdb2gmxShape(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HG11", "HG12"})
        && ContainsNone(siblings, {"HG13"});
}

}  // namespace


// ----------------------------------------------------------------------------
// Canonicality oracle
//
// Reads from AminoAcidType::atoms (the project-internal canonical chain
// inventory) plus the variant inventory derived per residue. A name is
// "canonical" when it appears either:
//   (a) in the chain residue's atom list for the AminoAcid type; OR
//   (b) in the per-variant overlay set (e.g. HID has HD1, HIE has HE2,
//       HIP has both); OR
//   (c) in the per-terminal cap atom inventory (H1/H2/H3 for charged
//       N-terminus; OXT for deprotonated C-terminus; HXT for protonated;
//       H1/H2 for neutral N-terminus).
//
// The oracle deliberately does NOT consult the substrate generator's
// AtomMechanicalIdentity tables; it consults the simpler chain-name
// inventory. The canonical names emitted by the generator are by
// construction also present in AminoAcidType::atoms when the variant
// is appropriate. (Two sources, same conclusion; the simpler is used
// here.)
// ----------------------------------------------------------------------------

bool NamingApplicator::IsCanonical(const NamingContext& ctx) const {
    if (ctx.input_name.empty()) return true;  // empty -> idempotent

    // (a) Chain residue inventory.
    const AminoAcidType& aatype = GetAminoAcidType(ctx.residue_type);
    for (const auto& a : aatype.atoms) {
        if (ctx.input_name == a.name) return true;
    }

    // (b) Per-variant overlay. Variants in AminoAcidType don't redeclare
    // the chain atoms; they just describe the variant's presence/absence
    // delta from canonical. The simplest model: HID has HD1, HIE has HE2,
    // HIP has HD1 and HE2, ASH has HD2, GLH has HE2, LYN has HZ2 and HZ3
    // (HZ1 absent), CYM has no HG, CYX has no HG. The chain inventory in
    // AminoAcidType already includes HD1, HE2, HG, HD1, HE2, HZ1/2/3, HG.
    // So (a) covers all canonical chain atoms across variants — what's
    // missing is to know that the *absent* atoms in a variant are still
    // legal-when-absent, which the oracle doesn't need to know (it only
    // says "yes, HZ2 is canonical for LYS"; the load path determines
    // which atoms actually appear).
    //
    // Special case: LYN canonical has HZ2 and HZ3 but NOT HZ1. The chain
    // inventory has HZ1, so HZ1 *is* a known LYS atom name. This is why
    // shift rules with sibling-aware predicates are required: the oracle
    // can't distinguish "HZ1 in canonical-charged-LYS" from "HZ1 in
    // non-canonical-LYN". Resolve() does that via map.empty + IsCanonical.

    // (c) Cap atoms per terminal_state.
    if (ctx.terminal_state == TerminalState::NtermCharged) {
        if (ctx.input_name == "H1" || ctx.input_name == "H2"
                || ctx.input_name == "H3") return true;
    } else if (ctx.terminal_state == TerminalState::NtermNeutral) {
        if (ctx.input_name == "H1" || ctx.input_name == "H2") return true;
    } else if (ctx.terminal_state == TerminalState::CtermDeprotonated) {
        if (ctx.input_name == "OXT") return true;
    } else if (ctx.terminal_state == TerminalState::CtermProtonated) {
        if (ctx.input_name == "OXT" || ctx.input_name == "HXT") return true;
    }

    // The N-terminal cap atom H1/H2/H3 might also appear on residues
    // we don't know are N-terminal (e.g. terminal_state defaulted to
    // Internal at load time). Whitelist them as canonical: the loader
    // pass will not be in a position to reject these even when
    // terminal_state is Internal. The substrate composition path takes
    // them as cap-only literals.
    if (ctx.input_name == "H1" || ctx.input_name == "H2"
            || ctx.input_name == "H3") return true;
    if (ctx.input_name == "OXT" || ctx.input_name == "HXT") return true;
    if (ctx.input_name == "H2N") return true;  // CHARMM-port NTERM literal

    return false;
}


// ----------------------------------------------------------------------------
// Collect: iterate every rule, accumulate matches.
// ----------------------------------------------------------------------------

std::vector<NamingApplication>
NamingApplicator::Collect(const NamingContext& ctx) const {
    std::vector<NamingApplication> applications;
    for (const NamingRule& rule : rules_) {
        if (rule.applies && rule.applies(ctx)) {
            applications.push_back({&rule, rule.output(ctx)});
        }
    }
    return applications;
}


// ----------------------------------------------------------------------------
// Resolve: turn the per-atom application map into a single canonical
// output. Body is the explicit per-case decision logic. Each branch
// cites the project decision authorising the choice; each branch is
// exercised by a property test.
// ----------------------------------------------------------------------------

std::string
NamingApplicator::Resolve(const std::vector<NamingApplication>& applications,
                          const NamingContext& ctx) const {
    // Branch 1: zero rules fired.
    //
    // Project decision: if the input is already canonical for this
    // (residue, variant, terminal), pass it through unchanged
    // (idempotent on canonical inputs). Otherwise fail loud — an
    // unrecognised atom name with no rule to canonicalise it is a
    // substrate gap that the user must address before composition can
    // succeed.
    //
    // Authority: spec/plan/naming-applicator-architecture-sketch-2026-05-06.md
    // §"Algorithm" — "Map empty AND input matches a known canonical
    // form for ctx: return input unchanged".
    if (applications.empty()) {
        if (IsCanonical(ctx)) return ctx.input_name;
        FailUnresolved(ctx, applications,
                       "no rule applies and input is not canonical");
    }

    // Branch 2: exactly one rule fired.
    //
    // Project decision: that rule's proposed output is the answer.
    // Trivially correct because no other rule contests the choice.
    //
    // Authority: ibid. §"Algorithm" — "Map has exactly one record:
    // return that record's proposed output".
    if (applications.size() == 1) {
        return applications.front().proposed_output;
    }

    // Branch 3: multiple rules fired.
    //
    // Project decision (3a): if every fired rule produces the SAME
    // proposed output, return that shared output. The rules agree;
    // there is no conflict. This typically arises when an atom is
    // canonical in two related sources (e.g. AmberFf14SBCanonical
    // pass-through + Markley1998 pass-through both return the input
    // unchanged for a canonical name).
    //
    // Authority: ibid. §"What's locked" — rule sets are preserved as
    // published; convergent agreement is the easiest case.
    {
        const std::string& first = applications.front().proposed_output;
        bool all_agree = true;
        for (size_t i = 1; i < applications.size(); ++i) {
            if (applications[i].proposed_output != first) {
                all_agree = false;
                break;
            }
        }
        if (all_agree) return first;
    }

    // Project decision (3b): no current branch covers this conflict
    // combination. Fail loud — the user introduces a new branch
    // (with a unit test) before the case can resolve.
    //
    // Authority: ibid. §"Algorithm" — "If a combination fires that
    // no branch covers: fail-loud (a new project decision is needed
    // before this case can be resolved)".
    FailUnresolved(ctx, applications,
                   "multiple rules fire with disagreeing outputs and no "
                   "documented branch in Resolve() covers the combination");
}


[[noreturn]] void
NamingApplicator::FailUnresolved(
        const NamingContext& ctx,
        const std::vector<NamingApplication>& applications,
        std::string_view reason) const {
    std::ostringstream map_str;
    for (const NamingApplication& app : applications) {
        map_str << "\n    " << NamingSourceName(app.rule->source) << "/"
                << app.rule->name << " -> '" << app.proposed_output << "'";
    }
    if (applications.empty()) map_str << "\n    (empty)";

    std::fprintf(stderr,
        "FATAL: NamingApplicator: atom '%s' in residue %s seq %d chain '%s' "
        "under source %s: %.*s. "
        "Map:%s. "
        "Context: variant_index=%d, terminal_state=%d. "
        "See spec/plan/naming-applicator-architecture-sketch-2026-05-06.md.\n",
        ctx.input_name.c_str(),
        GetAminoAcidType(ctx.residue_type).three_letter_code,
        ctx.residue_sequence_number,
        ctx.chain_id.c_str(),
        NamingSourceName(ctx.source),
        static_cast<int>(reason.size()), reason.data(),
        map_str.str().c_str(),
        ctx.variant_index,
        static_cast<int>(ctx.terminal_state));
    std::abort();
}


// ----------------------------------------------------------------------------
// Apply: single-atom canonicalisation. Collect + Resolve.
// ----------------------------------------------------------------------------

std::string NamingApplicator::Apply(const NamingContext& ctx) const {
    std::vector<NamingApplication> applications = Collect(ctx);
    return Resolve(applications, ctx);
}


// ----------------------------------------------------------------------------
// ApplyResidue: whole-residue convenience. Snapshots sibling input names
// ONCE for the residue, then iterates atoms with that snapshot in each
// per-atom NamingContext. This is essential for shift-pair rules: the
// rule's predicate must read the original sibling set, not a partially-
// rewritten one.
// ----------------------------------------------------------------------------

std::vector<std::string> NamingApplicator::ApplyResidue(
        const std::vector<std::string>& input_names,
        const std::vector<std::string>& parent_input_names,
        AminoAcid residue_type,
        int variant_index,
        TerminalState terminal_state,
        NamingSource source,
        int residue_sequence_number,
        std::string_view chain_id) const {
    // Snapshot: the original input names of every atom in the residue.
    // This is the set every per-atom rule predicate consults via
    // ctx.sibling_input_names. The snapshot is NOT updated as we
    // proceed atom-by-atom; the entire residue's predicates evaluate
    // against the same starting point.
    std::set<std::string> sibling_snapshot;
    for (const std::string& n : input_names) {
        if (!n.empty()) sibling_snapshot.insert(n);
    }

    std::vector<std::string> outputs;
    outputs.reserve(input_names.size());
    for (size_t k = 0; k < input_names.size(); ++k) {
        NamingContext ctx;
        ctx.source                  = source;
        ctx.input_name              = input_names[k];
        ctx.residue_type            = residue_type;
        ctx.variant_index           = variant_index;
        ctx.terminal_state          = terminal_state;
        ctx.sibling_input_names     = sibling_snapshot;
        ctx.parent_input_name       = (k < parent_input_names.size())
                                      ? parent_input_names[k]
                                      : std::string{};
        ctx.residue_sequence_number = residue_sequence_number;
        ctx.chain_id                = std::string(chain_id);
        outputs.push_back(Apply(ctx));
    }
    return outputs;
}


// ----------------------------------------------------------------------------
// InstallRules: every NamingRule is constructed here. Rules are tagged
// with the typed NamingSource that names what they actually represent,
// each carries a one-line rationale citing the source decision.
//
// Mapping from old (from_context, to_context) tags to new NamingSource:
//   (Charmm, Standard) for fleet_amber/{1P9J,1Z9B} TPR atom names
//      → NamingSource::Pdb2gmxAmberRtpDeviation (THIS IS THE RELABEL:
//        the old Charmm tag was overloaded to mean "pdb2gmx-AMBER-RTP
//        deviation in fleet_amber TPRs", NOT CHARMM force-field naming;
//        the new tag names what the rules actually represent).
//   (Standard, Charmm) for backbone H<->HN, OT1/OT2<->O/OXT
//      → NamingSource::CharmmLegacy (these rules historically lived in
//        the registry but are no longer fired by any active load path
//        post-2026-05-02 CHARMM retirement; retained for completeness
//        and any future legacy-CHARMM input).
//   (Charmm, Standard) for backbone HN->H, OT1/OT2->O/OXT (pre-Markley
//      Charmm-input collapse)
//      → NamingSource::CharmmLegacy as above. The fleet_amber TPRs that
//        flow through here were tagged as "Charmm" in the old code but
//        do not in fact carry these names (no HN, no OT1/OT2 in
//        fleet_amber topol.top).
//   (Standard, Amber) for LYN HZ shift, GLY HA collapse
//      → NamingSource::AmberFf14SBCanonical (these are AMBER ff14SB
//        canonical mappings, project decisions about what the canonical
//        target should be).
//
// Every rule's predicate now examines the sibling_input_names snapshot
// rather than firing unconditionally on a single (atom_name, residue,
// from_context, to_context) lookup. This makes the rules idempotent on
// canonical inputs (canonical sibling sets don't match the non-canonical
// pattern; the rule doesn't fire), removing the need for the
// RecanonicaliseAfterProtonation guard.
// ----------------------------------------------------------------------------

void NamingApplicator::InstallRules() {
    // ========================================================================
    // CHARMM legacy backbone atom-name rules (retired path; rules retained
    // for any future legacy input)
    // ========================================================================
    //
    // These rules fire only when ctx.source == CharmmLegacy. The active
    // load paths (PdbFileReader, FullSystemReader, OrcaRunLoader) do not
    // tag inputs with CharmmLegacy. The rules are kept so a future
    // legitimate CHARMM input path can be added by tagging its inputs
    // appropriately, without re-deriving the well-known H/HN, OT1/OT2,
    // O/OXT mappings.
    //
    // Reference: AMBER ff14SB cap convention encoded in the substrate
    // generator's kCapCterm tables; CHARMM36m residue topology.

    rules_.push_back(NamingRule{
        NamingSource::CharmmLegacy,
        "CharmmHnToCanonicalH",
        "CHARMM36m backbone amide HN -> AMBER ff14SB / IUPAC H",
        [](const NamingContext& c) {
            return c.source == NamingSource::CharmmLegacy
                && c.input_name == "HN";
        },
        [](const NamingContext&) { return std::string("H"); },
    });

    rules_.push_back(NamingRule{
        NamingSource::CharmmLegacy,
        "CharmmOt1ToCanonicalO",
        "CHARMM36m C-terminal OT1 (chain carbonyl O) -> AMBER ff14SB O",
        [](const NamingContext& c) {
            return c.source == NamingSource::CharmmLegacy
                && c.input_name == "OT1";
        },
        [](const NamingContext&) { return std::string("O"); },
    });

    rules_.push_back(NamingRule{
        NamingSource::CharmmLegacy,
        "CharmmOt2ToCanonicalOxt",
        "CHARMM36m C-terminal OT2 (carboxyl O) -> AMBER ff14SB OXT",
        [](const NamingContext& c) {
            return c.source == NamingSource::CharmmLegacy
                && c.input_name == "OT2";
        },
        [](const NamingContext&) { return std::string("OXT"); },
    });


    // ========================================================================
    // LYN protonation-variant H-on-NZ — sibling-aware shift to canonical
    // ========================================================================
    //
    // AMBER ff14SB LYN (deprotonated lysine, neutral NH2 amine on NZ)
    // has TWO H atoms on NZ canonically named HZ2 / HZ3 — preserving
    // the LYS NH3+ numbering convention with HZ1 absent (the proton
    // removed at deprotonation). Some upstream prep flows (notably
    // the 1Z9B fleet input.pdb) carry LYN-chemistry residues with
    // the H atoms named HZ1 / HZ2 (a pre-Markley-1998 numbering).
    //
    // SIBLING-AWARE: the rule fires only when ctx.sibling_input_names
    // contains HZ1+HZ2 but NOT HZ3 (the pre-Markley LYN signature).
    // This makes the rule IDEMPOTENT on:
    //   - canonical LYN (HZ2+HZ3 only): predicate evaluates false; the
    //     rule does not fire; the LynCanonicalHzPassThrough rule fires
    //     instead and returns the input unchanged.
    //   - canonical charged LYS (HZ1+HZ2+HZ3): predicate evaluates
    //     false; the LysAmmoniumHzPassThrough rule fires instead.
    //
    // Reference: AMBER ff14SB residue templates (data/ff14sb_params.dat
    // lines 433-453); Markley 1998 §2.1.1.

    rules_.push_back(NamingRule{
        NamingSource::AmberFf14SBCanonical,
        "LynHz1ToHz2_PreMarkleyShift",
        "LYN: HZ1 -> HZ2 when siblings have HZ1+HZ2 only (pre-Markley LYN)",
        [](const NamingContext& c) {
            return c.residue_type == AminoAcid::LYS
                && c.input_name == "HZ1"
                && IsLynPreMarkley(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HZ2"); },
    });

    rules_.push_back(NamingRule{
        NamingSource::AmberFf14SBCanonical,
        "LynHz2ToHz3_PreMarkleyShift",
        "LYN: HZ2 -> HZ3 when siblings have HZ1+HZ2 only (pre-Markley LYN)",
        [](const NamingContext& c) {
            return c.residue_type == AminoAcid::LYS
                && c.input_name == "HZ2"
                && IsLynPreMarkley(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HZ3"); },
    });

    // Pass-through for canonical LYN (HZ2+HZ3 in siblings, no HZ1).
    rules_.push_back(NamingRule{
        NamingSource::AmberFf14SBCanonical,
        "LynCanonicalHzPassThrough",
        "LYN: HZ2/HZ3 already canonical (siblings have HZ2+HZ3, no HZ1)",
        [](const NamingContext& c) {
            return c.residue_type == AminoAcid::LYS
                && (c.input_name == "HZ2" || c.input_name == "HZ3")
                && IsLynCanonical(c.sibling_input_names);
        },
        [](const NamingContext& c) { return c.input_name; },
    });

    // Pass-through for canonical charged LYS (HZ1+HZ2+HZ3 all present).
    rules_.push_back(NamingRule{
        NamingSource::AmberFf14SBCanonical,
        "LysAmmoniumHzPassThrough",
        "LYS: HZ1/HZ2/HZ3 canonical NH3+ when siblings have all three",
        [](const NamingContext& c) {
            return c.residue_type == AminoAcid::LYS
                && (c.input_name == "HZ1" || c.input_name == "HZ2"
                    || c.input_name == "HZ3")
                && IsLysAmmoniumCanonical(c.sibling_input_names);
        },
        [](const NamingContext& c) { return c.input_name; },
    });


    // ========================================================================
    // Glycine alpha methylene — HA <-> HA2/HA3 (sibling-aware)
    // ========================================================================
    //
    // Glycine has two prochiral alpha hydrogens; AMBER ff14SB +
    // IUPAC convention names them HA2 / HA3 (Markley 1998 §2.1.2).
    // Some fixtures collapse the methylene to a single "HA" name on
    // GLY (pre-Markley convention) when the file format only supplies
    // one of the two atoms. Map "HA" -> "HA2" so the substrate's
    // GLY chain table (which keys HA2 / HA3) finds a match.
    //
    // SIBLING-AWARE: the rule fires only when siblings contain "HA"
    // and lack HA2/HA3 (pre-Markley collapsed-methylene). On canonical
    // GLY (siblings contain HA2/HA3, no HA), predicate returns false
    // and the canonical rule fires instead.

    rules_.push_back(NamingRule{
        NamingSource::AmberFf14SBCanonical,
        "GlyHaToHa2_PreMarkley",
        "GLY: HA -> HA2 when siblings have HA only (pre-Markley collapse; Markley 1998 §2.1.2)",
        [](const NamingContext& c) {
            return c.residue_type == AminoAcid::GLY
                && c.input_name == "HA"
                && IsGlyHaPreMarkley(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HA2"); },
    });

    rules_.push_back(NamingRule{
        NamingSource::AmberFf14SBCanonical,
        "GlyHa1ToHa3_PreMarkley",
        "GLY: HA1 -> HA3 when siblings have HA1+HA2 (pre-Markley shift)",
        [](const NamingContext& c) {
            return c.residue_type == AminoAcid::GLY
                && c.input_name == "HA1"
                && IsGlyHa1Ha2PreMarkley(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HA3"); },
    });
    // Note: HA2 input with siblings {HA1,HA2} hits the empty-map
    // pass-through branch in Resolve() because HA2 is in GLY's
    // canonical chain inventory (AminoAcidType::atoms). No explicit
    // rule needed — the canonicality oracle covers this case.


    // ========================================================================
    // pdb2gmx-AMBER-RTP deviation rules — fleet-vetted load-time canonicalisation
    // ========================================================================
    //
    // When pdb2gmx writes an AMBER-ff14SB topology, side-chain methylene
    // and methyl-bearing-carbon atom names sometimes appear in older
    // AMBER-RTP form rather than the IUPAC/AMBER canonical form. The
    // 1P9J + 1Z9B fleet topol.top files carry these names (verified
    // inside their `#mol_X` rtp blocks). Without these rules the
    // canonical AMBER substrate's LookupBy fails on PRO HD1, LYS HD/HE,
    // ILE HD/CD, ILE HG1.
    //
    // THESE RULES FIRE ONLY WHEN ctx.source == Pdb2gmxAmberRtpDeviation,
    // AND siblings match the non-canonical pattern. Canonical inputs
    // (siblings already contain HD2+HD3 etc.) do not trigger these rules
    // — the predicate returns false. This makes the rules idempotent.
    //
    // Reference: spec/ChangesRequiredBeforeProductionH5Run.md;
    // h5-reader/notes/nmr_forensics/SUMMARY.md (empirical probe).

    // PRO δ-methylene: pdb2gmx HD1/HD2 → IUPAC HD2/HD3.
    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "ProHd1ToHd2_Pdb2gmxShift",
        "PRO: HD1 -> HD2 when siblings have HD1+HD2 (pdb2gmx-AMBER-RTP shift)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && c.residue_type == AminoAcid::PRO
                && c.input_name == "HD1"
                && IsProHdPdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HD2"); },
    });
    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "ProHd2ToHd3_Pdb2gmxShift",
        "PRO: HD2 -> HD3 when siblings have HD1+HD2 (pdb2gmx-AMBER-RTP shift)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && c.residue_type == AminoAcid::PRO
                && c.input_name == "HD2"
                && IsProHdPdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HD3"); },
    });

    // LYS δ-methylene: pdb2gmx HD1/HD2 → IUPAC HD2/HD3.
    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "LysHd1ToHd2_Pdb2gmxShift",
        "LYS: HD1 -> HD2 when siblings have HD1+HD2 (pdb2gmx-AMBER-RTP shift)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && c.residue_type == AminoAcid::LYS
                && c.input_name == "HD1"
                && IsLysHdPdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HD2"); },
    });
    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "LysHd2ToHd3_Pdb2gmxShift",
        "LYS: HD2 -> HD3 when siblings have HD1+HD2 (pdb2gmx-AMBER-RTP shift)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && c.residue_type == AminoAcid::LYS
                && c.input_name == "HD2"
                && IsLysHdPdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HD3"); },
    });

    // LYS ε-methylene: pdb2gmx HE1/HE2 → IUPAC HE2/HE3.
    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "LysHe1ToHe2_Pdb2gmxShift",
        "LYS: HE1 -> HE2 when siblings have HE1+HE2 (pdb2gmx-AMBER-RTP shift)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && c.residue_type == AminoAcid::LYS
                && c.input_name == "HE1"
                && IsLysHePdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HE2"); },
    });
    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "LysHe2ToHe3_Pdb2gmxShift",
        "LYS: HE2 -> HE3 when siblings have HE1+HE2 (pdb2gmx-AMBER-RTP shift)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && c.residue_type == AminoAcid::LYS
                && c.input_name == "HE2"
                && IsLysHePdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HE3"); },
    });

    // ARG δ-methylene: pdb2gmx HD1/HD2 → IUPAC HD2/HD3.
    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "ArgHd1ToHd2_Pdb2gmxShift",
        "ARG: HD1 -> HD2 when siblings have HD1+HD2 (pdb2gmx-AMBER-RTP shift)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && c.residue_type == AminoAcid::ARG
                && c.input_name == "HD1"
                && IsArgHdPdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HD2"); },
    });
    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "ArgHd2ToHd3_Pdb2gmxShift",
        "ARG: HD2 -> HD3 when siblings have HD1+HD2 (pdb2gmx-AMBER-RTP shift)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && c.residue_type == AminoAcid::ARG
                && c.input_name == "HD2"
                && IsArgHdPdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HD3"); },
    });

    // ILE δ-methyl: pdb2gmx HD1/HD2/HD3 → IUPAC HD11/HD12/HD13.
    // (3 H pseudoatom on a single methyl carbon; substrate collapses
    // their di_index to None after methyl detection, but the names
    // must canonicalise so LookupBy doesn't miss on the chain table.)
    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "IleHd1ToHd11_Pdb2gmxMethylShift",
        "ILE: HD1 -> HD11 when siblings have HD1+HD2+HD3 (pdb2gmx-AMBER-RTP methyl shift)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && c.residue_type == AminoAcid::ILE
                && c.input_name == "HD1"
                && IsIleHdMethylPdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HD11"); },
    });
    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "IleHd2ToHd12_Pdb2gmxMethylShift",
        "ILE: HD2 -> HD12 when siblings have HD1+HD2+HD3 (pdb2gmx-AMBER-RTP methyl shift)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && c.residue_type == AminoAcid::ILE
                && c.input_name == "HD2"
                && IsIleHdMethylPdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HD12"); },
    });
    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "IleHd3ToHd13_Pdb2gmxMethylShift",
        "ILE: HD3 -> HD13 when siblings have HD1+HD2+HD3 (pdb2gmx-AMBER-RTP methyl shift)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && c.residue_type == AminoAcid::ILE
                && c.input_name == "HD3"
                && IsIleHdMethylPdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HD13"); },
    });

    // ILE γ-carbon: pdb2gmx CD → IUPAC CD1.
    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "IleCdToCd1_Pdb2gmxShift",
        "ILE: CD -> CD1 when siblings have CD (no CD1) (pdb2gmx-AMBER-RTP shift)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && c.residue_type == AminoAcid::ILE
                && c.input_name == "CD"
                && IsIleCdPdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("CD1"); },
    });

    // ILE γ1-methylene: pdb2gmx HG11/HG12 → IUPAC HG12/HG13. The HG21,
    // HG22, HG23 names on the γ2-methyl already match canon — only HG1*
    // get rewritten.
    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "IleHg11ToHg12_Pdb2gmxShift",
        "ILE: HG11 -> HG12 when siblings have HG11+HG12 (pdb2gmx-AMBER-RTP shift)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && c.residue_type == AminoAcid::ILE
                && c.input_name == "HG11"
                && IsIleHg1Pdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HG12"); },
    });
    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "IleHg12ToHg13_Pdb2gmxShift",
        "ILE: HG12 -> HG13 when siblings have HG11+HG12 (pdb2gmx-AMBER-RTP shift)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && c.residue_type == AminoAcid::ILE
                && c.input_name == "HG12"
                && IsIleHg1Pdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HG13"); },
    });
}


NamingApplicator& GlobalNamingApplicator() {
    static NamingApplicator instance;
    return instance;
}

}  // namespace nmr
