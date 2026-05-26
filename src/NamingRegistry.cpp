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


NamingRegistry::NamingRegistry() {
    InitialiseStandardResidues();
    InitialiseAmberContext();
}


void NamingRegistry::InitialiseStandardResidues() {
    standard_residues_ = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    };

    for (const auto& name : standard_residues_)
        to_canonical_[name] = name;

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
        context_map_[{canonical, ToolContext::OpenMM, variant}] = amber_name;
    };

    add("HIS", "delta",   "HID");
    add("HIS", "epsilon", "HIE");
    add("HIS", "doubly",  "HIP");
    add("HIS", "",        "HIS");

    add("ASP", "protonated", "ASH");
    add("GLU", "protonated", "GLH");

    add("CYS", "disulfide",    "CYX");
    add("CYS", "deprotonated", "CYM");

    add("LYS", "deprotonated", "LYN");

    add("ARG", "deprotonated", "ARN");

    add("TYR", "deprotonated", "TYM");
}

bool NamingRegistry::IsKnownResidueName(const std::string& name) const {
    std::string upper = ToUpper(Trim(name));
    return to_canonical_.count(upper) > 0;
}

std::string NamingRegistry::ToCanonical(const std::string& name) const {
    std::string upper = ToUpper(Trim(name));
    auto it = to_canonical_.find(upper);
    if (it != to_canonical_.end())
        return it->second;
    return "";
}

std::string NamingRegistry::ResolveForTool(const std::string& canonical,
                                            ToolContext context,
                                            const std::string& variant) const {
    auto it = context_map_.find({canonical, context, variant});
    if (it != context_map_.end())
        return it->second;

    if (!variant.empty()) {
        it = context_map_.find({canonical, context, ""});
        if (it != context_map_.end())
            return it->second;
    }

    return canonical;
}

NamingRegistry& GlobalNamingRegistry() {
    static NamingRegistry instance;
    return instance;
}


NamingApplicator::NamingApplicator() {
    InstallRules();
}


// Custom rules for tests; Resolve() and the canonicality oracle are unchanged.
NamingApplicator::NamingApplicator(CustomRules custom)
    : rules_(std::move(custom.rules)) {
}


// Sibling-set predicates used by the rule InstallRules() function.

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

// LYN neutral amine canonical hydrogens are HZ2+HZ3.
bool IsLynPreMarkley(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HZ1", "HZ2"})
        && ContainsNone(siblings, {"HZ3"});
}

bool IsLynCanonical(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HZ2", "HZ3"})
        && ContainsNone(siblings, {"HZ1"});
}

bool IsLysAmmoniumCanonical(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HZ1", "HZ2", "HZ3"});
}

bool IsGlyHaPreMarkley(const std::set<std::string>& siblings) {
    return siblings.count("HA") != 0
        && ContainsNone(siblings, {"HA2", "HA3"});
}

bool IsGlyHa1Ha2PreMarkley(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HA1", "HA2"})
        && ContainsNone(siblings, {"HA3"});
}

bool IsProHdPdb2gmxShape(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HD1", "HD2"})
        && ContainsNone(siblings, {"HD3"});
}

bool IsLysHdPdb2gmxShape(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HD1", "HD2"})
        && ContainsNone(siblings, {"HD3"});
}

bool IsLysHePdb2gmxShape(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HE1", "HE2"})
        && ContainsNone(siblings, {"HE3"});
}

bool IsArgHdPdb2gmxShape(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HD1", "HD2"})
        && ContainsNone(siblings, {"HD3"});
}

// ILE pdb2gmx-RTP delta-methyl shape: HD1+HD2+HD3 + CD instead of canonical
// HD11+HD12+HD13 + CD1.
bool IsIleHdMethylPdb2gmxShape(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HD1", "HD2", "HD3"})
        && ContainsNone(siblings, {"HD11", "HD12", "HD13"});
}

bool IsIleCdPdb2gmxShape(const std::set<std::string>& siblings) {
    return siblings.count("CD") != 0
        && ContainsNone(siblings, {"CD1"});
}

bool IsIleHg1Pdb2gmxShape(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HG11", "HG12"})
        && ContainsNone(siblings, {"HG13"});
}

// beta-methylene pdb2gmx-RTP shape: HB1+HB2 instead of HB2+HB3.
bool IsBetaMethylenePdb2gmxShape(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HB1", "HB2"})
        && ContainsNone(siblings, {"HB3"});
}

// gamma-methylene pdb2gmx-RTP shape: HG1+HG2 instead of HG2+HG3.
bool IsGammaMethylenePdb2gmxShape(const std::set<std::string>& siblings) {
    return ContainsAll(siblings, {"HG1", "HG2"})
        && ContainsNone(siblings, {"HG3"});
}

// beta-methylene set: every residue with a CB-bearing methylene (CB has
// two H atoms). EXCLUDES:
//   - ALA (HB1/HB2/HB3 is canonical methyl on CB; not a methylene).
//   - GLY (no CB).
//   - VAL, ILE, THR (single beta H).
//
// gamma-methylene set: every residue with a CG-bearing methylene (CG has
// two H atoms named HG2/HG3 canonically). The set is exactly:
//   ARG, GLN, GLU, LYS, MET, PRO. Note: ILE has gamma1-methylene
//   (HG11/HG12 pdb2gmx to HG12/HG13 canonical) but ILE-specific rules
//   above already cover that case; ILE is not in this set.
bool IsBetaMethyleneResidue(AminoAcid aa) {
    switch (aa) {
        case AminoAcid::ARG:
        case AminoAcid::ASN:
        case AminoAcid::ASP:
        case AminoAcid::CYS:
        case AminoAcid::GLN:
        case AminoAcid::GLU:
        case AminoAcid::HIS:
        case AminoAcid::LEU:
        case AminoAcid::LYS:
        case AminoAcid::MET:
        case AminoAcid::PHE:
        case AminoAcid::PRO:
        case AminoAcid::SER:
        case AminoAcid::TRP:
        case AminoAcid::TYR:
            return true;
        default:
            return false;
    }
}

bool IsGammaMethyleneResidue(AminoAcid aa) {
    switch (aa) {
        case AminoAcid::ARG:
        case AminoAcid::GLN:
        case AminoAcid::GLU:
        case AminoAcid::LYS:
        case AminoAcid::MET:
        case AminoAcid::PRO:
            return true;
        default:
            return false;
    }
}

}  // namespace


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
// The oracle intentionally uses the chain-name inventory, not generated
// AtomMechanicalIdentity tables.

bool NamingApplicator::IsCanonical(const NamingContext& ctx) const {
    if (ctx.input_name.empty()) return true;  // empty -> idempotent

    // Some protonation variants delete atoms that are present in the
    // base AminoAcidType chain inventory. Under the resolved variant
    // (variant_index >= 0), the deleted atom is not canonical even
    // though the chain-inventory pass below would accept it. Without
    // this overlay, a chemistry mistake (e.g. CYX residue with HG
    // present) silently passes the oracle and substrate composition
    // looks up an atom that the resolved variant's substrate row does
    // not carry. The deletion overlay catches it before chain-pass.
    //
    //   CYS:  variant 0 = CYX (disulfide, no HG on Sγ)
    //         variant 1 = CYM (thiolate,  no HG on Sγ)
    //   LYS:  variant 0 = LYN (neutral amine; canonical hydrogens are
    //         HZ2+HZ3, no HZ1; Markley 1998 §2.1.1 numbering)
    //   TYR:  variant 0 = TYM (deprotonated phenolate, no HH on Oη)
    //   ARG:  variant 0 = ARN (deprotonated guanidinium, no HE on Nε)
    //
    // The unresolved-load state (variant_index == -1) remains permissive.
    if (ctx.variant_index >= 0) {
        switch (ctx.residue_type) {
            case AminoAcid::CYS:
                if ((ctx.variant_index == 0 /*CYX*/ ||
                     ctx.variant_index == 1 /*CYM*/) &&
                    ctx.input_name == "HG") {
                    return false;
                }
                break;
            case AminoAcid::LYS:
                if (ctx.variant_index == 0 /*LYN*/ &&
                    ctx.input_name == "HZ1") {
                    return false;
                }
                break;
            case AminoAcid::TYR:
                if (ctx.variant_index == 0 /*TYM*/ &&
                    ctx.input_name == "HH") {
                    return false;
                }
                break;
            case AminoAcid::ARG:
                if (ctx.variant_index == 0 /*ARN*/ &&
                    ctx.input_name == "HE") {
                    return false;
                }
                break;
            default: break;
        }
    }

    const AminoAcidType& aatype = GetAminoAcidType(ctx.residue_type);
    for (const auto& a : aatype.atoms) {
        if (ctx.input_name == a.name) return true;
    }

    // Per-variant atom-extension overlay. AminoAcidType::atoms is
    // populated for the maximally-protonated default state, with
    // variant-only atoms (HID's HD1, ASH's HD2, GLH's HE2) not included.
    // The substrate generator's per-variant tables encode these; we
    // hard-code the (residue, variant-extension) atom set here so the
    // oracle recognises canonical variant-only atoms without dragging
    // the substrate tables into runtime.
    //
    // HID adds HD1, HIE adds HE2 (already in chain for HIS), HIP adds HD1,
    // ASH adds HD2, GLH adds HE2. CYS variants (CYX, CYM), LYS LYN,
    // ARG ARN, TYR TYM either remove atoms (variant-deletion) or
    // re-purpose existing chain entries.
    //
    // variant_index == -1 accepts the union of extension atoms because
    // the tautomer may not be known yet. Resolved variants accept only
    // atoms belonging to that variant.
    const int v = ctx.variant_index;
    switch (ctx.residue_type) {
        case AminoAcid::HIS:
            if (ctx.input_name == "HD1") {
                // HIE(1) has only HE2 protonated.
                if (v == -1 || v == 0 || v == 2) return true;
            }
            break;
        case AminoAcid::ASP:
            if (ctx.input_name == "HD2") {
                if (v == -1 || v == 0) return true;
            }
            break;
        case AminoAcid::GLU:
            if (ctx.input_name == "HE2") {
                if (v == -1 || v == 0) return true;
            }
            break;
        default: break;
    }
    // Note: LYN canonical has HZ2 and HZ3, but not HZ1. The chain
    // inventory has HZ1, so HZ1 *is* a known LYS atom name; the oracle
    // can't distinguish "HZ1 in canonical-charged-LYS" from "HZ1 in
    // non-canonical-LYN". Resolve() does that via map.empty + IsCanonical:
    // LYN sibling-aware shift rules fire on the {HZ1,HZ2,no HZ3} pattern;
    // canonical sibling sets bypass the shift rules and hit the
    // canonicality oracle's chain-inventory match.

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

    // H1/H2/H3 may appear before terminal_state is known; substrate
    // composition handles them as cap-only literals.
    if (ctx.input_name == "H1" || ctx.input_name == "H2"
            || ctx.input_name == "H3") return true;
    if (ctx.input_name == "OXT" || ctx.input_name == "HXT") return true;

    return false;
}


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


std::string
NamingApplicator::Resolve(const std::vector<NamingApplication>& applications,
                          const NamingContext& ctx) const {
    if (applications.empty()) {
        if (IsCanonical(ctx)) return ctx.input_name;
        FailUnresolved(ctx, applications,
                       "no rule applies and input is not canonical");
    }

    if (applications.size() == 1) {
        return applications.front().proposed_output;
    }

    // Multiple agreeing rules are accepted; disagreeing rules abort below.
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
        "See spec/plan/bones/naming-applicator-architecture-sketch-2026-05-06.md.\n",
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


[[noreturn]] void
NamingApplicator::FailValidator(
        const NamingContext& ctx,
        const std::vector<NamingApplication>& applications,
        const std::string& chosen_output) const {
    std::ostringstream map_str;
    for (const NamingApplication& app : applications) {
        map_str << "\n    " << NamingSourceName(app.rule->source) << "/"
                << app.rule->name << " (" << app.rule->rationale
                << ") -> '" << app.proposed_output << "'";
    }
    if (applications.empty()) map_str << "\n    (empty)";

    std::fprintf(stderr,
        "FATAL: NamingApplicator post-resolution validator: rule produced "
        "non-canonical output for resolved chemistry context. "
        "Original input '%s' resolved to '%s' (BAD); residue %s seq %d "
        "chain '%s'; source %s; variant_index=%d; terminal_state=%d. "
        "Rules that fired:%s. "
        "Architectural contract (codex round-2, 2026-05-06): rules may "
        "recognise non-canonical INPUT (rules exist to repair) but their "
        "OUTPUT must be canonical for the resolved chemistry context; "
        "the canonicality oracle is the authority on canonical form. "
        "See spec/plan/bones/naming-applicator-architecture-sketch-2026-05-06.md "
        "and the 2026-05-06 codex round-2 review.\n",
        ctx.input_name.c_str(),
        chosen_output.c_str(),
        GetAminoAcidType(ctx.residue_type).three_letter_code,
        ctx.residue_sequence_number,
        ctx.chain_id.c_str(),
        NamingSourceName(ctx.source),
        ctx.variant_index,
        static_cast<int>(ctx.terminal_state),
        map_str.str().c_str());
    std::abort();
}


std::string NamingApplicator::Apply(const NamingContext& ctx) const {
    if (ctx.source == NamingSource::Unknown) {
        std::fprintf(stderr,
            "FATAL: NamingApplicator::Apply: ctx.source = "
            "NamingSource::Unknown is forbidden. Loaders must tag every "
            "atom with a real source (CifppPdbInput, "
            "Pdb2gmxAmberRtpDeviation, OrcaEcho, AmberFf14SBCanonical, "
            "etc.); Unknown indicates a loader bug. Atom '%s' in "
            "residue %s seq %d chain '%s'. "
            "See spec/plan/bones/naming-applicator-architecture-sketch-"
            "2026-05-06.md and codex Finding CC2.\n",
            ctx.input_name.c_str(),
            GetAminoAcidType(ctx.residue_type).three_letter_code,
            ctx.residue_sequence_number,
            ctx.chain_id.c_str());
        std::abort();
    }
    std::vector<NamingApplication> applications = Collect(ctx);
    std::string output = Resolve(applications, ctx);

    // Rules may repair non-canonical input, but the chosen output must
    // be canonical for the resolved chemistry context.
    NamingContext output_ctx = ctx;
    output_ctx.input_name = output;
    if (!IsCanonical(output_ctx)) {
        FailValidator(ctx, applications, output);
    }

    return output;
}


std::vector<std::string> NamingApplicator::ApplyResidue(
        const std::vector<std::string>& input_names,
        const std::vector<std::string>& parent_input_names,
        AminoAcid residue_type,
        int variant_index,
        TerminalState terminal_state,
        NamingSource source,
        int residue_sequence_number,
        std::string_view chain_id) const {
    // Every atom in the residue sees the original sibling-name set.
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


void NamingApplicator::InstallRules() {
    // AMBER ff14SB LYN (deprotonated lysine, neutral NH2 amine on NZ)
    // has two H atoms on NZ canonically named HZ2 / HZ3, preserving
    // the LYS NH3+ numbering convention with HZ1 absent (the proton
    // removed at deprotonation). Some upstream prep flows (notably
    // the 1Z9B fleet input.pdb) carry LYN-chemistry residues with
    // the H atoms named HZ1 / HZ2 (a pre-Markley-1998 numbering).
    //
    // Sibling predicates make the LYN and charged-LYS rules idempotent:
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

    // This rule represents charged LYS chemistry (NH3+, HZ1+HZ2+HZ3 on Nζ).
    // Under
    // resolved LYN (variant_index = 0), the residue is the neutral
    // amine NH2 form whose canonical hydrogens are HZ2+HZ3 only.
    rules_.push_back(NamingRule{
        NamingSource::AmberFf14SBCanonical,
        "LysAmmoniumHzPassThrough",
        "LYS: HZ1/HZ2/HZ3 canonical NH3+ when siblings have all three "
        "AND variant unresolved (default-charged-LYS); LYN gates this off",
        [](const NamingContext& c) {
            return c.residue_type == AminoAcid::LYS
                && c.variant_index < 0
                && (c.input_name == "HZ1" || c.input_name == "HZ2"
                    || c.input_name == "HZ3")
                && IsLysAmmoniumCanonical(c.sibling_input_names);
        },
        [](const NamingContext& c) { return c.input_name; },
    });

    // Glycine has two prochiral alpha hydrogens; AMBER ff14SB +
    // IUPAC convention names them HA2 / HA3 (Markley 1998 §2.1.2).
    // Some fixtures collapse the methylene to a single "HA" name on
    // GLY (pre-Markley convention) when the file format only supplies
    // one of the two atoms. Map "HA" -> "HA2" so the substrate's
    // GLY chain table (which keys HA2 / HA3) finds a match.
    //
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
    // When pdb2gmx writes an AMBER-ff14SB topology, side-chain methylene
    // and methyl-bearing-carbon atom names sometimes appear in older
    // AMBER-RTP form rather than the IUPAC/AMBER canonical form. The
    // 1P9J + 1Z9B fleet topol.top files carry these names (verified
    // inside their `#mol_X` rtp blocks). Without these rules the
    // canonical AMBER substrate's LookupBy fails on PRO HD1, LYS HD/HE,
    // ILE HD/CD, ILE HG1.
    //
    // These rules require ctx.source == Pdb2gmxAmberRtpDeviation and
    // a non-canonical sibling pattern, so canonical inputs do not fire.
    //
    // PRO delta-methylene: pdb2gmx HD1/HD2 -> IUPAC HD2/HD3.
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

    // LYS delta-methylene: pdb2gmx HD1/HD2 -> IUPAC HD2/HD3.
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

    // LYS epsilon-methylene: pdb2gmx HE1/HE2 -> IUPAC HE2/HE3.
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

    // ARG delta-methylene: pdb2gmx HD1/HD2 -> IUPAC HD2/HD3.
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

    // ILE delta-methyl: pdb2gmx HD1/HD2/HD3 -> IUPAC HD11/HD12/HD13.
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

    // ILE gamma-carbon: pdb2gmx CD -> IUPAC CD1.
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

    // ILE gamma1-methylene: pdb2gmx HG11/HG12 -> IUPAC HG12/HG13. The HG21,
    // HG22, HG23 names on the gamma2-methyl already match canon; only HG1*
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

    // Source rationale: AMBER pdb2gmx emits β/γ methylene H atoms in the
    // older AMBER RTP convention (HB1/HB2; HG1/HG2). Canonical AMBER ff14SB
    // (post-Markley 1998 §2.1.2 numbering) is HB2/HB3 and HG2/HG3.
    //
    // Residue-set scope:
    //   β-methylene: ARG, ASN, ASP, CYS, GLN, GLU, HIS, LEU, LYS, MET,
    //     PHE, PRO, SER, TRP, TYR (15 residues). ALA excluded
    //     (HB1/HB2/HB3 canonical methyl; not a methylene). VAL/ILE/THR
    //     excluded (single Hβ; not a methylene). GLY excluded (no β atom).
    //   γ-methylene: ARG, GLN, GLU, LYS, MET, PRO (6 residues). ILE
    //     γ1-methylene HG11/HG12 already covered by IleHg11/12_Pdb2gmxShift
    //     above; ILE not in this set.
    //
    // Reference: Markley et al. 1998 J. Biomol. NMR 12:1-23 §2.1.2.

    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "BetaMethyleneHb1ToHb2_Pdb2gmxShift",
        "β-methylene HB1 -> HB2 (Markley 1998 §2.1.2); applies to all "
        "residues with β-methylene side-chain (ARG, ASN, ASP, CYS, GLN, "
        "GLU, HIS, LEU, LYS, MET, PHE, PRO, SER, TRP, TYR); ALA excluded "
        "(HB1/HB2/HB3 canonical methyl); VAL/ILE/THR excluded (single Hβ); "
        "GLY excluded (no β atom)",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && IsBetaMethyleneResidue(c.residue_type)
                && c.input_name == "HB1"
                && IsBetaMethylenePdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HB2"); },
    });

    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "BetaMethyleneHb2ToHb3_Pdb2gmxShift",
        "β-methylene HB2 -> HB3 (Markley 1998 §2.1.2); residue-set + "
        "sibling-pattern predicate matches BetaMethyleneHb1ToHb2 above",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && IsBetaMethyleneResidue(c.residue_type)
                && c.input_name == "HB2"
                && IsBetaMethylenePdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HB3"); },
    });

    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "GammaMethyleneHg1ToHg2_Pdb2gmxShift",
        "γ-methylene HG1 -> HG2 (Markley 1998 §2.1.2); applies to "
        "residues with γ-methylene chain (ARG, GLN, GLU, LYS, MET, PRO); "
        "ILE γ1-methylene HG11/HG12 covered by ILE-specific rules above",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && IsGammaMethyleneResidue(c.residue_type)
                && c.input_name == "HG1"
                && IsGammaMethylenePdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HG2"); },
    });

    rules_.push_back(NamingRule{
        NamingSource::Pdb2gmxAmberRtpDeviation,
        "GammaMethyleneHg2ToHg3_Pdb2gmxShift",
        "γ-methylene HG2 -> HG3 (Markley 1998 §2.1.2); residue-set + "
        "sibling-pattern predicate matches GammaMethyleneHg1ToHg2 above",
        [](const NamingContext& c) {
            return c.source == NamingSource::Pdb2gmxAmberRtpDeviation
                && IsGammaMethyleneResidue(c.residue_type)
                && c.input_name == "HG2"
                && IsGammaMethylenePdb2gmxShape(c.sibling_input_names);
        },
        [](const NamingContext&) { return std::string("HG3"); },
    });
}


NamingApplicator& GlobalNamingApplicator() {
    static NamingApplicator instance;
    return instance;
}

}  // namespace nmr
