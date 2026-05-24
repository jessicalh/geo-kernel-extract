#include "AmberChargeResolver.h"
#include "AminoAcidType.h"
#include "AmberPreparedChargeSource.h"
#include "Atom.h"
#include "Residue.h"
#include "ChargeSource.h"
#include "ProteinBuildContext.h"
#include "OperationLog.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <utility>

namespace fs = std::filesystem;

namespace nmr {

const char* AmberPreparationPolicyName(AmberPreparationPolicy policy) {
    switch (policy) {
        case AmberPreparationPolicy::UseStockTermini:
            return "UseStockTermini";
        case AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants:
            return "UseCappedFragmentsForUnsupportedTerminalVariants";
        case AmberPreparationPolicy::FailOnUnsupportedTerminalVariants:
            return "FailOnUnsupportedTerminalVariants";
    }
    return "Unknown";
}

const char* AmberFlatTableCoverageKindName(AmberFlatTableCoverageKind kind) {
    switch (kind) {
        case AmberFlatTableCoverageKind::Satisfiable:
            return "Satisfiable";
        case AmberFlatTableCoverageKind::UnsupportedTerminalVariant:
            return "UnsupportedTerminalVariant";
        case AmberFlatTableCoverageKind::UnsupportedResidue:
            return "UnsupportedResidue";
        case AmberFlatTableCoverageKind::MissingAtomName:
            return "MissingAtomName";
    }
    return "Unknown";
}

const AmberFlatTableCoverageFailure&
AmberFlatTableCoverageVerdict::FirstFailure() const {
    return failures.front();
}

std::string AmberFlatTableCoverageVerdict::Detail() const {
    if (ok) return "Satisfiable";

    std::ostringstream oss;
    oss << "ff14SB flat table cannot represent " << failures.size()
        << " case" << (failures.size() == 1 ? "" : "s") << "; ";

    // First failure carries the prominent error text the existing tests
    // expect: terminal_token + ff_resname + "no canonical fallback".
    const auto& f0 = failures.front();
    switch (f0.kind) {
        case AmberFlatTableCoverageKind::UnsupportedTerminalVariant:
            oss << "no " << f0.terminal_token << " template for "
                << f0.ff_residue_name;
            break;
        case AmberFlatTableCoverageKind::UnsupportedResidue:
            oss << "no template for residue " << f0.ff_residue_name
                << " (" << f0.terminal_token << ")";
            break;
        case AmberFlatTableCoverageKind::MissingAtomName:
            oss << "no row for " << f0.terminal_token << " "
                << f0.ff_residue_name << " " << f0.atom_name;
            break;
        case AmberFlatTableCoverageKind::Satisfiable:
            break;
    }
    oss << "; no canonical fallback is allowed. Prepare this case with "
           "an AMBER-supported terminal template, cap/custom AMBER "
           "template, or prmtop-derived charges.";

    if (failures.size() > 1) {
        oss << " Additional cases: ";
        const size_t kMaxExtra = 8;
        const size_t n = std::min<size_t>(failures.size() - 1, kMaxExtra);
        for (size_t i = 0; i < n; ++i) {
            const auto& f = failures[i + 1];
            if (i > 0) oss << ", ";
            oss << f.terminal_token << " " << f.ff_residue_name;
            if (!f.atom_name.empty()) oss << " " << f.atom_name;
            oss << " (residue " << f.residue_sequence_number << ")";
        }
        if (failures.size() - 1 > kMaxExtra) {
            oss << " (+ " << (failures.size() - 1 - kMaxExtra) << " more)";
        }
    }

    return oss.str();
}

namespace {

// ---------------------------------------------------------------------------
// Flat-table key parser. Builds the set of valid (terminal, resname, atom)
// keys plus the set of residue prefixes (terminal+resname) so we can
// distinguish "residue not present at all" from "atom missing inside an
// otherwise-present residue template".
// ---------------------------------------------------------------------------

struct Ff14sbKeySet {
    std::unordered_set<std::string> atom_keys;          // "INTERNAL ALA CA"
    std::unordered_set<std::string> residue_prefixes;   // "INTERNAL ALA "
    bool loaded = false;
};

bool IsTerminalStateToken(const std::string& token) {
    return token == "INTERNAL" || token == "NTERM" ||
           token == "CTERM" || token == "NCTERM";
}

Ff14sbKeySet LoadKeySet(const std::string& path) {
    Ff14sbKeySet ks;
    std::ifstream in(path);
    if (!in.is_open()) {
        OperationLog::Error("AnalyzeFlatTableCoverage::LoadKeySet",
            "cannot open " + path);
        return ks;
    }

    std::string line;
    while (std::getline(in, line)) {
        auto comment = line.find('#');
        if (comment != std::string::npos) line = line.substr(0, comment);

        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        while (iss >> token) tokens.push_back(token);
        if (tokens.empty()) continue;

        if (IsTerminalStateToken(tokens[0])) {
            if (tokens.size() < 5) continue;
            try {
                double pb_radius = std::stod(tokens[4]);
                if (pb_radius <= 0.0) continue;
            } catch (...) {
                continue;
            }
            ks.atom_keys.insert(
                tokens[0] + " " + tokens[1] + " " + tokens[2]);
            ks.residue_prefixes.insert(tokens[0] + " " + tokens[1] + " ");
        } else {
            // Legacy table format: RESNAME ATOMNAME CHARGE LJ_EPSILON RADIUS
            if (tokens.size() < 5) continue;
            try {
                double pb_radius = std::stod(tokens[4]);
                if (pb_radius <= 0.0) continue;
            } catch (...) {
                continue;
            }
            ks.atom_keys.insert(
                "INTERNAL " + tokens[0] + " " + tokens[1]);
            ks.residue_prefixes.insert("INTERNAL " + tokens[0] + " ");
        }
    }

    ks.loaded = true;
    return ks;
}

bool ResidueNameAnyTerminalPresent(
        const Ff14sbKeySet& ks, const std::string& ff_resname) {
    for (const char* tok : {"INTERNAL ", "NTERM ", "CTERM ", "NCTERM "}) {
        if (ks.residue_prefixes.count(std::string(tok) + ff_resname + " ") > 0) {
            return true;
        }
    }
    return false;
}

std::string TerminalStateToken(ResidueTerminalState terminal_state) {
    switch (terminal_state) {
        case ResidueTerminalState::Internal:      return "INTERNAL";
        case ResidueTerminalState::NTerminus:     return "NTERM";
        case ResidueTerminalState::CTerminus:     return "CTERM";
        case ResidueTerminalState::NAndCTerminus: return "NCTERM";
        case ResidueTerminalState::Unknown:       return "UNKNOWN";
    }
    return "UNKNOWN";
}

// Resolve the AMBER residue name for the flat-table lookup. Mirrors
// ParamFileChargeSource's Ff14sbVariantResidueName so the verdict and
// the actual loader agree on what residue name they're checking.
bool Ff14sbResidueName(
        const Residue& res,
        std::string& name_out) {
    const AminoAcidType& aa_type = res.AminoAcidInfo();

    if (res.protonation_variant_index < 0) {
        if (res.type == AminoAcid::HIS) {
            name_out = "HIE";
        } else {
            name_out = ThreeLetterCodeForAminoAcid(res.type);
        }
        return true;
    }

    if (res.protonation_variant_index >=
            static_cast<int>(aa_type.variants.size())) {
        return false;
    }
    name_out = aa_type.variants[res.protonation_variant_index].name;
    return true;
}

std::vector<std::string> AtomNameCandidates(
        const std::string& atom_name, ResidueTerminalState terminal_state) {
    std::vector<std::string> candidates;
    candidates.push_back(atom_name);

    const bool n_terminal =
        terminal_state == ResidueTerminalState::NTerminus ||
        terminal_state == ResidueTerminalState::NAndCTerminus;
    if (n_terminal && (atom_name == "H" || atom_name == "HN")) {
        candidates.push_back("H1");
    }
    return candidates;
}

}  // namespace


// ============================================================================
// AnalyzeFlatTableCoverage
// ============================================================================

AmberFlatTableCoverageVerdict AnalyzeFlatTableCoverage(
        const Protein& protein,
        const std::string& flat_table_path) {

    AmberFlatTableCoverageVerdict verdict;

    Ff14sbKeySet ks = LoadKeySet(flat_table_path);
    if (!ks.loaded) {
        verdict.ok = false;
        AmberFlatTableCoverageFailure f;
        f.kind = AmberFlatTableCoverageKind::UnsupportedResidue;
        f.terminal_token = "UNKNOWN";
        f.ff_residue_name = "(table_unreadable)";
        verdict.failures.push_back(f);
        return verdict;
    }

    // Track residues already reported as missing-template/missing-residue
    // so we don't emit one failure per atom of the same residue.
    std::unordered_set<std::string> reported_residue_keys;

    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        const std::string terminal_token = TerminalStateToken(res.terminal_state);

        std::string ff_resname;
        if (!Ff14sbResidueName(res, ff_resname)) {
            // Invalid variant index → treat as unsupported residue. This is
            // a Protein-construction issue, not a table issue.
            verdict.ok = false;
            AmberFlatTableCoverageFailure f;
            f.kind = AmberFlatTableCoverageKind::UnsupportedResidue;
            f.terminal_token = terminal_token;
            f.ff_residue_name =
                ThreeLetterCodeForAminoAcid(res.type) + "(invalid_variant)";
            f.residue_sequence_number = res.sequence_number;
            f.chain_id = res.chain_id;
            verdict.failures.push_back(f);
            continue;
        }

        const std::string residue_key = terminal_token + " " + ff_resname;
        const std::string residue_prefix = residue_key + " ";

        const bool has_template_for_terminal =
            ks.residue_prefixes.count(residue_prefix) > 0;
        if (!has_template_for_terminal) {
            if (reported_residue_keys.insert(residue_key).second) {
                verdict.ok = false;
                AmberFlatTableCoverageFailure f;
                if (ResidueNameAnyTerminalPresent(ks, ff_resname)) {
                    f.kind = AmberFlatTableCoverageKind::UnsupportedTerminalVariant;
                } else {
                    f.kind = AmberFlatTableCoverageKind::UnsupportedResidue;
                }
                f.terminal_token = terminal_token;
                f.ff_residue_name = ff_resname;
                f.residue_sequence_number = res.sequence_number;
                f.chain_id = res.chain_id;
                verdict.failures.push_back(f);
            }
            continue;
        }

        // Template exists for this terminal_state + ff_resname. Check every
        // atom of the residue against (terminal_token, ff_resname, atom_name).
        for (size_t ai : res.atom_indices) {
            const Atom& atom = protein.AtomAt(ai);
            const auto candidates =
                AtomNameCandidates(atom.pdb_atom_name, res.terminal_state);
            bool atom_found = false;
            for (const auto& aname : candidates) {
                if (ks.atom_keys.count(
                        residue_key + " " + aname) > 0) {
                    atom_found = true;
                    break;
                }
            }
            if (!atom_found) {
                verdict.ok = false;
                AmberFlatTableCoverageFailure f;
                f.kind = AmberFlatTableCoverageKind::MissingAtomName;
                f.terminal_token = terminal_token;
                f.ff_residue_name = ff_resname;
                f.atom_name = atom.pdb_atom_name;
                f.residue_sequence_number = res.sequence_number;
                f.chain_id = res.chain_id;
                verdict.failures.push_back(f);
            }
        }
    }

    return verdict;
}


// ============================================================================
// ResolveAmberChargeSource
// ============================================================================

std::unique_ptr<ChargeSource> ResolveAmberChargeSource(
        const Protein& protein,
        const ProteinBuildContext& build_context,
        const AmberSourceConfig& config,
        std::string& error_out) {

    // Branch 1: upstream PRMTOP wins, always. Never re-tleap.
    if (!build_context.prmtop_path.empty()) {
        if (!fs::exists(build_context.prmtop_path)) {
            error_out =
                "ResolveAmberChargeSource: build_context.prmtop_path does "
                "not exist: " + build_context.prmtop_path;
            return nullptr;
        }
        return std::make_unique<PrmtopChargeSource>(
            build_context.prmtop_path, ForceField::Amber_ff14SB);
    }

    if (config.flat_table_path.empty()) {
        error_out = "ResolveAmberChargeSource: flat_table_path is empty";
        return nullptr;
    }

    // Branch 2: typed predicate over Protein state.
    auto verdict = AnalyzeFlatTableCoverage(protein, config.flat_table_path);
    if (verdict.Ok()) {
        return std::make_unique<ParamFileChargeSource>(config.flat_table_path);
    }

    // Branch 3: prepared topology via tleap.
    if (config.preparation_policy ==
            AmberPreparationPolicy::FailOnUnsupportedTerminalVariants) {
        error_out = "AMBER charge resolution failed under "
                    "FailOnUnsupportedTerminalVariants policy: " +
                    verdict.Detail();
        return nullptr;
    }

    // Step-4 scaffolding present: source constructs and exposes its
    // generated PDB / LEaP script for testing. Step 5 wires the actual
    // tleap invocation; until then LoadCharges returns a named
    // step-5-not-implemented error.
    return std::make_unique<AmberPreparedChargeSource>(
        protein, config.preparation_policy, std::move(verdict), config);
}

}  // namespace nmr
