#include "ChargeSource.h"
#include "AmberChargeResolver.h"
#include "Protein.h"
#include "RuntimeEnvironment.h"
#include "OperationLog.h"

#include <fstream>
#include <unordered_map>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <cmath>
#include <vector>

namespace fs = std::filesystem;

namespace nmr {

namespace {

struct ParamEntry {
    double partial_charge = 0.0;
    double pb_radius = 0.0;
};

bool IsTerminalStateToken(const std::string& token) {
    return token == "INTERNAL" || token == "NTERM" ||
           token == "CTERM" || token == "NCTERM";
}

std::string ParamKey(const std::string& terminal_state,
                     const std::string& residue_name,
                     const std::string& atom_name) {
    return terminal_state + " " + residue_name + " " + atom_name;
}

std::unordered_map<std::string, ParamEntry>
LoadFf14sbParamFile(const std::string& path) {
    std::unordered_map<std::string, ParamEntry> params;

    std::ifstream in(path);
    if (!in.is_open()) {
        OperationLog::Error("ParamFileChargeSource::LoadFf14sbParamFile",
            "cannot open " + path);
        return params;
    }

    std::string line;
    size_t legacy_rows = 0;
    size_t skipped_rows = 0;
    while (std::getline(in, line)) {
        auto comment = line.find('#');
        if (comment != std::string::npos) line = line.substr(0, comment);

        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        while (iss >> token) tokens.push_back(token);
        if (tokens.empty()) continue;

        try {
            if (IsTerminalStateToken(tokens[0])) {
                if (tokens.size() < 5) {
                    ++skipped_rows;
                    continue;
                }
                double charge = std::stod(tokens[3]);
                double pb_radius = std::stod(tokens[4]);
                if (pb_radius <= 0.0) {
                    ++skipped_rows;
                    continue;
                }
                params[ParamKey(tokens[0], tokens[1], tokens[2])] =
                    {charge, pb_radius};
            } else {
                // Legacy table format:
                //   RESNAME ATOMNAME CHARGE LJ_EPSILON RADIUS
                // Treat as INTERNAL only. Terminal residues must be present
                // explicitly in the regenerated table.
                if (tokens.size() < 5) {
                    ++skipped_rows;
                    continue;
                }
                double charge = std::stod(tokens[2]);
                double pb_radius = std::stod(tokens[4]);
                if (pb_radius <= 0.0) {
                    ++skipped_rows;
                    continue;
                }
                params[ParamKey("INTERNAL", tokens[0], tokens[1])] =
                    {charge, pb_radius};
                ++legacy_rows;
            }
        } catch (...) {
            ++skipped_rows;
        }
    }

    OperationLog::Info(LogCharges, "ParamFileChargeSource::LoadFf14sbParamFile",
        "loaded " + std::to_string(params.size()) + " entries from " + path +
        " legacy_rows=" + std::to_string(legacy_rows) +
        " skipped_rows=" + std::to_string(skipped_rows));

    return params;
}

bool Ff14sbVariantResidueName(
        AminoAcid aa,
        int protonation_variant_index,
        std::string& residue_name_out,
        std::string& error_out) {
    const AminoAcidType& aa_type = GetAminoAcidType(aa);

    if (protonation_variant_index < 0) {
        // ff14SB lacks a neutral generic HIS in the flat file path used here.
        if (aa == AminoAcid::HIS) {
            residue_name_out = "HIE";
        } else {
            residue_name_out = ThreeLetterCodeForAminoAcid(aa);
        }
        return true;
    }

    if (protonation_variant_index >=
            static_cast<int>(aa_type.variants.size())) {
        error_out = "invalid protonation variant index " +
                    std::to_string(protonation_variant_index) + " for " +
                    ThreeLetterCodeForAminoAcid(aa);
        return false;
    }

    residue_name_out = aa_type.variants[protonation_variant_index].name;
    return true;
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

std::vector<std::string> AtomNameCandidates(
        const std::string& atom_name,
        ResidueTerminalState terminal_state) {
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
// ParamFileChargeSource: ff14SB from the flat parameter file.
//
// The flat file is string-keyed, so string lookup is unavoidable inside this
// construction-boundary adapter. It produces only AtomChargeRadius rows for the
// prepared ForceFieldChargeTable; ChargeAssignmentResult is only a projection.
// ============================================================================

std::vector<AtomChargeRadius> ParamFileChargeSource::LoadCharges(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) const {

    // Single source of truth for "can the flat table cover this protein."
    // The verdict is the typed predicate; this loader only runs when the
    // verdict is satisfiable. The early-fail message preserves the
    // existing test contract (terminal_token + ff_resname + "no canonical
    // fallback" substrings).
    auto verdict = AnalyzeFlatTableCoverage(protein, path_);
    if (!verdict.Ok()) {
        error_out = verdict.Detail();
        OperationLog::Error("ParamFileChargeSource::LoadCharges", error_out);
        return {};
    }

    auto params = LoadFf14sbParamFile(path_);
    if (params.empty()) {
        error_out = "cannot load parameters from " + path_;
        return {};
    }

    std::vector<AtomChargeRadius> result(conf.AtomCount());
    size_t internal_matches = 0;
    size_t terminal_matches = 0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const Atom& identity = protein.AtomAt(ai);
        const Residue& res = protein.ResidueAt(identity.residue_index);

        // PDB LOADING BOUNDARY: variant residue name for param lookup.
        // The verdict already proved every (terminal, resname, atom)
        // triple resolves; here we only need the value.
        std::string ff_resname;
        std::string variant_error;
        if (!Ff14sbVariantResidueName(
                res.type, res.protonation_variant_index,
                ff_resname, variant_error)) {
            error_out = variant_error + " at residue " +
                        std::to_string(res.sequence_number);
            OperationLog::Error("ParamFileChargeSource::LoadCharges", error_out);
            return {};
        }
        const std::string terminal_token =
            TerminalStateToken(res.terminal_state);

        const auto atom_candidates =
            AtomNameCandidates(identity.pdb_atom_name, res.terminal_state);

        const ParamEntry* entry = nullptr;
        for (const auto& atom_name : atom_candidates) {
            auto it = params.find(ParamKey(terminal_token, ff_resname, atom_name));
            if (it != params.end()) {
                entry = &it->second;
                break;
            }
        }

        if (!entry) {
            // Verdict said satisfiable but per-atom lookup missed. This is
            // a contract violation between AnalyzeFlatTableCoverage and the
            // loader's parser — not a runtime data condition.
            fprintf(stderr,
                "FATAL: ParamFileChargeSource::LoadCharges: verdict claimed "
                "Satisfiable but lookup missed (%s %s %s). "
                "AnalyzeFlatTableCoverage and LoadFf14sbParamFile have "
                "diverged.\n",
                terminal_token.c_str(), ff_resname.c_str(),
                identity.pdb_atom_name.c_str());
            std::abort();
        }

        result[ai] = {entry->partial_charge, entry->pb_radius,
                      ChargeAssignmentStatus::Matched};
        if (terminal_token == "INTERNAL") ++internal_matches;
        else ++terminal_matches;
    }

    OperationLog::Info(LogCharges, "ParamFileChargeSource::LoadCharges",
        "internal_matches=" + std::to_string(internal_matches) +
        " terminal_matches=" + std::to_string(terminal_matches) +
        " atoms=" + std::to_string(result.size()));

    return result;
}


// ============================================================================
// GmxTprChargeSource: charges from GROMACS .tpr via `gmx dump`.
//
// PDB LOADING BOUNDARY: parses `gmx dump` text output for atom lines:
//   atom[  0]={type=0, ..., q=-3.00000e-01, ..., atomnumber=7}
//
// Atoms are in topology order which must match the Protein's atom order.
// The .tpr was built from the same protonated PDB we loaded, so atom
// count and order should match. Mismatch is a hard error.
// ============================================================================

std::vector<AtomChargeRadius> GmxTprChargeSource::LoadCharges(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) const {

    if (!fs::exists(tpr_)) {
        error_out = "tpr file not found: " + tpr_;
        return {};
    }
    if (!fs::exists(gmx_)) {
        error_out = "gmx binary not found: " + gmx_;
        return {};
    }

    // Run gmx dump, capture output in a guid-tagged temp file
    std::string tpr_stem = fs::path(tpr_).stem().string();
    std::string tmp_out = RuntimeEnvironment::TempFilePath(tpr_stem, "gmxdump.txt");
    std::string cmd = gmx_ + " dump -s " + tpr_ + " 2>/dev/null > " + tmp_out;
    int rc = std::system(cmd.c_str());
    if (rc != 0) {
        error_out = "gmx dump failed (rc=" + std::to_string(rc) + ") on " + tpr_;
        std::error_code ec;
        fs::remove(tmp_out, ec);
        return {};
    }

    // PDB LOADING BOUNDARY: parse atom lines from gmx dump output
    std::ifstream in(tmp_out);
    if (!in.is_open()) {
        error_out = "cannot read gmx dump output: " + tmp_out;
        return {};
    }

    // We need protein-only atoms. The tpr includes water and ions.
    // We extract ALL atom charges, then take the first N matching the
    // protein atom count. GROMACS puts protein atoms first in the topology.
    std::vector<AtomChargeRadius> all_atoms;
    std::string line;

    // Pattern: atom[  123]={..., q= 1.23456e-01, ..., atomnumber=  7}
    while (std::getline(in, line)) {
        auto q_pos = line.find("q=");
        auto an_pos = line.find("atomnumber=");
        if (q_pos == std::string::npos || an_pos == std::string::npos) continue;
        if (line.find("atom[") == std::string::npos) continue;

        // Extract charge
        double q = 0.0;
        std::sscanf(line.c_str() + q_pos + 2, "%lf", &q);

        all_atoms.push_back({
            q,
            kCompatibilityPlaceholderPbRadiusAngstrom,
            ChargeAssignmentStatus::PlaceholderPbRadius
        });
    }
    in.close();

    // Clean up
    std::error_code ec;
    fs::remove(tmp_out, ec);

    // The protein atoms are the first N atoms in the tpr topology.
    size_t n_protein = conf.AtomCount();
    if (all_atoms.size() < n_protein) {
        error_out = "tpr has " + std::to_string(all_atoms.size()) +
                    " atoms, protein needs " + std::to_string(n_protein);
        return {};
    }

    std::vector<AtomChargeRadius> result(all_atoms.begin(),
                                          all_atoms.begin() + n_protein);

    OperationLog::Info(LogCharges, "GmxTprChargeSource::LoadCharges",
        "loaded " + std::to_string(n_protein) + " charges from " + tpr_ +
        " (tpr total: " + std::to_string(all_atoms.size()) + ")");
    OperationLog::Warn("GmxTprChargeSource::LoadCharges",
        "using quarantined compatibility PB radii for GROMACS/CHARMM TPR path; "
        "TODO(charge-model): replace with defensible CHARMM/GROMACS PB radii");

    return result;
}


// ============================================================================
// PrmtopChargeSource: AMBER prmtop format parser.
//
// PDB LOADING BOUNDARY: reads Fortran-formatted text sections from the
// AMBER topology file. After parsing, everything is doubles.
//
// Format: sections delimited by %FLAG <NAME> / %FORMAT(<fmt>)
// Charges: %FLAG CHARGE, format 5E16.8, AMBER internal units (e * 18.2223)
// Radii: %FLAG RADII, format 5E16.8, PB radii in Angstroms
// ============================================================================

// Read a section of doubles from a prmtop file.
// Reads all values after %FLAG <flag_name> until the next %FLAG.
static std::vector<double> ReadPrmtopSection(const std::string& path,
                                              const std::string& flag_name) {
    std::vector<double> values;
    std::ifstream in(path);
    if (!in.is_open()) return values;

    std::string line;
    bool in_section = false;
    bool found_format = false;

    while (std::getline(in, line)) {
        if (line.find("%FLAG " + flag_name) != std::string::npos) {
            in_section = true;
            found_format = false;
            continue;
        }
        if (in_section && !found_format && line.find("%FORMAT") != std::string::npos) {
            found_format = true;
            continue;
        }
        if (in_section && found_format) {
            if (line.find("%FLAG") != std::string::npos) break;  // next section

            // Parse space-separated or fixed-width Fortran doubles
            // E16.8 format: 5 values per line, 16 chars each
            for (size_t pos = 0; pos + 15 < line.size(); pos += 16) {
                std::string token = line.substr(pos, 16);
                double val = 0.0;
                if (std::sscanf(token.c_str(), "%lf", &val) == 1) {
                    values.push_back(val);
                }
            }
            // Handle lines shorter than a full row
            if (line.size() > 0 && line.size() % 16 != 0) {
                size_t last_start = (line.size() / 16) * 16;
                if (last_start < line.size()) {
                    std::string token = line.substr(last_start);
                    double val = 0.0;
                    if (std::sscanf(token.c_str(), "%lf", &val) == 1) {
                        values.push_back(val);
                    }
                }
            }
        }
    }
    return values;
}


std::vector<AtomChargeRadius> PrmtopChargeSource::LoadCharges(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) const {
    (void)protein;

    if (!fs::exists(path_)) {
        error_out = "prmtop not found: " + path_;
        return {};
    }

    auto raw_charges = ReadPrmtopSection(path_, "CHARGE");
    auto raw_radii = ReadPrmtopSection(path_, "RADII");

    if (raw_charges.empty()) {
        error_out = "no CHARGE section in prmtop: " + path_;
        return {};
    }

    size_t n_protein = conf.AtomCount();

    if (raw_charges.size() < n_protein) {
        error_out = "prmtop has " + std::to_string(raw_charges.size()) +
                    " charges, protein needs " + std::to_string(n_protein);
        return {};
    }

    // RADII section may be absent in older prmtops
    if (raw_radii.size() < n_protein) {
        error_out = "prmtop has " + std::to_string(raw_radii.size()) +
                    " PB radii, protein needs " + std::to_string(n_protein);
        return {};
    }

    std::vector<AtomChargeRadius> result(n_protein);
    for (size_t i = 0; i < n_protein; ++i) {
        // Convert from AMBER internal units to elementary charges
        result[i].partial_charge = raw_charges[i] / AMBER_CHARGE_FACTOR;
        result[i].pb_radius = raw_radii[i];
        result[i].status = ChargeAssignmentStatus::Matched;
    }

    OperationLog::Info(LogCharges, "PrmtopChargeSource::LoadCharges",
        "loaded " + std::to_string(n_protein) + " charge/PB-radius rows from " + path_ +
        " (prmtop total: " + std::to_string(raw_charges.size()) + ")");

    return result;
}


// ============================================================================
// PreloadedChargeSource
// ============================================================================

std::vector<AtomChargeRadius> PreloadedChargeSource::LoadCharges(
        const Protein& /*protein*/,
        const ProteinConformation& conf,
        std::string& error_out) const {

    if (charges_.size() != conf.AtomCount()) {
        error_out = "PreloadedChargeSource: expected " +
                    std::to_string(conf.AtomCount()) + " atoms, have " +
                    std::to_string(charges_.size()) + " charges";
        return {};
    }
    return charges_;
}

}  // namespace nmr
