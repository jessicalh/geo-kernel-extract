#include "ChargeSource.h"
#include "Protein.h"
#include "ChargeAssignmentResult.h"
#include "RuntimeEnvironment.h"
#include "OperationLog.h"

#include <fstream>
#include <unordered_map>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <cmath>

namespace fs = std::filesystem;

namespace nmr {

// ============================================================================
// ParamFileChargeSource: ff14SB from the flat parameter file.
//
// This is the existing LoadParamFile + VariantResidueName logic, now
// behind the ChargeSource interface. The implementation delegates to
// ChargeAssignmentResult's existing param file parsing.
// ============================================================================

std::vector<AtomChargeRadius> ParamFileChargeSource::LoadCharges(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) const {

    // Use the existing LoadParamFile from ChargeAssignmentResult
    auto params = ChargeAssignmentResult::LoadParamFile(path_);
    if (params.empty()) {
        error_out = "cannot load parameters from " + path_;
        return {};
    }

    std::vector<AtomChargeRadius> result(conf.AtomCount());

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const Atom& identity = protein.AtomAt(ai);
        const Residue& res = protein.ResidueAt(identity.residue_index);

        // PDB LOADING BOUNDARY: variant residue name for param lookup
        std::string ff_resname = ChargeAssignmentResult::VariantResidueName(
            res.type, res.protonation_variant_index);

        std::string key = ff_resname + " " + identity.pdb_atom_name;
        auto it = params.find(key);

        if (it != params.end()) {
            result[ai] = {it->second.charge, it->second.radius};
        } else {
            // Fallback: standard residue name
            std::string std_name = ThreeLetterCodeForAminoAcid(res.type);
            std::string fallback_key = std_name + " " + identity.pdb_atom_name;
            auto fb = params.find(fallback_key);

            if (fb != params.end()) {
                result[ai] = {fb->second.charge, fb->second.radius};
            } else {
                // Element-based defaults (charge 0, reasonable radius)
                double r = 1.5;
                switch (identity.element) {
                    case Element::H: r = 0.6; break;
                    case Element::C: r = 1.7; break;
                    case Element::N: r = 1.625; break;
                    case Element::O: r = 1.48; break;
                    case Element::S: r = 1.782; break;
                    default: break;
                }
                result[ai] = {0.0, r};
            }
        }
    }

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

        // For radius: CHARMM36m uses sigma from LJ parameters.
        // The tpr dump doesn't directly give sigma per atom (it's in
        // the pair interaction matrix). Use element-based CHARMM defaults.
        int atomic_number = 0;
        std::sscanf(line.c_str() + an_pos + 11, "%d", &atomic_number);

        double radius = 1.5;
        switch (atomic_number) {
            case 1:  radius = 1.0;   break;  // H (CHARMM uses larger H than AMBER)
            case 6:  radius = 1.7;   break;  // C
            case 7:  radius = 1.625; break;  // N
            case 8:  radius = 1.48;  break;  // O
            case 16: radius = 1.782; break;  // S
        }

        all_atoms.push_back({q, radius});
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
    bool has_radii = raw_radii.size() >= n_protein;

    std::vector<AtomChargeRadius> result(n_protein);
    for (size_t i = 0; i < n_protein; ++i) {
        // Convert from AMBER internal units to elementary charges
        result[i].charge = raw_charges[i] / AMBER_CHARGE_FACTOR;
        result[i].radius = has_radii ? raw_radii[i] : 1.5;
    }

    OperationLog::Info(LogCharges, "PrmtopChargeSource::LoadCharges",
        "loaded " + std::to_string(n_protein) + " charges from " + path_ +
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
