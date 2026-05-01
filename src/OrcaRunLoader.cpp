#include "OrcaRunLoader.h"
#include "AminoAcidType.h"
#include "NamingRegistry.h"
#include "ChargeSource.h"
#include "AmberChargeResolver.h"
#include "ForceFieldChargeTable.h"
#include "OperationLog.h"

#include <fstream>
#include <sstream>
#include <filesystem>
#include <cstdio>

namespace fs = std::filesystem;

namespace nmr {

// Internal result type (pre-BuildResult). LoadWithPrmtop produces a Protein;
// the public BuildFromOrca wraps it with charges and net_charge.
struct OrcaLoadInternal {
    std::unique_ptr<Protein> protein;
    bool ok = false;
    std::string error;
};

// ============================================================================
// Prmtop section readers (PDB LOADING BOUNDARY)
// ============================================================================

static std::vector<std::string> ReadPrmtopStrings(
        const std::string& path, const std::string& flag_name, int field_width) {
    std::vector<std::string> values;
    std::ifstream in(path);
    if (!in.is_open()) return values;

    std::string line;
    bool in_section = false, found_format = false;

    while (std::getline(in, line)) {
        if (line.find("%FLAG " + flag_name) != std::string::npos) {
            in_section = true; found_format = false; continue;
        }
        if (in_section && !found_format && line.find("%FORMAT") != std::string::npos) {
            found_format = true; continue;
        }
        if (in_section && found_format) {
            if (line.find("%FLAG") != std::string::npos) break;
            for (size_t pos = 0; pos + static_cast<size_t>(field_width) <= line.size();
                 pos += field_width) {
                std::string s = line.substr(pos, field_width);
                size_t start = s.find_first_not_of(' ');
                size_t end = s.find_last_not_of(' ');
                if (start != std::string::npos)
                    values.push_back(s.substr(start, end - start + 1));
            }
        }
    }
    return values;
}

static std::vector<int> ReadPrmtopInts(
        const std::string& path, const std::string& flag_name) {
    std::vector<int> values;
    std::ifstream in(path);
    if (!in.is_open()) return values;

    std::string line;
    bool in_section = false, found_format = false;

    while (std::getline(in, line)) {
        if (line.find("%FLAG " + flag_name) != std::string::npos) {
            in_section = true; found_format = false; continue;
        }
        if (in_section && !found_format && line.find("%FORMAT") != std::string::npos) {
            found_format = true; continue;
        }
        if (in_section && found_format) {
            if (line.find("%FLAG") != std::string::npos) break;
            std::istringstream iss(line);
            int v;
            while (iss >> v) values.push_back(v);
        }
    }
    return values;
}


// ============================================================================
// XYZ reader (PDB LOADING BOUNDARY)
// ============================================================================

struct XyzAtom { std::string element; double x, y, z; };

static std::vector<XyzAtom> ReadXyz(const std::string& path) {
    std::vector<XyzAtom> atoms;
    std::ifstream in(path);
    if (!in.is_open()) return atoms;

    std::string line;
    if (!std::getline(in, line)) return atoms;
    int n = 0;
    std::sscanf(line.c_str(), "%d", &n);
    if (n <= 0) return atoms;

    std::getline(in, line);  // title

    atoms.reserve(n);
    while (std::getline(in, line) && static_cast<int>(atoms.size()) < n) {
        XyzAtom a;
        char elem[4] = {};
        if (std::sscanf(line.c_str(), " %3s %lf %lf %lf",
                         elem, &a.x, &a.y, &a.z) == 4) {
            a.element = elem;
            atoms.push_back(a);
        }
    }
    return atoms;
}

static Element ElementFromAtomicNumber(int an) {
    switch (an) {
        case 1:  return Element::H;
        case 6:  return Element::C;
        case 7:  return Element::N;
        case 8:  return Element::O;
        case 16: return Element::S;
        default: return Element::Unknown;
    }
}


// ============================================================================
// LoadOrcaRun: WITH prmtop path
// ============================================================================

static OrcaLoadInternal LoadWithPrmtop(const OrcaRunFiles& files,
                                      const std::vector<XyzAtom>& xyz) {
    OrcaLoadInternal result;

    auto atom_names = ReadPrmtopStrings(files.prmtop_path, "ATOM_NAME", 4);
    auto res_labels = ReadPrmtopStrings(files.prmtop_path, "RESIDUE_LABEL", 4);
    auto res_pointers = ReadPrmtopInts(files.prmtop_path, "RESIDUE_POINTER");
    auto atomic_numbers = ReadPrmtopInts(files.prmtop_path, "ATOMIC_NUMBER");

    if (atom_names.empty() || res_labels.empty() || res_pointers.empty()) {
        result.error = "incomplete prmtop: " + files.prmtop_path;
        return result;
    }

    size_t n_atoms = atom_names.size();
    size_t n_residues = res_labels.size();

    if (xyz.size() != n_atoms) {
        result.error = "XYZ has " + std::to_string(xyz.size()) +
                       " atoms, prmtop has " + std::to_string(n_atoms);
        return result;
    }

    auto& registry = GlobalNamingRegistry();

    // Build Protein using Protein's public API
    auto protein = std::make_unique<Protein>();

    // Build residue ranges (1-based → 0-based)
    struct ResRange { size_t start; size_t end; };
    std::vector<ResRange> ranges(n_residues);
    for (size_t ri = 0; ri < n_residues; ++ri) {
        ranges[ri].start = res_pointers[ri] - 1;
        ranges[ri].end = (ri + 1 < n_residues)
                         ? static_cast<size_t>(res_pointers[ri + 1] - 1)
                         : n_atoms;
    }

    // Add residues
    for (size_t ri = 0; ri < n_residues; ++ri) {
        Residue res;

        // PDB LOADING BOUNDARY: AMBER residue name → canonical
        std::string canonical = registry.ToCanonical(res_labels[ri]);
        if (canonical.empty()) {
            result.error = "unknown residue: " + res_labels[ri] +
                           " at position " + std::to_string(ri);
            return result;
        }
        res.type = AminoAcidFromThreeLetterCode(canonical);
        res.sequence_number = static_cast<int>(ri + 1);
        res.chain_id = "A";

        // Detect protonation variant from AMBER label
        if (res_labels[ri] != canonical) {
            const AminoAcidType& aat = GetAminoAcidType(res.type);
            for (int vi = 0; vi < static_cast<int>(aat.variants.size()); ++vi) {
                if (res_labels[ri] == aat.variants[vi].name) {
                    res.protonation_variant_index = vi;
                    break;
                }
            }
        }

        protein->AddResidue(std::move(res));
    }

    // Add atoms
    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto atom = std::make_unique<Atom>();
        atom->pdb_atom_name = atom_names[ai];

        if (ai < atomic_numbers.size())
            atom->element = ElementFromAtomicNumber(atomic_numbers[ai]);
        else
            atom->element = ElementFromSymbol(xyz[ai].element);

        // Find residue for this atom
        for (size_t ri = 0; ri < n_residues; ++ri) {
            if (ai >= ranges[ri].start && ai < ranges[ri].end) {
                atom->residue_index = ri;
                break;
            }
        }

        // H parent: in tleap ordering, H follows its parent heavy atom
        if (atom->element == Element::H && ai > 0) {
            for (size_t prev = ai; prev > 0; --prev) {
                size_t pi = prev - 1;
                if (pi < atomic_numbers.size() && atomic_numbers[pi] != 1) {
                    atom->parent_atom_index = pi;
                    break;
                }
            }
        }

        size_t idx = protein->AddAtom(std::move(atom));

        // Update residue atom_indices and backbone cache
        Residue& res = protein->MutableResidueAt(
            protein->AtomAt(idx).residue_index);
        res.atom_indices.push_back(idx);

        // Backbone index cache from atom name
        // PDB LOADING BOUNDARY: string name → typed backbone index
        const std::string& name = atom_names[ai];
        if (name == "N" && res.N == Residue::NONE) res.N = idx;
        else if (name == "CA" && res.CA == Residue::NONE) res.CA = idx;
        else if (name == "C" && res.C == Residue::NONE &&
                 protein->AtomAt(idx).element == Element::C &&
                 res.CA != Residue::NONE) res.C = idx;
        else if (name == "O" && res.O == Residue::NONE) res.O = idx;
        else if ((name == "H" || name == "HN" || name == "H1") &&
                 res.H == Residue::NONE) res.H = idx;
        else if ((name == "HA" || name == "HA2") &&
                 res.HA == Residue::NONE) res.HA = idx;
        else if (name == "CB" && res.CB == Residue::NONE) res.CB = idx;
    }

    // Build positions and create conformation
    std::vector<Vec3> positions;
    positions.reserve(n_atoms);
    for (const auto& a : xyz)
        positions.push_back(Vec3(a.x, a.y, a.z));

    // Set build context
    auto ctx = std::make_unique<ProteinBuildContext>();
    ctx->pdb_source = files.pdb_path;
    ctx->force_field = "ff14SB";
    ctx->protonation_tool = "tleap";
    ctx->prmtop_path = files.prmtop_path;
    ctx->tleap_script_path = files.tleap_script_path;
    protein->SetBuildContext(std::move(ctx));

    // Finalize: backbone indices, bonds, rings — same as PdbFileReader.
    // Must run BEFORE creating the conformation.
    protein->FinalizeConstruction(positions);

    // Create conformation with XYZ positions.
    // These are tleap-protonated AlphaFold structures — predictions, not crystals.
    protein->AddPrediction(std::move(positions), "AlphaFold+tleap");

    result.protein = std::move(protein);
    result.ok = true;
    return result;
}


// ============================================================================
// LoadOrcaRun: dispatcher
// ============================================================================

BuildResult BuildFromOrca(const OrcaRunFiles& files) {
    BuildResult result;

    OperationLog::Scope scope("BuildFromOrca",
        "pdb=" + files.pdb_path + " xyz=" + files.xyz_path);

    if (!fs::exists(files.xyz_path)) {
        result.error = "XYZ not found: " + files.xyz_path;
        return result;
    }

    auto xyz = ReadXyz(files.xyz_path);
    if (xyz.empty()) {
        result.error = "failed to read XYZ: " + files.xyz_path;
        return result;
    }

    // ORCA / mutant paths require an upstream PRMTOP. Missing or
    // unreadable PRMTOP is a hard load error — there is no fall-through
    // to the ff14SB flat table on the ORCA paths.
    if (files.prmtop_path.empty() || !fs::exists(files.prmtop_path)) {
        result.error = "no prmtop available for " + files.pdb_path +
                       ". --orca/--mutant require an upstream PRMTOP "
                       "(provide files.prmtop_path).";
        return result;
    }

    auto internal = LoadWithPrmtop(files, xyz);
    if (!internal.ok) {
        result.error = internal.error;
        return result;
    }
    result.protein = std::move(internal.protein);

    // Charges from prmtop, via the resolver (branch 1 short-circuit).
    AmberSourceConfig source_config;
    source_config.preparation_policy =
        AmberPreparationPolicy::FailOnUnsupportedTerminalVariants;
    std::string charge_err;
    result.charges = ResolveAmberChargeSource(
        *result.protein, result.protein->BuildContext(),
        source_config, charge_err);
    if (!result.charges) {
        result.error = charge_err;
        return result;
    }

    auto& conf = result.protein->Conformation();
    if (!result.protein->PrepareForceFieldCharges(*result.charges, conf, charge_err)) {
        result.error = "charge preparation failed: " + charge_err;
        return result;
    }
    double charge_sum = result.protein->ForceFieldCharges().TotalCharge();
    result.net_charge = static_cast<int>(
        charge_sum + (charge_sum > 0 ? 0.5 : -0.5));

    result.ok = true;

    OperationLog::Info(LogCharges, "BuildFromOrca",
        "loaded " + std::to_string(result.protein->AtomCount()) + " atoms, " +
        std::to_string(result.protein->ResidueCount()) + " residues, " +
        "net_charge=" + std::to_string(result.net_charge));

    return result;
}

}  // namespace nmr
