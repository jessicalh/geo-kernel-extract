#include "KamlProtonator.h"
#include "PropkaProtonator.h"  // reuse ApplyHendersonHasselbalch and WriteTempPdb pattern
#include "Protein.h"
#include "OperationLog.h"
#include "RuntimeEnvironment.h"

#include <fstream>
#include <sstream>
#include <filesystem>
#include <cstdlib>
#include <cstdio>
#include <cmath>

namespace fs = std::filesystem;

namespace nmr {

// ============================================================================
// PDB LOADING BOUNDARY: write temp PDB for KaML (heavy atoms only).
// Same logic as PropkaProtonator::WriteTempPdb.
// ============================================================================

std::string KamlProtonator::WriteTempPdb(
        const Protein& protein,
        const ProteinConformation& conf,
        const std::string& dir) {

    std::string path = dir + "/kaml_input.pdb";
    std::ofstream out(path);
    if (!out.is_open()) return "";

    int serial = 1;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        for (size_t ai : res.atom_indices) {
            const Atom& atom = protein.AtomAt(ai);
            if (atom.element == Element::H) continue;

            Vec3 pos = conf.AtomAt(ai).Position();

            char atomField[5];
            if (atom.pdb_atom_name.size() <= 3)
                std::snprintf(atomField, sizeof(atomField), " %-3s",
                              atom.pdb_atom_name.c_str());
            else
                std::snprintf(atomField, sizeof(atomField), "%-4s",
                              atom.pdb_atom_name.c_str());

            std::string resname = ThreeLetterCodeForAminoAcid(res.type);
            char chain = res.chain_id.empty() ? 'A' : res.chain_id[0];

            char line[82];
            std::snprintf(line, sizeof(line),
                "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                serial++, atomField, resname.c_str(), chain,
                res.sequence_number, pos.x(), pos.y(), pos.z());
            out << line;
        }
    }
    out << "END\n";
    return path;
}


// ============================================================================
// PDB LOADING BOUNDARY: parse KaML _pka.csv output.
// Format: "Residue, pKa\nASP-43,3.245\nHIS-24,6.789\n..."
// ============================================================================

std::vector<PkaResult> KamlProtonator::ParsePkaCsv(const std::string& path) {
    std::vector<PkaResult> results;
    std::ifstream in(path);
    if (!in.is_open()) return results;

    std::string line;
    // Skip header
    if (!std::getline(in, line)) return results;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        // Parse "ASP-43,3.245" or "HIS-24,6.789"
        auto comma = line.find(',');
        if (comma == std::string::npos) continue;

        std::string residue_part = line.substr(0, comma);
        std::string pka_str = line.substr(comma + 1);

        auto dash = residue_part.find('-');
        if (dash == std::string::npos) continue;

        PkaResult r;
        r.residue_type = residue_part.substr(0, dash);
        r.residue_number = 0;
        std::sscanf(residue_part.substr(dash + 1).c_str(), "%d", &r.residue_number);
        r.chain_id = "A";
        r.pKa = 0.0;
        std::sscanf(pka_str.c_str(), "%lf", &r.pKa);

        if (r.residue_number > 0 && !r.residue_type.empty()) {
            results.push_back(r);
        }
    }
    return results;
}


// ============================================================================
// PredictPka: call KaML, return raw pKa values.
// ============================================================================

std::vector<PkaResult> KamlProtonator::PredictPka(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) {

    OperationLog::Scope scope("KamlProtonator::PredictPka",
        "residues=" + std::to_string(protein.ResidueCount()));

    // TODO: KaML path should be a parameter when reintegrated into BuildFromPdb
    const std::string kaml_path = "KaML-CBtree.py";
    if (!fs::exists(kaml_path)) {
        error_out = "KaML not found at " + kaml_path;
        return {};
    }

    // KaML creates intermediate files in CWD — must run from temp dir
    std::string dir = RuntimeEnvironment::TempFilePath(
        "kaml", std::to_string(protein.ResidueCount()));
    fs::create_directories(dir);

    std::string pdb_path = WriteTempPdb(protein, conf, dir);
    if (pdb_path.empty()) {
        error_out = "cannot write temp PDB to " + dir;
        std::error_code ec;
        fs::remove_all(dir, ec);
        return {};
    }

    // KaML writes <stem>_pka.csv in CWD
    std::string pdb_stem = fs::path(pdb_path).stem().string();
    std::string cmd = "cd " + dir + " && python3 " + kaml_path +
                      " " + pdb_path + " 2>/dev/null";
    int rc = std::system(cmd.c_str());

    std::string pka_path = dir + "/" + pdb_stem + "_pka.csv";

    if (!fs::exists(pka_path)) {
        error_out = "KaML produced no output (rc=" + std::to_string(rc) +
                    "), expected " + pka_path;
        std::error_code ec;
        fs::remove_all(dir, ec);
        return {};
    }

    auto results = ParsePkaCsv(pka_path);

    // Clean up
    std::error_code ec;
    fs::remove_all(dir, ec);

    OperationLog::Info(LogProtonation, "KamlProtonator::PredictPka",
        std::to_string(results.size()) + " pKa predictions");

    return results;
}


// ============================================================================
// Protonate: predict pKa + Henderson-Hasselbalch.
// Same H-H logic as PropkaProtonator — the chemistry is the same,
// only the pKa source differs.
// ============================================================================

ProtonationResult KamlProtonator::Protonate(
        const Protein& protein,
        const ProteinConformation& conf,
        double pH) {

    OperationLog::Scope scope("KamlProtonator::Protonate",
        "pH=" + std::to_string(pH));

    ProtonationResult result;

    auto pkas = PredictPka(protein, conf, result.error);
    if (pkas.empty() && !result.error.empty()) {
        result.ok = false;
        return result;
    }

    // Apply Henderson-Hasselbalch — same logic as PROPKA
    // (reuse PropkaProtonator::ApplyHendersonHasselbalch if it were public,
    // but it's private. Duplicate the H-H logic here — it's the same chemistry.)
    ProtonationState state(
        "kaml_pH" + std::to_string(pH),
        pH,
        ProtonationTool::KaML,
        "KaML-CBTrees");

    for (const auto& pka : pkas) {
        // PDB LOADING BOUNDARY: match KaML's string output to our residues
        size_t res_idx = SIZE_MAX;
        for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
            const Residue& res = protein.ResidueAt(ri);
            if (res.sequence_number != pka.residue_number) continue;
            std::string code = ThreeLetterCodeForAminoAcid(res.type);
            if (code == pka.residue_type) { res_idx = ri; break; }
            if (res.type == AminoAcid::HIS &&
                (pka.residue_type == "HIS" || pka.residue_type == "HIE" ||
                 pka.residue_type == "HID" || pka.residue_type == "HIP")) {
                res_idx = ri; break;
            }
        }
        if (res_idx == SIZE_MAX) continue;

        const Residue& res = protein.ResidueAt(res_idx);
        bool protonated = pka.pKa > pH;

        ResidueProtonation decision;
        decision.residue_index = res_idx;
        decision.amino_acid = res.type;
        decision.pKa = pka.pKa;

        switch (res.type) {
            case AminoAcid::ASP:
                decision.variant_index = protonated ? 0 : -1;
                decision.is_charged = !protonated;
                break;
            case AminoAcid::GLU:
                decision.variant_index = protonated ? 0 : -1;
                decision.is_charged = !protonated;
                break;
            case AminoAcid::HIS:
                decision.variant_index = protonated ? 2 : 1;
                decision.is_charged = protonated;
                break;
            case AminoAcid::CYS:
                decision.variant_index = protonated ? -1 : 1;
                decision.is_charged = !protonated;
                break;
            case AminoAcid::LYS:
                decision.variant_index = protonated ? -1 : 0;
                decision.is_charged = protonated;
                break;
            case AminoAcid::ARG:
                decision.variant_index = protonated ? -1 : 0;
                decision.is_charged = protonated;
                break;
            case AminoAcid::TYR:
                decision.variant_index = protonated ? -1 : 0;
                decision.is_charged = !protonated;
                break;
            default:
                continue;
        }

        state.AddResidue(decision);
    }

    result.state = std::move(state);
    result.ok = true;

    OperationLog::Info(LogProtonation, "KamlProtonator::Protonate",
        result.state.Describe());

    return result;
}

}  // namespace nmr
