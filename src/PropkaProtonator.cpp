#include "PropkaProtonator.h"
#include "Protein.h"
#include "OperationLog.h"
#include "RuntimeEnvironment.h"
#include "AminoAcidType.h"

#include <fstream>
#include <sstream>
#include <filesystem>
#include <cstdlib>
#include <cmath>

namespace fs = std::filesystem;

namespace nmr {

// ============================================================================
// PDB LOADING BOUNDARY: write a temp PDB for PROPKA (heavy atoms only).
// PROPKA works on heavy atoms — it does not need or want hydrogens.
// Atom names and residue names are PDB-standard strings here.
// After parsing PROPKA output, we return to typed objects.
// ============================================================================

std::string PropkaProtonator::WriteTempPdb(
        const Protein& protein,
        const ProteinConformation& conf,
        const std::string& dir) {

    std::string path = dir + "/propka_input.pdb";
    std::ofstream out(path);
    if (!out.is_open()) return "";

    int serial = 1;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        for (size_t ai : res.atom_indices) {
            const Atom& atom = protein.AtomAt(ai);
            // Heavy atoms only
            if (atom.element == Element::H) continue;

            Vec3 pos = conf.AtomAt(ai).Position();

            // PDB LOADING BOUNDARY: format atom name for PDB columns 13-16
            char atomField[5];
            if (atom.pdb_atom_name.size() <= 3)
                std::snprintf(atomField, sizeof(atomField), " %-3s",
                              atom.pdb_atom_name.c_str());
            else
                std::snprintf(atomField, sizeof(atomField), "%-4s",
                              atom.pdb_atom_name.c_str());

            // PDB LOADING BOUNDARY: residue name from three-letter code
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
// PDB LOADING BOUNDARY: parse PROPKA .pka output.
// PROPKA writes lines like "   ASP  43 A     3.94" in the summary section.
// We accept ALL predictions — PROPKA is the authority for which residues
// are titratable and what their pKa is.
// ============================================================================

std::vector<PkaResult> PropkaProtonator::ParsePkaFile(const std::string& path) {
    std::vector<PkaResult> results;
    std::ifstream in(path);
    if (!in.is_open()) return results;

    bool in_summary = false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.find("SUMMARY OF THIS PREDICTION") != std::string::npos) {
            in_summary = true;
            continue;
        }
        if (in_summary && line.find("-----") != std::string::npos)
            continue;
        if (in_summary && line.empty()) break;

        if (in_summary && line.size() > 20) {
            // PDB LOADING BOUNDARY: parse "   ASP  43 A     3.94"
            std::istringstream iss(line);
            std::string res_type;
            int res_num;
            std::string chain;
            double pka;
            if (iss >> res_type >> res_num >> chain >> pka) {
                results.push_back({res_type, res_num, chain, pka});
            }
        }
    }
    return results;
}


// ============================================================================
// Henderson-Hasselbalch: pKa + pH → protonation decision.
//
// For acids (ASP, GLU, CYS, TYR): pKa > pH → protonated (neutral)
// For bases (HIS, LYS): pKa > pH → protonated (charged)
//
// PDB LOADING BOUNDARY: PROPKA output uses PDB residue names.
// We match to our typed Protein residues by sequence number and amino
// acid type, then store variant_index (typed), not strings.
// ============================================================================

ProtonationState PropkaProtonator::ApplyHendersonHasselbalch(
        const Protein& protein,
        const std::vector<PkaResult>& pkas,
        double pH) {

    ProtonationState state(
        "propka_pH" + std::to_string(pH),
        pH,
        ProtonationTool::PROPKA,
        "PROPKA 3.5.1");

    for (const auto& pka : pkas) {
        // PDB LOADING BOUNDARY: match PROPKA's string output to our residues.
        // After this block, everything is typed.
        size_t res_idx = SIZE_MAX;
        for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
            const Residue& res = protein.ResidueAt(ri);
            if (res.sequence_number != pka.residue_number) continue;

            std::string code = ThreeLetterCodeForAminoAcid(res.type);
            if (code == pka.residue_type) {
                res_idx = ri;
                break;
            }
            // HIS matches HIS/HID/HIE/HIP from PROPKA
            if (res.type == AminoAcid::HIS && pka.residue_type == "HIS") {
                res_idx = ri;
                break;
            }
        }
        if (res_idx == SIZE_MAX) continue;

        // Now typed: we have the Residue and its AminoAcidType
        const Residue& res = protein.ResidueAt(res_idx);
        bool protonated = pka.pKa > pH;

        ResidueProtonation decision;
        decision.residue_index = res_idx;
        decision.amino_acid = res.type;
        decision.pKa = pka.pKa;

        // Map H-H decision to variant_index using AminoAcidType::variants
        switch (res.type) {
            case AminoAcid::ASP:
                // protonated → ASH (variant 0), deprotonated → ASP (default, -1)
                decision.variant_index = protonated ? 0 : -1;
                decision.is_charged = !protonated;
                break;

            case AminoAcid::GLU:
                // protonated → GLH (variant 0), deprotonated → GLU (default, -1)
                decision.variant_index = protonated ? 0 : -1;
                decision.is_charged = !protonated;
                break;

            case AminoAcid::HIS:
                // protonated → HIP (variant 2, +1), deprotonated → HIE (variant 1, 0)
                // PROPKA cannot distinguish HID from HIE; default to HIE when neutral
                decision.variant_index = protonated ? 2 : 1;
                decision.is_charged = protonated;
                break;

            case AminoAcid::CYS:
                // protonated → CYS (default, -1), deprotonated → CYM (variant 1)
                // Note: CYX (disulfide, variant 0) is detected by bond, not pKa
                decision.variant_index = protonated ? -1 : 1;
                decision.is_charged = !protonated;
                break;

            case AminoAcid::LYS:
                // protonated → LYS (default, +1), deprotonated → LYN (variant 0, 0)
                decision.variant_index = protonated ? -1 : 0;
                decision.is_charged = protonated;
                break;

            case AminoAcid::ARG:
                // protonated → ARG (default, +1), deprotonated → ARN (variant 0, 0)
                // pKa ~12.5, almost always protonated at physiological pH
                decision.variant_index = protonated ? -1 : 0;
                decision.is_charged = protonated;
                break;

            case AminoAcid::TYR:
                // protonated → TYR (default), deprotonated → TYM (variant 0, -1)
                decision.variant_index = protonated ? -1 : 0;
                decision.is_charged = !protonated;
                break;

            default:
                continue;  // PROPKA may report N-term/C-term; skip
        }

        state.AddResidue(decision);
    }

    return state;
}


// ============================================================================
// PredictPka: call PROPKA, return raw pKa values.
// ============================================================================

std::vector<PkaResult> PropkaProtonator::PredictPka(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) {

    OperationLog::Scope scope("PropkaProtonator::PredictPka",
        "residues=" + std::to_string(protein.ResidueCount()));

    // Create unique temp directory
    std::string dir = RuntimeEnvironment::TempFilePath(
        "propka", std::to_string(protein.ResidueCount()));
    fs::create_directories(dir);

    // Write heavy-atom PDB
    std::string pdb_path = WriteTempPdb(protein, conf, dir);
    if (pdb_path.empty()) {
        error_out = "cannot write temp PDB to " + dir;
        std::error_code ec;
        fs::remove_all(dir, ec);
        return {};
    }

    // Call propka3
    // TODO: propka path should be a parameter when reintegrated into BuildFromPdb
    std::string cmd = "cd " + dir + " && propka3" +
                      " --quiet " + pdb_path + " 2>/dev/null";
    int rc = std::system(cmd.c_str());

    // PROPKA writes <stem>.pka in the working directory
    std::string stem = fs::path(pdb_path).stem().string();
    std::string pka_path = dir + "/" + stem + ".pka";

    if (!fs::exists(pka_path)) {
        error_out = "propka3 produced no output (rc=" + std::to_string(rc) +
                    "), expected " + pka_path;
        std::error_code ec;
        fs::remove_all(dir, ec);
        return {};
    }

    auto results = ParsePkaFile(pka_path);

    // Clean up
    std::error_code ec;
    fs::remove_all(dir, ec);

    OperationLog::Info(LogProtonation, "PropkaProtonator::PredictPka",
        std::to_string(results.size()) + " pKa predictions");

    return results;
}


// ============================================================================
// Protonate: the full pipeline — predict pKa, apply H-H, return state.
// ============================================================================

ProtonationResult PropkaProtonator::Protonate(
        const Protein& protein,
        const ProteinConformation& conf,
        double pH) {

    OperationLog::Scope scope("PropkaProtonator::Protonate",
        "pH=" + std::to_string(pH));

    ProtonationResult result;

    auto pkas = PredictPka(protein, conf, result.error);
    if (pkas.empty() && !result.error.empty()) {
        result.ok = false;
        return result;
    }

    result.state = ApplyHendersonHasselbalch(protein, pkas, pH);
    result.ok = true;

    OperationLog::Info(LogProtonation, "PropkaProtonator::Protonate",
        result.state.Describe());

    return result;
}

}  // namespace nmr
