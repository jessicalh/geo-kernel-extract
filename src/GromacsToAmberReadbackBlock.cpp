#include "GromacsToAmberReadbackBlock.h"

#include <fstream>
#include <sstream>
#include <string>

#include <nlohmann/json.hpp>

namespace nmr {

GromacsToAmberReadbackBlock ParseTopolTopReadback(
        const std::string& topol_top_path,
        std::string& error_out) {

    error_out.clear();
    GromacsToAmberReadbackBlock block;
    block.topol_top_path = topol_top_path;

    std::ifstream in(topol_top_path);
    if (!in.is_open()) {
        error_out = "ParseTopolTopReadback: failed to open " + topol_top_path;
        return block;
    }

    std::string line;
    while (std::getline(in, line)) {
        size_t i = line.find_first_not_of(" \t");
        if (i == std::string::npos || line[i] != ';') continue;
        ++i;
        while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) ++i;
        if (line.compare(i, 7, "residue") != 0) continue;
        i += 7;

        std::istringstream iss(line.substr(i));
        int seqid = 0;
        std::string name, rtp_kw, rtp, q_kw;
        double q = 0.0;
        if (!(iss >> seqid)) continue;
        if (!(iss >> name)) continue;
        if (!(iss >> rtp_kw) || rtp_kw != "rtp") continue;
        if (!(iss >> rtp)) continue;
        if (!(iss >> q_kw) || q_kw != "q") continue;
        if (!(iss >> q)) continue;

        if (seqid <= 0) continue;

        GromacsToAmberReadbackBlock::ResidueEntry entry;
        entry.tpr_name    = name;
        entry.rtp         = rtp;
        entry.source_line = line;
        entry.charge_q    = q;

        const std::string base = BaseFfPortNameFromGromacsRtp(rtp);
        entry.canonical_three = CanonicalThreeLetterFromGromacsRtp(rtp);
        if (!entry.canonical_three.empty()) {
            entry.aa = AminoAcidFromThreeLetterCode(entry.canonical_three);
            entry.variant_index =
                VariantIndexFromForceFieldName(entry.aa, base);
        }

        // Gaps stay default-constructed so consumers can detect missing rows.
        const size_t idx = static_cast<size_t>(seqid - 1);
        if (block.residues.size() <= idx) {
            block.residues.resize(idx + 1);
        }
        block.residues[idx] = std::move(entry);
    }

    if (block.residues.empty()) {
        error_out = "ParseTopolTopReadback: no rtp comment lines found in " +
                    topol_top_path;
        return block;
    }

    for (const auto& e : block.residues) {
        const std::string base = BaseFfPortNameFromGromacsRtp(e.rtp);
        if (!e.tpr_name.empty() && !base.empty() && e.tpr_name != base) {
            ++block.n_port_label_translations;
        }
        if (base == "CYX") {
            ++block.n_disulfide_residues;
        }
    }

    return block;
}


bool EmitGromacsToAmberReadbackBlockJson(
        const GromacsToAmberReadbackBlock& block,
        const std::string& output_path,
        std::string& error_out) {

    error_out.clear();

    // Preserve insertion order for stable audit output.
    nlohmann::ordered_json j;
    j["schema_version"]              = 1;
    j["topol_top_path"]              = block.topol_top_path;
    j["n_residues"]                  = block.residues.size();
    j["n_port_label_translations"]   = block.n_port_label_translations;
    j["n_disulfide_residues"]        = block.n_disulfide_residues;

    nlohmann::ordered_json residues = nlohmann::ordered_json::array();
    for (size_t i = 0; i < block.residues.size(); ++i) {
        const auto& e = block.residues[i];
        if (e.tpr_name.empty() && e.rtp.empty()) continue;  // skip gaps
        residues.push_back({
            {"index",           i},
            {"tpr_name",        e.tpr_name},
            {"rtp",             e.rtp},
            {"canonical_three", e.canonical_three},
            {"variant_index",   e.variant_index},
            {"charge_q",        e.charge_q},
        });
    }
    j["residues"] = std::move(residues);

    std::ofstream out(output_path);
    if (!out.is_open()) {
        error_out = "EmitGromacsToAmberReadbackBlockJson: failed to open " +
                    output_path;
        return false;
    }
    out << j.dump(2, ' ', false,
                  nlohmann::ordered_json::error_handler_t::replace) << "\n";
    if (!out.good()) {
        error_out = "EmitGromacsToAmberReadbackBlockJson: write failure on " +
                    output_path;
        return false;
    }
    return true;
}

}  // namespace nmr
