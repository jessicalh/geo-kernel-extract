#include "GromacsToAmberReadbackBlock.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

namespace nmr {

namespace {

// JSON string escape: backslash, double-quote, common control characters.
// Path strings and source-line strings are the only inputs; no Unicode
// escape needs (GROMACS topol.top is ASCII).
std::string EscapeJsonString(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 2);
    for (char c : s) {
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                    out += buf;
                } else {
                    out += c;
                }
        }
    }
    return out;
}

}  // namespace


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
        // Locate "; residue " (allow leading whitespace, allow no space after ';').
        size_t i = line.find_first_not_of(" \t");
        if (i == std::string::npos || line[i] != ';') continue;
        ++i;
        // Skip whitespace after ';'.
        while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) ++i;
        // Expect literal "residue" keyword.
        if (line.compare(i, 7, "residue") != 0) continue;
        i += 7;

        // Tokenise the rest: <seqid> <name> rtp <rtp> q <q>
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

        if (seqid <= 0) continue;  // defensive: malformed seqid

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

        // Place at 0-based index. Vector grows as needed; gaps stay
        // default-constructed (Unknown aa, empty strings) so consumers
        // can detect missing entries.
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

    // Audit counts.
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
    std::ofstream out(output_path);
    if (!out.is_open()) {
        error_out = "EmitGromacsToAmberReadbackBlockJson: failed to open " +
                    output_path;
        return false;
    }

    out << "{\n";
    out << "  \"schema_version\": 1,\n";
    out << "  \"topol_top_path\": \""
        << EscapeJsonString(block.topol_top_path) << "\",\n";
    out << "  \"n_residues\": " << block.residues.size() << ",\n";
    out << "  \"n_port_label_translations\": "
        << block.n_port_label_translations << ",\n";
    out << "  \"n_disulfide_residues\": "
        << block.n_disulfide_residues << ",\n";
    out << "  \"residues\": [\n";

    bool first = true;
    for (size_t i = 0; i < block.residues.size(); ++i) {
        const auto& e = block.residues[i];
        if (e.tpr_name.empty() && e.rtp.empty()) continue;  // skip gaps
        if (!first) out << ",\n";
        first = false;
        out << "    {";
        out << "\"index\": " << i;
        out << ", \"tpr_name\": \"" << EscapeJsonString(e.tpr_name) << "\"";
        out << ", \"rtp\": \"" << EscapeJsonString(e.rtp) << "\"";
        out << ", \"canonical_three\": \""
            << EscapeJsonString(e.canonical_three) << "\"";
        out << ", \"variant_index\": " << e.variant_index;
        out << ", \"charge_q\": " << e.charge_q;
        out << "}";
    }
    out << "\n  ]\n";
    out << "}\n";

    if (!out.good()) {
        error_out = "EmitGromacsToAmberReadbackBlockJson: write failure on " +
                    output_path;
        return false;
    }
    return true;
}

}  // namespace nmr
