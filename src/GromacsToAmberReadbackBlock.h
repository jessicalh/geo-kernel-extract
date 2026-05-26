#pragma once
//
// Parses GROMACS topol.top rtp comment lines and exposes the pdb2gmx
// residue decisions as typed residue facts for load-time use. The string
// fields are transient parse/audit data, not model state.
//

#include "AminoAcidType.h"
#include "Residue.h"

#include <cstddef>
#include <string>
#include <vector>

namespace nmr {

struct GromacsToAmberReadbackBlock {
    // One entry per topol.top "; residue N <name> rtp <rtp> q <q>" line.
    // Indexed by 0-based residue index (1-based topol seqid minus 1).
    // The vector covers all residues the topol.top declares; the
    // consumer (BuildProtein) reads only the protein-residue range and
    // ignores entries past the protein count (water/ion residues that
    // share the same rtp comment style).
    struct ResidueEntry {
        // Source-side strings used only during parse and audit emission.
        std::string tpr_name;        // GROMACS .name field, e.g. "HISH"
        std::string rtp;             // canonical AMBER rtp, e.g. "HIP" or "NHIP"
        std::string source_line;     // verbatim comment line, for audit JSON

        // Typed resolution consumed by FullSystemReader::BuildProtein.
        std::string canonical_three; // e.g. "HIS" (after stripping N/C prefix)
        AminoAcid   aa = AminoAcid::Unknown;
        int         variant_index = -1;     // -1 = no variant / canonical-charged form

        // Per-comment charge recorded by pdb2gmx; per-atom charges come
        // from the TPR via PreloadedChargeSource.
        double      charge_q = 0.0;
    };

    std::vector<ResidueEntry> residues;

    std::string topol_top_path;
    int n_port_label_translations = 0;  // count where tpr_name != base(rtp)
    int n_disulfide_residues = 0;       // count of CYX occurrences
};


// Parse the rtp comment lines from a topol.top file.
// On success: error_out empty, returns a populated block.
// On failure: error_out non-empty, returns a block with empty residues.
//
// Robust to whitespace variation: matches both "; residue   3 HISH rtp HIP q +1.0"
// and "; residue 3 HISH rtp HIP q +1.0" patterns. Lines that don't match the
// rtp-comment shape are skipped silently (e.g., section headers, [ atoms ]).
GromacsToAmberReadbackBlock ParseTopolTopReadback(
    const std::string& topol_top_path,
    std::string& error_out);


// Audit JSON only; runtime code does not consume the emitted file.
bool EmitGromacsToAmberReadbackBlockJson(
    const GromacsToAmberReadbackBlock& block,
    const std::string& output_path,
    std::string& error_out);

}  // namespace nmr
