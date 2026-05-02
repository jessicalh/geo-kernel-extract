#pragma once
//
// GromacsToAmberReadbackBlock: load-time compiler trace that reads back
// the chemistry decisions GROMACS pdb2gmx made (during prep) and applies
// them as typed facts on the existing object model. NOT live state.
//
// Lifecycle: the block lives on the loader's stack. ParseTopolTopReadback
// builds it from the topol.top rtp comment lines; FullSystemReader::
// BuildProtein consumes it during the per-residue loop to set typed fields
// (Residue.type, Residue.protonation_variant_index); EmitJson optionally
// writes an audit file alongside the topology; the block goes out of
// scope at the end of BuildProtein. Calculators NEVER see it.
//
// Why a block, not residue-name aliasing: GROMACS rewrites .name to FF-port
// labels (HISH/HISD/HISE/CYS-with-HG-stripped) but records the canonical
// AMBER chemistry decision (HID/HIE/HIP/CYX/...) in the topol.top "; residue
// N <name> rtp <rtp> q <q>" comment line. Reading the rtp directly is
// structurally correct: GROMACS is the chemistry authority, we read what
// it decided, we never re-decide. NamingRegistry expansion (alias HISH →
// HIS) was rejected because it hides what GROMACS decided and pollutes
// the stable naming vocabulary.
//
// This block does NOT add any new string field to model objects. After
// Apply, all chemistry information lives in typed slots:
//   Residue.type                       (AminoAcid enum)
//   Residue.protonation_variant_index  (int, per AminoAcidType.h contract)
//   Residue.terminal_state             (ResidueTerminalState enum;
//                                       set by Protein::ResolveResidueTerminalStates
//                                       from chain position, not from rtp;
//                                       cap-residue extensions are out of scope)
//   CovalentTopology bonds + categories (typed, set by Protein::FinalizeConstruction)
//   ForceFieldChargeTable               (typed)
//   LegacyAmberInvariants → LegacyAmberTopology (typed)
//
// Strings on the block (tpr_name, rtp, source_line) exist only for the
// duration of the load and the JSON audit emission. They are never copied
// onto Residue, Atom, or any persistent typed object.
//
// Design doc: spec/plan/gromacs-to-amber-readback-block-design-2026-05-02.md
// Companion memory: feedback_readback_block_is_a_compiler_trace
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
        // Source-side strings (transient, used only during parse + audit;
        // NEVER copied onto model objects).
        std::string tpr_name;        // GROMACS .name field, e.g. "HISH"
        std::string rtp;             // canonical AMBER rtp, e.g. "HIP" or "NHIP"
        std::string source_line;     // verbatim comment line, for audit JSON

        // Typed resolution (this is what BuildProtein consumes).
        std::string canonical_three; // e.g. "HIS" (after stripping N/C prefix)
        AminoAcid   aa = AminoAcid::Unknown;
        int         variant_index = -1;     // -1 = no variant / canonical-charged form

        // Per-comment charge as recorded by pdb2gmx. Cross-check signal only;
        // the authoritative per-atom charges come from the TPR via
        // PreloadedChargeSource.
        double      charge_q = 0.0;
    };

    std::vector<ResidueEntry> residues;

    // Audit fields.
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


// Emit the block as JSON alongside the topology. Returns true on success.
// JSON is for debugging / methods text / reviewer spot-check; calculators
// do NOT consume this file at runtime. Regeneratable from the authority
// files; emitting at extraction time is convenience, not load-bearing.
bool EmitGromacsToAmberReadbackBlockJson(
    const GromacsToAmberReadbackBlock& block,
    const std::string& output_path,
    std::string& error_out);

}  // namespace nmr
