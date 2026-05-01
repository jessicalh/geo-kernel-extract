#include "AmberLeapInput.h"
#include "AminoAcidType.h"
#include "Atom.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "Types.h"

#include <array>
#include <cstdio>
#include <ostream>
#include <set>
#include <unordered_map>

namespace nmr {
namespace amber_leap {

namespace {

// Resolve the AMBER residue name from typed state.
// Mirrors the lookup-name choice in ParamFileChargeSource so the same
// triple (terminal_state, ff_resname, atom_name) drives both the verdict
// and the generated PDB.
std::string AmberResidueNameFor(const Residue& res, bool is_disulfide_cys) {
    if (res.type == AminoAcid::CYS && is_disulfide_cys) {
        return "CYX";
    }

    if (res.protonation_variant_index < 0) {
        if (res.type == AminoAcid::HIS) return "HIE";
        return ThreeLetterCodeForAminoAcid(res.type);
    }

    const AminoAcidType& aa_type = res.AminoAcidInfo();
    if (res.protonation_variant_index >=
            static_cast<int>(aa_type.variants.size())) {
        return ThreeLetterCodeForAminoAcid(res.type);
    }
    return aa_type.variants[res.protonation_variant_index].name;
}

// PDB ATOM record atom-name column rule for 1-letter elements (the only
// elements present in protein chemistry from ff14SB):
//   length 4: occupy columns 13–16 fully.
//   length < 4: cols 14–16, col 13 is a space.
std::string FormatAtomNameColumn(const std::string& name) {
    char buf[5] = {' ', ' ', ' ', ' ', '\0'};
    if (name.size() >= 4) {
        for (int i = 0; i < 4; ++i) buf[i] = name[i];
    } else {
        for (size_t i = 0; i < name.size(); ++i) buf[i + 1] = name[i];
    }
    return std::string(buf, 4);
}

// ---------------------------------------------------------------------------
// Capping support (step 6).
//
// Under UseCappedFragmentsForUnsupportedTerminalVariants, terminal
// ASH / CYM / GLH / LYN — variants ff14SB ships only as INTERNAL — are
// modelled as INTERNAL by inserting an ACE acetyl cap at the N end
// and/or an NME N-methylamide cap at the C end. tleap then loads the
// original residue with its INTERNAL template; cap atoms are dropped
// during atom mapping (ResidueAmberMapping marks their PRMTOP indices
// NONE_FOR_CAP).
//
// Cap positions are physically reasonable but ideal (offset from the
// connecting backbone atom). Charges are fixed by ff14SB templates and
// independent of position; the only requirement is non-overlapping
// atoms tleap will accept.
//
// Atom names are AMBER library names exactly:
//   ACE: HH31 HH32 HH33 CH3 C O           (N-terminal cap)
//   NME: N H CH3 HH31 HH32 HH33           (C-terminal cap)
// ---------------------------------------------------------------------------

struct CapAtomTemplate {
    const char* name;
    Element     element;
    Vec3        offset;    // displacement from the connecting backbone atom
};

static const std::array<CapAtomTemplate, 6> kAceTemplate = {{
    // Negative-x from the original residue's N (i.e. before it in space).
    {"HH31", Element::H, Vec3(-3.5, +1.0, 0.0)},
    {"HH32", Element::H, Vec3(-3.5, -0.5, +0.87)},
    {"HH33", Element::H, Vec3(-3.5, -0.5, -0.87)},
    {"CH3",  Element::C, Vec3(-2.5, 0.0, 0.0)},
    {"C",    Element::C, Vec3(-1.33, 0.0, 0.0)},  // ~peptide bond length to N
    {"O",    Element::O, Vec3(-1.33, +1.23, 0.0)},
}};

static const std::array<CapAtomTemplate, 6> kNmeTemplate = {{
    // Positive-x from the original residue's C (i.e. after it in space).
    {"N",    Element::N, Vec3(+1.33, 0.0, 0.0)},  // ~peptide bond length from C
    {"H",    Element::H, Vec3(+1.33, +1.0, 0.0)},
    {"CH3",  Element::C, Vec3(+2.5, 0.0, 0.0)},
    {"HH31", Element::H, Vec3(+3.5, +1.0, 0.0)},
    {"HH32", Element::H, Vec3(+3.5, -0.5, +0.87)},
    {"HH33", Element::H, Vec3(+3.5, -0.5, -0.87)},
}};

bool VariantIsCappable(const Residue& res) {
    // Only ASH, CYM, GLH, LYN — the ff14SB INTERNAL variants — can be
    // rescued by capping. Other unsupported variants (TYM, ARN) lack
    // INTERNAL templates too; capping doesn't help, so we refuse.
    if (res.protonation_variant_index < 0) return false;
    switch (res.type) {
        case AminoAcid::ASP:
            return res.protonation_variant_index == 0;  // ASH
        case AminoAcid::CYS:
            return res.protonation_variant_index == 1;  // CYM
        case AminoAcid::GLU:
            return res.protonation_variant_index == 0;  // GLH
        case AminoAcid::LYS:
            return res.protonation_variant_index == 0;  // LYN
        default:
            return false;
    }
}

void WriteAtomRecord(std::ostream& out,
                     int serial,
                     const std::string& atom_name,
                     const std::string& res_name,
                     const std::string& chain_id,
                     int seq_num,
                     const std::string& insertion_code,
                     const Vec3& pos,
                     Element element) {
    char line[88];
    const std::string aname_col = FormatAtomNameColumn(atom_name);
    const char chain_char = chain_id.empty() ? 'A' : chain_id.front();
    const char icode_char =
        insertion_code.empty() ? ' ' : insertion_code.front();
    const std::string esym = SymbolForElement(element);

    std::snprintf(
        line, sizeof(line),
        "ATOM  %5d %4s %3s %1c%4d%1c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s",
        serial,
        aname_col.c_str(),
        res_name.c_str(),
        chain_char,
        seq_num,
        icode_char,
        pos.x(), pos.y(), pos.z(),
        1.00, 0.00,
        esym.c_str());
    out << line << "\n";
}

}  // namespace


// Build per-residue capping decisions from the verdict + policy.
// Returns sets of extractor residue indices that need ACE prepended
// (NTERM cap) and / or NME appended (CTERM cap).
struct CappingPlan {
    std::set<size_t> ace_before;   // extractor residue indices
    std::set<size_t> nme_after;    // extractor residue indices
    std::vector<std::string> applied_descriptions;  // "NTERM ASH residue 1 (chain A) capped with ACE"
};

CappingPlan BuildCappingPlan(const Protein& protein,
                              AmberPreparationPolicy policy,
                              const AmberFlatTableCoverageVerdict& verdict) {
    CappingPlan plan;
    if (policy != AmberPreparationPolicy::
            UseCappedFragmentsForUnsupportedTerminalVariants) {
        return plan;
    }

    for (const auto& f : verdict.failures) {
        if (f.kind != AmberFlatTableCoverageKind::UnsupportedTerminalVariant) {
            continue;
        }
        // Locate the residue from sequence_number + chain_id.
        for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
            const Residue& res = protein.ResidueAt(ri);
            if (res.sequence_number != f.residue_sequence_number) continue;
            if (res.chain_id != f.chain_id) continue;
            if (!VariantIsCappable(res)) continue;

            const bool n_end = (f.terminal_token == "NTERM" ||
                                f.terminal_token == "NCTERM");
            const bool c_end = (f.terminal_token == "CTERM" ||
                                f.terminal_token == "NCTERM");
            if (n_end) {
                plan.ace_before.insert(ri);
                plan.applied_descriptions.push_back(
                    "NTERM " + f.ff_residue_name + " residue " +
                    std::to_string(f.residue_sequence_number) +
                    " (chain " + f.chain_id + ") capped with ACE");
            }
            if (c_end) {
                plan.nme_after.insert(ri);
                plan.applied_descriptions.push_back(
                    "CTERM " + f.ff_residue_name + " residue " +
                    std::to_string(f.residue_sequence_number) +
                    " (chain " + f.chain_id + ") capped with NME");
            }
            break;
        }
    }
    return plan;
}

void EmitCapResidue(std::ostream& pdb_out,
                    int& serial,
                    const std::string& cap_resname,
                    const std::array<CapAtomTemplate, 6>& templ,
                    const Vec3& anchor_position,
                    const std::string& chain_id,
                    int seq_num,
                    const std::string& insertion_code) {
    for (const auto& t : templ) {
        const Vec3 pos = anchor_position + t.offset;
        WriteAtomRecord(pdb_out, serial, t.name, cap_resname,
                        chain_id, seq_num, insertion_code, pos, t.element);
        ++serial;
    }
}

void GenerateAmberPdb(const Protein& protein,
                      const ProteinConformation& conf,
                      AmberPreparationPolicy policy,
                      const AmberFlatTableCoverageVerdict& verdict,
                      std::ostream& pdb_out,
                      ResidueAmberMapping& map_out) {
    map_out.extractor_index_for_prmtop_residue.clear();
    map_out.extractor_index_for_prmtop_residue.reserve(
        protein.ResidueCount() + 8);

    // Disulfide detection by direct SG-SG distance check on typed CYS
    // residues. This is the AmberTools-tutorial methodology and does
    // not depend on CovalentTopology's bond perception.
    std::set<size_t> cyx_residues;
    auto disulfide_pairs = DetectDisulfides(protein, conf);
    for (const auto& pair : disulfide_pairs) {
        cyx_residues.insert(pair.first);
        cyx_residues.insert(pair.second);
    }

    // Step 6: capping plan — which residues get ACE / NME inserts
    // around them under UseCappedFragmentsForUnsupportedTerminalVariants.
    const CappingPlan capping = BuildCappingPlan(protein, policy, verdict);

    int serial = 1;
    std::string current_chain;
    bool any_residue_emitted = false;

    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);

        // TER between different chain_ids.
        if (any_residue_emitted && res.chain_id != current_chain) {
            pdb_out << "TER\n";
        }
        current_chain = res.chain_id;

        // If this residue gets an ACE cap before it, emit the cap
        // (its own residue, with a synthesized sequence number one less
        // than the host residue's).
        if (capping.ace_before.count(ri) > 0) {
            const size_t n_idx = res.N;
            const Vec3 anchor =
                (n_idx != Residue::NONE) ? conf.PositionAt(n_idx)
                                          : Vec3(0.0, 0.0, 0.0);
            EmitCapResidue(pdb_out, serial, "ACE", kAceTemplate,
                           anchor, res.chain_id,
                           res.sequence_number - 1, res.insertion_code);
            map_out.extractor_index_for_prmtop_residue.push_back(
                ResidueAmberMapping::NONE_FOR_CAP);
        }

        const bool is_disulfide_cys = cyx_residues.count(ri) > 0;
        const std::string ambr_name =
            AmberResidueNameFor(res, is_disulfide_cys);

        for (size_t ai : res.atom_indices) {
            const Atom& atom = protein.AtomAt(ai);
            const Vec3 pos = conf.PositionAt(ai);
            WriteAtomRecord(pdb_out, serial,
                            atom.pdb_atom_name, ambr_name,
                            res.chain_id, res.sequence_number,
                            res.insertion_code, pos, atom.element);
            ++serial;
        }
        map_out.extractor_index_for_prmtop_residue.push_back(ri);
        any_residue_emitted = true;

        // If this residue gets an NME cap after it, emit the cap
        // following its atoms.
        if (capping.nme_after.count(ri) > 0) {
            const size_t c_idx = res.C;
            const Vec3 anchor =
                (c_idx != Residue::NONE) ? conf.PositionAt(c_idx)
                                          : Vec3(0.0, 0.0, 0.0);
            EmitCapResidue(pdb_out, serial, "NME", kNmeTemplate,
                           anchor, res.chain_id,
                           res.sequence_number + 1, res.insertion_code);
            map_out.extractor_index_for_prmtop_residue.push_back(
                ResidueAmberMapping::NONE_FOR_CAP);
        }
    }

    if (any_residue_emitted) {
        pdb_out << "TER\n";
    }
    pdb_out << "END\n";
}


std::vector<std::pair<size_t, size_t>> DetectDisulfides(
        const Protein& protein,
        const ProteinConformation& conf,
        double max_ss_distance_angstroms) {
    // Collect SG atom index for every CYS residue (typed: residue.type
    // and atom.pdb_atom_name == "SG").
    std::vector<std::pair<size_t, size_t>> cys_sg;  // (residue_index, atom_index)
    cys_sg.reserve(64);
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.type != AminoAcid::CYS) continue;
        for (size_t ai : res.atom_indices) {
            const Atom& atom = protein.AtomAt(ai);
            if (atom.pdb_atom_name == "SG" && atom.element == Element::S) {
                cys_sg.emplace_back(ri, ai);
                break;
            }
        }
    }

    const double max_sq = max_ss_distance_angstroms * max_ss_distance_angstroms;
    std::vector<std::pair<size_t, size_t>> pairs;
    for (size_t i = 0; i < cys_sg.size(); ++i) {
        for (size_t j = i + 1; j < cys_sg.size(); ++j) {
            const Vec3 pa = conf.PositionAt(cys_sg[i].second);
            const Vec3 pb = conf.PositionAt(cys_sg[j].second);
            const double sq = (pa - pb).squaredNorm();
            if (sq <= max_sq) {
                pairs.emplace_back(cys_sg[i].first, cys_sg[j].first);
            }
        }
    }
    return pairs;
}


void GenerateLeapScript(const LeapScriptInputs& inputs,
                        std::ostream& script_out) {
    script_out << "source leaprc.protein.ff14SB\n";
    script_out << "set default PBRadii mbondi2\n";
    script_out << "mol = loadPdb " << inputs.pdb_path << "\n";
    for (const auto& pair : inputs.disulfide_residue_pairs_1based) {
        script_out << "bond mol." << pair.first << ".SG mol."
                   << pair.second << ".SG\n";
    }
    script_out << "saveamberparm mol " << inputs.prmtop_path
               << " " << inputs.inpcrd_path << "\n";
    script_out << "quit\n";
}

}  // namespace amber_leap
}  // namespace nmr
