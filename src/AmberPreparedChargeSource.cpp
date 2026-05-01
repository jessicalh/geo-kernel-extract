#include "AmberPreparedChargeSource.h"
#include "AminoAcidType.h"
#include "Atom.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinBuildContext.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RuntimeEnvironment.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace nmr {

AmberPreparedChargeSource::AmberPreparedChargeSource(
        const Protein& protein,
        AmberPreparationPolicy policy,
        AmberFlatTableCoverageVerdict reason,
        AmberSourceConfig config)
    : protein_(protein)
    , policy_(policy)
    , reason_(std::move(reason))
    , config_(std::move(config)) {}

std::string AmberPreparedChargeSource::Describe() const {
    std::ostringstream oss;
    oss << "AmberPreparedPrmtop:"
        << AmberPreparationPolicyName(policy_)
        << ":reason=" << reason_.Detail();

    // Methods-text-facing summary of capping decisions, when applicable.
    if (policy_ == AmberPreparationPolicy::
            UseCappedFragmentsForUnsupportedTerminalVariants) {
        bool first = true;
        for (const auto& f : reason_.failures) {
            if (f.kind != AmberFlatTableCoverageKind::
                    UnsupportedTerminalVariant) continue;
            // Only include cases that map to ACE/NME-capping.
            if (f.ff_residue_name != "ASH" && f.ff_residue_name != "CYM" &&
                    f.ff_residue_name != "GLH" && f.ff_residue_name != "LYN") {
                continue;
            }
            if (first) {
                oss << " caps=";
                first = false;
            } else {
                oss << ",";
            }
            const bool n_end = (f.terminal_token == "NTERM" ||
                                f.terminal_token == "NCTERM");
            const bool c_end = (f.terminal_token == "CTERM" ||
                                f.terminal_token == "NCTERM");
            if (n_end) oss << "ACE-";
            oss << f.ff_residue_name << "@" << f.chain_id
                << f.residue_sequence_number;
            if (c_end) oss << "-NME";
        }
    }
    return oss.str();
}

std::string AmberPreparedChargeSource::GeneratedPdb(
        const ProteinConformation& conf) const {
    std::ostringstream pdb_out;
    amber_leap::GenerateAmberPdb(protein_, conf, policy_, reason_,
                                  pdb_out, residue_mapping_);
    disulfide_extractor_pairs_ = amber_leap::DetectDisulfides(protein_, conf);
    mapping_built_ = true;
    return pdb_out.str();
}

std::string AmberPreparedChargeSource::GeneratedLeapScript(
        const std::string& pdb_path,
        const std::string& prmtop_path,
        const std::string& inpcrd_path) const {
    if (!mapping_built_) {
        std::fprintf(stderr,
            "FATAL: AmberPreparedChargeSource::GeneratedLeapScript called "
            "before GeneratedPdb populated the residue mapping.\n");
        std::abort();
    }
    amber_leap::LeapScriptInputs inputs;
    inputs.pdb_path = pdb_path;
    inputs.prmtop_path = prmtop_path;
    inputs.inpcrd_path = inpcrd_path;
    inputs.disulfide_residue_pairs_1based = DisulfidePairs1Based();

    std::ostringstream script_out;
    amber_leap::GenerateLeapScript(inputs, script_out);
    return script_out.str();
}

std::vector<std::pair<size_t, size_t>>
AmberPreparedChargeSource::DisulfidePairs1Based() const {
    // Project extractor residue indices through residue_mapping_ to get
    // 1-based PRMTOP indices. Skips any disulfide whose residues map
    // to NONE_FOR_CAP (cap policy doesn't apply to internal CYS pairs;
    // this is a safety condition only).
    std::unordered_map<size_t, size_t> ext_to_prmtop_1based;
    for (size_t prmtop_ri = 0;
         prmtop_ri < residue_mapping_.extractor_index_for_prmtop_residue.size();
         ++prmtop_ri) {
        size_t ext_ri =
            residue_mapping_.extractor_index_for_prmtop_residue[prmtop_ri];
        if (ext_ri == amber_leap::ResidueAmberMapping::NONE_FOR_CAP) continue;
        ext_to_prmtop_1based[ext_ri] = prmtop_ri + 1;
    }

    std::vector<std::pair<size_t, size_t>> pairs;
    for (const auto& ext_pair : disulfide_extractor_pairs_) {
        auto ita = ext_to_prmtop_1based.find(ext_pair.first);
        auto itb = ext_to_prmtop_1based.find(ext_pair.second);
        if (ita == ext_to_prmtop_1based.end() ||
                itb == ext_to_prmtop_1based.end()) continue;
        pairs.emplace_back(ita->second, itb->second);
    }
    return pairs;
}

// ============================================================================
// Local PRMTOP section readers. Duplicated from OrcaRunLoader.cpp /
// PrmtopChargeSource.cpp; unifying the three copies into a shared
// PrmtopParser is a queued post-slice cleanup.
// ============================================================================

namespace {

constexpr double kAmberChargeFactor = 18.2223;

std::vector<std::string> ReadPrmtopFlagLines(const std::string& path,
                                              const std::string& flag_name) {
    std::vector<std::string> lines;
    std::ifstream in(path);
    if (!in.is_open()) return lines;
    std::string line;
    bool in_flag = false;
    bool saw_format = false;
    while (std::getline(in, line)) {
        if (line.rfind("%FLAG ", 0) == 0) {
            if (in_flag) break;
            in_flag = (line.find("%FLAG " + flag_name) != std::string::npos);
            saw_format = false;
            continue;
        }
        if (!in_flag) continue;
        if (line.rfind("%FORMAT", 0) == 0) {
            saw_format = true;
            continue;
        }
        if (saw_format) lines.push_back(line);
    }
    return lines;
}

std::vector<double> ReadPrmtopDoubles(const std::string& path,
                                       const std::string& flag_name) {
    auto lines = ReadPrmtopFlagLines(path, flag_name);
    std::vector<double> values;
    for (const auto& line : lines) {
        std::istringstream iss(line);
        double v;
        while (iss >> v) values.push_back(v);
    }
    return values;
}

std::vector<int> ReadPrmtopInts(const std::string& path,
                                 const std::string& flag_name) {
    auto lines = ReadPrmtopFlagLines(path, flag_name);
    std::vector<int> values;
    for (const auto& line : lines) {
        std::istringstream iss(line);
        int v;
        while (iss >> v) values.push_back(v);
    }
    return values;
}

std::vector<std::string> ReadPrmtopStrings(const std::string& path,
                                            const std::string& flag_name,
                                            int field_width) {
    auto lines = ReadPrmtopFlagLines(path, flag_name);
    std::vector<std::string> values;
    for (const auto& line : lines) {
        for (size_t pos = 0;
             pos + static_cast<size_t>(field_width) <= line.size();
             pos += field_width) {
            std::string token = line.substr(pos, field_width);
            // Trim trailing spaces.
            while (!token.empty() && token.back() == ' ') token.pop_back();
            // Trim leading spaces.
            size_t lead = 0;
            while (lead < token.size() && token[lead] == ' ') ++lead;
            token = token.substr(lead);
            if (!token.empty()) values.push_back(token);
        }
    }
    return values;
}

struct PrmtopFields {
    std::vector<std::string> atom_names;
    std::vector<std::string> residue_labels;
    std::vector<int>         residue_pointers;   // 1-based per AMBER convention
    std::vector<double>      charges_amber;      // raw AMBER units; divide by 18.2223
    std::vector<double>      radii;
    std::string error;
    bool Ok() const { return error.empty(); }
};

PrmtopFields ReadPrmtopFields(const std::string& prmtop_path) {
    PrmtopFields p;
    p.atom_names       = ReadPrmtopStrings(prmtop_path, "ATOM_NAME", 4);
    p.residue_labels   = ReadPrmtopStrings(prmtop_path, "RESIDUE_LABEL", 4);
    p.residue_pointers = ReadPrmtopInts(prmtop_path, "RESIDUE_POINTER");
    p.charges_amber    = ReadPrmtopDoubles(prmtop_path, "CHARGE");
    p.radii            = ReadPrmtopDoubles(prmtop_path, "RADII");

    if (p.atom_names.empty())       p.error = "missing ATOM_NAME";
    else if (p.residue_labels.empty())  p.error = "missing RESIDUE_LABEL";
    else if (p.residue_pointers.empty()) p.error = "missing RESIDUE_POINTER";
    else if (p.charges_amber.empty()) p.error = "missing CHARGE";
    else if (p.radii.empty())          p.error = "missing RADII";
    else if (p.charges_amber.size() != p.atom_names.size())
        p.error = "CHARGE size mismatch ATOM_NAME";
    else if (p.radii.size() != p.atom_names.size())
        p.error = "RADII size mismatch ATOM_NAME";
    return p;
}

// PRMTOP atom index → 0-based PRMTOP residue index.
std::vector<size_t> BuildPrmtopResidueIndexForAtom(
        const PrmtopFields& p) {
    std::vector<size_t> result(p.atom_names.size(), 0);
    for (size_t r = 0; r < p.residue_pointers.size(); ++r) {
        size_t start = static_cast<size_t>(p.residue_pointers[r]) - 1;  // 1→0
        size_t stop  = (r + 1 < p.residue_pointers.size())
            ? static_cast<size_t>(p.residue_pointers[r + 1]) - 1
            : p.atom_names.size();
        for (size_t a = start; a < stop && a < p.atom_names.size(); ++a) {
            result[a] = r;
        }
    }
    return result;
}

}  // namespace


std::vector<AtomChargeRadius> AmberPreparedChargeSource::LoadCharges(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) const {
    OperationLog::Scope scope("AmberPreparedChargeSource::LoadCharges",
        std::string("policy=") + AmberPreparationPolicyName(policy_) +
        " atoms=" + std::to_string(conf.AtomCount()));

    // Hard preconditions: re-tleap protection. Violations indicate the
    // resolver dispatch upstream is broken; abort rather than degrade.
    if (protein.HasForceFieldCharges()) {
        std::fprintf(stderr,
            "FATAL: AmberPreparedChargeSource on Protein that already has "
            "ForceFieldCharges (re-tleap protection).\n");
        std::abort();
    }
    if (!protein.BuildContext().prmtop_path.empty()) {
        std::fprintf(stderr,
            "FATAL: AmberPreparedChargeSource when build_context.prmtop_path "
            "is set to %s (re-tleap protection).\n",
            protein.BuildContext().prmtop_path.c_str());
        std::abort();
    }

    // Resolve tleap binary: config.tleap_path overrides
    // RuntimeEnvironment::Tleap(); both paths must exist.
    std::string tleap_bin = config_.tleap_path;
    if (tleap_bin.empty()) tleap_bin = RuntimeEnvironment::Tleap();
    if (tleap_bin.empty() || !fs::exists(tleap_bin)) {
        error_out = "AmberPreparedChargeSource: no usable tleap binary "
                    "(config.tleap_path=\"" + config_.tleap_path + "\", "
                    "RuntimeEnvironment::Tleap()=\"" +
                    RuntimeEnvironment::Tleap() + "\"). Set AMBERHOME or "
                    "tleap in ~/.nmr_tools.toml.";
        return {};
    }

    // Resolve work directory: config.work_dir overrides RuntimeEnvironment.
    std::string work_dir = config_.work_dir;
    if (work_dir.empty()) {
        work_dir = RuntimeEnvironment::TempFilePath(
            "amber_prepared", "workdir");
    }
    std::error_code ec;
    fs::create_directories(work_dir, ec);
    if (ec) {
        error_out = "AmberPreparedChargeSource: cannot create work_dir " +
                    work_dir + ": " + ec.message();
        return {};
    }

    const std::string pdb_path     = work_dir + "/input.pdb";
    const std::string script_path  = work_dir + "/leap.in";
    const std::string prmtop_path  = work_dir + "/prep.prmtop";
    const std::string inpcrd_path  = work_dir + "/prep.inpcrd";
    const std::string log_path     = work_dir + "/tleap.log";

    // 1. Generate PDB + populate residue_mapping_ + disulfide_extractor_pairs_.
    {
        std::ofstream pdb_out(pdb_path);
        if (!pdb_out) {
            error_out = "AmberPreparedChargeSource: cannot write " + pdb_path;
            return {};
        }
        pdb_out << GeneratedPdb(conf);
    }

    // 2. Generate LEaP script.
    {
        std::ofstream script_out(script_path);
        if (!script_out) {
            error_out = "AmberPreparedChargeSource: cannot write " + script_path;
            return {};
        }
        script_out << GeneratedLeapScript(pdb_path, prmtop_path, inpcrd_path);
    }

    // 3. Run tleap. Stdout/stderr to log file.
    const std::string cmd = tleap_bin + " -f " + script_path +
                            " > " + log_path + " 2>&1";
    int rc = std::system(cmd.c_str());
    if (rc != 0 || !fs::exists(prmtop_path)) {
        error_out = "AmberPreparedChargeSource: tleap failed (rc=" +
                    std::to_string(rc) + "); see " + log_path;
        OperationLog::Error("AmberPreparedChargeSource::LoadCharges", error_out);
        return {};
    }

    // 4. Parse PRMTOP fields.
    auto prmtop = ReadPrmtopFields(prmtop_path);
    if (!prmtop.Ok()) {
        error_out = "AmberPreparedChargeSource: PRMTOP parse failed: " +
                    prmtop.error + " (file=" + prmtop_path + ")";
        return {};
    }

    // 5. Atom mapping: walk PRMTOP atoms, map to extractor atoms by
    //    (extractor residue index from residue_mapping_, atom name).
    if (residue_mapping_.extractor_index_for_prmtop_residue.size() !=
            prmtop.residue_labels.size()) {
        error_out = "AmberPreparedChargeSource: residue_mapping has " +
                    std::to_string(residue_mapping_
                        .extractor_index_for_prmtop_residue.size()) +
                    " entries but PRMTOP has " +
                    std::to_string(prmtop.residue_labels.size()) +
                    " residues";
        return {};
    }

    auto prmtop_res_for_atom = BuildPrmtopResidueIndexForAtom(prmtop);

    // Build extractor (residue_index, atom_name) → atom_index map.
    std::vector<std::map<std::string, size_t>>
        name_index_by_extractor_residue(protein.ResidueCount());
    for (size_t ai = 0; ai < protein.AtomCount(); ++ai) {
        const Atom& atom = protein.AtomAt(ai);
        auto& m = name_index_by_extractor_residue[atom.residue_index];
        if (m.count(atom.pdb_atom_name) > 0) {
            error_out = "AmberPreparedChargeSource: duplicate atom name " +
                        atom.pdb_atom_name + " in residue index " +
                        std::to_string(atom.residue_index);
            return {};
        }
        m[atom.pdb_atom_name] = ai;
    }

    const size_t n_extractor_atoms = protein.AtomCount();
    std::vector<AtomChargeRadius> result(n_extractor_atoms);
    std::vector<bool> seen(n_extractor_atoms, false);

    for (size_t prmtop_ai = 0; prmtop_ai < prmtop.atom_names.size(); ++prmtop_ai) {
        size_t prmtop_ri = prmtop_res_for_atom[prmtop_ai];
        size_t ext_ri = residue_mapping_
            .extractor_index_for_prmtop_residue[prmtop_ri];
        if (ext_ri == amber_leap::ResidueAmberMapping::NONE_FOR_CAP) {
            continue;  // ACE/NME/NHE inserted-cap atom — no extractor counterpart.
        }
        const std::string& aname = prmtop.atom_names[prmtop_ai];
        auto& m = name_index_by_extractor_residue[ext_ri];
        auto it = m.find(aname);
        if (it == m.end()) {
            error_out = "AmberPreparedChargeSource: PRMTOP atom \"" + aname +
                        "\" in residue \"" + prmtop.residue_labels[prmtop_ri] +
                        "\" (extractor residue index " +
                        std::to_string(ext_ri) +
                        ") has no counterpart in extractor protein";
            return {};
        }
        size_t ext_ai = it->second;
        if (seen[ext_ai]) {
            error_out = "AmberPreparedChargeSource: PRMTOP atom \"" + aname +
                        "\" mapped twice to extractor atom index " +
                        std::to_string(ext_ai);
            return {};
        }
        seen[ext_ai] = true;
        result[ext_ai].partial_charge =
            prmtop.charges_amber[prmtop_ai] / kAmberChargeFactor;
        result[ext_ai].pb_radius = prmtop.radii[prmtop_ai];
        result[ext_ai].status = ChargeAssignmentStatus::Matched;
    }

    // Every extractor atom must have been seen.
    for (size_t ai = 0; ai < n_extractor_atoms; ++ai) {
        if (seen[ai]) continue;
        const Atom& atom = protein.AtomAt(ai);
        const Residue& res = protein.ResidueAt(atom.residue_index);
        error_out = "AmberPreparedChargeSource: extractor atom \"" +
                    atom.pdb_atom_name + "\" in residue " +
                    std::to_string(res.sequence_number) +
                    " (extractor residue index " +
                    std::to_string(atom.residue_index) +
                    ") has no PRMTOP counterpart";
        return {};
    }

    OperationLog::Info(LogCharges, "AmberPreparedChargeSource::LoadCharges",
        "tleap_bin=" + tleap_bin + " work_dir=" + work_dir +
        " atoms=" + std::to_string(n_extractor_atoms) +
        " policy=" + AmberPreparationPolicyName(policy_));

    return result;
}

}  // namespace nmr
