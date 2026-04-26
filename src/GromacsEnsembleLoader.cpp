#include "GromacsEnsembleLoader.h"
#include "AminoAcidType.h"
#include "NamingRegistry.h"
#include "ChargeSource.h"
#include "OperationLog.h"

// GROMACS internal API — linked against libgromacs.so
#include "gromacs/fileio/tpxio.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/inputrec.h"

#include <fstream>
#include <filesystem>
#include <cstdio>
#include <set>
#include <regex>
#include <cmath>

namespace fs = std::filesystem;

namespace nmr {


// ============================================================================
// Amino acid name set — for identifying the protein molecule block
// ============================================================================

static const std::set<std::string> KNOWN_AMINO_ACIDS = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "HID","HIE","HIP","HSD","HSE","HSP","ASH","GLH","ASPP","GLUP",
    "CYX","CYM","CYS2","LYN","ARN","TYM","MSE"
};


// ============================================================================
// Element from GROMACS atomic number
// ============================================================================

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
// ensemble.json parser
// ============================================================================

struct PoseInfo {
    int pose = 0;
    int walker = 0;
    double time_ps = 0.0;
    double rmsd_nm = 0.0;
    double rg_nm = 0.0;
    double weight = 0.0;
};

static bool ParseEnsembleJson(const std::string& path,
                               std::vector<PoseInfo>& poses,
                               std::string& error) {
    std::ifstream in(path);
    if (!in.is_open()) {
        error = "cannot open ensemble.json: " + path;
        return false;
    }

    std::string line;
    // Full format (fleet production): pose, walker, time_ps, rmsd_nm, rg_nm, weight
    std::regex re_full(
        R"(\{\"pose\":\s*(\d+).*\"walker\":\s*(\d+).*\"time_ps\":\s*([0-9.eE+-]+))"
        R"(.*\"rmsd_nm\":\s*([0-9.eE+-]+).*\"rg_nm\":\s*([0-9.eE+-]+))"
        R"(.*\"weight\":\s*([0-9.eE+-]+))");
    // Minimal format (fes-sampler): pose, walker, time_ps, weight (no rmsd/rg)
    std::regex re_minimal(
        R"(\{\"pose\":\s*(\d+).*\"walker\":\s*(\d+).*\"time_ps\":\s*([0-9.eE+-]+))"
        R"(.*\"weight\":\s*([0-9.eE+-]+))");
    std::smatch m;

    while (std::getline(in, line)) {
        if (std::regex_search(line, m, re_full)) {
            PoseInfo pi;
            pi.pose = std::stoi(m[1].str());
            pi.walker = std::stoi(m[2].str());
            pi.time_ps = std::stod(m[3].str());
            pi.rmsd_nm = std::stod(m[4].str());
            pi.rg_nm = std::stod(m[5].str());
            pi.weight = std::stod(m[6].str());
            poses.push_back(pi);
        } else if (std::regex_search(line, m, re_minimal)) {
            PoseInfo pi;
            pi.pose = std::stoi(m[1].str());
            pi.walker = std::stoi(m[2].str());
            pi.time_ps = std::stod(m[3].str());
            pi.weight = std::stod(m[4].str());
            poses.push_back(pi);
        }
    }

    if (poses.empty()) {
        error = "no poses found in " + path;
        return false;
    }
    return true;
}


// ============================================================================
// PDB position reader — ATOM records only, in file order
// ============================================================================

static std::vector<Vec3> ReadPdbPositions(const std::string& path,
                                           size_t expected_count,
                                           std::string& error) {
    std::vector<Vec3> positions;
    std::ifstream in(path);
    if (!in.is_open()) {
        error = "cannot open PDB: " + path;
        return positions;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.size() < 54) continue;
        if (line.substr(0, 4) != "ATOM" && line.substr(0, 6) != "HETATM") continue;

        double x = std::stod(line.substr(30, 8));
        double y = std::stod(line.substr(38, 8));
        double z = std::stod(line.substr(46, 8));
        positions.emplace_back(x, y, z);

        if (positions.size() == expected_count) break;
    }

    if (positions.size() != expected_count) {
        error = path + ": expected " + std::to_string(expected_count) +
                " atoms, got " + std::to_string(positions.size());
        positions.clear();
    }
    return positions;
}


// ============================================================================
// HIS tautomer detection from CHARMM atom names within a residue
//
// GROMACS writes "HIS" for the pdb2gmx default tautomer and preserves
// "HSP" for doubly protonated. We cannot distinguish HSD from HSE by
// residue name alone. The atom names in the topology are authoritative:
//   HD1 present, HE2 absent → HSD (delta, variant 0)
//   HD1 absent,  HE2 present → HSE (epsilon, variant 1)
//   HD1 present, HE2 present → HSP (doubly, variant 2)
// ============================================================================

static int DetectHisVariantFromAtoms(const t_atoms& atoms, int res_start, int res_end) {
    bool has_HD1 = false, has_HE2 = false;
    for (int ai = res_start; ai < res_end; ++ai) {
        std::string name = *(atoms.atomname[ai]);
        if (name == "HD1") has_HD1 = true;
        if (name == "HE2") has_HE2 = true;
    }
    if (has_HD1 && has_HE2) return 2;  // HIP/HSP
    if (has_HD1)            return 0;  // HID/HSD
    if (has_HE2)            return 1;  // HIE/HSE
    return -1;
}


// ============================================================================
// CHARMM residue name → variant_index
// ============================================================================

static int VariantFromCharmmResidueName(const std::string& charmm_name,
                                         AminoAcid canonical_aa) {
    const AminoAcidType& aat = GetAminoAcidType(canonical_aa);
    if (!aat.is_titratable) return -1;

    for (int vi = 0; vi < static_cast<int>(aat.variants.size()); ++vi) {
        const auto& variant = aat.variants[vi];
        if (charmm_name == "HSP" && canonical_aa == AminoAcid::HIS &&
            std::string(variant.registry_key) == "doubly") return vi;
        if (charmm_name == "HSD" && canonical_aa == AminoAcid::HIS &&
            std::string(variant.registry_key) == "delta") return vi;
        if (charmm_name == "HSE" && canonical_aa == AminoAcid::HIS &&
            std::string(variant.registry_key) == "epsilon") return vi;
        if (charmm_name == "ASPP" && canonical_aa == AminoAcid::ASP &&
            std::string(variant.registry_key) == "protonated") return vi;
        if (charmm_name == "GLUP" && canonical_aa == AminoAcid::GLU &&
            std::string(variant.registry_key) == "protonated") return vi;
        if (charmm_name == "CYS2" && canonical_aa == AminoAcid::CYS &&
            std::string(variant.registry_key) == "disulfide") return vi;
    }

    return -1;
}


// ============================================================================
// Resolve a pose PDB filename in the poses directory.
//
// Try pose_%03d.pdb first (legacy naming). If not found, scan for any PDB
// whose name contains _%02d_ matching the pose number (fes-sampler naming).
// Returns empty string on failure.
// ============================================================================

static std::string ResolvePosePdb(const std::string& dir, int pose_num) {
    // Try legacy naming first
    char legacy[32];
    std::snprintf(legacy, sizeof(legacy), "pose_%03d.pdb", pose_num);
    std::string legacy_path = dir + "/" + legacy;
    if (fs::exists(legacy_path)) return legacy_path;

    // Scan for fes-sampler naming: {protein_id}_{NN}_{reason}_{weight}.pdb
    // Match _%02d_ in the first half of the filename (before the reason
    // field, which can also contain numbers like max_dev_phi_7).
    char tag[8];
    std::snprintf(tag, sizeof(tag), "_%02d_", pose_num);
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.path().extension() != ".pdb") continue;
        std::string name = entry.path().filename().string();
        auto pos = name.find(tag);
        if (pos != std::string::npos && pos < name.size() / 2)
            return entry.path().string();
    }
    return "";
}


// ============================================================================
// BuildProteinFromTpr: TPR topology → Protein + charges (no poses)
// ============================================================================

BuildResult BuildProteinFromTpr(const std::string& tpr_path,
                                const std::string& protein_id,
                                ForceField force_field) {
    BuildResult result;

    OperationLog::Scope scope("BuildProteinFromTpr", tpr_path);

    if (!fs::exists(tpr_path)) {
        result.error = "TPR not found: " + tpr_path;
        return result;
    }

    TpxFileHeader tpx = readTpxHeader(tpr_path.c_str(), true);
    t_inputrec ir;
    t_state state;
    gmx_mtop_t mtop;
    read_tpx_state(tpr_path.c_str(), tpx.bIr ? &ir : nullptr,
                   &state, tpx.bTop ? &mtop : nullptr);

    // Find protein moltype
    const gmx_moltype_t* protein_moltype = nullptr;
    for (const auto& molblock : mtop.molblock) {
        const auto& moltype = mtop.moltype[molblock.type];
        if (moltype.atoms.nres > 0) {
            std::string first_res = *(moltype.atoms.resinfo[0].name);
            if (KNOWN_AMINO_ACIDS.count(first_res) > 0) {
                protein_moltype = &moltype;
                break;
            }
        }
    }

    if (!protein_moltype) {
        result.error = "no protein molecule found in TPR: " + tpr_path;
        return result;
    }

    const t_atoms& tpr_atoms = protein_moltype->atoms;
    size_t n_protein_atoms = tpr_atoms.nr;
    size_t n_residues = tpr_atoms.nres;

    // Build Protein from TPR topology
    auto& registry = GlobalNamingRegistry();
    auto protein = std::make_unique<Protein>();

    struct ResRange { int start; int end; };
    std::vector<ResRange> ranges(n_residues);
    {
        std::vector<int> first_atom(n_residues, -1);
        for (int ai = 0; ai < tpr_atoms.nr; ++ai) {
            int ri = tpr_atoms.atom[ai].resind;
            if (ri >= 0 && ri < static_cast<int>(n_residues) && first_atom[ri] < 0)
                first_atom[ri] = ai;
        }
        for (size_t ri = 0; ri < n_residues; ++ri) {
            ranges[ri].start = first_atom[ri];
            ranges[ri].end = (ri + 1 < n_residues) ? first_atom[ri + 1]
                                                     : tpr_atoms.nr;
        }
    }

    // Add residues
    for (size_t ri = 0; ri < n_residues; ++ri) {
        Residue res;
        std::string charmm_name = *(tpr_atoms.resinfo[ri].name);
        std::string canonical = registry.ToCanonical(charmm_name);
        if (canonical.empty()) {
            result.error = "unknown residue '" + charmm_name +
                           "' at position " + std::to_string(ri) +
                           " in " + protein_id;
            return result;
        }
        res.type = AminoAcidFromThreeLetterCode(canonical);
        if (res.type == AminoAcid::Unknown) {
            result.error = "unrecognised amino acid '" + canonical +
                           "' from '" + charmm_name + "' at " + std::to_string(ri);
            return result;
        }
        res.sequence_number = tpr_atoms.resinfo[ri].nr;
        res.chain_id = "A";
        res.protonation_variant_index =
            VariantFromCharmmResidueName(charmm_name, res.type);
        if (res.type == AminoAcid::HIS && res.protonation_variant_index < 0) {
            res.protonation_variant_index =
                DetectHisVariantFromAtoms(tpr_atoms, ranges[ri].start, ranges[ri].end);
        }
        protein->AddResidue(std::move(res));
    }

    // Add atoms
    for (int ai = 0; ai < tpr_atoms.nr; ++ai) {
        auto atom = std::make_unique<Atom>();
        std::string charmm_atom_name = *(tpr_atoms.atomname[ai]);
        size_t ri = tpr_atoms.atom[ai].resind;
        std::string canonical_residue =
            ThreeLetterCodeForAminoAcid(protein->ResidueAt(ri).type);
        atom->iupac_name = registry.TranslateAtomName(
            charmm_atom_name, canonical_residue,
            ToolContext::Charmm, ToolContext::Standard);
        atom->element = ElementFromAtomicNumber(tpr_atoms.atom[ai].atomnumber);
        atom->residue_index = ri;
        size_t idx = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(ri).atom_indices.push_back(idx);
        const IupacAtomName& name = protein->AtomAt(idx).iupac_name;
        Residue& res_ref = protein->MutableResidueAt(ri);
        if      (name == "N"  && res_ref.N  == Residue::NONE) res_ref.N  = idx;
        else if (name == "CA" && res_ref.CA == Residue::NONE) res_ref.CA = idx;
        else if (name == "C"  && res_ref.C  == Residue::NONE &&
                 protein->AtomAt(idx).element == Element::C &&
                 res_ref.CA != Residue::NONE) res_ref.C = idx;
        else if (name == "O"  && res_ref.O  == Residue::NONE) res_ref.O  = idx;
        else if ((name == "H" || name == "HN") &&
                 res_ref.H == Residue::NONE) res_ref.H = idx;
        else if ((name == "HA" || name == "HA2") &&
                 res_ref.HA == Residue::NONE) res_ref.HA = idx;
        else if (name == "CB" && res_ref.CB == Residue::NONE) res_ref.CB = idx;
    }

    // Build context
    auto ctx = std::make_unique<ProteinBuildContext>();
    ctx->pdb_source = protein_id;
    ctx->force_field = ForceFieldName(force_field);
    ctx->protonation_tool = "pdb2gmx";
    protein->SetBuildContext(std::move(ctx));

    // Extract charges from TPR
    std::vector<AtomChargeRadius> charge_vec(n_protein_atoms);
    double charge_sum = 0.0;
    for (size_t ai = 0; ai < n_protein_atoms; ++ai) {
        charge_vec[ai].charge = tpr_atoms.atom[ai].q;
        charge_vec[ai].radius = 1.5;
        charge_sum += tpr_atoms.atom[ai].q;
    }
    result.net_charge = static_cast<int>(
        charge_sum + (charge_sum > 0 ? 0.5 : -0.5));
    result.charges = std::make_unique<PreloadedChargeSource>(
        std::move(charge_vec), force_field);

    result.protein = std::move(protein);
    result.ok = true;

    OperationLog::Info(LogCalcOther, "BuildProteinFromTpr",
        protein_id + ": " +
        std::to_string(n_protein_atoms) + " atoms, " +
        std::to_string(n_residues) + " residues");

    return result;
}


// ============================================================================
// LoadFleetEnsemble
// ============================================================================

BuildResult BuildFromGromacs(const FleetPaths& paths) {
    BuildResult result;

    OperationLog::Scope scope("BuildFromGromacs", paths.sampled_poses_dir);

    std::string protein_id = fs::path(paths.sampled_poses_dir).filename().string();

    // ---------------------------------------------------------------
    // 1. Read TPR via libgromacs
    // ---------------------------------------------------------------
    if (!fs::exists(paths.tpr_path)) {
        result.error = "TPR not found: " + paths.tpr_path;
        return result;
    }

    TpxFileHeader tpx = readTpxHeader(paths.tpr_path.c_str(), true);
    t_inputrec ir;
    t_state state;
    gmx_mtop_t mtop;
    read_tpx_state(paths.tpr_path.c_str(), tpx.bIr ? &ir : nullptr,
                   &state, tpx.bTop ? &mtop : nullptr);

    // Find the protein molecule block (first block with amino acid residues)
    const gmx_moltype_t* protein_moltype = nullptr;
    for (const auto& molblock : mtop.molblock) {
        const auto& moltype = mtop.moltype[molblock.type];
        if (moltype.atoms.nres > 0) {
            std::string first_res = *(moltype.atoms.resinfo[0].name);
            if (KNOWN_AMINO_ACIDS.count(first_res) > 0) {
                protein_moltype = &moltype;
                break;
            }
        }
    }

    if (!protein_moltype) {
        result.error = "no protein molecule found in TPR: " + paths.tpr_path;
        return result;
    }

    const t_atoms& tpr_atoms = protein_moltype->atoms;
    size_t n_protein_atoms = tpr_atoms.nr;
    size_t n_residues = tpr_atoms.nres;

    // ---------------------------------------------------------------
    // 2. Parse ensemble.json
    // ---------------------------------------------------------------
    std::string ensemble_path = paths.sampled_poses_dir + "/ensemble.json";
    std::vector<PoseInfo> poses;
    if (!ParseEnsembleJson(ensemble_path, poses, result.error))
        return result;

    // ---------------------------------------------------------------
    // 3. Build Protein from TPR topology
    // ---------------------------------------------------------------
    auto& registry = GlobalNamingRegistry();
    auto protein = std::make_unique<Protein>();

    // Build residue atom ranges
    struct ResRange { int start; int end; };
    std::vector<ResRange> ranges(n_residues);
    {
        // Walk atoms to find first atom of each residue
        std::vector<int> first_atom(n_residues, -1);
        for (int ai = 0; ai < tpr_atoms.nr; ++ai) {
            int ri = tpr_atoms.atom[ai].resind;
            if (ri >= 0 && ri < static_cast<int>(n_residues) && first_atom[ri] < 0)
                first_atom[ri] = ai;
        }
        for (size_t ri = 0; ri < n_residues; ++ri) {
            ranges[ri].start = first_atom[ri];
            ranges[ri].end = (ri + 1 < n_residues) ? first_atom[ri + 1]
                                                     : tpr_atoms.nr;
        }
    }

    // Add residues
    for (size_t ri = 0; ri < n_residues; ++ri) {
        Residue res;

        std::string charmm_name = *(tpr_atoms.resinfo[ri].name);
        std::string canonical = registry.ToCanonical(charmm_name);
        if (canonical.empty()) {
            result.error = "unknown residue '" + charmm_name +
                           "' at position " + std::to_string(ri) +
                           " in " + protein_id;
            return result;
        }

        res.type = AminoAcidFromThreeLetterCode(canonical);
        if (res.type == AminoAcid::Unknown) {
            result.error = "unrecognised amino acid '" + canonical +
                           "' from '" + charmm_name + "' at " + std::to_string(ri);
            return result;
        }

        res.sequence_number = tpr_atoms.resinfo[ri].nr;
        res.chain_id = "A";

        // Variant from CHARMM name if explicit
        res.protonation_variant_index =
            VariantFromCharmmResidueName(charmm_name, res.type);

        // For "HIS" (ambiguous), inspect atom names
        if (res.type == AminoAcid::HIS && res.protonation_variant_index < 0) {
            res.protonation_variant_index =
                DetectHisVariantFromAtoms(tpr_atoms, ranges[ri].start, ranges[ri].end);
        }

        protein->AddResidue(std::move(res));
    }

    // Add atoms
    for (int ai = 0; ai < tpr_atoms.nr; ++ai) {
        auto atom = std::make_unique<Atom>();

        std::string charmm_atom_name = *(tpr_atoms.atomname[ai]);
        size_t ri = tpr_atoms.atom[ai].resind;
        std::string canonical_residue =
            ThreeLetterCodeForAminoAcid(protein->ResidueAt(ri).type);

        // PDB LOADING BOUNDARY: CHARMM atom name → canonical
        atom->iupac_name = registry.TranslateAtomName(
            charmm_atom_name, canonical_residue,
            ToolContext::Charmm, ToolContext::Standard);

        atom->element = ElementFromAtomicNumber(tpr_atoms.atom[ai].atomnumber);
        atom->residue_index = ri;

        size_t idx = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(ri).atom_indices.push_back(idx);

        // Cache backbone indices from canonical atom name
        const IupacAtomName& name = protein->AtomAt(idx).iupac_name;
        Residue& res_ref = protein->MutableResidueAt(ri);
        if      (name == "N"  && res_ref.N  == Residue::NONE) res_ref.N  = idx;
        else if (name == "CA" && res_ref.CA == Residue::NONE) res_ref.CA = idx;
        else if (name == "C"  && res_ref.C  == Residue::NONE &&
                 protein->AtomAt(idx).element == Element::C &&
                 res_ref.CA != Residue::NONE) res_ref.C = idx;
        else if (name == "O"  && res_ref.O  == Residue::NONE) res_ref.O  = idx;
        else if ((name == "H" || name == "HN") &&
                 res_ref.H == Residue::NONE) res_ref.H = idx;
        else if ((name == "HA" || name == "HA2") &&
                 res_ref.HA == Residue::NONE) res_ref.HA = idx;
        else if (name == "CB" && res_ref.CB == Residue::NONE) res_ref.CB = idx;
    }

    // ---------------------------------------------------------------
    // 4. Read first pose positions and finalize construction
    // ---------------------------------------------------------------
    std::string first_pose_path = ResolvePosePdb(paths.sampled_poses_dir, poses[0].pose);
    if (first_pose_path.empty()) {
        result.error = "pose PDB not found for pose " +
                       std::to_string(poses[0].pose) + " in " +
                       paths.sampled_poses_dir;
        return result;
    }

    auto first_positions = ReadPdbPositions(first_pose_path, n_protein_atoms,
                                            result.error);
    if (first_positions.empty()) return result;

    // Build context
    auto ctx = std::make_unique<ProteinBuildContext>();
    ctx->pdb_source = protein_id;
    ctx->force_field = ForceFieldName(paths.force_field);
    ctx->protonation_tool = "pdb2gmx";
    protein->SetBuildContext(std::move(ctx));

    protein->FinalizeConstruction(first_positions);

    // ---------------------------------------------------------------
    // 5. Create MDFrameConformations for all poses
    // ---------------------------------------------------------------
    // Store pose names from PDB filenames when they carry descriptive
    // identity (fes-sampler output).  Legacy pose_NNN.pdb names are
    // not stored — RunFleet falls back to frame_%03zu for those.
    auto store_pose_name = [&](const std::string& pdb_path) {
        std::string stem = fs::path(pdb_path).stem().string();
        if (stem.rfind("pose_", 0) != 0)  // not legacy naming
            result.pose_names.push_back(stem);
    };
    store_pose_name(first_pose_path);

    for (const auto& pose : poses) {
        std::vector<Vec3> positions;

        if (pose.pose == poses[0].pose) {
            positions = first_positions;
        } else {
            std::string pose_path = ResolvePosePdb(
                paths.sampled_poses_dir, pose.pose);
            if (pose_path.empty()) {
                result.error = "pose PDB not found for pose " +
                               std::to_string(pose.pose) + " in " +
                               paths.sampled_poses_dir;
                return result;
            }
            store_pose_name(pose_path);
            positions = ReadPdbPositions(pose_path, n_protein_atoms, result.error);
            if (positions.empty()) return result;
        }

        protein->AddMDFrame(std::move(positions),
                            pose.walker, pose.time_ps, pose.weight,
                            pose.rmsd_nm, pose.rg_nm);
    }

    // ---------------------------------------------------------------
    // 6. Extract charges from TPR → PreloadedChargeSource
    // ---------------------------------------------------------------
    std::vector<AtomChargeRadius> charge_vec(n_protein_atoms);
    double charge_sum = 0.0;
    for (size_t ai = 0; ai < n_protein_atoms; ++ai) {
        charge_vec[ai].charge = tpr_atoms.atom[ai].q;
        charge_vec[ai].radius = 1.5;  // default PB radius
        charge_sum += tpr_atoms.atom[ai].q;
    }
    result.net_charge = static_cast<int>(
        charge_sum + (charge_sum > 0 ? 0.5 : -0.5));
    result.charges = std::make_unique<PreloadedChargeSource>(
        std::move(charge_vec), paths.force_field);

    result.protein = std::move(protein);
    result.ok = true;

    OperationLog::Info(LogCharges, "BuildFromGromacs",
        protein_id + ": " +
        std::to_string(n_protein_atoms) + " atoms, " +
        std::to_string(n_residues) + " residues, " +
        std::to_string(poses.size()) + " poses");

    return result;
}


// ============================================================================
// RunAllFrames: delegates to OperationRunner::RunEnsemble.
// ============================================================================

std::vector<RunResult> RunAllFrames(
        Protein& protein, const RunOptions& opts) {
    return OperationRunner::RunEnsemble(protein, opts);
}

}  // namespace nmr
