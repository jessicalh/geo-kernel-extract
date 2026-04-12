#include "GromacsProtein.h"
#include "OperationLog.h"
#include "Types.h"

#include <filesystem>
#include <fstream>
#include <set>

namespace fs = std::filesystem;
namespace nmr {


GromacsProtein::GromacsProtein() = default;
GromacsProtein::~GromacsProtein() = default;


// ── Build (fleet path — pre-extracted PDB poses) ─────────────────

bool GromacsProtein::Build(const FleetPaths& paths) {
    auto result = BuildFromGromacs(paths);
    if (!result.ok) {
        error_ = result.error;
        return false;
    }

    protein_ = std::move(result.protein);
    charges_ = std::move(result.charges);
    net_charge_ = result.net_charge;
    pose_names_ = std::move(result.pose_names);

    protein_id_ = fs::path(paths.sampled_poses_dir).filename().string();

    InitAtoms();

    OperationLog::Info(LogCalcOther, "GromacsProtein::Build",
        protein_id_ + ": " +
        std::to_string(protein_->AtomCount()) + " atoms, " +
        std::to_string(protein_->ResidueCount()) + " residues, " +
        std::to_string(protein_->ConformationCount()) + " poses loaded");

    return true;
}


// ── Build (trajectory path — TPR is authoritative) ───────────────

bool GromacsProtein::BuildFromTrajectory(const std::string& tpr_path) {

    // Full-system topology for frame splitting (protein/water/ion ranges)
    if (!sys_reader_.ReadTopology(tpr_path)) {
        error_ = "TPR topology: " + sys_reader_.error();
        return false;
    }

    // Protein construction from TPR (same builder as fleet path)
    protein_id_ = fs::path(tpr_path).stem().string();
    auto build = BuildProteinFromTpr(tpr_path, protein_id_);
    if (!build.Ok()) {
        error_ = "TPR protein: " + build.error;
        return false;
    }

    protein_ = std::move(build.protein);
    charges_ = std::move(build.charges);
    net_charge_ = build.net_charge;

    // NOTE: Protein is not finalized yet — FinalizeProtein() needs
    // positions from the first XTC frame for bond detection. Called
    // by GromacsFrameHandler after reading the first frame.

    OperationLog::Info(LogCalcOther, "GromacsProtein::BuildFromTrajectory",
        protein_id_ + ": " +
        std::to_string(protein_->AtomCount()) + " protein, " +
        std::to_string(sys_reader_.Topology().water_count) + " water, " +
        std::to_string(sys_reader_.Topology().ion_count) + " ions");

    return true;
}


// ── FinalizeProtein ──────────────────────────────────────────────

void GromacsProtein::FinalizeProtein(
        std::vector<Vec3> first_frame_positions, double time_ps) {
    // 1. Bond detection from first frame geometry
    protein_->FinalizeConstruction(first_frame_positions);

    // 2. Conformation 0 — permanent, lives in the Protein
    protein_->AddMDFrame(std::move(first_frame_positions),
                         /*walker=*/0, time_ps, /*weight=*/1.0,
                         /*rmsd_nm=*/0.0, /*rg_nm=*/0.0);

    // 3. Accumulators
    InitAtoms();

    OperationLog::Info(LogCalcOther, "GromacsProtein::FinalizeProtein",
        protein_id_ + ": " +
        std::to_string(protein_->AtomCount()) + " atoms, " +
        std::to_string(protein_->BondCount()) + " bonds, " +
        std::to_string(protein_->RingCount()) + " rings, " +
        "conf0 at t=" + std::to_string(time_ps) + "ps");
}


// ── InitAtoms ────────────────────────────────────────────────────

void GromacsProtein::InitAtoms() {
    atoms_.clear();
    atoms_.reserve(protein_->AtomCount());
    for (size_t i = 0; i < protein_->AtomCount(); ++i)
        atoms_.emplace_back(GromacsProteinAtom(i));
}


// ── AccumulateFrame ──────────────────────────────────────────────

void GromacsProtein::AccumulateFrame(
        const ProteinConformation& conf, int frame_idx) {
    const size_t N = atoms_.size();
    for (size_t ai = 0; ai < N; ++ai) {
        const auto& ca = conf.AtomAt(ai);
        auto& ga = atoms_[ai];

        ga.sasa.Update(ca.atom_sasa, frame_idx);
        ga.water_n_first.Update(
            static_cast<double>(ca.water_n_first), frame_idx);
        ga.water_n_second.Update(
            static_cast<double>(ca.water_n_second), frame_idx);
        ga.water_emag.Update(ca.water_efield.norm(), frame_idx);
        ga.water_emag_first.Update(ca.water_efield_first.norm(), frame_idx);
        ga.half_shell.Update(ca.half_shell_asymmetry, frame_idx);
        ga.dipole_cos.Update(ca.mean_water_dipole_cos, frame_idx);
        ga.nearest_ion_dist.Update(ca.nearest_ion_distance, frame_idx);

        Vec3 pos = ca.Position();
        ga.position_x.Update(pos.x(), frame_idx);
        ga.position_y.Update(pos.y(), frame_idx);
        ga.position_z.Update(pos.z(), frame_idx);

        if (ca.water_n_first == 0) ga.n_frames_dry++;
        if (ca.water_n_first >= 4) ga.n_frames_exposed++;
    }
}


// ── WriteCatalog ─────────────────────────────────────────────────

void GromacsProtein::WriteCatalog(const std::string& path) const {
    std::ofstream csv(path);
    csv << "atom_idx,element,resname,resnum,atom_name,"
        << "rmsf,"
        << "water_n_first_mean,water_n_first_std,"
        << "water_n_first_min,water_n_first_min_frame,"
        << "water_n_first_max,water_n_first_max_frame,"
        << "water_emag_mean,water_emag_std,"
        << "water_emag_min,water_emag_min_frame,"
        << "water_emag_max,water_emag_max_frame,"
        << "sasa_mean,sasa_std,"
        << "sasa_min,sasa_min_frame,"
        << "sasa_max,sasa_max_frame,"
        << "half_shell_mean,half_shell_std,"
        << "dipole_cos_mean,dipole_cos_std,"
        << "nearest_ion_mean,"
        << "n_frames_dry,n_frames_exposed,n_frames_total\n";

    for (size_t ai = 0; ai < atoms_.size(); ++ai) {
        const auto& ga = atoms_[ai];
        const auto& atom = protein_->AtomAt(ai);
        const auto& residue = protein_->ResidueAt(atom.residue_index);

        const char* elem = "?";
        switch (atom.element) {
            case Element::H: elem = "H"; break;
            case Element::C: elem = "C"; break;
            case Element::N: elem = "N"; break;
            case Element::O: elem = "O"; break;
            case Element::S: elem = "S"; break;
            default: break;
        }

        csv << ai << "," << elem << ","
            << GetAminoAcidType(residue.type).three_letter_code << ","
            << residue.sequence_number << ","
            << atom.pdb_atom_name << ","
            << ga.RMSF() << ","
            << ga.water_n_first.mean << "," << ga.water_n_first.Std() << ","
            << ga.water_n_first.min_val << "," << ga.water_n_first.min_frame << ","
            << ga.water_n_first.max_val << "," << ga.water_n_first.max_frame << ","
            << ga.water_emag.mean << "," << ga.water_emag.Std() << ","
            << ga.water_emag.min_val << "," << ga.water_emag.min_frame << ","
            << ga.water_emag.max_val << "," << ga.water_emag.max_frame << ","
            << ga.sasa.mean << "," << ga.sasa.Std() << ","
            << ga.sasa.min_val << "," << ga.sasa.min_frame << ","
            << ga.sasa.max_val << "," << ga.sasa.max_frame << ","
            << ga.half_shell.mean << "," << ga.half_shell.Std() << ","
            << ga.dipole_cos.mean << "," << ga.dipole_cos.Std() << ","
            << ga.nearest_ion_dist.mean << ","
            << ga.n_frames_dry << "," << ga.n_frames_exposed << ","
            << ga.water_n_first.count << "\n";
    }
}


// ── SelectFrames ─────────────────────────────────────────────────

std::vector<size_t> GromacsProtein::SelectFrames(size_t max_frames) const {
    std::set<size_t> selected;

    for (size_t ai = 0; ai < atoms_.size(); ++ai) {
        const auto& ga = atoms_[ai];
        if (ga.water_n_first.Std() > 0.5) {
            if (ga.water_n_first.min_frame >= 0)
                selected.insert(static_cast<size_t>(ga.water_n_first.min_frame));
            if (ga.water_n_first.max_frame >= 0)
                selected.insert(static_cast<size_t>(ga.water_n_first.max_frame));
        }
        if (ga.water_emag.Std() > 0.1) {
            if (ga.water_emag.min_frame >= 0)
                selected.insert(static_cast<size_t>(ga.water_emag.min_frame));
            if (ga.water_emag.max_frame >= 0)
                selected.insert(static_cast<size_t>(ga.water_emag.max_frame));
        }
    }

    std::vector<size_t> result(selected.begin(), selected.end());
    std::sort(result.begin(), result.end());
    if (result.size() > max_frames)
        result.resize(max_frames);
    return result;
}


}  // namespace nmr
