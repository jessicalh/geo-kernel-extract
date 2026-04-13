#include "GromacsProtein.h"
#include "DsspResult.h"
#include "OperationLog.h"
#include "Types.h"

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

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
    for (size_t i = 0; i < protein_->AtomCount(); ++i) {
        atoms_.emplace_back(GromacsProteinAtom(i));
        atoms_.back().n_bonded_neighbors =
            static_cast<int>(protein_->AtomAt(i).bond_indices.size());
    }

    // Initialize per-bond accumulators from covalent topology
    bonds_accum_.clear();
    bonds_accum_.reserve(protein_->BondCount());
    for (size_t bi = 0; bi < protein_->BondCount(); ++bi) {
        GromacsProteinBond gb;
        gb.atom_a = protein_->BondAt(bi).atom_index_a;
        gb.atom_b = protein_->BondAt(bi).atom_index_b;
        bonds_accum_.push_back(gb);
    }
}


// ── AccumulateFrame ──────────────────────────────────────────────

void GromacsProtein::AccumulateFrame(
        const ProteinConformation& conf, int frame_idx,
        double time_ps) {
    const size_t N = atoms_.size();
    const Protein& prot = conf.ProteinRef();

    // Store raw positions for movies and waveform training
    const auto& positions = conf.Positions();
    size_t base = frame_positions_.size();
    frame_positions_.resize(base + N * 3);
    for (size_t i = 0; i < N; ++i) {
        frame_positions_[base + i * 3 + 0] = positions[i].x();
        frame_positions_[base + i * 3 + 1] = positions[i].y();
        frame_positions_[base + i * 3 + 2] = positions[i].z();
    }
    frame_times_.push_back(time_ps);
    ++n_stored_frames_;

    // DSSP: per-residue data accessible if DSSP ran this frame
    const DsspResult* dssp = conf.HasResult<DsspResult>()
        ? &conf.Result<DsspResult>() : nullptr;

    for (size_t ai = 0; ai < N; ++ai) {
        const auto& ca = conf.AtomAt(ai);
        auto& ga = atoms_[ai];

        // ── Positions ──────────────────────────────────────────
        Vec3 pos = ca.Position();
        ga.position_x.Update(pos.x(), frame_idx);
        ga.position_y.Update(pos.y(), frame_idx);
        ga.position_z.Update(pos.z(), frame_idx);

        // ── Ring current (BiotSavart, HaighMallion, RingSusceptibility) ──
        double bs_t0 = ca.bs_shielding_contribution.T0;
        double bs_t2 = ca.bs_shielding_contribution.T2Magnitude();
        ga.bs_T0.Update(bs_t0, frame_idx);
        ga.bs_T2mag.Update(bs_t2, frame_idx);
        ga.bs_T0_delta.UpdateDelta(bs_t0, frame_idx);

        ga.hm_T0.Update(ca.hm_shielding_contribution.T0, frame_idx);
        ga.hm_T2mag.Update(ca.hm_shielding_contribution.T2Magnitude(), frame_idx);

        ga.rs_T0.Update(ca.ringchi_shielding_contribution.T0, frame_idx);
        ga.rs_T2mag.Update(ca.ringchi_shielding_contribution.T2Magnitude(), frame_idx);

        // ── Bond anisotropy (McConnell) ────────────────────────
        ga.mc_T0.Update(ca.mc_shielding_contribution.T0, frame_idx);
        ga.mc_T2mag.Update(ca.mc_shielding_contribution.T2Magnitude(), frame_idx);
        ga.mc_aromatic.Update(ca.mcconnell_aromatic_sum, frame_idx);
        ga.mc_backbone.Update(ca.mcconnell_co_sum + ca.mcconnell_cn_sum, frame_idx);

        // ── Pi-quadrupole ──────────────────────────────────────
        ga.pq_T0.Update(ca.piquad_shielding_contribution.T0, frame_idx);
        ga.pq_T2mag.Update(ca.piquad_shielding_contribution.T2Magnitude(), frame_idx);

        // ── Dispersion ─────────────────────────────────────────
        ga.disp_T0.Update(ca.disp_shielding_contribution.T0, frame_idx);
        ga.disp_T2mag.Update(ca.disp_shielding_contribution.T2Magnitude(), frame_idx);

        // ── H-bond ─────────────────────────────────────────────
        ga.hbond_inv_d3.Update(ca.hbond_inv_d3, frame_idx);
        ga.hbond_count.Update(
            static_cast<double>(ca.hbond_count_within_3_5A), frame_idx);

        // ── APBS solvated field ────────────────────────────────
        if (ca.apbs_efield.norm() > 0.0 || ca.apbs_efg_spherical.T0 != 0.0) {
            ga.apbs_emag.Update(ca.apbs_efield.norm(), frame_idx);
            ga.apbs_efg_T0.Update(ca.apbs_efg_spherical.T0, frame_idx);
            ga.apbs_efg_T2mag.Update(ca.apbs_efg_spherical.T2Magnitude(), frame_idx);
        }

        // ── AIMNet2 ────────────────────────────────────────────
        // aimnet2_charge is 0.0 if AIMNet2 didn't run — accumulate
        // conditionally to avoid poisoning the Welford with zeros.
        if (ca.aimnet2_charge != 0.0 || ca.aimnet2_EFG_total_spherical.T0 != 0.0) {
            ga.aimnet2_charge.Update(ca.aimnet2_charge, frame_idx);
            ga.aimnet2_charge_delta.UpdateDelta(ca.aimnet2_charge, frame_idx);
            ga.aimnet2_efg_T0.Update(
                ca.aimnet2_EFG_total_spherical.T0, frame_idx);
            ga.aimnet2_efg_T2mag.Update(
                ca.aimnet2_EFG_total_spherical.T2Magnitude(), frame_idx);
        }

        // ── SASA ───────────────────────────────────────────────
        ga.sasa.Update(ca.atom_sasa, frame_idx);
        ga.sasa_normal_x.Update(ca.sasa_normal.x(), frame_idx);
        ga.sasa_normal_y.Update(ca.sasa_normal.y(), frame_idx);
        ga.sasa_normal_z.Update(ca.sasa_normal.z(), frame_idx);
        ga.sasa_delta.UpdateDelta(ca.atom_sasa, frame_idx);

        // ── Water environment ──────────────────────────────────
        double nf = static_cast<double>(ca.water_n_first);
        ga.water_n_first.Update(nf, frame_idx);
        ga.water_n_first_delta.UpdateDelta(nf, frame_idx);
        ga.water_n_second.Update(
            static_cast<double>(ca.water_n_second), frame_idx);
        ga.water_emag.Update(ca.water_efield.norm(), frame_idx);
        ga.water_emag_first.Update(ca.water_efield_first.norm(), frame_idx);
        ga.water_efield_x.Update(ca.water_efield.x(), frame_idx);
        ga.water_efield_y.Update(ca.water_efield.y(), frame_idx);
        ga.water_efield_z.Update(ca.water_efield.z(), frame_idx);

        // ── Hydration geometry ─────────────────────────────────
        ga.half_shell.Update(ca.half_shell_asymmetry, frame_idx);
        ga.dipole_cos.Update(ca.mean_water_dipole_cos, frame_idx);
        if (ca.nearest_ion_distance < 1e30)  // infinity when no ion within cutoff
            ga.nearest_ion_dist.Update(ca.nearest_ion_distance, frame_idx);

        // ── DSSP dynamics ──────────────────────────────────────
        if (dssp) {
            size_t ri = prot.AtomAt(ai).residue_index;
            if (ri < dssp->AllResidues().size()) {
                const auto& dr = dssp->AllResidues()[ri];

                // Phi/psi (backbone dihedrals)
                ga.phi_cos.Update(std::cos(dr.phi), frame_idx);
                ga.psi_cos.Update(std::cos(dr.psi), frame_idx);

                // Secondary structure transitions
                double ss_val = static_cast<double>(dr.secondary_structure);
                ga.ss8_transitions.Update(ss_val, 0.5);

                // H-bond energy: best (most negative) of the 4 partners
                double best_e = 0.0;
                for (int h = 0; h < 2; ++h) {
                    if (dr.acceptors[h].energy < best_e)
                        best_e = dr.acceptors[h].energy;
                    if (dr.donors[h].energy < best_e)
                        best_e = dr.donors[h].energy;
                }
                ga.dssp_hbond_energy.Update(best_e, frame_idx);

                // Chi1-4: atom indices from Residue::chi[k], same as DsspResult
                const auto& res = prot.ResidueAt(ri);
                for (int k = 0; k < 4; ++k) {
                    if (!res.chi[k].Valid()) break;
                    Vec3 p0 = conf.PositionAt(res.chi[k].a[0]);
                    Vec3 p1 = conf.PositionAt(res.chi[k].a[1]);
                    Vec3 p2 = conf.PositionAt(res.chi[k].a[2]);
                    Vec3 p3 = conf.PositionAt(res.chi[k].a[3]);
                    Vec3 b1 = p1 - p0;
                    Vec3 b2 = p2 - p1;
                    Vec3 b3 = p3 - p2;
                    Vec3 n1 = b1.cross(b2);
                    Vec3 n2 = b2.cross(b3);
                    double n1n = n1.norm();
                    double n2n = n2.norm();
                    if (n1n > 1e-10 && n2n > 1e-10) {
                        n1 /= n1n;
                        n2 /= n2n;
                        double cos_chi = std::max(-1.0, std::min(1.0, n1.dot(n2)));
                        ga.chi_cos[k].Update(cos_chi, frame_idx);
                        ga.chi_transitions[k].Update(cos_chi, 0.5);
                    }
                }
            }
        }

        // ── Bond angles at this atom (for GNN) ─────────────────
        // Mean cos(angle) across all bond angles A-this-B
        const auto& bond_ids = prot.AtomAt(ai).bond_indices;
        if (bond_ids.size() >= 2) {
            double cos_sum = 0.0;
            int n_angles = 0;
            for (size_t bi = 0; bi < bond_ids.size(); ++bi) {
                const auto& b1 = prot.BondAt(bond_ids[bi]);
                size_t other1 = (b1.atom_index_a == ai) ? b1.atom_index_b : b1.atom_index_a;
                Vec3 d1 = (conf.PositionAt(other1) - pos).normalized();
                for (size_t bj = bi + 1; bj < bond_ids.size(); ++bj) {
                    const auto& b2 = prot.BondAt(bond_ids[bj]);
                    size_t other2 = (b2.atom_index_a == ai) ? b2.atom_index_b : b2.atom_index_a;
                    Vec3 d2 = (conf.PositionAt(other2) - pos).normalized();
                    cos_sum += d1.dot(d2);
                    ++n_angles;
                }
            }
            if (n_angles > 0)
                ga.mean_bond_angle_cos.Update(cos_sum / n_angles, frame_idx);
        }

        // ── Frame counting ─────────────────────────────────────
        if (ca.water_n_first == 0) ga.n_frames_dry++;
        if (ca.water_n_first >= 4) ga.n_frames_exposed++;
    }

    // ── Per-bond accumulation (bond lengths) ───────────────────
    for (size_t bi = 0; bi < bonds_accum_.size(); ++bi) {
        auto& gb = bonds_accum_[bi];
        double len = (positions[gb.atom_b] - positions[gb.atom_a]).norm();
        gb.length.Update(len, frame_idx);
        gb.length_delta.UpdateDelta(len, frame_idx);
    }
}


// ── WriteCatalog ─────────────────────────────────────────────────

void GromacsProtein::WriteCatalog(const std::string& path) const {
    std::ofstream csv(path);

    // Column list driven by AllWelfords() — single source of truth.
    // Add a Welford to AllWelfords() and it appears here automatically.
    auto columns = atoms_[0].AllWelfords();

    // Header: fixed identity columns + data-driven Welford columns
    csv << "atom_idx,element,resname,resnum,atom_name,n_bonded,rmsf,";
    for (const auto& col : columns)
        csv << col.name << "_mean," << col.name << "_std,";
    // Transition counters (not Welfords, appended separately)
    csv << "chi1_transitions,chi2_transitions,chi3_transitions,chi4_transitions,"
        << "ss8_transitions,"
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
            << ga.n_bonded_neighbors << ","
            << ga.RMSF() << ",";

        // Data-driven Welford columns
        auto atom_cols = ga.AllWelfords();
        for (const auto& col : atom_cols)
            csv << col.w->mean << "," << col.w->Std() << ",";

        // Transition counters
        for (int k = 0; k < 4; ++k)
            csv << ga.chi_transitions[k].transitions << ",";
        csv << ga.ss8_transitions.transitions << ",";

        // Frame counting
        csv << ga.n_frames_dry << "," << ga.n_frames_exposed << ","
            << ga.water_n_first.count << "\n";
    }
}


// ── WriteH5 ─────────────────────────────────────────────────────
//
// Master HDF5 file per protein.  Contains:
//   /positions          (T, N, 3)   float64  — raw xyz per frame
//   /frame_times        (T,)        float64  — time in ps
//   /rollup/mean        (N, K)      float64  — Welford means
//   /rollup/std         (N, K)      float64  — Welford stds
//   /rollup/names       (K,)        string   — column names
//   /bonds/length_mean  (B,)        float64
//   /bonds/length_std   (B,)        float64
//   /bonds/atom_a       (B,)        uint64
//   /bonds/atom_b       (B,)        uint64
//   /metadata           attributes  — protein_id, n_atoms, n_frames
//
// This is the training data for the waveform model, the analysis
// data for the ridge, and the movie data for VTK.

void GromacsProtein::WriteH5(const std::string& path) const {
    using namespace HighFive;
    File file(path, File::Truncate);

    const size_t N = atoms_.size();
    const size_t T = n_stored_frames_;
    const size_t B = bonds_accum_.size();

    // Metadata
    file.createAttribute("protein_id", protein_id_);
    file.createAttribute("n_atoms", N);
    file.createAttribute("n_frames", T);
    file.createAttribute("n_bonds", B);

    // /positions (T*N, 3) stored flat, reshaped by consumer to (T, N, 3)
    if (T > 0 && N > 0) {
        std::vector<std::vector<double>> pos_2d(T * N, std::vector<double>(3));
        for (size_t t = 0; t < T; ++t) {
            size_t base = t * N * 3;
            for (size_t i = 0; i < N; ++i) {
                pos_2d[t * N + i][0] = frame_positions_[base + i * 3 + 0];
                pos_2d[t * N + i][1] = frame_positions_[base + i * 3 + 1];
                pos_2d[t * N + i][2] = frame_positions_[base + i * 3 + 2];
            }
        }
        file.createDataSet("/positions", pos_2d);
        file.createAttribute("positions_shape_T", T);
        file.createAttribute("positions_shape_N", N);

        // /frame_times (T,)
        file.createDataSet("/frame_times", frame_times_);
    }

    // Rollup: data-driven from AllWelfords() — same source as WriteCatalog.
    // No switch statement, no manual indexing. Adding a Welford to
    // AllWelfords() makes it appear in the H5 automatically.
    auto columns = atoms_[0].AllWelfords();
    const size_t K = columns.size();

    std::vector<std::string> names;
    names.reserve(K);
    for (const auto& col : columns)
        names.push_back(col.name);

    std::vector<std::vector<double>> means(N, std::vector<double>(K));
    std::vector<std::vector<double>> stds(N, std::vector<double>(K));
    for (size_t i = 0; i < N; ++i) {
        auto atom_cols = atoms_[i].AllWelfords();
        for (size_t k = 0; k < K; ++k) {
            means[i][k] = atom_cols[k].w->mean;
            stds[i][k] = atom_cols[k].w->Std();
        }
    }

    auto rollup = file.createGroup("/rollup");
    rollup.createDataSet("mean", means);
    rollup.createDataSet("std", stds);
    rollup.createDataSet("names", names);

    // Bonds
    if (B > 0) {
        std::vector<double> bond_len_mean(B), bond_len_std(B);
        std::vector<size_t> bond_a(B), bond_b(B);
        for (size_t bi = 0; bi < B; ++bi) {
            bond_len_mean[bi] = bonds_accum_[bi].length.mean;
            bond_len_std[bi] = bonds_accum_[bi].length.Std();
            bond_a[bi] = bonds_accum_[bi].atom_a;
            bond_b[bi] = bonds_accum_[bi].atom_b;
        }
        auto bonds_grp = file.createGroup("/bonds");
        bonds_grp.createDataSet("length_mean", bond_len_mean);
        bonds_grp.createDataSet("length_std", bond_len_std);
        bonds_grp.createDataSet("atom_a", bond_a);
        bonds_grp.createDataSet("atom_b", bond_b);
    }

    OperationLog::Info(LogCalcOther, "GromacsProtein::WriteH5",
        path + ": " + std::to_string(T) + " frames, " +
        std::to_string(N) + " atoms, " + std::to_string(K) +
        " rollup columns, " + std::to_string(B) + " bonds");
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
