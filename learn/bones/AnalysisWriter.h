#pragma once
//
// AnalysisWriter: buffers per-frame ConformationAtom data and writes
// the analysis H5 file at the end of a trajectory analysis run.
//
// This is NOT a ConformationResult. It does not attach to conformations.
// It is a trajectory-level output writer that reads from conformations
// after all calculators have run.
//
// Usage:
//   AnalysisWriter writer(protein, protein_id);
//   writer.BeginTrajectory(n_frames_estimate, n_atoms);
//   for each frame:
//     OperationRunner::Run(conf, opts);
//     writer.HarvestFrame(conf, frame_idx, time_ps);
//   writer.WriteH5(path);
//
// All data is buffered in memory (float64). Write happens once at end.
// ~1.2 GB per protein for 600 frames × 335 atoms × 700 channels.
//
// See spec/ANALYSIS_TRAJECTORY_2026-04-14.md for the full schema.
//

#include "ProteinConformation.h"
#include "Protein.h"
#include "Bond.h"
#include "Ring.h"

#include <string>
#include <vector>
#include <cstddef>

namespace nmr {

class AnalysisWriter {
public:
    // Construct with protein reference (for topology metadata)
    // and human-readable ID.
    AnalysisWriter(const Protein& protein, const std::string& protein_id);

    // Allocate buffers for the expected number of frames and atoms.
    // stride is recorded as metadata (e.g. 2 = every other frame).
    void BeginTrajectory(size_t n_frames_estimate, size_t n_atoms,
                         size_t stride = 1);

    // Harvest one frame's data from a fully-computed conformation.
    // Called after OperationRunner::Run has attached all results.
    // frame_idx is the original XTC frame index (for the manifest).
    void HarvestFrame(const ProteinConformation& conf,
                      size_t frame_idx, double time_ps);

    // Write the complete analysis H5 to the given path.
    // Writes all buffered data in one shot. The file is either
    // complete or doesn't exist — no partial writes.
    void WriteH5(const std::string& path) const;

    // Write a PDB snapshot of the protein at the given conformation.
    // Uses PBC-fixed positions already on the conformation — no temp
    // files, no re-reading the trajectory. For external geometry
    // validation (MolProbity) independent of our tensor math.
    static void WritePdb(const Protein& protein,
                         const ProteinConformation& conf,
                         const std::string& path);

    // Number of frames harvested so far.
    size_t FrameCount() const { return n_frames_; }

private:
    const Protein& protein_;
    std::string protein_id_;
    size_t n_atoms_ = 0;
    size_t n_frames_ = 0;
    size_t stride_ = 1;

    // ── Frame metadata ──────────────────────────────────────────
    std::vector<double> frame_times_;       // (T,) ps
    std::vector<int32_t> frame_indices_;    // (T,) original XTC indices

    // ── Positions ───────────────────────────────────────────────
    // Flat buffer: T * N * 3, row-major (frame, atom, xyz)
    std::vector<double> positions_;

    // ── Ring current group ───────────────────────────────────────
    // Per-type arrays: T * N * 8 (T0) and T * N * 8 * 5 (T2)
    std::vector<double> bs_T0_per_type_;    // (T, N, 8)
    std::vector<double> bs_T2_per_type_;    // (T, N, 8, 5)
    std::vector<double> hm_T0_per_type_;    // (T, N, 8)
    std::vector<double> hm_T2_per_type_;    // (T, N, 8, 5)

    // SphericalTensor totals: T * N * 9  (T0 + T1[3] + T2[5])
    std::vector<double> bs_shielding_;      // (T, N, 9)
    std::vector<double> hm_shielding_;      // (T, N, 9)
    std::vector<double> rs_shielding_;      // (T, N, 9)  ringchi

    // Ring proximity scalars: T * N each
    std::vector<int16_t> n_rings_3A_;
    std::vector<int16_t> n_rings_5A_;
    std::vector<int16_t> n_rings_8A_;
    std::vector<int16_t> n_rings_12A_;
    std::vector<double> mean_ring_dist_;
    std::vector<double> nearest_ring_atom_;

    // Exponential-weighted sums
    std::vector<double> G_iso_exp_sum_;     // (T, N)
    std::vector<double> G_T2_exp_sum_;      // (T, N, 5)
    std::vector<double> G_iso_var_8A_;      // (T, N)

    // B-field vector
    std::vector<double> total_B_field_;     // (T, N, 3)

    // ── EFG group ───────────────────────────────────────────────
    std::vector<double> coulomb_efg_total_;       // (T, N, 9) SphericalTensor
    std::vector<double> coulomb_efg_backbone_;    // (T, N, 9)
    std::vector<double> coulomb_efg_aromatic_;    // (T, N, 9)
    std::vector<double> coulomb_E_total_;         // (T, N, 3)
    std::vector<double> coulomb_E_backbone_;      // (T, N, 3)
    std::vector<double> coulomb_E_sidechain_;     // (T, N, 3)
    std::vector<double> coulomb_E_aromatic_;      // (T, N, 3)
    std::vector<double> coulomb_E_solvent_;       // (T, N, 3)
    std::vector<double> coulomb_E_magnitude_;     // (T, N)
    std::vector<double> coulomb_E_bond_proj_;     // (T, N)
    std::vector<double> coulomb_E_bb_frac_;       // (T, N)
    std::vector<double> apbs_efg_;                // (T, N, 9)
    std::vector<double> apbs_efield_;             // (T, N, 3)
    std::vector<double> aimnet2_efg_total_;       // (T, N, 9)
    std::vector<double> aimnet2_efg_backbone_;    // (T, N, 9)
    std::vector<double> aimnet2_efg_aromatic_;    // (T, N, 9)
    std::vector<double> aimnet2_shielding_;       // (T, N, 9)

    // ── Bond anisotropy group ───────────────────────────────────
    std::vector<double> mc_shielding_;            // (T, N, 9)
    std::vector<double> mc_T2_backbone_;          // (T, N, 9)
    std::vector<double> mc_T2_sidechain_;         // (T, N, 9)
    std::vector<double> mc_T2_aromatic_;          // (T, N, 9)
    std::vector<double> mc_T2_CO_nearest_;        // (T, N, 9)
    std::vector<double> mc_T2_CN_nearest_;        // (T, N, 9)
    std::vector<double> mc_co_sum_;               // (T, N)
    std::vector<double> mc_cn_sum_;               // (T, N)
    std::vector<double> mc_sidechain_sum_;        // (T, N)
    std::vector<double> mc_aromatic_sum_;         // (T, N)
    std::vector<double> mc_co_nearest_;           // (T, N)
    std::vector<double> mc_nearest_CO_dist_;      // (T, N)
    std::vector<double> mc_nearest_CN_dist_;      // (T, N)
    std::vector<double> mc_nearest_CO_mid_;       // (T, N, 3)
    std::vector<double> mc_dir_nearest_CO_;       // (T, N, 3)

    // ── Quadrupole group ────────────────────────────────────────
    std::vector<double> pq_shielding_;            // (T, N, 9)
    std::vector<double> pq_T0_per_type_;          // (T, N, 8)
    std::vector<double> pq_T2_per_type_;          // (T, N, 8, 5)

    // ── Dispersion group ────────────────────────────────────────
    std::vector<double> disp_shielding_;          // (T, N, 9)
    std::vector<double> disp_T0_per_type_;        // (T, N, 8)
    std::vector<double> disp_T2_per_type_;        // (T, N, 8, 5)

    // ── H-bond group ────────────────────────────────────────────
    std::vector<double> hbond_shielding_;         // (T, N, 9)
    std::vector<double> hbond_nearest_sph_;       // (T, N, 9)
    std::vector<double> hbond_nearest_dist_;      // (T, N)
    std::vector<double> hbond_nearest_dir_;       // (T, N, 3)
    std::vector<double> hbond_inv_d3_;            // (T, N)
    std::vector<int16_t> hbond_count_;            // (T, N)
    std::vector<int8_t> hbond_is_donor_;          // (T, N)
    std::vector<int8_t> hbond_is_acceptor_;       // (T, N)
    std::vector<int8_t> hbond_is_backbone_;       // (T, N)

    // ── SASA group ──────────────────────────────────────────────
    std::vector<double> sasa_;                    // (T, N)
    std::vector<double> sasa_normal_;             // (T, N, 3)

    // ── Water group ─────────────────────────────────────────────
    std::vector<double> water_efield_;            // (T, N, 3)
    std::vector<double> water_efg_;               // (T, N, 9)
    std::vector<double> water_efield_first_;      // (T, N, 3)
    std::vector<double> water_efg_first_;         // (T, N, 9)
    std::vector<int16_t> water_n_first_;          // (T, N)
    std::vector<int16_t> water_n_second_;         // (T, N)
    std::vector<double> water_half_shell_;        // (T, N)
    std::vector<double> water_dipole_cos_;        // (T, N)
    std::vector<double> water_nearest_ion_dist_;  // (T, N)
    std::vector<double> water_nearest_ion_charge_;// (T, N)
    std::vector<double> water_dipole_vec_;        // (T, N, 3)
    std::vector<double> water_surface_normal_;    // (T, N, 3)
    std::vector<double> water_sasa_asym_;         // (T, N)
    std::vector<double> water_sasa_align_;        // (T, N)
    std::vector<double> water_sasa_cohere_;       // (T, N)
    std::vector<int16_t> water_sasa_first_n_;     // (T, N)

    // ── Charges group ───────────────────────────────────────────
    std::vector<double> aimnet2_charge_;           // (T, N)
    std::vector<double> eeq_charge_;              // (T, N)
    std::vector<double> eeq_cn_;                  // (T, N)

    // ── AIMNet2 embedding group ─────────────────────────────────
    std::vector<float> aimnet2_aim_;              // (T, N, 256) float32 native torch

    // ── Per-ring group (K=6 nearest rings) ──────────────────────
    static constexpr size_t K_RINGS = 6;
    // Geometry: (T, N, K, 6) — distance, rho, z, theta, cos_phi, sin_phi
    std::vector<double> per_ring_geom_;
    // Ring type: (T, N, K) — RingTypeIndex enum, -1 = no ring
    std::vector<int8_t> per_ring_type_;
    // Per-calculator T2: (T, N, K, 5) each
    std::vector<double> per_ring_bs_T2_;
    std::vector<double> per_ring_hm_T2_;
    std::vector<double> per_ring_chi_T2_;
    std::vector<double> per_ring_pq_T2_;
    std::vector<double> per_ring_hm_H_T2_;
    // Dispersion: scalar (T, N, K) + T2 (T, N, K, 5)
    std::vector<double> per_ring_disp_scalar_;
    std::vector<double> per_ring_disp_T2_;

    // ── Additional atom metadata (enrichment flags, graph, charges) ──
    // These are frame-invariant — written once from conformation 0.
    // Stored here to avoid re-reading conf0 at write time.
    // (Populated on first HarvestFrame call from conf0.)
    bool enrichment_captured_ = false;
    std::vector<int8_t> is_amide_H_;
    std::vector<int8_t> is_alpha_H_;
    std::vector<int8_t> is_methyl_;
    std::vector<int8_t> is_aromatic_H_;
    std::vector<int8_t> is_on_aromatic_residue_;
    std::vector<int8_t> is_hbond_donor_;
    std::vector<int8_t> is_hbond_acceptor_;
    std::vector<int8_t> parent_is_sp2_;
    std::vector<int32_t> graph_dist_N_;
    std::vector<int32_t> graph_dist_O_;
    std::vector<double> eneg_sum_1_;
    std::vector<double> eneg_sum_2_;
    std::vector<int32_t> n_pi_bonds_3_;
    std::vector<int32_t> bfs_to_nearest_ring_;
    std::vector<double> bfs_decay_;
    std::vector<double> partial_charge_;
    std::vector<double> vdw_radius_;

    // ── Per-frame ring geometry ─────────────────────────────────
    // (T, n_rings, 7) — center(3), normal(3), radius
    std::vector<double> ring_geometry_;

    // ── Coulomb shielding SphericalTensor ───────────────────────
    std::vector<double> coulomb_shielding_;        // (T, N, 9)

    // ── Bonded energy group (per-atom, from BondedEnergyResult) ─
    std::vector<double> bonded_bond_;          // (T, N)
    std::vector<double> bonded_angle_;         // (T, N)
    std::vector<double> bonded_ub_;            // (T, N)
    std::vector<double> bonded_proper_;        // (T, N)
    std::vector<double> bonded_improper_;      // (T, N)
    std::vector<double> bonded_cmap_;          // (T, N)
    std::vector<double> bonded_total_;         // (T, N)

    // ── Energy group (per-frame scalars + tensors from EDR) ────
    std::vector<double> energy_coulomb_sr_;     // (T,)
    std::vector<double> energy_coulomb_recip_;  // (T,)
    std::vector<double> energy_bond_;           // (T,)
    std::vector<double> energy_angle_;          // (T,)
    std::vector<double> energy_ub_;             // (T,)
    std::vector<double> energy_proper_dih_;     // (T,)
    std::vector<double> energy_improper_dih_;   // (T,)
    std::vector<double> energy_cmap_;           // (T,)
    std::vector<double> energy_lj_sr_;          // (T,)
    std::vector<double> energy_potential_;      // (T,)
    std::vector<double> energy_kinetic_;        // (T,)
    std::vector<double> energy_enthalpy_;       // (T,)
    std::vector<double> energy_temperature_;    // (T,)
    std::vector<double> energy_pressure_;       // (T,)
    std::vector<double> energy_volume_;         // (T,)
    std::vector<double> energy_density_;        // (T,)
    std::vector<double> energy_box_;            // (T, 3) — box_x, box_y, box_z
    std::vector<double> energy_virial_;         // (T, 9) — 3x3 tensor
    std::vector<double> energy_pres_tensor_;    // (T, 9) — 3x3 tensor
    std::vector<double> energy_T_protein_;      // (T,)
    std::vector<double> energy_T_non_protein_;  // (T,)

    // ── Dihedrals group (per-residue, T * R) ────────────────────
    std::vector<double> dih_phi_;          // (T, R) radians
    std::vector<double> dih_psi_;          // (T, R) radians
    std::vector<double> dih_omega_;        // (T, R) radians, CA-C-N-CA
    std::vector<double> dih_chi1_;         // (T, R) radians
    std::vector<double> dih_chi2_;         // (T, R) radians
    std::vector<double> dih_chi3_;         // (T, R) radians
    std::vector<double> dih_chi4_;         // (T, R) radians
    std::vector<double> dih_chi1_cos_;     // (T, R)
    std::vector<double> dih_chi1_sin_;     // (T, R)
    std::vector<double> dih_chi2_cos_;     // (T, R)
    std::vector<double> dih_chi2_sin_;     // (T, R)
    std::vector<double> dih_chi3_cos_;     // (T, R)
    std::vector<double> dih_chi3_sin_;     // (T, R)
    std::vector<double> dih_chi4_cos_;     // (T, R)
    std::vector<double> dih_chi4_sin_;     // (T, R)

    // ── DSSP group (per-residue, T * R) ─────────────────────────
    std::vector<int8_t> dssp_ss8_;         // (T, R) 8-class SS
    std::vector<double> dssp_hbond_energy_;// (T, R) strongest acceptor energy

    // ── Internal helpers (implementation in .cpp) ─────────────────

    // Append a SphericalTensor's 9 components to a flat buffer.
    static void AppendSpherical(std::vector<double>& buf,
                                const struct SphericalTensor& st);

    // Append a Vec3 to a flat buffer.
    static void AppendVec3(std::vector<double>& buf,
                           const Vec3& v);

    // Compute dihedral angle (radians) from four atom positions.
    // Returns NaN if degenerate (collinear atoms).
    static double Dihedral(const Vec3& p0, const Vec3& p1,
                           const Vec3& p2, const Vec3& p3);
};

}  // namespace nmr
