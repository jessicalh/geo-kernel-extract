#pragma once
//
// AnalysisFile: serialisable object model for the analysis H5 file format.
//
// Each nested struct corresponds to one HDF5 group.  All numeric data is
// stored in flat vectors whose logical shape is documented in comments
// using the notation (T, N, ...) where T = n_frames, N = n_atoms,
// R = n_residues.  Row-major (C order).
//
// ReadH5(path)  — populate from an existing file.
// WriteH5(path) — emit a new file from in-memory data.
//
// Round-trip guarantee: for files where aim is float32, every dataset
// is byte-identical after ReadH5 + WriteH5.  Files where aim was
// written as float64 (pre-fix 1CBH) will be normalised to float32 on
// read, so their aim group will differ.
//

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

struct AnalysisFile {

    // ── Dimensions (populated by ReadH5 from /meta attrs) ───────
    size_t n_frames    = 0;  // T
    size_t n_atoms     = 0;  // N
    size_t n_residues  = 0;  // R
    size_t n_bonds     = 0;  // B  (from /topology attrs)
    size_t n_rings     = 0;  //    (from /topology attrs)

    // ── /meta ───────────────────────────────────────────────────
    struct Meta {
        std::string protein_id;
        size_t stride = 0;
        std::vector<double>  frame_times;    // (T,)
        std::vector<int32_t> frame_indices;  // (T,)
    } meta;

    // ── /atoms — per-atom static properties ─────────────────────
    struct Atoms {
        std::vector<int32_t>     element;              // (N,) atomic number
        std::vector<int32_t>     residue_index;        // (N,)
        std::vector<std::string> atom_name;            // (N,) PDB name
        std::vector<int32_t>     atom_role;            // (N,) enum
        std::vector<int32_t>     hybridisation;        // (N,) enum
        std::vector<int32_t>     n_bonded;             // (N,)
        std::vector<int32_t>     graph_dist_ring;      // (N,)
        std::vector<int8_t>      is_backbone;          // (N,)
        std::vector<int8_t>      is_conjugated;        // (N,)
        // Enrichment flags (also static, written from conf 0)
        std::vector<int8_t>      is_amide_H;           // (N,)
        std::vector<int8_t>      is_alpha_H;           // (N,)
        std::vector<int8_t>      is_methyl;            // (N,)
        std::vector<int8_t>      is_aromatic_H;        // (N,)
        std::vector<int8_t>      is_on_aromatic_residue; // (N,)
        std::vector<int8_t>      is_hbond_donor;       // (N,)
        std::vector<int8_t>      is_hbond_acceptor;    // (N,)
        std::vector<int8_t>      parent_is_sp2;        // (N,)
        std::vector<int32_t>     graph_dist_N;         // (N,)
        std::vector<int32_t>     graph_dist_O;         // (N,)
        std::vector<double>      eneg_sum_1;           // (N,)
        std::vector<double>      eneg_sum_2;           // (N,)
        std::vector<int32_t>     n_pi_bonds_3;         // (N,)
        std::vector<int32_t>     bfs_to_nearest_ring;  // (N,)
        std::vector<double>      bfs_decay;            // (N,)
        std::vector<double>      partial_charge;       // (N,)
        std::vector<double>      vdw_radius;           // (N,)
    } atoms;

    // ── /residues — per-residue static properties ───────────────
    struct Residues {
        std::vector<std::string> residue_name;    // (R,)
        std::vector<int32_t>     residue_number;  // (R,)
        std::vector<std::string> chain_id;        // (R,)
    } residues;

    // ── /topology — bond graph + ring membership ────────────────
    struct Topology {
        std::vector<int32_t> bond_atoms;          // (B, 2) flat
        std::vector<int32_t> bond_category;       // (B,)
        std::vector<int32_t> bond_order;          // (B,)
        std::vector<int32_t> parent_atom_index;   // (N,) H→heavy, -1 if not H
        std::vector<int32_t> ring_type;           // (n_rings,)
        std::vector<int32_t> ring_residue;        // (n_rings,)
        std::vector<int32_t> ring_fused_partner;  // (n_rings,)
        std::vector<int32_t> ring_offsets;        // (n_rings+1,) CSR
        std::vector<int32_t> ring_atom_indices;   // (sum of ring sizes,) flat
    } topology;

    // ── /positions ──────────────────────────────────────────────
    struct Positions {
        std::vector<double> xyz;  // (T, N, 3) flat
    } positions;

    // ── /ring_current ──────────────────────────────────────────
    struct RingCurrent {
        std::vector<double>  bs_T0_per_type;   // (T, N, 8)
        std::vector<double>  bs_T2_per_type;   // (T, N, 8, 5)
        std::vector<double>  hm_T0_per_type;   // (T, N, 8)
        std::vector<double>  hm_T2_per_type;   // (T, N, 8, 5)
        std::vector<double>  bs_shielding;     // (T, N, 9) SphericalTensor
        std::vector<double>  hm_shielding;     // (T, N, 9) SphericalTensor
        std::vector<double>  rs_shielding;     // (T, N, 9) SphericalTensor
        std::vector<int16_t> n_rings_3A;       // (T, N)
        std::vector<int16_t> n_rings_5A;       // (T, N)
        std::vector<int16_t> n_rings_8A;       // (T, N)
        std::vector<int16_t> n_rings_12A;      // (T, N)
        std::vector<double>  mean_ring_dist;   // (T, N)
        std::vector<double>  nearest_ring_atom;// (T, N)
        std::vector<double>  G_iso_exp_sum;    // (T, N)
        std::vector<double>  G_T2_exp_sum;     // (T, N, 5)
        std::vector<double>  G_iso_var_8A;     // (T, N)
        std::vector<double>  total_B_field;    // (T, N, 3)
    } ring_current;

    // ── /efg ────────────────────────────────────────────────────
    struct EFG {
        std::vector<double> coulomb_total;      // (T, N, 9) SphericalTensor
        std::vector<double> coulomb_backbone;   // (T, N, 9) SphericalTensor
        std::vector<double> coulomb_aromatic;   // (T, N, 9) SphericalTensor
        std::vector<double> coulomb_shielding;  // (T, N, 9) SphericalTensor
        std::vector<double> E_total;            // (T, N, 3)
        std::vector<double> E_backbone;         // (T, N, 3)
        std::vector<double> E_sidechain;        // (T, N, 3)
        std::vector<double> E_aromatic;         // (T, N, 3)
        std::vector<double> E_solvent;          // (T, N, 3)
        std::vector<double> E_magnitude;        // (T, N)
        std::vector<double> E_bond_proj;        // (T, N)
        std::vector<double> E_backbone_frac;    // (T, N)
        std::vector<double> apbs_efg;           // (T, N, 9) SphericalTensor
        std::vector<double> apbs_efield;        // (T, N, 3)
        std::vector<double> aimnet2_total;      // (T, N, 9) SphericalTensor
        std::vector<double> aimnet2_backbone;   // (T, N, 9) SphericalTensor
        std::vector<double> aimnet2_aromatic;   // (T, N, 9) SphericalTensor
        std::vector<double> aimnet2_shielding;  // (T, N, 9) SphericalTensor
    } efg;

    // ── /bond_aniso ─────────────────────────────────────────────
    struct BondAniso {
        std::vector<double> mc_shielding;       // (T, N, 9) SphericalTensor
        std::vector<double> T2_backbone;        // (T, N, 9) SphericalTensor
        std::vector<double> T2_sidechain;       // (T, N, 9) SphericalTensor
        std::vector<double> T2_aromatic;        // (T, N, 9) SphericalTensor
        std::vector<double> T2_CO_nearest;      // (T, N, 9) SphericalTensor
        std::vector<double> T2_CN_nearest;      // (T, N, 9) SphericalTensor
        std::vector<double> co_sum;             // (T, N)
        std::vector<double> cn_sum;             // (T, N)
        std::vector<double> sidechain_sum;      // (T, N)
        std::vector<double> aromatic_sum;       // (T, N)
        std::vector<double> co_nearest;         // (T, N)
        std::vector<double> nearest_CO_dist;    // (T, N)
        std::vector<double> nearest_CN_dist;    // (T, N)
        std::vector<double> nearest_CO_midpoint;// (T, N, 3)
        std::vector<double> dir_nearest_CO;     // (T, N, 3)
    } bond_aniso;

    // ── /quadrupole ─────────────────────────────────────────────
    struct Quadrupole {
        std::vector<double> pq_shielding;       // (T, N, 9) SphericalTensor
        std::vector<double> pq_T0_per_type;     // (T, N, 8)
        std::vector<double> pq_T2_per_type;     // (T, N, 8, 5)
    } quadrupole;

    // ── /dispersion ─────────────────────────────────────────────
    struct Dispersion {
        std::vector<double> disp_shielding;     // (T, N, 9) SphericalTensor
        std::vector<double> disp_T0_per_type;   // (T, N, 8)
        std::vector<double> disp_T2_per_type;   // (T, N, 8, 5)
    } dispersion;

    // ── /hbond ──────────────────────────────────────────────────
    struct HBond {
        std::vector<double>  hbond_shielding;   // (T, N, 9) SphericalTensor
        std::vector<double>  nearest_spherical; // (T, N, 9) SphericalTensor
        std::vector<double>  nearest_dist;      // (T, N)
        std::vector<double>  nearest_dir;       // (T, N, 3)
        std::vector<double>  inv_d3;            // (T, N)
        std::vector<int16_t> count_3_5A;        // (T, N)
        std::vector<int8_t>  is_donor;          // (T, N)
        std::vector<int8_t>  is_acceptor;       // (T, N)
        std::vector<int8_t>  is_backbone;       // (T, N)
    } hbond;

    // ── /sasa ───────────────────────────────────────────────────
    struct SASA {
        std::vector<double> sasa;     // (T, N)
        std::vector<double> normal;   // (T, N, 3)
    } sasa;

    // ── /water ──────────────────────────────────────────────────
    struct Water {
        std::vector<double>  efield;            // (T, N, 3)
        std::vector<double>  efg;               // (T, N, 9) SphericalTensor
        std::vector<double>  efield_first;      // (T, N, 3)
        std::vector<double>  efg_first;         // (T, N, 9) SphericalTensor
        std::vector<int16_t> n_first;           // (T, N)
        std::vector<int16_t> n_second;          // (T, N)
        std::vector<double>  half_shell_asymmetry; // (T, N)
        std::vector<double>  dipole_cos;        // (T, N)
        std::vector<double>  nearest_ion_dist;  // (T, N)
        std::vector<double>  nearest_ion_charge;// (T, N)
        std::vector<double>  dipole_vector;     // (T, N, 3)
        std::vector<double>  surface_normal;    // (T, N, 3)
        std::vector<double>  sasa_asymmetry;    // (T, N)
        std::vector<double>  sasa_dipole_align; // (T, N)
        std::vector<double>  sasa_dipole_cohere;// (T, N)
        std::vector<int16_t> sasa_first_shell_n;// (T, N)
    } water;

    // ── /charges ────────────────────────────────────────────────
    struct Charges {
        std::vector<double> aimnet2_charge;  // (T, N)
        std::vector<double> eeq_charge;      // (T, N)
        std::vector<double> eeq_cn;          // (T, N)
    } charges;

    // ── /aimnet2_embedding ──────────────────────────────────────
    static constexpr size_t AIM_DIMS = 256;
    struct AimNet2Embedding {
        std::vector<float> aim;  // (T, N, 256) float32
        bool was_float64 = false; // true if source file had float64
    } aimnet2_embedding;

    // ── /per_ring (K=6 nearest rings per atom) ──────────────────
    static constexpr size_t K_RINGS = 6;
    struct PerRing {
        std::vector<double>      geometry;         // (T, N, K, 6)
        std::vector<std::string> geometry_fields;  // (6,)
        std::vector<int8_t>      ring_type;        // (T, N, K)
        std::vector<double>      bs_T2;            // (T, N, K, 5)
        std::vector<double>      hm_T2;            // (T, N, K, 5)
        std::vector<double>      chi_T2;           // (T, N, K, 5)
        std::vector<double>      pq_T2;            // (T, N, K, 5)
        std::vector<double>      hm_H_T2;          // (T, N, K, 5)
        std::vector<double>      disp_scalar;      // (T, N, K)
        std::vector<double>      disp_T2;          // (T, N, K, 5)
    } per_ring;

    // ── /ring_geometry ──────────────────────────────────────────
    struct RingGeometry {
        std::vector<double>      data;    // (T, n_rings, 7)
        std::vector<std::string> fields;  // (7,)
    } ring_geometry;

    // ── /bonded_energy ──────────────────────────────────────────
    struct BondedEnergy {
        std::vector<double> bond;         // (T, N)
        std::vector<double> angle;        // (T, N)
        std::vector<double> urey_bradley; // (T, N)
        std::vector<double> proper_dih;   // (T, N)
        std::vector<double> improper_dih; // (T, N)
        std::vector<double> cmap;         // (T, N)
        std::vector<double> total;        // (T, N)
    } bonded_energy;

    // ── /energy (per-frame system-level from GROMACS EDR) ───────
    struct Energy {
        std::vector<double> coulomb_sr;     // (T,)
        std::vector<double> coulomb_recip;  // (T,)
        std::vector<double> bond;           // (T,)
        std::vector<double> angle;          // (T,)
        std::vector<double> urey_bradley;   // (T,)
        std::vector<double> proper_dih;     // (T,)
        std::vector<double> improper_dih;   // (T,)
        std::vector<double> cmap_dih;       // (T,)
        std::vector<double> lj_sr;          // (T,)
        std::vector<double> potential;      // (T,)
        std::vector<double> kinetic;        // (T,)
        std::vector<double> enthalpy;       // (T,)
        std::vector<double> temperature;    // (T,)
        std::vector<double> pressure;       // (T,)
        std::vector<double> volume;         // (T,)
        std::vector<double> density;        // (T,)
        std::vector<double> box;            // (T, 3)
        std::vector<double> virial;         // (T, 9)
        std::vector<std::string> virial_layout; // (9,) "XX","XY",...
        std::vector<double> pressure_tensor;// (T, 9)
        std::vector<double> T_protein;      // (T,)
        std::vector<double> T_non_protein;  // (T,)
    } energy;

    // ── /dihedrals (per-residue, per-frame) ─────────────────────
    struct Dihedrals {
        std::vector<double> phi;       // (T, R)
        std::vector<double> psi;       // (T, R)
        std::vector<double> omega;     // (T, R)
        std::vector<double> chi1;      // (T, R)
        std::vector<double> chi2;      // (T, R)
        std::vector<double> chi3;      // (T, R)
        std::vector<double> chi4;      // (T, R)
        std::vector<double> chi1_cos;  // (T, R)
        std::vector<double> chi1_sin;  // (T, R)
        std::vector<double> chi2_cos;  // (T, R)
        std::vector<double> chi2_sin;  // (T, R)
        std::vector<double> chi3_cos;  // (T, R)
        std::vector<double> chi3_sin;  // (T, R)
        std::vector<double> chi4_cos;  // (T, R)
        std::vector<double> chi4_sin;  // (T, R)
    } dihedrals;

    // ── /dssp (per-residue, per-frame) ──────────────────────────
    struct DSSP {
        std::vector<int8_t> ss8;            // (T, R)
        std::vector<double> hbond_energy;   // (T, R)
    } dssp;

    // ── I/O ─────────────────────────────────────────────────────
    void ReadH5(const std::string& path);
    void WriteH5(const std::string& path) const;
};
