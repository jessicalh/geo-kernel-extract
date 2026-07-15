#pragma once
//
// ConformationAtom: per-atom computed data within a ProteinConformation.
//
// All computed fields are public and mutable. One ConformationResult
// instance per type is enforced; by the convention that each result
// writes only its own fields, each field then has a single writer.
//
// Identity (element, bonds, residue) is not stored here — look it up on
// the parallel Protein atom (same index).
//

#include "Types.h"
#include <vector>
#include <array>
#include <cstdint>

namespace nmr {

// AIMNet2 aim embedding dimensionality.
// Defined here (not AIMNet2Result.h) because ConformationAtom needs it
// for the array extent, and ConformationAtom.h must not include torch headers.
static constexpr size_t AIMNET2_AIM_DIMS = 256;
static constexpr size_t AIMNET2_AIM_PROJECTION_DIMS = 32;

class ProteinConformation;

// Per-atom, per-ring structured result shared by ring calculators.
struct RingNeighbourhood {
    size_t ring_index = 0;
    RingTypeIndex ring_type = RingTypeIndex::PheBenzene;
    double distance_to_center = 0.0;
    Vec3 direction_to_center = Vec3::Zero();
    double rho = 0.0;
    double z = 0.0;
    double theta = 0.0;
    Mat3 G_tensor = Mat3::Zero();
    SphericalTensor G_spherical;
    Vec3 B_field = Vec3::Zero();
    Vec3 B_cylindrical = Vec3::Zero();
    Mat3 hm_H_tensor = Mat3::Zero();           // raw surface integral H (symmetric, traceless, pure T2)
    SphericalTensor hm_H_spherical;             // Decompose(H) — T0=0, T1=0, all info in T2
    Vec3 hm_B_field = Vec3::Zero();             // effective B-field V = H . n
    Mat3 hm_G_tensor = Mat3::Zero();            // full shielding kernel G = -n ⊗ V (rank-1)
    SphericalTensor hm_G_spherical;             // Decompose(G) — T0, T1, T2 all non-zero
    double quad_scalar = 0.0;
    // PiQuadrupoleLocalTensorResult-owned companion payload.  The raw
    // tensor and its decomposition are expressed in the deterministic
    // ring-local frame whose axes are stored alongside them.  `evaluated`
    // distinguishes a geometrically invalid row (evaluated=true,
    // valid=false, NaN payload) from a row appended after that result ran.
    Mat3 piquad_local_tensor = Mat3::Zero();
    SphericalTensor piquad_local_spherical;
    Vec3 piquad_local_x_axis = Vec3::Zero();
    Vec3 piquad_local_y_axis = Vec3::Zero();
    Vec3 piquad_local_z_axis = Vec3::Zero();
    double piquad_local_distance = 0.0;
    double piquad_local_cos_theta = 0.0;
    bool piquad_local_evaluated = false;
    bool piquad_local_valid = false;
    double chi_scalar = 0.0;
    double disp_scalar = 0.0;
    int disp_contacts = 0;
    double cos_phi = 1.0;   // azimuthal: cos of in-plane angle to vertex 0
    double sin_phi = 0.0;   // azimuthal: sin of in-plane angle to vertex 0
};

// Per-atom, per-bond structured result written by McConnellResult.
struct BondNeighbourhood {
    size_t bond_index = 0;
    BondCategory bond_category = BondCategory::Unknown;
    double distance_to_midpoint = 0.0;
    Vec3 direction_to_midpoint = Vec3::Zero();
    Mat3 dipolar_tensor = Mat3::Zero();
    SphericalTensor dipolar_spherical;
    double mcconnell_scalar = 0.0;
};

// Per-atom spatial neighbour
struct AtomNeighbour {
    size_t atom_index = 0;
    double distance = 0.0;
    Vec3 direction = Vec3::Zero();
};

// Per-atom MOPAC quantum bond order to another atom
struct MopacBondNeighbour {
    size_t other_atom = 0;              // atom index in protein
    double wiberg_order = 0.0;          // continuous bond order (0.01–3.0)
    size_t topology_bond_index = SIZE_MAX; // into protein.Bonds(), SIZE_MAX if no covalent bond
};


class ConformationAtom {
    friend class ProteinConformation;
public:
    Vec3 Position() const { return position_; }

    // === Enrichment properties (set by EnrichmentResult) ===
    AtomRole role = AtomRole::Unknown;
    Hybridisation hybridisation = Hybridisation::Unassigned;
    bool is_backbone = false;
    bool is_amide_H = false;
    bool is_alpha_H = false;
    bool is_methyl = false;
    bool is_aromatic_H = false;
    bool is_on_aromatic_residue = false;
    bool is_hbond_donor = false;
    bool is_hbond_acceptor = false;
    bool parent_is_sp2 = false;

    // === Charges and PB radii (ChargeAssignmentResult) ===
    double partial_charge = 0.0;
    double pb_radius = 0.0;

    // === MOPAC semiempirical results (MopacResult) ===
    double mopac_charge = 0.0;            // Coulson charge, legacy F15.6 projection (e)
    double mopac_s_pop = 0.0;             // s-orbital population
    double mopac_p_pop = 0.0;             // p-orbital population
    double mopac_valency = 0.0;           // sum of legacy compact-row Wiberg orders
    std::vector<MopacBondNeighbour> mopac_bond_neighbours;  // legacy compact rows, F6.3 order

    // === Spatial neighbourhood (SpatialIndexResult) ===
    std::vector<AtomNeighbour> spatial_neighbours;

    // === Ring neighbourhood (BiotSavartResult et al.) ===
    std::vector<RingNeighbourhood> ring_neighbours;

    // === Bond neighbourhood (McConnellResult) ===
    std::vector<BondNeighbourhood> bond_neighbours;

    // === Ring current totals (BiotSavartResult, HaighMallionResult) ===
    Vec3 total_B_field = Vec3::Zero();
    Mat3 total_G_tensor = Mat3::Zero();
    SphericalTensor total_G_spherical;
    std::array<double, 8> per_type_G_T0_sum = {};
    std::array<std::array<double, 3>, 8> per_type_G_T1_sum = {};
    std::array<std::array<double, 5>, 8> per_type_G_T2_sum = {};
    std::array<double, 8> per_type_hm_T0_sum = {};
    std::array<std::array<double, 3>, 8> per_type_hm_T1_sum = {};
    std::array<std::array<double, 5>, 8> per_type_hm_T2_sum = {};
    SphericalTensor hm_shielding_contribution;
    int n_rings_within_3A = 0;
    int n_rings_within_5A = 0;
    int n_rings_within_8A = 0;
    int n_rings_within_12A = 0;
    double mean_ring_distance = 0.0;
    double nearest_ring_atom_distance = 0.0;
    double G_iso_exp_sum = 0.0;
    std::array<double, 5> G_T2_exp_sum = {};
    double G_iso_var_8A = 0.0;
    SphericalTensor bs_shielding_contribution;

    // === Bond anisotropy totals (McConnellResult) ===
    // Forward surface: 10 source categories x {fixed, MOPAC bond-order}
    // channels, each packed as a full even SphericalTensor.
    std::array<std::array<SphericalTensor, kMcConnellChannelCount>,
               kMcConnellSourceCategoryCount> mcconnell_source_tensors = {};
    SphericalTensor mcconnell_peptide_co_rhombic;

    // Counts used by the near-field audit. The geometric filter can reject
    // near contacts; this preserves both accepted and rejected danger-zone
    // counts without making the pair list itself part of the atom payload.
    int mcconnell_near_field_accepted_lt3A = 0;
    int mcconnell_near_field_rejected_lt3A = 0;

    // Legacy projections retained for trajectory/Welford consumers. Values
    // are populated from the new fixed/BO channel model.
    double mcconnell_co_sum = 0.0;
    double mcconnell_cn_sum = 0.0;
    double mcconnell_sidechain_sum = 0.0;
    double mcconnell_aromatic_sum = 0.0;
    double mcconnell_co_nearest = 0.0;
    Vec3 nearest_CO_midpoint = Vec3::Zero();
    double nearest_CO_dist = 0.0;
    double nearest_CN_dist = 0.0;
    size_t nearest_CO_bond_index = SIZE_MAX;
    size_t nearest_CN_bond_index = SIZE_MAX;
    SphericalTensor T2_CO_nearest;
    SphericalTensor T2_CN_nearest;
    SphericalTensor T2_backbone_total;
    SphericalTensor T2_sidechain_total;
    SphericalTensor T2_aromatic_total;
    Vec3 dir_nearest_CO = Vec3::Zero();
    SphericalTensor mc_shielding_contribution;

    // === MOPAC bond-order-weighted anisotropy (MopacMcConnellResult) ===
    // Same kernel as McConnellResult, each bond weighted by MOPAC Wiberg order.
    double mopac_mc_co_sum = 0.0;
    double mopac_mc_cn_sum = 0.0;
    double mopac_mc_sidechain_sum = 0.0;
    double mopac_mc_aromatic_sum = 0.0;
    double mopac_mc_co_nearest = 0.0;
    double mopac_mc_nearest_CO_dist = 0.0;
    double mopac_mc_nearest_CN_dist = 0.0;
    SphericalTensor mopac_mc_T2_CO_nearest;
    SphericalTensor mopac_mc_T2_CN_nearest;
    SphericalTensor mopac_mc_T2_backbone_total;
    SphericalTensor mopac_mc_T2_sidechain_total;
    SphericalTensor mopac_mc_T2_aromatic_total;
    SphericalTensor mopac_mc_shielding_contribution;

    // === Coulomb field totals (CoulombResult) ===
    Vec3 coulomb_E_total = Vec3::Zero();
    Vec3 coulomb_E_backbone = Vec3::Zero();
    Vec3 coulomb_E_sidechain = Vec3::Zero();
    Vec3 coulomb_E_aromatic = Vec3::Zero();
    Mat3 coulomb_EFG_total = Mat3::Zero();
    SphericalTensor coulomb_EFG_total_spherical;
    Mat3 coulomb_EFG_backbone = Mat3::Zero();
    SphericalTensor coulomb_EFG_backbone_spherical;
    Mat3 coulomb_EFG_sidechain = Mat3::Zero();
    SphericalTensor coulomb_EFG_sidechain_spherical;
    Mat3 coulomb_EFG_aromatic = Mat3::Zero();
    SphericalTensor coulomb_EFG_aromatic_spherical;
    double coulomb_E_magnitude = 0.0;
    double coulomb_E_bond_proj = 0.0;
    double coulomb_E_backbone_frac = 0.0;  // projection of E_bb along E_total dir (V/A)
    double aromatic_E_magnitude = 0.0;
    double aromatic_E_bond_proj = 0.0;
    int aromatic_n_sidechain_atoms = 0;
    SphericalTensor coulomb_shielding_contribution;

    // === MOPAC Coulomb field totals (MopacCoulombResult) ===
    // Same kernel as CoulombResult but with MOPAC QM charges.
    // Units: V/A (E-field), V/A^2 (EFG).
    Vec3 mopac_coulomb_E_total = Vec3::Zero();
    Vec3 mopac_coulomb_E_backbone = Vec3::Zero();
    Vec3 mopac_coulomb_E_sidechain = Vec3::Zero();
    Vec3 mopac_coulomb_E_aromatic = Vec3::Zero();
    Mat3 mopac_coulomb_EFG_total = Mat3::Zero();
    SphericalTensor mopac_coulomb_EFG_total_spherical;
    Mat3 mopac_coulomb_EFG_backbone = Mat3::Zero();
    SphericalTensor mopac_coulomb_EFG_backbone_spherical;
    Mat3 mopac_coulomb_EFG_sidechain = Mat3::Zero();
    SphericalTensor mopac_coulomb_EFG_sidechain_spherical;
    Mat3 mopac_coulomb_EFG_aromatic = Mat3::Zero();
    SphericalTensor mopac_coulomb_EFG_aromatic_spherical;
    double mopac_coulomb_E_magnitude = 0.0;
    double mopac_coulomb_E_bond_proj = 0.0;
    double mopac_coulomb_E_backbone_frac = 0.0;
    double mopac_coulomb_aromatic_E_magnitude = 0.0;
    double mopac_coulomb_aromatic_E_bond_proj = 0.0;
    int mopac_coulomb_aromatic_n_sidechain_atoms = 0;
    std::uint8_t mopac_coulomb_efield_clamp_mask = 0;
    double mopac_coulomb_efield_clamp_scale = 1.0;
    SphericalTensor mopac_coulomb_shielding_contribution;

    // === APBS reaction fields (ApbsFieldResult) ===
    // Canonical fields are total PB minus the homogeneous-vacuum reference.
    // Units: V/A (E-field), V/A^2 (EFG), using the configured temperature.
    Vec3 apbs_efield = Vec3::Zero();
    Mat3 apbs_efg = Mat3::Zero();
    SphericalTensor apbs_efg_spherical;
    double apbs_phi = 0.0;  // reaction potential, V
    std::uint8_t apbs_efield_clamp_mask = 0;
    double apbs_efield_clamp_scale = 1.0;
    // Bit mask for the retained legacy APBS nonfinite sanitizers:
    // bit 0 reaction E, bit 1 reaction EFG, bit 2 total E, bit 3 total EFG.
    // Every set bit is accompanied by an OperationLog warning.
    std::uint8_t apbs_nonfinite_sanitizer_mask = 0;
    Vec3 apbs_efield_total_diagnostic = Vec3::Zero();
    Mat3 apbs_efg_total_diagnostic = Mat3::Zero();
    SphericalTensor apbs_efg_total_diagnostic_spherical;

    // === H-bond properties (HBondResult) ===
    double hbond_nearest_dist = 0.0;
    Vec3 hbond_nearest_dir = Vec3::Zero();
    double hbond_inv_d3 = 0.0;
    bool hbond_is_backbone = false;
    int hbond_count_within_3_5A = 0;
    double hbond_mcconnell_scalar = 0.0;  // Σ (3cos²θ−1)/r³ over contributing H-bonds
    bool hbond_is_donor = false;
    bool hbond_is_acceptor = false;

    // === Per-type PiQuadrupole accumulation (PiQuadrupoleResult) ===
    std::array<double, 8> per_type_pq_scalar_sum = {};           // (3cos²θ-1)/r⁴ per ring type

    // === Per-type Dispersion accumulation (DispersionResult) ===
    std::array<double, 8> per_type_disp_scalar_sum = {};           // 1/r⁶ per ring type

    // === Graph topology (MolecularGraphResult) ===
    int graph_dist_ring = -1;
    int graph_dist_N = -1;
    int graph_dist_O = -1;
    double eneg_sum_1 = 0.0;
    double eneg_sum_2 = 0.0;
    int n_pi_bonds_3 = 0;
    bool is_conjugated = false;
    int bfs_to_nearest_ring_atom = -1;
    int nearest_ring_atom_index = -1;
    double bfs_decay = 0.0;

    // === ORCA DFT shielding (OrcaShieldingResult) ===
    // Per-conformation: THIS protein's DFT shielding at this atom.
    // WT and mutant are separate Proteins with separate conformations.
    // Comparison is done by MutationDeltaResult, not here.
    Mat3 orca_shielding_total = Mat3::Zero();
    SphericalTensor orca_shielding_total_spherical;
    Mat3 orca_shielding_diamagnetic = Mat3::Zero();
    SphericalTensor orca_shielding_diamagnetic_spherical;
    Mat3 orca_shielding_paramagnetic = Mat3::Zero();
    SphericalTensor orca_shielding_paramagnetic_spherical;
    bool has_orca_shielding = false;

    // === Larsen H-bond contributions (LarsenHBondShieldingResult) ===
    //
    // Per Larsen 2015 Eqs. 4-5: Δσ_HB^i (amide donor) and Δσ_HαB^i (Hα donor)
    // are H-bond contribution terms read from DFT grid lookups against
    // the 6 (donor × acceptor) ProCS15 archives. Each H-bond pair
    // contributes a 1° term (donor-residue effect) and a 2° term
    // (acceptor's residue i+1, for NMA acceptor only). Per-atom-type
    // dispatch follows Larsen Table 2; see LarsenContribDispatch in
    // src/LarsenHBondShieldingResult.h.
    //
    // Methods accumulate (feedback_methods_accumulate): these fields
    // coexist with HBondResult's kernel-form output for the
    // amide-H/backbone-O subset. Per-atom-type residuals between the
    // two are themselves thesis-reportable.
    //
    // larsen_hbond_shielding_tensor is the SUM over all contribution
    // classes that apply at this atom (1°HB + 2°HB + 1°HαB + 2°HαB) per
    // Larsen Table 2. Per-class fields hold each contribution separately
    // for ML feature stratification. Tensors are in protein lab frame
    // (already rotated from canonical donor frame via
    // RotateTensorToProteinLabFrame).
    //
    // Pattern 11 (PATTERNS.md): every tensor is stored as BOTH Mat3
    // AND SphericalTensor — every Mat3 below has a `*_spherical`
    // companion. Downstream consumers never decompose at point of use.
    Mat3            larsen_hbond_shielding_tensor    = Mat3::Zero();
    SphericalTensor larsen_hbond_shielding_spherical;
    Mat3            larsen_hbond_1pHB_tensor         = Mat3::Zero();
    SphericalTensor larsen_hbond_1pHB_spherical;
    Mat3            larsen_hbond_2pHB_tensor         = Mat3::Zero();
    SphericalTensor larsen_hbond_2pHB_spherical;
    Mat3            larsen_hbond_1pHaB_tensor        = Mat3::Zero();
    SphericalTensor larsen_hbond_1pHaB_spherical;
    Mat3            larsen_hbond_2pHaB_tensor        = Mat3::Zero();
    SphericalTensor larsen_hbond_2pHaB_spherical;
    // Cβ diagnostic — Larsen Table 2 says Cβ gets NO HB contribution;
    // we compute and emit it anyway to verify the parser→loader→
    // rotation pipeline produces near-zero where the physics expects
    // it (reality check per feedback_methods_accumulate).
    Mat3            larsen_hbond_diagnostic_CB       = Mat3::Zero();
    SphericalTensor larsen_hbond_diagnostic_CB_spherical;
    // Water term: 2.07 ppm isotropic on amide H atoms that received
    // ZERO H-bond pair contributions (Larsen Δσ_w, NMA-water complex
    // value). Zero for non-HN atoms and for HN atoms with any pair.
    double larsen_hbond_water_term = 0.0;
    // Pair count contributing to this atom (across all 4 classes).
    int  larsen_hbond_n_pairs = 0;
    // True iff any of the 8 trilinear corner cells in any grid lookup
    // serving this atom was an imputed (nearest-neighbour-filled) bin.
    bool larsen_hbond_any_corner_imputed = false;

    // === Prediction fields ===
    double predicted_T0 = 0.0;
    std::array<double, 5> predicted_T2 = {};
    double confidence = 0.0;
    HeuristicTier tier = HeuristicTier::SILENT;

    // === AIMNet2 neural network results (AIMNet2Result) ===
    // Hirshfeld charge from the loaded AIMNet2 model (elementary charges)
    double aimnet2_charge = 0.0;
    // Learned electronic structure embedding (256 dims, geometry-dependent).
    // float32: native torch precision. No upshift to double.
    std::array<float, AIMNET2_AIM_DIMS> aimnet2_aim = {};
    // Fixed, committed element-specific projection of aimnet2_aim. Written
    // once by AIMNet2Result::Compute; every emitter only reads this field.
    std::array<float, AIMNET2_AIM_PROJECTION_DIMS> aimnet2_aim_projection = {};
    double aimnet2_energy_mlp = 0.0;
    double aimnet2_energy_shifted_local = 0.0;
    double aimnet2_d3_e_disp_atom = 0.0;
    double aimnet2_d3_cn = 0.0;
    std::array<double, 3> aimnet2_d3_c6_stats = {};
    // Coulomb E/EFG from AIMNet2 charges — same kernel as CoulombResult.
    Vec3 aimnet2_E_total = Vec3::Zero();
    Vec3 aimnet2_E_backbone = Vec3::Zero();
    Vec3 aimnet2_E_sidechain = Vec3::Zero();
    Vec3 aimnet2_E_aromatic = Vec3::Zero();
    Mat3 aimnet2_EFG_total = Mat3::Zero();
    SphericalTensor aimnet2_EFG_total_spherical;
    Mat3 aimnet2_EFG_backbone = Mat3::Zero();
    SphericalTensor aimnet2_EFG_backbone_spherical;
    Mat3 aimnet2_EFG_sidechain = Mat3::Zero();
    SphericalTensor aimnet2_EFG_sidechain_spherical;
    Mat3 aimnet2_EFG_aromatic = Mat3::Zero();
    SphericalTensor aimnet2_EFG_aromatic_spherical;
    SphericalTensor aimnet2_shielding_contribution;

    // Charge-response gradient via autograd (AIMNet2ChargeResponseGradientResult).
    // Vector is dL/d(r_i) where L = sum_j q_j^2 over non-sentinel atoms;
    // scalar is its L2 norm.
    Vec3 aimnet2_charge_response_gradient_vector = Vec3::Zero();
    double aimnet2_charge_response_gradient_scalar = 0.0;

    // === Planar geometry (PlanarGeometryResult) ===
    // Non-negative out-of-plane displacement magnitude (Å). NaN means
    // non-applicable or invalid/degenerate; pyramidalization_valid
    // disambiguates real zero from absence, and center_type stores the
    // PlanarGroupKind enum value from the typed semantic substrate.
    double pyramidalization = 0.0;
    std::int8_t pyramidalization_valid = 0;
    std::int8_t pyramidalization_center_type = 0;

    // === Solvent-accessible surface area (SasaResult) ===
    double atom_sasa = 0.0;  // Shrake-Rupley SASA (A^2)
    Vec3 sasa_normal = Vec3::Zero();  // outward surface normal from non-occluded test points

    // === Explicit solvent fields (WaterFieldResult) ===
    // Electric field at this atom from water charges within cutoff (V/A)
    Vec3 water_efield = Vec3::Zero();
    std::int8_t water_efield_clamp_mask = 0;
    double water_efield_clamp_scale = 1.0;
    // Electric field gradient from water (V/A^2)
    Mat3 water_efg = Mat3::Zero();
    SphericalTensor water_efg_spherical;
    // First hydration shell (water O within 3.5A): E-field and EFG
    Vec3 water_efield_first = Vec3::Zero();
    std::int8_t water_efield_first_clamp_mask = 0;
    double water_efield_first_clamp_scale = 1.0;
    Mat3 water_efg_first = Mat3::Zero();
    SphericalTensor water_efg_first_spherical;
    // Shell counts
    int water_n_first = 0;    // water O within 3.5A
    int water_n_second = 0;   // water O within 3.5-5.5A

    // === Hydration shell geometry (HydrationShellResult) ===
    double half_shell_asymmetry = 0.0;  // fraction exposed vs buried
    double mean_water_dipole_cos = 0.0; // water orientation order parameter
    double nearest_ion_distance = std::numeric_limits<double>::infinity();  // distance to closest ion (A), inf = none within cutoff
    double nearest_ion_charge = 0.0;    // charge of nearest ion (e)

    // === Hydration geometry — SASA-normal reference frame (HydrationGeometryResult) ===
    Vec3 water_dipole_vector = Vec3::Zero();   // net first-shell water dipole (Debye-like, unnormalised)
    Vec3 water_surface_normal = Vec3::Zero();  // copy of sasa_normal for this block
    double sasa_half_shell_asymmetry = 0.0;    // exposed/total using SASA normal (not COM)
    double sasa_dipole_alignment = 0.0;        // cos(net dipole, SASA normal)
    double sasa_dipole_coherence = 0.0;        // legacy mean net dipole |Σ dᵢ| / n, e_Angstrom
    double sasa_dipole_order_parameter = 0.0;  // |Σ dᵢ| / Σ|dᵢ|, dimensionless [0, 1]
    int sasa_first_shell_count = 0;            // first-shell water O count

    // === EEQ charges (EeqResult) ===
    // project-local QEq/EEQ-style charge-equilibration model with error-function coordination number, CN-dependent electronegativity shift, Gaussian self term, and Ohno-Klopman off-diagonal kernel; parameters are in-repo/project-local and are not a validated dftd4/multicharge port.
    double eeq_charge = 0.0;  // geometry-dependent partial charge (elementary charges)
    double eeq_cn = 0.0;      // coordination number used to compute eeq_charge
    double eeq_chi_eff = 0.0; // CN-shifted electronegativity (atomic units)
    double eeq_eta = 0.0;     // element chemical hardness (atomic units)
    double eeq_self_hardness_diag = 0.0; // eta + Gaussian self term

    // === EEQ-charge Coulomb fields (EeqCoulombResult) ===
    Vec3 eeq_coulomb_E_total = Vec3::Zero();
    Vec3 eeq_coulomb_E_backbone = Vec3::Zero();
    Vec3 eeq_coulomb_E_sidechain = Vec3::Zero();
    Vec3 eeq_coulomb_E_aromatic = Vec3::Zero();
    Mat3 eeq_coulomb_EFG_total = Mat3::Zero();
    SphericalTensor eeq_coulomb_EFG_total_spherical;
    Mat3 eeq_coulomb_EFG_backbone = Mat3::Zero();
    SphericalTensor eeq_coulomb_EFG_backbone_spherical;
    Mat3 eeq_coulomb_EFG_sidechain = Mat3::Zero();
    SphericalTensor eeq_coulomb_EFG_sidechain_spherical;
    Mat3 eeq_coulomb_EFG_aromatic = Mat3::Zero();
    SphericalTensor eeq_coulomb_EFG_aromatic_spherical;
    double eeq_coulomb_E_magnitude = 0.0;
    double eeq_coulomb_E_bond_proj = 0.0;
    double eeq_coulomb_E_backbone_frac = 0.0;
    double eeq_coulomb_aromatic_E_magnitude = 0.0;
    double eeq_coulomb_aromatic_E_bond_proj = 0.0;
    int eeq_coulomb_aromatic_n_sidechain_atoms = 0;
    SphericalTensor eeq_coulomb_shielding_contribution;

    // === DemoResult fields (Pass 0) ===
    double demo_nearest_ring_distance = 0.0;
    Vec3 demo_nearest_ring_direction = Vec3::Zero();

private:
    explicit ConformationAtom(Vec3 pos) : position_(pos) {}
    const Vec3 position_;
};

}  // namespace nmr
