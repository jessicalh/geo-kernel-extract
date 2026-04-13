#pragma once
//
// ConformationAtom: per-atom computed data within a ProteinConformation.
//
// Private constructor -- only ProteinConformation can create these.
// Position is const after construction. All computed fields are public
// and default-initialised. The singleton guarantee means each field
// has exactly one writer (the ConformationResult that owns it).
//
// Identity (element, bonds, residue) goes through the Protein back-pointer,
// NOT duplicated here.
//
// ALL fields from OBJECT_MODEL.md are present, even if this pass
// does not populate them.
//

#include "Types.h"
#include <vector>
#include <array>

namespace nmr {

// AIMNet2 aim embedding dimensionality.
// Defined here (not AIMNet2Result.h) because ConformationAtom needs it
// for the array extent, and ConformationAtom.h must not include torch headers.
static constexpr size_t AIMNET2_AIM_DIMS = 256;

class ProteinConformation;

// Per-atom, per-ring structured result (placeholder for future ring calculators)
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
    Mat3 quad_tensor = Mat3::Zero();
    SphericalTensor quad_spherical;
    double quad_scalar = 0.0;
    Mat3 chi_tensor = Mat3::Zero();
    SphericalTensor chi_spherical;
    double chi_scalar = 0.0;
    Mat3 disp_tensor = Mat3::Zero();
    SphericalTensor disp_spherical;
    double disp_scalar = 0.0;
    int disp_contacts = 0;
    double gaussian_density = 0.0;
    double cos_phi = 1.0;   // azimuthal: cos of in-plane angle to vertex 0
    double sin_phi = 0.0;   // azimuthal: sin of in-plane angle to vertex 0
};

// Per-atom, per-bond structured result (placeholder for McConnell)
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

    // === Charges and radii (ChargeAssignmentResult) ===
    double partial_charge = 0.0;
    double vdw_radius = 0.0;

    // === MOPAC semiempirical results (MopacResult) ===
    double mopac_charge = 0.0;            // Mulliken charge (elementary charges)
    double mopac_s_pop = 0.0;             // s-orbital population
    double mopac_p_pop = 0.0;             // p-orbital population
    double mopac_valency = 0.0;           // sum of Wiberg bond orders (CSC diagonal)
    std::vector<MopacBondNeighbour> mopac_bond_neighbours;  // sorted descending by order

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
    std::array<std::array<double, 5>, 8> per_type_G_T2_sum = {};
    std::array<double, 8> per_type_hm_T0_sum = {};
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
    double mcconnell_co_sum = 0.0;
    double mcconnell_cn_sum = 0.0;
    double mcconnell_sidechain_sum = 0.0;
    double mcconnell_aromatic_sum = 0.0;
    double mcconnell_co_nearest = 0.0;
    Vec3 nearest_CO_midpoint = Vec3::Zero();
    double nearest_CO_dist = 0.0;
    double nearest_CN_dist = 0.0;
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
    Mat3 coulomb_EFG_aromatic = Mat3::Zero();
    SphericalTensor coulomb_EFG_aromatic_spherical;
    Vec3 coulomb_E_solvent = Vec3::Zero();
    Mat3 coulomb_EFG_solvent = Mat3::Zero();
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
    Mat3 mopac_coulomb_EFG_aromatic = Mat3::Zero();
    SphericalTensor mopac_coulomb_EFG_aromatic_spherical;
    double mopac_coulomb_E_magnitude = 0.0;
    double mopac_coulomb_E_bond_proj = 0.0;
    double mopac_coulomb_E_backbone_frac = 0.0;
    SphericalTensor mopac_coulomb_shielding_contribution;

    // === APBS solvated fields (ApbsFieldResult) ===
    // Units: V/A (E-field), V/A^2 (EFG). Converted from APBS kT/(e*A)
    // by KT_OVER_E_298K. Same units as CoulombResult for direct comparison.
    Vec3 apbs_efield = Vec3::Zero();
    Mat3 apbs_efg = Mat3::Zero();
    SphericalTensor apbs_efg_spherical;

    // === H-bond properties (HBondResult) ===
    double hbond_nearest_dist = 0.0;
    Vec3 hbond_nearest_dir = Vec3::Zero();
    Mat3 hbond_nearest_tensor = Mat3::Zero();
    SphericalTensor hbond_nearest_spherical;
    double hbond_inv_d3 = 0.0;
    bool hbond_is_backbone = false;
    int hbond_count_within_3_5A = 0;
    bool hbond_is_donor = false;
    bool hbond_is_acceptor = false;
    SphericalTensor hbond_shielding_contribution;

    // === Ring-based shielding contributions ===
    SphericalTensor piquad_shielding_contribution;
    SphericalTensor ringchi_shielding_contribution;
    SphericalTensor disp_shielding_contribution;

    // === Per-type PiQuadrupole accumulation (PiQuadrupoleResult) ===
    std::array<double, 8> per_type_pq_scalar_sum = {};           // (3cos²θ-1)/r⁴ per ring type
    std::array<std::array<double, 5>, 8> per_type_pq_T2_sum = {}; // EFG T2 per ring type

    // === Per-type Dispersion accumulation (DispersionResult) ===
    std::array<double, 8> per_type_disp_scalar_sum = {};           // 1/r⁶ per ring type
    std::array<std::array<double, 5>, 8> per_type_disp_T2_sum = {}; // disp T2 per ring type

    // === Graph topology (MolecularGraphResult) ===
    int graph_dist_ring = -1;
    int graph_dist_N = -1;
    int graph_dist_O = -1;
    double eneg_sum_1 = 0.0;
    double eneg_sum_2 = 0.0;
    int n_pi_bonds_3 = 0;
    bool is_conjugated = false;
    int bfs_to_nearest_ring_atom = -1;
    double bfs_decay = 0.0;

    // === ORCA DFT shielding (OrcaShieldingResult) ===
    // Per-conformation: THIS protein's DFT shielding at this atom.
    // WT and mutant are separate Proteins with separate conformations.
    // Comparison is done by MutantProteinConformationComparison, not here.
    Mat3 orca_shielding_total = Mat3::Zero();
    SphericalTensor orca_shielding_total_spherical;
    Mat3 orca_shielding_diamagnetic = Mat3::Zero();
    SphericalTensor orca_shielding_diamagnetic_spherical;
    Mat3 orca_shielding_paramagnetic = Mat3::Zero();
    SphericalTensor orca_shielding_paramagnetic_spherical;
    bool has_orca_shielding = false;

    // === Predictions (PredictionResult) ===
    double predicted_T0 = 0.0;
    std::array<double, 5> predicted_T2 = {};
    double confidence = 0.0;
    HeuristicTier tier = HeuristicTier::SILENT;

    // === AIMNet2 neural network results (AIMNet2Result) ===
    // Hirshfeld charge from AIMNet2 wB97M model (elementary charges)
    double aimnet2_charge = 0.0;
    // Learned electronic structure embedding (256 dims, geometry-dependent)
    std::array<double, AIMNET2_AIM_DIMS> aimnet2_aim = {};
    // Coulomb EFG from AIMNet2 charges — same kernel as CoulombResult
    Mat3 aimnet2_EFG_total = Mat3::Zero();
    SphericalTensor aimnet2_EFG_total_spherical;
    Mat3 aimnet2_EFG_backbone = Mat3::Zero();
    SphericalTensor aimnet2_EFG_backbone_spherical;
    Mat3 aimnet2_EFG_aromatic = Mat3::Zero();
    SphericalTensor aimnet2_EFG_aromatic_spherical;
    SphericalTensor aimnet2_shielding_contribution;

    // === Solvent-accessible surface area (SasaResult) ===
    double atom_sasa = 0.0;  // Shrake-Rupley SASA (A^2)
    Vec3 sasa_normal = Vec3::Zero();  // outward surface normal from non-occluded test points

    // === Explicit solvent fields (WaterFieldResult) ===
    // Electric field at this atom from water charges within cutoff (V/A)
    Vec3 water_efield = Vec3::Zero();
    // Electric field gradient from water (V/A^2)
    Mat3 water_efg = Mat3::Zero();
    SphericalTensor water_efg_spherical;
    // First hydration shell (water O within 3.5A): E-field and EFG
    Vec3 water_efield_first = Vec3::Zero();
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
    double sasa_dipole_coherence = 0.0;        // |Σ dᵢ| / n — ordered vs random
    int sasa_first_shell_count = 0;            // first-shell water O count

    // === EEQ charges (EeqResult — Caldeweyher 2019) ===
    double eeq_charge = 0.0;  // geometry-dependent partial charge (elementary charges)
    double eeq_cn = 0.0;      // coordination number used to compute eeq_charge

    // === DemoResult fields (Pass 0) ===
    double demo_nearest_ring_distance = 0.0;
    Vec3 demo_nearest_ring_direction = Vec3::Zero();

private:
    explicit ConformationAtom(Vec3 pos) : position_(pos) {}
    const Vec3 position_;
};

}  // namespace nmr
