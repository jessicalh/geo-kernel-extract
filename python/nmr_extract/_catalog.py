"""Producer SDK format contract for NMR-extract NPY files.

This is the SDK source of truth for producer-owned arrays on disk. If the
nmr_extract producer adds a new array, add an entry here. Reader-owned
rediscover/statistics sidecars live in h5-reader/src/rediscover/
ReaderOutputCatalog.h.

Generated from the WriteFeatures methods in src/*Result.cpp and
src/ConformationResult.cpp.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

import numpy as np

from ._tensors import (
    ShieldingTensor,
    EFGTensor,
    PositionField,
    VectorField,
    MagneticVectorField,
    PerRingTypeT0,
    PerRingTypeT1,
    PerRingTypeT2,
    PerBondCategoryT2,
    RingCounts,
    McConnellNearFieldCounts,
    McConnellScalars,
    CoulombScalars,
    HBondScalars,
    DsspScalars,
    MopacScalars,
    MopacGlobal,
    MopacAtomPopulations,
    MopacAtomicOrbitalPopulations,
    MopacAtomicOrbitalPopulationTotals,
    MopacUniqueBondOrders,
    MopacTopologyBondOrdersFull,
    BondOrders,
    DeltaScalars,
    DeltaAPBS,
    DeltaRingProximity,
    AIMNet2Charges,
    AIMNet2AimEmbedding,
    AIMNet2ChargeResponseGradient,
)
from ._ring import RingContributions, RingGeometry, RingPairGeometry


ALLOWED_NATIVE_AXES = frozenset({
    "atom",
    "residue",
    "aromatic_ring",
    "saturated_ring",
    "ring",
    "ring_contribution_pair",
    "bond",
    "ring_membership",
    "mutation_match_pair",
    "protein",
    "mopac_bond_neighbor_pair",
    "mopac_unique_pair",
    "atom_neighbor_pair",
    "ring_pair",
    "larsen_hbond_pair",
    "hbond_pair",
    "larsen_sidechain_donor_pair",
    "protein_water_hbond_pair",
    "sidechain_co_source",
    "mopac_atomic_orbital_row",
})


@dataclass(frozen=True)
class ArraySpec:
    """Metadata for one NPY file produced by the C++ extractor.

    Topology-sidecar metadata (2026-05-13, resolves
    TENTATIVE_OUTSTANDING_ISSUES OI-016):

    * ``native_axis`` declares which axis this array's rows are
      indexed along: ``atom`` / ``residue`` / ``aromatic_ring`` /
      ``saturated_ring`` / ``ring`` / ``ring_contribution_pair`` /
      ``bond`` / ``ring_membership`` / ``mutation_match_pair`` /
      ``larsen_hbond_pair`` / ``protein`` (plus the other sparse axes in
      ``ALLOWED_NATIVE_AXES``).
      R / Python analysis must read this metadata column rather than infer axis
      from filename.

    * ``irreps`` is legacy/non-authoritative for raw project-native
      tensors. ``tensor_basis``, ``tensor_component_order``,
      ``tensor_frame``, ``structural_zero_components``, and
      ``e3nn_export`` describe the producer payload. Consumers that need
      e3nn must call the explicit conversion APIs in ``_tensors.py``.

    * ``units`` carries the SI / NMR-standard unit string for
      consumers needing to compare values across calculators
      (``"ppm"``, ``"V/A"``, ``"V/A^2"``, ``"e"``, ``"kJ/mol"``,
      ``"Å"``, ``"radians"``, ``"degrees"``). Empty for
      dimensionless / categorical / structured-dtype arrays.

    * ``sign_convention`` documents the explicit physical sign
      convention where one applies (notably ring-current shielding
      kernels: ``"σ_ab = -dB_a^sec / dB_{0,b}"``). Empty otherwise.

    * ``tensor_rank`` is the rank of one row's tensor representation:
      0 (scalar / vector-of-scalars), 1 (Vec3), 2 (Mat3 / SphericalTensor).

    * ``parity`` is ``"even"`` / ``"odd"`` under spatial inversion, or
      ``"mixed"`` when a serialized payload has no single homogeneous
      improper-transform law (including mixed columns and categorical
      outcomes conditioned on chirality);
      shieldings (rank 2 even), axial B fields, and most scalars are
      even. Polar vector fields like E are odd. ``mixed`` is reserved for a
      single serialized payload whose columns/local components have genuinely
      different improper-transform laws and therefore cannot truthfully carry
      one homogeneous parity.

    * ``coordinate_frame`` names the frame in which serialized components
      are expressed. ``transformation`` is the authoritative physical law,
      including affine positions and column-level laws for mixed tables.
      ``validity`` names masks, NaN rules, physical-zero rules, and known
      legacy missing-value sentinels. These fields deliberately complement,
      rather than reinterpret, the legacy ``irreps`` / ``parity`` fields.

    * ``mechanism`` is a thesis-narrative grouping over physics:
      ``ring_current`` / ``ring_efg`` / ``ring_dispersion`` /
      ``bond_anisotropy`` / ``electrostatic_efg`` / ``hbond_kernel`` /
      ``hbond_grid`` / ``secondary_structure`` / ``solvation`` /
      ``charges`` / ``quantum_reference`` / ``mutation_delta`` /
      ``topology`` / ``gromacs_runtime`` / ``geometry``. R analysis
      should read this column rather than regex over feature names
      (see learn/R/stage1_bmrb_dimension_independence.R — current
      offender per the codex sidecar contract).
    """
    stem: str                   # filename without .npy
    group: str                  # logical group (e.g. "biot_savart")
    wrapper: type               # Python wrapper class
    cols: Optional[int]         # expected last dim (None = variable or 1D)
    required: bool              # False for optional calculators
    description: str
    # Topology-sidecar metadata. Defaults are minimal-information so a
    # legacy entry without explicit values is still loadable; populated
    # entries override.
    native_axis: str = "atom"
    irreps: str = ""
    units: str = ""
    sign_convention: str = ""
    tensor_rank: int = 0
    parity: str = "even"
    mechanism: str = "metadata"
    tensor_basis: str = ""
    tensor_component_order: str = ""
    tensor_frame: str = ""
    structural_zero_components: str = ""
    e3nn_export: str = ""
    # Producer-side physical scaling contract for geometry-only kernels.
    # Empty when an array is already in its final physical units.
    scaling_contract: str = ""
    # Directional-output freeze contract. Empty only for arrays whose scalar
    # or categorical behaviour needs no directional qualification.
    coordinate_frame: str = ""
    transformation: str = ""
    validity: str = ""


# fmt: off
# Common metadata strings — name once, use across entries of the same physics.
_SHIELD_IRREPS = "0e + 1e + 2e"
_SHIELD_SIGN   = "σ_ab = -dB_a^sec/dB_{0,b}"
_FULL9_BASIS = "project_native_full9_spherical_tensor_v1"
_T2_BASIS = "project_native_t2_isometric_real_tesseral_v1"
_FULL9_ORDER = "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
_T2_ORDER = "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
_EFG_STRUCTURAL_ZEROS = "T0,T1_x,T1_y,T1_z"
_E3NN_EXPORT = (
    "raw project tensor; call to_e3nn()/to_e3nn_T2() or "
    "project_t2_to_e3nn() before using e3nn Irreps"
)
_BS_INTENSITY_SCALING = (
    "unit-current geometry kernel; multiply each per-ring-type channel by "
    "extraction_manifest.json::feature_metadata.ring_current."
    "resolved_intensity_nA_per_T in the declared ring_type_order and sum "
    "over ring type to obtain physical ppm shielding"
)
_T2_TENSOR_METADATA = dict(
    tensor_basis=_T2_BASIS,
    tensor_component_order=_T2_ORDER,
    tensor_frame="conformation_cartesian_xyz",
    structural_zero_components=_EFG_STRUCTURAL_ZEROS,
    e3nn_export=_E3NN_EXPORT,
)

CATALOG: dict[str, ArraySpec] = {s.stem: s for s in [
    # ── Identity (ConformationResult.cpp) ────────────────────────
    ArraySpec("pos",              "identity",   PositionField,     3,    True,  "Atom positions (A)",
              native_axis="atom", units="Å", parity="odd", mechanism="topology"),
    ArraySpec("element",          "identity",   np.ndarray,        None, True,  "Atomic number (int32)",
              native_axis="atom", mechanism="topology"),
    ArraySpec("residue_index",    "identity",   np.ndarray,        None, True,  "Residue index (int32)",
              native_axis="atom", mechanism="topology"),
    ArraySpec("residue_type",     "identity",   np.ndarray,        None, True,  "Residue type enum (int32)",
              native_axis="atom", mechanism="topology"),
    ArraySpec("ring_contributions","identity",  RingContributions, 40,   True,  "Per-(atom,ring) pair contributions",
              native_axis="ring_contribution_pair", irreps=_SHIELD_IRREPS, units="",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="ring_current"),
    ArraySpec("ring_direction_to_center", "identity", VectorField, 3, False, "Sparse per-(atom,ring) rows aligned to ring_contributions: RingNeighbourhood.direction_to_center vector",
              native_axis="ring_contribution_pair", irreps="1o", tensor_rank=1, parity="odd", mechanism="geometry"),
    ArraySpec("ring_geometry",    "identity",   RingGeometry,      10,   True,  "Per-ring geometry reference",
              native_axis="aromatic_ring", units="Å", mechanism="topology"),
    ArraySpec("ring_pair_geometry", "identity", RingPairGeometry, 13, True, "All i<j aromatic-ring pair geometry [ring/residue/type ids, center distance, normal relations, signed normal offsets, in-plane slip, fused flag]",
              native_axis="ring_pair", units="mixed_index_A_dimensionless", mechanism="geometry"),

    # ── Covalent/local geometry (GeometryResult and
    #    LocalBackboneGeometryResult) ─────────────────────────────
    ArraySpec("bond_length", "geometry", np.ndarray, None, False, "Per-topology-bond length in bonds.npy row order",
              native_axis="bond", units="Å", mechanism="geometry"),
    ArraySpec("bond_direction", "geometry", VectorField, 3, False, "Unit vector from bonds.atom_index_a to bonds.atom_index_b in conformation Cartesian xyz",
              native_axis="bond", irreps="1o", tensor_rank=1, parity="odd", mechanism="geometry"),
    ArraySpec("bond_geometry_valid", "geometry", np.ndarray, None, False, "Per-bond uint8 validity: finite non-zero bond vector",
              native_axis="bond", mechanism="geometry"),
    ArraySpec("tau_N_CA_C", "local_backbone_geometry", np.ndarray, None, False, "Residue-local N-CA-C valence angle",
              native_axis="residue", units="radians", mechanism="geometry"),
    ArraySpec("tau_N_CA_C_valid", "local_backbone_geometry", np.ndarray, None, False, "uint8 validity for tau_N_CA_C",
              native_axis="residue", mechanism="geometry"),
    ArraySpec("angle_N_CA_CB", "local_backbone_geometry", np.ndarray, None, False, "Residue-local N-CA-CB valence angle",
              native_axis="residue", units="radians", mechanism="geometry"),
    ArraySpec("angle_N_CA_CB_valid", "local_backbone_geometry", np.ndarray, None, False, "uint8 validity for angle_N_CA_CB",
              native_axis="residue", mechanism="geometry"),
    ArraySpec("angle_CB_CA_C", "local_backbone_geometry", np.ndarray, None, False, "Residue-local CB-CA-C valence angle",
              native_axis="residue", units="radians", mechanism="geometry"),
    ArraySpec("angle_CB_CA_C_valid", "local_backbone_geometry", np.ndarray, None, False, "uint8 validity for angle_CB_CA_C",
              native_axis="residue", mechanism="geometry"),
    ArraySpec("angle_Cprev_N_CA", "local_backbone_geometry", np.ndarray, None, False, "Covalent-backbone predecessor C to N-CA valence angle",
              native_axis="residue", units="radians", mechanism="geometry"),
    ArraySpec("angle_Cprev_N_CA_valid", "local_backbone_geometry", np.ndarray, None, False, "uint8 validity for angle_Cprev_N_CA",
              native_axis="residue", mechanism="geometry"),
    ArraySpec("angle_CA_C_Nnext", "local_backbone_geometry", np.ndarray, None, False, "CA-C to covalent-backbone successor N valence angle",
              native_axis="residue", units="radians", mechanism="geometry"),
    ArraySpec("angle_CA_C_Nnext_valid", "local_backbone_geometry", np.ndarray, None, False, "uint8 validity for angle_CA_C_Nnext",
              native_axis="residue", mechanism="geometry"),
    ArraySpec("cb_deviation", "local_backbone_geometry", np.ndarray, None, False, "Magnitude of observed_CB - established ideal_CB",
              native_axis="residue", units="Å", mechanism="geometry"),
    ArraySpec("cb_deviation_valid", "local_backbone_geometry", np.ndarray, None, False, "uint8 validity for cb_deviation",
              native_axis="residue", mechanism="geometry"),
    ArraySpec("cb_residual_vector", "local_backbone_geometry", VectorField, 3, False, "observed_CB - established ideal_CB in conformation Cartesian xyz (not a local frame)",
              native_axis="residue", irreps="1o", units="Å", tensor_rank=1, parity="odd", mechanism="geometry"),
    ArraySpec("cb_residual_vector_valid", "local_backbone_geometry", np.ndarray, None, False, "uint8 validity for all three cb_residual_vector components",
              native_axis="residue", mechanism="geometry"),

    # ── Enrichment (EnrichmentResult.cpp) ───────────────────────────
    ArraySpec("enrichment_role",          "enrichment", np.ndarray, None, False, "AtomRole enum per atom (int32)",
              mechanism="topology"),
    ArraySpec("enrichment_hybridisation", "enrichment", np.ndarray, None, False, "Hybridisation enum per atom (int32)",
              mechanism="topology"),
    ArraySpec("enrichment_flags",         "enrichment", np.ndarray, 8,    False, "Enrichment boolean flags as int8 columns: is_backbone, is_amide_H, is_alpha_H, is_methyl, is_aromatic_H, is_hbond_donor, is_hbond_acceptor, is_on_aromatic_residue",
              mechanism="topology"),
    ArraySpec("enrichment_parent_is_sp2", "enrichment", np.ndarray, None, False, "Per-atom uint8 flag: hydrogen parent is typed sp2/aromatic",
              mechanism="topology"),
    ArraySpec("semantic_polar_h_kind",     "enrichment", np.ndarray, None, False, "Raw typed PolarHKind enum per atom (uint8)",
              mechanism="topology"),
    ArraySpec("semantic_planar_group_kind","enrichment", np.ndarray, None, False, "Raw typed PlanarGroupKind enum per atom (uint8)",
              mechanism="topology"),
    ArraySpec("semantic_formal_charge",    "enrichment", np.ndarray, None, False, "Raw typed formal charge per atom (int8)",
              mechanism="topology"),
    ArraySpec("semantic_ring_position",    "enrichment", np.ndarray, None, False, "Raw canonical primary RingPositionLabel per atom (uint8); complete secondary/tertiary membership remains in atoms_category_info",
              mechanism="topology"),
    ArraySpec("semantic_locant",           "enrichment", np.ndarray, None, False, "Raw typed Locant enum per atom (uint8)",
              mechanism="topology"),
    ArraySpec("enrichment_donor_class",    "enrichment", np.ndarray, None, False, "Typed donor class per atom: exact PolarHKind projection (uint8)",
              mechanism="topology"),
    ArraySpec("enrichment_acceptor_class", "enrichment", np.ndarray, None, False, "Typed acceptor projection: 0 none, 1 backbone carbonyl, 2 sidechain amide carbonyl, 3 carboxylate, 4 hydroxyl/oxide, 5 unprotonated ring N, 6 neutral other N/O/S",
              mechanism="topology"),
    ArraySpec("enrichment_hybridisation_class", "enrichment", np.ndarray, None, False, "Semantic hybridisation projection using Hybridisation numeric codes: 0 sp, 1 sp2, 2 sp3, 3 unassigned",
              mechanism="topology"),

    # ── Spatial index (SpatialIndexResult.cpp) ────────────────────
    ArraySpec("spatial_neighbors", "spatial_index", np.ndarray, 6, False, "Directed per-atom neighbour rows [atom_i, atom_j, unit_dx, unit_dy, unit_dz, distance_A], no self rows, 15 A cutoff",
              native_axis="atom_neighbor_pair", units="mixed_index_unit_vector_A", mechanism="geometry"),

    # ── Molecular graph (MolecularGraphResult.cpp) ────────────────
    ArraySpec("molecular_graph_int", "molecular_graph", np.ndarray, 6, False, "Raw graph integers [distance_to_ring, distance_to_N, distance_to_O, pi_bonds_within_3, is_conjugated, nearest_ring_atom_index]",
              mechanism="topology"),
    ArraySpec("molecular_graph_float", "molecular_graph", np.ndarray, 3, False, "Graph floating features [electronegativity_sum_1, electronegativity_sum_2, exp(-ring_graph_distance/4)]",
              mechanism="topology"),
    ArraySpec("molecular_graph", "molecular_graph", np.ndarray, 9, False, "Compatibility graph table combining molecular_graph_int and molecular_graph_float",
              mechanism="topology"),

    # ── Force-field charge assignment (ChargeAssignmentResult.cpp) ──
    ArraySpec("ff_partial_charge", "charge_assignment", np.ndarray, None, False, "Force-field partial charge per atom (e)",
              units="e", mechanism="charges"),
    ArraySpec("ff_pb_radius",      "charge_assignment", np.ndarray, None, False, "Force-field Poisson-Boltzmann radius per atom (Å)",
              units="Å", mechanism="charges"),

    # ── Biot-Savart (BiotSavartResult.cpp) ───────────────────────
    ArraySpec("bs_shielding",     "biot_savart", ShieldingTensor,  9,    True,  "BS ring current shielding",
              irreps=_SHIELD_IRREPS, units="ppm_T_per_nA", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="ring_current",
              scaling_contract="unscaled compatibility sum; use bs_per_type channels with the manifest-resolved intensity vector for physical ppm"),
    ArraySpec("bs_per_type_T0",   "biot_savart", PerRingTypeT0,    8,    True,  "BS isotropic per ring type",
              irreps="0e", units="ppm_T_per_nA", mechanism="ring_current", scaling_contract=_BS_INTENSITY_SCALING),
    ArraySpec("bs_per_type_T1",   "biot_savart", PerRingTypeT1,    24,   True,  "BS T1 per ring type",
              irreps="1e", units="ppm_T_per_nA", tensor_rank=1, parity="even", mechanism="ring_current", scaling_contract=_BS_INTENSITY_SCALING),
    ArraySpec("bs_per_type_T2",   "biot_savart", PerRingTypeT2,    40,   True,  "BS T2 per ring type",
              irreps="2e", units="ppm_T_per_nA", tensor_rank=2, mechanism="ring_current", scaling_contract=_BS_INTENSITY_SCALING),
    ArraySpec("bs_total_B",       "biot_savart", MagneticVectorField, 3, True,  "BS total B-field vector",
              irreps="1e", units="T", tensor_rank=1, parity="even", mechanism="ring_current"),
    ArraySpec("bs_ring_B_field",  "biot_savart", MagneticVectorField, 3, True,  "BS per-(atom,ring) B-field vector",
              native_axis="ring_contribution_pair", irreps="1e", units="T", tensor_rank=1, parity="even", mechanism="ring_current"),
    ArraySpec("bs_ring_B_cylindrical", "biot_savart", MagneticVectorField, 3, True, "BS per-(atom,ring) B-field in ring cylindrical frame",
              native_axis="ring_contribution_pair", irreps="1e", units="T", tensor_rank=1, parity="even", mechanism="ring_current"),
    ArraySpec("bs_ring_counts",   "biot_savart", RingCounts,       4,    True,  "Ring proximity counts (3/5/8/12 A)",
              mechanism="ring_current"),

    # ── Haigh-Mallion (HaighMallionResult.cpp) ───────────────────
    ArraySpec("hm_shielding",     "haigh_mallion", ShieldingTensor, 9,   True,  "HM ring current shielding",
              irreps=_SHIELD_IRREPS, units="Angstrom^-1", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="ring_current"),
    ArraySpec("hm_per_type_T0",   "haigh_mallion", PerRingTypeT0,   8,   True,  "HM isotropic per ring type",
              irreps="0e", units="Angstrom^-1", mechanism="ring_current"),
    ArraySpec("hm_per_type_T1",   "haigh_mallion", PerRingTypeT1,   24,  True,  "HM T1 per ring type",
              irreps="1e", units="Angstrom^-1", tensor_rank=1, parity="even", mechanism="ring_current"),
    ArraySpec("hm_per_type_T2",   "haigh_mallion", PerRingTypeT2,   40,  True,  "HM T2 per ring type",
              irreps="2e", units="Angstrom^-1", tensor_rank=2, mechanism="ring_current"),
    ArraySpec("hm_ring_B_field",  "haigh_mallion", MagneticVectorField, 3, True, "HM per-(atom,ring) effective B-field vector",
              native_axis="ring_contribution_pair", irreps="1e", units="Angstrom^-1", tensor_rank=1, parity="even", mechanism="ring_current"),

    # ── Pi-Quadrupole (PiQuadrupoleResult.cpp) ───────────────────
    ArraySpec("pq_per_type_T0",   "pi_quadrupole", PerRingTypeT0,   8,   True,  "PQ Buckingham A-term scalar per ring type",
              irreps="0e", units="Angstrom^-4", mechanism="ring_efg"),
    ArraySpec("piquad_axial_scalar_per_type_T0", "pi_quadrupole", PerRingTypeT0, 8, False, "Preferred additive alias of pq_per_type_T0; optional so pre-alias extractions load through the legacy required name",
              irreps="0e", units="Angstrom^-4", mechanism="ring_efg"),
    ArraySpec("piquad_quad_scalar", "pi_quadrupole", np.ndarray, None, False, "Sparse per-(atom,ring) rows aligned to ring_contributions: real computed RingNeighbourhood.quad_scalar, not the derived geometry scalar in ring_contributions column 7",
              native_axis="ring_contribution_pair", irreps="0e", units="Angstrom^-4", mechanism="ring_efg"),
    ArraySpec("piquad_local_tensor", "pi_quadrupole", ShieldingTensor, 9, True, "Raw symmetric-traceless pi-quadrupole geometry tensor in the deterministic ring-local frame",
              native_axis="ring_contribution_pair", irreps="2e", units="Angstrom^-5", tensor_rank=2, mechanism="ring_efg",
              tensor_basis=_FULL9_BASIS, tensor_component_order=_FULL9_ORDER, tensor_frame="ring_local_vertex0_gauge", structural_zero_components="T0,T1_x,T1_y,T1_z", e3nn_export=_E3NN_EXPORT),
    ArraySpec("piquad_local_T2", "pi_quadrupole", np.ndarray, 5, True, "T2 convenience view of piquad_local_tensor in project-native real tesseral order",
              native_axis="ring_contribution_pair", irreps="2e", units="Angstrom^-5", tensor_rank=2, mechanism="ring_efg",
              tensor_basis=_T2_BASIS, tensor_component_order=_T2_ORDER,
              tensor_frame="ring_local_vertex0_gauge",
              structural_zero_components=_EFG_STRUCTURAL_ZEROS,
              e3nn_export=_E3NN_EXPORT),
    ArraySpec("piquad_local_frame", "pi_quadrupole", np.ndarray, 9, True, "Deterministic ring-local frame columns [x_axis,y_axis,z_axis], vertex-0 gauge; NaN when invalid",
              native_axis="ring_contribution_pair", mechanism="geometry"),
    ArraySpec("piquad_local_geometry", "pi_quadrupole", np.ndarray, 8, True, "Row-aligned local tensor provenance [atom, ring, type, distance_A, cos_theta, existing_axial_scalar, tensor_valid, aromatic_only]",
              native_axis="ring_contribution_pair", units="mixed_index_A_dimensionless", mechanism="ring_efg"),

    # ── Ring susceptibility (RingSusceptibilityResult.cpp) ────────
    ArraySpec("ringchi_scalar", "ring_susceptibility", np.ndarray, None, True, "Sparse computed ring-susceptibility geometry scalar aligned to ring_contributions",
              native_axis="ring_contribution_pair", irreps="0e", units="Angstrom^-3", mechanism="ring_current"),
    ArraySpec("ringchi_per_type_T0", "ring_susceptibility", PerRingTypeT0, 8, True, "Dense per-atom ring-susceptibility scalar grouped by aromatic RingTypeIndex",
              irreps="0e", units="Angstrom^-3", mechanism="ring_current"),

    # ── Dispersion (DispersionResult.cpp) ────────────────────────
    ArraySpec("disp_per_type_T0", "dispersion", PerRingTypeT0,     8,    True,  "Deprecated name for unit-coefficient aromatic R^-6 proximity per ring type; not D3/D4 energy",
              irreps="0e", units="Angstrom^-6", mechanism="ring_dispersion"),
    ArraySpec("aromatic_r6_proximity_per_type_T0", "dispersion", PerRingTypeT0, 8, False, "Canonical name for unit-coefficient aromatic R^-6 proximity per ring type; not D3/D4 energy",
              irreps="0e", units="Angstrom^-6", mechanism="ring_dispersion"),

    # ── McConnell (McConnellResult.cpp) ──────────────────────────
    # Forward schema: project-native packed SphericalTensor
    # [T0,T1_x,T1_y,T1_z,T2_m-2..+2] per source category and
    # source-strength channel. This is not an e3nn basis; source metadata
    # lives in extraction_manifest.json::feature_metadata.mcconnell.
    ArraySpec("mc_peptide_co_fixed",      "mcconnell", ShieldingTensor, 9, True, "McConnell peptide C=O fixed source response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_peptide_co_bo",         "mcconnell", ShieldingTensor, 9, True, "McConnell peptide C=O Wiberg BO response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_peptide_co_rhombic",    "mcconnell", ShieldingTensor, 9, True, "McConnell peptide C=O additive rhombic source delta",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_peptide_cn_fixed",      "mcconnell", ShieldingTensor, 9, True, "McConnell peptide C-N fixed source response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_peptide_cn_bo",         "mcconnell", ShieldingTensor, 9, True, "McConnell peptide C-N Wiberg BO response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_backbone_other_fixed",  "mcconnell", ShieldingTensor, 9, True, "McConnell strict-backbone other fixed source response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_backbone_other_bo",     "mcconnell", ShieldingTensor, 9, True, "McConnell strict-backbone other Wiberg BO response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_sidechain_co_fixed",    "mcconnell", ShieldingTensor, 9, True, "McConnell sidechain C=O fixed source response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_sidechain_co_bo",       "mcconnell", ShieldingTensor, 9, True, "McConnell sidechain C=O Wiberg BO response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_sidechain_other_fixed", "mcconnell", ShieldingTensor, 9, True, "McConnell sidechain other fixed source response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_sidechain_other_bo",    "mcconnell", ShieldingTensor, 9, True, "McConnell sidechain other Wiberg BO response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_disulfide_fixed",       "mcconnell", ShieldingTensor, 9, True, "McConnell disulfide fixed source response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_disulfide_bo",          "mcconnell", ShieldingTensor, 9, True, "McConnell disulfide Wiberg BO response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_aromatic_zeroed_fixed", "mcconnell", ShieldingTensor, 9, True, "McConnell aromatic fixed source response zeroed to avoid BS/HM double-count",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_aromatic_zeroed_bo",    "mcconnell", ShieldingTensor, 9, True, "McConnell aromatic Wiberg BO response zeroed to avoid BS/HM double-count",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_backbone_xh_fixed",      "mcconnell", ShieldingTensor, 9, True, "McConnell backbone C/N/O-H fixed source response, classified from typed heavy-atom semantics",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_backbone_xh_bo",         "mcconnell", ShieldingTensor, 9, True, "McConnell backbone C/N/O-H Wiberg BO response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_sidechain_xh_fixed",     "mcconnell", ShieldingTensor, 9, True, "McConnell sidechain C/N/O-H fixed source response, classified from typed heavy-atom semantics",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_sidechain_xh_bo",        "mcconnell", ShieldingTensor, 9, True, "McConnell sidechain C/N/O-H Wiberg BO response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_s_h_fixed",              "mcconnell", ShieldingTensor, 9, True, "McConnell S-H fixed source response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_s_h_bo",                 "mcconnell", ShieldingTensor, 9, True, "McConnell S-H Wiberg BO response",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_nearfield_counts",      "mcconnell", McConnellNearFieldCounts, 2, False, "McConnell near-field accepted/rejected source-target pair counts below 3 A", units="count", mechanism="bond_anisotropy"),
    ArraySpec("mc_nearest_co_dir",        "mcconnell", VectorField, 3, False, "Nearest accepted peptide C=O source direction per atom from ConformationAtom::dir_nearest_CO",
              irreps="1o", units="", tensor_rank=1, parity="odd", mechanism="bond_anisotropy"),
    ArraySpec("mc_nearest_co_midpoint",   "mcconnell", PositionField, 3, False, "Nearest accepted peptide C=O source midpoint per atom from ConformationAtom::nearest_CO_midpoint",
              units="Å", parity="odd", mechanism="bond_anisotropy"),
    ArraySpec("mc_nearest_co_T2",         "mcconnell", ShieldingTensor, 9, False, "Nearest accepted peptide C=O response per atom from ConformationAtom::T2_CO_nearest, packed [T0, T1x, T1y, T1z, T2_m-2..+2]",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_nearest_cn_T2",         "mcconnell", ShieldingTensor, 9, False, "Nearest accepted peptide C-N response per atom from ConformationAtom::T2_CN_nearest, packed [T0, T1x, T1y, T1z, T2_m-2..+2]",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),

    # Typed side-chain carbonyl inventory + the canonical McConnell
    # SidechainCO response projected onto atom rows.  The source tables use a
    # filtered axis because peptide carbonyls and other bond categories are
    # deliberately absent rather than represented as placeholder rows.
    ArraySpec("sidechain_co_source_bonds", "sidechain_carbonyl_anisotropy", np.ndarray, 8, True, "Typed side-chain C=O source rows [bond_index, carbon_atom, oxygen_atom, residue_index, BondCategory, PlanarGroupKind, oxygen_semantic_class, source_valid]",
              native_axis="sidechain_co_source", mechanism="bond_anisotropy"),
    ArraySpec("sidechain_co_frame", "sidechain_carbonyl_anisotropy", np.ndarray, 12, True, "Side-chain C=O local source frame [bond_midpoint_A, carbon_to_oxygen_x, in_plane_y, plane_normal_z]",
              native_axis="sidechain_co_source", units="mixed_A_dimensionless", mechanism="bond_anisotropy"),
    ArraySpec("sidechain_co_frame_quality", "sidechain_carbonyl_anisotropy", np.ndarray, 4, True, "Side-chain C=O frame audit [bond_length_A, orthogonality_error, raw_normal_norm_A2, frame_valid]",
              native_axis="sidechain_co_source", units="mixed", mechanism="bond_anisotropy"),
    ArraySpec("sidechain_co_fixed_T2", "sidechain_carbonyl_anisotropy", ShieldingTensor, 9, True, "Canonical fixed-strength McConnell SidechainCO response, full project-native SphericalTensor pack",
              native_axis="atom", irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("sidechain_co_bo_T2", "sidechain_carbonyl_anisotropy", ShieldingTensor, 9, True, "MOPAC bond-order-weighted McConnell SidechainCO response, full project-native SphericalTensor pack; NaN rows when MOPAC is absent",
              native_axis="atom", irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("sidechain_co_scalar_audit", "sidechain_carbonyl_anisotropy", np.ndarray, 4, True, "Per-atom audit [fixed_T2_norm, bo_T2_norm, accepted_sidechain_CO_source_count, nearest_sidechain_CO_midpoint_distance_A]",
              native_axis="atom", units="mixed", mechanism="bond_anisotropy"),

    # Legacy McConnell arrays retained as optional/deprecated wrappers for
    # reading old extraction directories; new C++ emits 20 category/channel
    # tensors plus the PeptideCO rhombic audit tensor above.
    ArraySpec("mc_shielding",     "mcconnell_legacy", ShieldingTensor,    9,    False,  "Legacy McConnell aggregate shielding", irreps=_SHIELD_IRREPS, units="Angstrom^-3",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_category_T2",   "mcconnell_legacy", PerBondCategoryT2,  25,   False,  "Legacy McConnell T2 per old bond category", irreps="2e", units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_scalars",       "mcconnell_legacy", McConnellScalars,   6,    False,  "Legacy McConnell scalar sums + distances", mechanism="bond_anisotropy"),

    # ── Coulomb (CoulombResult.cpp) — emitted alongside canonical APBS
    # in production. It supplies the vacuum field plus direct aliases of the
    # APBS reaction field, and also feeds the MOPAC-vs-FF14SB probe. ──
    ArraySpec("coulomb_efg",            "coulomb", ShieldingTensor, 9,   False, "Coulomb bare total EFG (full 9-pack; T0/T1 structural zeros)",
              irreps=_SHIELD_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),
    ArraySpec("coulomb_efg_t2",         "coulomb", EFGTensor,       5,   False, "Coulomb bare total EFG T2-only companion copied from coulomb_efg columns 4:9",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("coulomb_E",              "coulomb", VectorField,     3,   False, "Coulomb total E-field",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_E_backbone",     "coulomb", VectorField,     3,   False, "Coulomb E-field backbone",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_E_sidechain",    "coulomb", VectorField,     3,   False, "Coulomb E-field sidechain",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_E_aromatic",     "coulomb", VectorField,     3,   False, "Coulomb E-field aromatic",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_efg_backbone",   "coulomb", EFGTensor,       5,   False, "Coulomb EFG backbone (T2 only, symmetric-traceless)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("coulomb_efg_sidechain",  "coulomb", EFGTensor,       5,   False, "Coulomb EFG sidechain (T2 only, symmetric-traceless)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("coulomb_efg_aromatic",   "coulomb", EFGTensor,       5,   False, "Coulomb EFG aromatic (T2 only, symmetric-traceless)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("coulomb_scalars",        "coulomb", CoulombScalars,  4,   False, "Coulomb E-field scalars",
              mechanism="electrostatic_efg"),
    ArraySpec("coulomb_aromatic_E_proj", "coulomb", np.ndarray,     None, False, "Coulomb aromatic E-field parent-to-H projection; NaN for non-H or parentless atoms",
              irreps="0e", units="V/A", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_aromatic_n_src",  "coulomb", np.ndarray,     None, False, "Count of sidechain aromatic source atoms contributing to the Coulomb aromatic field (int32)", units="count", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_E_solvent",       "coulomb", VectorField,     3,   False, "APBS reaction-field alias: canonical APBS E = total PB minus homogeneous-vacuum reference",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_efg_solvent",     "coulomb", EFGTensor,       5,   False, "APBS reaction-field alias: canonical APBS EFG T2 = total PB minus homogeneous-vacuum reference",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("coulomb_shielding",      "coulomb_legacy", ShieldingTensor, 9, False, "Legacy name for Coulomb bare total EFG", irreps=_SHIELD_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),

    # ── H-Bond (HBondResult.cpp) ─────────────────────────────────
    ArraySpec("hbond_scalars",    "hbond", HBondScalars,           4,    True,  "H-bond scalars (nearest_dist, 1/r^3, count, McConnell angular scalar Σ)",
              mechanism="hbond_kernel"),
    ArraySpec("hbond_nearest_dir", "hbond", VectorField,            3,    False, "Nearest accepted H-bond direction per atom from the explicit donor H source point to the target atom",
              irreps="1o", tensor_rank=1, parity="odd", mechanism="hbond_kernel"),
    ArraySpec("hbond_flags", "hbond", np.ndarray,                   3,    False, "H-bond boolean flags as int8 columns: hbond_is_backbone, hbond_is_donor, hbond_is_acceptor",
              mechanism="hbond_kernel"),
    ArraySpec("hbond_pairs_index", "hbond", np.ndarray,             6,    False, "Accepted backbone H-bond rows [donor_residue, donor_N_atom, donor_H_atom, acceptor_residue, acceptor_O_atom, backbone_sequence_separation]",
              native_axis="hbond_pair", mechanism="hbond_kernel"),
    ArraySpec("hbond_pairs_geometry", "hbond", np.ndarray,          5,    False, "Row-aligned accepted backbone H-bond geometry [H_O_distance_A, N_H_O_angle_rad_at_H, H_to_O_unit_x, H_to_O_unit_y, H_to_O_unit_z]",
              native_axis="hbond_pair", units="mixed_A_radians_dimensionless", mechanism="hbond_kernel"),
    ArraySpec("hbond_pairs_angle_valid", "hbond", np.ndarray,    None,    False, "Row-aligned uint8 validity for the N-H...O angle",
              native_axis="hbond_pair", mechanism="hbond_kernel"),

    # ── DSSP (DsspResult.cpp) ────────────────────────────────────
    ArraySpec("dssp_observed",    "dssp", np.ndarray,              None, False, "DSSP observation mask per atom: int8 1 when the parent residue mapped to a libdssp row, 0 otherwise",
              mechanism="secondary_structure"),
    ArraySpec("dssp_backbone",    "dssp", DsspScalars,             5,    True,  "DSSP backbone geometry, columns [phi_IUPAC_rad, psi_IUPAC_rad, sasa_A2, ss_helix, ss_sheet]: phi/psi are IUPAC-signed radians (-libdssp), SASA is A^2, helix/sheet are dimensionless flags; phi/psi/SASA are NaN for unobserved residues",
              units="radians_for_phi_psi; Angstrom^2_for_sasa; dimensionless_for_flags", mechanism="secondary_structure"),
    ArraySpec("dssp_ss8",         "dssp", np.ndarray,              8,    False, "DSSP 8-class SS one-hot (H/G/I/E/B/T/S/C)",
              mechanism="secondary_structure"),
    ArraySpec("dssp_ppii",        "dssp", np.ndarray,              None, False, "DSSP PPII flag: int8 1=PPII, 0=observed non-PPII, -1=no DSSP observation",
              mechanism="secondary_structure"),
    ArraySpec("dssp_hbond_energy","dssp", np.ndarray,              4,    False, "DSSP H-bond energies (acc0/acc1/don0/don1)",
              units="kcal/mol", mechanism="secondary_structure"),
    ArraySpec("dssp_chi",         "dssp", np.ndarray,              12,   False, "Chi1-4 cos/sin/exists (4 x 3 cols)",
              mechanism="secondary_structure"),
    ArraySpec("dssp_torsion_angle", "dssp", np.ndarray,             6,   False, "Residue-axis signed IUPAC torsion angles [phi, psi, chi1, chi2, chi3, chi4]; undefined entries are NaN",
              native_axis="residue", units="radians", mechanism="secondary_structure"),
    ArraySpec("dssp_torsion_sin", "dssp", np.ndarray,               6,   False, "Sine of [phi, psi, chi1, chi2, chi3, chi4]; undefined entries are NaN",
              native_axis="residue", mechanism="secondary_structure"),
    ArraySpec("dssp_torsion_cos", "dssp", np.ndarray,               6,   False, "Cosine of [phi, psi, chi1, chi2, chi3, chi4]; undefined entries are NaN",
              native_axis="residue", mechanism="secondary_structure"),
    ArraySpec("dssp_torsion_valid", "dssp", np.ndarray,             6,   False, "Residue-axis uint8 validity for [phi, psi, chi1, chi2, chi3, chi4]",
              native_axis="residue", mechanism="secondary_structure"),
    ArraySpec("dssp_hbond_partner_residue_index", "dssp", np.ndarray, 4, False, "Atom-broadcast DSSP H-bond partner residue indices in energy-slot order [acceptor0, acceptor1, donor0, donor1]; -1 means no partner",
              native_axis="atom", mechanism="secondary_structure"),

    # ── SASA (SasaResult.cpp) ───────────────────────────────────
    ArraySpec("atom_sasa",        "sasa", np.ndarray,              None, False, "Per-atom Shrake-Rupley SASA (A^2)",
              units="Å^2", mechanism="solvation"),
    ArraySpec("sasa_normal",      "sasa", VectorField,             3,    False, "SASA outward surface normal (unit vector)",
              irreps="1o", tensor_rank=1, parity="odd", mechanism="solvation"),
    ArraySpec("atom_sasa_fraction", "sasa", np.ndarray,            None, False, "Normalized per-atom Shrake-Rupley exposure fraction, atom_sasa / (4*pi*(Bondi_vdW + probe)^2)",
              mechanism="solvation"),

    # ── Explicit water (WaterFieldResult.cpp) ───────────────────
    ArraySpec("water_efield",       "water_field", VectorField,    3,    False, "Water E-field total (V/A)",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="solvation"),
    ArraySpec("water_hbond_candidates", "water_hbond_geometry", np.ndarray, 16, False, "Raw protein-water H-bond candidates [protein atom/residue, water, mode, typed protein role, O/H/heavy distances, D-H...A angle, water O/H coordinates, pass flag]",
              native_axis="protein_water_hbond_pair", units="mixed_ids_A_degrees", mechanism="hbond_kernel"),
    ArraySpec("water_hbond_counts", "water_hbond_geometry", np.ndarray, 6, False, "Per-protein-atom water H-bond counts [water-donor candidates/pass, protein-donor candidates/pass, nearest water index, nearest passing mode]",
              mechanism="hbond_kernel"),
    ArraySpec("water_hbond_nearest", "water_hbond_geometry", np.ndarray, 8, False, "Per-protein-atom nearest candidate [water-O distance, H-acceptor distance, angle, mode, water index, pass, candidate count, pass count]",
              units="mixed_A_degrees", mechanism="hbond_kernel"),
    ArraySpec("water_efield_first", "water_field", VectorField,    3,    False, "Water E-field first shell <3.5A (V/A)",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="solvation"),
    ArraySpec("water_efg",          "water_field", EFGTensor,      5,    False, "Water EFG total (T2 only, symmetric-traceless)",
              units="V/A^2", tensor_rank=2, mechanism="solvation", **_T2_TENSOR_METADATA),
    ArraySpec("water_efg_first",    "water_field", EFGTensor,      5,    False, "Water EFG first shell (T2 only, symmetric-traceless)",
              units="V/A^2", tensor_rank=2, mechanism="solvation", **_T2_TENSOR_METADATA),
    ArraySpec("water_shell_counts", "water_field", np.ndarray,     2,    False, "Water shell counts [n_first, n_second]",
              mechanism="solvation"),
    ArraySpec("water_efield_clamp_mask", "water_field", np.ndarray, None, False, "Water total E-field clamp mask per atom: int8 0/1",
              mechanism="solvation"),
    ArraySpec("water_efield_clamp_scale", "water_field", np.ndarray, None, False, "Water total E-field clamp scale per atom",
              units="dimensionless", mechanism="solvation"),
    ArraySpec("water_efield_first_clamp_mask", "water_field", np.ndarray, None, False, "Water first-shell E-field clamp mask per atom: int8 0/1",
              mechanism="solvation"),
    ArraySpec("water_efield_first_clamp_scale", "water_field", np.ndarray, None, False, "Water first-shell E-field clamp scale per atom",
              units="dimensionless", mechanism="solvation"),

    # ── Hydration shell (HydrationShellResult.cpp) ──────────────
    ArraySpec("hydration_shell",    "hydration",   np.ndarray,     4,    False, "Hydration geometry [asymmetry, dipole_cos, ion_dist, ion_charge]",
              mechanism="solvation"),

    # ── Hydration geometry — SASA-normal (HydrationGeometryResult.cpp) ─
    ArraySpec("water_polarization", "water_polarization", np.ndarray, 10, False, "Water polarisation [dipole(3), normal(3), asym, align, mean_net_dipole_eA legacy coherence, count]",
              mechanism="solvation"),

    # ── Project-local QEq/EEQ-style charges (EeqResult.cpp) ──────
    # project-local QEq/EEQ-style charge-equilibration model with error-function coordination number, CN-dependent electronegativity shift, Gaussian self term, and Ohno-Klopman off-diagonal kernel; parameters are in-repo/project-local and are not a validated dftd4/multicharge port.
    ArraySpec("eeq_charges",        "eeq",         np.ndarray,     None, False, "Project-local QEq/EEQ-style geometry-dependent charges (e)",
              units="e", mechanism="charges"),
    ArraySpec("eeq_cn",             "eeq",         np.ndarray,     None, False, "EEQ coordination number",
              mechanism="charges"),
    ArraySpec("eeq_chi_eff",        "eeq",         np.ndarray,     None, False, "Project-local QEq/EEQ-style CN-shifted effective electronegativity in atomic units",
              units="Hartree", mechanism="charges"),
    ArraySpec("eeq_hardness",       "eeq",         np.ndarray,     2,    False, "Project-local QEq/EEQ-style hardness diagnostics [eta, self_hardness_diagonal] in atomic units",
              units="Hartree", mechanism="charges"),

    # EEQ-charge Coulomb fields (EeqCoulombResult.cpp). This is the FF
    # Coulomb shape family evaluated with geometry-dependent EEQ charges.
    ArraySpec("eeq_coulomb_efg",           "eeq_coulomb", ShieldingTensor, 9, False, "EEQ-charge Coulomb bare total EFG (full 9-pack; T0/T1 structural zeros)",
              irreps=_SHIELD_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),
    ArraySpec("eeq_coulomb_E",             "eeq_coulomb", VectorField, 3, False, "EEQ-charge Coulomb total E-field",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("eeq_coulomb_E_backbone",    "eeq_coulomb", VectorField, 3, False, "EEQ-charge Coulomb E-field backbone",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("eeq_coulomb_E_sidechain",   "eeq_coulomb", VectorField, 3, False, "EEQ-charge Coulomb E-field sidechain",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("eeq_coulomb_E_aromatic",    "eeq_coulomb", VectorField, 3, False, "EEQ-charge Coulomb E-field aromatic",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("eeq_coulomb_efg_backbone",  "eeq_coulomb", EFGTensor, 5, False, "EEQ-charge Coulomb EFG backbone (T2 only)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("eeq_coulomb_efg_sidechain", "eeq_coulomb", EFGTensor, 5, False, "EEQ-charge Coulomb EFG sidechain (T2 only)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("eeq_coulomb_efg_aromatic",  "eeq_coulomb", EFGTensor, 5, False, "EEQ-charge Coulomb EFG aromatic (T2 only)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("eeq_coulomb_scalars",       "eeq_coulomb", CoulombScalars, 4, False, "EEQ-charge Coulomb scalars [E_magnitude, H-only E_bond_proj, E_backbone_frac, aromatic_E_magnitude]",
              mechanism="electrostatic_efg"),
    ArraySpec("eeq_coulomb_aromatic_E_proj", "eeq_coulomb", np.ndarray, None, False, "EEQ-charge Coulomb aromatic E-field parent-to-H projection; NaN for non-H or parentless atoms",
              irreps="0e", units="V/A", mechanism="electrostatic_efg"),
    ArraySpec("eeq_coulomb_aromatic_n_src", "eeq_coulomb", np.ndarray, None, False, "Count of sidechain aromatic source atoms contributing to the EEQ-charge Coulomb aromatic field (int32)",
              units="count", mechanism="electrostatic_efg"),

    # ── GROMACS energy (GromacsEnergyResult.cpp) ────────────────
    ArraySpec("gromacs_energy",     "gromacs",     np.ndarray,     43,   False, "Per-frame energy (43 cols: electrostatic 3, bonded 6, VdW 3, thermo 8, box 3, virial 9, pressure tensor 9, T_group 2)",
              native_axis="protein", units="kJ/mol", mechanism="gromacs_runtime"),

    # ── Bonded energy (BondedEnergyResult.cpp) ─────────────────
    ArraySpec("bonded_energy",      "bonded",      np.ndarray,      7,   False, "Per-atom bonded energy (bond,angle,UB,proper,improper,CMAP,total) kJ/mol",
              units="kJ/mol", mechanism="gromacs_runtime"),

    # ── MOPAC core (MopacResult.cpp) ─────────────────────────────
    ArraySpec("mopac_charges",    "mopac_core", np.ndarray,        None, False, "MOPAC Mulliken charges",
              units="e", mechanism="charges"),
    ArraySpec("mopac_scalars",    "mopac_core", MopacScalars,      4,    False, "MOPAC per-atom scalars",
              mechanism="charges"),
    ArraySpec("mopac_bond_orders","mopac_core", BondOrders,        3,    False, "MOPAC Wiberg bond orders: sparse rows [atom_i, atom_j, order]",
              native_axis="bond", mechanism="charges"),
    ArraySpec("mopac_bond_neighbors", "mopac_core", np.ndarray,    4,    False, "Directed per-atom MOPAC bond-neighbor sparse rows [atom_i, atom_j, order, topology_bond_index]; rows are sorted descending by order within atom_i, topology_bond_index is -1 when absent",
              native_axis="mopac_bond_neighbor_pair", units="dimensionless", mechanism="charges"),
    ArraySpec("mopac_global",     "mopac_core", MopacGlobal,       4,    False, "MOPAC graph-level scalars",
              native_axis="protein", mechanism="charges"),
    ArraySpec("mopac_atom_populations", "mopac_core", MopacAtomPopulations, 12, False, "MOPAC per-atom charge, density, shell populations, dipole contribution, and valencies", native_axis="atom", units="mixed", mechanism="charges"),
    ArraySpec("mopac_atomic_orbital_populations", "mopac_core", MopacAtomicOrbitalPopulations, 9, False, "Frame-dependent printed MOPAC atomic-orbital electron populations s/px/py/pz/d; diagnostic table, not model-facing invariant scalars", native_axis="atom", units="electron", mechanism="charges"),
    ArraySpec("mopac_atomic_orbital_population_totals", "mopac_core", MopacAtomicOrbitalPopulationTotals, 3, False, "Invariant MOPAC AO shell totals [s_total,p_total,d_total] derived from printed atomic-orbital populations", native_axis="atom", units="electron", mechanism="charges"),
    ArraySpec("mopac_bond_valencies", "mopac_core", np.ndarray,   None, False, "MOPAC bond-order diagonal valencies, not recomputed", native_axis="atom", units="dimensionless", mechanism="charges"),
    ArraySpec("mopac_bond_orders_unique", "mopac_core", MopacUniqueBondOrders, 8, False, "Deterministic symmetric projection over printed MOPAC bond rows", native_axis="mopac_unique_pair", units="dimensionless", mechanism="charges"),
    ArraySpec("mopac_topology_bond_orders_full", "mopac_core", MopacTopologyBondOrdersFull, 8, False, "Topology-bond bridge with present flag and absence reason id", native_axis="bond", units="dimensionless", mechanism="charges"),

    # ── MOPAC Coulomb (MopacCoulombResult.cpp) ───────────────────
    ArraySpec("mopac_coulomb_efg",           "mopac_coulomb", ShieldingTensor, 9,  False, "MOPAC Coulomb bare total EFG (full 9-pack; T0/T1 structural zeros)",
              irreps=_SHIELD_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_E",             "mopac_coulomb", VectorField,     3,  False, "MOPAC Coulomb E-field",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_E_backbone",    "mopac_coulomb", VectorField,     3,  False, "MOPAC Coulomb E-field backbone",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_E_sidechain",   "mopac_coulomb", VectorField,     3,  False, "MOPAC Coulomb E-field sidechain",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_E_aromatic",    "mopac_coulomb", VectorField,     3,  False, "MOPAC Coulomb E-field aromatic",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_efg_backbone",  "mopac_coulomb", EFGTensor,       5,  False, "MOPAC Coulomb EFG backbone (T2 only)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("mopac_coulomb_efg_sidechain", "mopac_coulomb", EFGTensor,       5,  False, "MOPAC Coulomb EFG sidechain (T2 only)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("mopac_coulomb_efg_aromatic",  "mopac_coulomb", EFGTensor,       5,  False, "MOPAC Coulomb EFG aromatic (T2 only)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("mopac_coulomb_scalars",       "mopac_coulomb", CoulombScalars,  4,  False, "MOPAC Coulomb scalars",
              mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_shielding",     "mopac_coulomb_legacy", ShieldingTensor, 9, False, "Legacy name for MOPAC Coulomb bare total EFG", irreps=_SHIELD_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),

    # ── Legacy MOPAC McConnell (replaced by McConnell BO channel) ─
    ArraySpec("mopac_mc_shielding",    "mopac_mcconnell_legacy", ShieldingTensor,   9,  False, "Legacy MOPAC McConnell aggregate shielding", irreps=_SHIELD_IRREPS, units="Angstrom^-3",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mopac_mc_category_T2",  "mopac_mcconnell_legacy", PerBondCategoryT2, 25, False, "Legacy MOPAC McConnell T2 per category", irreps="2e", units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mopac_mc_scalars",      "mopac_mcconnell_legacy", McConnellScalars,  6,  False, "Legacy MOPAC McConnell scalars", mechanism="bond_anisotropy"),

    # ── APBS (ApbsFieldResult.cpp) ───────────────────────────────
    ArraySpec("apbs_E",           "apbs", VectorField,             3,    False, "APBS canonical reaction E-field: total PB minus homogeneous-vacuum reference",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("apbs_efg",         "apbs", EFGTensor,               5,    False, "APBS canonical reaction EFG: total PB minus homogeneous-vacuum reference (T2 only, symmetric-traceless)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("apbs_phi",         "apbs", np.ndarray,              None, False, "APBS canonical reaction potential: total PB minus homogeneous-vacuum reference",
              units="V", mechanism="electrostatic_efg"),
    ArraySpec("apbs_E_clamp_mask", "apbs", np.ndarray,             None, False, "APBS canonical reaction E-field clamp mask (uint8 0/1)",
              units="dimensionless", mechanism="electrostatic_efg"),
    ArraySpec("apbs_E_clamp_scale", "apbs", np.ndarray,            None, False, "Scale applied to canonical APBS reaction E; 1.0 when unclamped",
              units="dimensionless", mechanism="electrostatic_efg"),
    ArraySpec("apbs_nonfinite_sanitizer_mask", "apbs", np.ndarray, None, False, "APBS finite-sanitizer bit mask per atom: bit0 reaction E, bit1 reaction EFG, bit2 total E, bit3 total EFG",
              units="dimensionless", mechanism="electrostatic_efg"),
    ArraySpec("apbs_E_total_diagnostic", "apbs", VectorField,      3,    False, "APBS raw total-PB E-field diagnostic, finite-sanitized and unclamped",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("apbs_efg_total_diagnostic", "apbs", EFGTensor,      5,    False, "APBS raw total-PB EFG diagnostic T2, finite-sanitized and unclamped",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),

    # ── Orca DFT (OrcaShieldingResult.cpp) ───────────────────────
    ArraySpec("orca_total",       "orca", ShieldingTensor,         9,    False, "Orca DFT total shielding",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="quantum_reference"),
    ArraySpec("orca_diamagnetic", "orca", ShieldingTensor,         9,    False, "Orca DFT diamagnetic",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="quantum_reference"),
    ArraySpec("orca_paramagnetic","orca", ShieldingTensor,         9,    False, "Orca DFT paramagnetic",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="quantum_reference"),

    # ── Mutation delta (MutationDeltaResult.cpp) ─────────────────
    ArraySpec("delta_shielding",       "delta", ShieldingTensor,       9,    False, "WT-ALA shielding delta (total)",
              native_axis="mutation_match_pair", irreps=_SHIELD_IRREPS, units="ppm",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="mutation_delta"),
    ArraySpec("delta_scalars",         "delta", DeltaScalars,          6,    False, "Delta metadata + match info",
              native_axis="mutation_match_pair", mechanism="mutation_delta"),
    ArraySpec("delta_graph",           "delta", np.ndarray,             5,    False, "Graph deltas [matched, has_graph_delta, delta_graph_dist_ring, delta_bfs_decay, delta_is_conjugated]",
              native_axis="mutation_match_pair", mechanism="mutation_delta"),
    ArraySpec("delta_apbs",            "delta", DeltaAPBS,             12,   False, "APBS delta_E(3) + legacy full-9 EFG envelope; only columns 7:12 are physical EFG T2",
              native_axis="mutation_match_pair", mechanism="mutation_delta",
              tensor_basis=_T2_BASIS,
              tensor_component_order="delta_E_x,delta_E_y,delta_E_z,T0_compat_zero,T1_x_compat_zero,T1_y_compat_zero,T1_z_compat_zero,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2",
              tensor_frame="conformation_cartesian_xyz",
              structural_zero_components=_EFG_STRUCTURAL_ZEROS,
              e3nn_export="raw project tensor; call DeltaAPBS.delta_efg_t2.to_e3nn() before using e3nn Irreps"),
    ArraySpec("delta_ring_proximity",  "delta", DeltaRingProximity,    None, False, "Removed ring geometry (variable cols)",
              native_axis="mutation_match_pair", units="Å", mechanism="mutation_delta"),
    # DFT shielding component decomposition: WT side, mut side, deltas;
    # diamagnetic and paramagnetic. sigma_total = sigma_dia + sigma_para;
    # the existing delta_shielding satisfies that identity at ORCA's
    # output precision (~1e-3 ppm). Stratifies mutation shifts by
    # physical mechanism.
    ArraySpec("wt_shielding_diamagnetic",     "delta", ShieldingTensor, 9, False, "WT diamagnetic shielding (matched, by WT atom row)",
              native_axis="mutation_match_pair", irreps=_SHIELD_IRREPS, units="ppm",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="mutation_delta"),
    ArraySpec("wt_shielding_paramagnetic",    "delta", ShieldingTensor, 9, False, "WT paramagnetic shielding (matched, by WT atom row)",
              native_axis="mutation_match_pair", irreps=_SHIELD_IRREPS, units="ppm",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="mutation_delta"),
    ArraySpec("mut_shielding_diamagnetic",    "delta", ShieldingTensor, 9, False, "mut diamagnetic shielding (matched, by WT atom row)",
              native_axis="mutation_match_pair", irreps=_SHIELD_IRREPS, units="ppm",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="mutation_delta"),
    ArraySpec("mut_shielding_paramagnetic",   "delta", ShieldingTensor, 9, False, "mut paramagnetic shielding (matched, by WT atom row)",
              native_axis="mutation_match_pair", irreps=_SHIELD_IRREPS, units="ppm",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="mutation_delta"),
    ArraySpec("delta_shielding_diamagnetic",  "delta", ShieldingTensor, 9, False, "WT - mut diamagnetic shielding delta",
              native_axis="mutation_match_pair", irreps=_SHIELD_IRREPS, units="ppm",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="mutation_delta"),
    ArraySpec("delta_shielding_paramagnetic", "delta", ShieldingTensor, 9, False, "WT - mut paramagnetic shielding delta",
              native_axis="mutation_match_pair", irreps=_SHIELD_IRREPS, units="ppm",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="mutation_delta"),

    # ── Per-atom invariant categorical record (CategoryInfoProjection.cpp) ──
    # Structured-dtype NPY (~31 fields). One-shot per protein, NOT per-
    # conformation. Wrapper is np.ndarray here (delivered as the raw
    # structured array) and `load()` in _protein.py wraps it as
    # CategoryInfo — a circular import would result if the catalog
    # referenced the wrapper class directly.
    ArraySpec("atoms_category_info",   "identity", np.ndarray, None, False, "Per-atom invariant categorical record (structured dtype)",
              native_axis="atom", mechanism="topology"),

    # ── AIMNet2 (AIMNet2Result.cpp) ─────────────────────────────
    # Required in production output per project_aimnet2_contract_20260426.
    # The 2026-04-26 contract was articulated in the memory entry but
    # the required=True flag was not landed in the catalog at the time;
    # landed 2026-05-04 alongside AIMNet2 wire-in to smoke tests.
    ArraySpec("aimnet2_charges",             "aimnet2", AIMNet2Charges,            None, True,  "AIMNet2 Hirshfeld charges",
              units="e", mechanism="charges"),
    ArraySpec("aimnet2_aim",                 "aimnet2", AIMNet2AimEmbedding,       256,  True,  "AIMNet2 256-dim electronic embedding",
              mechanism="charges"),
    ArraySpec("aimnet2_efg",                 "aimnet2", EFGTensor,                 5,    True,  "AIMNet2 Coulomb EFG total (T2 only)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("aimnet2_efg_aromatic",        "aimnet2", EFGTensor,                 5,    True,  "AIMNet2 Coulomb EFG aromatic (T2 only)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("aimnet2_efg_backbone",        "aimnet2", EFGTensor,                 5,    True,  "AIMNet2 Coulomb EFG backbone (T2 only)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("aimnet2_efg_sidechain",       "aimnet2", EFGTensor,                 5,    True,  "AIMNet2 Coulomb EFG sidechain (T2 only)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("aimnet2_E",                   "aimnet2", VectorField,               3,    True,  "AIMNet2 charge-derived total E-field",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("aimnet2_E_backbone",          "aimnet2", VectorField,               3,    True,  "AIMNet2 charge-derived E-field backbone",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("aimnet2_E_sidechain",         "aimnet2", VectorField,               3,    True,  "AIMNet2 charge-derived E-field sidechain",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("aimnet2_E_aromatic",          "aimnet2", VectorField,               3,    True,  "AIMNet2 charge-derived E-field aromatic",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("aimnet2_energy_mlp",          "aimnet2", np.ndarray,                None, True, "AIMNet2 per-atom energy after the energy_mlp head",
              units="eV", mechanism="charges"),
    ArraySpec("aimnet2_energy_shifted_local", "aimnet2", np.ndarray,               None, True, "AIMNet2 per-atom shifted local energy after the atomic_shift head",
              units="eV", mechanism="charges"),
    ArraySpec("aimnet2_energy_terms",        "aimnet2", np.ndarray,                6,    True, "AIMNet2 protein-level energy terms [local_sum, e_lrcoulomb, e_dftd3, total, conditioned_net_charge, neutral_conditioning_flag]",
              native_axis="protein", units="mixed:eV[0:4],e[4],dimensionless[5]", mechanism="charges"),
    ArraySpec("aimnet2_d3_e_disp_atom",      "aimnet2", np.ndarray,                None, True, "AIMNet2 D3 per-atom dispersion-energy increment",
              units="eV", mechanism="charges"),
    ArraySpec("aimnet2_d3_cn",               "aimnet2", np.ndarray,                None, True, "AIMNet2 D3 coordination number",
              units="dimensionless", mechanism="charges"),
    ArraySpec("aimnet2_d3_c6_stats",         "aimnet2", np.ndarray,                3,    True, "AIMNet2 D3 C6 statistics [sum, mean, max] over valid long-range neighbours",
              units="Hartree*bohr^6", mechanism="charges"),
    ArraySpec("aimnet2_aim_projection",      "aimnet2", np.ndarray,                32,   True, "AIMNet2 fixed 32-d projection of raw aim; basis splitmix64_0xA17E20260708_achlioptas_32x256_element_HCNOS",
              units="dimensionless", mechanism="charges"),

    # ── AIMNet2 charge-response gradient (AIMNet2ChargeResponseGradientResult.cpp)
    # Always-on after the --aimnet2 model is loaded (per the
    # 2026-05-09 promotion of Amendment 2026-05-08(b) from a test flag
    # to standard non-trajectory pipeline; trajectory mode unchanged).
    # Vector is dL/d(r_i) where L = sum_j q_j^2 over non-sentinel atoms
    # (charge-conservation makes sum(q) gradient near-zero, so the L2
    # of charges is the cheapest single-pass objective with non-trivial
    # gradient). Scalar is the L2 norm of the vector. required=False
    # because old extraction outputs (pre-2026-05-09) do not include
    # them; new outputs always do when AIMNet2 is loaded. Renamed from
    # "polarisability" 2026-05-20 (commit 58594f5) — the emission is
    # ∂(Σq²)/∂r, NOT a Buckingham α = ∂μ/∂E.
    ArraySpec("aimnet2_charge_response_gradient",        "aimnet2", AIMNet2ChargeResponseGradient, 3,    False, "AIMNet2 charge-response proxy/diagnostic: d(sum q_j^2)/dR_i; NOT Buckingham polarizability",
              irreps="1o", units="e^2/Å", tensor_rank=1, parity="odd", mechanism="charges"),
    ArraySpec("aimnet2_charge_response_gradient_scalar", "aimnet2", np.ndarray,            None, False, "AIMNet2 charge-response proxy/diagnostic scalar: L2 norm of d(sum q_j^2)/dR_i; NOT Buckingham polarizability",
              units="e^2/Å", mechanism="charges"),

    # ── Planar geometry (PlanarGeometryResult.cpp) ───────────────────
    # Per Amendment 2026-05-08(a). Conformation-only quantities derived
    # from positions; runs whenever the substrate (LegacyAmber
    # AtomSemanticTable) is populated. Twelve NPYs with three different
    # axes (per-atom, per-residue, per-ring) — the reader's wrapper
    # picks shape from the catalog cols field.
    #
    # Conventions: pyramidalization is a non-negative magnitude with
    # explicit valid/type arrays; omega in radians; Cremer-Pople θ in
    # degrees [0, 360).
    # required=False because pre-2026-05-09 extractions don't carry them.
    ArraySpec("pyramidalization",  "planar_geometry", np.ndarray, None, False, "Per-atom non-negative sp2 out-of-plane displacement magnitude (Å); NaN for non-applicable or invalid/degenerate rows; see pyramidalization_valid and pyramidalization_center_type",
              units="Å", mechanism="geometry"),
    ArraySpec("pyramidalization_valid", "planar_geometry", np.ndarray, None, False, "Per-atom int8 mask: 1 when pyramidalization is a finite computed scalar, 0 for non-applicable or invalid/degenerate rows",
              native_axis="atom", mechanism="geometry"),
    ArraySpec("pyramidalization_center_type", "planar_geometry", np.ndarray, None, False, "Per-atom int8 PlanarGroupKind code: 0 None, 1 PeptideAmide, 2 SidechainAmide, 3 Guanidinium, 4 Imidazole, 5 Aromatic6Ring, 6 Aromatic5Ring, 7 Carboxylate, 8 AromaticHydroxyl, 9 AromaticOxide",
              native_axis="atom", mechanism="geometry"),
    ArraySpec("omega_actual",      "planar_geometry", np.ndarray, None, False, "Per-residue ω (Cα-C-N-Cα to next), radians; NaN only when the covalent backbone successor omega is undefined",
              native_axis="residue", units="radians", mechanism="geometry"),
    ArraySpec("omega_deviation",   "planar_geometry", np.ndarray, None, False, "Per-residue WrapPi(ω - π) in [-π, π]; emitted for every well-defined covalent backbone successor including X→Pro; use omega_is_xpro to analyse X→Pro separately",
              native_axis="residue", units="radians", mechanism="geometry"),
    ArraySpec("aromatic_chi2",     "planar_geometry", np.ndarray, None, False, "Per-aromatic-ring χ₂ (parent residue, radians); ring-flip observable per Akke-Weininger 2023",
              native_axis="aromatic_ring", units="radians", mechanism="geometry"),
    ArraySpec("pucker_Q",          "planar_geometry", np.ndarray, None, False, "Per-saturated-ring Cremer-Pople puckering amplitude (Å); 5-rings only",
              native_axis="saturated_ring", units="Å", mechanism="geometry"),
    ArraySpec("pucker_theta",      "planar_geometry", np.ndarray, None, False, "Per-saturated-ring Cremer-Pople phase angle (degrees, [0, 360))",
              native_axis="saturated_ring", units="degrees", mechanism="geometry"),
    ArraySpec("omega_is_xpro",     "planar_geometry", np.ndarray, None, False, "Per-residue mask: 1 where the covalent backbone successor is Pro; X→Pro cis/trans is real signal and is not NaN-filled",
              native_axis="residue", mechanism="geometry"),
    ArraySpec("omega_sin",         "planar_geometry", np.ndarray, None, False, "Sine of the signed IUPAC peptide omega angle; NaN when undefined",
              native_axis="residue", mechanism="geometry"),
    ArraySpec("omega_cos",         "planar_geometry", np.ndarray, None, False, "Cosine of the signed IUPAC peptide omega angle; NaN when undefined",
              native_axis="residue", mechanism="geometry"),
    ArraySpec("omega_valid",       "planar_geometry", np.ndarray, None, False, "Per-residue uint8 validity for omega_actual/omega_sin/omega_cos",
              native_axis="residue", mechanism="geometry"),

    # ── Tripeptide DFT shielding ────────────────────────────────────
    # ProCS15 (Larsen 2015) tripeptide DFT lookup. σ_BB^i emitted on
    # backbone N/CA/C/O/H/HA and central-residue sidechain atoms per
    # the typed-identity-matched LarsenResidue model. Δσ_BB^{i±1}
    # neighbor correction emitted at the central residue's atoms per
    # Larsen Eq 3 cap-side reading. required=False because the
    # tensorcs15 DB is not available on every host.
    ArraySpec("tripeptide_bb_shielding",          "tripeptide", ShieldingTensor, 9,    False, "σ_BB^i — Mat3 (ppm) from typed-identity match against Larsen 2015 AXA tripeptide DFT row",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="quantum_reference"),
    ArraySpec("tripeptide_bb_residual_vec",       "tripeptide", VectorField,     3,    False, "σ_BB^i match residual: aligned_dft - protein position; Vec3 ML feature (magnitude + direction)",
              irreps="1o", units="Å", tensor_rank=1, parity="odd", mechanism="quantum_reference"),
    ArraySpec("tripeptide_bb_match_distance",     "tripeptide", np.ndarray,      None, False, "σ_BB^i match distance (Å) — magnitude of residual_vec",
              units="Å", mechanism="quantum_reference"),
    ArraySpec("tripeptide_bb_method_tag",         "tripeptide", np.ndarray,      None, False, "DFT method discriminator: 0=none, 1=OPBE Gaussian (Larsen), 2=PBE ORCA (project SER regen)",
              mechanism="quantum_reference"),
    ArraySpec("tripeptide_bb_match_atoms",        "tripeptide", np.ndarray,      5,    False, "σ_BB^i atom-match metadata columns [residue_index, has_match, matched_dft_atom_idx, match_distance_A, method_tag]",
              native_axis="atom", mechanism="quantum_reference"),
    ArraySpec("tripeptide_neighbor_shielding",    "tripeptide", ShieldingTensor, 9,    False, "Δσ_BB^{i±1} — neighbour correction at residue i from i±1 cap reads (Larsen 2015 Eq 3)",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="quantum_reference"),
    ArraySpec("tripeptide_neighbor_shielding_prev", "tripeptide", ShieldingTensor, 9,  False, "Δσ_BB^{i-1} prev-direction neighbour correction tensor before prev+next summation",
              native_axis="atom", irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="quantum_reference"),
    ArraySpec("tripeptide_neighbor_shielding_next", "tripeptide", ShieldingTensor, 9,  False, "Δσ_BB^{i+1} next-direction neighbour correction tensor before prev+next summation",
              native_axis="atom", irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="quantum_reference"),
    ArraySpec("tripeptide_neighbor_residual_vec_prev", "tripeptide", VectorField, 3,   False, "Δσ_BB^{i-1} match residual at the C-term ALA cap of (i-1)'s tripeptide; Vec3, NaN where i-1 direction had no contribution",
              irreps="1o", units="Å", tensor_rank=1, parity="odd", mechanism="quantum_reference"),
    ArraySpec("tripeptide_neighbor_residual_vec_next", "tripeptide", VectorField, 3,   False, "Δσ_BB^{i+1} match residual at the N-term ALA cap of (i+1)'s tripeptide; Vec3, NaN where i+1 direction had no contribution",
              irreps="1o", units="Å", tensor_rank=1, parity="odd", mechanism="quantum_reference"),
    ArraySpec("tripeptide_neighbor_reference",    "tripeptide", np.ndarray,      5,    False, "AAA reference metadata columns [aaa_calc_id, aaa_frame_type_code, aaa_phi_db_deg, aaa_psi_db_deg, any_mixed_method_flag]",
              native_axis="protein", mechanism="quantum_reference"),
    ArraySpec("tripeptide_bb_diagnostics", "tripeptide", np.ndarray, 28, False, "Per-residue central tripeptide lookup diagnostics including exact status (0 miss, 1 ok, 2 unsupported HID, 3 unsupported HIE, 4 perception failed), natural/used chi depth, dropped-chi distance, actual/target/database angles, residue and HIS variant",
              native_axis="residue", units="mixed_ids_A_degrees", mechanism="quantum_reference"),
    ArraySpec("tripeptide_neighbor_diagnostics", "tripeptide", np.ndarray, 59, False, "Per-central-residue previous/next tripeptide lookup diagnostics (28 columns each) plus AAA frame, mixed-method, and any-match flags; status codes 0 miss, 1 ok, 2 unsupported HID, 3 unsupported HIE, 4 perception failed",
              native_axis="residue", units="mixed_ids_A_degrees", mechanism="quantum_reference"),

    # ────────────────────────────────────────────────────────────────
    # Larsen H-bond shielding contributions
    # (src/LarsenHBondShieldingResult.cpp). Direct DFT-grid lookup
    # against Larsen 2015 ProCS15 6-archive scan (NMA|ALA donor ×
    # NMA|HOMe|Acetate acceptor). Donor enumeration is spatial:
    # backbone amide H (HN) and any α-hydrogen (HA, plus GLY HA2/HA3)
    # within 4.2 Å of any acceptor O classified as one of
    # BackboneCarbonyl / SidechainCarbonyl / Hydroxyl / Carboxylate.
    # No DSSP — Larsen's framework is geometric (θ ≥ 90° is the gate),
    # and the spatial sweep IS the H-bond finder.
    #
    # Per-class Mat3 contributions are decomposed per Pattern 11 into
    # SphericalTensor (T0 + T1 + T2 = 9 columns) and accumulated to
    # ConformationAtom. Methods accumulate side-by-side with the
    # scalar-geometry HBondResult — see feedback_methods_accumulate memory
    # entry.
    #
    # STRUCTURE NOTE (per audit 2026-05-13):
    #   `larsen_hbond_shielding` is the ELEMENT-WISE SUM of the four
    #   per-class tensors, i.e. an exact linear combination of the
    #   per-class breakdown. That is a property for the analysis to
    #   handle at modeling time (exact aggregates are grouped via the
    #   relaimpo family / redundancy grouping), not a reason to pre-drop
    #   any array here. The per-class breakdown carries distinct physics
    #   per Larsen Table 1.
    #
    #   `larsen_hbond_diagnostic_CB_shielding` is a parser-pipeline
    #   reality check (Larsen Table 2 says Cβ gets NO H-bond term; we
    #   emit anyway to verify the rotation pipeline produces what the
    #   physics predicts).
    #
    #   `larsen_hbond_count` is per-atom metadata (number of pairs the
    #   atom received contributions from), correlated with shielding
    #   magnitude.
    #
    #   `larsen_hbond_water_term` is the Larsen Δσ_w = 2.07 ppm isotropic
    #   offset on solvent-exposed amide Hs.
    # ────────────────────────────────────────────────────────────────
    ArraySpec("larsen_hbond_shielding",                  "larsen_hbond", ShieldingTensor, 9, False, "Σ Larsen H-bond contributions across all four Table 2 classes (1°HB + 2°HB + 1°HαB + 2°HαB) — ppm, lab frame. Structurally = sum of the four per-class arrays; NOT a feature.", irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_1pHB_shielding",             "larsen_hbond", ShieldingTensor, 9, False, "Δσ_1°HB per Larsen 2015 Table 2 — primary amide-H donor contribution applied to donor residue i atoms",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_2pHB_shielding",             "larsen_hbond", ShieldingTensor, 9, False, "Δσ_2°HB per Larsen 2015 Table 2 — secondary amide-H donor contribution. N/Cα/Hα/HN apply to acceptor residue j+1; C' applies to acceptor's OWN residue j",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_1pHaB_shielding",            "larsen_hbond", ShieldingTensor, 9, False, "Δσ_1°HαB per Larsen 2015 Table 2 — primary Hα donor contribution applied to donor residue i atoms",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_2pHaB_shielding",            "larsen_hbond", ShieldingTensor, 9, False, "Δσ_2°HαB per Larsen 2015 Table 2 — secondary Hα donor contribution. N/Cα/Hα/HN apply to acceptor residue j+1; C' applies to acceptor's OWN residue j",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_diagnostic_CB_shielding",    "larsen_hbond", ShieldingTensor, 9, False, "Cβ diagnostic — Larsen Table 2 says Cβ gets NO contribution; emitted as parser→loader→frame-rotation reality check. NOT a feature.", irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_water_term",                 "larsen_hbond", np.ndarray,      None, False, "Δσ_w = 2.07 ppm isotropic on amide H atoms with NO valid, symmetry-selected geometric H-bond candidate (θ ≥ 90° in 4.2 Å); finite valid grid misses still confirm an H-bond, while invalid frames and filtered carboxylate siblings do not", units="ppm", mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_count",                      "larsen_hbond", np.ndarray,      None, False, "Per-atom count of H-bond pairs that contributed under any of the four Table 2 classes; metadata, NOT a feature.", mechanism="hbond_grid"),
    ArraySpec("larsen_corner_imputed",                   "larsen_hbond", np.ndarray,      None, False, "Per-atom int8 flag: 1 iff any Larsen H-bond grid lookup corner serving this atom was imputed", mechanism="hbond_grid"),
    ArraySpec("larsen_imputed_pair_count",               "larsen_hbond", np.ndarray,      None, False, "Per-atom count of Table2 shielding pairs with at least one imputed interpolation corner; imputed_policy=emitted_unmasked", mechanism="hbond_grid"),
    ArraySpec("larsen_sidechain_carbonyl_pair_count",    "larsen_hbond", np.ndarray,      None, False, "Per-atom count of contributing SidechainCarbonyl acceptor pairs; SidechainCarbonyl uses BackboneCarbonyl/NMA archive approximation", mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_pairs_index",                "larsen_hbond", np.ndarray,      16, False, "Raw per-candidate Larsen integer provenance [donor_atom, acceptor_atom, donor_residue, acceptor_residue, donor_class, acceptor_class, disposition, donor_anchor_atom, donor_third_atom, acceptor_C_atom, acceptor_third_atom, target_mask_1pHB, target_mask_2pHB, target_mask_1pHaB, target_mask_2pHaB, target_mask_diagnostic_CB]; disposition 0=missing_frame_atoms, 1=theta_out_of_range, 2=grid_miss, 3=success, 4=invalid_frame, 5=carboxylate_symmetry_filtered",
              native_axis="larsen_hbond_pair", mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_pairs_geometry",             "larsen_hbond", np.ndarray,       6, False, "Per-candidate geometry [r_A, theta_deg, rho_deg, any_corner_imputed, imputed_corner_count, frame_valid]; frame_valid requires finite query geometry and a non-degenerate donor rotation frame; imputed_policy=emitted_unmasked",
              native_axis="larsen_hbond_pair", units="mixed_A_degrees", mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_pairs_isotropic",            "larsen_hbond", np.ndarray,       6, False, "Per-candidate isotropic terms [1pHB,2pHB,1pHaB,2pHaB,diagnostic_CB,Table2_total] ppm",
              native_axis="larsen_hbond_pair", units="ppm", mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_pairs",                      "larsen_hbond", np.ndarray,      28, False, "Compatibility Larsen pair table: 16 index + 6 geometry + 6 isotropic columns",
              native_axis="larsen_hbond_pair", units="mixed", mechanism="hbond_grid"),
    ArraySpec("larsen_sidechain_donor_atoms", "larsen_sidechain_donor_audit", np.ndarray, 6, False, "Per-atom audit [is_sidechain_polar_H, PolarHKind, parent, residue, acceptor_candidates_3p5A, geometry_pass_count]",
              mechanism="hbond_grid"),
    ArraySpec("larsen_sidechain_donor_candidates", "larsen_sidechain_donor_audit", np.ndarray, 13, False, "Raw sidechain polar-H donor candidate geometry; audit-only, modeled_by_larsen_table2 is always false",
              native_axis="larsen_sidechain_donor_pair", units="mixed_ids_A_degrees", mechanism="hbond_grid"),

    # ────────────────────────────────────────────────────────────────
    # Topology sidecar (TopologySidecar.cpp) — additive 2026-05-13.
    # Three structured-NPY projections of LegacyAmberTopology +
    # RingTopology. Emitted alongside atoms_category_info.npy whenever
    # the protein has a populated typed substrate.
    # ────────────────────────────────────────────────────────────────
    ArraySpec("residues",         "topology", np.ndarray, None, True,  "Per-residue record: residue_index, chain_id, residue_number, insertion_code, residue_type, amber/iupac/one-letter names, protonation_variant_index, terminal_state, prev/next links, atom_count, is_proline/aromatic/titratable, has_amide_h",
              native_axis="residue", mechanism="topology"),
    ArraySpec("bonds",            "topology", np.ndarray, None, True,  "Per-bond record: bond_index, atom_index_a/b, bond_order, bond_category, is_rotatable, is_aromatic, is_peptide, is_backbone",
              native_axis="bond", mechanism="topology"),
    ArraySpec("rings",            "topology", np.ndarray, None, True,  "Per-ring record: ring_id, ring_kind (aromatic|saturated), ring_type_index, atom_count, native_axis_index, parent_residue_index, parent_residue_number, fused_partner_ring_id",
              native_axis="ring", mechanism="topology"),
    ArraySpec("ring_membership",  "topology", np.ndarray, None, True,  "Per (ring, ring-vertex-atom) row: ring_id, atom_index, ring_atom_order, is_vertex, is_substituent",
              native_axis="ring_membership", mechanism="topology"),
]}
# fmt: on


def _with_project_tensor_metadata(spec: ArraySpec) -> ArraySpec:
    if spec.wrapper is ShieldingTensor:
        structural_zeros = spec.structural_zero_components
        if spec.stem in {
            "coulomb_efg",
            "eeq_coulomb_efg",
            "mopac_coulomb_efg",
            "mopac_coulomb_shielding",
        }:
            structural_zeros = _EFG_STRUCTURAL_ZEROS
        return replace(
            spec,
            tensor_basis=spec.tensor_basis or _FULL9_BASIS,
            tensor_component_order=spec.tensor_component_order or _FULL9_ORDER,
            tensor_frame=spec.tensor_frame or "conformation_cartesian_xyz",
            structural_zero_components=structural_zeros,
            e3nn_export=spec.e3nn_export or _E3NN_EXPORT,
        )
    if spec.wrapper in (PerRingTypeT2, PerBondCategoryT2):
        return replace(
            spec,
            tensor_basis=spec.tensor_basis or _T2_BASIS,
            tensor_component_order=spec.tensor_component_order or _T2_ORDER,
            tensor_frame=spec.tensor_frame or "conformation_cartesian_xyz",
            e3nn_export=spec.e3nn_export or _E3NN_EXPORT,
        )
    return spec


CATALOG = {
    stem: _with_project_tensor_metadata(spec)
    for stem, spec in CATALOG.items()
}


# Directional-output freeze contract.  Keep this additive table separate from
# the historical declarations above: wrappers, legacy irreps, and loader
# compatibility stay untouched while the serialized physical law is explicit.
_CARTESIAN_FRAME = "conformation_cartesian_xyz"
_INTRINSIC_FRAME = "intrinsic_geometry"
_POLAR_VECTOR = "polar_vector: v'=R v"
_AXIAL_VECTOR = "axial_vector: a'=det(R) R a"
_AFFINE_POSITION = "cartesian_position: p'=R p+t"
_EVEN_RANK2 = "even_rank2: T'=R T R^T"
_NATIVE_T2 = (
    "even_rank2_native_T2: reconstruct Cartesian T, apply T'=R T R^T, "
    "then decompose in project-native T2 basis"
)


def _set_contract(stems, *, coordinate_frame, transformation, validity,
                  **metadata):
    """Apply one exact producer contract to existing catalog entries."""
    for stem in stems:
        if stem not in CATALOG:
            raise RuntimeError(
                f"directional contract references unknown NPY stem {stem!r}")
        CATALOG[stem] = replace(
            CATALOG[stem],
            coordinate_frame=coordinate_frame,
            transformation=transformation,
            validity=validity,
            **metadata,
        )


# Writers serialize every MutationDeltaResult table on the WT atom axis, not
# on a compact matched-pair axis.  MOPAC AO rows are likewise their own parsed
# row axis: the writer does not guarantee K == atom_count.
for _stem in (
    "delta_shielding", "delta_scalars", "delta_graph", "delta_apbs",
    "delta_ring_proximity", "wt_shielding_diamagnetic",
    "wt_shielding_paramagnetic", "mut_shielding_diamagnetic",
    "mut_shielding_paramagnetic", "delta_shielding_diamagnetic",
    "delta_shielding_paramagnetic",
):
    CATALOG[_stem] = replace(CATALOG[_stem], native_axis="atom")
for _stem in (
    "mopac_atomic_orbital_populations",
    "mopac_atomic_orbital_population_totals",
):
    CATALOG[_stem] = replace(
        CATALOG[_stem], native_axis="mopac_atomic_orbital_row")


# Cartesian positions and homogeneous vectors.
_set_contract(
    ("pos",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=_AFFINE_POSITION,
    validity="one row per atom; no separate position-validity mask")
_set_contract(
    ("bond_direction",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=_POLAR_VECTOR,
    validity="bond_geometry_valid.npy; invalid/zero-length bonds serialize zero")
_set_contract(
    ("bond_geometry_valid",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant uint8 validity for the row-aligned bond "
        "direction; translation invariant"
    ),
    validity=(
        "1 iff the bond displacement is finite and nonzero; 0 gates the "
        "legacy zero in bond_direction.npy"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("ring_direction_to_center",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=_POLAR_VECTOR,
    validity=(
        "row-aligned with ring_contributions.npy; center-to-atom unit vector "
        "despite the legacy stem; coincident center serializes zero without a mask"
    ))
_set_contract(
    ("mc_nearest_co_dir",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=_POLAR_VECTOR,
    validity=(
        "NaN triplet when nearest_CO_dist is NO_DATA_SENTINEL; finite only for "
        "an accepted peptide C=O source"
    ))
_set_contract(
    ("mc_nearest_co_midpoint",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=_AFFINE_POSITION,
    validity=(
        "NaN triplet when nearest_CO_dist is NO_DATA_SENTINEL; finite only for "
        "an accepted peptide C=O source"
    ))
_set_contract(
    ("hbond_nearest_dir",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=_POLAR_VECTOR,
    validity="hbond_flags.npy column 0 identifies an accepted nearest source; otherwise zero")
_set_contract(
    ("hbond_flags",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant int8 atom classifications for fixed typed "
        "topology and accepted-pair set; translation invariant"
    ),
    validity=(
        "columns are is_backbone,is_donor,is_acceptor; column 0 gates the "
        "legacy zero in hbond_nearest_dir.npy"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("hbond_pairs_index",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant sparse donor/acceptor residue and atom identity "
        "for a fixed accepted-pair set; deterministic row order"
    ),
    validity=(
        "row identity for hbond_pairs_geometry.npy and "
        "hbond_pairs_angle_valid.npy; zero rows means no accepted pair"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("hbond_pairs_angle_valid",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant uint8 angle-availability mask, row-aligned with "
        "hbond_pairs_index.npy"
    ),
    validity=(
        "1 iff the N-H...O angle in hbond_pairs_geometry.npy is valid; "
        "distance and direction columns have their own producer semantics"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("cb_deviation_valid", "cb_residual_vector_valid"),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant uint8 availability for the corresponding "
        "chiral-L-CB output; translation invariant"
    ),
    validity=(
        "1 iff the corresponding cb_deviation.npy scalar or all three "
        "cb_residual_vector.npy components are finite; 0 gates their NaNs"
    ),
    irreps="0e", parity="even", tensor_rank=0)

_set_contract(
    ("coulomb_E", "coulomb_E_backbone", "coulomb_E_sidechain",
     "coulomb_E_aromatic"),
    coordinate_frame=_CARTESIAN_FRAME, transformation=_POLAR_VECTOR,
    validity=(
        "force-field Coulomb producer replaces non-finite values with zero; "
        "no sanitizer mask (legacy unavailable/zero ambiguity); all four "
        "vectors are jointly scaled when total-E exceeds the configured clamp, "
        "with provenance only in GeometryChoice"
    ))
_set_contract(
    ("coulomb_E_solvent",),
    coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "continuum polar_vector APBS alias: v'=R v; translation invariant. "
        "The live axis-aligned finite-difference APBS source has no exact "
        "O(3) law; transformed production reruns use the recorded 1.8e-2 "
        "V/A absolute + 5e-2 relative finite-grid envelope"
    ),
    validity=(
        "optional clamped APBS reaction-field alias; absent without APBS; "
        "consult apbs_nonfinite_sanitizer_mask.npy plus "
        "apbs_E_clamp_mask.npy/apbs_E_clamp_scale.npy when present"
    ))
_set_contract(
    ("mopac_coulomb_E", "mopac_coulomb_E_backbone",
     "mopac_coulomb_E_sidechain", "mopac_coulomb_E_aromatic"),
    coordinate_frame=_CARTESIAN_FRAME, transformation=_POLAR_VECTOR,
    validity=(
        "MOPAC-charge Coulomb producer replaces non-finite values with zero; "
        "no sanitizer mask (legacy unavailable/zero ambiguity); all four "
        "vectors are jointly scaled when total-E exceeds the configured clamp, "
        "with provenance only in GeometryChoice"
    ))
_set_contract(
    ("eeq_coulomb_E", "eeq_coulomb_E_backbone",
     "eeq_coulomb_E_sidechain", "eeq_coulomb_E_aromatic"),
    coordinate_frame=_CARTESIAN_FRAME, transformation=_POLAR_VECTOR,
    validity=(
        "whole EEQ Coulomb result is absent if a non-finite field is produced; "
        "all four vectors are jointly scaled when total-E exceeds the "
        "configured clamp, with provenance only in GeometryChoice"
    ))
_set_contract(
    ("apbs_E", "apbs_E_total_diagnostic"),
    coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "continuum polar_vector: v'=R v; translation invariant. The live "
        "axis-aligned finite-difference APBS solve has no exact O(3) law; "
        "transformed production reruns use the recorded 1.8e-2 V/A absolute "
        "+ 5e-2 relative finite-grid envelope"
    ),
    validity=(
        "apbs_nonfinite_sanitizer_mask.npy bits 0/2 respectively; canonical "
        "reaction field also has apbs_E_clamp_mask.npy and apbs_E_clamp_scale.npy"
    ))
_set_contract(
    ("apbs_phi",),
    coordinate_frame="lab_fixed_apbs_finite_difference_grid",
    transformation=(
        "continuum rotation-invariant electrostatic potential; the live "
        "axis-aligned finite-difference APBS solve has no exact O(3) law and "
        "is only approximately invariant within the recorded 1.2e-3 V "
        "covariance-test envelope"
    ),
    validity=(
        "whole ApbsFieldResult is absent if APBS fails or any reaction "
        "potential is non-finite; potential rows have no zero sanitizer"
    ),
    irreps="", parity="mixed", tensor_rank=0)
_set_contract(
    ("apbs_E_clamp_mask",),
    coordinate_frame="lab_fixed_apbs_finite_difference_grid",
    transformation=(
        "continuum rotation/translation/reflection-invariant scalar "
        "threshold diagnostic derived from |E|; the live axis-aligned "
        "finite-difference APBS solve has no exact O(3) law"
    ),
    validity=(
        "uint8 0/1; 1 iff the canonical reaction E-field row was "
        "magnitude-clamped; whole ApbsFieldResult is absent if APBS fails"
    ),
    irreps="", parity="mixed", tensor_rank=0)
_set_contract(
    ("apbs_E_clamp_scale",),
    coordinate_frame="lab_fixed_apbs_finite_difference_grid",
    transformation=(
        "continuum rotation/translation/reflection-invariant scalar derived "
        "from |E| and the configured clamp threshold; the live axis-aligned "
        "finite-difference APBS solve has no exact O(3) law"
    ),
    validity=(
        "finite in (0,1] when the APBS result exists; 1 when unclamped; "
        "whole ApbsFieldResult is absent if APBS fails"
    ),
    irreps="", parity="mixed", tensor_rank=0)
_set_contract(
    ("apbs_nonfinite_sanitizer_mask",),
    coordinate_frame="lab_fixed_apbs_finite_difference_grid",
    transformation=(
        "continuum rotation/translation/reflection-invariant finite-value "
        "diagnostic; the live axis-aligned finite-difference APBS solve has "
        "no exact O(3) outcome law"
    ),
    validity=(
        "uint8 bit field: bit0 reaction E, bit1 reaction EFG, bit2 total E, "
        "bit3 total EFG; a set bit means that derivative was sanitized, while "
        "zero means no sanitizer fired; whole result is absent if APBS fails"
    ),
    irreps="", parity="mixed", tensor_rank=0)
_set_contract(
    ("water_efield", "water_efield_first"),
    coordinate_frame=_CARTESIAN_FRAME, transformation=_POLAR_VECTOR,
    validity="whole WaterFieldResult is absent without solvent; clamp provenance is emitted separately")
_set_contract(
    ("aimnet2_E", "aimnet2_E_backbone", "aimnet2_E_sidechain",
     "aimnet2_E_aromatic", "aimnet2_charge_response_gradient"),
    coordinate_frame=_CARTESIAN_FRAME, transformation=_POLAR_VECTOR,
    validity="AIMNet2 calculation is required as a whole; no manufactured per-row fallback")
_set_contract(
    ("aimnet2_charge_response_gradient_scalar",),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact rotation/translation/reflection-invariant scalar: per-atom "
        "Euclidean norm of the polar charge-response gradient"
    ),
    validity=(
        "whole AIMNet2 charge-response-gradient result is absent if the "
        "model/gradient calculation is unavailable; no manufactured zero"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("tripeptide_bb_residual_vec", "tripeptide_neighbor_residual_vec_prev",
     "tripeptide_neighbor_residual_vec_next"),
    coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "polar_vector under proper rotations: v'=R v; lookup/alignment is "
        "L-amino-acid chirality-conditioned and has no improper-transform contract"
    ),
    validity="NaN where the corresponding typed DFT match/direction is absent")
_set_contract(
    ("tripeptide_bb_match_distance",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "residual magnitude invariant under proper rotations; typed "
        "L-amino-acid lookup/proper-Kabsch alignment has no improper-transform "
        "contract against the unchanged chiral DFT source"
    ),
    validity="NaN where tripeptide_bb_method_tag.npy is zero",
    parity="mixed")
_set_contract(
    ("tripeptide_bb_match_atoms",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "all lookup-derived columns are invariant under proper rotations; "
        "the chiral lookup/proper-Kabsch match outcome (cols1:5) has no "
        "improper-transform contract against the unchanged DFT source"
    ),
    validity=(
        "column1 is has_match; unmatched col2 DFT index uses legacy zero, "
        "col3 distance is NaN, and col4 method tag is zero"
    ), parity="mixed")
_set_contract(
    ("tripeptide_bb_method_tag",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "proper-rotation invariant lookup/method outcome; no improper-transform "
        "contract for the chiral typed lookup/proper-Kabsch alignment"
    ),
    validity="zero means no matched DFT source; 1/2 identify the matched method",
    parity="mixed")
_set_contract(
    ("tripeptide_neighbor_reference",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "AAA source-reference cols0:4 are proper-rotation invariant; col4 "
        "any-mixed-method is lookup-outcome-dependent and has no improper-transform "
        "contract for the chiral neighbor lookup"
    ),
    validity="protein-level row; NaN reference angles mean the AAA source is unavailable",
    parity="mixed")
_set_contract(
    ("atom_sasa", "atom_sasa_fraction"),
    coordinate_frame="lab_fixed_fibonacci_sampling_grid",
    transformation=(
        "continuum rotation-invariant scalar; live finite lab-fixed Fibonacci "
        "estimator has no exact O(3) law and is only approximately invariant "
        "within the recorded covariance-test envelope"
    ),
    validity=(
        "whole SasaResult is absent when sasa_n_points is invalid; ordinary "
        "supported-element rows are finite and have no per-row validity mask"
    ), irreps="", parity="mixed", tensor_rank=0)
_set_contract(
    ("sasa_normal",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "outward_polar_vector: v'=R v in the continuum limit; finite fixed "
        "lab-axis Fibonacci sampling has no exact O(3) law and is only "
        "approximately rotation-covariant within the recorded test envelope"
    ),
    validity=(
        "zero for fully buried atoms or near-cancelling exposed samples; no "
        "separate normal-validity mask"
    ), irreps="", parity="mixed", wrapper=np.ndarray, e3nn_export="")
_set_contract(
    ("bs_total_B", "bs_ring_B_field", "hm_ring_B_field"),
    coordinate_frame=_CARTESIAN_FRAME, transformation=_AXIAL_VECTOR,
    validity=(
        "dense totals use physical zero when no accepted source; sparse ring "
        "fields are row-aligned with ring_contributions.npy"
    ))


# Full Cartesian rank-2 tensors and project-native T2 projections.
_set_contract(
    ("bs_shielding", "hm_shielding"),
    coordinate_frame=_CARTESIAN_FRAME, transformation=_EVEN_RANK2,
    validity="physical zero when no accepted aromatic-ring source")
_set_contract(
    ("bs_per_type_T1", "hm_per_type_T1"),
    coordinate_frame=_CARTESIAN_FRAME, transformation=_AXIAL_VECTOR,
    validity=(
        "eight contiguous Cartesian Levi-Civita-dual T1 blocks in RingTypeIndex "
        "order; zero block means no accepted source of that type"
    ))
_set_contract(
    ("bs_per_type_T2", "hm_per_type_T2"),
    coordinate_frame=_CARTESIAN_FRAME, transformation=_NATIVE_T2,
    validity=(
        "eight contiguous native-T2 blocks in RingTypeIndex order; omitted "
        "T0/T1 are separate signals, not structural zeros"
    ))

_MCCONNELL_FULL9 = (
    "mc_peptide_co_fixed", "mc_peptide_co_bo", "mc_peptide_co_rhombic",
    "mc_peptide_cn_fixed", "mc_peptide_cn_bo",
    "mc_backbone_other_fixed", "mc_backbone_other_bo",
    "mc_sidechain_co_fixed", "mc_sidechain_co_bo",
    "mc_sidechain_other_fixed", "mc_sidechain_other_bo",
    "mc_disulfide_fixed", "mc_disulfide_bo",
    "mc_aromatic_zeroed_fixed", "mc_aromatic_zeroed_bo",
    "mc_backbone_xh_fixed", "mc_backbone_xh_bo",
    "mc_sidechain_xh_fixed", "mc_sidechain_xh_bo",
    "mc_s_h_fixed", "mc_s_h_bo", "mc_nearest_co_T2",
    "mc_nearest_cn_T2", "sidechain_co_fixed_T2", "sidechain_co_bo_T2",
)
_set_contract(
    _MCCONNELL_FULL9, coordinate_frame=_CARTESIAN_FRAME,
    transformation=_EVEN_RANK2,
    validity=(
        "fixed channels use physical zero for no accepted source; legacy BO "
        "channels in McConnellResult use zero when MOPAC is unavailable; "
        "sidechain_co_bo_T2 instead uses NaN when MOPAC is unavailable; "
        "nearest CO/CN tensors use NaN when the corresponding nearest distance "
        "is NO_DATA_SENTINEL"
    ))
for _stem in ("mc_aromatic_zeroed_fixed", "mc_aromatic_zeroed_bo"):
    CATALOG[_stem] = replace(
        CATALOG[_stem],
        structural_zero_components=(
            "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
        ),
    )

_set_contract(
    ("coulomb_efg",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=_EVEN_RANK2,
    validity=(
        "T0/T1 structural zeros; non-finite producer values become zero "
        "without a sanitizer mask"
    ))
_set_contract(
    ("mopac_coulomb_efg",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=_EVEN_RANK2,
    validity=(
        "T0/T1 structural zeros; non-finite producer values become zero "
        "without a sanitizer mask"
    ))
_set_contract(
    ("eeq_coulomb_efg",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=_EVEN_RANK2,
    validity="T0/T1 structural zeros; whole result absent on non-finite output")

_EFG_T2 = (
    "coulomb_efg_t2", "coulomb_efg_backbone", "coulomb_efg_sidechain",
    "coulomb_efg_aromatic", "coulomb_efg_solvent",
    "mopac_coulomb_efg_backbone", "mopac_coulomb_efg_sidechain",
    "mopac_coulomb_efg_aromatic", "eeq_coulomb_efg_backbone",
    "eeq_coulomb_efg_sidechain", "eeq_coulomb_efg_aromatic",
    "apbs_efg", "apbs_efg_total_diagnostic", "water_efg",
    "water_efg_first", "aimnet2_efg", "aimnet2_efg_aromatic",
    "aimnet2_efg_backbone", "aimnet2_efg_sidechain",
)
_set_contract(
    _EFG_T2, coordinate_frame=_CARTESIAN_FRAME,
    transformation=_NATIVE_T2,
    validity=(
        "symmetric-traceless source: T0 and Cartesian-Levi-Civita T1 are "
        "structural zeros; calculator-specific absence/clamp/sanitizer arrays apply"
    ))
for _stem in (
    "coulomb_efg_t2", "coulomb_efg_backbone",
    "coulomb_efg_sidechain", "coulomb_efg_aromatic",
):
    CATALOG[_stem] = replace(
        CATALOG[_stem],
        validity=(
            "symmetric-traceless source (T0/T1 structural zeros); non-finite "
            "force-field producer values become zero without a sanitizer mask"
        ),
    )
CATALOG["coulomb_efg_solvent"] = replace(
    CATALOG["coulomb_efg_solvent"],
    transformation=(
        "continuum even_rank2_native_T2 APBS alias: reconstruct Cartesian T, "
        "apply T'=R T R^T, then decompose in project-native T2 basis. The "
        "live finite-difference APBS source has no exact O(3) law; transformed "
        "production reruns use the recorded 4e-2 V/A^2 absolute + 5e-2 "
        "relative finite-grid envelope"
    ),
    validity=(
        "symmetric-traceless APBS reaction-field alias (T0/T1 structural "
        "zeros); absent without APBS; sanitizer provenance is "
        "apbs_nonfinite_sanitizer_mask.npy bit1"
    ),
)
for _stem in (
    "mopac_coulomb_efg_backbone", "mopac_coulomb_efg_sidechain",
    "mopac_coulomb_efg_aromatic",
):
    CATALOG[_stem] = replace(
        CATALOG[_stem],
        validity=(
            "symmetric-traceless source (T0/T1 structural zeros); non-finite "
            "MOPAC-charge producer values become zero without a sanitizer mask"
        ),
    )
for _stem in (
    "eeq_coulomb_efg_backbone", "eeq_coulomb_efg_sidechain",
    "eeq_coulomb_efg_aromatic",
):
    CATALOG[_stem] = replace(
        CATALOG[_stem],
        validity=(
            "symmetric-traceless source (T0/T1 structural zeros); whole EEQ "
            "Coulomb result is absent if any output is non-finite"
        ),
    )
CATALOG["apbs_efg"] = replace(
    CATALOG["apbs_efg"],
    transformation=(
        "continuum even_rank2_native_T2: reconstruct Cartesian T, apply "
        "T'=R T R^T, then decompose in project-native T2 basis. The live "
        "axis-aligned finite-difference APBS solve has no exact O(3) law; "
        "transformed production reruns use the recorded 4e-2 V/A^2 absolute "
        "+ 5e-2 relative finite-grid envelope"
    ),
    validity=(
        "symmetric-traceless reaction EFG (T0/T1 structural zeros); "
        "apbs_nonfinite_sanitizer_mask.npy bit1 identifies sanitized rows"
    ),
)
CATALOG["apbs_efg_total_diagnostic"] = replace(
    CATALOG["apbs_efg_total_diagnostic"],
    transformation=(
        "continuum even_rank2_native_T2: reconstruct Cartesian T, apply "
        "T'=R T R^T, then decompose in project-native T2 basis. The live "
        "axis-aligned finite-difference APBS solve has no exact O(3) law; "
        "transformed production reruns use the recorded 4e-2 V/A^2 absolute "
        "+ 5e-2 relative finite-grid envelope"
    ),
    validity=(
        "symmetric-traceless total diagnostic EFG (T0/T1 structural zeros); "
        "apbs_nonfinite_sanitizer_mask.npy bit3 identifies sanitized rows"
    ),
)
CATALOG["water_efg"] = replace(
    CATALOG["water_efg"],
    validity=(
        "symmetric-traceless water EFG (T0/T1 structural zeros); whole "
        "WaterFieldResult is absent when solvent is unavailable"
    ),
)
CATALOG["water_efg_first"] = replace(
    CATALOG["water_efg_first"],
    validity=(
        "symmetric-traceless first-shell EFG (T0/T1 structural zeros); "
        "physical zero when no water lies in the first shell; whole result "
        "is absent when solvent is unavailable"
    ),
)
for _stem in (
    "aimnet2_efg", "aimnet2_efg_aromatic",
    "aimnet2_efg_backbone", "aimnet2_efg_sidechain",
):
    CATALOG[_stem] = replace(
        CATALOG[_stem],
        validity=(
            "symmetric-traceless source (T0/T1 structural zeros); AIMNet2 "
            "model/result is required as a whole with no per-row fallback"
        ),
    )

_set_contract(
    ("orca_total", "orca_diamagnetic", "orca_paramagnetic"),
    coordinate_frame="orca_output_cartesian_xyz",
    transformation=(
        "even_rank2 in ORCA output frame: T'=R T R^T; parser-only producer "
        "does not validate or rotate ORCA coordinates into the conformation frame"
    ),
    validity="whole optional ORCA result absent when the parsed shielding source is unavailable",
    tensor_frame="orca_output_cartesian_xyz")

_TRIPEPTIDE_TENSORS = (
    "tripeptide_bb_shielding", "tripeptide_neighbor_shielding",
    "tripeptide_neighbor_shielding_prev",
    "tripeptide_neighbor_shielding_next",
)
_set_contract(
    _TRIPEPTIDE_TENSORS, coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "even_rank2 under proper rotations: T'=R T R^T; typed tripeptide "
        "lookup/Kabsch alignment is L-amino-acid chirality-conditioned and has "
        "no improper-transform contract"
    ),
    validity="NaN where the corresponding typed DFT match/direction is absent")

# The central lookup is only covariant for proper rotations because the DFT
# source is chiral.  Raw arrays preserve the PackFull9/xyz payload without
# advertising an O(3) wrapper or e3nn parity that the calculator cannot obey.
for _stem in ("tripeptide_bb_shielding", "tripeptide_bb_residual_vec"):
    CATALOG[_stem] = replace(
        CATALOG[_stem], wrapper=np.ndarray, irreps="", parity="mixed",
        e3nn_export="")

# The neighbour lookup uses the same chiral DFT source and proper-Kabsch
# alignment.  Preserve its raw PackFull9/xyz values, but do not attach a
# homogeneous O(3) parity or e3nn export contract.
for _stem in (
    "tripeptide_neighbor_shielding",
    "tripeptide_neighbor_shielding_prev",
    "tripeptide_neighbor_shielding_next",
    "tripeptide_neighbor_residual_vec_prev",
    "tripeptide_neighbor_residual_vec_next",
):
    CATALOG[_stem] = replace(
        CATALOG[_stem], wrapper=np.ndarray, irreps="", parity="mixed",
        e3nn_export="")

_LARSEN_TENSORS = (
    "larsen_hbond_shielding", "larsen_hbond_1pHB_shielding",
    "larsen_hbond_2pHB_shielding", "larsen_hbond_1pHaB_shielding",
    "larsen_hbond_2pHaB_shielding",
    "larsen_hbond_diagnostic_CB_shielding",
)
_set_contract(
    _LARSEN_TENSORS, coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "even_rank2 under proper rotations: T'=R T R^T; signed-rho DFT-grid "
        "lookup is chirality-conditioned and has no improper-transform contract"
    ),
    validity=(
        "physical zero for no contributing pair; whole optional group absent "
        "without Larsen grids; imputed corners are emitted unmasked"
    ))

# Signed-rho grid lookup is chiral and therefore has only an SO(3) contract.
# Keep the packed tensor representation explicit, but expose it as raw data so
# neither the SDK nor e3nn metadata invents an unsupported improper law.
for _stem in _LARSEN_TENSORS:
    CATALOG[_stem] = replace(
        CATALOG[_stem], wrapper=np.ndarray, irreps="", parity="mixed",
        e3nn_export="")

_MUTATION_TENSORS = (
    "delta_shielding", "wt_shielding_diamagnetic",
    "wt_shielding_paramagnetic", "mut_shielding_diamagnetic",
    "mut_shielding_paramagnetic", "delta_shielding_diamagnetic",
    "delta_shielding_paramagnetic",
)
_set_contract(
    _MUTATION_TENSORS,
    coordinate_frame="shared_wt_mut_orca_cartesian_xyz",
    transformation=(
        "even_rank2 only when WT and mutant ORCA tensors share one Cartesian "
        "frame: T'=R T R^T; MutationDeltaResult performs no frame alignment"
    ),
    validity="WT-atom rows; unmatched atoms serialize zero and delta_scalars.npy column 0 is the match mask",
    tensor_frame="shared_wt_mut_orca_cartesian_xyz")


# Local-frame and mixed-column payloads.
_set_contract(
    ("ring_contributions",), coordinate_frame="mixed_intrinsic_and_conformation_cartesian_xyz",
    transformation=(
        "mixed blocks: cols0:5,6:9,36:40 invariant; col5 pseudoscalar "
        "z'=det(R)z; col6 is the acute atan2(rho,abs(z)) angle and remains "
        "invariant; cols9:18 BS full even-rank2; cols18:27 HM symmetric-"
        "traceless even-rank2 with structural-zero T0/T1; cols27:36 HM full "
        "even-rank2; translation invariant"
    ),
    validity=(
        "sparse rows are atom/ring identity-aligned; azimuth cols38:40 are "
        "NaN for an invalid vertex-0 gauge; tensor blocks otherwise use producer zeros"
    ), irreps="", parity="mixed", tensor_rank=0,
    tensor_basis=(
        "mixed row; tensor blocks cols9:18,18:27,27:36 each use "
        + _FULL9_BASIS
    ),
    tensor_component_order=(
        "within each tensor block cols9:18,18:27,27:36: " + _FULL9_ORDER
    ),
    tensor_frame=(
        "tensor blocks cols9:36 conformation_cartesian_xyz; other columns "
        "intrinsic identifiers/scalars"
    ),
    structural_zero_components=(
        "cols18:22 (HM symmetric-traceless H: T0,T1_x,T1_y,T1_z)"
    ), e3nn_export="")
_set_contract(
    ("bs_ring_B_cylindrical",), coordinate_frame="ring_cylindrical_components",
    transformation=(
        "local components [B_rho,structural_zero_B_phi,B_z]: invariant under "
        "proper rotations; under improper transforms B_rho is pseudoscalar, "
        "B_phi remains structural zero, and B_z is invariant"
    ),
    validity=(
        "row-aligned with ring_contributions.npy; B_phi is structural zero; "
        "rho-axis degeneracy serializes B_rho=0 without a separate mask"
    ), structural_zero_components="B_phi", irreps="", parity="mixed",
    tensor_rank=0, wrapper=np.ndarray)
_set_contract(
    ("ring_geometry",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "mixed blocks: cols0:3 invariant ids; cols3:6 affine position "
        "p'=R p+t; cols6:9 axial normal a'=det(R)R a for an ordinary "
        "nondegenerate ring; col9 invariant radius. A collinear/collapsed "
        ">=3-member ring has an underdetermined SVD normal with no closed law"
    ),
    validity=(
        "no explicit normal-valid mask; Ring::ComputeGeometry unconditionally "
        "uses SVD V.col(2) for >=3 vertices, so collinear/collapsed geometry "
        "can carry a lab-basis-dependent unit normal; zero is only the "
        "<3-member/default case"
    ),
    parity="mixed")
_set_contract(
    ("ring_pair_geometry",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "mixed intrinsic scalars: cols0:9,11:13 invariant; signed normal "
        "offsets cols9:11 are pseudoscalars s'=det(R)s; translation invariant"
    ),
    validity="one deterministic i<j row per aromatic-ring pair; no separate validity mask",
    parity="mixed")
_set_contract(
    ("spatial_neighbors",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "mixed sparse row: cols0:2 invariant atom ids; cols2:5 polar unit "
        "vector v'=R v; col5 invariant distance; translation invariant"
    ),
    validity="zero-distance/self rows are omitted; retained row identities are deterministic",
    parity="mixed")
_set_contract(
    ("hbond_pairs_geometry",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "mixed sparse row: cols0:2 invariant distance/angle; cols2:5 polar "
        "H-to-O unit vector v'=R v; translation invariant"
    ),
    validity="hbond_pairs_angle_valid.npy applies to the angle column",
    parity="mixed")
_set_contract(
    ("cb_deviation",), coordinate_frame="intrinsic_chiral_lookup",
    transformation=(
        "rotation-invariant under proper rotations; ideal-L-CB construction "
        "is chirality-conditioned and has no improper-transform contract"
    ),
    validity="cb_deviation_valid.npy; invalid/non-applicable rows are NaN",
    irreps="", parity="mixed", tensor_rank=0)
_set_contract(
    ("cb_residual_vector",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "polar displacement under proper rotations: v'=R v; the ideal-L-CB "
        "construction mixes bond vectors with a cross product and therefore "
        "has no single improper-transform parity"
    ),
    validity="cb_residual_vector_valid.npy; invalid/non-applicable rows are NaN",
    irreps="", parity="mixed", wrapper=np.ndarray)

_set_contract(
    ("piquad_quad_scalar",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "rotation_invariant scalar under proper and improper transforms: "
        "q' = q for q=(3*cos(theta)^2-1)/r^4; translation invariant"
    ),
    validity=(
        "row-aligned with ring_contributions.npy; no dedicated validity mask; "
        "zero is ambiguous between the physical magic-angle value and the "
        "legacy/default value for a row not accepted by PiQuadrupoleResult"
    ))

_set_contract(
    ("piquad_local_tensor", "piquad_local_T2"),
    coordinate_frame="ring_local_vertex0_gauge",
    transformation=(
        "ring-local symmetric-traceless even-rank2: coefficients are invariant "
        "under a global proper rotation rerun; under an improper transform the "
        "local Cartesian matrix obeys T_local'=P T_local P with P=diag(1,1,-1)"
    ),
    validity=(
        "NaN for invalid local frame/singular evaluation; "
        "piquad_local_geometry.npy column 6 is combined tensor-evaluation "
        "validity (valid frame and nonsingular distance)"
    ), irreps="", parity="mixed", e3nn_export="", wrapper=np.ndarray)
_set_contract(
    ("piquad_local_frame",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "mixed basis columns [x,y,z]: x and y are polar v'=R v; z is axial "
        "a'=det(R)R a; vertex-0 gauge"
    ),
    validity=(
        "all nine components are NaN unless the combined tensor evaluation "
        "succeeds; piquad_local_geometry.npy column 6 is the combined validity"
    ), parity="mixed")
_set_contract(
    ("piquad_local_geometry",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "mixed row: ids/distance/axial scalar/flags invariant; col4 cos_theta "
        "is pseudoscalar s'=det(R)s; translation invariant"
    ),
    validity=(
        "column 6 is combined tensor-evaluation validity; distance is computed "
        "before frame construction and cos_theta before vertex-0 frame "
        "construction, so either may remain finite when column 6 is zero"
    ), parity="mixed")
_set_contract(
    ("sidechain_co_frame",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "mixed frame [origin,x,y,z]: cols0:3 affine position p'=R p+t; "
        "x/y polar; z cross-product axial"
    ),
    validity=(
        "origin cols0:3 are assigned before frame validation and may remain "
        "finite on an invalid row; axes cols3:12 are NaN when invalid; "
        "sidechain_co_frame_quality.npy column 3 is frame_valid"
    ), parity="mixed")
_set_contract(
    ("water_polarization",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "mixed row: cols0:3 net-water-dipole polar vector; cols3:6 outward "
        "SASA normal is continuum polar but live finite-grid non-covariant; "
        "cols6:8 are continuum invariant but inherit the live finite-grid "
        "normal; cols8:10 are exact O(3)-invariant scalars"
    ),
    validity=(
        "whole result absent without solvent; no-shell dipole/scalars use "
        "physical zero; copied SASA normal inherits its zero/finite-grid semantics"
    ), parity="mixed")
_set_contract(
    ("water_hbond_candidates",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "mixed sparse row: cols0:9,15 invariant; cols9:12 water-O affine "
        "position p'=R p+t; cols12:15 selected-water-H affine position"
    ),
    validity="selected-H position cols12:15 are NaN in mode 2 where water is the acceptor",
    parity="mixed")

_set_contract(
    ("delta_apbs",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "mixed WT-atom row: cols0:3 polar delta-E; cols3:7 compatibility "
        "structural zeros (T0/T1); cols7:12 native T2 even-rank2"
    ),
    validity="unmatched WT atoms serialize zero; delta_scalars.npy column 0 is the match mask",
    parity="mixed")
_set_contract(
    ("delta_ring_proximity",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "repeated 6-column removed-ring blocks: distance/rho/factor/decay "
        "invariant; z pseudoscalar; theta maps to pi-theta under an improper transform"
    ),
    validity="unmatched WT atoms serialize zero; delta_scalars.npy column 0 is the match mask",
    parity="mixed")


# External-frame mixed diagnostics.
_set_contract(
    ("mopac_global",), coordinate_frame="mopac_output_cartesian_xyz",
    transformation=(
        "mixed graph row: col0 heat of formation invariant; cols1:4 molecular "
        "dipole polar vector. For a charged system the dipole is origin- and "
        "translation-dependent according to the MOPAC source convention"
    ),
    validity="missing parsed MOPAC dipole block currently serializes zero without an availability mask",
    parity="mixed")
_set_contract(
    ("mopac_atom_populations",), coordinate_frame="mopac_output_cartesian_xyz",
    transformation=(
        "mixed atom row: intended cols6:9 per-atom polar dipole components; "
        "other columns are scalar populations/valencies"
    ),
    validity=(
        "live parser never populates dipole cols6:10, so they serialize NaN; "
        "no per-column availability mask"
    ), parity="mixed")
_set_contract(
    ("mopac_atomic_orbital_populations",),
    coordinate_frame="mopac_output_atomic_orbital_axes",
    transformation=(
        "frame-dependent diagonal s/px/py/pz/d populations are not a closed "
        "O(3) representation; omitted AO density coherences prevent rotation"
    ),
    validity="parsed K-row diagnostic; missing printed cells are NaN and atom identity is not serialized",
    parity="mixed")
_set_contract(
    ("mopac_atomic_orbital_population_totals",),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation="rotation_invariant shell totals [s_total,p_total,d_total]",
    validity="parsed K-row diagnostic; NaN propagates from missing printed AO populations")
_set_contract(
    ("gromacs_energy",), coordinate_frame="gromacs_simulation_cartesian_xyz",
    transformation=(
        "mixed graph row: cols0:23 and41:43 scalar diagnostics; cols23:32 "
        "virial and cols32:41 pressure are full even-rank2 row-major tensors "
        "T'=R T R^T; box lengths cols20:23 are axis-tied cell diagnostics"
    ),
    validity="external EDR data are not recomputed from conformation coordinates; unavailable EDR terms are NaN",
    units=(
        "mixed: cols0:16 kJ/mol; col16 K; col17 bar; col18 nm^3; "
        "col19 kg/m^3; cols20:23 nm; cols23:32 kJ/mol virial; "
        "cols32:41 bar; cols41:43 K"
    ), parity="mixed")


# Reflection-sensitive scalar/diagnostic channels.  Each entry states whether
# its production rerun is O(3)-invariant or has an improper-transform law.
_set_contract(
    ("dssp_observed",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact rotation/translation/reflection-invariant atom-broadcast "
        "libdssp residue-observation mask for a fixed typed topology"
    ),
    validity="int8 1 means observed/mapped; 0 means every DSSP payload for the atom is unavailable")
_set_contract(
    ("dssp_backbone",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "mixed atom-broadcast row: phi/psi cols0:2 are signed dihedral "
        "pseudoscalars (wrapped sign reversal); col2 SASA and helix/sheet "
        "cols3:5 have physical O(3)-invariant laws. The live temporary PDB "
        "rounds coordinates to 0.001 Angstrom, so transformed production "
        "reruns can show bounded scalar drift or cross a category boundary"
    ),
    validity="phi/psi/SASA NaN for an unobserved DSSP residue; dssp_observed.npy is the mask",
    parity="mixed")
_set_contract(
    ("dssp_ss8",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "physical O(3)-invariant libdssp eight-class one-hot after the PPII-"
        "to-coil collapse; transformed production reruns can cross a category "
        "boundary because the temporary PDB rounds coordinates to 0.001 "
        "Angstrom"
    ),
    validity=(
        "all-zero row for an unobserved DSSP residue; dssp_observed.npy is "
        "the observation mask"
    ), parity="even")
_set_contract(
    ("dssp_ppii",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "proper-rotation invariant chirality-conditioned categorical flag; "
        "no homogeneous improper-transform law because reflection reverses "
        "the signed phi/psi region used by libdssp's PPII classifier; 0.001-"
        "Angstrom temporary-PDB rounding can affect a boundary value"
    ),
    validity="int8 -1 means no DSSP observation; otherwise 0/1 is the recomputed class",
    parity="mixed")
_set_contract(
    ("dssp_chi",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "four [cos,sin,exists] blocks: cos/exists invariant; sin is "
        "pseudoscalar and changes sign under an improper transform"
    ),
    validity="undefined chi blocks serialize [0,0,0]; exists is the block-local mask",
    parity="mixed")
_set_contract(
    ("dssp_torsion_angle",), coordinate_frame=_INTRINSIC_FRAME,
    transformation="signed_dihedral_pseudoscalars: angle'=det(R) angle with canonical wrapping",
    validity="NaN where undefined; dssp_torsion_valid.npy is the elementwise mask",
    irreps="0o", parity="odd")
_set_contract(
    ("dssp_torsion_sin",), coordinate_frame=_INTRINSIC_FRAME,
    transformation="pseudoscalar_sines: sin'=det(R) sin",
    validity="NaN where undefined; dssp_torsion_valid.npy is the elementwise mask",
    irreps="0o", parity="odd")
_set_contract(
    ("dssp_torsion_cos",), coordinate_frame=_INTRINSIC_FRAME,
    transformation="rotation_invariant cosines of signed dihedrals",
    validity="NaN where undefined; dssp_torsion_valid.npy is the elementwise mask")
_set_contract(
    ("dssp_torsion_valid",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact rotation/translation/reflection-invariant elementwise "
        "availability mask for ordinary nondegenerate torsions"
    ),
    validity="uint8 1 means the corresponding angle/sine/cosine is finite; 0 means NaN")
_set_contract(
    ("dssp_hbond_energy",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "continuum rotation/translation/reflection-invariant DSSP "
        "electrostatic H-bond energies; transformed production reruns use a "
        "5e-3 kcal/mol absolute envelope because the temporary PDB rounds "
        "coordinates to 0.001 Angstrom"
    ),
    validity=(
        "dssp_observed.npy is the observation mask; use the corresponding "
        "dssp_hbond_partner_residue_index.npy slot to distinguish a mapped "
        "partner from the legacy zero no-partner value"
    ))
_set_contract(
    ("dssp_hbond_partner_residue_index",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "continuum O(3)-invariant residue identity of the distance-derived "
        "DSSP H-bond partner; a transformed production rerun can select a "
        "different boundary partner after 0.001-Angstrom temporary-PDB "
        "rounding"
    ),
    validity=(
        "nonnegative values are protein residue indices; -1 means no mapped "
        "partner or no DSSP observation; dssp_observed.npy is the observation mask"
    ))
_set_contract(
    ("omega_actual", "omega_deviation"), coordinate_frame=_INTRINSIC_FRAME,
    transformation="signed_dihedral_pseudoscalar with canonical wrapping under improper transforms",
    validity="NaN where undefined; omega_valid.npy is the mask",
    irreps="0o", parity="odd")
_set_contract(
    ("omega_sin",), coordinate_frame=_INTRINSIC_FRAME,
    transformation="pseudoscalar_sine: sin'=det(R) sin",
    validity="NaN where undefined; omega_valid.npy is the mask",
    irreps="0o", parity="odd")
_set_contract(
    ("omega_cos",), coordinate_frame=_INTRINSIC_FRAME,
    transformation="rotation_invariant cosine of signed omega",
    validity="NaN where undefined; omega_valid.npy is the mask")
_set_contract(
    ("aromatic_chi2",), coordinate_frame=_INTRINSIC_FRAME,
    transformation="signed_dihedral_pseudoscalar with canonical wrapping under improper transforms",
    validity="NaN when the parent residue has no valid chi2 atom quartet",
    irreps="0o", parity="odd")
_set_contract(
    ("omega_valid",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant uint8 availability for the corresponding "
        "signed peptide omega; translation invariant"
    ),
    validity=(
        "1 iff omega_actual.npy, omega_sin.npy and omega_cos.npy are finite "
        "for the residue; 0 gates their NaNs"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("pucker_Q",), coordinate_frame=_INTRINSIC_FRAME,
    transformation="rotation_invariant Cremer-Pople amplitude",
    validity=(
        "NaN for non-five-membered rings or mean-plane failure; Q remains "
        "finite (including Q=0) when only the oriented phase is degenerate"
    ))
_set_contract(
    ("pucker_theta",), coordinate_frame=_INTRINSIC_FRAME,
    transformation="oriented Cremer-Pople phase: theta'=(theta+180 degrees) mod 360 under an improper transform",
    validity="NaN for non-five-membered rings or sub-amplitude/plane degeneracy",
    parity="mixed")

_set_contract(
    ("tripeptide_bb_diagnostics",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "mixed proper-rotation invariants with signed actual torsions: cols10:12 "
        "and14:18 reverse sign (with wrapping) under an improper transform; "
        "lookup outcomes remain chirality-conditioned"
    ),
    validity="NaN fields are governed by the embedded status/has-match columns 1/2",
    parity="mixed")
_set_contract(
    ("tripeptide_neighbor_diagnostics",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "two 28-column side blocks: actual phi/psi at base+11:13 and actual "
        "chi at base+15:19 are signed; lookup outcomes are chirality-conditioned"
    ),
    validity="NaN fields are governed by each side block's status/has-match columns",
    parity="mixed")

_set_contract(
    ("larsen_hbond_pairs_geometry",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "mixed row: cols0:2 r/theta are exact O(3) invariants; col2 signed "
        "rho is a wrapped pseudoscalar and changes sign under an improper "
        "transform; cols3:5 imputation provenance has no improper-transform "
        "law; col5 frame_valid is exactly O(3)-invariant"
    ),
    validity=(
        "col5 is donor-frame validity; r/theta/rho cols0:3 are computed before "
        "donor-frame validation and may remain finite when col5 is zero; "
        "imputed provenance is cols3:5"
    ), parity="mixed")
_set_contract(
    ("larsen_hbond_pairs_isotropic",), coordinate_frame="larsen_signed_rho_grid",
    transformation=(
        "proper-rotation invariant grid values; no improper-transform parity "
        "because lookup is conditioned on signed rho"
    ),
    validity="NaN for non-successful pair dispositions; imputed successful values are emitted unmasked",
    parity="mixed")
_set_contract(
    ("larsen_hbond_count", "larsen_sidechain_carbonyl_pair_count"),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact rotation/translation/reflection-invariant count: successful "
        "pair membership is gated by topology, distance, and theta; rho is "
        "periodic and only selects tensor values and imputation provenance"
    ),
    validity="per-atom count; zero is a real no-contribution count",
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("larsen_hbond_water_term",),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact rotation/translation/reflection-invariant scalar: the water "
        "gate uses only topology, distance, and theta; both successful and "
        "out-of-r-grid geometric pairs suppress the term"
    ),
    validity=(
        "2.07 ppm on geometrically unpaired amide H atoms and physical zero "
        "elsewhere; whole optional group absent without Larsen grids"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("larsen_corner_imputed", "larsen_imputed_pair_count"),
    coordinate_frame="larsen_signed_rho_grid",
    transformation=(
        "proper-rotation invariant imputation provenance; no improper-transform "
        "law because reflection changes which signed-rho validity-mask corners "
        "serve the otherwise successful periodic query"
    ),
    validity="per-atom imputation provenance; zero is a real no-imputation count",
    parity="mixed")
_set_contract(
    ("larsen_hbond_pairs_index",),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "all integer identity, class, disposition, anchor, and target-mask "
        "columns are exactly rotation/translation/reflection invariant; rho is "
        "periodic and affects tensor values/imputation provenance, not success"
    ),
    validity="col6 disposition is the authoritative row-status code",
    parity="even")
_set_contract(
    ("larsen_hbond_pairs",), coordinate_frame="mixed_larsen_pair_contract",
    transformation=(
        "compatibility table: cols0:16 index and cols16:18 r/theta are exact "
        "O(3) invariants; col18 signed-rho pseudoscalar; cols19:21 imputation "
        "provenance and cols22:28 signed-rho-conditioned isotropic values "
        "have no improper parity; col21 frame_valid is exactly O(3)-invariant"
    ),
    validity=(
        "disposition is col6 and frame_valid is col21; geometry cols16:19 may "
        "remain finite when the donor frame is invalid, while unavailable "
        "lookup isotropics are NaN"
    ), parity="mixed")
_set_contract(
    ("larsen_sidechain_donor_candidates",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "mixed audit row: cols0:10 identity/distance/angle and cols11:13 "
        "status are exact O(3) invariants; col10 signed rho is a wrapped "
        "pseudoscalar and changes sign under an improper transform"
    ),
    validity="rho is NaN when acceptor frame atoms are unavailable; audit-only and never modeled by Table 2",
    parity="mixed")


del _stem
