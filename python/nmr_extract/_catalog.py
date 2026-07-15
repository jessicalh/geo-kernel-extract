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
    RingCounts,
    McConnellNearFieldCounts,
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
    "mopac_compact_pair",
    "mopac_unique_pair",
    "mopac_csc_entry",
    "mopac_csc_pointer",
    "mopac_lmo",
    "mopac_occupied_lmo",
    "mopac_virtual_lmo",
    "mopac_occupied_lmo_atom_entry",
    "mopac_virtual_lmo_atom_entry",
    "mopac_occupied_lmo_coefficient",
    "mopac_virtual_lmo_coefficient",
    "mopac_icocc_storage",
    "mopac_icvir_storage",
    "mopac_cocc_storage",
    "mopac_cvir_storage",
    "atom_neighbor_pair",
    "ring_pair",
    "larsen_hbond_pair",
    "hbond_pair",
    "larsen_sidechain_donor_pair",
    "protein_water_hbond_pair",
    "sidechain_co_source",
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
      ``"Å"``, ``"radians"``, ``"degrees"``). Dimensionless,
      categorical, index, mask and structured/mixed layouts are named
      explicitly rather than encoded as an ambiguous empty string.

    * ``sign_convention`` documents the explicit physical sign
      convention where one applies (notably ring-current shielding
      kernels: ``"σ_ab = -dB_a^sec / dB_{0,b}"``). Empty otherwise.

    * ``tensor_rank`` is the rank of one row's tensor representation:
      0 (scalar / vector-of-scalars), 1 (Vec3), 2 (Mat3 / SphericalTensor).

    * ``parity`` is ``"even"`` / ``"odd"`` under spatial inversion, or
      the legacy sentinel ``"mixed"`` when the producer boundary cannot
      certify one unconditional homogeneous improper-transform law.  That
      includes genuinely mixed-column payloads, chiral lookup outcomes, and
      outputs of a caller-selected learned model whose parity is unchecked;
      ``"mixed"`` does not by itself assert that both even and odd columns
      are present.
      shieldings (rank 2 even), axial B fields, and most scalars are
      even. Polar vector fields like E are odd. Consumers must read
      ``transformation`` for the reason a ``mixed`` row lacks an
      unconditional parity contract.

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
    # Empty when the array is already the final declared quantity or no
    # additional physical scale is part of the producer contract.
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
_BS_FIELD_SCALING = (
    "diagnostic magnetic field evaluated at exactly 1 nA source current; the "
    "producer does not apply the ring-type physical current. Sparse rows can "
    "be combined with ring identity and the manifest-resolved intensities; "
    "the dense all-types sum has lost that separation, so use the per-type "
    "BS channels for physical ppm shielding"
)
_MOPAC_DIRECT_SOURCE = (
    "Source: diskless PM7/MOZYME/1SCF through the pinned generic-Release "
    "libmopac.so.2 (MOPAC 23.2.5, commit "
    "052691223d19935a89f0fe18cd12301bd83e4201); single point "
    "(natom_move=0) at the prepared eight-decimal geometry, exact topology "
    "charge, vacuum epsilon=1, tolerance=1.0 and max_time=7200 s; "
    "convergence is not loosened. The mozyme_scf API runs in a crash-"
    "contained worker and no MOPAC text output is parsed."
)
_MOPAC_RESULT_ABSENCE = (
    "MopacResult is optional as a whole: worker/API/dimension/non-finite "
    "failure means every MOPAC core/direct NPY is absent, never zero-filled."
)
_MOPAC_COULOMB_SOURCE = (
    " Source: all-pairs Coulomb kernel over the legacy F15.6 projection of "
    "diskless-libmopac PM7 Coulson charges (not the full-precision charge "
    "NPY), with aromatic-ring membership taking precedence over cached "
    "backbone N/CA/C/O/H/HA/CB and all other sources classed sidechain."
)
_PIQUAD_SCALING = (
    "geometry-only kernel; the producer applies no physical quadrupole "
    "prefactor. extraction_manifest.json declares the multiplier as "
    "deferred_learnable"
)
_HM_SCALING = (
    "unscaled Angstrom^-1 Haigh-Mallion surface-integral geometry kernel; "
    "the producer applies no ring-current intensity or conversion to ppm. "
    "Any physical shielding scale is downstream/model-defined"
)
_MCCONNELL_SCALING = (
    "unscaled Angstrom^-3 McConnell geometry response with unit susceptibility "
    "shape; no physical susceptibility magnitude or conversion to ppm is "
    "applied. extraction_manifest.json::feature_metadata.mcconnell declares "
    "the axial scale learned; the peptide-C=O principal-shape ratios are "
    "pinned, not the emitted response magnitude"
)
_RINGCHI_SCALING = (
    "geometry-only unit-susceptibility (3*cos(theta)^2-1)/r^3 kernel; "
    "the producer applies no ring-susceptibility coefficient or conversion "
    "to ppm"
)
_DISPERSION_SCALING = (
    "geometry-only unit-C6 switched S(r)/r^6 kernel; the producer applies "
    "no C6 coefficient or conversion to a dispersion energy"
)
_RING_UNION_SCALING = (
    "mixed calculator row: BS cols9:18 are unit-current kernels governed by "
    "the ring-current manifest; HM cols18:36 are unscaled Angstrom^-1 "
    "surface-integral kernels; col7 is a unit-susceptibility Angstrom^-3 "
    "geometry scalar; col36 is a unit-C6 Angstrom^-6 geometry scalar. The "
    "producer applies none of the HM/susceptibility/C6 physical magnitudes"
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
    ArraySpec("ring_contributions","identity",  RingContributions, 40,   True,  "Per-(atom,ring) pair geometry and unscaled calculator-kernel contributions",
              native_axis="ring_contribution_pair", irreps=_SHIELD_IRREPS, units="",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="ring_current",
              scaling_contract=_RING_UNION_SCALING),
    ArraySpec("ring_direction_to_center", "identity", VectorField, 3, False, "Sparse per-(atom,ring) rows aligned to ring_contributions: RingNeighbourhood.direction_to_center vector",
              native_axis="ring_contribution_pair", irreps="1o", tensor_rank=1, parity="odd", mechanism="geometry"),
    ArraySpec("ring_geometry",    "identity",   RingGeometry,      10,   True,  "Per-aromatic-ring geometry [ring/type/residue ids, center_A, unit normal, radius_A]",
              native_axis="aromatic_ring", units="mixed_index_A_dimensionless", mechanism="topology"),
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
    ArraySpec("enrichment_flags",         "enrichment", np.ndarray, 8,    False, "Enrichment boolean flags as int8 columns: is_backbone, is_amide_H, is_alpha_H, is_methyl, is_aromatic_H, donor candidate (H with N/O parent), coarse acceptor candidate (any N/O heavy atom, not chemical acceptor typing), is_on_aromatic_residue",
              mechanism="topology"),
    ArraySpec("enrichment_parent_is_sp2", "enrichment", np.ndarray, None, False, "Per-atom uint8 compatibility flag: for H only, parent is an aromatic-ring member or the residue's cached backbone C/N; this is not a general typed-hybridisation lookup",
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
    ArraySpec("enrichment_acceptor_class", "enrichment", np.ndarray, None, False, "Typed nonpositive-formal-charge acceptor projection: 0 none, 1 backbone carbonyl O, 2 sidechain amide O, 3 carboxylate O, 4 O in aromatic hydroxyl/oxide/PlanarGroup::None, 5 unprotonated ring N, 6 N with PlanarGroup::None or S",
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
    ArraySpec("ff_partial_charge", "charge_assignment", np.ndarray, None, False, "Prepared per-atom partial charge from the configured typed ChargeSource (parameter file, AMBER/runtime-tleap prmtop, or caller-preloaded table); not hard-wired to one force field",
              units="e", mechanism="charges"),
    ArraySpec("ff_pb_radius",      "charge_assignment", np.ndarray, None, False, "Prepared per-atom Poisson-Boltzmann radius from the same configured typed ChargeSource table",
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
    ArraySpec("bs_total_B",       "biot_savart", MagneticVectorField, 3, True,  "Dense unit-current (1 nA) BS geometry field summed over modeled ring types; diagnostic TrpPerimeter rings are excluded",
              irreps="1e", units="T", tensor_rank=1, parity="even", mechanism="ring_current", scaling_contract=_BS_FIELD_SCALING),
    ArraySpec("bs_ring_B_field",  "biot_savart", MagneticVectorField, 3, True,  "Unit-current (1 nA) BS field on union (atom,ring) rows; diagnostic rows may exist and a zero is calculator-absence-or-physical-zero without a BS-specific mask",
              native_axis="ring_contribution_pair", irreps="1e", units="T", tensor_rank=1, parity="even", mechanism="ring_current", scaling_contract=_BS_FIELD_SCALING),
    ArraySpec("bs_ring_B_cylindrical", "biot_savart", MagneticVectorField, 3, True, "Unit-current (1 nA) BS per-(atom,ring) B-field in the ring cylindrical frame, on the union ring_contributions row axis",
              native_axis="ring_contribution_pair", irreps="1e", units="T", tensor_rank=1, parity="even", mechanism="ring_current", scaling_contract=_BS_FIELD_SCALING),
    ArraySpec("bs_ring_counts",   "biot_savart", RingCounts,       4,    True,  "Ring proximity counts (3/5/8/12 A)",
              mechanism="ring_current"),

    # ── Haigh-Mallion (HaighMallionResult.cpp) ───────────────────
    ArraySpec("hm_shielding",     "haigh_mallion", ShieldingTensor, 9,   True,  "Unscaled HM ring-current surface-integral geometry kernel G",
              irreps=_SHIELD_IRREPS, units="Angstrom^-1", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="ring_current", scaling_contract=_HM_SCALING),
    ArraySpec("hm_per_type_T0",   "haigh_mallion", PerRingTypeT0,   8,   True,  "Unscaled HM geometry-kernel T0 per ring type",
              irreps="0e", units="Angstrom^-1", mechanism="ring_current", scaling_contract=_HM_SCALING),
    ArraySpec("hm_per_type_T1",   "haigh_mallion", PerRingTypeT1,   24,  True,  "Unscaled HM geometry-kernel T1 per ring type",
              irreps="1e", units="Angstrom^-1", tensor_rank=1, parity="even", mechanism="ring_current", scaling_contract=_HM_SCALING),
    ArraySpec("hm_per_type_T2",   "haigh_mallion", PerRingTypeT2,   40,  True,  "Unscaled HM geometry-kernel T2 per ring type",
              irreps="2e", units="Angstrom^-1", tensor_rank=2, mechanism="ring_current", scaling_contract=_HM_SCALING),
    ArraySpec("hm_ring_B_field",  "haigh_mallion", MagneticVectorField, 3, True, "Unscaled HM per-(atom,ring) effective surface-integral field vector",
              native_axis="ring_contribution_pair", irreps="1e", units="Angstrom^-1", tensor_rank=1, parity="even", mechanism="ring_current", scaling_contract=_HM_SCALING),

    # ── Pi-Quadrupole (PiQuadrupoleResult.cpp) ───────────────────
    ArraySpec("pq_per_type_T0",   "pi_quadrupole", PerRingTypeT0,   8,   True,  "Geometry-only Buckingham A-term kernel sum (3 cos^2(theta)-1)/r^4 per aromatic ring type; not a physically scaled quadrupole response",
              irreps="0e", units="Angstrom^-4", mechanism="ring_efg", scaling_contract=_PIQUAD_SCALING),
    ArraySpec("piquad_axial_scalar_per_type_T0", "pi_quadrupole", PerRingTypeT0, 8, False, "Additive alias of pq_per_type_T0: the same geometry-only per-type kernel, optional for pre-alias extractions",
              irreps="0e", units="Angstrom^-4", mechanism="ring_efg", scaling_contract=_PIQUAD_SCALING),
    ArraySpec("piquad_quad_scalar", "pi_quadrupole", np.ndarray, None, False, "Sparse geometry-only (3 cos^2(theta)-1)/r^4 values aligned to ring_contributions; a default zero can mean PiQuadrupoleResult did not accept that union row",
              native_axis="ring_contribution_pair", irreps="0e", units="Angstrom^-4", mechanism="ring_efg", scaling_contract=_PIQUAD_SCALING),
    ArraySpec("piquad_local_tensor", "pi_quadrupole", ShieldingTensor, 9, True, "Raw geometry-only symmetric-traceless A^-5 tensor in the deterministic ring-local frame; no physical quadrupole prefactor is applied",
              native_axis="ring_contribution_pair", irreps="2e", units="Angstrom^-5", tensor_rank=2, mechanism="ring_efg",
              tensor_basis=_FULL9_BASIS, tensor_component_order=_FULL9_ORDER, tensor_frame="ring_local_vertex0_gauge", structural_zero_components="T0,T1_x,T1_y,T1_z", e3nn_export=_E3NN_EXPORT, scaling_contract=_PIQUAD_SCALING),
    ArraySpec("piquad_local_T2", "pi_quadrupole", np.ndarray, 5, True, "Project-native T2 view of the same geometry-only local A^-5 tensor; not e3nn-ready and not physically prefactored",
              native_axis="ring_contribution_pair", irreps="2e", units="Angstrom^-5", tensor_rank=2, mechanism="ring_efg",
              tensor_basis=_T2_BASIS, tensor_component_order=_T2_ORDER,
              tensor_frame="ring_local_vertex0_gauge",
              structural_zero_components=_EFG_STRUCTURAL_ZEROS,
              e3nn_export=_E3NN_EXPORT, scaling_contract=_PIQUAD_SCALING),
    ArraySpec("piquad_local_frame", "pi_quadrupole", np.ndarray, 9, True, "Deterministic ring-local frame columns [x_axis,y_axis,z_axis], vertex-0 gauge; NaN when invalid",
              native_axis="ring_contribution_pair", mechanism="geometry"),
    ArraySpec("piquad_local_geometry", "pi_quadrupole", np.ndarray, 8, True, "Row-aligned local-tensor provenance [atom, ring, type, distance_A, cos_theta, union-row axial scalar, local_tensor_valid, aromatic_only]; col5 zero is not an availability flag",
              native_axis="ring_contribution_pair",
              units="mixed:index[0:3],Å[3],dimensionless[4],Å^-4[5],mask[6:8]",
              mechanism="ring_efg",
              scaling_contract=_PIQUAD_SCALING + "; applies to col5 only; identity, distance, cosine and mask columns are not scaled"),

    # ── Ring susceptibility (RingSusceptibilityResult.cpp) ────────
    ArraySpec("ringchi_scalar", "ring_susceptibility", np.ndarray, None, True, "Sparse ring-susceptibility geometry scalar aligned to the union ring_contributions rows; zero is ambiguous between a physical zero and a row not accepted by RingSusceptibilityResult",
              native_axis="ring_contribution_pair", irreps="0e", units="Angstrom^-3", mechanism="ring_current", scaling_contract=_RINGCHI_SCALING),
    ArraySpec("ringchi_per_type_T0", "ring_susceptibility", PerRingTypeT0, 8, True, "Dense per-atom sums of accepted ring-susceptibility geometry scalars by aromatic RingTypeIndex; zero means an empty/zero accepted sum",
              irreps="0e", units="Angstrom^-3", mechanism="ring_current", scaling_contract=_RINGCHI_SCALING),

    # ── Dispersion (DispersionResult.cpp) ────────────────────────
    ArraySpec("disp_per_type_T0", "dispersion", PerRingTypeT0,     8,    True,  "Deprecated name for the per-type sum over accepted aromatic-ring vertices of the configured switched S(r)/r^6 kernel; ring vertices and atoms bonded to any vertex are excluded; not D3/D4 energy",
              irreps="0e", units="Angstrom^-6", mechanism="ring_dispersion", scaling_contract=_DISPERSION_SCALING),
    ArraySpec("aromatic_r6_proximity_per_type_T0", "dispersion", PerRingTypeT0, 8, False, "Canonical alias: per-type sum over accepted aromatic-ring vertices of configured switched S(r)/r^6, with vertex/bonded exclusions; not D3/D4 energy",
              irreps="0e", units="Angstrom^-6", mechanism="ring_dispersion", scaling_contract=_DISPERSION_SCALING),

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
    ArraySpec("mc_nearest_co_T2",         "mcconnell", ShieldingTensor, 9, False, "Nearest accepted peptide C=O axial ComputePairKernel response per atom (not the additive channel's PeptideCO rhombic source-shape response), packed in project Full9 order",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_nearest_cn_T2",         "mcconnell", ShieldingTensor, 9, False, "Nearest accepted peptide C-N axial ComputePairKernel response per atom, packed in project Full9 order",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_nearest_co_bond_index", "mcconnell", np.ndarray, None, False, "Topology bond index of the geometrically nearest accepted peptide C=O source per atom; -1 means no accepted source",
              irreps="0e", units="index", mechanism="bond_anisotropy"),
    ArraySpec("mc_nearest_cn_bond_index", "mcconnell", np.ndarray, None, False, "Topology bond index of the geometrically nearest accepted peptide C-N source per atom; -1 means no accepted source",
              irreps="0e", units="index", mechanism="bond_anisotropy"),

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
    ArraySpec("sidechain_co_fixed_T2", "sidechain_carbonyl_anisotropy", ShieldingTensor, 9, True, "Canonical fixed-unit-weight McConnell SidechainCO geometry response, full project-native SphericalTensor pack",
              native_axis="atom", irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("sidechain_co_bo_T2", "sidechain_carbonyl_anisotropy", ShieldingTensor, 9, True, "MOPAC bond-order-weighted McConnell SidechainCO geometry response, full project-native SphericalTensor pack; NaN rows when MOPAC is absent",
              native_axis="atom", irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("sidechain_co_scalar_audit", "sidechain_carbonyl_anisotropy", np.ndarray, 4, True, "Per-atom audit [fixed_T2_norm, bo_T2_norm, accepted_sidechain_CO_source_count, nearest_sidechain_CO_midpoint_distance_A]",
              native_axis="atom", units="mixed", mechanism="bond_anisotropy"),

    # ── Coulomb (CoulombResult.cpp) — vacuum force-field-charge field. ──
    ArraySpec("coulomb_efg",            "coulomb", ShieldingTensor, 9,   False, "Coulomb bare total EFG full 9-pack; T0/T1 are structural zeros (bitwise-symmetric source, explicit traceless projection)",
              irreps=_SHIELD_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),
    ArraySpec("coulomb_efg_t2",         "coulomb", EFGTensor,       5,   False, "Coulomb bare total EFG T2-only companion copied from coulomb_efg columns 4:9",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("coulomb_E",              "coulomb", VectorField,     3,   False, "Coulomb total E-field",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_E_backbone",     "coulomb", VectorField,     3,   False, "Force-field-charge Coulomb field from non-aromatic cached backbone N/CA/C/O/H/HA/CB sources (CB is included)",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_E_sidechain",    "coulomb", VectorField,     3,   False, "Force-field-charge Coulomb field from sources that are neither aromatic-ring members nor cached backbone N/CA/C/O/H/HA/CB",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_E_aromatic",     "coulomb", VectorField,     3,   False, "Force-field-charge Coulomb field from aromatic-ring-member sources; aromatic membership takes precedence over backbone classification",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_efg_backbone",   "coulomb", EFGTensor,       5,   False, "Force-field-charge EFG T2 from non-aromatic cached backbone N/CA/C/O/H/HA/CB sources (CB included)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("coulomb_efg_sidechain",  "coulomb", EFGTensor,       5,   False, "Force-field-charge EFG T2 from sources neither aromatic nor cached backbone",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("coulomb_efg_aromatic",   "coulomb", EFGTensor,       5,   False, "Force-field-charge EFG T2 from aromatic-ring members, with aromatic precedence over backbone classification",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("coulomb_scalars",        "coulomb", CoulombScalars,  4,   False, "Force-field Coulomb scalars [|E_total|, parent-to-H E_total projection (NaN when unavailable), signed E_backbone projection onto E_total direction (historical *_frac name; zero at near-zero |E_total|), |E_aromatic|]",
              units="V/A", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_aromatic_E_proj", "coulomb", np.ndarray,     None, False, "Coulomb aromatic E-field parent-to-H projection; NaN for non-H or parentless atoms",
              irreps="0e", units="V/A", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_aromatic_n_src",  "coulomb", np.ndarray,     None, False, "Count of sidechain aromatic source atoms contributing to the Coulomb aromatic field (int32)", units="count", mechanism="electrostatic_efg"),
    # ── H-Bond (HBondResult.cpp) ─────────────────────────────────
    ArraySpec("hbond_scalars",    "hbond", HBondScalars,           4,    True,  "H-bond summary [nearest accepted H-to-target distance_A, its inverse-cube_A^-3, count of accepted source evaluations within configured hbond_counting_radius (default 3.5 A), sum over every accepted evaluation of (3cos^2(theta)-1)/r^3_A^-3]; legacy zero distance/inverse-cube means no accepted source when hbond_flags col0 is zero",
              units="mixed_A_A^-3_count_A^-3", mechanism="hbond_kernel"),
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
    ArraySpec("hydration_shell",    "hydration",   np.ndarray,     4,    False, "Explicit-solvent per-atom [centroid-proxy exposed first-shell-water fraction, mean water-dipole/atom-to-water cosine, nearest-ion distance within configured cutoff, nearest-ion charge]; no shell gives 0,0 and no ion gives +inf,0",
              units="mixed_dimensionless_A_e", mechanism="solvation"),

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
    ArraySpec("eeq_coulomb_E_backbone",    "eeq_coulomb", VectorField, 3, False, "EEQ-charge field from non-aromatic cached backbone N/CA/C/O/H/HA/CB sources (CB is included)",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("eeq_coulomb_E_sidechain",   "eeq_coulomb", VectorField, 3, False, "EEQ-charge field from sources that are neither aromatic-ring members nor cached backbone N/CA/C/O/H/HA/CB",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("eeq_coulomb_E_aromatic",    "eeq_coulomb", VectorField, 3, False, "EEQ-charge field from aromatic-ring-member sources; aromatic membership takes precedence over backbone classification",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("eeq_coulomb_efg_backbone",  "eeq_coulomb", EFGTensor, 5, False, "EEQ-charge EFG T2 from non-aromatic cached backbone N/CA/C/O/H/HA/CB sources (CB included)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("eeq_coulomb_efg_sidechain", "eeq_coulomb", EFGTensor, 5, False, "EEQ-charge EFG T2 from sources neither aromatic nor cached backbone",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("eeq_coulomb_efg_aromatic",  "eeq_coulomb", EFGTensor, 5, False, "EEQ-charge EFG T2 from aromatic-ring members, with aromatic precedence over backbone classification",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("eeq_coulomb_scalars",       "eeq_coulomb", CoulombScalars, 4, False, "EEQ-charge Coulomb scalars [|E_total|, parent-to-H E_total projection (NaN when unavailable), signed E_backbone projection onto E_total direction (historical *_frac name; zero at near-zero |E_total|), |E_aromatic|]",
              units="V/A", mechanism="electrostatic_efg"),
    ArraySpec("eeq_coulomb_aromatic_E_proj", "eeq_coulomb", np.ndarray, None, False, "EEQ-charge Coulomb aromatic E-field parent-to-H projection; NaN for non-H or parentless atoms",
              irreps="0e", units="V/A", mechanism="electrostatic_efg"),
    ArraySpec("eeq_coulomb_aromatic_n_src", "eeq_coulomb", np.ndarray, None, False, "Count of sidechain aromatic source atoms contributing to the EEQ-charge Coulomb aromatic field (int32)",
              units="count", mechanism="electrostatic_efg"),

    # ── GROMACS energy (GromacsEnergyResult.cpp) ────────────────
    ArraySpec("gromacs_energy",     "gromacs",     np.ndarray,     43,   False, "Per-frame energy (43 cols: electrostatic 3, bonded 6, VdW 3, thermo 8, box 3, virial 9, pressure tensor 9, T_group 2)",
              native_axis="protein", units="kJ/mol", mechanism="gromacs_runtime"),

    # ── Bonded energy (BondedEnergyResult.cpp) ─────────────────
    ArraySpec("bonded_energy",      "bonded",      np.ndarray,      7,   False, "Per-atom equal-share attribution of locally evaluated supplied BondedParameters using CHARMM36m-style forms [bond, angle, Urey-Bradley, proper, improper, CMAP, total]; not a GROMACS EDR per-atom decomposition",
              units="kJ/mol", parity="mixed", mechanism="bonded_parameters",
              coordinate_frame="intrinsic_signed_valence_geometry",
              transformation="translation/proper-rotation invariant local bonded energies; no unconditional improper-transform law because signed dihedrals feed supplied phase offsets and an arbitrary supplied CMAP grid",
              validity="physical equal-share zero when no supplied bonded term touches an atom; total is the sum of columns0:6"),

    # ── Legacy trajectory Welford NPY projections ───────────────
    # The production trajectory CLI writes the richer H5 groups. These three
    # NPYs remain live through TrajectoryProtein::WriteFeatures and therefore
    # still belong in the producer format contract.
    ArraySpec("hm_welford", "trajectory_rollup", np.ndarray, 11, False,
              "(N,11) unscaled Haigh-Mallion geometry-kernel trajectory rollup [T0 mean/std/min/max, |T2| mean/std/min/max, consecutive-frame T0-delta mean/std, n_frames]; source is hm_shielding_contribution and no directional component is serialized",
              native_axis="atom", irreps="0e", units="mixed_A^-1_count", mechanism="ring_current",
              coordinate_frame="intrinsic_statistics", transformation="exact O(3)-invariant scalar statistics of T0 and the project-native T2 norm; frame count invariant",
              validity="optional output of TrajectoryProtein::WriteFeatures when HmWelfordTrajectoryResult is attached; with zero primary samples mean/std are NaN and min/max retain +/-infinity sentinels; delta mean/std columns8:10 are NaN with fewer than two frames, while zero is physical only after at least one actual delta sample",
              scaling_contract=_HM_SCALING + "; scale signed and norm statistics according to their semantics; frame count is unscaled"),
    ArraySpec("mc_welford", "trajectory_rollup", np.ndarray, 11, False,
              "(N,11) unscaled McConnell unit-susceptibility geometry-response trajectory rollup [T0 mean/std/min/max, |T2| mean/std/min/max, consecutive-frame T0-delta mean/std, n_frames]; source is mc_shielding_contribution and no directional component is serialized",
              native_axis="atom", irreps="0e", units="mixed_A^-3_count", mechanism="bond_anisotropy",
              coordinate_frame="intrinsic_statistics", transformation="exact O(3)-invariant scalar statistics of T0 and the project-native T2 norm; frame count invariant",
              validity="optional output of TrajectoryProtein::WriteFeatures when McConnellWelfordTrajectoryResult is attached; with zero primary samples mean/std are NaN and min/max retain +/-infinity sentinels; delta mean/std columns8:10 are NaN with fewer than two frames, while zero is physical only after at least one actual delta sample",
              scaling_contract=_MCCONNELL_SCALING + "; scale signed and norm statistics according to their semantics; frame count is unscaled"),
    ArraySpec("sasa_welford", "trajectory_rollup", np.ndarray, 7, False,
              "(N,7) SASA trajectory rollup [SASA mean/std/min/max, consecutive-frame SASA-delta mean/std, n_frames] from the finite lab-fixed Fibonacci surface estimator",
              native_axis="atom", units="mixed_A^2_count", parity="mixed", mechanism="solvation",
              coordinate_frame="lab_fixed_fibonacci_surface_grid", transformation="cols0-5 inherit the finite sampled SASA estimator and have no exact O(3) law; col6 frame count is exact invariant",
              validity="optional output of TrajectoryProtein::WriteFeatures when SasaWelfordTrajectoryResult is attached; with zero primary samples mean/std are NaN and min/max retain +/-infinity sentinels; delta mean/std columns4:6 are NaN with fewer than two frames, while zero is physical only after at least one actual delta sample"),

    # ── MOPAC direct API (MopacResult.cpp) ───────────────────────
    # Structural N/A for this SCF configuration: mopac_properties
    # coord_update/coord_deriv have zero length because natom_move=0;
    # freq/disp are unavailable because this is not a vibration job;
    # lattice_update/lattice_deriv have zero length and stress is unavailable
    # because nlattice=nlattice_move=0. nerror/error_msg are control-plane:
    # any nonzero nerror hard-aborts before emission, so they are not data NPYs.
    ArraySpec("mopac_charges", "mopac_core", np.ndarray, None, False,
              "(N,) legacy printed-format charge projection reconstructed from mopac_properties.charge at F15.6 (six decimal places). The unquantized binary64 API values are in mopac_charges_full_precision; that name does not claim physical exactness. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="e", mechanism="charges"),
    ArraySpec("mopac_scalars", "mopac_core", MopacScalars, 4, False,
              "(N,4) legacy printed-format projection [charge,s population,p population,compact-bond valency]: charge is F15.6, s/p are F12.5, and valency is the per-atom sum of F6.3 orders greater than 0.01 among MOPAC's sorted first six compact bond-row entries, all reconstructed from the charge, atom_ao_density, and bond CSC structs. Atoms without p AOs carry the mathematical p-population zero here, matching the old scalar surface. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="mixed_e_dimensionless", mechanism="charges"),
    ArraySpec("mopac_bond_orders", "mopac_core", BondOrders, 3, False,
              "(C,3) legacy compact-table rows [atom_i,atom_j,order], reconstructed from mopac_properties' CSC bond fields by reproducing MOPAC's sorted first-six entries per atom, F6.3 (three-decimal) formatting, strict >0.01 retention, and unordered-pair maximum. C can be smaller than the retained API-unique pair count U; the unquantized set is in mopac_bond_orders_full_precision. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_compact_pair", units="mixed_index_dimensionless", mechanism="charges"),
    ArraySpec("mopac_bond_neighbors", "mopac_core", np.ndarray, 4, False,
              "(2C,4) directed read-back of the legacy compact-pair projection [atom_i,atom_j,F6.3 order,topology_bond_index], sorted by decreasing order within atom_i; the pair set comes from MOPAC's reconstructed first-six/>0.01 compact rows, and topology index -1 means no covalent-topology bond. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_bond_neighbor_pair", units="mixed_index_dimensionless", mechanism="charges"),
    ArraySpec("mopac_global", "mopac_core", MopacGlobal, 4, False,
              "(4,) legacy printed-format graph vector [heat_kcal_mol,dipole_x,dipole_y,dipole_z], reconstructed from mopac_properties.heat/dipole with heat at five decimal places and each dipole component at three. Unquantized binary64 API values are emitted separately; that split does not claim physical exactness. " + _MOPAC_DIRECT_SOURCE,
              native_axis="protein", units="mixed_kcal_per_mol_Debye", mechanism="quantum_reference"),
    ArraySpec("mopac_atom_populations", "mopac_core", MopacAtomPopulations, 12, False,
              "(N,12) legacy printed-format columns [charge,electron trace,s,p,d,N/A f,N/A per-atom dipole x/y/z/total,CSC-diagonal valency,compact-bond valency], reconstructed from structs. Precisions are 6,4,5,5,5 decimal places respectively for charge/trace/s/p/d, three for the CSC diagonal, and compact-bond valency is the sum of retained first-six F6.3 orders. p/d are NaN without a live p/d shell; the f and per-atom-dipole text-only columns are NaN, never zero-filled. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="mixed_e_dimensionless", mechanism="quantum_reference"),
    ArraySpec("mopac_atomic_orbital_populations", "mopac_core", MopacAtomicOrbitalPopulations, 9, False,
              "(N,9) legacy F10.5 (five-decimal) AO-population projection, reconstructed from the diagonal of mopac_properties.atom_ao_density in MOPAC s/px/py/pz/d order only through each atom's live ao_orbitals width; non-existent per-atom AOs are NaN/N/A. The unquantized globally padded API blocks remain in mopac_atom_ao_density. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="electron", mechanism="quantum_reference"),
    ArraySpec("mopac_atomic_orbital_population_totals", "mopac_core", MopacAtomicOrbitalPopulationTotals, 3, False,
              "(N,3) legacy shell totals [s,p,d], formed exactly as before by summing the individually F10.5-projected AO entries; the resulting sums are not a second quantization of the separately printed 3F12.5 shell table. p/d are NaN where that atom has no live p/d shell. Unquantized binary64 API-derived sums are emitted separately as mopac_atom_{s,p,d}_population. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="electron", mechanism="quantum_reference"),
    ArraySpec("mopac_bond_valencies", "mopac_core", np.ndarray, None, False,
              "(N,) legacy three-decimal printed-format projection of the diagonal valencies in mopac_properties' CSC Wiberg matrix, not recomputed from off-diagonal orders. Unquantized binary64 API diagonal values are in mopac_bond_valencies_full_precision; that name does not claim physical exactness. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="dimensionless", mechanism="charges"),
    ArraySpec("mopac_bond_orders_unique", "mopac_core", MopacUniqueBondOrders, 8, False,
              "(U,8) legacy-compatible API-unique a<b projection [a,b,max_order,mean_order,NaN,NaN,NaN,topology_bond_index]. Each of the two symmetric CSC entries is first put through the F6.3 compatibility projection, then max/mean are formed as the retired text parser did; U is the complete retained API sparse pair set, independent of compact first-six visibility. Columns 4, 5, and 6 were ALLBONDS text-only count/entry indices and are explicitly NaN rather than fabricated; col7 is NaN when the retained pair has no covalent-topology bond. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_unique_pair", units="mixed_index_dimensionless", mechanism="charges"),
    ArraySpec("mopac_topology_bond_orders_full", "mopac_core", MopacTopologyBondOrdersFull, 8, False,
              "(B,8) legacy-compatible topology bridge [bond_index,a,b,three-decimal API-pair order,present,API_unique_pair_index,API_sparse_absence_flag,NaN]. Order and unique-pair index are NaN when the topology pair is absent from mopac_properties' sparse CSC set; the final ALLBONDS printed-entry-count column is text-only and explicitly NaN. This API-pair subset is independent of compact first-six visibility. " + _MOPAC_DIRECT_SOURCE,
              native_axis="bond", units="mixed_index_dimensionless", mechanism="charges"),

    ArraySpec("mopac_charges_full_precision", "mopac_direct", np.ndarray, None, False,
              "(N,) unquantized binary64 net atomic Coulson charges copied from mopac_properties.charge without legacy decimal projection; full_precision denotes that compatibility split, not physical exactness. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="e", mechanism="charges"),
    ArraySpec("mopac_bond_orders_full_precision", "mopac_direct", BondOrders, 3, False,
              "(U,3) all retained libmopac sparse-API unique a<b rows [atom_i,atom_j,Wiberg order]; values are unquantized binary64 API values (the meaning of full_precision here), with no compact first-six/F6.3 projection. Pairs omitted by libmopac's strict >0.01 sparse finalizer have no row and are not measured zeros. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_unique_pair", units="mixed_index_dimensionless", mechanism="charges"),
    ArraySpec("mopac_bond_valencies_full_precision", "mopac_direct", np.ndarray, None, False,
              "(N,) unquantized binary64 diagonal valencies copied from mopac_properties.bond_order's CSC diagonal entries without legacy decimal projection; full_precision is not a claim of physical exactness. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="dimensionless", mechanism="charges"),

    ArraySpec("mopac_heat_kcal_mol", "mopac_direct", np.ndarray, None, False,
              "(1,) unquantized binary64 PM7 heat of formation copied from mopac_properties.heat; no claim of physical exactness. " + _MOPAC_DIRECT_SOURCE,
              native_axis="protein", units="kcal/mol", mechanism="quantum_reference"),
    ArraySpec("mopac_dipole_debye", "mopac_direct", np.ndarray, 3, False,
              "(3,) unquantized binary64 molecular electric dipole xyz copied from mopac_properties.dipole; no claim of physical exactness. " + _MOPAC_DIRECT_SOURCE,
              native_axis="protein", units="Debye", mechanism="quantum_reference"),
    ArraySpec("mopac_dipole_point_charge_debye", "mopac_direct", np.ndarray, 3, False,
              "(3,) unquantized binary64 point-charge contribution xyz copied from mopac_properties.dipole_point_charge; no claim of physical exactness. " + _MOPAC_DIRECT_SOURCE,
              native_axis="protein", units="Debye", mechanism="quantum_reference"),
    ArraySpec("mopac_dipole_hybridization_debye", "mopac_direct", np.ndarray, 3, False,
              "(3,) unquantized binary64 one-centre hybridization contribution xyz copied from mopac_properties.dipole_hybridization; no claim of physical exactness. " + _MOPAC_DIRECT_SOURCE,
              native_axis="protein", units="Debye", mechanism="quantum_reference"),
    ArraySpec("mopac_bond_index", "mopac_direct", np.ndarray, None, False,
              "(N+1,) int32 zero-based CSC offsets from mopac_properties.bond_index; final value E sizes the directed bond arrays. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_csc_pointer", units="index", mechanism="quantum_reference"),
    ArraySpec("mopac_bond_atom", "mopac_direct", np.ndarray, None, False,
              "(E,) int32 zero-based CSC row/neighbor atom indices from mopac_properties.bond_atom, parallel to mopac_bond_order and directed density. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_csc_entry", units="index", mechanism="quantum_reference"),
    ArraySpec("mopac_bond_order", "mopac_direct", np.ndarray, None, False,
              "(E,) retained directed sparse-CSC values from mopac_properties.bond_order: diagonal valencies and strict >0.01 off-diagonal Wiberg indices; an omitted pair is not a measured zero. Values are unquantized API binary64, not a claim of physical exactness. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_csc_entry", units="dimensionless", mechanism="quantum_reference"),
    ArraySpec("mopac_ao_max_orbitals", "mopac_direct", np.ndarray, None, False,
              "(1,) int32 AO block width W from mopac_properties.ao_max_orbitals. " + _MOPAC_DIRECT_SOURCE,
              native_axis="protein", units="count", mechanism="quantum_reference"),
    ArraySpec("mopac_ao_orbitals_per_atom", "mopac_direct", np.ndarray, None, False,
              "(N,) int32 live AO width per atom from mozyme_state.iorbs, cross-checked exactly equal to mopac_properties.ao_orbitals so both populated API views are represented. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="count", mechanism="quantum_reference"),
    ArraySpec("mopac_atom_ao_density", "mopac_direct", np.ndarray, None, False,
              "(N,W,W) atom-centred s/p/d AO-basis Coulson population matrices from mopac_properties.atom_ao_density, converted from Fortran [W,W,N] to C order. These are reducible AO-basis matrices, not Cartesian rank-2/e3nn tensors; global-W cells beyond each atom's live width are structural N/A zero padding. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="electron", mechanism="quantum_reference"),
    ArraySpec("mopac_atomic_orbital_populations_full_precision", "mopac_direct", np.ndarray, None, False,
              "(N,W) unquantized binary64 diagonal of each mopac_properties.atom_ao_density block in MOPAC s/px/py/pz/d AO order. Diagonal-only components are not closed under rotation because coherences are omitted. Global-W cells beyond live width are structural N/A zero padding; the unsuffixed NPY is the F10.5 (N,9) compatibility projection. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="electron", mechanism="quantum_reference"),
    ArraySpec("mopac_bond_ao_density_directed", "mopac_direct", np.ndarray, None, False,
              "(E,W,W) Coulson AO-basis blocks parallel to the retained sparse CSC bond_atom/bond_order arrays, including atom-diagonal entries and retained off-diagonal forward/reverse entries. They are left/right AO-basis matrices, not Cartesian rank-2/e3nn tensors; cells beyond live endpoint widths are structural padding. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_csc_entry", units="electron", mechanism="quantum_reference"),
    ArraySpec("mopac_bond_density_pairs", "mopac_direct", np.ndarray, 2, False,
              "(U,2) int32 zero-based a<b atom pairs indexing the probe-compatible mopac_bond_ao_density view. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_unique_pair", units="index", mechanism="quantum_reference"),
    ArraySpec("mopac_bond_ao_density", "mopac_direct", np.ndarray, None, False,
              "(U,W,W) a<b Coulson AO-basis blocks selected from the retained sparse directed field; squared live-block norms reconstruct retained paired Wiberg indices. Missing sparse pairs have no row, and padded cells outside live endpoint widths are structural N/A zeros. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_unique_pair", units="electron", mechanism="quantum_reference"),
    ArraySpec("mopac_atom_electron_population", "mopac_direct", np.ndarray, None, False,
              "(N,) trace of each atom AO-density block. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="electron", mechanism="quantum_reference"),
    ArraySpec("mopac_atom_s_population", "mopac_direct", np.ndarray, None, False,
              "(N,) atomic s-shell Coulson population from AO diagonal entry 0. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="electron", mechanism="quantum_reference"),
    ArraySpec("mopac_atom_p_population", "mopac_direct", np.ndarray, None, False,
              "(N,) atomic p-shell Coulson population summed over px/py/pz; mathematical zero where no p AO exists. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="electron", mechanism="quantum_reference"),
    ArraySpec("mopac_atom_d_population", "mopac_direct", np.ndarray, None, False,
              "(N,) atomic d-shell Coulson population summed over MOPAC's five d AOs; mathematical zero where no d AO exists. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="electron", mechanism="quantum_reference"),
    ArraySpec("mopac_lewis_bond_count", "mopac_direct", np.ndarray, None, False,
              "(N,) int32 MOZYME Lewis-bond count from mozyme_state.nbonds. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="count", mechanism="quantum_reference"),
    ArraySpec("mopac_lewis_bond_atoms", "mopac_direct", np.ndarray, 9, False,
              "(N,9) int32 zero-based MOZYME Lewis-bond atoms from mozyme_state.ibonds; unused slots are explicit -1. " + _MOPAC_DIRECT_SOURCE,
              native_axis="atom", units="index", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_energy_levels", "mopac_direct", np.ndarray, None, False,
              "(O+V,) localized molecular-orbital energies from mopac_properties.lmo_energy, occupied first then virtual; these are not canonical HOMO/LUMO energies. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_lmo", units="eV", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_occupied_atom_counts", "mopac_direct", np.ndarray, None, False,
              "(O,) int32 occupied-LMO atom counts from mozyme_state.ncf. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_occupied_lmo", units="count", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_occupied_atoms", "mopac_direct", np.ndarray, None, False,
              "(sum ncf,) int32 packed zero-based occupied-LMO atom lists selected from mozyme_state.icocc using native libmopac offsets. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_occupied_lmo_atom_entry", units="index", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_occupied_coefficients", "mopac_direct", np.ndarray, None, False,
              "(sum live AO widths,) packed occupied-LMO AO coefficients selected from mozyme_state.cocc using native offsets and listed-atom AO widths. Orbital sign, ordering and localization gauge are not fixed across rotated/re-run SCF calculations, so rows have no closed e3nn/O(3) law. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_occupied_lmo_coefficient", units="dimensionless", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_virtual_atom_counts", "mopac_direct", np.ndarray, None, False,
              "(V,) int32 virtual-LMO atom counts from mozyme_state.nce; MOPAC notes that the virtual set is not re-localized. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_virtual_lmo", units="count", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_virtual_atoms", "mopac_direct", np.ndarray, None, False,
              "(sum nce,) int32 packed zero-based virtual-LMO atom lists selected from mozyme_state.icvir using native libmopac offsets. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_virtual_lmo_atom_entry", units="index", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_virtual_coefficients", "mopac_direct", np.ndarray, None, False,
              "(sum live AO widths,) packed virtual-LMO AO coefficients selected from mozyme_state.cvir using native offsets and listed-atom AO widths. Orbital sign, ordering and localization gauge are not fixed across rotated/re-run SCF calculations, so rows have no closed e3nn/O(3) law. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_virtual_lmo_coefficient", units="dimensionless", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_occupied_atom_offsets_native", "mopac_direct", np.ndarray, None, False,
              "(O,) int32 exact mopac_properties.lmo_occupied_atom_offset values into native icocc storage; not offsets into the compact view. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_occupied_lmo", units="index", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_virtual_atom_offsets_native", "mopac_direct", np.ndarray, None, False,
              "(V,) int32 exact mopac_properties.lmo_virtual_atom_offset values into native icvir storage; not compact-view offsets. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_virtual_lmo", units="index", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_occupied_coefficient_offsets_native", "mopac_direct", np.ndarray, None, False,
              "(O,) int32 exact mopac_properties.lmo_occupied_coefficient_offset values into native cocc storage. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_occupied_lmo", units="index", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_virtual_coefficient_offsets_native", "mopac_direct", np.ndarray, None, False,
              "(V,) int32 exact mopac_properties.lmo_virtual_coefficient_offset values into native cvir storage. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_virtual_lmo", units="index", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_occupied_atom_storage_native", "mopac_direct", np.ndarray, None, False,
              "(icocc_dim,) int32 exact one-based mozyme_state.icocc native storage; only offset+count slices are live scientific lists and remaining allocated capacity is retained restart state. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_icocc_storage", units="MOPAC_one_based_index", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_virtual_atom_storage_native", "mopac_direct", np.ndarray, None, False,
              "(icvir_dim,) int32 exact one-based mozyme_state.icvir native storage; only offset+count slices are live and remaining allocated capacity is retained restart state. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_icvir_storage", units="MOPAC_one_based_index", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_occupied_coefficient_storage_native", "mopac_direct", np.ndarray, None, False,
              "(cocc_dim,) exact mozyme_state.cocc native occupied-LMO coefficient storage; offsets identify live slices. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_cocc_storage", units="dimensionless", mechanism="quantum_reference"),
    ArraySpec("mopac_lmo_virtual_coefficient_storage_native", "mopac_direct", np.ndarray, None, False,
              "(cvir_dim,) exact mozyme_state.cvir native virtual-LMO coefficient storage; offsets identify live slices. " + _MOPAC_DIRECT_SOURCE,
              native_axis="mopac_cvir_storage", units="dimensionless", mechanism="quantum_reference"),
    ArraySpec("mopac_mozyme_state_dimensions", "mopac_direct", np.ndarray, 7, False,
              "(7,) int32 exact mozyme_state scalar fields [numat,noccupied,nvirtual,icocc_dim,icvir_dim,cocc_dim,cvir_dim]. " + _MOPAC_DIRECT_SOURCE,
              native_axis="protein", units="count", mechanism="quantum_reference"),

    # ── MOPAC Coulomb (MopacCoulombResult.cpp) ───────────────────
    ArraySpec("mopac_coulomb_efg",           "mopac_coulomb", ShieldingTensor, 9,  False, "MOPAC-charge bare total EFG in project Full9; finite analytic rows have structural-zero T0/T1, but the unmasked post-projection elementwise nonfinite sanitizer can break T0." + _MOPAC_COULOMB_SOURCE,
              irreps=_SHIELD_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_E",             "mopac_coulomb", VectorField,     3,  False, "MOPAC-charge all-pairs total E-field." + _MOPAC_COULOMB_SOURCE,
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_E_backbone",    "mopac_coulomb", VectorField,     3,  False, "MOPAC-charge field from non-aromatic cached backbone N/CA/C/O/H/HA/CB sources." + _MOPAC_COULOMB_SOURCE,
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_E_sidechain",   "mopac_coulomb", VectorField,     3,  False, "MOPAC-charge field from sources neither aromatic nor cached backbone." + _MOPAC_COULOMB_SOURCE,
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_E_aromatic",    "mopac_coulomb", VectorField,     3,  False, "MOPAC-charge field from aromatic-ring-member sources." + _MOPAC_COULOMB_SOURCE,
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_efg_backbone",  "mopac_coulomb", EFGTensor,       5,  False, "MOPAC-charge EFG native T2 from non-aromatic cached backbone N/CA/C/O/H/HA/CB sources." + _MOPAC_COULOMB_SOURCE,
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("mopac_coulomb_efg_sidechain", "mopac_coulomb", EFGTensor,       5,  False, "MOPAC-charge EFG native T2 from sources neither aromatic nor cached backbone." + _MOPAC_COULOMB_SOURCE,
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("mopac_coulomb_efg_aromatic",  "mopac_coulomb", EFGTensor,       5,  False, "MOPAC-charge EFG native T2 from aromatic-ring-member sources." + _MOPAC_COULOMB_SOURCE,
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("mopac_coulomb_scalars",       "mopac_coulomb", CoulombScalars,  4,  False, "MOPAC-charge scalars [|E_total|, parent-to-H E_total projection (NaN when unavailable), signed E_backbone projection onto E_total direction (historical *_frac name; zero at near-zero |E_total|), |E_aromatic|]." + _MOPAC_COULOMB_SOURCE,
              units="V/A", mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_aromatic_E_proj", "mopac_coulomb", np.ndarray, None, False, "Signed parent-to-H bond-axis projection of the MOPAC-charge aromatic E-field." + _MOPAC_COULOMB_SOURCE,
              irreps="0e", units="V/A", mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_aromatic_n_src", "mopac_coulomb", np.ndarray, None, False, "Count of contributing aromatic sources outside the cached backbone set in the MOPAC-charge all-pairs sum." + _MOPAC_COULOMB_SOURCE,
              irreps="0e", units="count", mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_E_clamp_mask", "mopac_coulomb", np.ndarray, None, False, "Per-atom uint8 mask identifying joint MOPAC Coulomb E-vector rescaling by the configured magnitude clamp.",
              irreps="0e", units="dimensionless", mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_E_clamp_scale", "mopac_coulomb", np.ndarray, None, False, "Per-atom scale applied jointly to the four MOPAC Coulomb E vectors; 1.0 when unclamped.",
              irreps="0e", units="dimensionless", mechanism="electrostatic_efg"),

    # ── MOPAC McConnell BO projection (MopacMcConnellResult.cpp) ───
    ArraySpec("mopac_mc_co_sum", "mopac_mcconnell", np.ndarray, None, False, "T0 scalar of the summed MOPAC Wiberg-weighted PeptideCO response produced by McConnellResult.",
              irreps="0e", units="Angstrom^-3", mechanism="bond_anisotropy", scaling_contract=_MCCONNELL_SCALING),
    ArraySpec("mopac_mc_cn_sum", "mopac_mcconnell", np.ndarray, None, False, "T0 scalar of the summed MOPAC Wiberg-weighted PeptideCN response produced by McConnellResult.",
              irreps="0e", units="Angstrom^-3", mechanism="bond_anisotropy", scaling_contract=_MCCONNELL_SCALING),
    ArraySpec("mopac_mc_sidechain_sum", "mopac_mcconnell", np.ndarray, None, False, "T0 scalar summed across the MOPAC Wiberg-weighted sidechain response categories produced by McConnellResult.",
              irreps="0e", units="Angstrom^-3", mechanism="bond_anisotropy", scaling_contract=_MCCONNELL_SCALING),
    ArraySpec("mopac_mc_aromatic_sum", "mopac_mcconnell", np.ndarray, None, False, "MOPAC Wiberg-weighted aromatic scalar channel, structurally zeroed by McConnellResult to avoid ring-current double counting.",
              irreps="0e", units="Angstrom^-3", mechanism="bond_anisotropy", scaling_contract=_MCCONNELL_SCALING),
    ArraySpec("mopac_mc_co_nearest", "mopac_mcconnell", np.ndarray, None, False, "MOPAC Wiberg-weighted scalar response of the geometrically nearest accepted PeptideCO source; zero when that source's post-floor bond order is zero.",
              irreps="0e", units="Angstrom^-3", mechanism="bond_anisotropy", scaling_contract=_MCCONNELL_SCALING),
    ArraySpec("mopac_mc_nearest_co_dist", "mopac_mcconnell", np.ndarray, None, False, "Distance to the geometrically nearest accepted PeptideCO source, independent of its MOPAC Wiberg bond order.",
              irreps="0e", units="Å", mechanism="bond_anisotropy"),
    ArraySpec("mopac_mc_nearest_cn_dist", "mopac_mcconnell", np.ndarray, None, False, "Distance to the geometrically nearest accepted PeptideCN source, independent of its MOPAC Wiberg bond order.",
              irreps="0e", units="Å", mechanism="bond_anisotropy"),
    ArraySpec("mopac_mc_nearest_co_T2", "mopac_mcconnell", ShieldingTensor, 9, False, "Full project-native SphericalTensor read-back of the MOPAC Wiberg-weighted response for the geometrically nearest accepted PeptideCO source.",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy", scaling_contract=_MCCONNELL_SCALING),
    ArraySpec("mopac_mc_nearest_cn_T2", "mopac_mcconnell", ShieldingTensor, 9, False, "Full project-native SphericalTensor read-back of the MOPAC Wiberg-weighted response for the geometrically nearest accepted PeptideCN source.",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy", scaling_contract=_MCCONNELL_SCALING),
    ArraySpec("mopac_mc_backbone_total", "mopac_mcconnell", ShieldingTensor, 9, False, "Full project-native SphericalTensor read-back of the summed MOPAC Wiberg-weighted backbone response.",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy", scaling_contract=_MCCONNELL_SCALING),
    ArraySpec("mopac_mc_sidechain_total", "mopac_mcconnell", ShieldingTensor, 9, False, "Full project-native SphericalTensor read-back of the summed MOPAC Wiberg-weighted sidechain response.",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy", scaling_contract=_MCCONNELL_SCALING),
    ArraySpec("mopac_mc_aromatic_total", "mopac_mcconnell", ShieldingTensor, 9, False, "Full project-native SphericalTensor read-back of the MOPAC Wiberg-weighted aromatic response, structurally zeroed by McConnellResult.",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy", scaling_contract=_MCCONNELL_SCALING),
    ArraySpec("mopac_mc_shielding", "mopac_mcconnell", ShieldingTensor, 9, False, "Full project-native SphericalTensor sum of all MOPAC Wiberg-weighted McConnell categories; unscaled geometry response, not ppm.",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", sign_convention="producer response sign; physical shielding sign depends on the downstream susceptibility scale", tensor_rank=2, mechanism="bond_anisotropy", scaling_contract=_MCCONNELL_SCALING),

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
    ArraySpec("delta_shielding",       "delta", ShieldingTensor,       9,    False, "WT-minus-mutant total ORCA shielding delta on the dense WT atom axis",
              native_axis="atom", irreps=_SHIELD_IRREPS, units="ppm",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="mutation_delta"),
    ArraySpec("delta_scalars",         "delta", DeltaScalars,          6,    False, "Dense WT-atom rows [matched, WT-mutant shielding T0, nearest removed-ring distance, configured-ChargeSource charge delta, legacy F15.6 MOPAC charge delta, typed WT-to-mutant match distance]; matched rows with no removed ring use 99 Å; MOPAC-delta zero is ambiguous between unavailable MOPAC and a true zero",
              native_axis="atom", units="mixed_mask_ppm_A_e_e_A", mechanism="mutation_delta"),
    ArraySpec("mutation_atom_map",     "delta", np.ndarray,             None, False, "Dense WT-atom rows containing the corresponding mutant atom index in the separately emitted ala/ calcset; -1 means unmatched",
              native_axis="atom", irreps="0e", units="index", mechanism="mutation_delta"),
    ArraySpec("delta_graph",           "delta", np.ndarray,             5,    False, "Graph deltas [matched, has_graph_delta, delta_graph_dist_ring, delta_bfs_decay, delta_is_conjugated]",
              native_axis="atom", mechanism="mutation_delta"),
    ArraySpec("delta_apbs",            "delta", DeltaAPBS,             12,   False, "Optional APBS delta_E(3) + legacy full-9 EFG envelope; only columns 7:12 are physical EFG T2, and the whole file is absent when the mutation calculation has no paired APBS delta",
              native_axis="atom", mechanism="mutation_delta",
              tensor_basis=_T2_BASIS,
              tensor_component_order="delta_E_x,delta_E_y,delta_E_z,T0_compat_zero,T1_x_compat_zero,T1_y_compat_zero,T1_z_compat_zero,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2",
              tensor_frame="conformation_cartesian_xyz",
              structural_zero_components=_EFG_STRUCTURAL_ZEROS,
              e3nn_export="raw project tensor; call DeltaAPBS.delta_efg_t2.to_e3nn() before using e3nn Irreps"),
    ArraySpec("delta_ring_proximity",  "delta", DeltaRingProximity,    None, False, "Optional dense WT-atom rows with one six-column block per removed WT ring [distance_A, signed z_A, rho_A, theta_rad, (3cos^2(theta)-1)/r^3_A^-3, exp(-distance/L)]; the whole file is absent without a mutation site that removes at least one WT ring",
              native_axis="atom", units="mixed_A_radians_A^-3_dimensionless", mechanism="mutation_delta"),
    # DFT shielding component decomposition: WT side, mut side, deltas;
    # diamagnetic and paramagnetic. sigma_total = sigma_dia + sigma_para;
    # the existing delta_shielding satisfies that identity at ORCA's
    # output precision (~1e-3 ppm). Stratifies mutation shifts by
    # physical mechanism.
    ArraySpec("wt_shielding_diamagnetic",     "delta", ShieldingTensor, 9, False, "WT diamagnetic shielding (matched, by WT atom row)",
              native_axis="atom", irreps=_SHIELD_IRREPS, units="ppm",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="mutation_delta"),
    ArraySpec("wt_shielding_paramagnetic",    "delta", ShieldingTensor, 9, False, "WT paramagnetic shielding (matched, by WT atom row)",
              native_axis="atom", irreps=_SHIELD_IRREPS, units="ppm",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="mutation_delta"),
    ArraySpec("mut_shielding_diamagnetic",    "delta", ShieldingTensor, 9, False, "mut diamagnetic shielding (matched, by WT atom row)",
              native_axis="atom", irreps=_SHIELD_IRREPS, units="ppm",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="mutation_delta"),
    ArraySpec("mut_shielding_paramagnetic",   "delta", ShieldingTensor, 9, False, "mut paramagnetic shielding (matched, by WT atom row)",
              native_axis="atom", irreps=_SHIELD_IRREPS, units="ppm",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="mutation_delta"),
    ArraySpec("delta_shielding_diamagnetic",  "delta", ShieldingTensor, 9, False, "WT - mut diamagnetic shielding delta",
              native_axis="atom", irreps=_SHIELD_IRREPS, units="ppm",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="mutation_delta"),
    ArraySpec("delta_shielding_paramagnetic", "delta", ShieldingTensor, 9, False, "WT - mut paramagnetic shielding delta",
              native_axis="atom", irreps=_SHIELD_IRREPS, units="ppm",
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
    ArraySpec("aimnet2_charges",             "aimnet2", AIMNet2Charges,            None, True,  "Per-atom `charges` output of the loaded AIMNet2 TorchScript model (intended Hirshfeld charges for a conforming model); float64 NPY widens the model result and does not imply binary64 model evaluation",
              units="e", mechanism="charges"),
    ArraySpec("aimnet2_aim",                 "aimnet2", AIMNet2AimEmbedding,       256,  True,  "Raw float32 (N,256) `aim` output of the loaded AIMNet2 TorchScript model; 256 opaque learned channels, not Cartesian components or producer-certified e3nn irreps",
              mechanism="charges"),
    ArraySpec("aimnet2_efg",                 "aimnet2", EFGTensor,                 5,    True,  "AIMNet2 Coulomb EFG total (T2 only)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("aimnet2_efg_aromatic",        "aimnet2", EFGTensor,                 5,    True,  "AIMNet2-charge EFG T2 from aromatic-ring members, with aromatic precedence over backbone classification",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("aimnet2_efg_backbone",        "aimnet2", EFGTensor,                 5,    True,  "AIMNet2-charge EFG T2 from non-aromatic cached backbone N/CA/C/O/H/HA/CB sources (CB included)",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("aimnet2_efg_sidechain",       "aimnet2", EFGTensor,                 5,    True,  "AIMNet2-charge EFG T2 from sources neither aromatic nor cached backbone",
              units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg", **_T2_TENSOR_METADATA),
    ArraySpec("aimnet2_E",                   "aimnet2", VectorField,               3,    True,  "AIMNet2 charge-derived total E-field",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("aimnet2_E_backbone",          "aimnet2", VectorField,               3,    True,  "AIMNet2-charge field from non-aromatic cached backbone N/CA/C/O/H/HA/CB sources (CB included)",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("aimnet2_E_sidechain",         "aimnet2", VectorField,               3,    True,  "AIMNet2-charge field from sources neither aromatic nor cached backbone",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("aimnet2_E_aromatic",          "aimnet2", VectorField,               3,    True,  "AIMNet2-charge field from aromatic-ring members, with aromatic precedence over backbone classification",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("aimnet2_energy_mlp",          "aimnet2", np.ndarray,                None, True, "Per-atom energy after replaying the loaded model's energy_mlp head; float64 NPY widens the Torch result",
              units="eV", mechanism="charges"),
    ArraySpec("aimnet2_energy_shifted_local", "aimnet2", np.ndarray,               None, True, "Per-atom local energy after replaying energy_mlp then atomic_shift; float64 NPY widens the Torch result",
              units="eV", mechanism="charges"),
    ArraySpec("aimnet2_energy_terms",        "aimnet2", np.ndarray,                6,    True, "One protein row [sum after atomic_sum, lrcoulomb increment, dftd3 increment, final energy, configured model-conditioning net charge, neutral-conditioning 0/1 flag]",
              native_axis="protein", units="mixed:eV[0:4],e[4],dimensionless[5]", mechanism="charges"),
    ArraySpec("aimnet2_d3_e_disp_atom",      "aimnet2", np.ndarray,                None, True, "Per-atom D3 dispersion-energy increment recomputed from the loaded model's private dftd3 state and long-range neighbour table",
              units="eV", mechanism="charges"),
    ArraySpec("aimnet2_d3_cn",               "aimnet2", np.ndarray,                None, True, "Per-atom D3 coordination number recomputed from the loaded model's private dftd3 state and long-range neighbour table",
              units="dimensionless", mechanism="charges"),
    ArraySpec("aimnet2_d3_c6_stats",         "aimnet2", np.ndarray,                3,    True, "Per-atom D3 C6 [sum,mean,max] over valid long-range neighbour slots; exactly [0,0,0] when no slot is valid",
              units="Hartree*bohr^6", mechanism="charges"),
    ArraySpec("aimnet2_aim_projection",      "aimnet2", np.ndarray,                32,   True, "Float32 fixed element-conditioned 32-d linear projection of raw aim, computed once in Compute; basis splitmix64_0xA17E20260708_achlioptas_32x256_element_HCNOS; channels remain opaque learned features",
              units="dimensionless", mechanism="charges"),

    # ── AIMNet2 charge-response gradient (AIMNet2ChargeResponseGradientResult.cpp)
    # Always-on after the --aimnet2 model is loaded (per the
    # 2026-05-09 promotion of Amendment 2026-05-08(b) from a test flag
    # to standard non-trajectory pipeline; trajectory mode unchanged).
    # Vector is dL/d(r_i) where L = sum_j q_j^2 over non-sentinel atoms
    # (for a conforming model, charge conservation makes the sum(q) gradient
    # near-zero, so the L2 of charges is the cheapest single-pass objective
    # with non-trivial gradient). Scalar is the L2 norm. required=False
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
    # Four structured-NPY projections of the Protein topology. The writer
    # attempts them for every protein; empty semantic/ring content is carried
    # by the structured rows/sentinels rather than suppressing the sidecar.
    # ────────────────────────────────────────────────────────────────
    ArraySpec("residues",         "topology", np.ndarray, None, True,  "Per-residue record: residue_index, chain_id, residue_number, insertion_code, residue_type, amber/iupac/one-letter names, protonation_variant_index, terminal_state, prev/next links, atom_count, is_proline/aromatic/titratable, has_amide_h, is_xpro_context",
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
    if spec.wrapper is PerRingTypeT2:
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


def _set_units(stems, units):
    """Set an explicit unit/layout label on existing catalog entries."""
    for stem in stems:
        if stem not in CATALOG:
            raise RuntimeError(f"unit contract references unknown NPY stem {stem!r}")
        CATALOG[stem] = replace(CATALOG[stem], units=units)


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
_LOCAL_BACKBONE_ANGLES = (
    "tau_N_CA_C", "angle_N_CA_CB", "angle_CB_CA_C",
    "angle_Cprev_N_CA", "angle_CA_C_Nnext",
)
_set_contract(
    _LOCAL_BACKBONE_ANGLES, coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant valence angle in radians; translation invariant"
    ),
    validity=(
        "the same-stem _valid.npy uint8 array gates independently computed "
        "values; invalid/non-applicable rows are NaN"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    tuple(f"{stem}_valid" for stem in _LOCAL_BACKBONE_ANGLES),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant uint8 availability for the corresponding "
        "residue-local valence angle"
    ),
    validity="1 iff the corresponding angle is finite; 0 gates its NaN",
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
    ("mc_nearest_co_bond_index", "mc_nearest_cn_bond_index"),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)- and translation-invariant topology bond identity selected "
        "solely by accepted source-midpoint geometry"
    ),
    validity=(
        "-1 means no accepted source of that category; every nonnegative value "
        "indexes bonds.npy and the same row of mopac_topology_bond_orders_full.npy"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("hbond_nearest_dir",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=_POLAR_VECTOR,
    validity="hbond_flags.npy column 0 identifies an accepted nearest source; otherwise zero")
_set_contract(
    ("hbond_scalars",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "four exact O(3)-invariant scalars for a fixed typed accepted-pair set; "
        "translation invariant"
    ),
    validity=(
        "hbond_flags.npy column 0 gates nearest-distance and inverse-cube "
        "legacy zeros; count/angular-sum zero is a physical empty/zero sum"
    ),
    irreps="0e", parity="even", tensor_rank=0)
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
    ("hydration_shell",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant scalar row for a jointly transformed protein, "
        "waters and ions; translation invariant"
    ),
    validity=(
        "whole NPY is absent when explicit solvent is unavailable; no "
        "first-shell waters gives cols0:2=0; no finite-dipole sample gives "
        "col1=0; no ion inside hydration_ion_cutoff gives cols2:4=[+inf,0]"
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
    coordinate_frame=_CARTESIAN_FRAME, transformation=(
        "polar_vector: v'=R v"
    ),
    validity=(
        "force-field Coulomb producer replaces non-finite values with zero; "
        "no sanitizer mask (legacy unavailable/zero ambiguity); all four "
        "vectors are jointly scaled when total-E exceeds the configured clamp, "
        "with provenance only in GeometryChoice"
    ), irreps="1o", parity="odd", tensor_rank=1)
_set_contract(
    ("mopac_coulomb_E", "mopac_coulomb_E_backbone",
     "mopac_coulomb_E_sidechain", "mopac_coulomb_E_aromatic"),
    coordinate_frame=_CARTESIAN_FRAME, transformation=(
        "polar_vector: v'=R v"
    ),
    validity=(
        "MOPAC-charge Coulomb producer replaces non-finite values with zero; "
        "no sanitizer mask (physical-zero/sanitized-zero ambiguity); all four "
        "vectors are jointly scaled when total-E exceeds the configured clamp, "
        "with provenance in mopac_coulomb_E_clamp_mask.npy and "
        "mopac_coulomb_E_clamp_scale.npy; all arrays are absent when the "
        "prerequisite MopacResult is absent"
    ), irreps="1o", parity="odd", tensor_rank=1)
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
    ("coulomb_scalars",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "O(3)-invariant magnitudes and signed vector projections"
    ),
    validity=(
        "col1 is NaN without a valid parent-to-H direction; col2 is a signed "
        "V/A projection (not a fraction) and is zero at near-zero total E; "
        "field sanitizer/clamp has no NPY mask"
    ),
    irreps="", parity="even", tensor_rank=0)
_set_contract(
    ("mopac_coulomb_scalars",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "O(3)-invariant magnitudes and signed projections"
    ),
    validity=(
        "col1 is NaN without a valid parent-to-H direction; col2 is a signed "
        "V/A projection (not a fraction) and is zero at near-zero total E; "
        "absent with MopacResult; consult mopac_coulomb_E_clamp_mask.npy and "
        "mopac_coulomb_E_clamp_scale.npy for field rescaling; the legacy "
        "nonfinite sanitizer remains unmasked"
    ),
    irreps="", parity="even", tensor_rank=0)

_set_contract(
    ("coulomb_aromatic_E_proj",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "O(3)-invariant scalar parent-bond projection"
    ),
    validity=(
        "NaN for non-H or parentless atoms; otherwise physical signed V/A "
        "projection"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("mopac_coulomb_aromatic_E_proj",),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation="exact O(3)-invariant scalar parent-bond projection",
    validity=(
        "NaN for non-H or parentless atoms; otherwise the signed projection "
        "of the post-clamp aromatic MOPAC Coulomb field; absent with MopacResult"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("mopac_coulomb_aromatic_n_src",),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant contributing-source count for fixed topology, "
        "charges, filters and charge floor"
    ),
    validity=(
        "zero is a physical empty count; includes only aromatic sources outside "
        "the cached backbone set that survive the source filters and charge floor; "
        "absent with MopacResult"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("mopac_coulomb_E_clamp_mask", "mopac_coulomb_E_clamp_scale"),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant clamp diagnostic for the analytic all-pairs "
        "MOPAC Coulomb field; translation invariant"
    ),
    validity=(
        "mask is uint8 1 iff all four E vectors were jointly rescaled; scale "
        "is clamp/original_total_magnitude when masked and exactly 1 otherwise; "
        "absent with MopacResult"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("eeq_coulomb_scalars",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant magnitudes and signed vector projections for "
        "the serialized EEQ-charge fields"
    ),
    validity=(
        "col1 is NaN without a valid parent-to-H direction; col2 is a signed "
        "V/A projection (not a fraction) and is zero at near-zero total E; "
        "the whole result is absent on non-finite fields"
    ),
    irreps="0e", parity="even", tensor_rank=0)
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
    ("water_shell_counts",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant first/second-shell water counts derived from "
        "protein-water distances; translation invariant"
    ),
    validity=(
        "columns are n_first,n_second; zero is the physical empty-shell "
        "count; whole WaterFieldResult is absent without solvent"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("water_efield_clamp_mask", "water_efield_first_clamp_mask"),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant int8 threshold diagnostic derived from the "
        "corresponding water E-field magnitude; translation invariant"
    ),
    validity=(
        "0/1; 1 iff the corresponding total/first-shell E row was clamped; "
        "whole WaterFieldResult is absent without solvent"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("water_efield_clamp_scale", "water_efield_first_clamp_scale"),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant scalar scale derived from the corresponding "
        "water E-field magnitude and configured threshold; translation invariant"
    ),
    validity=(
        "finite in (0,1] when the result exists; 1 means unclamped; whole "
        "WaterFieldResult is absent without solvent"
    ),
    irreps="0e", parity="even", tensor_rank=0)
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
        "dense total uses physical zero when no accepted modeled source and "
        "excludes diagnostic TrpPerimeter rings; sparse fields use the union "
        "ring_contributions.npy axis, where calculator-default zero is "
        "ambiguous with a physical zero and no per-calculator mask is emitted"
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
for _stem in _MCCONNELL_FULL9:
    CATALOG[_stem] = replace(
        CATALOG[_stem],
        description=(
            CATALOG[_stem].description +
            "; unscaled unit-susceptibility geometry response, not ppm"
        ),
        scaling_contract=_MCCONNELL_SCALING,
    )
CATALOG["sidechain_co_scalar_audit"] = replace(
    CATALOG["sidechain_co_scalar_audit"],
    description=(
        CATALOG["sidechain_co_scalar_audit"].description +
        "; cols0:2 are norms of unscaled unit-susceptibility responses"
    ),
    scaling_contract=(
        _MCCONNELL_SCALING +
        "; cols0:2 scale by the magnitude of the downstream response scale; "
        "count and distance columns are unscaled"
    ),
)
for _stem in ("mc_aromatic_zeroed_fixed", "mc_aromatic_zeroed_bo"):
    CATALOG[_stem] = replace(
        CATALOG[_stem],
        structural_zero_components=(
            "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
        ),
    )

_MOPAC_MCCONNELL_FULL9 = (
    "mopac_mc_nearest_co_T2", "mopac_mc_nearest_cn_T2",
    "mopac_mc_backbone_total", "mopac_mc_sidechain_total",
    "mopac_mc_aromatic_total", "mopac_mc_shielding",
)
_set_contract(
    _MOPAC_MCCONNELL_FULL9,
    coordinate_frame=_CARTESIAN_FRAME,
    transformation=_EVEN_RANK2,
    validity=(
        "whole family is absent unless MopacMcConnellResult attaches; nearest "
        "CO/CN tensors are NaN when the corresponding mopac_mc_nearest_*_dist "
        "is NO_DATA_SENTINEL; a valid nearest tensor is physical zero when the "
        "geometry-selected bond's post-floor BO is zero; aggregate zero means "
        "an empty/zero BO sum"
    ))
for _stem in _MOPAC_MCCONNELL_FULL9:
    CATALOG[_stem] = replace(
        CATALOG[_stem], scaling_contract=_MCCONNELL_SCALING)
CATALOG["mopac_mc_aromatic_total"] = replace(
    CATALOG["mopac_mc_aromatic_total"],
    structural_zero_components=(
        "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
    ))

_set_contract(
    ("mopac_mc_co_sum", "mopac_mc_cn_sum", "mopac_mc_sidechain_sum",
     "mopac_mc_aromatic_sum", "mopac_mc_co_nearest"),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant scalar from the MOPAC Wiberg-weighted "
        "McConnell response; translation invariant"
    ),
    validity=(
        "whole family is absent unless MopacMcConnellResult attaches; sums use "
        "physical zero for empty/zero BO channels; co_nearest is gated by "
        "mopac_mc_nearest_co_dist.npy and is physical zero when the geometry-"
        "selected bond's post-floor BO is zero; aromatic_sum is structurally zero"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("mopac_mc_nearest_co_dist", "mopac_mc_nearest_cn_dist"),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant source-midpoint distance; translation invariant"
    ),
    validity=(
        "NO_DATA_SENTINEL means no geometrically accepted source of that category; "
        "selection is independent of MOPAC bond order and matches the corresponding "
        "mc_nearest_*_bond_index.npy; whole family is absent unless "
        "MopacMcConnellResult attaches"
    ), irreps="0e", parity="even", tensor_rank=0)

for _stem in (
    "mc_peptide_co_bo", "mc_peptide_cn_bo", "mc_backbone_other_bo",
    "mc_sidechain_co_bo", "mc_sidechain_other_bo", "mc_disulfide_bo",
    "mc_backbone_xh_bo", "mc_sidechain_xh_bo", "mc_s_h_bo",
    "sidechain_co_bo_T2",
):
    CATALOG[_stem] = replace(
        CATALOG[_stem],
        transformation=_EVEN_RANK2,
        irreps=_SHIELD_IRREPS,
    )

_set_contract(
    ("coulomb_efg",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "even_rank2: T'=R T R^T"
    ),
    validity=(
        "symmetric-traceless source: T0/T1 are structural zeros"
    ), irreps=_SHIELD_IRREPS)
_set_contract(
    ("mopac_coulomb_efg",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "even_rank2: T'=R T R^T"
    ),
    validity=(
        "symmetric-traceless source: T0/T1 are structural zeros; absent when "
        "the prerequisite MopacResult is absent"
    ), irreps=_SHIELD_IRREPS)
_set_contract(
    ("eeq_coulomb_efg",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=_EVEN_RANK2,
    validity="T0/T1 structural zeros; whole result absent on non-finite output")

_EFG_T2 = (
    "coulomb_efg_t2", "coulomb_efg_backbone", "coulomb_efg_sidechain",
    "coulomb_efg_aromatic",
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
        transformation=(
            "even_rank2_native_T2: reconstruct Cartesian T, apply "
            "T'=R T R^T, then decompose in project-native T2 basis"
        ),
        validity=(
            "finite symmetric-traceless analytic source; only native T2 is "
            "serialized"
        ),
        structural_zero_components=(
            "T0,T1_x,T1_y,T1_z"
        ),
        e3nn_export=(
            "finite analytic rows only: raw project-native T2; convert "
            "explicitly before e3nn use. Exceptional sanitized rows are not "
            "a closed O(3) representation and have no validity mask"
        ),
    )
CATALOG["coulomb_efg"] = replace(
    CATALOG["coulomb_efg"],
    structural_zero_components="T0,T1_x,T1_y,T1_z",
    irreps=_SHIELD_IRREPS,
)
for _stem in (
    "mopac_coulomb_efg_backbone", "mopac_coulomb_efg_sidechain",
    "mopac_coulomb_efg_aromatic",
):
    CATALOG[_stem] = replace(
        CATALOG[_stem],
        transformation=(
            "even_rank2_native_T2: reconstruct Cartesian T, apply "
            "T'=R T R^T, then decompose in project-native T2 basis"
        ),
        validity=(
            "finite symmetric-traceless analytic source; only native T2 is "
            "serialized; absent when the prerequisite MopacResult is absent"
        ),
        structural_zero_components=(
            "T0,T1_x,T1_y,T1_z"
        ),
        e3nn_export=(
            "finite analytic rows only: raw project-native T2; convert "
            "explicitly before e3nn use. Exceptional sanitized rows are not "
            "a closed O(3) representation and have no validity mask"
        ),
    )
CATALOG["mopac_coulomb_efg"] = replace(
    CATALOG["mopac_coulomb_efg"],
    structural_zero_components=(
        "T0,T1_x,T1_y,T1_z"
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

# AIMNet2 loads the caller-selected TorchScript file.  The producer checks its
# runtime interface, dimensions and finiteness, but does not prove that an
# arbitrary compatible model is O(3)-equivariant.  Preserve the conditional
# physical laws below without upgrading them to unconditional rerun promises.
_AIMNET_MODEL_SCALAR = (
    "O(3)- and translation-invariant scalar"
)
_AIMNET_MAIN_VALIDITY = (
    "finite-validated output of the single AIMNet2Result Compute; current CLI "
    "extraction requires that result and aborts on model/interface/non-finite "
    "failure rather than fabricating a row"
)
_set_contract(
    ("aimnet2_charges", "aimnet2_energy_mlp",
     "aimnet2_energy_shifted_local", "aimnet2_d3_e_disp_atom",
     "aimnet2_d3_cn"),
    coordinate_frame="loaded_aimnet2_model_scalar_channel",
    transformation=_AIMNET_MODEL_SCALAR,
    validity=(
        _AIMNET_MAIN_VALIDITY + "; zero is a finite model/calculator value, "
        "not an absence sentinel; float64 NPY storage widens a Torch result "
        "and does not add numerical precision"
    ),
    irreps="", parity="even", tensor_rank=0)
_set_contract(
    ("aimnet2_energy_terms",),
    coordinate_frame="loaded_aimnet2_model_and_configuration_scalars",
    transformation=(
        "mixed scalar row: cols0:4 have the loaded-model conditional scalar "
        "law; col4 is the configured integer net-charge conditioning value "
        "and col5 its 0/1 neutral-conditioning flag"
    ),
    validity=(
        _AIMNET_MAIN_VALIDITY + "; cols0:4 are finite-validated, col4 may "
        "legitimately be zero, and col5 is exactly 0 or 1"
    ),
    irreps="", parity="even", tensor_rank=0)
_set_contract(
    ("aimnet2_d3_c6_stats",),
    coordinate_frame="loaded_aimnet2_d3_scalar_channels",
    transformation=_AIMNET_MODEL_SCALAR,
    validity=(
        _AIMNET_MAIN_VALIDITY + "; columns are [sum,mean,max] over valid "
        "long-range neighbour slots and [0,0,0] is the explicit structural "
        "empty-neighbour result; float64 NPY storage widens Torch values"
    ),
    irreps="", parity="even", tensor_rank=0)
_set_contract(
    ("aimnet2_aim", "aimnet2_aim_projection"),
    coordinate_frame="loaded_aimnet2_learned_channel_basis",
    transformation=(
        "opaque learned-channel vector with no producer-guaranteed Cartesian, "
        "O(3), parity or e3nn transformation law; do not rotate or reinterpret "
        "the channel axis. The 32-channel array is only a fixed element-"
        "conditioned linear projection of the same opaque 256 channels"
    ),
    validity=(
        _AIMNET_MAIN_VALIDITY + "; every emitted channel is finite float32; "
        "there is no per-channel availability mask and numeric zero is not "
        "an absence marker"
    ),
    irreps="", parity="even", tensor_rank=0, e3nn_export="")

_AIMNET_CONDITIONAL_VECTOR = (
    "conditional on the same emitted scalar AIMNet2 charge vector and source "
    "membership, the analytic Coulomb field is polar (v'=R v); a full "
    "producer rerun reevaluates an arbitrary caller-selected float32 .jpt, "
    "whose scalar invariance is not certified, so there is no unconditional "
    "rotated-rerun O(3) law"
)
_set_contract(
    ("aimnet2_E", "aimnet2_E_backbone", "aimnet2_E_sidechain",
     "aimnet2_E_aromatic"),
    coordinate_frame=_CARTESIAN_FRAME,
    transformation="polar_vector: v'=R v",
    validity=(
        "finite analytic V/A vectors; exact zero can mean an empty accepted "
        "source partition or cancellation and has no separate mask; any "
        "non-finite E/EFG rejects the whole AIMNet2Result"
    ),
    irreps="1o", parity="odd", tensor_rank=1)
for _stem in (
    "aimnet2_efg", "aimnet2_efg_aromatic",
    "aimnet2_efg_backbone", "aimnet2_efg_sidechain",
):
    CATALOG[_stem] = replace(
        CATALOG[_stem],
        coordinate_frame=_CARTESIAN_FRAME,
        transformation=(
            "even_rank2_native_T2: reconstruct Cartesian T, apply "
            "T'=R T R^T, then decompose in project-native T2 basis"
        ),
        validity=(
            "finite symmetric-traceless analytic source with omitted T0 and "
            "Cartesian-Levi-Civita T1 structural zeros; exact zero T2 can mean "
            "an empty source partition or cancellation and has no separate "
            "mask; any non-finite E/EFG rejects the whole AIMNet2Result"
        ),
        parity="even",
    )

_set_contract(
    ("aimnet2_charge_response_gradient",),
    coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "polar_vector: v'=R v"
    ),
    validity=(
        "separate AIMNet2ChargeResponseGradientResult: any missing autograd "
        "path or non-finite component rejects that result/current extraction; "
        "no manufactured zero. Catalog required=False retains pre-promotion "
        "file compatibility, although current successful CLI runs emit it"
    ),
    irreps="1o", parity="odd", tensor_rank=1, e3nn_export="")
_set_contract(
    ("aimnet2_charge_response_gradient_scalar",),
    coordinate_frame="intrinsic_norm_of_aimnet2_gradient",
    transformation=(
        "exact Euclidean norm of the emitted three Cartesian components; "
        "O(3)-invariant"
    ),
    validity=(
        "separate gradient result; finite nonnegative norm, including a real "
        "zero; absent rather than zero-filled on gradient failure. Catalog "
        "required=False retains pre-promotion file compatibility"
    ),
    irreps="0e", parity="even", tensor_rank=0)

_set_contract(
    ("orca_total", "orca_diamagnetic", "orca_paramagnetic"),
    coordinate_frame="orca_output_cartesian_xyz",
    transformation=(
        "even_rank2 in ORCA output frame: T'=R T R^T; parser-only producer "
        "does not validate or rotate ORCA coordinates into the conformation frame"
    ),
    validity="whole optional ORCA result absent when the parsed shielding source is unavailable",
    tensor_frame="orca_output_cartesian_xyz")

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
        "unit-current 1 nA row aligned with the union ring_contributions.npy "
        "axis; a whole default-zero row is ambiguous between no accepted BS "
        "source and physical zero because there is no BS-specific mask; B_phi "
        "is structural zero, and rho-axis degeneracy serializes B_rho=0 "
        "without a separate mask"
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
        "construction, so either may remain finite when column 6 is zero; "
        "column 6 does not validate col5, whose default zero is ambiguous"
    ), parity="mixed")
_set_contract(
    ("ringchi_scalar",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant (3*cos(theta)^2-1)/r^3 geometry scalar; "
        "translation invariant"
    ),
    validity=(
        "row-aligned with the union ring_contributions.npy axis; no dedicated "
        "RingSusceptibilityResult mask, so zero is physical-zero-or-unaccepted-row"
    ))
_set_contract(
    ("ringchi_per_type_T0",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "eight exact O(3)-invariant accepted-scalar sums in aromatic "
        "RingTypeIndex order"
    ),
    validity="zero is a physical empty/zero sum; no per-source availability mask")
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
    ("sidechain_co_source_bonds",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant sparse typed bond/atom/residue identity and "
        "semantic-class row; deterministic source order"
    ),
    validity=(
        "column7 is source_valid; the row axis is shared with "
        "sidechain_co_frame.npy and sidechain_co_frame_quality.npy"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("sidechain_co_frame_quality",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant bond length, orthogonality error, raw normal "
        "norm and frame-valid flag, row-aligned with sidechain_co_source_bonds.npy"
    ),
    validity=(
        "column3 is frame_valid and gates the axes in sidechain_co_frame.npy; "
        "earlier quality columns may remain finite on an invalid frame"
    ),
    irreps="0e", parity="even", tensor_rank=0)
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
    ("water_hbond_counts",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant per-protein-atom candidate/pass counts, nearest "
        "water identity and nearest passing mode for a fixed candidate set"
    ),
    validity=(
        "columns0:4 are candidate/pass counts and nearest-water index; "
        "column5 is nearest passing mode; producer sentinels distinguish no candidate"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("water_hbond_nearest",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant mixed scalar/identity row: distances, angle, "
        "mode, water identity, pass flag and candidate/pass counts"
    ),
    validity=(
        "nearest-candidate fields use the producer's no-candidate sentinels; "
        "candidate/pass counts in columns6 and7 expose row availability"
    ),
    irreps="0e", parity="even", tensor_rank=0)

_set_contract(
    ("delta_apbs",), coordinate_frame=_CARTESIAN_FRAME,
    transformation=(
        "mixed WT-atom row: cols0:3 polar delta-E; cols3:7 compatibility "
        "structural zeros (T0/T1); cols7:12 native T2 even-rank2"
    ),
    validity=(
        "whole NPY absent when MutationDeltaResult has no paired APBS delta; "
        "when present, unmatched WT atoms serialize zero and "
        "delta_scalars.npy column0 is the match mask"
    ),
    parity="mixed")
_set_contract(
    ("delta_scalars",), coordinate_frame="shared_wt_mut_intrinsic",
    transformation=(
        "exact O(3)-invariant WT-atom row under a joint rigid transform of "
        "WT and mutant: match mask, shielding T0 delta, removed-ring distance, "
        "charge deltas and typed-match distance"
    ),
    validity=(
        "column0 is the match mask for every delta payload; unmatched rows "
        "are legacy all-zero compatibility rows; on a matched row column4 "
        "MOPAC charge delta is unmasked zero both when MOPAC is unavailable "
        "and when the true compatibility-charge delta is zero"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("mutation_atom_map",), coordinate_frame="shared_wt_mut_intrinsic",
    transformation=(
        "exact O(3)- and translation-invariant typed WT-to-mutant atom "
        "identity"
    ),
    validity=(
        "-1 means unmatched; every nonnegative value indexes the atom axis "
        "of the separately emitted ala/ calcset"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("delta_graph",), coordinate_frame="shared_wt_mut_intrinsic",
    transformation=(
        "exact O(3)-invariant WT-atom graph row: match mask, graph-source "
        "mask, ring-distance delta, BFS decay delta and conjugation delta"
    ),
    validity=(
        "column0 is matched and column1 is has_graph_delta; unmatched or "
        "graph-unavailable fields use documented compatibility zeros"
    ),
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("delta_ring_proximity",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "repeated 6-column removed-ring blocks: distance/rho/factor/decay "
        "invariant; z pseudoscalar; theta maps to pi-theta under an improper transform"
    ),
    validity=(
        "whole NPY absent unless at least one mutation site removes a WT "
        "ring; when present, unmatched WT atoms serialize zero and "
        "delta_scalars.npy column0 is the match mask"
    ),
    parity="mixed")


# External-frame mixed diagnostics.
_set_contract(
    ("mopac_global",), coordinate_frame="mopac_output_cartesian_xyz",
    transformation=(
        "mixed graph row: col0 heat of formation invariant; cols1:4 molecular "
        "dipole polar vector. For a charged system the dipole is origin- and "
        "translation-dependent according to the MOPAC input-origin convention; "
        "componentwise F3 compatibility projection and independent MOZYME/SCF "
        "reruns prevent exact arbitrary-rotation covariance"
    ),
    validity=(
        "heat/dipole are finite-validated direct API values before F5/F3 "
        "compatibility quantization; failure omits the whole MopacResult"
    ),
    parity="mixed")
_set_contract(
    ("mopac_atom_populations",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "scalar population/valency row; compatibility cols5 and6:10 are NaN "
        "placeholders with no physical vector or transform meaning. Physical "
        "populations are scalars for a fixed electronic solution, but decimal "
        "projection and an independent MOZYME rerun are not promised exactly "
        "O(3)-invariant numerically"
    ),
    validity=(
        "one atom-aligned row; p/d are NaN without a live shell, and f plus "
        "per-atom-dipole text-only compatibility columns are always NaN"
    ), parity="even")
_set_contract(
    ("mopac_atomic_orbital_populations",),
    coordinate_frame="mopac_output_atomic_orbital_axes",
    transformation=(
        "frame-dependent diagonal s/px/py/pz/d populations are not a closed "
        "O(3) representation; omitted AO density coherences prevent rotation"
    ),
    validity=(
        "exactly one reconstructed row per atom; entries beyond the atom's "
        "live AO width are NaN, never a different unlabelled K-row axis"
    ),
    parity="even")
_set_contract(
    ("mopac_atomic_orbital_population_totals",),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "physical shell-trace scalars [s_total,p_total,d_total] for a fixed "
        "electronic solution; the producer first applies F10.5 to each AO "
        "diagonal before summing, and independent MOZYME reruns are not "
        "promised exact numerical O(3) invariance"
    ),
    validity=(
        "exactly one reconstructed row per atom; p/d are NaN when that shell "
        "does not exist, rather than a structural zero"
    ))
_set_contract(
    ("mopac_dipole_debye", "mopac_dipole_point_charge_debye"),
    coordinate_frame="mopac_output_cartesian_xyz",
    transformation=(
        "Cartesian polar vector v'=R v; for net-charged systems the total and "
        "point-charge dipoles are origin/translation dependent, and independent "
        "MOZYME/SCF reruns are not promised exact numerical covariance"
    ),
    validity="finite-validated direct API vector; whole MopacResult absent on failure",
    irreps="1o", parity="odd", tensor_rank=1)
_set_contract(
    ("mopac_dipole_hybridization_debye",),
    coordinate_frame="mopac_output_cartesian_xyz",
    transformation=(
        "one-centre Cartesian polar vector v'=R v and translation invariant; "
        "independent MOZYME/SCF reruns are not promised exact numerical covariance"
    ),
    validity="finite-validated direct API vector; whole MopacResult absent on failure",
    irreps="1o", parity="odd", tensor_rank=1)
_set_contract(
    ("mopac_atom_ao_density",),
    coordinate_frame="mopac_atom_centered_ao_basis",
    transformation=(
        "reducible AO-basis matrix: under an ideal fixed-shell rotation "
        "P'=D_AO(R) P D_AO(R)^T; not a Cartesian rank-2/e3nn tensor and not "
        "an exact independent-MOZYME-rerun covariance promise"
    ),
    validity=(
        "mopac_ao_orbitals_per_atom gives live widths; global-W cells beyond "
        "each live block are structural N/A zeros"
    ),
    irreps="", parity="mixed", tensor_rank=0)
_set_contract(
    ("mopac_bond_ao_density_directed", "mopac_bond_ao_density"),
    coordinate_frame="mopac_endpoint_ao_bases",
    transformation=(
        "reducible interatomic AO block: ideally P_ab'=D_a(R) P_ab D_b(R)^T; "
        "not a Cartesian rank-2/e3nn tensor and not an exact independent-"
        "MOZYME-rerun covariance promise"
    ),
    validity=(
        "only retained sparse API pairs have rows; endpoint cells beyond live "
        "AO widths are structural N/A zeros, and an omitted pair is not zero"
    ),
    irreps="", parity="mixed", tensor_rank=0)
_set_contract(
    ("mopac_atomic_orbital_populations_full_precision",),
    coordinate_frame="mopac_output_atomic_orbital_axes",
    transformation=(
        "diagonal s/px/py/pz/d AO populations are not a closed O(3) "
        "representation because off-diagonal coherences are omitted"
    ),
    validity=(
        "unquantized API binary64 values; global-W cells beyond live width are "
        "structural N/A zeros, not absent-shell measurements"
    ),
    parity="even")
_set_contract(
    ("mopac_atom_electron_population", "mopac_atom_s_population",
     "mopac_atom_p_population", "mopac_atom_d_population"),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "physical AO-block or complete-shell trace scalar for a fixed "
        "electronic solution; unquantized binary64 serialization removes the "
        "compatibility decimal projection but independent MOZYME reruns are "
        "not promised exact numerical O(3) invariance"
    ),
    validity=(
        "direct unquantized API-derived scalar; p/d use mathematical zero when "
        "the shell is absent, as explicitly distinct from compatibility NaN"
    ),
    irreps="", parity="even", tensor_rank=0)
_set_contract(
    ("mopac_lmo_occupied_coefficients", "mopac_lmo_virtual_coefficients"),
    coordinate_frame="mopac_localized_orbital_ao_gauge",
    transformation=(
        "no fixed rowwise O(3)/e3nn law: AO components rotate within shells, "
        "while localized orbitals also have sign, ordering and localization-gauge freedom"
    ),
    validity=(
        "packed live slices selected by the companion atom/support counts and "
        "native offsets; whole MopacResult absent on API/state validation failure"
    ),
    irreps="", parity="mixed", tensor_rank=0)

# Remaining diskless-MOPAC scalar, sparse-support, and native-state surfaces.
# Their ideal physical quantity may be invariant, but the producer reruns a
# localized SCF after eight-decimal coordinate preparation.  Do not turn that
# conditional physical law into an exact rowwise covariance promise.
_set_contract(
    ("mopac_charges", "mopac_scalars", "mopac_bond_valencies",
     "mopac_charges_full_precision", "mopac_bond_valencies_full_precision",
     "mopac_heat_kcal_mol"),
    coordinate_frame="intrinsic_mopac_electronic_solution",
    transformation=(
        "physical scalar(s) for a fixed PM7/MOZYME electronic solution under "
        "rigid O(3)/translation; independent SCF/localization reruns on the "
        "eight-decimal prepared coordinates are not promised exact numerical "
        "invariance, and compatibility decimal projection can amplify drift"
    ),
    validity=(
        "finite-validated calculator values; compatibility zero is a real "
        "serialized value, except sparse compact-valency zero only establishes "
        "that no retained first-six projected order contributed"
    ), irreps="", parity="even", tensor_rank=0)
_set_contract(
    ("mopac_bond_orders", "mopac_bond_neighbors",
     "mopac_bond_orders_unique", "mopac_bond_orders_full_precision",
     "mopac_bond_index", "mopac_bond_atom", "mopac_bond_order",
     "mopac_bond_density_pairs"),
    coordinate_frame="mopac_atom_identity_and_sparse_support",
    transformation=(
        "indices/pointers are discrete identities and Wiberg values are "
        "scalars conditional on identical retained support; an independent "
        "rerun can cross strict >0.01, F6.3, first-six or near-tie boundaries, "
        "changing row membership/order, so there is no unconditional rowwise "
        "O(3) correspondence"
    ),
    validity=(
        "omitted sparse pair means no row, never a measured zero; raw CSC "
        "pointer/atom/order companions must be interpreted together. Compact "
        "C rows add first-six/F6.3 censoring and mopac_bond_orders row order is "
        "unordered-map iteration order, not a sorted semantic axis"
    ), irreps="", parity="even", tensor_rank=0)
_set_contract(
    ("mopac_topology_bond_orders_full",),
    coordinate_frame="topology_bond_identity",
    transformation=(
        "fixed topology-bond indices and masks are identities; a present "
        "Wiberg order is scalar, but the present flag can change if an "
        "independent SCF rerun crosses the retained API >0.01 boundary"
    ),
    validity=(
        "columns3 and5 are NaN and col6 is one when the topology pair is "
        "absent from retained API support; col7 is always retired-artifact NaN"
    ), irreps="", parity="even", tensor_rank=0)
_set_contract(
    ("mopac_ao_max_orbitals", "mopac_ao_orbitals_per_atom"),
    coordinate_frame="mopac_atom_centered_ao_schema",
    transformation=(
        "discrete AO-basis extents, not Cartesian components; invariant for "
        "the same typed atoms and pinned libmopac parameterization"
    ),
    validity=(
        "global W is validated in [1,64]; every live atom width is in [1,W] "
        "and mozyme_state.iorbs must equal mopac_properties.ao_orbitals"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("mopac_lewis_bond_count", "mopac_lewis_bond_atoms"),
    coordinate_frame="mopac_mozyme_lewis_state",
    transformation=(
        "discrete localized-Lewis assignment, not geometric components; atom "
        "indices are identities conditional on the same chosen MOZYME state, "
        "but independent localization reruns have no unconditional rowwise "
        "O(3) correspondence"
    ),
    validity=(
        "counts are hard-validated in [0,9]; the first count slots are valid "
        "zero-based atoms and unused slots are structural -1 sentinels"
    ), irreps="", parity="even", tensor_rank=0)
_set_contract(
    ("mopac_lmo_energy_levels",),
    coordinate_frame="mopac_localized_orbital_gauge",
    transformation=(
        "energy is scalar only conditional on a matched localized orbital; "
        "LMO sign/order/localization/support are not fixed under rerun, so the "
        "array has no stable rowwise O(3) correspondence and is not a "
        "canonical HOMO/LUMO spectrum"
    ),
    validity=(
        "occupied rows precede virtual rows according to the companion state "
        "dimensions; every emitted energy is finite-validated"
    ), irreps="", parity="even", tensor_rank=0)
_set_contract(
    ("mopac_lmo_occupied_atom_counts", "mopac_lmo_occupied_atoms",
     "mopac_lmo_virtual_atom_counts", "mopac_lmo_virtual_atoms",
     "mopac_lmo_occupied_atom_offsets_native",
     "mopac_lmo_virtual_atom_offsets_native",
     "mopac_lmo_occupied_coefficient_offsets_native",
     "mopac_lmo_virtual_coefficient_offsets_native",
     "mopac_lmo_occupied_atom_storage_native",
     "mopac_lmo_virtual_atom_storage_native"),
    coordinate_frame="mopac_localized_orbital_gauge",
    transformation=(
        "discrete support/index state, not geometric components; meaningful "
        "only conditional on the same LMO ordering and localization, which an "
        "independent rerun is not promised to preserve rowwise"
    ),
    validity=(
        "counts and native offsets identify live support slices; compact atom "
        "lists are zero-based, native atom storage remains one-based, and "
        "allocated tail capacity is restart state rather than live support"
    ), irreps="", parity="even", tensor_rank=0)
_set_contract(
    ("mopac_lmo_occupied_coefficient_storage_native",
     "mopac_lmo_virtual_coefficient_storage_native"),
    coordinate_frame="mopac_localized_orbital_ao_gauge",
    transformation=(
        "no fixed rowwise O(3)/e3nn law: AO coefficients rotate within shells "
        "while localized orbitals have sign, ordering and localization-gauge "
        "freedom"
    ),
    validity=(
        "only companion native-offset-defined slices are live coefficients; "
        "remaining allocated capacity is native restart state"
    ), irreps="", parity="mixed", tensor_rank=0)
_set_contract(
    ("mopac_mozyme_state_dimensions",),
    coordinate_frame="mopac_native_state_schema",
    transformation=(
        "nonnegative electronic-state and allocation extents, not geometric "
        "components; values are conditional on the same returned MOZYME state "
        "and pinned build"
    ),
    validity=(
        "numat equals N; occupied/virtual counts and four storage extents are "
        "the authority for interpreting every LMO companion array"
    ), irreps="", parity="even", tensor_rank=0)

# The diskless worker/API result is atomic at the NPY boundary. Preserve each
# array's field-level rule above, then append the shared whole-result absence
# rule to every core/direct MOPAC payload.
for _stem, _spec in tuple(CATALOG.items()):
    if _spec.group in {"mopac_core", "mopac_direct"}:
        _validity = _spec.validity.rstrip("; ")
        CATALOG[_stem] = replace(
            _spec,
            validity=(
                f"{_validity}; {_MOPAC_RESULT_ABSENCE}"
                if _validity else _MOPAC_RESULT_ABSENCE
            ),
        )
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
_set_contract(
    ("larsen_sidechain_donor_atoms",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant per-atom typed donor identity, parent/residue "
        "identity and distance/angle-gated candidate/pass counts"
    ),
    validity=(
        "column0 identifies a typed sidechain polar H; remaining identity "
        "and count fields use the producer's non-donor sentinels"
    ),
    irreps="0e", parity="even", tensor_rank=0)


# Scalar, categorical, and structured records still need an explicit rigid-
# transform and absence contract: blank metadata is not a useful learning
# boundary even when no directional component is serialized.
_set_contract(
    ("element", "residue_index", "residue_type"),
    coordinate_frame="intrinsic_topology",
    transformation=(
        "exact O(3)- and translation-invariant per-atom integer identity; no "
        "coordinate component is serialized"
    ),
    validity="required dense atom axis with exactly one identity value per atom",
    irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("bond_length",), coordinate_frame=_INTRINSIC_FRAME,
    transformation="exact O(3)-invariant bond length; translation invariant",
    validity=(
        "bond_geometry_valid.npy gates zero/non-finite/degenerate bond rows; "
        "a valid zero length is not produced"
    ), irreps="0e", parity="even", tensor_rank=0)

_set_contract(
    ("enrichment_role", "enrichment_hybridisation", "enrichment_flags",
     "enrichment_parent_is_sp2"),
    coordinate_frame="intrinsic_topology",
    transformation=(
        "exact O(3)- and translation-invariant categorical/flag projection "
        "from atom, residue, ring and cached-backbone identity"
    ),
    validity=(
        "always computed by EnrichmentResult; Unknown/Unassigned/false are "
        "real compatibility categories, not missing-value sentinels"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("semantic_polar_h_kind", "semantic_planar_group_kind",
     "semantic_formal_charge", "semantic_ring_position", "semantic_locant",
     "enrichment_donor_class", "enrichment_acceptor_class",
     "enrichment_hybridisation_class"),
    coordinate_frame="intrinsic_topology",
    transformation=(
        "exact O(3)- and translation-invariant typed-semantic enum/charge "
        "projection; no directional component is serialized"
    ),
    validity=(
        "when AtomSemantic is unavailable the producer leaves zero defaults "
        "(hybridisation class 3), without a per-NPY mask; consult manifest "
        "has_atom_semantic, so zero can be absence or a real enum/charge value"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("molecular_graph_int", "molecular_graph_float", "molecular_graph"),
    coordinate_frame="intrinsic_topology",
    transformation=(
        "exact O(3)- and translation-invariant bond-graph distances, counts, "
        "element-weighted sums and decay features"
    ),
    validity=(
        "unreachable graph distances and nearest-ring atom use -1; unreachable "
        "ring decay uses 0, while neighbour sums/counts can be physical zero"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("ff_partial_charge", "ff_pb_radius"),
    coordinate_frame="intrinsic_charge_table",
    transformation=(
        "exact O(3)- and translation-invariant scalar table values from the "
        "configured typed ChargeSource"
    ),
    validity=(
        "whole ChargeAssignmentResult is absent if table preparation/row count "
        "fails; per-row status is not serialized, so a missing-charge numeric "
        "zero and placeholder PB radius 1.5 A are unmasked"
    ), irreps="0e", parity="even", tensor_rank=0)

_set_contract(
    ("bs_per_type_T0", "hm_per_type_T0"),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "eight contiguous exact O(3)-invariant T0 scalar sums in "
        "RingTypeIndex order; translation invariant"
    ),
    validity=(
        "zero block is an empty/cancelled accepted sum; RingTypeIndex 4 "
        "TrpPerimeter is structural zero because canonical BS/HM totals "
        "exclude the diagnostic perimeter ring"
    ), irreps="0e", parity="even", tensor_rank=0)
for _stem in ("bs_per_type_T1", "bs_per_type_T2", "hm_per_type_T1",
              "hm_per_type_T2"):
    CATALOG[_stem] = replace(
        CATALOG[_stem],
        validity=(
            CATALOG[_stem].validity + "; RingTypeIndex 4 TrpPerimeter block "
            "is structural zero because canonical BS/HM totals exclude it"
        ),
    )
_set_contract(
    ("pq_per_type_T0", "piquad_axial_scalar_per_type_T0"),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant per-type sum of (3cos^2(theta)-1)/r^4; "
        "translation invariant"
    ),
    validity=(
        "the two NPYs are additive aliases when both exist; zero is an empty/"
        "cancelled accepted geometry sum, and no physical prefactor is applied"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("disp_per_type_T0", "aromatic_r6_proximity_per_type_T0"),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant per-type switched-distance-kernel sum; "
        "translation invariant"
    ),
    validity=(
        "deprecated/canonical aliases of the same values; zero means no "
        "accepted vertex contribution after vertex/bonded exclusions"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("bs_ring_counts",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)- and translation-invariant counts in the 3/5/8/12 A shells"
    ),
    validity=(
        "zero means no BiotSavart-accepted ring row in that shell; counts are "
        "not over every topology ring"
    ), irreps="0e", parity="even", tensor_rank=0)

_set_contract(
    ("mc_nearfield_counts",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)- and translation-invariant accepted/rejected pair counts "
        "for source-midpoint distance below 3 A"
    ),
    validity=(
        "columns are accepted then filter-rejected queried pairs after the "
        "configured XH skip; zero is a real no-pair count"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("sidechain_co_scalar_audit",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "cols0:2 are project-T2 norms and col2 count/col3 distance are scalar "
        "invariants conditional on fixed Wiberg weights; an independent MOPAC "
        "rerun is not promised exact numerical covariance"
    ),
    validity=(
        "col0 zero can be empty/cancelled; col1 is NaN without MOPAC; col2 "
        "zero and col3 NaN mean no accepted SidechainCO source"
    ), irreps="", parity="even", tensor_rank=0)

_set_contract(
    ("coulomb_aromatic_n_src",), coordinate_frame=_INTRINSIC_FRAME,
    transformation="exact O(3)- and translation-invariant source count",
    validity=(
        "zero means no non-backbone aromatic source within cutoff passed the "
        "filters and charge floor"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("eeq_charges", "eeq_cn", "eeq_chi_eff", "eeq_hardness"),
    coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)- and translation-invariant element/distance-derived EEQ "
        "scalar channel(s)"
    ),
    validity=(
        "whole EeqResult is absent for empty input or a failed/non-finite "
        "solve; successful rows are finite and have no per-row mask"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("eeq_coulomb_aromatic_E_proj",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)- and translation-invariant signed projection of aromatic "
        "E onto the parent-to-H direction"
    ),
    validity=(
        "NaN without a valid H-parent direction; whole EEQ Coulomb result is "
        "absent on a non-finite field"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("eeq_coulomb_aromatic_n_src",), coordinate_frame=_INTRINSIC_FRAME,
    transformation="exact O(3)- and translation-invariant source count",
    validity=(
        "zero means no qualifying non-backbone aromatic source within cutoff "
        "passed the filters and charge floor"
    ), irreps="0e", parity="even", tensor_rank=0)

_set_contract(
    ("pyramidalization",), coordinate_frame=_INTRINSIC_FRAME,
    transformation=(
        "exact O(3)-invariant nonnegative magnitude of the signed out-of-plane "
        "distance; translation invariant"
    ),
    validity=(
        "NaN for non-applicable or invalid/degenerate rows; "
        "pyramidalization_valid.npy is the authoritative mask; whole trio "
        "absent when typed semantics are unavailable"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("pyramidalization_valid",), coordinate_frame=_INTRINSIC_FRAME,
    transformation="exact O(3)- and translation-invariant availability mask",
    validity=(
        "1 iff pyramidalization is finite; 0 combines non-applicability and "
        "invalid/degenerate geometry"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("pyramidalization_center_type",), coordinate_frame="intrinsic_topology",
    transformation="exact O(3)- and translation-invariant PlanarGroupKind enum",
    validity=(
        "nonzero marks typed applicability even if per-frame geometry is "
        "invalid; whole trio absent when typed semantics are unavailable"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("omega_is_xpro",), coordinate_frame="intrinsic_topology",
    transformation="exact O(3)- and translation-invariant categorical mask",
    validity=(
        "1 only when a covalent Pro successor and both CA atoms are cached; 0 "
        "conflates non-XPro with missing successor/CA, so consult omega_valid "
        "and residues.is_xpro_context; whole result absent without semantics"
    ), irreps="0e", parity="even", tensor_rank=0)

_set_contract(
    ("atoms_category_info",), coordinate_frame="intrinsic_topology",
    transformation=(
        "exact O(3)- and translation-invariant structured categorical/index/"
        "string record"
    ),
    validity=(
        "semantic bytes are zero-filled when AtomSemantic is absent, so enum/"
        "equivalence zero can be absence or real zero; force-field type may be "
        "empty and IUPAC/BMRB names may fall back to AMBER with provenance"
    ), irreps="0e", parity="even", tensor_rank=0)
_set_contract(
    ("residues", "bonds", "rings", "ring_membership"),
    coordinate_frame="intrinsic_topology",
    transformation=(
        "exact O(3)- and translation-invariant structured identity/enum/flag "
        "record; no Cartesian component is serialized"
    ),
    validity=(
        "required sidecar written even for empty semantic/ring content; -1 "
        "marks absent residue links/fused partners. Ring membership contains "
        "cyclic vertices only, so is_vertex=1 and is_substituent=0 are structural"
    ), irreps="0e", parity="even", tensor_rank=0)


# Unit/layout labels are explicit even for identity, enum, mask, count and
# mixed structured records.  This keeps an empty string from being mistaken
# for "dimensionless" at the learning boundary.
_set_units(("element",), "atomic_number")
_set_units(("residue_index",), "index")
_set_units(("residue_type",), "enum:AminoAcid")
_set_units(("ring_contributions",), (
    "mixed:index_or_enum[0:3],Å[3:6],radians[6],Å^-3[7],"
    "dimensionless[8],ppm_T_per_nA[9:18],Å^-1[18:36],Å^-6[36],"
    "count[37],dimensionless[38:40]"
))
_set_units(("ring_direction_to_center", "bond_direction",
            "mc_nearest_co_dir", "hbond_nearest_dir", "sasa_normal"),
           "dimensionless_unit_vector")
_set_units(("bond_geometry_valid", "tau_N_CA_C_valid",
            "angle_N_CA_CB_valid", "angle_CB_CA_C_valid",
            "angle_Cprev_N_CA_valid", "angle_CA_C_Nnext_valid",
            "cb_deviation_valid", "cb_residual_vector_valid",
            "enrichment_parent_is_sp2", "pyramidalization_valid",
            "omega_is_xpro", "omega_valid", "dssp_observed",
            "dssp_torsion_valid", "hbond_pairs_angle_valid"), "mask")
_set_units(("enrichment_role",), "enum:AtomRole")
_set_units(("enrichment_hybridisation",), "enum:Hybridisation")
_set_units(("enrichment_flags",), "mask_columns")
_set_units(("semantic_polar_h_kind", "enrichment_donor_class"),
           "enum:PolarHKind")
_set_units(("semantic_planar_group_kind",), "enum:PlanarGroupKind")
_set_units(("semantic_formal_charge",), "e")
_set_units(("semantic_ring_position",), "enum:RingPositionLabel")
_set_units(("semantic_locant",), "enum:Locant")
_set_units(("enrichment_acceptor_class",), "enum:acceptor_class")
_set_units(("enrichment_hybridisation_class",),
           "enum:hybridisation_class")
_set_units(("molecular_graph_int",),
           "mixed:graph_hops,count,mask,index")
_set_units(("molecular_graph_float",), "dimensionless")
_set_units(("molecular_graph",),
           "mixed:graph_hops,count,mask,index,Pauling_scale,dimensionless")
_set_units(("bs_ring_counts", "water_shell_counts", "larsen_hbond_count",
            "larsen_corner_imputed", "larsen_imputed_pair_count",
            "larsen_sidechain_carbonyl_pair_count"), "count")
_set_units(("piquad_local_frame",), "dimensionless_orthonormal_frame")
_set_units(("sidechain_co_source_bonds",),
           "mixed:index[0:4],enum[4:7],mask[7]")
_set_units(("sidechain_co_frame_quality",),
           "mixed:Å[0],dimensionless[1],Å^2[2],mask[3]")
_set_units(("sidechain_co_scalar_audit",),
           "mixed:Å^-3[0:2],count[2],Å[3]")
_set_units(("hbond_flags",), "mask_columns")
_set_units(("hbond_pairs_index",),
           "mixed:index[0:5],residue_separation_count[5]")
_set_units(("dssp_ss8",), "one_hot_class_mask")
_set_units(("dssp_ppii",), "enum:ternary_observation_class")
_set_units(("dssp_chi",),
           "dimensionless:cos_sin_exists_blocks")
_set_units(("dssp_torsion_sin", "dssp_torsion_cos", "omega_sin",
            "omega_cos", "atom_sasa_fraction", "eeq_cn"), "dimensionless")
_set_units(("dssp_hbond_partner_residue_index",), "index")
_set_units(("water_hbond_counts",),
           "mixed:count[0:4],index[4],enum[5]")
_set_units(("water_efield_clamp_mask", "water_efield_first_clamp_mask"),
           "mask")
_set_units(("water_polarization",),
           "mixed:e*Å[0:3],dimensionless[3:8],e*Å[8],count[9]")
_set_units(("delta_graph",),
           "mixed:mask[0:2],graph_hops[2],dimensionless[3:5]")
_set_units(("delta_apbs",),
           "mixed:V/Å[0:3],V/Å^2[3:12]")
_set_units(("atoms_category_info",),
           "structured:index,atomic_number,string,enum,mask,e")
_set_units(("aimnet2_aim", "aimnet2_aim_projection"),
           "model_native_learned_units_unspecified")
_set_units(("pyramidalization_center_type",), "enum:PlanarGroupKind")
_set_units(("larsen_hbond_pairs_index",),
           "mixed:index[0:6,7:11],enum[6],mask[11:16]")
_set_units(("larsen_hbond_pairs",),
           "mixed:index[0:6,7:11],enum[6],mask[11:16,19,21],"
           "Å[16],degrees[17:19],count[20],ppm[22:28]")
_set_units(("larsen_sidechain_donor_atoms",),
           "mixed:mask[0],enum[1],index[2:4],count[4:6]")
_set_units(("residues",), "structured:index,string,enum,count,mask")
_set_units(("bonds",), "structured:index,enum,mask")
_set_units(("rings",), "structured:index,enum,count")
_set_units(("ring_membership",), "structured:index,order,mask")
_set_units(("larsen_corner_imputed",), "mask")

# Positions are affine rank-1 coordinate fields even though their parity is
# not a homogeneous translation law.  Keep tensor_rank consistent with the
# ArraySpec definition used by downstream schema generation.
for _stem in ("pos", "mc_nearest_co_midpoint"):
    CATALOG[_stem] = replace(CATALOG[_stem], tensor_rank=1)


del _stem
