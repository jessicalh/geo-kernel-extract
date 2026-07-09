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
from ._ring import RingContributions, RingGeometry


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
      ``protein``.
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

    * ``parity`` is ``"even"`` / ``"odd"`` under spatial inversion;
      shieldings (rank 2 even), axial B fields, and most scalars are
      even. Polar vector fields like E are odd.

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
_T2_TENSOR_METADATA = dict(
    tensor_basis=_T2_BASIS,
    tensor_component_order=_T2_ORDER,
    tensor_frame="conformation_cartesian_xyz",
    structural_zero_components=_EFG_STRUCTURAL_ZEROS,
    e3nn_export=_E3NN_EXPORT,
)

CATALOG: dict[str, ArraySpec] = {s.stem: s for s in [
    # ── Identity (ConformationResult.cpp) ────────────────────────
    ArraySpec("pos",              "identity",   VectorField,       3,    True,  "Atom positions (A)",
              native_axis="atom", irreps="1o", units="Å", tensor_rank=1, parity="odd", mechanism="topology"),
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

    # ── Enrichment (EnrichmentResult.cpp) ───────────────────────────
    ArraySpec("enrichment_role",          "enrichment", np.ndarray, None, False, "AtomRole enum per atom (int32)",
              mechanism="topology"),
    ArraySpec("enrichment_hybridisation", "enrichment", np.ndarray, None, False, "Hybridisation enum per atom (int32)",
              mechanism="topology"),
    ArraySpec("enrichment_flags",         "enrichment", np.ndarray, 8,    False, "Enrichment boolean flags as int8 columns: is_backbone, is_amide_H, is_alpha_H, is_methyl, is_aromatic_H, is_hbond_donor, is_hbond_acceptor, is_on_aromatic_residue",
              mechanism="topology"),

    # ── Force-field charge assignment (ChargeAssignmentResult.cpp) ──
    ArraySpec("ff_partial_charge", "charge_assignment", np.ndarray, None, False, "Force-field partial charge per atom (e)",
              units="e", mechanism="charges"),
    ArraySpec("ff_pb_radius",      "charge_assignment", np.ndarray, None, False, "Force-field Poisson-Boltzmann radius per atom (Å)",
              units="Å", mechanism="charges"),

    # ── Biot-Savart (BiotSavartResult.cpp) ───────────────────────
    ArraySpec("bs_shielding",     "biot_savart", ShieldingTensor,  9,    True,  "BS ring current shielding",
              irreps=_SHIELD_IRREPS, units="ppm_T_per_nA", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="ring_current"),
    ArraySpec("bs_per_type_T0",   "biot_savart", PerRingTypeT0,    8,    True,  "BS isotropic per ring type",
              irreps="0e", units="ppm_T_per_nA", mechanism="ring_current"),
    ArraySpec("bs_per_type_T1",   "biot_savart", PerRingTypeT1,    24,   True,  "BS T1 per ring type",
              irreps="1e", units="ppm_T_per_nA", tensor_rank=1, parity="even", mechanism="ring_current"),
    ArraySpec("bs_per_type_T2",   "biot_savart", PerRingTypeT2,    40,   True,  "BS T2 per ring type",
              irreps="2e", units="ppm_T_per_nA", tensor_rank=2, mechanism="ring_current"),
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
    ArraySpec("piquad_quad_scalar", "pi_quadrupole", np.ndarray, None, False, "Sparse per-(atom,ring) rows aligned to ring_contributions: real computed RingNeighbourhood.quad_scalar, not the derived geometry scalar in ring_contributions column 7",
              native_axis="ring_contribution_pair", irreps="0e", units="Angstrom^-4", mechanism="ring_efg"),

    # ── Dispersion (DispersionResult.cpp) ────────────────────────
    ArraySpec("disp_per_type_T0", "dispersion", PerRingTypeT0,     8,    True,  "Dispersion scalar per ring type",
              irreps="0e", units="Angstrom^-6", mechanism="ring_dispersion"),

    # ── McConnell (McConnellResult.cpp) ──────────────────────────
    # Forward schema: packed SphericalTensor
    # [0e,1e_x,1e_y,1e_z,2e_m-2..+2] per source category and
    # source-strength channel. Source model metadata
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
    ArraySpec("mc_nearfield_counts",      "mcconnell", McConnellNearFieldCounts, 2, False, "McConnell near-field accepted/rejected source-target pair counts below 3 A", units="count", mechanism="bond_anisotropy"),
    ArraySpec("mc_nearest_co_dir",        "mcconnell", VectorField, 3, False, "Nearest accepted peptide C=O source direction per atom from ConformationAtom::dir_nearest_CO",
              irreps="1o", units="", tensor_rank=1, parity="odd", mechanism="bond_anisotropy"),
    ArraySpec("mc_nearest_co_midpoint",   "mcconnell", VectorField, 3, False, "Nearest accepted peptide C=O source midpoint per atom from ConformationAtom::nearest_CO_midpoint",
              irreps="1o", units="Å", tensor_rank=1, parity="odd", mechanism="bond_anisotropy"),
    ArraySpec("mc_nearest_co_T2",         "mcconnell", ShieldingTensor, 9, False, "Nearest accepted peptide C=O response per atom from ConformationAtom::T2_CO_nearest, packed [T0, T1x, T1y, T1z, T2_m-2..+2]",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_nearest_cn_T2",         "mcconnell", ShieldingTensor, 9, False, "Nearest accepted peptide C-N response per atom from ConformationAtom::T2_CN_nearest, packed [T0, T1x, T1y, T1z, T2_m-2..+2]",
              irreps=_SHIELD_IRREPS, units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    # Legacy McConnell arrays retained as optional/deprecated wrappers for
    # reading old extraction directories; new C++ emits the 15 tensor arrays above.
    ArraySpec("mc_shielding",     "mcconnell_legacy", ShieldingTensor,    9,    False,  "Legacy McConnell aggregate shielding", irreps=_SHIELD_IRREPS, units="Angstrom^-3",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_category_T2",   "mcconnell_legacy", PerBondCategoryT2,  25,   False,  "Legacy McConnell T2 per old bond category", irreps="2e", units="Angstrom^-3", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_scalars",       "mcconnell_legacy", McConnellScalars,   6,    False,  "Legacy McConnell scalar sums + distances", mechanism="bond_anisotropy"),

    # ── Coulomb (CoulombResult.cpp) — optional; retired from production
    # (APBS is canonical), so present only in the FullFatFrameExtraction
    # trajectory (--mopac), where it feeds the MOPAC-vs-FF14SB probe. ──
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
    ArraySpec("coulomb_aromatic_E_proj", "coulomb", np.ndarray,     None, False, "Coulomb aromatic E-field projection along the primary bond direction",
              irreps="0e", units="V/A", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_aromatic_n_src",  "coulomb", np.ndarray,     None, False, "Count of sidechain aromatic source atoms contributing to the Coulomb aromatic field (int32)", units="count", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_shielding",      "coulomb_legacy", ShieldingTensor, 9, False, "Legacy name for Coulomb bare total EFG", irreps=_SHIELD_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),

    # ── H-Bond (HBondResult.cpp) ─────────────────────────────────
    ArraySpec("hbond_scalars",    "hbond", HBondScalars,           4,    True,  "H-bond scalars (nearest_dist, 1/r^3, count, McConnell angular scalar Σ)",
              mechanism="hbond_kernel"),
    ArraySpec("hbond_nearest_dir", "hbond", VectorField,            3,    False, "Nearest accepted H-bond direction per atom from H-bond midpoint to target atom",
              irreps="1o", tensor_rank=1, parity="odd", mechanism="hbond_kernel"),
    ArraySpec("hbond_flags", "hbond", np.ndarray,                   3,    False, "H-bond boolean flags as int8 columns: hbond_is_backbone, hbond_is_donor, hbond_is_acceptor",
              mechanism="hbond_kernel"),

    # ── DSSP (DsspResult.cpp) ────────────────────────────────────
    ArraySpec("dssp_observed",    "dssp", np.ndarray,              None, False, "DSSP observation mask per atom: int8 1 when the parent residue mapped to a libdssp row, 0 otherwise",
              mechanism="secondary_structure"),
    ArraySpec("dssp_backbone",    "dssp", DsspScalars,             5,    True,  "DSSP backbone geometry, columns [phi_IUPAC_rad, psi_IUPAC_rad, sasa_A2, ss_helix, ss_sheet]: phi/psi are IUPAC-signed radians (-libdssp), SASA is A^2, helix/sheet are dimensionless flags; phi/psi/SASA are NaN for unobserved residues",
              units="radians_for_phi_psi; Angstrom^2_for_sasa; dimensionless_for_flags", mechanism="secondary_structure"),
    ArraySpec("dssp_ss8",         "dssp", np.ndarray,              8,    False, "DSSP 8-class SS one-hot (H/G/I/E/B/T/S/C)",
              mechanism="secondary_structure"),
    ArraySpec("dssp_hbond_energy","dssp", np.ndarray,              4,    False, "DSSP H-bond energies (acc0/acc1/don0/don1)",
              units="kcal/mol", mechanism="secondary_structure"),
    ArraySpec("dssp_chi",         "dssp", np.ndarray,              12,   False, "Chi1-4 cos/sin/exists (4 x 3 cols)",
              mechanism="secondary_structure"),

    # ── SASA (SasaResult.cpp) ───────────────────────────────────
    ArraySpec("atom_sasa",        "sasa", np.ndarray,              None, False, "Per-atom Shrake-Rupley SASA (A^2)",
              units="Å^2", mechanism="solvation"),
    ArraySpec("sasa_normal",      "sasa", VectorField,             3,    False, "SASA outward surface normal (unit vector)",
              irreps="1o", tensor_rank=1, parity="odd", mechanism="solvation"),

    # ── Explicit water (WaterFieldResult.cpp) ───────────────────
    ArraySpec("water_efield",       "water_field", VectorField,    3,    False, "Water E-field total (V/A)",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="solvation"),
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

    # ── EEQ charges (EeqResult.cpp — Caldeweyher 2019) ─────────
    ArraySpec("eeq_charges",        "eeq",         np.ndarray,     None, False, "EEQ geometry-dependent charges (e)",
              units="e", mechanism="charges"),
    ArraySpec("eeq_cn",             "eeq",         np.ndarray,     None, False, "EEQ coordination number",
              mechanism="charges"),

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
    ArraySpec("apbs_E",           "apbs", VectorField,             3,    False, "APBS solvated E-field",
              irreps="1o", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("apbs_efg",         "apbs", EFGTensor,               5,    False, "APBS solvated EFG (T2 only, symmetric-traceless)",
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
    ArraySpec("aimnet2_charge_response_gradient",        "aimnet2", AIMNet2ChargeResponseGradient, 3,    False, "AIMNet2 per-atom charge-response gradient (d(sum q_j^2)/d(r_i))",
              irreps="1o", units="e^2/Å", tensor_rank=1, parity="odd", mechanism="charges"),
    ArraySpec("aimnet2_charge_response_gradient_scalar", "aimnet2", np.ndarray,            None, False, "AIMNet2 per-atom charge-response gradient scalar (L2 norm of vector)",
              units="e^2/Å", mechanism="charges"),

    # ── Planar geometry (PlanarGeometryResult.cpp) ───────────────────
    # Per Amendment 2026-05-08(a). Conformation-only quantities derived
    # from positions; runs whenever the substrate (LegacyAmber
    # AtomSemanticTable) is populated. Nine NPYs with three different
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
    ArraySpec("tripeptide_neighbor_shielding",    "tripeptide", ShieldingTensor, 9,    False, "Δσ_BB^{i±1} — neighbour correction at residue i from i±1 cap reads (Larsen 2015 Eq 3)",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="quantum_reference"),
    ArraySpec("tripeptide_neighbor_residual_vec_prev", "tripeptide", VectorField, 3,   False, "Δσ_BB^{i-1} match residual at the C-term ALA cap of (i-1)'s tripeptide; Vec3, NaN where i-1 direction had no contribution",
              irreps="1o", units="Å", tensor_rank=1, parity="odd", mechanism="quantum_reference"),
    ArraySpec("tripeptide_neighbor_residual_vec_next", "tripeptide", VectorField, 3,   False, "Δσ_BB^{i+1} match residual at the N-term ALA cap of (i+1)'s tripeptide; Vec3, NaN where i+1 direction had no contribution",
              irreps="1o", units="Å", tensor_rank=1, parity="odd", mechanism="quantum_reference"),

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
    ArraySpec("larsen_hbond_water_term",                 "larsen_hbond", np.ndarray,      None, False, "Δσ_w = 2.07 ppm isotropic on amide H atoms with NO geometric H-bond candidate (θ ≥ 90° in 4.2 Å); proxies the NMA+water complex Larsen scanned for solvent-exposed amides", units="ppm", mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_count",                      "larsen_hbond", np.ndarray,      None, False, "Per-atom count of H-bond pairs that contributed under any of the four Table 2 classes; metadata, NOT a feature.", mechanism="hbond_grid"),
    ArraySpec("larsen_corner_imputed",                   "larsen_hbond", np.ndarray,      None, False, "Per-atom int8 flag: 1 iff any Larsen H-bond grid lookup corner serving this atom was imputed", mechanism="hbond_grid"),

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
