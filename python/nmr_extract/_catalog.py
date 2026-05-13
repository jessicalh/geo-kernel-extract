"""Format contract: every NPY file the C++ extractor can produce.

This is the single source of truth for what's on disk. If the C++ side
adds a new array, add an entry here. Nothing else in the SDK needs to
change — the loader picks it up automatically.

Generated from the WriteFeatures methods in src/*Result.cpp and
src/ConformationResult.cpp.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from ._tensors import (
    ShieldingTensor,
    EFGTensor,
    VectorField,
    PerRingTypeT0,
    PerRingTypeT2,
    PerBondCategoryT2,
    RingCounts,
    McConnellScalars,
    CoulombScalars,
    HBondScalars,
    DsspScalars,
    MopacScalars,
    MopacGlobal,
    BondOrders,
    DeltaScalars,
    DeltaAPBS,
    DeltaRingProximity,
    AIMNet2Charges,
    AIMNet2AimEmbedding,
    AIMNet2Polarisability,
)
from ._ring import RingContributions, RingGeometry


@dataclass(frozen=True)
class ArraySpec:
    """Metadata for one NPY file produced by the C++ extractor.

    `is_feature` is the ML-eligibility flag: True means this array is
    intended as input to downstream ridge / model fitting; False means
    it is metadata (counts, water-term-as-isotropic-offset, parser-
    pipeline diagnostics) or an aggregate that would introduce a
    structural linear dependence with other arrays already in the
    feature set (e.g. a `total = sum(per-class)` shielding). Loaders
    in `learn/` should filter by `is_feature=True` when assembling
    the regression design matrix. Default True so existing arrays
    keep their current behaviour.

    Topology-sidecar metadata (2026-05-13, resolves
    TENTATIVE_OUTSTANDING_ISSUES OI-016):

    * ``native_axis`` declares which axis this array's rows are
      indexed along: ``atom`` / ``residue`` / ``aromatic_ring`` /
      ``saturated_ring`` / ``ring`` / ``ring_contribution_pair`` /
      ``bond`` / ``ring_membership`` / ``mutation_match_pair`` /
      ``protein``. R / Python analysis must read this metadata
      column rather than infer axis from filename.

    * ``irreps`` carries the e3nn-style irrep decomposition for
      tensor-valued arrays: ``"0e + 1e + 2e"`` for a 9-component
      SphericalTensor (T0 + T1 + T2), ``"1e"`` for a Vec3,
      ``"0e + 2e"`` for an EFG (symmetric traceless), ``""`` for
      scalar / structured / categorical arrays where irreps do not
      apply.

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
      shieldings (rank 2 even) and most scalars are even, vector
      fields like E and B are odd.

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
    is_feature: bool = True
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


# fmt: off
# Common metadata strings — name once, use across entries of the same physics.
_SHIELD_IRREPS = "0e + 1e + 2e"
_SHIELD_SIGN   = "σ_ab = -dB_a^sec/dB_{0,b}"
_EFG_IRREPS    = "0e + 2e"   # symmetric traceless

CATALOG: dict[str, ArraySpec] = {s.stem: s for s in [
    # ── Identity (ConformationResult.cpp) ────────────────────────
    ArraySpec("pos",              "identity",   VectorField,       3,    True,  "Atom positions (A)",
              native_axis="atom", irreps="1e", units="Å", tensor_rank=1, parity="odd", mechanism="topology"),
    ArraySpec("element",          "identity",   np.ndarray,        None, True,  "Atomic number (int32)",
              native_axis="atom", mechanism="topology"),
    ArraySpec("residue_index",    "identity",   np.ndarray,        None, True,  "Residue index (int32)",
              native_axis="atom", mechanism="topology"),
    ArraySpec("residue_type",     "identity",   np.ndarray,        None, True,  "Residue type enum (int32)",
              native_axis="atom", mechanism="topology"),
    ArraySpec("ring_contributions","identity",  RingContributions, 59,   True,  "Per-(atom,ring) pair contributions",
              native_axis="ring_contribution_pair", irreps=_SHIELD_IRREPS, units="ppm",
              sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="ring_current"),
    ArraySpec("ring_geometry",    "identity",   RingGeometry,      10,   True,  "Per-ring geometry reference",
              native_axis="aromatic_ring", units="Å", mechanism="topology"),

    # ── Biot-Savart (BiotSavartResult.cpp) ───────────────────────
    ArraySpec("bs_shielding",     "biot_savart", ShieldingTensor,  9,    True,  "BS ring current shielding",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="ring_current"),
    ArraySpec("bs_per_type_T0",   "biot_savart", PerRingTypeT0,    8,    True,  "BS isotropic per ring type",
              irreps="0e", units="ppm", mechanism="ring_current"),
    ArraySpec("bs_per_type_T2",   "biot_savart", PerRingTypeT2,    40,   True,  "BS T2 per ring type",
              irreps="2e", units="ppm", tensor_rank=2, mechanism="ring_current"),
    ArraySpec("bs_total_B",       "biot_savart", VectorField,      3,    True,  "BS total B-field vector",
              irreps="1e", units="T", tensor_rank=1, parity="odd", mechanism="ring_current"),
    ArraySpec("bs_ring_counts",   "biot_savart", RingCounts,       4,    True,  "Ring proximity counts (3/5/8/12 A)",
              mechanism="ring_current"),

    # ── Haigh-Mallion (HaighMallionResult.cpp) ───────────────────
    ArraySpec("hm_shielding",     "haigh_mallion", ShieldingTensor, 9,   True,  "HM ring current shielding",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="ring_current"),
    ArraySpec("hm_per_type_T0",   "haigh_mallion", PerRingTypeT0,   8,   True,  "HM isotropic per ring type",
              irreps="0e", units="ppm", mechanism="ring_current"),
    ArraySpec("hm_per_type_T2",   "haigh_mallion", PerRingTypeT2,   40,  True,  "HM T2 per ring type",
              irreps="2e", units="ppm", tensor_rank=2, mechanism="ring_current"),

    # ── Pi-Quadrupole (PiQuadrupoleResult.cpp) ───────────────────
    ArraySpec("pq_shielding",     "pi_quadrupole", ShieldingTensor, 9,   True,  "PQ EFG shielding",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="ring_efg"),
    ArraySpec("pq_per_type_T0",   "pi_quadrupole", PerRingTypeT0,   8,   True,  "PQ isotropic per ring type",
              irreps="0e", units="ppm", mechanism="ring_efg"),
    ArraySpec("pq_per_type_T2",   "pi_quadrupole", PerRingTypeT2,   40,  True,  "PQ T2 per ring type",
              irreps="2e", units="ppm", tensor_rank=2, mechanism="ring_efg"),

    # ── Dispersion (DispersionResult.cpp) ────────────────────────
    ArraySpec("disp_shielding",   "dispersion", ShieldingTensor,   9,    True,  "Dispersion shielding (1/r^6)",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="ring_dispersion"),
    ArraySpec("disp_per_type_T0", "dispersion", PerRingTypeT0,     8,    True,  "Dispersion scalar per ring type",
              irreps="0e", units="ppm", mechanism="ring_dispersion"),
    ArraySpec("disp_per_type_T2", "dispersion", PerRingTypeT2,     40,   True,  "Dispersion T2 per ring type",
              irreps="2e", units="ppm", tensor_rank=2, mechanism="ring_dispersion"),

    # ── Ring Susceptibility (RingSusceptibilityResult.cpp) ───────
    ArraySpec("ringchi_shielding","ring_susceptibility", ShieldingTensor, 9, True, "Ring susceptibility shielding",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="ring_current"),

    # ── McConnell (McConnellResult.cpp) ──────────────────────────
    ArraySpec("mc_shielding",     "mcconnell", ShieldingTensor,    9,    True,  "McConnell bond anisotropy shielding",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_category_T2",   "mcconnell", PerBondCategoryT2,  25,   True,  "McConnell T2 per bond category",
              irreps="2e", units="ppm", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mc_scalars",       "mcconnell", McConnellScalars,   6,    True,  "McConnell scalar sums + distances",
              mechanism="bond_anisotropy"),

    # ── Coulomb (CoulombResult.cpp) — optional via --no-coulomb ──
    ArraySpec("coulomb_shielding",      "coulomb", ShieldingTensor, 9,   False, "Coulomb E-field shielding",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="electrostatic_efg"),
    ArraySpec("coulomb_E",              "coulomb", VectorField,     3,   False, "Coulomb total E-field",
              irreps="1e", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("coulomb_efg_backbone",   "coulomb", EFGTensor,       9,   False, "Coulomb EFG backbone",
              irreps=_EFG_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),
    ArraySpec("coulomb_efg_aromatic",   "coulomb", EFGTensor,       9,   False, "Coulomb EFG aromatic",
              irreps=_EFG_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),
    ArraySpec("coulomb_scalars",        "coulomb", CoulombScalars,  4,   False, "Coulomb E-field scalars",
              mechanism="electrostatic_efg"),

    # ── H-Bond (HBondResult.cpp) ─────────────────────────────────
    ArraySpec("hbond_shielding",  "hbond", ShieldingTensor,        9,    True,  "H-bond shielding",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="hbond_kernel"),
    ArraySpec("hbond_scalars",    "hbond", HBondScalars,           3,    True,  "H-bond scalars",
              mechanism="hbond_kernel"),

    # ── DSSP (DsspResult.cpp) ────────────────────────────────────
    ArraySpec("dssp_backbone",    "dssp", DsspScalars,             5,    True,  "DSSP backbone geometry",
              units="radians", mechanism="secondary_structure"),
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
              irreps="1e", tensor_rank=1, parity="odd", mechanism="solvation"),

    # ── Explicit water (WaterFieldResult.cpp) ───────────────────
    ArraySpec("water_efield",       "water_field", VectorField,    3,    False, "Water E-field total (V/A)",
              irreps="1e", units="V/A", tensor_rank=1, parity="odd", mechanism="solvation"),
    ArraySpec("water_efield_first", "water_field", VectorField,    3,    False, "Water E-field first shell <3.5A (V/A)",
              irreps="1e", units="V/A", tensor_rank=1, parity="odd", mechanism="solvation"),
    ArraySpec("water_efg",          "water_field", EFGTensor,      9,    False, "Water EFG total (SphericalTensor)",
              irreps=_EFG_IRREPS, units="V/A^2", tensor_rank=2, mechanism="solvation"),
    ArraySpec("water_efg_first",    "water_field", EFGTensor,      9,    False, "Water EFG first shell (SphericalTensor)",
              irreps=_EFG_IRREPS, units="V/A^2", tensor_rank=2, mechanism="solvation"),
    ArraySpec("water_shell_counts", "water_field", np.ndarray,     2,    False, "Water shell counts [n_first, n_second]",
              mechanism="solvation"),

    # ── Hydration shell (HydrationShellResult.cpp) ──────────────
    ArraySpec("hydration_shell",    "hydration",   np.ndarray,     4,    False, "Hydration geometry [asymmetry, dipole_cos, ion_dist, ion_charge]",
              mechanism="solvation"),

    # ── Hydration geometry — SASA-normal (HydrationGeometryResult.cpp) ─
    ArraySpec("water_polarization", "water_polarization", np.ndarray, 10, False, "Water polarisation [dipole(3), normal(3), asym, align, coher, count]",
              mechanism="solvation"),

    # ── EEQ charges (EeqResult.cpp — Caldeweyher 2019) ─────────
    ArraySpec("eeq_charges",        "eeq",         np.ndarray,     None, False, "EEQ geometry-dependent charges (e)",
              units="e", mechanism="charges"),
    ArraySpec("eeq_cn",             "eeq",         np.ndarray,     None, False, "EEQ coordination number",
              mechanism="charges"),

    # ── GROMACS energy (GromacsEnergyResult.cpp) ────────────────
    ArraySpec("gromacs_energy",     "gromacs",     np.ndarray,     42,   False, "Per-frame energy (42 cols: electrostatic, bonded, VdW, thermo, box, virial, pressure tensor, T_group)",
              native_axis="protein", units="kJ/mol", mechanism="gromacs_runtime"),

    # ── Bonded energy (BondedEnergyResult.cpp) ─────────────────
    ArraySpec("bonded_energy",      "bonded",      np.ndarray,      7,   False, "Per-atom bonded energy (bond,angle,UB,proper,improper,CMAP,total) kJ/mol",
              units="kJ/mol", mechanism="gromacs_runtime"),

    # ── MOPAC core (MopacResult.cpp) ─────────────────────────────
    ArraySpec("mopac_charges",    "mopac_core", np.ndarray,        None, False, "MOPAC Mulliken charges",
              units="e", mechanism="charges"),
    ArraySpec("mopac_scalars",    "mopac_core", MopacScalars,      4,    False, "MOPAC per-atom scalars",
              mechanism="charges"),
    ArraySpec("mopac_bond_orders","mopac_core", BondOrders,        3,    False, "MOPAC Wiberg bond orders",
              mechanism="charges"),
    ArraySpec("mopac_global",     "mopac_core", MopacGlobal,       4,    False, "MOPAC graph-level scalars",
              native_axis="protein", mechanism="charges"),

    # ── MOPAC Coulomb (MopacCoulombResult.cpp) ───────────────────
    ArraySpec("mopac_coulomb_shielding",     "mopac_coulomb", ShieldingTensor, 9,  False, "MOPAC Coulomb shielding",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_E",             "mopac_coulomb", VectorField,     3,  False, "MOPAC Coulomb E-field",
              irreps="1e", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_efg_backbone",  "mopac_coulomb", EFGTensor,       9,  False, "MOPAC Coulomb EFG backbone",
              irreps=_EFG_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_efg_aromatic",  "mopac_coulomb", EFGTensor,       9,  False, "MOPAC Coulomb EFG aromatic",
              irreps=_EFG_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),
    ArraySpec("mopac_coulomb_scalars",       "mopac_coulomb", CoulombScalars,  4,  False, "MOPAC Coulomb scalars",
              mechanism="electrostatic_efg"),

    # ── MOPAC McConnell (MopacMcConnellResult.cpp) ───────────────
    ArraySpec("mopac_mc_shielding",    "mopac_mcconnell", ShieldingTensor,   9,  False, "MOPAC McConnell shielding",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mopac_mc_category_T2",  "mopac_mcconnell", PerBondCategoryT2, 25, False, "MOPAC McConnell T2 per category",
              irreps="2e", units="ppm", tensor_rank=2, mechanism="bond_anisotropy"),
    ArraySpec("mopac_mc_scalars",      "mopac_mcconnell", McConnellScalars,  6,  False, "MOPAC McConnell scalars",
              mechanism="bond_anisotropy"),

    # ── APBS (ApbsFieldResult.cpp) ───────────────────────────────
    ArraySpec("apbs_E",           "apbs", VectorField,             3,    False, "APBS solvated E-field",
              irreps="1e", units="V/A", tensor_rank=1, parity="odd", mechanism="electrostatic_efg"),
    ArraySpec("apbs_efg",         "apbs", EFGTensor,               9,    False, "APBS solvated EFG",
              irreps=_EFG_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),

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
    ArraySpec("delta_apbs",            "delta", DeltaAPBS,             12,   False, "APBS delta E + EFG",
              native_axis="mutation_match_pair", mechanism="mutation_delta"),
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
    ArraySpec("aimnet2_efg",                 "aimnet2", EFGTensor,                 9,    True,  "AIMNet2 Coulomb EFG total",
              irreps=_EFG_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),
    ArraySpec("aimnet2_efg_aromatic",        "aimnet2", EFGTensor,                 9,    True,  "AIMNet2 Coulomb EFG aromatic",
              irreps=_EFG_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),
    ArraySpec("aimnet2_efg_backbone",        "aimnet2", EFGTensor,                 9,    True,  "AIMNet2 Coulomb EFG backbone",
              irreps=_EFG_IRREPS, units="V/A^2", tensor_rank=2, mechanism="electrostatic_efg"),

    # ── AIMNet2 polarisability (AIMNet2PolarisabilityResult.cpp) ─────
    # Always-on after the JobSpec --aimnet2 model is loaded (per the
    # 2026-05-09 promotion of Amendment 2026-05-08(b) from a test flag
    # to standard non-trajectory pipeline; trajectory mode unchanged).
    # Vector is dL/d(r_i) where L = sum_j q_j^2 over non-sentinel atoms
    # (charge-conservation makes sum(q) gradient near-zero, so the L2
    # of charges is the cheapest single-pass objective with non-trivial
    # gradient). Scalar is the L2 norm of the vector. required=False
    # because old extraction outputs (pre-2026-05-09) do not include
    # them; new outputs always do when AIMNet2 is loaded.
    ArraySpec("aimnet2_polarisability",        "aimnet2", AIMNet2Polarisability, 3,    False, "AIMNet2 per-atom polarisability gradient (d(sum q_j^2)/d(r_i))",
              irreps="1e", units="e^2/Å", tensor_rank=1, parity="odd", mechanism="charges"),
    ArraySpec("aimnet2_polarisability_scalar", "aimnet2", np.ndarray,            None, False, "AIMNet2 per-atom polarisability scalar (L2 norm of gradient)",
              units="e^2/Å", mechanism="charges"),

    # ── Planar geometry (PlanarGeometryResult.cpp) ───────────────────
    # Per Amendment 2026-05-08(a). Conformation-only quantities derived
    # from positions; runs whenever the substrate (LegacyAmber
    # AtomSemanticTable) is populated. Six NPYs with three different
    # axes (per-atom, per-residue, per-ring) — the reader's wrapper
    # picks shape from the catalog cols field.
    #
    # Conventions: pyramidalization sign by improper-dihedral right-
    # hand rule; omega in radians; Cremer-Pople θ in degrees [0, 360).
    # required=False because pre-2026-05-09 extractions don't carry them.
    ArraySpec("pyramidalization",  "planar_geometry", np.ndarray, None, False, "Per-atom signed sp2 out-of-plane displacement (Å); 0 for non-planar atoms",
              units="Å", mechanism="geometry"),
    ArraySpec("omega_actual",      "planar_geometry", np.ndarray, None, False, "Per-residue ω (Cα-C-N-Cα to next), radians; NaN at C-term and X-Pro",
              native_axis="residue", units="radians", mechanism="geometry"),
    ArraySpec("omega_deviation",   "planar_geometry", np.ndarray, None, False, "Per-residue ω - π wrapped to (-π, π]; NaN where omega_actual is NaN",
              native_axis="residue", units="radians", mechanism="geometry"),
    ArraySpec("aromatic_chi2",     "planar_geometry", np.ndarray, None, False, "Per-aromatic-ring χ₂ (parent residue, radians); ring-flip observable per Akke-Weininger 2023",
              native_axis="aromatic_ring", units="radians", mechanism="geometry"),
    ArraySpec("pucker_Q",          "planar_geometry", np.ndarray, None, False, "Per-saturated-ring Cremer-Pople puckering amplitude (Å); 5-rings only",
              native_axis="saturated_ring", units="Å", mechanism="geometry"),
    ArraySpec("pucker_theta",      "planar_geometry", np.ndarray, None, False, "Per-saturated-ring Cremer-Pople phase angle (degrees, [0, 360))",
              native_axis="saturated_ring", units="degrees", mechanism="geometry"),
    ArraySpec("omega_is_xpro",     "planar_geometry", np.ndarray, None, False, "Per-residue mask: 1 where the bond into i+1 is X→Pro (cis/trans isomerism is real signal there, not 'non-planar amide' deviation)",
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
              irreps="1e", units="Å", tensor_rank=1, parity="odd", mechanism="quantum_reference"),
    ArraySpec("tripeptide_bb_match_distance",     "tripeptide", np.ndarray,      None, False, "σ_BB^i match distance (Å) — magnitude of residual_vec",
              units="Å", mechanism="quantum_reference"),
    ArraySpec("tripeptide_bb_method_tag",         "tripeptide", np.ndarray,      None, False, "DFT method discriminator: 0=none, 1=OPBE Gaussian (Larsen), 2=PBE ORCA (project SER regen)",
              mechanism="quantum_reference"),
    ArraySpec("tripeptide_neighbor_shielding",    "tripeptide", ShieldingTensor, 9,    False, "Δσ_BB^{i±1} — neighbour correction at residue i from i±1 cap reads (Larsen 2015 Eq 3)",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="quantum_reference"),
    ArraySpec("tripeptide_neighbor_residual_vec_prev", "tripeptide", VectorField, 3,   False, "Δσ_BB^{i-1} match residual at the C-term ALA cap of (i-1)'s tripeptide; Vec3, NaN where i-1 direction had no contribution",
              irreps="1e", units="Å", tensor_rank=1, parity="odd", mechanism="quantum_reference"),
    ArraySpec("tripeptide_neighbor_residual_vec_next", "tripeptide", VectorField, 3,   False, "Δσ_BB^{i+1} match residual at the N-term ALA cap of (i+1)'s tripeptide; Vec3, NaN where i+1 direction had no contribution",
              irreps="1e", units="Å", tensor_rank=1, parity="odd", mechanism="quantum_reference"),

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
    # kernel-form HBondResult — see feedback_methods_accumulate memory
    # entry.
    #
    # ML ELIGIBILITY NOTE (per audit 2026-05-13):
    #   `larsen_hbond_shielding` is the ELEMENT-WISE SUM of the four
    #   per-class tensors. Including it AND the per-class breakdown in
    #   the same regression introduces a structural linear dependence
    #   that ridge will silently distribute coefficients across. The
    #   per-class breakdown is the more informative form (each class
    #   carries distinct physics per Larsen Table 1) so the aggregate
    #   is marked `is_feature=False`. Pick one form per regression.
    #
    #   `larsen_hbond_diagnostic_CB_shielding` is a parser-pipeline
    #   reality check (Larsen Table 2 says Cβ gets NO H-bond term; we
    #   emit anyway to verify the rotation pipeline produces what the
    #   physics predicts). NOT a feature — `is_feature=False`.
    #
    #   `larsen_hbond_count` is per-atom metadata (number of pairs the
    #   atom received contributions from). Strongly correlated with
    #   shielding magnitude; not a clean physics feature on its own.
    #   `is_feature=False`.
    #
    #   `larsen_hbond_water_term` IS a real Larsen term (Δσ_w = 2.07
    #   ppm isotropic on solvent-exposed amide Hs) and stays
    #   `is_feature=True` even though it's binary 0-or-2.07.
    # ────────────────────────────────────────────────────────────────
    ArraySpec("larsen_hbond_shielding",                  "larsen_hbond", ShieldingTensor, 9, False, "Σ Larsen H-bond contributions across all four Table 2 classes (1°HB + 2°HB + 1°HαB + 2°HαB) — ppm, lab frame. Structurally = sum of the four per-class arrays; NOT a feature.",
              is_feature=False, irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_1pHB_shielding",             "larsen_hbond", ShieldingTensor, 9, False, "Δσ_1°HB per Larsen 2015 Table 2 — primary amide-H donor contribution applied to donor residue i atoms",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_2pHB_shielding",             "larsen_hbond", ShieldingTensor, 9, False, "Δσ_2°HB per Larsen 2015 Table 2 — secondary amide-H donor contribution. N/Cα/Hα/HN apply to acceptor residue j+1; C' applies to acceptor's OWN residue j",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_1pHaB_shielding",            "larsen_hbond", ShieldingTensor, 9, False, "Δσ_1°HαB per Larsen 2015 Table 2 — primary Hα donor contribution applied to donor residue i atoms",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_2pHaB_shielding",            "larsen_hbond", ShieldingTensor, 9, False, "Δσ_2°HαB per Larsen 2015 Table 2 — secondary Hα donor contribution. N/Cα/Hα/HN apply to acceptor residue j+1; C' applies to acceptor's OWN residue j",
              irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_diagnostic_CB_shielding",    "larsen_hbond", ShieldingTensor, 9, False, "Cβ diagnostic — Larsen Table 2 says Cβ gets NO contribution; emitted as parser→loader→frame-rotation reality check. NOT a feature.",
              is_feature=False, irreps=_SHIELD_IRREPS, units="ppm", sign_convention=_SHIELD_SIGN, tensor_rank=2, mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_water_term",                 "larsen_hbond", np.ndarray,      None, False, "Δσ_w = 2.07 ppm isotropic on amide H atoms with NO geometric H-bond candidate (θ ≥ 90° in 4.2 Å); proxies the NMA+water complex Larsen scanned for solvent-exposed amides",
              units="ppm", mechanism="hbond_grid"),
    ArraySpec("larsen_hbond_count",                      "larsen_hbond", np.ndarray,      None, False, "Per-atom count of H-bond pairs that contributed under any of the four Table 2 classes; metadata, NOT a feature.",
              is_feature=False, mechanism="hbond_grid"),

    # ────────────────────────────────────────────────────────────────
    # Topology sidecar (TopologySidecar.cpp) — additive 2026-05-13.
    # Three structured-NPY projections of LegacyAmberTopology +
    # RingTopology. Emitted alongside atoms_category_info.npy whenever
    # the protein has a populated typed substrate.
    # ────────────────────────────────────────────────────────────────
    ArraySpec("bonds",            "topology", np.ndarray, None, True,  "Per-bond record: bond_index, atom_index_a/b, bond_order, bond_category, is_rotatable, is_aromatic, is_peptide, is_backbone",
              native_axis="bond", is_feature=False, mechanism="topology"),
    ArraySpec("rings",            "topology", np.ndarray, None, True,  "Per-ring record: ring_id, ring_kind (aromatic|saturated), ring_type_index, atom_count, native_axis_index, parent_residue_index, parent_residue_number, fused_partner_ring_id",
              native_axis="ring", is_feature=False, mechanism="topology"),
    ArraySpec("ring_membership",  "topology", np.ndarray, None, True,  "Per (ring, ring-vertex-atom) row: ring_id, atom_index, ring_atom_order, is_vertex, is_substituent",
              native_axis="ring_membership", is_feature=False, mechanism="topology"),
]}
# fmt: on


def feature_specs() -> dict[str, ArraySpec]:
    """Return only the catalog entries marked as ML features.

    Loaders in `learn/` should call this when assembling the regression
    design matrix so that metadata arrays (counts, diagnostic CB,
    structural-sum aggregates) are not silently handed to ridge. See
    the ArraySpec docstring for the `is_feature` semantics.
    """
    return {stem: spec for stem, spec in CATALOG.items() if spec.is_feature}
