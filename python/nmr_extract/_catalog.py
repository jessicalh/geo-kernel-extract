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
    AIMNet2ChargeSensitivity,
)
from ._ring import RingContributions, RingGeometry


@dataclass(frozen=True)
class ArraySpec:
    """Metadata for one NPY file produced by the C++ extractor."""
    stem: str                   # filename without .npy
    group: str                  # logical group (e.g. "biot_savart")
    wrapper: type               # Python wrapper class
    cols: Optional[int]         # expected last dim (None = variable or 1D)
    required: bool              # False for optional calculators
    description: str


# fmt: off
CATALOG: dict[str, ArraySpec] = {s.stem: s for s in [
    # ── Identity (ConformationResult.cpp) ────────────────────────
    ArraySpec("pos",                "identity", VectorField,       3,    True,  "Atom positions (A)"),
    ArraySpec("element",            "identity", np.ndarray,        None, True,  "Atomic number (int32)"),
    ArraySpec("residue_index",      "identity", np.ndarray,        None, True,  "Residue index (int32)"),
    ArraySpec("residue_type",       "identity", np.ndarray,        None, True,  "Residue type enum (int32)"),
    # 2026-04-26: rich per-atom topology + role exports for post-hoc
    # classification by IUPAC topology rules (Markley 1998).
    #
    # Marked required=False today only because the SDK's golden fixtures
    # under tests/data/sdk_geo_only/ predate this surface and would all
    # fail loading if these were promoted to required. The extractor
    # always emits them after FinalizeConstruction + EnrichmentResult,
    # which every job spec includes — so once the fixtures are
    # regenerated, flip required to True. Tracked: regenerate fixtures
    # then promote.
    ArraySpec("atom_topology",      "identity", np.ndarray,        14,   False, "AtomTopology (locant, branch, diast, prochiral, planar, pseudo, ring_pos, polar_h_kind, stamped, ring_count, chi[0..3]) (int32)"),
    ArraySpec("atom_relationships", "identity", np.ndarray,        2,    False, "[partner_atom_index, parent_atom_index] (int32, -1 = none)"),
    ArraySpec("atom_role",          "identity", np.ndarray,        None, False, "AtomRole enum / nmr_class (int32)"),
    ArraySpec("atom_flags",         "identity", np.ndarray,        10,   False, "[hybridisation, is_backbone, is_amide_H, is_alpha_H, is_methyl, is_aromatic_H, is_on_aromatic_residue, is_hbond_donor, is_hbond_acceptor, parent_is_sp2] (int32)"),
    ArraySpec("residue_context",    "identity", np.ndarray,        4,    False, "[prev_residue_type, next_residue_type, is_n_terminal, is_c_terminal] (int32)"),
    ArraySpec("ring_contributions", "identity", RingContributions, 59,   True,  "Per-(atom,ring) pair contributions"),
    ArraySpec("ring_geometry",      "identity", RingGeometry,      10,   True,  "Per-ring geometry reference"),

    # ── Biot-Savart (BiotSavartResult.cpp) ───────────────────────
    ArraySpec("bs_shielding",     "biot_savart", ShieldingTensor,  9,    True,  "BS ring current shielding"),
    ArraySpec("bs_per_type_T0",   "biot_savart", PerRingTypeT0,    8,    True,  "BS isotropic per ring type"),
    ArraySpec("bs_per_type_T2",   "biot_savart", PerRingTypeT2,    40,   True,  "BS T2 per ring type"),
    ArraySpec("bs_total_B",       "biot_savart", VectorField,      3,    True,  "BS total B-field vector"),
    ArraySpec("bs_ring_counts",   "biot_savart", RingCounts,       4,    True,  "Ring proximity counts (3/5/8/12 A)"),

    # ── Haigh-Mallion (HaighMallionResult.cpp) ───────────────────
    ArraySpec("hm_shielding",     "haigh_mallion", ShieldingTensor, 9,   True,  "HM ring current shielding"),
    ArraySpec("hm_per_type_T0",   "haigh_mallion", PerRingTypeT0,   8,   True,  "HM isotropic per ring type"),
    ArraySpec("hm_per_type_T2",   "haigh_mallion", PerRingTypeT2,   40,  True,  "HM T2 per ring type"),

    # ── Pi-Quadrupole (PiQuadrupoleResult.cpp) ───────────────────
    ArraySpec("pq_shielding",     "pi_quadrupole", ShieldingTensor, 9,   True,  "PQ EFG shielding"),
    ArraySpec("pq_per_type_T0",   "pi_quadrupole", PerRingTypeT0,   8,   True,  "PQ isotropic per ring type"),
    ArraySpec("pq_per_type_T2",   "pi_quadrupole", PerRingTypeT2,   40,  True,  "PQ T2 per ring type"),

    # ── Dispersion (DispersionResult.cpp) ────────────────────────
    ArraySpec("disp_shielding",   "dispersion", ShieldingTensor,   9,    True,  "Dispersion shielding (1/r^6)"),
    ArraySpec("disp_per_type_T0", "dispersion", PerRingTypeT0,     8,    True,  "Dispersion scalar per ring type"),
    ArraySpec("disp_per_type_T2", "dispersion", PerRingTypeT2,     40,   True,  "Dispersion T2 per ring type"),

    # ── Ring Susceptibility (RingSusceptibilityResult.cpp) ───────
    ArraySpec("ringchi_shielding","ring_susceptibility", ShieldingTensor, 9, True, "Ring susceptibility shielding"),

    # ── McConnell (McConnellResult.cpp) ──────────────────────────
    ArraySpec("mc_shielding",     "mcconnell", ShieldingTensor,    9,    True,  "McConnell bond anisotropy shielding"),
    ArraySpec("mc_category_T2",   "mcconnell", PerBondCategoryT2,  25,   True,  "McConnell T2 per bond category"),
    ArraySpec("mc_scalars",       "mcconnell", McConnellScalars,   6,    True,  "McConnell scalar sums + distances"),

    # ── Coulomb (CoulombResult.cpp) — optional via --no-coulomb ──
    ArraySpec("coulomb_shielding",      "coulomb", ShieldingTensor, 9,   False, "Coulomb E-field shielding"),
    ArraySpec("coulomb_E",              "coulomb", VectorField,     3,   False, "Coulomb total E-field"),
    ArraySpec("coulomb_efg_backbone",   "coulomb", EFGTensor,       9,   False, "Coulomb EFG backbone"),
    ArraySpec("coulomb_efg_aromatic",   "coulomb", EFGTensor,       9,   False, "Coulomb EFG aromatic"),
    ArraySpec("coulomb_scalars",        "coulomb", CoulombScalars,  4,   False, "Coulomb E-field scalars"),

    # ── H-Bond (HBondResult.cpp) ─────────────────────────────────
    ArraySpec("hbond_shielding",  "hbond", ShieldingTensor,        9,    True,  "H-bond shielding"),
    ArraySpec("hbond_scalars",    "hbond", HBondScalars,           3,    True,  "H-bond scalars"),

    # ── DSSP (DsspResult.cpp) ────────────────────────────────────
    ArraySpec("dssp_backbone",    "dssp", DsspScalars,             5,    True,  "DSSP backbone geometry"),
    ArraySpec("dssp_ss8",         "dssp", np.ndarray,              8,    False, "DSSP 8-class SS one-hot (H/G/I/E/B/T/S/C)"),
    ArraySpec("dssp_hbond_energy","dssp", np.ndarray,              4,    False, "DSSP H-bond energies (acc0/acc1/don0/don1)"),
    ArraySpec("dssp_chi",         "dssp", np.ndarray,              12,   False, "Chi1-4 cos/sin/exists (4 x 3 cols)"),

    # ── SASA (SasaResult.cpp) ───────────────────────────────────
    ArraySpec("atom_sasa",        "sasa", np.ndarray,              None, False, "Per-atom Shrake-Rupley SASA (A^2)"),
    ArraySpec("sasa_normal",      "sasa", VectorField,             3,    False, "SASA outward surface normal (unit vector)"),

    # ── Explicit water (WaterFieldResult.cpp) ───────────────────
    ArraySpec("water_efield",       "water_field", VectorField,    3,    False, "Water E-field total (V/A)"),
    ArraySpec("water_efield_first", "water_field", VectorField,    3,    False, "Water E-field first shell <3.5A (V/A)"),
    ArraySpec("water_efg",          "water_field", EFGTensor,      9,    False, "Water EFG total (SphericalTensor)"),
    ArraySpec("water_efg_first",    "water_field", EFGTensor,      9,    False, "Water EFG first shell (SphericalTensor)"),
    ArraySpec("water_shell_counts", "water_field", np.ndarray,     2,    False, "Water shell counts [n_first, n_second]"),

    # ── Hydration shell (HydrationShellResult.cpp) ──────────────
    ArraySpec("hydration_shell",    "hydration",   np.ndarray,     4,    False, "Hydration geometry [asymmetry, dipole_cos, ion_dist, ion_charge]"),

    # ── Hydration geometry — SASA-normal (HydrationGeometryResult.cpp) ─
    ArraySpec("water_polarization", "water_polarization", np.ndarray, 10, False, "Water polarisation [dipole(3), normal(3), asym, align, coher, count]"),

    # ── EEQ charges (EeqResult.cpp — Caldeweyher 2019) ─────────
    ArraySpec("eeq_charges",        "eeq",         np.ndarray,     None, False, "EEQ geometry-dependent charges (e)"),
    ArraySpec("eeq_cn",             "eeq",         np.ndarray,     None, False, "EEQ coordination number"),

    # ── GROMACS energy (GromacsEnergyResult.cpp) ────────────────
    ArraySpec("gromacs_energy",     "gromacs",     np.ndarray,     42,   False, "Per-frame energy (42 cols: electrostatic, bonded, VdW, thermo, box, virial, pressure tensor, T_group)"),

    # ── Bonded energy (BondedEnergyResult.cpp) ─────────────────
    ArraySpec("bonded_energy",      "bonded",      np.ndarray,      7,   False, "Per-atom bonded energy (bond,angle,UB,proper,improper,CMAP,total) kJ/mol"),

    # ── MOPAC core (MopacResult.cpp) ─────────────────────────────
    ArraySpec("mopac_charges",    "mopac_core", np.ndarray,        None, False, "MOPAC Mulliken charges"),
    ArraySpec("mopac_scalars",    "mopac_core", MopacScalars,      4,    False, "MOPAC per-atom scalars"),
    ArraySpec("mopac_bond_orders","mopac_core", BondOrders,        3,    False, "MOPAC Wiberg bond orders"),
    ArraySpec("mopac_global",     "mopac_core", MopacGlobal,       4,    False, "MOPAC graph-level scalars"),

    # ── MOPAC Coulomb (MopacCoulombResult.cpp) ───────────────────
    ArraySpec("mopac_coulomb_shielding",     "mopac_coulomb", ShieldingTensor, 9,  False, "MOPAC Coulomb shielding"),
    ArraySpec("mopac_coulomb_E",             "mopac_coulomb", VectorField,     3,  False, "MOPAC Coulomb E-field"),
    ArraySpec("mopac_coulomb_efg_backbone",  "mopac_coulomb", EFGTensor,       9,  False, "MOPAC Coulomb EFG backbone"),
    ArraySpec("mopac_coulomb_efg_aromatic",  "mopac_coulomb", EFGTensor,       9,  False, "MOPAC Coulomb EFG aromatic"),
    ArraySpec("mopac_coulomb_scalars",       "mopac_coulomb", CoulombScalars,  4,  False, "MOPAC Coulomb scalars"),

    # ── MOPAC McConnell (MopacMcConnellResult.cpp) ───────────────
    ArraySpec("mopac_mc_shielding",    "mopac_mcconnell", ShieldingTensor,   9,  False, "MOPAC McConnell shielding"),
    ArraySpec("mopac_mc_category_T2",  "mopac_mcconnell", PerBondCategoryT2, 25, False, "MOPAC McConnell T2 per category"),
    ArraySpec("mopac_mc_scalars",      "mopac_mcconnell", McConnellScalars,  6,  False, "MOPAC McConnell scalars"),

    # ── APBS (ApbsFieldResult.cpp) ───────────────────────────────
    ArraySpec("apbs_E",           "apbs", VectorField,             3,    False, "APBS solvated E-field"),
    ArraySpec("apbs_efg",         "apbs", EFGTensor,               9,    False, "APBS solvated EFG"),

    # ── Orca DFT (OrcaShieldingResult.cpp) ───────────────────────
    ArraySpec("orca_total",       "orca", ShieldingTensor,         9,    False, "Orca DFT total shielding"),
    ArraySpec("orca_diamagnetic", "orca", ShieldingTensor,         9,    False, "Orca DFT diamagnetic"),
    ArraySpec("orca_paramagnetic","orca", ShieldingTensor,         9,    False, "Orca DFT paramagnetic"),

    # ── Mutation delta (MutationDeltaResult.cpp) ─────────────────
    # Total-channel delta — preserved name (downstream consumers).
    ArraySpec("delta_shielding",                "delta", ShieldingTensor,    9,    False, "WT-ALA shielding delta (total channel)"),
    # 2026-04-26: per-channel deltas + WT/mutant raw shielding for residual
    # analysis. Paramagnetic isolates aromatic-ring + heavy-atom effects;
    # diamagnetic is the local-density baseline.
    ArraySpec("delta_shielding_diamagnetic",    "delta", ShieldingTensor,    9,    False, "WT-ALA shielding delta (diamagnetic channel)"),
    ArraySpec("delta_shielding_paramagnetic",   "delta", ShieldingTensor,    9,    False, "WT-ALA shielding delta (paramagnetic channel)"),
    ArraySpec("wt_shielding_total",             "delta", ShieldingTensor,    9,    False, "WT raw shielding (total channel)"),
    ArraySpec("wt_shielding_diamagnetic",       "delta", ShieldingTensor,    9,    False, "WT raw shielding (diamagnetic channel)"),
    ArraySpec("wt_shielding_paramagnetic",      "delta", ShieldingTensor,    9,    False, "WT raw shielding (paramagnetic channel)"),
    ArraySpec("mut_shielding_total",            "delta", ShieldingTensor,    9,    False, "Mutant raw shielding (total channel)"),
    ArraySpec("mut_shielding_diamagnetic",      "delta", ShieldingTensor,    9,    False, "Mutant raw shielding (diamagnetic channel)"),
    ArraySpec("mut_shielding_paramagnetic",     "delta", ShieldingTensor,    9,    False, "Mutant raw shielding (paramagnetic channel)"),
    ArraySpec("delta_scalars",                  "delta", DeltaScalars,       6,    False, "Delta metadata + match info"),
    ArraySpec("delta_apbs",                     "delta", DeltaAPBS,          12,   False, "APBS delta E + EFG"),
    ArraySpec("delta_ring_proximity",           "delta", DeltaRingProximity, None, False, "Removed ring geometry (variable cols)"),
    # 2026-04-26: per-WT-atom mutant index + per-mutation-site removed
    # rings. mutation_match is N rows (same as atom NPYs), -1 sentinel
    # for unmatched. mutation_removed_rings is R rows where R = total
    # rings removed across all mutation sites; absent when no aromatic
    # mutations occurred.
    ArraySpec("mutation_match",                 "delta", np.ndarray,     None, False, "WT-to-mutant atom mapping (int32, -1 unmatched)"),
    ArraySpec("mutation_removed_rings",         "delta", np.ndarray,     3,    False, "[residue_index, ring_index, ring_type] per removed ring (int32)"),

    # ── AIMNet2 (AIMNet2Result.cpp) ─────────────────────────────
    ArraySpec("aimnet2_charges",             "aimnet2", AIMNet2Charges,            None, False, "AIMNet2 Hirshfeld charges"),
    ArraySpec("aimnet2_aim",                 "aimnet2", AIMNet2AimEmbedding,       256,  False, "AIMNet2 256-dim electronic embedding"),
    ArraySpec("aimnet2_efg",                 "aimnet2", EFGTensor,                 9,    False, "AIMNet2 Coulomb EFG total"),
    ArraySpec("aimnet2_efg_aromatic",        "aimnet2", EFGTensor,                 9,    False, "AIMNet2 Coulomb EFG aromatic"),
    ArraySpec("aimnet2_efg_backbone",        "aimnet2", EFGTensor,                 9,    False, "AIMNet2 Coulomb EFG backbone"),
    # aimnet2_charge_sensitivity: REMOVED from calculator output (2026-04-12).
    # Perturbation approach deleted. Charge sensitivity is now computed
    # by GromacsFrameHandler as ensemble charge variance (Welford on
    # aimnet2_charge across frames) and optionally autograd on selected
    # frames. The SDK entry is kept for backward compatibility with old
    # extractions that have the file.
    ArraySpec("aimnet2_charge_sensitivity",  "aimnet2", AIMNet2ChargeSensitivity,  None, False, "LEGACY: per-atom charge sensitivity (old perturbation method)"),
]}
# fmt: on
