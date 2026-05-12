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
    ArraySpec("pos",              "identity",   VectorField,       3,    True,  "Atom positions (A)"),
    ArraySpec("element",          "identity",   np.ndarray,        None, True,  "Atomic number (int32)"),
    ArraySpec("residue_index",    "identity",   np.ndarray,        None, True,  "Residue index (int32)"),
    ArraySpec("residue_type",     "identity",   np.ndarray,        None, True,  "Residue type enum (int32)"),
    ArraySpec("ring_contributions","identity",  RingContributions, 59,   True,  "Per-(atom,ring) pair contributions"),
    ArraySpec("ring_geometry",    "identity",   RingGeometry,      10,   True,  "Per-ring geometry reference"),

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
    ArraySpec("delta_shielding",       "delta", ShieldingTensor,       9,    False, "WT-ALA shielding delta (total)"),
    ArraySpec("delta_scalars",         "delta", DeltaScalars,          6,    False, "Delta metadata + match info"),
    ArraySpec("delta_apbs",            "delta", DeltaAPBS,             12,   False, "APBS delta E + EFG"),
    ArraySpec("delta_ring_proximity",  "delta", DeltaRingProximity,    None, False, "Removed ring geometry (variable cols)"),
    # DFT shielding component decomposition: WT side, mut side, deltas;
    # diamagnetic and paramagnetic. sigma_total = sigma_dia + sigma_para;
    # the existing delta_shielding satisfies that identity at ORCA's
    # output precision (~1e-3 ppm). Stratifies mutation shifts by
    # physical mechanism.
    ArraySpec("wt_shielding_diamagnetic",     "delta", ShieldingTensor, 9, False, "WT diamagnetic shielding (matched, by WT atom row)"),
    ArraySpec("wt_shielding_paramagnetic",    "delta", ShieldingTensor, 9, False, "WT paramagnetic shielding (matched, by WT atom row)"),
    ArraySpec("mut_shielding_diamagnetic",    "delta", ShieldingTensor, 9, False, "mut diamagnetic shielding (matched, by WT atom row)"),
    ArraySpec("mut_shielding_paramagnetic",   "delta", ShieldingTensor, 9, False, "mut paramagnetic shielding (matched, by WT atom row)"),
    ArraySpec("delta_shielding_diamagnetic",  "delta", ShieldingTensor, 9, False, "WT - mut diamagnetic shielding delta"),
    ArraySpec("delta_shielding_paramagnetic", "delta", ShieldingTensor, 9, False, "WT - mut paramagnetic shielding delta"),

    # ── Per-atom invariant categorical record (CategoryInfoProjection.cpp) ──
    # Structured-dtype NPY (~31 fields). One-shot per protein, NOT per-
    # conformation. Wrapper is np.ndarray here (delivered as the raw
    # structured array) and `load()` in _protein.py wraps it as
    # CategoryInfo — a circular import would result if the catalog
    # referenced the wrapper class directly.
    ArraySpec("atoms_category_info",   "identity", np.ndarray, None, False, "Per-atom invariant categorical record (structured dtype)"),

    # ── AIMNet2 (AIMNet2Result.cpp) ─────────────────────────────
    # Required in production output per project_aimnet2_contract_20260426.
    # The 2026-04-26 contract was articulated in the memory entry but
    # the required=True flag was not landed in the catalog at the time;
    # landed 2026-05-04 alongside AIMNet2 wire-in to smoke tests.
    ArraySpec("aimnet2_charges",             "aimnet2", AIMNet2Charges,            None, True,  "AIMNet2 Hirshfeld charges"),
    ArraySpec("aimnet2_aim",                 "aimnet2", AIMNet2AimEmbedding,       256,  True,  "AIMNet2 256-dim electronic embedding"),
    ArraySpec("aimnet2_efg",                 "aimnet2", EFGTensor,                 9,    True,  "AIMNet2 Coulomb EFG total"),
    ArraySpec("aimnet2_efg_aromatic",        "aimnet2", EFGTensor,                 9,    True,  "AIMNet2 Coulomb EFG aromatic"),
    ArraySpec("aimnet2_efg_backbone",        "aimnet2", EFGTensor,                 9,    True,  "AIMNet2 Coulomb EFG backbone"),

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
    ArraySpec("aimnet2_polarisability",        "aimnet2", AIMNet2Polarisability, 3,    False, "AIMNet2 per-atom polarisability gradient (d(sum q_j^2)/d(r_i))"),
    ArraySpec("aimnet2_polarisability_scalar", "aimnet2", np.ndarray,            None, False, "AIMNet2 per-atom polarisability scalar (L2 norm of gradient)"),

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
    ArraySpec("pyramidalization",  "planar_geometry", np.ndarray, None, False, "Per-atom signed sp2 out-of-plane displacement (Å); 0 for non-planar atoms"),
    ArraySpec("omega_actual",      "planar_geometry", np.ndarray, None, False, "Per-residue ω (Cα-C-N-Cα to next), radians; NaN at C-term and X-Pro"),
    ArraySpec("omega_deviation",   "planar_geometry", np.ndarray, None, False, "Per-residue ω - π wrapped to (-π, π]; NaN where omega_actual is NaN"),
    ArraySpec("aromatic_chi2",     "planar_geometry", np.ndarray, None, False, "Per-aromatic-ring χ₂ (parent residue, radians); ring-flip observable per Akke-Weininger 2023"),
    ArraySpec("pucker_Q",          "planar_geometry", np.ndarray, None, False, "Per-saturated-ring Cremer-Pople puckering amplitude (Å); 5-rings only"),
    ArraySpec("pucker_theta",      "planar_geometry", np.ndarray, None, False, "Per-saturated-ring Cremer-Pople phase angle (degrees, [0, 360))"),
    ArraySpec("omega_is_xpro",     "planar_geometry", np.ndarray, None, False, "Per-residue mask: 1 where the bond into i+1 is X→Pro (cis/trans isomerism is real signal there, not 'non-planar amide' deviation)"),

    # ── Tripeptide DFT shielding ────────────────────────────────────
    # ProCS15 (Larsen 2015) tripeptide DFT lookup. σ_BB^i emitted on
    # backbone N/CA/C/O/H/HA and central-residue sidechain atoms per
    # the typed-identity-matched LarsenResidue model. Δσ_BB^{i±1}
    # neighbor correction emitted at the central residue's atoms per
    # Larsen Eq 3 cap-side reading. required=False because the
    # tensorcs15 DB is not available on every host.
    ArraySpec("tripeptide_bb_shielding",          "tripeptide", ShieldingTensor, 9,    False, "σ_BB^i — Mat3 (ppm) from typed-identity match against Larsen 2015 AXA tripeptide DFT row"),
    ArraySpec("tripeptide_bb_residual_vec",       "tripeptide", VectorField,     3,    False, "σ_BB^i match residual: aligned_dft - protein position; Vec3 ML feature (magnitude + direction)"),
    ArraySpec("tripeptide_bb_match_distance",     "tripeptide", np.ndarray,      None, False, "σ_BB^i match distance (Å) — magnitude of residual_vec"),
    ArraySpec("tripeptide_bb_method_tag",         "tripeptide", np.ndarray,      None, False, "DFT method discriminator: 0=none, 1=OPBE Gaussian (Larsen), 2=PBE ORCA (project SER regen)"),
    ArraySpec("tripeptide_neighbor_shielding",    "tripeptide", ShieldingTensor, 9,    False, "Δσ_BB^{i±1} — neighbour correction at residue i from i±1 cap reads (Larsen 2015 Eq 3)"),
    ArraySpec("tripeptide_neighbor_residual_vec_prev", "tripeptide", VectorField, 3,   False, "Δσ_BB^{i-1} match residual at the C-term ALA cap of (i-1)'s tripeptide; Vec3, NaN where i-1 direction had no contribution"),
    ArraySpec("tripeptide_neighbor_residual_vec_next", "tripeptide", VectorField, 3,   False, "Δσ_BB^{i+1} match residual at the N-term ALA cap of (i+1)'s tripeptide; Vec3, NaN where i+1 direction had no contribution"),

    # ────────────────────────────────────────────────────────────────
    # Larsen H-bond shielding contributions
    # (src/LarsenHBondShieldingResult.cpp). Phase 1 emits amide-H /
    # backbone-O subset (DSSP-resolved); Hα donors + sidechain
    # acceptors land in Phase 2. Methods accumulate side-by-side with
    # the kernel-form HBondResult — see feedback_methods_accumulate
    # memory entry. NPY layout is SphericalTensor-packed
    # (T0 + T1 + T2 = 9 columns) per HBondResult convention.
    # ────────────────────────────────────────────────────────────────
    ArraySpec("larsen_hbond_shielding",                  "larsen_hbond", ShieldingTensor, 9, False, "Σ Larsen H-bond contributions across all four Table 2 classes (1°HB + 2°HB + 1°HαB + 2°HαB) — ppm, lab frame"),
    ArraySpec("larsen_hbond_1pHB_shielding",             "larsen_hbond", ShieldingTensor, 9, False, "Δσ_1°HB per Larsen 2015 — primary amide-H donor contribution on donor residue i atoms"),
    ArraySpec("larsen_hbond_2pHB_shielding",             "larsen_hbond", ShieldingTensor, 9, False, "Δσ_2°HB per Larsen 2015 — secondary amide-H donor contribution on acceptor's residue i+1 atoms"),
    ArraySpec("larsen_hbond_1pHaB_shielding",            "larsen_hbond", ShieldingTensor, 9, False, "Δσ_1°HαB per Larsen 2015 — primary Hα donor contribution (Phase 2)"),
    ArraySpec("larsen_hbond_2pHaB_shielding",            "larsen_hbond", ShieldingTensor, 9, False, "Δσ_2°HαB per Larsen 2015 — secondary Hα donor contribution (Phase 2)"),
    ArraySpec("larsen_hbond_diagnostic_CB_shielding",    "larsen_hbond", ShieldingTensor, 9, False, "Cβ diagnostic readout — Larsen Table 2 says Cβ gets NO contribution; we emit anyway as reality check, should be near-zero"),
    ArraySpec("larsen_hbond_water_term",                 "larsen_hbond", np.ndarray,      None, False, "Δσ_w = 2.07 ppm isotropic on amide H atoms that DSSP detected as solvent-exposed (no H-bond partner)"),
    ArraySpec("larsen_hbond_count",                      "larsen_hbond", np.ndarray,      None, False, "Per-atom count of H-bond pairs that contributed under any of the four Table 2 classes (diagnostic CB does NOT count here)"),
]}
# fmt: on
