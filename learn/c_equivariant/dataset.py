"""
Calibration dataset using typed protein features.

Loads --mutant extraction output via load_protein(), uses C++ matched_mask
instead of Python atom matching, assembles kernel T2s and scalar features
from WT environment directly.

The delta framing means:
  - Ring kernels (BS, HM, PQ, Disp per-type T2) ARE the delta —
    ALA has no aromatic rings, so WT value = perturbation.
  - Non-ring kernels (MC, Coulomb, HBond etc) provide environment context
    for the scalar MLP to weight ring kernels correctly.
  - Target is DFT delta T2 (delta.shielding.T2), computed in C++.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein import load_protein, list_proteins
from features import RingType, BondCategory


# ── Kernel layout ────────────────────────────────────────────────────
#
#  0-7    BS per-ring-type T2   (8 ring types)
#  8-15   HM per-ring-type T2
# 16-23   Disp per-ring-type T2
# 24-28   MC bond-category T2   (5 categories)
# 29-33   MopacMC bond-category T2  (5 categories)
# 34      MC total T2
# 35      Coulomb total T2
# 36      HBond total T2
# 37      RingSusc total T2
# 38      MopacCoulomb total T2
# 39      MopacMC total T2
# --- EFG T2 kernels (anisotropic electric environment) ---
# 40      Coulomb EFG backbone T2
# 41      Coulomb EFG aromatic T2
# 42      MopacCoulomb EFG backbone T2
# 43      MopacCoulomb EFG aromatic T2
# 44      APBS EFG T2
# 45      Delta APBS EFG T2 (WT-ALA electrostatic change)
#
# PQ (Pople-Karplus) dropped: r < 0.03 across all ring types in analysis.

N_RING_TYPE_KERNELS = 24  # 3 calculators (BS, HM, Disp) × 8 ring types
N_CATEGORY_KERNELS = 10   # 2 calculators × 5 bond categories
N_TOTAL_KERNELS = 6       # MC, Coulomb, HBond, RingSusc, MopacCoulomb, MopacMC
N_EFG_KERNELS = 6         # EFG T2: backbone/aromatic × ff14SB/MOPAC + APBS + delta APBS
N_KERNELS = N_RING_TYPE_KERNELS + N_CATEGORY_KERNELS + N_TOTAL_KERNELS + N_EFG_KERNELS  # 46


def _build_features(protein_dir: Path):
    """Load one protein and return raw arrays.

    Returns None if protein lacks delta or has too few matched atoms.
    Kernel T2s are per-protein normalized (angular structure preserved,
    magnitude stripped).  Scalar normalization happens after stacking.
    """
    try:
        p = load_protein(protein_dir)
    except (FileNotFoundError, ValueError) as e:
        return None, str(e)

    if p.delta is None:
        return None, "no delta"

    mask = p.delta.scalars.matched_mask
    idx = np.where(mask)[0]
    if len(idx) < 10:
        return None, f"only {len(idx)} matched atoms"

    M = len(idx)

    # ── Target: DFT delta T2 from C++ ────────────────────────────
    target = p.delta.shielding.T2[idx]  # (M, 5)

    # ── Distance to nearest removed ring (for stratified evaluation) ──
    ring_dist = p.delta.scalars.nearest_removed_ring_dist[idx]  # (M,)

    # ── Kernel T2s ─────────────────────────────────────────────
    kernels = np.zeros((M, N_KERNELS, 5), dtype=np.float64)

    # Ring-type T2: BS, HM, Disp (PQ dropped — near-zero correlation)
    for t in range(8):
        kernels[:, t, :]      = p.biot_savart.per_type_T2.for_ring_type(RingType(t))[idx]
        kernels[:, 8 + t, :]  = p.haigh_mallion.per_type_T2.for_ring_type(RingType(t))[idx]
        kernels[:, 16 + t, :] = p.dispersion.per_type_T2.for_ring_type(RingType(t))[idx]

    for c in range(5):
        kernels[:, 24 + c, :] = p.mcconnell.category_T2.for_category(BondCategory(c))[idx]
    if p.mopac:
        for c in range(5):
            kernels[:, 29 + c, :] = p.mopac.mcconnell.category_T2.for_category(BondCategory(c))[idx]

    kernels[:, 34, :] = p.mcconnell.shielding.T2[idx]
    kernels[:, 35, :] = p.coulomb.shielding.T2[idx]
    kernels[:, 36, :] = p.hbond.shielding.T2[idx]
    kernels[:, 37, :] = p.ring_susceptibility.T2[idx]
    if p.mopac:
        kernels[:, 38, :] = p.mopac.coulomb.shielding.T2[idx]
        kernels[:, 39, :] = p.mopac.mcconnell.shielding.T2[idx]

    # EFG T2 kernels — anisotropic electric environment
    kernels[:, 40, :] = p.coulomb.efg_backbone.T2[idx]
    kernels[:, 41, :] = p.coulomb.efg_aromatic.T2[idx]
    if p.mopac:
        kernels[:, 42, :] = p.mopac.coulomb.efg_backbone.T2[idx]
        kernels[:, 43, :] = p.mopac.coulomb.efg_aromatic.T2[idx]
    if p.apbs:
        kernels[:, 44, :] = p.apbs.efg.T2[idx]
    if p.delta and p.delta.apbs is not None:
        kernels[:, 45, :] = p.delta.apbs.delta_efg.T2[idx]

    # Per-protein kernel normalization: isolates angular structure,
    # strips magnitude (which varies with ring count/distance/protein size).
    # Scalar divisor → preserves equivariance.
    # The normalization factors are preserved as scalar features so the MLP
    # knows the scale it was stripped of — bridging per-protein and global.
    kernel_scales = np.zeros(N_KERNELS, dtype=np.float64)
    for k in range(N_KERNELS):
        kstd = kernels[:, k, :].std()
        if kstd > 1e-10:
            kernels[:, k, :] /= kstd
            kernel_scales[k] = kstd

    # Broadcast per-protein kernel scales to all atoms in this protein
    s_kernel_scales = np.tile(kernel_scales, (M, 1))  # (M, N_KERNELS)

    # ── Scalar features ──────────────────────────────────────────
    # Element one-hot (H=1, C=6, N=7, O=8, S=16)
    elem = p.element[idx]
    s_elem = np.zeros((M, 5), dtype=np.float64)
    s_elem[:, 0] = (elem == 1).astype(float)
    s_elem[:, 1] = (elem == 6).astype(float)
    s_elem[:, 2] = (elem == 7).astype(float)
    s_elem[:, 3] = (elem == 8).astype(float)
    s_elem[:, 4] = (elem == 16).astype(float)

    # Residue type one-hot (20 amino acids)
    s_restype = np.zeros((M, 20), dtype=np.float64)
    for aa in range(20):
        s_restype[:, aa] = (p.residue_type[idx] == aa).astype(float)

    # Ring proximity counts (4)
    s_ring = p.biot_savart.ring_counts.data[idx] / 10.0

    # McConnell scalars (6)
    s_mc = p.mcconnell.scalars.data[idx].copy()
    s_mc[:, 4] = np.where(s_mc[:, 4] > 90, 0.0, 1.0 / (s_mc[:, 4] + 1e-3))
    s_mc[:, 5] = np.where(s_mc[:, 5] > 90, 0.0, 1.0 / (s_mc[:, 5] + 1e-3))

    # Coulomb scalars (4)
    s_coul = p.coulomb.scalars.data[idx]

    # HBond scalars (3)
    s_hb = p.hbond.scalars.data[idx].copy()
    s_hb[:, 0] = np.where(s_hb[:, 0] > 90, 0.0, 1.0 / (s_hb[:, 0] + 1e-3))

    # Delta metadata (4)
    s_delta = np.zeros((M, 4), dtype=np.float64)
    ds = p.delta.scalars
    s_delta[:, 0] = np.where(ds.nearest_removed_ring_dist[idx] > 90, 0.0,
                             1.0 / (ds.nearest_removed_ring_dist[idx] + 1e-3))
    s_delta[:, 1] = ds.delta_partial_charge[idx]
    s_delta[:, 2] = ds.delta_mopac_charge[idx]
    s_delta[:, 3] = ds.match_distance[idx]

    # MOPAC charge (1)
    s_mopac = np.zeros((M, 1), dtype=np.float64)
    if p.mopac:
        s_mopac[:, 0] = p.mopac.core.charges[idx]

    # T1 magnitudes (10)
    t1_mags = np.zeros((M, 10), dtype=np.float64)
    t1_mags[:, 0] = np.linalg.norm(p.biot_savart.shielding.T1[idx], axis=-1)
    t1_mags[:, 1] = np.linalg.norm(p.haigh_mallion.shielding.T1[idx], axis=-1)
    t1_mags[:, 2] = np.linalg.norm(p.pople_karplus.shielding.T1[idx], axis=-1)
    t1_mags[:, 3] = np.linalg.norm(p.dispersion.shielding.T1[idx], axis=-1)
    t1_mags[:, 4] = np.linalg.norm(p.mcconnell.shielding.T1[idx], axis=-1)
    t1_mags[:, 5] = np.linalg.norm(p.coulomb.shielding.T1[idx], axis=-1)
    t1_mags[:, 6] = np.linalg.norm(p.hbond.shielding.T1[idx], axis=-1)
    t1_mags[:, 7] = np.linalg.norm(p.ring_susceptibility.T1[idx], axis=-1)
    if p.mopac:
        t1_mags[:, 8] = np.linalg.norm(p.mopac.coulomb.shielding.T1[idx], axis=-1)
        t1_mags[:, 9] = np.linalg.norm(p.mopac.mcconnell.shielding.T1[idx], axis=-1)

    # Ring proximity detail — top-6 nearest removed rings (62% of proteins
    # have >3 rings; top-3 was discarding most of the ring environment).
    # Per ring: [mcconnell_factor, exp_decay, 1/dist, z, rho, theta] + n_rings
    TOP_K_RINGS = 6
    RING_COLS = 6
    s_ringprox = np.zeros((M, TOP_K_RINGS * RING_COLS + 1), dtype=np.float64)
    if p.delta.ring_proximity is not None:
        nr = p.delta.ring_proximity.n_removed_rings
        s_ringprox[:, -1] = nr
        dists_all = np.column_stack(
            [p.delta.ring_proximity.distance(i)[idx] for i in range(nr)]
        )
        for atom_j in range(M):
            order = np.argsort(dists_all[atom_j])
            for ki, ri in enumerate(order[:TOP_K_RINGS]):
                base = ki * RING_COLS
                aj = idx[atom_j]
                d = p.delta.ring_proximity.distance(ri)[aj]
                s_ringprox[atom_j, base]     = p.delta.ring_proximity.mcconnell_factor(ri)[aj]
                s_ringprox[atom_j, base + 1] = p.delta.ring_proximity.exp_decay(ri)[aj]
                s_ringprox[atom_j, base + 2] = 1.0 / (d + 1e-3)
                s_ringprox[atom_j, base + 3] = p.delta.ring_proximity.z(ri)[aj]
                s_ringprox[atom_j, base + 4] = p.delta.ring_proximity.rho(ri)[aj]
                s_ringprox[atom_j, base + 5] = p.delta.ring_proximity.theta(ri)[aj]

    # Per-atom bond order features (4)
    s_bond = np.zeros((M, 4), dtype=np.float64)
    if p.mopac:
        bo_dense = p.mopac.core.bond_orders.to_dense(p.n_atoms)
        per_atom = bo_dense[idx]
        s_bond[:, 0] = per_atom.max(axis=1)                     # max bond order
        s_bond[:, 1] = (per_atom > 0.01).sum(axis=1).astype(float)  # n_bonds
        s_bond[:, 2] = per_atom.sum(axis=1)                     # total bond order
        s_bond[:, 3] = (per_atom > 1.2).sum(axis=1).astype(float)   # n_aromatic (pi)

    # Mutation identity — which aromatic residue was removed (4 one-hot)
    # Inferred from which ring-type kernels are non-zero in this protein.
    # PHE=has PHE rings, TYR=has TYR rings, TRP=has TRP rings, HIS=has HIS/HID/HIE.
    s_mutid = np.zeros((M, 4), dtype=np.float64)
    # Check which ring types have nonzero BS per-type T0
    bs_t0 = p.biot_savart.per_type_T0.data[idx]  # (M, 8)
    protein_ring_active = np.abs(bs_t0).sum(axis=0) > 1e-10  # (8,)
    if protein_ring_active[RingType.PHE]:
        s_mutid[:, 0] = 1.0
    if protein_ring_active[RingType.TYR]:
        s_mutid[:, 1] = 1.0
    if any(protein_ring_active[t] for t in [RingType.TRP_benzene,
           RingType.TRP_pyrrole, RingType.TRP_perimeter]):
        s_mutid[:, 2] = 1.0
    if any(protein_ring_active[t] for t in [RingType.HIS,
           RingType.HID, RingType.HIE]):
        s_mutid[:, 3] = 1.0

    # Per-ring-type T0 from all calculators — how strongly each calculator
    # sees each ring type at this atom (scalar context for weighting T2 kernels)
    s_t0_bs = p.biot_savart.per_type_T0.data[idx]       # (M, 8)
    s_t0_hm = p.haigh_mallion.per_type_T0.data[idx]     # (M, 8)
    s_t0_disp = p.dispersion.per_type_T0.data[idx]      # (M, 8)

    # MOPAC electronic structure detail
    s_mopac_elec = np.zeros((M, 2), dtype=np.float64)
    if p.mopac:
        s_mopac_elec[:, 0] = p.mopac.core.scalars.s_pop[idx]
        s_mopac_elec[:, 1] = p.mopac.core.scalars.p_pop[idx]

    # T0 (isotropic) magnitudes per calculator — large T0 with small T2 means
    # strong isotropic effect with weak angular structure, different physics.
    s_t0_mags = np.zeros((M, 6), dtype=np.float64)
    s_t0_mags[:, 0] = np.abs(p.mcconnell.shielding.isotropic[idx])
    s_t0_mags[:, 1] = np.abs(p.coulomb.shielding.isotropic[idx])
    s_t0_mags[:, 2] = np.abs(p.biot_savart.shielding.isotropic[idx])
    s_t0_mags[:, 3] = np.abs(p.hbond.shielding.isotropic[idx])
    if p.mopac:
        s_t0_mags[:, 4] = np.abs(p.mopac.coulomb.shielding.isotropic[idx])
        s_t0_mags[:, 5] = np.abs(p.mopac.mcconnell.shielding.isotropic[idx])

    # MOPAC Coulomb scalars — parallel to s_coul but with PM7 charges
    s_mopac_coul = np.zeros((M, 4), dtype=np.float64)
    if p.mopac:
        s_mopac_coul = p.mopac.coulomb.scalars.data[idx]

    # MOPAC McConnell scalars — parallel to s_mc but with Wiberg bond orders
    s_mopac_mc = np.zeros((M, 6), dtype=np.float64)
    if p.mopac:
        s_mopac_mc = p.mopac.mcconnell.scalars.data[idx].copy()
        s_mopac_mc[:, 4] = np.where(s_mopac_mc[:, 4] > 90, 0.0,
                                     1.0 / (s_mopac_mc[:, 4] + 1e-3))
        s_mopac_mc[:, 5] = np.where(s_mopac_mc[:, 5] > 90, 0.0,
                                     1.0 / (s_mopac_mc[:, 5] + 1e-3))

    scalars = np.concatenate([
        s_elem,       # 5
        s_restype,    # 20
        s_ring,       # 4
        s_mc,         # 6
        s_coul,       # 4
        s_hb,         # 3
        s_delta,      # 4
        s_mopac,      # 1
        t1_mags,      # 10
        s_ringprox,   # 37
        s_bond,       # 4
        s_mutid,      # 4
        s_t0_bs,      # 8
        s_t0_hm,      # 8
        s_t0_disp,    # 8
        s_mopac_elec, # 2
        s_t0_mags,    # 6
        s_mopac_coul,    # 4
        s_mopac_mc,      # 6
        s_kernel_scales, # 46 — per-protein norm factors (magnitude the kernels were stripped of)
    ], axis=1)           # total: 190

    return {
        "scalars": scalars,
        "kernels": kernels,
        "target": target,
        "ring_dist": ring_dist,
    }, None


N_SCALAR_FEATURES = 190

# Columns that are one-hot / categorical — skip z-score for these
_ONEHOT_COLS = list(range(25)) + list(range(98, 102))  # 5 elem + 20 restype + 4 mutid


@dataclass
class NormStats:
    """Normalization statistics computed from training set, shared with val."""
    target_std: float
    scalar_mean: np.ndarray   # (N_SCALAR_FEATURES,)
    scalar_std: np.ndarray    # (N_SCALAR_FEATURES,)
    gate_threshold: np.ndarray  # (N_KERNELS,) — per-kernel median magnitude


class CalibrationDataset(Dataset):
    """WT-ALA delta T2 dataset from typed protein features.

    Target: delta_shielding T2 (5 components) from MutationDeltaResult.
    Inputs: 46 kernel T2s (L=2) + 138 scalar features (L=0).

    Normalization:
      - Kernels: per-protein per-kernel std (isolates angular structure)
      - Scalars: z-score (skip one-hot columns)
      - Targets: divide by target_std

    Pass norm_stats from training set to val set so both use the same scale.
    """

    def __init__(self, protein_ids: list[str], features_dir: Path,
                 norm_stats: NormStats | None = None):
        raw_scalars = []
        raw_kernels = []
        raw_targets = []
        raw_ring_dist = []
        protein_ids_loaded = []
        protein_boundaries = []  # (start, end) per protein for per-protein R²
        loaded = 0
        skipped = 0
        offset = 0

        for pid in protein_ids:
            result, err = _build_features(features_dir / pid)
            if result is None:
                skipped += 1
                if err:
                    print(f"  skip {pid}: {err}")
                continue
            n = len(result["target"])
            raw_scalars.append(result["scalars"])
            raw_kernels.append(result["kernels"])
            raw_targets.append(result["target"])
            raw_ring_dist.append(result["ring_dist"])
            protein_ids_loaded.append(pid)
            protein_boundaries.append((offset, offset + n))
            offset += n
            loaded += 1

        print(f"  loaded {loaded} proteins, skipped {skipped}")

        if not raw_scalars:
            raise ValueError("No proteins loaded — check features_dir and extraction")

        self.protein_ids = protein_ids_loaded
        self.protein_boundaries = protein_boundaries

        scalars_np = np.vstack(raw_scalars)
        kernels_np = np.vstack(raw_kernels)
        targets_np = np.vstack(raw_targets)
        self.ring_dist = torch.tensor(
            np.concatenate(raw_ring_dist), dtype=torch.float32)

        # Kernels already per-protein normalized in _build_features —
        # no global kernel normalization needed.

        if norm_stats is None:
            # Training set: compute stats
            target_std = float(targets_np.std())

            scalar_mean = scalars_np.mean(axis=0)
            scalar_std = scalars_np.std(axis=0)
            scalar_mean[_ONEHOT_COLS] = 0.0
            scalar_std[_ONEHOT_COLS] = 1.0
            scalar_std = np.maximum(scalar_std, 1e-8)

            # Per-kernel gate threshold: median magnitude of nonzero atoms.
            # Each kernel's noise floor is set by its own distribution —
            # whisper kernels get a low threshold, shouters get a high one.
            kernel_mags = np.linalg.norm(kernels_np, axis=-1)  # (N, K)
            gate_threshold = np.ones(N_KERNELS)
            for k in range(N_KERNELS):
                nonzero = kernel_mags[:, k][kernel_mags[:, k] > 1e-8]
                if len(nonzero) > 100:
                    gate_threshold[k] = float(np.median(nonzero))

            self.norm_stats = NormStats(
                target_std=target_std,
                scalar_mean=scalar_mean,
                scalar_std=scalar_std,
                gate_threshold=gate_threshold,
            )
        else:
            self.norm_stats = norm_stats

        ns = self.norm_stats

        # Apply scalar z-score
        scalars_np = (scalars_np - ns.scalar_mean) / ns.scalar_std

        # Apply target scaling
        targets_np = targets_np / ns.target_std

        self.scalars = torch.tensor(scalars_np, dtype=torch.float32)
        self.kernels = torch.tensor(kernels_np, dtype=torch.float32)
        self.targets = torch.tensor(targets_np, dtype=torch.float32)
        self.target_std = torch.tensor(ns.target_std, dtype=torch.float32)

    def to(self, device):
        self.scalars = self.scalars.to(device)
        self.kernels = self.kernels.to(device)
        self.targets = self.targets.to(device)
        self.target_std = self.target_std.to(device)
        self.ring_dist = self.ring_dist.to(device)
        return self

    def __len__(self):
        return len(self.scalars)

    def __getitem__(self, idx):
        return self.scalars[idx], self.kernels[idx], self.targets[idx]


# ── Evaluation ───────────────────────────────────────────────────────

# Distance bands for stratified R² (Angstroms from nearest removed ring)
DIST_BANDS = [(0, 4), (4, 8), (8, 12), (12, 999)]
DIST_LABELS = ["0-4A", "4-8A", "8-12A", "12+A"]


def compute_r2(pred: torch.Tensor, tgt: torch.Tensor,
               ring_dist: torch.Tensor | None = None) -> dict:
    """Per-component, overall, and distance-stratified R² for T2 predictions.

    Args:
        pred: (N, 5) predicted T2 in ppm
        tgt:  (N, 5) target T2 in ppm
        ring_dist: (N,) distance to nearest removed ring (for stratification)
    """
    result = {}

    # Per-component R²
    for i, m in enumerate([-2, -1, 0, 1, 2]):
        p, t = pred[:, i], tgt[:, i]
        ss_res = torch.sum((t - p) ** 2).item()
        ss_tot = torch.sum((t - t.mean()) ** 2).item()
        result[f"r2_m{m:+d}"] = 1.0 - ss_res / max(ss_tot, 1e-12)

    result["r2_mean"] = np.mean([result[f"r2_m{m:+d}"] for m in [-2, -1, 0, 1, 2]])

    # Flat R² (per-component mean subtracted — proper for T2)
    ss_res_all = torch.sum((tgt - pred) ** 2).item()
    ss_tot_all = torch.sum((tgt - tgt.mean(dim=0, keepdim=True)) ** 2).item()
    result["r2_flat"] = 1.0 - ss_res_all / max(ss_tot_all, 1e-12)

    # RMSE per atom
    per_atom_err = torch.sqrt(torch.sum((tgt - pred) ** 2, dim=1))
    result["rmse_per_atom_ppm"] = per_atom_err.mean().item()

    # Distance-stratified R² — where the model actually matters
    if ring_dist is not None:
        for (lo, hi), label in zip(DIST_BANDS, DIST_LABELS):
            mask = (ring_dist >= lo) & (ring_dist < hi)
            n = mask.sum().item()
            if n < 10:
                result[f"r2_{label}"] = float("nan")
                result[f"n_{label}"] = n
                continue
            p_band = pred[mask]
            t_band = tgt[mask]
            ss_res = torch.sum((t_band - p_band) ** 2).item()
            ss_tot = torch.sum(
                (t_band - t_band.mean(dim=0, keepdim=True)) ** 2).item()
            result[f"r2_{label}"] = 1.0 - ss_res / max(ss_tot, 1e-12)
            result[f"n_{label}"] = n

    return result


def compute_naive_baselines(ds: CalibrationDataset) -> dict:
    """Physics baselines that require no learning.

    Kernels are per-protein normalized (unit variance per kernel per protein),
    targets are scaled by target_std.  The baselines test angular structure:
    does the kernel T2 *direction* match the DFT delta T2 direction?

    A positive R² here means the raw kernels already point the right way —
    the mixer just needs to refine weights.
    """
    result = {}
    # Work in scaled space (same as model training)
    tgt = ds.targets  # (N, 5) scaled
    kernels = ds.kernels  # (N, 48, 5) per-protein normalized

    ss_tot = torch.sum((tgt - tgt.mean(dim=0, keepdim=True)) ** 2).item()

    # Baseline 1: unweighted sum of ring kernels only (0-31)
    ring_sum = kernels[:, :32, :].sum(dim=1)
    ss_res = torch.sum((tgt - ring_sum) ** 2).item()
    result["r2_ring_sum"] = 1.0 - ss_res / max(ss_tot, 1e-12)

    # Baseline 2: unweighted sum of ALL kernels
    all_sum = kernels.sum(dim=1)
    ss_res = torch.sum((tgt - all_sum) ** 2).item()
    result["r2_all_sum"] = 1.0 - ss_res / max(ss_tot, 1e-12)

    # Baseline 3: just BS ring kernels (0-7) — the dominant calculator
    bs_sum = kernels[:, :8, :].sum(dim=1)
    ss_res = torch.sum((tgt - bs_sum) ** 2).item()
    result["r2_bs_only"] = 1.0 - ss_res / max(ss_tot, 1e-12)

    # Baseline 4: optimal single scalar per kernel (ridge, no MLP)
    # Fits w_k for each kernel: pred = sum_k w_k * kernel_k
    # This is the Level A result — what linear mixing achieves
    X = kernels.reshape(kernels.shape[0], -1)  # (N, 48*5)
    y = tgt  # (N, 5)
    try:
        lam = 1e-3
        XtX = X.T @ X + lam * torch.eye(X.shape[1], device=X.device)
        Xty = X.T @ y
        w = torch.linalg.solve(XtX, Xty)
        ridge_pred = X @ w
        ss_res = torch.sum((tgt - ridge_pred) ** 2).item()
        result["r2_ridge"] = 1.0 - ss_res / max(ss_tot, 1e-12)
    except Exception:
        result["r2_ridge"] = float("nan")

    # Distance-stratified for ring_sum
    if ds.ring_dist is not None:
        for (lo, hi), label in zip(DIST_BANDS, DIST_LABELS):
            mask = (ds.ring_dist >= lo) & (ds.ring_dist < hi)
            n = mask.sum().item()
            if n < 10:
                result[f"ring_sum_{label}"] = float("nan")
                continue
            t_band = tgt[mask]
            p_band = ring_sum[mask]
            ss_res = torch.sum((t_band - p_band) ** 2).item()
            ss_tot_b = torch.sum(
                (t_band - t_band.mean(dim=0, keepdim=True)) ** 2).item()
            result[f"ring_sum_{label}"] = 1.0 - ss_res / max(ss_tot_b, 1e-12)

    return result


def compute_per_protein_r2(pred: torch.Tensor, tgt: torch.Tensor,
                           boundaries: list[tuple[int, int]],
                           protein_ids: list[str]) -> dict[str, float]:
    """R² per protein — identifies outliers and bad DFT."""
    result = {}
    for pid, (lo, hi) in zip(protein_ids, boundaries):
        p, t = pred[lo:hi], tgt[lo:hi]
        if hi - lo < 10:
            result[pid] = float("nan")
            continue
        ss_res = torch.sum((t - p) ** 2).item()
        ss_tot = torch.sum((t - t.mean(dim=0, keepdim=True)) ** 2).item()
        result[pid] = 1.0 - ss_res / max(ss_tot, 1e-12)
    return result
