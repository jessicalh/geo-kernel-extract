"""
Calibration dataset: WT-ALA delta T2 from typed protein features.

Uses nmr_extract SDK for all data access.  The delta framing means:
  - Ring kernels ARE the delta (ALA has no rings, so WT value = perturbation)
  - Non-ring kernels provide environment context for the MLP
  - Target is DFT delta T2 from MutationDeltaResult

Kernel layout (46 T2 vectors, L=2):
    0–7    BS per-ring-type         8–15  HM per-ring-type
   16–23   Disp per-ring-type      24–28  MC bond-category
   29–33   MopacMC bond-category   34–39  calculator totals
   40–45   EFG T2 (backbone/aromatic × ff14SB/MOPAC, APBS, delta APBS)
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "python"))
from nmr_extract import load, RingType, BondCategory

N_KERNELS = 46
N_SCALAR_FEATURES = 190

# Columns that are one-hot / categorical — skip z-score
_ELEMENT_COLS = list(range(5))            # 0–4
_RESTYPE_COLS = list(range(5, 25))        # 5–24
_MUTID_COLS   = list(range(98, 102))      # 98–101
_ONEHOT_COLS  = _ELEMENT_COLS + _RESTYPE_COLS + _MUTID_COLS

_ELEMENTS = [1, 6, 7, 8, 16]  # H, C, N, O, S
_TOP_K_RINGS = 6


# ── Kernel assembly ─────────────────────────────────────────────────

def _assemble_kernels(p, idx):
    """Build (M, 46, 5) kernel T2 array from protein features."""
    M = len(idx)
    K = np.zeros((M, N_KERNELS, 5), dtype=np.float64)

    # Ring-type T2: BS(0–7), HM(8–15), Disp(16–23)
    K[:, 0:8, :]   = p.biot_savart.per_type_T2.as_block()[idx]
    K[:, 8:16, :]  = p.haigh_mallion.per_type_T2.as_block()[idx]
    K[:, 16:24, :] = p.dispersion.per_type_T2.as_block()[idx]

    # Bond-category T2: MC(24–28), MopacMC(29–33)
    K[:, 24:29, :] = p.mcconnell.category_T2.as_block()[idx]
    if p.mopac:
        K[:, 29:34, :] = p.mopac.mcconnell.category_T2.as_block()[idx]

    # Calculator totals (34–39)
    K[:, 34, :] = p.mcconnell.shielding.T2[idx]
    K[:, 35, :] = p.coulomb.shielding.T2[idx]
    K[:, 36, :] = p.hbond.shielding.T2[idx]
    K[:, 37, :] = p.ring_susceptibility.T2[idx]
    if p.mopac:
        K[:, 38, :] = p.mopac.coulomb.shielding.T2[idx]
        K[:, 39, :] = p.mopac.mcconnell.shielding.T2[idx]

    # EFG T2 (40–45)
    K[:, 40, :] = p.coulomb.efg_backbone.T2[idx]
    K[:, 41, :] = p.coulomb.efg_aromatic.T2[idx]
    if p.mopac:
        K[:, 42, :] = p.mopac.coulomb.efg_backbone.T2[idx]
        K[:, 43, :] = p.mopac.coulomb.efg_aromatic.T2[idx]
    if p.apbs:
        K[:, 44, :] = p.apbs.efg.T2[idx]
    if p.delta and p.delta.apbs is not None:
        K[:, 45, :] = p.delta.apbs.delta_efg.T2[idx]

    return K


def _normalize_kernels(kernels):
    """Per-protein normalization.  Returns (normalized_kernels, scale_per_kernel)."""
    scales = np.zeros(N_KERNELS, dtype=np.float64)
    for k in range(N_KERNELS):
        s = kernels[:, k, :].std()
        if s > 1e-10:
            kernels[:, k, :] /= s
            scales[k] = s
    return kernels, scales


# ── Scalar assembly ──────────────────────────────────────────────────

def _one_hot(values, n_classes):
    """One-hot encode integer values."""
    M = len(values)
    oh = np.zeros((M, n_classes), dtype=np.float64)
    for i, v in enumerate(values):
        if 0 <= v < n_classes:
            oh[i, v] = 1.0
    return oh


def _sentinel_to_inv(arr, sentinel=90.0):
    """Convert sentinel-flagged distance to 1/d, zero for sentinels."""
    return np.where(arr > sentinel, 0.0, 1.0 / (arr + 1e-3))


def _assemble_scalars(p, idx, kernel_scales):
    """Build (M, N_SCALAR_FEATURES) scalar array.

    Groups are assembled in order matching _ONEHOT_COLS indices.
    """
    M = len(idx)
    blocks = []

    # ── Identity (29) ────────────────────────────────────────────
    # Element one-hot (5)
    blocks.append(_one_hot([_ELEMENTS.index(e) if e in _ELEMENTS else -1
                            for e in p.element[idx]], 5))
    # Residue type one-hot (20)
    blocks.append(_one_hot(p.residue_type[idx], 20))
    # Ring proximity counts (4)
    blocks.append(p.biot_savart.ring_counts.data[idx] / 10.0)

    # ── Calculator scalars (17) ──────────────────────────────────
    # McConnell (6): sums + 1/nearest_dist
    mc = p.mcconnell.scalars.data[idx].copy()
    mc[:, 4] = _sentinel_to_inv(mc[:, 4])
    mc[:, 5] = _sentinel_to_inv(mc[:, 5])
    blocks.append(mc)

    # Coulomb (4)
    blocks.append(p.coulomb.scalars.data[idx])

    # HBond (3): nearest_dist → 1/d
    hb = p.hbond.scalars.data[idx].copy()
    hb[:, 0] = _sentinel_to_inv(hb[:, 0])
    blocks.append(hb)

    # ── Delta metadata (4) ───────────────────────────────────────
    ds = p.delta.scalars
    delta = np.column_stack([
        _sentinel_to_inv(ds.nearest_removed_ring_dist[idx]),
        ds.delta_partial_charge[idx],
        ds.delta_mopac_charge[idx],
        ds.match_distance[idx],
    ])
    blocks.append(delta)

    # ── MOPAC charge (1) ─────────────────────────────────────────
    mopac_q = np.zeros((M, 1), dtype=np.float64)
    if p.mopac:
        mopac_q[:, 0] = p.mopac.core.charges[idx]
    blocks.append(mopac_q)

    # ── T1 magnitudes (10) ───────────────────────────────────────
    shieldings = [
        p.biot_savart.shielding, p.haigh_mallion.shielding,
        p.pi_quadrupole.shielding, p.dispersion.shielding,
        p.mcconnell.shielding, p.coulomb.shielding,
        p.hbond.shielding, p.ring_susceptibility,
    ]
    t1 = np.zeros((M, 10), dtype=np.float64)
    for i, s in enumerate(shieldings):
        t1[:, i] = np.linalg.norm(s.T1[idx], axis=-1)
    if p.mopac:
        t1[:, 8] = np.linalg.norm(p.mopac.coulomb.shielding.T1[idx], axis=-1)
        t1[:, 9] = np.linalg.norm(p.mopac.mcconnell.shielding.T1[idx], axis=-1)
    blocks.append(t1)

    # ── Per-ring geometry: top-6 nearest rings (37) ──────────────
    ringprox = np.zeros((M, _TOP_K_RINGS * 6 + 1), dtype=np.float64)
    if p.ring_contributions is not None:
        rc = p.ring_contributions
        for j in range(M):
            aj = idx[j]
            atom_rc = rc.for_atom(aj)
            if atom_rc.n_pairs == 0:
                continue
            order = np.argsort(atom_rc.distance)
            for ki, ri in enumerate(order[:_TOP_K_RINGS]):
                base = ki * 6
                ringprox[j, base]     = atom_rc.mcconnell_factor[ri]
                ringprox[j, base + 1] = atom_rc.exp_decay[ri]
                ringprox[j, base + 2] = 1.0 / (atom_rc.distance[ri] + 1e-3)
                ringprox[j, base + 3] = atom_rc.z[ri]
                ringprox[j, base + 4] = atom_rc.rho[ri]
                ringprox[j, base + 5] = atom_rc.theta[ri]
            ringprox[j, -1] = atom_rc.n_pairs
    elif p.delta and p.delta.ring_proximity is not None:
        # Fallback for old extractions without ring_contributions
        rp = p.delta.ring_proximity
        nr = rp.n_removed_rings
        ringprox[:, -1] = nr
        dists_all = np.column_stack([rp.distance(i)[idx] for i in range(nr)])
        for j in range(M):
            order = np.argsort(dists_all[j])
            for ki, ri in enumerate(order[:_TOP_K_RINGS]):
                base = ki * 6
                aj = idx[j]
                ringprox[j, base]     = rp.mcconnell_factor(ri)[aj]
                ringprox[j, base + 1] = rp.exp_decay(ri)[aj]
                ringprox[j, base + 2] = 1.0 / (rp.distance(ri)[aj] + 1e-3)
                ringprox[j, base + 3] = rp.z(ri)[aj]
                ringprox[j, base + 4] = rp.rho(ri)[aj]
                ringprox[j, base + 5] = rp.theta(ri)[aj]
    blocks.append(ringprox)

    # ── Bond orders (4) ──────────────────────────────────────────
    bond = np.zeros((M, 4), dtype=np.float64)
    if p.mopac:
        bo = p.mopac.core.bond_orders.to_dense(p.n_atoms)[idx]
        bond[:, 0] = bo.max(axis=1)
        bond[:, 1] = (bo > 0.01).sum(axis=1).astype(float)
        bond[:, 2] = bo.sum(axis=1)
        bond[:, 3] = (bo > 1.2).sum(axis=1).astype(float)
    blocks.append(bond)

    # ── Mutation identity (4 one-hot) ────────────────────────────
    mutid = np.zeros((M, 4), dtype=np.float64)
    bs_t0 = p.biot_savart.per_type_T0.data[idx]
    active = np.abs(bs_t0).sum(axis=0) > 1e-10
    if active[RingType.PHE]:       mutid[:, 0] = 1.0
    if active[RingType.TYR]:       mutid[:, 1] = 1.0
    if any(active[t] for t in [RingType.TRP_benzene,
           RingType.TRP_pyrrole, RingType.TRP_perimeter]):
                                   mutid[:, 2] = 1.0
    if any(active[t] for t in [RingType.HIS, RingType.HID, RingType.HIE]):
                                   mutid[:, 3] = 1.0
    blocks.append(mutid)

    # ── Per-ring-type T0: isotropic context (24) ─────────────────
    blocks.append(p.biot_savart.per_type_T0.data[idx])    # 8
    blocks.append(p.haigh_mallion.per_type_T0.data[idx])  # 8
    blocks.append(p.dispersion.per_type_T0.data[idx])     # 8

    # ── MOPAC electronic structure (2) ───────────────────────────
    mopac_elec = np.zeros((M, 2), dtype=np.float64)
    if p.mopac:
        mopac_elec[:, 0] = p.mopac.core.scalars.s_pop[idx]
        mopac_elec[:, 1] = p.mopac.core.scalars.p_pop[idx]
    blocks.append(mopac_elec)

    # ── T0 isotropic magnitudes (6) ──────────────────────────────
    t0 = np.zeros((M, 6), dtype=np.float64)
    t0[:, 0] = np.abs(p.mcconnell.shielding.isotropic[idx])
    t0[:, 1] = np.abs(p.coulomb.shielding.isotropic[idx])
    t0[:, 2] = np.abs(p.biot_savart.shielding.isotropic[idx])
    t0[:, 3] = np.abs(p.hbond.shielding.isotropic[idx])
    if p.mopac:
        t0[:, 4] = np.abs(p.mopac.coulomb.shielding.isotropic[idx])
        t0[:, 5] = np.abs(p.mopac.mcconnell.shielding.isotropic[idx])
    blocks.append(t0)

    # ── MOPAC Coulomb scalars (4) ────────────────────────────────
    mopac_coul = np.zeros((M, 4), dtype=np.float64)
    if p.mopac:
        mopac_coul = p.mopac.coulomb.scalars.data[idx]
    blocks.append(mopac_coul)

    # ── MOPAC McConnell scalars (6) ──────────────────────────────
    mopac_mc = np.zeros((M, 6), dtype=np.float64)
    if p.mopac:
        mopac_mc = p.mopac.mcconnell.scalars.data[idx].copy()
        mopac_mc[:, 4] = _sentinel_to_inv(mopac_mc[:, 4])
        mopac_mc[:, 5] = _sentinel_to_inv(mopac_mc[:, 5])
    blocks.append(mopac_mc)

    # ── Kernel normalization scales (46) ─────────────────────────
    # Per-protein constants: how strong each kernel was before normalization.
    # Bridges per-protein and global by telling the MLP what was stripped.
    blocks.append(np.tile(kernel_scales, (M, 1)))

    scalars = np.concatenate(blocks, axis=1)
    assert scalars.shape[1] == N_SCALAR_FEATURES, \
        f"Expected {N_SCALAR_FEATURES} scalars, got {scalars.shape[1]}"
    return scalars


# ── Per-protein feature builder ──────────────────────────────────────

def _build_features(protein_dir: Path):
    """Load one protein → (kernels, scalars, target, ring_dist) or None."""
    try:
        p = load(protein_dir)
    except (FileNotFoundError, ValueError) as e:
        return None, str(e)

    if p.delta is None:
        return None, "no delta"

    idx = np.where(p.delta.scalars.matched_mask)[0]
    if len(idx) < 10:
        return None, f"only {len(idx)} matched atoms"

    kernels = _assemble_kernels(p, idx)
    kernels, kernel_scales = _normalize_kernels(kernels)
    scalars = _assemble_scalars(p, idx, kernel_scales)
    target = p.delta.shielding.T2[idx]
    ring_dist = p.delta.scalars.nearest_removed_ring_dist[idx]

    return {"scalars": scalars, "kernels": kernels,
            "target": target, "ring_dist": ring_dist}, None


# ── Dataset ──────────────────────────────────────────────────────────

@dataclass
class NormStats:
    """Normalization statistics computed from training set, shared with val."""
    target_std: float
    scalar_mean: np.ndarray     # (N_SCALAR_FEATURES,)
    scalar_std: np.ndarray      # (N_SCALAR_FEATURES,)
    gate_threshold: np.ndarray  # (N_KERNELS,)


class CalibrationDataset(Dataset):
    """WT-ALA delta T2 dataset.

    Inputs: 46 kernel T2s (L=2) + 190 scalar features (L=0).
    Target: delta_shielding T2 (5 components).

    Normalization:
      - Kernels: per-protein std (angular structure preserved)
      - Scalars: z-score (train stats shared to val)
      - Targets: divide by target_std
    """

    def __init__(self, protein_ids: list[str], features_dir: Path,
                 norm_stats: NormStats | None = None):
        all_s, all_k, all_t, all_d = [], [], [], []
        protein_ids_loaded = []
        protein_boundaries = []
        offset = 0

        for pid in protein_ids:
            result, err = _build_features(features_dir / pid)
            if result is None:
                if err:
                    print(f"  skip {pid}: {err}")
                continue
            n = len(result["target"])
            all_s.append(result["scalars"])
            all_k.append(result["kernels"])
            all_t.append(result["target"])
            all_d.append(result["ring_dist"])
            protein_ids_loaded.append(pid)
            protein_boundaries.append((offset, offset + n))
            offset += n

        print(f"  loaded {len(protein_ids_loaded)} proteins, "
              f"skipped {len(protein_ids) - len(protein_ids_loaded)}")
        if not all_s:
            raise ValueError("No proteins loaded")

        self.protein_ids = protein_ids_loaded
        self.protein_boundaries = protein_boundaries

        scalars_np = np.vstack(all_s)
        kernels_np = np.vstack(all_k)
        targets_np = np.vstack(all_t)
        self.ring_dist = torch.tensor(np.concatenate(all_d), dtype=torch.float32)

        if norm_stats is None:
            target_std = float(targets_np.std())
            scalar_mean = scalars_np.mean(axis=0)
            scalar_std = np.maximum(scalars_np.std(axis=0), 1e-8)
            scalar_mean[_ONEHOT_COLS] = 0.0
            scalar_std[_ONEHOT_COLS] = 1.0

            # Per-kernel gate threshold: median magnitude of active atoms
            kernel_mags = np.linalg.norm(kernels_np, axis=-1)
            gate_threshold = np.ones(N_KERNELS)
            for k in range(N_KERNELS):
                nonzero = kernel_mags[:, k][kernel_mags[:, k] > 1e-8]
                if len(nonzero) > 100:
                    gate_threshold[k] = float(np.median(nonzero))

            self.norm_stats = NormStats(target_std, scalar_mean,
                                        scalar_std, gate_threshold)
        else:
            self.norm_stats = norm_stats

        ns = self.norm_stats
        scalars_np = (scalars_np - ns.scalar_mean) / ns.scalar_std
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

DIST_BANDS = [(0, 4), (4, 8), (8, 12), (12, 999)]
DIST_LABELS = ["0-4A", "4-8A", "8-12A", "12+A"]


def compute_r2(pred, tgt, ring_dist=None):
    """Per-component, overall, and distance-stratified R²."""
    result = {}
    for i, m in enumerate([-2, -1, 0, 1, 2]):
        ss_res = torch.sum((tgt[:, i] - pred[:, i]) ** 2).item()
        ss_tot = torch.sum((tgt[:, i] - tgt[:, i].mean()) ** 2).item()
        result[f"r2_m{m:+d}"] = 1.0 - ss_res / max(ss_tot, 1e-12)
    result["r2_mean"] = np.mean([result[f"r2_m{m:+d}"] for m in [-2, -1, 0, 1, 2]])

    ss_res = torch.sum((tgt - pred) ** 2).item()
    ss_tot = torch.sum((tgt - tgt.mean(dim=0, keepdim=True)) ** 2).item()
    result["r2_flat"] = 1.0 - ss_res / max(ss_tot, 1e-12)
    result["rmse_per_atom_ppm"] = torch.sqrt(
        torch.sum((tgt - pred) ** 2, dim=1)).mean().item()

    if ring_dist is not None:
        for (lo, hi), label in zip(DIST_BANDS, DIST_LABELS):
            mask = (ring_dist >= lo) & (ring_dist < hi)
            n = mask.sum().item()
            if n < 10:
                result[f"r2_{label}"] = float("nan")
                result[f"n_{label}"] = n
                continue
            p_b, t_b = pred[mask], tgt[mask]
            ss_res = torch.sum((t_b - p_b) ** 2).item()
            ss_tot = torch.sum((t_b - t_b.mean(dim=0, keepdim=True)) ** 2).item()
            result[f"r2_{label}"] = 1.0 - ss_res / max(ss_tot, 1e-12)
            result[f"n_{label}"] = n
    return result


def compute_naive_baselines(ds):
    """Physics baselines requiring no learning."""
    tgt = ds.targets
    kernels = ds.kernels
    ss_tot = torch.sum((tgt - tgt.mean(dim=0, keepdim=True)) ** 2).item()
    result = {}

    for name, slc in [("bs_only", slice(0, 8)),
                       ("ring_sum", slice(0, 24)),
                       ("all_sum", slice(None))]:
        pred = kernels[:, slc, :].sum(dim=1)
        ss_res = torch.sum((tgt - pred) ** 2).item()
        result[f"r2_{name}"] = 1.0 - ss_res / max(ss_tot, 1e-12)

    # Ridge baseline
    X = kernels.reshape(len(kernels), -1)
    try:
        lam = 1e-3
        XtX = X.T @ X + lam * torch.eye(X.shape[1], device=X.device)
        w = torch.linalg.solve(XtX, X.T @ tgt)
        ss_res = torch.sum((tgt - X @ w) ** 2).item()
        result["r2_ridge"] = 1.0 - ss_res / max(ss_tot, 1e-12)
    except Exception:
        result["r2_ridge"] = float("nan")

    if ds.ring_dist is not None:
        ring_sum = kernels[:, :24, :].sum(dim=1)
        for (lo, hi), label in zip(DIST_BANDS, DIST_LABELS):
            mask = (ds.ring_dist >= lo) & (ds.ring_dist < hi)
            n = mask.sum().item()
            if n < 10:
                result[f"ring_sum_{label}"] = float("nan")
                continue
            t_b, p_b = tgt[mask], ring_sum[mask]
            ss_res = torch.sum((t_b - p_b) ** 2).item()
            ss_tot_b = torch.sum((t_b - t_b.mean(dim=0, keepdim=True)) ** 2).item()
            result[f"ring_sum_{label}"] = 1.0 - ss_res / max(ss_tot_b, 1e-12)

    return result


def compute_per_protein_r2(pred, tgt, boundaries, protein_ids):
    """R² per protein."""
    result = {}
    for pid, (lo, hi) in zip(protein_ids, boundaries):
        if hi - lo < 10:
            result[pid] = float("nan")
            continue
        p, t = pred[lo:hi], tgt[lo:hi]
        ss_res = torch.sum((t - p) ** 2).item()
        ss_tot = torch.sum((t - t.mean(dim=0, keepdim=True)) ** 2).item()
        result[pid] = 1.0 - ss_res / max(ss_tot, 1e-12)
    return result
