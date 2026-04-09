#!/usr/bin/env python3
"""Tool 4: Per-ring-type R² ablation.

Separates "HIE geometric kernels are wrong" from "HIE DFT contrast is
unreliable" by computing ridge R² on ring-type-stratified atom subsets.

Two modes:
  1. Full-kernel R² per stratum: how well do ALL kernels explain atoms
     near each ring type?
  2. Self-fit R²: for ring type RT, use ONLY the 4 RT-specific kernels
     (BS_RT, HM_RT, Disp_RT, PQ_RT) to explain atoms near RT.
     This isolates kernel quality from DFT quality.

Outputs in {output_dir}/ablation/:
    ring_type_r2.csv   — per-stratum R² with kernel family slicing
    self_fit_r2.csv    — per-ring-type self-fit (4 kernels only)
    atom_counts.csv    — atom counts and target magnitudes per stratum

Usage:
    cd learn/src
    python -m secondary ablation --config calibration.toml
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .loader import (
    iter_proteins, setup_sdk, ridge_fit,
    RING_TYPE_NAMES, N_RING_TYPES,
)
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels


def run(cfg: Config, max_proteins: int = 0):
    out_dir = cfg.secondary.output_dir / "ablation"
    out_dir.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    strata = cfg.secondary.strata
    lam = cfg.secondary.ridge_lambda

    # Per-stratum accumulation
    strata_kernels: dict[str, list] = {s: [] for s in strata}
    strata_targets: dict[str, list] = {s: [] for s in strata}

    # Per-ring-type accumulation (for self-fit)
    # Atoms that see ring type i → their kernels + targets
    rt_kernels: dict[int, list] = {i: [] for i in range(N_RING_TYPES)}
    rt_targets: dict[int, list] = {i: [] for i in range(N_RING_TYPES)}

    print("Loading proteins...")
    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx
        kernels = assemble_kernels(p, idx, layout)  # (M, K, 5)
        target = p.delta.shielding.T2[idx]           # (M, 5)

        for j in range(len(idx)):
            mask = int(rec.atom_ring_masks[j])
            from .loader import assign_stratum
            s = assign_stratum(mask, strata)
            if s:
                strata_kernels[s].append(kernels[j])
                strata_targets[s].append(target[j])
            if s != "all":
                strata_kernels["all"].append(kernels[j])
                strata_targets["all"].append(target[j])

            # Per-ring-type: atom sees ring type i
            for i in range(N_RING_TYPES):
                if mask & (1 << i):
                    rt_kernels[i].append(kernels[j])
                    rt_targets[i].append(target[j])

    print(f"  Loaded {n_loaded} proteins")

    # ── 1. Per-stratum R² with kernel family slicing ────────────────
    r2_rows = []
    for s in strata:
        K_list = strata_kernels[s]
        T_list = strata_targets[s]
        if not K_list:
            continue
        K = np.array(K_list)
        T = np.array(T_list)
        N_s = len(T)

        row = {"stratum": s, "n_atoms": N_s}

        # Full kernel set
        row["r2_all"] = _safe_ridge(K, T, layout.n_kernels, lam)

        # Ring-type kernels only (first ring_type_end)
        row["r2_ring"] = _safe_ridge(
            K[:, :layout.ring_type_end, :], T, layout.ring_type_end, lam)

        # Bond-category kernels only
        n_bond = layout.bond_cat_end - layout.ring_type_end
        row["r2_bond"] = _safe_ridge(
            K[:, layout.ring_type_end:layout.bond_cat_end, :], T, n_bond, lam)

        # Calculator totals only
        n_total = layout.total_end - layout.bond_cat_end
        row["r2_total"] = _safe_ridge(
            K[:, layout.bond_cat_end:layout.total_end, :], T, n_total, lam)

        # EFG only
        n_efg = layout.efg_end - layout.total_end
        row["r2_efg"] = _safe_ridge(
            K[:, layout.total_end:layout.efg_end, :], T, n_efg, lam)

        # Per-ring individual kernels only
        n_perring = layout.n_kernels - layout.efg_end
        row["r2_perring"] = _safe_ridge(
            K[:, layout.efg_end:, :], T, n_perring, lam)

        r2_rows.append(row)

    if r2_rows:
        _write_csv(out_dir / "ring_type_r2.csv", r2_rows)
        print("\n  Per-stratum R²:")
        for r in r2_rows:
            print(f"    {r['stratum']:15s}  n={r['n_atoms']:6d}  "
                  f"R²_all={r['r2_all']:.3f}  R²_ring={r['r2_ring']:.3f}  "
                  f"R²_efg={r['r2_efg']:.3f}")

    # ── 2. Per-ring-type self-fit ───────────────────────────────────
    # For ring type RT: use only BS_RT, HM_RT, Disp_RT, PQ_RT
    # These are at known positions in the ring_type block:
    #   BS:   columns 0..7   → RT is column RT
    #   HM:   columns 8..15  → RT is column 8 + RT
    #   Disp: columns 16..23 → RT is column 16 + RT
    #   PQ:   columns 24..31 → RT is column 24 + RT
    self_fit_rows = []
    for rt_idx in range(N_RING_TYPES):
        K_list = rt_kernels[rt_idx]
        T_list = rt_targets[rt_idx]
        if not K_list:
            self_fit_rows.append({
                "ring_type": RING_TYPE_NAMES[rt_idx],
                "n_atoms": 0,
                "r2_self_4kernel": float("nan"),
                "r2_all_kernels": float("nan"),
                "target_mag_mean": float("nan"),
            })
            continue

        K = np.array(K_list)
        T = np.array(T_list)
        N_rt = len(T)

        # Extract the 4 ring-type-specific kernels
        rt_cols = [rt_idx, 8 + rt_idx, 16 + rt_idx, 24 + rt_idx]
        # Only include if within ring_type_end
        rt_cols = [c for c in rt_cols if c < layout.ring_type_end]
        K_self = K[:, rt_cols, :]

        r2_self = _safe_ridge(K_self, T, len(rt_cols), lam)
        r2_all = _safe_ridge(K, T, layout.n_kernels, lam)
        target_mag = np.linalg.norm(T, axis=1).mean()

        self_fit_rows.append({
            "ring_type": RING_TYPE_NAMES[rt_idx],
            "n_atoms": N_rt,
            "r2_self_4kernel": float(r2_self),
            "r2_all_kernels": float(r2_all),
            "target_mag_mean": float(target_mag),
        })

    if self_fit_rows:
        _write_csv(out_dir / "self_fit_r2.csv", self_fit_rows)
        print("\n  Per-ring-type self-fit:")
        for r in self_fit_rows:
            print(f"    {r['ring_type']:15s}  n={r['n_atoms']:6d}  "
                  f"R²_self={r['r2_self_4kernel']:.3f}  "
                  f"R²_all={r['r2_all_kernels']:.3f}  "
                  f"T2_mag={r['target_mag_mean']:.3f}")

    # ── 3. Atom counts ──────────────────────────────────────────────
    count_rows = []
    for s in strata:
        T_list = strata_targets[s]
        if not T_list:
            continue
        T = np.array(T_list)
        count_rows.append({
            "stratum": s,
            "n_atoms": len(T),
            "target_mag_mean": float(np.linalg.norm(T, axis=1).mean()),
            "target_mag_std": float(np.linalg.norm(T, axis=1).std()),
        })
    if count_rows:
        _write_csv(out_dir / "atom_counts.csv", count_rows)


def _safe_ridge(K: np.ndarray, T: np.ndarray, n_kernels: int,
                lam: float) -> float:
    """Ridge fit, returning NaN if underdetermined."""
    N = len(T)
    n_features = n_kernels * 5
    if N < max(20, n_features // 5):
        return float("nan")
    X = K.reshape(N, -1)
    try:
        _, r2 = ridge_fit(X, T, lam)
        return float(r2)
    except np.linalg.LinAlgError:
        return float("nan")


def _write_csv(path: Path, rows: list[dict]):
    fields = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  Wrote {path} ({len(rows)} rows)")
