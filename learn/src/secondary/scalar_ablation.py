#!/usr/bin/env python3
"""Tool 5: Scalar group ablation — what drives the linear→equivariant gap.

Two analyses:
  1. Interaction ridge: for each scalar group, add scalar × kernel
     interaction terms to the base ridge and measure R² improvement.
     Shows which scalars help even in a linear model.
  2. Drop-one-group: for each scalar group, remove it from the full
     scalar set and measure R² of the base ridge on the remaining
     kernel-only features (doesn't require retraining the MLP —
     measures how much each scalar group contributes via interaction).

The base ridge uses the top-K kernels from forward selection to keep
feature counts manageable.

Outputs in {output_dir}/scalar_ablation/:
    interaction_r2.csv     — R² with each scalar group as interactions
    group_summary.csv      — scalar groups with widths and ridge impact
    eigenspectrum.csv      — interaction feature space dimensionality

Usage:
    cd learn/src
    python -m secondary scalar_ablation --config calibration.toml
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .loader import (
    iter_proteins, setup_sdk, ridge_fit,
    ELEMENT_NAMES, RING_TYPE_NAMES,
)
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels, normalize_kernels
from mutation_set.scalars import assemble_scalars, ScalarLayout


def _forward_select_top_k(kernels: np.ndarray, target: np.ndarray,
                           layout: KernelLayout, lam: float,
                           k: int) -> list[int]:
    """Greedy forward selection of top-k kernels by ridge R²."""
    N = len(target)
    remaining = set(range(layout.n_kernels))
    selected = []
    for _ in range(k):
        best_idx, best_r2 = -1, -np.inf
        for ki in remaining:
            trial = selected + [ki]
            X = kernels[:, trial, :].reshape(N, -1)
            _, r2 = ridge_fit(X, target, lam)
            if r2 > best_r2:
                best_r2, best_idx = r2, ki
        if best_idx < 0:
            break
        selected.append(best_idx)
        remaining.remove(best_idx)
    return selected


def run(cfg: Config, max_proteins: int = 0):
    out_dir = cfg.secondary.output_dir / "scalar_ablation"
    out_dir.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    lam = cfg.secondary.ridge_lambda
    TOP_K = 10  # top kernels for interaction analysis

    # ── Load all data ───────────────────────────────────────────────
    print("Loading proteins...")
    all_kernels, all_targets, all_scalars = [], [], []
    scalar_layout = None
    n_loaded = 0

    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        p = rec.protein
        idx = rec.matched_idx
        kernels = assemble_kernels(p, idx, layout)
        kernels, kernel_scales = normalize_kernels(
            kernels, cfg.normalization.kernel_std_floor)
        scalars, sl = assemble_scalars(p, idx, kernel_scales, cfg)
        target = p.delta.shielding.T2[idx]

        all_kernels.append(kernels)
        all_targets.append(target)
        all_scalars.append(scalars)
        if scalar_layout is None:
            scalar_layout = sl
        n_loaded += 1

    print(f"  Loaded {n_loaded} proteins")
    if not all_kernels:
        print("No data.")
        return

    K = np.vstack(all_kernels)   # (N, n_kernels, 5)
    T = np.vstack(all_targets)   # (N, 5)
    S = np.vstack(all_scalars)   # (N, n_scalars)
    N = len(T)
    print(f"  {N} atoms, {layout.n_kernels} kernels, "
          f"{scalar_layout.total} scalars")

    # ── Forward select top-K kernels ────────────────────────────────
    print(f"\n  Forward selecting top {TOP_K} kernels...")
    top_k = _forward_select_top_k(K, T, layout, lam, TOP_K)
    top_names = [layout.names[k] for k in top_k]
    print(f"  Selected: {', '.join(top_names)}")

    K_top = K[:, top_k, :]  # (N, TOP_K, 5)
    X_base = K_top.reshape(N, -1)  # (N, TOP_K * 5)
    _, r2_base = ridge_fit(X_base, T, lam)
    print(f"  Base ridge R² (top-{TOP_K} kernels): {r2_base:.4f}")

    # ── Full kernel ridge ───────────────────────────────────────────
    X_full = K.reshape(N, -1)
    _, r2_full = ridge_fit(X_full, T, lam)
    print(f"  Full ridge R² (all {layout.n_kernels} kernels): {r2_full:.4f}")

    # ── Interaction analysis ────────────────────────────────────────
    # For each scalar group: create interaction features
    # X_interact = [X_base, s_1 * X_base, s_2 * X_base, ...]
    # where s_i are the scalars in that group.
    # Then ridge fit → R² with interaction terms.

    print(f"\n{'=' * 60}")
    print("SCALAR GROUP INTERACTION ANALYSIS")
    print(f"{'=' * 60}")
    print(f"  Base: top-{TOP_K} kernels → {X_base.shape[1]} features → "
          f"R²={r2_base:.4f}")
    print(f"  For each scalar group: add scalar × kernel interactions\n")

    rows = []
    offset = 0
    for gi, (gname, gwidth) in enumerate(
            zip(scalar_layout.names, scalar_layout.widths)):
        s_group = S[:, offset:offset + gwidth]  # (N, gwidth)
        offset += gwidth

        # Build interaction features: each scalar × each kernel feature
        # Shape: (N, gwidth * TOP_K * 5)
        interactions = []
        for si in range(gwidth):
            interactions.append(s_group[:, si:si+1] * X_base)
        X_interact = np.concatenate([X_base] + interactions, axis=1)

        n_features = X_interact.shape[1]
        if N < max(n_features * 2, 100):
            r2_with = float("nan")
        else:
            try:
                _, r2_with = ridge_fit(X_interact, T, lam)
            except np.linalg.LinAlgError:
                r2_with = float("nan")

        delta = r2_with - r2_base if not np.isnan(r2_with) else float("nan")

        rows.append({
            "group_name": gname,
            "group_width": gwidth,
            "n_interaction_features": n_features - X_base.shape[1],
            "r2_base": float(r2_base),
            "r2_with_interactions": float(r2_with),
            "r2_delta": float(delta),
        })

        bar = "#" * min(int(max(delta, 0) * 200), 40) if not np.isnan(delta) else ""
        print(f"  {gname:25s}  width={gwidth:3d}  "
              f"R²={r2_with:.4f}  Δ={delta:+.4f}  {bar}")

    # Sort by impact
    rows.sort(key=lambda r: r["r2_delta"] if not np.isnan(r["r2_delta"]) else -999,
              reverse=True)

    _write_csv(out_dir / "interaction_r2.csv", rows)

    # ── Summary with cumulative ─────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("CUMULATIVE INTERACTION (groups added in order of impact)")
    print(f"{'=' * 60}\n")

    cum_features = [X_base]
    cum_r2 = r2_base
    cum_rows = [{"step": 0, "group": "base_kernels",
                 "r2": float(r2_base), "delta": 0.0}]

    offset_map = {}
    o = 0
    for gn, gw in zip(scalar_layout.names, scalar_layout.widths):
        offset_map[gn] = (o, o + gw)
        o += gw

    for r in rows:
        gname = r["group_name"]
        if np.isnan(r["r2_delta"]) or r["r2_delta"] < 0.001:
            continue
        lo, hi = offset_map[gname]
        s_group = S[:, lo:hi]
        interactions = []
        for si in range(hi - lo):
            interactions.append(s_group[:, si:si+1] * X_base)
        cum_features.extend(interactions)
        X_cum = np.concatenate(cum_features, axis=1)
        if N < X_cum.shape[1] * 2:
            break
        try:
            _, r2_cum = ridge_fit(X_cum, T, lam)
        except np.linalg.LinAlgError:
            break

        step = len(cum_rows)
        delta_cum = r2_cum - cum_r2
        cum_rows.append({"step": step, "group": gname,
                         "r2": float(r2_cum), "delta": float(delta_cum)})
        cum_r2 = r2_cum
        print(f"  +{gname:25s}  R²={r2_cum:.4f}  (+{delta_cum:.4f})")

    _write_csv(out_dir / "group_summary.csv", cum_rows)

    # ── Eigenspectrum of interaction space ───────────────────────────
    # How many effective dimensions does the full interaction space have?
    print(f"\n  Computing eigenspectrum of full interaction space...")
    all_interactions = [X_base]
    for gi, (gname, gwidth) in enumerate(
            zip(scalar_layout.names, scalar_layout.widths)):
        lo, hi = offset_map[gname]
        s_group = S[:, lo:hi]
        for si in range(gwidth):
            all_interactions.append(s_group[:, si:si+1] * X_base)

    X_all_interact = np.concatenate(all_interactions, axis=1)
    print(f"  Full interaction features: {X_all_interact.shape[1]}")

    # Covariance eigenspectrum (feature × feature)
    if X_all_interact.shape[1] < 5000:
        cov = X_all_interact.T @ X_all_interact / N
        eigvals = np.linalg.eigvalsh(cov)[::-1]
        eigvals = np.maximum(eigvals, 0)
        total_var = eigvals.sum()
        cumvar = np.cumsum(eigvals) / total_var if total_var > 0 else np.zeros_like(eigvals)

        eigen_rows = []
        for rank, (ev, cv) in enumerate(zip(eigvals[:50], cumvar[:50])):
            eigen_rows.append({"rank": rank + 1, "eigenvalue": float(ev),
                               "cumulative_variance": float(cv)})
        _write_csv(out_dir / "eigenspectrum.csv", eigen_rows)

        dim90 = int(np.searchsorted(cumvar, 0.90)) + 1
        dim95 = int(np.searchsorted(cumvar, 0.95)) + 1
        print(f"  Effective dim: {dim90} (90%), {dim95} (95%)")
    else:
        print(f"  Skipping eigenspectrum (too many features: "
              f"{X_all_interact.shape[1]})")

    # ── Ridge on full interactions ──────────────────────────────────
    if N > X_all_interact.shape[1] * 2:
        try:
            _, r2_all_interact = ridge_fit(X_all_interact, T, lam)
            print(f"\n  Full interaction ridge R²: {r2_all_interact:.4f}")
            print(f"  (Base: {r2_base:.4f}, Full kernels: {r2_full:.4f}, "
                  f"Equivariant target: ~0.61)")
        except np.linalg.LinAlgError:
            print("  Full interaction ridge failed")


def _write_csv(path: Path, rows: list[dict]):
    fields = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  Wrote {path} ({len(rows)} rows)")
