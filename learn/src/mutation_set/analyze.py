#!/usr/bin/env python3
"""Cold analytical diagnostic — no learning, pure numpy.

Examines the raw relationship between per-protein-normalized kernel T2
and DFT delta T2 across all proteins.  All thresholds from config.

Usage:
    cd learn/src
    python -m mutation_set.analyze --config calibration.toml
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from .config import load_config
from .dataset import list_proteins, _build_one
from .kernels import KernelLayout


def _corr(a, b):
    if len(a) < 3:
        return float("nan")
    a, b = a - a.mean(), b - b.mean()
    denom = np.sqrt((a ** 2).sum() * (b ** 2).sum())
    return float((a * b).sum() / denom) if denom > 1e-12 else 0.0


def _r2(pred, tgt):
    ss_res = np.sum((tgt - pred) ** 2)
    ss_tot = np.sum((tgt - tgt.mean(axis=0, keepdims=True)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0


def analyze(cfg):
    layout = KernelLayout.from_config(cfg)
    features_dir = cfg.paths.features
    proteins = list_proteins(features_dir)
    if cfg.data.max_proteins > 0:
        proteins = proteins[:cfg.data.max_proteins]

    print(f"Loading {len(proteins)} proteins...")
    all_data = []
    for pid in proteins:
        result, err = _build_one(features_dir / pid, cfg, layout)
        if result is not None:
            all_data.append({
                "kernels": result["kernels"],
                "target": result["target"],
                "ring_dist": result["ring_dist"],
                "pid": pid,
            })
    print(f"Loaded {len(all_data)} proteins with delta\n")
    if not all_data:
        return

    kernels = np.vstack([d["kernels"] for d in all_data])
    target = np.vstack([d["target"] for d in all_data])
    ring_dist = np.concatenate([d["ring_dist"] for d in all_data])
    N = len(target)
    lam = cfg.analysis.ridge_lambda

    print(f"Atoms: {N}, Kernels: {layout.n_kernels}")
    print(f"Target T2 std: {target.std():.4f} ppm\n")

    # ── 1. Per-kernel correlation ────────────────────────────────
    print("=" * 60)
    print("1. PER-KERNEL CORRELATION WITH DFT DELTA T2")
    print("-" * 60)
    corrs = [(layout.names[k], _corr(kernels[:, k, :].ravel(), target.ravel()))
             for k in range(layout.n_kernels)]
    corrs.sort(key=lambda x: abs(x[1]), reverse=True)
    for name, r in corrs:
        bar = "#" * int(abs(r) * 50)
        sign = "+" if r > 0 else "-"
        print(f"  {name:30s}  r={r:+.4f}  {sign}{bar}")

    # ── 2. Per-distance ridge ────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("2. SIGNAL BY DISTANCE")
    print("-" * 60)
    for lo, hi in cfg.analysis.distance_bands:
        label = f"{int(lo)}-{int(hi)}A" if hi < 900 else f"{int(lo)}+A"
        mask = (ring_dist >= lo) & (ring_dist < hi)
        n = mask.sum()
        if n < cfg.analysis.min_atoms_per_band:
            print(f"  {label:6s}  n={n:6d}  (too few)")
            continue
        X = kernels[mask].reshape(n, -1)
        y = target[mask]
        try:
            XtX = X.T @ X + lam * np.eye(X.shape[1])
            w = np.linalg.solve(XtX, X.T @ y)
            r2 = _r2(X @ w, y)
        except np.linalg.LinAlgError:
            r2 = float("nan")
        mag = np.linalg.norm(y, axis=1)
        print(f"  {label:6s}  n={n:6d}  T2_mag={mag.mean():.3f} ppm  "
              f"ridge_R²={r2:.3f}")

    # ── 3. Per-protein ridge ─────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("3. PER-PROTEIN RIDGE R²")
    print("-" * 60)
    protein_r2s = []
    for d in all_data:
        k, t = d["kernels"], d["target"]
        M = len(t)
        if M < 20:
            continue
        X = k.reshape(M, -1)
        try:
            XtX = X.T @ X + lam * np.eye(X.shape[1])
            w = np.linalg.solve(XtX, X.T @ t)
            r2 = _r2(X @ w, t)
        except np.linalg.LinAlgError:
            r2 = float("nan")
        protein_r2s.append((d["pid"], r2, M))

    protein_r2s.sort(key=lambda x: x[1])
    vals = [r for _, r, _ in protein_r2s if not np.isnan(r)]
    print(f"  Proteins: {len(protein_r2s)}")
    print(f"  Ridge R²: median={np.median(vals):.3f}  "
          f"mean={np.mean(vals):.3f}")
    print(f"  Worst 3:")
    for pid, r2, m in protein_r2s[:3]:
        print(f"    {pid:20s}  R²={r2:.3f}  atoms={m}")
    print(f"  Best 3:")
    for pid, r2, m in protein_r2s[-3:]:
        print(f"    {pid:20s}  R²={r2:.3f}  atoms={m}")

    # ── 4. Greedy forward selection ──────────────────────────────
    print(f"\n{'=' * 60}")
    print("4. GREEDY FORWARD KERNEL SELECTION")
    print("-" * 60)
    remaining = set(range(layout.n_kernels))
    selected = []
    cumulative_r2 = []
    for step in range(min(layout.n_kernels, cfg.analysis.max_forward_steps)):
        best_k, best_r2 = -1, -np.inf
        for k in remaining:
            trial = selected + [k]
            X = kernels[:, trial, :].reshape(N, -1)
            try:
                XtX = X.T @ X + lam * np.eye(X.shape[1])
                w = np.linalg.solve(XtX, X.T @ target)
                r2 = _r2(X @ w, target)
            except np.linalg.LinAlgError:
                r2 = -np.inf
            if r2 > best_r2:
                best_r2, best_k = r2, k
        if best_k < 0:
            break
        selected.append(best_k)
        remaining.remove(best_k)
        cumulative_r2.append(best_r2)
        delta = best_r2 - (cumulative_r2[-2] if len(cumulative_r2) > 1 else 0.0)
        print(f"  {step+1:2d}. +{layout.names[best_k]:30s}  "
              f"R²={best_r2:.4f}  (+{delta:.4f})")

    # ── 5. Residual subspace ─────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("5. RESIDUAL SUBSPACE PROJECTION")
    print("-" * 60)
    X_all = kernels.reshape(N, -1)
    try:
        XtX = X_all.T @ X_all + lam * np.eye(X_all.shape[1])
        w_full = np.linalg.solve(XtX, X_all.T @ target)
        r2_ridge = _r2(X_all @ w_full, target)

        # OLS (no regularisation) vs ridge to measure in-span residual
        ols_pred = np.zeros_like(target)
        for comp in range(5):
            w_ols, _, _, _ = np.linalg.lstsq(X_all, target[:, comp], rcond=1e-8)
            ols_pred[:, comp] = X_all @ w_ols
        r2_ols = _r2(ols_pred, target)

        in_span = (r2_ols - r2_ridge) / (1 - r2_ridge) if r2_ridge < 1.0 else 0.0
        print(f"  Ridge R² (λ={lam}):    {r2_ridge:.4f}")
        print(f"  OLS R² (no reg):       {r2_ols:.4f}")
        print(f"  In-span (OLS-ridge):   {r2_ols - r2_ridge:.4f} "
              f"({in_span:.1%} of unexplained)")
        print(f"  Out-of-span (1-OLS):   {1 - r2_ols:.4f}")
    except np.linalg.LinAlgError:
        print("  Ridge fit failed")

    print(f"\n{'=' * 60}")
    print("Done.")


def main():
    parser = argparse.ArgumentParser(description="Cold analytical diagnostic")
    parser.add_argument("--config", type=str, default="calibration.toml")
    parser.add_argument("--max-proteins", type=int, default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    if args.max_proteins is not None:
        from dataclasses import replace
        cfg = replace(cfg, data=replace(cfg.data, max_proteins=args.max_proteins))
    sys.path.insert(0, str(cfg.paths.sdk))
    analyze(cfg)


if __name__ == "__main__":
    main()
