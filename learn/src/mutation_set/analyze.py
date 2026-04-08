#!/usr/bin/env python3
"""Cold analytical diagnostic — no learning, pure numpy.

Sections:
    1. Per-kernel correlation with DFT delta T2
    2. Per-element breakdown (target magnitude + best kernel)
    3. Signal by distance to nearest removed ring
    4. Per-protein ridge R² (optimal linear mixing)
    5. Greedy forward kernel selection
    6. Leave-one-out kernel ablation
    7. Residual subspace projection (OLS vs ridge + by element/distance)

All thresholds from config.  Kernel names from KernelLayout.

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

ELEMENT_NAMES = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S"}


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


def _ridge_fit(X, y, lam):
    """Ridge regression, returns (predictions, R²)."""
    XtX = X.T @ X + lam * np.eye(X.shape[1])
    w = np.linalg.solve(XtX, X.T @ y)
    pred = X @ w
    return pred, _r2(pred, y)


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
                "element": result["element"],
                "pid": pid,
            })
    print(f"Loaded {len(all_data)} proteins with delta\n")
    if not all_data:
        return

    kernels = np.vstack([d["kernels"] for d in all_data])
    target = np.vstack([d["target"] for d in all_data])
    ring_dist = np.concatenate([d["ring_dist"] for d in all_data])
    element = np.concatenate([d["element"] for d in all_data])
    N = len(target)
    n_kernels = layout.n_kernels
    lam = cfg.analysis.ridge_lambda
    min_band = cfg.analysis.min_atoms_per_band

    print(f"Atoms: {N}, Kernels: {n_kernels}")
    print(f"Target T2 std: {target.std():.4f} ppm")
    print(f"Target T2 mean magnitude: {np.linalg.norm(target, axis=1).mean():.4f} ppm\n")

    # ── 1. Per-kernel correlation ────────────────────────────────
    print("=" * 60)
    print("1. PER-KERNEL CORRELATION WITH DFT DELTA T2")
    print("-" * 60)
    corrs = [(layout.names[k], _corr(kernels[:, k, :].ravel(), target.ravel()))
             for k in range(n_kernels)]
    corrs.sort(key=lambda x: abs(x[1]), reverse=True)
    for name, r in corrs:
        bar = "#" * int(abs(r) * 50)
        sign = "+" if r > 0 else "-"
        active = kernels[:, layout.names.index(name), :].std() > 1e-10
        flag = "" if active else " [ZERO]"
        print(f"  {name:30s}  r={r:+.4f}  {sign}{bar}{flag}")
    top5 = corrs[:5]
    print(f"\n  Top 5: {', '.join(f'{n}({r:+.3f})' for n, r in top5)}")

    # ── 2. Per-element breakdown ─────────────────────────────────
    print(f"\n{'=' * 60}")
    print("2. PER-ELEMENT TARGET MAGNITUDE AND BEST KERNEL")
    print("-" * 60)
    for elem_id, elem_name in sorted(ELEMENT_NAMES.items()):
        mask = element == elem_id
        n = mask.sum()
        if n < 5:
            continue
        t = target[mask]
        mag = np.linalg.norm(t, axis=1)
        best_r, best_k = 0, ""
        for k in range(n_kernels):
            kdat = kernels[mask, k, :]
            if kdat.std() < 1e-10:
                continue
            r = _corr(kdat.ravel(), t.ravel())
            if abs(r) > abs(best_r):
                best_r, best_k = r, layout.names[k]
        print(f"  {elem_name:2s}  n={n:6d}  T2_mag={mag.mean():.3f}+/-{mag.std():.3f} ppm  "
              f"best={best_k}(r={best_r:+.3f})")

    # ── 3. Signal by distance ────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("3. SIGNAL BY DISTANCE TO NEAREST REMOVED RING")
    print("-" * 60)
    for lo, hi in cfg.analysis.distance_bands:
        label = f"{int(lo)}-{int(hi)}A" if hi < 900 else f"{int(lo)}+A"
        mask = (ring_dist >= lo) & (ring_dist < hi)
        n = mask.sum()
        if n < min_band:
            print(f"  {label:6s}  n={n:6d}  (too few)")
            continue
        X = kernels[mask].reshape(n, -1)
        y = target[mask]
        try:
            _, r2 = _ridge_fit(X, y, lam)
        except np.linalg.LinAlgError:
            r2 = float("nan")
        mag = np.linalg.norm(y, axis=1)
        print(f"  {label:6s}  n={n:6d}  T2_mag={mag.mean():.3f}+/-{mag.std():.3f} ppm  "
              f"ridge_R²={r2:.3f}")

    # ── 4. Per-protein ridge ─────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("4. PER-PROTEIN RIDGE R²")
    print("-" * 60)
    protein_r2s = []
    for d in all_data:
        k, t = d["kernels"], d["target"]
        M = len(t)
        if M < 20:
            continue
        X = k.reshape(M, -1)
        try:
            _, r2 = _ridge_fit(X, t, lam)
        except np.linalg.LinAlgError:
            r2 = float("nan")
        protein_r2s.append((d["pid"], r2, M))

    protein_r2s.sort(key=lambda x: x[1])
    vals = [r for _, r, _ in protein_r2s if not np.isnan(r)]
    print(f"  Proteins: {len(protein_r2s)}")
    print(f"  Ridge R²: median={np.median(vals):.3f}  "
          f"mean={np.mean(vals):.3f}  std={np.std(vals):.3f}")
    print(f"  Range: [{min(vals):.3f}, {max(vals):.3f}]")
    print(f"\n  Worst 5:")
    for pid, r2, m in protein_r2s[:5]:
        print(f"    {pid:20s}  R²={r2:.3f}  atoms={m}")
    print(f"  Best 5:")
    for pid, r2, m in protein_r2s[-5:]:
        print(f"    {pid:20s}  R²={r2:.3f}  atoms={m}")

    # ── 5. Greedy forward kernel selection ────────────────────────
    print(f"\n{'=' * 60}")
    print("5. GREEDY FORWARD KERNEL SELECTION")
    print("-" * 60)
    remaining = set(range(n_kernels))
    selected = []
    cumulative_r2 = []
    for step in range(min(n_kernels, cfg.analysis.max_forward_steps)):
        best_k, best_r2 = -1, -np.inf
        for k in remaining:
            trial = selected + [k]
            X = kernels[:, trial, :].reshape(N, -1)
            try:
                _, r2 = _ridge_fit(X, target, lam)
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

    if len(cumulative_r2) >= 5:
        print(f"\n  First 5 reach R²={cumulative_r2[4]:.4f}, "
              f"all {len(selected)} reach R²={cumulative_r2[-1]:.4f}")

    # ── 6. Leave-one-out kernel ablation ─────────────────────────
    print(f"\n{'=' * 60}")
    print("6. LEAVE-ONE-OUT KERNEL ABLATION")
    print("-" * 60)
    X_all = kernels.reshape(N, -1)
    try:
        _, r2_full = _ridge_fit(X_all, target, lam)
    except np.linalg.LinAlgError:
        r2_full = float("nan")
    print(f"  Full ridge R² ({n_kernels} kernels): {r2_full:.4f}\n")

    ablations = []
    for k in range(n_kernels):
        keep = [j for j in range(n_kernels) if j != k]
        X = kernels[:, keep, :].reshape(N, -1)
        try:
            _, r2 = _ridge_fit(X, target, lam)
        except np.linalg.LinAlgError:
            r2 = float("nan")
        ablations.append((layout.names[k], r2_full - r2, r2))

    ablations.sort(key=lambda x: x[1], reverse=True)
    for name, drop, r2_without in ablations:
        bar = "#" * min(int(drop * 500), 50) if drop > 1e-5 else ""
        print(f"  {name:30s}  drop={drop:+.4f}  R²_without={r2_without:.4f}  {bar}")

    unique = [n for n, d, _ in ablations if d > 0.005]
    redundant = [n for n, d, _ in ablations if d < 0.0005]
    print(f"\n  Unique (drop > 0.005):    {', '.join(unique) if unique else 'none'}")
    print(f"  Redundant (drop < 0.0005): {len(redundant)} kernels")

    # ── 7. Residual subspace projection ──────────────────────────
    print(f"\n{'=' * 60}")
    print("7. RESIDUAL SUBSPACE PROJECTION")
    print("-" * 60)
    try:
        ridge_pred, _ = _ridge_fit(X_all, target, lam)
        residual = target - ridge_pred

        # OLS vs ridge
        ols_pred = np.zeros_like(target)
        for comp in range(5):
            w_ols, _, _, _ = np.linalg.lstsq(X_all, target[:, comp], rcond=1e-8)
            ols_pred[:, comp] = X_all @ w_ols
        r2_ols = _r2(ols_pred, target)

        in_span = (r2_ols - r2_full) / (1 - r2_full) if r2_full < 1.0 else 0.0
        print(f"  Ridge R² (λ={lam}):    {r2_full:.4f}")
        print(f"  OLS R² (no reg):       {r2_ols:.4f}")
        print(f"  In-span (OLS-ridge):   {r2_ols - r2_full:.4f} "
              f"({in_span:.1%} of unexplained)")
        print(f"  Out-of-span (1-OLS):   {1 - r2_ols:.4f}")

        if in_span > 0.15:
            print(f"\n  >> Ridge leaving {in_span:.0%} of residual on the table."
                  f"\n     MLP environment-dependent weights should recover this.")
        else:
            print(f"\n  >> Regularisation costs little. The missing {1-r2_ols:.0%}"
                  f" is truly out-of-span.")

        # Residual by element
        print(f"\n  Residual by element:")
        for elem_id, elem_name in sorted(ELEMENT_NAMES.items()):
            mask = element == elem_id
            n = mask.sum()
            if n < 5:
                continue
            res = residual[mask]
            print(f"    {elem_name:2s}  n={n:6d}  "
                  f"residual_mag={np.linalg.norm(res, axis=1).mean():.3f} ppm  "
                  f"bias={res.mean():.4f}")

        # Residual by distance
        print(f"\n  Residual by distance:")
        for lo, hi in cfg.analysis.distance_bands:
            label = f"{int(lo)}-{int(hi)}A" if hi < 900 else f"{int(lo)}+A"
            mask = (ring_dist >= lo) & (ring_dist < hi)
            n = mask.sum()
            if n < min_band:
                continue
            res = residual[mask]
            print(f"    {label:6s}  n={n:6d}  "
                  f"residual_mag={np.linalg.norm(res, axis=1).mean():.3f} ppm  "
                  f"bias={res.mean():.4f}")

        # Per-kernel alignment with residual
        print(f"\n  Per-kernel alignment with residual (mean |cos|):")
        alignments = []
        for k in range(n_kernels):
            kv = kernels[:, k, :]
            k_norm = np.linalg.norm(kv, axis=1, keepdims=True)
            r_norm = np.linalg.norm(residual, axis=1, keepdims=True)
            valid = ((k_norm > 1e-10) & (r_norm > 1e-10)).ravel()
            if valid.sum() < 10:
                alignments.append((layout.names[k], 0.0, 0))
                continue
            cos = np.sum(kv[valid] * residual[valid], axis=1) / (
                k_norm[valid].ravel() * r_norm[valid].ravel())
            alignments.append((layout.names[k], np.mean(np.abs(cos)), int(valid.sum())))

        alignments.sort(key=lambda x: x[1], reverse=True)
        for name, mean_cos, nv in alignments[:15]:
            bar = "#" * int(mean_cos * 50)
            print(f"    {name:30s}  |cos|={mean_cos:.3f}  n={nv:6d}  {bar}")

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
