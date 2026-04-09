#!/usr/bin/env python3
"""Clean per-element ridge calibration on 55 core kernels.

No MLP, no gating, no scalars. The ridge coefficient vector IS the
table of calibrated physical constants. Output: one table per element,
sorted by |coefficient|, with cumulative R² showing how many kernels
you actually need.

Usage:
    cd learn/src
    python3 -c "
    from mutation_set.config import load_config
    import sys; cfg = load_config('calibration.toml')
    sys.path.insert(0, str(cfg.paths.sdk))
    from actual_physics.clean_calibration import run
    run(cfg)
    "
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from secondary.loader import iter_proteins, setup_sdk, ELEMENT_NAMES
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout


# The 55 core kernels: ring-type + bond-category + totals + EFG.
# No per-ring residuals, no RBF, no angular basis.
def _core_kernel_mask(layout: KernelLayout) -> list[int]:
    """Indices of the 55 core kernels."""
    return list(range(layout.efg_end))


def _ridge_coeffs(X, y, lam=1.0):
    """Ridge regression returning coefficients, predictions, R²."""
    XtX = X.T @ X + lam * np.eye(X.shape[1])
    w = np.linalg.solve(XtX, X.T @ y)
    pred = X @ w
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean(axis=0)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
    return w, pred, r2


def _forward_r2(X, y, lam, n_kernels, kernel_names):
    """Greedy forward selection returning ordered list of
    (kernel_name, marginal_delta, cumulative_R²)."""
    selected = []
    available = set(range(n_kernels))
    results = []
    prev_r2 = 0.0

    for _ in range(min(n_kernels, 20)):
        best_k, best_r2 = -1, -1e30
        for k in available:
            trial = selected + [k]
            Xt = X[:, trial, :].reshape(X.shape[0], -1)
            _, _, r2 = _ridge_coeffs(Xt, y, lam)
            if r2 > best_r2:
                best_k, best_r2 = k, r2
        if best_k < 0:
            break
        selected.append(best_k)
        available.discard(best_k)
        delta = best_r2 - prev_r2
        results.append((kernel_names[best_k], delta, best_r2))
        prev_r2 = best_r2
    return results


def run(cfg: Config, max_proteins: int = 0):
    out = Path("output/actual_physics/clean_calibration")
    out.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    core = _core_kernel_mask(layout)
    core_names = [layout.names[i] for i in core]
    n_core = len(core)
    lam = getattr(cfg.secondary, 'ridge_lambda', 1.0)

    print(f"Core kernels: {n_core} (of {layout.n_kernels} total)")
    print(f"  ring_type: 0-{layout.ring_type_end}")
    print(f"  bond_cat:  {layout.ring_type_end}-{layout.bond_cat_end}")
    print(f"  total:     {layout.bond_cat_end}-{layout.total_end}")
    print(f"  efg:       {layout.total_end}-{layout.efg_end}")

    # Accumulate
    ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O"}
    elem_K = {e: [] for e in ELEMENTS}
    elem_Y = {e: [] for e in ELEMENTS}
    all_K, all_Y = [], []

    print("Loading proteins...")
    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx

        # Assemble kernels manually for just core (reuse the full assembly
        # but slice to core indices)
        from mutation_set.kernels import assemble_kernels
        kernels_full = assemble_kernels(p, idx, layout)
        kernels = kernels_full[:, core, :]  # (M, 55, 5)
        target = p.delta.shielding.T2[idx]  # (M, 5)

        for j in range(len(idx)):
            e = int(rec.element[j])
            if e in elem_K:
                elem_K[e].append(kernels[j])
                elem_Y[e].append(target[j])
            all_K.append(kernels[j])
            all_Y.append(target[j])

        if n_loaded % 50 == 0:
            print(f"  {n_loaded} proteins")

    print(f"Loaded {n_loaded} proteins")

    # ── Per-element ridge ────────────────────────────────────────────

    print(f"\n{'='*70}")
    print(f"  Per-Element Ridge Calibration — {n_core} Core Kernels")
    print(f"{'='*70}")

    summary = {}

    for e in list(ELEMENTS.keys()) + [0]:
        if e == 0:
            K_arr = np.array(all_K)
            Y_arr = np.array(all_Y)
            name = "all"
        else:
            if len(elem_K[e]) < 100:
                continue
            K_arr = np.array(elem_K[e])
            Y_arr = np.array(elem_Y[e])
            name = ELEMENTS[e]

        N = K_arr.shape[0]
        X = K_arr.reshape(N, -1)  # (N, 55*5)
        w, pred, r2 = _ridge_coeffs(X, Y_arr, lam)

        # Compute per-kernel "importance": norm of the 5 coefficients
        # for each kernel (each kernel has 5 columns in the X matrix)
        kernel_importance = np.zeros(n_core)
        for ki in range(n_core):
            w_k = w[ki*5:(ki+1)*5, :]  # (5, 5) block
            kernel_importance[ki] = np.linalg.norm(w_k)

        # Sort by importance
        order = np.argsort(-kernel_importance)

        # Forward selection (top 20)
        fwd = _forward_r2(K_arr, Y_arr, lam, n_core, core_names)

        print(f"\n  {name} ({N:,d} atoms): ridge R² = {r2:.4f}")
        print(f"  Top 10 by |coefficient|:")
        for rank, ki in enumerate(order[:10], 1):
            print(f"    {rank:2d}. {core_names[ki]:25s}  "
                  f"|w|={kernel_importance[ki]:.4f}")

        print(f"  Forward selection (cumulative R²):")
        for rank, (kn, delta, cum) in enumerate(fwd[:10], 1):
            print(f"    {rank:2d}. {kn:25s}  +{delta:.4f}  R²={cum:.4f}")

        # How many kernels to reach 90% of full R²?
        threshold_90 = r2 * 0.90
        n_for_90 = n_core
        for rank, (_, _, cum) in enumerate(fwd, 1):
            if cum >= threshold_90:
                n_for_90 = rank
                break

        print(f"  Kernels for 90% of R²: {n_for_90}")

        summary[name] = {
            "n_atoms": N,
            "r2": round(r2, 4),
            "n_for_90pct": n_for_90,
            "top10_coeff": [(core_names[ki], round(float(kernel_importance[ki]), 5))
                           for ki in order[:10]],
            "forward_selection": [(kn, round(delta, 4), round(cum, 4))
                                 for kn, delta, cum in fwd],
        }

        # Write per-element coefficient table
        with open(out / f"coefficients_{name}.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["rank", "kernel", "importance", "fwd_delta", "fwd_cumulative_r2"])
            for rank in range(n_core):
                ki = order[rank]
                fwd_delta = ""
                fwd_cum = ""
                for fkn, fd, fc in fwd:
                    if fkn == core_names[ki]:
                        fwd_delta = f"{fd:.4f}"
                        fwd_cum = f"{fc:.4f}"
                        break
                writer.writerow([
                    rank + 1, core_names[ki],
                    f"{kernel_importance[ki]:.5f}",
                    fwd_delta, fwd_cum
                ])

    # ── Comparison table ─────────────────────────────────────────────

    print(f"\n{'='*70}")
    print(f"  Comparison: 55-kernel ridge vs 121-kernel ridge vs gated MLP")
    print(f"{'='*70}")
    print(f"  {'Element':8s} {'55-core':>10s} {'n_for_90%':>10s}")
    for name in ["H", "C", "N", "O", "all"]:
        if name in summary:
            s = summary[name]
            print(f"  {name:8s} {s['r2']:10.4f} {s['n_for_90pct']:10d}")

    weighted = sum(
        summary[ELEMENTS[e]]["r2"] * summary[ELEMENTS[e]]["n_atoms"]
        for e in ELEMENTS if ELEMENTS[e] in summary
    ) / sum(
        summary[ELEMENTS[e]]["n_atoms"]
        for e in ELEMENTS if ELEMENTS[e] in summary
    )
    print(f"  {'weighted':8s} {weighted:10.4f}")

    summary["weighted_r2"] = round(weighted, 4)

    with open(out / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nResults in {out}/")
