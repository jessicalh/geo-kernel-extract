#!/usr/bin/env python3
"""Two complementary views of the accumulation basis.

View 1: Raw kernel space.  3 predictive dimensions in PCA.
  - What are those 3 directions?  Project the ridge coefficients
    into PCA space and see which raw-PCA dimensions carry prediction.
  - Per-kernel noise gating: mask low-magnitude entries, re-compute
    eigenspectrum.  Does the gated spectrum sharpen?

View 2: Ridge-projected space.  The calibrated coefficients define
  per-element linear maps from kernel T2 to predicted T2.
  - What is the effective rank of the coefficient matrix?
  - How much variance in the raw kernel space is predictive vs noise?
  - Per-kernel contribution to predicted variance.

No magic distances.  Kernel magnitudes sort themselves.

Usage:
    cd learn/src
    python3 -c "
    from mutation_set.config import load_config
    import sys; cfg = load_config('calibration.toml')
    sys.path.insert(0, str(cfg.paths.sdk))
    from actual_physics.accumulation_basis import run
    run(cfg)
    "
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from secondary.loader import iter_proteins, ELEMENT_NAMES, RING_TYPE_NAMES
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels


_MOPAC_APBS_PREFIXES = (
    "MopacMC_", "MopacCoulomb_total", "MopacMC_total",
    "MopacEFG_", "APBS_", "DeltaAPBS_",
)

def _is_mopac(name):
    return any(name.startswith(p) or name == p for p in _MOPAC_APBS_PREFIXES)

def _geo_indices(layout):
    return [i for i in range(layout.efg_end) if not _is_mopac(layout.names[i])]


def _ridge_coeffs(X, y, lam):
    """Returns (w, pred, R²).  X: (N, D), y: (N, 5)."""
    XtX = X.T @ X + lam * np.eye(X.shape[1])
    w = np.linalg.solve(XtX, X.T @ y)
    pred = X @ w
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean(axis=0)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
    return w, pred, r2


def run(cfg: Config, max_proteins: int = 0):
    out = Path("output/actual_physics/accumulation_basis")
    out.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    geo_idx = _geo_indices(layout)
    geo_names = [layout.names[i] for i in geo_idx]
    n_geo = len(geo_idx)
    lam = getattr(cfg.secondary, 'ridge_lambda', 1e-2)

    ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O"}

    elem_data = {e: {"K": [], "Y": []} for e in ELEMENTS}

    print(f"Geometry-only kernels: {n_geo}")
    print("Loading proteins...")
    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx
        kernels = assemble_kernels(p, idx, layout)[:, geo_idx, :]
        target = p.delta.shielding.T2[idx]
        for j in range(len(idx)):
            e = int(rec.element[j])
            if e in elem_data:
                elem_data[e]["K"].append(kernels[j])
                elem_data[e]["Y"].append(target[j])
        if n_loaded % 50 == 0:
            print(f"  {n_loaded} proteins")
    print(f"Loaded {n_loaded} proteins\n")

    results = {}

    for e, ename in ELEMENTS.items():
        d = elem_data[e]
        if len(d["K"]) < 200:
            continue

        K = np.array(d["K"])  # (N, 44, 5)
        Y = np.array(d["Y"])  # (N, 5)
        N = K.shape[0]

        print("=" * 70)
        print(f"  {ename} ({N:,d} atoms)")
        print("=" * 70)

        er = {"n_atoms": N}

        # ══════════════════════════════════════════════════════════
        # VIEW 1: RAW KERNEL SPACE
        # ══════════════════════════════════════════════════════════

        X = K.reshape(N, -1)  # (N, 220)

        # --- 1a. Raw PCA ---
        mu = X.mean(axis=0)
        Xc = X - mu
        U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
        total_var = (S ** 2).sum()
        frac = S ** 2 / total_var

        print(f"\n  RAW PCA eigenspectrum (top 15):")
        print(f"    {'PC':>4s}  {'var%':>8s}  {'cum%':>8s}")
        cum = 0
        for k in range(min(15, len(frac))):
            cum += frac[k]
            bar = "#" * int(frac[k] * 200)
            print(f"    {k+1:4d}  {frac[k]*100:7.2f}%  {cum*100:7.2f}%  {bar}")

        # --- 1b. Ridge fit on raw kernels ---
        w_raw, pred_raw, r2_raw = _ridge_coeffs(X, Y, lam)
        print(f"\n  Raw ridge: R² = {r2_raw:.4f}")

        # --- 1c. Project ridge weights into PCA space ---
        # w_raw is (220, 5).  Vt is (220, 220) orthonormal rows.
        # w_pca = Vt @ w_raw gives the ridge weight in PCA coordinates.
        w_pca = Vt @ w_raw  # (min(N,220), 5)
        # Importance of each PC for prediction: ||w_pca[k,:]||
        pc_pred_importance = np.linalg.norm(w_pca, axis=1)
        # Normalize
        total_imp = pc_pred_importance.sum()
        pc_pred_frac = pc_pred_importance / total_imp if total_imp > 0 else pc_pred_importance

        print(f"\n  Ridge weight distribution across PCA components:")
        print(f"    {'PC':>4s}  {'var%':>8s}  {'ridge%':>8s}  {'ratio':>8s}")
        for k in range(min(15, len(frac))):
            ratio = pc_pred_frac[k] / frac[k] if frac[k] > 1e-10 else float('inf')
            print(f"    {k+1:4d}  {frac[k]*100:7.2f}%  {pc_pred_frac[k]*100:7.2f}%  {ratio:8.2f}")

        # How many PCs hold 90% of ridge weight?
        sorted_imp = np.sort(pc_pred_frac)[::-1]
        cum_imp = np.cumsum(sorted_imp)
        n_for_90_ridge = int(np.searchsorted(cum_imp, 0.90)) + 1
        print(f"    PCs for 90% of ridge weight: {n_for_90_ridge}")

        er["raw_pca"] = {
            "top10_var_frac": [round(float(frac[k]), 4) for k in range(min(10, len(frac)))],
            "top10_ridge_frac": [round(float(pc_pred_frac[k]), 4) for k in range(min(10, len(frac)))],
            "pcs_for_90pct_ridge": n_for_90_ridge,
            "r2_raw": round(r2_raw, 4),
        }

        # --- 1d. Per-kernel noise gating ---
        # Compute per-kernel magnitude distribution
        # magnitude = ||K[:, ki, :]|| for each atom and kernel
        K_mag = np.linalg.norm(K, axis=2)  # (N, 44)

        print(f"\n  Per-kernel magnitude distribution:")
        print(f"    {'kernel':>25s}  {'median':>8s}  {'p25':>8s}  {'p75':>8s}  "
              f"{'frac>0':>8s}")
        mag_stats = []
        for ki in range(n_geo):
            mags = K_mag[:, ki]
            nonzero = (mags > 1e-15).mean()
            if nonzero > 0.01:
                med = np.median(mags[mags > 1e-15])
                p25 = np.percentile(mags[mags > 1e-15], 25)
                p75 = np.percentile(mags[mags > 1e-15], 75)
            else:
                med = p25 = p75 = 0
            mag_stats.append((geo_names[ki], med, p25, p75, nonzero))

        # Sort by median magnitude (highest signal first)
        mag_stats.sort(key=lambda x: -x[1])
        for name, med, p25, p75, frac_nz in mag_stats[:15]:
            print(f"    {name:>25s}  {med:8.5f}  {p25:8.5f}  {p75:8.5f}  "
                  f"{frac_nz:7.1%}")

        # Gate: zero out entries below per-kernel 25th percentile
        K_gated = K.copy()
        for ki in range(n_geo):
            mags = K_mag[:, ki]
            nz = mags > 1e-15
            if nz.sum() > 100:
                threshold = np.percentile(mags[nz], 25)
                gate = mags[:, np.newaxis] > threshold  # (N, 1) broadcast
                # Soft gate: magnitude / (magnitude + threshold)
                soft_gate = mags / (mags + threshold + 1e-30)
                K_gated[:, ki, :] *= soft_gate[:, np.newaxis]

        X_gated = K_gated.reshape(N, -1)

        # Eigenspectrum after gating
        Xc_g = X_gated - X_gated.mean(axis=0)
        _, S_g, _ = np.linalg.svd(Xc_g, full_matrices=False)
        total_g = (S_g ** 2).sum()
        frac_g = S_g ** 2 / total_g if total_g > 0 else S_g ** 2

        # Ridge on gated
        _, _, r2_gated = _ridge_coeffs(X_gated, Y, lam)

        print(f"\n  After per-kernel soft gating:")
        print(f"    Ridge R²: raw={r2_raw:.4f}, gated={r2_gated:.4f}")
        print(f"    Eigenspectrum (top 10):")
        cum_g = 0
        for k in range(10):
            cum_g += frac_g[k]
            print(f"      PC{k+1}: {frac_g[k]*100:.2f}%  cum={cum_g*100:.2f}%")

        er["gated"] = {
            "r2_raw": round(r2_raw, 4),
            "r2_gated": round(r2_gated, 4),
            "top10_var_frac": [round(float(frac_g[k]), 4) for k in range(10)],
        }

        # ══════════════════════════════════════════════════════════
        # VIEW 2: RIDGE-PROJECTED SPACE
        # ══════════════════════════════════════════════════════════

        # The ridge coefficient matrix w_raw: (220, 5) = (44*5, 5)
        # Reshape to (44, 5, 5): for each kernel, a 5x5 matrix mapping
        # kernel T2 components to predicted T2 components.
        W = w_raw.reshape(n_geo, 5, 5)  # (44, 5, 5)

        # Per-kernel predictive importance: Frobenius norm of W block
        kernel_importance = np.array([np.linalg.norm(W[ki]) for ki in range(n_geo)])
        order = np.argsort(-kernel_importance)

        print(f"\n  RIDGE COEFFICIENT STRUCTURE:")
        print(f"    Per-kernel importance (||W_k||_F):")
        print(f"    {'rank':>4s}  {'kernel':>25s}  {'||W||':>8s}  {'cum%':>8s}")
        total_w = kernel_importance.sum()
        cum_w = 0
        for rank, ki in enumerate(order[:20], 1):
            w_frac = kernel_importance[ki] / total_w
            cum_w += w_frac
            print(f"    {rank:4d}  {geo_names[ki]:>25s}  "
                  f"{kernel_importance[ki]:8.5f}  {cum_w*100:7.2f}%")

        # How many kernels for 90% of ridge weight?
        cum_sorted = np.cumsum(np.sort(kernel_importance)[::-1]) / total_w
        n_kernels_90 = int(np.searchsorted(cum_sorted, 0.90)) + 1
        print(f"    Kernels for 90% of coefficient weight: {n_kernels_90}")

        er["ridge_structure"] = {
            "top15_kernels": [(geo_names[order[i]], round(float(kernel_importance[order[i]]), 5))
                              for i in range(min(15, n_geo))],
            "kernels_for_90pct": n_kernels_90,
        }

        # --- 2b. Effective rank of the coefficient matrix ---
        # W viewed as a linear operator from (44×5) → 5.
        # SVD of w_raw: (220, 5)
        _, S_w, _ = np.linalg.svd(w_raw, full_matrices=False)
        S_w_frac = S_w / S_w.sum()
        print(f"\n    Coefficient matrix SVD (effective rank):")
        for k in range(min(5, len(S_w))):
            print(f"      SV{k+1}: {S_w[k]:.5f} ({S_w_frac[k]*100:.1f}%)")

        # Effective rank = exp(entropy of normalized singular values)
        S_norm = S_w / S_w.sum()
        eff_rank = np.exp(-np.sum(S_norm * np.log(S_norm + 1e-30)))
        print(f"      Effective rank: {eff_rank:.2f}")

        er["coeff_svd"] = {
            "singular_values": [round(float(s), 5) for s in S_w],
            "effective_rank": round(eff_rank, 2),
        }

        # --- 2c. Predicted variance decomposition ---
        # For each kernel, how much of the predicted variance does it contribute?
        # pred = sum_k K[:, k, :] @ W[k]
        # Var(pred) = sum_k sum_l Cov(K_k @ W_k, K_l @ W_l)
        # Simpler: per-kernel contribution = Var(K[:, k, :] @ W[k])
        pred_var_total = np.var(pred_raw, axis=0).sum()
        per_kernel_pred_var = np.zeros(n_geo)
        for ki in range(n_geo):
            pred_k = K[:, ki, :] @ W[ki]  # (N, 5)
            per_kernel_pred_var[ki] = np.var(pred_k, axis=0).sum()

        pv_order = np.argsort(-per_kernel_pred_var)
        print(f"\n    Predicted variance by kernel:")
        print(f"    {'rank':>4s}  {'kernel':>25s}  {'pred_var%':>10s}  {'cum%':>8s}")
        cum_pv = 0
        for rank, ki in enumerate(pv_order[:15], 1):
            frac_pv = per_kernel_pred_var[ki] / pred_var_total if pred_var_total > 0 else 0
            cum_pv += frac_pv
            print(f"    {rank:4d}  {geo_names[ki]:>25s}  "
                  f"{frac_pv*100:9.2f}%  {cum_pv*100:7.2f}%")

        er["pred_var_by_kernel"] = [
            (geo_names[pv_order[i]], round(float(per_kernel_pred_var[pv_order[i]] / pred_var_total), 4))
            for i in range(min(15, n_geo))
        ]

        # --- 2d. Cross-kernel prediction correlation ---
        # Which pairs of kernels produce correlated predictions?
        # This reveals which kernels are functionally redundant for prediction.
        pred_per_k = np.zeros((N, n_geo, 5))
        for ki in range(n_geo):
            pred_per_k[:, ki, :] = K[:, ki, :] @ W[ki]

        # Flatten to (N, 44*5) and compute correlation matrix
        ppk_flat = pred_per_k.reshape(N, -1)
        # Only compute for top-10 by predicted variance (44x44 is a lot)
        top10 = pv_order[:10]
        pred_top10 = np.zeros((N, 10, 5))
        for i, ki in enumerate(top10):
            pred_top10[:, i, :] = pred_per_k[:, ki, :]
        pt_flat = pred_top10.reshape(N, -1)
        corr = np.corrcoef(pt_flat.T)  # (50, 50)

        # Block correlation: mean |corr| between kernels (5x5 blocks)
        print(f"\n    Cross-kernel prediction correlation (top 10):")
        print(f"    {'':>14s}", end="")
        for i in range(min(10, len(top10))):
            print(f"  {geo_names[top10[i]][:6]:>6s}", end="")
        print()
        for i in range(min(10, len(top10))):
            print(f"    {geo_names[top10[i]]:>14s}", end="")
            for j in range(min(10, len(top10))):
                # Mean absolute correlation between 5x5 blocks
                block = corr[i*5:(i+1)*5, j*5:(j+1)*5]
                mean_abs = np.mean(np.abs(block))
                print(f"  {mean_abs:6.3f}", end="")
            print()

        results[ename] = er

    with open(out / "accumulation_basis.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults in {out}/")
