#!/usr/bin/env python3
"""Effective dimensionality of the geometry-only kernel space.

Two questions:
  1. Is the high dimensionality after normalization real signal or
     noise in the far eigenvalues?  Test: PCA-then-ridge — project
     onto top-k PCs, fit ridge, see where R² plateaus.

  2. Does the effective dimensionality depend on distance to the
     nearest ring?  Test: stratify by distance band, compute
     eigenspectrum and PCA-ridge per band.

Both analyses on geometry-only kernels, per element.

Usage:
    cd learn/src
    python3 -c "
    from mutation_set.config import load_config
    import sys; cfg = load_config('calibration.toml')
    sys.path.insert(0, str(cfg.paths.sdk))
    from actual_physics.dimensionality_test import run
    run(cfg)
    "
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from secondary.loader import iter_proteins, ELEMENT_NAMES
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels, normalize_kernels


# ── Reuse partition from geometry_only_basis ─────────────────────────

_MOPAC_APBS_PREFIXES = (
    "MopacMC_", "MopacCoulomb_total", "MopacMC_total",
    "MopacEFG_", "APBS_", "DeltaAPBS_",
)

def _is_mopac_apbs(name):
    for prefix in _MOPAC_APBS_PREFIXES:
        if name.startswith(prefix) or name == prefix:
            return True
    return False

def _geo_indices(layout):
    return [i for i in range(layout.efg_end)
            if not _is_mopac_apbs(layout.names[i])]


def _ridge_r2(X, y, lam=1e-2):
    if X.shape[0] < X.shape[1] + 10:
        return float('nan')
    XtX = X.T @ X + lam * np.eye(X.shape[1])
    w = np.linalg.solve(XtX, X.T @ y)
    pred = X @ w
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean(axis=0)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0


def _pca_ridge_curve(K_arr, Y_arr, lam, max_k=50):
    """Project onto top-k PCs, fit ridge, return [(k, R²)].

    Uses train/val split (80/20) to detect overfitting."""
    N = K_arr.shape[0]
    X = K_arr.reshape(N, -1)
    # Center
    mu = X.mean(axis=0)
    Xc = X - mu

    # PCA via SVD (more stable for wide matrices)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)

    # Train/val split (deterministic)
    n_val = max(50, N // 5)
    n_train = N - n_val
    idx = np.arange(N)
    np.random.seed(42)
    np.random.shuffle(idx)
    train, val = idx[:n_train], idx[n_train:]

    results = []
    for k in range(1, min(max_k + 1, len(S) + 1)):
        # Project onto top-k components
        V_k = Vt[:k, :].T  # (D, k)
        Z_train = Xc[train] @ V_k  # (n_train, k)
        Z_val = Xc[val] @ V_k      # (n_val, k)

        # Ridge on projections
        Y_train = Y_arr[train]
        Y_val = Y_arr[val]

        XtX = Z_train.T @ Z_train + lam * np.eye(k)
        w = np.linalg.solve(XtX, Z_train.T @ Y_train)

        pred_val = Z_val @ w
        ss_res = np.sum((Y_val - pred_val) ** 2)
        ss_tot = np.sum((Y_val - Y_val.mean(axis=0)) ** 2)
        r2_val = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0

        # Also train R² for comparison (detect overfitting)
        pred_train = Z_train @ w
        ss_res_t = np.sum((Y_train - pred_train) ** 2)
        ss_tot_t = np.sum((Y_train - Y_train.mean(axis=0)) ** 2)
        r2_train = 1.0 - ss_res_t / ss_tot_t if ss_tot_t > 1e-12 else 0.0

        results.append((k, r2_train, r2_val))

    return results


def run(cfg: Config, max_proteins: int = 0):
    out = Path("output/actual_physics/dimensionality")
    out.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    geo_idx = _geo_indices(layout)
    geo_names = [layout.names[i] for i in geo_idx]
    n_geo = len(geo_idx)
    lam = getattr(cfg.secondary, 'ridge_lambda', 1e-2)

    ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O"}
    DISTANCE_BANDS = [(0, 4), (4, 8), (8, 12), (12, 999)]
    BAND_NAMES = ["0-4A", "4-8A", "8-12A", "12+A"]

    # Accumulators
    elem_data = {e: {"K_norm": [], "K_raw": [], "Y": [], "dist": []}
                 for e in ELEMENTS}

    print(f"Geometry-only kernels: {n_geo}")
    print("Loading proteins...")
    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx
        M = len(idx)

        kernels_full = assemble_kernels(p, idx, layout)[:, :layout.efg_end, :]
        K_raw = kernels_full[:, geo_idx, :]

        K_norm_copy = K_raw.copy()
        K_norm, _ = normalize_kernels(K_norm_copy, cfg.normalization.kernel_std_floor)

        target = p.delta.shielding.T2[idx]
        dist = rec.ring_dist  # distance to nearest removed ring

        for j in range(M):
            e = int(rec.element[j])
            if e not in elem_data:
                continue
            elem_data[e]["K_norm"].append(K_norm[j])
            elem_data[e]["K_raw"].append(K_raw[j])
            elem_data[e]["Y"].append(target[j])
            elem_data[e]["dist"].append(dist[j])

        if n_loaded % 50 == 0:
            print(f"  {n_loaded} proteins")

    print(f"Loaded {n_loaded} proteins\n")

    results = {}

    for e, ename in ELEMENTS.items():
        d = elem_data[e]
        if len(d["K_norm"]) < 200:
            continue

        K_norm = np.array(d["K_norm"])   # (N, 44, 5)
        K_raw = np.array(d["K_raw"])
        Y = np.array(d["Y"])
        dist = np.array(d["dist"])
        N = K_norm.shape[0]

        print("=" * 70)
        print(f"  {ename} ({N:,d} atoms)")
        print("=" * 70)

        er = {"n_atoms": N}

        # ── 1. PCA-then-ridge: where does R² plateau? ────────────

        print(f"\n  PCA-then-ridge (normalized, 80/20 split):")
        curve_norm = _pca_ridge_curve(K_norm, Y, lam, max_k=50)
        print(f"    {'k':>4s}  {'train':>8s}  {'val':>8s}  {'gap':>8s}")
        milestones = [1, 2, 3, 5, 8, 10, 15, 20, 30, 40, 50]
        for k, r2t, r2v in curve_norm:
            if k in milestones:
                print(f"    {k:4d}  {r2t:8.4f}  {r2v:8.4f}  {r2t-r2v:+8.4f}")

        # Find where val R² plateaus (less than 0.005 improvement per step)
        plateau_k = len(curve_norm)
        for i in range(2, len(curve_norm)):
            k, _, r2v = curve_norm[i]
            _, _, r2v_prev = curve_norm[i-1]
            if r2v - r2v_prev < 0.002 and all(
                curve_norm[j][2] - curve_norm[j-1][2] < 0.002
                for j in range(i, min(i+3, len(curve_norm)))
            ):
                plateau_k = k
                break
        print(f"    Plateau at k ≈ {plateau_k} "
              f"(val R² = {curve_norm[plateau_k-1][2]:.4f})")

        er["pca_ridge"] = {
            "curve": [(k, round(r2t, 4), round(r2v, 4))
                      for k, r2t, r2v in curve_norm],
            "plateau_k": plateau_k,
        }

        # Same on RAW (unnormalized) kernels
        print(f"\n  PCA-then-ridge (raw/unnormalized):")
        curve_raw = _pca_ridge_curve(K_raw, Y, lam, max_k=50)
        for k, r2t, r2v in curve_raw:
            if k in milestones:
                print(f"    {k:4d}  {r2t:8.4f}  {r2v:8.4f}  {r2t-r2v:+8.4f}")

        plateau_raw = len(curve_raw)
        for i in range(2, len(curve_raw)):
            k, _, r2v = curve_raw[i]
            _, _, r2v_prev = curve_raw[i-1]
            if r2v - r2v_prev < 0.002 and all(
                curve_raw[j][2] - curve_raw[j-1][2] < 0.002
                for j in range(i, min(i+3, len(curve_raw)))
            ):
                plateau_raw = k
                break
        print(f"    Plateau at k ≈ {plateau_raw}")
        er["pca_ridge_raw_plateau"] = plateau_raw

        # ── 2. Distance-stratified eigenspectrum ──────────────────

        print(f"\n  Distance-stratified analysis:")
        band_results = {}
        for (lo, hi), bname in zip(DISTANCE_BANDS, BAND_NAMES):
            mask = (dist >= lo) & (dist < hi) & np.isfinite(dist)
            n_band = mask.sum()
            if n_band < 100:
                print(f"    {bname:8s} ({n_band:5d} atoms): too few")
                continue

            K_band = K_raw[mask]  # use RAW here — distance is the organising variable
            Y_band = Y[mask]

            X_flat = K_band.reshape(n_band, -1)
            cov = np.cov(X_flat, rowvar=False)
            eigvals = np.linalg.eigvalsh(cov)[::-1]
            eigvals = np.maximum(eigvals, 0)
            total = eigvals.sum()
            if total < 1e-15:
                continue
            cum = np.cumsum(eigvals) / total

            dim90 = int(np.searchsorted(cum, 0.90)) + 1
            dim95 = int(np.searchsorted(cum, 0.95)) + 1
            dim98 = int(np.searchsorted(cum, 0.98)) + 1

            # Ridge R² in this band
            r2_band = _ridge_r2(X_flat, Y_band, lam)

            # PCA-ridge curve (abbreviated, up to 30)
            if n_band >= 200:
                curve_band = _pca_ridge_curve(K_band, Y_band, lam, max_k=30)
                # Find effective predictive dims
                pred_plateau = len(curve_band)
                for i in range(2, len(curve_band)):
                    _, _, r2v = curve_band[i]
                    _, _, r2v_prev = curve_band[i-1]
                    if r2v - r2v_prev < 0.003 and i + 2 < len(curve_band) and all(
                        curve_band[j][2] - curve_band[j-1][2] < 0.003
                        for j in range(i, min(i+3, len(curve_band)))
                    ):
                        pred_plateau = curve_band[i][0]
                        break
            else:
                pred_plateau = -1

            # Top eigenvalue concentration
            top5_frac = sum(eigvals[:5]) / total
            top10_frac = sum(eigvals[:10]) / total

            print(f"    {bname:8s} ({n_band:5d} atoms): "
                  f"dims@98%={dim98:3d}  top5={top5_frac:.3f}  "
                  f"top10={top10_frac:.3f}  "
                  f"R²={r2_band:.4f}  pred_dims≈{pred_plateau}")

            band_results[bname] = {
                "n_atoms": int(n_band),
                "dim90": dim90, "dim95": dim95, "dim98": dim98,
                "top5_frac": round(top5_frac, 4),
                "top10_frac": round(top10_frac, 4),
                "ridge_r2": round(r2_band, 4),
                "predictive_dims": pred_plateau,
            }

        er["distance_bands"] = band_results

        # ── 3. Eigenspectrum: raw vs normalized comparison ────────

        X_raw = K_raw.reshape(N, -1)
        X_norm = K_norm.reshape(N, -1)

        cov_raw = np.cov(X_raw, rowvar=False)
        cov_norm = np.cov(X_norm, rowvar=False)

        eig_raw = np.linalg.eigvalsh(cov_raw)[::-1]
        eig_norm = np.linalg.eigvalsh(cov_norm)[::-1]
        eig_raw = np.maximum(eig_raw, 0)
        eig_norm = np.maximum(eig_norm, 0)
        tot_raw = eig_raw.sum()
        tot_norm = eig_norm.sum()

        print(f"\n  Eigenspectrum comparison (raw vs normalized):")
        print(f"    {'rank':>4s}  {'raw':>8s}  {'norm':>8s}  {'raw_cum':>8s}  {'norm_cum':>8s}")
        cum_raw = 0
        cum_norm = 0
        for k in range(min(20, len(eig_raw))):
            fr = eig_raw[k] / tot_raw if tot_raw > 0 else 0
            fn = eig_norm[k] / tot_norm if tot_norm > 0 else 0
            cum_raw += fr
            cum_norm += fn
            if k < 15 or k + 1 in [20]:
                print(f"    {k+1:4d}  {fr:8.4f}  {fn:8.4f}  {cum_raw:8.4f}  {cum_norm:8.4f}")

        er["eigenspectrum_comparison"] = {
            "raw_top10_cumulative": round(sum(eig_raw[:10]) / tot_raw, 4) if tot_raw > 0 else 0,
            "norm_top10_cumulative": round(sum(eig_norm[:10]) / tot_norm, 4) if tot_norm > 0 else 0,
        }

        results[ename] = er

    # ── Cross-element summary ───────────────────────────────────────

    print(f"\n{'=' * 70}")
    print(f"  Summary: Predictive Dimensionality")
    print(f"{'=' * 70}")
    print(f"  {'elem':>5s}  {'norm_plateau':>12s}  {'raw_plateau':>12s}  "
          f"{'near_dims':>10s}  {'far_dims':>10s}")
    for ename in ["H", "C", "N", "O"]:
        if ename not in results:
            continue
        r = results[ename]
        bands = r.get("distance_bands", {})
        near = bands.get("0-4A", {}).get("predictive_dims", "?")
        far = bands.get("8-12A", {}).get("predictive_dims", "?")
        print(f"  {ename:>5s}  {r['pca_ridge']['plateau_k']:>12d}  "
              f"{r['pca_ridge_raw_plateau']:>12d}  "
              f"{str(near):>10s}  {str(far):>10s}")

    with open(out / "dimensionality.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults in {out}/")
