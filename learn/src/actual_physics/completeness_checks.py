#!/usr/bin/env python3
"""Completeness checks for Stage 1 thesis.

1. Nonlinear check: kernel ridge (RBF) vs linear ridge per element.
2. Leave-protein-out cross-validation.
3. Coefficient bootstrap stability.
4. Residual structure: per-protein R² distribution.

Usage:
    cd learn/src
    python3 -c "
    from mutation_set.config import load_config
    import sys; cfg = load_config('calibration.toml')
    sys.path.insert(0, str(cfg.paths.sdk))
    from actual_physics.completeness_checks import run
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


ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O"}


def ridge_r2(X, y, lam=1e-2):
    XtX = X.T @ X + lam * np.eye(X.shape[1])
    w = np.linalg.solve(XtX, X.T @ y)
    pred = X @ w
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean(axis=0)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0


def ridge_fit(X, y, lam=1e-2):
    XtX = X.T @ X + lam * np.eye(X.shape[1])
    w = np.linalg.solve(XtX, X.T @ y)
    return w


def run(cfg: Config, max_proteins: int = 0):
    out = Path("output/actual_physics/completeness")
    out.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    n_core = layout.efg_end
    lam = getattr(cfg.secondary, 'ridge_lambda', 1e-2)

    # ── Load data per protein (need protein grouping for LPOCV) ──
    protein_data = []  # list of (pid, element_array, K_norm, Y)

    print("Loading proteins...")
    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx

        kernels = assemble_kernels(p, idx, layout)[:, :n_core, :]
        K_norm = kernels.copy()
        K_norm, _ = normalize_kernels(K_norm, cfg.normalization.kernel_std_floor)
        target = p.delta.shielding.T2[idx]

        protein_data.append({
            "pid": rec.protein_id if hasattr(rec, 'protein_id') else str(n_loaded),
            "elem": rec.element,
            "K": K_norm,
            "Y": target,
        })

        if n_loaded % 100 == 0:
            print(f"  {n_loaded} proteins")

    print(f"Loaded {n_loaded} proteins\n")

    results = {}

    for e, en in ELEMENTS.items():
        # Collect per-protein data for this element
        prot_blocks = []
        for pd in protein_data:
            mask = pd["elem"] == e
            if mask.sum() < 5:
                continue
            prot_blocks.append({
                "pid": pd["pid"],
                "K": pd["K"][mask],
                "Y": pd["Y"][mask],
                "n": int(mask.sum()),
            })

        if len(prot_blocks) < 20:
            continue

        total_atoms = sum(pb["n"] for pb in prot_blocks)
        print("=" * 70)
        print(f"  {en}: {total_atoms:,d} atoms across {len(prot_blocks)} proteins")
        print("=" * 70)

        # Assemble full arrays
        K_all = np.vstack([pb["K"] for pb in prot_blocks])
        Y_all = np.vstack([pb["Y"] for pb in prot_blocks])
        X_all = K_all.reshape(total_atoms, -1)

        er = {"n_atoms": total_atoms, "n_proteins": len(prot_blocks)}

        # ── 1. Leave-protein-out CV ──────────────────────────────
        print(f"\n  Leave-protein-out cross-validation...")
        n_folds = min(len(prot_blocks), 50)  # max 50 folds for speed
        fold_size = len(prot_blocks) // n_folds
        np.random.seed(42)
        prot_order = np.random.permutation(len(prot_blocks))

        lpocv_r2s = []
        for fold in range(n_folds):
            test_idx = prot_order[fold * fold_size:(fold + 1) * fold_size]
            train_idx = np.setdiff1d(np.arange(len(prot_blocks)), test_idx)

            K_train = np.vstack([prot_blocks[i]["K"] for i in train_idx])
            Y_train = np.vstack([prot_blocks[i]["Y"] for i in train_idx])
            K_test = np.vstack([prot_blocks[i]["K"] for i in test_idx])
            Y_test = np.vstack([prot_blocks[i]["Y"] for i in test_idx])

            X_train = K_train.reshape(K_train.shape[0], -1)
            X_test = K_test.reshape(K_test.shape[0], -1)

            w = ridge_fit(X_train, Y_train, lam)
            pred = X_test @ w
            ss_res = np.sum((Y_test - pred) ** 2)
            ss_tot = np.sum((Y_test - Y_test.mean(axis=0)) ** 2)
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
            lpocv_r2s.append(r2)

        full_r2 = ridge_r2(X_all, Y_all, lam)
        cv_mean = np.mean(lpocv_r2s)
        cv_std = np.std(lpocv_r2s)
        cv_min = np.min(lpocv_r2s)
        cv_max = np.max(lpocv_r2s)

        print(f"    Full-data R²:    {full_r2:.4f}")
        print(f"    LPOCV R²:        {cv_mean:.4f} ± {cv_std:.4f}")
        print(f"    LPOCV range:     [{cv_min:.4f}, {cv_max:.4f}]")
        print(f"    Overfit gap:     {full_r2 - cv_mean:+.4f}")

        er["lpocv"] = {
            "full_r2": round(full_r2, 4),
            "cv_mean": round(cv_mean, 4),
            "cv_std": round(cv_std, 4),
            "cv_min": round(cv_min, 4),
            "cv_max": round(cv_max, 4),
            "n_folds": n_folds,
        }

        # ── 2. Per-protein R² distribution ───────────────────────
        print(f"\n  Per-protein R² distribution...")
        per_prot_r2 = []
        for pb in prot_blocks:
            X_p = pb["K"].reshape(pb["n"], -1)
            # Use coefficients from full fit
            w_full = ridge_fit(X_all, Y_all, lam)
            pred_p = X_p @ w_full
            ss_res = np.sum((pb["Y"] - pred_p) ** 2)
            ss_tot = np.sum((pb["Y"] - pb["Y"].mean(axis=0)) ** 2)
            r2_p = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
            per_prot_r2.append(r2_p)

        per_prot_r2 = np.array(per_prot_r2)
        print(f"    Median:  {np.median(per_prot_r2):.4f}")
        print(f"    Mean:    {np.mean(per_prot_r2):.4f}")
        print(f"    Std:     {np.std(per_prot_r2):.4f}")
        print(f"    Range:   [{np.min(per_prot_r2):.4f}, {np.max(per_prot_r2):.4f}]")
        print(f"    <0:      {(per_prot_r2 < 0).sum()} proteins")
        print(f"    >0.5:    {(per_prot_r2 > 0.5).sum()} proteins")

        er["per_protein_r2"] = {
            "median": round(np.median(per_prot_r2), 4),
            "mean": round(np.mean(per_prot_r2), 4),
            "std": round(np.std(per_prot_r2), 4),
            "min": round(float(np.min(per_prot_r2)), 4),
            "max": round(float(np.max(per_prot_r2)), 4),
            "n_negative": int((per_prot_r2 < 0).sum()),
            "n_above_05": int((per_prot_r2 > 0.5).sum()),
        }

        # ── 3. Coefficient bootstrap ────────────────────────────
        print(f"\n  Coefficient bootstrap (100 resamples of proteins)...")
        n_boot = 100
        boot_coeffs = []
        np.random.seed(42)
        for b in range(n_boot):
            boot_idx = np.random.choice(len(prot_blocks), size=len(prot_blocks), replace=True)
            K_boot = np.vstack([prot_blocks[i]["K"] for i in boot_idx])
            Y_boot = np.vstack([prot_blocks[i]["Y"] for i in boot_idx])
            X_boot = K_boot.reshape(K_boot.shape[0], -1)
            w_boot = ridge_fit(X_boot, Y_boot, lam)
            boot_coeffs.append(w_boot)

        boot_coeffs = np.array(boot_coeffs)  # (100, n_core*5, 5)
        # Per-kernel importance: Frobenius norm of the 5x5 weight block
        boot_importance = np.zeros((n_boot, n_core))
        for ki in range(n_core):
            block = boot_coeffs[:, ki*5:(ki+1)*5, :]  # (100, 5, 5)
            boot_importance[:, ki] = np.linalg.norm(block.reshape(n_boot, -1), axis=1)

        # Full-data importance
        w_full_reshaped = w_full.reshape(n_core, 5, 5)
        full_importance = np.array([np.linalg.norm(w_full_reshaped[ki]) for ki in range(n_core)])
        order = np.argsort(-full_importance)

        print(f"    {'rank':>4s}  {'kernel':>25s}  {'importance':>10s}  {'boot_mean':>10s}  {'boot_std':>10s}  {'CI95':>20s}  {'stable':>6s}")
        n_stable = 0
        for rank, ki in enumerate(order[:20], 1):
            imp = full_importance[ki]
            bm = np.mean(boot_importance[:, ki])
            bs = np.std(boot_importance[:, ki])
            lo = np.percentile(boot_importance[:, ki], 2.5)
            hi = np.percentile(boot_importance[:, ki], 97.5)
            stable = "yes" if lo > 0.01 * hi else "no"
            if lo > 0.01 * hi:
                n_stable += 1
            print(f"    {rank:4d}  {layout.names[ki]:>25s}  {imp:10.5f}  {bm:10.5f}  {bs:10.5f}  [{lo:.4f}, {hi:.4f}]  {stable:>6s}")

        er["bootstrap"] = {
            "n_resamples": n_boot,
            "n_stable_top20": n_stable,
        }

        # ── 4. Nonlinear check (subsample for speed) ────────────
        print(f"\n  Nonlinear check (random forest on 50K subsample)...")
        n_sub = min(50000, total_atoms)
        sub_idx = np.random.choice(total_atoms, size=n_sub, replace=False)
        X_sub = X_all[sub_idx]
        Y_sub = Y_all[sub_idx]

        # Train/test split
        n_test = n_sub // 5
        X_train_nl, X_test_nl = X_sub[n_test:], X_sub[:n_test]
        Y_train_nl, Y_test_nl = Y_sub[n_test:], Y_sub[:n_test]

        # Linear ridge
        w_lin = ridge_fit(X_train_nl, Y_train_nl, lam)
        pred_lin = X_test_nl @ w_lin
        ss_res_lin = np.sum((Y_test_nl - pred_lin) ** 2)
        ss_tot_nl = np.sum((Y_test_nl - Y_test_nl.mean(0)) ** 2)
        r2_lin = 1 - ss_res_lin / ss_tot_nl

        # Random forest (per T2 component)
        try:
            from sklearn.ensemble import RandomForestRegressor
            r2_rf_components = []
            for comp in range(5):
                rf = RandomForestRegressor(n_estimators=50, max_depth=15,
                                           n_jobs=-1, random_state=42)
                rf.fit(X_train_nl, Y_train_nl[:, comp])
                pred_rf = rf.predict(X_test_nl)
                ss_res_rf = np.sum((Y_test_nl[:, comp] - pred_rf) ** 2)
                ss_tot_comp = np.sum((Y_test_nl[:, comp] - Y_test_nl[:, comp].mean()) ** 2)
                r2_rf_components.append(1 - ss_res_rf / ss_tot_comp if ss_tot_comp > 1e-12 else 0)
            r2_rf = np.mean(r2_rf_components)

            print(f"    Linear ridge (test): {r2_lin:.4f}")
            print(f"    Random forest (test, mean over 5 T2 components): {r2_rf:.4f}")
            print(f"    RF - Ridge: {r2_rf - r2_lin:+.4f}")
            if r2_rf > r2_lin + 0.01:
                print(f"    *** Nonlinear signal detected ({r2_rf - r2_lin:+.3f}) ***")
            else:
                print(f"    No significant nonlinear signal.")

            er["nonlinear"] = {
                "ridge_test": round(r2_lin, 4),
                "rf_test": round(r2_rf, 4),
                "rf_per_component": [round(r, 4) for r in r2_rf_components],
                "nonlinear_signal": bool(r2_rf > r2_lin + 0.01),
            }
        except ImportError:
            print(f"    sklearn not available — skipping RF")
            er["nonlinear"] = {"error": "sklearn not available"}

        results[en] = er
        print()

    # ── Save ─────────────────────────────────────────────────────
    with open(out / "completeness_checks.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"Results in {out}/")
