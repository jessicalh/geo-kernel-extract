#!/usr/bin/env python3
"""
Level A: Linear kernel regression.

The physical model:
    DFT_T2[m] = Σ_k  w_k × G_k_T2[m]     for m = -2,...,+2

Each calculator contributes ONE scalar weight to ALL 5 T2 components.
This is equivariant by construction — scalar × L=2 = L=2.

The learned weights ARE the physical parameters:
    w_BS_PHE  →  ring current intensity for PHE  (nA × PPM_FACTOR)
    w_MC_CO   →  Δχ for C=O bonds  (Å³ × ppm)
    w_Coulomb →  Buckingham γ  (ppm / (V/Å²))

The residual after this fit is the T2 that classical kernels cannot explain.
"""

import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))
from load import (load_conformation, list_proteins, T2, T0,
                  FEATURES_DIR, RING_TYPES)


def build_features_and_targets(protein_ids):
    """Build the stacked regression matrix.

    Returns:
        X: (N_total * 5, K) — each row is one T2 component of all K kernels for one atom
        y: (N_total * 5,) — corresponding DFT T2 component
        meta: list of (protein_id, atom_index, t2_component) per row
        feature_names: list of K kernel names
    """
    # Feature columns: per-ring-type T2 for ring calculators,
    # per-category T2 for McConnell, total T2 for Coulomb/HBond/RingSusc
    feature_names = []
    # BS: 8 ring types
    for rt in RING_TYPES:
        feature_names.append(f"BS_{rt}")
    # HM: 8 ring types
    for rt in RING_TYPES:
        feature_names.append(f"HM_{rt}")
    # PiQuad: 8 ring types
    for rt in RING_TYPES:
        feature_names.append(f"PQ_{rt}")
    # Dispersion: 8 ring types
    for rt in RING_TYPES:
        feature_names.append(f"Disp_{rt}")
    # McConnell: 5 bond categories
    mc_cats = ["CO", "CN", "BB_other", "SC_CO", "Aromatic"]
    for cat in mc_cats:
        feature_names.append(f"MC_{cat}")
    # Coulomb EFG total
    feature_names.append("Coulomb")
    # RingSusc total
    feature_names.append("RingSusc")
    # HBond total
    feature_names.append("HBond")

    K = len(feature_names)

    all_X = []
    all_y = []
    all_meta = []

    for pid in protein_ids:
        wt = load_conformation(pid, "wt")
        if wt is None or wt.orca is None:
            continue

        N = wt.n_atoms
        dft_t2 = T2(wt.orca)  # (N, 5)

        # Build per-atom feature matrix: (N, K, 5)
        # Each feature k gives a T2 vector (5 components) per atom
        features = np.zeros((N, K, 5))

        # BS per-type: (N, 40) = 8 types × 5 components
        for t in range(8):
            features[:, t, :] = wt.bs_type_T2[:, t*5:(t+1)*5]
        # HM per-type
        for t in range(8):
            features[:, 8+t, :] = wt.hm_type_T2[:, t*5:(t+1)*5]
        # PQ per-type
        for t in range(8):
            features[:, 16+t, :] = wt.pq_type_T2[:, t*5:(t+1)*5]
        # Disp per-type
        for t in range(8):
            features[:, 24+t, :] = wt.disp_type_T2[:, t*5:(t+1)*5]
        # McConnell categories: (N, 25) = 5 cats × 5 components
        for c in range(5):
            features[:, 32+c, :] = wt.mc_category_T2[:, c*5:(c+1)*5]
        # Coulomb total T2
        features[:, 37, :] = T2(wt.coulomb)
        # RingSusc total T2
        features[:, 38, :] = T2(wt.ringchi)
        # HBond total T2
        features[:, 39, :] = T2(wt.hbond)

        # Stack: for equivariant regression, each T2 component is a separate observation
        # with the SAME weight. X shape: (N*5, K), y shape: (N*5,)
        for m in range(5):
            all_X.append(features[:, :, m])  # (N, K)
            all_y.append(dft_t2[:, m])        # (N,)
            for i in range(N):
                all_meta.append((pid, i, m))

    X = np.vstack(all_X)
    y = np.concatenate(all_y)

    return X, y, all_meta, feature_names


def fit_and_report(protein_ids, alpha=1.0):
    """Fit ridge regression and report results."""
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import cross_val_score

    print(f"Building features from {len(protein_ids)} proteins...")
    X, y, meta, names = build_features_and_targets(protein_ids)
    K = len(names)
    n_atoms = X.shape[0] // 5

    print(f"  Atoms: {n_atoms}")
    print(f"  Observations: {X.shape[0]} (atoms × 5 T2 components)")
    print(f"  Features: {K}")
    print()

    # Feature scale report (what the regression is working with)
    print("Feature scales (std of each column):")
    stds = X.std(axis=0)
    for i, name in enumerate(names):
        if stds[i] > 1e-10:
            print(f"  {name:15s}: std={stds[i]:.6f}")
    print()

    # Fit
    model = Ridge(alpha=alpha, fit_intercept=False)
    model.fit(X, y)

    y_pred = model.predict(X)
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot

    # Per-T2-component R²
    for m in range(5):
        idx = np.arange(m, len(y), 5)
        ym = y[idx]
        yp = y_pred[idx]
        r2m = 1.0 - np.sum((ym - yp)**2) / np.sum((ym - ym.mean())**2)
        print(f"  T2 component m={m-2:+d}: R² = {r2m:.4f}")

    print(f"\n  Overall T2 R² = {r2:.4f}  (alpha={alpha})")
    print()

    # Learned weights = physical parameters
    print("Learned weights (= physical parameters in kernel units → ppm):")
    print(f"  {'Feature':15s} {'Weight':>12s} {'|Weight|×std':>12s}")
    print("  " + "-" * 42)
    contributions = np.abs(model.coef_) * stds
    order = np.argsort(contributions)[::-1]
    for i in order:
        if stds[i] < 1e-10:
            continue
        print(f"  {names[i]:15s} {model.coef_[i]:+12.4f} {contributions[i]:12.4f}")
    print()

    # Residual analysis
    resid = y - y_pred
    resid_per_atom = resid.reshape(5, n_atoms).T  # (n_atoms, 5) back to T2 shape
    # Wait, the stacking interleaves atoms and components differently
    # Let me reshape correctly
    resid_by_m = []
    for m in range(5):
        idx = np.arange(m * n_atoms, (m+1) * n_atoms)
        resid_by_m.append(resid[idx])
    resid_per_atom = np.column_stack(resid_by_m)  # (n_atoms, 5)

    resid_norms = np.sqrt(np.sum(resid_per_atom**2, axis=1))
    dft_t2_norms_flat = []
    for m in range(5):
        idx = np.arange(m * n_atoms, (m+1) * n_atoms)
        dft_t2_norms_flat.append(y[idx])
    dft_per_atom = np.column_stack(dft_t2_norms_flat)
    dft_norms = np.sqrt(np.sum(dft_per_atom**2, axis=1))

    print(f"Residual |T2| after fit:")
    print(f"  Mean:   {resid_norms.mean():.4f} ppm")
    print(f"  Median: {np.median(resid_norms):.4f} ppm")
    print(f"  Max:    {resid_norms.max():.4f} ppm")
    print(f"  DFT |T2| mean: {dft_norms.mean():.4f} ppm")
    print(f"  Fraction explained: {1.0 - resid_norms.mean()/dft_norms.mean():.4f}")
    print()

    return model, names, X, y, meta


def main():
    proteins = list_proteins()
    if not proteins:
        print("No complete pairs extracted yet.")
        return

    print(f"{'='*60}")
    print(f"LEVEL A: LINEAR KERNEL REGRESSION")
    print(f"{'='*60}\n")

    # Use all available proteins (train on everything for now)
    model, names, X, y, meta = fit_and_report(proteins, alpha=1.0)

    # Also try different regularisation
    print(f"\n{'='*60}")
    print(f"REGULARISATION SWEEP")
    print(f"{'='*60}\n")
    for alpha in [0.01, 0.1, 1.0, 10.0, 100.0]:
        from sklearn.linear_model import Ridge
        m = Ridge(alpha=alpha, fit_intercept=False)
        m.fit(X, y)
        yp = m.predict(X)
        r2 = 1.0 - np.sum((y - yp)**2) / np.sum((y - y.mean())**2)
        print(f"  alpha={alpha:8.2f}  R²={r2:.4f}  max|w|={np.abs(m.coef_).max():.4f}")


if __name__ == "__main__":
    main()
