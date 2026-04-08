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
from load import (load_pair, list_proteins, T2, T0,
                  FEATURES_DIR, RING_TYPES)
from delta import match_atoms


def build_features_and_targets(protein_ids):
    """Build the stacked regression matrix from WT-ALA deltas.

    The mutation delta isolates the ring effect. Target is delta DFT T2,
    features are delta kernel T2s.

    Returns:
        X: (M_total * 5, K) — each row is one T2 component of all K kernels for one matched atom
        y: (M_total * 5,) — corresponding delta DFT T2 component
        meta: list of (protein_id, atom_index, t2_component) per row
        feature_names: list of K kernel names
    """
    feature_names = []
    for rt in RING_TYPES:
        feature_names.append(f"BS_{rt}")
    for rt in RING_TYPES:
        feature_names.append(f"HM_{rt}")
    for rt in RING_TYPES:
        feature_names.append(f"PQ_{rt}")
    for rt in RING_TYPES:
        feature_names.append(f"Disp_{rt}")
    mc_cats = ["CO", "CN", "BB_other", "SC_CO", "Aromatic"]
    for cat in mc_cats:
        feature_names.append(f"MC_{cat}")
    feature_names.append("Coulomb")
    feature_names.append("RingSusc")
    feature_names.append("HBond")
    feature_names.append("MopacCoulomb")
    for cat in mc_cats:
        feature_names.append(f"MopacMC_{cat}")

    K = len(feature_names)

    all_X = []
    all_y = []
    all_meta = []

    for pid in protein_ids:
        wt, ala = load_pair(pid)
        if wt is None or ala is None:
            continue
        if wt.orca is None or ala.orca is None:
            continue

        wt_idx, ala_idx, _ = match_atoms(
            wt.pos, wt.element, ala.pos, ala.element
        )
        if len(wt_idx) < 10:
            continue

        M = len(wt_idx)
        dft_t2 = T2(wt.orca[wt_idx] - ala.orca[ala_idx])  # (M, 5)

        features = np.zeros((M, K, 5))

        # Delta per-type T2s
        d_bs = wt.bs_type_T2[wt_idx] - ala.bs_type_T2[ala_idx]
        for t in range(8):
            features[:, t, :] = d_bs[:, t*5:(t+1)*5]
        d_hm = wt.hm_type_T2[wt_idx] - ala.hm_type_T2[ala_idx]
        for t in range(8):
            features[:, 8+t, :] = d_hm[:, t*5:(t+1)*5]
        d_pq = wt.pq_type_T2[wt_idx] - ala.pq_type_T2[ala_idx]
        for t in range(8):
            features[:, 16+t, :] = d_pq[:, t*5:(t+1)*5]
        d_disp = wt.disp_type_T2[wt_idx] - ala.disp_type_T2[ala_idx]
        for t in range(8):
            features[:, 24+t, :] = d_disp[:, t*5:(t+1)*5]
        d_mc = wt.mc_category_T2[wt_idx] - ala.mc_category_T2[ala_idx]
        for c in range(5):
            features[:, 32+c, :] = d_mc[:, c*5:(c+1)*5]
        features[:, 37, :] = T2(wt.coulomb[wt_idx] - ala.coulomb[ala_idx])
        features[:, 38, :] = T2(wt.ringchi[wt_idx] - ala.ringchi[ala_idx])
        features[:, 39, :] = T2(wt.hbond[wt_idx] - ala.hbond[ala_idx])
        features[:, 40, :] = T2(wt.mopac_coulomb[wt_idx] - ala.mopac_coulomb[ala_idx])
        d_mopac_mc = wt.mopac_mc_category_T2[wt_idx] - ala.mopac_mc_category_T2[ala_idx]
        for c in range(5):
            features[:, 41+c, :] = d_mopac_mc[:, c*5:(c+1)*5]

        for m in range(5):
            all_X.append(features[:, :, m])
            all_y.append(dft_t2[:, m])
            for i in range(M):
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
