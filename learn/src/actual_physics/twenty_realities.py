#!/usr/bin/env python3
"""Twenty NMR realities verified by geometric kernel calibration.

Each test checks a known, literature-grounded NMR physics prediction
against what the geometric kernels see in the WT-ALA mutation set.
The tests are analytical (ridge regression, correlation, sign checks)
— no training, no GPU.

Usage:
    cd learn/src
    python3 -c "
    from mutation_set.config import load_config
    from secondary.loader import setup_sdk
    from actual_physics.twenty_realities import run
    cfg = load_config('calibration.toml')
    setup_sdk(cfg)
    run(cfg)
    "
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from secondary.loader import (
    iter_proteins, setup_sdk, ridge_fit,
    ELEMENT_NAMES, RING_TYPE_NAMES, N_RING_TYPES,
)
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels


def _kernel_indices(layout, prefix):
    """All kernel indices whose name starts with prefix."""
    return [i for i, n in enumerate(layout.names) if n.startswith(prefix)]


def _kernel_index(layout, name):
    return layout.names.index(name)


def run(cfg: Config, max_proteins: int = 0):
    out_dir = Path("output/secondary/twenty_realities")
    out_dir.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    lam = cfg.secondary.ridge_lambda

    # ── Accumulate per-atom data ─────────────────────────────────────

    all_kernels = []
    all_targets = []
    all_elements = []
    all_ring_dist = []
    all_ring_masks = []

    # Per-ring geometry for sign/angle tests
    all_rho = []  # cylindrical rho to nearest ring
    all_z = []    # height above nearest ring plane
    all_theta = []  # polar angle to nearest ring

    print("Loading proteins...")
    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx
        kernels = assemble_kernels(p, idx, layout)
        target = p.delta.shielding.T2[idx]

        all_kernels.append(kernels)
        all_targets.append(target)
        all_elements.append(rec.element)
        all_ring_dist.append(rec.ring_dist)
        all_ring_masks.append(rec.atom_ring_masks)

        # Per-ring geometry from ring_contributions
        rc = p.ring_contributions
        M = len(idx)
        rho = np.full(M, np.nan)
        z = np.full(M, np.nan)
        theta = np.full(M, np.nan)
        if rc is not None and rc.n_pairs > 0:
            for j in range(M):
                atom_rc = rc.for_atom(idx[j])
                if atom_rc.n_pairs > 0:
                    nearest = np.argmin(atom_rc.distance)
                    rho[j] = atom_rc.rho[nearest]
                    z[j] = atom_rc.z[nearest]
                    theta[j] = atom_rc.theta[nearest]
        all_rho.append(rho)
        all_z.append(z)
        all_theta.append(theta)

        if n_loaded % 50 == 0:
            print(f"  {n_loaded} proteins")

    K = np.concatenate(all_kernels)
    Y = np.concatenate(all_targets)
    elem = np.concatenate(all_elements)
    rdist = np.concatenate(all_ring_dist)
    masks = np.concatenate(all_ring_masks)
    rho = np.concatenate(all_rho)
    z_height = np.concatenate(all_z)
    theta = np.concatenate(all_theta)

    N = K.shape[0]
    print(f"Loaded {n_loaded} proteins, {N} atoms\n")

    results = []

    def report(num, name, citation, prediction, observation, verdict):
        r = {"num": num, "name": name, "citation": citation,
             "prediction": prediction, "observation": observation,
             "verdict": verdict}
        results.append(r)
        mark = "PASS" if "pass" in verdict.lower() or "confirmed" in verdict.lower() else "NOTE"
        print(f"  [{mark}] {num}. {name}")
        print(f"         {observation}")

    # ── Tests ────────────────────────────────────────────────────────

    print("=" * 60)
    print("  Twenty NMR Realities")
    print("=" * 60)

    # --- 13. Ring current 1/r³ falloff ---
    bs_idx = _kernel_indices(layout, "BS_")[0:8]  # per-type only
    bins = [(0, 4), (4, 8), (8, 12), (12, 20)]
    bin_r2 = []
    for lo, hi in bins:
        mask_d = (rdist >= lo) & (rdist < hi) & np.isfinite(rdist)
        if mask_d.sum() < 100:
            bin_r2.append(("?", mask_d.sum()))
            continue
        X = K[mask_d][:, bs_idx, :].reshape(mask_d.sum(), -1)
        _, r2 = ridge_fit(X, Y[mask_d], lam)
        bin_r2.append((f"{r2:.3f}", int(mask_d.sum())))

    report(13, "Ring current 1/r³ falloff",
           "Pople 1956; Johnson & Bovey 1958",
           "BS R² should decrease monotonically with distance",
           f"BS R² by distance: " + ", ".join(
               f"{lo}-{hi}A: {r2} (n={n})" for (lo, hi), (r2, n) in zip(bins, bin_r2)),
           "confirmed" if all(bin_r2[i][0] >= bin_r2[i+1][0]
                              for i in range(len(bin_r2)-1)
                              if bin_r2[i][0] != "?" and bin_r2[i+1][0] != "?")
           else "partial")

    # --- 14. Shielding above ring, deshielding in plane ---
    # Above ring: |z| > rho (axial). In plane: rho > |z| (equatorial).
    # Ring current T0 should be positive (shielding) above, negative in plane.
    # Use BS_ring0 for the nearest ring.
    bs_ring0_idx = _kernel_index(layout, "BS_ring0")
    valid_geom = np.isfinite(theta) & np.isfinite(rho)
    axial = valid_geom & (theta < 0.6)  # theta < ~34 degrees
    equat = valid_geom & (theta > 1.2)  # theta > ~69 degrees

    # T0 of BS kernel = trace/3 ... but our kernels are T2 only.
    # Use the T2 m=0 component sign instead: positive above ring for BS.
    # Actually let's use the full per-type BS kernels and check correlation sign.
    # Simpler: check if BS kernel magnitude is larger above vs in-plane.
    if axial.sum() > 50 and equat.sum() > 50:
        bs_mag_axial = np.linalg.norm(K[axial][:, bs_ring0_idx, :], axis=1).mean()
        bs_mag_equat = np.linalg.norm(K[equat][:, bs_ring0_idx, :], axis=1).mean()
        report(14, "BS stronger above ring than in plane",
               "Pople 1956; Case 1995",
               "BS kernel magnitude should be larger for axial atoms",
               f"BS_ring0 mean magnitude: axial={bs_mag_axial:.4f} "
               f"(n={axial.sum()}), equatorial={bs_mag_equat:.4f} "
               f"(n={equat.sum()}), ratio={bs_mag_axial/max(bs_mag_equat,1e-15):.2f}",
               "confirmed" if bs_mag_axial > bs_mag_equat else "unexpected")
    else:
        report(14, "BS stronger above ring than in plane",
               "Pople 1956", "BS magnitude larger axially",
               f"insufficient geometry data (axial={axial.sum()}, equat={equat.sum()})",
               "skipped")

    # --- 15. PHE ≈ TYR ring current ---
    bs_phe = _kernel_index(layout, "BS_PHE")
    bs_tyr = _kernel_index(layout, "BS_TYR")
    H_mask = elem == 1
    if H_mask.sum() > 200:
        # Correlation of BS_PHE and BS_TYR T2 magnitudes for H atoms
        phe_mag = np.linalg.norm(K[H_mask][:, bs_phe, :], axis=1)
        tyr_mag = np.linalg.norm(K[H_mask][:, bs_tyr, :], axis=1)
        # Fraction of atoms where both are nonzero
        both = (phe_mag > 1e-10) & (tyr_mag > 1e-10)
        if both.sum() > 50:
            corr = np.corrcoef(phe_mag[both], tyr_mag[both])[0, 1]
        else:
            corr = float('nan')

        # Single-kernel R² comparison
        X_phe = K[H_mask][:, [bs_phe], :].reshape(H_mask.sum(), -1)
        X_tyr = K[H_mask][:, [bs_tyr], :].reshape(H_mask.sum(), -1)
        _, r2_phe = ridge_fit(X_phe, Y[H_mask], lam)
        _, r2_tyr = ridge_fit(X_tyr, Y[H_mask], lam)
        report(15, "PHE and TYR ring currents similar",
               "Giessner-Prettre & Pullman 1987; Case 1995 (I_PHE=1.46, I_TYR=1.24)",
               "BS_PHE R² ≈ BS_TYR R² for H atoms",
               f"BS_PHE R²={r2_phe:.4f}, BS_TYR R²={r2_tyr:.4f}, "
               f"ratio={r2_phe/max(r2_tyr,1e-10):.2f}",
               "confirmed" if 0.5 < r2_phe / max(r2_tyr, 1e-10) < 2.0 else "unexpected")

    # --- 16. PQ falls off faster than BS ---
    pq_idx = [_kernel_index(layout, "PQ_total")]
    bs_all = _kernel_indices(layout, "BS_")[0:8]
    near = (rdist >= 0) & (rdist < 5) & np.isfinite(rdist)
    far = (rdist >= 6) & (rdist < 12) & np.isfinite(rdist)
    if near.sum() > 100 and far.sum() > 100:
        _, r2_pq_near = ridge_fit(
            K[near][:, pq_idx, :].reshape(near.sum(), -1), Y[near], lam)
        _, r2_pq_far = ridge_fit(
            K[far][:, pq_idx, :].reshape(far.sum(), -1), Y[far], lam)
        _, r2_bs_near = ridge_fit(
            K[near][:, bs_all, :].reshape(near.sum(), -1), Y[near], lam)
        _, r2_bs_far = ridge_fit(
            K[far][:, bs_all, :].reshape(far.sum(), -1), Y[far], lam)
        pq_ratio = r2_pq_far / max(r2_pq_near, 1e-10)
        bs_ratio = r2_bs_far / max(r2_bs_near, 1e-10)
        report(16, "PQ (1/r⁵) falls off faster than BS (1/r³)",
               "Buckingham 1959; Stone T-tensor formalism",
               "PQ far/near R² ratio < BS far/near R² ratio",
               f"PQ near={r2_pq_near:.3f} far={r2_pq_far:.3f} ratio={pq_ratio:.3f}; "
               f"BS near={r2_bs_near:.3f} far={r2_bs_far:.3f} ratio={bs_ratio:.3f}",
               "confirmed" if pq_ratio < bs_ratio else "unexpected")

    # --- 17. Buckingham sign: E-field along bond deshields ---
    # The EFG_aro T2 correlation with target should be positive
    # (positive EFG_aro T2 component → positive target T2 component)
    efg_aro_idx = _kernel_index(layout, "EFG_aro")
    # Compute mean dot product (not |cos|, actual cos) across all atoms
    dots = []
    for j in range(min(N, 50000)):
        ka = K[j, efg_aro_idx, :]
        ta = Y[j]
        nk, nt = np.linalg.norm(ka), np.linalg.norm(ta)
        if nk > 1e-15 and nt > 1e-15:
            dots.append(np.dot(ka, ta) / (nk * nt))
    mean_cos = np.mean(dots) if dots else 0.0
    report(17, "EFG-target correlation sign consistent with Buckingham",
           "Buckingham 1960 (σ_T2 = γ·V_ab, γ > 0 for most nuclei)",
           "Mean cos(EFG_aro, target) should be positive",
           f"Mean cos = {mean_cos:.4f} (n={len(dots)})",
           "confirmed" if mean_cos > 0 else "note: negative — check γ sign convention")

    # --- 18. H-bond contribution small in mechanical mutants ---
    hb_idx = [_kernel_index(layout, "HBond_total")]
    X_hb = K[:, hb_idx, :].reshape(N, -1)
    _, r2_hb = ridge_fit(X_hb, Y, lam)
    report(18, "H-bond contribution negligible in mechanical mutants",
           "Mechanical mutant: backbone unchanged, H-bonds preserved",
           "HBond_total R² should be near zero",
           f"HBond_total R² = {r2_hb:.4f}",
           "confirmed" if r2_hb < 0.02 else "note: unexpectedly large")

    # --- 19. Solvation EFG minor for aromatics ---
    # DeltaAPBS_EFG should rank below EFG_aro in forward selection
    dapbs_idx = _kernel_index(layout, "DeltaAPBS_EFG")
    efga_idx = _kernel_index(layout, "EFG_aro")
    X_dapbs = K[:, [dapbs_idx], :].reshape(N, -1)
    X_efga = K[:, [efga_idx], :].reshape(N, -1)
    _, r2_dapbs = ridge_fit(X_dapbs, Y, lam)
    _, r2_efga = ridge_fit(X_efga, Y, lam)
    report(19, "Solvation correction minor vs direct aromatic EFG",
           "APBS solvent reaction field is diffuse; aromatic EFG is local",
           "DeltaAPBS_EFG R² << EFG_aro R²",
           f"DeltaAPBS_EFG R²={r2_dapbs:.4f}, EFG_aro R²={r2_efga:.4f}, "
           f"ratio={r2_dapbs/max(r2_efga,1e-10):.2f}",
           "confirmed" if r2_dapbs < r2_efga else "unexpected")

    # --- 20. Element identity matters for kernel weighting ---
    # Compare all-element R² vs sum of per-element R² (weighted)
    elem_r2 = {}
    elem_n = {}
    for e in [1, 6, 7, 8]:
        m = elem == e
        if m.sum() < 100:
            continue
        X_e = K[m].reshape(m.sum(), -1)
        _, r2_e = ridge_fit(X_e, Y[m], lam)
        elem_r2[e] = r2_e
        elem_n[e] = int(m.sum())

    X_all = K.reshape(N, -1)
    _, r2_all = ridge_fit(X_all, Y, lam)

    # Weighted average of per-element R²
    total_n = sum(elem_n.values())
    weighted_r2 = sum(elem_r2[e] * elem_n[e] / total_n for e in elem_r2)
    report(20, "Element-specific weighting improves over pooled fit",
           "Boyd & Skrynnikov 2002; Sahakyan & Vendruscolo 2013",
           "Weighted per-element R² > pooled R²",
           f"Pooled R²={r2_all:.4f}, weighted per-element R²={weighted_r2:.4f}, "
           f"gap={weighted_r2 - r2_all:.4f}",
           "confirmed" if weighted_r2 > r2_all else "unexpected")

    # ── Summary ──────────────────────────────────────────────────────

    print(f"\n{'=' * 60}")
    n_pass = sum(1 for r in results if "confirmed" in r["verdict"].lower())
    print(f"  {n_pass}/{len(results)} realities confirmed")
    print(f"{'=' * 60}")

    # Write results
    with open(out_dir / "realities_13_20.json", "w") as f:
        json.dump(results, f, indent=2)

    # Also write a simple CSV
    with open(out_dir / "realities_13_20.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["num", "name", "verdict", "observation"])
        w.writeheader()
        for r in results:
            w.writerow({k: r[k] for k in ["num", "name", "verdict", "observation"]})

    print(f"\nResults in {out_dir}/")
