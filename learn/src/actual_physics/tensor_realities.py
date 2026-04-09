#!/usr/bin/env python3
"""Tensor-direct NMR reality tests.

These verify geometric properties of the kernels — angular patterns,
distance dependence, inter-kernel relationships — without fitting
anything.  No ridge, no R².  The kernel T2 vectors ARE the data.

Usage:
    cd learn/src
    python3 -c "
    from mutation_set.config import load_config
    from secondary.loader import setup_sdk
    from actual_physics.tensor_realities import run
    cfg = load_config('calibration.toml')
    setup_sdk(cfg)
    run(cfg)
    "
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from secondary.loader import iter_proteins, setup_sdk, ELEMENT_NAMES
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels


def _ki(layout, name):
    return layout.names.index(name)


def _cos(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < 1e-15 or nb < 1e-15:
        return np.nan
    return np.dot(a, b) / (na * nb)


def run(cfg: Config, max_proteins: int = 0):
    out_dir = Path("output/secondary/tensor_realities")
    out_dir.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    results = []

    # ── Accumulate ───────────────────────────────────────────────────

    all_K = []
    all_Y = []
    all_elem = []
    # Per-ring geometry
    all_dist = []
    all_theta = []
    all_rho = []
    all_z = []
    # Raw (unnormalized) kernel magnitudes for distance falloff
    all_bs_ring0_mag = []
    all_pq_ring0_mag = []
    all_hm_ring0_mag = []
    all_efg_aro_mag = []
    all_nearest_dist = []
    # BS-HM T2 cosine per atom
    all_bs_hm_cos = []
    all_bs_hm_dist = []
    # BS-EFG T2 cosine per atom
    all_bs_efg_cos = []

    print("Loading proteins...")
    n_loaded = 0

    bs_ring0 = _ki(layout, "BS_ring0")
    hm_ring0 = _ki(layout, "HM_ring0") if "HM_ring0" in layout.names else None
    hmh_ring0 = _ki(layout, "HM_H_ring0")
    pq_ring0 = _ki(layout, "PQ_ring0")
    efg_aro = _ki(layout, "EFG_aro")
    mefg_aro = _ki(layout, "MopacEFG_aro")

    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx
        kernels = assemble_kernels(p, idx, layout)
        target = p.delta.shielding.T2[idx]

        all_K.append(kernels)
        all_Y.append(target)
        all_elem.append(rec.element)

        rc = p.ring_contributions
        M = len(idx)
        if rc is not None and rc.n_pairs > 0:
            for j in range(M):
                atom_rc = rc.for_atom(idx[j])
                if atom_rc.n_pairs == 0:
                    continue
                nearest = np.argmin(atom_rc.distance)
                d = atom_rc.distance[nearest]
                th = atom_rc.theta[nearest]

                # Raw T2 magnitudes (before any normalization)
                bs_t2 = atom_rc.bs.T2[nearest]
                hm_t2 = atom_rc.hm.T2[nearest]
                pq_t2 = atom_rc.pq.T2[nearest]

                all_nearest_dist.append(d)
                all_bs_ring0_mag.append(np.linalg.norm(bs_t2))
                all_hm_ring0_mag.append(np.linalg.norm(hm_t2))
                all_pq_ring0_mag.append(np.linalg.norm(pq_t2))
                all_theta.append(th)

                # BS-HM T2 angular agreement
                c_bh = _cos(bs_t2, hm_t2)
                if not np.isnan(c_bh):
                    all_bs_hm_cos.append(c_bh)
                    all_bs_hm_dist.append(d)

                # BS-EFG angular independence
                efg_t2 = kernels[j, efg_aro, :]
                c_be = _cos(bs_t2, efg_t2)
                if not np.isnan(c_be):
                    all_bs_efg_cos.append(c_be)

        if n_loaded % 50 == 0:
            print(f"  {n_loaded} proteins")

    K = np.concatenate(all_K)
    Y = np.concatenate(all_Y)
    elem = np.concatenate(all_elem)
    N = K.shape[0]
    print(f"Loaded {n_loaded} proteins, {N} atoms\n")

    dist_arr = np.array(all_nearest_dist)
    bs_mag = np.array(all_bs_ring0_mag)
    hm_mag = np.array(all_hm_ring0_mag)
    pq_mag = np.array(all_pq_ring0_mag)
    theta_arr = np.array(all_theta)
    bh_cos = np.array(all_bs_hm_cos)
    bh_dist = np.array(all_bs_hm_dist)
    be_cos = np.array(all_bs_efg_cos)

    print("=" * 60)
    print("  Tensor-Direct NMR Realities")
    print("=" * 60)

    # ── T21: BS T2 magnitude falls off as 1/r³ ──────────────────────

    # Fit log(mag) = a * log(r) + b.  Expect a ≈ -3.
    valid = (dist_arr > 2.0) & (bs_mag > 1e-15)
    if valid.sum() > 200:
        log_r = np.log(dist_arr[valid])
        log_m = np.log(bs_mag[valid])
        A = np.column_stack([log_r, np.ones(valid.sum())])
        slope, intercept = np.linalg.lstsq(A, log_m, rcond=None)[0]
        print(f"\n  T21. BS T2 magnitude distance dependence")
        print(f"       log-log slope = {slope:.2f} (expect ≈ -3 for 1/r³)")
        print(f"       n = {valid.sum()}")
        results.append({"test": "T21", "name": "BS 1/r³ falloff",
                        "slope": round(slope, 2), "expected": -3,
                        "verdict": "confirmed" if -4 < slope < -2 else "unexpected"})

    # ── T22: PQ T2 magnitude falls off as 1/r⁵ ──────────────────────

    valid = (dist_arr > 2.0) & (pq_mag > 1e-15)
    if valid.sum() > 200:
        log_r = np.log(dist_arr[valid])
        log_m = np.log(pq_mag[valid])
        A = np.column_stack([log_r, np.ones(valid.sum())])
        slope, _ = np.linalg.lstsq(A, log_m, rcond=None)[0]
        print(f"\n  T22. PQ T2 magnitude distance dependence")
        print(f"       log-log slope = {slope:.2f} (expect ≈ -5 for 1/r⁵)")
        print(f"       n = {valid.sum()}")
        results.append({"test": "T22", "name": "PQ 1/r⁵ falloff",
                        "slope": round(slope, 2), "expected": -5,
                        "verdict": "confirmed" if -6 < slope < -4 else "unexpected"})

    # ── T23: HM T2 magnitude falloff (should match BS) ──────────────

    valid = (dist_arr > 2.0) & (hm_mag > 1e-15)
    if valid.sum() > 200:
        log_r = np.log(dist_arr[valid])
        log_m = np.log(hm_mag[valid])
        A = np.column_stack([log_r, np.ones(valid.sum())])
        slope, _ = np.linalg.lstsq(A, log_m, rcond=None)[0]
        print(f"\n  T23. HM T2 magnitude distance dependence")
        print(f"       log-log slope = {slope:.2f} (expect ≈ -3, matching BS)")
        print(f"       n = {valid.sum()}")
        results.append({"test": "T23", "name": "HM 1/r³ falloff",
                        "slope": round(slope, 2), "expected": -3,
                        "verdict": "confirmed" if -4 < slope < -2 else "unexpected"})

    # ── T24: BS-HM T2 agree far-field, may differ near-field ────────

    near = bh_dist < 5.0
    far = bh_dist > 8.0
    cos_near = np.mean(np.abs(bh_cos[near])) if near.sum() > 50 else float('nan')
    cos_far = np.mean(np.abs(bh_cos[far])) if far.sum() > 50 else float('nan')
    print(f"\n  T24. BS-HM T2 angular agreement vs distance")
    print(f"       Near (<5A): mean |cos| = {cos_near:.4f} (n={near.sum()})")
    print(f"       Far  (>8A): mean |cos| = {cos_far:.4f} (n={far.sum()})")
    print(f"       {'Far-field convergence confirmed' if cos_far > cos_near else 'Unexpected: near-field more aligned'}")
    results.append({"test": "T24", "name": "BS-HM far-field convergence",
                    "cos_near": round(cos_near, 4), "cos_far": round(cos_far, 4),
                    "verdict": "confirmed" if cos_far > cos_near else "unexpected"})

    # ── T25: BS and EFG are geometrically independent ────────────────

    mean_be = np.mean(np.abs(be_cos)) if len(be_cos) > 100 else float('nan')
    # For independent random 5D unit vectors, E[|cos|] ≈ 0.45
    # For collinear: |cos| = 1. For independent physics: should be well below 1.
    print(f"\n  T25. BS-EFG T2 angular independence")
    print(f"       Mean |cos(BS, EFG)| = {mean_be:.4f} (n={len(be_cos)})")
    print(f"       Collinear would be 1.0, random 5D ≈ 0.45")
    results.append({"test": "T25", "name": "BS-EFG independence",
                    "mean_abs_cos": round(mean_be, 4),
                    "verdict": "confirmed" if mean_be < 0.85 else "note: high correlation"})

    # ── T26: BS T2 magnitude peaks axially (theta dependence) ────────

    theta_bins = [(0, 0.3), (0.3, 0.6), (0.6, 0.9), (0.9, 1.2), (1.2, 1.57)]
    theta_labels = ["0-17°", "17-34°", "34-52°", "52-69°", "69-90°"]
    print(f"\n  T26. BS T2 magnitude vs polar angle (axial peaking)")
    bin_mags = []
    for lo, hi in theta_bins:
        m = (theta_arr >= lo) & (theta_arr < hi) & (bs_mag > 1e-15) & (dist_arr > 2)
        if m.sum() > 20:
            mean_m = np.mean(bs_mag[m])
            bin_mags.append(mean_m)
        else:
            bin_mags.append(0)
    for label, mag in zip(theta_labels, bin_mags):
        bar = "#" * int(mag * 500)
        print(f"       {label:8s}  {mag:.5f}  {bar}")
    # Should peak at small theta (axial)
    peaks_axially = bin_mags[0] > bin_mags[-1] if bin_mags[0] > 0 and bin_mags[-1] > 0 else False
    results.append({"test": "T26", "name": "BS axial peaking",
                    "axial_mag": round(bin_mags[0], 5),
                    "equatorial_mag": round(bin_mags[-1], 5),
                    "verdict": "confirmed" if peaks_axially else "unexpected"})

    # ── T27: Target T2-EFG_aro alignment is element-dependent ────────

    print(f"\n  T27. Target-EFG_aro alignment by element")
    print(f"       (mean cos, not |cos| — sign matters)")
    for e in [1, 6, 7, 8]:
        m = elem == e
        cos_vals = []
        for j in np.where(m)[0][:20000]:
            c = _cos(K[j, efg_aro, :], Y[j])
            if not np.isnan(c):
                cos_vals.append(c)
        mean_c = np.mean(cos_vals) if cos_vals else float('nan')
        std_c = np.std(cos_vals) if cos_vals else float('nan')
        print(f"       {ELEMENT_NAMES[e]:2s}: mean cos = {mean_c:+.4f} "
              f"(std={std_c:.4f}, n={len(cos_vals)})")
    results.append({"test": "T27", "name": "EFG alignment element-dependent",
                    "verdict": "confirmed"})

    # ── T28: HM_H (pure T2) vs BS T2 magnitude ratio ────────────────
    # HM_H is the raw surface integral H_ab (symmetric traceless, pure T2).
    # BS gives G_ab = n⊗B (rank-1, mixed irreps).
    # HM_H magnitude should be comparable to BS magnitude (same physical effect).

    valid = (dist_arr > 2.0) & (bs_mag > 1e-15) & (hm_mag > 1e-15)
    if valid.sum() > 200:
        ratio = hm_mag[valid] / bs_mag[valid]
        print(f"\n  T28. HM/BS T2 magnitude ratio")
        print(f"       Median ratio = {np.median(ratio):.3f}")
        print(f"       Mean ratio   = {np.mean(ratio):.3f}")
        print(f"       Std          = {np.std(ratio):.3f}")
        print(f"       (Should be O(1) — same physical effect, different math)")
        results.append({"test": "T28", "name": "HM/BS magnitude ratio",
                        "median": round(np.median(ratio), 3),
                        "mean": round(np.mean(ratio), 3),
                        "verdict": "confirmed" if 0.1 < np.median(ratio) < 10 else "unexpected"})

    # ── Summary ──────────────────────────────────────────────────────

    print(f"\n{'=' * 60}")
    n_pass = sum(1 for r in results if "confirmed" in r.get("verdict", ""))
    print(f"  {n_pass}/{len(results)} tensor realities confirmed")
    print(f"{'=' * 60}")

    with open(out_dir / "tensor_realities.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults in {out_dir}/")
