#!/usr/bin/env python3
"""
Status calculator: where are we?

Run this at the start of any session to see the current state of
extraction, data quality, and training baselines. It exercises the
full pipeline (load → match → delta → Ridge fit → T2 diagnostics)
so any breakage surfaces immediately.

Usage:
    python learn/status.py
"""

import numpy as np
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from load import load_pair, list_proteins, T2, T0, T1, RING_TYPES
from delta import match_atoms


def section(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}\n")


def main():
    t0_time = time.time()

    # ============================================================
    # 1. EXTRACTION STATE
    # ============================================================
    section("1. EXTRACTION STATE")

    features_dir = Path(__file__).parent / "features" / "FirstExtraction"
    proteins = list_proteins(features_dir)
    print(f"Complete WT+ALA pairs: {len(proteins)}")

    # Count total proteins with at least WT
    wt_only = sum(1 for d in features_dir.iterdir()
                  if d.is_dir() and (d / "wt" / "pos.npy").exists())
    print(f"WT extracted (any): {wt_only}")

    # Check array count on first protein
    if proteins:
        sample = features_dir / proteins[0] / "wt"
        n_arrays = len(list(sample.glob("*.npy")))
        print(f"Arrays per conformation: {n_arrays} (expect 46)")
        has_mopac = (sample / "mopac_charges.npy").exists()
        has_orca = (sample / "orca_total.npy").exists()
        print(f"MOPAC present: {has_mopac}  ORCA DFT present: {has_orca}")

    # Check extraction log
    log_path = features_dir / "extract_log.jsonl"
    if log_path.exists():
        import json
        lines = log_path.read_text().strip().splitlines()
        entries = [json.loads(l) for l in lines]
        ok = sum(1 for e in entries if e["ok"])
        fail = sum(1 for e in entries if not e["ok"])
        print(f"Extraction log: {ok} OK, {fail} failed, {len(entries)} total")
        if fail > 0:
            print("  FAILURES:")
            for e in entries:
                if not e["ok"]:
                    print(f"    {e['protein_id']}/{e['variant']}: {e['error'][:80]}")

    if len(proteins) < 5:
        print("\nNot enough data for analysis. Waiting for extraction.")
        return

    # ============================================================
    # 2. DELTA T2 — THE TARGET
    # ============================================================
    section("2. DELTA T2 (WT-ALA mutation delta = ring effect)")

    delta_dft_t2_norms = []
    delta_dft_t0_vals = []
    n_matched_total = 0
    valid_proteins = []

    for pid in proteins:
        wt, ala = load_pair(pid, features_dir)
        if wt is None or ala is None or wt.orca is None or ala.orca is None:
            continue
        wt_idx, ala_idx, _ = match_atoms(wt.pos, wt.element, ala.pos, ala.element)
        if len(wt_idx) < 10:
            continue
        valid_proteins.append(pid)
        n_matched_total += len(wt_idx)

        delta_dft = wt.orca[wt_idx] - ala.orca[ala_idx]
        t2 = T2(delta_dft)
        t0 = T0(delta_dft).ravel()
        delta_dft_t2_norms.extend(np.sqrt((t2**2).sum(axis=1)).tolist())
        delta_dft_t0_vals.extend(np.abs(t0).tolist())

    print(f"Valid proteins (WT+ALA+DFT, >=10 matched atoms): {len(valid_proteins)}")
    print(f"Total matched heavy atoms: {n_matched_total}")
    print(f"Delta DFT |T2| mean: {np.mean(delta_dft_t2_norms):.4f} ppm")
    print(f"Delta DFT |T0| mean: {np.mean(delta_dft_t0_vals):.4f} ppm")
    print(f"  (T2 is what our kernels explain. T0 is harder — needs magnitude parameters.)")

    # ============================================================
    # 3. RIDGE BASELINE (Level A) ON DELTAS
    # ============================================================
    section("3. RIDGE BASELINE (Level A) on delta T2")

    from sklearn.linear_model import Ridge

    mc_cats = ["CO", "CN", "BB_other", "SC_CO", "Aromatic"]
    feature_names = []
    for rt in RING_TYPES:
        feature_names.append(f"BS_{rt}")
    for rt in RING_TYPES:
        feature_names.append(f"HM_{rt}")
    for rt in RING_TYPES:
        feature_names.append(f"PQ_{rt}")
    for rt in RING_TYPES:
        feature_names.append(f"Disp_{rt}")
    for cat in mc_cats:
        feature_names.append(f"MC_{cat}")
    feature_names += ["Coulomb", "RingSusc", "HBond", "MopacCoulomb"]
    for cat in mc_cats:
        feature_names.append(f"MopacMC_{cat}")
    K = len(feature_names)

    all_X = []
    all_y = []

    for pid in valid_proteins:
        wt, ala = load_pair(pid, features_dir)
        wt_idx, ala_idx, _ = match_atoms(wt.pos, wt.element, ala.pos, ala.element)
        if len(wt_idx) < 10:
            continue
        M = len(wt_idx)
        dft_t2 = T2(wt.orca[wt_idx] - ala.orca[ala_idx])

        features = np.zeros((M, K, 5))
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

    X = np.vstack(all_X)
    y = np.concatenate(all_y)
    n_atoms = X.shape[0] // 5

    print(f"Matched atoms: {n_atoms}")
    print(f"Observations: {X.shape[0]} (atoms x 5 T2 components)")
    print(f"Features: {K}")

    model = Ridge(alpha=1.0, fit_intercept=False)
    model.fit(X, y)
    y_pred = model.predict(X)
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot

    # Per-component R²
    for m in range(5):
        idx = np.arange(m * n_atoms, (m+1) * n_atoms)
        ym, yp = y[idx], y_pred[idx]
        r2m = 1.0 - np.sum((ym - yp)**2) / np.sum((ym - ym.mean())**2)
        print(f"  T2 component m={m-2:+d}: R² = {r2m:.4f}")
    print(f"\n  OVERALL DELTA T2 R² = {r2:.4f}  (Ridge alpha=1.0)")
    print(f"  (This is the linear floor. Equivariant model should beat this.)")

    # Top contributors
    stds = X.std(axis=0)
    contributions = np.abs(model.coef_) * stds
    order = np.argsort(contributions)[::-1]
    print(f"\n  Top 10 contributors (|weight| x feature std):")
    for rank, i in enumerate(order[:10]):
        if stds[i] < 1e-10:
            continue
        print(f"    {rank+1:2d}. {feature_names[i]:18s} w={model.coef_[i]:+8.2f}  "
              f"|w|*std={contributions[i]:.4f}")

    # ============================================================
    # 4. T2 INDEPENDENCE (pairwise calculator |cos|)
    # ============================================================
    section("4. T2 INDEPENDENCE (pairwise |cos| on WT, random=0.36)")

    def cos5d(a, b):
        na = np.sqrt((a**2).sum(axis=1))
        nb = np.sqrt((b**2).sum(axis=1))
        mask = (na > 1e-12) & (nb > 1e-12)
        if mask.sum() == 0:
            return 0.0
        dots = (a[mask] * b[mask]).sum(axis=1) / (na[mask] * nb[mask])
        return np.mean(np.abs(dots))

    # Accumulate from a few proteins
    calc_t2s = {name: [] for name in
                ["BS", "HM", "MC", "Coul", "MopCoul", "MopMC", "RSusc", "HBond"]}
    for pid in valid_proteins[:50]:
        wt, _ = load_pair(pid, features_dir)
        if wt is None:
            continue
        calc_t2s["BS"].append(T2(wt.bs))
        calc_t2s["HM"].append(T2(wt.hm))
        calc_t2s["MC"].append(T2(wt.mc))
        calc_t2s["Coul"].append(T2(wt.coulomb))
        calc_t2s["MopCoul"].append(T2(wt.mopac_coulomb))
        calc_t2s["MopMC"].append(T2(wt.mopac_mc))
        calc_t2s["RSusc"].append(T2(wt.ringchi))
        calc_t2s["HBond"].append(T2(wt.hbond))

    for name in calc_t2s:
        calc_t2s[name] = np.vstack(calc_t2s[name])

    names = list(calc_t2s.keys())
    print(f"  {'':12s}", end="")
    for n in names:
        print(f" {n:>7s}", end="")
    print()
    for i, ni in enumerate(names):
        print(f"  {ni:12s}", end="")
        for j, nj in enumerate(names):
            if j <= i:
                c = cos5d(calc_t2s[ni], calc_t2s[nj])
                print(f" {c:7.3f}", end="")
            else:
                print(f"       ", end="")
        print()
    print(f"\n  Key pairs:")
    print(f"    Coulomb vs MopacCoulomb: {cos5d(calc_t2s['Coul'], calc_t2s['MopCoul']):.3f}"
          f"  (charge polarisation effect)")
    print(f"    MC vs MopacMC:           {cos5d(calc_t2s['MC'], calc_t2s['MopMC']):.3f}"
          f"  (bond-order weighting effect)")
    print(f"    BS vs HM:                {cos5d(calc_t2s['BS'], calc_t2s['HM']):.3f}"
          f"  (same physics, different approximation)")

    # ============================================================
    # 5. EQUIVARIANT MODEL STATUS
    # ============================================================
    section("5. EQUIVARIANT MODEL RUNS")

    runs_dir = Path(__file__).parent / "runs"
    if runs_dir.exists():
        import json
        for run_dir in sorted(runs_dir.iterdir()):
            if not run_dir.is_dir():
                continue
            header = run_dir / "header.json"
            summary = run_dir / "summary.json"
            if header.exists():
                h = json.loads(header.read_text())
                print(f"  {run_dir.name}:")
                print(f"    proteins={h.get('n_proteins', '?')}  "
                      f"epochs={h.get('epochs', '?')}  "
                      f"lr={h.get('lr', '?')}  "
                      f"warmup={h.get('warmup_epochs', '?')}")
                if summary.exists():
                    s = json.loads(summary.read_text())
                    print(f"    Train R²={s.get('train_r2', '?')}  "
                          f"Val R²={s.get('val_r2', '?')}  "
                          f"Val RMSE={s.get('val_rmse_ppm', '?')} ppm")
                    print(f"    Best epoch: {s.get('best_epoch', '?')}")
                if h.get("notes"):
                    print(f"    Notes: {h['notes']}")
    else:
        print("  No runs yet.")

    # ============================================================
    # 6. SUMMARY
    # ============================================================
    section("6. SUMMARY")

    print(f"  Proteins extracted:     {len(proteins)} complete pairs")
    print(f"  Valid for training:     {len(valid_proteins)} (WT+ALA+DFT, >=10 matched)")
    print(f"  Ridge delta T2 R²:     {r2:.4f}")
    print(f"  Delta |T2| scale:      {np.mean(delta_dft_t2_norms):.4f} ppm")
    print(f"  Kernels (K):           {K}")
    print(f"  MOPAC present:         yes (calculators 9+10)")

    elapsed = time.time() - t0_time
    print(f"\n  Status computed in {elapsed:.1f}s")

    # Save timestamped snapshot
    from datetime import datetime
    import json
    snapshots_dir = Path(__file__).parent / "snapshots"
    snapshots_dir.mkdir(exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    snapshot = {
        "timestamp": datetime.now().isoformat(),
        "n_proteins_extracted": len(proteins),
        "n_valid_proteins": len(valid_proteins),
        "n_matched_atoms": n_matched_total,
        "delta_dft_t2_mean_ppm": round(np.mean(delta_dft_t2_norms), 4),
        "delta_dft_t0_mean_ppm": round(np.mean(delta_dft_t0_vals), 4),
        "ridge_r2": round(r2, 4),
        "ridge_alpha": 1.0,
        "n_kernels": K,
        "top_contributors": [
            {"name": feature_names[i], "weight": round(float(model.coef_[i]), 4),
             "contribution": round(float(contributions[i]), 4)}
            for i in order[:10] if stds[i] > 1e-10
        ],
        "independence": {
            "Coulomb_vs_MopacCoulomb": round(cos5d(calc_t2s['Coul'], calc_t2s['MopCoul']), 4),
            "MC_vs_MopacMC": round(cos5d(calc_t2s['MC'], calc_t2s['MopMC']), 4),
            "BS_vs_HM": round(cos5d(calc_t2s['BS'], calc_t2s['HM']), 4),
        },
    }
    snap_path = snapshots_dir / f"status_{ts}.json"
    with open(snap_path, "w") as f:
        json.dump(snapshot, f, indent=2)
    print(f"\n  Snapshot saved: {snap_path.name}")

    # Known issues / next steps
    print(f"\n  NOTES:")
    print(f"    - Equivariant model scale fix (2026-04-06): weights/sqrt(n_kernels)")
    print(f"      moved R² from -1.27 to +0.48. Still improving at epoch 500.")
    print(f"    - Round 1 results (K=40, pre-MOPAC) are lost. Cannot compare directly.")
    print(f"    - bones/ contains superseded code (mopac_extract.py, ring_geometry.py).")
    print(f"    - Extraction running ~6h total. Check with:")
    print(f"      wc -l learn/features/FirstExtraction/extract_log.jsonl")


if __name__ == "__main__":
    main()
