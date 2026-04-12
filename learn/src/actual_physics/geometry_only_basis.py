#!/usr/bin/env python3
"""Do the 3 effective kernel dimensions survive in geometry-only mode?

Full analysis: per-element, per-physics-group, with fair normalization,
forward selection, and the critical ff14SB-vs-MOPAC EFG comparison.

In geometry-only mode (no MOPAC/APBS), 44 of 55 core kernels survive.
The 11 dropped: MopacMC (5 bond categories), MopacCoulomb_total,
MopacMC_total, MopacEFG_bb, MopacEFG_aro, APBS_EFG, DeltaAPBS_EFG.

Usage:
    cd learn/src
    python3 -c "
    from mutation_set.config import load_config
    import sys; cfg = load_config('calibration.toml')
    sys.path.insert(0, str(cfg.paths.sdk))
    from actual_physics.geometry_only_basis import run
    run(cfg)
    "
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from secondary.loader import (
    iter_proteins, setup_sdk, ELEMENT_NAMES, RING_TYPE_NAMES,
    ridge_fit,
)
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels, normalize_kernels


# ── Classify kernels ─────────────────────────────────────────────────

_MOPAC_APBS_PREFIXES = (
    "MopacMC_", "MopacCoulomb_total", "MopacMC_total",
    "MopacEFG_", "APBS_", "DeltaAPBS_",
)

def _is_mopac_apbs(name: str) -> bool:
    for prefix in _MOPAC_APBS_PREFIXES:
        if name.startswith(prefix) or name == prefix:
            return True
    return False


def _partition_core(layout: KernelLayout):
    """Return (geo_indices, mopac_indices) into the core kernel range."""
    geo, mopac = [], []
    for i in range(layout.efg_end):
        (mopac if _is_mopac_apbs(layout.names[i]) else geo).append(i)
    return geo, mopac


def _classify_geo_kernels(layout, geo_idx):
    """Physics groups within geometry-only kernels.
    Returns dict mapping group name → list of positions within geo_idx."""
    groups = {
        "ring_current": [],  # BS, HM per type + RingSusc_total
        "bond_aniso": [],    # MC categories + MC_total + HBond_total
        "ff14sb_efg": [],    # Coulomb_total, EFG_bb, EFG_aro
        "quadrupole": [],    # PQ per type + PQ_total
        "dispersion": [],    # Disp per type
    }
    for pos, ki in enumerate(geo_idx):
        name = layout.names[ki]
        if name.startswith("BS_") or name.startswith("HM_") or name == "RingSusc_total":
            groups["ring_current"].append(pos)
        elif name.startswith("MC_") or name == "HBond_total":
            groups["bond_aniso"].append(pos)
        elif name.startswith("PQ_") or name == "PQ_total":
            groups["quadrupole"].append(pos)
        elif name.startswith("Disp_"):
            groups["dispersion"].append(pos)
        elif name in ("Coulomb_total", "EFG_bb", "EFG_aro"):
            groups["ff14sb_efg"].append(pos)
        else:
            print(f"  WARNING: unclassified geo kernel: {name}")
    return groups


def _ridge(X, y, lam=1.0):
    XtX = X.T @ X + lam * np.eye(X.shape[1])
    w = np.linalg.solve(XtX, X.T @ y)
    pred = X @ w
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean(axis=0)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
    return w, r2


def _forward_selection(K_arr, Y_arr, names, lam, max_steps=20):
    """Greedy forward selection on (N, n_kernels, 5) kernel array.
    Returns list of (name, marginal_delta, cumulative_R²)."""
    n_kernels = K_arr.shape[1]
    selected = []
    available = set(range(n_kernels))
    results = []
    prev_r2 = 0.0
    for _ in range(min(max_steps, n_kernels)):
        best_k, best_r2 = -1, -1e30
        for k in available:
            trial = selected + [k]
            Xt = K_arr[:, trial, :].reshape(K_arr.shape[0], -1)
            _, r2 = _ridge(Xt, Y_arr, lam)
            if r2 > best_r2:
                best_k, best_r2 = k, r2
        if best_k < 0:
            break
        selected.append(best_k)
        available.discard(best_k)
        results.append((names[best_k], best_r2 - prev_r2, best_r2))
        prev_r2 = best_r2
    return results


def _effective_dims(cov, thresholds=(0.90, 0.95, 0.98)):
    eigvals = np.linalg.eigvalsh(cov)[::-1]
    eigvals = np.maximum(eigvals, 0)
    total = eigvals.sum()
    if total < 1e-15:
        return {t: 0 for t in thresholds}, eigvals
    cumulative = np.cumsum(eigvals) / total
    dims = {}
    for t in thresholds:
        dims[t] = int(np.searchsorted(cumulative, t)) + 1
    return dims, eigvals


def run(cfg: Config, max_proteins: int = 0):
    out = Path("output/actual_physics/geometry_only_basis")
    out.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    geo_idx, mopac_idx = _partition_core(layout)
    geo_names = [layout.names[i] for i in geo_idx]
    mopac_names = [layout.names[i] for i in mopac_idx]
    n_geo = len(geo_idx)
    n_core = layout.efg_end
    lam = getattr(cfg.secondary, 'ridge_lambda', 1e-2)

    # Physics groups within geo-only
    geo_groups = _classify_geo_kernels(layout, geo_idx)

    # Also identify the specific EFG kernels for the head-to-head
    efg_ff14sb_pos = [pos for pos, ki in enumerate(geo_idx)
                      if layout.names[ki] in ("EFG_bb", "EFG_aro")]
    efg_mopac_core = [i for i in range(n_core)
                      if layout.names[i] in ("MopacEFG_bb", "MopacEFG_aro")]

    print(f"Core kernels: {n_core} (geo-only: {n_geo}, MOPAC/APBS: {len(mopac_idx)})")
    print(f"\nPhysics groups in geo-only ({n_geo} kernels):")
    for gname, gpos in geo_groups.items():
        print(f"  {gname:15s}: {len(gpos):2d} kernels")

    ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O"}

    # Per-element accumulators
    elem_data = {e: {
        "K_full": [], "K_geo": [], "Y": [],
        "scales_full": [], "scales_geo": [],
        "mut_type": [],
        # For EFG head-to-head
        "efg_ff14sb": [], "efg_mopac": [],
    } for e in ELEMENTS}

    print("\nLoading proteins...")
    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx
        M = len(idx)

        # Full 55-kernel assembly
        kernels_full = assemble_kernels(p, idx, layout)[:, :n_core, :]

        # Per-protein normalization: full and geo-only independently
        kf = kernels_full.copy()
        kf_norm, scales_full = normalize_kernels(kf, cfg.normalization.kernel_std_floor)

        kg = kernels_full[:, geo_idx, :].copy()
        kg_norm, scales_geo = normalize_kernels(kg, cfg.normalization.kernel_std_floor)

        target = p.delta.shielding.T2[idx]

        # Mutation type
        rg = p.ring_geometry
        if rg is not None and rg.n_rings > 0:
            type_counts = np.bincount(rg.ring_type.astype(int),
                                       minlength=len(RING_TYPE_NAMES))
            mut_type = int(np.argmax(type_counts))
        else:
            mut_type = -1

        # EFG head-to-head: raw (unnormalized) for isolation
        efg_ff = kernels_full[:, [geo_idx[p] for p in efg_ff14sb_pos], :]
        efg_mop = kernels_full[:, efg_mopac_core, :] if efg_mopac_core else None

        for j in range(M):
            e = int(rec.element[j])
            if e not in elem_data:
                continue
            elem_data[e]["K_full"].append(kf_norm[j])
            elem_data[e]["K_geo"].append(kg_norm[j])
            elem_data[e]["Y"].append(target[j])
            elem_data[e]["scales_full"].append(scales_full)
            elem_data[e]["scales_geo"].append(scales_geo)
            elem_data[e]["mut_type"].append(mut_type)
            elem_data[e]["efg_ff14sb"].append(efg_ff[j])
            if efg_mop is not None:
                elem_data[e]["efg_mopac"].append(efg_mop[j])

        if n_loaded % 50 == 0:
            print(f"  {n_loaded} proteins")

    print(f"Loaded {n_loaded} proteins\n")

    # ══════════════════════════════════════════════════════════════════
    results = {}

    for e, ename in ELEMENTS.items():
        d = elem_data[e]
        if len(d["K_geo"]) < 200:
            continue

        K_full = np.array(d["K_full"])      # (N, 55, 5) normalized
        K_geo = np.array(d["K_geo"])         # (N, 44, 5) normalized
        Y = np.array(d["Y"])                 # (N, 5)
        scales_full = np.array(d["scales_full"])  # (N, 55)
        scales_geo = np.array(d["scales_geo"])    # (N, 44)
        mut_types = np.array(d["mut_type"])
        N = K_geo.shape[0]

        print("=" * 70)
        print(f"  {ename} ({N:,d} atoms)")
        print("=" * 70)

        er = {"n_atoms": N}

        # ── 1. Progressive calibration: geo-only vs full ─────────

        print(f"\n  Progressive calibration:")
        print(f"  {'level':<40s}  {'full':>8s}  {'geo-only':>8s}  {'gap':>8s}")

        # Level 0: normalized kernels only
        X_f0 = K_full.reshape(N, -1)
        X_g0 = K_geo.reshape(N, -1)
        _, r2_f0 = _ridge(X_f0, Y, lam)
        _, r2_g0 = _ridge(X_g0, Y, lam)
        print(f"  {'Normalized kernels':<40s}  {r2_f0:8.4f}  {r2_g0:8.4f}  {r2_f0-r2_g0:+8.4f}")

        # Level 1: + kernel scale factors as interaction
        X_f1_extra = np.zeros((N, n_core * 5))
        for ki in range(n_core):
            X_f1_extra[:, ki*5:(ki+1)*5] = K_full[:, ki, :] * scales_full[:, ki:ki+1]
        X_f1 = np.hstack([X_f0, X_f1_extra])

        X_g1_extra = np.zeros((N, n_geo * 5))
        for ki in range(n_geo):
            X_g1_extra[:, ki*5:(ki+1)*5] = K_geo[:, ki, :] * scales_geo[:, ki:ki+1]
        X_g1 = np.hstack([X_g0, X_g1_extra])

        _, r2_f1 = _ridge(X_f1, Y, lam)
        _, r2_g1 = _ridge(X_g1, Y, lam)
        print(f"  {'+ kernel scales':<40s}  {r2_f1:8.4f}  {r2_g1:8.4f}  {r2_f1-r2_g1:+8.4f}")

        # Level 2: + mutation type as one-hot interaction
        unique_types = sorted(set(t for t in mut_types if t >= 0))
        n_types = len(unique_types)
        type_map = {t: i for i, t in enumerate(unique_types)}

        X_f2_mut = np.zeros((N, n_types * n_core * 5))
        X_g2_mut = np.zeros((N, n_types * n_geo * 5))
        for j in range(N):
            mt = mut_types[j]
            if mt < 0:
                continue
            ti = type_map.get(mt, -1)
            if ti < 0:
                continue
            X_f2_mut[j, ti*n_core*5:(ti+1)*n_core*5] = X_f0[j]
            X_g2_mut[j, ti*n_geo*5:(ti+1)*n_geo*5] = X_g0[j]

        X_f2 = np.hstack([X_f0, X_f2_mut])
        X_g2 = np.hstack([X_g0, X_g2_mut])
        _, r2_f2 = _ridge(X_f2, Y, lam)
        _, r2_g2 = _ridge(X_g2, Y, lam)
        print(f"  {'+ mutation type':<40s}  {r2_f2:8.4f}  {r2_g2:8.4f}  {r2_f2-r2_g2:+8.4f}")

        # Fair set: normalized + scales + mutation type
        X_f_fair = np.hstack([X_f0, X_f1_extra, X_f2_mut])
        X_g_fair = np.hstack([X_g0, X_g1_extra, X_g2_mut])
        _, r2_f_fair = _ridge(X_f_fair, Y, lam)
        _, r2_g_fair = _ridge(X_g_fair, Y, lam)
        print(f"  {'Fair set (norm+scales+mut)':<40s}  {r2_f_fair:8.4f}  {r2_g_fair:8.4f}  {r2_f_fair-r2_g_fair:+8.4f}")

        er["progressive"] = {
            "normalized": {"full": round(r2_f0, 4), "geo": round(r2_g0, 4)},
            "plus_scales": {"full": round(r2_f1, 4), "geo": round(r2_g1, 4)},
            "plus_mut": {"full": round(r2_f2, 4), "geo": round(r2_g2, 4)},
            "fair": {"full": round(r2_f_fair, 4), "geo": round(r2_g_fair, 4)},
        }

        # ── 2. Physics-group decomposition (geo-only, normalized) ─

        print(f"\n  Physics groups (geo-only, normalized):")
        group_r2 = {}
        for gname, gpos in geo_groups.items():
            if not gpos:
                continue
            X_g = K_geo[:, gpos, :].reshape(N, -1)
            _, r2_g = _ridge(X_g, Y, lam)
            gnames = [geo_names[p] for p in gpos]
            print(f"    {gname:15s} ({len(gpos):2d} kernels): R² = {r2_g:.4f}")
            group_r2[gname] = round(r2_g, 4)
        er["group_r2"] = group_r2

        # ── 3. EFG head-to-head: ff14SB vs MOPAC (raw, isolated) ──

        efg_ff = np.array(d["efg_ff14sb"])   # (N, 2, 5)
        X_ff = efg_ff.reshape(N, -1)
        _, r2_ff = _ridge(X_ff, Y, lam)

        if d["efg_mopac"]:
            efg_mop = np.array(d["efg_mopac"])  # (N, 2, 5)
            X_mop = efg_mop.reshape(N, -1)
            _, r2_mop = _ridge(X_mop, Y, lam)

            # Combined: ff14SB + MOPAC EFG together
            X_both = np.hstack([X_ff, X_mop])
            _, r2_both = _ridge(X_both, Y, lam)

            # Cosine similarity between ff14SB and MOPAC EFG per atom
            cos_vals = []
            for j in range(min(N, 30000)):
                ff_vec = efg_ff[j].ravel()
                mop_vec = efg_mop[j].ravel()
                nf, nm = np.linalg.norm(ff_vec), np.linalg.norm(mop_vec)
                if nf > 1e-15 and nm > 1e-15:
                    cos_vals.append(np.dot(ff_vec, mop_vec) / (nf * nm))
            mean_cos = np.mean(cos_vals) if cos_vals else float('nan')

            print(f"\n  EFG head-to-head (raw, bb+aro only):")
            print(f"    ff14SB Coulomb EFG:   R² = {r2_ff:.4f}")
            print(f"    MOPAC Mulliken EFG:   R² = {r2_mop:.4f}")
            print(f"    Both combined:        R² = {r2_both:.4f}")
            print(f"    ff14SB-MOPAC cos:     {mean_cos:.4f} (1.0 = identical direction)")
            print(f"    MOPAC advantage:      {r2_mop - r2_ff:+.4f}")

            er["efg_head_to_head"] = {
                "ff14sb": round(r2_ff, 4),
                "mopac": round(r2_mop, 4),
                "combined": round(r2_both, 4),
                "cosine": round(mean_cos, 4),
            }
        else:
            print(f"\n  EFG (ff14SB only, no MOPAC data):")
            print(f"    ff14SB Coulomb EFG:   R² = {r2_ff:.4f}")
            er["efg_head_to_head"] = {"ff14sb": round(r2_ff, 4)}

        # ── 4. Eigenspectrum (normalized geo-only) ────────────────

        X_geo_flat = K_geo.reshape(N, -1)
        cov_geo = np.cov(X_geo_flat, rowvar=False)
        dims_geo, eigvals_geo = _effective_dims(cov_geo)
        total_var = eigvals_geo.sum()

        print(f"\n  Eigenspectrum (normalized geo-only):")
        print(f"    Effective dims: 90%={dims_geo[0.90]}, "
              f"95%={dims_geo[0.95]}, 98%={dims_geo[0.98]}")

        # Check for the 5-component block structure
        # If physics sources are distinct, eigenvalues come in blocks of 5
        print(f"    Top 15 eigenvalues (fraction):")
        fracs = eigvals_geo / total_var if total_var > 0 else eigvals_geo
        for k in range(min(15, len(fracs))):
            bar = "#" * int(fracs[k] * 200)
            block = k // 5 + 1
            print(f"      {k+1:2d} (block {block}): {fracs[k]:.4f}  {bar}")

        # Block variance: sum eigenvalues 1-5, 6-10, 11-15
        block_var = []
        for b in range(3):
            bv = sum(fracs[b*5 + i] for i in range(5) if b*5 + i < len(fracs))
            block_var.append(bv)
        print(f"    Block sums: source1={block_var[0]:.3f}, "
              f"source2={block_var[1]:.3f}, source3={block_var[2]:.3f}")

        er["eigenspectrum"] = {
            "dims": {str(k): v for k, v in dims_geo.items()},
            "top15_frac": [round(float(fracs[k]), 4) for k in range(min(15, len(fracs)))],
            "block_sums": [round(bv, 4) for bv in block_var],
        }

        # ── 5. Forward selection (geo-only, normalized) ───────────

        print(f"\n  Forward selection (geo-only, normalized, top 15):")
        fwd = _forward_selection(K_geo, Y, geo_names, lam, max_steps=15)
        fwd_data = []
        for rank, (kn, delta, cum) in enumerate(fwd, 1):
            print(f"    {rank:2d}. {kn:25s}  +{delta:.4f}  R²={cum:.4f}")
            fwd_data.append({"kernel": kn, "delta": round(delta, 4),
                             "cumulative": round(cum, 4)})
        er["forward_selection"] = fwd_data

        # How many geo-only kernels to reach 90% of the full fair R²?
        if r2_f_fair > 0:
            target_90 = r2_f_fair * 0.90
            n_for_90 = len(fwd)
            for rank, (_, _, cum) in enumerate(fwd, 1):
                if cum >= target_90:
                    n_for_90 = rank
                    break
            print(f"    Kernels for 90% of full-fair R² ({target_90:.3f}): {n_for_90}")
            er["n_for_90pct_of_full"] = n_for_90

        # ── 6. What MOPAC adds beyond geo-only (residual analysis) ─

        # Fit geo-only, compute residual, then fit MOPAC kernels to residual
        w_geo, r2_geo_base = _ridge(X_g0, Y, lam)
        pred_geo = X_g0 @ w_geo
        residual = Y - pred_geo

        if d["efg_mopac"]:
            # Can MOPAC EFG explain the geo-only residual?
            X_mop_flat = np.array(d["efg_mopac"]).reshape(N, -1)
            _, r2_mop_on_resid = _ridge(X_mop_flat, residual, lam)

            # Full MOPAC set (all 11 kernels) on residual
            K_mopac_only = np.array(d["K_full"])[:, mopac_idx, :]
            X_mop_all = K_mopac_only.reshape(N, -1)
            _, r2_all_mop_on_resid = _ridge(X_mop_all, residual, lam)

            resid_mag = np.sqrt(np.mean(residual ** 2))
            target_mag = np.sqrt(np.mean(Y ** 2))

            print(f"\n  Residual analysis (what MOPAC adds):")
            print(f"    Geo-only residual RMS:         {resid_mag:.4f} "
                  f"({resid_mag/target_mag*100:.1f}% of target)")
            print(f"    MOPAC EFG R² on residual:      {r2_mop_on_resid:.4f}")
            print(f"    All 11 MOPAC R² on residual:   {r2_all_mop_on_resid:.4f}")

            er["residual"] = {
                "resid_rms": round(resid_mag, 4),
                "target_rms": round(target_mag, 4),
                "mopac_efg_on_resid": round(r2_mop_on_resid, 4),
                "all_mopac_on_resid": round(r2_all_mop_on_resid, 4),
            }

        results[ename] = er

    # ── Summary ─────────────────────────────────────────────────────

    print(f"\n{'=' * 70}")
    print(f"  Summary: Fair-Set R² (norm + scales + mutation type)")
    print(f"{'=' * 70}")
    print(f"  {'elem':>5s}  {'full':>8s}  {'geo-only':>8s}  {'gap':>8s}  "
          f"{'ff14sb_efg':>10s}  {'mopac_efg':>10s}")
    for ename in ["H", "C", "N", "O"]:
        if ename not in results:
            continue
        r = results[ename]
        p = r["progressive"]
        efg = r.get("efg_head_to_head", {})
        print(f"  {ename:>5s}  {p['fair']['full']:>8.4f}  {p['fair']['geo']:>8.4f}  "
              f"{p['fair']['full']-p['fair']['geo']:>+8.4f}  "
              f"{efg.get('ff14sb', 0):>10.4f}  "
              f"{efg.get('mopac', 0):>10.4f}")

    # Weighted fair R²
    total_n = sum(results[ELEMENTS[e]]["n_atoms"] for e in ELEMENTS
                  if ELEMENTS[e] in results)
    for label, key in [("full", "full"), ("geo-only", "geo")]:
        w_r2 = sum(results[ELEMENTS[e]]["progressive"]["fair"][key]
                   * results[ELEMENTS[e]]["n_atoms"] / total_n
                   for e in ELEMENTS if ELEMENTS[e] in results)
        print(f"  {'weighted':>5s}  {label + ':':>8s} {w_r2:.4f}")

    with open(out / "geometry_only_basis.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults in {out}/")
