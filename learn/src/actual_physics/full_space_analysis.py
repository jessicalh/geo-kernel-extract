#!/usr/bin/env python3
"""Full kernel space analysis for Stage 1 thesis.

Produces the complete characterisation of the 55-kernel T2 space:
per-element, raw vs normalised, group-to-group and within-group
cosine matrices, per-kernel and per-group R², forward selection.

Outputs JSON + CSV for R graphics.  All numbers reproducible from
calibration.toml config.

Usage:
    cd learn/src
    python3 -c "
    from mutation_set.config import load_config
    import sys; cfg = load_config('calibration.toml')
    sys.path.insert(0, str(cfg.paths.sdk))
    from actual_physics.full_space_analysis import run
    run(cfg)
    "
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from secondary.loader import iter_proteins, ELEMENT_NAMES
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels, normalize_kernels


# ── Physics group classification ─────────────────────────────────

def classify_kernels(layout: KernelLayout) -> dict[int, str]:
    """Map kernel index → physics group name."""
    n_core = layout.efg_end
    groups = {}
    for i in range(n_core):
        n = layout.names[i]
        if n.startswith("BS_") or n.startswith("HM_") or n == "RingSusc_total":
            groups[i] = "ring_current"
        elif n.startswith("Disp_"):
            groups[i] = "dispersion"
        elif n.startswith("PQ_") or n == "PQ_total":
            groups[i] = "quadrupole"
        elif n.startswith("MC_") or n == "HBond_total":
            groups[i] = "bond_aniso"
        elif n.startswith("MopacMC_") or n == "MopacMC_total":
            groups[i] = "mopac_bond"
        elif n in ("Coulomb_total", "EFG_bb", "EFG_aro"):
            groups[i] = "ff14sb_efg"
        elif n in ("MopacEFG_bb", "MopacEFG_aro", "MopacCoulomb_total"):
            groups[i] = "mopac_efg"
        elif n in ("APBS_EFG", "DeltaAPBS_EFG"):
            groups[i] = "solvation"
        else:
            groups[i] = "other"
    return groups


GROUP_ORDER = [
    "ring_current", "ff14sb_efg", "mopac_efg", "bond_aniso",
    "mopac_bond", "quadrupole", "dispersion", "solvation",
]

ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O"}


# ── Core functions ───────────────────────────────────────────────

def cosine_matrix(K: np.ndarray) -> np.ndarray:
    """(N, n_kernels, 5) → (n_kernels, n_kernels) mean |cos| matrix."""
    nk = K.shape[1]
    mags = np.linalg.norm(K, axis=2)
    C = np.full((nk, nk), np.nan)
    for i in range(nk):
        for j in range(i, nk):
            mi, mj = mags[:, i], mags[:, j]
            mask = (mi > 1e-15) & (mj > 1e-15)
            if mask.sum() < 50:
                continue
            dots = np.sum(K[mask, i, :] * K[mask, j, :], axis=1) / (mi[mask] * mj[mask])
            C[i, j] = C[j, i] = float(np.mean(np.abs(dots)))
    return C


def kernel_target_cosines(K: np.ndarray, Y: np.ndarray) -> list[float]:
    """Per-kernel mean |cos| with target.  Returns list of length n_kernels."""
    nk = K.shape[1]
    mags_t = np.linalg.norm(Y, axis=1)
    result = []
    for ki in range(nk):
        mags_k = np.linalg.norm(K[:, ki, :], axis=1)
        mask = (mags_k > 1e-15) & (mags_t > 1e-15)
        if mask.sum() < 50:
            result.append(float('nan'))
            continue
        dots = np.sum(K[mask, ki, :] * Y[mask], axis=1) / (mags_k[mask] * mags_t[mask])
        result.append(float(np.mean(np.abs(dots))))
    return result


def ridge_r2(X: np.ndarray, y: np.ndarray, lam: float = 1e-2) -> float:
    """Ridge regression R²."""
    if X.shape[0] < X.shape[1] + 10:
        return float('nan')
    XtX = X.T @ X + lam * np.eye(X.shape[1])
    w = np.linalg.solve(XtX, X.T @ y)
    pred = X @ w
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean(axis=0)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0


def forward_selection(K: np.ndarray, Y: np.ndarray, names: list[str],
                      groups: dict[int, str], lam: float,
                      max_steps: int = 20) -> list[dict]:
    """Greedy forward selection.  Returns list of step dicts."""
    N, nk, _ = K.shape
    selected = []
    available = set(range(nk))
    prev_r2 = 0.0
    results = []
    for _ in range(min(max_steps, nk)):
        best_k, best_r2 = -1, -1e30
        for k in available:
            trial = selected + [k]
            Xt = K[:, trial, :].reshape(N, -1)
            r2 = ridge_r2(Xt, Y, lam)
            if r2 > best_r2:
                best_k, best_r2 = k, r2
        if best_k < 0:
            break
        selected.append(best_k)
        available.discard(best_k)
        delta = best_r2 - prev_r2
        results.append({
            "rank": len(selected),
            "kernel": names[best_k],
            "group": groups.get(best_k, "other"),
            "delta": round(delta, 5),
            "cumulative": round(best_r2, 5),
        })
        prev_r2 = best_r2
    return results


def group_r2(K: np.ndarray, Y: np.ndarray, groups: dict[int, str],
             lam: float) -> dict[str, float]:
    """R² for each physics group (all kernels in group, simultaneously)."""
    N = K.shape[0]
    result = {}
    for g in GROUP_ORDER:
        gidx = [i for i, grp in groups.items() if grp == g]
        if not gidx:
            continue
        Xg = K[:, gidx, :].reshape(N, -1)
        result[g] = round(ridge_r2(Xg, Y, lam), 5)
    # All
    result["ALL"] = round(ridge_r2(K.reshape(N, -1), Y, lam), 5)
    return result


def group_cosine_summary(C: np.ndarray, groups: dict[int, str]) -> dict:
    """Group-to-group mean |cos| from the full cosine matrix."""
    result = {}
    for g1 in GROUP_ORDER:
        idx1 = [i for i, grp in groups.items() if grp == g1]
        if not idx1:
            continue
        row = {}
        for g2 in GROUP_ORDER:
            idx2 = [i for i, grp in groups.items() if grp == g2]
            if not idx2:
                continue
            if g1 == g2:
                pairs = [(i, j) for i in idx1 for j in idx2 if i < j]
            else:
                pairs = [(i, j) for i in idx1 for j in idx2]
            vals = [C[i, j] for i, j in pairs if not np.isnan(C[i, j])]
            row[g2] = round(np.mean(vals), 4) if vals else None
        result[g1] = row
    return result


# ── Main ─────────────────────────────────────────────────────────

def run(cfg: Config, max_proteins: int = 0):
    out = Path("output/actual_physics/full_space")
    out.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    n_core = layout.efg_end
    groups = classify_kernels(layout)
    names = layout.names[:n_core]
    lam = getattr(cfg.secondary, 'ridge_lambda', 1e-2)

    # ── Load data ────────────────────────────────────────────────
    elem_data = {e: {"K_raw": [], "K_norm": [], "Y": []} for e in ELEMENTS}

    print(f"Core kernels: {n_core}")
    print("Loading proteins...")
    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx

        kernels_raw = assemble_kernels(p, idx, layout)[:, :n_core, :]
        kernels_norm = kernels_raw.copy()
        kernels_norm, _ = normalize_kernels(kernels_norm, cfg.normalization.kernel_std_floor)
        target = p.delta.shielding.T2[idx]
        elem = rec.element

        for j in range(len(idx)):
            e = int(elem[j])
            if e in elem_data:
                elem_data[e]["K_raw"].append(kernels_raw[j])
                elem_data[e]["K_norm"].append(kernels_norm[j])
                elem_data[e]["Y"].append(target[j])

        if n_loaded % 100 == 0:
            print(f"  {n_loaded} proteins")

    print(f"Loaded {n_loaded} proteins\n")

    # ── Per-element analysis ─────────────────────────────────────
    all_results = {}

    for e, en in ELEMENTS.items():
        d = elem_data[e]
        N = len(d["K_raw"])
        if N < 500:
            continue

        K_raw = np.array(d["K_raw"])
        K_norm = np.array(d["K_norm"])
        Y = np.array(d["Y"])

        print("=" * 80)
        print(f"  {en}: {N:,d} atoms")
        print("=" * 80)

        er = {"n_atoms": N}

        # ── Cosine matrices (raw only — norm is identical) ───────
        print("  Computing cosine matrix...")
        C = cosine_matrix(K_raw)

        # Full per-kernel cosine matrix as CSV
        csv_path = out / f"cosine_matrix_{en}.csv"
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow([""] + names)
            for i in range(n_core):
                row = [names[i]] + [f"{C[i,j]:.4f}" if not np.isnan(C[i,j]) else ""
                                     for j in range(n_core)]
                w.writerow(row)
        print(f"  Saved {csv_path}")

        er["group_cosines"] = group_cosine_summary(C, groups)

        # ── Kernel → target cosine ───────────────────────────────
        kt_raw = kernel_target_cosines(K_raw, Y)
        kt_norm = kernel_target_cosines(K_norm, Y)

        er["kernel_target_cos"] = {
            names[i]: {"raw": round(kt_raw[i], 4) if not np.isnan(kt_raw[i]) else None,
                       "norm": round(kt_norm[i], 4) if not np.isnan(kt_norm[i]) else None,
                       "group": groups[i]}
            for i in range(n_core)
        }

        # ── Per-kernel R² ────────────────────────────────────────
        kr2_raw = {}
        kr2_norm = {}
        for ki in range(n_core):
            kr2_raw[names[ki]] = round(ridge_r2(K_raw[:, ki, :], Y, lam), 5)
            kr2_norm[names[ki]] = round(ridge_r2(K_norm[:, ki, :], Y, lam), 5)
        er["per_kernel_r2"] = {"raw": kr2_raw, "norm": kr2_norm}

        # ── Per-group R² ─────────────────────────────────────────
        gr2_raw = group_r2(K_raw, Y, groups, lam)
        gr2_norm = group_r2(K_norm, Y, groups, lam)
        er["per_group_r2"] = {"raw": gr2_raw, "norm": gr2_norm}

        print(f"\n  Per-group R² (raw / norm):")
        for g in GROUP_ORDER:
            if g in gr2_raw:
                print(f"    {g:>14s}: {gr2_raw[g]:.4f} / {gr2_norm[g]:.4f}  "
                      f"({gr2_norm[g] - gr2_raw[g]:+.4f})")
        print(f"    {'ALL':>14s}: {gr2_raw['ALL']:.4f} / {gr2_norm['ALL']:.4f}  "
              f"({gr2_norm['ALL'] - gr2_raw['ALL']:+.4f})")

        # ── Forward selection (both raw and norm) ────────────────
        print(f"\n  Forward selection (raw, top 20):")
        fwd_raw = forward_selection(K_raw, Y, names, groups, lam, 20)
        for step in fwd_raw:
            print(f"    {step['rank']:2d}. {step['kernel']:>25s}  "
                  f"{step['delta']:+.4f}  {step['cumulative']:.4f}  {step['group']}")

        print(f"\n  Forward selection (normalised, top 20):")
        fwd_norm = forward_selection(K_norm, Y, names, groups, lam, 20)
        for step in fwd_norm:
            print(f"    {step['rank']:2d}. {step['kernel']:>25s}  "
                  f"{step['delta']:+.4f}  {step['cumulative']:.4f}  {step['group']}")

        er["forward_selection"] = {"raw": fwd_raw, "norm": fwd_norm}

        # ── Within-group cosine detail ───────────────────────────
        print(f"\n  Within-group cosine structure:")
        within = {}
        for g in GROUP_ORDER:
            gidx = [i for i, grp in groups.items() if grp == g]
            if len(gidx) < 2:
                continue
            pairs = []
            for ii, i in enumerate(gidx):
                for j in gidx[ii + 1:]:
                    if not np.isnan(C[i, j]):
                        pairs.append({
                            "k1": names[i], "k2": names[j],
                            "cos": round(C[i, j], 4),
                        })
            if not pairs:
                continue
            cos_vals = [p["cos"] for p in pairs]
            within[g] = {
                "n_kernels": len(gidx),
                "n_pairs": len(pairs),
                "mean_cos": round(np.mean(cos_vals), 4),
                "min_cos": round(min(cos_vals), 4),
                "max_cos": round(max(cos_vals), 4),
                "std_cos": round(np.std(cos_vals), 4),
                "pairs": sorted(pairs, key=lambda x: -x["cos"]),
            }
            print(f"    {g:>14s} ({len(gidx)} kernels, {len(pairs)} pairs): "
                  f"mean={np.mean(cos_vals):.3f}  "
                  f"range=[{min(cos_vals):.3f}, {max(cos_vals):.3f}]  "
                  f"std={np.std(cos_vals):.3f}")

            # Top 3 most aligned and most independent pairs
            sorted_pairs = sorted(pairs, key=lambda x: -x["cos"])
            if len(sorted_pairs) > 3:
                print(f"      Most aligned:     {sorted_pairs[0]['k1']:>20s} ↔ {sorted_pairs[0]['k2']:<20s} cos={sorted_pairs[0]['cos']:.3f}")
                print(f"      Most independent: {sorted_pairs[-1]['k1']:>20s} ↔ {sorted_pairs[-1]['k2']:<20s} cos={sorted_pairs[-1]['cos']:.3f}")

        er["within_group"] = within

        # ── Export per-kernel CSV for R ──────────────────────────
        csv_path = out / f"per_kernel_{en}.csv"
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["kernel", "group", "cos_target_raw", "cos_target_norm",
                        "r2_raw", "r2_norm"])
            for ki in range(n_core):
                w.writerow([
                    names[ki], groups[ki],
                    f"{kt_raw[ki]:.5f}" if not np.isnan(kt_raw[ki]) else "",
                    f"{kt_norm[ki]:.5f}" if not np.isnan(kt_norm[ki]) else "",
                    f"{kr2_raw[names[ki]]:.5f}",
                    f"{kr2_norm[names[ki]]:.5f}",
                ])
        print(f"\n  Saved {csv_path}")

        all_results[en] = er
        print()

    # ── Save JSON ────────────────────────────────────────────────
    json_path = out / "full_space_analysis.json"
    with open(json_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"Results in {out}/")
    print(f"JSON: {json_path}")
    print(f"CSVs: cosine_matrix_{{H,C,N,O}}.csv, per_kernel_{{H,C,N,O}}.csv")
