#!/usr/bin/env python3
"""Per-element kernel physics decomposition.

Tests the literature prediction (Boyd & Skrynnikov 2002; Sahakyan &
Vendruscolo 2013) that ring current dominates the L=2 shielding
perturbation for H atoms, while the electric field gradient dominates
for heavy atoms (C, N, O).

For each element, computes ridge R² using physics-based kernel groups:
  1. ring_current — BS + HM per ring type, RingSusc total, per-ring
                    BS/HM/Chi/HM_H (magnetic mechanism)
  2. efg          — Coulomb EFG backbone + aromatic, MOPAC EFG, APBS,
                    DeltaAPBS, Coulomb/MOPAC totals (Buckingham mechanism)
  3. bond_aniso   — McConnell + MOPAC McConnell categories and totals,
                    HBond total (bond magnetic anisotropy)
  4. quadrupole   — PQ per ring type, PQ total, per-ring PQ
                    (pi-electron quadrupole, electrostatic)
  5. dispersion   — Disp per ring type, per-ring DispChi
                    (van der Waals)

Also computes per-kernel correlations stratified by element.

Outputs in {output_dir}/element_physics/:
    element_group_r2.csv   — R² per element × kernel group
    element_kernel_corr.csv — per-kernel mean |cos| with target per element
    element_counts.csv     — atom counts per element
    element_forward.csv    — greedy forward selection per element

Usage:
    cd learn/src
    python -m secondary element_physics --config calibration.toml
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .loader import (
    iter_proteins, setup_sdk, ridge_fit,
    ELEMENT_NAMES,
)
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels


# ── Physics-based kernel groups ──────────────────────────────────────

def _classify_kernels(layout: KernelLayout) -> dict[str, list[int]]:
    """Map each kernel to a physics group by name pattern."""
    groups: dict[str, list[int]] = {
        "ring_current": [],
        "efg": [],
        "bond_aniso": [],
        "quadrupole": [],
        "dispersion": [],
    }

    for i, name in enumerate(layout.names):
        if name.startswith("BS_") or name.startswith("HM_"):
            if "ring" in name:
                # Per-ring BS/HM — ring current
                groups["ring_current"].append(i)
            else:
                # Per-type BS/HM — ring current
                groups["ring_current"].append(i)
        elif name.startswith("Chi_") or name.startswith("HM_H_"):
            groups["ring_current"].append(i)
        elif name == "RingSusc_total":
            groups["ring_current"].append(i)
        elif name.startswith("RbfBs") or name.startswith("AngBs"):
            # Radial/angular BS partitions — still ring current
            groups["ring_current"].append(i)
        elif name.startswith("EFG_") or name.startswith("MopacEFG_"):
            groups["efg"].append(i)
        elif name.startswith("APBS_") or name.startswith("DeltaAPBS_"):
            groups["efg"].append(i)
        elif name == "Coulomb_total" or name == "MopacCoulomb_total":
            groups["efg"].append(i)
        elif name.startswith("MC_") or name.startswith("MopacMC_"):
            groups["bond_aniso"].append(i)
        elif name == "MC_total" or name == "MopacMC_total":
            groups["bond_aniso"].append(i)
        elif name == "HBond_total":
            groups["bond_aniso"].append(i)
        elif name.startswith("PQ_"):
            if "ring" in name:
                groups["quadrupole"].append(i)
            else:
                groups["quadrupole"].append(i)
        elif name == "PQ_total":
            groups["quadrupole"].append(i)
        elif name.startswith("Disp_") or name.startswith("DispChi_"):
            groups["dispersion"].append(i)
        else:
            # Shouldn't happen — flag it
            print(f"  WARNING: unclassified kernel {i}: {name}")

    return groups


def _cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """Cosine similarity between two 5D T2 vectors."""
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na < 1e-15 or nb < 1e-15:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def _forward_selection(X: np.ndarray, y: np.ndarray,
                       names: list[str], lam: float,
                       max_steps: int = 20) -> list[tuple[str, float, float]]:
    """Greedy forward selection: add the kernel that most improves R².

    Returns list of (kernel_name, marginal_delta, cumulative_R²).
    """
    n_kernels = X.shape[1]
    selected: list[int] = []
    available = set(range(n_kernels))
    results = []
    prev_r2 = 0.0

    for step in range(min(max_steps, n_kernels)):
        best_k, best_r2 = -1, -1e30
        for k in available:
            trial = selected + [k]
            X_trial = X[:, trial, :].reshape(X.shape[0], -1)
            _, r2 = ridge_fit(X_trial, y, lam)
            if r2 > best_r2:
                best_k, best_r2 = k, r2

        if best_k < 0:
            break
        selected.append(best_k)
        available.discard(best_k)
        delta = best_r2 - prev_r2
        results.append((names[best_k], delta, best_r2))
        prev_r2 = best_r2

    return results


def run(cfg: Config, max_proteins: int = 0):
    out_dir = cfg.secondary.output_dir / "element_physics"
    out_dir.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    lam = cfg.secondary.ridge_lambda
    groups = _classify_kernels(layout)

    print(f"Kernel layout: {layout.n_kernels} kernels")
    for gname, indices in groups.items():
        print(f"  {gname}: {len(indices)} kernels")

    # Element accumulators
    ELEMENTS = [1, 6, 7, 8, 16]
    elem_kernels: dict[int, list] = {e: [] for e in ELEMENTS}
    elem_targets: dict[int, list] = {e: [] for e in ELEMENTS}
    all_kernels: list = []
    all_targets: list = []

    print("Loading proteins...")
    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx
        kernels = assemble_kernels(p, idx, layout)  # (M, K, 5)
        target = p.delta.shielding.T2[idx]           # (M, 5)

        for j in range(len(idx)):
            elem = int(rec.element[j])
            if elem in elem_kernels:
                elem_kernels[elem].append(kernels[j])
                elem_targets[elem].append(target[j])
            all_kernels.append(kernels[j])
            all_targets.append(target[j])

        if n_loaded % 100 == 0:
            total = sum(len(v) for v in elem_kernels.values())
            print(f"  {n_loaded} proteins, {total} atoms")

    print(f"Loaded {n_loaded} proteins")
    for elem in ELEMENTS:
        print(f"  {ELEMENT_NAMES.get(elem, '?')}: {len(elem_kernels[elem])} atoms")

    # ── Per-element × per-group ridge R² ─────────────────────────────

    print("\nComputing per-element R² by physics group...")

    rows_r2 = []
    for elem in ELEMENTS + [0]:  # 0 = "all"
        if elem == 0:
            K_arr = np.array(all_kernels)
            y_arr = np.array(all_targets)
            ename = "all"
        else:
            if len(elem_kernels[elem]) < 50:
                continue
            K_arr = np.array(elem_kernels[elem])
            y_arr = np.array(elem_targets[elem])
            ename = ELEMENT_NAMES[elem]

        n_atoms = K_arr.shape[0]

        # All kernels
        X_all = K_arr.reshape(n_atoms, -1)
        _, r2_all = ridge_fit(X_all, y_arr, lam)

        row = {"element": ename, "n_atoms": n_atoms, "all": f"{r2_all:.4f}"}

        # Each physics group
        for gname, gidx in groups.items():
            if len(gidx) == 0:
                row[gname] = "0.0000"
                continue
            X_g = K_arr[:, gidx, :].reshape(n_atoms, -1)
            _, r2_g = ridge_fit(X_g, y_arr, lam)
            row[gname] = f"{r2_g:.4f}"

        # BS-only (just the 8 per-type BS kernels — the Biot-Savart wire model)
        bs_idx = [i for i, n in enumerate(layout.names) if n.startswith("BS_") and "ring" not in n]
        if bs_idx:
            X_bs = K_arr[:, bs_idx, :].reshape(n_atoms, -1)
            _, r2_bs = ridge_fit(X_bs, y_arr, lam)
            row["bs_only"] = f"{r2_bs:.4f}"

        # HM-only (just the 8 per-type HM kernels — the surface integral)
        hm_idx = [i for i, n in enumerate(layout.names) if n.startswith("HM_") and "ring" not in n]
        if hm_idx:
            X_hm = K_arr[:, hm_idx, :].reshape(n_atoms, -1)
            _, r2_hm = ridge_fit(X_hm, y_arr, lam)
            row["hm_only"] = f"{r2_hm:.4f}"

        # EFG aromatic only (MopacEFG_aro — the single dominant kernel from workbench)
        efg_aro_idx = [layout.index("MopacEFG_aro")]
        X_ea = K_arr[:, efg_aro_idx, :].reshape(n_atoms, -1)
        _, r2_ea = ridge_fit(X_ea, y_arr, lam)
        row["mopac_efg_aro"] = f"{r2_ea:.4f}"

        rows_r2.append(row)
        print(f"  {ename:3s} ({n_atoms:7d}): "
              f"ring_current={row['ring_current']}  efg={row['efg']}  "
              f"bond={row['bond_aniso']}  quad={row['quadrupole']}  "
              f"BS={row.get('bs_only','?')}  HM={row.get('hm_only','?')}  "
              f"EFG_aro={row['mopac_efg_aro']}  all={row['all']}")

    fieldnames = ["element", "n_atoms", "ring_current", "efg", "bond_aniso",
                  "quadrupole", "dispersion", "bs_only", "hm_only",
                  "mopac_efg_aro", "all"]
    with open(out_dir / "element_group_r2.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows_r2)

    # ── Per-kernel mean |cosine| with target, stratified by element ───

    print("\nComputing per-kernel correlations by element...")

    rows_corr = []
    for ki, kname in enumerate(layout.names):
        row = {"kernel": kname}
        for elem in ELEMENTS + [0]:
            if elem == 0:
                K_arr = np.array(all_kernels)
                y_arr = np.array(all_targets)
                ename = "all"
            else:
                if len(elem_kernels[elem]) < 50:
                    continue
                K_arr = np.array(elem_kernels[elem])
                y_arr = np.array(elem_targets[elem])
                ename = ELEMENT_NAMES[elem]

            # Mean |cosine similarity| between kernel T2 and target T2
            cos_vals = []
            for j in range(K_arr.shape[0]):
                c = _cosine_similarity(K_arr[j, ki, :], y_arr[j])
                cos_vals.append(abs(c))
            row[f"cos_{ename}"] = f"{np.mean(cos_vals):.4f}" if cos_vals else "0.0000"
        rows_corr.append(row)

    corr_fields = ["kernel"] + [f"cos_{ELEMENT_NAMES.get(e, 'all')}" for e in ELEMENTS + [0]]
    with open(out_dir / "element_kernel_corr.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=corr_fields)
        w.writeheader()
        w.writerows(rows_corr)

    # ── Forward selection per element ────────────────────────────────

    print("\nForward selection per element (top 15)...")

    for elem in ELEMENTS + [0]:
        if elem == 0:
            K_arr = np.array(all_kernels)
            y_arr = np.array(all_targets)
            ename = "all"
        else:
            if len(elem_kernels[elem]) < 200:
                continue
            K_arr = np.array(elem_kernels[elem])
            y_arr = np.array(elem_targets[elem])
            ename = ELEMENT_NAMES[elem]

        results = _forward_selection(K_arr, y_arr, layout.names, lam, max_steps=15)

        print(f"\n  {ename} ({K_arr.shape[0]} atoms):")
        rows_fwd = []
        for rank, (kname, delta, cum_r2) in enumerate(results, 1):
            print(f"    {rank:2d}. {kname:25s}  +{delta:.4f}  R²={cum_r2:.4f}")
            rows_fwd.append({"rank": rank, "kernel": kname,
                             "delta": f"{delta:.4f}", "r2": f"{cum_r2:.4f}"})

        with open(out_dir / f"forward_{ename}.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["rank", "kernel", "delta", "r2"])
            w.writeheader()
            w.writerows(rows_fwd)

    # ── Atom counts ──────────────────────────────────────────────────

    with open(out_dir / "element_counts.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["element", "n_atoms"])
        w.writeheader()
        for elem in ELEMENTS:
            w.writerow({"element": ELEMENT_NAMES.get(elem, str(elem)),
                        "n_atoms": len(elem_kernels[elem])})

    print(f"\nResults in {out_dir}/")
