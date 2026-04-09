#!/usr/bin/env python3
"""Physics-grounded calibration: per-protein normalized kernels +
fair categorical scalars + per-element ridge.

The "fair" scalars are physical identifiers, not features:
  - element (H/C/N/O) — handled by stratification
  - mutation_type (which ring was removed) — ring current I is type-dependent
  - kernel_scales (what per-protein normalization stripped) — bridges angular/magnitude

Optional continuous physics scalars tested as interaction terms:
  - mopac_valence (sum of Wiberg bond orders) — electronic environment
  - nearest_bond_order — bond anisotropy depends on this

Usage:
    cd learn/src
    python3 -c "
    from mutation_set.config import load_config
    import sys; cfg = load_config('calibration.toml')
    sys.path.insert(0, str(cfg.paths.sdk))
    from actual_physics.physics_calibration import run
    run(cfg)
    "
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from secondary.loader import iter_proteins, setup_sdk, ELEMENT_NAMES, RING_TYPE_NAMES
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels, normalize_kernels


def _core_indices(layout):
    return list(range(layout.efg_end))


def _ridge(X, y, lam=1.0):
    XtX = X.T @ X + lam * np.eye(X.shape[1])
    w = np.linalg.solve(XtX, X.T @ y)
    pred = X @ w
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean(axis=0)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
    return w, r2


def run(cfg: Config, max_proteins: int = 0):
    out = Path("output/actual_physics/physics_calibration")
    out.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    core = _core_indices(layout)
    core_names = [layout.names[i] for i in core]
    n_core = len(core)
    lam = getattr(cfg.secondary, 'ridge_lambda', 1.0)

    ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O"}

    # Accumulators per element
    elem_data = {e: {"K": [], "Y": [], "scales": [], "mut_type": [],
                      "valence": [], "bond_order": [], "dipole_mag": []}
                 for e in ELEMENTS}

    print(f"Core kernels: {n_core}")
    print("Loading proteins with per-protein normalization...")

    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx
        M = len(idx)

        # Assemble and normalize per-protein
        kernels_full = assemble_kernels(p, idx, layout)
        kernels_norm, scales = normalize_kernels(
            kernels_full.copy(), cfg.normalization.kernel_std_floor)
        K = kernels_norm[:, core, :]

        target = p.delta.shielding.T2[idx]

        # Core scales for the 55 kernels
        core_scales = scales[core]

        # Mutation type: which ring types exist in this protein
        rg = p.ring_geometry
        if rg is not None and rg.n_rings > 0:
            # Dominant ring type (most rings of this type)
            type_counts = np.bincount(rg.ring_type.astype(int),
                                       minlength=len(RING_TYPE_NAMES))
            mut_type = int(np.argmax(type_counts))
        else:
            mut_type = -1

        # MOPAC valence (sum of Wiberg bond orders per atom)
        valence = np.zeros(M)
        bond_order = np.zeros(M)
        dipole_mag = np.zeros(M)
        if p.mopac is not None:
            if hasattr(p.mopac, 'valence'):
                valence = p.mopac.valence[idx]
            if hasattr(p.mopac, 'bond_order_max'):
                bond_order = p.mopac.bond_order_max[idx]
            # Molecular dipole magnitude (scalar, same for all atoms in protein)
            if hasattr(p.mopac, 'dipole') and p.mopac.dipole is not None:
                dip = p.mopac.dipole
                if hasattr(dip, 'magnitude'):
                    dipole_mag[:] = dip.magnitude
                elif isinstance(dip, np.ndarray) and len(dip) >= 3:
                    dipole_mag[:] = np.linalg.norm(dip[:3])

        for j in range(M):
            e = int(rec.element[j])
            if e not in elem_data:
                continue
            elem_data[e]["K"].append(K[j])
            elem_data[e]["Y"].append(target[j])
            elem_data[e]["scales"].append(core_scales)
            elem_data[e]["mut_type"].append(mut_type)
            elem_data[e]["valence"].append(valence[j] if j < len(valence) else 0.0)
            elem_data[e]["bond_order"].append(bond_order[j] if j < len(bond_order) else 0.0)
            elem_data[e]["dipole_mag"].append(dipole_mag[j] if j < len(dipole_mag) else 0.0)

        if n_loaded % 50 == 0:
            print(f"  {n_loaded} proteins")

    print(f"Loaded {n_loaded} proteins\n")

    # ── Progressive calibration per element ──────────────────────────

    print(f"{'='*70}")
    print(f"  Progressive calibration: what helps?")
    print(f"{'='*70}")

    results = {}

    for e, ename in ELEMENTS.items():
        d = elem_data[e]
        if len(d["K"]) < 200:
            continue

        K_arr = np.array(d["K"])          # (N, 55, 5)
        Y_arr = np.array(d["Y"])          # (N, 5)
        scales = np.array(d["scales"])    # (N, 55)
        mut_types = np.array(d["mut_type"])
        valence = np.array(d["valence"])
        bond_order = np.array(d["bond_order"])

        N = K_arr.shape[0]
        X_base = K_arr.reshape(N, -1)     # (N, 275)

        # Level 0: raw kernels (per-protein normalized)
        _, r2_base = _ridge(X_base, Y_arr, lam)

        # Level 1: + kernel scale factors as interaction
        # scale_k × kernel_k = un-normalized kernel contribution
        # This lets ridge learn scale-dependent weights
        X_scales = np.zeros((N, n_core * 5))
        for ki in range(n_core):
            X_scales[:, ki*5:(ki+1)*5] = K_arr[:, ki, :] * scales[:, ki:ki+1]
        X_l1 = np.hstack([X_base, X_scales])
        _, r2_scales = _ridge(X_l1, Y_arr, lam)

        # Level 2: + mutation type as one-hot interaction
        unique_types = sorted(set(mut_types))
        unique_types = [t for t in unique_types if t >= 0]
        n_types = len(unique_types)
        type_map = {t: i for i, t in enumerate(unique_types)}

        X_mut = np.zeros((N, n_types * n_core * 5))
        for j in range(N):
            mt = mut_types[j]
            if mt < 0:
                continue
            ti = type_map.get(mt, -1)
            if ti < 0:
                continue
            offset = ti * n_core * 5
            X_mut[j, offset:offset + n_core * 5] = X_base[j]
        X_l2 = np.hstack([X_base, X_mut])
        _, r2_mut = _ridge(X_l2, Y_arr, lam)

        # Level 3: + valence as kernel interaction
        # valence × kernel = electronic-environment-modulated weight
        if np.std(valence) > 1e-8:
            v_z = (valence - valence.mean()) / max(valence.std(), 1e-8)
            X_val = X_base * v_z[:, np.newaxis]
            X_l3 = np.hstack([X_base, X_val])
            _, r2_valence = _ridge(X_l3, Y_arr, lam)
        else:
            r2_valence = r2_base

        # Level 4: + bond order as kernel interaction
        if np.std(bond_order) > 1e-8:
            bo_z = (bond_order - bond_order.mean()) / max(bond_order.std(), 1e-8)
            X_bo = X_base * bo_z[:, np.newaxis]
            X_l4 = np.hstack([X_base, X_bo])
            _, r2_bond = _ridge(X_l4, Y_arr, lam)
        else:
            r2_bond = r2_base

        # Level 5: + dipole magnitude as kernel interaction
        dipole = np.array(d["dipole_mag"])
        if np.std(dipole) > 1e-8:
            dip_z = (dipole - dipole.mean()) / max(dipole.std(), 1e-8)
            X_dip = X_base * dip_z[:, np.newaxis]
            X_l5 = np.hstack([X_base, X_dip])
            _, r2_dipole = _ridge(X_l5, Y_arr, lam)
        else:
            r2_dipole = r2_base

        # Fair set: base + scales + mutation type
        X_fair = np.hstack([X_base, X_scales, X_mut])
        _, r2_fair = _ridge(X_fair, Y_arr, lam)

        # Full physics set: base + scales + mut + valence + bond + dipole
        extras = [X_base, X_scales, X_mut]
        if np.std(valence) > 1e-8:
            v_z = (valence - valence.mean()) / max(valence.std(), 1e-8)
            extras.append(X_base * v_z[:, np.newaxis])
        if np.std(bond_order) > 1e-8:
            bo_z = (bond_order - bond_order.mean()) / max(bond_order.std(), 1e-8)
            extras.append(X_base * bo_z[:, np.newaxis])
        if np.std(dipole) > 1e-8:
            dip_z = (dipole - dipole.mean()) / max(dipole.std(), 1e-8)
            extras.append(X_base * dip_z[:, np.newaxis])
        X_full = np.hstack(extras)
        _, r2_full = _ridge(X_full, Y_arr, lam)

        print(f"\n  {ename} ({N:,d} atoms):")
        print(f"    Base (55 normalized kernels):         R² = {r2_base:.4f}")
        print(f"    + kernel scales:                      R² = {r2_scales:.4f}  (+{r2_scales-r2_base:.4f})")
        print(f"    + mutation type:                      R² = {r2_mut:.4f}  (+{r2_mut-r2_base:.4f})")
        print(f"    + valence interaction:                R² = {r2_valence:.4f}  (+{r2_valence-r2_base:.4f})")
        print(f"    + bond order interaction:             R² = {r2_bond:.4f}  (+{r2_bond-r2_base:.4f})")
        print(f"    + dipole interaction:                 R² = {r2_dipole:.4f}  (+{r2_dipole-r2_base:.4f})")
        print(f"    Fair set (base+scales+mut):           R² = {r2_fair:.4f}  (+{r2_fair-r2_base:.4f})")
        print(f"    Full physics set:                     R² = {r2_full:.4f}  (+{r2_full-r2_base:.4f})")

        results[ename] = {
            "n_atoms": N,
            "base": round(r2_base, 4),
            "plus_scales": round(r2_scales, 4),
            "plus_mut_type": round(r2_mut, 4),
            "plus_valence": round(r2_valence, 4),
            "plus_bond_order": round(r2_bond, 4),
            "plus_dipole": round(r2_dipole, 4),
            "fair_set": round(r2_fair, 4),
            "full_physics": round(r2_full, 4),
        }

    # Weighted
    total_n = sum(results[ELEMENTS[e]]["n_atoms"] for e in ELEMENTS if ELEMENTS[e] in results)
    for level in ["base", "fair_set"]:
        w_r2 = sum(results[ELEMENTS[e]][level] * results[ELEMENTS[e]]["n_atoms"] / total_n
                    for e in ELEMENTS if ELEMENTS[e] in results)
        print(f"\n  Weighted {level}: {w_r2:.4f}")
        results[f"weighted_{level}"] = round(w_r2, 4)

    with open(out / "progressive_calibration.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults in {out}/")
