#!/usr/bin/env python3
"""Tool 2: Ring-type conditional analysis.

What does the kernel vector look like near HIE vs PHE?  Does the
MOPAC/classical divergence concentrate at HIE?

Outputs in {output_dir}/ring_strata/:
    protein_tags.csv      — per-protein ring type inventory
    atom_tags.csv         — per-atom ring type exposure flags
    kernel_by_stratum.csv — mean T2 magnitude per kernel per stratum
    ridge_by_stratum.csv  — ridge R² per stratum

Usage:
    cd learn/src
    python -m secondary ring_strata --config calibration.toml
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .loader import (
    iter_proteins, assign_stratum, setup_sdk,
    ridge_fit, nearest_ring_type,
    ELEMENT_NAMES, RING_TYPE_NAMES, N_RING_TYPES,
)
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels, normalize_kernels


def run(cfg: Config, max_proteins: int = 0):
    out_dir = cfg.secondary.output_dir / "ring_strata"
    out_dir.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    strata = cfg.secondary.strata
    lam = cfg.secondary.ridge_lambda

    # Accumulate
    protein_rows = []
    atom_rows = []
    # Per-stratum: kernels (N, K, 5) and targets (N, 5)
    strata_kernels: dict[str, list] = {s: [] for s in strata}
    strata_targets: dict[str, list] = {s: [] for s in strata}

    print("Loading proteins...")
    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx
        M = len(idx)

        # Protein-level tags
        ring_counts = [0] * N_RING_TYPES
        rg = p.ring_geometry
        if rg is not None:
            for rt in rg.ring_type.astype(int):
                ring_counts[rt] += 1

        protein_rows.append({
            "protein_id": rec.pid,
            "n_matched_atoms": M,
            **{f"has_{RING_TYPE_NAMES[i]}": int(i in rec.protein_ring_types)
               for i in range(N_RING_TYPES)},
            **{f"n_rings_{RING_TYPE_NAMES[i]}": ring_counts[i]
               for i in range(N_RING_TYPES)},
        })

        # Assemble kernels for this protein
        kernels = assemble_kernels(p, idx, layout)  # (M, K, 5)
        target = p.delta.shielding.T2[idx]           # (M, 5)

        # Atom-level tags + stratum assignment
        for j in range(M):
            mask = int(rec.atom_ring_masks[j])
            nrt = nearest_ring_type(p, int(idx[j]))
            stratum = assign_stratum(mask, strata)

            atom_rows.append({
                "protein_id": rec.pid,
                "atom_index": int(idx[j]),
                "element": ELEMENT_NAMES.get(int(rec.element[j]),
                                              str(int(rec.element[j]))),
                "ring_dist": float(rec.ring_dist[j]),
                **{f"sees_{RING_TYPE_NAMES[i]}": int(bool(mask & (1 << i)))
                   for i in range(N_RING_TYPES)},
                "nearest_ring_type": (RING_TYPE_NAMES[nrt]
                                      if nrt is not None else "none"),
                "stratum": stratum if stratum else "unclassified",
            })

            # Accumulate into strata
            if stratum:
                strata_kernels[stratum].append(kernels[j])
                strata_targets[stratum].append(target[j])
            # "all" always gets everything
            if stratum != "all":
                strata_kernels["all"].append(kernels[j])
                strata_targets["all"].append(target[j])

    print(f"  Loaded {n_loaded} proteins, {len(atom_rows)} atoms")

    # Write protein tags
    if protein_rows:
        _write_csv(out_dir / "protein_tags.csv", protein_rows)

    # Write atom tags
    if atom_rows:
        _write_csv(out_dir / "atom_tags.csv", atom_rows)

    # Kernel magnitude by stratum
    kernel_summary = []
    for s in strata:
        if not strata_kernels[s]:
            continue
        K = np.array(strata_kernels[s])  # (N_s, n_kernels, 5)
        for k in range(layout.n_kernels):
            mags = np.linalg.norm(K[:, k, :], axis=1)
            kernel_summary.append({
                "stratum": s,
                "kernel_name": layout.names[k],
                "kernel_index": k,
                "n_atoms": len(K),
                "t2_mag_mean": float(mags.mean()),
                "t2_mag_std": float(mags.std()),
                "t2_mag_median": float(np.median(mags)),
                "frac_nonzero": float((mags > 1e-10).mean()),
            })
    if kernel_summary:
        _write_csv(out_dir / "kernel_by_stratum.csv", kernel_summary)

    # Ridge R² by stratum
    ridge_rows = []
    for s in strata:
        if not strata_kernels[s]:
            continue
        K = np.array(strata_kernels[s])
        T = np.array(strata_targets[s])
        N_s = len(T)

        # Need enough atoms for ridge regression
        n_features = layout.n_kernels * 5
        if N_s < max(20, n_features // 5):
            ridge_rows.append({
                "stratum": s,
                "n_atoms": N_s,
                "ridge_r2_all": float("nan"),
                "ridge_r2_ring": float("nan"),
                "ridge_r2_bond": float("nan"),
                "ridge_r2_efg": float("nan"),
            })
            continue

        X_all = K.reshape(N_s, -1)
        _, r2_all = ridge_fit(X_all, T, lam)

        # Kernel family slices
        X_ring = K[:, :layout.ring_type_end, :].reshape(N_s, -1)
        _, r2_ring = ridge_fit(X_ring, T, lam)

        X_bond = K[:, layout.ring_type_end:layout.bond_cat_end, :].reshape(N_s, -1)
        _, r2_bond = ridge_fit(X_bond, T, lam) if X_bond.shape[1] > 0 else (None, float("nan"))

        X_efg = K[:, layout.total_end:layout.efg_end, :].reshape(N_s, -1)
        _, r2_efg = ridge_fit(X_efg, T, lam) if X_efg.shape[1] > 0 else (None, float("nan"))

        ridge_rows.append({
            "stratum": s,
            "n_atoms": N_s,
            "ridge_r2_all": float(r2_all),
            "ridge_r2_ring": float(r2_ring),
            "ridge_r2_bond": float(r2_bond),
            "ridge_r2_efg": float(r2_efg),
        })

    if ridge_rows:
        _write_csv(out_dir / "ridge_by_stratum.csv", ridge_rows)

    # Print summary
    print(f"\n  Stratum summary:")
    for r in ridge_rows:
        print(f"    {r['stratum']:15s}  n={r['n_atoms']:6d}  "
              f"R²={r['ridge_r2_all']:.3f}")


def _write_csv(path: Path, rows: list[dict]):
    fields = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  Wrote {path} ({len(rows)} rows)")
