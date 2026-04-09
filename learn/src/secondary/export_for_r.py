#!/usr/bin/env python3
"""Export kernel matrix + scalars + target + metadata for R analysis.

One Python run, unlimited R experiments.  Writes compact binary files
(numpy .npz) that R reads via reticulate, plus CSV metadata.

Outputs in {output_dir}/r_export/:
    kernels.npz        — (N, K, 5) float32 kernel T2 matrix
    target.npz         — (N, 5) float32 DFT delta T2
    scalars.npz        — (N, S) float32 scalar features
    kernel_scales.npz  — (N, K) float32 pre-normalization stds
    metadata.csv       — per-atom: protein_id, atom_index, element,
                         ring_dist, ring_type_mask, nearest_ring_type
    protein_meta.csv   — per-protein: pid, n_atoms, ring types present
    kernel_names.csv   — kernel index → name mapping
    scalar_layout.csv  — scalar block names and widths
    config.json        — kernel count, scalar count, ridge lambda, etc.

Usage:
    cd learn/src
    python -m secondary export_for_r --config calibration.toml
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from .loader import (
    iter_proteins, nearest_ring_type, setup_sdk,
    ELEMENT_NAMES, RING_TYPE_NAMES, N_RING_TYPES,
)
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels, normalize_kernels
from mutation_set.scalars import assemble_scalars


def run(cfg: Config, max_proteins: int = 0):
    out_dir = cfg.secondary.output_dir / "r_export"
    out_dir.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)

    all_kernels = []     # (N, K, 5) unnormalized
    all_norm_k = []      # (N, K, 5) per-protein normalized
    all_scales = []      # (K,) per protein
    all_target = []      # (N, 5)
    all_scalars = []     # (N, S)
    meta_rows = []       # per-atom metadata
    protein_rows = []    # per-protein metadata
    scalar_layout = None

    print("Loading proteins for R export...")
    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        p = rec.protein
        idx = rec.matched_idx
        M = len(idx)

        kernels_raw = assemble_kernels(p, idx, layout)
        kernels_norm, k_scales = normalize_kernels(
            kernels_raw.copy(), cfg.normalization.kernel_std_floor)
        scalars, sl = assemble_scalars(p, idx, k_scales, cfg)
        target = p.delta.shielding.T2[idx]

        all_kernels.append(kernels_raw)
        all_norm_k.append(kernels_norm)
        all_scales.append(np.tile(k_scales, (M, 1)))
        all_target.append(target)
        all_scalars.append(scalars)

        if scalar_layout is None:
            scalar_layout = sl

        # Per-protein metadata
        protein_rows.append({
            "protein_id": rec.pid,
            "n_atoms": M,
            "has_mopac": int(rec.has_mopac),
            **{f"has_{RING_TYPE_NAMES[i]}": int(i in rec.protein_ring_types)
               for i in range(N_RING_TYPES)},
        })

        # Per-atom metadata
        for j in range(M):
            nrt = nearest_ring_type(p, int(idx[j]))
            meta_rows.append({
                "protein_id": rec.pid,
                "atom_index": int(idx[j]),
                "element": ELEMENT_NAMES.get(int(rec.element[j]),
                                              str(int(rec.element[j]))),
                "ring_dist": round(float(rec.ring_dist[j]), 3),
                "ring_type_mask": int(rec.atom_ring_masks[j]),
                "nearest_ring_type": (RING_TYPE_NAMES[nrt]
                                      if nrt is not None else "none"),
            })

        n_loaded += 1
        if n_loaded % 100 == 0:
            print(f"  {n_loaded} proteins...")

    print(f"  Loaded {n_loaded} proteins, {len(meta_rows)} atoms")

    if not all_kernels:
        print("No data.")
        return

    # Stack and save
    K_raw = np.vstack(all_kernels).astype(np.float32)
    K_norm = np.vstack(all_norm_k).astype(np.float32)
    T = np.vstack(all_target).astype(np.float32)
    S = np.vstack(all_scalars).astype(np.float32)
    KS = np.vstack(all_scales).astype(np.float32)

    np.savez_compressed(out_dir / "kernels_raw.npz", data=K_raw)
    np.savez_compressed(out_dir / "kernels_norm.npz", data=K_norm)
    np.savez_compressed(out_dir / "target.npz", data=T)
    np.savez_compressed(out_dir / "scalars.npz", data=S)
    np.savez_compressed(out_dir / "kernel_scales.npz", data=KS)

    print(f"  Kernels: {K_raw.shape}  Target: {T.shape}  "
          f"Scalars: {S.shape}")

    # CSV metadata
    _write_csv(out_dir / "metadata.csv", meta_rows)
    _write_csv(out_dir / "protein_meta.csv", protein_rows)

    # Kernel names
    kn_rows = [{"index": i, "name": n, "group": _kernel_group(i, layout)}
               for i, n in enumerate(layout.names)]
    _write_csv(out_dir / "kernel_names.csv", kn_rows)

    # Scalar layout
    sl_rows = []
    offset = 0
    for name, width, cat in zip(scalar_layout.names, scalar_layout.widths,
                                 scalar_layout.categorical_mask):
        sl_rows.append({"name": name, "offset": offset, "width": width,
                        "categorical": int(cat)})
        offset += width
    _write_csv(out_dir / "scalar_layout.csv", sl_rows)

    # Config summary
    config = {
        "n_atoms": len(meta_rows),
        "n_proteins": n_loaded,
        "n_kernels": layout.n_kernels,
        "n_scalars": scalar_layout.total,
        "ridge_lambda": cfg.secondary.ridge_lambda,
        "kernel_groups": {
            "ring_type": [0, layout.ring_type_end],
            "bond_cat": [layout.ring_type_end, layout.bond_cat_end],
            "total": [layout.bond_cat_end, layout.total_end],
            "efg": [layout.total_end, layout.efg_end],
            "per_ring": [layout.efg_end, layout.n_kernels],
        },
        "per_ring_k": layout.per_ring_k,
        "per_ring_calcs": layout.per_ring_calcs,
        "features_dir": str(cfg.paths.features),
    }
    with open(out_dir / "config.json", "w") as f:
        json.dump(config, f, indent=2)

    print(f"\n  R export complete: {out_dir}/")
    print(f"  Load in R with: source('load_export.R')")


def _kernel_group(idx: int, layout: KernelLayout) -> str:
    if idx < layout.ring_type_end:
        return "ring_type"
    elif idx < layout.bond_cat_end:
        return "bond_cat"
    elif idx < layout.total_end:
        return "total"
    elif idx < layout.efg_end:
        return "efg"
    else:
        return "per_ring"


def _write_csv(path: Path, rows: list[dict]):
    fields = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  Wrote {path} ({len(rows)} rows)")
