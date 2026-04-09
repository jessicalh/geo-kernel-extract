#!/usr/bin/env python3
"""Tool 3: Inter-kernel T2 structure (DFT-free).

Pairwise T2 cosine similarity between all 91 kernels, eigenspectrum
of the kernel covariance, stratified by ring-type contact.  No DFT
delta needed — this is about the geometry of the kernel space itself.

Outputs in {output_dir}/kernel_structure/:
    cosine_matrix.csv  — long format: kernel_i, kernel_j, cosine_sim, stratum
    eigenspectrum.csv  — eigenvalue rank + cumulative variance per stratum

Usage:
    cd learn/src
    python -m secondary kernel_structure --config calibration.toml
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .loader import (
    iter_proteins, assign_stratum, setup_sdk,
    RING_TYPE_NAMES,
)
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels


def _pairwise_cosine(K_flat: np.ndarray) -> np.ndarray:
    """Cosine similarity matrix from (n_kernels, N*5) matrix."""
    norms = np.linalg.norm(K_flat, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-12)
    normed = K_flat / norms
    return normed @ normed.T


def _eigenspectrum(K_flat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Eigenvalues and cumulative variance explained from (n_kernels, N*5)."""
    # Kernel covariance: (n_kernels, n_kernels)
    cov = K_flat @ K_flat.T / max(K_flat.shape[1], 1)
    eigvals = np.linalg.eigvalsh(cov)[::-1]  # descending
    eigvals = np.maximum(eigvals, 0)
    total = eigvals.sum()
    cumvar = np.cumsum(eigvals) / total if total > 1e-12 else np.zeros_like(eigvals)
    return eigvals, cumvar


def run(cfg: Config, max_proteins: int = 0):
    out_dir = cfg.secondary.output_dir / "kernel_structure"
    out_dir.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    strata = cfg.secondary.strata

    # Accumulate kernels and stratum labels
    strata_data: dict[str, list] = {s: [] for s in strata}

    print("Loading proteins...")
    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx
        kernels = assemble_kernels(p, idx, layout)  # (M, K, 5)

        for j in range(len(idx)):
            mask = int(rec.atom_ring_masks[j])
            s = assign_stratum(mask, strata)
            if s:
                strata_data[s].append(kernels[j])
            if s != "all":
                strata_data["all"].append(kernels[j])

    print(f"  Loaded {n_loaded} proteins")

    # Compute per stratum
    cosine_rows = []
    eigen_rows = []

    for s in strata:
        if not strata_data[s]:
            continue
        K = np.array(strata_data[s])  # (N_s, n_kernels, 5)
        N_s = len(K)
        n_k = layout.n_kernels
        print(f"  Stratum {s}: {N_s} atoms")

        # Reshape to (n_kernels, N_s * 5) for pairwise comparison
        K_flat = K.transpose(1, 0, 2).reshape(n_k, -1)  # (K, N*5)

        # Cosine matrix
        cos_mat = _pairwise_cosine(K_flat)
        for i in range(n_k):
            for j in range(i, n_k):
                cosine_rows.append({
                    "kernel_i": layout.names[i],
                    "kernel_j": layout.names[j],
                    "cosine_sim": float(cos_mat[i, j]),
                    "stratum": s,
                })

        # Eigenspectrum
        eigvals, cumvar = _eigenspectrum(K_flat)
        for rank, (ev, cv) in enumerate(zip(eigvals, cumvar)):
            eigen_rows.append({
                "stratum": s,
                "rank": rank + 1,
                "eigenvalue": float(ev),
                "cumulative_variance": float(cv),
            })

        # Summary: effective dimensionality (90% variance)
        dim90 = int(np.searchsorted(cumvar, 0.90)) + 1
        dim95 = int(np.searchsorted(cumvar, 0.95)) + 1
        print(f"    Effective dim: {dim90} (90% var), {dim95} (95% var)")

    if cosine_rows:
        _write_csv(out_dir / "cosine_matrix.csv", cosine_rows)
    if eigen_rows:
        _write_csv(out_dir / "eigenspectrum.csv", eigen_rows)


def _write_csv(path: Path, rows: list[dict]):
    fields = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  Wrote {path} ({len(rows)} rows)")
