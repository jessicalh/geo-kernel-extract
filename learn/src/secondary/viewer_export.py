#!/usr/bin/env python3
"""Tool 6: Viewer annotation export.

Produces a JSON file per protein with per-atom secondary analysis
annotations, keyed by atom_index.  Designed for the Qt/VTK viewer
to overlay as colour maps on the molecule.

Each atom carries its annotation values PLUS reality-check fields
(element, residue_index, position) so the viewer can verify
referential integrity at load time — the viewer built this protein,
so it knows what atom 42 should be.

Annotations exported:
  - divergence_cosine:   classical vs MOPAC T2 cosine (coulomb pair)
  - divergence_mag:      |T2(classical) - T2(MOPAC)| (ppm)
  - stratum:             ring-type stratum label
  - nearest_ring_type:   ring type name of nearest ring
  - self_fit_residual:   T2 magnitude of residual after 4-kernel self-fit
  - target_t2_mag:       DFT delta T2 magnitude (what we're trying to explain)
  - ridge_residual_mag:  T2 magnitude of residual after full 91-kernel ridge

Outputs in {output_dir}/viewer/:
  {protein_id}.json — one file per protein

Usage:
    cd learn/src
    python -m secondary viewer_export --config calibration.toml --protein A0A7C5FAR6
    python -m secondary viewer_export --config calibration.toml --max-proteins 10
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .loader import (
    iter_proteins, assign_stratum, nearest_ring_type, setup_sdk,
    ridge_fit, ELEMENT_NAMES, RING_TYPE_NAMES, N_RING_TYPES,
)
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels


def _cosine_5d(a: np.ndarray, b: np.ndarray) -> float:
    """5D cosine similarity between two vectors."""
    dot = float(np.dot(a, b))
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na < 1e-12 or nb < 1e-12:
        return float("nan")
    return dot / (na * nb)


def _resolve(protein, dotpath: str):
    obj = protein
    for part in dotpath.split("."):
        obj = getattr(obj, part, None)
        if obj is None:
            return None
    return obj


def _self_fit_residual(kernels: np.ndarray, target: np.ndarray,
                       atom_ring_masks: np.ndarray, layout: KernelLayout,
                       lam: float) -> np.ndarray:
    """Per-atom T2 residual magnitude from self-fit (4 ring-type kernels).

    For each atom, uses only the 4 kernels of its nearest ring type.
    Returns (M,) array of residual magnitudes.
    """
    M = len(target)
    residual_mag = np.full(M, np.nan)

    # Group atoms by nearest ring type
    for rt_idx in range(N_RING_TYPES):
        rt_mask = np.array([(m & (1 << rt_idx)) != 0 for m in atom_ring_masks])
        if rt_mask.sum() < 20:
            continue

        rt_cols = [rt_idx, 8 + rt_idx, 16 + rt_idx, 24 + rt_idx]
        rt_cols = [c for c in rt_cols if c < layout.ring_type_end]
        K_rt = kernels[rt_mask][:, rt_cols, :]
        T_rt = target[rt_mask]
        n = int(rt_mask.sum())

        X = K_rt.reshape(n, -1)
        if n < max(20, X.shape[1]):
            continue
        try:
            pred, _ = ridge_fit(X, T_rt, lam)
            resid = T_rt - pred
            mags = np.linalg.norm(resid, axis=1)
            residual_mag[rt_mask] = mags
        except np.linalg.LinAlgError:
            pass

    return residual_mag


def run(cfg: Config, max_proteins: int = 0, protein_filter: str = ""):
    out_dir = cfg.secondary.output_dir / "viewer"
    out_dir.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    lam = cfg.secondary.ridge_lambda

    print("Exporting viewer annotations...")
    n_exported = 0

    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        if protein_filter and rec.pid != protein_filter:
            continue

        p = rec.protein
        idx = rec.matched_idx
        M = len(idx)

        # Assemble kernels and target
        kernels = assemble_kernels(p, idx, layout)
        target = p.delta.shielding.T2[idx]

        # Full ridge residual (magnitude + T2 components for butterfly rendering)
        X_all = kernels.reshape(M, -1)
        ridge_residual_mag = np.full(M, np.nan)
        ridge_residual_t2 = np.zeros((M, 5))
        if M > max(20, X_all.shape[1] // 5):
            try:
                pred, _ = ridge_fit(X_all, target, lam)
                ridge_residual_t2 = target - pred
                ridge_residual_mag = np.linalg.norm(ridge_residual_t2, axis=1)
            except np.linalg.LinAlgError:
                pass

        # Self-fit residual
        self_residual = _self_fit_residual(
            kernels, target, rec.atom_ring_masks, layout, lam)

        # Divergence (coulomb pair)
        c_obj = _resolve(p, "coulomb.shielding")
        m_obj = _resolve(p, "mopac.coulomb.shielding")
        has_mopac = c_obj is not None and m_obj is not None

        # Build per-atom annotations
        atoms = {}
        for j in range(M):
            ai = int(idx[j])
            pos = p.pos[ai]
            elem = int(rec.element[j])
            nrt = nearest_ring_type(p, ai)

            annotation = {
                # Reality checks — viewer verifies these match its own data
                "element": ELEMENT_NAMES.get(elem, str(elem)),
                "element_number": elem,
                "residue_index": int(p.residue_index[ai]),
                "position": [round(float(pos[0]), 3),
                             round(float(pos[1]), 3),
                             round(float(pos[2]), 3)],

                # Analysis results
                "stratum": assign_stratum(int(rec.atom_ring_masks[j]),
                                          cfg.secondary.strata) or "unclassified",
                "nearest_ring_type": (RING_TYPE_NAMES[nrt]
                                      if nrt is not None else "none"),
                "ring_dist": round(float(rec.ring_dist[j]), 2),
                "target_t2_mag": round(float(np.linalg.norm(target[j])), 4),
                "ridge_residual_mag": round(float(ridge_residual_mag[j]), 4),
                "ridge_residual_t2": [round(float(v), 5)
                                      for v in ridge_residual_t2[j]],
                "target_t2": [round(float(v), 5) for v in target[j]],
                "self_fit_residual": round(float(self_residual[j]), 4),
            }

            if has_mopac:
                c_t2 = c_obj.T2[ai]
                m_t2 = m_obj.T2[ai]
                diff = c_t2 - m_t2
                annotation["divergence_cosine"] = round(
                    _cosine_5d(c_t2, m_t2), 4)
                annotation["divergence_mag"] = round(
                    float(np.linalg.norm(diff)), 4)

            atoms[str(ai)] = annotation

        # Write JSON
        output = {
            "protein_id": rec.pid,
            "n_matched_atoms": M,
            "n_total_atoms": int(p.n_atoms),
            "ring_types_present": sorted(
                RING_TYPE_NAMES[rt] for rt in rec.protein_ring_types),
            "ridge_lambda": lam,
            "n_kernels": layout.n_kernels,
            "reality_check": {
                "description": "element, residue_index, and position are "
                               "from the extraction. The viewer should "
                               "verify these match its own protein data "
                               "before applying annotations.",
                "position_tolerance_angstrom": 0.01,
            },
            "annotations": atoms,
        }

        path = out_dir / f"{rec.pid}.json"
        with open(path, "w") as f:
            json.dump(output, f, indent=2, allow_nan=True)

        n_exported += 1
        if n_exported % 50 == 0 or protein_filter:
            print(f"  {rec.pid}: {M} atoms, "
                  f"rings={output['ring_types_present']}")

    print(f"\n  Exported {n_exported} proteins to {out_dir}/")
