#!/usr/bin/env python3
"""Export per-atom physics annotation for the Qt/VTK viewer.

Produces a JSON file per protein with per-atom fields that the viewer
can map to colour, glyph size, or tensor ellipsoid overlays.

Fields per atom:
    pos            — [x, y, z] Angstroms
    element        — atomic number
    element_name   — H, C, N, O, S
    dist_to_ring   — distance to nearest removed ring (Angstroms)
    target_T2      — [5] DFT delta T2 components (m=-2..+2)
    predicted_T2   — [5] ridge-predicted T2 from 55 core kernels
    residual_T2    — [5] target - predicted
    target_mag     — |target_T2|
    residual_mag   — |residual_T2|
    fit_quality    — 1 - |residual|/|target| (1=perfect, 0=no signal)
    dominant_group — "ring_current", "efg", "bond_aniso", "quadrupole"
    dominant_kernel — name of the single kernel with highest |cos(k, target)|
    efg_cos        — cos(EFG_aro, target) — alignment with dominant mechanism
    nearest_ring_type — name of nearest ring type

Usage:
    cd learn/src
    python3 -m viewer_export.export_physics A0A3A5VK63
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from secondary.loader import ELEMENT_NAMES, RING_TYPE_NAMES
from mutation_set.config import load_config, Config
from mutation_set.kernels import KernelLayout, assemble_kernels, normalize_kernels


PHYSICS_GROUPS = {
    "ring_current": lambda n: (n.startswith("BS_") or n.startswith("HM_") or
                               n == "RingSusc_total" or n.startswith("Chi_") or
                               n.startswith("HM_H_")),
    "efg": lambda n: (n.startswith("EFG_") or n.startswith("MopacEFG_") or
                      n.startswith("APBS_") or n.startswith("DeltaAPBS_") or
                      n in ("Coulomb_total", "MopacCoulomb_total")),
    "bond_aniso": lambda n: (n.startswith("MC_") or n.startswith("MopacMC_") or
                             n in ("MC_total", "MopacMC_total", "HBond_total")),
    "quadrupole": lambda n: n.startswith("PQ_") or n == "PQ_total",
}


def _cos(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < 1e-15 or nb < 1e-15:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def _classify_kernel(name):
    for group, pred in PHYSICS_GROUPS.items():
        if pred(name):
            return group
    return "other"


def export_protein(cfg: Config, pid: str, out_dir: Path):
    """Export one protein's physics annotation as JSON."""
    from nmr_extract import load

    features = cfg.paths.features
    p = load(features / pid)
    if p.delta is None:
        print(f"  {pid}: no delta data")
        return None

    layout = KernelLayout.from_config(cfg)
    core = list(range(layout.efg_end))
    core_names = [layout.names[i] for i in core]

    idx = np.where(p.delta.scalars.matched_mask)[0]
    M = len(idx)

    # Assemble and normalize
    kernels_full = assemble_kernels(p, idx, layout)
    kernels_norm, scales = normalize_kernels(
        kernels_full.copy(), cfg.normalization.kernel_std_floor)
    K = kernels_norm[:, core, :]
    target = p.delta.shielding.T2[idx]

    # Ridge fit on this protein (per-protein ridge)
    X = K.reshape(M, -1)
    XtX = X.T @ X + 1.0 * np.eye(X.shape[1])
    w = np.linalg.solve(XtX, X.T @ target)
    predicted = X @ w

    # Ring geometry
    rc = p.ring_contributions
    pos = p.pos.data[idx]

    atoms = []
    for j in range(M):
        t = target[j]
        pred = predicted[j]
        res = t - pred
        t_mag = float(np.linalg.norm(t))
        r_mag = float(np.linalg.norm(res))

        # Fit quality
        fq = 1.0 - r_mag / t_mag if t_mag > 1e-10 else 0.0

        # Dominant group: which physics group has highest sum |cos(k, target)|
        group_scores = {g: 0.0 for g in PHYSICS_GROUPS}
        best_kernel = ""
        best_cos = 0.0
        for ki, kname in enumerate(core_names):
            c = abs(_cos(K[j, ki, :], t))
            group = _classify_kernel(kname)
            if group in group_scores:
                group_scores[group] += c
            if c > best_cos:
                best_cos = c
                best_kernel = kname

        dominant = max(group_scores, key=group_scores.get)

        # EFG alignment
        efg_aro_ki = core_names.index("EFG_aro") if "EFG_aro" in core_names else -1
        efg_cos = _cos(K[j, efg_aro_ki, :], t) if efg_aro_ki >= 0 else 0.0

        # Distance and ring type
        dist = float(p.delta.scalars.nearest_removed_ring_dist[idx[j]])
        nearest_rt = ""
        if rc is not None and rc.n_pairs > 0:
            atom_rc = rc.for_atom(idx[j])
            if atom_rc.n_pairs > 0:
                nearest = np.argmin(atom_rc.distance)
                rt_idx = int(atom_rc.ring_type[nearest])
                if rt_idx < len(RING_TYPE_NAMES):
                    nearest_rt = RING_TYPE_NAMES[rt_idx]

        elem = int(p.element[idx[j]])

        atoms.append({
            "atom_index": int(idx[j]),
            "pos": [round(float(pos[j, k]), 3) for k in range(3)],
            "element": elem,
            "element_name": ELEMENT_NAMES.get(elem, str(elem)),
            "dist_to_ring": round(dist, 2),
            "target_T2": [round(float(t[k]), 5) for k in range(5)],
            "predicted_T2": [round(float(pred[k]), 5) for k in range(5)],
            "residual_T2": [round(float(res[k]), 5) for k in range(5)],
            "target_mag": round(t_mag, 4),
            "residual_mag": round(r_mag, 4),
            "fit_quality": round(fq, 3),
            "dominant_group": dominant,
            "dominant_kernel": best_kernel,
            "efg_cos": round(efg_cos, 4),
            "nearest_ring_type": nearest_rt,
        })

    result = {
        "protein": pid,
        "n_atoms": M,
        "n_kernels": len(core_names),
        "kernel_names": core_names,
        "ring_types_present": [RING_TYPE_NAMES[int(t)]
                               for t in sorted(set(int(t) for t in p.ring_geometry.ring_type))
                               if int(t) < len(RING_TYPE_NAMES)]
                              if p.ring_geometry is not None else [],
        "atoms": atoms,
    }

    out_path = out_dir / f"{pid}_physics.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"  {pid}: {M} atoms, {len(result['ring_types_present'])} ring types → {out_path}")
    return result


def main():
    if len(sys.argv) < 2:
        print("Usage: python -m viewer_export.export_physics <PROTEIN_ID> [PROTEIN_ID2 ...]")
        sys.exit(1)

    cfg = load_config("calibration.toml")
    sys.path.insert(0, str(cfg.paths.sdk))

    out_dir = Path("output/viewer_export")
    out_dir.mkdir(parents=True, exist_ok=True)

    for pid in sys.argv[1:]:
        export_protein(cfg, pid, out_dir)


if __name__ == "__main__":
    main()
