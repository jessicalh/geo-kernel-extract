#!/usr/bin/env python3
"""Tool 1: Classical vs MOPAC T2 divergence.

Where do ff14SB charges and MOPAC Mulliken charges produce different
T2 angular patterns?  The delta between classical and MOPAC pairs is
a polarisation map — no DFT ground truth needed.

4 pairs compared:
    1. coulomb.shielding       vs  mopac.coulomb.shielding
    2. mcconnell.shielding     vs  mopac.mcconnell.shielding
    3. coulomb.efg_backbone    vs  mopac.coulomb.efg_backbone
    4. coulomb.efg_aromatic    vs  mopac.coulomb.efg_aromatic

Output CSVs in {output_dir}/divergence/:
    atom_divergence.csv     — per-atom, per-pair T2 difference
    summary_by_element.csv  — grouped means/stds
    summary_by_ring_type.csv — grouped by nearest ring type

Usage:
    cd learn/src
    python -m secondary divergence --config calibration.toml
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .loader import (
    iter_proteins, nearest_ring_type, setup_sdk,
    ELEMENT_NAMES, RING_TYPE_NAMES, N_RING_TYPES,
)
from mutation_set.config import Config

# Classical dotpath → MOPAC dotpath, display name
PAIRS = [
    ("coulomb.shielding",    "mopac.coulomb.shielding",    "coulomb_shielding"),
    ("mcconnell.shielding",  "mopac.mcconnell.shielding",  "mcconnell_shielding"),
    ("coulomb.efg_backbone", "mopac.coulomb.efg_backbone",  "efg_backbone"),
    ("coulomb.efg_aromatic", "mopac.coulomb.efg_aromatic",  "efg_aromatic"),
]


def _resolve(protein, dotpath: str):
    """Navigate p.mopac.coulomb.shielding etc."""
    obj = protein
    for part in dotpath.split("."):
        obj = getattr(obj, part, None)
        if obj is None:
            return None
    return obj


def _cosine_sim_5d(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Per-row 5D cosine similarity.  Returns (M,) array."""
    dot = np.sum(a * b, axis=1)
    na = np.linalg.norm(a, axis=1)
    nb = np.linalg.norm(b, axis=1)
    denom = na * nb
    valid = denom > 1e-12
    result = np.full(len(a), np.nan)
    result[valid] = dot[valid] / denom[valid]
    return result


def run(cfg: Config, max_proteins: int = 0):
    out_dir = cfg.secondary.output_dir / "divergence"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Accumulate rows
    atom_rows = []
    n_loaded = 0
    n_skipped_mopac = 0

    print("Loading proteins (MOPAC required)...")
    for rec in iter_proteins(cfg, require_mopac=True,
                             max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx
        M = len(idx)

        for classical_path, mopac_path, pair_name in PAIRS:
            c_obj = _resolve(p, classical_path)
            m_obj = _resolve(p, mopac_path)
            if c_obj is None or m_obj is None:
                continue

            c_t2 = c_obj.T2[idx]  # (M, 5)
            m_t2 = m_obj.T2[idx]  # (M, 5)
            diff = c_t2 - m_t2
            diff_mag = np.linalg.norm(diff, axis=1)
            cos = _cosine_sim_5d(c_t2, m_t2)
            c_mag = np.linalg.norm(c_t2, axis=1)
            m_mag = np.linalg.norm(m_t2, axis=1)
            mag_ratio = np.where(c_mag > 1e-12, m_mag / c_mag, np.nan)

            for j in range(M):
                nrt = nearest_ring_type(p, int(idx[j]))
                atom_rows.append({
                    "protein_id": rec.pid,
                    "atom_index": int(idx[j]),
                    "element": ELEMENT_NAMES.get(int(rec.element[j]),
                                                  str(int(rec.element[j]))),
                    "ring_dist": float(rec.ring_dist[j]),
                    "ring_type_mask": int(rec.atom_ring_masks[j]),
                    "nearest_ring_type": (RING_TYPE_NAMES[nrt]
                                          if nrt is not None else "none"),
                    "pair_name": pair_name,
                    "diff_m2": float(diff[j, 0]),
                    "diff_m1": float(diff[j, 1]),
                    "diff_0":  float(diff[j, 2]),
                    "diff_p1": float(diff[j, 3]),
                    "diff_p2": float(diff[j, 4]),
                    "diff_mag": float(diff_mag[j]),
                    "cosine_sim": float(cos[j]),
                    "mag_ratio": float(mag_ratio[j]),
                })

    print(f"  Loaded {n_loaded} proteins with MOPAC")
    print(f"  {len(atom_rows)} atom × pair rows")

    if not atom_rows:
        print("No data — nothing to write.")
        return

    # Write atom-level CSV
    fields = list(atom_rows[0].keys())
    atom_csv = out_dir / "atom_divergence.csv"
    with open(atom_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(atom_rows)
    print(f"  Wrote {atom_csv} ({len(atom_rows)} rows)")

    # Summary by element
    _write_summary(atom_rows, "element", out_dir / "summary_by_element.csv")

    # Summary by nearest ring type
    _write_summary(atom_rows, "nearest_ring_type",
                   out_dir / "summary_by_ring_type.csv")


def _write_summary(rows: list[dict], group_key: str, path: Path):
    """Group rows by (group_key, pair_name), compute mean/std of divergence."""
    from collections import defaultdict
    groups: dict[tuple, list] = defaultdict(list)
    for r in rows:
        key = (r[group_key], r["pair_name"])
        groups[key].append(r)

    summary = []
    for (gval, pair), members in sorted(groups.items()):
        mags = [m["diff_mag"] for m in members]
        coss = [m["cosine_sim"] for m in members if not np.isnan(m["cosine_sim"])]
        rats = [m["mag_ratio"] for m in members if not np.isnan(m["mag_ratio"])]
        summary.append({
            group_key: gval,
            "pair_name": pair,
            "n_atoms": len(members),
            "diff_mag_mean": float(np.mean(mags)),
            "diff_mag_std": float(np.std(mags)),
            "cosine_mean": float(np.mean(coss)) if coss else float("nan"),
            "cosine_std": float(np.std(coss)) if coss else float("nan"),
            "mag_ratio_mean": float(np.mean(rats)) if rats else float("nan"),
            "mag_ratio_std": float(np.std(rats)) if rats else float("nan"),
        })

    fields = list(summary[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(summary)
    print(f"  Wrote {path} ({len(summary)} rows)")
