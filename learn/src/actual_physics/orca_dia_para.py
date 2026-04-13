#!/usr/bin/env python3
"""Parse orca dia/para tensors from raw nmr.out files, compute WT-ALA delta.

The extractor writes orca_diamagnetic.npy and orca_paramagnetic.npy as
absolute WT values.  This script parses both WT and ALA .out files and
computes delta_dia = WT_dia - ALA_dia, delta_para = WT_para - ALA_para
per matched atom.

Usage:
    cd learn/src
    python3 -c "
    from mutation_set.config import load_config
    import sys; cfg = load_config('calibration.toml')
    sys.path.insert(0, str(cfg.paths.sdk))
    from actual_physics.orca_dia_para import run
    run(cfg)
    "
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from secondary.loader import iter_proteins, ELEMENT_NAMES
from mutation_set.config import Config

CALIBRATION = Path("/shared/2026Thesis/nmr-shielding/calibration")


_ELEMENT_MAP = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16, "P": 15, "F": 9, "CL": 17}


def _parse_xyz(xyz_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Parse positions and element numbers from .xyz file.
    Returns (positions (N,3), elements (N,))."""
    lines = xyz_path.read_text().strip().split('\n')
    n_atoms = int(lines[0].strip())
    positions = []
    elements = []
    for line in lines[2:2 + n_atoms]:
        parts = line.split()
        elements.append(_ELEMENT_MAP.get(parts[0].upper(), 0))
        positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(positions), np.array(elements, dtype=int)


def parse_orca_dia_para(nmr_out_path: Path) -> tuple[np.ndarray, np.ndarray] | None:
    """Parse dia and para 3x3 tensors from an ORCA nmr.out file.

    Returns (dia, para) each (N, 3, 3) or None if parsing fails.
    """
    text = nmr_out_path.read_text(errors='replace')

    dia_tensors = []
    para_tensors = []

    # Split by nucleus blocks
    # Pattern: " Nucleus   0N :" then "Diamagnetic contribution..." then 3 lines of matrix
    nuc_pattern = re.compile(r'Nucleus\s+\d+\w+\s*:')
    dia_pattern = re.compile(r'Diamagnetic contribution to the shielding tensor')
    para_pattern = re.compile(r'Paramagnetic contribution to the shielding tensor')

    lines = text.split('\n')
    i = 0
    while i < len(lines):
        if nuc_pattern.search(lines[i]):
            # Find dia
            dia_mat = None
            para_mat = None
            j = i + 1
            while j < min(i + 30, len(lines)):
                if dia_pattern.search(lines[j]):
                    try:
                        rows = []
                        for k in range(1, 4):
                            vals = [float(x) for x in lines[j + k].split()]
                            if len(vals) == 3:
                                rows.append(vals)
                        if len(rows) == 3:
                            dia_mat = np.array(rows)
                    except (ValueError, IndexError):
                        pass
                if para_pattern.search(lines[j]):
                    try:
                        rows = []
                        for k in range(1, 4):
                            vals = [float(x) for x in lines[j + k].split()]
                            if len(vals) == 3:
                                rows.append(vals)
                        if len(rows) == 3:
                            para_mat = np.array(rows)
                    except (ValueError, IndexError):
                        pass
                j += 1

            if dia_mat is not None and para_mat is not None:
                dia_tensors.append(dia_mat)
                para_tensors.append(para_mat)
            i = j
        else:
            i += 1

    if not dia_tensors:
        return None

    return np.array(dia_tensors), np.array(para_tensors)


def mat33_to_t2(M):
    """(N, 3, 3) → (N, 5) T2 in isometric normalization (matches C++ SphericalTensor::Decompose).

    T2[0] = sqrt(2) * Sxy        (m = -2)
    T2[1] = sqrt(2) * Syz        (m = -1)
    T2[2] = sqrt(3/2) * Szz      (m =  0)
    T2[3] = sqrt(2) * Sxz        (m = +1)
    T2[4] = (Sxx - Syy) / sqrt(2) (m = +2)

    where S is the traceless symmetric part of M.
    """
    S = 0.5 * (M + M.transpose(0, 2, 1))
    tr = np.trace(S, axis1=1, axis2=2)
    S_tl = S - (tr / 3.0)[:, np.newaxis, np.newaxis] * np.eye(3)

    SQRT2 = np.sqrt(2.0)
    SQRT3_2 = np.sqrt(3.0 / 2.0)

    Sxx = S_tl[:, 0, 0]
    Syy = S_tl[:, 1, 1]
    Szz = S_tl[:, 2, 2]
    Sxy = S_tl[:, 0, 1]
    Sxz = S_tl[:, 0, 2]
    Syz = S_tl[:, 1, 2]

    return np.column_stack([
        SQRT2 * Sxy,              # m = -2
        SQRT2 * Syz,              # m = -1
        SQRT3_2 * Szz,            # m =  0
        SQRT2 * Sxz,              # m = +1
        (Sxx - Syy) / SQRT2,      # m = +2
    ])


def cos_sim(a, b):
    na = np.linalg.norm(a, axis=1)
    nb = np.linalg.norm(b, axis=1)
    mask = (na > 1e-15) & (nb > 1e-15)
    if mask.sum() < 10:
        return float('nan'), 0
    dots = np.sum(a[mask] * b[mask], axis=1) / (na[mask] * nb[mask])
    return float(np.mean(np.abs(dots))), int(mask.sum())


def ridge_r2(X, y, lam=1e-2):
    if X.shape[0] < X.shape[1] + 10:
        return float('nan')
    XtX = X.T @ X + lam * np.eye(X.shape[1])
    w = np.linalg.solve(XtX, X.T @ y)
    p = X @ w
    ss_r = np.sum((y - p) ** 2)
    ss_t = np.sum((y - y.mean(0)) ** 2)
    return 1 - ss_r / ss_t if ss_t > 1e-12 else 0


def run(cfg: Config, max_proteins: int = 0):
    out = Path("output/actual_physics/dia_para")
    out.mkdir(parents=True, exist_ok=True)

    features = Path(cfg.paths.features)
    ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O"}

    data = {e: {"delta_dia_t2": [], "delta_para_t2": [],
                "delta_total_t2": [], "bs_t2": [], "efg_aro_t2": []}
            for e in ELEMENTS}

    import nmr_extract

    proteins = sorted(d.name for d in features.iterdir()
                      if d.is_dir() and (d / "pos.npy").exists()
                      and (d / "ring_contributions.npy").exists())
    if max_proteins > 0:
        proteins = proteins[:max_proteins]

    n_loaded = 0
    n_parsed = 0
    n_failed = 0

    for pid in proteins:
        wt_out = CALIBRATION / pid / f"{pid}_WT_nmr.out"
        ala_out = CALIBRATION / pid / f"{pid}_ALA_nmr.out"

        if not wt_out.exists() or not ala_out.exists():
            continue

        wt_result = parse_orca_dia_para(wt_out)
        ala_result = parse_orca_dia_para(ala_out)
        if wt_result is None or ala_result is None:
            n_failed += 1
            continue

        wt_dia, wt_para = wt_result
        ala_dia, ala_para = ala_result

        try:
            p = nmr_extract.load(features / pid)
        except Exception:
            continue

        if p.delta is None or p.delta.shielding is None:
            continue

        idx = np.where(p.delta.scalars.matched_mask)[0]
        if len(idx) < 10:
            continue

        # Match WT and ALA atoms by 3D position (same method as C++).
        # Parse ALA positions from .xyz file.
        n_wt = wt_dia.shape[0]
        n_ala = ala_dia.shape[0]
        n_atoms = p.element.shape[0]

        if n_wt != n_atoms:
            n_failed += 1
            continue

        # Read positions from xyz files
        ala_xyz = CALIBRATION / pid / f"{pid}_ALA.xyz"
        if not ala_xyz.exists():
            n_failed += 1
            continue

        try:
            ala_pos, ala_elem = _parse_xyz(ala_xyz)
        except Exception:
            n_failed += 1
            continue

        if len(ala_pos) != n_ala:
            n_failed += 1
            continue

        wt_pos = p.pos.data  # (n_wt, 3)
        wt_elem = p.element   # (n_wt,)

        # For each matched WT atom, find nearest ALA atom of SAME ELEMENT.
        delta_9 = np.load(features / pid / "delta_shielding.npy")[idx]  # (M, 9)
        delta_33 = delta_9.reshape(-1, 3, 3)

        ala_indices = np.full(len(idx), -1, dtype=int)
        used = set()
        for j in range(len(idx)):
            wt_j = idx[j]
            e = int(wt_elem[wt_j])
            # Only consider ALA atoms of same element
            same_elem = np.where(ala_elem == e)[0]
            if len(same_elem) == 0:
                continue
            dists = np.linalg.norm(ala_pos[same_elem] - wt_pos[wt_j], axis=1)
            # Mask used
            mask = np.ones(len(same_elem), dtype=bool)
            for k, ai in enumerate(same_elem):
                if ai in used:
                    mask[k] = False
                    dists[k] = 1e30
            best_local = int(np.argmin(dists))
            best = same_elem[best_local]
            if dists[best_local] < 0.1:  # 0.1 Å — mechanical mutant, positions identical
                ala_indices[j] = best
                used.add(best)

        valid = ala_indices >= 0
        if valid.sum() < 10:
            n_failed += 1
            continue

        valid_idx = np.where(valid)[0]
        ala_j = ala_indices[valid_idx]

        delta_dia_33 = wt_dia[idx[valid_idx]] - ala_dia[ala_j]
        delta_para_33 = wt_para[idx[valid_idx]] - ala_para[ala_j]

        delta_dia_t2 = mat33_to_t2(delta_dia_33)
        delta_para_t2 = mat33_to_t2(delta_para_33)

        # The C++ delta is delta_shielding = WT_total - ALA_total (full 3x3).
        # Our delta_dia + delta_para should reproduce this exactly if matching is right.
        # Compare against the stored delta as full 3x3, not just T2.
        recon_33 = delta_dia_33 + delta_para_33
        recon_err = np.abs(recon_33 - delta_33[valid_idx]).max()

        delta_target_t2 = p.delta.shielding.T2[idx[valid_idx]]

        elem = p.element[idx[valid_idx]]
        bs_t2 = p.biot_savart.shielding.T2[idx[valid_idx]]
        efg_t2 = p.coulomb.efg_aromatic.T2[idx[valid_idx]]

        for j in range(len(valid_idx)):
            e = int(elem[j])
            if e not in data:
                continue
            data[e]["delta_dia_t2"].append(delta_dia_t2[j])
            data[e]["delta_para_t2"].append(delta_para_t2[j])
            data[e]["delta_total_t2"].append(delta_target_t2[j])
            data[e]["bs_t2"].append(bs_t2[j])
            data[e]["efg_aro_t2"].append(efg_t2[j])

        n_parsed += 1
        n_loaded += 1
        if n_loaded % 50 == 0:
            print(f"  {n_loaded} proteins (parsed {n_parsed}, failed {n_failed}, "
                  f"last recon err {recon_err:.2e})")

    print(f"\nLoaded {n_loaded}, parsed dia/para for {n_parsed}, failed {n_failed}\n")

    # ── Results ──────────────────────────────────────────────────────

    print("Atom counts:")
    for e, en in ELEMENTS.items():
        print(f"  {en}: {len(data[e]['delta_dia_t2']):,d}")

    # Verification: dia + para = total
    print("\nVerification: |delta_dia + delta_para - delta_total| max per element:")
    for e, en in ELEMENTS.items():
        d = data[e]
        if len(d["delta_dia_t2"]) < 10:
            continue
        dd = np.array(d["delta_dia_t2"])
        dp = np.array(d["delta_para_t2"])
        dt = np.array(d["delta_total_t2"])
        err = np.abs(dd + dp - dt).max()
        print(f"  {en}: {err:.2e}")

    print("\n" + "=" * 80)
    print("  DIA/PARA DELTA: T2 magnitude and variability per element")
    print("=" * 80)
    print(f"  {'elem':>5s}  {'n':>7s}  {'|Δdia|':>10s}  {'|Δpara|':>10s}  {'|Δtot|':>10s}  "
          f"{'dia_std':>10s}  {'para_std':>10s}  {'para/dia':>10s}")
    for e, en in ELEMENTS.items():
        d = data[e]
        n = len(d["delta_dia_t2"])
        if n < 100:
            continue
        dd_mag = np.linalg.norm(d["delta_dia_t2"], axis=1)
        dp_mag = np.linalg.norm(d["delta_para_t2"], axis=1)
        dt_mag = np.linalg.norm(d["delta_total_t2"], axis=1)
        print(f"  {en:>5s}  {n:7d}  {dd_mag.mean():10.4f}  {dp_mag.mean():10.4f}  "
              f"{dt_mag.mean():10.4f}  {dd_mag.std():10.4f}  {dp_mag.std():10.4f}  "
              f"{dp_mag.std()/dd_mag.std() if dd_mag.std() > 1e-10 else 0:10.2f}")

    print("\n" + "=" * 80)
    print("  KERNEL → DELTA DIA/PARA: cosine alignment")
    print("=" * 80)
    print(f"  {'elem':>5s}  {'BS→Δdia':>9s}  {'BS→Δpara':>9s}  {'EFG→Δdia':>9s}  {'EFG→Δpara':>9s}  {'BS→Δtot':>9s}  {'EFG→Δtot':>9s}")
    for e, en in ELEMENTS.items():
        d = data[e]
        if len(d["delta_dia_t2"]) < 100:
            continue
        bs = np.array(d["bs_t2"])
        efg = np.array(d["efg_aro_t2"])
        dd = np.array(d["delta_dia_t2"])
        dp = np.array(d["delta_para_t2"])
        dt = np.array(d["delta_total_t2"])
        print(f"  {en:>5s}  {cos_sim(bs,dd)[0]:9.4f}  {cos_sim(bs,dp)[0]:9.4f}  "
              f"{cos_sim(efg,dd)[0]:9.4f}  {cos_sim(efg,dp)[0]:9.4f}  "
              f"{cos_sim(bs,dt)[0]:9.4f}  {cos_sim(efg,dt)[0]:9.4f}")

    print("\n" + "=" * 80)
    print("  KERNEL → DELTA DIA/PARA: ridge R²")
    print("=" * 80)
    print(f"  {'elem':>5s}  {'BS→Δdia':>9s}  {'BS→Δpara':>9s}  {'EFG→Δdia':>9s}  {'EFG→Δpara':>9s}  {'BS→Δtot':>9s}  {'EFG→Δtot':>9s}")
    for e, en in ELEMENTS.items():
        d = data[e]
        if len(d["delta_dia_t2"]) < 200:
            continue
        bs = np.array(d["bs_t2"])
        efg = np.array(d["efg_aro_t2"])
        dd = np.array(d["delta_dia_t2"])
        dp = np.array(d["delta_para_t2"])
        dt = np.array(d["delta_total_t2"])
        print(f"  {en:>5s}  {ridge_r2(bs,dd):9.4f}  {ridge_r2(bs,dp):9.4f}  "
              f"{ridge_r2(efg,dd):9.4f}  {ridge_r2(efg,dp):9.4f}  "
              f"{ridge_r2(bs,dt):9.4f}  {ridge_r2(efg,dt):9.4f}")

    # Save results
    results = {}
    for e, en in ELEMENTS.items():
        d = data[e]
        n = len(d["delta_dia_t2"])
        if n < 100:
            continue
        results[en] = {"n_atoms": n}

    with open(out / "dia_para_delta.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults in {out}/")
