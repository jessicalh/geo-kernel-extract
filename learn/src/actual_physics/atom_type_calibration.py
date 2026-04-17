#!/usr/bin/env python3
"""Per-atom-type ridge calibration.

Stratifies C, N, O by AMBER atom name from the prmtop.  Same ridge
as clean_calibration and physics_calibration but within atom-type
sub-groups.

Carbon:  CA (alpha), C (carbonyl), CB (beta), C_side (CG/CD/CE/CZ/...)
Nitrogen: N (backbone amide), N_side (ND/NE/NZ/NH/...)
Oxygen:  O (backbone carbonyl), O_side (OG/OH/OD/OE/...)
Hydrogen: kept as one group (already R²=0.92).

Usage:
    cd learn/src
    python3 -c "
    from mutation_set.config import load_config
    import sys; cfg = load_config('calibration.toml')
    sys.path.insert(0, str(cfg.paths.sdk))
    from actual_physics.atom_type_calibration import run
    run(cfg)
    "
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from secondary.loader import iter_proteins, setup_sdk, RING_TYPE_NAMES
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels, normalize_kernels


# ── prmtop parser ───────────────────────────────────────────────────

def parse_atom_names(prmtop_path: Path) -> list[str]:
    """Read ATOM_NAME from an AMBER prmtop file."""
    names = []
    in_section = False
    with open(prmtop_path) as f:
        for line in f:
            if line.startswith('%FLAG ATOM_NAME'):
                in_section = True
                continue
            if in_section and line.startswith('%FORMAT'):
                continue
            if in_section:
                if line.startswith('%FLAG'):
                    break
                # 20a4 format: 20 names, each 4 chars
                for i in range(0, len(line.rstrip('\n')), 4):
                    name = line[i:i+4].strip()
                    if name:
                        names.append(name)
    return names


# ── atom type classification ────────────────────────────────────────

CARBON_TYPES = {
    'CA': 'CA',
    'C':  'C=O',
    'CB': 'CB',
}

def classify_carbon(atom_name: str) -> str:
    return CARBON_TYPES.get(atom_name, 'C_side')

def classify_nitrogen(atom_name: str) -> str:
    return 'N_bb' if atom_name == 'N' else 'N_side'

def classify_oxygen(atom_name: str) -> str:
    return 'O_bb' if atom_name == 'O' else 'O_side'


# ── ridge ───────────────────────────────────────────────────────────

def _ridge(X, y, lam=1.0):
    XtX = X.T @ X + lam * np.eye(X.shape[1])
    w = np.linalg.solve(XtX, X.T @ y)
    pred = X @ w
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean(axis=0)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
    return r2


# ── main ────────────────────────────────────────────────────────────

def run(cfg: Config, max_proteins: int = 0):
    out = Path("output/actual_physics/atom_type_calibration")
    out.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    core = list(range(layout.efg_end))
    n_core = len(core)
    lam = getattr(cfg.secondary, 'ridge_lambda', 1.0)

    # Raw calibration directory (prmtops)
    calibration_dir = Path(cfg.paths.features).parent.parent

    # Accumulators
    accum = {}

    def ensure(label):
        if label not in accum:
            accum[label] = {"K_raw": [], "K_norm": [], "Y": [], "scales": [],
                            "mut_type": []}

    print(f"Core kernels: {n_core}")
    print(f"Calibration dir: {calibration_dir}")
    print("Loading proteins...")

    n_loaded = 0
    n_no_prmtop = 0

    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        pid = rec.pid
        p = rec.protein
        idx = rec.matched_idx
        M = len(idx)

        # Parse prmtop for atom names
        prmtop = calibration_dir / pid / f"{pid}_WT.prmtop"
        if not prmtop.exists():
            n_no_prmtop += 1
            continue

        atom_names = parse_atom_names(prmtop)
        if len(atom_names) < max(idx) + 1:
            continue

        n_loaded += 1

        # Assemble kernels
        kernels_full = assemble_kernels(p, idx, layout)
        K_raw = kernels_full[:, core, :]
        kernels_norm, scales = normalize_kernels(
            kernels_full.copy(), cfg.normalization.kernel_std_floor)
        K_norm = kernels_norm[:, core, :]
        core_scales = scales[core]

        target = p.delta.shielding.T2[idx]

        # Mutation type: dominant ring type in this protein
        rg = p.ring_geometry
        if rg is not None and rg.n_rings > 0:
            type_counts = np.bincount(rg.ring_type.astype(int),
                                       minlength=len(RING_TYPE_NAMES))
            mut_type = int(np.argmax(type_counts))
        else:
            mut_type = -1

        for j in range(M):
            ai = int(idx[j])
            e = int(rec.element[j])
            aname = atom_names[ai]

            if e == 1:
                label = 'H'
            elif e == 6:
                label = f'C_{classify_carbon(aname)}'
            elif e == 7:
                label = f'N_{classify_nitrogen(aname)}'
            elif e == 8:
                label = f'O_{classify_oxygen(aname)}'
            else:
                continue

            ensure(label)
            accum[label]["K_raw"].append(K_raw[j])
            accum[label]["K_norm"].append(K_norm[j])
            accum[label]["Y"].append(target[j])
            accum[label]["scales"].append(core_scales)
            accum[label]["mut_type"].append(mut_type)

        if n_loaded % 50 == 0:
            print(f"  {n_loaded} proteins")

    print(f"Loaded {n_loaded} proteins ({n_no_prmtop} skipped, no prmtop)")

    # ── Ridge per atom type ─────────────────────────────────────────

    # Ordering for display
    order = ['H',
             'C_CA', 'C_C=O', 'C_CB', 'C_C_side',
             'N_N_bb', 'N_N_side',
             'O_O_bb', 'O_O_side']

    results = {}

    print(f"\n{'='*70}")
    print(f"  Atom-Type Ridge: raw / norm / +scales / +mut / fair")
    print(f"{'='*70}")

    def run_ridge_suite(label, K_raw, K_norm, Y, scales_arr, mut_types):
        N = K_raw.shape[0]

        # Raw ridge
        X_raw = K_raw.reshape(N, -1)
        r2_raw = _ridge(X_raw, Y, lam)

        # Normalised base
        X_norm = K_norm.reshape(N, -1)
        r2_norm = _ridge(X_norm, Y, lam)

        # + scales
        X_scales = np.zeros((N, n_core * 5))
        for ki in range(n_core):
            X_scales[:, ki*5:(ki+1)*5] = K_norm[:, ki, :] * scales_arr[:, ki:ki+1]
        r2_scales = _ridge(np.hstack([X_norm, X_scales]), Y, lam)

        # + mutation type (one-hot interaction)
        unique_types = sorted(set(int(t) for t in mut_types if t >= 0))
        n_types = len(unique_types)
        type_map = {t: i for i, t in enumerate(unique_types)}
        X_mut = np.zeros((N, n_types * n_core * 5))
        for j in range(N):
            mt = int(mut_types[j])
            if mt < 0:
                continue
            ti = type_map.get(mt, -1)
            if ti < 0:
                continue
            offset = ti * n_core * 5
            X_mut[j, offset:offset + n_core * 5] = X_norm[j]
        r2_mut = _ridge(np.hstack([X_norm, X_mut]), Y, lam)

        # Fair set: base + scales + mut
        r2_fair = _ridge(np.hstack([X_norm, X_scales, X_mut]), Y, lam)

        print(f"  {label:12s}  n={N:>8,d}  raw={r2_raw:.4f}  "
              f"norm={r2_norm:.4f}  +scl={r2_scales:.4f}  "
              f"+mut={r2_mut:.4f}  fair={r2_fair:.4f}")

        results[label] = {
            "n_atoms": N,
            "r2_raw": round(r2_raw, 4),
            "base": round(r2_norm, 4),
            "plus_scales": round(r2_scales, 4),
            "plus_mut": round(r2_mut, 4),
            "fair_set": round(r2_fair, 4),
        }

    for label in order:
        if label not in accum or len(accum[label]["K_raw"]) < 50:
            continue
        d = accum[label]
        run_ridge_suite(
            label,
            np.array(d["K_raw"]),
            np.array(d["K_norm"]),
            np.array(d["Y"]),
            np.array(d["scales"]),
            np.array(d["mut_type"]),
        )

    # Element totals for comparison
    for elem_prefix, elem_name in [('C_', 'C_all'), ('N_', 'N_all'), ('O_', 'O_all')]:
        K_raw_all, K_norm_all, Y_all, S_all, M_all = [], [], [], [], []
        for label in order:
            if label.startswith(elem_prefix) and label in accum:
                K_raw_all.extend(accum[label]["K_raw"])
                K_norm_all.extend(accum[label]["K_norm"])
                Y_all.extend(accum[label]["Y"])
                S_all.extend(accum[label]["scales"])
                M_all.extend(accum[label]["mut_type"])
        if len(K_raw_all) < 50:
            continue
        run_ridge_suite(
            elem_name,
            np.array(K_raw_all),
            np.array(K_norm_all),
            np.array(Y_all),
            np.array(S_all),
            np.array(M_all),
        )

    with open(out / "atom_type_results.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults in {out}/")
