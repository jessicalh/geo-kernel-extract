#!/usr/bin/env python3
"""Export per-atom tensor data for R plotting.

Produces CSVs with the raw geometric quantities R needs to draw
every figure in the 28 realities document.
"""

from __future__ import annotations
from pathlib import Path
import csv
import numpy as np

from secondary.loader import iter_proteins, setup_sdk, ELEMENT_NAMES
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout, assemble_kernels


def _cos(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < 1e-15 or nb < 1e-15:
        return np.nan
    return float(np.dot(a, b) / (na * nb))


def run(cfg: Config, max_proteins: int = 0):
    out = Path("output/actual_physics")
    out.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)
    efg_aro = layout.names.index("EFG_aro")
    mefg_aro = layout.names.index("MopacEFG_aro")
    bs_ring0 = layout.names.index("BS_ring0")

    # Per-type BS/HM indices
    bs_type = [i for i, n in enumerate(layout.names)
               if n.startswith("BS_") and "ring" not in n]
    hm_type = [i for i, n in enumerate(layout.names)
               if n.startswith("HM_") and "ring" not in n and "HM_H" not in n]

    f_atom = open(out / "atom_tensor_data.csv", "w", newline="")
    w_atom = csv.writer(f_atom)
    w_atom.writerow([
        "protein", "element", "dist", "theta", "rho", "z",
        "bs_mag", "hm_mag", "pq_mag",
        "bs_hm_cos", "bs_efg_cos",
        "efg_target_cos", "mefg_target_cos",
        "target_mag",
        "bs_type_r2_contrib",  # dot(bs_types, target) proxy
    ])

    # Forward selection data already exported by element_physics.py.
    # Ring-type self-fit: export from ablation data if available.

    print("Exporting per-atom tensor data...")
    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx
        kernels = assemble_kernels(p, idx, layout)
        target = p.delta.shielding.T2[idx]
        M = len(idx)

        rc = p.ring_contributions
        for j in range(M):
            el = int(rec.element[j])
            el_name = ELEMENT_NAMES.get(el, str(el))

            # Ring geometry
            dist = theta = rho = z_h = np.nan
            bs_raw_mag = hm_raw_mag = pq_raw_mag = np.nan
            bs_hm_cos = np.nan

            if rc is not None and rc.n_pairs > 0:
                atom_rc = rc.for_atom(idx[j])
                if atom_rc.n_pairs > 0:
                    nearest = np.argmin(atom_rc.distance)
                    dist = float(atom_rc.distance[nearest])
                    theta = float(atom_rc.theta[nearest])
                    rho = float(atom_rc.rho[nearest])
                    z_h = float(atom_rc.z[nearest])

                    bs_t2 = atom_rc.bs.T2[nearest]
                    hm_t2 = atom_rc.hm.T2[nearest]
                    pq_t2 = atom_rc.pq.T2[nearest]

                    bs_raw_mag = float(np.linalg.norm(bs_t2))
                    hm_raw_mag = float(np.linalg.norm(hm_t2))
                    pq_raw_mag = float(np.linalg.norm(pq_t2))
                    bs_hm_cos = _cos(bs_t2, hm_t2)

            # Kernel-target cosines (from assembled, possibly normalized)
            efg_cos = _cos(kernels[j, efg_aro, :], target[j])
            mefg_cos = _cos(kernels[j, mefg_aro, :], target[j])

            # BS-EFG independence
            bs_efg_c = _cos(kernels[j, bs_ring0, :], kernels[j, efg_aro, :])

            target_mag = float(np.linalg.norm(target[j]))

            # BS type contribution proxy: sum of |cos(bs_k, target)|
            bs_contrib = 0.0
            for ki in bs_type:
                c = _cos(kernels[j, ki, :], target[j])
                if not np.isnan(c):
                    bs_contrib += abs(c)

            w_atom.writerow([
                rec.pid, el_name,
                f"{dist:.3f}" if np.isfinite(dist) else "",
                f"{theta:.4f}" if np.isfinite(theta) else "",
                f"{rho:.3f}" if np.isfinite(rho) else "",
                f"{z_h:.3f}" if np.isfinite(z_h) else "",
                f"{bs_raw_mag:.6f}" if np.isfinite(bs_raw_mag) else "",
                f"{hm_raw_mag:.6f}" if np.isfinite(hm_raw_mag) else "",
                f"{pq_raw_mag:.6f}" if np.isfinite(pq_raw_mag) else "",
                f"{bs_hm_cos:.5f}" if np.isfinite(bs_hm_cos) else "",
                f"{bs_efg_c:.5f}" if np.isfinite(bs_efg_c) else "",
                f"{efg_cos:.5f}" if np.isfinite(efg_cos) else "",
                f"{mefg_cos:.5f}" if np.isfinite(mefg_cos) else "",
                f"{target_mag:.6f}",
                f"{bs_contrib:.4f}",
            ])

        if n_loaded % 50 == 0:
            print(f"  {n_loaded} proteins")

    f_atom.close()
    print(f"Exported {n_loaded} proteins to {out}/atom_tensor_data.csv")
