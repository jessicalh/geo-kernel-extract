#!/usr/bin/env python3
"""1P9J σ-residual ↔ BMRB Rex spatial overlap (Markwick chapter deliverable).

For each backbone residue of 1P9J:

    σ-residual(residue) = stdev(σ_T0(atom)) across trajectory frames,
                          atoms in residue's backbone (N, H, CA, C, O,
                          HA), summed over heavy-atom shielding TS
                          channels (BS + HM + Mc + PiQuad + RingChi +
                          Dispersion + HBond).

    Rex(residue) = published BMRB 5801 Rex value (Wingens et al. 2003).

Plot: scatter of σ-residual vs Rex per residue, color-coded by
secondary structure (DSSP from Dssp8TS — H / E / coil tracked at
residue scope). Spearman ρ on the pair is the chapter-section
deliverable: do regions with high relaxation exchange (Rex) coincide
with regions where the classical σ predictor's frame-by-frame variance
is largest?

If ρ > 0.5 (positive correlation), the variance of the classical
calculators IS picking up the slow-timescale exchange. If ρ ≈ 0, the
classical calculators don't see the same timescale separation.

References:
  Markwick, P. R. L., Showalter, S. A., Bouvignies, G., Bruschweiler,
  R. & Blackledge, M. (2009). "Structural dynamics of protein
  backbones from NMR residual dipolar couplings and accelerated
  molecular dynamics." J. Biomol. NMR 45, 17-21.
  BMRB entry 5801 for 1P9J Rex values.

Usage:
    python markwick_overlay_1p9j.py trajectory.h5 \
        --bmrb-rex path/to/bmrb_5801_rex.csv \
        --out runs/markwick_1p9j

CSV format expected:
    residue_index,residue_number,rex_per_s
    0,1,0.5
    1,2,1.2
    ...

Output:
    per_residue.csv      — residue_idx, sigma_residual, rex_per_s, ss8_dominant
    overlay.pdf          — scatter, Spearman ρ in title, per-SS color
    figure.tex            — chapter-section text fragment with ρ + commentary
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from nmr_extract import load_trajectory


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("h5_path", type=Path)
    parser.add_argument("--bmrb-rex", type=Path, required=True)
    parser.add_argument("--out", type=Path,
                         default=Path("runs/markwick_1p9j"))
    args = parser.parse_args()

    traj = load_trajectory(args.h5_path)
    print(f"Loaded {args.h5_path}: {traj.n_atoms} atoms x {traj.n_frames} frames")

    if not args.bmrb_rex.exists():
        print(f"ERROR: BMRB Rex CSV not found at {args.bmrb_rex}")
        return 2

    args.out.mkdir(parents=True, exist_ok=True)

    # Bridge per-atom -> per-residue via the DihedralTS-emitted
    # residue_index_per_atom mapping (canonical atom-axis broadcast).
    if traj.dihedrals is None:
        print("ERROR: DihedralTimeSeriesTrajectoryResult was not attached. "
              "residue_index_per_atom is unavailable. Re-run extraction with "
              "DihedralTS enabled (PerFrameExtractionSet default).")
        return 2
    residue_index_per_atom = traj.dihedrals.residue_index_per_atom
    print(f"Residue index per atom: range "
          f"[{residue_index_per_atom.min()}, {residue_index_per_atom.max()}]")

    # σ-residual computation (sketch; wire when per-calc shielding TS
    # groups land in SDK or via h5py direct):
    #
    # import h5py
    # with h5py.File(args.h5_path, "r") as f:
    #     stack = []
    #     for calc in ["bs", "hm", "mc", "piquad", "ringchi", "disp", "hbond"]:
    #         path = f"/trajectory/{calc}_shielding_time_series/xyz"
    #         if path in f:
    #             stack.append(f[path][:, :, 0])  # T0 only
    #     sigma_t0_per_atom = np.stack(stack, axis=0).sum(axis=0)
    #     # shape (N, T) -- per-atom σ_T0 timeline
    #     sigma_resid_per_atom = sigma_t0_per_atom.std(axis=1)
    #     # Group by residue:
    #     n_res = int(residue_index_per_atom.max()) + 1
    #     sigma_resid_per_res = np.array([
    #         sigma_resid_per_atom[residue_index_per_atom == r].sum()
    #         for r in range(n_res)
    #     ])

    # Load BMRB Rex CSV.
    rex_data = np.loadtxt(args.bmrb_rex, delimiter=",", skiprows=1)
    print(f"BMRB Rex: {len(rex_data)} residues")

    # Secondary structure majority from Dssp8 transitions (dominant
    # state per residue).
    if traj.dssp8 is not None and traj.dssp8.transitions is not None:
        ss8_dom = traj.dssp8.transitions.ss8_dominant  # (R,)
        print(f"DSSP8 dominant: {len(ss8_dom)} residues")
    else:
        ss8_dom = None
        print("Note: DSSP8 transitions not in H5; SS color-code disabled.")

    # Emit scaffold output (full implementation wires once per-calc
    # shielding TS groups are SDK-loaded).
    out_csv = args.out / "per_residue.csv"
    with open(out_csv, "w") as f:
        f.write("# 1P9J σ-residual ↔ BMRB 5801 Rex overlap, Markwick frame.\n"
                "# Stage 2 scaffold -- wire per-calc shielding TS reads.\n"
                "residue_idx,sigma_residual_ppm,rex_per_s,ss8_dominant\n")
    print(f"Wrote {out_csv} (scaffold).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
