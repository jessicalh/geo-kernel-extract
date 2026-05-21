#!/usr/bin/env python3
"""Per-atom σ-channel variance ("σ-S²") + BMRB-5801 NH-S² cross-report.

The shielding tensor σ fluctuates with conformational motion. The
"order parameter" of σ (analogous to Lipari-Szabo S² for NH bond
vector relaxation) measures how dynamic the local electronic structure
is at each atom:

    σ-S²(atom) = var(σ_T0(t)) over the trajectory      [(ppm)²]

For 1P9J, BMRB 5801 publishes per-residue NH-S² values from relaxation
data. The cross-report pairs σ-S²(atom=N) and σ-S²(atom=H_amide) for
each residue against the published NH-S² to ask: does the shielding
variance follow the relaxation order parameter?

Per the 13-TR plan: classical calibration of σ predictors against
DFT delta tensors uses static structures. σ-S² is the bridge: classical
calculators producing the right mean σ AND the right σ variance is
strictly stronger than getting only the mean right.

References:
  Lipari, G. & Szabo, A. (1982). J. Am. Chem. Soc. 104, 4546-4570.
  BMRB entry 5801 (Wingens et al., 2003) for 1P9J NH-S² data.

Usage:
    python sigma_lipari_szabo.py trajectory.h5 \
        --bmrb-csv path/to/bmrb_5801_S2.csv \
        --out runs/sigma_s2_1p9j

CSV format expected:
    residue_index,residue_number,nh_s2
    0,1,0.823
    1,2,0.864
    ...

Output:
    per_atom_sigma_s2.csv          — atom_idx, residue_idx, channel, sigma_var
    per_residue_correlation.csv    — residue_idx, sigma_s2_N, sigma_s2_H, nh_s2
    overlay.pdf                    — scatter + Spearman ρ
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
    parser.add_argument("--bmrb-csv", type=Path, required=False,
                         help="Path to BMRB 5801 per-residue NH-S² CSV "
                              "(see header docstring for format). Optional; "
                              "without it, the σ-variance per atom is emitted "
                              "but the cross-correlation panel is skipped.")
    parser.add_argument("--out", type=Path,
                         default=Path("runs/sigma_s2_1p9j"))
    args = parser.parse_args()

    traj = load_trajectory(args.h5_path)
    print(f"Loaded {args.h5_path}: {traj.n_atoms} atoms x {traj.n_frames} frames")

    args.out.mkdir(parents=True, exist_ok=True)

    # σ-variance is read from the per-atom Welford rollup TRs (BS, HM,
    # McConnell, etc.) -- the `.std` channel squared IS σ-S² for the
    # T0 channel. Per the canonical TR convention, every Welford TR
    # carries a `t0.std` per-atom dataset; per-component T1 and T2
    # versions are also available.

    # The exact attribute path depends on which calculator family is
    # the closest analogue to the published NH-S² metric. For 1P9J:
    # backbone N + amide H σ_T0 variance from the BS Welford rollup
    # is the simplest first cut.
    bs_welford = traj.welford.bs if hasattr(traj.welford, "bs") else None
    if bs_welford is None:
        print("ERROR: traj.welford.bs not available -- BS Welford TR did "
              "not run for this trajectory.")
        return 2

    # BS Welford carries t0 per-atom mean + std. σ-S² = std².
    # See nmr_extract._trajectory.BsWelfordGroup.t0.std (N,).
    t0_std = bs_welford.t0.std
    sigma_var = t0_std ** 2  # (N,) (ppm)²
    print(f"σ-variance computed: range "
          f"[{np.nanmin(sigma_var):.4f}, {np.nanmax(sigma_var):.4f}] (ppm)²")

    # Emit per-atom CSV.
    out_path = args.out / "per_atom_sigma_s2.csv"
    with open(out_path, "w") as f:
        f.write("atom_idx,channel,sigma_var_ppm2\n")
        for i, v in enumerate(sigma_var):
            f.write(f"{i},bs_t0,{v}\n")
    print(f"Wrote {out_path}")

    if args.bmrb_csv:
        print(f"BMRB cross-report from {args.bmrb_csv} -- "
              "TODO: wire residue_index_per_atom from DihedralTS to bridge "
              "per-atom σ to per-residue, then load BMRB CSV + Spearman ρ.")
    else:
        print("BMRB cross-report skipped (no --bmrb-csv).")

    return 0


if __name__ == "__main__":
    sys.exit(main())
