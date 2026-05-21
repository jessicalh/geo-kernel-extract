#!/usr/bin/env python3
"""Essential dynamics (SVD / PCA) on per-frame σ-tensor stack.

Stages per-frame shielding tensors across calculators into a single
(T, N, K) array (T frames, N atoms, K=9 SphericalTensor channels per
calc × ~6 calcs ≈ 54 channels at full strength) and computes the SVD
to test whether the conformational ensemble's σ-output is low-rank.

If the leading 5-10 singular values explain > 90% of the variance,
the conformational σ-manifold is effectively low-dimensional and the
calibration regression's degrees of freedom are over-counted. If the
spectrum is flat, the ensemble explores genuinely high-dimensional
σ-space and full per-channel fits are warranted.

Per the 13-TR plan: this is one of the gates for low-rank ansatz
ML model framings. Empirical answer determines whether to ship the
"low-rank σ-manifold" claim in the thesis chapter.

Method:
  1. Reshape per-frame σ tensors to a (T, N*K) matrix X.
  2. Center per-column (per-(atom, channel) column gets its time-mean
     subtracted, so SVD captures FLUCTUATIONS not means).
  3. SVD(X) -> eigenvalues λ_i. Cumulative variance = sum(λ_i²) /
     sum(λ_total²).

Usage:
    python sigma_essential_dynamics.py trajectory.h5 \
        --calcs bs,hm,mc,piquad,ringchi,disp,hbond \
        --out runs/sigma_ed

Output:
    eigenvalues.csv               — index, λ_i, cumulative_variance_fraction
    leading_eigenvector_top.csv   — top-10 (atom, channel) contributors per
                                     leading singular vector (k=1..5)
    spectrum.pdf                  — log-log scree plot
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
    parser.add_argument("--calcs", default="bs,hm,mc,piquad,ringchi,disp",
                         help="Comma-separated list of shielding-TS prefixes "
                              "to stack. Each becomes K=9 channels per atom.")
    parser.add_argument("--out", type=Path, default=Path("runs/sigma_ed"))
    args = parser.parse_args()

    traj = load_trajectory(args.h5_path)
    print(f"Loaded {args.h5_path}: {traj.n_atoms} atoms x {traj.n_frames} frames")

    # Each shielding TS exists as `/trajectory/<name>_shielding_time_series/`
    # with shape (N, T, 9). Currently the SDK doesn't expose a typed
    # group per shielding TS; the cleanest read is via h5py direct
    # access for the stacking. SDK extension (per-calc shielding TS
    # group wrapper) is a Stage-2 SDK followup.
    print("NOTE: this script reads the per-calc shielding TS groups via "
          "h5py directly until the SDK exposes per-group wrappers.")
    print("  Calcs to stack: " + args.calcs)

    args.out.mkdir(parents=True, exist_ok=True)

    # Sketch of compute (to be wired against the real h5 once stage-2
    # output is on disk):
    #
    # import h5py
    # with h5py.File(args.h5_path, "r") as f:
    #     blocks = []
    #     for calc in args.calcs.split(","):
    #         path = f"/trajectory/{calc}_shielding_time_series"
    #         if path in f:
    #             # (N, T, 9) -> (T, N*9)
    #             g = f[path]
    #             xyz = g["xyz"][:]
    #             N, T, K = xyz.shape
    #             blocks.append(xyz.transpose(1, 0, 2).reshape(T, N * K))
    # X = np.concatenate(blocks, axis=1)
    # X_centered = X - X.mean(axis=0, keepdims=True)
    # U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    # cum_var = np.cumsum(S**2) / np.sum(S**2)

    # Placeholder output so the directory layout is visible.
    eigen_path = args.out / "eigenvalues.csv"
    with open(eigen_path, "w") as f:
        f.write("# Stage-2 essential dynamics scaffold.\n"
                "# Wire to the per-calc shielding TS groups once "
                "trajectory.h5 is on disk.\n"
                "index,singular_value,cumulative_variance_fraction\n")
    print(f"Wrote {eigen_path} (scaffold).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
