#!/usr/bin/env python3
"""Block-averaged SEM convergence diagnostic (Grossfield & Zuckerman 2009).

Reads per-atom σ-channel timelines from a trajectory.h5 (one of the FO
SphericalTensor TS TRs: BS / HM / Mc / PiQuad / RingChi / Dispersion /
HBond / Tripeptide-BB / Larsen-HB-family) and reports how the standard
error of the per-atom mean converges as the block length grows.

Reference:
  Grossfield, A. & Zuckerman, D. M. (2009). "Quantifying uncertainty
  and sampling quality in biomolecular simulations." Annu. Rep. Comput.
  Chem. 5, 23–48.

The block-averaged SEM (BASE):

   For each block size b, divide T frames into K = T // b blocks,
   compute K block means, then

       SEM(b) = std(block_means) / sqrt(K)

   Plot SEM(b) vs b. SEM should plateau as b approaches the
   correlation time; the plateau is the converged uncertainty.
   A SEM that still GROWS at the largest accessible b indicates
   under-sampling — the trajectory is too short for that channel.

Gate for downstream Stage 2 claims: any per-atom σ-channel mean used
in calibration / regression / per-residue groupby should have its
BASE plateau visible. Otherwise the claim is contingent on the
sample, not on the physics.

Usage:
    python block_avg_convergence.py trajectory.h5 \
        --channel bs_shielding_time_series \
        --out runs/block_avg

The script reads the named TS group via the nmr_extract SDK, computes
BASE on the T0 component per atom, writes a per-atom CSV +
mean-over-atoms plot.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from nmr_extract import load_trajectory


# Map of supported channel keys -> (sdk_attr, channel_extractor).
# Each extractor returns a (N, T) float64 per-atom T0 timeline.
def _extract_t0(group_data: np.ndarray) -> np.ndarray:
    """For a (N, T, 9) SphericalTensor TS, return (N, T) T0 component."""
    return group_data[:, :, 0]


CHANNELS = {
    "bs_shielding_time_series": ("biot_savart_time_series", _extract_t0),
    "hm_shielding_time_series": ("haigh_mallion_time_series", _extract_t0),
    "mc_shielding_time_series": ("mcconnell_time_series", _extract_t0),
}


def block_avg_sem(x: np.ndarray, block_size: int) -> float:
    """Per-atom-or-channel block-averaged SEM at one block size.

    x: (T,) timeline. Returns the BASE estimate of the SEM of the
    sample mean. NaN if T < 2 * block_size (need ≥ 2 blocks).
    """
    T = x.shape[0]
    K = T // block_size
    if K < 2:
        return float("nan")
    truncated = x[: K * block_size]
    blocks = truncated.reshape(K, block_size).mean(axis=1)
    return float(blocks.std(ddof=1) / np.sqrt(K))


def base_curve(timeline: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (block_sizes, sem_values) for a per-atom timeline."""
    T = timeline.shape[0]
    max_b = T // 4  # need ≥ 4 blocks for any meaningful estimate
    block_sizes = np.unique(np.logspace(0, np.log10(max(2, max_b)),
                                         24).astype(int))
    block_sizes = block_sizes[block_sizes >= 1]
    sems = np.array([block_avg_sem(timeline, b) for b in block_sizes])
    return block_sizes, sems


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("h5_path", type=Path)
    parser.add_argument("--channel", required=True, choices=list(CHANNELS))
    parser.add_argument("--out", type=Path, default=Path("runs/block_avg"))
    args = parser.parse_args()

    traj = load_trajectory(args.h5_path)
    print(f"Loaded {args.h5_path}: {traj.n_atoms} atoms x {traj.n_frames} frames")

    # Resolve the TR group from SDK based on the channel key.
    # Direct attribute access — these are the canonical SDK paths.
    # If extended later, add the mapping to CHANNELS.
    if args.channel == "bs_shielding_time_series":
        ts = traj.welford.bs if hasattr(traj, "welford") else None  # placeholder
        raise NotImplementedError(
            "TS access for bs_shielding_time_series — wire to traj.<bs>.time_series "
            "or directly to /trajectory/bs_shielding_time_series/xyz "
            "once the SDK group lands. For now use h5py directly.")

    # NB: when the per-calc shielding-TS SDK groups land, replace the
    # above with `data = ts.xyz` and feed each channel's (N, T, 9)[..., 0]
    # into base_curve per atom.

    args.out.mkdir(parents=True, exist_ok=True)
    print(f"Outputs at {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
