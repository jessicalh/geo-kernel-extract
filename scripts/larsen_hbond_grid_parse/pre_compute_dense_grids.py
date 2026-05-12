"""Pre-compute dense (r, theta, rho) tensor grids from scattered Larsen H-bond
DFT archives.

Input: the 6 scattered NPZ archives at
``data/larsen_hbond_grids/<STEM>_grid.npz`` produced by
``parse_larsen_hbond_grids.py``. Each holds (N,) arrays of measured (r, theta,
rho) plus (N, 3, 3) Delta-sigma tensors per readout atom on the donor and (for
NMA-acceptor archives) on the acceptor.

Output: ``data/larsen_hbond_grids/<STEM>_dense.npz`` — one regular dense grid
per archive, ready for trivial C++ trilinear lookup. Per archive:

* ``r_axis``, ``theta_axis``, ``rho_axis``  — 1D float64 axes.
* ``donor_<X>`` (Nr, Ntheta, Nrho, 3, 3) — Delta-sigma tensor per readout.
* ``acceptor_<X>`` (Nr, Ntheta, Nrho, 3, 3) — for NMA-acceptor archives only.

Pipeline per archive:

1. Load NPZ + meta JSON.
2. Derive nominal axes:
     - r:    0.125 A step (NMA donor) or 0.200 A step (ALA donor).
     - theta: 10 deg step.
     - rho:   15 deg step, periodic at +-180 deg.
   The actual scattered geometry shows a 0.0001 A or 0.0001 deg dithering plus
   small rho scatter (max ~7.5 deg from nominal). We snap each scattered point
   to its nearest nominal bin.
3. Average tensors within each bin. Most bins contain exactly one or two
   scattered points; averaging tolerates the dithering.
4. Detect "holes" — nominal bins with no scattered measurement. Fill via
   nearest-neighbour in the (r, theta, rho) regular grid (treating rho as
   periodic in distance).
5. Build cubic spline (``RegularGridInterpolator(method='cubic')``) on the
   binned grid, with rho extended by one period on each side so the spline is
   periodic in rho.
6. Evaluate the cubic spline on a denser regular grid:
     - r:    5x density (0.04 A on NMA, 0.040 A on ALA -- 1/5 of step).
     - theta: 2x density (5 deg).
     - rho:   3x density (5 deg).
7. Crop the dense rho axis back to [-180, 180) and emit.

Acceptance check at the end: for each archive, sample 50 random scattered
points, query the dense interpolator at the binned (r, theta, rho), and
report the max/median residual against the binned tensor. The dense grid
cannot be more accurate than the binning, but it should match at sample
positions to within ~1% relative.

Author: Jessica Larsen H-bond pipeline, 2026-05-11.
Only depends on numpy + scipy.
"""

from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import RegularGridInterpolator


# -------------------------------------------------------------------------
# Constants
# -------------------------------------------------------------------------

DATA_DIR = "/shared/2026Thesis/nmr-shielding/data/larsen_hbond_grids"

ARCHIVE_STEMS = ["ALANMA", "NMANMA", "ALACOH", "ALACOO", "NMACOH", "NMACOO"]

# Per-archive nominal grid steps (read directly from the actual data spacing).
# r step depends on donor kind. theta and rho are uniform across archives.
NOMINAL_R_STEP = {
    "ALANMA": 0.200,
    "ALACOH": 0.200,
    "ALACOO": 0.200,
    "NMANMA": 0.125,
    "NMACOH": 0.125,
    "NMACOO": 0.125,
}
NOMINAL_THETA_STEP = 10.0   # degrees
NOMINAL_RHO_STEP = 15.0     # degrees

# Density multipliers for the dense output grid.
DENSE_R_MULT = 5     # 5x on r
DENSE_THETA_MULT = 2  # 2x on theta
DENSE_RHO_MULT = 3   # 3x on rho


# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------


def build_nominal_axis(values: np.ndarray, step: float) -> np.ndarray:
    """Build the nominal regular axis spanning the scattered data range.

    Anchors to ``round(min(values) / step) * step`` so every scattered value
    snaps cleanly to a bin centre. The axis is dense enough to cover the full
    observed range with one bin past either side as safety.
    """
    v_min = float(values.min())
    v_max = float(values.max())
    # Snap min/max to the nearest nominal multiple of step.
    k_min = int(np.round(v_min / step))
    k_max = int(np.round(v_max / step))
    axis = np.arange(k_min, k_max + 1, dtype=np.int64).astype(np.float64) * step
    return axis


def bin_to_axis(values: np.ndarray, axis: np.ndarray) -> np.ndarray:
    """Return the index along ``axis`` (uniform spacing) nearest each value."""
    if len(axis) == 1:
        return np.zeros_like(values, dtype=np.int64)
    step = axis[1] - axis[0]
    idx = np.round((values - axis[0]) / step).astype(np.int64)
    idx = np.clip(idx, 0, len(axis) - 1)
    return idx


def fill_holes_nearest(
    grid: np.ndarray, valid: np.ndarray, r_axis: np.ndarray,
    theta_axis: np.ndarray, rho_axis: np.ndarray
) -> np.ndarray:
    """Fill missing bins by nearest valid neighbour in (r, theta, rho).

    ``grid`` has shape (Nr, Ntheta, Nrho, 3, 3); ``valid`` has shape
    (Nr, Ntheta, Nrho) and is True where a measurement exists. rho is treated
    as periodic with period 360 deg for distance computation.
    """
    filled = grid.copy()
    if valid.all():
        return filled

    nr, nth, nrh = valid.shape
    # Cartesian coordinates of every bin (scale-free; weighted equally).
    # Use index-space distance for simplicity since axes are uniform.
    rr, tt, pp = np.meshgrid(
        np.arange(nr), np.arange(nth), np.arange(nrh), indexing="ij"
    )
    valid_idx = np.argwhere(valid)            # (Nv, 3) of (ir, ith, irh)
    missing_idx = np.argwhere(~valid)          # (Nm, 3)

    if valid_idx.size == 0:
        raise RuntimeError("All bins empty -- cannot fill holes.")

    # For each missing bin, find the nearest valid bin in index space, treating
    # rho as periodic with period nrh. With nrh ~ 24 this is cheap.
    for mi, (ir, ith, irh) in enumerate(missing_idx):
        dr = valid_idx[:, 0] - ir
        dt = valid_idx[:, 1] - ith
        dp = valid_idx[:, 2] - irh
        # Periodic distance in rho.
        dp = np.minimum(np.abs(dp), nrh - np.abs(dp))
        dist2 = dr * dr + dt * dt + dp * dp
        best = int(np.argmin(dist2))
        bi = tuple(valid_idx[best])
        filled[ir, ith, irh] = grid[bi]
    return filled


def cubic_interp_periodic_rho(
    grid: np.ndarray, r_axis: np.ndarray, theta_axis: np.ndarray,
    rho_axis: np.ndarray, dense_r: np.ndarray, dense_theta: np.ndarray,
    dense_rho: np.ndarray
) -> np.ndarray:
    """Cubic-spline interpolate ``grid`` (Nr, Ntheta, Nrho, 3, 3) onto the
    dense grid, with rho treated as periodic.

    Implementation: extend grid by 2 periods on either side of rho (so the
    spline has wrap context), build interpolator on the extended axis, evaluate
    at dense points, return.
    """
    nrh = len(rho_axis)
    rho_step = rho_axis[1] - rho_axis[0]
    # Extend by 2 planes on each side; spline only needs 2 for a clean cubic.
    pad = 2
    ext_rho = np.concatenate([
        rho_axis[-pad:] - 360.0,
        rho_axis,
        rho_axis[:pad] + 360.0,
    ])
    ext_grid = np.concatenate([
        grid[:, :, -pad:, ...],
        grid,
        grid[:, :, :pad, ...],
    ], axis=2)

    # Build interpolator. RegularGridInterpolator with method='cubic' uses
    # cubic spline along each axis. Inputs are the axis vectors plus the
    # n-d value array.
    interp = RegularGridInterpolator(
        (r_axis, theta_axis, ext_rho),
        ext_grid,
        method="cubic",
        bounds_error=False,
        fill_value=None,  # extrapolate at edges
    )

    # Dense query points.
    # Reshape dense_rho into [-180, 180) range explicitly.
    rr, tt, pp = np.meshgrid(dense_r, dense_theta, dense_rho, indexing="ij")
    pts = np.stack([rr.ravel(), tt.ravel(), pp.ravel()], axis=-1)
    vals = interp(pts)
    out = vals.reshape(len(dense_r), len(dense_theta), len(dense_rho), 3, 3)
    return out


# -------------------------------------------------------------------------
# Per-archive processing
# -------------------------------------------------------------------------


@dataclass
class ArchiveSummary:
    stem: str
    n_scattered: int
    n_bins_filled: int
    n_bins_total: int
    n_holes: int
    dense_shape: tuple
    output_size_mb: float
    readout_residuals: dict  # readout_name -> (max_rel, median_rel)


def process_archive(stem: str) -> ArchiveSummary:
    """Process one scattered archive into a dense one."""
    in_npz_path = os.path.join(DATA_DIR, f"{stem}_grid.npz")
    in_meta_path = os.path.join(DATA_DIR, f"{stem}_meta.json")
    out_path = os.path.join(DATA_DIR, f"{stem}_dense.npz")

    if not os.path.isfile(in_npz_path):
        raise FileNotFoundError(in_npz_path)
    if not os.path.isfile(in_meta_path):
        raise FileNotFoundError(in_meta_path)

    npz = np.load(in_npz_path)
    with open(in_meta_path) as f:
        meta = json.load(f)

    r = npz["r"].astype(np.float64)
    theta = npz["theta"].astype(np.float64)
    rho = npz["rho"].astype(np.float64)
    n_scat = len(r)

    # 1. Build nominal axes.
    r_step = NOMINAL_R_STEP[stem]
    r_axis = build_nominal_axis(r, r_step)
    theta_axis = build_nominal_axis(theta, NOMINAL_THETA_STEP)
    # rho is periodic. Build axis covering [-180, 180) at the nominal step.
    rho_axis = np.arange(-180.0, 180.0, NOMINAL_RHO_STEP, dtype=np.float64)
    # snap rho to nearest 15 deg, wrap into [-180, 180).
    rho_snap = ((rho + 180.0) % 360.0) - 180.0

    # 2. Bin each scattered point.
    ir = bin_to_axis(r, r_axis)
    ith = bin_to_axis(theta, theta_axis)
    # For rho, periodic snap to nearest bin centre.
    irh = np.round((rho_snap - rho_axis[0]) / NOMINAL_RHO_STEP).astype(np.int64)
    irh = irh % len(rho_axis)

    nr, nth, nrh = len(r_axis), len(theta_axis), len(rho_axis)
    bin_sum = np.zeros((nr, nth, nrh, 3, 3), dtype=np.float64)
    bin_count = np.zeros((nr, nth, nrh), dtype=np.int64)

    # Identify readout atoms present in archive.
    donor_names = [k.split("_", 1)[1] for k in npz.files if k.startswith("donor_")]
    acceptor_names = [k.split("_", 1)[1] for k in npz.files if k.startswith("acceptor_")]

    # 3. Accumulate per-readout-atom sums + counts in one pass.
    # We need separate (bin_sum, bin_count) per readout, since they share the
    # same bin map. Build a dict.
    bins_per_readout: dict[str, np.ndarray] = {}  # full name -> bin_sum array
    counts = np.zeros((nr, nth, nrh), dtype=np.int64)
    # First, count once (shared across readouts since geometry is identical).
    np.add.at(counts, (ir, ith, irh), 1)

    # Iterate as (prefix, atom) pairs so donor and acceptor with the same
    # short name (e.g. both "N") are kept in separate dict entries.
    typed_readouts = [("donor", n) for n in donor_names] + \
                     [("acceptor", n) for n in acceptor_names]
    for prefix, name in typed_readouts:
        full_key = f"{prefix}_{name}"
        tensors = npz[full_key].astype(np.float64)   # (N, 3, 3)
        accum = np.zeros((nr, nth, nrh, 3, 3), dtype=np.float64)
        np.add.at(accum, (ir, ith, irh), tensors)
        bins_per_readout[full_key] = accum

    # 4. Detect holes and convert sums to means.
    valid = counts > 0
    n_holes = int((~valid).sum())
    n_bins_filled = int(valid.sum())
    n_bins_total = int(valid.size)

    # Per-readout binned-mean tensors. Holes get NaN first then filled later.
    binned_per_readout: dict[str, np.ndarray] = {}
    for full_key, accum in bins_per_readout.items():
        mean = np.where(
            valid[..., None, None],
            accum / np.maximum(counts, 1)[..., None, None],
            np.nan,
        )
        # Fill nearest-neighbour for holes (cubic spline can't tolerate NaN).
        if n_holes > 0:
            mean = fill_holes_nearest(mean, valid, r_axis, theta_axis, rho_axis)
        binned_per_readout[full_key] = mean

    # 5/6. Build dense axes and cubic-spline-interpolate.
    dense_r = np.linspace(
        r_axis[0], r_axis[-1],
        (nr - 1) * DENSE_R_MULT + 1
    )
    dense_theta = np.linspace(
        theta_axis[0], theta_axis[-1],
        (nth - 1) * DENSE_THETA_MULT + 1
    )
    # Dense rho still covers [-180, 180) exactly. Construct it explicitly.
    dense_rho_step = NOMINAL_RHO_STEP / DENSE_RHO_MULT
    dense_rho = np.arange(-180.0, 180.0, dense_rho_step, dtype=np.float64)

    dense_per_readout: dict[str, np.ndarray] = {}
    for full_key, binned in binned_per_readout.items():
        dense_per_readout[full_key] = cubic_interp_periodic_rho(
            binned, r_axis, theta_axis, rho_axis,
            dense_r, dense_theta, dense_rho,
        )

    # 7. Emit dense NPZ. We store dense_per_readout as float32 to halve disk
    # while keeping plenty of precision (binning + cubic are far above the
    # float32 ULP).
    out_dict: dict[str, np.ndarray] = {
        "r_axis": dense_r.astype(np.float64),
        "theta_axis": dense_theta.astype(np.float64),
        "rho_axis": dense_rho.astype(np.float64),
    }
    for full_key, vals in dense_per_readout.items():
        out_dict[full_key] = vals.astype(np.float32)

    np.savez_compressed(out_path, **out_dict)
    output_size_mb = os.path.getsize(out_path) / (1024.0 * 1024.0)

    # 8. Accuracy check: sample 50 scattered points, compute dense
    # interpolation at their binned bin centre, compare to the binned mean
    # tensor. Both should match modulo cubic spline smoothing.
    rng = np.random.default_rng(seed=0)
    sample_n = min(50, n_scat)
    sample_idx = rng.choice(n_scat, size=sample_n, replace=False)
    residuals: dict[str, tuple[float, float]] = {}

    for full_key, dense_vals in dense_per_readout.items():
        rel_residuals = []
        for s in sample_idx:
            ri = ir[s]
            ti = ith[s]
            phi_i = irh[s]
            bin_centre_r = r_axis[ri]
            bin_centre_t = theta_axis[ti]
            bin_centre_p = rho_axis[phi_i]
            interp = RegularGridInterpolator(
                (dense_r, dense_theta, dense_rho),
                dense_vals.astype(np.float64),
                method="linear",  # trilinear, matches future C++ lookup
                bounds_error=False,
                fill_value=None,
            )
            # Query at bin centre (which exists in the dense grid because the
            # dense grid is a superset of the binned grid).
            queried = interp(
                np.array([[bin_centre_r, bin_centre_t, bin_centre_p]])
            )[0]
            binned_val = binned_per_readout[full_key][ri, ti, phi_i]
            # Relative residual norm, with epsilon to avoid div0.
            denom = np.linalg.norm(binned_val) + 1e-9
            rel = np.linalg.norm(queried - binned_val) / denom
            rel_residuals.append(rel)
        rel_residuals = np.asarray(rel_residuals)
        residuals[full_key] = (float(rel_residuals.max()),
                                float(np.median(rel_residuals)))

    dense_shape = (len(dense_r), len(dense_theta), len(dense_rho))
    return ArchiveSummary(
        stem=stem,
        n_scattered=n_scat,
        n_bins_filled=n_bins_filled,
        n_bins_total=n_bins_total,
        n_holes=n_holes,
        dense_shape=dense_shape,
        output_size_mb=output_size_mb,
        readout_residuals=residuals,
    )


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------


def main() -> int:
    summaries: list[ArchiveSummary] = []
    total_mb = 0.0
    for stem in ARCHIVE_STEMS:
        print(f"\n=== {stem} ===", flush=True)
        try:
            summary = process_archive(stem)
        except Exception as exc:  # noqa: BLE001
            print(f"  FAIL: {exc}", flush=True)
            return 1
        summaries.append(summary)
        total_mb += summary.output_size_mb
        print(f"  scattered points: {summary.n_scattered}")
        print(f"  binned grid:      {summary.n_bins_filled}/"
              f"{summary.n_bins_total} bins filled "
              f"({summary.n_holes} holes nearest-filled)")
        print(f"  dense grid:       {summary.dense_shape[0]} x "
              f"{summary.dense_shape[1]} x {summary.dense_shape[2]} = "
              f"{np.prod(summary.dense_shape):,} points")
        print(f"  output size:      {summary.output_size_mb:.2f} MB")
        print(f"  accuracy at sample bin centres (median / max rel):")
        for full_key, (rmax, rmed) in summary.readout_residuals.items():
            print(f"    {full_key:<16}: {rmed:.2e} / {rmax:.2e}")

    print(f"\n=== Total ===")
    print(f"  Total output disk size: {total_mb:.2f} MB")
    if total_mb >= 200.0:
        print(f"  WARNING: exceeds 200 MB budget.")
    else:
        print(f"  Within 200 MB budget.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
