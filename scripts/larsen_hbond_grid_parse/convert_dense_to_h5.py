#!/usr/bin/env python3
"""Convert dense NPZ grids to HDF5 for C++ consumption via HighFive.

Reads `data/larsen_hbond_grids/<STEM>_dense.npz` (output of
`pre_compute_dense_grids.py`) and emits `<STEM>_dense.h5` with the
same arrays, organized as HighFive-friendly datasets.

The C++ side (`src/LarsenHBondGrid.cpp`) consumes the .h5 via HighFive
(vendored in `extern/HighFive/`). NPZ is kept around for Python
inspection; .h5 is the C++-facing artefact.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


ARCHIVE_STEMS = ["NMANMA", "NMACOH", "NMACOO", "ALANMA", "ALACOH", "ALACOO"]


def convert(npz_path: Path, h5_path: Path) -> None:
    data = np.load(npz_path)
    with h5py.File(h5_path, "w") as h5:
        # Axes
        for axis_name in ("r_axis", "theta_axis", "rho_axis"):
            h5.create_dataset(axis_name, data=data[axis_name].astype(np.float64))
        # Donor/acceptor readouts (5D float32, shape (Nr, Nθ, Nρ, 3, 3))
        for key in data.files:
            if key.startswith(("donor_", "acceptor_")):
                # Store as float32 to match the source; C++ side casts to double
                # at runtime if needed (Eigen Mat3 is double precision).
                h5.create_dataset(
                    key, data=data[key].astype(np.float32),
                    compression="gzip", compression_opts=4,
                )
        # Stem (for self-identification)
        h5.attrs["stem"] = h5_path.stem.replace("_dense", "")
        h5.attrs["nr"] = int(data["r_axis"].shape[0])
        h5.attrs["ntheta"] = int(data["theta_axis"].shape[0])
        h5.attrs["nrho"] = int(data["rho_axis"].shape[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--grids-dir",
        default=Path("data/larsen_hbond_grids"),
        type=Path,
    )
    args = ap.parse_args()
    for stem in ARCHIVE_STEMS:
        npz = args.grids_dir / f"{stem}_dense.npz"
        h5 = args.grids_dir / f"{stem}_dense.h5"
        if not npz.exists():
            print(f"  [{stem}] SKIP: {npz} missing")
            continue
        convert(npz, h5)
        size_in = npz.stat().st_size / 1e6
        size_out = h5.stat().st_size / 1e6
        print(f"  [{stem}] {npz.name} ({size_in:.1f} MB) → {h5.name} ({size_out:.1f} MB)")


if __name__ == "__main__":
    main()
