"""Trajectory data reader — loads the H5 master file from a GROMACS scan.

    from nmr_extract import load_trajectory
    traj = load_trajectory("output/trajectory.h5")

    traj.positions          # (T, N, 3) float64 — per-frame xyz
    traj.frame_times        # (T,) float64 — time in ps
    traj.n_atoms            # int
    traj.n_frames           # int

    # Rollup: per-atom Welford statistics across frames
    traj.rollup.bs_T0.mean  # (N,) ring current isotropic mean
    traj.rollup.bs_T0.std   # (N,) ring current isotropic std
    traj.rollup.names       # list of 43 column names

    # Bonds: per-bond length statistics
    traj.bonds.atom_a       # (B,) atom indices
    traj.bonds.atom_b       # (B,) atom indices
    traj.bonds.length_mean  # (B,) mean bond length (A)
    traj.bonds.length_std   # (B,) bond length std
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np


@dataclass(frozen=True)
class WelfordColumn:
    """Mean and std for one accumulated quantity across all atoms."""
    name: str
    mean: np.ndarray    # (N,)
    std: np.ndarray     # (N,)


class TrajectoryRollup:
    """Named access to per-atom Welford statistics.

    Access by attribute: rollup.bs_T0.mean, rollup.aimnet2_charge.std
    Access by index: rollup[0].mean (same as rollup.bs_T0.mean)
    """

    def __init__(self, names: list[str],
                 means: np.ndarray, stds: np.ndarray):
        self._names = names
        self._means = means   # (N, K)
        self._stds = stds     # (N, K)
        self._by_name: dict[str, int] = {n: i for i, n in enumerate(names)}

    @property
    def names(self) -> list[str]:
        return list(self._names)

    @property
    def n_columns(self) -> int:
        return len(self._names)

    def __getattr__(self, name: str) -> WelfordColumn:
        if name.startswith("_"):
            raise AttributeError(name)
        if name not in self._by_name:
            raise AttributeError(
                f"No rollup column '{name}'. "
                f"Available: {', '.join(self._names[:10])}...")
        k = self._by_name[name]
        return WelfordColumn(
            name=name,
            mean=self._means[:, k],
            std=self._stds[:, k],
        )

    def __getitem__(self, k: int) -> WelfordColumn:
        return WelfordColumn(
            name=self._names[k],
            mean=self._means[:, k],
            std=self._stds[:, k],
        )

    def as_block(self, columns: Optional[list[str]] = None) -> np.ndarray:
        """Return (N, K) mean array for selected columns (default: all)."""
        if columns is None:
            return self._means.copy()
        idxs = [self._by_name[c] for c in columns]
        return self._means[:, idxs]


@dataclass(frozen=True)
class BondRollup:
    """Per-bond length statistics across trajectory frames."""
    atom_a: np.ndarray       # (B,) uint64
    atom_b: np.ndarray       # (B,) uint64
    length_mean: np.ndarray  # (B,) float64
    length_std: np.ndarray   # (B,) float64

    @property
    def n_bonds(self) -> int:
        return len(self.atom_a)


@dataclass
class TrajectoryData:
    """All trajectory-level data for one protein from the H5 master file.

    This is the training data for the waveform model, the analysis
    data for the ridge, and the movie data for VTK.
    """
    protein_id: str
    n_atoms: int
    n_frames: int

    # Per-frame positions: (T, N, 3) for VTK movies + waveform input
    positions: np.ndarray       # (T, N, 3) float64
    frame_times: np.ndarray     # (T,) float64

    # Per-atom rollup statistics (43 Welford columns)
    rollup: TrajectoryRollup

    # Per-bond rollup statistics
    bonds: Optional[BondRollup]


def load_trajectory(path: str | Path) -> TrajectoryData:
    """Load a trajectory H5 master file.

    The file is written by GromacsProtein::WriteH5 after a trajectory scan.
    """
    import h5py

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Trajectory H5 not found: {path}")

    with h5py.File(path, "r") as f:
        protein_id = f.attrs.get("protein_id", path.stem)
        n_atoms = int(f.attrs["n_atoms"])
        n_frames = int(f.attrs["n_frames"])

        # Positions: stored as (T*N, 3), reshape to (T, N, 3)
        pos_raw = f["positions"][:]          # (T*N, 3)
        T_stored = int(f.attrs.get("positions_shape_T", n_frames))
        N_stored = int(f.attrs.get("positions_shape_N", n_atoms))
        positions = pos_raw.reshape(T_stored, N_stored, 3)

        frame_times = f["frame_times"][:]    # (T,)

        # Rollup
        means = f["rollup/mean"][:]          # (N, K)
        stds = f["rollup/std"][:]            # (N, K)
        # HDF5 stores strings as bytes — decode
        names_raw = f["rollup/names"][:]
        names = [n.decode() if isinstance(n, bytes) else n for n in names_raw]

        # Bonds (optional)
        bonds = None
        if "bonds" in f:
            bonds = BondRollup(
                atom_a=f["bonds/atom_a"][:],
                atom_b=f["bonds/atom_b"][:],
                length_mean=f["bonds/length_mean"][:],
                length_std=f["bonds/length_std"][:],
            )

    return TrajectoryData(
        protein_id=protein_id if isinstance(protein_id, str) else protein_id.decode(),
        n_atoms=n_atoms,
        n_frames=n_frames,
        positions=positions,
        frame_times=frame_times,
        rollup=TrajectoryRollup(names, means, stds),
        bonds=bonds,
    )
