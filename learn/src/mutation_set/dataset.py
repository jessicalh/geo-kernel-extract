"""CalibrationDataset: thin orchestrator over kernels and scalars.

Loads proteins via nmr_extract SDK, assembles kernels and scalars
through their respective modules, applies normalization.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset

from nmr_extract import load

from .config import Config
from .kernels import KernelLayout, assemble_kernels, normalize_kernels
from .scalars import ScalarLayout, assemble_scalars


@dataclass
class NormStats:
    """Normalization statistics from training set, shared with val."""
    target_std: float
    scalar_mean: np.ndarray
    scalar_std: np.ndarray
    gate_threshold: np.ndarray  # per-kernel median magnitude
    scalar_layout: ScalarLayout


def list_proteins(features_dir: Path) -> list[str]:
    """Protein IDs with complete extraction (including ring_contributions)."""
    return sorted(
        d.name for d in features_dir.iterdir()
        if d.is_dir()
        and (d / "pos.npy").exists()
        and (d / "ring_contributions.npy").exists()
    )


def _build_one(protein_dir: Path, cfg: Config, layout: KernelLayout):
    """Load one protein → dict of arrays, or (None, error)."""
    try:
        p = load(protein_dir)
    except (FileNotFoundError, ValueError) as e:
        return None, str(e)

    if p.delta is None:
        return None, "no delta"

    idx = np.where(p.delta.scalars.matched_mask)[0]
    if len(idx) < cfg.data.min_matched_atoms:
        return None, f"only {len(idx)} matched atoms"

    kernels = assemble_kernels(p, idx, layout)
    kernels, kernel_scales = normalize_kernels(
        kernels, cfg.normalization.kernel_std_floor)
    scalars, scalar_layout = assemble_scalars(p, idx, kernel_scales, cfg)

    return {
        "scalars": scalars,
        "kernels": kernels,
        "target": p.delta.shielding.T2[idx],
        "ring_dist": p.delta.scalars.nearest_removed_ring_dist[idx],
        "scalar_layout": scalar_layout,
    }, None


class CalibrationDataset(Dataset):
    """WT-ALA delta T2 dataset.

    Normalization:
      - Kernels: per-protein std (angular structure preserved)
      - Scalars: z-score (train stats shared to val)
      - Targets: divide by target_std
    """

    def __init__(self, protein_ids: list[str], features_dir: Path,
                 cfg: Config, layout: KernelLayout,
                 norm_stats: NormStats | None = None):
        all_s, all_k, all_t, all_d = [], [], [], []
        self.protein_ids = []
        self.protein_boundaries = []
        scalar_layout = None
        offset = 0

        for pid in protein_ids:
            result, err = _build_one(features_dir / pid, cfg, layout)
            if result is None:
                if err:
                    print(f"  skip {pid}: {err}")
                continue
            n = len(result["target"])
            all_s.append(result["scalars"])
            all_k.append(result["kernels"])
            all_t.append(result["target"])
            all_d.append(result["ring_dist"])
            self.protein_ids.append(pid)
            self.protein_boundaries.append((offset, offset + n))
            offset += n
            if scalar_layout is None:
                scalar_layout = result["scalar_layout"]

        print(f"  loaded {len(self.protein_ids)} proteins, "
              f"skipped {len(protein_ids) - len(self.protein_ids)}")
        if not all_s:
            raise ValueError("No proteins loaded")

        scalars_np = np.vstack(all_s)
        kernels_np = np.vstack(all_k)
        targets_np = np.vstack(all_t)
        self.ring_dist = torch.tensor(np.concatenate(all_d), dtype=torch.float32)

        if norm_stats is None:
            target_std = float(targets_np.std())
            scalar_mean = scalars_np.mean(axis=0)
            scalar_std = np.maximum(scalars_np.std(axis=0),
                                    cfg.normalization.scalar_std_floor)

            # Skip z-score for categorical columns
            cat_cols = scalar_layout.categorical_columns
            scalar_mean[cat_cols] = 0.0
            scalar_std[cat_cols] = 1.0

            # Per-kernel gate threshold from training data
            kernel_mags = np.linalg.norm(kernels_np, axis=-1)
            gate_threshold = np.ones(layout.n_kernels)
            for k in range(layout.n_kernels):
                nonzero = kernel_mags[:, k][kernel_mags[:, k] > 1e-8]
                if len(nonzero) > cfg.normalization.gate_threshold_min_atoms:
                    gate_threshold[k] = float(np.median(nonzero))

            self.norm_stats = NormStats(
                target_std=target_std,
                scalar_mean=scalar_mean,
                scalar_std=scalar_std,
                gate_threshold=gate_threshold,
                scalar_layout=scalar_layout,
            )
        else:
            self.norm_stats = norm_stats

        ns = self.norm_stats
        scalars_np = (scalars_np - ns.scalar_mean) / ns.scalar_std
        targets_np = targets_np / ns.target_std

        self.scalars = torch.tensor(scalars_np, dtype=torch.float32)
        self.kernels = torch.tensor(kernels_np, dtype=torch.float32)
        self.targets = torch.tensor(targets_np, dtype=torch.float32)
        self.target_std = torch.tensor(ns.target_std, dtype=torch.float32)

    def to(self, device):
        self.scalars = self.scalars.to(device)
        self.kernels = self.kernels.to(device)
        self.targets = self.targets.to(device)
        self.target_std = self.target_std.to(device)
        self.ring_dist = self.ring_dist.to(device)
        return self

    def __len__(self):
        return len(self.scalars)

    def __getitem__(self, idx):
        return self.scalars[idx], self.kernels[idx], self.targets[idx]
