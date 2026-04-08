"""Evaluation metrics: R^2, baselines, per-protein breakdown.

All distance bands and thresholds come from Config — no magic numbers.
"""

from __future__ import annotations

import numpy as np
import torch

from .config import Config


def compute_r2(pred: torch.Tensor, tgt: torch.Tensor,
               ring_dist: torch.Tensor | None, cfg: Config) -> dict:
    """Per-component, overall, and distance-stratified R^2."""
    result = {}

    for i, m in enumerate([-2, -1, 0, 1, 2]):
        ss_res = torch.sum((tgt[:, i] - pred[:, i]) ** 2).item()
        ss_tot = torch.sum((tgt[:, i] - tgt[:, i].mean()) ** 2).item()
        result[f"r2_m{m:+d}"] = 1.0 - ss_res / max(ss_tot, 1e-12)
    result["r2_mean"] = np.mean([result[f"r2_m{m:+d}"] for m in [-2, -1, 0, 1, 2]])

    ss_res = torch.sum((tgt - pred) ** 2).item()
    ss_tot = torch.sum((tgt - tgt.mean(dim=0, keepdim=True)) ** 2).item()
    result["r2_flat"] = 1.0 - ss_res / max(ss_tot, 1e-12)
    result["rmse_per_atom_ppm"] = torch.sqrt(
        torch.sum((tgt - pred) ** 2, dim=1)).mean().item()

    if ring_dist is not None:
        for lo, hi in cfg.analysis.distance_bands:
            label = f"{int(lo)}-{int(hi)}A" if hi < 900 else f"{int(lo)}+A"
            mask = (ring_dist >= lo) & (ring_dist < hi)
            n = mask.sum().item()
            if n < cfg.analysis.min_atoms_per_band:
                result[f"r2_{label}"] = float("nan")
                result[f"n_{label}"] = n
                continue
            p_b, t_b = pred[mask], tgt[mask]
            ss_res = torch.sum((t_b - p_b) ** 2).item()
            ss_tot = torch.sum((t_b - t_b.mean(dim=0, keepdim=True)) ** 2).item()
            result[f"r2_{label}"] = 1.0 - ss_res / max(ss_tot, 1e-12)
            result[f"n_{label}"] = n

    return result


def compute_naive_baselines(ds, cfg: Config) -> dict:
    """Physics baselines requiring no learning."""
    tgt = ds.targets
    kernels = ds.kernels
    ss_tot = torch.sum((tgt - tgt.mean(dim=0, keepdim=True)) ** 2).item()
    result = {}

    for name, slc in [("bs_only", slice(0, 8)),
                       ("ring_sum", slice(0, 24)),
                       ("all_sum", slice(None))]:
        pred = kernels[:, slc, :].sum(dim=1)
        ss_res = torch.sum((tgt - pred) ** 2).item()
        result[f"r2_{name}"] = 1.0 - ss_res / max(ss_tot, 1e-12)

    X = kernels.reshape(len(kernels), -1)
    try:
        lam = cfg.analysis.ridge_lambda
        XtX = X.T @ X + lam * torch.eye(X.shape[1], device=X.device)
        w = torch.linalg.solve(XtX, X.T @ tgt)
        ss_res = torch.sum((tgt - X @ w) ** 2).item()
        result["r2_ridge"] = 1.0 - ss_res / max(ss_tot, 1e-12)
    except Exception:
        result["r2_ridge"] = float("nan")

    return result


def compute_per_protein_r2(pred: torch.Tensor, tgt: torch.Tensor,
                           boundaries: list[tuple[int, int]],
                           protein_ids: list[str]) -> dict[str, float]:
    """R^2 per protein — identifies outliers and bad DFT."""
    result = {}
    for pid, (lo, hi) in zip(protein_ids, boundaries):
        if hi - lo < 10:
            result[pid] = float("nan")
            continue
        p, t = pred[lo:hi], tgt[lo:hi]
        ss_res = torch.sum((t - p) ** 2).item()
        ss_tot = torch.sum((t - t.mean(dim=0, keepdim=True)) ** 2).item()
        result[pid] = 1.0 - ss_res / max(ss_tot, 1e-12)
    return result
