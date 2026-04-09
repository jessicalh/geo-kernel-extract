#!/usr/bin/env python3
"""Probe MLP hidden layer dimensionality.

Loads a trained model, runs the scalar→hidden mapping on all atoms,
captures hidden activations, and computes eigenspectrum.  Answers:
how many of the hidden dimensions does the model actually use?

Also reports per-kernel weight statistics: which kernels get weight,
which are suppressed, which flip sign across atoms.

Usage:
    cd learn/src
    python -m secondary hidden_probe --config calibration.toml --run ../runs/azimuthal_scalars_only_50
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch

from .loader import setup_sdk
from mutation_set.config import load_config, Config
from mutation_set.model import ShieldingT2Model
from mutation_set.dataset import CalibrationDataset
from mutation_set.kernels import KernelLayout


def run(cfg: Config, run_dir: str, max_proteins: int = 0):
    run_path = Path(run_dir)
    header = json.loads((run_path / "header.json").read_text())

    n_scalars = header["n_scalar_features"]
    n_kernels = header["n_kernels"]
    kernel_names = header["kernel_names"]

    print(f"Model: {n_kernels} kernels, {n_scalars} scalars, "
          f"hidden={cfg.model.hidden_mix}")

    # Load model
    model = ShieldingT2Model(cfg, n_scalars, n_kernels)
    state = torch.load(run_path / "best_model.pt",
                       map_location="cpu", weights_only=True)
    model.load_state_dict(state)
    model.eval()

    # Load data matching the model's training config.
    # Override kernel layout to match what the model was trained with.
    from dataclasses import replace as dc_replace
    train_calcs = header.get("per_ring_calcs",
                              ["bs", "hm", "chi", "pq", "hm_H", "disp_chi"])
    cfg_matched = dc_replace(cfg, kernels=dc_replace(
        cfg.kernels, per_ring_calculators=train_calcs))
    layout = KernelLayout.from_config(cfg_matched)
    assert layout.n_kernels == n_kernels, \
        f"Layout {layout.n_kernels} != model {n_kernels}"

    from mutation_set.dataset import list_proteins, _build_one
    proteins = list_proteins(cfg.paths.features)
    if max_proteins > 0:
        proteins = proteins[:max_proteins]
    elif header.get("n_proteins", 0) > 0:
        proteins = proteins[:header["n_proteins"]]

    all_s, all_k = [], []
    for pid in proteins:
        result, err = _build_one(cfg.paths.features / pid, cfg_matched, layout)
        if result is not None:
            all_s.append(result["scalars"])
            all_k.append(result["kernels"])

    scalars = torch.tensor(np.vstack(all_s), dtype=torch.float32)
    kernels = torch.tensor(np.vstack(all_k), dtype=torch.float32)
    N = len(scalars)
    print(f"Loaded: {N} atoms")

    # Hook hidden layers
    hidden_outputs = {}

    def make_hook(name):
        def hook(module, inp, out):
            hidden_outputs[name] = out.detach()
        return hook

    # FullyConnectedNet layers: linear, act, linear, act, linear
    # We want the activations after each nonlinearity
    hooks = []
    for i, layer in enumerate(model.mixing.mlp):
        h = layer.register_forward_hook(make_hook(f"layer_{i}"))
        hooks.append(h)

    # Forward pass
    with torch.no_grad():
        relative = model.mixing.mlp(scalars)
        magnitude = kernels.norm(dim=-1)
        gate = magnitude / (magnitude + model.mixing.gate_threshold)
        import math
        weights = (relative * gate / math.sqrt(n_kernels)).numpy()
        relative = relative.numpy()

    for h in hooks:
        h.remove()

    # Analyse each captured layer
    print(f"\n{'=' * 60}")
    print("HIDDEN LAYER DIMENSIONALITY")
    print(f"{'=' * 60}")

    for name in sorted(hidden_outputs.keys()):
        act = hidden_outputs[name].numpy()
        if act.ndim != 2 or act.shape[1] < 2:
            continue
        _analyse(act, f"Layer {name} ({act.shape[1]} units)")

    _analyse(weights, f"Kernel weights, gated ({weights.shape[1]} kernels)")
    _analyse(relative, f"Relative weights, pre-gate ({relative.shape[1]} kernels)")

    # Per-kernel weight stats
    print(f"\n{'=' * 60}")
    print("PER-KERNEL WEIGHT STATISTICS (gated)")
    print(f"{'=' * 60}\n")
    print(f"  {'Kernel':25s}  {'mean|w|':>8s}  {'std|w|':>8s}  "
          f"{'%pos':>6s}  {'%active':>8s}")
    print(f"  {'-'*65}")

    stats = []
    for k in range(n_kernels):
        w = weights[:, k]
        stats.append((
            kernel_names[k],
            float(np.mean(np.abs(w))),
            float(np.std(np.abs(w))),
            float((w > 0).mean()),
            float((np.abs(w) > 0.001).mean()),
        ))

    stats.sort(key=lambda x: x[1], reverse=True)
    for name, ma, sa, fp, fa in stats[:25]:
        print(f"  {name:25s}  {ma:8.4f}  {sa:8.4f}  {fp:6.1%}  {fa:8.1%}")

    print(f"\n  Kernels with mean|w| > 0.01: "
          f"{sum(1 for _, ma, *_ in stats if ma > 0.01)}")
    print(f"  Kernels with mean|w| > 0.001: "
          f"{sum(1 for _, ma, *_ in stats if ma > 0.001)}")
    print(f"  Kernels effectively dead (mean|w| < 0.0001): "
          f"{sum(1 for _, ma, *_ in stats if ma < 0.0001)}")


def _analyse(mat: np.ndarray, label: str):
    """Eigenspectrum of an (N, D) matrix."""
    mat_c = mat - mat.mean(axis=0, keepdims=True)
    cov = mat_c.T @ mat_c / len(mat_c)
    eig = np.linalg.eigvalsh(cov)[::-1]
    eig = np.maximum(eig, 0)
    total = eig.sum()
    if total < 1e-12:
        print(f"\n  {label}: all zeros")
        return
    cumvar = np.cumsum(eig) / total
    dim90 = int(np.searchsorted(cumvar, 0.90)) + 1
    dim95 = int(np.searchsorted(cumvar, 0.95)) + 1
    dim99 = int(np.searchsorted(cumvar, 0.99)) + 1
    active = int((eig > eig[0] * 0.01).sum())

    print(f"\n  {label}:")
    print(f"    Total dims: {mat.shape[1]}")
    print(f"    Active (>1% of max eig): {active}")
    print(f"    90% variance: {dim90} dims")
    print(f"    95% variance: {dim95} dims")
    print(f"    99% variance: {dim99} dims")
