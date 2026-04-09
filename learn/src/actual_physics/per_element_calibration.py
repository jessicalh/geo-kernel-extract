#!/usr/bin/env python3
"""Per-element calibration: train tiny models per element, extract
weight vectors and gating thresholds, compare to pooled model.

The key finding: per-element models with hidden=4 should match or
exceed the pooled hidden=64 model, because the physics IS element-
dependent and the pooled model can't learn that from scalar features.

Also extracts the weight vector from any existing trained model and
interprets it as a table of calibrated physical constants.

Outputs:
    per_element_summary.json — R² per element (ridge, model, noise floor)
    weight_vectors.csv       — per-kernel weights by element model
    gate_thresholds.csv      — gating thresholds (learned effective ranges)
    noise_floor.csv          — naive baselines per element

Usage:
    cd learn/src
    python3 -c "
    from mutation_set.config import load_config
    from actual_physics.per_element_calibration import run
    import sys; sys.path.insert(0, str(load_config('calibration.toml').paths.sdk))
    run(load_config('calibration.toml'))
    "
"""

from __future__ import annotations

import json
import csv
import math
from pathlib import Path
from dataclasses import replace

import numpy as np
import torch
import torch.nn as nn

from mutation_set.config import load_config, Config
from mutation_set.dataset import CalibrationDataset, list_proteins
from mutation_set.kernels import KernelLayout
from mutation_set.model import KernelMixingHead


def _ridge_r2(X, y, lam=1.0):
    """Ridge R² on numpy arrays. X: (N, D), y: (N, 5)."""
    XtX = X.T @ X + lam * np.eye(X.shape[1])
    w = np.linalg.solve(XtX, X.T @ y)
    pred = X @ w
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean(axis=0)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0


def _train_mixing_head(scalars, kernels, targets, gate_thresh,
                       n_kernels, hidden=4, epochs=300, lr=1e-3):
    """Train a small KernelMixingHead, return model and val R²."""
    N = scalars.shape[0]
    n_val = max(50, N // 5)
    n_train = N - n_val

    # Shuffle
    perm = torch.randperm(N)
    s_train, s_val = scalars[perm[:n_train]], scalars[perm[n_train:]]
    k_train, k_val = kernels[perm[:n_train]], kernels[perm[n_train:]]
    t_train, t_val = targets[perm[:n_train]], targets[perm[n_train:]]

    n_scalars = scalars.shape[1]
    model = KernelMixingHead(n_scalars, n_kernels, hidden)
    model.set_gate_thresholds(gate_thresh)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    s_train, k_train, t_train = s_train.to(device), k_train.to(device), t_train.to(device)
    s_val, k_val, t_val = s_val.to(device), k_val.to(device), t_val.to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)
    bs = min(512, n_train)

    best_val_loss = float("inf")
    best_state = None

    for epoch in range(epochs):
        model.train()
        idx = torch.randperm(n_train, device=device)
        for i in range(0, n_train, bs):
            batch = idx[i:i+bs]
            pred = model(s_train[batch], k_train[batch])
            loss = nn.functional.mse_loss(pred, t_train[batch])
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        scheduler.step()

        model.eval()
        with torch.no_grad():
            val_pred = model(s_val, k_val)
            val_loss = nn.functional.mse_loss(val_pred, t_val).item()
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}

    # Restore best and compute R²
    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        val_pred = model(s_val, k_val)
    ss_res = ((val_pred - t_val) ** 2).sum().item()
    ss_tot = ((t_val - t_val.mean(dim=0)) ** 2).sum().item()
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0

    return model.cpu(), r2


def _extract_mean_weights(model, scalars, kernels):
    """Run all atoms through model, capture mean gated weight vector."""
    model.eval()
    device = next(model.parameters()).device
    with torch.no_grad():
        relative = model.mlp(scalars.to(device))
        magnitude = kernels.to(device).norm(dim=-1)
        gate = magnitude / (magnitude + model.gate_threshold.to(device))
        weights = relative * gate / math.sqrt(model.n_kernels)
    return weights.cpu().numpy()


def run(cfg: Config, max_proteins: int = 0):
    out_dir = Path("output/actual_physics/calibration")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Build dataset (same as training pipeline)
    features_dir = cfg.paths.features
    proteins = list_proteins(features_dir)
    if max_proteins > 0:
        proteins = proteins[:max_proteins]

    np.random.seed(cfg.data.seed)
    np.random.shuffle(proteins)

    layout = KernelLayout.from_config(cfg)
    ds = CalibrationDataset(proteins, features_dir, cfg, layout)

    print(f"Dataset: {len(ds)} atoms, {layout.n_kernels} kernels, "
          f"{ds.scalars.shape[1]} scalars")

    # Element masks
    # Element is not directly stored in dataset — need to reconstruct
    # from the protein loading. Let's load element from the proteins.
    print("Reconstructing element labels...")
    from nmr_extract import load as nmr_load
    all_elements = []
    for pid in ds.protein_ids:
        p = nmr_load(features_dir / pid)
        idx = np.where(p.delta.scalars.matched_mask)[0]
        all_elements.append(p.element[idx])
    elements = np.concatenate(all_elements)
    assert len(elements) == len(ds), f"Element count mismatch: {len(elements)} vs {len(ds)}"

    ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O"}
    elem_masks = {e: elements == e for e in ELEMENTS}

    # ── Noise floor and ridge per element ────────────────────────────

    print("\nNoise floor and ridge per element:")
    noise_rows = []
    for e, name in ELEMENTS.items():
        m = elem_masks[e]
        n = m.sum()
        if n < 100:
            continue

        k_np = ds.kernels[m].numpy()
        t_np = ds.targets[m].numpy()

        # Naive baseline: unweighted kernel sum
        k_sum = k_np.sum(axis=1)  # (N, 5)
        ss_res = np.sum((t_np - k_sum) ** 2)
        ss_tot = np.sum((t_np - t_np.mean(axis=0)) ** 2)
        r2_sum = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0

        # Ridge
        X_flat = k_np.reshape(n, -1)
        r2_ridge = _ridge_r2(X_flat, t_np, lam=1.0)

        # Target magnitude stats
        tgt_mag = np.linalg.norm(t_np, axis=1)

        noise_rows.append({
            "element": name, "n_atoms": int(n),
            "r2_naive_sum": round(r2_sum, 4),
            "r2_ridge": round(r2_ridge, 4),
            "target_mean_mag": round(tgt_mag.mean(), 4),
            "target_std_mag": round(tgt_mag.std(), 4),
        })
        print(f"  {name}: n={n:6d}  naive={r2_sum:.4f}  ridge={r2_ridge:.4f}  "
              f"target_mag={tgt_mag.mean():.3f}+/-{tgt_mag.std():.3f}")

    # Pooled
    k_np = ds.kernels.numpy()
    t_np = ds.targets.numpy()
    r2_ridge_pooled = _ridge_r2(k_np.reshape(len(ds), -1), t_np, lam=1.0)
    print(f"  ALL: n={len(ds):6d}  ridge={r2_ridge_pooled:.4f}")

    with open(out_dir / "noise_floor.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(noise_rows[0].keys()))
        w.writeheader()
        w.writerows(noise_rows)

    # ── Per-element models (hidden=4) ────────────────────────────────

    print("\nTraining per-element models (hidden=4)...")
    gate_thresh = torch.tensor(ds.norm_stats.gate_threshold, dtype=torch.float32)

    element_models = {}
    element_r2 = {}
    element_weights = {}

    for e, name in ELEMENTS.items():
        m = elem_masks[e]
        n = int(m.sum())
        if n < 200:
            print(f"  {name}: skipping (n={n})")
            continue

        print(f"  {name} (n={n})...", end=" ", flush=True)
        model, r2 = _train_mixing_head(
            ds.scalars[m], ds.kernels[m], ds.targets[m],
            gate_thresh, layout.n_kernels, hidden=4, epochs=300)
        print(f"val R²={r2:.4f}")

        element_models[name] = model
        element_r2[name] = r2

        # Extract mean weight vector
        weights = _extract_mean_weights(model, ds.scalars[m], ds.kernels[m])
        mean_w = weights.mean(axis=0)
        element_weights[name] = mean_w

    # Pooled model (hidden=4 for fair comparison)
    print(f"  ALL (n={len(ds)})...", end=" ", flush=True)
    model_all, r2_all = _train_mixing_head(
        ds.scalars, ds.kernels, ds.targets,
        gate_thresh, layout.n_kernels, hidden=4, epochs=300)
    print(f"val R²={r2_all:.4f}")
    weights_all = _extract_mean_weights(model_all, ds.scalars, ds.kernels)
    element_weights["all"] = weights_all.mean(axis=0)
    element_r2["all"] = r2_all

    # Weighted per-element R²
    total_n = sum(elem_masks[e].sum() for e in ELEMENTS if ELEMENTS[e] in element_r2)
    weighted_r2 = sum(
        element_r2[ELEMENTS[e]] * elem_masks[e].sum() / total_n
        for e in ELEMENTS if ELEMENTS[e] in element_r2
    )

    print(f"\n{'='*60}")
    print(f"  Pooled model (hidden=4):     R² = {r2_all:.4f}")
    print(f"  Weighted per-element:        R² = {weighted_r2:.4f}")
    print(f"  Gap:                              {weighted_r2 - r2_all:.4f}")
    print(f"{'='*60}")

    # ── Weight vectors ───────────────────────────────────────────────

    print("\nExporting weight vectors...")
    with open(out_dir / "weight_vectors.csv", "w", newline="") as f:
        fieldnames = ["kernel"] + list(element_weights.keys())
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for ki, kname in enumerate(layout.names):
            row = {"kernel": kname}
            for ename, wvec in element_weights.items():
                row[ename] = f"{wvec[ki]:.6f}"
            w.writerow(row)

    # ── Gate thresholds ──────────────────────────────────────────────

    print("Exporting gate thresholds...")
    with open(out_dir / "gate_thresholds.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["kernel", "threshold", "interpretation"])
        for ki, kname in enumerate(layout.names):
            thresh = float(gate_thresh[ki])
            # Low threshold = sensitive at low magnitude = long range
            # High threshold = only active at high magnitude = short range
            interp = "long-range" if thresh < np.median(gate_thresh.numpy()) else "short-range"
            w.writerow([kname, f"{thresh:.6f}", interp])

    # ── Summary ──────────────────────────────────────────────────────

    summary = {
        "per_element_r2": {k: round(v, 4) for k, v in element_r2.items()},
        "weighted_per_element_r2": round(weighted_r2, 4),
        "pooled_r2": round(r2_all, 4),
        "gap": round(weighted_r2 - r2_all, 4),
        "ridge_per_element": {r["element"]: r["r2_ridge"] for r in noise_rows},
        "ridge_pooled": round(r2_ridge_pooled, 4),
        "n_atoms_per_element": {ELEMENTS[e]: int(elem_masks[e].sum()) for e in ELEMENTS},
    }

    with open(out_dir / "per_element_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nResults in {out_dir}/")
    return summary
