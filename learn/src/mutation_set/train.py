#!/usr/bin/env python3
"""Training loop for mutation-set calibration.

GPU-resident: all data on device, no DataLoader overhead.
All parameters from calibration.toml, with CLI overrides.

Usage:
    cd learn/src
    python -m mutation_set.train --config calibration.toml
    python -m mutation_set.train --config calibration.toml --max-proteins 40 --epochs 50
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

from .config import load_config
from .dataset import CalibrationDataset, NormStats, list_proteins
from .evaluate import compute_r2, compute_naive_baselines, compute_per_protein_r2
from .kernels import KernelLayout
from .model import ShieldingT2Model


def build_datasets(cfg):
    """Build train/val split from config."""
    features_dir = cfg.paths.features
    proteins = list_proteins(features_dir)
    if cfg.data.max_proteins > 0:
        proteins = proteins[:cfg.data.max_proteins]

    np.random.seed(cfg.data.seed)
    np.random.shuffle(proteins)
    n_val = max(2, int(len(proteins) * cfg.data.val_fraction))

    layout = KernelLayout.from_config(cfg)
    train_ds = CalibrationDataset(
        proteins[n_val:], features_dir, cfg, layout)
    val_ds = CalibrationDataset(
        proteins[:n_val], features_dir, cfg, layout,
        norm_stats=train_ds.norm_stats)
    return train_ds, val_ds, layout, len(proteins)


def train(cfg, run_name: str, notes: str = ""):
    runs_dir = cfg.paths.runs
    runs_dir.mkdir(parents=True, exist_ok=True)
    this_run = runs_dir / run_name
    this_run.mkdir(exist_ok=True)

    train_ds, val_ds, layout, n_proteins = build_datasets(cfg)

    header = {
        "run_name": run_name,
        "started": datetime.now().isoformat(),
        "n_proteins": n_proteins,
        "n_train_atoms": len(train_ds),
        "n_val_atoms": len(val_ds),
        "target_std_ppm": train_ds.target_std.item(),
        "n_kernels": layout.n_kernels,
        "kernel_names": layout.names,
        "n_scalar_features": train_ds.scalars.shape[1],
        "scalar_blocks": list(zip(
            train_ds.norm_stats.scalar_layout.names,
            train_ds.norm_stats.scalar_layout.widths)),
        "epochs": cfg.training.epochs,
        "batch_size": cfg.training.batch_size,
        "lr": cfg.training.lr,
        "use_correction": cfg.model.use_correction,
        "notes": notes,
    }

    print(f"Run: {run_name}")
    print(f"Train: {len(train_ds)} atoms, Val: {len(val_ds)} atoms")
    print(f"Target std: {train_ds.target_std:.4f} ppm")
    print(f"Kernels: {layout.n_kernels}, Scalars: {train_ds.scalars.shape[1]}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    train_ds.to(device)
    val_ds.to(device)
    print(f"Device: {device}")

    # Baselines
    print(f"\n{'─' * 60}")
    for name, ds in [("train", train_ds), ("val", val_ds)]:
        bl = compute_naive_baselines(ds, cfg, layout)
        print(f"  {name}: ridge={bl['r2_ridge']:.4f}  "
              f"ring_sum={bl['r2_ring_sum']:.4f}  "
              f"bs_only={bl['r2_bs_only']:.4f}")
        header[f"baseline_{name}"] = {
            k: round(v, 4) if isinstance(v, float) else v
            for k, v in bl.items()}
    print(f"{'─' * 60}\n")

    # Model
    n_scalars = train_ds.scalars.shape[1]
    model = ShieldingT2Model(cfg, n_scalars, layout.n_kernels)
    model.mixing.set_gate_thresholds(
        torch.tensor(train_ds.norm_stats.gate_threshold, dtype=torch.float32))
    model = model.to(device)

    print(f"Model: {model.parameter_count()} params")
    header["n_parameters"] = model.parameter_count()

    with open(this_run / "header.json", "w") as f:
        json.dump(header, f, indent=2)

    optimizer = torch.optim.Adam(
        model.parameters(), lr=cfg.training.lr,
        weight_decay=cfg.training.weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=cfg.training.epochs)

    log_path = this_run / "epochs.jsonl"
    best_val_loss = float("inf")
    best_epoch = 0
    N_train = len(train_ds)
    bs = cfg.training.batch_size

    for epoch in range(cfg.training.epochs):
        model.train()
        perm = torch.randperm(N_train, device=device)
        train_loss = 0.0
        for i in range(0, N_train, bs):
            idx = perm[i:i + bs]
            pred = model(train_ds.scalars[idx], train_ds.kernels[idx])
            loss = nn.functional.mse_loss(pred, train_ds.targets[idx])
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * len(idx)
        train_loss /= N_train

        model.eval()
        with torch.no_grad():
            val_pred = model(val_ds.scalars, val_ds.kernels)
            val_loss = nn.functional.mse_loss(val_pred, val_ds.targets).item()
        scheduler.step()

        is_best = val_loss < best_val_loss
        if is_best:
            best_val_loss = val_loss
            best_epoch = epoch + 1
            torch.save(model.state_dict(), this_run / "best_model.pt")

        train_ppm = np.sqrt(train_loss) * train_ds.target_std.item()
        val_ppm = np.sqrt(val_loss) * train_ds.target_std.item()

        entry = {
            "epoch": epoch + 1,
            "train_rmse_ppm": round(train_ppm, 4),
            "val_rmse_ppm": round(val_ppm, 4),
            "lr": scheduler.get_last_lr()[0],
            "best": is_best,
        }
        if cfg.model.use_correction and hasattr(model, "correction_scale"):
            entry["correction_scale"] = round(model.correction_scale.item(), 6)

        with open(log_path, "a") as f:
            f.write(json.dumps(entry) + "\n")

        if (epoch + 1) % 10 == 0 or epoch == 0 or is_best:
            corr_str = ""
            if cfg.model.use_correction and hasattr(model, "correction_scale"):
                corr_str = f"  cs={model.correction_scale.item():.4f}"
            marker = " *" if is_best else ""
            print(f"  epoch {epoch+1:4d}  train={train_ppm:.4f}  "
                  f"val={val_ppm:.4f}  lr={scheduler.get_last_lr()[0]:.1e}"
                  f"{corr_str}{marker}")

    # ── Final evaluation ─────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print(f"Best val at epoch {best_epoch}")
    model.load_state_dict(torch.load(this_run / "best_model.pt",
                                     weights_only=True))
    model.eval()
    summary = {"finished": datetime.now().isoformat(), "best_epoch": best_epoch}

    for split_name, ds in [("Train", train_ds), ("Val", val_ds)]:
        with torch.no_grad():
            pred = model(ds.scalars, ds.kernels) * ds.target_std
        tgt = ds.targets * ds.target_std
        r2 = compute_r2(pred, tgt, ds.ring_dist, cfg)

        print(f"\n{split_name}:")
        print(f"  R² (mean): {r2['r2_mean']:.4f}  RMSE: {r2['rmse_per_atom_ppm']:.4f} ppm")
        for m in [-2, -1, 0, 1, 2]:
            print(f"    m={m:+d}: {r2[f'r2_m{m:+d}']:.3f}", end="")
        print()

        pp = compute_per_protein_r2(pred, tgt, ds.protein_boundaries,
                                    ds.protein_ids)
        vals = [v for v in pp.values() if not np.isnan(v)]
        if vals:
            print(f"  Per-protein: median={np.median(vals):.3f}  "
                  f"min={min(vals):.3f}  max={max(vals):.3f}")

        lk = split_name.lower()
        for key, val in r2.items():
            summary[f"{lk}_{key}"] = round(val, 4) if isinstance(val, float) else val
        summary[f"{lk}_per_protein_r2"] = {k: round(v, 4) for k, v in pp.items()}

    if cfg.model.use_correction:
        print(f"\nCorrection scale: {model.correction_scale.item():.4f}")

    with open(this_run / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nRun logged to {this_run}/")


def main():
    parser = argparse.ArgumentParser(description="Mutation-set calibration training")
    parser.add_argument("--config", type=str, default="calibration.toml")
    parser.add_argument("--run-name", type=str,
                        default=datetime.now().strftime("run_%Y%m%d_%H%M%S"))
    parser.add_argument("--notes", type=str, default="")
    # CLI overrides for quick experiments
    parser.add_argument("--epochs", type=int, default=None)
    parser.add_argument("--lr", type=float, default=None)
    parser.add_argument("--max-proteins", type=int, default=None)
    parser.add_argument("--correction", action="store_true", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)

    # Apply CLI overrides by rebuilding the config (frozen dataclass)
    if args.epochs or args.lr or args.max_proteins is not None or args.correction is not None:
        from dataclasses import replace
        if args.epochs:
            cfg = replace(cfg, training=replace(cfg.training, epochs=args.epochs))
        if args.lr:
            cfg = replace(cfg, training=replace(cfg.training, lr=args.lr))
        if args.max_proteins is not None:
            cfg = replace(cfg, data=replace(cfg.data, max_proteins=args.max_proteins))
        if args.correction is not None:
            cfg = replace(cfg, model=replace(cfg.model, use_correction=args.correction))

    # Ensure SDK is importable
    sys.path.insert(0, str(cfg.paths.sdk))

    train(cfg, args.run_name, args.notes)


if __name__ == "__main__":
    main()
