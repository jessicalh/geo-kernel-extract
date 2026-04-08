#!/usr/bin/env python3
"""
Training loop for Level C equivariant calibration model.

Uses CalibrationDataset (typed protein loader, C++ matched_mask,
48 kernel T2s, 68 scalar features) against calibration/features/.

Usage:
    python learn/c_equivariant/train.py --run CalibrationExtractionTest
    python learn/c_equivariant/train.py --run CalibrationExtractionTest --correction
    python learn/c_equivariant/train.py --run CalibrationExtractionTest --distance-weighted
    python learn/c_equivariant/train.py --run CalibrationExtractionTest --epochs 500
"""

import argparse
import json as json_mod
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

sys.path.insert(0, str(Path(__file__).parent.parent))
from c_equivariant.dataset import (
    CalibrationDataset, N_KERNELS, N_SCALAR_FEATURES,
    compute_r2, compute_naive_baselines, compute_per_protein_r2,
)
from c_equivariant.model import make_model
def list_proteins(run_dir: Path) -> list[str]:
    """Protein IDs with complete extracted features in a run directory."""
    return sorted(d.name for d in run_dir.iterdir()
                  if d.is_dir() and (d / "pos.npy").exists()
                  and (d / "ring_contributions.npy").exists())

REPO = Path(__file__).resolve().parent.parent.parent
DEFAULT_FEATURES = REPO / "calibration" / "features"

# Distance weighting: exp(-d/tau) focuses loss on atoms near removed rings
DIST_WEIGHT_TAU = 8.0  # Angstroms


def build_datasets(run_name: str, features_base: Path = DEFAULT_FEATURES):
    """Build train/val datasets from a named extraction run.

    Val set uses train set's normalization stats so they share the same scale.
    """
    features_dir = features_base / run_name
    if not features_dir.is_dir():
        print(f"ERROR: {features_dir} not found. Run extract.py first.",
              file=sys.stderr)
        sys.exit(1)

    proteins = list_proteins(features_dir)
    np.random.seed(42)
    np.random.shuffle(proteins)
    n_val = max(2, len(proteins) // 5)
    train_ds = CalibrationDataset(proteins[n_val:], features_dir)
    val_ds = CalibrationDataset(proteins[:n_val], features_dir,
                                norm_stats=train_ds.norm_stats)
    return train_ds, val_ds, len(proteins)


def train(args):
    run_dir = Path(__file__).parent.parent / "runs"
    run_dir.mkdir(exist_ok=True)
    run_name = args.run_name or datetime.now().strftime("run_%Y%m%d_%H%M%S")
    this_run = run_dir / run_name
    this_run.mkdir(exist_ok=True)

    train_ds, val_ds, n_proteins = build_datasets(args.extraction_run)

    if args.shuffle_targets:
        print("*** SANITY CHECK: targets randomly permuted ***")
        perm = torch.randperm(len(train_ds))
        train_ds.targets = train_ds.targets[perm]
        train_ds.ring_dist = train_ds.ring_dist[perm]
        perm_v = torch.randperm(len(val_ds))
        val_ds.targets = val_ds.targets[perm_v]
        val_ds.ring_dist = val_ds.ring_dist[perm_v]

    header = {
        "run_name": run_name,
        "extraction_run": args.extraction_run,
        "started": datetime.now().isoformat(),
        "n_proteins": n_proteins,
        "n_train_atoms": len(train_ds),
        "n_val_atoms": len(val_ds),
        "target_std_ppm": train_ds.target_std.item(),
        "n_kernels": N_KERNELS,
        "n_scalar_features": train_ds.scalars.shape[1],
        "epochs": args.epochs,
        "batch_size": args.batch_size,
        "lr": args.lr,
        "use_correction": args.correction,
        "distance_weighted": args.distance_weighted,
        "shuffle_targets": args.shuffle_targets,
        "notes": args.notes or "",
    }

    print(f"Run: {run_name}")
    print(f"Train atoms: {len(train_ds)}, Val atoms: {len(val_ds)}")
    print(f"Target std: {train_ds.target_std:.4f} ppm")
    print(f"Proteins: {n_proteins}")

    active = sum(1 for k in range(N_KERNELS)
                 if train_ds.kernels[:, k, :].std() > 1e-10)
    print(f"Active kernels: {active}/{N_KERNELS}")
    header["active_kernels"] = active

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    train_ds.to(device)
    val_ds.to(device)
    print(f"Device: {device}")

    # ── Naive physics baselines (before any learning) ────────────
    print(f"\n{'─' * 60}")
    print("Physics baselines (no learning):")
    for name, ds in [("train", train_ds), ("val", val_ds)]:
        bl = compute_naive_baselines(ds)
        print(f"  {name}:")
        print(f"    BS only:      R²={bl['r2_bs_only']:.4f}")
        print(f"    Ring sum:     R²={bl['r2_ring_sum']:.4f}")
        print(f"    All sum:      R²={bl['r2_all_sum']:.4f}")
        if "ring_sum_0-4A" in bl:
            print(f"    Ring sum by dist: ", end="")
            for label in ["0-4A", "4-8A", "8-12A", "12+A"]:
                v = bl.get(f"ring_sum_{label}", float("nan"))
                print(f"{label}:{v:.3f}  ", end="")
            print()
        header[f"baseline_{name}"] = {k: round(v, 4) if isinstance(v, float) else v
                                      for k, v in bl.items()}
    print(f"{'─' * 60}\n")

    # ── Distance weights for loss ────────────────────────────────
    dist_weights_train = None
    dist_weights_val = None
    if args.distance_weighted:
        print(f"Distance-weighted loss: tau={DIST_WEIGHT_TAU} A")
        dist_weights_train = torch.exp(-train_ds.ring_dist / DIST_WEIGHT_TAU)
        dist_weights_val = torch.exp(-val_ds.ring_dist / DIST_WEIGHT_TAU)
        # Normalize so mean weight = 1 (doesn't change effective LR)
        dist_weights_train /= dist_weights_train.mean()
        dist_weights_val /= dist_weights_val.mean()

    model = make_model(
        n_scalar_features=N_SCALAR_FEATURES,
        n_kernels=N_KERNELS,
        use_correction=args.correction,
    )
    # Set per-kernel gate thresholds from training data
    model.mixing.set_gate_thresholds(
        torch.tensor(train_ds.norm_stats.gate_threshold, dtype=torch.float32))
    model = model.to(device)

    print(f"Model: {'mixing + correction' if args.correction else 'mixing only'}")
    print(f"Parameters: {model.parameter_count()}\n")
    header["n_parameters"] = model.parameter_count()

    with open(this_run / "header.json", "w") as f:
        json_mod.dump(header, f, indent=2)

    log_path = this_run / "epochs.jsonl"
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=args.epochs)

    best_val_loss = float("inf")
    best_epoch = 0

    # All data already on GPU — skip DataLoader, shuffle indices on device
    N_train = len(train_ds)
    N_val = len(val_ds)
    bs = args.batch_size

    for epoch in range(args.epochs):
        model.train()
        perm = torch.randperm(N_train, device=device)
        train_loss = 0.0
        for i in range(0, N_train, bs):
            idx = perm[i:i+bs]
            s = train_ds.scalars[idx]
            k = train_ds.kernels[idx]
            tgt = train_ds.targets[idx]
            pred = model(s, k)
            if dist_weights_train is not None:
                w = dist_weights_train[idx]
                per_atom = ((pred - tgt) ** 2).mean(dim=1)
                loss = (per_atom * w).mean()
            else:
                loss = nn.functional.mse_loss(pred, tgt)
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
        if args.correction and hasattr(model, "correction_scale"):
            entry["correction_scale"] = round(model.correction_scale.item(), 6)
        with open(log_path, "a") as f:
            f.write(json_mod.dumps(entry) + "\n")

        if (epoch + 1) % 10 == 0 or epoch == 0 or is_best:
            marker = " *" if is_best else ""
            corr_str = ""
            if args.correction and hasattr(model, "correction_scale"):
                corr_str = f"  cs={model.correction_scale.item():.4f}"
            print(f"  epoch {epoch+1:4d}  train={train_ppm:.4f} ppm  "
                  f"val={val_ppm:.4f} ppm  lr={scheduler.get_last_lr()[0]:.1e}"
                  f"{corr_str}{marker}")

    # ── Final evaluation ─────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print(f"Best val at epoch {best_epoch}")
    model.load_state_dict(torch.load(this_run / "best_model.pt",
                                     weights_only=True))
    model.eval()

    summary = {"finished": datetime.now().isoformat(), "best_epoch": best_epoch}

    for split_name, ds in [("Train", train_ds), ("Val", val_ds)]:
        all_pred, all_tgt = [], []
        with torch.no_grad():
            for s, k, tgt in DataLoader(ds, batch_size=4096):
                all_pred.append(model(s, k))
                all_tgt.append(tgt)
        pred = torch.cat(all_pred) * ds.target_std
        tgt = torch.cat(all_tgt) * ds.target_std
        r2 = compute_r2(pred, tgt, ring_dist=ds.ring_dist)

        print(f"\n{split_name}:")
        print(f"  R² (mean of components): {r2['r2_mean']:.4f}")
        print(f"  R² (flat):               {r2['r2_flat']:.4f}")
        print(f"  RMSE per atom:           {r2['rmse_per_atom_ppm']:.4f} ppm")
        print(f"  Per component:  ", end="")
        for m in [-2, -1, 0, 1, 2]:
            print(f"m={m:+d}:{r2[f'r2_m{m:+d}']:.3f}  ", end="")
        print()
        print(f"  By distance:    ", end="")
        for label in ["0-4A", "4-8A", "8-12A", "12+A"]:
            n = r2.get(f"n_{label}", 0)
            v = r2.get(f"r2_{label}", float("nan"))
            print(f"{label}:{v:.3f}({n})  ", end="")
        print()

        # Per-protein R²
        pp_r2 = compute_per_protein_r2(pred, tgt, ds.protein_boundaries,
                                       ds.protein_ids)
        vals = [v for v in pp_r2.values() if not np.isnan(v)]
        if vals:
            print(f"  Per-protein R²: median={np.median(vals):.3f}  "
                  f"min={min(vals):.3f}  max={max(vals):.3f}")
            worst = sorted(pp_r2.items(), key=lambda x: x[1])[:3]
            print(f"  Worst proteins: {', '.join(f'{p}({v:.3f})' for p, v in worst)}")

        lk = split_name.lower()
        for key, val in r2.items():
            summary[f"{lk}_{key}"] = round(val, 4) if isinstance(val, float) else val
        summary[f"{lk}_per_protein_r2"] = {k: round(v, 4) for k, v in pp_r2.items()}

    if args.correction:
        print(f"\nCorrection scale: {model.correction_scale.item():.4f}")

    # ── Save kernel diagnostics for inspection ─────────────────────
    with torch.no_grad():
        sample_s = train_ds.scalars[:1000]
        sample_k = train_ds.kernels[:1000]
        # MLP relative weights (before gating)
        sample_relative = model.mixing.mlp(sample_s)  # (1000, N_KERNELS)
        summary["kernel_weight_mean"] = sample_relative.mean(dim=0).cpu().tolist()
        summary["kernel_weight_std"] = sample_relative.std(dim=0).cpu().tolist()
        # Kernel self-reported magnitudes (the gate signal)
        sample_magnitude = sample_k.norm(dim=-1)  # (1000, N_KERNELS)
        summary["kernel_magnitude_mean"] = sample_magnitude.mean(dim=0).cpu().tolist()
        summary["kernel_magnitude_std"] = sample_magnitude.std(dim=0).cpu().tolist()
        # Gated weights (what actually drives prediction)
        sample_gated = sample_relative * sample_magnitude
        summary["kernel_gated_mean"] = sample_gated.mean(dim=0).cpu().tolist()
        summary["kernel_gated_std"] = sample_gated.std(dim=0).cpu().tolist()

    with open(this_run / "summary.json", "w") as f:
        json_mod.dump(summary, f, indent=2)
    print(f"\nRun logged to {this_run}/")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--extraction-run", type=str,
                        default="CalibrationExtractionTest",
                        help="extraction run name under calibration/features/")
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--correction", action="store_true",
                        help="use equivariant correction head")
    parser.add_argument("--distance-weighted", action="store_true",
                        help="DEPRECATED: weight loss by exp(-d/8A). "
                             "Superseded by kernel self-gating in the mixing head.")
    parser.add_argument("--run-name", type=str, default=None,
                        help="name for this run (default: timestamped)")
    parser.add_argument("--notes", type=str, default=None)
    parser.add_argument("--shuffle-targets", action="store_true",
                        help="permute targets as sanity check (R² should go negative)")
    args = parser.parse_args()
    train(args)


if __name__ == "__main__":
    main()
