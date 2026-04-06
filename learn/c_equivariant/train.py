#!/usr/bin/env python3
"""
Training loop for Level C equivariant model.

Loads extracted features, builds datasets, trains, and evaluates.
Designed to run on whatever proteins are available (works with 10 or 700).

Usage:
    python learn/c_equivariant/train.py                  # mixing only
    python learn/c_equivariant/train.py --correction      # full model
    python learn/c_equivariant/train.py --epochs 100      # more epochs
"""

import argparse
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).parent.parent))
from load import load_conformation, load_pair, list_proteins, T2, RING_TYPES
from delta import match_atoms
from c_equivariant.model import make_model


class T2Dataset(Dataset):
    """WT-ALA delta T2 dataset from extracted features.

    The mutation delta isolates the ring effect: everything that isn't
    ring-related cancels between WT and ALA. The target is delta DFT T2,
    the features are delta kernel T2s.

    Normalization: ALL scaling is global (across the full dataset),
    never per-protein. Train set stats can be passed to val set via
    norm_stats to ensure consistent scaling.
    """

    def __init__(self, protein_ids, norm_stats=None):
        raw_scalars = []
        raw_kernels = []
        raw_mc = []
        raw_coulomb = []
        raw_bs = []
        raw_mopac_coulomb = []
        raw_targets = []

        skipped = 0
        for pid in protein_ids:
            wt, ala = load_pair(pid)
            if wt is None or ala is None:
                skipped += 1
                continue
            if wt.orca is None or ala.orca is None:
                skipped += 1
                continue

            # Match heavy atoms by position
            wt_idx, ala_idx, _ = match_atoms(
                wt.pos, wt.element, ala.pos, ala.element
            )
            if len(wt_idx) < 10:
                skipped += 1
                continue

            M = len(wt_idx)  # matched atoms

            # Delta DFT T2 = WT - ALA (effect of the ring being present)
            dft_t2 = T2(wt.orca[wt_idx] - ala.orca[ala_idx])  # (M, 5)

            # Scalar features from WT (ring proximity, element — these describe
            # the atom's environment WITH the ring present)
            # 13 features: 3 element one-hots (C/N/O — no H in heavy-atom matching),
            # 4 ring proximity counts, 6 per-protein-normalized kernel magnitudes
            s = np.zeros((M, 13))
            s[:, 0] = (wt.element[wt_idx] == 6).astype(float)
            s[:, 1] = (wt.element[wt_idx] == 7).astype(float)
            s[:, 2] = (wt.element[wt_idx] == 8).astype(float)
            s[:, 3:7] = wt.ring_counts[wt_idx] / 10.0
            # Magnitude scalars: per-protein normalized (matching clean backup approach)
            def _pnorm(v):
                """Per-protein normalization: divide by own std."""
                st = v.std()
                return v / (st + 1e-8)
            s[:, 7] = _pnorm(np.sqrt(np.sum(T2(wt.mc[wt_idx] - ala.mc[ala_idx])**2, axis=1)))
            s[:, 8] = _pnorm(np.sqrt(np.sum(T2(wt.coulomb[wt_idx] - ala.coulomb[ala_idx])**2, axis=1)))
            s[:, 9] = _pnorm(np.sqrt(np.sum(T2(wt.bs[wt_idx] - ala.bs[ala_idx])**2, axis=1)))
            s[:, 10] = _pnorm(np.sqrt(np.sum(T2(wt.hbond[wt_idx] - ala.hbond[ala_idx])**2, axis=1)))
            s[:, 11] = _pnorm(np.sqrt(np.sum(T2(wt.mopac_coulomb[wt_idx] - ala.mopac_coulomb[ala_idx])**2, axis=1)))
            s[:, 12] = _pnorm(np.sqrt(np.sum(T2(wt.mopac_mc[wt_idx] - ala.mopac_mc[ala_idx])**2, axis=1)))

            # Delta kernel T2s: 46 kernels × 5 components
            kernels = np.zeros((M, 46, 5))
            # BS per-type delta
            delta_bs_type = wt.bs_type_T2[wt_idx] - ala.bs_type_T2[ala_idx]
            for t in range(8):
                kernels[:, t, :] = delta_bs_type[:, t*5:(t+1)*5]
            # HM per-type delta
            delta_hm_type = wt.hm_type_T2[wt_idx] - ala.hm_type_T2[ala_idx]
            for t in range(8):
                kernels[:, 8+t, :] = delta_hm_type[:, t*5:(t+1)*5]
            # PQ per-type delta
            delta_pq_type = wt.pq_type_T2[wt_idx] - ala.pq_type_T2[ala_idx]
            for t in range(8):
                kernels[:, 16+t, :] = delta_pq_type[:, t*5:(t+1)*5]
            # Disp per-type delta
            delta_disp_type = wt.disp_type_T2[wt_idx] - ala.disp_type_T2[ala_idx]
            for t in range(8):
                kernels[:, 24+t, :] = delta_disp_type[:, t*5:(t+1)*5]
            # McConnell categories delta
            delta_mc_cat = wt.mc_category_T2[wt_idx] - ala.mc_category_T2[ala_idx]
            for c in range(5):
                kernels[:, 32+c, :] = delta_mc_cat[:, c*5:(c+1)*5]
            # Coulomb total delta
            kernels[:, 37, :] = T2(wt.coulomb[wt_idx] - ala.coulomb[ala_idx])
            # RingSusc total delta
            kernels[:, 38, :] = T2(wt.ringchi[wt_idx] - ala.ringchi[ala_idx])
            # HBond total delta
            kernels[:, 39, :] = T2(wt.hbond[wt_idx] - ala.hbond[ala_idx])
            # MopacCoulomb total delta
            kernels[:, 40, :] = T2(wt.mopac_coulomb[wt_idx] - ala.mopac_coulomb[ala_idx])
            # MopacMC categories delta
            delta_mopac_mc_cat = wt.mopac_mc_category_T2[wt_idx] - ala.mopac_mc_category_T2[ala_idx]
            for c in range(5):
                kernels[:, 41+c, :] = delta_mopac_mc_cat[:, c*5:(c+1)*5]

            # Per-protein kernel normalization (matching clean backup approach)
            # Each protein's kernels normalized by own std before stacking.
            # MLP sees relative kernel patterns within a protein, not absolute magnitudes.
            for k in range(46):
                kstd = kernels[:, k, :].std()
                if kstd > 1e-10:
                    kernels[:, k, :] /= kstd

            raw_scalars.append(s)
            raw_kernels.append(kernels)
            raw_mc.append(T2(wt.mc[wt_idx] - ala.mc[ala_idx]))
            raw_coulomb.append(T2(wt.coulomb[wt_idx] - ala.coulomb[ala_idx]))
            raw_bs.append(T2(wt.bs[wt_idx] - ala.bs[ala_idx]))
            raw_mopac_coulomb.append(T2(wt.mopac_coulomb[wt_idx] - ala.mopac_coulomb[ala_idx]))
            raw_targets.append(dft_t2)

        if skipped:
            print(f"  (skipped {skipped} proteins: missing data or <10 matched atoms)")

        self.scalars = torch.tensor(np.vstack(raw_scalars), dtype=torch.float32)
        self.kernels = torch.tensor(np.vstack(raw_kernels), dtype=torch.float32)
        self.mc_t2 = torch.tensor(np.vstack(raw_mc), dtype=torch.float32)
        self.coulomb_t2 = torch.tensor(np.vstack(raw_coulomb), dtype=torch.float32)
        self.bs_t2 = torch.tensor(np.vstack(raw_bs), dtype=torch.float32)
        self.mopac_coulomb_t2 = torch.tensor(np.vstack(raw_mopac_coulomb), dtype=torch.float32)
        self.targets = torch.tensor(np.vstack(raw_targets), dtype=torch.float32)

        # Target scaling (for training stability)
        self.target_std = self.targets.std()
        self.targets_scaled = self.targets / self.target_std

        # Scale individual T2 inputs for correction head
        self.mc_std = self.mc_t2.std() + 1e-8
        self.mc_t2 = self.mc_t2 / self.mc_std
        self.coul_std = self.coulomb_t2.std() + 1e-8
        self.coulomb_t2 = self.coulomb_t2 / self.coul_std
        self.bs_std = self.bs_t2.std() + 1e-8
        self.bs_t2 = self.bs_t2 / self.bs_std
        self.mopac_coul_std = self.mopac_coulomb_t2.std() + 1e-8
        self.mopac_coulomb_t2 = self.mopac_coulomb_t2 / self.mopac_coul_std

        # Kernels: already per-protein normalized before stacking.
        # Scalars: element one-hots are 0/1 (fine as categorical),
        # ring counts /10, magnitudes per-protein normalized.

    def norm_stats(self):
        """Return normalization stats (for logging, not for val set)."""
        return {
            'target_std': self.target_std.item(),
        }

    def to(self, device):
        """Move entire dataset to device (GPU). Returns self."""
        self.scalars = self.scalars.to(device)
        self.kernels = self.kernels.to(device)
        self.mc_t2 = self.mc_t2.to(device)
        self.coulomb_t2 = self.coulomb_t2.to(device)
        self.bs_t2 = self.bs_t2.to(device)
        self.mopac_coulomb_t2 = self.mopac_coulomb_t2.to(device)
        self.targets = self.targets.to(device)
        self.targets_scaled = self.targets_scaled.to(device)
        self.target_std = self.target_std.to(device)
        return self

    def __len__(self):
        return len(self.scalars)

    def __getitem__(self, idx):
        return (self.scalars[idx], self.kernels[idx],
                self.mc_t2[idx], self.coulomb_t2[idx], self.bs_t2[idx],
                self.mopac_coulomb_t2[idx],
                self.targets_scaled[idx])


def build_datasets():
    """Build train/val datasets from .npy files. Each split normalizes independently."""
    proteins = list_proteins()
    np.random.seed(42)
    np.random.shuffle(proteins)
    n_val = max(2, len(proteins) // 5)
    train_ds = T2Dataset(proteins[n_val:])
    val_ds = T2Dataset(proteins[:n_val])
    return train_ds, val_ds


def compute_r2(model, dataset, target_std, device):
    """Compute R² on a dataset. Returns (r2, rmse_ppm)."""
    loader = DataLoader(dataset, batch_size=4096)
    all_pred = []
    all_tgt = []
    with torch.no_grad():
        for batch in loader:
            s, k, mc, co, bs, mopac_co, tgt = batch
            all_pred.append(model(s, k, mc, co, bs, mopac_co))
            all_tgt.append(tgt)
    pred = torch.cat(all_pred) * target_std
    tgt = torch.cat(all_tgt) * target_std
    ss_res = torch.sum((tgt - pred) ** 2).item()
    ss_tot = torch.sum((tgt - tgt.mean()) ** 2).item()
    r2 = 1.0 - ss_res / ss_tot
    rmse = np.sqrt(ss_res / pred.numel())  # per scalar element, not per atom
    return r2, rmse


def train(args):
    import json as json_mod
    from datetime import datetime

    # Run output directory
    run_dir = Path(__file__).parent.parent / "runs"
    run_dir.mkdir(exist_ok=True)
    run_name = args.run_name or datetime.now().strftime("run_%Y%m%d_%H%M%S")
    this_run = run_dir / run_name
    this_run.mkdir(exist_ok=True)

    train_ds, val_ds = build_datasets()

    n_proteins = len(list_proteins())
    header = {
        "run_name": run_name,
        "started": datetime.now().isoformat(),
        "n_proteins": n_proteins,
        "n_train_atoms": len(train_ds),
        "n_val_atoms": len(val_ds),
        "target_std_ppm": train_ds.target_std.item(),
        "n_kernels": 46,
        "n_scalar_features": 13,
        "epochs": args.epochs,
        "warmup_epochs": args.warmup_epochs,
        "batch_size": args.batch_size,
        "lr": args.lr,
        "use_correction": args.correction,
        "notes": args.notes or "",
    }

    print(f"Run: {run_name}")
    print(f"Train atoms: {len(train_ds)}, Val atoms: {len(val_ds)}")
    print(f"Target std: {train_ds.target_std:.4f} ppm")
    print(f"Proteins: {n_proteins}")

    # Report kernel activity (per-protein normalized, so check across stacked data)
    active = sum(1 for k in range(46) if train_ds.kernels[:, k, :].std() > 1e-10)
    print(f"Active kernels: {active}/46")
    header["active_kernels"] = active

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    train_ds.to(device)
    val_ds.to(device)
    print(f"Datasets on {device}\n")

    train_loader = DataLoader(train_ds, batch_size=args.batch_size,
                              shuffle=True, drop_last=True)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size)

    model = make_model(
        n_scalar_features=13, n_kernels=46,
        use_correction=args.correction
    ).to(device)

    print(f"Model: {'mixing + correction' if args.correction else 'mixing only'}")
    print(f"Parameters: {model.parameter_count()}")
    print(f"Device: {device}\n")
    header["n_parameters"] = model.parameter_count()

    # Write header
    with open(this_run / "header.json", "w") as f:
        json_mod.dump(header, f, indent=2)

    # Epoch log — one JSON line per epoch
    log_path = this_run / "epochs.jsonl"

    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    best_val_loss = float("inf")
    best_epoch = 0
    for epoch in range(args.epochs):
        # Train
        model.train()
        train_loss = 0.0
        n_train = 0
        for batch in train_loader:
            s, k, mc, co, bs, mopac_co, tgt = batch
            pred = model(s, k, mc, co, bs, mopac_co)
            loss = nn.functional.mse_loss(pred, tgt)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * len(s)
            n_train += len(s)
        train_loss /= n_train

        # Validate
        model.eval()
        val_loss = 0.0
        n_val_atoms = 0
        with torch.no_grad():
            for batch in val_loader:
                s, k, mc, co, bs, mopac_co, tgt = batch
                pred = model(s, k, mc, co, bs, mopac_co)
                loss = nn.functional.mse_loss(pred, tgt)
                val_loss += loss.item() * len(s)
                n_val_atoms += len(s)
        val_loss /= n_val_atoms

        scheduler.step()

        is_best = val_loss < best_val_loss
        if is_best:
            best_val_loss = val_loss
            best_epoch = epoch + 1
            torch.save(model.state_dict(), this_run / "best_model.pt")

        train_ppm = np.sqrt(train_loss) * train_ds.target_std.item()
        val_ppm = np.sqrt(val_loss) * train_ds.target_std.item()

        # Log every epoch
        entry = {
            "epoch": epoch + 1,
            "train_rmse_ppm": round(train_ppm, 4),
            "val_rmse_ppm": round(val_ppm, 4),
            "lr": scheduler.get_last_lr()[0],
            "best": is_best,
        }
        with open(log_path, "a") as f:
            f.write(json_mod.dumps(entry) + "\n")

        if (epoch + 1) % 10 == 0 or epoch == 0 or is_best:
            marker = " *" if is_best else ""
            print(f"  epoch {epoch+1:4d}  train={train_ppm:.4f} ppm  "
                  f"val={val_ppm:.4f} ppm  lr={scheduler.get_last_lr()[0]:.1e}{marker}")

    # Final evaluation — load best model
    print(f"\n{'='*50}")
    print(f"Best val at epoch {best_epoch}")
    model.load_state_dict(torch.load(this_run / "best_model.pt", weights_only=True))
    model.eval()
    target_std = train_ds.target_std

    r2_train, rmse_train = compute_r2(model, train_ds, target_std, device)
    r2_val, rmse_val = compute_r2(model, val_ds, target_std, device)

    print(f"Train R²: {r2_train:.4f}  RMSE: {rmse_train:.4f} ppm")
    print(f"Val R²:   {r2_val:.4f}  RMSE: {rmse_val:.4f} ppm")
    if args.correction:
        cs = model.correction_scale.item()
        print(f"Correction scale: {cs:.4f}")

    # Write final summary
    summary = {
        "finished": datetime.now().isoformat(),
        "best_epoch": best_epoch,
        "train_r2": round(r2_train, 4),
        "val_r2": round(r2_val, 4),
        "train_rmse_ppm": round(rmse_train, 4),
        "val_rmse_ppm": round(rmse_val, 4),
    }
    with open(this_run / "summary.json", "w") as f:
        json_mod.dump(summary, f, indent=2)

    print(f"\nRun logged to {this_run}/")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--warmup-epochs", type=int, default=0)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--correction", action="store_true",
                        help="use equivariant correction head (Level C full)")
    parser.add_argument("--no-cache", action="store_true",
                        help="rebuild dataset from .npy files (ignore cache)")
    parser.add_argument("--run-name", type=str, default=None,
                        help="name for this run (default: timestamped)")
    parser.add_argument("--notes", type=str, default=None,
                        help="free-text notes saved with the run")
    args = parser.parse_args()
    train(args)


if __name__ == "__main__":
    main()
