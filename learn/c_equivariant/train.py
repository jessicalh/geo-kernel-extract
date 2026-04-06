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
from load import load_conformation, list_proteins, T2, RING_TYPES
from c_equivariant.model import make_model


class T2Dataset(Dataset):
    """Per-atom T2 dataset from extracted features.

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
        raw_targets = []

        for pid in protein_ids:
            wt = load_conformation(pid, "wt")
            if wt is None or wt.orca is None:
                continue

            N = wt.n_atoms
            dft_t2 = T2(wt.orca)  # (N, 5)

            # Scalar features (L=0) — raw values, no per-protein normalization
            s = np.zeros((N, 12))
            s[:, 0] = (wt.element == 1).astype(float)
            s[:, 1] = (wt.element == 6).astype(float)
            s[:, 2] = (wt.element == 7).astype(float)
            s[:, 3] = (wt.element == 8).astype(float)
            s[:, 4:8] = wt.ring_counts / 10.0
            # Raw magnitudes — normalized globally after stacking
            s[:, 8] = np.sqrt(np.sum(T2(wt.mc)**2, axis=1))
            s[:, 9] = np.sqrt(np.sum(T2(wt.coulomb)**2, axis=1))
            s[:, 10] = np.sqrt(np.sum(T2(wt.bs)**2, axis=1))
            s[:, 11] = np.sqrt(np.sum(T2(wt.hbond)**2, axis=1))

            # Kernel T2s: 40 kernels × 5 components — RAW, no per-protein scaling
            kernels = np.zeros((N, 40, 5))
            for t in range(8):
                kernels[:, t, :] = wt.bs_type_T2[:, t*5:(t+1)*5]
            for t in range(8):
                kernels[:, 8+t, :] = wt.hm_type_T2[:, t*5:(t+1)*5]
            for t in range(8):
                kernels[:, 16+t, :] = wt.pq_type_T2[:, t*5:(t+1)*5]
            for t in range(8):
                kernels[:, 24+t, :] = wt.disp_type_T2[:, t*5:(t+1)*5]
            for c in range(5):
                kernels[:, 32+c, :] = wt.mc_category_T2[:, c*5:(c+1)*5]
            kernels[:, 37, :] = T2(wt.coulomb)
            kernels[:, 38, :] = T2(wt.ringchi)
            kernels[:, 39, :] = T2(wt.hbond)

            raw_scalars.append(s)
            raw_kernels.append(kernels)
            raw_mc.append(T2(wt.mc))
            raw_coulomb.append(T2(wt.coulomb))
            raw_bs.append(T2(wt.bs))
            raw_targets.append(dft_t2)

        self.scalars = torch.tensor(np.vstack(raw_scalars), dtype=torch.float32)
        self.kernels = torch.tensor(np.vstack(raw_kernels), dtype=torch.float32)
        self.mc_t2 = torch.tensor(np.vstack(raw_mc), dtype=torch.float32)
        self.coulomb_t2 = torch.tensor(np.vstack(raw_coulomb), dtype=torch.float32)
        self.bs_t2 = torch.tensor(np.vstack(raw_bs), dtype=torch.float32)
        self.targets = torch.tensor(np.vstack(raw_targets), dtype=torch.float32)

        # Compute or reuse normalization stats
        if norm_stats is None:
            # GLOBAL normalization across entire dataset
            self.kernel_stds = torch.zeros(40)
            for k in range(40):
                self.kernel_stds[k] = self.kernels[:, k, :].std()

            self.scalar_stds = torch.zeros(12)
            for i in range(12):
                self.scalar_stds[i] = self.scalars[:, i].std()

            self.target_std = self.targets.std()
            self.mc_std = self.mc_t2.std() + 1e-8
            self.coul_std = self.coulomb_t2.std() + 1e-8
            self.bs_std = self.bs_t2.std() + 1e-8
        else:
            # Reuse train set stats (for val set)
            self.kernel_stds = norm_stats['kernel_stds']
            self.scalar_stds = norm_stats['scalar_stds']
            self.target_std = norm_stats['target_std']
            self.mc_std = norm_stats['mc_std']
            self.coul_std = norm_stats['coul_std']
            self.bs_std = norm_stats['bs_std']

        # Apply global normalization
        for k in range(40):
            if self.kernel_stds[k] > 1e-10:
                self.kernels[:, k, :] /= self.kernel_stds[k]

        for i in [8, 9, 10, 11]:  # magnitude scalars only
            if self.scalar_stds[i] > 1e-10:
                self.scalars[:, i] /= self.scalar_stds[i]

        self.targets_scaled = self.targets / self.target_std
        self.mc_t2 = self.mc_t2 / self.mc_std
        self.coulomb_t2 = self.coulomb_t2 / self.coul_std
        self.bs_t2 = self.bs_t2 / self.bs_std

    def norm_stats(self):
        """Return normalization stats for passing to val set."""
        return {
            'kernel_stds': self.kernel_stds,
            'scalar_stds': self.scalar_stds,
            'target_std': self.target_std,
            'mc_std': self.mc_std,
            'coul_std': self.coul_std,
            'bs_std': self.bs_std,
        }

    def to(self, device):
        """Move entire dataset to device (GPU). Returns self."""
        self.scalars = self.scalars.to(device)
        self.kernels = self.kernels.to(device)
        self.mc_t2 = self.mc_t2.to(device)
        self.coulomb_t2 = self.coulomb_t2.to(device)
        self.bs_t2 = self.bs_t2.to(device)
        self.targets = self.targets.to(device)
        self.targets_scaled = self.targets_scaled.to(device)
        self.target_std = self.target_std.to(device)
        return self

    def __len__(self):
        return len(self.scalars)

    def __getitem__(self, idx):
        return (self.scalars[idx], self.kernels[idx],
                self.mc_t2[idx], self.coulomb_t2[idx], self.bs_t2[idx],
                self.targets_scaled[idx])


def load_or_build_datasets(cache_path=None):
    """Load cached dataset or build from .npy files."""
    if cache_path and Path(cache_path).exists():
        t0 = time.time()
        c = torch.load(cache_path, weights_only=False)
        train_ds = T2Dataset.__new__(T2Dataset)
        train_ds.scalars = c['train_scalars']
        train_ds.kernels = c['train_kernels']
        train_ds.mc_t2 = c['train_mc']
        train_ds.coulomb_t2 = c['train_coulomb']
        train_ds.bs_t2 = c['train_bs']
        train_ds.targets_scaled = c['train_targets_scaled']
        train_ds.targets = c['train_targets']
        train_ds.target_std = c['train_target_std']
        stats = c['norm_stats']
        train_ds.kernel_stds = stats['kernel_stds']
        train_ds.scalar_stds = stats['scalar_stds']
        train_ds.mc_std = stats['mc_std']
        train_ds.coul_std = stats['coul_std']
        train_ds.bs_std = stats['bs_std']

        val_ds = T2Dataset.__new__(T2Dataset)
        val_ds.scalars = c['val_scalars']
        val_ds.kernels = c['val_kernels']
        val_ds.mc_t2 = c['val_mc']
        val_ds.coulomb_t2 = c['val_coulomb']
        val_ds.bs_t2 = c['val_bs']
        val_ds.targets_scaled = c['val_targets_scaled']
        val_ds.targets = c['val_targets']
        val_ds.target_std = c['train_target_std']
        print(f"Loaded cached dataset in {time.time()-t0:.1f}s")
        return train_ds, val_ds

    proteins = list_proteins()
    np.random.seed(42)
    np.random.shuffle(proteins)
    n_val = max(2, len(proteins) // 5)
    train_ds = T2Dataset(proteins[n_val:])
    val_ds = T2Dataset(proteins[:n_val], norm_stats=train_ds.norm_stats())
    return train_ds, val_ds


def train(args):
    cache_path = Path(__file__).parent.parent / "features" / "FirstExtraction" / "cached_dataset.pt"
    train_ds, val_ds = load_or_build_datasets(cache_path if not args.no_cache else None)

    print(f"Train atoms: {len(train_ds)}, Val atoms: {len(val_ds)}")
    print(f"Target std: {train_ds.target_std:.4f} ppm")

    # Report kernel scales for sanity
    active = [(k, train_ds.kernel_stds[k].item()) for k in range(40)
              if train_ds.kernel_stds[k] > 1e-10]
    print(f"Active kernels: {len(active)}/40")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    train_ds.to(device)
    val_ds.to(device)
    print(f"Datasets on {device}\n")

    train_loader = DataLoader(train_ds, batch_size=args.batch_size,
                              shuffle=True, drop_last=True)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size)

    model = make_model(
        n_scalar_features=12, n_kernels=40,
        use_correction=args.correction
    ).to(device)

    print(f"Model: {'mixing + correction' if args.correction else 'mixing only'}")
    print(f"Parameters: {model.parameter_count()}")
    print(f"Device: {device}\n")

    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(
        optimizer, T_0=50, T_mult=2  # restarts at 50, 150, 350
    )

    best_val_loss = float("inf")
    for epoch in range(args.epochs):
        # Train
        model.train()
        train_loss = 0.0
        n_train = 0
        for batch in train_loader:
            s, k, mc, co, bs, tgt = batch
            pred = model(s, k, mc, co, bs)
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
                s, k, mc, co, bs, tgt = batch
                pred = model(s, k, mc, co, bs)
                loss = nn.functional.mse_loss(pred, tgt)
                val_loss += loss.item() * len(s)
                n_val_atoms += len(s)
        val_loss /= n_val_atoms

        scheduler.step()

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            marker = " *"
        else:
            marker = ""

        if (epoch + 1) % 10 == 0 or epoch == 0:
            train_ppm = np.sqrt(train_loss) * train_ds.target_std.item()
            val_ppm = np.sqrt(val_loss) * train_ds.target_std.item()  # SAME std
            print(f"  epoch {epoch+1:4d}  train={train_ppm:.4f} ppm  "
                  f"val={val_ppm:.4f} ppm  lr={scheduler.get_last_lr()[0]:.1e}{marker}")

    # Final evaluation
    print(f"\n{'='*50}")
    model.eval()
    target_std = train_ds.target_std
    with torch.no_grad():
        all_pred = []
        all_tgt = []
        for batch in DataLoader(train_ds, batch_size=4096):
            s, k, mc, co, bs, tgt = batch
            all_pred.append(model(s, k, mc, co, bs))
            all_tgt.append(tgt)
        pred = torch.cat(all_pred) * target_std
        tgt = torch.cat(all_tgt) * target_std

        ss_res = torch.sum((tgt - pred) ** 2).item()
        ss_tot = torch.sum((tgt - tgt.mean()) ** 2).item()
        r2_train = 1.0 - ss_res / ss_tot

        all_pred = []
        all_tgt = []
        for batch in DataLoader(val_ds, batch_size=4096):
            s, k, mc, co, bs, tgt = batch
            all_pred.append(model(s, k, mc, co, bs))
            all_tgt.append(tgt)
        pred = torch.cat(all_pred) * target_std  # SAME std as train
        tgt = torch.cat(all_tgt) * target_std

        ss_res = torch.sum((tgt - pred) ** 2).item()
        ss_tot = torch.sum((tgt - tgt.mean()) ** 2).item()
        r2_val = 1.0 - ss_res / ss_tot

    print(f"Train R²: {r2_train:.4f}")
    print(f"Val R²:   {r2_val:.4f}")
    print(f"(Level A linear baseline: R² ≈ 0.175)")
    if args.correction:
        cs = model.correction_scale.item()
        print(f"Correction scale: {cs:.4f} (fraction of L=2 beyond linear mixing)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--batch-size", type=int, default=2048)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--correction", action="store_true",
                        help="use equivariant correction head (Level C full)")
    parser.add_argument("--no-cache", action="store_true",
                        help="rebuild dataset from .npy files (ignore cache)")
    args = parser.parse_args()
    train(args)


if __name__ == "__main__":
    main()
