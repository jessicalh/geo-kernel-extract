#!/usr/bin/env python3
"""Shared split, centering, and normalization helpers for e3nn T2 fitters."""

import numpy as np


TEST_FRAME_FRACTION = 0.20
FRAME_SPLIT_SEED = 0
MIN_FRAME_SPLIT_FRAMES = 5
DEFAULT_SPLIT = "blocked"
PURGE_NEIGHBOUR_FRAMES = 1


def centred_by_train_atom(x, group_atom_idx, train_mask):
    """Center each atom's rows using only that atom's training rows."""
    x = np.asarray(x, dtype=float)
    group_atom_idx = np.asarray(group_atom_idx)
    train_mask = np.asarray(train_mask, dtype=bool)
    out = np.full_like(x, np.nan, dtype=float)
    for atom in np.unique(group_atom_idx):
        m = group_atom_idx == atom
        train_atom = m & train_mask
        if train_atom.sum() == 0:
            continue
        out[m] = x[m] - x[train_atom].mean(axis=0, keepdims=True)
    return out


def normalised_by_train_rows(x, train_mask, eps=1e-8):
    """Normalize columns from training rows only; never fit scalers on held-out rows."""
    x = np.asarray(x, dtype=np.float32)
    train_mask = np.asarray(train_mask, dtype=bool)
    if x.ndim != 2:
        raise ValueError(f"expected a 2D feature matrix, got shape {x.shape}")
    if train_mask.shape[0] != x.shape[0]:
        raise ValueError(f"train mask length {train_mask.shape[0]} != rows {x.shape[0]}")

    source = x[train_mask]
    if source.shape[0] == 0:
        mean = np.zeros(x.shape[1], dtype=np.float32)
        scale = np.ones(x.shape[1], dtype=np.float32)
    else:
        mean = source.mean(axis=0)
        std = source.std(axis=0)
        scale = np.where(std > eps, std, 1.0).astype(np.float32)
    return (x - mean) / scale, mean, scale


def make_split_masks(row_frames, strategy=DEFAULT_SPLIT, seed=FRAME_SPLIT_SEED,
                     purge_frames=PURGE_NEIGHBOUR_FRAMES):
    frames = np.sort(np.unique(row_frames))
    empty = np.zeros(len(row_frames), dtype=bool)
    if len(frames) < MIN_FRAME_SPLIT_FRAMES:
        return empty.copy(), empty.copy(), {
            "split_strategy": strategy,
            "test_frames": 0,
            "purged_train_frames": 0,
            "cross_split_lag1_pairs": 0,
        }
    n_test = max(1, int(TEST_FRAME_FRACTION * len(frames)))
    rng = np.random.default_rng(seed)

    if strategy == "random":
        test_frames = set(rng.choice(frames, n_test, replace=False))
        train_frames = set(frames) - test_frames
        purged = set()
    elif strategy == "blocked":
        start = int(rng.integers(0, len(frames) - n_test + 1))
        stop = start + n_test
        test_frames = set(frames[start:stop])
        purge_lo = max(0, start - max(0, purge_frames))
        purge_hi = min(len(frames), stop + max(0, purge_frames))
        purged = set(frames[purge_lo:start]) | set(frames[stop:purge_hi])
        train_frames = set(frames) - test_frames - purged
    else:
        raise ValueError(f"unknown split strategy {strategy!r}")

    train = np.array([f in train_frames for f in row_frames], dtype=bool)
    test = np.array([f in test_frames for f in row_frames], dtype=bool)

    frame_to_split = {
        f: ("test" if f in test_frames else "train" if f in train_frames else "purged")
        for f in frames
    }
    cross = 0
    for a, b in zip(frames[:-1], frames[1:]):
        sa = frame_to_split[a]
        sb = frame_to_split[b]
        if {sa, sb} == {"train", "test"}:
            cross += 1
    return train, test, {
        "split_strategy": strategy,
        "test_frames": int(len(test_frames)),
        "purged_train_frames": int(len(purged)),
        "cross_split_lag1_pairs": int(cross),
    }
