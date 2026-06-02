#!/usr/bin/env python3
"""Layer-3 joint ensemble fit on the corrected rediscover substrate.

Reads only emitted CSV/NPY payloads. No H5, no ORCA, no C++ producer calls, no
kernel reconstruction, and no tensor projection recompute. Local T2 sidecars and
emitted McConnell source tensors are mapped once through the frozen
change_of_basis.get_C() matrix.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import change_of_basis as cob


FV_HN = {1, 2}
FV_N = {4, 5}
FV_CA = {6}
FV_C = {7}
FV_O = {8}
FV_HA = {9}

STRATA_ORDER = ["N", "CA", "C", "O", "HN", "HA"]
PHYSICS_BLOCK_ORDER = ["ring", "charge", "mcconnell", "field_buckingham", "efg"]
AIMNET_BLOCK_ORDER = ["aimnet2_charge", "aimnet2_crg", "aimnet2_embedding"]
MC_CATEGORIES = {
    0: "PeptideCO",
    1: "PeptideCN",
    3: "SidechainCO",
}
TEST_FRAME_FRACTION = 0.20
THIN_ATOMS = 10


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--substrate-dir",
        default="/tmp/rediscover-corrected-backbone-snapshot-1p9j",
        help="root emitted substrate directory; subdirs can be overridden individually",
    )
    ap.add_argument("--broad-dir", default=None)
    ap.add_argument("--buckingham-dir", default=None)
    ap.add_argument("--efg-dir", default=None)
    ap.add_argument("--aimnet-dir", default=None)
    ap.add_argument(
        "--out-dir",
        default=None,
        help="artifact output directory; defaults to <substrate-dir>/analysis/ensemble_layer3",
    )
    ap.add_argument(
        "--report-md",
        default="src/rediscover/ENSEMBLE_LAYER3.md",
        help="markdown report path, relative to cwd unless absolute",
    )
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--ridge-alpha", type=float, default=10.0)
    ap.add_argument("--embedding-components", type=int, default=32)
    ap.add_argument("--within-split", choices=("random", "blocked"), default="random")
    ap.add_argument("--test-frame-fraction", type=float, default=TEST_FRAME_FRACTION)
    ap.add_argument("--mc-chunksize", type=int, default=500_000)
    return ap.parse_args()


def stratum_of(frame_variant: object) -> str | None:
    fv = int(frame_variant)
    if fv in FV_N:
        return "N"
    if fv in FV_CA:
        return "CA"
    if fv in FV_C:
        return "C"
    if fv in FV_O:
        return "O"
    if fv in FV_HN:
        return "HN"
    if fv in FV_HA:
        return "HA"
    return None


def require_file(path: Path) -> Path:
    if not path.exists():
        raise SystemExit(f"FATAL: required emitted substrate missing: {path}")
    return path


def require_columns(df: pd.DataFrame, cols: Iterable[str], label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(f"FATAL: {label} missing required columns: {missing}")


def finite_fmt(x: object, digits: int = 4) -> str:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(v):
        return "nan"
    av = abs(v)
    if av != 0.0 and (av < 1.0e-3 or av >= 1.0e4):
        return f"{v:.{digits}e}"
    return f"{v:.{digits}f}"


def r2_score(pred: np.ndarray, target: np.ndarray) -> float:
    pred = np.asarray(pred, dtype=float)
    target = np.asarray(target, dtype=float)
    if pred.ndim == 1:
        pred = pred.reshape(-1, 1)
    if target.ndim == 1:
        target = target.reshape(-1, 1)
    ok = np.isfinite(pred).all(axis=1) & np.isfinite(target).all(axis=1)
    if ok.sum() < 2:
        return math.nan
    p = pred[ok]
    y = target[ok]
    den = float(((y - y.mean(axis=0, keepdims=True)) ** 2).sum())
    if den <= 0.0:
        return math.nan
    return float(1.0 - ((y - p) ** 2).sum() / den)


def atom_means(values: np.ndarray, atoms: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values.reshape(-1, 1)
    labels = np.sort(np.unique(atoms))
    means = []
    counts = []
    for atom in labels:
        m = atoms == atom
        means.append(np.nanmean(values[m], axis=0))
        counts.append(int(m.sum()))
    return labels, np.asarray(means, dtype=float), np.asarray(counts, dtype=float)


def centred_by_train_atom(values: np.ndarray, atoms: np.ndarray, train: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values.reshape(-1, 1)
    out = np.full_like(values, np.nan, dtype=float)
    for atom in np.unique(atoms):
        m = atoms == atom
        mt = m & train
        if mt.sum() == 0:
            continue
        out[m] = values[m] - np.nanmean(values[mt], axis=0, keepdims=True)
    return out


def variance_shares(target: np.ndarray, atoms: np.ndarray) -> tuple[float, float]:
    y = np.asarray(target, dtype=float)
    if y.ndim == 1:
        y = y.reshape(-1, 1)
    ok = np.isfinite(y).all(axis=1)
    y = y[ok]
    a = np.asarray(atoms)[ok]
    if len(y) < 3:
        return math.nan, math.nan
    grand = y.mean(axis=0, keepdims=True)
    total = float(((y - grand) ** 2).sum())
    if total <= 0.0:
        return math.nan, math.nan
    between = 0.0
    within = 0.0
    for atom in np.unique(a):
        m = a == atom
        mean = y[m].mean(axis=0, keepdims=True)
        between += float(m.sum() * ((mean - grand) ** 2).sum())
        within += float(((y[m] - mean) ** 2).sum())
    return between / total, within / total


def effective_n_components(target: np.ndarray, atoms: np.ndarray, frames: np.ndarray) -> dict[str, float]:
    y = np.asarray(target, dtype=float)
    if y.ndim == 1:
        y = y.reshape(-1, 1)
    neff = np.zeros(y.shape[1], dtype=float)
    rhos = []
    for atom in np.unique(atoms):
        idx = np.flatnonzero(atoms == atom)
        idx = idx[np.argsort(frames[idx])]
        for c in range(y.shape[1]):
            vals = y[idx, c]
            vals = vals[np.isfinite(vals)]
            n = len(vals)
            if n < 3:
                neff[c] += float(n)
                continue
            yd = vals - vals.mean()
            den = float(np.dot(yd, yd))
            if den <= 0.0:
                rho = 0.0
                n_eff = float(n)
            else:
                rho = float(np.dot(yd[:-1], yd[1:]) / den)
                rho = max(min(rho, 0.999), -0.999)
                n_eff = n * (1.0 - rho) / (1.0 + rho)
                n_eff = min(max(n_eff, 1.0), float(n))
            neff[c] += n_eff
            rhos.append(rho)
    return {
        "within_N_eff_min": float(np.min(neff)) if len(neff) else math.nan,
        "within_N_eff_median": float(np.median(neff)) if len(neff) else math.nan,
        "within_N_eff_max": float(np.max(neff)) if len(neff) else math.nan,
        "median_lag1_rho": float(np.median(rhos)) if rhos else math.nan,
    }


def split_masks(
    frames: np.ndarray,
    strategy: str,
    seed: int,
    test_fraction: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    frames = np.asarray(frames)
    unique = np.sort(np.unique(frames))
    if len(unique) < 5:
        z = np.zeros(len(frames), dtype=bool)
        return z.copy(), z.copy(), {"test_frames": 0, "cross_split_lag1_pairs": 0}
    n_test = max(1, int(test_fraction * len(unique)))
    rng = np.random.default_rng(seed)
    if strategy == "random":
        test_frames = set(rng.choice(unique, n_test, replace=False))
    else:
        start = int(rng.integers(0, len(unique) - n_test + 1))
        test_frames = set(unique[start:start + n_test])
    train_frames = set(unique) - test_frames
    train = np.array([f in train_frames for f in frames], dtype=bool)
    test = np.array([f in test_frames for f in frames], dtype=bool)
    frame_split = {f: "test" if f in test_frames else "train" for f in unique}
    cross = 0
    for a, b in zip(unique[:-1], unique[1:]):
        if {frame_split[a], frame_split[b]} == {"train", "test"}:
            cross += 1
    return train, test, {"test_frames": int(len(test_frames)), "cross_split_lag1_pairs": int(cross)}


def fit_ridge_predict(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    alpha: float,
) -> np.ndarray:
    x_train = np.asarray(x_train, dtype=float)
    y_train = np.asarray(y_train, dtype=float)
    x_test = np.asarray(x_test, dtype=float)
    if y_train.ndim == 1:
        y_train = y_train.reshape(-1, 1)
    if x_train.ndim == 1:
        x_train = x_train.reshape(-1, 1)
    if x_test.ndim == 1:
        x_test = x_test.reshape(-1, 1)

    y_ok = np.isfinite(y_train).all(axis=1)
    x_ok = np.isfinite(x_train).any(axis=1)
    ok = y_ok & x_ok
    out_dim = y_train.shape[1]
    x_train = x_train[ok]
    y_train = y_train[ok]
    if len(x_train) < 2:
        return np.full((len(x_test), out_dim), np.nan)

    mean = np.nanmean(x_train, axis=0)
    mean[~np.isfinite(mean)] = 0.0
    x_fill = np.where(np.isfinite(x_train), x_train, mean)
    std = x_fill.std(axis=0)
    std[~np.isfinite(std) | (std <= 0.0)] = 1.0
    xz = (x_fill - mean) / std
    xt = np.where(np.isfinite(x_test), x_test, mean)
    xtz = (xt - mean) / std

    design = np.column_stack([np.ones(len(xz)), xz])
    design_t = np.column_stack([np.ones(len(xtz)), xtz])
    penalty = np.eye(design.shape[1]) * max(0.0, float(alpha))
    penalty[0, 0] = 0.0
    beta = np.linalg.solve(design.T @ design + penalty, design.T @ y_train)
    return design_t @ beta


def jackknife_metric_values(pred: np.ndarray, target: np.ndarray, groups: np.ndarray) -> np.ndarray:
    groups = np.asarray(groups)
    labels = np.unique(groups)
    if len(labels) < 3:
        return np.asarray([], dtype=float)
    vals = []
    for label in labels:
        keep = groups != label
        if keep.sum() < 3:
            continue
        vals.append(r2_score(np.asarray(pred)[keep], np.asarray(target)[keep]))
    vals = np.asarray(vals, dtype=float)
    return vals[np.isfinite(vals)]


def jackknife_se_from_values(vals: np.ndarray) -> float:
    vals = np.asarray(vals, dtype=float)
    if len(vals) < 3:
        return math.nan
    mean = float(vals.mean())
    return float(math.sqrt((len(vals) - 1) / len(vals) * np.sum((vals - mean) ** 2)))


def jackknife_metric_se(pred: np.ndarray, target: np.ndarray, groups: np.ndarray) -> float:
    return jackknife_se_from_values(jackknife_metric_values(pred, target, groups))


def jackknife_difference_se(
    pred: np.ndarray,
    target: np.ndarray,
    groups: np.ndarray,
    drop_pred: np.ndarray,
    drop_target: np.ndarray,
    drop_groups: np.ndarray,
) -> float:
    pred = np.asarray(pred, dtype=float)
    target = np.asarray(target, dtype=float)
    groups = np.asarray(groups)
    drop_pred = np.asarray(drop_pred, dtype=float)
    drop_target = np.asarray(drop_target, dtype=float)
    drop_groups = np.asarray(drop_groups)
    if (
        len(pred) != len(drop_pred)
        or target.shape != drop_target.shape
        or len(groups) != len(drop_groups)
        or not np.array_equal(groups, drop_groups)
        or not np.allclose(target, drop_target, equal_nan=True)
    ):
        return math.nan
    labels = np.unique(groups)
    vals = []
    for label in labels:
        keep = groups != label
        if keep.sum() < 3:
            continue
        vals.append(
            r2_score(pred[keep], target[keep])
            - r2_score(drop_pred[keep], drop_target[keep])
        )
    vals = np.asarray(vals, dtype=float)
    return jackknife_se_from_values(vals[np.isfinite(vals)])


def evaluate_between(x: np.ndarray, y: np.ndarray, atoms: np.ndarray, alpha: float) -> dict[str, object]:
    labels, xbar, _ = atom_means(x, atoms)
    _, ybar, _ = atom_means(y, atoms)
    pred = np.full_like(ybar, np.nan, dtype=float)
    for i in range(len(labels)):
        train = np.ones(len(labels), dtype=bool)
        train[i] = False
        pred[i] = fit_ridge_predict(xbar[train], ybar[train], xbar[i:i + 1], alpha)[0]
    return {
        "R2": r2_score(pred, ybar),
        "jackknife_se": jackknife_metric_se(pred, ybar, labels),
        "pred": pred,
        "target": ybar,
        "groups": labels,
    }


def evaluate_within(
    x: np.ndarray,
    y: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    alpha: float,
    split_strategy: str,
    seed: int,
    test_fraction: float,
) -> dict[str, object]:
    train, test, split = split_masks(frames, split_strategy, seed, test_fraction)
    xw = centred_by_train_atom(x, atoms, train)
    yw = centred_by_train_atom(y, atoms, train)
    fit = train & np.isfinite(yw).all(axis=1) & np.isfinite(xw).any(axis=1)
    score = test & np.isfinite(yw).all(axis=1) & np.isfinite(xw).any(axis=1)
    if fit.sum() < 3 or score.sum() < 3:
        pred = np.full((score.sum(), yw.shape[1]), np.nan)
    else:
        pred = fit_ridge_predict(xw[fit], yw[fit], xw[score], alpha)
    target = yw[score]
    groups = np.asarray(atoms)[score]
    return {
        "R2": r2_score(pred, target),
        "jackknife_se": jackknife_metric_se(pred, target, groups),
        "pred": pred,
        "target": target,
        "groups": groups,
        **split,
    }


def reduce_embedding(embedding: np.ndarray, n_components: int) -> tuple[np.ndarray, dict[str, object]]:
    embedding = np.asarray(embedding, dtype=float)
    original_dims = embedding.shape[1]
    if n_components <= 0 or n_components >= original_dims:
        return embedding, {
            "method": "identity_full_embedding",
            "original_dims": int(original_dims),
            "used_dims": int(original_dims),
            "explained_variance_ratio": None,
        }
    mean = np.nanmean(embedding, axis=0)
    mean[~np.isfinite(mean)] = 0.0
    x = np.where(np.isfinite(embedding), embedding, mean) - mean
    cov = (x.T @ x) / max(1, x.shape[0] - 1)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    k = min(int(n_components), original_dims)
    total = float(np.sum(vals[vals > 0.0]))
    ratio = float(np.sum(vals[:k]) / total) if total > 0.0 else math.nan
    return x @ vecs[:, :k], {
        "method": "global_unsupervised_pca_on_emitted_embedding",
        "original_dims": int(original_dims),
        "used_dims": int(k),
        "explained_variance_ratio": ratio,
    }


def t2_sidecar(directory: Path, name: str, C: np.ndarray) -> np.ndarray:
    return np.load(require_file(directory / name)) @ C.T


def source_sum_valid_mcconnell_categories(
    broad_dir: Path,
    n_rows: int,
    C: np.ndarray,
    chunksize: int,
) -> tuple[np.ndarray, dict[str, object]]:
    source_path = require_file(broad_dir / "broad_backbone_sources.csv")
    t2_cols = [f"mc_lit_T2_local_{i}" for i in range(5)]
    usecols = ["row_id", "bond_category", "mc_lit_kernel_present", "mc_source_is_self_or_bonded", *t2_cols]
    out = np.zeros((n_rows, len(MC_CATEGORIES), 5), dtype=float)
    cat_to_col = {cat: i for i, cat in enumerate(MC_CATEGORIES)}
    counts = {MC_CATEGORIES[cat]: 0 for cat in MC_CATEGORIES}
    filtered_self_or_bonded = {MC_CATEGORIES[cat]: 0 for cat in MC_CATEGORIES}
    source_rows = 0
    present_rows = 0
    for chunk in pd.read_csv(source_path, usecols=usecols, chunksize=chunksize):
        source_rows += len(chunk)
        present = chunk["mc_lit_kernel_present"].to_numpy(int) == 1
        valid = chunk["mc_source_is_self_or_bonded"].to_numpy(int) == 0
        present_rows += int(present.sum())
        categories = chunk["bond_category"].to_numpy(int)
        for cat, col in cat_to_col.items():
            cat_mask = present & (categories == cat)
            filtered_self_or_bonded[MC_CATEGORIES[cat]] += int((cat_mask & ~valid).sum())
            mask = cat_mask & valid
            if not mask.any():
                continue
            rows = chunk.loc[mask, "row_id"].to_numpy(int)
            if rows.min(initial=0) < 0 or rows.max(initial=0) >= n_rows:
                raise SystemExit(f"FATAL: source row_id outside aggregate range 0..{n_rows - 1}")
            vals = chunk.loc[mask, t2_cols].to_numpy(float)
            counts[MC_CATEGORIES[cat]] += int(mask.sum())
            for k in range(5):
                np.add.at(out[:, col, k], rows, vals[:, k])
    flat = out.reshape(-1, 5) @ C.T
    out = flat.reshape(n_rows, len(MC_CATEGORIES), 5)
    return out, {
        "source_csv": str(source_path),
        "source_rows": int(source_rows),
        "present_source_rows": int(present_rows),
        "valid_category_counts": counts,
        "filtered_self_or_bonded_counts": filtered_self_or_bonded,
        "category_order": [MC_CATEGORIES[cat] for cat in MC_CATEGORIES],
    }


def load_joined(args: argparse.Namespace, C: np.ndarray) -> tuple[pd.DataFrame, dict[str, np.ndarray], dict[str, object]]:
    substrate_dir = Path(args.substrate_dir)
    broad_dir = Path(args.broad_dir) if args.broad_dir else substrate_dir / "broad_backbone"
    buckingham_dir = Path(args.buckingham_dir) if args.buckingham_dir else substrate_dir / "buckingham_efield"
    efg_dir = Path(args.efg_dir) if args.efg_dir else substrate_dir / "efg"
    aimnet_dir = Path(args.aimnet_dir) if args.aimnet_dir else substrate_dir / "aimnet2_features"

    broad = pd.read_csv(require_file(broad_dir / "broad_backbone_aggregated.csv"))
    broad["broad_row_id"] = np.arange(len(broad), dtype=int)
    broad["stratum"] = [stratum_of(v) for v in broad["frame_variant"]]
    require_columns(
        broad,
        [
            "atom_index",
            "h5_row",
            "frame_variant",
            "frame_valid",
            "dft_present",
            "dft_sigma_iso",
            "dft_local_frame_valid",
            "ring_literature_kernel_present",
            "charge_literature_kernel_present",
            "mc_lit_valid_kernel_present",
            "field_z",
            "field_mag",
            "bond_cutoff_A",
            "mc_source_near_field_ratio",
            *[f"mc_lit_T2_local_valid_{i}" for i in range(5)],
        ],
        "broad backbone aggregate",
    )

    buck = pd.read_csv(require_file(buckingham_dir / "buckingham_efield_aggregated.csv"))
    require_columns(
        buck,
        ["atom_index", "h5_row", "apbs_efield_present", "E_proj", "E_mag"],
        "buckingham efield aggregate",
    )
    buck = buck[["atom_index", "h5_row", "apbs_efield_present", "E_proj", "E_mag"]].rename(
        columns={
            "apbs_efield_present": "buckingham_present",
            "E_proj": "buckingham_E_proj",
            "E_mag": "buckingham_E_mag",
        }
    )

    efg = pd.read_csv(require_file(efg_dir / "efg_aggregated.csv"))
    require_columns(
        efg,
        ["row_id", "atom_index", "h5_row", "apbs_efg_present", "efg_T2_magnitude"],
        "APBS EFG aggregate",
    )
    efg = efg[["row_id", "atom_index", "h5_row", "apbs_efg_present", "efg_T2_magnitude"]].rename(
        columns={"row_id": "efg_row_id", "apbs_efg_present": "efg_present"}
    )

    aim = pd.read_csv(require_file(aimnet_dir / "aimnet2_features_aggregated.csv"))
    require_columns(
        aim,
        [
            "row_id",
            "atom_index",
            "h5_row",
            "aimnet2_charge_present",
            "aimnet2_charge",
            "aimnet2_charge_response_gradient_present",
            "aimnet2_charge_response_gradient_scalar",
            "aimnet2_charge_response_gradient_local_x",
            "aimnet2_charge_response_gradient_local_y",
            "aimnet2_charge_response_gradient_local_z",
            "aimnet2_embedding_present",
            "aimnet2_embedding_dim",
        ],
        "AIMNet2 aggregate",
    )
    aim = aim.rename(columns={"row_id": "aimnet2_row_id"})
    aim = aim[
        [
            "atom_index",
            "h5_row",
            "aimnet2_row_id",
            "aimnet2_charge_present",
            "aimnet2_charge",
            "aimnet2_charge_response_gradient_present",
            "aimnet2_charge_response_gradient_scalar",
            "aimnet2_charge_response_gradient_local_x",
            "aimnet2_charge_response_gradient_local_y",
            "aimnet2_charge_response_gradient_local_z",
            "aimnet2_embedding_present",
            "aimnet2_embedding_dim",
        ]
    ]

    df = broad.merge(buck, on=["atom_index", "h5_row"], how="left", sort=False)
    df = df.merge(efg, on=["atom_index", "h5_row"], how="left", sort=False)
    df = df.merge(aim, on=["atom_index", "h5_row"], how="left", sort=False)
    if len(df) != len(broad):
        raise SystemExit(f"FATAL: join changed row count: broad={len(broad)} joined={len(df)}")
    if df["aimnet2_row_id"].isna().any():
        raise SystemExit("FATAL: AIMNet2 feature export does not cover every broad row")
    if df["efg_row_id"].isna().any():
        raise SystemExit("FATAL: APBS EFG export does not cover every broad row")

    row_idx = df["broad_row_id"].to_numpy(int)
    efg_idx = df["efg_row_id"].to_numpy(int)
    aim_idx = df["aimnet2_row_id"].to_numpy(int)
    df["field_mag_sq"] = df["field_mag"].to_numpy(float) ** 2
    df["buckingham_E_mag_sq"] = df["buckingham_E_mag"].to_numpy(float) ** 2

    mc_cat, mc_audit = source_sum_valid_mcconnell_categories(broad_dir, len(df), C, args.mc_chunksize)
    mc_valid_cols = df[[f"mc_lit_T2_local_valid_{i}" for i in range(5)]].to_numpy(float) @ C.T
    mc_total = mc_cat.sum(axis=1)
    mc_agg_diff = float(np.nanmax(np.abs(mc_total - mc_valid_cols)))

    arrays = {
        "target_T2": t2_sidecar(broad_dir, "broad_backbone_aggregated_target_local_T2.npy", C)[row_idx],
        "ring_T2": t2_sidecar(broad_dir, "broad_backbone_aggregated_ring_literature_kernel_T2.npy", C)[row_idx],
        "charge_T2": t2_sidecar(broad_dir, "broad_backbone_aggregated_charge_literature_kernel_T2.npy", C)[row_idx],
        "mcconnell_cat_T2": mc_cat[row_idx],
        "efg_T2": t2_sidecar(efg_dir, "efg_feature_T2.npy", C)[efg_idx],
        "embedding": np.load(require_file(aimnet_dir / "aimnet2_features_embedding.npy"))[aim_idx],
        "crg_local": np.load(require_file(aimnet_dir / "aimnet2_features_charge_response_gradient_local.npy"))[aim_idx],
        "crg_scalar": np.load(require_file(aimnet_dir / "aimnet2_features_charge_response_gradient_scalar.npy"))[
            aim_idx
        ].reshape(-1, 1),
    }
    if arrays["embedding"].shape[1] != 256:
        raise SystemExit(f"FATAL: AIMNet2 embedding sidecar has shape {arrays['embedding'].shape}")

    audit = {
        "substrate_dir": str(substrate_dir),
        "broad_dir": str(broad_dir),
        "buckingham_dir": str(buckingham_dir),
        "efg_dir": str(efg_dir),
        "aimnet_dir": str(aimnet_dir),
        "broad_rows": int(len(df)),
        "bond_cutoff_A_unique": sorted({float(x) for x in df["bond_cutoff_A"].dropna().unique()}),
        "mc_source_near_field_ratio_unique": sorted(
            {float(x) for x in df["mc_source_near_field_ratio"].dropna().unique()}
        ),
        "mcconnell_per_category": mc_audit,
        "mcconnell_category_sum_vs_valid_aggregate_max_abs_diff_e3nn": mc_agg_diff,
        "python_reads": [
            "emitted CSV/NPY substrate",
            "broad_backbone_sources.csv emitted McConnell T2 columns only",
        ],
        "python_does_not": ["open trajectory.h5", "run ORCA", "call C++ emitter", "recompute kernels"],
    }
    return df, arrays, audit


def block_features(df: pd.DataFrame, arrays: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], dict[str, list[str]]]:
    mc = arrays["mcconnell_cat_T2"].reshape(len(df), len(MC_CATEGORIES) * 5)
    mc_labels = []
    for cat in MC_CATEGORIES:
        name = MC_CATEGORIES[cat]
        mc_labels += [f"mcconnell:{name}:T2_{i}" for i in range(5)]
    blocks = {
        "ring": arrays["ring_T2"],
        "charge": arrays["charge_T2"],
        "mcconnell": mc,
        "field_buckingham": df[
            ["field_z", "field_mag_sq", "buckingham_E_proj", "buckingham_E_mag_sq"]
        ].to_numpy(float),
        "efg": arrays["efg_T2"],
        "aimnet2_charge": df[["aimnet2_charge"]].to_numpy(float),
        "aimnet2_crg": np.column_stack([arrays["crg_local"], arrays["crg_scalar"]]),
        "aimnet2_embedding": arrays["embedding"],
    }
    labels = {
        "ring": [f"ring:literature_T2_{i}" for i in range(5)],
        "charge": [f"charge:q_over_r3_T2_{i}" for i in range(5)],
        "mcconnell": mc_labels,
        "field_buckingham": [
            "field:ff14sb_E_z",
            "field:ff14sb_absE_sq",
            "buckingham:apbs_E_proj",
            "buckingham:apbs_absE_sq",
        ],
        "efg": [f"efg:apbs_local_T2_{i}" for i in range(5)],
        "aimnet2_charge": ["aimnet2:charge"],
        "aimnet2_crg": [
            "aimnet2:charge_response_gradient_local_x",
            "aimnet2:charge_response_gradient_local_y",
            "aimnet2:charge_response_gradient_local_z",
            "aimnet2:charge_response_gradient_scalar",
        ],
        "aimnet2_embedding": [f"aimnet2:embedding_pc_{i:03d}" for i in range(arrays["embedding"].shape[1])],
    }
    return blocks, labels


def concat_blocks(
    blocks: dict[str, np.ndarray],
    labels: dict[str, list[str]],
    names: list[str],
) -> tuple[np.ndarray, list[str]]:
    return np.column_stack([blocks[name] for name in names]), [label for name in names for label in labels[name]]


def result_row(
    target: str,
    stratum: str,
    feature_scope: str,
    x: np.ndarray,
    y: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    alpha: float,
    split_strategy: str,
    seed: int,
    test_fraction: float,
) -> dict[str, object]:
    shares = variance_shares(y, atoms)
    neff = effective_n_components(y, atoms, frames)
    between = evaluate_between(x, y, atoms, alpha)
    within = evaluate_within(x, y, atoms, frames, alpha, split_strategy, seed, test_fraction)
    return {
        "target": target,
        "stratum": stratum,
        "feature_scope": feature_scope,
        "rows": int(len(y)),
        "n_atoms_between": int(len(np.unique(atoms))),
        "n_features": int(x.shape[1]),
        "ridge_alpha": float(alpha),
        "variance_share_between": shares[0],
        "variance_share_within": shares[1],
        **neff,
        "thin_flag": "THIN" if len(np.unique(atoms)) < THIN_ATOMS else "",
        "between_LOAO_R2": between["R2"],
        "between_LOAO_R2_jackknife_se": between["jackknife_se"],
        "within_frame_R2": within["R2"],
        "within_frame_R2_jackknife_se": within["jackknife_se"],
        "within_split_strategy": split_strategy,
        "test_frames": int(within.get("test_frames", 0)),
        "cross_split_lag1_pairs": int(within.get("cross_split_lag1_pairs", 0)),
        "_between_pred": between["pred"],
        "_between_target": between["target"],
        "_between_groups": between["groups"],
        "_within_pred": within["pred"],
        "_within_target": within["target"],
        "_within_groups": within["groups"],
    }


def strip_private(row: dict[str, object]) -> dict[str, object]:
    return {k: v for k, v in row.items() if not k.startswith("_")}


def evaluate_drop_one(
    full: dict[str, object],
    target: str,
    stratum: str,
    feature_scope: str,
    block_name: str,
    kept_blocks: list[str],
    blocks: dict[str, np.ndarray],
    labels: dict[str, list[str]],
    y: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    alpha: float,
    split_strategy: str,
    seed: int,
    test_fraction: float,
) -> dict[str, object]:
    x_drop, _ = concat_blocks(blocks, labels, kept_blocks)
    drop = result_row(
        target,
        stratum,
        f"{feature_scope}_drop_{block_name}",
        x_drop,
        y,
        atoms,
        frames,
        alpha,
        split_strategy,
        seed,
        test_fraction,
    )
    between_delta = float(full["between_LOAO_R2"]) - float(drop["between_LOAO_R2"])
    within_delta = float(full["within_frame_R2"]) - float(drop["within_frame_R2"])
    return {
        "target": target,
        "stratum": stratum,
        "feature_scope": feature_scope,
        "block": block_name,
        "rows": full["rows"],
        "n_atoms_between": full["n_atoms_between"],
        "ridge_alpha": alpha,
        "full_between_LOAO_R2": full["between_LOAO_R2"],
        "drop_between_LOAO_R2": drop["between_LOAO_R2"],
        "delta_between_LOAO_R2": between_delta,
        "delta_between_LOAO_R2_jackknife_se": jackknife_difference_se(
            full["_between_pred"],
            full["_between_target"],
            full["_between_groups"],
            drop["_between_pred"],
            drop["_between_target"],
            drop["_between_groups"],
        ),
        "full_within_frame_R2": full["within_frame_R2"],
        "drop_within_frame_R2": drop["within_frame_R2"],
        "delta_within_frame_R2": within_delta,
        "delta_within_frame_R2_jackknife_se": jackknife_difference_se(
            full["_within_pred"],
            full["_within_target"],
            full["_within_groups"],
            drop["_within_pred"],
            drop["_within_target"],
            drop["_within_groups"],
        ),
    }


def add_contribution_shares(contrib: pd.DataFrame) -> pd.DataFrame:
    if contrib.empty:
        return contrib
    out = contrib.copy()
    out["positive_delta_between"] = np.maximum(out["delta_between_LOAO_R2"].to_numpy(float), 0.0)
    out["positive_delta_within"] = np.maximum(out["delta_within_frame_R2"].to_numpy(float), 0.0)
    sums = out.groupby(["target", "stratum", "feature_scope"], sort=False)[
        ["positive_delta_between", "positive_delta_within"]
    ].transform("sum")
    out["share_between_positive_drop"] = np.divide(
        out["positive_delta_between"],
        sums["positive_delta_between"],
        out=np.zeros(len(out), dtype=float),
        where=sums["positive_delta_between"].to_numpy(float) > 0.0,
    )
    out["share_within_positive_drop"] = np.divide(
        out["positive_delta_within"],
        sums["positive_delta_within"],
        out=np.zeros(len(out), dtype=float),
        where=sums["positive_delta_within"].to_numpy(float) > 0.0,
    )
    return out


def verdict_from_r2(value: float, se: float) -> str:
    if not np.isfinite(value):
        return "can't-make-it-work-for-now"
    if value >= 0.25 and (not np.isfinite(se) or value > se):
        return "form-recovered-scale-fitted"
    if value >= 0.10 and (not np.isfinite(se) or value > 0.5 * se):
        return "form-recovered-scale-fitted"
    return "can't-make-it-work-for-now"


def mc_threshold_positive(delta: float, se: float) -> bool:
    return bool(np.isfinite(delta) and delta > 0.02 and (not np.isfinite(se) or delta > se))


def mc_verdict_bucket(delta: float, se: float) -> str:
    if mc_threshold_positive(delta, se):
        return "form-recovered-scale-fitted"
    return "can't-make-it-work-for-now"


def mc_add_read(delta: float, se: float) -> str:
    if np.isfinite(delta) and delta > 0.02 and (not np.isfinite(se) or delta > se):
        return "threshold-positive"
    if np.isfinite(delta) and delta > 0.0:
        return "weak-positive"
    return "no-positive-drop-one-lift"


def pivot_metric(df: pd.DataFrame, feature_scope: str, metric: str) -> pd.DataFrame:
    sub = df[df["feature_scope"] == feature_scope][["target", "stratum", metric]]
    return sub.rename(columns={metric: feature_scope})


def md_table(df: pd.DataFrame, cols: list[str]) -> list[str]:
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in df[cols].iterrows():
        vals = []
        for col in cols:
            v = row[col]
            if isinstance(v, (float, np.floating)):
                vals.append(finite_fmt(v))
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    return lines


def write_report(
    results: pd.DataFrame,
    contributions: pd.DataFrame,
    audit: dict[str, object],
    run: dict[str, object],
    path: Path,
) -> None:
    phys = results[results["feature_scope"] == "physics_only"].copy()
    plus = results[results["feature_scope"] == "physics_plus_aimnet2"].copy()
    joined = phys.merge(
        plus,
        on=["target", "stratum"],
        suffixes=("_physics", "_plus_aimnet2"),
        how="inner",
    )
    joined["aimnet2_gap_between_R2"] = (
        joined["between_LOAO_R2_plus_aimnet2"] - joined["between_LOAO_R2_physics"]
    )
    joined["aimnet2_gap_within_R2"] = (
        joined["within_frame_R2_plus_aimnet2"] - joined["within_frame_R2_physics"]
    )
    joined["physics_verdict"] = [
        verdict_from_r2(r2, se)
        for r2, se in zip(joined["between_LOAO_R2_physics"], joined["between_LOAO_R2_jackknife_se_physics"])
    ]
    joined["plus_aimnet2_verdict"] = [
        verdict_from_r2(r2, se)
        for r2, se in zip(
            joined["between_LOAO_R2_plus_aimnet2"],
            joined["between_LOAO_R2_jackknife_se_plus_aimnet2"],
        )
    ]

    physics_contrib = contributions[contributions["feature_scope"] == "physics_only"].copy()
    top = physics_contrib.copy()
    top["rank_value"] = top["positive_delta_between"]
    idx = top.groupby(["target", "stratum"])["rank_value"].idxmax()
    top = top.loc[idx, ["target", "stratum", "block", "delta_between_LOAO_R2", "share_between_positive_drop"]]
    top = top.rename(
        columns={
            "block": "top_physics_block_between",
            "delta_between_LOAO_R2": "top_delta_between_R2",
            "share_between_positive_drop": "top_share_between",
        }
    )
    mc = physics_contrib[physics_contrib["block"] == "mcconnell"][
        [
            "target",
            "stratum",
            "delta_between_LOAO_R2",
            "delta_between_LOAO_R2_jackknife_se",
            "share_between_positive_drop",
        ]
    ].copy()
    mc["mcconnell_joint_add_read"] = [
        mc_add_read(delta, se)
        for delta, se in zip(mc["delta_between_LOAO_R2"], mc["delta_between_LOAO_R2_jackknife_se"])
    ]
    mc["mcconnell_verdict_bucket"] = [
        mc_verdict_bucket(delta, se)
        for delta, se in zip(mc["delta_between_LOAO_R2"], mc["delta_between_LOAO_R2_jackknife_se"])
    ]
    mc = mc.rename(
        columns={
            "delta_between_LOAO_R2": "mcconnell_delta_between_R2",
            "delta_between_LOAO_R2_jackknife_se": "mcconnell_delta_between_R2_se",
            "share_between_positive_drop": "mcconnell_share_between",
        }
    )
    summary = joined.merge(top, on=["target", "stratum"], how="left").merge(
        mc, on=["target", "stratum"], how="left"
    )
    summary = summary.sort_values(
        ["target", "stratum"],
        key=lambda s: s.map({name: i for i, name in enumerate(STRATA_ORDER)}).fillna(s)
        if s.name == "stratum"
        else s,
    )

    share_rows = []
    for (target, stratum), group in physics_contrib.groupby(["target", "stratum"], sort=False):
        row = {"target": target, "stratum": stratum}
        for block in PHYSICS_BLOCK_ORDER:
            b = group[group["block"] == block]
            if b.empty:
                row[block] = "nan"
            else:
                first = b.iloc[0]
                row[block] = (
                    f"{finite_fmt(first['share_between_positive_drop'])} "
                    f"(dR2 {finite_fmt(first['delta_between_LOAO_R2'])})"
                )
        share_rows.append(row)
    share_df = pd.DataFrame(share_rows)
    if not share_df.empty:
        share_df = share_df.sort_values(
            ["target", "stratum"],
            key=lambda s: s.map({name: i for i, name in enumerate(STRATA_ORDER)}).fillna(s)
            if s.name == "stratum"
            else s,
        )

    best_physics = joined.loc[joined["between_LOAO_R2_physics"].idxmax()]
    best_plus = joined.loc[joined["between_LOAO_R2_plus_aimnet2"].idxmax()]
    positive_physics = int((joined["between_LOAO_R2_physics"] > 0.0).sum())
    positive_plus = int((joined["between_LOAO_R2_plus_aimnet2"] > 0.0).sum())
    mc_add_count = int((mc["mcconnell_joint_add_read"] == "threshold-positive").sum())

    lines = [
        "# Ensemble Layer 3",
        "",
        "Joint per-stratum ridge on the corrected emitted substrate. This is a within-protein "
        "signal-capture diagnostic, not a population prediction claim.",
        "",
        "## Lead",
        "",
        (
            f"- Physics-only between-atom LOAO R2 is positive in {positive_physics}/12 target-stratum "
            f"rows. Best physics-only row: {best_physics['target']} {best_physics['stratum']} "
            f"R2={finite_fmt(best_physics['between_LOAO_R2_physics'])}."
        ),
        (
            f"- Physics + AIMNet2 is positive in {positive_plus}/12 rows. Best ceiling row: "
            f"{best_plus['target']} {best_plus['stratum']} "
            f"R2={finite_fmt(best_plus['between_LOAO_R2_plus_aimnet2'])}."
        ),
        (
            f"- McConnell's joint drop-one result: {mc_add_count} threshold-positive rows "
            "across the 12 between-led tests. "
            "This is the fair joint test; standalone McConnell remains a weaker read."
        ),
        "- Verdict bucket for the layer-3 physics ensemble: form-recovered-scale-fitted where "
        "between LOAO R2 is positive and stable; otherwise can't-make-it-work-for-now.",
        "",
        "## Joint LOAO R2",
        "",
        "Headline is between-atom LOAO R2. Within-frame R2 is a train/test frame diagnostic after "
        "train-only per-atom centering. Static sigma_iso is between-led.",
        "",
    ]

    joint_cols = [
        "target",
        "stratum",
        "n_atoms_between_physics",
        "variance_share_between_physics",
        "between_LOAO_R2_physics",
        "between_LOAO_R2_jackknife_se_physics",
        "between_LOAO_R2_plus_aimnet2",
        "between_LOAO_R2_jackknife_se_plus_aimnet2",
        "aimnet2_gap_between_R2",
        "within_frame_R2_physics",
        "within_frame_R2_plus_aimnet2",
        "within_N_eff_median_physics",
        "median_lag1_rho_physics",
        "thin_flag_physics",
        "physics_verdict",
        "plus_aimnet2_verdict",
    ]
    lines.extend(md_table(summary, joint_cols))
    lines.extend(
        [
            "",
            "## Physics Contribution Shares",
            "",
            "Shares are drop-one positive-delta shares within the physics-only joint fit. "
            "Each cell is `share (dR2 delta)`, so negative or redundant blocks can have zero share.",
            "",
        ]
    )
    if not share_df.empty:
        lines.extend(md_table(share_df, ["target", "stratum", *PHYSICS_BLOCK_ORDER]))
    lines.extend(
        [
            "",
            "## McConnell Verdict",
            "",
            "Drop-one values ask whether valid-source per-category McConnell improves the already-joint "
            "physics fit. The form-recovered-scale-fitted bucket requires positive dR2 above 0.02 and "
            "above the delete-atom jackknife SE.",
            "",
        ]
    )
    lines.extend(
        md_table(
            mc.assign(stratum_order=mc["stratum"].map({name: i for i, name in enumerate(STRATA_ORDER)}))
            .sort_values(["target", "stratum_order"])
            .reset_index(drop=True),
            [
                "target",
                "stratum",
                "mcconnell_delta_between_R2",
                "mcconnell_delta_between_R2_se",
                "mcconnell_share_between",
                "mcconnell_verdict_bucket",
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Honest Bottom Line",
            "",
            "The recovered kernels jointly capture real backbone signal only where the between-atom "
            "LOAO rows are positive with tolerable jackknife uncertainty. Ring, charge, field/Buckingham, "
            "and EFG compete in the joint basis; McConnell gets its fair valid-source per-category test "
            "as a component rather than as a standalone total-target fit. AIMNet2 remains a ceiling: "
            "embedding PCs, AIMNet2 charge, and CRG are learnable upper-bound features, not recovered laws.",
            "",
            "Ring and charge remain visible but are not universal joint drop-one leaders in this "
            "one-protein layer-3 fit. Charge is clearest on HA/HN T2-style rows, ring appears in HA T2 "
            "and CA sigma_iso, while McConnell, field/Buckingham, and EFG carry several between-led "
            "deltas once all kernels compete in the same ridge.",
            "",
            "This is n=1 protein evidence. The report is per stratum, correlate-not-match, and does not "
            "make population claims.",
            "",
            "## Self-Audit",
            "",
            f"- Substrate: `{audit['substrate_dir']}`.",
            f"- Broad bond cutoff unique values: `{audit['bond_cutoff_A_unique']}`.",
            f"- McConnell near-field ratio unique values: `{audit['mc_source_near_field_ratio_unique']}`.",
            f"- Valid McConnell category counts: `{audit['mcconnell_per_category']['valid_category_counts']}`.",
            "- Python read emitted CSV/NPY features only; it did not open `trajectory.h5`, run ORCA, "
            "call the C++ emitter, or recompute kernels.",
            f"- Frozen basis audit `|C.T C - I|max`: {finite_fmt(run['change_of_basis_orthogonality_max'], 6)}.",
            f"- Ridge alpha: {finite_fmt(run['ridge_alpha'])}; embedding PCs: "
            f"{run['embedding_transform']['used_dims']} of {run['embedding_transform']['original_dims']}.",
            f"- Artifacts: `{run['out_dir']}`.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    substrate_dir = Path(args.substrate_dir)
    out_dir = Path(args.out_dir) if args.out_dir else substrate_dir / "analysis" / "ensemble_layer3"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_md = Path(args.report_md)

    C = cob.get_C()
    orth = float(np.abs(C.T @ C - np.eye(5)).max())
    print(f"change_of_basis.get_C() |C^T C - I|max={orth:.3e}")
    df, arrays, audit = load_joined(args, C)
    arrays["embedding"], embedding_transform = reduce_embedding(arrays["embedding"], args.embedding_components)

    all_blocks, all_labels = block_features(df, arrays)
    feature_label_payload = {name: all_labels[name] for name in all_labels}
    results: list[dict[str, object]] = []
    contributions: list[dict[str, object]] = []

    for stratum in STRATA_ORDER:
        base_mask = (
            (df["stratum"].to_numpy() == stratum)
            & (df["dft_present"].to_numpy(int) == 1)
            & (df["frame_valid"].to_numpy(int) == 1)
            & (df["dft_local_frame_valid"].to_numpy(int) == 1)
            & (df["ring_literature_kernel_present"].to_numpy(int) == 1)
            & (df["charge_literature_kernel_present"].to_numpy(int) == 1)
            & (df["mc_lit_valid_kernel_present"].to_numpy(int) == 1)
            & (df["buckingham_present"].fillna(0).to_numpy(int) == 1)
            & (df["efg_present"].fillna(0).to_numpy(int) == 1)
            & (df["aimnet2_charge_present"].to_numpy(int) == 1)
            & (df["aimnet2_charge_response_gradient_present"].to_numpy(int) == 1)
            & (df["aimnet2_embedding_present"].to_numpy(int) == 1)
        )
        idx = np.flatnonzero(base_mask)
        if len(idx) < 4:
            continue
        d = df.iloc[idx].reset_index(drop=True)
        atoms = d["atom_index"].to_numpy(int)
        frames = d["h5_row"].to_numpy(int)
        local_blocks = {name: value[idx] for name, value in all_blocks.items()}

        scopes = {
            "physics_only": PHYSICS_BLOCK_ORDER,
            "physics_plus_aimnet2": PHYSICS_BLOCK_ORDER + AIMNET_BLOCK_ORDER,
        }
        targets = {
            "T2": arrays["target_T2"][idx],
            "sigma_iso": d["dft_sigma_iso"].to_numpy(float).reshape(-1, 1),
        }
        for target_name, y in targets.items():
            finite_target = np.isfinite(y).all(axis=1)
            if finite_target.sum() < 4:
                continue
            y_fit = y[finite_target]
            atoms_fit = atoms[finite_target]
            frames_fit = frames[finite_target]
            scoped_full: dict[str, dict[str, object]] = {}
            scoped_blocks: dict[str, dict[str, np.ndarray]] = {}
            for scope_name, block_names in scopes.items():
                local = {name: local_blocks[name][finite_target] for name in block_names}
                x, labels = concat_blocks(local, all_labels, block_names)
                full = result_row(
                    target_name,
                    stratum,
                    scope_name,
                    x,
                    y_fit,
                    atoms_fit,
                    frames_fit,
                    args.ridge_alpha,
                    args.within_split,
                    args.seed,
                    args.test_frame_fraction,
                )
                full["feature_blocks"] = "+".join(block_names)
                results.append(full)
                scoped_full[scope_name] = full
                scoped_blocks[scope_name] = local
                for block in block_names:
                    kept = [name for name in block_names if name != block]
                    contributions.append(
                        evaluate_drop_one(
                            full,
                            target_name,
                            stratum,
                            scope_name,
                            block,
                            kept,
                            local,
                            all_labels,
                            y_fit,
                            atoms_fit,
                            frames_fit,
                            args.ridge_alpha,
                            args.within_split,
                            args.seed,
                            args.test_frame_fraction,
                        )
                    )

    if not results:
        raise SystemExit("FATAL: no layer-3 results produced")

    res = pd.DataFrame([strip_private(row) for row in results])
    res["stratum_order"] = res["stratum"].map({name: i for i, name in enumerate(STRATA_ORDER)})
    res["scope_order"] = res["feature_scope"].map({"physics_only": 0, "physics_plus_aimnet2": 1})
    res = res.sort_values(["target", "stratum_order", "scope_order"]).reset_index(drop=True)
    contrib = add_contribution_shares(pd.DataFrame(contributions))
    contrib["stratum_order"] = contrib["stratum"].map({name: i for i, name in enumerate(STRATA_ORDER)})
    contrib = contrib.sort_values(["target", "stratum_order", "feature_scope", "block"]).reset_index(drop=True)

    run = {
        "substrate_dir": str(substrate_dir),
        "out_dir": str(out_dir),
        "report_md": str(report_md),
        "seed": int(args.seed),
        "ridge_alpha": float(args.ridge_alpha),
        "within_split": args.within_split,
        "test_frame_fraction": float(args.test_frame_fraction),
        "embedding_transform": embedding_transform,
        "feature_blocks": {
            "physics": PHYSICS_BLOCK_ORDER,
            "aimnet2_ceiling": AIMNET_BLOCK_ORDER,
        },
        "mcconnell_source_mode": "valid_source_per_category",
        "change_of_basis_orthogonality_max": orth,
        "result_rows": int(len(res)),
        "contribution_rows": int(len(contrib)),
    }

    result_path = out_dir / "ensemble_layer3_results.csv"
    contrib_path = out_dir / "ensemble_layer3_contributions.csv"
    run_path = out_dir / "ensemble_layer3_run.json"
    audit_path = out_dir / "self_audit.json"
    labels_path = out_dir / "feature_labels.json"
    res.drop(columns=["stratum_order", "scope_order"]).to_csv(result_path, index=False)
    contrib.drop(columns=["stratum_order"]).to_csv(contrib_path, index=False)
    run_path.write_text(json.dumps(run, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    labels_path.write_text(json.dumps(feature_label_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(res, contrib, audit, run, report_md)

    print(f"wrote {result_path}")
    print(f"wrote {contrib_path}")
    print(f"wrote {run_path}")
    print(f"wrote {audit_path}")
    print(f"wrote {labels_path}")
    print(f"wrote {report_md}")
    print("\ntarget,stratum,scope,between_R2,within_R2")
    for _, row in res.iterrows():
        print(
            f"{row['target']},{row['stratum']},{row['feature_scope']},"
            f"{finite_fmt(row['between_LOAO_R2'])},{finite_fmt(row['within_frame_R2'])}"
        )


if __name__ == "__main__":
    main()
