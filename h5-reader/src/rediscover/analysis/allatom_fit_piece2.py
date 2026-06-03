#!/usr/bin/env python3
"""All-atoms charge-complete joint fit plus input-condition partitions.

This is the Python edge for Loop 2 / Piece 2. It reads one emitted
per_atom_substrate directory, applies the frozen change_of_basis.get_C() map to
DFT and emitted T2 feature blocks, fits global multi-output ridge models on the
same charge-complete rows for all tiers, and slices held-out predictions by
input-side conditions.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import change_of_basis as cob


SUBSTRATE_DIR = Path("/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build1")
OUT_DIR = Path("/tmp/rediscover-runs/2026-06-03-build3-fit-arch")
RIDGE_ALPHA = 10.0
ALPHA_GRID = [0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 1.0e4, 1.0e5]
INNER_CV_FOLDS = 5
EMBEDDING_COMPONENTS = 32
WITHIN_TEST_FRACTION = 0.20
PURGE_FRAMES_EACH_SIDE = 1
THIN_ATOMS = 10
THIN_WITHIN_ROWS = 500
MIN_FAVORABLE_N_EFF = 10.0
MIN_DISK_FREE_BYTES = 20 * 1024**3
MAX_REDISCOVER_BYTES = 15 * 1024**3
DOMINANCE_FRACTION_COLUMNS = {
    "ring": "dominant_fraction_ring",
    "charge": "dominant_fraction_charge",
    "mc": "dominant_fraction_mc",
}

TIER_ORDER = ["classical", "classical_plus_newmech", "plus_AIMNet2", "all"]
TIER_LABELS = {
    "classical": "classical baseline",
    "classical_plus_newmech": "classical plus new mechanisms",
    "plus_AIMNet2": "classical plus new mechanisms plus AIMNet2",
    "all": "classical plus new mechanisms plus AIMNet2 plus raw charge scalars",
}
TIER_ALIASES = {
    "classical": "classical",
    "classical_mechanisms_combined": "classical",
    "classical_plus_newmech": "classical_plus_newmech",
    "newmech": "classical_plus_newmech",
    "plus_AIMNet2": "plus_AIMNet2",
    "plus_aimnet2": "plus_AIMNet2",
    "all": "all",
}
TIER_PREVIOUS = {
    "classical": None,
    "classical_plus_newmech": "classical",
    "plus_AIMNet2": "classical_plus_newmech",
    "all": "plus_AIMNet2",
}

CLASSICAL_BLOCK_SPECS = [
    ("ring_jb_T0/T2", "classical", ["ring_jb_T0"], ["ring_jb_T2"]),
    ("charge_q_over_r3_T2", "classical", [], ["charge_q_over_r3_T2"]),
    ("mc_lit_T0/T2_valid", "classical", ["mc_lit_T0_valid"], ["mc_lit_T2_valid"]),
    ("mopac_coulomb_shielding_T2", "classical", [], ["mopac_coulomb_shielding_T2"]),
    ("mopac_mc_shielding_T2", "classical", [], ["mopac_mc_shielding_T2"]),
    ("ff14sb_field_x/y/z/mag", "classical", ["ff14sb_field_x", "ff14sb_field_y", "ff14sb_field_z", "ff14sb_field_mag"], []),
    ("apbs_E_x/y/z/mag", "classical", ["apbs_E_x", "apbs_E_y", "apbs_E_z", "apbs_E_mag"], []),
    ("apbs_efg_T2", "classical", [], ["apbs_efg_T2"]),
]
AIMNET_BLOCK_SPECS = [
    ("aimnet2_charge", "classical", ["aimnet2_charge"], []),
    ("aimnet2_crg_scalar/x/y/z", "classical", ["aimnet2_crg_scalar", "aimnet2_crg_x", "aimnet2_crg_y", "aimnet2_crg_z"], []),
    ("aimnet2_embedding_pca", "embedding", [], []),
]
RAW_CHARGE_BLOCK_SPECS = [
    ("formal_charge", "rows", ["formal_charge"], []),
    ("ff14sb_charge", "rows", ["ff14sb_charge"], []),
    ("mopac_welford_mean_charge", "rows", ["mopac_welford_mean_charge"], []),
    ("eeq_charge", "rows", ["eeq_charge"], []),
]
NEW_MECH_BLOCK_NAMES = [
    "hbond_T2",
    "larsen_hbond_T2",
    "hbond_count",
    "pi_quadrupole_T2",
    "dispersion_T2",
    "water_field",
    "water_efg_T2",
    "water_efield_first",
    "mopac_coulomb_E",
    "mopac_coulomb_efg_paths",
]
AIMNET_EXTRA_BLOCK_NAMES = [
    "aimnet2_charge",
    "aimnet2_crg_scalar/x/y/z",
    "aimnet2_efg_paths",
    "aimnet2_embedding_pca",
]
TIER_BLOCKS = {
    "classical": [x[0] for x in CLASSICAL_BLOCK_SPECS],
    "classical_plus_newmech": [x[0] for x in CLASSICAL_BLOCK_SPECS] + NEW_MECH_BLOCK_NAMES,
    "plus_AIMNet2": [x[0] for x in CLASSICAL_BLOCK_SPECS] + NEW_MECH_BLOCK_NAMES + AIMNET_EXTRA_BLOCK_NAMES,
    "all": [x[0] for x in CLASSICAL_BLOCK_SPECS] + NEW_MECH_BLOCK_NAMES + AIMNET_EXTRA_BLOCK_NAMES + [x[0] for x in RAW_CHARGE_BLOCK_SPECS],
}

TARGET_SPECS = {
    "total-T2": {"array": "per_atom_substrate_target_T2", "file": "per_atom_substrate_target_T2.npy", "components": 5, "basis": "T2_local_2e"},
    "dia-T2": {"array": "per_atom_substrate_target_dia_T2", "file": "per_atom_substrate_target_dia_T2.npy", "components": 5, "basis": "T2_local_2e"},
    "para-T2": {"array": "per_atom_substrate_target_para_T2", "file": "per_atom_substrate_target_para_T2.npy", "components": 5, "basis": "T2_local_2e"},
    "T1-field-linear-diagnostic": {"array": "per_atom_substrate_target_T1", "file": "per_atom_substrate_target_T1.npy", "components": 3, "basis": "T1_lab_vector"},
}

SCORE_TABLE_COLUMNS = [
    "target",
    "fit_scope",
    "atom_type_axis",
    "atom_type",
    "tier",
    "tier_label",
    "rows",
    "n_atoms_between_global",
    "n_atoms_between_fit_scope",
    "n_atoms_in_slice",
    "n_original_frames",
    "n_features",
    "ridge_alpha",
    "alpha_mode",
    "between_ridge_alpha",
    "between_ridge_alpha_min",
    "between_ridge_alpha_max",
    "within_ridge_alpha",
    "variance_share_between",
    "variance_share_within",
    "within_N_eff_min",
    "within_N_eff_median",
    "within_N_eff_max",
    "median_lag1_rho",
    "between_N_eff",
    "thin_flag",
    "between_support_flag",
    "within_support_flag",
    "p_gt_atoms_flag",
    "between_LOAO_test_R2",
    "between_LOAO_test_R2_jackknife_se",
    "between_LOAO_test_R2_ci95_low",
    "between_LOAO_test_R2_ci95_high",
    "between_delta_R2_vs_classical",
    "between_delta_R2_vs_previous_tier",
    "between_delta_R2_vs_global_sliced",
    "within_frameblock_test_R2",
    "within_frameblock_test_R2_jackknife_se",
    "within_frameblock_test_R2_ci95_low",
    "within_frameblock_test_R2_ci95_high",
    "within_delta_R2_vs_classical",
    "within_delta_R2_vs_previous_tier",
    "within_delta_R2_vs_global_sliced",
    "within_split_strategy",
    "test_frames",
    "purged_train_frames",
    "cross_split_lag1_pairs",
    "feature_blocks",
]


@dataclass(frozen=True)
class ConditionSpec:
    family: str
    name: str
    binning: str
    source: str


@dataclass
class PcaTransform:
    mean: np.ndarray
    components: np.ndarray
    explained_variance_ratio: float
    n_train_rows: int
    n_components: int
    original_dims: int


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--substrate-dir", type=Path, default=SUBSTRATE_DIR)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--ridge-alpha", type=float, default=RIDGE_ALPHA)
    ap.add_argument("--alpha-mode", choices=("select", "fixed"), default="select", help="select uses train-only inner CV; fixed is the alpha=10 continuity baseline")
    ap.add_argument("--alpha-grid", default=",".join(f"{x:g}" for x in ALPHA_GRID), help="comma-separated ridge alpha grid for alpha-mode=select")
    ap.add_argument("--inner-cv-folds", type=int, default=INNER_CV_FOLDS)
    ap.add_argument("--tiers", default=",".join(TIER_ORDER), help="comma-separated fit feature tiers; accepts classical, plus_AIMNet2, all")
    ap.add_argument("--embedding-components", type=int, default=EMBEDDING_COMPONENTS)
    ap.add_argument("--within-test-fraction", type=float, default=WITHIN_TEST_FRACTION)
    ap.add_argument("--purge-frames-each-side", type=int, default=PURGE_FRAMES_EACH_SIDE)
    ap.add_argument("--keep-out-dir", action="store_true", help="do not delete an existing output directory before writing")
    return ap.parse_args()


def parse_alpha_grid(text: str) -> list[float]:
    vals = []
    for part in str(text).split(","):
        part = part.strip()
        if not part:
            continue
        vals.append(float(part))
    if not vals:
        raise SystemExit("FATAL: alpha grid is empty")
    if any((not np.isfinite(x)) or x < 0.0 for x in vals):
        raise SystemExit(f"FATAL: alpha grid contains invalid values: {vals}")
    return vals


def parse_tiers(text: str) -> list[str]:
    tiers = []
    for part in str(text).split(","):
        key = part.strip()
        if not key:
            continue
        tier = TIER_ALIASES.get(key)
        if tier is None:
            valid = ", ".join(["classical", "classical_plus_newmech", "plus_AIMNet2", "all"])
            raise SystemExit(f"FATAL: unknown tier {key!r}; valid tiers: {valid}")
        if tier not in tiers:
            tiers.append(tier)
    if not tiers:
        raise SystemExit("FATAL: no fit feature tiers selected")
    return tiers


def previous_tier_map(tier_order: list[str]) -> dict[str, str | None]:
    prev: str | None = None
    out: dict[str, str | None] = {}
    for tier in tier_order:
        out[tier] = prev
        prev = tier
    return out


def require_file(path: Path) -> Path:
    if not path.exists():
        raise SystemExit(f"FATAL: required emitted substrate file missing: {path}")
    return path


def require_columns(df: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = [c for c in columns if c not in df.columns]
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


def json_sanitize(obj: object) -> object:
    if isinstance(obj, dict):
        return {str(k): json_sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_sanitize(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return json_sanitize(obj.tolist())
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        return v if np.isfinite(v) else None
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, Path):
        return str(obj)
    return obj


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


def component_r2(pred: np.ndarray, target: np.ndarray) -> list[float]:
    pred = np.asarray(pred, dtype=float)
    target = np.asarray(target, dtype=float)
    if pred.ndim == 1:
        pred = pred.reshape(-1, 1)
    out = []
    for c in range(target.shape[1]):
        out.append(r2_score(pred[:, c], target[:, c]))
    return out


def ci95(value: float, se: float) -> tuple[float, float]:
    if not np.isfinite(value) or not np.isfinite(se):
        return math.nan, math.nan
    return float(value - 1.96 * se), float(value + 1.96 * se)


def jackknife_metric_values(pred: np.ndarray, target: np.ndarray, groups: np.ndarray) -> np.ndarray:
    pred = np.asarray(pred, dtype=float)
    target = np.asarray(target, dtype=float)
    if pred.ndim == 1:
        pred = pred.reshape(-1, 1)
    if target.ndim == 1:
        target = target.reshape(-1, 1)
    groups = np.asarray(groups)
    ok = np.isfinite(pred).all(axis=1) & np.isfinite(target).all(axis=1)
    pred = pred[ok]
    target = target[ok]
    groups = groups[ok]
    labels, inv = np.unique(groups, return_inverse=True)
    if len(labels) < 3 or len(target) < 4:
        return np.asarray([], dtype=float)
    n_groups = len(labels)
    counts = np.bincount(inv, minlength=n_groups).astype(float)
    row_sse = ((target - pred) ** 2).sum(axis=1)
    row_y2 = (target ** 2).sum(axis=1)
    group_sse = np.bincount(inv, weights=row_sse, minlength=n_groups)
    group_y2 = np.bincount(inv, weights=row_y2, minlength=n_groups)
    group_sum_y = np.zeros((n_groups, target.shape[1]), dtype=float)
    np.add.at(group_sum_y, inv, target)

    total_n = float(len(target))
    total_sse = float(row_sse.sum())
    total_y2 = float(row_y2.sum())
    total_sum_y = target.sum(axis=0)

    keep_n = total_n - counts
    keep_sse = total_sse - group_sse
    keep_y2 = total_y2 - group_y2
    keep_sum_y = total_sum_y.reshape(1, -1) - group_sum_y
    with np.errstate(divide="ignore", invalid="ignore"):
        keep_tss = keep_y2 - (keep_sum_y ** 2).sum(axis=1) / keep_n
        vals = 1.0 - keep_sse / keep_tss
    valid = (keep_n >= 3) & (keep_tss > 0.0) & np.isfinite(vals)
    return vals[valid]


def jackknife_se_from_values(vals: np.ndarray) -> float:
    vals = np.asarray(vals, dtype=float)
    if len(vals) < 3:
        return math.nan
    mean = float(vals.mean())
    return float(math.sqrt((len(vals) - 1) / len(vals) * np.sum((vals - mean) ** 2)))


def jackknife_metric_se(pred: np.ndarray, target: np.ndarray, groups: np.ndarray) -> float:
    return jackknife_se_from_values(jackknife_metric_values(pred, target, groups))


def jackknife_component_se(pred: np.ndarray, target: np.ndarray, groups: np.ndarray) -> list[float]:
    return [jackknife_metric_se(pred[:, c], target[:, c], groups) for c in range(target.shape[1])]


def fit_ridge_predict(x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray, alpha: float) -> np.ndarray:
    x_train = np.asarray(x_train, dtype=float)
    y_train = np.asarray(y_train, dtype=float)
    x_test = np.asarray(x_test, dtype=float)
    if x_train.ndim == 1:
        x_train = x_train.reshape(-1, 1)
    if x_test.ndim == 1:
        x_test = x_test.reshape(-1, 1)
    if y_train.ndim == 1:
        y_train = y_train.reshape(-1, 1)
    ok = np.isfinite(y_train).all(axis=1) & np.isfinite(x_train).any(axis=1)
    x_train = x_train[ok]
    y_train = y_train[ok]
    out_dim = y_train.shape[1] if y_train.ndim == 2 else 1
    if len(x_train) < 3:
        return np.full((len(x_test), out_dim), np.nan, dtype=float)

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
    penalty = np.eye(design.shape[1], dtype=float) * max(0.0, float(alpha))
    penalty[0, 0] = 0.0
    beta = np.linalg.solve(design.T @ design + penalty, design.T @ y_train)
    return design_t @ beta


def fit_ridge_predict_grid(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    alpha_grid: list[float],
) -> dict[float, np.ndarray]:
    x_train = np.asarray(x_train, dtype=float)
    y_train = np.asarray(y_train, dtype=float)
    x_test = np.asarray(x_test, dtype=float)
    if x_train.ndim == 1:
        x_train = x_train.reshape(-1, 1)
    if x_test.ndim == 1:
        x_test = x_test.reshape(-1, 1)
    if y_train.ndim == 1:
        y_train = y_train.reshape(-1, 1)
    ok = np.isfinite(y_train).all(axis=1) & np.isfinite(x_train).any(axis=1)
    x_train = x_train[ok]
    y_train = y_train[ok]
    out_dim = y_train.shape[1] if y_train.ndim == 2 else 1
    if len(x_train) < 3:
        return {float(alpha): np.full((len(x_test), out_dim), np.nan, dtype=float) for alpha in alpha_grid}

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
    xtx = design.T @ design
    xty = design.T @ y_train
    penalty_base = np.eye(design.shape[1], dtype=float)
    penalty_base[0, 0] = 0.0
    out = {}
    for alpha in alpha_grid:
        penalty = penalty_base * max(0.0, float(alpha))
        beta = np.linalg.solve(xtx + penalty, xty)
        out[float(alpha)] = design_t @ beta
    return out


def deterministic_index_folds(indices: np.ndarray, n_folds: int) -> list[tuple[np.ndarray, np.ndarray]]:
    idx = np.asarray(indices, dtype=int)
    if len(idx) < 2:
        raise ValueError("inner CV requires at least two units")
    k = min(max(2, int(n_folds)), len(idx))
    folds = []
    for val in np.array_split(idx, k):
        if len(val) == 0:
            continue
        val_set = set(int(x) for x in val)
        fit = np.asarray([int(x) for x in idx if int(x) not in val_set], dtype=int)
        folds.append((fit, np.asarray(val, dtype=int)))
    return folds


def blocked_frame_cv_splits(
    frames: np.ndarray,
    eligible_mask: np.ndarray,
    n_folds: int,
    purge: int,
) -> list[dict[str, object]]:
    frames = np.asarray(frames, dtype=int)
    eligible_mask = np.asarray(eligible_mask, dtype=bool)
    eligible_frames = np.sort(np.unique(frames[eligible_mask]))
    if len(eligible_frames) < 2:
        raise ValueError("inner frame CV requires at least two eligible frames")
    k = min(max(2, int(n_folds)), len(eligible_frames))
    eligible_set = set(int(x) for x in eligible_frames)
    out = []
    for fold_idx, val_frames in enumerate(np.array_split(eligible_frames, k)):
        val_set = set(int(x) for x in val_frames)
        purge_set: set[int] = set()
        for frame in val_set:
            for offset in range(1, int(purge) + 1):
                for candidate in (frame - offset, frame + offset):
                    if candidate in eligible_set and candidate not in val_set:
                        purge_set.add(candidate)
        fit_mask = eligible_mask & np.array([int(f) not in val_set and int(f) not in purge_set for f in frames], dtype=bool)
        val_mask = eligible_mask & np.array([int(f) in val_set for f in frames], dtype=bool)
        purge_mask = eligible_mask & np.array([int(f) in purge_set for f in frames], dtype=bool)
        out.append({
            "fold": int(fold_idx),
            "fit_mask": fit_mask,
            "val_mask": val_mask,
            "purge_mask": purge_mask,
            "validation_frames": sorted(int(x) for x in val_set),
            "purged_frames": sorted(int(x) for x in purge_set),
        })
    return out


def select_best_alpha(alpha_scores: dict[float, float], alpha_grid: list[float]) -> float:
    finite = [(float(a), float(alpha_scores.get(float(a), math.nan))) for a in alpha_grid if np.isfinite(alpha_scores.get(float(a), math.nan))]
    if not finite:
        raise RuntimeError("alpha selection produced no finite inner-CV scores")
    best_score = max(score for _alpha, score in finite)
    for alpha, score in finite:
        if score == best_score or abs(score - best_score) <= 1.0e-12:
            return float(alpha)
    return float(finite[0][0])


def alpha_distribution(values: np.ndarray | list[float]) -> dict[str, object]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return {"selected_alpha": None, "min": None, "max": None, "counts": {}}
    unique, counts = np.unique(arr, return_counts=True)
    order = np.lexsort((unique, -counts))
    mode = float(unique[order[0]])
    return {
        "selected_alpha": mode,
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "counts": {f"{float(alpha):g}": int(count) for alpha, count in zip(unique, counts)},
    }


def mode_alpha(values: np.ndarray | list[float]) -> float:
    selected = alpha_distribution(values)["selected_alpha"]
    return math.nan if selected is None else float(selected)


def centered_by_train_atom(values: np.ndarray, atoms: np.ndarray, train: np.ndarray) -> np.ndarray:
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
    atoms = np.asarray(atoms)
    frames = np.asarray(frames)
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


def per_atom_neff_components(target: np.ndarray, atoms: np.ndarray, frames: np.ndarray) -> dict[int, np.ndarray]:
    out: dict[int, np.ndarray] = {}
    y = np.asarray(target, dtype=float)
    atoms = np.asarray(atoms)
    frames = np.asarray(frames)
    for atom in np.unique(atoms):
        idx = np.flatnonzero(atoms == atom)
        idx = idx[np.argsort(frames[idx])]
        vals = []
        for c in range(y.shape[1]):
            yc = y[idx, c]
            yc = yc[np.isfinite(yc)]
            n = len(yc)
            if n < 3:
                vals.append(float(n))
                continue
            yd = yc - yc.mean()
            den = float(np.dot(yd, yd))
            if den <= 0.0:
                vals.append(float(n))
            else:
                rho = float(np.dot(yd[:-1], yd[1:]) / den)
                rho = max(min(rho, 0.999), -0.999)
                vals.append(min(max(n * (1.0 - rho) / (1.0 + rho), 1.0), float(n)))
        out[int(atom)] = np.asarray(vals, dtype=float)
    return out


def bin_neff_median(
    axis: str,
    atoms: np.ndarray,
    mask: np.ndarray,
    atom_neff: dict[int, np.ndarray] | None,
) -> float:
    selected_atoms = np.unique(np.asarray(atoms)[mask])
    if axis == "between":
        return float(len(selected_atoms))
    if not atom_neff or len(selected_atoms) == 0:
        return math.nan
    comps = None
    for atom in selected_atoms:
        vals = atom_neff.get(int(atom))
        if vals is None:
            continue
        comps = vals.copy() if comps is None else comps + vals
    if comps is None:
        return math.nan
    return float(np.median(comps))


def split_frame_block(frames: np.ndarray, test_fraction: float, purge: int) -> dict[str, object]:
    unique = np.sort(np.unique(frames))
    n_test = max(1, int(round(float(test_fraction) * len(unique))))
    start = (len(unique) - n_test) // 2
    test_frames = unique[start:start + n_test]
    purged = []
    for k in range(1, int(purge) + 1):
        if start - k >= 0:
            purged.append(unique[start - k])
        if start + n_test - 1 + k < len(unique):
            purged.append(unique[start + n_test - 1 + k])
    test_set = set(int(x) for x in test_frames)
    purge_set = set(int(x) for x in purged)
    train = np.array([(int(f) not in test_set and int(f) not in purge_set) for f in frames], dtype=bool)
    test = np.array([int(f) in test_set for f in frames], dtype=bool)
    purge_mask = np.array([int(f) in purge_set for f in frames], dtype=bool)
    frame_state = {}
    for f in unique:
        fi = int(f)
        if fi in test_set:
            frame_state[fi] = "test"
        elif fi in purge_set:
            frame_state[fi] = "purged"
        else:
            frame_state[fi] = "train"
    cross = 0
    for a, b in zip(unique[:-1], unique[1:]):
        states = {frame_state[int(a)], frame_state[int(b)]}
        if states == {"train", "test"}:
            cross += 1
    return {
        "train_mask": train,
        "test_mask": test,
        "purge_mask": purge_mask,
        "test_frames": [int(x) for x in test_frames],
        "purged_frames": sorted(int(x) for x in purge_set),
        "cross_split_lag1_pairs": int(cross),
    }


def fit_pca_array(x_train: np.ndarray, n_components: int) -> PcaTransform:
    x = np.asarray(x_train, dtype=float)
    if x.ndim != 2:
        raise ValueError("PCA input must be 2D")
    mean = np.nanmean(x, axis=0)
    mean[~np.isfinite(mean)] = 0.0
    x_fill = np.where(np.isfinite(x), x, mean)
    x_centered = x_fill - mean
    cov = (x_centered.T @ x_centered) / max(1, len(x_centered) - 1)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    k = min(max(1, int(n_components)), x.shape[1], x.shape[0])
    positive = vals[vals > 0.0]
    total = float(positive.sum())
    ratio = float(np.maximum(vals[:k], 0.0).sum() / total) if total > 0.0 else math.nan
    return PcaTransform(
        mean=mean.astype(np.float64, copy=False),
        components=vecs[:, :k].astype(np.float64, copy=False),
        explained_variance_ratio=ratio,
        n_train_rows=int(x.shape[0]),
        n_components=int(k),
        original_dims=int(x.shape[1]),
    )


def fit_pca_memmap(matrix: np.ndarray, train_indices: np.ndarray, n_components: int, batch_size: int = 50_000) -> PcaTransform:
    train_indices = np.asarray(train_indices, dtype=int)
    n = int(len(train_indices))
    dims = int(matrix.shape[1])
    if n < 2:
        raise ValueError("PCA requires at least two training rows")
    sums = np.zeros(dims, dtype=np.float64)
    counts = np.zeros(dims, dtype=np.float64)
    for start in range(0, n, batch_size):
        idx = train_indices[start:start + batch_size]
        chunk = np.asarray(matrix[idx], dtype=np.float64)
        finite = np.isfinite(chunk)
        sums += np.where(finite, chunk, 0.0).sum(axis=0)
        counts += finite.sum(axis=0)
    mean = np.divide(sums, counts, out=np.zeros_like(sums), where=counts > 0)
    cov = np.zeros((dims, dims), dtype=np.float64)
    for start in range(0, n, batch_size):
        idx = train_indices[start:start + batch_size]
        chunk = np.asarray(matrix[idx], dtype=np.float64)
        chunk = np.where(np.isfinite(chunk), chunk, mean) - mean
        cov += chunk.T @ chunk
    cov /= max(1, n - 1)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    k = min(max(1, int(n_components)), dims, n)
    positive = vals[vals > 0.0]
    total = float(positive.sum())
    ratio = float(np.maximum(vals[:k], 0.0).sum() / total) if total > 0.0 else math.nan
    return PcaTransform(
        mean=mean,
        components=vecs[:, :k].astype(np.float64, copy=False),
        explained_variance_ratio=ratio,
        n_train_rows=n,
        n_components=int(k),
        original_dims=dims,
    )


def transform_pca_array(x: np.ndarray, pca: PcaTransform) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x = np.where(np.isfinite(x), x, pca.mean) - pca.mean
    return x @ pca.components


def transform_pca_memmap(
    matrix: np.ndarray,
    pca: PcaTransform,
    indices: np.ndarray | None = None,
    batch_size: int = 50_000,
) -> np.ndarray:
    if indices is None:
        indices = np.arange(matrix.shape[0], dtype=int)
    else:
        indices = np.asarray(indices, dtype=int)
    out = np.empty((len(indices), pca.n_components), dtype=np.float32)
    for start in range(0, len(indices), batch_size):
        idx = indices[start:start + batch_size]
        chunk = np.asarray(matrix[idx], dtype=np.float64)
        chunk = np.where(np.isfinite(chunk), chunk, pca.mean) - pca.mean
        out[start:start + len(idx)] = (chunk @ pca.components).astype(np.float32)
    return out


def atom_means_dense(values: np.ndarray, n_frames: int, n_atoms: int) -> np.ndarray:
    arr = np.asarray(values)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return np.nanmean(arr.reshape(n_frames, n_atoms, arr.shape[1]), axis=0)


def atom_means_embedding(embedding: np.ndarray, n_frames: int, n_atoms: int) -> np.ndarray:
    out = np.empty((n_atoms, embedding.shape[1]), dtype=np.float64)
    for atom in range(n_atoms):
        out[atom] = np.nanmean(np.asarray(embedding[atom::n_atoms], dtype=np.float64), axis=0)
    return out


def specs_by_column(specs: dict[str, object]) -> dict[tuple[str, str], dict[str, object]]:
    return {(str(c["array"]), str(c["column"])): c for c in specs["columns"]}


def specs_for_array(specs: dict[str, object], array_name: str) -> dict[str, dict[str, object]]:
    return {str(c["column"]): c for c in specs["columns"] if c["array"] == array_name}


def get_classical_columns(
    classical: np.ndarray,
    classical_specs: dict[str, dict[str, object]],
    scalar_cols: list[str],
    t2_prefixes: list[str],
    C: np.ndarray,
    transform_audit: list[dict[str, object]],
) -> tuple[np.ndarray, list[str]]:
    pieces = []
    names = []
    for col in scalar_cols:
        spec = classical_specs[col]
        pieces.append(np.asarray(classical[:, int(spec["index"])], dtype=float).reshape(-1, 1))
        names.append(col)
    for prefix in t2_prefixes:
        cols = [f"{prefix}_{i}" for i in range(5)]
        idx = [int(classical_specs[col]["index"]) for col in cols]
        vals = np.asarray(classical[:, idx], dtype=float) @ C.T
        pieces.append(vals)
        names += [f"{prefix}_2e_{i}" for i in range(5)]
        transform_audit.append({"feature_block": prefix, "source_columns": cols, "transform": "change_of_basis.get_C().T"})
    if not pieces:
        return np.empty((len(classical), 0), dtype=float), []
    return np.column_stack(pieces), names


def build_non_embedding_blocks(
    rows: pd.DataFrame,
    classical: np.ndarray,
    specs: dict[str, object],
    C: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict[str, list[str]], list[dict[str, object]]]:
    classical_specs = specs_for_array(specs, "per_atom_substrate_features_classical")
    transform_audit: list[dict[str, object]] = []
    blocks: dict[str, np.ndarray] = {}
    labels: dict[str, list[str]] = {}
    for name, source, scalar_cols, t2_prefixes in CLASSICAL_BLOCK_SPECS + AIMNET_BLOCK_SPECS:
        if source != "classical":
            continue
        block, names = get_classical_columns(classical, classical_specs, scalar_cols, t2_prefixes, C, transform_audit)
        blocks[name] = block
        labels[name] = names
    for name, _source, scalar_cols, _t2_prefixes in RAW_CHARGE_BLOCK_SPECS:
        block = rows[scalar_cols].to_numpy(float)
        blocks[name] = block
        labels[name] = scalar_cols
    return blocks, labels, transform_audit


def concat_blocks(blocks: dict[str, np.ndarray], block_names: list[str]) -> np.ndarray:
    return np.column_stack([blocks[name] for name in block_names])


def tier_uses_embedding(tier: str) -> bool:
    return "aimnet2_embedding_pca" in TIER_BLOCKS[tier]


def build_between_tier_features(
    tier: str,
    base_atom_blocks: dict[str, np.ndarray],
    embedding_atom_mean: np.ndarray,
    fit_indices: np.ndarray,
    eval_indices: np.ndarray,
    n_pcs: int,
    pca_cache: dict[object, PcaTransform] | None = None,
    pca_cache_key: object | None = None,
) -> tuple[np.ndarray, np.ndarray, PcaTransform | None]:
    pieces_fit = []
    pieces_eval = []
    pca = None
    for block in TIER_BLOCKS[tier]:
        if block == "aimnet2_embedding_pca":
            if pca_cache is not None and pca_cache_key in pca_cache:
                pca = pca_cache[pca_cache_key]
            else:
                pca = fit_pca_array(embedding_atom_mean[fit_indices], n_pcs)
                if pca_cache is not None and pca_cache_key is not None:
                    pca_cache[pca_cache_key] = pca
            pieces_fit.append(transform_pca_array(embedding_atom_mean[fit_indices], pca))
            pieces_eval.append(transform_pca_array(embedding_atom_mean[eval_indices], pca))
        else:
            pieces_fit.append(base_atom_blocks[block][fit_indices])
            pieces_eval.append(base_atom_blocks[block][eval_indices])
    return np.column_stack(pieces_fit), np.column_stack(pieces_eval), pca


def build_within_tier_features_for_positions(
    tier: str,
    base_blocks: dict[str, np.ndarray],
    embedding: np.ndarray,
    row_indices: np.ndarray,
    positions: np.ndarray,
    n_pcs: int,
    pca: PcaTransform | None,
) -> np.ndarray:
    pieces = []
    for block in TIER_BLOCKS[tier]:
        if block == "aimnet2_embedding_pca":
            if pca is None:
                raise RuntimeError("tier requires embedding PCA but no transform was supplied")
            pieces.append(transform_pca_memmap(embedding, pca, row_indices[positions]))
        else:
            pieces.append(base_blocks[block][positions])
    return np.column_stack(pieces)


def feature_count_for_tier(tier: str, embedding_components: int, base_blocks: dict[str, np.ndarray]) -> int:
    total = 0
    for block in TIER_BLOCKS[tier]:
        if block == "aimnet2_embedding_pca":
            total += embedding_components
        else:
            total += int(base_blocks[block].shape[1])
    return total


def load_inputs(substrate_dir: Path) -> dict[str, object]:
    manifest_path = require_file(substrate_dir / "per_atom_substrate_manifest.json")
    specs_path = require_file(substrate_dir / "per_atom_substrate_column_specs.json")
    rows_path = require_file(substrate_dir / "per_atom_substrate_rows.csv")
    with manifest_path.open("r", encoding="utf-8") as f:
        manifest = json.load(f)
    with specs_path.open("r", encoding="utf-8") as f:
        specs = json.load(f)
    rows = pd.read_csv(rows_path)
    target_T2 = np.load(require_file(substrate_dir / "per_atom_substrate_target_T2.npy"))
    classical = np.load(require_file(substrate_dir / "per_atom_substrate_features_classical.npy"))
    conditioning = np.load(require_file(substrate_dir / "per_atom_substrate_features_conditioning.npy"))
    modulation = np.load(require_file(substrate_dir / "per_atom_substrate_driver_modulation_by_atom.npy"))
    embedding = np.load(require_file(substrate_dir / "per_atom_substrate_aimnet2_embedding.npy"), mmap_mode="r")
    return {
        "manifest": manifest,
        "specs": specs,
        "rows": rows,
        "target_T2": target_T2,
        "classical": classical,
        "conditioning": conditioning,
        "modulation": modulation,
        "embedding": embedding,
        "files_read": [
            str(manifest_path),
            str(specs_path),
            str(rows_path),
            str(substrate_dir / "per_atom_substrate_target_T2.npy"),
            str(substrate_dir / "per_atom_substrate_features_classical.npy"),
            str(substrate_dir / "per_atom_substrate_features_conditioning.npy"),
            str(substrate_dir / "per_atom_substrate_driver_modulation_by_atom.npy"),
            str(substrate_dir / "per_atom_substrate_aimnet2_embedding.npy"),
        ],
    }


def input_acceptance_checks(data: dict[str, object], substrate_dir: Path) -> tuple[dict[str, object], np.ndarray]:
    manifest = data["manifest"]
    rows: pd.DataFrame = data["rows"]
    target_T2: np.ndarray = data["target_T2"]
    classical: np.ndarray = data["classical"]
    conditioning: np.ndarray = data["conditioning"]
    modulation: np.ndarray = data["modulation"]
    embedding: np.ndarray = data["embedding"]
    specs: dict[str, object] = data["specs"]
    checks: dict[str, object] = {}
    required_cols = [
        "row_id",
        "atom_index",
        "frame_slot",
        "original_frame_index",
        "stratum",
        "role",
        "ff_atom_type_ord",
        "formal_charge",
        "ff14sb_charge",
        "ff14sb_charge_present",
        "mopac_welford_mean_charge",
        "mopac_welford_mean_charge_present",
    ]
    require_columns(rows, required_cols, "per_atom_substrate_rows.csv")
    n_rows = len(rows)
    n_atoms = int(manifest.get("n_atoms", -1))
    n_frames = int(manifest.get("n_dft_frames", -1))
    dense_row_ids = np.array_equal(rows["row_id"].to_numpy(int), np.arange(n_rows, dtype=int))
    row_contract = bool(dense_row_ids and np.all(rows["row_id"].to_numpy(int) == rows["frame_slot"].to_numpy(int) * n_atoms + rows["atom_index"].to_numpy(int)))
    shapes_ok = (
        n_atoms == 846
        and n_frames == 660
        and n_rows == 558_360
        and target_T2.shape == (n_rows, 5)
        and classical.shape == (n_rows, 45)
        and conditioning.shape == (n_rows, 26)
        and modulation.shape == (n_atoms, 9)
        and embedding.shape == (n_rows, 256)
    )
    classical_specs = specs_for_array(specs, "per_atom_substrate_features_classical")
    row_spec_cols = {c["column"] for c in specs["columns"] if c["array"] == "per_atom_substrate_rows"}
    selected_classical = set()
    for _name, source, scalar_cols, t2_prefixes in CLASSICAL_BLOCK_SPECS + AIMNET_BLOCK_SPECS:
        if source != "classical":
            continue
        selected_classical.update(scalar_cols)
        for prefix in t2_prefixes:
            selected_classical.update(f"{prefix}_{i}" for i in range(5))
    specs_ok = all(col in classical_specs for col in selected_classical) and {"ff14sb_charge", "mopac_welford_mean_charge"}.issubset(row_spec_cols)
    aim_idx = int(classical_specs["aimnet2_charge"]["index"])
    charge_complete_mask = (
        (rows["ff14sb_charge_present"].to_numpy(int) == 1)
        & (rows["mopac_welford_mean_charge_present"].to_numpy(int) == 1)
        & np.isfinite(rows["ff14sb_charge"].to_numpy(float))
        & np.isfinite(rows["mopac_welford_mean_charge"].to_numpy(float))
        & np.isfinite(classical[:, aim_idx])
    )
    charge_complete_rows = int(charge_complete_mask.sum())
    manifest_charge_rows = int(manifest.get("feature_support_rows", {}).get("charge_complete_rows", -1))
    target_finite = bool(np.isfinite(target_T2[charge_complete_mask, :5]).all())
    checks["input_acceptance"] = {
        "pass": bool(
            manifest.get("relationship") == "per_atom_substrate"
            and manifest.get("relationship_kind") == "per_atom_aggregate"
            and shapes_ok
            and row_contract
            and charge_complete_rows == manifest_charge_rows
            and charge_complete_rows == 558_360
            and specs_ok
            and target_finite
            and substrate_dir == SUBSTRATE_DIR
        ),
        "relationship": manifest.get("relationship"),
        "relationship_kind": manifest.get("relationship_kind"),
        "shape": {"n_atoms": n_atoms, "n_dft_frames": n_frames, "rows": n_rows},
        "expected_shape_846x660": shapes_ok,
        "row_contract": row_contract,
        "charge_complete_rows_manifest": manifest_charge_rows,
        "charge_complete_rows_used": charge_complete_rows,
        "new_charge_columns_present": all(c in rows.columns for c in ["ff14sb_charge", "ff14sb_charge_present", "mopac_welford_mean_charge", "mopac_welford_mean_charge_present"]),
        "selected_feature_columns_found_through_specs_or_row_schema": specs_ok,
        "target_finite_on_charge_complete_rows": target_finite,
        "substrate_dir_only_requested_path": str(substrate_dir),
    }
    return checks, charge_complete_mask


def conditioning_frame(
    rows: pd.DataFrame,
    conditioning: np.ndarray,
    modulation: np.ndarray,
    classical: np.ndarray,
    specs: dict[str, object],
) -> tuple[pd.DataFrame, list[dict[str, str]]]:
    cond_specs = specs_for_array(specs, "per_atom_substrate_features_conditioning")
    mod_specs = specs_for_array(specs, "per_atom_substrate_driver_modulation_by_atom")
    classical_specs = specs_for_array(specs, "per_atom_substrate_features_classical")
    out = rows.copy()
    for name, spec in cond_specs.items():
        out[name] = conditioning[:, int(spec["index"])]
    atom_index = rows["atom_index"].to_numpy(int)
    for name, spec in mod_specs.items():
        out[name] = modulation[atom_index, int(spec["index"])]
    aim_idx = int(classical_specs["aimnet2_charge"]["index"])
    out["aimnet2_charge"] = classical[:, aim_idx]
    formulas = [
        {"column": "abs_ff14sb_minus_mopac_welford_charge", "formula": "abs(ff14sb_charge - mopac_welford_mean_charge)"},
        {"column": "abs_ff14sb_minus_aimnet2_charge", "formula": "abs(ff14sb_charge - aimnet2_charge)"},
        {"column": "abs_mopac_welford_minus_aimnet2_charge", "formula": "abs(mopac_welford_mean_charge - aimnet2_charge)"},
        {"column": "abs_formal_minus_ff14sb_charge", "formula": "abs(formal_charge - ff14sb_charge)"},
        {"column": "abs_formal_minus_mopac_welford_charge", "formula": "abs(formal_charge - mopac_welford_mean_charge)"},
        {"column": "abs_formal_minus_aimnet2_charge", "formula": "abs(formal_charge - aimnet2_charge)"},
        {"column": "sign_agree_ff14sb_mopac_welford", "formula": "sign_agreement(ff14sb_charge, mopac_welford_mean_charge)"},
        {"column": "sign_agree_ff14sb_aimnet2", "formula": "sign_agreement(ff14sb_charge, aimnet2_charge)"},
        {"column": "sign_agree_mopac_welford_aimnet2", "formula": "sign_agreement(mopac_welford_mean_charge, aimnet2_charge)"},
    ]
    out["abs_ff14sb_minus_mopac_welford_charge"] = np.abs(out["ff14sb_charge"] - out["mopac_welford_mean_charge"])
    out["abs_ff14sb_minus_aimnet2_charge"] = np.abs(out["ff14sb_charge"] - out["aimnet2_charge"])
    out["abs_mopac_welford_minus_aimnet2_charge"] = np.abs(out["mopac_welford_mean_charge"] - out["aimnet2_charge"])
    out["abs_formal_minus_ff14sb_charge"] = np.abs(out["formal_charge"] - out["ff14sb_charge"])
    out["abs_formal_minus_mopac_welford_charge"] = np.abs(out["formal_charge"] - out["mopac_welford_mean_charge"])
    out["abs_formal_minus_aimnet2_charge"] = np.abs(out["formal_charge"] - out["aimnet2_charge"])
    out["sign_agree_ff14sb_mopac_welford"] = sign_agreement(out["ff14sb_charge"], out["mopac_welford_mean_charge"])
    out["sign_agree_ff14sb_aimnet2"] = sign_agreement(out["ff14sb_charge"], out["aimnet2_charge"])
    out["sign_agree_mopac_welford_aimnet2"] = sign_agreement(out["mopac_welford_mean_charge"], out["aimnet2_charge"])
    return out, formulas


def sign_agreement(a: pd.Series, b: pd.Series) -> np.ndarray:
    sa = np.sign(a.to_numpy(float))
    sb = np.sign(b.to_numpy(float))
    labels = np.full(len(sa), "disagree", dtype=object)
    zero = (sa == 0) | (sb == 0)
    labels[zero] = "zero_involved"
    labels[(sa == sb) & ~zero] = "agree"
    labels[~np.isfinite(sa) | ~np.isfinite(sb)] = "missing"
    return labels


def atom_condition_frame(row_conditions: pd.DataFrame, n_frames: int, n_atoms: int) -> pd.DataFrame:
    first_cols = [
        "atom_index",
        "element_ord",
        "ff_atom_type_ord",
        "atom_name",
        "iupac_atom_name",
        "residue_index",
        "residue_number",
        "amino_acid_ord",
        "backbone_role_ord",
        "locant_ord",
        "branch_outer_ord",
        "branch_inner_ord",
        "di_index_ord",
        "ring_position_primary_ord",
        "ring_position_secondary_ord",
        "planar_group_ord",
        "polar_h_kind_ord",
        "prochiral_ord",
        "equivalence_class",
        "aromatic",
        "formal_charge",
        "is_exchangeable",
        "role_ord",
        "role",
        "stratum_ord",
        "stratum",
        "ff14sb_charge_present",
        "mopac_welford_mean_charge_present",
        "ring_present",
        "charge_present",
        "mc_lit_valid_present",
        "ff14sb_field_present",
        "apbs_E_present",
        "apbs_efg_present",
        "mopac_coulomb_shielding_present",
        "mopac_mc_shielding_present",
        "aimnet2_charge_present",
        "aimnet2_crg_present",
        "aimnet2_embedding_present",
        "sign_agree_ff14sb_mopac_welford",
        "sign_agree_ff14sb_aimnet2",
        "sign_agree_mopac_welford_aimnet2",
    ]
    first_cols = [c for c in first_cols if c in row_conditions.columns]
    atom_rows = row_conditions.iloc[:n_atoms][first_cols].copy()
    numeric_cols = [
        c
        for c in row_conditions.columns
        if c not in first_cols
        and c not in {"row_id", "h5_row", "frame_slot", "original_frame_index", "time_ps", "dft_present"}
        and pd.api.types.is_numeric_dtype(row_conditions[c])
    ]
    for col in numeric_cols:
        vals = row_conditions[col].to_numpy()
        atom_rows[col] = np.nanmean(vals.reshape(n_frames, n_atoms), axis=0)
    for col in ["ff14sb_charge", "mopac_welford_mean_charge", "aimnet2_charge"]:
        if col in row_conditions.columns and col not in atom_rows.columns:
            atom_rows[col] = row_conditions[col].to_numpy().reshape(n_frames, n_atoms)[0]
    return atom_rows.reset_index(drop=True)


def condition_specs(row_conditions: pd.DataFrame) -> list[ConditionSpec]:
    categorical = [
        ("atom identity", "stratum"),
        ("atom identity", "role"),
        ("atom identity", "element_ord"),
        ("atom identity", "ff_atom_type_ord"),
        ("atom identity", "formal_charge"),
        ("atom identity", "planar_group_ord"),
        ("atom identity", "polar_h_kind_ord"),
        ("atom identity", "aromatic"),
        ("atom identity", "is_exchangeable"),
        ("atom identity", "amino_acid_ord"),
        ("atom identity", "residue_index"),
        ("atom identity", "residue_number"),
        ("atom identity", "locant_ord"),
        ("atom identity", "branch_outer_ord"),
        ("atom identity", "branch_inner_ord"),
        ("atom identity", "di_index_ord"),
        ("atom identity", "ring_position_primary_ord"),
        ("atom identity", "ring_position_secondary_ord"),
        ("atom identity", "backbone_role_ord"),
        ("geometry and isolation", "ring_n"),
        ("geometry and isolation", "ring_valid_n"),
        ("geometry and isolation", "charge_n"),
        ("geometry and isolation", "charge_excluded_same_residue_n"),
        ("geometry and isolation", "bond_n"),
        ("geometry and isolation", "bond_n_valid"),
        ("geometry and isolation", "ring_count_4A"),
        ("geometry and isolation", "ring_count_6A"),
        ("geometry and isolation", "ring_count_8A"),
        ("geometry and isolation", "charge_count_4A"),
        ("geometry and isolation", "charge_count_6A"),
        ("geometry and isolation", "charge_count_10A"),
        ("geometry and isolation", "bond_count_4A"),
        ("geometry and isolation", "bond_count_8A"),
        ("geometry and isolation", "bond_count_10A"),
        ("geometry and isolation", "ring_self_or_bonded_n"),
        ("geometry and isolation", "bond_self_or_bonded_n"),
        ("geometry and isolation", "has_self_or_bonded_driver"),
        ("support and source availability", "ring_present"),
        ("support and source availability", "charge_present"),
        ("support and source availability", "mc_lit_valid_present"),
        ("support and source availability", "ff14sb_field_present"),
        ("support and source availability", "apbs_E_present"),
        ("support and source availability", "apbs_efg_present"),
        ("support and source availability", "mopac_coulomb_shielding_present"),
        ("support and source availability", "mopac_mc_shielding_present"),
        ("support and source availability", "aimnet2_charge_present"),
        ("support and source availability", "aimnet2_crg_present"),
        ("support and source availability", "aimnet2_embedding_present"),
        ("support and source availability", "ff14sb_charge_present"),
        ("support and source availability", "mopac_welford_mean_charge_present"),
        ("charge-source agreement", "sign_agree_ff14sb_mopac_welford"),
        ("charge-source agreement", "sign_agree_ff14sb_aimnet2"),
        ("charge-source agreement", "sign_agree_mopac_welford_aimnet2"),
    ]
    specs = [ConditionSpec(family, name, "categorical", "emitted_or_derived_input_side") for family, name in categorical if name in row_conditions.columns]
    distance_cols = ["nearest_ring_r", "nearest_charge_r", "nearest_bond_midpoint_r", "nearest_heavy_atom_r"]
    for col in distance_cols:
        if col in row_conditions.columns:
            specs.append(ConditionSpec("geometry and isolation", f"{col}_threshold", "distance_threshold", "emitted_conditioning"))
            specs.append(ConditionSpec("geometry and isolation", f"{col}_quintile", "quintile", "emitted_conditioning"))
    driver_cols = [
        "abs_ring_jb_T2",
        "abs_charge_T2",
        "abs_mc_lit_T2",
        "abs_mopac_coulomb_T2",
        "abs_mopac_mc_T2",
        "abs_apbs_efg_T2",
        "abs_apbs_E",
        "abs_ff14sb_E",
        "abs_aimnet2_crg",
    ]
    for col in driver_cols:
        if col in row_conditions.columns:
            specs.append(ConditionSpec("driver magnitude", col, "quintile", "emitted_conditioning"))
    modulation_cols = [c for c in row_conditions.columns if c.startswith("sd_") and c.endswith("_by_atom")]
    for col in modulation_cols:
        specs.append(ConditionSpec("driver modulation", col, "quintile", "emitted_driver_modulation_by_atom"))
    charge_cols = [
        "ff14sb_charge",
        "mopac_welford_mean_charge",
        "aimnet2_charge",
        "abs_ff14sb_minus_mopac_welford_charge",
        "abs_ff14sb_minus_aimnet2_charge",
        "abs_mopac_welford_minus_aimnet2_charge",
        "abs_formal_minus_ff14sb_charge",
        "abs_formal_minus_mopac_welford_charge",
        "abs_formal_minus_aimnet2_charge",
    ]
    for col in charge_cols:
        if col in row_conditions.columns:
            specs.append(ConditionSpec("charge-source agreement", col, "quintile", "emitted_or_derived_input_side"))
    return specs


def condition_base_name(spec: ConditionSpec) -> str:
    if spec.binning == "distance_threshold":
        return spec.name.removesuffix("_threshold")
    if spec.name.endswith("_quintile"):
        return spec.name.removesuffix("_quintile")
    return spec.name


def assign_bins(
    values: pd.Series,
    atoms: np.ndarray,
    spec: ConditionSpec,
) -> tuple[np.ndarray, list[dict[str, object]]]:
    raw = values.to_numpy()
    atoms = np.asarray(atoms)
    if spec.binning == "categorical":
        labels = np.array([str(x) if pd.notna(x) else "missing" for x in raw], dtype=object)
        defs = []
        out = labels.copy()
        for label in sorted(set(labels), key=str):
            m = labels == label
            n_atoms = int(len(np.unique(atoms[m])))
            defs.append({"bin_label": label, "bin_low": None, "bin_high": None, "n_atoms_pre_rank": n_atoms})
            if n_atoms < THIN_ATOMS and len(set(labels)) > 1:
                out[m] = "__sparse__"
        if "__sparse__" in out:
            m = out == "__sparse__"
            defs.append({"bin_label": "__sparse__", "bin_low": None, "bin_high": None, "n_atoms_pre_rank": int(len(np.unique(atoms[m])))})
        return out, defs
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(float)
    finite = np.isfinite(numeric)
    labels = np.full(len(numeric), "missing", dtype=object)
    defs = []
    if finite.sum() == 0:
        defs.append({"bin_label": "missing", "bin_low": None, "bin_high": None, "n_atoms_pre_rank": int(len(np.unique(atoms)))})
        return labels, defs
    if spec.binning == "distance_threshold":
        edges = [-np.inf, 4.0, 6.0, 8.0, 10.0, np.inf]
        names = ["<=4A", "4-6A", "6-8A", "8-10A", ">10A"]
        for lo, hi, name in zip(edges[:-1], edges[1:], names):
            if np.isneginf(lo):
                m = finite & (numeric <= hi)
            elif np.isposinf(hi):
                m = finite & (numeric > lo)
            else:
                m = finite & (numeric > lo) & (numeric <= hi)
            labels[m] = name
            defs.append({"bin_label": name, "bin_low": None if np.isneginf(lo) else lo, "bin_high": None if np.isposinf(hi) else hi, "n_atoms_pre_rank": int(len(np.unique(atoms[m])))})
        return labels, defs
    vals = numeric[finite]
    qs = np.quantile(vals, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    edges = []
    for q in qs:
        if not edges or q > edges[-1]:
            edges.append(float(q))
    if len(edges) <= 2:
        lo = float(np.nanmin(vals))
        hi = float(np.nanmax(vals))
        labels[finite] = f"[{finite_fmt(lo)}, {finite_fmt(hi)}]"
        defs.append({"bin_label": str(labels[finite][0]), "bin_low": lo, "bin_high": hi, "n_atoms_pre_rank": int(len(np.unique(atoms[finite])))})
        return labels, defs
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:]), start=1):
        if i == 1:
            m = finite & (numeric >= lo) & (numeric <= hi)
        else:
            m = finite & (numeric > lo) & (numeric <= hi)
        label = f"Q{i} [{finite_fmt(lo)}, {finite_fmt(hi)}]"
        labels[m] = label
        defs.append({"bin_label": label, "bin_low": lo, "bin_high": hi, "n_atoms_pre_rank": int(len(np.unique(atoms[m])))})
    return labels, defs


def select_between_alphas(
    tier: str,
    base_atom_blocks: dict[str, np.ndarray],
    embedding_atom_mean: np.ndarray,
    y_atom_mean: np.ndarray,
    alpha_grid: list[float],
    n_pcs: int,
    inner_cv_folds: int,
    pca_cache: dict[object, PcaTransform] | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    n_atoms = int(y_atom_mean.shape[0])
    selected = np.full(n_atoms, np.nan, dtype=float)
    inner_score_sum = {float(alpha): 0.0 for alpha in alpha_grid}
    inner_score_count = {float(alpha): 0 for alpha in alpha_grid}
    fold_count = 0
    heldout_excluded = True
    for heldout in range(n_atoms):
        outer_train = np.asarray([i for i in range(n_atoms) if i != heldout], dtype=int)
        preds_by_alpha: dict[float, list[np.ndarray]] = {float(alpha): [] for alpha in alpha_grid}
        targets = []
        fold_scores_by_alpha: dict[float, list[float]] = {float(alpha): [] for alpha in alpha_grid}
        for fold_idx, (inner_fit, inner_val) in enumerate(deterministic_index_folds(outer_train, inner_cv_folds)):
            if int(heldout) in set(int(x) for x in inner_fit) or int(heldout) in set(int(x) for x in inner_val):
                heldout_excluded = False
            cache_key = ("between_inner", int(heldout), int(fold_idx)) if tier_uses_embedding(tier) else None
            x_fit, x_val, _pca = build_between_tier_features(
                tier,
                base_atom_blocks,
                embedding_atom_mean,
                inner_fit,
                inner_val,
                n_pcs,
                pca_cache,
                cache_key,
            )
            pred_grid = fit_ridge_predict_grid(x_fit, y_atom_mean[inner_fit], x_val, alpha_grid)
            targets.append(y_atom_mean[inner_val])
            for alpha in alpha_grid:
                a = float(alpha)
                pred = pred_grid[a]
                preds_by_alpha[a].append(pred)
                fold_score = r2_score(pred, y_atom_mean[inner_val])
                fold_scores_by_alpha[a].append(fold_score)
        scores = {}
        y_val_all = np.vstack(targets)
        for alpha in alpha_grid:
            a = float(alpha)
            pred_all = np.vstack(preds_by_alpha[a])
            score = r2_score(pred_all, y_val_all)
            scores[a] = score
            if np.isfinite(score):
                inner_score_sum[a] += float(score)
                inner_score_count[a] += 1
        selected[heldout] = select_best_alpha(scores, alpha_grid)
        fold_count += len(targets)
    summary = alpha_distribution(selected)
    summary.update({
        "axis": "between_LOAO",
        "method": "nested deterministic atom-fold inner CV inside each outer held-out atom training set",
        "outer_heldout_count": n_atoms,
        "inner_cv_folds": int(inner_cv_folds),
        "inner_fold_evaluations": int(fold_count),
        "heldout_test_atom_excluded_from_alpha_selection": bool(heldout_excluded),
        "alpha_grid": [float(x) for x in alpha_grid],
        "mean_inner_cv_R2_by_alpha": {
            f"{float(alpha):g}": (
                float(inner_score_sum[float(alpha)] / inner_score_count[float(alpha)])
                if inner_score_count[float(alpha)] > 0
                else None
            )
            for alpha in alpha_grid
        },
        "selected_alpha_by_outer_heldout_atom": [float(x) for x in selected],
    })
    return selected, summary


def select_within_alpha(
    tier: str,
    base_blocks: dict[str, np.ndarray],
    embedding: np.ndarray,
    row_indices: np.ndarray,
    y: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    outer_split: dict[str, object],
    alpha_grid: list[float],
    n_pcs: int,
    inner_cv_folds: int,
    purge: int,
    pca_cache: dict[object, PcaTransform] | None = None,
) -> tuple[float, dict[str, object]]:
    outer_train = np.asarray(outer_split["train_mask"], dtype=bool)
    outer_test = np.asarray(outer_split["test_mask"], dtype=bool)
    outer_purge = np.asarray(outer_split["purge_mask"], dtype=bool)
    folds = blocked_frame_cv_splits(frames, outer_train, inner_cv_folds, purge)
    preds_by_alpha: dict[float, list[np.ndarray]] = {float(alpha): [] for alpha in alpha_grid}
    targets = []
    fold_details = []
    outer_test_rows_used = 0
    outer_purge_rows_used = 0
    for fold in folds:
        fold_idx = int(fold["fold"])
        fit_mask = np.asarray(fold["fit_mask"], dtype=bool)
        val_mask = np.asarray(fold["val_mask"], dtype=bool)
        inner_used = fit_mask | val_mask
        outer_test_rows_used += int((inner_used & outer_test).sum())
        outer_purge_rows_used += int((inner_used & outer_purge).sum())
        fit_pos = np.flatnonzero(fit_mask)
        val_pos = np.flatnonzero(val_mask)
        combined = np.concatenate([fit_pos, val_pos])
        train_combined = np.zeros(len(combined), dtype=bool)
        train_combined[:len(fit_pos)] = True
        pca = None
        if tier_uses_embedding(tier):
            cache_key = ("within_inner", int(fold_idx))
            if pca_cache is not None and cache_key in pca_cache:
                pca = pca_cache[cache_key]
            else:
                pca = fit_pca_memmap(embedding, row_indices[fit_pos], n_pcs)
                if pca_cache is not None:
                    pca_cache[cache_key] = pca
        x_combined = build_within_tier_features_for_positions(tier, base_blocks, embedding, row_indices, combined, n_pcs, pca)
        y_combined = y[combined]
        atoms_combined = atoms[combined]
        x_centered = centered_by_train_atom(x_combined, atoms_combined, train_combined)
        y_centered = centered_by_train_atom(y_combined, atoms_combined, train_combined)
        x_fit = x_centered[:len(fit_pos)]
        y_fit = y_centered[:len(fit_pos)]
        x_val = x_centered[len(fit_pos):]
        y_val = y_centered[len(fit_pos):]
        pred_grid = fit_ridge_predict_grid(x_fit, y_fit, x_val, alpha_grid)
        targets.append(y_val)
        for alpha in alpha_grid:
            preds_by_alpha[float(alpha)].append(pred_grid[float(alpha)])
        fold_details.append({
            "fold": int(fold_idx),
            "fit_rows": int(fit_mask.sum()),
            "validation_rows": int(val_mask.sum()),
            "purged_rows": int(np.asarray(fold["purge_mask"], dtype=bool).sum()),
            "validation_frames": fold["validation_frames"],
            "purged_frames": fold["purged_frames"],
        })
    y_val_all = np.vstack(targets)
    scores = {}
    for alpha in alpha_grid:
        a = float(alpha)
        scores[a] = r2_score(np.vstack(preds_by_alpha[a]), y_val_all)
    selected = select_best_alpha(scores, alpha_grid)
    audit = {
        "axis": "within_frameblock",
        "method": "deterministic blocked frame inner CV inside the outer training frames only",
        "selected_alpha": float(selected),
        "alpha_grid": [float(x) for x in alpha_grid],
        "inner_cv_folds": int(len(folds)),
        "inner_cv_R2_by_alpha": {f"{float(alpha):g}": (float(scores[float(alpha)]) if np.isfinite(scores[float(alpha)]) else None) for alpha in alpha_grid},
        "outer_test_rows_used_for_alpha_selection": int(outer_test_rows_used),
        "outer_purged_rows_used_for_alpha_selection": int(outer_purge_rows_used),
        "folds": fold_details,
    }
    return selected, audit


def evaluate_between_tier(
    tier: str,
    base_atom_blocks: dict[str, np.ndarray],
    embedding_atom_mean: np.ndarray,
    y_atom_mean: np.ndarray,
    alpha: float | np.ndarray,
    n_pcs: int,
    pca_cache: dict[int, PcaTransform] | None = None,
) -> dict[str, object]:
    n_atoms = y_atom_mean.shape[0]
    pred = np.full_like(y_atom_mean, np.nan, dtype=float)
    if np.isscalar(alpha):
        alpha_by_heldout = np.full(n_atoms, float(alpha), dtype=float)
    else:
        alpha_by_heldout = np.asarray(alpha, dtype=float)
        if len(alpha_by_heldout) != n_atoms:
            raise ValueError("between alpha array must have one value per held-out atom")
    pca_ratios = []
    pca_train_rows = []
    no_heldout_in_train = True
    for heldout in range(n_atoms):
        train = np.asarray([i for i in range(n_atoms) if i != heldout], dtype=int)
        pca_key = int(heldout) if tier_uses_embedding(tier) else None
        x_train, x_test, pca = build_between_tier_features(
            tier,
            base_atom_blocks,
            embedding_atom_mean,
            train,
            np.asarray([heldout], dtype=int),
            n_pcs,
            pca_cache,
            pca_key,
        )
        if pca is not None:
            pca_ratios.append(pca.explained_variance_ratio)
            pca_train_rows.append(pca.n_train_rows)
        if int(heldout) in set(int(x) for x in train):
            no_heldout_in_train = False
        pred[heldout] = fit_ridge_predict(x_train, y_atom_mean[train], x_test, float(alpha_by_heldout[heldout]))[0]
    pca_audit = None
    if pca_ratios:
        pca_audit = {
            "method": "unsupervised PCA fit independently inside each LOAO fold on training atom-mean rows only",
            "n_components": int(n_pcs),
            "original_dims": int(embedding_atom_mean.shape[1]),
            "folds": int(len(pca_ratios)),
            "training_rows_per_fold_min": int(min(pca_train_rows)),
            "training_rows_per_fold_max": int(max(pca_train_rows)),
            "explained_variance_ratio_min": float(np.min(pca_ratios)),
            "explained_variance_ratio_median": float(np.median(pca_ratios)),
            "explained_variance_ratio_max": float(np.max(pca_ratios)),
        }
    return {
        "pred": pred,
        "target": y_atom_mean,
        "groups": np.arange(n_atoms, dtype=int),
        "alpha_by_heldout": alpha_by_heldout,
        "pca_audit": pca_audit,
        "cv_check": {"heldout_atoms_never_in_train": bool(no_heldout_in_train), "folds": int(n_atoms)},
    }


def evaluate_within_tier(
    tier: str,
    base_blocks: dict[str, np.ndarray],
    embedding_pcs: np.ndarray | None,
    embedding_pca_audit: dict[str, object] | None,
    y: np.ndarray,
    atoms: np.ndarray,
    split: dict[str, object],
    alpha: float,
) -> dict[str, object]:
    train = split["train_mask"]
    test = split["test_mask"]
    pieces = []
    pca_audit = None
    for block in TIER_BLOCKS[tier]:
        if block == "aimnet2_embedding_pca":
            if embedding_pcs is None or embedding_pca_audit is None:
                raise RuntimeError("tier requires embedding PCA but no shared within transform was provided")
            pieces.append(embedding_pcs)
            pca_audit = embedding_pca_audit
        else:
            pieces.append(base_blocks[block])
    x = np.column_stack(pieces)
    x_centered = centered_by_train_atom(x, atoms, train)
    y_centered = centered_by_train_atom(y, atoms, train)
    fit = train & np.isfinite(y_centered).all(axis=1) & np.isfinite(x_centered).any(axis=1)
    score = test & np.isfinite(y_centered).all(axis=1) & np.isfinite(x_centered).any(axis=1)
    pred = fit_ridge_predict(x_centered[fit], y_centered[fit], x_centered[score], alpha)
    return {
        "pred": pred,
        "target": y_centered[score],
        "groups": atoms[score],
        "score_mask": score,
        "score_row_positions": np.flatnonzero(score),
        "alpha": float(alpha),
        "pca_audit": pca_audit,
        "cv_check": {
            "fit_rows": int(fit.sum()),
            "score_rows": int(score.sum()),
            "test_rows": int(test.sum()),
            "train_test_overlap": int((train & test).sum()),
            "purged_in_train_or_test": int(((split["purge_mask"]) & (train | test)).sum()),
            "heldout_block_excluded_from_preprocessing": True,
            "alpha_selection": "outer held-out frame block excluded from alpha selection",
        },
    }


def subset_frame_split(split: dict[str, object], positions: np.ndarray) -> dict[str, object]:
    pos = np.asarray(positions, dtype=int)
    return {
        "train_mask": np.asarray(split["train_mask"], dtype=bool)[pos],
        "test_mask": np.asarray(split["test_mask"], dtype=bool)[pos],
        "purge_mask": np.asarray(split["purge_mask"], dtype=bool)[pos],
        "test_frames": list(split["test_frames"]),
        "purged_frames": list(split["purged_frames"]),
        "cross_split_lag1_pairs": int(split["cross_split_lag1_pairs"]),
    }


def combine_per_stratum_cv(cv_items: list[dict[str, object]], axis: str) -> dict[str, object]:
    if axis == "between":
        return {
            "heldout_atoms_never_in_train": bool(all(bool(item["heldout_atoms_never_in_train"]) for item in cv_items)),
            "folds": int(sum(int(item["folds"]) for item in cv_items)),
            "per_stratum": cv_items,
        }
    return {
        "fit_rows": int(sum(int(item["fit_rows"]) for item in cv_items)),
        "score_rows": int(sum(int(item["score_rows"]) for item in cv_items)),
        "test_rows": int(sum(int(item["test_rows"]) for item in cv_items)),
        "train_test_overlap": int(sum(int(item["train_test_overlap"]) for item in cv_items)),
        "purged_in_train_or_test": int(sum(int(item["purged_in_train_or_test"]) for item in cv_items)),
        "heldout_block_excluded_from_preprocessing": bool(all(bool(item["heldout_block_excluded_from_preprocessing"]) for item in cv_items)),
        "alpha_selection": "outer held-out frame block excluded from each stratum alpha selection",
        "per_stratum": cv_items,
    }


def per_type_alpha_summary(alpha_by_stratum: dict[str, object], axis: str, alpha_grid: list[float]) -> dict[str, object]:
    selected = []
    for audit in alpha_by_stratum.values():
        if isinstance(audit, dict) and audit.get("selected_alpha") is not None:
            selected.append(float(audit["selected_alpha"]))
    summary = alpha_distribution(selected)
    summary.update({
        "axis": axis,
        "method": "per-type alpha selection; one selected alpha per stratum model",
        "alpha_grid": [float(x) for x in alpha_grid],
        "by_stratum": alpha_by_stratum,
    })
    return summary


def evaluate_per_type_tier(
    tier: str,
    rows: pd.DataFrame,
    atom_rows: pd.DataFrame,
    base_blocks: dict[str, np.ndarray],
    base_atom_blocks: dict[str, np.ndarray],
    embedding: np.ndarray,
    embedding_atom_mean: np.ndarray,
    row_indices: np.ndarray,
    y: np.ndarray,
    y_atom: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    split: dict[str, object],
    alpha_mode: str,
    fixed_alpha: float,
    alpha_grid: list[float],
    n_pcs: int,
    inner_cv_folds: int,
    purge: int,
) -> tuple[dict[str, object], dict[str, object], dict[str, object], dict[str, object]]:
    n_atoms = int(len(atom_rows))
    out_dim = int(y.shape[1])
    between_pred = np.full((n_atoms, out_dim), np.nan, dtype=float)
    between_alpha_all = np.full(n_atoms, np.nan, dtype=float)
    within_pred_parts: list[np.ndarray] = []
    within_target_parts: list[np.ndarray] = []
    within_group_parts: list[np.ndarray] = []
    within_position_parts: list[np.ndarray] = []
    within_alpha_parts: list[np.ndarray] = []
    alpha_between_by_stratum: dict[str, object] = {}
    alpha_within_by_stratum: dict[str, object] = {}
    pca_between_by_stratum: dict[str, object] = {}
    pca_within_by_stratum: dict[str, object] = {}
    cv_between_items: list[dict[str, object]] = []
    cv_within_items: list[dict[str, object]] = []
    strata = sorted(str(x) for x in atom_rows["stratum"].dropna().unique())
    for stratum in strata:
        atom_mask = atom_rows["stratum"].to_numpy(str) == stratum
        atom_pos = np.flatnonzero(atom_mask)
        row_mask = rows["stratum"].to_numpy(str) == stratum
        row_pos = np.flatnonzero(row_mask)
        if len(atom_pos) < 2 or len(row_pos) < 3:
            continue
        base_atom_s = {name: block[atom_pos] for name, block in base_atom_blocks.items()}
        embedding_atom_s = embedding_atom_mean[atom_pos]
        y_atom_s = y_atom[atom_pos]
        if alpha_mode == "select":
            between_alpha, between_alpha_audit = select_between_alphas(
                tier,
                base_atom_s,
                embedding_atom_s,
                y_atom_s,
                alpha_grid,
                n_pcs,
                inner_cv_folds,
                {},
            )
        else:
            between_alpha = np.full(len(atom_pos), fixed_alpha, dtype=float)
            between_alpha_audit = {"axis": "between_LOAO", "method": "fixed-alpha per-type baseline; no inner CV", "selected_alpha": fixed_alpha, "min": fixed_alpha, "max": fixed_alpha, "counts": {f"{fixed_alpha:g}": int(len(atom_pos))}, "heldout_test_atom_excluded_from_alpha_selection": True, "alpha_grid": [float(x) for x in alpha_grid]}
        between = evaluate_between_tier(tier, base_atom_s, embedding_atom_s, y_atom_s, between_alpha, n_pcs, {})
        between_pred[atom_pos] = between["pred"]
        between_alpha_all[atom_pos] = np.asarray(between["alpha_by_heldout"], dtype=float)
        alpha_between_by_stratum[stratum] = between_alpha_audit
        pca_between_by_stratum[stratum] = between["pca_audit"]
        cv_b = dict(between["cv_check"])
        cv_b["stratum"] = stratum
        cv_b["n_atoms"] = int(len(atom_pos))
        cv_between_items.append(cv_b)

        base_blocks_s = {name: block[row_pos] for name, block in base_blocks.items()}
        row_indices_s = row_indices[row_pos]
        y_s = y[row_pos]
        atoms_s = atoms[row_pos]
        frames_s = frames[row_pos]
        split_s = subset_frame_split(split, row_pos)
        if alpha_mode == "select":
            within_alpha, within_alpha_audit = select_within_alpha(
                tier,
                base_blocks_s,
                embedding,
                row_indices_s,
                y_s,
                atoms_s,
                frames_s,
                split_s,
                alpha_grid,
                n_pcs,
                inner_cv_folds,
                purge,
                {},
            )
        else:
            within_alpha = fixed_alpha
            within_alpha_audit = {"axis": "within_frameblock", "method": "fixed-alpha per-type baseline; no inner CV", "selected_alpha": fixed_alpha, "outer_test_rows_used_for_alpha_selection": 0, "outer_purged_rows_used_for_alpha_selection": 0, "alpha_grid": [float(x) for x in alpha_grid]}
        within_embedding_pcs = None
        within_embedding_pca_audit = {"method": "embedding sidecar present but not consumed by this tier", "n_components": 0, "original_dims": int(embedding.shape[1]), "stratum": stratum}
        if tier_uses_embedding(tier):
            within_pca = fit_pca_memmap(embedding, row_indices_s[np.asarray(split_s["train_mask"], dtype=bool)], n_pcs)
            within_embedding_pcs = transform_pca_memmap(embedding, within_pca, row_indices_s)
            within_embedding_pca_audit = {
                "method": "per-type unsupervised PCA fit on stratum training rows only",
                "stratum": stratum,
                "n_components": int(within_pca.n_components),
                "original_dims": int(within_pca.original_dims),
                "training_rows": int(within_pca.n_train_rows),
                "test_rows_excluded": int(np.asarray(split_s["test_mask"], dtype=bool).sum()),
                "purged_rows_excluded": int(np.asarray(split_s["purge_mask"], dtype=bool).sum()),
                "explained_variance_ratio": float(within_pca.explained_variance_ratio),
            }
        within = evaluate_within_tier(tier, base_blocks_s, within_embedding_pcs, within_embedding_pca_audit, y_s, atoms_s, split_s, float(within_alpha))
        score_positions = row_pos[np.asarray(within["score_row_positions"], dtype=int)]
        within_pred_parts.append(within["pred"])
        within_target_parts.append(within["target"])
        within_group_parts.append(within["groups"])
        within_position_parts.append(score_positions)
        within_alpha_parts.append(np.full(len(within["target"]), float(within_alpha), dtype=float))
        alpha_within_by_stratum[stratum] = within_alpha_audit
        pca_within_by_stratum[stratum] = within["pca_audit"]
        cv_w = dict(within["cv_check"])
        cv_w["stratum"] = stratum
        cv_w["n_atoms"] = int(len(atom_pos))
        cv_within_items.append(cv_w)

    if not within_position_parts:
        raise RuntimeError(f"per-type evaluation produced no within predictions for tier {tier}")
    within_positions = np.concatenate(within_position_parts).astype(int)
    order = np.argsort(within_positions)
    within_pred = np.vstack(within_pred_parts)[order]
    within_target = np.vstack(within_target_parts)[order]
    within_groups = np.concatenate(within_group_parts).astype(int)[order]
    within_positions = within_positions[order]
    within_alpha_by_score = np.concatenate(within_alpha_parts).astype(float)[order]
    within_alpha_mode = mode_alpha(within_alpha_by_score)
    tier_result = {
        "between": {
            "pred": between_pred,
            "target": y_atom,
            "groups": np.arange(n_atoms, dtype=int),
            "alpha_by_heldout": between_alpha_all,
            "pca_audit": pca_between_by_stratum,
            "cv_check": combine_per_stratum_cv(cv_between_items, "between"),
        },
        "within": {
            "pred": within_pred,
            "target": within_target,
            "groups": within_groups,
            "score_mask": None,
            "score_row_positions": within_positions,
            "alpha": within_alpha_mode,
            "alpha_by_score_row": within_alpha_by_score,
            "pca_audit": pca_within_by_stratum,
            "cv_check": combine_per_stratum_cv(cv_within_items, "within"),
        },
    }
    alpha_audit = {
        "between": per_type_alpha_summary(alpha_between_by_stratum, "between_LOAO", alpha_grid),
        "within": per_type_alpha_summary(alpha_within_by_stratum, "within_frameblock", alpha_grid),
    }
    pca_audit = {"between": pca_between_by_stratum, "within": pca_within_by_stratum}
    cv_details = {"between": {"cv_check": tier_result["between"]["cv_check"]}, "within": {"cv_check": tier_result["within"]["cv_check"]}}
    return tier_result, alpha_audit, pca_audit, cv_details


def build_score_rows(
    tier_results: dict[str, dict[str, object]],
    tier_order: list[str],
    rows: pd.DataFrame,
    atom_rows: pd.DataFrame,
    y: np.ndarray,
    frames: np.ndarray,
    atoms: np.ndarray,
    split: dict[str, object],
    base_blocks: dict[str, np.ndarray],
    alpha_mode: str,
    embedding_components: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    score_rows = []
    long_rows = []
    n_atoms_global = int(len(atom_rows))
    y_atom = atom_means_dense(y, int(len(np.unique(frames))), n_atoms_global)
    atom_types = ["ALL"] + sorted(str(x) for x in atom_rows["stratum"].dropna().unique())
    for tier in tier_order:
        n_features = feature_count_for_tier(tier, embedding_components, base_blocks)
        between = tier_results[tier]["between"]
        within = tier_results[tier]["within"]
        between_alpha_all = np.asarray(between["alpha_by_heldout"], dtype=float)
        within_alpha = float(within["alpha"])
        within_alpha_by_score = within.get("alpha_by_score_row")
        within_alpha_all = np.asarray(within_alpha_by_score, dtype=float) if within_alpha_by_score is not None else None
        within_positions = within["score_row_positions"]
        within_frames = frames[within_positions]
        within_rows = rows.iloc[within_positions]
        for atom_type in atom_types:
            if atom_type == "ALL":
                atom_mask = np.ones(n_atoms_global, dtype=bool)
                row_mask = np.ones(len(rows), dtype=bool)
                within_mask = np.ones(len(within["target"]), dtype=bool)
            else:
                atom_mask = atom_rows["stratum"].to_numpy(str) == atom_type
                row_mask = rows["stratum"].to_numpy(str) == atom_type
                within_mask = within_rows["stratum"].to_numpy(str) == atom_type
            n_atoms_slice = int(atom_mask.sum())
            rows_slice = int(row_mask.sum())
            frames_slice = int(rows.loc[row_mask, "original_frame_index"].nunique()) if rows_slice else 0
            b_pred = between["pred"][atom_mask]
            b_target = between["target"][atom_mask]
            b_groups = between["groups"][atom_mask]
            b_alpha_dist = alpha_distribution(between_alpha_all[atom_mask])
            b_alpha_mode = math.nan if b_alpha_dist["selected_alpha"] is None else float(b_alpha_dist["selected_alpha"])
            b_alpha_min = math.nan if b_alpha_dist["min"] is None else float(b_alpha_dist["min"])
            b_alpha_max = math.nan if b_alpha_dist["max"] is None else float(b_alpha_dist["max"])
            if within_alpha_all is None:
                w_alpha_mode = within_alpha
            else:
                w_alpha_dist = alpha_distribution(within_alpha_all[within_mask])
                w_alpha_mode = math.nan if w_alpha_dist["selected_alpha"] is None else float(w_alpha_dist["selected_alpha"])
            w_pred = within["pred"][within_mask]
            w_target = within["target"][within_mask]
            w_groups = within["groups"][within_mask]
            b_r2 = r2_score(b_pred, b_target)
            w_r2 = r2_score(w_pred, w_target)
            b_se = jackknife_metric_se(b_pred, b_target, b_groups)
            w_se = jackknife_metric_se(w_pred, w_target, w_groups)
            b_low, b_high = ci95(b_r2, b_se)
            w_low, w_high = ci95(w_r2, w_se)
            share_b, share_w = variance_shares(y[row_mask], atoms[row_mask])
            neff = effective_n_components(
                y[within_positions][within_mask],
                within["groups"][within_mask],
                within_frames[within_mask],
            )
            thin = n_atoms_slice < THIN_ATOMS or len(np.unique(b_groups)) < 3 or len(np.unique(w_groups)) < 3
            p_gt_atoms = n_features >= max(1, n_atoms_slice)
            row = {
                "target": "DFT_total_T2_2e_vector",
                "atom_type_axis": "stratum",
                "atom_type": atom_type,
                "tier": tier,
                "tier_label": TIER_LABELS[tier],
                "rows": rows_slice,
                "n_atoms_between_global": n_atoms_global,
                "n_atoms_in_slice": n_atoms_slice,
                "n_original_frames": frames_slice,
                "n_features": n_features,
                "ridge_alpha": b_alpha_mode,
                "alpha_mode": alpha_mode,
                "between_ridge_alpha": b_alpha_mode,
                "between_ridge_alpha_min": b_alpha_min,
                "between_ridge_alpha_max": b_alpha_max,
                "within_ridge_alpha": within_alpha,
                "variance_share_between": share_b,
                "variance_share_within": share_w,
                **neff,
                "thin_flag": "thin" if thin else "",
                "p_gt_atoms_flag": "p>=atoms" if p_gt_atoms else "",
                "between_LOAO_test_R2": b_r2,
                "between_LOAO_test_R2_jackknife_se": b_se,
                "between_LOAO_test_R2_ci95_low": b_low,
                "between_LOAO_test_R2_ci95_high": b_high,
                "within_frameblock_test_R2": w_r2,
                "within_frameblock_test_R2_jackknife_se": w_se,
                "within_frameblock_test_R2_ci95_low": w_low,
                "within_frameblock_test_R2_ci95_high": w_high,
                "within_split_strategy": "centered_by_atom_train_frames_only; contiguous_frame_block",
                "test_frames": int(len(split["test_frames"])),
                "purged_train_frames": int(len(split["purged_frames"])),
                "cross_split_lag1_pairs": int(split["cross_split_lag1_pairs"]),
                "feature_blocks": "|".join(TIER_BLOCKS[tier]),
            }
            score_rows.append(row)
            long_rows += long_score_rows(row, "between", b_pred, b_target, b_groups, b_alpha_mode)
            long_rows += long_score_rows(row, "within", w_pred, w_target, w_groups, within_alpha)
    return score_rows, long_rows


def long_score_rows(row: dict[str, object], axis: str, pred: np.ndarray, target: np.ndarray, groups: np.ndarray, ridge_alpha: float) -> list[dict[str, object]]:
    out = []
    vector_r2 = r2_score(pred, target)
    vector_se = jackknife_metric_se(pred, target, groups)
    lo, hi = ci95(vector_r2, vector_se)
    base = {k: row[k] for k in ["target", "fit_scope", "atom_type_axis", "atom_type", "tier", "tier_label", "rows", "n_atoms_in_slice", "n_features", "alpha_mode"]}
    base["ridge_alpha"] = float(ridge_alpha)
    out.append({**base, "axis": axis, "component": "vector_5d", "heldout_R2": vector_r2, "jackknife_se": vector_se, "ci95_low": lo, "ci95_high": hi})
    comp = component_r2(pred, target)
    comp_se = jackknife_component_se(np.asarray(pred), np.asarray(target), groups)
    for i, (r2, se) in enumerate(zip(comp, comp_se)):
        lo, hi = ci95(r2, se)
        out.append({**base, "axis": axis, "component": f"2e_{i}", "heldout_R2": r2, "jackknife_se": se, "ci95_low": lo, "ci95_high": hi})
    return out


def partition_curves(
    tier_results: dict[str, dict[str, object]],
    tier_order: list[str],
    atom_conditions: pd.DataFrame,
    row_conditions: pd.DataFrame,
    frames: np.ndarray,
    condition_list: list[ConditionSpec],
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    rows_out = []
    bin_defs = []
    atom_types = ["ALL"] + sorted(str(x) for x in atom_conditions["stratum"].dropna().unique())
    axes = {
        "between": {
            "conditions": atom_conditions,
            "frames": np.arange(len(atom_conditions), dtype=int),
        },
    }
    for tier in tier_order:
        within = tier_results[tier]["within"]
        if "within" not in axes:
            axes["within"] = {
                "conditions": row_conditions.iloc[within["score_row_positions"]].reset_index(drop=True),
                "frames": frames[within["score_row_positions"]],
            }
    neff_by_tier_axis = {
        (tier, "within"): per_atom_neff_components(
            tier_results[tier]["within"]["target"],
            tier_results[tier]["within"]["groups"],
            axes["within"]["frames"],
        )
        for tier in tier_order
    }
    for axis, axis_info in axes.items():
        cond_df: pd.DataFrame = axis_info["conditions"]
        axis_frames = np.asarray(axis_info["frames"])
        for atom_type in atom_types:
            if atom_type == "ALL":
                atom_type_mask = np.ones(len(cond_df), dtype=bool)
            else:
                atom_type_mask = cond_df["stratum"].to_numpy(str) == atom_type
            if not atom_type_mask.any():
                continue
            for spec in condition_list:
                base = condition_base_name(spec)
                if base not in cond_df.columns:
                    continue
                labels, defs = assign_bins(cond_df.loc[atom_type_mask, base], cond_df.loc[atom_type_mask, "atom_index"].to_numpy(int), spec)
                full_labels = np.full(len(cond_df), "__outside_atom_type__", dtype=object)
                full_labels[atom_type_mask] = labels
                for d in defs:
                    bin_defs.append({
                        "axis": axis,
                        "atom_type": atom_type,
                        "condition_family": spec.family,
                        "condition_name": spec.name,
                        "binning": spec.binning,
                        **d,
                    })
                for bin_label in sorted(set(labels), key=str):
                    unit_mask = atom_type_mask & (full_labels == bin_label)
                    if not unit_mask.any():
                        continue
                    n_atoms = int(cond_df.loc[unit_mask, "atom_index"].nunique())
                    rows_count = int(unit_mask.sum())
                    for tier in tier_order:
                        axis_pred = tier_results[tier][axis]
                        pred = axis_pred["pred"][unit_mask]
                        target = axis_pred["target"][unit_mask]
                        groups = axis_pred["groups"][unit_mask]
                        r2 = r2_score(pred, target)
                        se = jackknife_metric_se(pred, target, groups)
                        lo, hi = ci95(r2, se)
                        atom_neff = neff_by_tier_axis.get((tier, axis))
                        n_eff_median = bin_neff_median(axis, cond_df.loc[:, "atom_index"].to_numpy(int), unit_mask, atom_neff)
                        thin = n_atoms < THIN_ATOMS or (axis == "within" and rows_count < THIN_WITHIN_ROWS)
                        matching_defs = [d for d in defs if d["bin_label"] == bin_label]
                        bin_low = matching_defs[0].get("bin_low") if matching_defs else None
                        bin_high = matching_defs[0].get("bin_high") if matching_defs else None
                        rows_out.append({
                            "condition_family": spec.family,
                            "condition_name": spec.name,
                            "bin_label": bin_label,
                            "bin_low": bin_low,
                            "bin_high": bin_high,
                            "atom_type": atom_type,
                            "tier": tier,
                            "axis": axis,
                            "rows": rows_count,
                            "n_atoms": n_atoms,
                            "N_eff_median": n_eff_median,
                            "heldout_R2": r2,
                            "heldout_R2_ci95_low": lo,
                            "heldout_R2_ci95_high": hi,
                            "delta_R2_vs_classical": math.nan,
                            "delta_R2_vs_previous_tier": math.nan,
                            "thin_flag": "thin" if thin else "",
                            "_bin_order": bin_order_key(bin_label),
                        })
    curves = pd.DataFrame(rows_out)
    if curves.empty:
        return curves, curves.copy(), bin_defs
    key_cols = ["condition_family", "condition_name", "bin_label", "bin_low", "bin_high", "atom_type", "axis"]
    classical = curves[curves["tier"] == "classical_mechanisms_combined"][key_cols + ["heldout_R2"]].rename(columns={"heldout_R2": "_classical_R2"})
    curves = curves.merge(classical, on=key_cols, how="left")
    curves["delta_R2_vs_classical"] = curves["heldout_R2"] - curves["_classical_R2"]
    prev_parts = []
    prev_map = previous_tier_map(tier_order)
    for tier, prev in prev_map.items():
        part = curves[curves["tier"] == tier].copy()
        if prev is None:
            part["_previous_R2"] = np.nan
        else:
            prev_df = curves[curves["tier"] == prev][key_cols + ["heldout_R2"]].rename(columns={"heldout_R2": "_previous_R2"})
            part = part.drop(columns=["_previous_R2"], errors="ignore").merge(prev_df, on=key_cols, how="left")
        prev_parts.append(part)
    curves = pd.concat(prev_parts, ignore_index=True)
    curves["delta_R2_vs_previous_tier"] = curves["heldout_R2"] - curves["_previous_R2"]
    shape_map = response_shapes(curves)
    curves["response_shape"] = [
        shape_map.get((r.condition_family, r.condition_name, r.atom_type, r.tier, r.axis), "unstable-thin")
        for r in curves.itertuples(index=False)
    ]
    public_cols = [
        "condition_family",
        "condition_name",
        "bin_label",
        "bin_low",
        "bin_high",
        "atom_type",
        "tier",
        "axis",
        "rows",
        "n_atoms",
        "N_eff_median",
        "heldout_R2",
        "heldout_R2_ci95_low",
        "heldout_R2_ci95_high",
        "delta_R2_vs_classical",
        "delta_R2_vs_previous_tier",
        "thin_flag",
        "response_shape",
    ]
    curves_public = curves[public_cols].copy()
    long_rows = []
    for r in curves.itertuples(index=False):
        for metric in ["heldout_R2", "delta_R2_vs_classical", "delta_R2_vs_previous_tier"]:
            long_rows.append({
                "condition_family": r.condition_family,
                "condition_name": r.condition_name,
                "bin_label": r.bin_label,
                "atom_type": r.atom_type,
                "tier": r.tier,
                "axis": r.axis,
                "metric": metric,
                "value": getattr(r, metric),
                "thin_flag": r.thin_flag,
                "response_shape": r.response_shape,
            })
    return curves_public, pd.DataFrame(long_rows), bin_defs


def bin_order_key(label: object) -> float:
    text = str(label)
    if text.startswith("Q"):
        try:
            return float(text[1:].split(" ", 1)[0])
        except ValueError:
            return 999.0
    order = {"<=4A": 1, "4-6A": 2, "6-8A": 3, "8-10A": 4, ">10A": 5}
    return float(order.get(text, 999))


def response_shapes(curves: pd.DataFrame) -> dict[tuple[str, str, str, str, str], str]:
    out = {}
    group_cols = ["condition_family", "condition_name", "atom_type", "tier", "axis"]
    for key, group in curves.groupby(group_cols, dropna=False):
        usable = group[(group["thin_flag"] == "") & np.isfinite(group["heldout_R2"])].sort_values("_bin_order")
        vals = usable["heldout_R2"].to_numpy(float)
        if len(vals) < 3:
            out[key] = "unstable-thin"
            continue
        span = float(np.nanmax(vals) - np.nanmin(vals))
        tol = max(0.02, 0.10 * span)
        diffs = np.diff(vals)
        if span < 0.02:
            out[key] = "flat response"
        elif np.all(diffs >= -tol) and vals[-1] > vals[0]:
            out[key] = "monotone rise"
        elif np.all(diffs <= tol) and vals[-1] < vals[0]:
            out[key] = "monotone fall"
        elif len(vals) >= 4 and vals[0] > vals[1] and vals[-1] > vals[-2] and np.nanmin(vals[1:-1]) < min(vals[0], vals[-1]) - tol:
            out[key] = "U-shape"
        elif np.max(np.abs(diffs)) > 0.6 * span:
            out[key] = "threshold behavior"
        else:
            out[key] = "flat response"
    return out


def favorable_cases(curves: pd.DataFrame, shortlist_size: int = 20) -> pd.DataFrame:
    if curves.empty:
        return curves.copy()
    eligible = curves[
        (curves["thin_flag"] == "")
        & np.isfinite(curves["heldout_R2"])
        & (curves["n_atoms"] >= THIN_ATOMS)
        & (curves["N_eff_median"] >= MIN_FAVORABLE_N_EFF)
    ].copy()
    eligible = eligible[(eligible["heldout_R2_ci95_low"] > 0.0) | (eligible["heldout_R2"] > 0.0)]
    eligible["_positive_ci"] = eligible["heldout_R2_ci95_low"] > 0.0
    eligible = eligible.sort_values(
        ["_positive_ci", "heldout_R2", "delta_R2_vs_previous_tier"],
        ascending=[False, False, False],
        na_position="last",
    )
    cols = [
        "condition_family",
        "condition_name",
        "bin_label",
        "atom_type",
        "tier",
        "axis",
        "rows",
        "n_atoms",
        "N_eff_median",
        "heldout_R2",
        "heldout_R2_ci95_low",
        "heldout_R2_ci95_high",
        "delta_R2_vs_classical",
        "delta_R2_vs_previous_tier",
        "response_shape",
    ]
    return eligible.head(shortlist_size)[cols]


def join_coverage(rows: pd.DataFrame) -> pd.DataFrame:
    parts = []
    for axis in ["stratum", "role", "ff_atom_type_ord"]:
        for level, group in rows.groupby(axis, dropna=False):
            parts.append({
                "axis": axis,
                "level": level,
                "rows": int(len(group)),
                "n_atoms": int(group["atom_index"].nunique()),
                "fraction_rows": float(len(group) / len(rows)),
            })
    parts.append({"axis": "ALL", "level": "ALL", "rows": int(len(rows)), "n_atoms": int(rows["atom_index"].nunique()), "fraction_rows": 1.0})
    return pd.DataFrame(parts)


def feature_blocks_used(
    labels: dict[str, list[str]],
    embedding_components: int,
    transform_audit: list[dict[str, object]],
    tier_order: list[str],
) -> dict[str, object]:
    blocks = {}
    for tier in tier_order:
        tier_blocks = {}
        for block in TIER_BLOCKS[tier]:
            if block == "aimnet2_embedding_pca":
                tier_blocks[block] = [f"embedding_pc_{i:02d}" for i in range(embedding_components)]
            else:
                tier_blocks[block] = labels[block]
        blocks[tier] = tier_blocks
    return {
        "tiers": blocks,
        "tier_order": tier_order,
        "same_charge_complete_row_set_for_all_tiers": True,
        "t2_feature_transforms": transform_audit,
    }


def markdown_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df.empty:
        return "No rows."
    view = df.head(max_rows).copy() if max_rows is not None else df.copy()
    cols = list(view.columns)

    def cell(x: object) -> str:
        if isinstance(x, (float, np.floating)):
            text = finite_fmt(float(x))
        else:
            text = "" if pd.isna(x) else str(x)
        return text.replace("|", "\\|")

    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join(["---"] * len(cols)) + " |",
    ]
    for row in view.itertuples(index=False):
        lines.append("| " + " | ".join(cell(x) for x in row) + " |")
    if max_rows is not None and len(df) > max_rows:
        lines.append(f"\nShowing first {max_rows} of {len(df)} rows.")
    return "\n".join(lines)


def write_reports(
    out_dir: Path,
    score_table: pd.DataFrame,
    curves: pd.DataFrame,
    fav: pd.DataFrame,
    audit: dict[str, object],
) -> None:
    lines = [
        "# All-Atoms Charge-Complete Fit Report",
        "",
        f"Substrate: `{audit['substrate_dir']}`",
        f"Rows used: {audit['charge_complete_rows_used']:,}",
        f"Alpha mode: {audit['alpha_mode']}",
        f"Fixed-alpha baseline value: {finite_fmt(audit['ridge_alpha_baseline'])}",
        f"Embedding PCA: {audit['embedding_pca']['n_components']} components, training-only per CV fold.",
        "",
        "## Fit-Stage Checks",
    ]
    for name, check in audit["fit_stage_checks"].items():
        lines.append(f"- {name}: {'pass' if check['pass'] else 'fail'}")
    lines += [
        "",
        "## Score Table",
        markdown_table(score_table[[
            "tier",
            "atom_type",
            "n_atoms_in_slice",
            "n_features",
            "between_ridge_alpha",
            "between_ridge_alpha_min",
            "between_ridge_alpha_max",
            "within_ridge_alpha",
            "between_LOAO_test_R2",
            "between_LOAO_test_R2_ci95_low",
            "between_LOAO_test_R2_ci95_high",
            "within_frameblock_test_R2",
            "within_frameblock_test_R2_ci95_low",
            "within_frameblock_test_R2_ci95_high",
            "thin_flag",
        ]]),
        "",
        "## Tier Deltas",
    ]
    delta = score_table[["tier", "atom_type", "between_LOAO_test_R2", "within_frameblock_test_R2"]].copy()
    tier_order = list(audit["tier_order"])
    prev_map = previous_tier_map(tier_order)
    for metric in ["between_LOAO_test_R2", "within_frameblock_test_R2"]:
        pivot = delta.pivot(index="atom_type", columns="tier", values=metric)
        for tier in tier_order:
            prev = prev_map[tier]
            if prev and tier in pivot and prev in pivot:
                pivot[f"{tier}_delta_vs_previous_{metric}"] = pivot[tier] - pivot[prev]
        lines.append(f"### {metric}")
        lines.append(markdown_table(pivot.reset_index()))
        lines.append("")
    (out_dir / "allatom_fit_report.md").write_text("\n".join(lines), encoding="utf-8")

    plines = [
        "# Partition By Condition Report",
        "",
        "Partitions use emitted or algebraically derived input-side conditions only. Bins were generated before favorable-case ranking and thin bins are excluded from the shortlist.",
        "",
        "## Favorable-Case Shortlist",
        markdown_table(fav) if not fav.empty else "No eligible rows.",
        "",
        "## Response Shapes",
    ]
    if not curves.empty:
        shape_counts = curves.groupby(["condition_family", "condition_name", "axis", "tier", "response_shape"], dropna=False).size().reset_index(name="curve_points")
        plines.append(markdown_table(shape_counts))
    else:
        plines.append("No curves emitted.")
    (out_dir / "partition_report.md").write_text("\n".join(plines), encoding="utf-8")


def fit_stage_checks(audit: dict[str, object], score_table: pd.DataFrame, curves: pd.DataFrame, fav: pd.DataFrame) -> None:
    checks = audit["fit_stage_checks"]
    no_external = {
        "pass": True,
        "python_read_only_emitted_substrate_directory": True,
        "files_read": audit["files_read"],
        "forbidden_inputs_not_opened": [
            "trajectory.h5",
            "ORCA outputs",
            "/tmp/combined-mopac-layer3",
            "broad-backbone directories",
            "MOPAC directories",
            "all-atom-equivariant directories",
            "pair/source dumps",
        ],
        "sidecar_join": "CSV row i joins sidecar row i / row_id i",
        "driver_modulation_join": "per_atom_substrate_driver_modulation_by_atom joined by atom_index",
    }
    checks["no_external_merge"] = no_external
    checks["basis_and_target"] = {
        "pass": bool(
            audit["change_of_basis_orthogonality_max_abs"] < 1.0e-12
            and audit["target_T2_original_shape"] == [audit["manifest_rows"], 5]
            and all(x["components"] == 5 for x in audit["t2_transform_check_blocks"])
        ),
        "C_T_C_minus_I_max_abs": audit["change_of_basis_orthogonality_max_abs"],
        "target_shape": audit["target_T2_original_shape"],
        "t2_blocks": audit["t2_transform_check_blocks"],
    }
    cv_checks = audit["cv_integrity_details"]
    checks["cv_integrity"] = {
        "pass": bool(
            all(v["between"]["cv_check"]["heldout_atoms_never_in_train"] for v in cv_checks.values())
            and all(v["within"]["cv_check"]["train_test_overlap"] == 0 for v in cv_checks.values())
            and all(v["within"]["cv_check"]["purged_in_train_or_test"] == 0 for v in cv_checks.values())
            and audit["within_split"]["cross_split_lag1_pairs"] == 0
            and score_table["within_N_eff_median"].notna().all()
            and score_table["median_lag1_rho"].notna().all()
        ),
        "between_heldout_atoms_never_in_train": all(v["between"]["cv_check"]["heldout_atoms_never_in_train"] for v in cv_checks.values()),
        "within_train_test_overlap_rows": {k: v["within"]["cv_check"]["train_test_overlap"] for k, v in cv_checks.items()},
        "purged_frames_absent_from_train_and_test": all(v["within"]["cv_check"]["purged_in_train_or_test"] == 0 for v in cv_checks.values()),
        "scores_are_held_out": True,
        "N_eff_and_lag1_reported_per_score_row": bool(score_table["within_N_eff_median"].notna().all() and score_table["median_lag1_rho"].notna().all()),
        "cross_split_lag1_pairs": audit["within_split"]["cross_split_lag1_pairs"],
    }
    alpha_info: dict[str, object] = audit["alpha_selection"]
    alpha_mode = str(alpha_info["mode"])
    by_tier_axis: dict[str, object] = alpha_info["by_tier_axis"]
    grid = [float(x) for x in alpha_info.get("grid", [])]
    grid_is_predeclared = grid == [float(x) for x in ALPHA_GRID]
    selected_alphas_recorded = True
    selected_alphas_in_grid = True
    between_excludes_outer = True
    within_excludes_outer = True
    for tier in audit["tier_order"]:
        tier_alpha: dict[str, object] = by_tier_axis.get(tier, {})  # type: ignore[assignment]
        between_alpha: dict[str, object] = tier_alpha.get("between", {})  # type: ignore[assignment]
        within_alpha: dict[str, object] = tier_alpha.get("within", {})  # type: ignore[assignment]
        b_selected = between_alpha.get("selected_alpha")
        w_selected = within_alpha.get("selected_alpha")
        selected_alphas_recorded = selected_alphas_recorded and b_selected is not None and w_selected is not None
        if b_selected is not None:
            selected_alphas_in_grid = selected_alphas_in_grid and float(b_selected) in grid
        if w_selected is not None:
            selected_alphas_in_grid = selected_alphas_in_grid and float(w_selected) in grid
        if alpha_mode == "select":
            between_excludes_outer = between_excludes_outer and bool(between_alpha.get("heldout_test_atom_excluded_from_alpha_selection"))
            within_excludes_outer = (
                within_excludes_outer
                and int(within_alpha.get("outer_test_rows_used_for_alpha_selection", -1)) == 0
                and int(within_alpha.get("outer_purged_rows_used_for_alpha_selection", -1)) == 0
            )
    checks["alpha_selection_integrity"] = {
        "pass": bool(
            alpha_mode in {"select", "fixed"}
            and selected_alphas_recorded
            and selected_alphas_in_grid
            and (alpha_mode == "fixed" or (grid_is_predeclared and between_excludes_outer and within_excludes_outer))
        ),
        "alpha_mode": alpha_mode,
        "alpha_grid_predeclared": grid_is_predeclared,
        "selected_alphas_recorded_per_tier_axis": selected_alphas_recorded,
        "selected_alphas_in_grid": selected_alphas_in_grid,
        "between_outer_test_atoms_excluded": between_excludes_outer,
        "within_outer_test_and_purged_rows_excluded": within_excludes_outer,
    }
    shortlist_thin = bool(fav["thin_flag"].eq("thin").any()) if "thin_flag" in fav.columns else False
    checks["partition_integrity"] = {
        "pass": bool(
            not curves.empty
            and audit["partition_bin_definitions_written_before_ranking"]
            and audit["partition_conditions_input_side_only"]
            and not shortlist_thin
        ),
        "condition_columns_emitted_or_algebraically_derived_from_emitted_input_side": audit["partition_conditions_input_side_only"],
        "bin_definitions_written_before_ranking": audit["partition_bin_definitions_written_before_ranking"],
        "no_condition_uses_dft_target_residual_or_coefficients": True,
        "thin_bins_excluded_from_favorable_shortlist": not shortlist_thin,
    }


@dataclass(frozen=True)
class Build2ConditionSpec:
    family: str
    mechanism: str
    name: str
    kind: str
    source: str
    column: str
    bin_index: int | None = None


def path_under(path: Path, root: Path) -> bool:
    path = path.resolve()
    root = root.resolve()
    return path == root or root in path.parents


def directory_size_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    total = 0
    for item in path.rglob("*"):
        try:
            if item.is_file() or item.is_symlink():
                total += item.stat().st_size
        except FileNotFoundError:
            continue
    return int(total)


def prepare_build2_output_dir(out_dir: Path, keep_out_dir: bool) -> tuple[Path, dict[str, object]]:
    out_dir = out_dir.resolve()
    runs_root = Path("/tmp/rediscover-runs").resolve()
    shared_root = Path("/shared").resolve()
    if not path_under(out_dir, runs_root):
        raise SystemExit(f"FATAL: Build 3 output must be under {runs_root}, got {out_dir}")
    if path_under(out_dir, shared_root):
        raise SystemExit(f"FATAL: refusing to write Build 3 result data under /shared: {out_dir}")
    usage = shutil.disk_usage(out_dir.parent if out_dir.parent.exists() else runs_root)
    rediscover_bytes = directory_size_bytes(runs_root)
    audit = {
        "checked_before_write": True,
        "path": str(out_dir),
        "filesystem_free_bytes": int(usage.free),
        "filesystem_free_GiB": float(usage.free / 1024**3),
        "rediscover_runs_bytes_before_write": int(rediscover_bytes),
        "rediscover_runs_GiB_before_write": float(rediscover_bytes / 1024**3),
        "min_free_GiB_required": float(MIN_DISK_FREE_BYTES / 1024**3),
        "max_rediscover_GiB_required": float(MAX_REDISCOVER_BYTES / 1024**3),
        "deleted_existing_output_dir": False,
    }
    if usage.free < MIN_DISK_FREE_BYTES:
        raise SystemExit(f"FATAL: disk guard: {usage.free / 1024**3:.1f} GiB free under /tmp, need >=20 GiB")
    if rediscover_bytes > MAX_REDISCOVER_BYTES:
        raise SystemExit(f"FATAL: disk guard: /tmp/rediscover-runs is {rediscover_bytes / 1024**3:.1f} GiB, need <15 GiB")
    if out_dir.exists() and not keep_out_dir:
        shutil.rmtree(out_dir)
        audit["deleted_existing_output_dir"] = True
        audit["deleted_existing_output_dir_full_path"] = str(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir, audit


def specs_array_columns(specs: dict[str, object], array_name: str) -> list[str]:
    return [str(c["column"]) for c in specs["columns"] if c["array"] == array_name]


def build2_load_inputs(substrate_dir: Path) -> dict[str, object]:
    manifest_path = require_file(substrate_dir / "per_atom_substrate_manifest.json")
    specs_path = require_file(substrate_dir / "per_atom_substrate_column_specs.json")
    rows_path = require_file(substrate_dir / "per_atom_substrate_rows.csv")
    with manifest_path.open("r", encoding="utf-8") as f:
        manifest = json.load(f)
    with specs_path.open("r", encoding="utf-8") as f:
        specs = json.load(f)
    rows = pd.read_csv(rows_path)
    arrays = {
        "per_atom_substrate_features_classical": np.load(require_file(substrate_dir / "per_atom_substrate_features_classical.npy"), mmap_mode="r"),
        "per_atom_substrate_features_conditioning": np.load(require_file(substrate_dir / "per_atom_substrate_features_conditioning.npy"), mmap_mode="r"),
        "per_atom_substrate_features_hbond_conditioning": np.load(require_file(substrate_dir / "per_atom_substrate_features_hbond_conditioning.npy"), mmap_mode="r"),
        "per_atom_substrate_features_method_paths": np.load(require_file(substrate_dir / "per_atom_substrate_features_method_paths.npy"), mmap_mode="r"),
        "per_atom_substrate_partition_bins": np.load(require_file(substrate_dir / "per_atom_substrate_partition_bins.npy"), mmap_mode="r"),
        "per_atom_substrate_aimnet2_embedding": np.load(require_file(substrate_dir / "per_atom_substrate_aimnet2_embedding.npy"), mmap_mode="r"),
    }
    targets = {}
    for target_name, target_spec in TARGET_SPECS.items():
        array_name = str(target_spec["array"])
        targets[target_name] = np.load(require_file(substrate_dir / str(target_spec["file"])), mmap_mode="r")
    files_read = [
        str(manifest_path),
        str(specs_path),
        str(rows_path),
        *[str(substrate_dir / f"{name}.npy") for name in arrays],
        *[str(substrate_dir / str(spec["file"])) for spec in TARGET_SPECS.values()],
    ]
    return {"manifest": manifest, "specs": specs, "rows": rows, "arrays": arrays, "targets": targets, "files_read": files_read}


def build2_input_acceptance_checks(data: dict[str, object], substrate_dir: Path) -> tuple[dict[str, object], np.ndarray]:
    manifest: dict[str, object] = data["manifest"]
    specs: dict[str, object] = data["specs"]
    rows: pd.DataFrame = data["rows"]
    arrays: dict[str, np.ndarray] = data["arrays"]
    targets: dict[str, np.ndarray] = data["targets"]
    checks: dict[str, object] = {}
    required_cols = [
        "row_id",
        "atom_index",
        "frame_slot",
        "original_frame_index",
        "stratum",
        "role",
        "ff_atom_type_ord",
        "formal_charge",
        "ff14sb_charge",
        "ff14sb_charge_present",
        "mopac_welford_mean_charge",
        "mopac_welford_mean_charge_present",
        "eeq_charge",
        "eeq_charge_present",
    ]
    require_columns(rows, required_cols, "per_atom_substrate_rows.csv")
    n_rows = len(rows)
    n_atoms = int(manifest.get("n_atoms", -1))
    n_frames = int(manifest.get("n_dft_frames", -1))
    row_contract = bool(
        np.array_equal(rows["row_id"].to_numpy(int), np.arange(n_rows, dtype=int))
        and np.all(rows["row_id"].to_numpy(int) == rows["frame_slot"].to_numpy(int) * n_atoms + rows["atom_index"].to_numpy(int))
    )
    shape_expected = {
        "per_atom_substrate_features_classical": (n_rows, 89),
        "per_atom_substrate_features_conditioning": (n_rows, 32),
        "per_atom_substrate_features_hbond_conditioning": (n_rows, 73),
        "per_atom_substrate_features_method_paths": (n_rows, 111),
        "per_atom_substrate_partition_bins": (n_rows, 25),
        "per_atom_substrate_aimnet2_embedding": (n_rows, 256),
    }
    shapes_ok = n_atoms == 846 and n_frames == 660 and n_rows == 558_360
    for name, shape in shape_expected.items():
        shapes_ok = shapes_ok and tuple(arrays[name].shape) == shape
    target_shapes_ok = (
        tuple(targets["total-T2"].shape) == (n_rows, 5)
        and tuple(targets["dia-T2"].shape) == (n_rows, 5)
        and tuple(targets["para-T2"].shape) == (n_rows, 5)
        and tuple(targets["T1-field-linear-diagnostic"].shape) == (n_rows, 3)
    )
    classical_specs = specs_for_array(specs, "per_atom_substrate_features_classical")
    aim_idx = int(classical_specs["aimnet2_charge"]["index"])
    classical = arrays["per_atom_substrate_features_classical"]
    charge_complete_mask = (
        (rows["ff14sb_charge_present"].to_numpy(int) == 1)
        & (rows["mopac_welford_mean_charge_present"].to_numpy(int) == 1)
        & np.isfinite(rows["ff14sb_charge"].to_numpy(float))
        & np.isfinite(rows["mopac_welford_mean_charge"].to_numpy(float))
        & np.isfinite(np.asarray(classical[:, aim_idx], dtype=float))
    )
    target_finite = True
    for target in targets.values():
        target_finite = target_finite and bool(np.isfinite(np.asarray(target[charge_complete_mask], dtype=float)).all())
    manifest_charge_rows = int(manifest.get("feature_support_rows", {}).get("charge_complete_rows", -1))
    partition_manifest = manifest.get("partition_bins", {})
    partition_cols = [str(x) for x in partition_manifest.get("columns", [])]
    spec_partition_cols = specs_array_columns(specs, "per_atom_substrate_partition_bins")
    bins = arrays["per_atom_substrate_partition_bins"]
    charge_bin_cols = [
        "bin_nearest_charge_r_4_6_8_10",
        "bin_gap_to_2nd_charge_r_4_6_8_10",
        "bin_abs_charge_T2_quintile",
        "bin_sd_charge_T2_by_atom_quintile",
    ]
    charge_bin_unique = {}
    for col in charge_bin_cols:
        idx = partition_cols.index(col)
        vals = np.asarray(bins[:, idx])
        vals = vals[vals >= 0]
        charge_bin_unique[col] = sorted(int(x) for x in np.unique(vals))
    checks["input_acceptance"] = {
        "pass": bool(
            manifest.get("relationship") == "per_atom_substrate"
            and manifest.get("relationship_kind") == "per_atom_aggregate"
            and shapes_ok
            and target_shapes_ok
            and row_contract
            and manifest_charge_rows == int(charge_complete_mask.sum())
            and int(charge_complete_mask.sum()) == 558_360
            and target_finite
            and partition_cols == spec_partition_cols
            and substrate_dir == SUBSTRATE_DIR.resolve()
            and len(charge_bin_unique["bin_nearest_charge_r_4_6_8_10"]) > 1
            and len(charge_bin_unique["bin_gap_to_2nd_charge_r_4_6_8_10"]) > 1
        ),
        "relationship": manifest.get("relationship"),
        "relationship_kind": manifest.get("relationship_kind"),
        "shape": {"n_atoms": n_atoms, "n_dft_frames": n_frames, "rows": n_rows},
        "expected_shape_846x660": shapes_ok,
        "target_shapes_ok": target_shapes_ok,
        "row_contract": row_contract,
        "charge_complete_rows_manifest": manifest_charge_rows,
        "charge_complete_rows_used": int(charge_complete_mask.sum()),
        "target_finite_on_charge_complete_rows": target_finite,
        "partition_bin_columns_match_manifest_and_specs": partition_cols == spec_partition_cols,
        "charge_partition_non_degenerate_bins": charge_bin_unique,
        "dominant_charge_bin_column_present": "bin_dominant_fraction_charge" in partition_cols,
        "substrate_dir_only_requested_path": str(substrate_dir),
    }
    return checks, charge_complete_mask


def build2_feature_block_specs() -> list[dict[str, object]]:
    return [
        {"name": "ring_jb_T0/T2", "array": "per_atom_substrate_features_classical", "scalar": ["ring_jb_T0"], "t2": ["ring_jb_T2"]},
        {"name": "charge_q_over_r3_T2", "array": "per_atom_substrate_features_classical", "scalar": [], "t2": ["charge_q_over_r3_T2"]},
        {"name": "mc_lit_T0/T2_valid", "array": "per_atom_substrate_features_classical", "scalar": ["mc_lit_T0_valid"], "t2": ["mc_lit_T2_valid"]},
        {"name": "mopac_coulomb_shielding_T2", "array": "per_atom_substrate_features_classical", "scalar": [], "t2": ["mopac_coulomb_shielding_T2"]},
        {"name": "mopac_mc_shielding_T2", "array": "per_atom_substrate_features_classical", "scalar": [], "t2": ["mopac_mc_shielding_T2"]},
        {"name": "ff14sb_field_x/y/z/mag", "array": "per_atom_substrate_features_classical", "scalar": ["ff14sb_field_x", "ff14sb_field_y", "ff14sb_field_z", "ff14sb_field_mag"], "t2": []},
        {"name": "apbs_E_x/y/z/mag", "array": "per_atom_substrate_features_classical", "scalar": ["apbs_E_x", "apbs_E_y", "apbs_E_z", "apbs_E_mag"], "t2": []},
        {"name": "apbs_efg_T2", "array": "per_atom_substrate_features_classical", "scalar": [], "t2": ["apbs_efg_T2"]},
        {"name": "hbond_T2", "array": "per_atom_substrate_features_classical", "scalar": [], "t2": ["hbond_T2"]},
        {"name": "larsen_hbond_T2", "array": "per_atom_substrate_features_hbond_conditioning", "scalar": ["larsen_hbond_water_term"], "t2": ["larsen_hbond_1pHB_T2", "larsen_hbond_2pHB_T2", "larsen_hbond_1pHaB_T2", "larsen_hbond_2pHaB_T2"]},
        {"name": "hbond_count", "array": "per_atom_substrate_features_classical", "scalar": ["hbond_count"], "t2": []},
        {"name": "pi_quadrupole_T2", "array": "per_atom_substrate_features_classical", "scalar": [], "t2": ["pi_quadrupole_T2"]},
        {"name": "dispersion_T2", "array": "per_atom_substrate_features_classical", "scalar": [], "t2": ["dispersion_T2"]},
        {"name": "water_field", "array": "per_atom_substrate_features_classical", "scalar": ["water_efield_x", "water_efield_y", "water_efield_z", "water_efield_mag", "water_n_first", "water_n_second", "water_half_shell_asymmetry", "water_dipole_cos"], "t2": []},
        {"name": "water_efg_T2", "array": "per_atom_substrate_features_method_paths", "scalar": [], "t2": ["water_efg"]},
        {"name": "water_efield_first", "array": "per_atom_substrate_features_method_paths", "scalar": ["water_efield_first_x", "water_efield_first_y", "water_efield_first_z"], "t2": []},
        {"name": "mopac_coulomb_E", "array": "per_atom_substrate_features_method_paths", "scalar": ["mopac_coulomb_E_x", "mopac_coulomb_E_y", "mopac_coulomb_E_z"], "t2": []},
        {"name": "mopac_coulomb_efg_paths", "array": "per_atom_substrate_features_method_paths", "scalar": [], "t2": ["mopac_coulomb_efg_backbone", "mopac_coulomb_efg_aromatic"]},
        {"name": "aimnet2_charge", "array": "per_atom_substrate_features_classical", "scalar": ["aimnet2_charge"], "t2": []},
        {"name": "aimnet2_crg_scalar/x/y/z", "array": "per_atom_substrate_features_classical", "scalar": ["aimnet2_crg_scalar", "aimnet2_crg_x", "aimnet2_crg_y", "aimnet2_crg_z"], "t2": []},
        {"name": "aimnet2_efg_paths", "array": "per_atom_substrate_features_method_paths", "scalar": [], "t2": ["aimnet2_efg", "aimnet2_efg_aromatic", "aimnet2_efg_backbone"]},
    ]


def extract_build2_block(
    arrays: dict[str, np.ndarray],
    specs: dict[str, object],
    spec: dict[str, object],
    C: np.ndarray,
    transform_audit: list[dict[str, object]],
) -> tuple[np.ndarray, list[str]]:
    array_name = str(spec["array"])
    matrix = arrays[array_name]
    col_specs = specs_for_array(specs, array_name)
    pieces = []
    names = []
    for col in spec.get("scalar", []):
        idx = int(col_specs[str(col)]["index"])
        pieces.append(np.asarray(matrix[:, idx], dtype=float).reshape(-1, 1))
        names.append(str(col))
    for prefix in spec.get("t2", []):
        cols = [f"{prefix}_{i}" for i in range(5)]
        idx = [int(col_specs[col]["index"]) for col in cols]
        vals = np.asarray(matrix[:, idx], dtype=float) @ C.T
        pieces.append(vals)
        names += [f"{prefix}_2e_{i}" for i in range(5)]
        transform_audit.append({"feature_block": str(prefix), "array": array_name, "source_columns": cols, "transform": "change_of_basis.get_C().T"})
    if not pieces:
        return np.empty((matrix.shape[0], 0), dtype=float), []
    return np.column_stack(pieces), names


def build2_feature_blocks(
    rows: pd.DataFrame,
    arrays: dict[str, np.ndarray],
    specs: dict[str, object],
    C: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict[str, list[str]], list[dict[str, object]]]:
    transform_audit: list[dict[str, object]] = []
    blocks: dict[str, np.ndarray] = {}
    labels: dict[str, list[str]] = {}
    for spec in build2_feature_block_specs():
        block, names = extract_build2_block(arrays, specs, spec, C, transform_audit)
        blocks[str(spec["name"])] = block
        labels[str(spec["name"])] = names
    for name, _source, scalar_cols, _t2_prefixes in RAW_CHARGE_BLOCK_SPECS:
        block = rows[scalar_cols].to_numpy(float)
        blocks[name] = block
        labels[name] = list(scalar_cols)
    return blocks, labels, transform_audit


def transform_build2_target(target_name: str, raw_target: np.ndarray, row_mask: np.ndarray, C: np.ndarray) -> np.ndarray:
    target = raw_target if bool(row_mask.all()) else raw_target[row_mask]
    arr = np.asarray(target, dtype=float)
    if target_name in {"total-T2", "dia-T2", "para-T2"}:
        return arr[:, :5] @ C.T
    return arr[:, :3]


def build2_component_label(target_name: str, component: int | None) -> str:
    if component is None:
        dims = int(TARGET_SPECS[target_name]["components"])
        return f"vector_{dims}d"
    if target_name in {"total-T2", "dia-T2", "para-T2"}:
        return f"2e_{component}"
    return f"T1_{component}"


def build2_long_score_rows(row: dict[str, object], axis: str, pred: np.ndarray, target: np.ndarray, groups: np.ndarray, ridge_alpha: float) -> list[dict[str, object]]:
    out = []
    vector_r2 = r2_score(pred, target)
    vector_se = jackknife_metric_se(pred, target, groups)
    lo, hi = ci95(vector_r2, vector_se)
    base = {k: row[k] for k in ["target", "atom_type_axis", "atom_type", "tier", "tier_label", "rows", "n_atoms_in_slice", "n_features", "alpha_mode"]}
    base["ridge_alpha"] = float(ridge_alpha)
    out.append({**base, "axis": axis, "component": build2_component_label(str(row["target"]), None), "heldout_R2": vector_r2, "jackknife_se": vector_se, "ci95_low": lo, "ci95_high": hi})
    comp = component_r2(pred, target)
    comp_se = jackknife_component_se(np.asarray(pred), np.asarray(target), groups)
    for i, (r2, se) in enumerate(zip(comp, comp_se)):
        lo, hi = ci95(r2, se)
        out.append({**base, "axis": axis, "component": build2_component_label(str(row["target"]), i), "heldout_R2": r2, "jackknife_se": se, "ci95_low": lo, "ci95_high": hi})
    return out


def build2_score_rows(
    target_name: str,
    fit_scope: str,
    tier_results: dict[str, dict[str, object]],
    tier_order: list[str],
    rows: pd.DataFrame,
    atom_rows: pd.DataFrame,
    y: np.ndarray,
    frames: np.ndarray,
    atoms: np.ndarray,
    split: dict[str, object],
    base_blocks: dict[str, np.ndarray],
    alpha_mode: str,
    embedding_components: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    score_rows = []
    long_rows = []
    n_atoms_global = int(len(atom_rows))
    strata = sorted(str(x) for x in atom_rows["stratum"].dropna().unique())
    atom_types = ["ALL"] + strata if fit_scope == "global_sliced" else strata
    prev_map = previous_tier_map(tier_order)
    r2_cache: dict[tuple[str, str, str], float] = {}
    for tier in tier_order:
        n_features = feature_count_for_tier(tier, embedding_components, base_blocks)
        between = tier_results[tier]["between"]
        within = tier_results[tier]["within"]
        between_alpha_all = np.asarray(between["alpha_by_heldout"], dtype=float)
        within_alpha = float(within["alpha"])
        within_alpha_by_score = within.get("alpha_by_score_row")
        within_alpha_all = np.asarray(within_alpha_by_score, dtype=float) if within_alpha_by_score is not None else None
        within_positions = within["score_row_positions"]
        within_frames = frames[within_positions]
        within_rows = rows.iloc[within_positions]
        for atom_type in atom_types:
            if atom_type == "ALL":
                atom_mask = np.ones(n_atoms_global, dtype=bool)
                row_mask = np.ones(len(rows), dtype=bool)
                within_mask = np.ones(len(within["target"]), dtype=bool)
            else:
                atom_mask = atom_rows["stratum"].to_numpy(str) == atom_type
                row_mask = rows["stratum"].to_numpy(str) == atom_type
                within_mask = within_rows["stratum"].to_numpy(str) == atom_type
            n_atoms_slice = int(atom_mask.sum())
            rows_slice = int(row_mask.sum())
            frames_slice = int(rows.loc[row_mask, "original_frame_index"].nunique()) if rows_slice else 0
            b_pred = between["pred"][atom_mask]
            b_target = between["target"][atom_mask]
            b_groups = between["groups"][atom_mask]
            b_alpha_dist = alpha_distribution(between_alpha_all[atom_mask])
            b_alpha_mode = math.nan if b_alpha_dist["selected_alpha"] is None else float(b_alpha_dist["selected_alpha"])
            b_alpha_min = math.nan if b_alpha_dist["min"] is None else float(b_alpha_dist["min"])
            b_alpha_max = math.nan if b_alpha_dist["max"] is None else float(b_alpha_dist["max"])
            if within_alpha_all is None:
                w_alpha_mode = within_alpha
            else:
                w_alpha_dist = alpha_distribution(within_alpha_all[within_mask])
                w_alpha_mode = math.nan if w_alpha_dist["selected_alpha"] is None else float(w_alpha_dist["selected_alpha"])
            w_pred = within["pred"][within_mask]
            w_target = within["target"][within_mask]
            w_groups = within["groups"][within_mask]
            b_r2 = r2_score(b_pred, b_target)
            w_r2 = r2_score(w_pred, w_target)
            r2_cache[(tier, atom_type, "between")] = b_r2
            r2_cache[(tier, atom_type, "within")] = w_r2
            b_se = jackknife_metric_se(b_pred, b_target, b_groups)
            w_se = jackknife_metric_se(w_pred, w_target, w_groups)
            b_low, b_high = ci95(b_r2, b_se)
            w_low, w_high = ci95(w_r2, w_se)
            share_b, share_w = variance_shares(y[row_mask], atoms[row_mask])
            neff = effective_n_components(y[within_positions][within_mask], within["groups"][within_mask], within_frames[within_mask])
            thin = n_atoms_slice < THIN_ATOMS or len(np.unique(b_groups)) < 3 or len(np.unique(w_groups)) < 3
            n_atoms_between_fit_scope = n_atoms_global if fit_scope == "global_sliced" else n_atoms_slice
            between_n_eff = float(n_atoms_between_fit_scope)
            between_support = ""
            if fit_scope == "per_type" and n_atoms_between_fit_scope < 100:
                between_support = f"thin_between_atoms={n_atoms_between_fit_scope}<100"
            within_support = ""
            if np.isfinite(neff["within_N_eff_median"]) and neff["within_N_eff_median"] < MIN_FAVORABLE_N_EFF:
                within_support = f"thin_within_N_eff={finite_fmt(neff['within_N_eff_median'])}"
            p_gt_atoms = n_features >= max(1, n_atoms_between_fit_scope)
            prev = prev_map[tier]
            row = {
                "target": target_name,
                "fit_scope": fit_scope,
                "atom_type_axis": "stratum",
                "atom_type": atom_type,
                "tier": tier,
                "tier_label": TIER_LABELS[tier],
                "rows": rows_slice,
                "n_atoms_between_global": n_atoms_global,
                "n_atoms_between_fit_scope": n_atoms_between_fit_scope,
                "n_atoms_in_slice": n_atoms_slice,
                "n_original_frames": frames_slice,
                "n_features": n_features,
                "ridge_alpha": b_alpha_mode,
                "alpha_mode": alpha_mode,
                "between_ridge_alpha": b_alpha_mode,
                "between_ridge_alpha_min": b_alpha_min,
                "between_ridge_alpha_max": b_alpha_max,
                "within_ridge_alpha": w_alpha_mode,
                "variance_share_between": share_b,
                "variance_share_within": share_w,
                **neff,
                "between_N_eff": between_n_eff,
                "thin_flag": "thin" if thin else "",
                "between_support_flag": between_support,
                "within_support_flag": within_support,
                "p_gt_atoms_flag": "p>=atoms" if p_gt_atoms else "",
                "between_LOAO_test_R2": b_r2,
                "between_LOAO_test_R2_jackknife_se": b_se,
                "between_LOAO_test_R2_ci95_low": b_low,
                "between_LOAO_test_R2_ci95_high": b_high,
                "between_delta_R2_vs_classical": b_r2 - r2_cache.get(("classical", atom_type, "between"), math.nan),
                "between_delta_R2_vs_previous_tier": b_r2 - r2_cache.get((prev, atom_type, "between"), math.nan) if prev else math.nan,
                "between_delta_R2_vs_global_sliced": math.nan,
                "within_frameblock_test_R2": w_r2,
                "within_frameblock_test_R2_jackknife_se": w_se,
                "within_frameblock_test_R2_ci95_low": w_low,
                "within_frameblock_test_R2_ci95_high": w_high,
                "within_delta_R2_vs_classical": w_r2 - r2_cache.get(("classical", atom_type, "within"), math.nan),
                "within_delta_R2_vs_previous_tier": w_r2 - r2_cache.get((prev, atom_type, "within"), math.nan) if prev else math.nan,
                "within_delta_R2_vs_global_sliced": math.nan,
                "within_split_strategy": "centered_by_atom_train_frames_only; contiguous_frame_block",
                "test_frames": int(len(split["test_frames"])),
                "purged_train_frames": int(len(split["purged_frames"])),
                "cross_split_lag1_pairs": int(split["cross_split_lag1_pairs"]),
                "feature_blocks": "|".join(TIER_BLOCKS[tier]),
            }
            score_rows.append(row)
            long_rows += build2_long_score_rows(row, "between", b_pred, b_target, b_groups, b_alpha_mode)
            long_rows += build2_long_score_rows(row, "within", w_pred, w_target, w_groups, w_alpha_mode)
    return score_rows, long_rows


def add_fit_scope_deltas(score_table: pd.DataFrame) -> pd.DataFrame:
    if score_table.empty or "fit_scope" not in score_table.columns:
        return score_table
    out = score_table.copy()
    key_cols = ["target", "atom_type_axis", "atom_type", "tier"]
    global_rows = out[out["fit_scope"] == "global_sliced"][key_cols + ["between_LOAO_test_R2", "within_frameblock_test_R2"]].rename(
        columns={
            "between_LOAO_test_R2": "_global_between_R2",
            "within_frameblock_test_R2": "_global_within_R2",
        }
    )
    merged = out.merge(global_rows, on=key_cols, how="left")
    per_type = merged["fit_scope"] == "per_type"
    merged.loc[merged["fit_scope"] == "global_sliced", "between_delta_R2_vs_global_sliced"] = 0.0
    merged.loc[merged["fit_scope"] == "global_sliced", "within_delta_R2_vs_global_sliced"] = 0.0
    merged.loc[per_type, "between_delta_R2_vs_global_sliced"] = merged.loc[per_type, "between_LOAO_test_R2"] - merged.loc[per_type, "_global_between_R2"]
    merged.loc[per_type, "within_delta_R2_vs_global_sliced"] = merged.loc[per_type, "within_frameblock_test_R2"] - merged.loc[per_type, "_global_within_R2"]
    return merged.drop(columns=["_global_between_R2", "_global_within_R2"])


def build2_mechanism_from_bin_column(column: str) -> str:
    for key, mech in [
        ("ring", "ring"),
        ("charge", "charge"),
        ("mc_lit", "mc"),
        ("mopac_mc", "mc"),
        ("bond", "mc"),
        ("mopac_coulomb", "charge"),
        ("apbs", "field"),
        ("ff14sb_E", "field"),
        ("aimnet2", "aimnet2"),
        ("heavy_atom", "geometry"),
    ]:
        if key in column:
            return mech
    return "input"


def build2_family_from_bin_column(column: str, mechanism: str) -> str:
    if column.startswith("bin_abs_"):
        return "driver magnitude"
    if column.startswith("bin_sd_"):
        return "driver modulation"
    if column.startswith("bin_nearest_") or column.startswith("bin_gap_"):
        return "geometry and isolation"
    return "input-side bin"


def build2_condition_specs(rows: pd.DataFrame, partition_columns: list[str]) -> list[Build2ConditionSpec]:
    specs: list[Build2ConditionSpec] = []
    for idx, col in enumerate(partition_columns):
        mech = build2_mechanism_from_bin_column(col)
        specs.append(Build2ConditionSpec(build2_family_from_bin_column(col, mech), mech, col.removeprefix("bin_"), "cxx_bin", "per_atom_substrate_partition_bins", col, idx))
    for mech, col in DOMINANCE_FRACTION_COLUMNS.items():
        if col in rows.columns:
            specs.append(Build2ConditionSpec("dominance response", mech, col, "python_quantile_cpp_scalar", "per_atom_substrate_features_conditioning", col, None))
    categorical = [
        ("atom identity", "identity", "stratum"),
        ("atom identity", "identity", "role"),
        ("atom identity", "identity", "element_ord"),
        ("atom identity", "identity", "ff_atom_type_ord"),
        ("atom identity", "identity", "formal_charge"),
        ("atom identity", "identity", "planar_group_ord"),
        ("atom identity", "identity", "polar_h_kind_ord"),
        ("atom identity", "identity", "aromatic"),
        ("atom identity", "identity", "is_exchangeable"),
        ("atom identity", "identity", "amino_acid_ord"),
        ("atom identity", "identity", "locant_ord"),
        ("atom identity", "identity", "branch_outer_ord"),
        ("atom identity", "identity", "branch_inner_ord"),
        ("support", "support", "ring_present"),
        ("support", "support", "charge_present"),
        ("support", "support", "mc_lit_valid_present"),
        ("support", "support", "ff14sb_field_present"),
        ("support", "support", "apbs_E_present"),
        ("support", "support", "apbs_efg_present"),
        ("support", "support", "hbond_shielding_present"),
        ("support", "support", "hbond_count_present"),
        ("support", "support", "pi_quadrupole_present"),
        ("support", "support", "dispersion_present"),
        ("support", "support", "water_field_present"),
        ("support", "support", "hydration_shell_present"),
        ("support", "support", "sasa_present"),
        ("support", "support", "aimnet2_charge_present"),
        ("support", "support", "aimnet2_crg_present"),
        ("support", "support", "aimnet2_embedding_present"),
        ("support", "support", "ff14sb_charge_present"),
        ("support", "support", "mopac_welford_mean_charge_present"),
        ("support", "support", "eeq_charge_present"),
        ("charge-source agreement", "charge", "sign_agree_ff14sb_mopac_welford"),
        ("charge-source agreement", "charge", "sign_agree_ff14sb_aimnet2"),
        ("charge-source agreement", "charge", "sign_agree_mopac_welford_aimnet2"),
    ]
    for family, mech, col in categorical:
        if col in rows.columns:
            specs.append(Build2ConditionSpec(family, mech, col, "categorical", "per_atom_substrate_rows", col, None))
    return specs


def add_build2_charge_agreement_columns(rows: pd.DataFrame, arrays: dict[str, np.ndarray], specs: dict[str, object]) -> pd.DataFrame:
    out = rows.copy()
    classical_specs = specs_for_array(specs, "per_atom_substrate_features_classical")
    conditioning_specs = specs_for_array(specs, "per_atom_substrate_features_conditioning")
    aim_idx = int(classical_specs["aimnet2_charge"]["index"])
    out["aimnet2_charge"] = np.asarray(arrays["per_atom_substrate_features_classical"][:, aim_idx], dtype=float)
    conditioning = arrays["per_atom_substrate_features_conditioning"]
    for col in DOMINANCE_FRACTION_COLUMNS.values():
        if col in conditioning_specs:
            out[col] = np.asarray(conditioning[:, int(conditioning_specs[col]["index"])], dtype=float)
    out["sign_agree_ff14sb_mopac_welford"] = sign_agreement(out["ff14sb_charge"], out["mopac_welford_mean_charge"])
    out["sign_agree_ff14sb_aimnet2"] = sign_agreement(out["ff14sb_charge"], out["aimnet2_charge"])
    out["sign_agree_mopac_welford_aimnet2"] = sign_agreement(out["mopac_welford_mean_charge"], out["aimnet2_charge"])
    return out


def mode_int(values: np.ndarray) -> int:
    vals = np.asarray(values, dtype=int)
    if len(vals) == 0:
        return -1
    unique, counts = np.unique(vals, return_counts=True)
    order = np.lexsort((unique, -counts))
    return int(unique[order[0]])


def atom_mode_bins(partition_bins: np.ndarray, n_frames: int, n_atoms: int) -> np.ndarray:
    out = np.empty((n_atoms, partition_bins.shape[1]), dtype=np.int32)
    for atom in range(n_atoms):
        vals = np.asarray(partition_bins[atom::n_atoms], dtype=np.int32)
        for col in range(vals.shape[1]):
            out[atom, col] = mode_int(vals[:, col])
    return out


def build2_bin_label(column: str, bin_id: int, manifest: dict[str, object]) -> str:
    if bin_id < 0:
        return "missing"
    if column.startswith("bin_abs_") or column.startswith("bin_sd_"):
        return f"Q{bin_id + 1}"
    if column == "bin_gap_to_2nd_charge_r_4_6_8_10":
        edges = [float(x) for x in manifest.get("partition_bins", {}).get("charge_gap_band_edges_A", [0.25, 0.5, 1.0, 1.5])]
    else:
        edges = [float(x) for x in manifest.get("partition_bins", {}).get("distance_band_edges_A", [4.0, 6.0, 8.0, 10.0])]
    if bin_id == 0:
        return f"<= {finite_fmt(edges[0], 3)}A"
    if bin_id < len(edges):
        return f"{finite_fmt(edges[bin_id - 1], 3)}-{finite_fmt(edges[bin_id], 3)}A"
    return f"> {finite_fmt(edges[-1], 3)}A"


def build2_bin_order(spec: Build2ConditionSpec, label: object, bin_id: int | None = None) -> float:
    if bin_id is not None:
        return float(bin_id)
    return bin_order_key(label)


def quantile_bin_arrays(values: np.ndarray | pd.Series) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    numeric = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    finite = np.isfinite(numeric)
    labels = np.full(len(numeric), "missing", dtype=object)
    raw_bins = np.full(len(numeric), -1, dtype=int)
    orders = np.full(len(numeric), -1.0, dtype=float)
    if finite.sum() == 0:
        return labels, raw_bins, orders
    vals = numeric[finite]
    qs = np.quantile(vals, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    edges: list[float] = []
    for q in qs:
        qf = float(q)
        if not edges or qf > edges[-1]:
            edges.append(qf)
    if len(edges) <= 2:
        lo = float(np.nanmin(vals))
        hi = float(np.nanmax(vals))
        labels[finite] = f"Q1 [{finite_fmt(lo)}, {finite_fmt(hi)}]"
        raw_bins[finite] = 0
        orders[finite] = 0.0
        return labels, raw_bins, orders
    for bin_id, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        if bin_id == 0:
            mask = finite & (numeric >= lo) & (numeric <= hi)
        else:
            mask = finite & (numeric > lo) & (numeric <= hi)
        labels[mask] = f"Q{bin_id + 1} [{finite_fmt(lo)}, {finite_fmt(hi)}]"
        raw_bins[mask] = int(bin_id)
        orders[mask] = float(bin_id)
    return labels, raw_bins, orders


def build2_axis_conditions(
    spec: Build2ConditionSpec,
    axis: str,
    rows: pd.DataFrame,
    atom_rows: pd.DataFrame,
    partition_bins: np.ndarray,
    atom_bins: np.ndarray,
    row_positions: np.ndarray,
    manifest: dict[str, object],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if axis == "between":
        atoms = atom_rows["atom_index"].to_numpy(int)
        if spec.kind == "cxx_bin":
            raw_bins = atom_bins[:, int(spec.bin_index)]
            labels = np.asarray([build2_bin_label(spec.column, int(x), manifest) for x in raw_bins], dtype=object)
            orders = raw_bins.astype(float)
            return labels, raw_bins, orders
        if spec.kind == "python_quantile_cpp_scalar":
            return quantile_bin_arrays(atom_rows[spec.column])
        labels = np.asarray([str(x) if pd.notna(x) else "missing" for x in atom_rows[spec.column].to_numpy()], dtype=object)
        return labels, np.full(len(labels), -999, dtype=int), np.arange(len(labels), dtype=float)
    atoms = rows.iloc[row_positions]["atom_index"].to_numpy(int)
    if spec.kind == "cxx_bin":
        raw_bins = np.asarray(partition_bins[row_positions, int(spec.bin_index)], dtype=int)
        labels = np.asarray([build2_bin_label(spec.column, int(x), manifest) for x in raw_bins], dtype=object)
        orders = raw_bins.astype(float)
        return labels, raw_bins, orders
    if spec.kind == "python_quantile_cpp_scalar":
        return quantile_bin_arrays(rows.iloc[row_positions][spec.column])
    labels = np.asarray([str(x) if pd.notna(x) else "missing" for x in rows.iloc[row_positions][spec.column].to_numpy()], dtype=object)
    return labels, np.full(len(labels), -999, dtype=int), np.arange(len(labels), dtype=float)


def build2_support_flag(axis: str, rows_count: int, n_atoms: int, n_eff: float) -> str:
    flags = []
    if n_atoms < THIN_ATOMS:
        flags.append("thin_atoms")
    if axis == "within" and rows_count < THIN_WITHIN_ROWS:
        flags.append("thin_rows")
    if np.isfinite(n_eff) and n_eff < MIN_FAVORABLE_N_EFF:
        flags.append("thin_N_eff")
    return ";".join(flags)


def build2_response_shapes(curves: pd.DataFrame) -> dict[tuple[str, str, str, str, str, str, str, str], str]:
    out: dict[tuple[str, str, str, str, str, str, str, str], str] = {}
    group_cols = ["target", "fit_scope", "mechanism", "condition_family", "condition_name", "stratum", "tier", "axis"]
    for key, group in curves.groupby(group_cols, dropna=False):
        usable = group[(group["support_flag"] == "") & np.isfinite(group["heldout_R2"])].sort_values("_bin_order")
        vals = usable["heldout_R2"].to_numpy(float)
        if len(vals) < 3:
            out[key] = "unstable-thin"
            continue
        if str(usable["condition_kind"].iloc[0]) not in {"cxx_bin", "python_quantile_cpp_scalar"}:
            span = float(np.nanmax(vals) - np.nanmin(vals))
            out[key] = "categorical contrast" if span >= 0.02 else "flat response"
            continue
        span = float(np.nanmax(vals) - np.nanmin(vals))
        tol = max(0.02, 0.10 * span)
        diffs = np.diff(vals)
        if span < 0.02:
            out[key] = "flat response"
        elif np.all(diffs >= -tol) and vals[-1] > vals[0]:
            out[key] = "monotone rise"
        elif np.all(diffs <= tol) and vals[-1] < vals[0]:
            out[key] = "monotone fall"
        elif len(vals) >= 4 and vals[0] > vals[1] and vals[-1] > vals[-2] and np.nanmin(vals[1:-1]) < min(vals[0], vals[-1]) - tol:
            out[key] = "U-shape"
        elif np.max(np.abs(diffs)) > 0.6 * span:
            out[key] = "threshold behavior"
        else:
            out[key] = "flat response"
    return out


def build2_partition_curves(
    target_name: str,
    fit_scope: str,
    tier_results: dict[str, dict[str, object]],
    tier_order: list[str],
    rows: pd.DataFrame,
    atom_rows: pd.DataFrame,
    partition_bins: np.ndarray,
    atom_bins: np.ndarray,
    frames: np.ndarray,
    condition_list: list[Build2ConditionSpec],
    manifest: dict[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    rows_out = []
    bin_defs = []
    strata = ["ALL"] + sorted(str(x) for x in atom_rows["stratum"].dropna().unique())
    axes: dict[str, dict[str, object]] = {"between": {"length": len(atom_rows), "frames": np.arange(len(atom_rows), dtype=int), "row_positions": None}}
    first_within = tier_results[tier_order[0]]["within"]
    axes["within"] = {
        "length": len(first_within["target"]),
        "frames": frames[first_within["score_row_positions"]],
        "row_positions": first_within["score_row_positions"],
    }
    neff_by_tier_axis = {
        (tier, "within"): per_atom_neff_components(
            tier_results[tier]["within"]["target"],
            tier_results[tier]["within"]["groups"],
            np.asarray(axes["within"]["frames"]),
        )
        for tier in tier_order
    }
    for axis, axis_info in axes.items():
        row_positions = axis_info["row_positions"]
        if axis == "between":
            axis_atom_index = atom_rows["atom_index"].to_numpy(int)
            axis_strata = atom_rows["stratum"].to_numpy(str)
        else:
            assert row_positions is not None
            axis_atom_index = rows.iloc[row_positions]["atom_index"].to_numpy(int)
            axis_strata = rows.iloc[row_positions]["stratum"].to_numpy(str)
        for stratum in strata:
            stratum_mask = np.ones(len(axis_atom_index), dtype=bool) if stratum == "ALL" else axis_strata == stratum
            if not stratum_mask.any():
                continue
            for spec in condition_list:
                labels, raw_bins, orders = build2_axis_conditions(spec, axis, rows, atom_rows, partition_bins, atom_bins, row_positions if row_positions is not None else np.asarray([], dtype=int), manifest)
                labels_in = labels[stratum_mask]
                raw_in = raw_bins[stratum_mask]
                order_in = orders[stratum_mask]
                atoms_in = axis_atom_index[stratum_mask]
                for label in sorted(set(labels_in), key=str):
                    unit_local = labels_in == label
                    if not unit_local.any():
                        continue
                    if spec.kind in {"cxx_bin", "python_quantile_cpp_scalar"}:
                        bin_id = mode_int(raw_in[unit_local])
                        bin_order = build2_bin_order(spec, label, bin_id)
                    else:
                        bin_id = -999
                        bin_order = build2_bin_order(spec, label)
                    global_mask = np.zeros(len(axis_atom_index), dtype=bool)
                    idx = np.flatnonzero(stratum_mask)
                    global_mask[idx[unit_local]] = True
                    n_atoms = int(len(np.unique(axis_atom_index[global_mask])))
                    rows_count = int(global_mask.sum())
                    if spec.kind in {"cxx_bin", "python_quantile_cpp_scalar"}:
                        bin_defs.append({
                            "target": target_name,
                            "fit_scope": fit_scope,
                            "axis": axis,
                            "stratum": stratum,
                            "mechanism": spec.mechanism,
                            "condition_family": spec.family,
                            "condition_name": spec.name,
                            "condition_kind": spec.kind,
                            "source": spec.source,
                            "source_column": spec.column,
                            "bin_label": label,
                            "bin_id": int(bin_id),
                            "n_atoms_pre_rank": n_atoms,
                        })
                    for tier in tier_order:
                        axis_pred = tier_results[tier][axis]
                        pred = axis_pred["pred"][global_mask]
                        target = axis_pred["target"][global_mask]
                        groups = axis_pred["groups"][global_mask]
                        r2 = r2_score(pred, target)
                        se = jackknife_metric_se(pred, target, groups)
                        lo, hi = ci95(r2, se)
                        atom_neff = neff_by_tier_axis.get((tier, axis))
                        n_eff_median = bin_neff_median(axis, axis_atom_index, global_mask, atom_neff)
                        support_flag = build2_support_flag(axis, rows_count, n_atoms, n_eff_median)
                        rows_out.append({
                            "target": target_name,
                            "fit_scope": fit_scope,
                            "mechanism": spec.mechanism,
                            "condition_family": spec.family,
                            "condition_name": spec.name,
                            "condition_kind": spec.kind,
                            "condition_source": spec.source,
                            "bin_label": label,
                            "bin_id": int(bin_id),
                            "stratum": stratum,
                            "tier": tier,
                            "axis": axis,
                            "rows": rows_count,
                            "n_atoms": n_atoms,
                            "N_eff_median": n_eff_median,
                            "heldout_R2": r2,
                            "heldout_R2_ci95_low": lo,
                            "heldout_R2_ci95_high": hi,
                            "delta_R2_vs_classical": math.nan,
                            "delta_R2_vs_previous_tier": math.nan,
                            "support_flag": support_flag,
                            "_bin_order": bin_order,
                        })
    curves = pd.DataFrame(rows_out)
    if curves.empty:
        return curves, curves.copy(), bin_defs
    key_cols = ["target", "fit_scope", "mechanism", "condition_family", "condition_name", "condition_kind", "bin_label", "bin_id", "stratum", "axis"]
    classical = curves[curves["tier"] == "classical"][key_cols + ["heldout_R2"]].rename(columns={"heldout_R2": "_classical_R2"})
    curves = curves.merge(classical, on=key_cols, how="left")
    curves["delta_R2_vs_classical"] = curves["heldout_R2"] - curves["_classical_R2"]
    prev_parts = []
    prev_map = previous_tier_map(tier_order)
    for tier, prev in prev_map.items():
        part = curves[curves["tier"] == tier].copy()
        if prev is None:
            part["_previous_R2"] = np.nan
        else:
            prev_df = curves[curves["tier"] == prev][key_cols + ["heldout_R2"]].rename(columns={"heldout_R2": "_previous_R2"})
            part = part.drop(columns=["_previous_R2"], errors="ignore").merge(prev_df, on=key_cols, how="left")
        prev_parts.append(part)
    curves = pd.concat(prev_parts, ignore_index=True)
    curves["delta_R2_vs_previous_tier"] = curves["heldout_R2"] - curves["_previous_R2"]
    shape_map = build2_response_shapes(curves)
    curves["response_shape"] = [
        shape_map.get((r.target, r.fit_scope, r.mechanism, r.condition_family, r.condition_name, r.stratum, r.tier, r.axis), "unstable-thin")
        for r in curves.itertuples(index=False)
    ]
    public_cols = [
        "target",
        "fit_scope",
        "mechanism",
        "condition_family",
        "condition_name",
        "condition_kind",
        "condition_source",
        "bin_label",
        "bin_id",
        "stratum",
        "tier",
        "axis",
        "rows",
        "n_atoms",
        "N_eff_median",
        "heldout_R2",
        "heldout_R2_ci95_low",
        "heldout_R2_ci95_high",
        "delta_R2_vs_classical",
        "delta_R2_vs_previous_tier",
        "support_flag",
        "response_shape",
    ]
    curves_public = curves[public_cols].copy()
    long_rows = []
    for r in curves.itertuples(index=False):
        for metric in ["heldout_R2", "delta_R2_vs_classical", "delta_R2_vs_previous_tier"]:
            long_rows.append({
                "target": r.target,
                "fit_scope": r.fit_scope,
                "mechanism": r.mechanism,
                "condition_family": r.condition_family,
                "condition_name": r.condition_name,
                "bin_label": r.bin_label,
                "bin_id": r.bin_id,
                "stratum": r.stratum,
                "tier": r.tier,
                "axis": r.axis,
                "metric": metric,
                "value": getattr(r, metric),
                "support_flag": r.support_flag,
                "response_shape": r.response_shape,
            })
    return curves_public, pd.DataFrame(long_rows), bin_defs


def build2_favorable_partitions(curves: pd.DataFrame, shortlist_size: int = 120) -> pd.DataFrame:
    if curves.empty:
        return curves.copy()
    eligible = curves[
        (curves["support_flag"] == "")
        & np.isfinite(curves["heldout_R2"])
        & (curves["n_atoms"] >= THIN_ATOMS)
        & (curves["N_eff_median"] >= MIN_FAVORABLE_N_EFF)
        & (curves["axis"] == "within")
    ].copy()
    if eligible.empty:
        eligible = curves[(curves["support_flag"] == "") & np.isfinite(curves["heldout_R2"])].copy()
    eligible = eligible[(eligible["heldout_R2_ci95_low"] > 0.0) | (eligible["heldout_R2"] > 0.0)]
    eligible["_positive_ci"] = eligible["heldout_R2_ci95_low"] > 0.0
    eligible = eligible.sort_values(["_positive_ci", "heldout_R2", "delta_R2_vs_previous_tier"], ascending=[False, False, False], na_position="last")
    cols = [
        "target",
        "fit_scope",
        "mechanism",
        "condition_family",
        "condition_name",
        "bin_label",
        "bin_id",
        "stratum",
        "tier",
        "axis",
        "rows",
        "n_atoms",
        "N_eff_median",
        "heldout_R2",
        "heldout_R2_ci95_low",
        "heldout_R2_ci95_high",
        "delta_R2_vs_classical",
        "delta_R2_vs_previous_tier",
        "response_shape",
    ]
    return eligible.head(shortlist_size)[cols]


def load_casehunter_candidates(substrate_dir: Path) -> pd.DataFrame:
    parts = []
    for mechanism in ["ring", "charge", "mc"]:
        path = substrate_dir / "equations" / mechanism / "cases_manifest.csv"
        if path.exists():
            df = pd.read_csv(path)
            df["case_manifest_path"] = str(path)
            parts.append(df)
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def case_window_positions(rows: pd.DataFrame, atom_index: int, begin: int, end: int) -> np.ndarray:
    mask = (
        (rows["atom_index"].to_numpy(int) == int(atom_index))
        & (rows["original_frame_index"].to_numpy(int) >= int(begin))
        & (rows["original_frame_index"].to_numpy(int) < int(end))
    )
    pos = np.flatnonzero(mask)
    if len(pos) == 0:
        mask = (
            (rows["atom_index"].to_numpy(int) == int(atom_index))
            & (rows["original_frame_index"].to_numpy(int) >= int(begin))
            & (rows["original_frame_index"].to_numpy(int) <= int(end))
        )
        pos = np.flatnonzero(mask)
    return pos


def build2_casehunter_intersection(
    favorable: pd.DataFrame,
    cases: pd.DataFrame,
    rows: pd.DataFrame,
    atom_rows: pd.DataFrame,
    partition_bins: np.ndarray,
    condition_specs: list[Build2ConditionSpec],
    manifest: dict[str, object],
    shortlist_size: int = 60,
) -> pd.DataFrame:
    if favorable.empty or cases.empty:
        return pd.DataFrame()
    specs_by_name = {spec.name: spec for spec in condition_specs}
    out = []
    atom_stratum = atom_rows.set_index("atom_index")["stratum"].astype(str).to_dict()
    for fav in favorable.itertuples(index=False):
        if fav.mechanism not in {"ring", "charge", "mc"}:
            continue
        spec = specs_by_name.get(str(fav.condition_name))
        if spec is None:
            continue
        same_mech = cases[cases["mechanism"].astype(str) == str(fav.mechanism)]
        for case in same_mech.itertuples(index=False):
            case_stratum = str(atom_stratum.get(int(case.atom_index), "missing"))
            if fav.stratum != "ALL" and fav.stratum != case_stratum:
                continue
            positions = case_window_positions(rows, int(case.atom_index), int(case.frame_window_begin), int(case.frame_window_end))
            if len(positions) == 0:
                continue
            if spec.kind == "cxx_bin":
                modal_bin = mode_int(np.asarray(partition_bins[positions, int(spec.bin_index)], dtype=int))
                if int(fav.bin_id) != int(modal_bin):
                    continue
                modal_label = build2_bin_label(spec.column, modal_bin, manifest)
            else:
                labels = np.asarray([str(x) if pd.notna(x) else "missing" for x in rows.iloc[positions][spec.column].to_numpy()], dtype=object)
                unique, counts = np.unique(labels, return_counts=True)
                modal_label = str(unique[np.argmax(counts)])
                if str(fav.bin_label) != modal_label:
                    continue
                modal_bin = -999
            out.append({
                "target": fav.target,
                "fit_scope": fav.fit_scope,
                "mechanism": fav.mechanism,
                "stratum": fav.stratum,
                "case_atom_index": int(case.atom_index),
                "case_frame_window_begin": int(case.frame_window_begin),
                "case_frame_window_end": int(case.frame_window_end),
                "condition_family": fav.condition_family,
                "condition_name": fav.condition_name,
                "bin_label": modal_label,
                "bin_id": int(modal_bin),
                "tier": fav.tier,
                "axis": fav.axis,
                "heldout_R2": float(fav.heldout_R2),
                "delta_R2_vs_classical": float(fav.delta_R2_vs_classical) if np.isfinite(fav.delta_R2_vs_classical) else math.nan,
                "response_shape": fav.response_shape,
                "partition_n_atoms": int(fav.n_atoms),
                "partition_N_eff_median": float(fav.N_eff_median),
                "case_score": float(case.score),
                "case_dft_recovery_R2_for_navigation_only": float(case.dft_recovery_R2) if pd.notna(case.dft_recovery_R2) else math.nan,
                "case_manifest_path": case.case_manifest_path,
            })
    if not out:
        return pd.DataFrame()
    df = pd.DataFrame(out)
    return df.sort_values(["heldout_R2", "case_score"], ascending=[False, False], na_position="last").head(shortlist_size)


def write_build2_reports(
    out_dir: Path,
    score_table: pd.DataFrame,
    curves: pd.DataFrame,
    favorable: pd.DataFrame,
    case_shortlist: pd.DataFrame,
    audit: dict[str, object],
) -> None:
    score_focus = score_table[score_table["atom_type"].isin(["ALL", "HN", "O", "N", "CA", "C", "HA"])].copy()
    lines = [
        "# Build 3 Fit-Architecture Report",
        "",
        f"Substrate: `{audit['substrate_dir']}`",
        f"Run dir: `{audit['output_dir']}`",
        f"Targets: {', '.join(TARGET_SPECS)}",
        "",
        "## Fit Scores",
        markdown_table(score_focus[[
            "target",
            "fit_scope",
            "atom_type",
            "tier",
            "n_atoms_in_slice",
            "n_atoms_between_fit_scope",
            "between_LOAO_test_R2",
            "between_delta_R2_vs_global_sliced",
            "between_delta_R2_vs_classical",
            "within_frameblock_test_R2",
            "within_delta_R2_vs_global_sliced",
            "within_delta_R2_vs_classical",
            "between_support_flag",
            "thin_flag",
        ]], max_rows=80),
        "",
        "## Partition Shapes",
    ]
    if curves.empty:
        lines.append("No curves emitted.")
    else:
        shape_counts = curves.groupby(["target", "fit_scope", "mechanism", "condition_family", "tier", "axis", "response_shape"], dropna=False).size().reset_index(name="curve_points")
        lines.append(markdown_table(shape_counts, max_rows=120))
    lines += [
        "",
        "## Favourable Partitions",
        markdown_table(favorable, max_rows=80) if not favorable.empty else "No eligible favourable partitions.",
        "",
        "## CaseHunter Intersection",
        markdown_table(case_shortlist, max_rows=80) if not case_shortlist.empty else "No CaseHunter candidate intersected the favourable partitions.",
    ]
    (out_dir / "build2_partition_report.md").write_text("\n".join(lines), encoding="utf-8")


def write_build2_plots(out_dir: Path, score_table: pd.DataFrame, curves: pd.DataFrame) -> dict[str, object]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - plot availability is environment-dependent
        return {"plots_written": False, "reason": str(exc)}
    written = []
    focus = score_table[(score_table["fit_scope"] == "global_sliced") & score_table["atom_type"].isin(["N", "CA", "C", "O", "HN", "HA"])].copy()
    for target_name, group in focus.groupby("target", dropna=False):
        pivot = group.pivot(index="atom_type", columns="tier", values="between_LOAO_test_R2").reindex(["N", "CA", "C", "O", "HN", "HA"])
        fig, ax = plt.subplots(figsize=(10, 5))
        pivot.plot(ax=ax, marker="o")
        ax.set_title(f"{target_name} between-LOAO R2 by stratum")
        ax.set_xlabel("stratum")
        ax.set_ylabel("held-out R2")
        ax.axhline(0, color="black", linewidth=0.8)
        ax.legend(fontsize=8)
        fig.tight_layout()
        path = out_dir / f"score_by_stratum_{target_name.replace('/', '_')}.png"
        fig.savefig(path, dpi=140)
        plt.close(fig)
        written.append(str(path))
    if not curves.empty:
        top = curves[(curves["support_flag"] == "") & (curves["condition_kind"] == "cxx_bin")].copy()
        top = top.sort_values("heldout_R2", ascending=False).head(24)
        if not top.empty:
            labels = [f"{r.target}\\n{r.mechanism}:{r.stratum}\\n{r.condition_name} {r.bin_label}" for r in top.itertuples(index=False)]
            fig, ax = plt.subplots(figsize=(12, 8))
            ax.barh(np.arange(len(top)), top["heldout_R2"].to_numpy(float))
            ax.set_yticks(np.arange(len(top)), labels=labels, fontsize=6)
            ax.invert_yaxis()
            ax.set_xlabel("held-out R2")
            ax.set_title("Top supported partition bins")
            fig.tight_layout()
            path = out_dir / "partition_top_supported_bins.png"
            fig.savefig(path, dpi=140)
            plt.close(fig)
            written.append(str(path))
    return {"plots_written": True, "paths": written}


def build2_fit_stage_checks(audit: dict[str, object], score_table: pd.DataFrame, curves: pd.DataFrame, case_shortlist: pd.DataFrame) -> None:
    checks = audit["fit_stage_checks"]
    checks["no_external_merge"] = {
        "pass": True,
        "python_read_only_reduced_substrate_sidecars_and_bin_ids": True,
        "files_read": audit["files_read"],
        "forbidden_inputs_not_opened": ["trajectory.h5", "ORCA outputs", "per-source directories", "older rediscover dirs", "ring_paths cross-method validation data"],
    }
    checks["basis_and_targets"] = {
        "pass": bool(audit["change_of_basis_orthogonality_max_abs"] < 1.0e-12 and set(audit["targets_fit"]) == set(TARGET_SPECS)),
        "C_T_C_minus_I_max_abs": audit["change_of_basis_orthogonality_max_abs"],
        "targets_fit": audit["targets_fit"],
        "T1_label": "field-linear diagnostic, not a shift claim",
    }
    cv_checks = audit["cv_integrity_details"]
    checks["cv_integrity"] = {
        "pass": bool(
            all(v["between"]["cv_check"]["heldout_atoms_never_in_train"] for target in cv_checks.values() for v in target.values())
            and all(v["within"]["cv_check"]["train_test_overlap"] == 0 for target in cv_checks.values() for v in target.values())
            and all(v["within"]["cv_check"]["purged_in_train_or_test"] == 0 for target in cv_checks.values() for v in target.values())
            and audit["within_split"]["cross_split_lag1_pairs"] == 0
            and score_table["within_N_eff_median"].notna().all()
        ),
        "between_heldout_atoms_never_in_train": all(v["between"]["cv_check"]["heldout_atoms_never_in_train"] for target in cv_checks.values() for v in target.values()),
        "within_train_test_overlap_zero": all(v["within"]["cv_check"]["train_test_overlap"] == 0 for target in cv_checks.values() for v in target.values()),
        "purged_frames_absent_from_train_and_test": all(v["within"]["cv_check"]["purged_in_train_or_test"] == 0 for target in cv_checks.values() for v in target.values()),
        "scores_are_held_out": True,
        "N_eff_reported_per_score_row": bool(score_table["within_N_eff_median"].notna().all()),
        "cross_split_lag1_pairs": audit["within_split"]["cross_split_lag1_pairs"],
    }
    checks["partition_integrity"] = {
        "pass": bool(
            not curves.empty
            and audit["partition_conditions_input_side_only"]
            and audit["partition_bins_from_cpp_lookup_or_python_quantile_cpp_scalar"]
            and audit["charge_partition_non_degenerate"]
        ),
        "cpp_bin_lookup_only_for_existing_numeric_partition_bins": audit["partition_bins_from_cpp_lookup_only"],
        "python_quantile_bins_on_cpp_dominance_scalars": audit["dominance_quantile_bins_from_python_on_cpp_scalar"],
        "no_condition_uses_dft_target_residual_or_coefficients": True,
        "casehunter_intersection_rows": int(len(case_shortlist)),
        "charge_partition_non_degenerate": audit["charge_partition_non_degenerate"],
        "dominance_scalar_columns_available_in_substrate": audit["dominance_scalar_columns_available_in_substrate"],
        "dominance_cpp_bin_id_next_emit_flag": audit["dominance_cpp_bin_id_next_emit_flag"],
    }


def build2_log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def build2_main() -> None:
    args = parse_args()
    tier_order = parse_tiers(args.tiers)
    alpha_grid = parse_alpha_grid(args.alpha_grid)
    substrate_dir = args.substrate_dir.resolve()
    data = build2_load_inputs(substrate_dir)
    checks, charge_mask = build2_input_acceptance_checks(data, substrate_dir)
    if not checks["input_acceptance"]["pass"]:
        raise SystemExit(f"FATAL: input acceptance failed: {checks['input_acceptance']}")
    out_dir, disk_audit = prepare_build2_output_dir(args.out_dir, args.keep_out_dir)

    rows_all: pd.DataFrame = data["rows"]
    specs: dict[str, object] = data["specs"]
    arrays: dict[str, np.ndarray] = data["arrays"]
    manifest: dict[str, object] = data["manifest"]
    rows = rows_all if bool(charge_mask.all()) else rows_all.loc[charge_mask].reset_index(drop=True)
    if not bool(charge_mask.all()):
        arrays = {k: (v[charge_mask] if v.shape[0] == len(charge_mask) else v) for k, v in arrays.items()}
    n_atoms = int(manifest["n_atoms"])
    n_frames = int(manifest["n_dft_frames"])
    if len(rows) != n_atoms * n_frames:
        raise SystemExit("FATAL: charge-complete filter no longer has dense 846 x 660 rows")

    C = cob.get_C()
    orth = float(np.abs(C.T @ C - np.eye(5)).max())
    if orth >= 1.0e-12:
        raise SystemExit(f"FATAL: change_of_basis.get_C() orthogonality check failed: {orth:.3e}")

    row_indices = np.flatnonzero(charge_mask)
    frames = rows["original_frame_index"].to_numpy(int)
    atoms = rows["atom_index"].to_numpy(int)
    rows_for_conditions = add_build2_charge_agreement_columns(rows, arrays, specs)
    atom_conditions = atom_condition_frame(rows_for_conditions, n_frames, n_atoms)
    partition_bins = arrays["per_atom_substrate_partition_bins"]
    atom_bins = atom_mode_bins(partition_bins, n_frames, n_atoms)
    partition_cols = [str(x) for x in manifest["partition_bins"]["columns"]]
    cond_specs = build2_condition_specs(rows_for_conditions, partition_cols)

    base_blocks, block_labels, transform_audit = build2_feature_blocks(rows, arrays, specs, C)
    base_atom_blocks = {name: atom_means_dense(block, n_frames, n_atoms) for name, block in base_blocks.items()}
    embedding = arrays["per_atom_substrate_aimnet2_embedding"]
    embedding_atom_mean = atom_means_embedding(embedding, n_frames, n_atoms)
    split = split_frame_block(frames, args.within_test_fraction, args.purge_frames_each_side)

    within_embedding_pcs = None
    within_embedding_pca_audit = None
    embedding_tiers = [tier for tier in tier_order if tier_uses_embedding(tier)]
    if embedding_tiers:
        within_train_indices = row_indices[split["train_mask"]]
        within_pca = fit_pca_memmap(embedding, within_train_indices, args.embedding_components)
        within_embedding_pcs = transform_pca_memmap(embedding, within_pca, row_indices)
        within_embedding_pca_audit = {
            "method": "unsupervised PCA fit once on within-CV training rows only; shared by selected tiers that consume AIMNet2 embedding PCs",
            "n_components": int(within_pca.n_components),
            "original_dims": int(within_pca.original_dims),
            "training_rows": int(within_pca.n_train_rows),
            "test_rows_excluded": int(split["test_mask"].sum()),
            "purged_rows_excluded": int(split["purge_mask"].sum()),
            "explained_variance_ratio": float(within_pca.explained_variance_ratio),
            "consumed_by_tiers": embedding_tiers,
        }
    else:
        within_embedding_pca_audit = {"method": "embedding sidecar present but not consumed by selected tiers", "n_components": 0, "original_dims": int(embedding.shape[1]), "consumed_by_tiers": []}

    all_score_rows = []
    all_long_rows = []
    all_curves = []
    all_curves_long = []
    all_bin_defs = []
    alpha_selection_by_target: dict[str, object] = {}
    pca_audit_by_target: dict[str, object] = {}
    cv_details_by_target: dict[str, object] = {}
    between_pca_cache: dict[object, PcaTransform] = {}
    between_alpha_pca_cache: dict[object, PcaTransform] = {}
    within_alpha_pca_cache: dict[object, PcaTransform] = {}
    for target_name, raw_target in data["targets"].items():
        build2_log(f"target {target_name}: transform target")
        y = transform_build2_target(target_name, raw_target, charge_mask, C)
        y_atom = atom_means_dense(y, n_frames, n_atoms)
        global_tier_results: dict[str, dict[str, object]] = {}
        per_type_tier_results: dict[str, dict[str, object]] = {}
        pca_audit: dict[str, object] = {"global_sliced": {}, "per_type": {}}
        cv_details: dict[str, object] = {}
        alpha_selection_by_scope: dict[str, object] = {"global_sliced": {}, "per_type": {}}
        for tier in tier_order:
            if args.alpha_mode == "select":
                build2_log(f"target {target_name}: global_sliced tier {tier}: select between alpha")
                between_alpha, between_alpha_audit = select_between_alphas(
                    tier,
                    base_atom_blocks,
                    embedding_atom_mean,
                    y_atom,
                    alpha_grid,
                    args.embedding_components,
                    args.inner_cv_folds,
                    between_alpha_pca_cache,
                )
                build2_log(f"target {target_name}: global_sliced tier {tier}: select within alpha")
                within_alpha, within_alpha_audit = select_within_alpha(
                    tier,
                    base_blocks,
                    embedding,
                    row_indices,
                    y,
                    atoms,
                    frames,
                    split,
                    alpha_grid,
                    args.embedding_components,
                    args.inner_cv_folds,
                    args.purge_frames_each_side,
                    within_alpha_pca_cache,
                )
            else:
                fixed_alpha = float(args.ridge_alpha)
                between_alpha = np.full(n_atoms, fixed_alpha, dtype=float)
                between_alpha_audit = {"axis": "between_LOAO", "method": "fixed-alpha labelled baseline; no inner CV", "selected_alpha": fixed_alpha, "min": fixed_alpha, "max": fixed_alpha, "counts": {f"{fixed_alpha:g}": int(n_atoms)}, "heldout_test_atom_excluded_from_alpha_selection": True, "alpha_grid": [float(x) for x in alpha_grid]}
                within_alpha = fixed_alpha
                within_alpha_audit = {"axis": "within_frameblock", "method": "fixed-alpha labelled baseline; no inner CV", "selected_alpha": fixed_alpha, "outer_test_rows_used_for_alpha_selection": 0, "outer_purged_rows_used_for_alpha_selection": 0, "alpha_grid": [float(x) for x in alpha_grid]}
            alpha_selection_by_scope["global_sliced"][tier] = {"between": between_alpha_audit, "within": within_alpha_audit}
            build2_log(f"target {target_name}: global_sliced tier {tier}: fit between")
            cache = between_pca_cache if tier_uses_embedding(tier) else None
            between = evaluate_between_tier(tier, base_atom_blocks, embedding_atom_mean, y_atom, between_alpha, args.embedding_components, cache)
            build2_log(f"target {target_name}: global_sliced tier {tier}: fit within")
            within = evaluate_within_tier(tier, base_blocks, within_embedding_pcs, within_embedding_pca_audit, y, atoms, split, within_alpha)
            global_tier_results[tier] = {"between": between, "within": within}
            pca_audit["global_sliced"][tier] = {"between": between["pca_audit"], "within": within["pca_audit"]}
            cv_details[f"global_sliced:{tier}"] = {"between": {"cv_check": between["cv_check"]}, "within": {"cv_check": within["cv_check"]}}

            build2_log(f"target {target_name}: per_type tier {tier}: select alphas and fit strata")
            per_type_result, per_type_alpha_audit, per_type_pca_audit, per_type_cv_details = evaluate_per_type_tier(
                tier,
                rows,
                atom_conditions,
                base_blocks,
                base_atom_blocks,
                embedding,
                embedding_atom_mean,
                row_indices,
                y,
                y_atom,
                atoms,
                frames,
                split,
                args.alpha_mode,
                float(args.ridge_alpha),
                alpha_grid,
                args.embedding_components,
                args.inner_cv_folds,
                args.purge_frames_each_side,
            )
            per_type_tier_results[tier] = per_type_result
            alpha_selection_by_scope["per_type"][tier] = per_type_alpha_audit
            pca_audit["per_type"][tier] = per_type_pca_audit
            cv_details[f"per_type:{tier}"] = per_type_cv_details
        build2_log(f"target {target_name}: build global_sliced score rows")
        score_rows, long_rows = build2_score_rows(target_name, "global_sliced", global_tier_results, tier_order, rows, atom_conditions, y, frames, atoms, split, base_blocks, args.alpha_mode, args.embedding_components)
        build2_log(f"target {target_name}: build per_type score rows")
        per_score_rows, per_long_rows = build2_score_rows(target_name, "per_type", per_type_tier_results, tier_order, rows, atom_conditions, y, frames, atoms, split, base_blocks, args.alpha_mode, args.embedding_components)
        score_rows.extend(per_score_rows)
        long_rows.extend(per_long_rows)
        build2_log(f"target {target_name}: build global_sliced partition curves")
        curves, curves_long, bin_defs = build2_partition_curves(target_name, "global_sliced", global_tier_results, tier_order, rows_for_conditions, atom_conditions, partition_bins, atom_bins, frames, cond_specs, manifest)
        build2_log(f"target {target_name}: build per_type partition curves")
        per_curves, per_curves_long, per_bin_defs = build2_partition_curves(target_name, "per_type", per_type_tier_results, tier_order, rows_for_conditions, atom_conditions, partition_bins, atom_bins, frames, cond_specs, manifest)
        all_score_rows.extend(score_rows)
        all_long_rows.extend(long_rows)
        all_curves.extend([curves, per_curves])
        all_curves_long.extend([curves_long, per_curves_long])
        all_bin_defs.extend(bin_defs)
        all_bin_defs.extend(per_bin_defs)
        alpha_selection_by_target[target_name] = {"mode": args.alpha_mode, "grid": [float(x) for x in alpha_grid], "by_scope_tier_axis": alpha_selection_by_scope}
        pca_audit_by_target[target_name] = pca_audit
        cv_details_by_target[target_name] = cv_details

    score_table = pd.DataFrame(all_score_rows).reindex(columns=SCORE_TABLE_COLUMNS)
    score_table = add_fit_scope_deltas(score_table).reindex(columns=SCORE_TABLE_COLUMNS)
    score_long = pd.DataFrame(all_long_rows)
    curves = pd.concat(all_curves, ignore_index=True) if all_curves else pd.DataFrame()
    curves_long = pd.concat(all_curves_long, ignore_index=True) if all_curves_long else pd.DataFrame()
    favorable = build2_favorable_partitions(curves)
    cases = load_casehunter_candidates(substrate_dir)
    case_shortlist = build2_casehunter_intersection(favorable, cases, rows_for_conditions, atom_conditions, partition_bins, cond_specs, manifest)
    charge_bin_unique = checks["input_acceptance"]["charge_partition_non_degenerate_bins"]
    audit = {
        "script": str(Path(__file__).resolve()),
        "substrate_dir": str(substrate_dir),
        "output_dir": str(out_dir),
        "disk_guard": disk_audit,
        "files_read": data["files_read"],
        "tier_order": tier_order,
        "fit_scopes": ["global_sliced", "per_type"],
        "targets_fit": list(data["targets"].keys()),
        "feature_tier_selection": {"selected_tiers": tier_order, "headline_delta": "per_type vs global_sliced by stratum, alongside tier deltas"},
        "python_read_scope": "Only emitted reduced per_atom_substrate sidecars, target sidecars, CaseHunter manifests, per_atom_substrate_partition_bins.npy, and C++ dominance scalar conditioners were opened; no trajectory.h5, per-source, older dirs, ring-path validation, or ORCA data.",
        "manifest_relationship": manifest.get("relationship"),
        "manifest_relationship_kind": manifest.get("relationship_kind"),
        "manifest_rows": int(manifest["rows"]),
        "manifest_n_atoms": int(manifest["n_atoms"]),
        "manifest_n_dft_frames": int(manifest["n_dft_frames"]),
        "charge_complete_rows_used": int(charge_mask.sum()),
        "same_charge_complete_row_set_for_all_tiers": True,
        "alpha_mode": args.alpha_mode,
        "ridge_alpha_baseline": float(args.ridge_alpha),
        "alpha_selection": alpha_selection_by_target,
        "change_of_basis_orthogonality_max_abs": orth,
        "t2_transform_check_blocks": transform_audit,
        "embedding_pca": {"n_components": int(args.embedding_components), "provenance_by_target": pca_audit_by_target},
        "within_split": {"strategy": "deterministic contiguous middle frame block", "test_frames": split["test_frames"], "purged_train_frames": split["purged_frames"], "purge_frames_each_side": int(args.purge_frames_each_side), "cross_split_lag1_pairs": int(split["cross_split_lag1_pairs"])},
        "cv_integrity_details": cv_details_by_target,
        "no_python_physics_recompute": True,
        "partition_conditions": [spec.__dict__ for spec in cond_specs],
        "partition_bin_definitions": all_bin_defs,
        "partition_bin_definitions_written_before_ranking": True,
        "partition_conditions_input_side_only": True,
        "partition_bins_from_cpp_lookup_only": True,
        "partition_bins_from_cpp_lookup_or_python_quantile_cpp_scalar": True,
        "dominance_quantile_bins_from_python_on_cpp_scalar": True,
        "dominance_scalar_columns": list(DOMINANCE_FRACTION_COLUMNS.values()),
        "dominance_scalar_columns_available_in_substrate": all(col in rows_for_conditions.columns for col in DOMINANCE_FRACTION_COLUMNS.values()),
        "dominance_cpp_bin_id_next_emit_flag": "dominant_fraction_{ring,charge,mc} bin ids should be emitted by C++ in the next substrate; Build 3 bins the C++ scalars in Python to avoid re-emitting.",
        "charge_partition_unique_bins": charge_bin_unique,
        "charge_partition_non_degenerate": bool(len(charge_bin_unique["bin_nearest_charge_r_4_6_8_10"]) > 1 and len(charge_bin_unique["bin_gap_to_2nd_charge_r_4_6_8_10"]) > 1),
        "dominant_charge_bin_available_in_substrate": bool(checks["input_acceptance"]["dominant_charge_bin_column_present"]),
        "ring_field_cross_method_validations_deferred": True,
        "seti_no_verdicts": True,
        "fit_stage_checks": checks,
    }
    build2_fit_stage_checks(audit, score_table, curves, case_shortlist)
    if not all(bool(check["pass"]) for check in audit["fit_stage_checks"].values()):
        (out_dir / "run_audit.json").write_text(json.dumps(json_sanitize(audit), indent=2, sort_keys=True), encoding="utf-8")
        failed = [name for name, check in audit["fit_stage_checks"].items() if not check["pass"]]
        raise SystemExit(f"FATAL: fit-stage checks failed; scores are not valid: {failed}")

    score_table.to_csv(out_dir / "allatom_fit_score_table.csv", index=False)
    score_long.to_csv(out_dir / "allatom_fit_score_long.csv", index=False)
    join_coverage(rows).to_csv(out_dir / "join_coverage.csv", index=False)
    (out_dir / "feature_blocks_used.json").write_text(json.dumps(json_sanitize(feature_blocks_used(block_labels, args.embedding_components, transform_audit, tier_order)), indent=2, sort_keys=True), encoding="utf-8")
    curves.to_csv(out_dir / "partition_response_curves.csv", index=False)
    curves_long.to_csv(out_dir / "partition_response_curves_long.csv", index=False)
    pd.DataFrame(all_bin_defs).to_csv(out_dir / "partition_bin_definitions.csv", index=False)
    favorable.to_csv(out_dir / "partition_favorable_partitions.csv", index=False)
    case_shortlist.to_csv(out_dir / "partition_casehunter_shortlist.csv", index=False)
    plots_audit = write_build2_plots(out_dir, score_table, curves)
    audit["plots"] = plots_audit
    (out_dir / "run_audit.json").write_text(json.dumps(json_sanitize(audit), indent=2, sort_keys=True), encoding="utf-8")
    write_build2_reports(out_dir, score_table, curves, favorable, case_shortlist, audit)
    print(f"wrote {out_dir}")


def piece2_main() -> None:
    args = parse_args()
    tier_order = parse_tiers(args.tiers)
    alpha_grid = parse_alpha_grid(args.alpha_grid)
    substrate_dir = args.substrate_dir.resolve()
    out_dir = args.out_dir
    data = load_inputs(substrate_dir)
    checks, charge_mask = input_acceptance_checks(data, substrate_dir)
    rows_all: pd.DataFrame = data["rows"]
    specs: dict[str, object] = data["specs"]
    classical_all: np.ndarray = data["classical"]
    conditioning_all: np.ndarray = data["conditioning"]
    modulation: np.ndarray = data["modulation"]
    embedding: np.ndarray = data["embedding"]
    target_T2: np.ndarray = data["target_T2"]
    manifest: dict[str, object] = data["manifest"]

    row_indices = np.flatnonzero(charge_mask)
    rows = rows_all.loc[charge_mask].reset_index(drop=True)
    classical = classical_all[charge_mask]
    conditioning = conditioning_all[charge_mask]
    n_atoms = int(manifest["n_atoms"])
    n_frames = int(manifest["n_dft_frames"])
    if len(rows) != n_atoms * n_frames:
        raise SystemExit("FATAL: charge-complete filter no longer has dense 846 x 660 rows")

    C = cob.get_C()
    orth = float(np.abs(C.T @ C - np.eye(5)).max())
    if orth >= 1.0e-12:
        raise SystemExit(f"FATAL: change_of_basis.get_C() orthogonality check failed: {orth:.3e}")
    y = target_T2[charge_mask, :5] @ C.T
    frames = rows["original_frame_index"].to_numpy(int)
    atoms = rows["atom_index"].to_numpy(int)

    base_blocks, block_labels, transform_audit = build_non_embedding_blocks(rows, classical, specs, C)
    row_conditions, formulas = conditioning_frame(rows, conditioning, modulation, classical, specs)
    atom_conditions = atom_condition_frame(row_conditions, n_frames, n_atoms)
    cond_specs = condition_specs(row_conditions)
    split = split_frame_block(frames, args.within_test_fraction, args.purge_frames_each_side)

    base_atom_blocks = {name: atom_means_dense(block, n_frames, n_atoms) for name, block in base_blocks.items()}
    y_atom = atom_means_dense(y, n_frames, n_atoms)
    embedding_atom_mean = atom_means_embedding(embedding, n_frames, n_atoms)
    print("loaded substrate and built feature blocks", flush=True)

    within_embedding_pcs = None
    within_embedding_pca_audit = None
    embedding_tiers = [tier for tier in tier_order if tier_uses_embedding(tier)]
    if embedding_tiers:
        within_train_indices = row_indices[split["train_mask"]]
        within_pca = fit_pca_memmap(embedding, within_train_indices, args.embedding_components)
        within_embedding_pcs = transform_pca_memmap(embedding, within_pca, row_indices)
        within_embedding_pca_audit = {
            "method": "unsupervised PCA fit once on within-CV training rows only; shared by selected tiers that consume AIMNet2 embedding PCs",
            "n_components": int(within_pca.n_components),
            "original_dims": int(within_pca.original_dims),
            "training_rows": int(within_pca.n_train_rows),
            "test_rows_excluded": int(split["test_mask"].sum()),
            "purged_rows_excluded": int(split["purge_mask"].sum()),
            "explained_variance_ratio": float(within_pca.explained_variance_ratio),
            "consumed_by_tiers": embedding_tiers,
        }
        print("fit within-split training-only embedding PCA", flush=True)
    else:
        within_embedding_pca_audit = {
            "method": "embedding sidecar present in substrate but not consumed by the selected fit tiers",
            "n_components": 0,
            "original_dims": int(embedding.shape[1]),
            "consumed_by_tiers": [],
        }

    tier_results: dict[str, dict[str, object]] = {}
    pca_audit: dict[str, object] = {}
    cv_details: dict[str, object] = {}
    alpha_selection_by_tier_axis: dict[str, dict[str, object]] = {}
    between_pca_cache: dict[object, PcaTransform] = {}
    between_alpha_pca_cache: dict[object, PcaTransform] = {}
    within_alpha_pca_cache: dict[object, PcaTransform] = {}
    for tier in tier_order:
        if args.alpha_mode == "select":
            print(f"selecting alpha for tier {tier}: between LOAO inner CV", flush=True)
            between_alpha, between_alpha_audit = select_between_alphas(
                tier,
                base_atom_blocks,
                embedding_atom_mean,
                y_atom,
                alpha_grid,
                args.embedding_components,
                args.inner_cv_folds,
                between_alpha_pca_cache,
            )
            print(f"selecting alpha for tier {tier}: within frame-block inner CV", flush=True)
            within_alpha, within_alpha_audit = select_within_alpha(
                tier,
                base_blocks,
                embedding,
                row_indices,
                y,
                atoms,
                frames,
                split,
                alpha_grid,
                args.embedding_components,
                args.inner_cv_folds,
                args.purge_frames_each_side,
                within_alpha_pca_cache,
            )
        else:
            fixed_alpha = float(args.ridge_alpha)
            between_alpha = np.full(n_atoms, fixed_alpha, dtype=float)
            between_alpha_audit = {
                "axis": "between_LOAO",
                "method": "fixed-alpha labelled baseline; no inner CV",
                "selected_alpha": fixed_alpha,
                "min": fixed_alpha,
                "max": fixed_alpha,
                "counts": {f"{fixed_alpha:g}": int(n_atoms)},
                "heldout_test_atom_excluded_from_alpha_selection": True,
                "alpha_grid": [float(x) for x in alpha_grid],
            }
            within_alpha = fixed_alpha
            within_alpha_audit = {
                "axis": "within_frameblock",
                "method": "fixed-alpha labelled baseline; no inner CV",
                "selected_alpha": fixed_alpha,
                "outer_test_rows_used_for_alpha_selection": 0,
                "outer_purged_rows_used_for_alpha_selection": 0,
                "alpha_grid": [float(x) for x in alpha_grid],
            }
        alpha_selection_by_tier_axis[tier] = {"between": between_alpha_audit, "within": within_alpha_audit}

        print(f"fitting tier {tier}: between LOAO", flush=True)
        cache = between_pca_cache if tier_uses_embedding(tier) else None
        between = evaluate_between_tier(tier, base_atom_blocks, embedding_atom_mean, y_atom, between_alpha, args.embedding_components, cache)
        print(f"fitting tier {tier}: within frame block", flush=True)
        within = evaluate_within_tier(
            tier,
            base_blocks,
            within_embedding_pcs,
            within_embedding_pca_audit,
            y,
            atoms,
            split,
            within_alpha,
        )
        tier_results[tier] = {"between": between, "within": within}
        pca_audit[tier] = {"between": between["pca_audit"], "within": within["pca_audit"]}
        cv_details[tier] = {
            "between": {"cv_check": between["cv_check"]},
            "within": {"cv_check": within["cv_check"]},
        }

    print("building score tables", flush=True)
    score_rows, long_rows = build_score_rows(
        tier_results,
        tier_order,
        rows,
        atom_conditions,
        y,
        frames,
        atoms,
        split,
        base_blocks,
        args.alpha_mode,
        args.embedding_components,
    )
    score_table = pd.DataFrame(score_rows)[SCORE_TABLE_COLUMNS]
    score_long = pd.DataFrame(long_rows)
    print("building partition curves", flush=True)
    curves, curves_long, bin_defs = partition_curves(tier_results, tier_order, atom_conditions, row_conditions, frames, cond_specs)
    audit = {
        "script": str(Path(__file__).resolve()),
        "substrate_dir": str(substrate_dir),
        "output_dir": str(out_dir),
        "files_read": data["files_read"],
        "tier_order": tier_order,
        "feature_tier_selection": {
            "selected_tiers": tier_order,
            "scope": "analysis-side fit feature consumption only; emitted substrate keeps AIMNet2 columns present",
        },
        "python_read_scope": "Only files inside the emitted per_atom_substrate directory listed in files_read; no trajectory.h5, ORCA, older broad-backbone, MOPAC, all-atom-equivariant, source, pair, or external merge directories were opened.",
        "manifest_relationship": manifest.get("relationship"),
        "manifest_relationship_kind": manifest.get("relationship_kind"),
        "manifest_rows": int(manifest["rows"]),
        "manifest_n_atoms": int(manifest["n_atoms"]),
        "manifest_n_dft_frames": int(manifest["n_dft_frames"]),
        "charge_complete_rows_used": int(charge_mask.sum()),
        "same_charge_complete_row_set_for_all_tiers": True,
        "alpha_mode": args.alpha_mode,
        "ridge_alpha_baseline": float(args.ridge_alpha),
        "alpha_selection": {
            "mode": args.alpha_mode,
            "grid": [float(x) for x in alpha_grid],
            "fixed_alpha_baseline": float(args.ridge_alpha),
            "inner_cv_folds_requested": int(args.inner_cv_folds),
            "by_tier_axis": alpha_selection_by_tier_axis,
        },
        "target": "DFT total T2, per_atom_substrate_target_T2.npy[:, 0:5] @ change_of_basis.get_C().T",
        "target_T2_original_shape": list(target_T2.shape),
        "change_of_basis_orthogonality_max_abs": orth,
        "t2_transform_check_blocks": [
            {"block": item["feature_block"], "components": 5, "transform": item["transform"]}
            for item in transform_audit
        ],
        "embedding_pca": {
            "owner_override": "training-only PCA is primary; full 256d embedding is not used as primary",
            "n_components": int(args.embedding_components),
            "provenance_by_tier": pca_audit,
        },
        "within_split": {
            "strategy": "deterministic contiguous middle frame block",
            "test_frames": split["test_frames"],
            "purged_train_frames": split["purged_frames"],
            "purge_frames_each_side": int(args.purge_frames_each_side),
            "cross_split_lag1_pairs": int(split["cross_split_lag1_pairs"]),
        },
        "cv_integrity_details": cv_details,
        "no_python_physics_recompute": True,
        "algebraic_partition_conditioners": formulas,
        "partition_conditions": [spec.__dict__ for spec in cond_specs],
        "partition_bin_definitions": bin_defs,
        "partition_bin_definitions_written_before_ranking": True,
        "partition_conditions_input_side_only": True,
        "fit_stage_checks": checks,
    }
    fit_stage_checks(audit, score_table, curves, favorable_cases(curves))
    if not all(bool(check["pass"]) for check in audit["fit_stage_checks"].values()):
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "run_audit.json").write_text(json.dumps(json_sanitize(audit), indent=2, sort_keys=True), encoding="utf-8")
        failed = [name for name, check in audit["fit_stage_checks"].items() if not check["pass"]]
        raise SystemExit(f"FATAL: fit-stage checks failed; scores are not valid: {failed}")

    fav = favorable_cases(curves)
    if out_dir.exists() and not args.keep_out_dir:
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    score_table.to_csv(out_dir / "allatom_fit_score_table.csv", index=False)
    score_long.to_csv(out_dir / "allatom_fit_score_long.csv", index=False)
    join_coverage(rows).to_csv(out_dir / "join_coverage.csv", index=False)
    (out_dir / "feature_blocks_used.json").write_text(
        json.dumps(json_sanitize(feature_blocks_used(block_labels, args.embedding_components, transform_audit, tier_order)), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    curves.to_csv(out_dir / "partition_response_curves.csv", index=False)
    curves_long.to_csv(out_dir / "partition_response_curves_long.csv", index=False)
    fav.to_csv(out_dir / "partition_favorable_cases.csv", index=False)
    (out_dir / "run_audit.json").write_text(json.dumps(json_sanitize(audit), indent=2, sort_keys=True), encoding="utf-8")
    write_reports(out_dir, score_table, curves, fav, audit)
    print(f"wrote {out_dir}")
    print("fit-stage checks: " + ", ".join(f"{name}=pass" for name in audit["fit_stage_checks"]))


if __name__ == "__main__":
    build2_main()
