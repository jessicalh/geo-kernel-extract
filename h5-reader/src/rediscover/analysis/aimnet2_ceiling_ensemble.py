#!/usr/bin/env python3
"""AIMNet2 ceiling and ensemble diagnostic for rediscover substrates.

Reads only emitted CSV/NPY payloads. No H5, no kernel reconstruction, no field
or projection recompute. T2 targets and T2 mechanism sidecars are mapped through
the frozen library->e3nn 2e basis from change_of_basis.get_C().

The AIMNet2 embedding rows are reported as a learnable ceiling diagnostic, not
as a recovered physical law. The AIMNet2 charge-response-gradient is labelled
CRG throughout; it is not a polarizability.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

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
TEST_FRAME_FRACTION = 0.20

PHYSICS_SCALAR_COLS = [
    "ring_sum_dipolar",
    "bond_sum_dipolar",
    "field_local_x",
    "field_local_y",
    "field_local_z",
    "field_z",
    "field_mag_sq",
    "buckingham_E_local_x",
    "buckingham_E_local_y",
    "buckingham_E_local_z",
    "buckingham_E_proj",
    "buckingham_E_mag_sq",
]

FEATURE_SET_COMPONENTS = {
    "ensemble_physics_only": (),
    "ensemble_plus_crg": ("crg",),
    "ensemble_plus_embedding": ("embedding",),
    "ensemble_plus_aimnet2_charge": ("charge",),
    "ensemble_plus_crg_embedding": ("crg", "embedding"),
    "ensemble_plus_aimnet2": ("charge", "crg", "embedding"),
}

FEATURE_SET_ORDER = list(FEATURE_SET_COMPONENTS)
REQUIRED_FEATURE_SETS = {
    "ensemble_physics_only",
    "ensemble_plus_crg",
    "ensemble_plus_embedding",
    "ensemble_plus_aimnet2_charge",
    "ensemble_plus_aimnet2",
}


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--broad-dir", default="/tmp/rdc-broad-backbone-axes")
    ap.add_argument("--buckingham-dir", default="/tmp/rdc-buckingham-capstone")
    ap.add_argument("--aimnet-dir", default="/tmp/rdc-aimnet2-features")
    ap.add_argument("--out-dir", default="/tmp/rediscover-aimnet2-ensemble")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--physics-alpha", type=float, default=1.0e-6)
    ap.add_argument("--aimnet-alpha", type=float, default=10.0)
    ap.add_argument("--embedding-components", type=int, default=32,
                    help="Unsupervised PCA components from the emitted 256-d embedding used in fits; 0 keeps all 256.")
    ap.add_argument("--within-split", choices=["random", "blocked"], default="random")
    ap.add_argument("--test-frame-fraction", type=float, default=TEST_FRAME_FRACTION)
    return ap.parse_args()


def stratum_of(frame_variant):
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


def finite_fmt(x):
    if x is None or not np.isfinite(x):
        return "nan"
    return f"{x:.4f}"


def r2_score(pred, target):
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


def manual_r2(pred, target):
    pred = np.asarray(pred, dtype=float)
    target = np.asarray(target, dtype=float)
    if pred.ndim == 1:
        pred = pred.reshape(-1, 1)
    if target.ndim == 1:
        target = target.reshape(-1, 1)
    ok = np.isfinite(pred).all(axis=1) & np.isfinite(target).all(axis=1)
    y = target[ok]
    p = pred[ok]
    return 1.0 - float(((y - p) ** 2).sum()) / float(((y - y.mean(axis=0, keepdims=True)) ** 2).sum())


def manual_between_loao_r2_from_features(x, y, atoms, alpha):
    labels, xbar, _ = atom_means(x, atoms)
    _, ybar, _ = atom_means(y, atoms)
    pred = np.full_like(ybar, np.nan, dtype=float)
    for i in range(len(labels)):
        train = np.ones(len(labels), dtype=bool)
        train[i] = False
        pred[i] = fit_ridge_predict(xbar[train], ybar[train], xbar[i:i + 1], alpha)[0]
    return manual_r2(pred, ybar)


def jackknife_metric_values(pred, target, groups):
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


def jackknife_se_from_values(vals):
    vals = np.asarray(vals, dtype=float)
    if len(vals) < 3:
        return math.nan
    mean = float(vals.mean())
    return float(math.sqrt((len(vals) - 1) / len(vals) * np.sum((vals - mean) ** 2)))


def jackknife_metric_se(pred, target, groups):
    return jackknife_se_from_values(jackknife_metric_values(pred, target, groups))


def jackknife_difference_se(pred, target, groups, base_pred, base_target, base_groups):
    pred = np.asarray(pred, dtype=float)
    target = np.asarray(target, dtype=float)
    groups = np.asarray(groups)
    base_pred = np.asarray(base_pred, dtype=float)
    base_target = np.asarray(base_target, dtype=float)
    base_groups = np.asarray(base_groups)
    if (
        len(pred) != len(base_pred)
        or target.shape != base_target.shape
        or len(groups) != len(base_groups)
        or not np.array_equal(groups, base_groups)
        or not np.allclose(target, base_target, equal_nan=True)
    ):
        return math.nan
    labels = np.unique(groups)
    if len(labels) < 3:
        return math.nan
    vals = []
    for label in labels:
        keep = groups != label
        if keep.sum() < 3:
            continue
        vals.append(
            r2_score(pred[keep], target[keep])
            - r2_score(base_pred[keep], base_target[keep])
        )
    vals = np.asarray(vals, dtype=float)
    return jackknife_se_from_values(vals[np.isfinite(vals)])


def atom_means(values, atoms):
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


def centred_by_train_atom(values, atoms, train):
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


def variance_shares(target, atoms):
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
    between = 0.0
    within = 0.0
    for atom in np.unique(a):
        m = a == atom
        mean = y[m].mean(axis=0, keepdims=True)
        between += float(m.sum() * ((mean - grand) ** 2).sum())
        within += float(((y[m] - mean) ** 2).sum())
    if total <= 0.0:
        return math.nan, math.nan
    return between / total, within / total


def split_masks(frames, strategy, seed, test_fraction):
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
    return train, test, {"test_frames": len(test_frames), "cross_split_lag1_pairs": cross}


def fit_ridge_predict(x_train, y_train, x_test, alpha):
    x_train = np.asarray(x_train, dtype=float)
    y_train = np.asarray(y_train, dtype=float)
    x_test = np.asarray(x_test, dtype=float)
    if y_train.ndim == 1:
        y_train = y_train.reshape(-1, 1)

    y_ok = np.isfinite(y_train).all(axis=1)
    x_ok = np.isfinite(x_train).any(axis=1)
    ok = y_ok & x_ok
    x_train = x_train[ok]
    y_train = y_train[ok]
    if len(x_train) < 2:
        return np.full((len(x_test), y_train.shape[1] if y_train.ndim == 2 else 1), np.nan)

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
    if alpha > 0.0:
        penalty = np.eye(design.shape[1]) * alpha
        penalty[0, 0] = 0.0
        beta = np.linalg.solve(design.T @ design + penalty, design.T @ y_train)
    else:
        beta = np.linalg.lstsq(design, y_train, rcond=None)[0]
    return design_t @ beta


def evaluate_between(x, y, atoms, alpha):
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


def evaluate_within(x, y, atoms, frames, alpha, split_strategy, seed, test_fraction):
    train, test, split = split_masks(frames, split_strategy, seed, test_fraction)
    xw = centred_by_train_atom(x, atoms, train)
    yw = centred_by_train_atom(y, atoms, train)
    fit = train & np.isfinite(yw).all(axis=1) & np.isfinite(xw).any(axis=1)
    score = test & np.isfinite(yw).all(axis=1) & np.isfinite(xw).any(axis=1)
    if fit.sum() < 3 or score.sum() < 3:
        pred = np.full((score.sum(), yw.shape[1]), np.nan)
        target = yw[score]
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


def require_file(path):
    path = Path(path)
    if not path.exists():
        raise SystemExit(f"FATAL: required emitted substrate missing: {path}")
    return path


def t2_sidecar(directory, name, C):
    return np.load(require_file(Path(directory) / name)) @ C.T


def reduce_embedding(embedding, n_components):
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
    k = min(n_components, original_dims)
    total = float(np.sum(vals[vals > 0.0]))
    ratio = float(np.sum(vals[:k]) / total) if total > 0.0 else math.nan
    return x @ vecs[:, :k], {
        "method": "global_unsupervised_pca_on_emitted_embedding",
        "original_dims": int(original_dims),
        "used_dims": int(k),
        "explained_variance_ratio": ratio,
    }


def load_joined(broad_dir, buckingham_dir, aimnet_dir, C):
    broad_dir = Path(broad_dir)
    buckingham_dir = Path(buckingham_dir)
    aimnet_dir = Path(aimnet_dir)

    broad = pd.read_csv(require_file(broad_dir / "broad_backbone_aggregated.csv"))
    broad["broad_row_id"] = np.arange(len(broad), dtype=int)
    broad["stratum"] = [stratum_of(v) for v in broad["frame_variant"]]

    buck = pd.read_csv(require_file(buckingham_dir / "buckingham_efield_aggregated.csv"))
    buck_cols = [
        "atom_index", "h5_row", "apbs_efield_present",
        "efield_local_x", "efield_local_y", "efield_local_z", "E_proj", "E_mag",
    ]
    buck = buck[buck_cols].rename(columns={
        "apbs_efield_present": "buckingham_present",
        "efield_local_x": "buckingham_E_local_x",
        "efield_local_y": "buckingham_E_local_y",
        "efield_local_z": "buckingham_E_local_z",
        "E_proj": "buckingham_E_proj",
        "E_mag": "buckingham_E_mag",
    })

    aim = pd.read_csv(require_file(aimnet_dir / "aimnet2_features_aggregated.csv"))
    aim = aim.rename(columns={"row_id": "aimnet2_row_id"})
    aim_cols = [
        "atom_index", "h5_row", "aimnet2_row_id",
        "aimnet2_charge_present", "aimnet2_charge",
        "aimnet2_charge_response_gradient_present",
        "aimnet2_charge_response_gradient_scalar",
        "aimnet2_charge_response_gradient_local_x",
        "aimnet2_charge_response_gradient_local_y",
        "aimnet2_charge_response_gradient_local_z",
        "aimnet2_embedding_present", "aimnet2_embedding_dim",
    ]
    aim = aim[aim_cols]

    df = broad.merge(buck, on=["atom_index", "h5_row"], how="left", sort=False)
    df = df.merge(aim, on=["atom_index", "h5_row"], how="left", sort=False)
    if len(df) != len(broad):
        raise SystemExit(f"FATAL: join changed row count: broad={len(broad)} joined={len(df)}")
    if df["aimnet2_row_id"].isna().any():
        raise SystemExit("FATAL: AIMNet2 feature export does not cover every broad_backbone row")

    field_mag = df["field_mag"].to_numpy(float)
    buck_mag = df["buckingham_E_mag"].to_numpy(float)
    df["field_mag_sq"] = field_mag * field_mag
    df["buckingham_E_mag_sq"] = buck_mag * buck_mag

    row_idx = df["broad_row_id"].to_numpy(int)
    aim_idx = df["aimnet2_row_id"].to_numpy(int)

    arrays = {
        "target_T2": t2_sidecar(broad_dir, "broad_backbone_aggregated_target_local_T2.npy", C)[row_idx],
        "ring_T2": t2_sidecar(broad_dir, "broad_backbone_aggregated_ring_literature_kernel_T2.npy", C)[row_idx],
        "bond_T2": t2_sidecar(broad_dir, "broad_backbone_aggregated_bond_literature_kernel_T2.npy", C)[row_idx],
        "charge_T2": t2_sidecar(broad_dir, "broad_backbone_aggregated_charge_literature_kernel_T2.npy", C)[row_idx],
        "embedding": np.load(require_file(aimnet_dir / "aimnet2_features_embedding.npy"))[aim_idx],
        "crg_local": np.load(require_file(aimnet_dir / "aimnet2_features_charge_response_gradient_local.npy"))[aim_idx],
        "crg_scalar": np.load(require_file(aimnet_dir / "aimnet2_features_charge_response_gradient_scalar.npy"))[aim_idx].reshape(-1, 1),
    }
    if arrays["embedding"].shape[1] != 256:
        raise SystemExit(f"FATAL: AIMNet2 embedding sidecar has shape {arrays['embedding'].shape}, expected second dim 256")
    return df, arrays


def physics_blocks(df, arrays):
    parts = [
        df[PHYSICS_SCALAR_COLS].to_numpy(float),
        arrays["ring_T2"],
        arrays["bond_T2"],
        arrays["charge_T2"],
    ]
    labels = [f"physics:{c}" for c in PHYSICS_SCALAR_COLS]
    labels += [f"physics:ring_T2_{i}" for i in range(5)]
    labels += [f"physics:mcconnell_T2_{i}" for i in range(5)]
    labels += [f"physics:charge_field_gradient_T2_{i}" for i in range(5)]
    return parts, labels


def feature_blocks(df, arrays, feature_set):
    if feature_set not in FEATURE_SET_COMPONENTS:
        raise ValueError(f"unknown feature_set {feature_set}")
    parts, labels = physics_blocks(df, arrays)
    components = FEATURE_SET_COMPONENTS[feature_set]
    if "charge" in components:
        parts.append(df[["aimnet2_charge"]].to_numpy(float))
        labels.append("aimnet2:hirshfeld_charge")
    if "crg" in components:
        parts.extend([arrays["crg_local"], arrays["crg_scalar"]])
        labels += [f"aimnet2:charge_response_gradient_local_{c}" for c in "xyz"]
        labels.append("aimnet2:charge_response_gradient_scalar")
    if "embedding" in components:
        parts.append(arrays["embedding"].astype(float))
        labels += [f"aimnet2:embedding_pc_{i:03d}" for i in range(arrays["embedding"].shape[1])]
    return np.column_stack(parts), labels


def result_row(target_name, stratum, feature_set, x, y, atoms, frames, alpha,
               split_strategy, seed, test_fraction, feature_count, labels):
    shares = variance_shares(y, atoms)
    between = evaluate_between(x, y, atoms, alpha)
    within = evaluate_within(x, y, atoms, frames, alpha, split_strategy, seed, test_fraction)
    return {
        "target": target_name,
        "stratum": stratum,
        "feature_set": feature_set,
        "feature_components": "+".join(FEATURE_SET_COMPONENTS[feature_set]) or "physics",
        "rows": int(len(y)),
        "n_atoms_between": int(len(np.unique(atoms))),
        "n_features": int(feature_count),
        "ridge_alpha": float(alpha),
        "variance_share_between": shares[0],
        "variance_share_within": shares[1],
        "between_LOAO_R2": between["R2"],
        "between_LOAO_R2_jackknife_se": between["jackknife_se"],
        "within_frame_R2": within["R2"],
        "within_frame_R2_jackknife_se": within["jackknife_se"],
        "within_split_strategy": split_strategy,
        "test_frames": int(within.get("test_frames", 0)),
        "cross_split_lag1_pairs": int(within.get("cross_split_lag1_pairs", 0)),
        "feature_labels": labels,
        "_between_pred": between["pred"],
        "_between_target": between["target"],
        "_between_groups": between["groups"],
        "_within_pred": within["pred"],
        "_within_target": within["target"],
        "_within_groups": within["groups"],
    }


def strip_private(row):
    return {k: v for k, v in row.items() if not k.startswith("_") and k != "feature_labels"}


def add_physics_lifts(rows):
    base_rows = {
        (row["target"], row["stratum"]): row
        for row in rows
        if row["feature_set"] == "ensemble_physics_only"
    }
    for row in rows:
        base = base_rows.get((row["target"], row["stratum"]))
        if base is None:
            row["gap_between_R2_vs_physics"] = math.nan
            row["gap_within_R2_vs_physics"] = math.nan
            row["gap_between_R2_jackknife_se"] = math.nan
            row["gap_within_R2_jackknife_se"] = math.nan
            continue
        row["gap_between_R2_vs_physics"] = row["between_LOAO_R2"] - base["between_LOAO_R2"]
        row["gap_within_R2_vs_physics"] = row["within_frame_R2"] - base["within_frame_R2"]
        if row["feature_set"] == "ensemble_physics_only":
            row["gap_between_R2_jackknife_se"] = 0.0
            row["gap_within_R2_jackknife_se"] = 0.0
        else:
            row["gap_between_R2_jackknife_se"] = jackknife_difference_se(
                row["_between_pred"], row["_between_target"], row["_between_groups"],
                base["_between_pred"], base["_between_target"], base["_between_groups"],
            )
            row["gap_within_R2_jackknife_se"] = jackknife_difference_se(
                row["_within_pred"], row["_within_target"], row["_within_groups"],
                base["_within_pred"], base["_within_target"], base["_within_groups"],
            )


def safe_div(num, den):
    if den is None or not np.isfinite(den) or abs(den) <= 1.0e-12:
        return math.nan
    return num / den


def dominant_lift(row, prefix):
    lifts = {
        "CRG": row[f"{prefix}_crg_lift"],
        "embedding": row[f"{prefix}_embedding_lift"],
        "charge": row[f"{prefix}_charge_lift"],
    }
    finite = {k: v for k, v in lifts.items() if np.isfinite(v)}
    if not finite:
        return "nan"
    name, value = max(finite.items(), key=lambda kv: kv[1])
    if value <= 0.0:
        return "none_positive"
    return name


def build_decomposition(out):
    rows = []
    by_set = {name: name for name in FEATURE_SET_COMPONENTS}
    for (target, stratum), group in out.groupby(["target", "stratum"], sort=False):
        keyed = {row["feature_set"]: row for _, row in group.iterrows()}
        if not all(name in keyed for name in by_set):
            continue
        base = keyed["ensemble_physics_only"]
        row = {
            "target": target,
            "stratum": stratum,
            "n_atoms_between": int(base["n_atoms_between"]),
            "rows": int(base["rows"]),
            "variance_share_between": float(base["variance_share_between"]),
            "variance_share_within": float(base["variance_share_within"]),
        }
        for prefix, metric in [("between", "between_LOAO_R2"), ("within", "within_frame_R2")]:
            physics = float(keyed["ensemble_physics_only"][metric])
            crg = float(keyed["ensemble_plus_crg"][metric]) - physics
            embedding = float(keyed["ensemble_plus_embedding"][metric]) - physics
            charge = float(keyed["ensemble_plus_aimnet2_charge"][metric]) - physics
            crg_embedding = float(keyed["ensemble_plus_crg_embedding"][metric]) - physics
            combined = float(keyed["ensemble_plus_aimnet2"][metric]) - physics
            isolated_sum = crg + embedding + charge
            crg_embedding_separate_sum = crg + embedding
            crg_embedding_overlap = crg_embedding_separate_sum - crg_embedding
            isolated_overlap = isolated_sum - combined
            row.update({
                f"{prefix}_physics_R2": physics,
                f"{prefix}_crg_lift": crg,
                f"{prefix}_embedding_lift": embedding,
                f"{prefix}_charge_lift": charge,
                f"{prefix}_crg_embedding_lift": crg_embedding,
                f"{prefix}_combined_lift": combined,
                f"{prefix}_isolated_sum_lift": isolated_sum,
                f"{prefix}_crg_embedding_overlap": crg_embedding_overlap,
                f"{prefix}_crg_embedding_overlap_fraction_of_separate_sum": safe_div(
                    crg_embedding_overlap, crg_embedding_separate_sum
                ),
                f"{prefix}_isolated_overlap_vs_combined": isolated_overlap,
                f"{prefix}_crg_share_of_combined_lift": safe_div(crg, combined),
                f"{prefix}_embedding_share_of_combined_lift": safe_div(embedding, combined),
                f"{prefix}_charge_share_of_combined_lift": safe_div(charge, combined),
            })
        row["dominant_between_lift"] = dominant_lift(row, "between")
        row["dominant_within_lift"] = dominant_lift(row, "within")
        rows.append(row)
    return pd.DataFrame(rows)


def add_full_lift_fractions(out):
    full = out[out["feature_set"] == "ensemble_plus_aimnet2"][
        ["target", "stratum", "gap_between_R2_vs_physics", "gap_within_R2_vs_physics"]
    ].rename(columns={
        "gap_between_R2_vs_physics": "aimnet2_full_between_lift",
        "gap_within_R2_vs_physics": "aimnet2_full_within_lift",
    })
    out = out.merge(full, on=["target", "stratum"], how="left")
    out["between_lift_fraction_of_aimnet2_full"] = [
        safe_div(num, den)
        for num, den in zip(out["gap_between_R2_vs_physics"], out["aimnet2_full_between_lift"])
    ]
    out["within_lift_fraction_of_aimnet2_full"] = [
        safe_div(num, den)
        for num, den in zip(out["gap_within_R2_vs_physics"], out["aimnet2_full_within_lift"])
    ]
    return out


def write_markdown(df, path):
    cols = [
        "target", "stratum", "feature_set", "feature_components", "n_atoms_between", "n_features",
        "variance_share_between", "variance_share_within",
        "between_LOAO_R2", "between_LOAO_R2_jackknife_se",
        "within_frame_R2", "within_frame_R2_jackknife_se",
        "gap_between_R2_vs_physics", "gap_between_R2_jackknife_se",
        "gap_within_R2_vs_physics", "gap_within_R2_jackknife_se",
        "between_lift_fraction_of_aimnet2_full", "within_lift_fraction_of_aimnet2_full",
    ]
    lines = [
        "# AIMNet2 ceiling and ensemble",
        "",
        "AIMNet2 embedding columns are interpreted as a learnable ceiling diagnostic, not a recovered law.",
        "AIMNet2 CRG means charge-response-gradient, not polarizability.",
        "",
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join(["---"] * len(cols)) + " |",
    ]
    for _, row in df[cols].iterrows():
        vals = []
        for col in cols:
            v = row[col]
            if isinstance(v, float):
                vals.append("nan" if not np.isfinite(v) else f"{v:.4g}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_decomposition_markdown(df, path):
    cols = [
        "target", "stratum", "n_atoms_between",
        "between_crg_lift", "between_embedding_lift", "between_charge_lift",
        "between_combined_lift", "between_crg_embedding_overlap",
        "within_crg_lift", "within_embedding_lift", "within_charge_lift",
        "within_combined_lift", "within_crg_embedding_overlap",
        "dominant_between_lift", "dominant_within_lift",
    ]
    lines = [
        "# AIMNet2 isolated lift decomposition",
        "",
        "Lifts are delta R2 against ensemble_physics_only, per stratum and target.",
        "Positive CRG+embedding overlap means the isolated CRG and embedding lifts double-count shared signal; negative overlap means the joint CRG+embedding fit exceeds their isolated sum.",
        "",
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join(["---"] * len(cols)) + " |",
    ]
    for _, row in df[cols].iterrows():
        vals = []
        for col in cols:
            v = row[col]
            if isinstance(v, float):
                vals.append("nan" if not np.isfinite(v) else f"{v:.4g}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    C = cob.get_C()
    orth = float(np.abs(C.T @ C - np.eye(5)).max())
    print(f"change_of_basis.get_C() |C^T C - I|max={orth:.3e}")

    df, arrays = load_joined(args.broad_dir, args.buckingham_dir, args.aimnet_dir, C)
    arrays["embedding"], embedding_transform = reduce_embedding(arrays["embedding"], args.embedding_components)
    rows = []
    audit = None
    feature_labels = {}

    for stratum in STRATA_ORDER:
        mask = (
            (df["stratum"].to_numpy() == stratum)
            & (df["dft_present"].to_numpy(int) == 1)
            & (df["frame_valid"].to_numpy(int) == 1)
            & (df["dft_local_frame_valid"].to_numpy(int) == 1)
        )
        idx = np.flatnonzero(mask)
        if len(idx) < 4:
            continue
        d = df.iloc[idx].reset_index(drop=True)
        atoms = d["atom_index"].to_numpy(int)
        frames = d["h5_row"].to_numpy(int)

        local_arrays = {k: v[idx] for k, v in arrays.items()}
        feature_matrices = {}
        for feature_set in FEATURE_SET_ORDER:
            x, labels = feature_blocks(d, local_arrays, feature_set)
            feature_matrices[feature_set] = (x, labels)
            feature_labels[feature_set] = labels

        targets = {
            "sigma_iso": d["dft_sigma_iso"].to_numpy(float).reshape(-1, 1),
            "T2": arrays["target_T2"][idx],
        }
        for target_name, y in targets.items():
            for feature_set in FEATURE_SET_ORDER:
                x, labels = feature_matrices[feature_set]
                alpha = args.physics_alpha if feature_set == "ensemble_physics_only" else args.aimnet_alpha
                row = result_row(
                    target_name, stratum, feature_set, x, y, atoms, frames, alpha,
                    args.within_split, args.seed, args.test_frame_fraction,
                    x.shape[1], labels,
                )
                if audit is None and target_name == "sigma_iso" and stratum == "CA" and feature_set == "ensemble_plus_crg":
                    direct = manual_r2(row["_between_pred"], row["_between_target"])
                    emitted_direct = manual_between_loao_r2_from_features(x, y, atoms, alpha)
                    audit = {
                        "row": {
                            "target": target_name,
                            "stratum": stratum,
                            "feature_set": feature_set,
                        },
                        "source": "manual_R2_from_CA_sigma_iso_LOAO_predictions_fitted_from_emitted_physics_plus_charge_response_gradient_features",
                        "crg_feature_labels": [
                            label for label in labels if "charge_response_gradient" in label
                        ],
                        "reported_between_LOAO_R2": row["between_LOAO_R2"],
                        "manual_between_LOAO_R2_from_predictions": direct,
                        "manual_between_LOAO_R2_from_emitted_features": emitted_direct,
                        "abs_diff_predictions": abs(direct - row["between_LOAO_R2"]),
                        "abs_diff_emitted_features": abs(emitted_direct - row["between_LOAO_R2"]),
                    }
                rows.append(row)

    missing = REQUIRED_FEATURE_SETS - {row["feature_set"] for row in rows}
    if missing:
        raise SystemExit(f"FATAL: required feature sets missing: {sorted(missing)}")
    add_physics_lifts(rows)
    out = pd.DataFrame([strip_private(r) for r in rows])
    if out.empty:
        raise SystemExit("FATAL: no fitted rows produced")
    out["feature_set_order"] = out["feature_set"].map({name: i for i, name in enumerate(FEATURE_SET_ORDER)})
    out["stratum_order"] = out["stratum"].map({name: i for i, name in enumerate(STRATA_ORDER)})
    out = add_full_lift_fractions(out)
    out = out.sort_values(["target", "stratum_order", "feature_set_order"]).reset_index(drop=True)
    decomposition = build_decomposition(out)
    decomposition["stratum_order"] = decomposition["stratum"].map({name: i for i, name in enumerate(STRATA_ORDER)})
    decomposition = decomposition.sort_values(["target", "stratum_order"]).reset_index(drop=True)

    csv_path = out_dir / "aimnet2_ceiling_ensemble.csv"
    md_path = out_dir / "aimnet2_ceiling_ensemble.md"
    decomp_csv_path = out_dir / "aimnet2_lift_decomposition.csv"
    decomp_md_path = out_dir / "aimnet2_lift_decomposition.md"
    json_path = out_dir / "aimnet2_ceiling_ensemble_run.json"
    audit_path = out_dir / "self_audit.json"
    labels_path = out_dir / "feature_labels.json"
    out.drop(columns=["feature_set_order", "stratum_order"]).to_csv(csv_path, index=False)
    decomposition.drop(columns=["stratum_order"]).to_csv(decomp_csv_path, index=False)
    write_markdown(out, md_path)
    write_decomposition_markdown(decomposition, decomp_md_path)
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    labels_path.write_text(json.dumps(feature_labels, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    json_path.write_text(
        json.dumps(
            {
                "broad_dir": args.broad_dir,
                "buckingham_dir": args.buckingham_dir,
                "aimnet_dir": args.aimnet_dir,
                "seed": args.seed,
                "within_split": args.within_split,
                "test_frame_fraction": args.test_frame_fraction,
                "physics_alpha": args.physics_alpha,
                "aimnet_alpha": args.aimnet_alpha,
                "embedding_dims_emitted": 256,
                "embedding_transform": embedding_transform,
                "crg_label": "aimnet2_charge_response_gradient_not_polarizability",
                "embedding_interpretation": "learnable_ceiling_not_physical_law",
                "feature_set_components": FEATURE_SET_COMPONENTS,
                "required_feature_sets": sorted(REQUIRED_FEATURE_SETS),
                "crg_embedding_overlap_definition": "crg_lift + embedding_lift - crg_embedding_lift",
                "change_of_basis_orthogonality_max": orth,
                "rows": int(len(out)),
                "decomposition_rows": int(len(decomposition)),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"wrote {csv_path}")
    print(f"wrote {md_path}")
    print(f"wrote {decomp_csv_path}")
    print(f"wrote {decomp_md_path}")
    print(f"wrote {json_path}")
    print(f"wrote {audit_path}")
    print("\ntarget,stratum,feature_set,between_R2,within_R2,gap_between,gap_within")
    for _, row in out.iterrows():
        print(
            f"{row['target']},{row['stratum']},{row['feature_set']},"
            f"{finite_fmt(row['between_LOAO_R2'])},{finite_fmt(row['within_frame_R2'])},"
            f"{finite_fmt(row['gap_between_R2_vs_physics'])},"
            f"{finite_fmt(row['gap_within_R2_vs_physics'])}"
        )
    print("\ntarget,stratum,between_crg,between_embedding,between_charge,between_combined,between_crg_embedding_overlap,within_crg,within_embedding,within_charge,within_combined,within_crg_embedding_overlap")
    for _, row in decomposition.iterrows():
        print(
            f"{row['target']},{row['stratum']},"
            f"{finite_fmt(row['between_crg_lift'])},{finite_fmt(row['between_embedding_lift'])},"
            f"{finite_fmt(row['between_charge_lift'])},{finite_fmt(row['between_combined_lift'])},"
            f"{finite_fmt(row['between_crg_embedding_overlap'])},"
            f"{finite_fmt(row['within_crg_lift'])},{finite_fmt(row['within_embedding_lift'])},"
            f"{finite_fmt(row['within_charge_lift'])},{finite_fmt(row['within_combined_lift'])},"
            f"{finite_fmt(row['within_crg_embedding_overlap'])}"
        )


if __name__ == "__main__":
    main()
