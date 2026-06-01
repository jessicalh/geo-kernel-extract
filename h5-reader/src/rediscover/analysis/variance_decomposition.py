#!/usr/bin/env python3
"""Between/within variance decomposition for rediscover sidecar substrates.

Reads only emitted CSV/NPY payloads. No H5, no kernel reconstruction, no field
or projection recompute. Tensor features and targets are mapped through the
frozen library->e3nn 2e basis from change_of_basis.get_C().
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
FRAME_SPLIT_SEED = 0
THIN_ATOMS = 10


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--broad-dir", default="/tmp/rdc-broad-backbone-axes")
    ap.add_argument("--buckingham-dir", default="/tmp/rdc-buckingham")
    ap.add_argument("--efg-dir", default="/tmp/rdc-efg")
    ap.add_argument("--out-dir", default="/tmp/rediscover-variance-decomposition")
    ap.add_argument("--seed", type=int, default=FRAME_SPLIT_SEED)
    ap.add_argument("--efg-purge-frames", type=int, default=1)
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
    return f"{x:+.4f}"


def corrcoef(a, b):
    a = np.asarray(a, dtype=float).reshape(-1)
    b = np.asarray(b, dtype=float).reshape(-1)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return math.nan
    a = a[ok] - a[ok].mean()
    b = b[ok] - b[ok].mean()
    den = math.sqrt(float(np.dot(a, a) * np.dot(b, b)))
    if den <= 0.0:
        return math.nan
    return float(np.dot(a, b) / den)


def r2_score(pred, target):
    pred = np.asarray(pred, dtype=float)
    target = np.asarray(target, dtype=float)
    ok = np.isfinite(pred) & np.isfinite(target)
    if pred.ndim == 2:
        ok = np.isfinite(pred).all(axis=1) & np.isfinite(target).all(axis=1)
    if ok.sum() < 2:
        return math.nan
    p = pred[ok]
    y = target[ok]
    den = float(((y - y.mean(axis=0, keepdims=True)) ** 2).sum())
    if den <= 0.0:
        return math.nan
    return float(1.0 - ((y - p) ** 2).sum() / den)


def mean_row_cosine(pred, target):
    pred = np.asarray(pred, dtype=float)
    target = np.asarray(target, dtype=float)
    dot = np.sum(pred * target, axis=1)
    den = np.linalg.norm(pred, axis=1) * np.linalg.norm(target, axis=1)
    ok = np.isfinite(dot) & np.isfinite(den) & (den > 0.0)
    if ok.sum() < 1:
        return math.nan
    return float(np.mean(dot[ok] / den[ok]))


def split_masks(frames, strategy, seed, purge_frames=0):
    frames = np.asarray(frames)
    unique = np.sort(np.unique(frames))
    if len(unique) < 5:
        z = np.zeros(len(frames), dtype=bool)
        return z.copy(), z.copy(), {
            "test_frames": 0,
            "purged_train_frames": 0,
            "cross_split_lag1_pairs": 0,
        }
    n_test = max(1, int(TEST_FRAME_FRACTION * len(unique)))
    rng = np.random.default_rng(seed)

    if strategy == "random":
        test_frames = set(rng.choice(unique, n_test, replace=False))
        purged = set()
        train_frames = set(unique) - test_frames
    elif strategy == "blocked":
        start = int(rng.integers(0, len(unique) - n_test + 1))
        stop = start + n_test
        test_frames = set(unique[start:stop])
        purge_lo = max(0, start - max(0, purge_frames))
        purge_hi = min(len(unique), stop + max(0, purge_frames))
        purged = set(unique[purge_lo:start]) | set(unique[stop:purge_hi])
        train_frames = set(unique) - test_frames - purged
    else:
        raise ValueError(strategy)

    train = np.array([f in train_frames for f in frames], dtype=bool)
    test = np.array([f in test_frames for f in frames], dtype=bool)
    frame_split = {
        f: "test" if f in test_frames else "train" if f in train_frames else "purged"
        for f in unique
    }
    cross = 0
    for a, b in zip(unique[:-1], unique[1:]):
        if {frame_split[a], frame_split[b]} == {"train", "test"}:
            cross += 1
    return train, test, {
        "test_frames": int(len(test_frames)),
        "purged_train_frames": int(len(purged)),
        "cross_split_lag1_pairs": int(cross),
    }


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
        out[m] = values[m] - values[mt].mean(axis=0, keepdims=True)
    return out


def atom_means(values, atoms):
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values.reshape(-1, 1)
    labels = np.sort(np.unique(atoms))
    means = []
    counts = []
    for atom in labels:
        m = atoms == atom
        means.append(values[m].mean(axis=0))
        counts.append(int(m.sum()))
    return labels, np.asarray(means), np.asarray(counts, dtype=float)


def variance_shares(target, atoms):
    y = np.asarray(target, dtype=float)
    if y.ndim == 1:
        y = y.reshape(-1, 1)
    grand = y.mean(axis=0, keepdims=True)
    total = float(((y - grand) ** 2).sum())
    between = 0.0
    within = 0.0
    for atom in np.unique(atoms):
        m = atoms == atom
        mean = y[m].mean(axis=0, keepdims=True)
        between += float(m.sum() * ((mean - grand) ** 2).sum())
        within += float(((y[m] - mean) ** 2).sum())
    if total <= 0.0:
        return math.nan, math.nan
    return between / total, within / total


def effective_n_components(target, atoms, frames):
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


def fit_beta(x, y):
    return np.linalg.lstsq(x, y, rcond=None)[0]


def scalar_between_loao(x, y, atoms):
    labels, xbar, _ = atom_means(x, atoms)
    _, ybar, _ = atom_means(y, atoms)
    xbar_i = np.column_stack([np.ones(len(labels)), xbar])
    pred = np.full((len(labels), 1), np.nan)
    for i in range(len(labels)):
        train = np.ones(len(labels), dtype=bool)
        train[i] = False
        if train.sum() < xbar_i.shape[1]:
            continue
        beta = fit_beta(xbar_i[train], ybar[train])
        pred[i] = xbar_i[i] @ beta
    return r2_score(pred, ybar)


def scalar_within_split(x, y, atoms, frames, strategy, seed, purge_frames=0):
    train, test, split = split_masks(frames, strategy, seed, purge_frames)
    xw = centred_by_train_atom(x, atoms, train)
    yw = centred_by_train_atom(y, atoms, train)
    fit = train & np.isfinite(xw).all(axis=1) & np.isfinite(yw).all(axis=1)
    score = test & np.isfinite(xw).all(axis=1) & np.isfinite(yw).all(axis=1)
    if fit.sum() < xw.shape[1] or score.sum() < 2:
        return math.nan, split
    beta = fit_beta(xw[fit], yw[fit])
    pred = xw[score] @ beta
    return r2_score(pred, yw[score]), split


def tensor_gamma(x, y):
    den = float((x * x).sum())
    if den <= 0.0:
        return 0.0
    return float((x * y).sum() / den)


def tensor_metrics(pred, target):
    return {
        "R2": r2_score(pred, target),
        "absT2_r": corrcoef(np.linalg.norm(pred, axis=1), np.linalg.norm(target, axis=1)),
        "component_r": corrcoef(pred.reshape(-1), target.reshape(-1)),
        "mean_row_cosine": mean_row_cosine(pred, target),
    }


def tensor_between_loao(x, y, atoms):
    labels, xbar, _ = atom_means(x, atoms)
    _, ybar, _ = atom_means(y, atoms)
    pred = np.full_like(ybar, np.nan)
    for i in range(len(labels)):
        train = np.ones(len(labels), dtype=bool)
        train[i] = False
        if train.sum() < 2:
            continue
        gamma = tensor_gamma(xbar[train], ybar[train])
        pred[i] = gamma * xbar[i]
    return tensor_metrics(pred, ybar)


def tensor_within_split(x, y, atoms, frames, strategy, seed, purge_frames=0):
    train, test, split = split_masks(frames, strategy, seed, purge_frames)
    xw = centred_by_train_atom(x, atoms, train)
    yw = centred_by_train_atom(y, atoms, train)
    fit = train & np.isfinite(xw).all(axis=1) & np.isfinite(yw).all(axis=1)
    score = test & np.isfinite(xw).all(axis=1) & np.isfinite(yw).all(axis=1)
    if fit.sum() < 2 or score.sum() < 2:
        return {
            "R2": math.nan,
            "absT2_r": math.nan,
            "component_r": math.nan,
            "mean_row_cosine": math.nan,
        }, split
    gamma = tensor_gamma(xw[fit], yw[fit])
    return tensor_metrics(gamma * xw[score], yw[score]), split


def base_row(mechanism, target, stratum, treatment, n_atoms, rows, shares, neff, split):
    return {
        "mechanism": mechanism,
        "target": target,
        "stratum": stratum,
        "treatment_label": treatment,
        "rows": int(rows),
        "n_atoms_between": int(n_atoms),
        "thin_flag": "THIN" if n_atoms < THIN_ATOMS else "",
        "variance_share_between": shares[0],
        "variance_share_within": shares[1],
        **neff,
        "split_strategy": split.get("strategy", ""),
        "test_frames": split.get("test_frames", 0),
        "purged_train_frames": split.get("purged_train_frames", 0),
        "cross_split_lag1_pairs": split.get("cross_split_lag1_pairs", 0),
    }


def scalar_result(df, mechanism, stratum, feature_cols, treatment, strategy, seed, purge_frames=0):
    x = df[feature_cols].to_numpy(float)
    if "field_mag_sq" in feature_cols and "field_mag_sq" not in df:
        raise RuntimeError("field_mag_sq missing")
    y = df["dft_sigma_iso"].to_numpy(float).reshape(-1, 1)
    atoms = df["atom_index"].to_numpy(int)
    frames = df["h5_row"].to_numpy(int)
    shares = variance_shares(y, atoms)
    neff = effective_n_components(y, atoms, frames)
    between = scalar_between_loao(x, y, atoms)
    within, split = scalar_within_split(x, y, atoms, frames, strategy, seed, purge_frames)
    row = base_row(mechanism, "sigma_iso", stratum, treatment,
                   len(np.unique(atoms)), len(df), shares, neff,
                   {"strategy": strategy, **split})
    row.update({
        "between_LOAO_R2": between,
        "between_LOAO_absT2_r": math.nan,
        "between_LOAO_component_r": math.nan,
        "within_frame_R2": within,
        "within_frame_absT2_r": math.nan,
        "within_frame_component_r": math.nan,
        "within_frame_row_cosine": math.nan,
    })
    return row


def tensor_result(df, x, y, mechanism, stratum, treatment, strategy, seed, purge_frames=0):
    atoms = df["atom_index"].to_numpy(int)
    frames = df["h5_row"].to_numpy(int)
    shares = variance_shares(y, atoms)
    neff = effective_n_components(y, atoms, frames)
    between = tensor_between_loao(x, y, atoms)
    within, split = tensor_within_split(x, y, atoms, frames, strategy, seed, purge_frames)
    row = base_row(mechanism, "T2", stratum, treatment,
                   len(np.unique(atoms)), len(df), shares, neff,
                   {"strategy": strategy, **split})
    row.update({
        "between_LOAO_R2": between["R2"],
        "between_LOAO_absT2_r": between["absT2_r"],
        "between_LOAO_component_r": between["component_r"],
        "between_LOAO_row_cosine": between["mean_row_cosine"],
        "within_frame_R2": within["R2"],
        "within_frame_absT2_r": within["absT2_r"],
        "within_frame_component_r": within["component_r"],
        "within_frame_row_cosine": within["mean_row_cosine"],
    })
    return row


def load_broad(broad_dir, C, seed):
    broad_dir = Path(broad_dir)
    df = pd.read_csv(broad_dir / "broad_backbone_aggregated.csv")
    df["stratum"] = [stratum_of(v) for v in df["frame_variant"]]
    df["field_mag_sq"] = df["field_mag"].to_numpy(float) ** 2
    rows = []

    scalar_specs = [
        ("ring_current", ["ring_sum_dipolar"], "magnetic, solvation-insensitive"),
        ("bond_anisotropy_scalar", ["bond_sum_dipolar"], "magnetic, solvation-insensitive"),
        ("ff14sb_charge_field", ["field_z", "field_mag_sq"],
         "FF14SB vacuum-Coulomb proxy (treatment-mismatched to CPCM)"),
    ]
    for stratum in STRATA_ORDER:
        base = df[
            (df["stratum"] == stratum)
            & (df["dft_present"] == 1)
            & (df["frame_valid"] == 1)
            & np.isfinite(df["dft_sigma_iso"].to_numpy(float))
        ].copy()
        for mechanism, cols, treatment in scalar_specs:
            mask = np.isfinite(base[cols].to_numpy(float)).all(axis=1)
            d = base.loc[mask].reset_index(drop=True)
            if len(d) >= 4 and d["atom_index"].nunique() >= 2:
                rows.append(scalar_result(d, mechanism, stratum, cols, treatment,
                                          "random", seed))

    target = np.load(broad_dir / "broad_backbone_aggregated_target_local_T2.npy") @ C.T
    tensor_specs = [
        ("broad_total_T2", "literature_kernel_present",
         "broad_backbone_aggregated_literature_kernel_T2.npy",
         "mixed emitted sum (not a ranked composite score)"),
        ("ring_current_T2", "ring_literature_kernel_present",
         "broad_backbone_aggregated_ring_literature_kernel_T2.npy",
         "magnetic, solvation-insensitive"),
        ("bond_anisotropy_T2", "bond_literature_kernel_present",
         "broad_backbone_aggregated_bond_literature_kernel_T2.npy",
         "magnetic, solvation-insensitive"),
        ("ff14sb_charge_EFG_T2", "charge_literature_kernel_present",
         "broad_backbone_aggregated_charge_literature_kernel_T2.npy",
         "FF14SB vacuum-Coulomb proxy (treatment-mismatched to CPCM)"),
    ]
    for mechanism, present_col, npy, treatment in tensor_specs:
        feature = np.load(broad_dir / npy) @ C.T
        present = df[present_col].to_numpy(int) == 1
        finite = np.isfinite(feature).all(axis=1) & np.isfinite(target).all(axis=1)
        ok = (
            present
            & finite
            & (df["dft_present"].to_numpy(int) == 1)
            & (df["dft_local_frame_valid"].to_numpy(int) == 1)
            & df["stratum"].notna().to_numpy()
        )
        for stratum in STRATA_ORDER:
            idx = np.flatnonzero(ok & (df["stratum"].to_numpy() == stratum))
            if len(idx) < 4:
                continue
            d = df.iloc[idx].reset_index(drop=True)
            if d["atom_index"].nunique() < 2:
                continue
            rows.append(tensor_result(d, feature[idx], target[idx], mechanism,
                                      stratum, treatment, "random", seed))
    return rows


def load_buckingham(buckingham_dir, seed):
    buckingham_dir = Path(buckingham_dir)
    df = pd.read_csv(buckingham_dir / "buckingham_efield_aggregated.csv")
    df["stratum"] = [stratum_of(v) for v in df["frame_variant"]]
    df["E_mag_sq"] = df["E_mag"].to_numpy(float) ** 2
    rows = []
    for stratum in STRATA_ORDER:
        d = df[
            (df["stratum"] == stratum)
            & (df["dft_present"] == 1)
            & (df["apbs_efield_present"] == 1)
            & (df["frame_valid"] == 1)
        ].copy()
        cols = ["E_proj", "E_mag_sq"]
        mask = np.isfinite(d[cols + ["dft_sigma_iso"]].to_numpy(float)).all(axis=1)
        d = d.loc[mask].reset_index(drop=True)
        if len(d) >= 4 and d["atom_index"].nunique() >= 2:
            rows.append(scalar_result(
                d, "apbs_efield_buckingham", stratum, cols,
                "PB-continuum of FF charges (treatment-mismatched to CPCM)",
                "random", seed,
            ))
    return rows


def load_efg(efg_dir, C, seed, purge_frames):
    efg_dir = Path(efg_dir)
    df = pd.read_csv(efg_dir / "efg_aggregated.csv")
    df["stratum"] = [stratum_of(v) for v in df["frame_variant"]]
    feature = np.load(efg_dir / "efg_feature_T2.npy") @ C.T
    target = np.load(efg_dir / "efg_target_T2.npy") @ C.T
    ok = (
        (df["dft_present"].to_numpy(int) == 1)
        & (df["apbs_efg_present"].to_numpy(int) == 1)
        & np.isfinite(df["efg_T2_magnitude"].to_numpy(float))
        & np.isfinite(feature).all(axis=1)
        & np.isfinite(target).all(axis=1)
        & df["stratum"].notna().to_numpy()
    )
    rows = []
    treatment = "PB-continuum of FF charges (treatment-mismatched to CPCM)"
    for stratum in STRATA_ORDER:
        idx = np.flatnonzero(ok & (df["stratum"].to_numpy() == stratum))
        if len(idx) < 4:
            continue
        d = df.iloc[idx].reset_index(drop=True)
        if d["atom_index"].nunique() < 2:
            continue
        rows.append(tensor_result(d, feature[idx], target[idx],
                                  "apbs_efg_T2_old_random", stratum, treatment,
                                  "random", seed))
        rows.append(tensor_result(d, feature[idx], target[idx],
                                  "apbs_efg_T2", stratum, treatment,
                                  "blocked", seed, purge_frames))
    return rows


def write_markdown(df, path):
    cols = [
        "mechanism", "target", "stratum", "n_atoms_between",
        "variance_share_between", "variance_share_within",
        "within_N_eff_median", "median_lag1_rho",
        "between_LOAO_R2", "between_LOAO_absT2_r",
        "within_frame_R2", "within_frame_absT2_r",
        "split_strategy", "cross_split_lag1_pairs", "thin_flag",
    ]
    lines = ["# Rediscover variance decomposition", ""]
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
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

    rows = []
    rows.extend(load_broad(args.broad_dir, C, args.seed))
    rows.extend(load_buckingham(args.buckingham_dir, args.seed))
    rows.extend(load_efg(args.efg_dir, C, args.seed, args.efg_purge_frames))
    df = pd.DataFrame(rows)
    df = df.sort_values(["target", "mechanism", "stratum", "split_strategy"]).reset_index(drop=True)

    csv_path = out_dir / "variance_decomposition.csv"
    md_path = out_dir / "variance_decomposition.md"
    json_path = out_dir / "variance_decomposition_run.json"
    df.to_csv(csv_path, index=False)
    write_markdown(df, md_path)
    json_path.write_text(
        json.dumps(
            {
                "broad_dir": args.broad_dir,
                "buckingham_dir": args.buckingham_dir,
                "efg_dir": args.efg_dir,
                "seed": args.seed,
                "efg_purge_frames": args.efg_purge_frames,
                "change_of_basis_orthogonality_max": orth,
                "rows": int(len(df)),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"wrote {csv_path}")
    print(f"wrote {md_path}")
    print(f"wrote {json_path}")
    print("\nmechanism,target,stratum,between_share,within_share,N_eff_med,between,within,split,cross_lag1")
    for _, row in df.iterrows():
        between = row["between_LOAO_R2"]
        if row["target"] == "T2":
            between = row["between_LOAO_absT2_r"]
        within = row["within_frame_R2"]
        if row["target"] == "T2":
            within = row["within_frame_absT2_r"]
        print(
            f"{row['mechanism']},{row['target']},{row['stratum']},"
            f"{row['variance_share_between']:.4g},{row['variance_share_within']:.4g},"
            f"{row['within_N_eff_median']:.4g},{finite_fmt(between)},"
            f"{finite_fmt(within)},{row['split_strategy']},"
            f"{int(row['cross_split_lag1_pairs'])}"
        )


if __name__ == "__main__":
    main()
