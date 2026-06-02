#!/usr/bin/env python3
"""Per-category McConnell Delta-chi calibration from emitted CSV/NPY substrate.

This is a read-only Python consumer. It reads the already-emitted broad backbone
aggregate/source CSVs and the emitted local DFT T2 sidecar. It does not open H5,
run ORCA, rebuild tensors from geometry, or call the C++ emitter.

The broad source CSV carries one emitted McConnell tensor per source row plus its
`bond_category`. This script source-sums those emitted tensor components by
aggregate `row_id` and category, then removes the exact C++ Williamson-Asakura
scale to express the features as geometric tensor coefficients. No coordinate
formula is evaluated in Python.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


STRATA = ("HN", "N", "CA", "C", "O", "HA")
T2_COLS = [f"mc_lit_T2_local_{i}" for i in range(5)]
AGG_ALL_T2_COLS = [f"mc_lit_T2_local_{i}" for i in range(5)]
AGG_VALID_T2_COLS = [f"mc_lit_T2_local_valid_{i}" for i in range(5)]

AVOGADRO = 6.02214076e23
TENSOR_PREF = 1.0e24 / AVOGADRO
SCALAR_PREF = TENSOR_PREF / 3.0

CATEGORIES = {
    0: {"name": "PeptideCO", "q_wa": 2.41, "short": "CO"},
    1: {"name": "PeptideCN", "q_wa": -5.42, "short": "CN"},
    3: {"name": "SidechainCO", "q_wa": 2.41, "short": "SC"},
}
CATEGORY_IDS = tuple(CATEGORIES)
CATEGORY_NAMES = tuple(CATEGORIES[c]["name"] for c in CATEGORY_IDS)
AROMATIC_CATEGORY = 4

LITERATURE = {
    "PeptideCO": {
        "WA": 2.41,
        "Abraham aliphatic amide": 6.34,
        "ApSimon": 12.65,
        "Schneider": 14.45,
    },
    "PeptideCN": {
        "WA": -5.42,
        "Abraham aliphatic amide": -14.25,
        "ApSimon": -3.61,
        "Schneider": -7.23,
    },
    "SidechainCO": {
        "WA peptide-like": 2.41,
        "Abraham aliphatic amide": 6.34,
        "ApSimon carbonyl": 12.65,
        "Schneider carbonyl": 14.45,
    },
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out-dir",
        default="/tmp/rediscover-mc-lit-fresh",
        help="directory containing broad_backbone_aggregated.csv and broad_backbone_sources.csv",
    )
    ap.add_argument(
        "--decirc-csv",
        default="/tmp/rediscover-mc-lit-decirc/mcconnell_literature_decirc.csv",
        help="single-gamma literature-scaled CSV used for the Step-1 implied Delta-chi table",
    )
    ap.add_argument(
        "--artifact-dir",
        default="/tmp/rediscover-mc-dchi-calibration",
        help="directory for generated CSV/JSON artifacts",
    )
    ap.add_argument(
        "--report-md",
        default="src/rediscover/analysis/MCCONNELL_DCHI_CALIBRATION.md",
        help="markdown report path",
    )
    ap.add_argument(
        "--chunksize",
        type=int,
        default=500_000,
        help="source CSV rows per chunk while source-summing emitted tensors",
    )
    ap.add_argument(
        "--thin-atom-neff",
        type=float,
        default=5.0,
        help="flag strata with atom signal participation below this value",
    )
    ap.add_argument(
        "--mc-source-mode",
        choices=("valid", "all"),
        default="valid",
        help="calibrate against producer-valid McConnell source rows, or legacy all-source rows",
    )
    return ap.parse_args()


def require_columns(df: pd.DataFrame, cols: Iterable[str], label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise RuntimeError(f"{label} missing required columns: {missing}")


def finite_rows(*arrays: np.ndarray) -> np.ndarray:
    mask = None
    for arr in arrays:
        a = np.asarray(arr, dtype=float)
        ok = np.isfinite(a).all(axis=tuple(range(1, a.ndim))) if a.ndim > 1 else np.isfinite(a)
        mask = ok if mask is None else (mask & ok)
    if mask is None:
        raise RuntimeError("finite_rows requires at least one array")
    return mask


def corrcoef(a: np.ndarray, b: np.ndarray) -> float:
    aa = np.asarray(a, dtype=float).reshape(-1)
    bb = np.asarray(b, dtype=float).reshape(-1)
    ok = np.isfinite(aa) & np.isfinite(bb)
    if ok.sum() < 3:
        return math.nan
    xa = aa[ok] - aa[ok].mean()
    xb = bb[ok] - bb[ok].mean()
    den = math.sqrt(float(np.dot(xa, xa) * np.dot(xb, xb)))
    if den <= 0.0:
        return math.nan
    return float(np.dot(xa, xb) / den)


def row_stratum(name: object) -> str | None:
    atom = str(name).strip().upper()
    if atom == "N":
        return "N"
    if atom == "CA":
        return "CA"
    if atom == "C":
        return "C"
    if atom == "O":
        return "O"
    if atom in {"H", "HN"}:
        return "HN"
    if atom.startswith("HA"):
        return "HA"
    return None


def literature_nearest(category: str, q: float) -> tuple[str, float]:
    if not np.isfinite(q):
        return "NA", math.nan
    vals = LITERATURE[category]
    name, val = min(vals.items(), key=lambda kv: abs(q - kv[1]))
    return name, val


def literature_position(category: str, q: float) -> str:
    if not np.isfinite(q):
        return "NA"
    vals = list(LITERATURE[category].values())
    lo = min(vals)
    hi = max(vals)
    if category.endswith("CO") and q < 0.0:
        return "sign-flip"
    if category == "PeptideCN" and q > 0.0:
        return "sign-flip"
    if lo <= q <= hi:
        return "inside spread"
    if q < lo:
        return "below spread"
    return "above spread"


def fmt(x: object, digits: int = 3) -> str:
    if x is None:
        return "NA"
    if isinstance(x, str):
        return x
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(v):
        return "NA"
    av = abs(v)
    if av != 0.0 and (av < 1.0e-4 or av >= 1.0e4):
        return f"{v:.{digits}e}"
    return f"{v:.{digits}f}"


def q_pair(q: object, se: object) -> str:
    return f"{fmt(q)} +/- {fmt(se)}"


def source_sum_weighted_tensors(
    source_path: Path,
    n_rows: int,
    chunksize: int,
    source_mode: str,
) -> tuple[np.ndarray, dict[str, object]]:
    usecols = ["row_id", "bond_category", "mc_lit_kernel_present", *T2_COLS]
    if source_mode == "valid":
        usecols.append("mc_source_is_self_or_bonded")
    weighted = np.zeros((n_rows, len(CATEGORY_IDS), 5), dtype=float)
    category_to_col = {cat: i for i, cat in enumerate(CATEGORY_IDS)}
    counts = {str(cat): 0 for cat in CATEGORY_IDS}
    filtered_counts = {str(cat): 0 for cat in CATEGORY_IDS}
    aromatic_rows = 0
    aromatic_abs_sum = 0.0
    total_rows = 0
    present_rows = 0
    filtered_present_rows = 0

    for chunk in pd.read_csv(source_path, usecols=usecols, chunksize=chunksize):
        require_columns(chunk, usecols, "broad_backbone_sources.csv")
        total_rows += len(chunk)
        present = chunk["mc_lit_kernel_present"].to_numpy(int) == 1
        present_rows += int(present.sum())
        if source_mode == "valid":
            producer_valid = chunk["mc_source_is_self_or_bonded"].to_numpy(int) == 0
        else:
            producer_valid = np.ones(len(chunk), dtype=bool)
        filtered_present_rows += int((present & ~producer_valid).sum())

        aromatic = present & (chunk["bond_category"].to_numpy(int) == AROMATIC_CATEGORY)
        aromatic_rows += int(aromatic.sum())
        if aromatic.any():
            vals = chunk.loc[aromatic, T2_COLS].to_numpy(float)
            aromatic_abs_sum += float(np.abs(vals).sum())

        for cat, col in category_to_col.items():
            cat_mask = present & (chunk["bond_category"].to_numpy(int) == cat)
            if source_mode == "valid":
                filtered_counts[str(cat)] += int((cat_mask & ~producer_valid).sum())
            mask = cat_mask & producer_valid
            if not mask.any():
                continue
            rows = chunk.loc[mask, "row_id"].to_numpy(int)
            if rows.min(initial=0) < 0 or rows.max(initial=0) >= n_rows:
                raise RuntimeError(f"source row_id outside aggregate range 0..{n_rows - 1}")
            vals = chunk.loc[mask, T2_COLS].to_numpy(float)
            counts[str(cat)] += int(mask.sum())
            for k in range(5):
                np.add.at(weighted[:, col, k], rows, vals[:, k])

    audit = {
        "source_rows": total_rows,
        "present_source_rows": present_rows,
        "source_mode": source_mode,
        "filtered_present_source_rows": filtered_present_rows,
        "category_source_rows": counts,
        "category_filtered_source_rows": filtered_counts,
        "aromatic_present_rows_excluded": aromatic_rows,
        "aromatic_abs_mc_lit_T2_sum": aromatic_abs_sum,
    }
    return weighted, audit


def deweight_to_geometry(weighted: np.ndarray) -> tuple[np.ndarray, dict[str, float]]:
    geom = np.empty_like(weighted)
    scales: dict[str, float] = {}
    for i, cat in enumerate(CATEGORY_IDS):
        info = CATEGORIES[cat]
        q = float(info["q_wa"])
        scale = -TENSOR_PREF * q / 3.0
        if abs(scale) <= 0.0:
            raise RuntimeError(f"cannot deweight zero-scale category {info['name']}")
        geom[:, i, :] = weighted[:, i, :] / scale
        scales[info["name"]] = scale
    return geom, scales


def center_by_atom(
    x: np.ndarray,
    y: np.ndarray,
    atoms: np.ndarray,
    valid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    out_x = np.full_like(x, np.nan, dtype=float)
    out_y = np.full_like(y, np.nan, dtype=float)
    for atom in np.sort(np.unique(atoms[valid])):
        m = valid & (atoms == atom)
        if m.sum() == 0:
            continue
        out_x[m] = x[m] - x[m].mean(axis=0, keepdims=True)
        out_y[m] = y[m] - y[m].mean(axis=0, keepdims=True)
    return out_x, out_y


def atom_means(
    x: np.ndarray,
    y: np.ndarray,
    atoms: np.ndarray,
    valid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    labels: list[int] = []
    xm: list[np.ndarray] = []
    ym: list[np.ndarray] = []
    base = valid & finite_rows(x, y)
    for atom in np.sort(np.unique(atoms[base])):
        m = base & (atoms == atom)
        if m.sum() == 0:
            continue
        labels.append(int(atom))
        xm.append(x[m].mean(axis=0))
        ym.append(y[m].mean(axis=0))
    if not labels:
        return np.array([], dtype=int), np.empty((0, x.shape[1], x.shape[2])), np.empty((0, y.shape[1]))
    return np.asarray(labels, dtype=int), np.stack(xm), np.vstack(ym)


def obs_matrix(x: np.ndarray, y: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ok_rows = mask & finite_rows(x, y)
    if ok_rows.sum() == 0:
        return np.empty((0, x.shape[1])), np.empty((0,))
    xx = x[ok_rows].transpose(0, 2, 1).reshape(-1, x.shape[1])
    yy = y[ok_rows].reshape(-1)
    ok = np.isfinite(xx).all(axis=1) & np.isfinite(yy)
    return xx[ok], yy[ok]


def fit_beta_no_intercept(x: np.ndarray, y: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, int]:
    xx, yy = obs_matrix(x, y, mask)
    if xx.shape[0] < x.shape[1]:
        return np.full(x.shape[1], math.nan), 0
    beta, _, rank, _ = np.linalg.lstsq(xx, yy, rcond=None)
    return beta.astype(float), int(rank)


def fit_beta_with_intercept(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    ok = finite_rows(x, y)
    if ok.sum() < x.shape[1]:
        return np.full(x.shape[1], math.nan), np.full(y.shape[1], math.nan), 0
    xx0 = x[ok]
    yy0 = y[ok]
    x_mean = xx0.mean(axis=0, keepdims=True)
    y_mean = yy0.mean(axis=0, keepdims=True)
    xx = (xx0 - x_mean).transpose(0, 2, 1).reshape(-1, x.shape[1])
    yy = (yy0 - y_mean).reshape(-1)
    finite = np.isfinite(xx).all(axis=1) & np.isfinite(yy)
    if finite.sum() < x.shape[1]:
        return np.full(x.shape[1], math.nan), np.full(y.shape[1], math.nan), 0
    beta, _, rank, _ = np.linalg.lstsq(xx[finite], yy[finite], rcond=None)
    intercept = y_mean.reshape(-1) - np.einsum("c,ck->k", beta, x_mean.reshape(x.shape[1], x.shape[2]))
    return beta.astype(float), intercept.astype(float), int(rank)


def predict(beta: np.ndarray, x: np.ndarray, intercept: np.ndarray | None = None) -> np.ndarray:
    pred = np.einsum("c,nck->nk", beta, x)
    if intercept is not None:
        pred = pred + intercept.reshape(1, -1)
    return pred


def r2_score(pred: np.ndarray, target: np.ndarray, mask: np.ndarray | None = None) -> float:
    p = np.asarray(pred, dtype=float)
    y = np.asarray(target, dtype=float)
    ok = finite_rows(p, y)
    if mask is not None:
        ok = ok & mask
    if ok.sum() < 3:
        return math.nan
    pp = p[ok]
    yy = y[ok]
    den = float(((yy - yy.mean(axis=0, keepdims=True)) ** 2).sum())
    if den <= 0.0:
        return math.nan
    return float(1.0 - ((yy - pp) ** 2).sum() / den)


def metric_values(pred: np.ndarray, target: np.ndarray, mask: np.ndarray | None = None) -> dict[str, float]:
    p = np.asarray(pred, dtype=float)
    y = np.asarray(target, dtype=float)
    ok = finite_rows(p, y)
    if mask is not None:
        ok = ok & mask
    if ok.sum() < 3:
        return {"component_r": math.nan, "absT2_r": math.nan, "R2": math.nan}
    return {
        "component_r": corrcoef(p[ok].reshape(-1), y[ok].reshape(-1)),
        "absT2_r": corrcoef(np.linalg.norm(p[ok], axis=1), np.linalg.norm(y[ok], axis=1)),
        "R2": r2_score(p, y, ok),
    }


def jackknife_se(vals: Iterable[np.ndarray]) -> np.ndarray:
    arr = np.asarray(list(vals), dtype=float)
    if arr.ndim != 2 or arr.shape[0] < 2:
        width = arr.shape[1] if arr.ndim == 2 else len(CATEGORY_IDS)
        return np.full(width, math.nan)
    good = np.isfinite(arr).all(axis=1)
    arr = arr[good]
    if arr.shape[0] < 2:
        return np.full(arr.shape[1] if arr.ndim == 2 else len(CATEGORY_IDS), math.nan)
    mean = arr.mean(axis=0, keepdims=True)
    return np.sqrt((arr.shape[0] - 1) / arr.shape[0] * ((arr - mean) ** 2).sum(axis=0))


def jk_se_within(x: np.ndarray, y: np.ndarray, atoms: np.ndarray, mask: np.ndarray) -> np.ndarray:
    vals = []
    labels = np.sort(np.unique(atoms[mask]))
    if labels.size <= len(CATEGORY_IDS):
        return np.full(len(CATEGORY_IDS), math.nan)
    for atom in labels:
        beta, _ = fit_beta_no_intercept(x, y, mask & (atoms != atom))
        vals.append(beta)
    return jackknife_se(vals)


def jk_se_between(x: np.ndarray, y: np.ndarray, labels: np.ndarray) -> np.ndarray:
    vals = []
    if labels.size <= len(CATEGORY_IDS):
        return np.full(len(CATEGORY_IDS), math.nan)
    for atom in labels:
        keep = labels != atom
        beta, _, _ = fit_beta_with_intercept(x[keep], y[keep])
        vals.append(beta)
    return jackknife_se(vals)


def loao_within_r2(x: np.ndarray, y: np.ndarray, atoms: np.ndarray, mask: np.ndarray) -> float:
    pred = np.full_like(y, np.nan, dtype=float)
    labels = np.sort(np.unique(atoms[mask]))
    if labels.size <= len(CATEGORY_IDS):
        return math.nan
    for atom in labels:
        train = mask & (atoms != atom)
        score = mask & (atoms == atom)
        beta, _ = fit_beta_no_intercept(x, y, train)
        if np.isfinite(beta).all():
            pred[score] = predict(beta, x[score])
    return r2_score(pred, y, mask)


def loao_between_r2(x: np.ndarray, y: np.ndarray, labels: np.ndarray) -> float:
    pred = np.full_like(y, np.nan, dtype=float)
    if labels.size <= len(CATEGORY_IDS):
        return math.nan
    for i, atom in enumerate(labels):
        keep = labels != atom
        beta, intercept, _ = fit_beta_with_intercept(x[keep], y[keep])
        if np.isfinite(beta).all() and np.isfinite(intercept).all():
            pred[i] = predict(beta, x[i : i + 1], intercept)[0]
    return r2_score(pred, y)


def signal_effective_atoms(x: np.ndarray, atoms: np.ndarray, mask: np.ndarray) -> tuple[float, int, int]:
    flat = x.reshape(x.shape[0], -1)
    weights = []
    active = 0
    for atom in np.sort(np.unique(atoms[mask])):
        m = mask & (atoms == atom)
        if m.sum() < 2:
            continue
        xc = flat[m] - flat[m].mean(axis=0, keepdims=True)
        w = float((xc * xc).sum())
        if w > 0.0:
            active += 1
            weights.append(w)
    w_arr = np.asarray(weights, dtype=float)
    if w_arr.size == 0 or float((w_arr * w_arr).sum()) <= 0.0:
        return math.nan, 0, active
    participation = float((w_arr.sum() * w_arr.sum()) / (w_arr * w_arr).sum())
    ordered = np.sort(w_arr)[::-1]
    cum = np.cumsum(ordered)
    top90 = int(np.searchsorted(cum / cum[-1], 0.90) + 1)
    return participation, top90, active


def atom_mean_signal_neff(x: np.ndarray) -> float:
    flat = x.reshape(x.shape[0], -1)
    if flat.shape[0] < 2:
        return math.nan
    xc = flat - flat.mean(axis=0, keepdims=True)
    w = np.sum(xc * xc, axis=1)
    w = w[w > 0.0]
    if w.size == 0 or float((w * w).sum()) <= 0.0:
        return math.nan
    return float((w.sum() * w.sum()) / (w * w).sum())


def ar1_effective_frames(x: np.ndarray, atoms: np.ndarray, frames: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
    flat = x.reshape(x.shape[0], -1)
    mag = np.linalg.norm(flat, axis=1)
    total = 0.0
    rhos = []
    for atom in np.sort(np.unique(atoms[mask])):
        m = mask & (atoms == atom)
        order = np.argsort(frames[m])
        vals = mag[m][order]
        vals = vals[np.isfinite(vals)]
        n = vals.size
        if n < 3:
            continue
        yd = vals - vals.mean()
        den = float(np.dot(yd, yd))
        if den <= 0.0:
            continue
        rho = float(np.dot(yd[:-1], yd[1:]) / den)
        rho = max(min(rho, 0.999), -0.999)
        neff = n * (1.0 - rho) / (1.0 + rho)
        neff = min(max(neff, 1.0), float(n))
        total += neff
        rhos.append(rho)
    return total, float(np.median(rhos)) if rhos else math.nan


def thin_label(neff: float, threshold: float) -> str:
    if not np.isfinite(neff):
        return "thin"
    return "thin" if neff < threshold else "ok"


def beta_to_q(beta: np.ndarray) -> np.ndarray:
    return -np.asarray(beta, dtype=float) / SCALAR_PREF


def summarize_stratum(
    rows: list[dict[str, object]],
    stratum: str,
    x: np.ndarray,
    y: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    valid: np.ndarray,
    thin_threshold: float,
) -> None:
    valid = valid & finite_rows(x, y)
    atom_neff, top90, active_atoms = signal_effective_atoms(x, atoms, valid)
    ar1_neff, ar1_lag1 = ar1_effective_frames(x, atoms, frames, valid)
    atoms_total = int(np.unique(atoms[valid]).size)

    xw, yw = center_by_atom(x, y, atoms, valid)
    score = valid & finite_rows(xw, yw)
    beta, rank = fit_beta_no_intercept(xw, yw, score)
    beta_se = jk_se_within(xw, yw, atoms, score)
    q = beta_to_q(beta)
    q_se = np.abs(beta_to_q(beta_se))
    pred = predict(beta, xw)
    metrics = metric_values(pred, yw, score)
    rows.append(
        make_row(
            stratum,
            "within",
            int(score.sum()),
            atoms_total,
            active_atoms,
            atom_neff,
            top90,
            math.nan,
            ar1_neff,
            ar1_lag1,
            thin_label(atom_neff, thin_threshold),
            beta,
            beta_se,
            q,
            q_se,
            rank,
            metrics,
            loao_within_r2(xw, yw, atoms, score),
            math.nan,
        )
    )

    labels, xm, ym = atom_means(x, y, atoms, valid)
    beta_b, intercept, rank_b = fit_beta_with_intercept(xm, ym)
    beta_b_se = jk_se_between(xm, ym, labels)
    q_b = beta_to_q(beta_b)
    q_b_se = np.abs(beta_to_q(beta_b_se))
    pred_b = predict(beta_b, xm, intercept)
    metrics_b = metric_values(pred_b, ym)
    rows.append(
        make_row(
            stratum,
            "between",
            int(labels.size),
            atoms_total,
            active_atoms,
            atom_neff,
            top90,
            atom_mean_signal_neff(xm),
            ar1_neff,
            ar1_lag1,
            thin_label(atom_neff, thin_threshold),
            beta_b,
            beta_b_se,
            q_b,
            q_b_se,
            rank_b,
            metrics_b,
            loao_between_r2(xm, ym, labels),
            float(np.linalg.norm(intercept)) if np.isfinite(intercept).all() else math.nan,
        )
    )


def make_row(
    stratum: str,
    axis: str,
    rows_scored: int,
    atoms_total: int,
    atoms_active: int,
    atom_signal_neff: float,
    top90: int,
    atom_mean_neff: float,
    ar1_neff: float,
    ar1_lag1: float,
    thin_flag: str,
    beta: np.ndarray,
    beta_se: np.ndarray,
    q: np.ndarray,
    q_se: np.ndarray,
    rank: int,
    metrics: dict[str, float],
    loao_r2: float,
    intercept_norm: float,
) -> dict[str, object]:
    row: dict[str, object] = {
        "stratum": stratum,
        "axis": axis,
        "rows_scored": rows_scored,
        "atoms_total": atoms_total,
        "atoms_active": atoms_active,
        "atom_signal_neff": atom_signal_neff,
        "atom_signal_top90_count": top90,
        "atom_mean_signal_neff": atom_mean_neff,
        "ar1_frame_neff": ar1_neff,
        "kernel_median_lag1": ar1_lag1,
        "thin_flag": thin_flag,
        "rank": rank,
        "component_r": metrics["component_r"],
        "absT2_r": metrics["absT2_r"],
        "R2": metrics["R2"],
        "LOAO_R2": loao_r2,
        "intercept_norm": intercept_norm,
    }
    for i, name in enumerate(CATEGORY_NAMES):
        row[f"beta_{name}"] = float(beta[i])
        row[f"beta_{name}_se"] = float(beta_se[i])
        row[f"q_{name}"] = float(q[i])
        row[f"q_{name}_se"] = float(q_se[i])
        nearest, nearest_val = literature_nearest(name, float(q[i]))
        row[f"nearest_{name}"] = nearest
        row[f"nearest_{name}_value"] = nearest_val
        row[f"spread_{name}"] = literature_position(name, float(q[i]))
    co = float(q[0])
    cn = float(q[1])
    row["CO_CN_ratio"] = co / cn if np.isfinite(co) and np.isfinite(cn) and abs(cn) > 0.0 else math.nan
    return row


def implied_rows(decirc_csv: Path) -> list[dict[str, object]]:
    df = pd.read_csv(decirc_csv)
    require_columns(df, ["stratum", "target", "axis", "gamma_lit", "gamma_lit_se"], "decirc CSV")
    out: list[dict[str, object]] = []
    dft = df[df["target"] == "T2"].copy()
    for _, row in dft.iterrows():
        gamma = float(row["gamma_lit"])
        se = float(row["gamma_lit_se"])
        qco = gamma * 2.41
        qcn = gamma * -5.42
        qco_se = abs(se * 2.41)
        qcn_se = abs(se * -5.42)
        co_near, co_val = literature_nearest("PeptideCO", qco)
        cn_near, cn_val = literature_nearest("PeptideCN", qcn)
        out.append(
            {
                "stratum": row["stratum"],
                "axis": row["axis"],
                "gamma_lit": gamma,
                "gamma_lit_se": se,
                "implied_CO": qco,
                "implied_CO_se": qco_se,
                "implied_CN": qcn,
                "implied_CN_se": qcn_se,
                "CO_nearest": co_near,
                "CO_nearest_value": co_val,
                "CO_position": literature_position("PeptideCO", qco),
                "CN_nearest": cn_near,
                "CN_nearest_value": cn_val,
                "CN_position": literature_position("PeptideCN", qcn),
            }
        )
    return out


def write_csv_rows(rows: list[dict[str, object]], path: Path) -> None:
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_report(
    path: Path,
    out_dir: Path,
    decirc_csv: Path,
    artifact_dir: Path,
    implied: list[dict[str, object]],
    fit_rows: list[dict[str, object]],
    audit: dict[str, object],
    source_mode: str,
) -> None:
    wa_ratio = 2.41 / -5.42
    lines = [
        "# McConnell Delta-Chi Calibration",
        "",
        "Status: ORCA-free layer-2 calibration from emitted substrate only. The calibrated Delta-chi values below are provisional DFT-calibrated coefficients with within-protein confidence, not a de-circularisation claim.",
        "",
        f"Input out-dir: `{out_dir}`",
        f"Single-gamma CSV: `{decirc_csv}`",
        f"Artifacts: `{artifact_dir}`",
        f"McConnell source mode: `{source_mode}`",
        "",
        "Relation used: the C++ McConnell source scales the local geometric tensor as `T2_ppm = -(1e24/N_A) * q / 3 * T2_geom`. With `q = Delta-chi/(10^-6 cm^3 mol^-1)`, `1e24/N_A = 1.660539067`, so `beta_geom = -0.553513022 * q` and `q_cal = -beta_geom / 0.553513022`.",
        "",
        "The broad source CSV does not carry aggregate per-category unweighted columns. It does carry emitted source tensor rows with `bond_category`; this report source-sums those emitted tensor components by `row_id` and category, optionally filters with the emitted C++ producer-valid flag, then removes the exact C++ WA scalar. It does not evaluate a distance, angle, tensor projection, H5 read, ORCA job, or emitter.",
        "",
        "## Calibrated Delta-Chi Lead",
        "",
        "| stratum | axis | q_CO | q_CN | q_sidechain_CO | CO/CN | absT2 r | R2 | LOAO R2 | N_eff atom | nearest CO / CN / SC |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in fit_rows:
        nearest = (
            f"{row['nearest_PeptideCO']} / "
            f"{row['nearest_PeptideCN']} / "
            f"{row['nearest_SidechainCO']}"
        )
        lines.append(
            "| {stratum} | {axis} | {co} | {cn} | {sc} | {ratio} | {mag} | {r2} | {loao} | {neff} | {nearest} |".format(
                stratum=row["stratum"],
                axis=row["axis"],
                co=q_pair(row["q_PeptideCO"], row["q_PeptideCO_se"]),
                cn=q_pair(row["q_PeptideCN"], row["q_PeptideCN_se"]),
                sc=q_pair(row["q_SidechainCO"], row["q_SidechainCO_se"]),
                ratio=fmt(row["CO_CN_ratio"]),
                mag=fmt(row["absT2_r"]),
                r2=fmt(row["R2"]),
                loao=fmt(row["LOAO_R2"]),
                neff=fmt(row["atom_signal_neff"]),
                nearest=nearest,
            )
        )

    lines.extend(
        [
            "",
            f"WA assumes CO/CN = {wa_ratio:.3f}. The fitted CO/CN ratios vary by stratum and axis; rows with sign flips or very large SEs should be treated as diagnostics rather than transferable constants.",
            "",
            "## Fit Diagnostics",
            "",
            "| stratum | axis | rows | atoms | active | rank | component r | absT2 r | R2 | LOAO R2 | AR1 N_eff | intercept norm | thin |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in fit_rows:
        lines.append(
            "| {stratum} | {axis} | {rows} | {atoms} | {active} | {rank} | {comp} | {mag} | {r2} | {loao} | {ar1} | {intercept} | {thin} |".format(
                stratum=row["stratum"],
                axis=row["axis"],
                rows=row["rows_scored"],
                atoms=row["atoms_total"],
                active=row["atoms_active"],
                rank=row["rank"],
                comp=fmt(row["component_r"]),
                mag=fmt(row["absT2_r"]),
                r2=fmt(row["R2"]),
                loao=fmt(row["LOAO_R2"]),
                ar1=fmt(row["ar1_frame_neff"]),
                intercept=fmt(row["intercept_norm"]),
                thin=row["thin_flag"],
            )
        )

    lines.extend(
        [
            "",
            "## Step 1 Implied Readout",
            "",
            "Arithmetic from the prior single-gamma WA run: `q_CO = gamma_lit * 2.41`, `q_CN = gamma_lit * -5.42`.",
            "",
            "| stratum | axis | gamma_lit | implied CO | CO spread | implied CN | CN spread |",
            "|---|---:|---:|---:|---|---:|---|",
        ]
    )
    for row in implied:
        lines.append(
            "| {stratum} | {axis} | {gamma} | {co} | {copos}; nearest {conear} | {cn} | {cnpos}; nearest {cnnear} |".format(
                stratum=row["stratum"],
                axis=row["axis"],
                gamma=q_pair(row["gamma_lit"], row["gamma_lit_se"]),
                co=q_pair(row["implied_CO"], row["implied_CO_se"]),
                copos=row["CO_position"],
                conear=row["CO_nearest"],
                cn=q_pair(row["implied_CN"], row["implied_CN_se"]),
                cnpos=row["CN_position"],
                cnnear=row["CN_nearest"],
            )
        )

    lines.extend(
        [
            "",
            "Step-1 disagreement flags: N between implies CO about +6.3 and CN about -14.1, closest to Abraham aliphatic amide. C within/between sign-flips both categories (CO about -5, CN about +10 to +12), consistent with the carbonyl-C near-field/self-bond warning rather than a transferable far-field Delta-chi. O between is far above the literature spread and weakly predictive out-of-sample. HN/HA within are low-scale diagnostics.",
            "",
            "## Literature Spread Used",
            "",
            "| category | WA | Abraham aliphatic amide | ApSimon | Schneider |",
            "|---|---:|---:|---:|---:|",
            "| PeptideCO | 2.41 | 6.34 | 12.65 | 14.45 |",
            "| PeptideCN | -5.42 | -14.25 | -3.61 | -7.23 |",
            "| SidechainCO | 2.41 | 6.34 | 12.65 | 14.45 |",
            "",
            "## Self Audit",
            "",
            "- Aggregated source-sum vs emitted aggregate `{cols}`: RMS={rms}, max_abs={mx}. This checks the category source-sum is the emitted McConnell tensor, not a rebuilt one.".format(
                cols=audit["aggregate_t2_columns"],
                rms=fmt(audit["source_sum_vs_aggregate_rms"], 8),
                mx=fmt(audit["source_sum_vs_aggregate_max_abs"], 8),
            ),
            "- Producer-valid filter: mode={mode}, filtered_present_source_rows={rows}, category_filtered_source_rows={cats}.".format(
                mode=source_mode,
                rows=audit["filtered_present_source_rows"],
                cats=audit["category_filtered_source_rows"],
            ),
            "- Aromatic category excluded: rows={rows}, emitted abs T2 sum={abs_sum}. RING owns the pi current.".format(
                rows=audit["aromatic_present_rows_excluded"],
                abs_sum=fmt(audit["aromatic_abs_mc_lit_T2_sum"], 8),
            ),
            "- Manual conversion example: {stratum} {axis} PeptideCO beta={beta}; q=-beta/{pref}={q}; reported q={reported}; abs_diff={diff}.".format(
                stratum=audit["manual_stratum"],
                axis=audit["manual_axis"],
                beta=fmt(audit["manual_beta"], 8),
                pref=fmt(SCALAR_PREF, 8),
                q=fmt(audit["manual_q"], 8),
                reported=fmt(audit["manual_reported_q"], 8),
                diff=fmt(audit["manual_abs_diff"], 8),
            ),
            "- Python boundary: no ORCA, no `trajectory.h5`, and no coordinate/tensor physics recompute in Python.",
            "",
            "## Rerun",
            "",
            "`python3 src/rediscover/analysis/mcconnell_dchi_calibration.py --out-dir <new-broad-emitted-dir> --decirc-csv <new-single-gamma-csv>`",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    artifact_dir = Path(args.artifact_dir)
    agg_path = out_dir / "broad_backbone_aggregated.csv"
    source_path = out_dir / "broad_backbone_sources.csv"
    target_path = out_dir / "broad_backbone_aggregated_target_local_T2.npy"
    decirc_csv = Path(args.decirc_csv)

    if not agg_path.exists():
        raise RuntimeError(f"missing aggregate CSV: {agg_path}")
    if not source_path.exists():
        raise RuntimeError(f"missing source CSV: {source_path}")
    if not target_path.exists():
        raise RuntimeError(f"missing local target T2 sidecar: {target_path}")
    if not decirc_csv.exists():
        raise RuntimeError(f"missing Step-1 decirc CSV: {decirc_csv}")

    agg_t2_cols = AGG_VALID_T2_COLS if args.mc_source_mode == "valid" else AGG_ALL_T2_COLS
    agg = pd.read_csv(
        agg_path,
        usecols=[
            "row_id",
            "atom_index",
            "atom_name",
            "h5_row",
            "dft_present",
            "dft_local_frame_valid",
            *agg_t2_cols,
        ],
    )
    require_columns(
        agg,
        ["row_id", "atom_index", "atom_name", "h5_row", "dft_present", "dft_local_frame_valid", *agg_t2_cols],
        "broad_backbone_aggregated.csv",
    )
    if not np.array_equal(agg["row_id"].to_numpy(int), np.arange(len(agg))):
        raise RuntimeError("aggregate row_id is not contiguous 0..N-1; source grouping needs explicit remap")
    target = np.load(target_path)
    if target.shape != (len(agg), 5):
        raise RuntimeError(f"{target_path} shape {target.shape} != ({len(agg)}, 5)")

    weighted, source_audit = source_sum_weighted_tensors(
        source_path, len(agg), args.chunksize, args.mc_source_mode
    )
    geom, scales = deweight_to_geometry(weighted)

    summed_weighted = weighted.sum(axis=1)
    agg_weighted = agg[agg_t2_cols].to_numpy(float)
    ok = finite_rows(summed_weighted, agg_weighted)
    diff = summed_weighted[ok] - agg_weighted[ok]

    agg = agg.copy()
    agg["stratum2"] = [row_stratum(v) for v in agg["atom_name"]]
    atoms = agg["atom_index"].to_numpy(int)
    frames = agg["h5_row"].to_numpy(int)
    present = (
        agg["dft_present"].to_numpy(int).astype(bool)
        & agg["dft_local_frame_valid"].to_numpy(int).astype(bool)
    )

    fit_rows: list[dict[str, object]] = []
    for stratum in STRATA:
        sm = (agg["stratum2"] == stratum).to_numpy(bool)
        if not sm.any():
            continue
        summarize_stratum(
            fit_rows,
            stratum,
            geom,
            target.astype(float),
            atoms,
            frames,
            present & sm,
            args.thin_atom_neff,
        )

    implied = implied_rows(decirc_csv)
    fit_csv = artifact_dir / "mcconnell_dchi_calibration.csv"
    implied_csv = artifact_dir / "mcconnell_dchi_implied_readout.csv"
    audit_path = artifact_dir / "mcconnell_dchi_calibration_audit.json"
    write_csv_rows(fit_rows, fit_csv)
    write_csv_rows(implied, implied_csv)

    manual_row = next(r for r in fit_rows if r["stratum"] == "N" and r["axis"] == "between")
    manual_beta = float(manual_row["beta_PeptideCO"])
    manual_q = -manual_beta / SCALAR_PREF
    audit = {
        "out_dir": str(out_dir),
        "aggregate_csv": str(agg_path),
        "source_csv": str(source_path),
        "target_local_T2": str(target_path),
        "decirc_csv": str(decirc_csv),
        "artifact_fit_csv": str(fit_csv),
        "artifact_implied_csv": str(implied_csv),
        "source_mode": args.mc_source_mode,
        "aggregate_t2_columns": agg_t2_cols,
        "category_ids": {str(k): v["name"] for k, v in CATEGORIES.items()},
        "wa_q": {v["name"]: v["q_wa"] for v in CATEGORIES.values()},
        "deweight_scales": scales,
        "tensor_prefactor_before_one_third": TENSOR_PREF,
        "scalar_prefactor_after_one_third": SCALAR_PREF,
        "source_sum_vs_aggregate_rms": float(math.sqrt(float(np.mean(diff * diff)))) if diff.size else math.nan,
        "source_sum_vs_aggregate_max_abs": float(np.max(np.abs(diff))) if diff.size else math.nan,
        "manual_stratum": "N",
        "manual_axis": "between",
        "manual_beta": manual_beta,
        "manual_q": manual_q,
        "manual_reported_q": float(manual_row["q_PeptideCO"]),
        "manual_abs_diff": abs(manual_q - float(manual_row["q_PeptideCO"])),
        "no_orca": True,
        "no_h5_read": True,
        "no_coordinate_physics_recompute": True,
        "aromatic_excluded": True,
    }
    audit.update(source_audit)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    write_report(Path(args.report_md), out_dir, decirc_csv, artifact_dir, implied, fit_rows,
                 audit, args.mc_source_mode)
    print(f"wrote {fit_csv}")
    print(f"wrote {implied_csv}")
    print(f"wrote {audit_path}")
    print(f"wrote {args.report_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
