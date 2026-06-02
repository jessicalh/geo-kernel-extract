#!/usr/bin/env python3
"""McConnell literature-scaled de-circularisation from emitted mc_lit columns.

Read-only consumer for the rediscover McConnell Δχ emit. Inputs are the
already-emitted broad_backbone aggregate CSV columns:

  * mc_lit_T0
  * mc_lit_T2_local_0..4

and the emitted local DFT T2 sidecar. Python does not open H5, rebuild a kernel,
project tensors, or apply per-category Δχ. It only groups, de-means, correlates,
and reports scale diagnostics per backbone stratum.
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


ALL_T2_COLS = [f"mc_lit_T2_local_{i}" for i in range(5)]
VALID_T2_COLS = [f"mc_lit_T2_local_valid_{i}" for i in range(5)]
STRATA = ("HN", "N", "CA", "C", "O", "HA")
GAMMA_DEN_EPS = 1.0e-18


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out-dir",
        default="/tmp/rediscover-mc-lit/composed",
        help="directory containing broad_backbone_aggregated.csv and sidecars",
    )
    ap.add_argument(
        "--artifact-dir",
        default="/tmp/rediscover-mc-lit-decirc",
        help="directory for CSV/JSON analysis artifacts",
    )
    ap.add_argument(
        "--report-md",
        default="src/rediscover/analysis/MCCONNELL_LITERATURE_DECIRC.md",
        help="markdown report path",
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
        help="use producer-valid McConnell source sum columns, or legacy all-source columns",
    )
    return ap.parse_args()


def require_columns(df: pd.DataFrame, cols: Iterable[str], label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise RuntimeError(f"{label} missing required columns: {missing}")


def as_matrix(a: np.ndarray) -> np.ndarray:
    arr = np.asarray(a, dtype=float)
    if arr.ndim == 1:
        return arr.reshape(-1, 1)
    return arr


def finite_rows(*arrays: np.ndarray) -> np.ndarray:
    mask = None
    for arr in arrays:
        a = np.asarray(arr, dtype=float)
        ok = np.isfinite(a).all(axis=1) if a.ndim > 1 else np.isfinite(a)
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
    if den <= GAMMA_DEN_EPS:
        return math.nan
    return float(np.dot(xa, xb) / den)


def r2_score(pred: np.ndarray, target: np.ndarray) -> float:
    pred = as_matrix(pred)
    target = as_matrix(target)
    ok = finite_rows(pred, target)
    if ok.sum() < 3:
        return math.nan
    p = pred[ok]
    y = target[ok]
    den = float(((y - y.mean(axis=0, keepdims=True)) ** 2).sum())
    if den <= 0.0:
        return math.nan
    return float(1.0 - ((y - p) ** 2).sum() / den)


def center_by_atom(
    x: np.ndarray,
    y: np.ndarray,
    atoms: np.ndarray,
    valid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    x = as_matrix(x)
    y = as_matrix(y)
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
    x = as_matrix(x)
    y = as_matrix(y)
    labels: list[int] = []
    xm: list[np.ndarray] = []
    ym: list[np.ndarray] = []
    for atom in np.sort(np.unique(atoms[valid])):
        m = valid & (atoms == atom) & finite_rows(x, y)
        if m.sum() == 0:
            continue
        labels.append(int(atom))
        xm.append(x[m].mean(axis=0))
        ym.append(y[m].mean(axis=0))
    if not labels:
        return np.array([], dtype=int), np.empty((0, x.shape[1])), np.empty((0, y.shape[1]))
    return np.asarray(labels, dtype=int), np.vstack(xm), np.vstack(ym)


def gamma_no_intercept(x: np.ndarray, y: np.ndarray, mask: np.ndarray) -> float:
    x = as_matrix(x)
    y = as_matrix(y)
    score = np.asarray(mask, dtype=bool) & finite_rows(x, y)
    if score.sum() < 3:
        return math.nan
    den = float((x[score] * x[score]).sum())
    if den <= 0.0:
        return math.nan
    return float((x[score] * y[score]).sum() / den)


def gamma_with_intercept(x: np.ndarray, y: np.ndarray) -> tuple[float, np.ndarray]:
    x = as_matrix(x)
    y = as_matrix(y)
    ok = finite_rows(x, y)
    if ok.sum() < 3:
        return math.nan, np.full(y.shape[1], math.nan)
    xx = x[ok]
    yy = y[ok]
    xm = xx.mean(axis=0, keepdims=True)
    ym = yy.mean(axis=0, keepdims=True)
    xc = xx - xm
    yc = yy - ym
    den = float((xc * xc).sum())
    if den <= GAMMA_DEN_EPS:
        return math.nan, np.full(y.shape[1], math.nan)
    gamma = float((xc * yc).sum() / den)
    intercept = (ym - gamma * xm).reshape(-1)
    return gamma, intercept


def jackknife_se(vals: Iterable[float]) -> float:
    arr = np.asarray(list(vals), dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < 2:
        return math.nan
    mean = float(arr.mean())
    return float(math.sqrt((arr.size - 1) / arr.size * ((arr - mean) ** 2).sum()))


def jk_se_within(x: np.ndarray, y: np.ndarray, atoms: np.ndarray, mask: np.ndarray) -> float:
    vals = []
    labels = np.sort(np.unique(atoms[mask]))
    if labels.size <= 3:
        return math.nan
    for atom in labels:
        keep = mask & (atoms != atom)
        vals.append(gamma_no_intercept(x, y, keep))
    return jackknife_se(vals)


def jk_se_between(x: np.ndarray, y: np.ndarray, atom_labels: np.ndarray) -> tuple[float, float]:
    gammas = []
    intercept0 = []
    if atom_labels.size <= 3:
        return math.nan, math.nan
    for atom in atom_labels:
        keep = atom_labels != atom
        gamma, intercept = gamma_with_intercept(x[keep], y[keep])
        gammas.append(gamma)
        intercept0.append(float(intercept[0]) if intercept.size else math.nan)
    return jackknife_se(gammas), jackknife_se(intercept0)


def loao_within_r2(x: np.ndarray, y: np.ndarray, atoms: np.ndarray, mask: np.ndarray) -> float:
    x = as_matrix(x)
    y = as_matrix(y)
    pred = np.full_like(y, np.nan, dtype=float)
    labels = np.sort(np.unique(atoms[mask]))
    if labels.size <= 3:
        return math.nan
    for atom in labels:
        train = mask & (atoms != atom)
        score = mask & (atoms == atom)
        gamma = gamma_no_intercept(x, y, train)
        if np.isfinite(gamma):
            pred[score] = gamma * x[score]
    return r2_score(pred[mask], y[mask])


def loao_between_r2(x: np.ndarray, y: np.ndarray, atom_labels: np.ndarray) -> float:
    x = as_matrix(x)
    y = as_matrix(y)
    pred = np.full_like(y, np.nan, dtype=float)
    if atom_labels.size <= 3:
        return math.nan
    for i, atom in enumerate(atom_labels):
        train = atom_labels != atom
        gamma, intercept = gamma_with_intercept(x[train], y[train])
        if np.isfinite(gamma):
            pred[i] = intercept.reshape(1, -1) + gamma * x[i]
    return r2_score(pred, y)


def signal_effective_atoms(x: np.ndarray, atoms: np.ndarray, mask: np.ndarray) -> tuple[float, int, int]:
    x = as_matrix(x)
    weights = []
    active = 0
    for atom in np.sort(np.unique(atoms[mask])):
        m = mask & (atoms == atom)
        if m.sum() < 2:
            continue
        xc = x[m] - x[m].mean(axis=0, keepdims=True)
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
    x = as_matrix(x)
    if x.shape[0] < 2:
        return math.nan
    xc = x - x.mean(axis=0, keepdims=True)
    w = np.sum(xc * xc, axis=1)
    w = w[w > 0.0]
    if w.size == 0 or float((w * w).sum()) <= 0.0:
        return math.nan
    return float((w.sum() * w.sum()) / (w * w).sum())


def ar1_effective_frames(x: np.ndarray, atoms: np.ndarray, frames: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
    x = as_matrix(x)
    mag = x[:, 0] if x.shape[1] == 1 else np.linalg.norm(x, axis=1)
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


def metric_values(target: str, x: np.ndarray, y: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    x = as_matrix(x)
    y = as_matrix(y)
    score = np.asarray(mask, dtype=bool) & finite_rows(x, y)
    if score.sum() < 3:
        return {"scalar_r": math.nan, "component_r": math.nan, "absT2_r": math.nan}
    if target == "T0":
        return {
            "scalar_r": corrcoef(x[score], y[score]),
            "component_r": math.nan,
            "absT2_r": math.nan,
        }
    return {
        "scalar_r": math.nan,
        "component_r": corrcoef(x[score].reshape(-1), y[score].reshape(-1)),
        "absT2_r": corrcoef(np.linalg.norm(x[score], axis=1), np.linalg.norm(y[score], axis=1)),
    }


def verdict(gamma: float, se: float) -> str:
    if not np.isfinite(gamma):
        return "insufficient"
    if np.isfinite(se) and abs(gamma - 1.0) <= 2.0 * se:
        return "zero-circularity recovered law"
    return "form-recovered-scale-fitted"


def thin_label(neff: float, threshold: float) -> str:
    if not np.isfinite(neff):
        return "thin"
    return "thin" if neff < threshold else "ok"


def axis_neff(row: dict[str, object]) -> float:
    if row["axis"] == "between":
        neff = float(row["atom_mean_signal_neff"])
        if np.isfinite(neff):
            return neff
        return float(row["rows_scored"])
    return float(row["atom_signal_neff"])


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


def summarize_axis(
    rows: list[dict[str, object]],
    stratum: str,
    target: str,
    x: np.ndarray,
    y: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    valid: np.ndarray,
    thin_threshold: float,
) -> None:
    x = as_matrix(x)
    y = as_matrix(y)
    valid = valid & finite_rows(x, y)
    atom_neff, atom_top90, active_atoms = signal_effective_atoms(x, atoms, valid)
    ar1_neff, ar1_lag1 = ar1_effective_frames(x, atoms, frames, valid)
    atoms_total = int(np.unique(atoms[valid]).size)

    labels, xm, ym = atom_means(x, y, atoms, valid)
    atom_mean_neff = atom_mean_signal_neff(xm)
    if target == "T0":
        gamma_b = math.nan
        intercept = np.full(ym.shape[1], math.nan)
        gamma_se = math.nan
        intercept_se = math.nan
        loao = math.nan
    else:
        gamma_b, intercept = gamma_with_intercept(xm, ym)
        gamma_se, intercept_se = jk_se_between(xm, ym, labels)
        loao = loao_between_r2(xm, ym, labels)
    metrics = metric_values(target, xm, ym, finite_rows(xm, ym))
    rows.append(
        {
            "stratum": stratum,
            "target": target,
            "axis": "between",
            "rows_scored": int(labels.size),
            "atoms_total": atoms_total,
            "atoms_active": int(labels.size),
            "atom_signal_neff": atom_neff,
            "atom_signal_top90_count": atom_top90,
            "atom_mean_signal_neff": atom_mean_neff,
            "ar1_frame_neff": ar1_neff,
            "kernel_median_lag1": ar1_lag1,
            "thin_flag": thin_label(atom_mean_neff if np.isfinite(atom_mean_neff) else float(labels.size), thin_threshold),
            "gamma_lit": gamma_b,
            "gamma_lit_se": gamma_se,
            "intercept": float(intercept[0]) if intercept.size == 1 else math.nan,
            "intercept_se": intercept_se if intercept.size == 1 else math.nan,
            "intercept_norm": float(np.linalg.norm(intercept)) if intercept.size > 1 else math.nan,
            "scalar_r": metrics["scalar_r"],
            "component_r": metrics["component_r"],
            "absT2_r": metrics["absT2_r"],
            "loao_scaled_R2": loao,
            "verdict_bucket": verdict(gamma_b, gamma_se),
        }
    )

    xw, yw = center_by_atom(x, y, atoms, valid)
    score = valid & finite_rows(xw, yw)
    if target == "T0":
        gamma = math.nan
        se = math.nan
        loao = math.nan
    else:
        gamma = gamma_no_intercept(xw, yw, score)
        se = jk_se_within(xw, yw, atoms, score)
        loao = loao_within_r2(xw, yw, atoms, score)
    metrics = metric_values(target, xw, yw, score)
    rows.append(
        {
            "stratum": stratum,
            "target": target,
            "axis": "within",
            "rows_scored": int(score.sum()),
            "atoms_total": atoms_total,
            "atoms_active": active_atoms,
            "atom_signal_neff": atom_neff,
            "atom_signal_top90_count": atom_top90,
            "atom_mean_signal_neff": math.nan,
            "ar1_frame_neff": ar1_neff,
            "kernel_median_lag1": ar1_lag1,
            "thin_flag": thin_label(atom_neff, thin_threshold),
            "gamma_lit": gamma,
            "gamma_lit_se": se,
            "intercept": math.nan,
            "intercept_se": math.nan,
            "intercept_norm": math.nan,
            "scalar_r": metrics["scalar_r"],
            "component_r": metrics["component_r"],
            "absT2_r": metrics["absT2_r"],
            "loao_scaled_R2": loao,
            "verdict_bucket": verdict(gamma, se),
        }
    )


def fmt(x: object, digits: int = 4) -> str:
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


def gamma_pair(row: dict[str, object]) -> str:
    return f"{fmt(row['gamma_lit'])} +/- {fmt(row['gamma_lit_se'])}"


def write_csv(rows: list[dict[str, object]], path: Path) -> None:
    fields = [
        "stratum",
        "target",
        "axis",
        "rows_scored",
        "atoms_total",
        "atoms_active",
        "atom_signal_neff",
        "atom_signal_top90_count",
        "atom_mean_signal_neff",
        "ar1_frame_neff",
        "kernel_median_lag1",
        "thin_flag",
        "gamma_lit",
        "gamma_lit_se",
        "intercept",
        "intercept_se",
        "intercept_norm",
        "scalar_r",
        "component_r",
        "absT2_r",
        "loao_scaled_R2",
        "verdict_bucket",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_report(
    rows: list[dict[str, object]],
    audit: dict[str, object],
    path: Path,
    csv_path: Path,
    audit_path: Path,
    out_dir: Path,
    source_mode: str,
    t2_cols: list[str],
    t0_col: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("# McConnell literature-scaled de-circularisation\n\n")
        f.write(
            f"Read-only analysis of emitted `{t0_col}` and "
            f"`{t2_cols[0]}`..`{t2_cols[-1]}` columns from "
            "`broad_backbone_aggregated.csv`. "
            "Python only correlates emitted columns against emitted DFT targets; "
            "it does not open `trajectory.h5`, apply Delta-chi, rebuild tensors, "
            "or call the frozen change-of-basis helper.\n\n"
        )
        f.write(f"McConnell source mode: `{source_mode}`\n\n")
        f.write(
            "Delta-chi values are provisional single-family Williamson-Asakura "
            "values from `src/rediscover/MCCONNELL_DCHI_LITERATURE.md`: peptide "
            "C=O +2.41, peptide C-N -5.42, sidechain C=O +2.41, aromatic 0 "
            "in 10^-6 cm^3/mol. Aromatic is zero because RING carries the pi "
            "current.\n\n"
        )
        f.write(f"Input out-dir: `{out_dir}`\n\n")
        f.write(f"CSV artifact: `{csv_path}`\n\n")
        f.write(f"Audit artifact: `{audit_path}`\n\n")
        f.write("## T2 static lead\n\n")
        f.write("Between-axis atom means are the default static read for this near-static McConnell mechanism.\n\n")
        f.write("| stratum | axis | rows | atoms | N_eff(atom means) | gamma_lit +/- SE | component_r | absT2_r | LOAO R2 | bucket |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---|\n")
        for row in rows:
            if row["target"] != "T2" or row["axis"] != "between":
                continue
            f.write(
                "| {stratum} | {axis} | {rows} | {atoms} | {neff} | {gamma} | {comp} | {mag} | {loao} | {bucket} |\n".format(
                    stratum=row["stratum"],
                    axis=row["axis"],
                    rows=row["rows_scored"],
                    atoms=row["atoms_total"],
                    neff=fmt(axis_neff(row)),
                    gamma=gamma_pair(row),
                    comp=fmt(row["component_r"]),
                    mag=fmt(row["absT2_r"]),
                    loao=fmt(row["loao_scaled_R2"]),
                    bucket=row["verdict_bucket"],
                )
            )
        f.write("\n## T2 dynamic diagnostics\n\n")
        f.write("Within-axis rows are per-atom de-meaned frame-modulation diagnostics, not the static headline.\n\n")
        f.write("| stratum | axis | rows | atoms | N_eff(atom modulation) | N_eff(AR1 frames) | gamma_lit +/- SE | component_r | absT2_r | LOAO R2 | bucket |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n")
        for row in rows:
            if row["target"] != "T2" or row["axis"] != "within":
                continue
            f.write(
                "| {stratum} | {axis} | {rows} | {atoms} | {neff} | {ar1} | {gamma} | {comp} | {mag} | {loao} | {bucket} |\n".format(
                    stratum=row["stratum"],
                    axis=row["axis"],
                    rows=row["rows_scored"],
                    atoms=row["atoms_total"],
                    neff=fmt(axis_neff(row)),
                    ar1=fmt(row["ar1_frame_neff"]),
                    gamma=gamma_pair(row),
                    comp=fmt(row["component_r"]),
                    mag=fmt(row["absT2_r"]),
                    loao=fmt(row["loao_scaled_R2"]),
                    bucket=row["verdict_bucket"],
                )
            )
        f.write("\n## T0 trace audit\n\n")
        f.write("| stratum | axis | rows | atoms | gamma_lit +/- SE | scalar_r | max_abs_mc_lit_T0_ppm | bucket |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---|\n")
        max_t0_by_stratum = audit.get("max_abs_mc_lit_T0_by_stratum", {})
        for row in rows:
            if row["target"] != "T0":
                continue
            f.write(
                "| {stratum} | {axis} | {rows} | {atoms} | {gamma} | {r} | {max_t0} | {bucket} |\n".format(
                    stratum=row["stratum"],
                    axis=row["axis"],
                    rows=row["rows_scored"],
                    atoms=row["atoms_total"],
                    gamma=gamma_pair(row),
                    r=fmt(row["scalar_r"]),
                    max_t0=fmt(max_t0_by_stratum.get(row["stratum"], math.nan), 6),
                    bucket=row["verdict_bucket"],
                )
            )
        f.write("\n## Notes\n\n")
        f.write("- `mc_lit_T0` is a trace audit channel; the emitted PCS tensor is traceless, so gamma is usually undefined.\n")
        f.write("- `gamma_lit` is a scale diagnostic only. Correlations are unfitted component and magnitude comparisons.\n")
        f.write("- Confidence is within-protein only: between rows use atom-count / atom-mean signal N_eff; within diagnostics additionally report AR(1)-deflated frame N_eff. Thin strata are flagged, not forced.\n")
        f.write("- Verdict buckets are mechanical diagnostics; interpretation is reserved for the lead.\n")


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    artifact_dir = Path(args.artifact_dir)
    agg_path = out_dir / "broad_backbone_aggregated.csv"
    target_path = out_dir / "broad_backbone_aggregated_target_local_T2.npy"
    if not agg_path.exists():
        raise RuntimeError(f"missing broad aggregate CSV: {agg_path}")
    if not target_path.exists():
        raise RuntimeError(f"missing local T2 target sidecar: {target_path}")

    if args.mc_source_mode == "valid":
        t2_cols = VALID_T2_COLS
        t0_col = "mc_lit_T0_valid"
    else:
        t2_cols = ALL_T2_COLS
        t0_col = "mc_lit_T0"

    df = pd.read_csv(agg_path)
    target_t2 = np.load(target_path)
    if target_t2.shape != (len(df), 5):
        raise RuntimeError(f"target local T2 shape {target_t2.shape} does not match {(len(df), 5)}")

    required = [
        "atom_index",
        "atom_name",
        "h5_row",
        "dft_present",
        "dft_local_frame_valid",
        "dft_sigma_iso",
        t0_col,
        *t2_cols,
    ]
    require_columns(df, required, "broad_backbone_aggregated.csv")
    df = df.copy()
    df["stratum"] = [row_stratum(v) for v in df["atom_name"]]

    rows: list[dict[str, object]] = []
    atoms = df["atom_index"].to_numpy(int)
    frames = df["h5_row"].to_numpy(int)
    present = (
        df["dft_present"].to_numpy(int).astype(bool)
        & df["dft_local_frame_valid"].to_numpy(int).astype(bool)
    )
    x0_all = df[t0_col].to_numpy(float).reshape(-1, 1)
    y0_all = df["dft_sigma_iso"].to_numpy(float).reshape(-1, 1)
    x2_all = df[t2_cols].to_numpy(float)
    y2_all = target_t2.astype(float)

    max_t0_by_stratum: dict[str, float] = {}
    for stratum in STRATA:
        sm = (df["stratum"] == stratum).to_numpy(bool)
        if sm.sum() == 0:
            continue
        valid0 = present & sm & finite_rows(x0_all, y0_all)
        valid2 = present & sm & finite_rows(x2_all, y2_all)
        if valid0.any():
            max_t0_by_stratum[stratum] = float(np.nanmax(np.abs(x0_all[valid0])))
        summarize_axis(
            rows,
            stratum,
            "T2",
            x2_all,
            y2_all,
            atoms,
            frames,
            valid2,
            args.thin_atom_neff,
        )
        summarize_axis(
            rows,
            stratum,
            "T0",
            x0_all,
            y0_all,
            atoms,
            frames,
            valid0,
            args.thin_atom_neff,
        )

    csv_path = artifact_dir / "mcconnell_literature_decirc.csv"
    audit_path = artifact_dir / "mcconnell_literature_decirc_audit.json"
    write_csv(rows, csv_path)
    audit = {
        "out_dir": str(out_dir),
        "rows": len(df),
        "strata_present": sorted(s for s in STRATA if (df["stratum"] == s).any()),
        "max_abs_mc_lit_T0_by_stratum": max_t0_by_stratum,
        "inputs": {
            "aggregate_csv": str(agg_path),
            "target_local_T2": str(target_path),
            "source_mode": args.mc_source_mode,
            "prediction_columns": [t0_col, *t2_cols],
        },
        "no_physics_recompute": True,
    }
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    write_report(rows, audit, Path(args.report_md), csv_path, audit_path, out_dir,
                 args.mc_source_mode, t2_cols, t0_col)
    print(f"wrote {csv_path}")
    print(f"wrote {audit_path}")
    print(f"wrote {args.report_md}")


if __name__ == "__main__":
    main()
