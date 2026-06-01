#!/usr/bin/env python3
"""Frame-to-frame differencing system-ID evidence for broad-backbone emit.

This is a read-only Python consumer of the emitted broad-backbone CSV/NPY
substrate. It subtracts emitted per-frame columns across DFT-sampled frame
windows, aligns emitted source rows with stable physical keys, and reports
evidence for the geometric perturbation probe:

    delta emitted geometry/features -> delta emitted shielding target

It only subtracts emitted values; no physical quantity is rebuilt here.

Outputs under --analysis-dir:
  atom_pairs.csv.gz
  aggregate_delta_correlations.csv
  aggregate_gate_summary.csv
  aggregate_decorrelation.csv
  source_delta_correlations.csv
  source_feature_gate.csv
  source_decorrelation.csv
  source_churn.csv
  differencing_system_id_report.md

The full source-pair table can be very large; write it only with
--write-source-pairs.
"""

from __future__ import annotations

import argparse
import gzip
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


FV_HN = {1, 2}
FV_N = {4, 5}
FV_CA = {6}
FV_C = {7}
FV_O = {8}
FV_HA = {9}

STRATA_ORDER = ["ALL", "N", "CA", "C", "O", "HN", "HA", "HA2", "HA3"]
MECHANISMS = ["ring", "bond", "charge"]

TARGET_LOCAL_T2 = [f"target_local_T2_{i}" for i in range(5)]

VECTOR_GROUPS = {
    "target_local_T2": TARGET_LOCAL_T2,
    "literature_kernel_T2": [f"literature_kernel_T2_{i}" for i in range(5)],
    "ring_literature_kernel_T2": [f"ring_literature_kernel_T2_{i}" for i in range(5)],
    "bond_literature_kernel_T2": [f"bond_literature_kernel_T2_{i}" for i in range(5)],
    "charge_literature_kernel_T2": [
        f"charge_literature_kernel_T2_{i}" for i in range(5)
    ],
    "field_local": ["field_local_x", "field_local_y", "field_local_z"],
    "mu_local": ["mu_local_x", "mu_local_y", "mu_local_z"],
}

AGG_FEATURES_BY_MECH = {
    "ring": [
        "ring_n",
        "ring_sum_dipolar",
        *VECTOR_GROUPS["ring_literature_kernel_T2"],
    ],
    "bond": [
        "bond_n",
        "bond_sum_dipolar",
        *VECTOR_GROUPS["bond_literature_kernel_T2"],
    ],
    "charge": [
        "charge_n",
        "field_local_x",
        "field_local_y",
        "field_local_z",
        "field_z",
        "field_mag",
        "mu_local_x",
        "mu_local_y",
        "mu_local_z",
        *VECTOR_GROUPS["charge_literature_kernel_T2"],
    ],
    "combined": [*VECTOR_GROUPS["literature_kernel_T2"]],
}

VECTOR_MECH = {
    "literature_kernel_T2": "combined",
    "ring_literature_kernel_T2": "ring",
    "bond_literature_kernel_T2": "bond",
    "charge_literature_kernel_T2": "charge",
    "field_local": "charge",
    "mu_local": "charge",
}

CUTOFF_COLUMN = {
    "ring": "ring_cutoff_A",
    "bond": "bond_cutoff_A",
    "charge": "charge_cutoff_A",
}

SOURCE_KEY_COLUMNS = {
    "ring": ["ring_index"],
    "bond": ["bond_index", "bond_atom_a", "bond_atom_b"],
    "charge": ["source_atom_index", "charge_source_code"],
}

SOURCE_GEOMETRY_COLUMNS = {
    "ring": [
        "disp_local_x",
        "disp_local_y",
        "disp_local_z",
        "r",
        "cos_theta",
        "dipolar",
        "source_normal_local_x",
        "source_normal_local_y",
        "source_normal_local_z",
    ],
    "bond": [
        "disp_local_x",
        "disp_local_y",
        "disp_local_z",
        "r",
        "cos_theta",
        "dipolar",
        "bond_axis_local_x",
        "bond_axis_local_y",
        "bond_axis_local_z",
    ],
    "charge": [
        "disp_local_x",
        "disp_local_y",
        "disp_local_z",
        "r",
        "cos_theta",
        "dipolar",
        "source_q_e",
    ],
}


class CorrAccumulator:
    def __init__(self) -> None:
        self.n = 0
        self.sx = 0.0
        self.sy = 0.0
        self.sxx = 0.0
        self.syy = 0.0
        self.sxy = 0.0

    def add(self, x, y) -> None:
        xa = np.asarray(x, dtype=np.float64)
        ya = np.asarray(y, dtype=np.float64)
        mask = np.isfinite(xa) & np.isfinite(ya)
        if not np.any(mask):
            return
        xa = xa[mask]
        ya = ya[mask]
        self.n += int(xa.size)
        self.sx += float(xa.sum())
        self.sy += float(ya.sum())
        self.sxx += float(np.dot(xa, xa))
        self.syy += float(np.dot(ya, ya))
        self.sxy += float(np.dot(xa, ya))

    def corr(self) -> float:
        if self.n < 3:
            return math.nan
        n = float(self.n)
        cov = self.sxy - self.sx * self.sy / n
        vx = self.sxx - self.sx * self.sx / n
        vy = self.syy - self.sy * self.sy / n
        if vx <= 0.0 or vy <= 0.0:
            return math.nan
        return float(cov / math.sqrt(vx * vy))


class LagAccumulator:
    def __init__(self) -> None:
        self.n = 0
        self.sx = 0.0
        self.sy = 0.0
        self.sxx = 0.0
        self.syy = 0.0
        self.sxy = 0.0
        self.sd = 0.0
        self.sdd = 0.0

    def add(self, x_t, x_tp1) -> None:
        x = np.asarray(x_t, dtype=np.float64)
        y = np.asarray(x_tp1, dtype=np.float64)
        mask = np.isfinite(x) & np.isfinite(y)
        if not np.any(mask):
            return
        x = x[mask]
        y = y[mask]
        d = y - x
        self.n += int(x.size)
        self.sx += float(x.sum())
        self.sy += float(y.sum())
        self.sxx += float(np.dot(x, x))
        self.syy += float(np.dot(y, y))
        self.sxy += float(np.dot(x, y))
        self.sd += float(d.sum())
        self.sdd += float(np.dot(d, d))

    def rho(self) -> float:
        if self.n < 3:
            return math.nan
        n = float(self.n)
        cov = self.sxy - self.sx * self.sy / n
        vx = self.sxx - self.sx * self.sx / n
        vy = self.syy - self.sy * self.sy / n
        if vx <= 0.0 or vy <= 0.0:
            return math.nan
        return float(cov / math.sqrt(vx * vy))

    def noise_cost(self) -> float:
        if self.n < 3:
            return math.nan
        n = float(self.n)
        vx = self.sxx - self.sx * self.sx / n
        vd = self.sdd - self.sd * self.sd / n
        if vx <= 0.0:
            return math.nan
        return float(vd / (2.0 * vx))


class MatrixAccumulator:
    def __init__(self, columns: list[str]) -> None:
        self.columns = columns
        self.n = 0
        k = len(columns)
        self.s = np.zeros(k, dtype=np.float64)
        self.ss = np.zeros((k, k), dtype=np.float64)

    def add(self, x) -> None:
        if len(self.columns) < 2:
            return
        a = np.asarray(x, dtype=np.float64)
        if a.ndim != 2 or a.shape[1] != len(self.columns):
            raise RuntimeError("matrix accumulator shape mismatch")
        mask = np.isfinite(a).all(axis=1)
        if not np.any(mask):
            return
        a = a[mask]
        self.n += int(a.shape[0])
        self.s += a.sum(axis=0)
        self.ss += a.T @ a

    def stats(self) -> dict[str, float | int]:
        return corr_matrix_stats_from_sums(self.n, self.s, self.ss, self.columns)


def parse_args():
    ap = argparse.ArgumentParser(
        description="Broad-backbone frame differencing system-ID evidence."
    )
    ap.add_argument(
        "out_dir",
        nargs="?",
        default="/tmp/rdc-broad-backbone-axes",
        help="directory containing broad_backbone_aggregated.csv and sources CSV",
    )
    ap.add_argument(
        "--analysis-dir",
        default=None,
        help="where to write differencing evidence; default: <out_dir>/differencing_system_id",
    )
    ap.add_argument("--windows", nargs="*", type=int, default=[1, 3, 5])
    ap.add_argument("--chunksize", type=int, default=1_000_000)
    ap.add_argument("--smooth-rho-min", type=float, default=0.30)
    ap.add_argument("--decorrelation-drop-min", type=float, default=0.05)
    ap.add_argument("--min-n", type=int, default=30)
    ap.add_argument("--write-source-pairs", action="store_true")
    ap.add_argument(
        "--skip-source",
        action="store_true",
        help="only build atom pairs and aggregate summaries",
    )
    return ap.parse_args()


def stratum_of(frame_variant, atom_name) -> str | None:
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
        nm = str(atom_name).strip().upper()
        if nm == "HA2":
            return "HA2"
        if nm == "HA3":
            return "HA3"
        return "HA"
    return None


def finite_pearson(x, y, min_n: int = 3) -> tuple[float, int]:
    xa = np.asarray(x, dtype=np.float64)
    ya = np.asarray(y, dtype=np.float64)
    mask = np.isfinite(xa) & np.isfinite(ya)
    n = int(mask.sum())
    if n < min_n:
        return math.nan, n
    xa = xa[mask]
    ya = ya[mask]
    xa = xa - xa.mean()
    ya = ya - ya.mean()
    denom = math.sqrt(float(np.dot(xa, xa) * np.dot(ya, ya)))
    if denom <= 0.0:
        return math.nan, n
    return float(np.dot(xa, ya) / denom), n


def abs_or_nan(x: float) -> float:
    return abs(x) if np.isfinite(x) else math.nan


def fmt_float(x: float, digits: int = 4) -> str:
    if x is None or not np.isfinite(x):
        return "nan"
    return f"{x:.{digits}g}"


def require_columns(df: pd.DataFrame, names: list[str], label: str) -> None:
    missing = [name for name in names if name not in df.columns]
    if missing:
        raise RuntimeError(f"{label} missing emitted columns: {', '.join(missing)}")


def add_vector_norms(df: pd.DataFrame) -> list[str]:
    created = []
    for prefix, cols in VECTOR_GROUPS.items():
        if all(c in df.columns for c in cols):
            name = f"levelnorm_{prefix}"
            df[name] = np.linalg.norm(df[cols].to_numpy(np.float64), axis=1)
            created.append(name)
    return created


def load_aggregate(out_dir: Path) -> pd.DataFrame:
    agg_path = out_dir / "broad_backbone_aggregated.csv"
    target_path = out_dir / "broad_backbone_aggregated_target_local_T2.npy"
    if not agg_path.exists():
        raise RuntimeError(f"missing aggregate CSV: {agg_path}")
    if not target_path.exists():
        raise RuntimeError(f"missing local T2 target sidecar: {target_path}")

    agg = pd.read_csv(agg_path)
    require_columns(
        agg,
        [
            "row_id",
            "atom_index",
            "residue_index",
            "residue_number",
            "atom_name",
            "frame_variant",
            "frame_valid",
            "h5_row",
            "original_index",
            "time_ps",
            "dft_present",
            "dft_sigma_iso",
            "dft_local_frame_valid",
            "charge_source",
        ],
        "aggregate CSV",
    )
    local_t2 = np.load(target_path)
    if local_t2.shape != (len(agg), 5):
        raise RuntimeError(
            f"local T2 sidecar shape {local_t2.shape} does not match {(len(agg), 5)}"
        )
    for i, col in enumerate(TARGET_LOCAL_T2):
        agg[col] = local_t2[:, i]
    agg["stratum"] = [
        stratum_of(fv, atom)
        for fv, atom in zip(agg["frame_variant"].to_numpy(), agg["atom_name"].to_numpy())
    ]
    codes, uniques = pd.factorize(agg["charge_source"].astype(str), sort=True)
    agg["charge_source_code"] = codes.astype(np.int64)
    agg.attrs["charge_source_levels"] = [str(x) for x in uniques]
    add_vector_norms(agg)
    return agg


def aggregate_feature_columns(agg: pd.DataFrame) -> list[str]:
    cols = []
    for values in AGG_FEATURES_BY_MECH.values():
        cols.extend(values)
    cols.extend(CUTOFF_COLUMN.values())
    cols.extend(["charge_source_code"])
    cols.extend([f"levelnorm_{name}" for name in VECTOR_GROUPS if name != "target_local_T2"])
    out = []
    seen = set()
    for col in cols:
        if col in agg.columns and col not in seen:
            out.append(col)
            seen.add(col)
    return out


def build_atom_pairs(agg: pd.DataFrame, windows: list[int]) -> pd.DataFrame:
    feature_cols = aggregate_feature_columns(agg)
    pair_value_cols = [
        "dft_sigma_iso",
        *TARGET_LOCAL_T2,
        "levelnorm_target_local_T2",
        *feature_cols,
    ]
    pair_value_cols = [c for c in pair_value_cols if c in agg.columns]

    keep = (
        (agg["dft_present"] == 1)
        & (agg["frame_valid"] == 1)
        & pd.notna(agg["stratum"])
    )
    base = agg.loc[keep].sort_values(["atom_index", "h5_row"]).copy()
    pieces = []
    pair_next = 0
    for window in sorted(set(windows)):
        if window <= 0:
            raise RuntimeError("windows must be positive DFT-sample offsets")
        for _, g in base.groupby("atom_index", sort=False):
            if len(g) <= window:
                continue
            left = g.iloc[:-window].reset_index(drop=True)
            right = g.iloc[window:].reset_index(drop=True)
            n = len(left)
            data = {
                "pair_id": np.arange(pair_next, pair_next + n, dtype=np.int64),
                "window": np.full(n, window, dtype=np.int64),
                "atom_index": left["atom_index"].to_numpy(np.int64),
                "residue_index": left["residue_index"].to_numpy(np.int64),
                "residue_number": left["residue_number"].to_numpy(np.int64),
                "atom_name": left["atom_name"].astype(str).to_numpy(),
                "stratum": left["stratum"].astype(str).to_numpy(),
                "frame_variant": left["frame_variant"].to_numpy(np.int64),
                "row_id_t": left["row_id"].to_numpy(np.int64),
                "row_id_tp1": right["row_id"].to_numpy(np.int64),
                "h5_row_t": left["h5_row"].to_numpy(np.int64),
                "h5_row_tp1": right["h5_row"].to_numpy(np.int64),
                "original_index_t": left["original_index"].to_numpy(np.int64),
                "original_index_tp1": right["original_index"].to_numpy(np.int64),
                "time_ps_t": left["time_ps"].to_numpy(np.float64),
                "time_ps_tp1": right["time_ps"].to_numpy(np.float64),
                "dt_ps": right["time_ps"].to_numpy(np.float64)
                - left["time_ps"].to_numpy(np.float64),
            }
            for col in pair_value_cols:
                lt = left[col].to_numpy(np.float64)
                rt = right[col].to_numpy(np.float64)
                data[f"{col}_t"] = lt
                data[f"{col}_tp1"] = rt
                data[f"d_{col}"] = rt - lt
            for prefix, cols in VECTOR_GROUPS.items():
                if all(f"d_{c}" in data for c in cols):
                    dcols = [f"d_{c}" for c in cols]
                    data[f"dnorm_{prefix}"] = np.linalg.norm(
                        np.column_stack([data[c] for c in dcols]), axis=1
                    )
            out = pd.DataFrame(data)
            pair_next += n
            pieces.append(out)
    if not pieces:
        raise RuntimeError("no DFT frame pairs could be built")
    return pd.concat(pieces, ignore_index=True)


def lag1_autocorr(x) -> float:
    arr = np.asarray(x, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size < 3:
        return math.nan
    y = arr - arr.mean()
    denom = float(np.dot(y, y))
    if denom <= 0.0:
        return math.nan
    return float(np.dot(y[:-1], y[1:]) / denom)


def level_rhos_by_atom_groups(groups: list[pd.DataFrame], col: str) -> list[tuple[float, int]]:
    out = []
    for g in groups:
        x = g[col].to_numpy(np.float64)
        rho = lag1_autocorr(x)
        if np.isfinite(rho):
            out.append((rho, int(np.isfinite(x).sum())))
    return out


def effective_n_from_rhos(rhos_and_n: list[tuple[float, int]]) -> float:
    total = 0.0
    for rho, n in rhos_and_n:
        rho = max(min(float(rho), 0.999), -0.999)
        neff = n * (1.0 - rho) / (1.0 + rho)
        total += min(max(neff, 1.0), float(n))
    return float(total)


def direct_noise_cost(level, delta) -> float:
    x = np.asarray(level, dtype=np.float64)
    d = np.asarray(delta, dtype=np.float64)
    mask = np.isfinite(x) & np.isfinite(d)
    if int(mask.sum()) < 3:
        return math.nan
    x = x[mask]
    d = d[mask]
    vx = float(np.var(x))
    if vx <= 0.0:
        return math.nan
    return float(np.var(d) / (2.0 * vx))


def median_window_noise_cost_by_atom_groups(groups: list[pd.DataFrame], col: str) -> float:
    level_col = f"{col}_t"
    delta_col = f"d_{col}"
    if not groups or level_col not in groups[0].columns or delta_col not in groups[0].columns:
        return math.nan
    values = []
    for g in groups:
        cost = direct_noise_cost(g[level_col], g[delta_col])
        if np.isfinite(cost):
            values.append(cost)
    return float(np.median(values)) if values else math.nan


def stratum_mask(df: pd.DataFrame, stratum: str) -> pd.Series:
    if stratum == "ALL":
        return pd.Series(True, index=df.index)
    return df["stratum"] == stratum


def available_strata(df: pd.DataFrame) -> list[str]:
    present = set(str(x) for x in df["stratum"].dropna().unique())
    return [s for s in STRATA_ORDER if s == "ALL" or s in present]


def aggregate_gate_summary(
    agg: pd.DataFrame,
    pairs: pd.DataFrame,
    windows: list[int],
    smooth_rho_min: float,
) -> pd.DataFrame:
    target_cols = ["dft_sigma_iso", "levelnorm_target_local_T2"]
    feature_cols = aggregate_feature_columns(agg)
    cols = [c for c in [*target_cols, *feature_cols] if c in agg.columns]
    rows = []
    for stratum in available_strata(agg):
        df_s = agg.loc[
            (agg["dft_present"] == 1)
            & (agg["frame_valid"] == 1)
            & pd.notna(agg["stratum"])
            & stratum_mask(agg, stratum)
        ]
        groups = [
            g
            for _, g in df_s.sort_values(["atom_index", "h5_row"]).groupby(
                "atom_index", sort=False
            )
        ]
        pair_s_all = pairs.loc[stratum_mask(pairs, stratum)]
        pair_groups_by_window = {}
        for window in sorted(set(windows)):
            pair_w = pair_s_all[pair_s_all["window"] == window]
            pair_groups_by_window[window] = [
                g for _, g in pair_w.groupby("atom_index", sort=False)
            ]
        for col in cols:
            rhos_n = level_rhos_by_atom_groups(groups, col)
            rhos = np.array([x[0] for x in rhos_n], dtype=np.float64)
            neff = effective_n_from_rhos(rhos_n)
            med_rho = float(np.nanmedian(rhos)) if rhos.size else math.nan
            for window in sorted(set(windows)):
                pair_s = pair_s_all[pair_s_all["window"] == window]
                d_col = f"d_{col}"
                if d_col not in pair_s.columns:
                    continue
                rows.append(
                    {
                        "scope": "target"
                        if col in target_cols
                        else "aggregate_feature",
                        "window": window,
                        "stratum": stratum,
                        "feature": col,
                        "n_level_rows": int(len(df_s)),
                        "n_pairs": int(len(pair_s)),
                        "n_series": int(len(rhos_n)),
                        "effective_n": neff,
                        "median_lag1_rho": med_rho,
                        "median_noise_approx_1_minus_rho": 1.0 - med_rho
                        if np.isfinite(med_rho)
                        else math.nan,
                        "median_atom_window_noise_cost": median_window_noise_cost_by_atom_groups(
                            pair_groups_by_window[window], col
                        ),
                        "smooth_pass": bool(
                            np.isfinite(med_rho) and med_rho >= smooth_rho_min
                        ),
                    }
                )
    return pd.DataFrame(rows)


def atom_pair_output_columns(pairs: pd.DataFrame) -> list[str]:
    base = [
        "pair_id",
        "window",
        "atom_index",
        "residue_index",
        "residue_number",
        "atom_name",
        "stratum",
        "frame_variant",
        "row_id_t",
        "row_id_tp1",
        "h5_row_t",
        "h5_row_tp1",
        "original_index_t",
        "original_index_tp1",
        "time_ps_t",
        "time_ps_tp1",
        "dt_ps",
        "dft_sigma_iso_t",
        "dft_sigma_iso_tp1",
        "d_dft_sigma_iso",
        "dnorm_target_local_T2",
        "d_levelnorm_target_local_T2",
    ]
    target = []
    for col in TARGET_LOCAL_T2:
        target.extend([f"{col}_t", f"{col}_tp1", f"d_{col}"])
    feature_delta = []
    for col in sorted({c for cols in AGG_FEATURES_BY_MECH.values() for c in cols}):
        feature_delta.append(f"d_{col}")
    for col in CUTOFF_COLUMN.values():
        feature_delta.extend([f"{col}_t", f"{col}_tp1", f"d_{col}"])
    for prefix in VECTOR_GROUPS:
        if prefix != "target_local_T2":
            feature_delta.append(f"dnorm_{prefix}")
            feature_delta.append(f"d_levelnorm_{prefix}")
    wanted = [*base, *target, *feature_delta]
    return [c for c in dict.fromkeys(wanted) if c in pairs.columns]


def corr_matrix_stats_from_sums(
    n: int, sums: np.ndarray, cross: np.ndarray, columns: list[str]
) -> dict[str, float | int]:
    if n < 3 or len(columns) < 2:
        return {
            "n_rows": int(n),
            "n_features_active": 0,
            "max_abs_corr": math.nan,
            "median_abs_corr": math.nan,
        }
    nf = float(n)
    cov = cross - np.outer(sums, sums) / nf
    var = np.diag(cov)
    active = var > 0.0
    if int(active.sum()) < 2:
        return {
            "n_rows": int(n),
            "n_features_active": int(active.sum()),
            "max_abs_corr": math.nan,
            "median_abs_corr": math.nan,
        }
    cov = cov[np.ix_(active, active)]
    std = np.sqrt(np.diag(cov))
    corr = cov / np.outer(std, std)
    off = np.abs(corr[np.triu_indices_from(corr, k=1)])
    off = off[np.isfinite(off)]
    if off.size == 0:
        max_abs = math.nan
        med_abs = math.nan
    else:
        max_abs = float(np.max(off))
        med_abs = float(np.median(off))
    return {
        "n_rows": int(n),
        "n_features_active": int(active.sum()),
        "max_abs_corr": max_abs,
        "median_abs_corr": med_abs,
    }


def corr_matrix_stats(x: np.ndarray, columns: list[str]) -> dict[str, float | int]:
    if x.ndim != 2:
        raise RuntimeError("expected matrix")
    if x.shape[1] != len(columns):
        raise RuntimeError("matrix column count mismatch")
    mask = np.isfinite(x).all(axis=1)
    if not np.any(mask):
        return corr_matrix_stats_from_sums(0, np.zeros(len(columns)), np.zeros((len(columns), len(columns))), columns)
    a = x[mask].astype(np.float64, copy=False)
    return corr_matrix_stats_from_sums(a.shape[0], a.sum(axis=0), a.T @ a, columns)


def aggregate_decorrelation(
    pairs: pd.DataFrame, windows: list[int], drop_min: float
) -> pd.DataFrame:
    rows = []
    for stratum in available_strata(pairs):
        pair_s_all = pairs.loc[stratum_mask(pairs, stratum)]
        for window in sorted(set(windows)):
            pair_s = pair_s_all[pair_s_all["window"] == window]
            for mech, cols in AGG_FEATURES_BY_MECH.items():
                cols = [
                    c
                    for c in cols
                    if f"{c}_t" in pair_s.columns and f"d_{c}" in pair_s.columns
                ]
                if len(cols) < 2:
                    continue
                level = pair_s[[f"{c}_t" for c in cols]].to_numpy(np.float64)
                diff = pair_s[[f"d_{c}" for c in cols]].to_numpy(np.float64)
                lev = corr_matrix_stats(level, cols)
                dif = corr_matrix_stats(diff, cols)
                lev_med = float(lev["median_abs_corr"])
                dif_med = float(dif["median_abs_corr"])
                rows.append(
                    {
                        "window": window,
                        "stratum": stratum,
                        "feature_set": mech,
                        "level_n_rows": lev["n_rows"],
                        "diff_n_rows": dif["n_rows"],
                        "level_n_features_active": lev["n_features_active"],
                        "diff_n_features_active": dif["n_features_active"],
                        "level_max_abs_corr": lev["max_abs_corr"],
                        "diff_max_abs_corr": dif["max_abs_corr"],
                        "level_median_abs_corr": lev_med,
                        "diff_median_abs_corr": dif_med,
                        "median_abs_corr_drop": lev_med - dif_med
                        if np.isfinite(lev_med) and np.isfinite(dif_med)
                        else math.nan,
                        "decorrelation_pass": bool(
                            np.isfinite(lev_med)
                            and np.isfinite(dif_med)
                            and (lev_med - dif_med) >= drop_min
                        ),
                    }
                )
    return pd.DataFrame(rows)


def aggregate_delta_correlations(
    pairs: pd.DataFrame, windows: list[int], min_n: int
) -> pd.DataFrame:
    rows = []
    target_specs = [
        ("dft_sigma_iso", "d_dft_sigma_iso"),
        ("target_local_T2_delta_norm", "dnorm_target_local_T2"),
    ]
    scalar_features = []
    for mech, cols in AGG_FEATURES_BY_MECH.items():
        for col in cols:
            if col in pairs.columns and f"d_{col}" in pairs.columns:
                scalar_features.append((mech, col, f"d_{col}", "column_delta"))
    for prefix, cols in VECTOR_GROUPS.items():
        if prefix == "target_local_T2":
            continue
        if f"dnorm_{prefix}" in pairs.columns:
            scalar_features.append(
                (
                    VECTOR_MECH.get(prefix, "combined"),
                    prefix,
                    f"dnorm_{prefix}",
                    "vector_delta_norm",
                )
            )

    for stratum in available_strata(pairs):
        pair_s_all = pairs.loc[stratum_mask(pairs, stratum)]
        for window in sorted(set(windows)):
            pair_s = pair_s_all[pair_s_all["window"] == window]
            atoms = int(pair_s["atom_index"].nunique())
            for target_name, target_col in target_specs:
                if target_col not in pair_s.columns:
                    continue
                y = pair_s[target_col].to_numpy(np.float64)
                for mech, feature, delta_col, metric in scalar_features:
                    if delta_col not in pair_s.columns:
                        continue
                    r, n = finite_pearson(
                        pair_s[delta_col].to_numpy(np.float64), y, min_n=min_n
                    )
                    rows.append(
                        {
                            "window": window,
                            "stratum": stratum,
                            "mechanism": mech,
                            "feature": feature,
                            "feature_metric": metric,
                            "target": target_name,
                            "target_metric": "delta",
                            "n_pairs": n,
                            "n_atoms": atoms,
                            "pearson_r": r,
                            "abs_pearson_r": abs_or_nan(r),
                        }
                    )

            if not all(f"d_{c}" in pair_s.columns for c in TARGET_LOCAL_T2):
                continue
            y_components = pair_s[[f"d_{c}" for c in TARGET_LOCAL_T2]].to_numpy(
                np.float64
            )
            for prefix, cols in VECTOR_GROUPS.items():
                if prefix == "target_local_T2" or len(cols) != 5:
                    continue
                if not all(f"d_{c}" in pair_s.columns for c in cols):
                    continue
                x_components = pair_s[[f"d_{c}" for c in cols]].to_numpy(np.float64)
                r, n = finite_pearson(
                    x_components.reshape(-1),
                    y_components.reshape(-1),
                    min_n=min_n,
                )
                rows.append(
                    {
                        "window": window,
                        "stratum": stratum,
                        "mechanism": VECTOR_MECH.get(prefix, "combined"),
                        "feature": prefix,
                        "feature_metric": "vector_components_flat",
                        "target": "target_local_T2",
                        "target_metric": "component_delta",
                        "n_pairs": n // 5,
                        "n_atoms": atoms,
                        "pearson_r": r,
                        "abs_pearson_r": abs_or_nan(r),
                    }
                )
    return pd.DataFrame(rows)


def read_sources_for_mechanism(
    source_path: Path,
    mechanism: str,
    chunksize: int,
    row_id_keep: set[int],
    charge_source_by_row: pd.Series,
) -> pd.DataFrame:
    usecols = ["row_id", "mechanism", *SOURCE_KEY_COLUMNS[mechanism], *SOURCE_GEOMETRY_COLUMNS[mechanism]]
    usecols = [c for c in dict.fromkeys(usecols) if c != "charge_source_code"]
    frames = []
    for chunk in pd.read_csv(source_path, usecols=usecols, chunksize=chunksize):
        c = chunk[(chunk["mechanism"] == mechanism) & (chunk["row_id"].isin(row_id_keep))]
        if c.empty:
            continue
        c = c.drop(columns=["mechanism"]).copy()
        if mechanism == "charge":
            c["charge_source_code"] = c["row_id"].map(charge_source_by_row).astype("int64")
        frames.append(c)
    if not frames:
        return pd.DataFrame(columns=[c for c in usecols if c != "mechanism"])
    out = pd.concat(frames, ignore_index=True)
    return out


def update_churn_counter(counter, key, present_both, t_only, tp1_only, pair_ids, n_total_pairs):
    values = counter[key]
    values["present_both"] += int(present_both)
    values["t_only"] += int(t_only)
    values["tp1_only"] += int(tp1_only)
    values["union_sources"] += int(present_both + t_only + tp1_only)
    values["pairs_with_any_source"] += int(pd.unique(pair_ids).size) if len(pair_ids) else 0
    values["total_atom_pairs"] += int(n_total_pairs)


def source_pair_output_frame(
    merged: pd.DataFrame,
    mechanism: str,
    geom_cols: list[str],
    target_cols: list[str],
    key_cols: list[str],
) -> pd.DataFrame:
    out = merged[["pair_id", "window", "stratum", "mechanism", *key_cols]].copy()
    out["present_t"] = merged["_merge"].isin(["left_only", "both"])
    out["present_tp1"] = merged["_merge"].isin(["right_only", "both"])
    if "row_id_t" in merged.columns:
        out["source_row_id_t"] = merged["row_id_t"]
    if "row_id_tp1" in merged.columns:
        out["source_row_id_tp1"] = merged["row_id_tp1"]
    for col in target_cols:
        out[col] = merged[col]
    for col in geom_cols:
        lt = f"{col}_t"
        rt = f"{col}_tp1"
        if lt in merged.columns:
            out[lt] = merged[lt]
        if rt in merged.columns:
            out[rt] = merged[rt]
        if lt in merged.columns and rt in merged.columns:
            out[f"d_{col}"] = merged[rt] - merged[lt]
    return out


def append_csv_gz(path: Path, df: pd.DataFrame, header: bool) -> None:
    mode = "wt" if header else "at"
    with gzip.open(path, mode, newline="") as f:
        df.to_csv(f, index=False, header=header)


def process_sources(
    out_dir: Path,
    analysis_dir: Path,
    pairs: pd.DataFrame,
    windows: list[int],
    chunksize: int,
    min_n: int,
    smooth_rho_min: float,
    drop_min: float,
    write_source_pairs: bool,
    charge_source_by_row: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    source_path = out_dir / "broad_backbone_sources.csv"
    if not source_path.exists():
        raise RuntimeError(f"missing source CSV: {source_path}")

    row_id_keep = set(pairs["row_id_t"].astype(int)).union(
        set(pairs["row_id_tp1"].astype(int))
    )
    corr_acc = defaultdict(CorrAccumulator)
    lag_acc = defaultdict(LagAccumulator)
    level_mat = {}
    diff_mat = {}
    churn = defaultdict(
        lambda: {
            "present_both": 0,
            "t_only": 0,
            "tp1_only": 0,
            "union_sources": 0,
            "pairs_with_any_source": 0,
            "total_atom_pairs": 0,
            "duplicates_t": 0,
            "duplicates_tp1": 0,
        }
    )
    source_pairs_path = analysis_dir / "source_pairs.csv.gz"
    wrote_source_header = False
    if write_source_pairs and source_pairs_path.exists():
        source_pairs_path.unlink()

    actual_strata = [s for s in available_strata(pairs) if s != "ALL"]
    for mechanism in MECHANISMS:
        print(f"[source] loading {mechanism} rows", flush=True)
        src = read_sources_for_mechanism(
            source_path, mechanism, chunksize, row_id_keep, charge_source_by_row
        )
        print(f"[source] {mechanism} rows kept: {len(src)}", flush=True)
        if src.empty:
            continue
        geom_cols = [c for c in SOURCE_GEOMETRY_COLUMNS[mechanism] if c in src.columns]
        key_cols = SOURCE_KEY_COLUMNS[mechanism]
        merge_cols = ["pair_id", *key_cols]
        cutoff_col = CUTOFF_COLUMN[mechanism]

        for window in sorted(set(windows)):
            pair_w = pairs[pairs["window"] == window]
            for stratum in actual_strata:
                pair_s = pair_w[pair_w["stratum"] == stratum]
                if pair_s.empty:
                    continue
                t_map = pd.Series(
                    pair_s["pair_id"].to_numpy(np.int64),
                    index=pair_s["row_id_t"].to_numpy(np.int64),
                )
                tp1_map = pd.Series(
                    pair_s["pair_id"].to_numpy(np.int64),
                    index=pair_s["row_id_tp1"].to_numpy(np.int64),
                )
                left = src[src["row_id"].isin(t_map.index)].copy()
                right = src[src["row_id"].isin(tp1_map.index)].copy()
                left["pair_id"] = left["row_id"].map(t_map).astype("int64")
                right["pair_id"] = right["row_id"].map(tp1_map).astype("int64")
                dup_t = int(left.duplicated(merge_cols).sum())
                dup_tp1 = int(right.duplicated(merge_cols).sum())
                if dup_t:
                    left = left.drop_duplicates(merge_cols, keep="first")
                if dup_tp1:
                    right = right.drop_duplicates(merge_cols, keep="first")
                merged = left.merge(
                    right,
                    on=merge_cols,
                    how="outer",
                    suffixes=("_t", "_tp1"),
                    indicator=True,
                )
                if merged.empty:
                    continue
                meta_cols = [
                    "pair_id",
                    "window",
                    "stratum",
                    "d_dft_sigma_iso",
                    "dnorm_target_local_T2",
                    f"{cutoff_col}_t",
                ]
                meta = pair_s[meta_cols].rename(columns={f"{cutoff_col}_t": "cutoff_A"})
                merged = merged.merge(meta, on="pair_id", how="left", validate="many_to_one")
                merged["mechanism"] = mechanism

                for label in (stratum, "ALL"):
                    for cutoff, g in merged.groupby("cutoff_A", dropna=False):
                        both = int((g["_merge"] == "both").sum())
                        t_only = int((g["_merge"] == "left_only").sum())
                        tp1_only = int((g["_merge"] == "right_only").sum())
                        ckey = (window, label, mechanism, cutoff)
                        update_churn_counter(
                            churn,
                            ckey,
                            both,
                            t_only,
                            tp1_only,
                            g["pair_id"].to_numpy(),
                            len(pair_s),
                        )
                        churn[ckey]["duplicates_t"] += dup_t
                        churn[ckey]["duplicates_tp1"] += dup_tp1

                both = merged[merged["_merge"] == "both"].copy()
                if both.empty:
                    continue
                targets = [
                    ("dft_sigma_iso", "d_dft_sigma_iso"),
                    ("target_local_T2_delta_norm", "dnorm_target_local_T2"),
                ]
                for label in (stratum, "ALL"):
                    g = both if label == stratum else both
                    for geom in geom_cols:
                        lt = f"{geom}_t"
                        rt = f"{geom}_tp1"
                        if lt not in g.columns or rt not in g.columns:
                            continue
                        dgeom = g[rt].to_numpy(np.float64) - g[lt].to_numpy(np.float64)
                        lag_acc[(window, label, mechanism, geom)].add(
                            g[lt].to_numpy(np.float64), g[rt].to_numpy(np.float64)
                        )
                        for target_name, target_col in targets:
                            corr_acc[
                                (window, label, mechanism, target_name, geom)
                            ].add(dgeom, g[target_col].to_numpy(np.float64))

                    mat_key = (window, label, mechanism)
                    if mat_key not in level_mat:
                        level_mat[mat_key] = MatrixAccumulator(geom_cols)
                        diff_mat[mat_key] = MatrixAccumulator(geom_cols)
                    level_mat[mat_key].add(
                        g[[f"{c}_t" for c in geom_cols]].to_numpy(np.float64)
                    )
                    diff_mat[mat_key].add(
                        (
                            g[[f"{c}_tp1" for c in geom_cols]].to_numpy(np.float64)
                            - g[[f"{c}_t" for c in geom_cols]].to_numpy(np.float64)
                        )
                    )

                if write_source_pairs:
                    frame = source_pair_output_frame(
                        merged,
                        mechanism,
                        geom_cols,
                        ["d_dft_sigma_iso", "dnorm_target_local_T2"],
                        key_cols,
                    )
                    append_csv_gz(source_pairs_path, frame, not wrote_source_header)
                    wrote_source_header = True

        del src

    corr_rows = []
    for (window, stratum, mechanism, target, geom), acc in sorted(corr_acc.items()):
        r = acc.corr()
        corr_rows.append(
            {
                "window": window,
                "stratum": stratum,
                "mechanism": mechanism,
                "target": target,
                "source_geometry": geom,
                "n_source_pairs": acc.n,
                "pearson_r": r,
                "abs_pearson_r": abs_or_nan(r),
                "min_n_pass": bool(acc.n >= min_n),
            }
        )
    lag_rows = []
    for (window, stratum, mechanism, geom), acc in sorted(lag_acc.items()):
        rho = acc.rho()
        lag_rows.append(
            {
                "window": window,
                "stratum": stratum,
                "mechanism": mechanism,
                "source_geometry": geom,
                "n_source_pairs": acc.n,
                "lag_corr_t_to_tp1": rho,
                "noise_cost": acc.noise_cost(),
                "smooth_pass": bool(np.isfinite(rho) and rho >= smooth_rho_min),
            }
        )
    decor_rows = []
    for key in sorted(level_mat):
        window, stratum, mechanism = key
        lev = level_mat[key].stats()
        dif = diff_mat[key].stats()
        lev_med = float(lev["median_abs_corr"])
        dif_med = float(dif["median_abs_corr"])
        decor_rows.append(
            {
                "window": window,
                "stratum": stratum,
                "mechanism": mechanism,
                "level_n_rows": lev["n_rows"],
                "diff_n_rows": dif["n_rows"],
                "level_n_features_active": lev["n_features_active"],
                "diff_n_features_active": dif["n_features_active"],
                "level_max_abs_corr": lev["max_abs_corr"],
                "diff_max_abs_corr": dif["max_abs_corr"],
                "level_median_abs_corr": lev_med,
                "diff_median_abs_corr": dif_med,
                "median_abs_corr_drop": lev_med - dif_med
                if np.isfinite(lev_med) and np.isfinite(dif_med)
                else math.nan,
                "decorrelation_pass": bool(
                    np.isfinite(lev_med)
                    and np.isfinite(dif_med)
                    and (lev_med - dif_med) >= drop_min
                ),
            }
        )
    churn_rows = []
    for (window, stratum, mechanism, cutoff), values in sorted(churn.items()):
        union = values["union_sources"]
        churn_count = values["t_only"] + values["tp1_only"]
        churn_rows.append(
            {
                "window": window,
                "stratum": stratum,
                "mechanism": mechanism,
                "cutoff_A": cutoff,
                **values,
                "churn_count": churn_count,
                "churn_rate": churn_count / union if union else math.nan,
                "entry_rate": values["tp1_only"] / union if union else math.nan,
                "exit_rate": values["t_only"] / union if union else math.nan,
            }
        )
    return (
        pd.DataFrame(corr_rows),
        pd.DataFrame(lag_rows),
        pd.DataFrame(decor_rows),
        pd.DataFrame(churn_rows),
    )


def markdown_table(df: pd.DataFrame, columns: list[str], max_rows: int = 12) -> str:
    if df.empty:
        return "_no rows_\n"
    d = df.loc[:, columns].head(max_rows).copy()
    for col in d.columns:
        if pd.api.types.is_float_dtype(d[col]):
            d[col] = [fmt_float(float(x)) for x in d[col]]
    rendered = d.astype(str)
    widths = {
        col: max(len(str(col)), *(len(x) for x in rendered[col].tolist()))
        for col in rendered.columns
    }
    header = "| " + " | ".join(str(col).ljust(widths[col]) for col in rendered.columns) + " |"
    sep = "| " + " | ".join("-" * widths[col] for col in rendered.columns) + " |"
    body = [
        "| "
        + " | ".join(str(row[col]).ljust(widths[col]) for col in rendered.columns)
        + " |"
        for _, row in rendered.iterrows()
    ]
    return "\n".join([header, sep, *body]) + "\n"


def write_report(
    path: Path,
    aggregate_corr: pd.DataFrame,
    aggregate_gate: pd.DataFrame,
    aggregate_decor: pd.DataFrame,
    source_corr: pd.DataFrame | None,
    source_gate: pd.DataFrame | None,
    source_decor: pd.DataFrame | None,
    source_churn: pd.DataFrame | None,
    pairs: pd.DataFrame,
    out_dir: Path,
) -> None:
    lines = []
    lines.append("# Differencing System-ID Evidence\n")
    lines.append(
        "Framing: geometric perturbation probe over emitted columns. "
        "The report gives correlations, gate diagnostics, effective N, and churn. "
        "It does not declare an overall system-ID verdict.\n"
    )
    lines.append(f"Substrate: `{out_dir}`\n")
    lines.append(
        f"Atom pairs: {len(pairs)} rows across windows "
        f"{sorted(pairs.window.unique().tolist())}; atoms={pairs.atom_index.nunique()}.\n"
    )

    lines.append("## Autocorr / Noise Gate\n")
    gate = aggregate_gate[
        (aggregate_gate["scope"] == "target")
        & (aggregate_gate["feature"].isin(["dft_sigma_iso", "levelnorm_target_local_T2"]))
    ].copy()
    gate = gate.sort_values(["window", "stratum", "feature"])
    lines.append(
        markdown_table(
            gate,
            [
                "window",
                "stratum",
                "feature",
                "n_pairs",
                "effective_n",
                "median_lag1_rho",
                "median_atom_window_noise_cost",
                "smooth_pass",
            ],
            max_rows=30,
        )
    )

    lines.append("## Aggregate Delta Correlations\n")
    top = aggregate_corr.dropna(subset=["abs_pearson_r"]).sort_values(
        "abs_pearson_r", ascending=False
    )
    lines.append(
        markdown_table(
            top,
            [
                "window",
                "stratum",
                "mechanism",
                "feature",
                "feature_metric",
                "target",
                "n_pairs",
                "pearson_r",
            ],
            max_rows=20,
        )
    )
    ca = aggregate_corr[
        (aggregate_corr["stratum"] == "CA")
        & np.isfinite(aggregate_corr["abs_pearson_r"])
    ].sort_values("abs_pearson_r", ascending=False)
    lines.append("### CA Aggregate Focus\n")
    lines.append(
        markdown_table(
            ca,
            [
                "window",
                "mechanism",
                "feature",
                "feature_metric",
                "target",
                "n_pairs",
                "pearson_r",
            ],
            max_rows=15,
        )
    )

    lines.append("## Aggregate Decorrelation\n")
    dec = (
        aggregate_decor.sort_values(["window", "stratum", "feature_set"])
        if not aggregate_decor.empty
        else aggregate_decor
    )
    lines.append(
        markdown_table(
            dec,
            [
                "window",
                "stratum",
                "feature_set",
                "level_median_abs_corr",
                "diff_median_abs_corr",
                "median_abs_corr_drop",
                "decorrelation_pass",
            ],
            max_rows=30,
        )
    )

    if source_corr is not None and not source_corr.empty:
        lines.append("## Aligned Source Delta Geometry Correlations\n")
        src_top = source_corr.dropna(subset=["abs_pearson_r"]).sort_values(
            "abs_pearson_r", ascending=False
        )
        lines.append(
            markdown_table(
                src_top,
                [
                    "window",
                    "stratum",
                    "mechanism",
                    "source_geometry",
                    "target",
                    "n_source_pairs",
                    "pearson_r",
                ],
                max_rows=20,
            )
        )
        src_ca = source_corr[
            (source_corr["stratum"] == "CA")
            & np.isfinite(source_corr["abs_pearson_r"])
        ].sort_values("abs_pearson_r", ascending=False)
        lines.append("### CA Source Focus\n")
        lines.append(
            markdown_table(
                src_ca,
                [
                    "window",
                    "mechanism",
                    "source_geometry",
                    "target",
                    "n_source_pairs",
                    "pearson_r",
                ],
                max_rows=15,
            )
        )

    if source_gate is not None and not source_gate.empty:
        lines.append("## Source Geometry Autocorr Gate\n")
        sg = source_gate.sort_values(["window", "stratum", "mechanism", "source_geometry"])
        lines.append(
            markdown_table(
                sg,
                [
                    "window",
                    "stratum",
                    "mechanism",
                    "source_geometry",
                    "n_source_pairs",
                    "lag_corr_t_to_tp1",
                    "noise_cost",
                    "smooth_pass",
                ],
                max_rows=30,
            )
        )

    if source_decor is not None and not source_decor.empty:
        lines.append("## Source Geometry Decorrelation\n")
        sd = source_decor.sort_values(["window", "stratum", "mechanism"])
        lines.append(
            markdown_table(
                sd,
                [
                    "window",
                    "stratum",
                    "mechanism",
                    "level_median_abs_corr",
                    "diff_median_abs_corr",
                    "median_abs_corr_drop",
                    "decorrelation_pass",
                ],
                max_rows=30,
            )
        )

    if source_churn is not None and not source_churn.empty:
        lines.append("## Source Churn\n")
        churn = source_churn.sort_values(["window", "stratum", "mechanism"])
        lines.append(
            markdown_table(
                churn,
                [
                    "window",
                    "stratum",
                    "mechanism",
                    "cutoff_A",
                    "present_both",
                    "t_only",
                    "tp1_only",
                    "churn_rate",
                ],
                max_rows=36,
            )
        )

    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    analysis_dir = (
        Path(args.analysis_dir)
        if args.analysis_dir is not None
        else out_dir / "differencing_system_id"
    )
    analysis_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] aggregate from {out_dir}", flush=True)
    agg = load_aggregate(out_dir)
    required_features = sorted(
        set(c for cols in AGG_FEATURES_BY_MECH.values() for c in cols)
        | set(CUTOFF_COLUMN.values())
    )
    require_columns(agg, required_features, "aggregate CSV")

    print(f"[pairs] building windows {args.windows}", flush=True)
    pairs = build_atom_pairs(agg, args.windows)
    atom_pairs_path = analysis_dir / "atom_pairs.csv.gz"
    print(f"[pairs] writing compact atom pair table to {atom_pairs_path}", flush=True)
    pairs.loc[:, atom_pair_output_columns(pairs)].to_csv(
        atom_pairs_path, index=False, compression="gzip"
    )
    print(f"[pairs] wrote {atom_pairs_path} ({len(pairs)} rows in memory)", flush=True)

    print("[aggregate] gate/decorrelation/correlation summaries", flush=True)
    aggregate_gate = aggregate_gate_summary(
        agg, pairs, args.windows, args.smooth_rho_min
    )
    aggregate_decor = aggregate_decorrelation(
        pairs, args.windows, args.decorrelation_drop_min
    )
    aggregate_corr = aggregate_delta_correlations(pairs, args.windows, args.min_n)
    aggregate_gate.to_csv(analysis_dir / "aggregate_gate_summary.csv", index=False)
    aggregate_decor.to_csv(analysis_dir / "aggregate_decorrelation.csv", index=False)
    aggregate_corr.to_csv(analysis_dir / "aggregate_delta_correlations.csv", index=False)

    source_corr = source_gate = source_decor = source_churn = None
    if not args.skip_source:
        charge_source_by_row = pd.Series(
            agg["charge_source_code"].to_numpy(np.int64),
            index=agg["row_id"].to_numpy(np.int64),
        )
        print("[source] aligning stable source keys", flush=True)
        source_corr, source_gate, source_decor, source_churn = process_sources(
            out_dir,
            analysis_dir,
            pairs,
            args.windows,
            args.chunksize,
            args.min_n,
            args.smooth_rho_min,
            args.decorrelation_drop_min,
            args.write_source_pairs,
            charge_source_by_row,
        )
        source_corr.to_csv(analysis_dir / "source_delta_correlations.csv", index=False)
        source_gate.to_csv(analysis_dir / "source_feature_gate.csv", index=False)
        source_decor.to_csv(analysis_dir / "source_decorrelation.csv", index=False)
        source_churn.to_csv(analysis_dir / "source_churn.csv", index=False)

    report_path = analysis_dir / "differencing_system_id_report.md"
    write_report(
        report_path,
        aggregate_corr,
        aggregate_gate,
        aggregate_decor,
        source_corr,
        source_gate,
        source_decor,
        source_churn,
        pairs,
        out_dir,
    )
    print(f"[report] wrote {report_path}", flush=True)
    print(
        "[done] evidence only: correlations, gate labels, effective N, and churn; no verdict emitted",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
