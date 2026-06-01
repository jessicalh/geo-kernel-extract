#!/usr/bin/env python3
"""Region-partitioned Johnson-Bovey recovery on aromatic ring-facing H.

Read-only Python consumer:
  * reads the C++-emitted ring_current source substrate;
  * groups emitted source columns by atom/frame and geometry region;
  * correlates those emitted features against emitted DFT targets.

It does not open trajectory.h5, evaluate the Johnson-Bovey field, compute
`(3cos^2-1)/r^3`, or project tensors. The J-B/BS source T0/T2 columns and the
DFT local T2 sidecar are emitted by the C++ spine.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Region:
    name: str
    description: str
    mask: np.ndarray


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out-dir",
        default="/tmp/rediscover-jb-parity/composed",
        help="directory containing C++-emitted ring_current CSV/NPY files",
    )
    ap.add_argument(
        "--artifact-dir",
        default="/tmp/rediscover-jb-region-recovery",
        help="directory for CSV/JSON artifacts",
    )
    ap.add_argument(
        "--report-md",
        default="src/rediscover/analysis/JOHNSON_BOVEY_REGION_RECOVERY.md",
        help="markdown report path",
    )
    ap.add_argument("--in-plane-abs-cos-max", type=float, default=0.20)
    ap.add_argument("--axial-abs-cos-min", type=float, default=0.70)
    ap.add_argument("--far-r-min", type=float, default=4.0)
    return ap.parse_args()


def require_columns(df: pd.DataFrame, cols: list[str], label: str) -> None:
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
        a = np.asarray(arr)
        ok = np.isfinite(a).all(axis=1) if a.ndim > 1 else np.isfinite(a)
        mask = ok if mask is None else (mask & ok)
    if mask is None:
        raise RuntimeError("finite_rows requires at least one array")
    return mask


def centre_by_atom(x: np.ndarray, y: np.ndarray, atoms: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = as_matrix(x)
    y = as_matrix(y)
    atoms = np.asarray(atoms, dtype=int)
    xw = np.full_like(x, np.nan, dtype=float)
    yw = np.full_like(y, np.nan, dtype=float)
    for atom in np.sort(np.unique(atoms)):
        m = atoms == atom
        if m.sum() == 0:
            continue
        xw[m] = x[m] - x[m].mean(axis=0, keepdims=True)
        yw[m] = y[m] - y[m].mean(axis=0, keepdims=True)
    return xw, yw


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


def gamma_fit(x: np.ndarray, y: np.ndarray, mask: np.ndarray | None = None) -> float:
    x = as_matrix(x)
    y = as_matrix(y)
    if mask is None:
        mask = finite_rows(x, y)
    else:
        mask = np.asarray(mask, dtype=bool) & finite_rows(x, y)
    if mask.sum() < 3:
        return math.nan
    den = float((x[mask] * x[mask]).sum())
    if den <= 0.0:
        return math.nan
    return float((x[mask] * y[mask]).sum() / den)


def loao_scaled_r2(x: np.ndarray, y: np.ndarray, atoms: np.ndarray, mask: np.ndarray) -> float:
    x = as_matrix(x)
    y = as_matrix(y)
    pred = np.full_like(y, np.nan, dtype=float)
    labels = np.sort(np.unique(atoms[mask]))
    if labels.size <= 3:
        return math.nan
    for atom in labels:
        score = mask & (atoms == atom)
        train = mask & (atoms != atom)
        gamma = gamma_fit(x, y, train)
        if np.isfinite(gamma):
            pred[score] = gamma * x[score]
    return r2_score(pred[mask], y[mask])


def metric_value(kind: str, x: np.ndarray, y: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    x = as_matrix(x)
    y = as_matrix(y)
    mask = np.asarray(mask, dtype=bool) & finite_rows(x, y)
    if mask.sum() < 3:
        return {
            "pearson_r": math.nan,
            "pearson_R2": math.nan,
            "component_r": math.nan,
            "component_R2": math.nan,
            "absT2_r": math.nan,
        }
    if kind == "T0":
        r = corrcoef(x[mask], y[mask])
        return {
            "pearson_r": r,
            "pearson_R2": r * r if np.isfinite(r) else math.nan,
            "component_r": math.nan,
            "component_R2": math.nan,
            "absT2_r": math.nan,
        }
    comp = corrcoef(x[mask].reshape(-1), y[mask].reshape(-1))
    mag = corrcoef(np.linalg.norm(x[mask], axis=1), np.linalg.norm(y[mask], axis=1))
    return {
        "pearson_r": math.nan,
        "pearson_R2": math.nan,
        "component_r": comp,
        "component_R2": comp * comp if np.isfinite(comp) else math.nan,
        "absT2_r": mag,
    }


def jackknife_metric_se(
    kind: str,
    metric: str,
    x: np.ndarray,
    y: np.ndarray,
    atoms: np.ndarray,
    mask: np.ndarray,
) -> float:
    vals = []
    labels = np.sort(np.unique(atoms[mask]))
    if labels.size <= 3:
        return math.nan
    for atom in labels:
        keep = mask & (atoms != atom)
        vals.append(metric_value(kind, x, y, keep)[metric])
    arr = np.asarray(vals, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < 2:
        return math.nan
    mean = float(arr.mean())
    return float(math.sqrt((arr.size - 1) / arr.size * ((arr - mean) ** 2).sum()))


def jackknife_gamma_se(x: np.ndarray, y: np.ndarray, atoms: np.ndarray, mask: np.ndarray) -> float:
    vals = []
    labels = np.sort(np.unique(atoms[mask]))
    if labels.size <= 3:
        return math.nan
    for atom in labels:
        vals.append(gamma_fit(x, y, mask & (atoms != atom)))
    arr = np.asarray(vals, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < 2:
        return math.nan
    mean = float(arr.mean())
    return float(math.sqrt((arr.size - 1) / arr.size * ((arr - mean) ** 2).sum()))


def signal_effective_n(x: np.ndarray, atoms: np.ndarray, mask: np.ndarray) -> tuple[float, int]:
    x = as_matrix(x)
    weights = []
    for atom in np.sort(np.unique(atoms[mask])):
        m = mask & (atoms == atom)
        if m.sum() < 2:
            continue
        xc = x[m] - x[m].mean(axis=0, keepdims=True)
        weights.append(float((xc * xc).sum()))
    w = np.asarray(weights, dtype=float)
    if w.size == 0 or float((w * w).sum()) <= 0.0:
        return math.nan, 0
    participation = float((w.sum() * w.sum()) / (w * w).sum())
    ordered = np.sort(w)[::-1]
    cum = np.cumsum(ordered)
    top90 = int(np.searchsorted(cum / cum[-1], 0.90) + 1)
    return participation, top90


def frame_effective_n(y: np.ndarray, atoms: np.ndarray, frames: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
    y = as_matrix(y)
    mag = y[:, 0] if y.shape[1] == 1 else np.linalg.norm(y, axis=1)
    total = 0.0
    rhos = []
    for atom in np.sort(np.unique(atoms[mask])):
        m = mask & (atoms == atom)
        order = np.argsort(frames[m])
        vals = mag[m][order]
        vals = vals[np.isfinite(vals)]
        n = vals.size
        if n < 3:
            total += float(n)
            continue
        yd = vals - vals.mean()
        den = float(np.dot(yd, yd))
        if den <= 0.0:
            rho = 0.0
        else:
            rho = float(np.dot(yd[:-1], yd[1:]) / den)
            rho = max(min(rho, 0.999), -0.999)
        neff = n * (1.0 - rho) / (1.0 + rho)
        neff = min(max(neff, 1.0), float(n))
        rhos.append(rho)
        total += neff
    return total, float(np.median(rhos)) if rhos else math.nan


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


def aggregate_region(
    df: pd.DataFrame,
    regions: list[Region],
    n_groups: int,
) -> dict[tuple[str, str, str], np.ndarray]:
    out: dict[tuple[str, str, str], np.ndarray] = {}
    for region in regions:
        m = region.mask
        part = df.loc[m, ["gid", "dipolar_3cos2m1_over_r3", "jb_T0", *[f"jb_T2_local_{i}" for i in range(5)]]]
        counts = np.zeros(n_groups, dtype=float)
        if len(part):
            count_series = part.groupby("gid").size()
            counts[count_series.index.to_numpy(int)] = count_series.to_numpy(float)
        out[(region.name, "source_count", "count")] = counts.reshape(-1, 1)

        point = np.zeros(n_groups, dtype=float)
        if len(part):
            s = part.groupby("gid")["dipolar_3cos2m1_over_r3"].sum()
            point[s.index.to_numpy(int)] = s.to_numpy(float)
        out[(region.name, "point_dipole", "T0")] = point.reshape(-1, 1)

        jb0 = np.zeros(n_groups, dtype=float)
        if len(part):
            s = part.groupby("gid")["jb_T0"].sum()
            jb0[s.index.to_numpy(int)] = s.to_numpy(float)
        out[(region.name, "johnson_bovey", "T0")] = jb0.reshape(-1, 1)

        jb2 = np.zeros((n_groups, 5), dtype=float)
        if len(part):
            s2 = part.groupby("gid")[[f"jb_T2_local_{i}" for i in range(5)]].sum()
            jb2[s2.index.to_numpy(int)] = s2.to_numpy(float)
        out[(region.name, "johnson_bovey", "T2")] = jb2
    return out


def summarize(
    rows: list[dict[str, object]],
    region: Region,
    model: str,
    target_kind: str,
    x: np.ndarray,
    y: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    valid: np.ndarray,
    source_counts: np.ndarray,
) -> None:
    x = as_matrix(x)
    y = as_matrix(y)
    valid = valid & finite_rows(x, y)
    xw, yw = centre_by_atom(x, y, atoms)
    score = valid & finite_rows(xw, yw)
    metrics = metric_value(target_kind, xw, yw, score)
    gamma = gamma_fit(xw, yw, score)
    gamma_se = jackknife_gamma_se(xw, yw, atoms, score)
    metric_for_se = "pearson_R2" if target_kind == "T0" else "component_r"
    metric_se = jackknife_metric_se(target_kind, metric_for_se, xw, yw, atoms, score)
    abs_se = (
        jackknife_metric_se(target_kind, "absT2_r", xw, yw, atoms, score)
        if target_kind == "T2"
        else math.nan
    )
    atom_neff, atom_top90 = signal_effective_n(xw, atoms, score)
    frame_neff, lag1 = frame_effective_n(yw, atoms, frames, score)
    rows.append(
        {
            "stratum": "aromatic_H_ring_current",
            "region": region.name,
            "region_definition": region.description,
            "model": model,
            "target": target_kind,
            "rows_scored": int(score.sum()),
            "atoms_total": int(np.unique(atoms[score]).size),
            "source_rows": int(np.nansum(source_counts[score])),
            "groups_with_sources": int(np.sum(source_counts[score] > 0)),
            "atom_signal_neff": atom_neff,
            "atom_signal_top90_count": atom_top90,
            "frame_effective_n_sum": frame_neff,
            "target_median_lag1": lag1,
            "pearson_r": metrics["pearson_r"],
            "pearson_R2": metrics["pearson_R2"],
            "pearson_R2_jackknife_se": metric_se if target_kind == "T0" else math.nan,
            "component_r": metrics["component_r"],
            "component_r_jackknife_se": metric_se if target_kind == "T2" else math.nan,
            "component_R2": metrics["component_R2"],
            "absT2_r": metrics["absT2_r"],
            "absT2_r_jackknife_se": abs_se,
            "gamma": gamma,
            "gamma_jackknife_se": gamma_se,
            "loao_gamma_scaled_R2": loao_scaled_r2(xw, yw, atoms, score),
            "notes": (
                "point-dipole is emitted scalar sum_dipolar; T2 is not defined for that scalar"
                if model == "point_dipole"
                else "J-B/BS source kernel emitted by C++ with fixed literature G-P intensity"
            ),
        }
    )


def write_csv(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "stratum",
        "region",
        "region_definition",
        "model",
        "target",
        "rows_scored",
        "atoms_total",
        "source_rows",
        "groups_with_sources",
        "atom_signal_neff",
        "atom_signal_top90_count",
        "frame_effective_n_sum",
        "target_median_lag1",
        "pearson_r",
        "pearson_R2",
        "pearson_R2_jackknife_se",
        "component_r",
        "component_r_jackknife_se",
        "component_R2",
        "absT2_r",
        "absT2_r_jackknife_se",
        "gamma",
        "gamma_jackknife_se",
        "loao_gamma_scaled_R2",
        "notes",
    ]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def write_report(rows: list[dict[str, object]], audit: dict[str, object], path: Path, csv_path: Path, out_dir: Path) -> None:
    lines = [
        "# Johnson-Bovey Region Recovery",
        "",
        f"Substrate: `{out_dir}`",
        f"CSV artifact: `{csv_path}`",
        "",
        "Read-only analysis of the aromatic ring-facing H stratum. Point-dipole rows use the emitted scalar `sum_dipolar` source column; J-B rows use the emitted source-level `jb_T0` / `jb_T2_local_*` columns computed in C++ from the Johnson-Bovey double-loop with literature G-P intensities.",
        "",
        "Geometry regions are parameterized and source-level: `all_valid` excludes self/bonded rings; `in_plane_band` is `|cos_theta| <= 0.2`; `near_field_r_lt4` is `r < 4 A`; `axial_far` is `r >= 4 A` and `|cos_theta| >= 0.7`. The H5 `ring_in_plane_angle` is retained as emitted geometry but is an in-plane azimuth, so region membership uses `cos_theta` and `r`.",
        "",
        "All scores are within-atom centered. SEs are delete-atom jackknife within this one protein. `loao_gamma_scaled_R2` fits only a scalar gamma on all other atoms before scoring each held-out atom; direct correlation columns are unfitted.",
        "",
        "## Results",
        "",
    ]
    table_cols = [
        "region",
        "model",
        "target",
        "atoms_total",
        "atom_signal_neff",
        "source_rows",
        "pearson_r",
        "pearson_R2",
        "pearson_R2_jackknife_se",
        "component_r",
        "component_r_jackknife_se",
        "absT2_r",
        "absT2_r_jackknife_se",
        "gamma",
        "gamma_jackknife_se",
        "loao_gamma_scaled_R2",
    ]
    lines.append("| " + " | ".join(table_cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(table_cols)) + " |")
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(c)) for c in table_cols) + " |")

    by_key = {(r["region"], r["model"], r["target"]): r for r in rows}
    ip_point = by_key.get(("in_plane_band", "point_dipole", "T0"), {})
    ip_jb0 = by_key.get(("in_plane_band", "johnson_bovey", "T0"), {})
    ip_jb2 = by_key.get(("in_plane_band", "johnson_bovey", "T2"), {})
    nf_point = by_key.get(("near_field_r_lt4", "point_dipole", "T0"), {})
    nf_jb0 = by_key.get(("near_field_r_lt4", "johnson_bovey", "T0"), {})
    nf_jb2 = by_key.get(("near_field_r_lt4", "johnson_bovey", "T2"), {})
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"- In-plane band: point-dipole T0 R2={fmt(ip_point.get('pearson_R2'))} +/- {fmt(ip_point.get('pearson_R2_jackknife_se'))}; fixed J-B T0 R2={fmt(ip_jb0.get('pearson_R2'))} +/- {fmt(ip_jb0.get('pearson_R2_jackknife_se'))}; fixed J-B T2 component r={fmt(ip_jb2.get('component_r'))} +/- {fmt(ip_jb2.get('component_r_jackknife_se'))}. This does not show J-B recovering the in-plane DFT modulation that the scalar point-dipole misses.",
            f"- r < 4 A near-field diagnostic: point-dipole T0 R2={fmt(nf_point.get('pearson_R2'))} +/- {fmt(nf_point.get('pearson_R2_jackknife_se'))}; fixed J-B T0 R2={fmt(nf_jb0.get('pearson_R2'))} +/- {fmt(nf_jb0.get('pearson_R2_jackknife_se'))}; fixed J-B T2 component r={fmt(nf_jb2.get('component_r'))} +/- {fmt(nf_jb2.get('component_r_jackknife_se'))}, |T2| r={fmt(nf_jb2.get('absT2_r'))} +/- {fmt(nf_jb2.get('absT2_r_jackknife_se'))}. This region is thin, with atom-signal N_eff={fmt(nf_jb2.get('atom_signal_neff'))}.",
            "- The all-valid rows still carry the main ring signal, especially J-B T2 magnitude, but the in-plane source partition alone is near-null against DFT in this single-protein substrate. Verdict remains reserved for the lead.",
        ]
    )

    lines.extend(
        [
            "",
            "## Self Audit",
            "",
            f"- Source-level emitted unit-current J-B T0 valid-sum vs emitted H5 `bare_T0`: r={fmt(audit['jb_valid_vs_bare_T0_r'])}, RMS={fmt(audit['jb_valid_vs_bare_T0_rms'])} ppm_T_per_nA, max_abs={fmt(audit['jb_valid_vs_bare_T0_max_abs'])} ppm_T_per_nA.",
            f"- Manual emitted-row sum check for one in-plane atom/frame: atom={audit['manual_atom_index']}, h5_row={audit['manual_h5_row']}, source_rows={audit['manual_source_rows']}, manual_jb_T0={fmt(audit['manual_jb_T0'], 8)}, grouped_jb_T0={fmt(audit['grouped_jb_T0'], 8)}, abs_diff={fmt(audit['manual_abs_diff'], 8)}.",
            "- Python did not open trajectory.h5, evaluate Biot-Savart/J-B, compute `(3cos^2-1)/r^3`, or project tensors; it summed emitted source columns and read emitted target sidecars.",
            "- Point-dipole T2 is not reported because `sum_dipolar` is a scalar T0 proxy, not an emitted tensor law.",
            "",
            "## De-Circularised Read",
            "",
            "The de-circularised test is the direct J-B correlation/gamma rows: the kernel uses fixed literature G-P intensities in C++ and no fitted scale for the reported correlation. Gamma and LOAO scaled R2 are diagnostics, not the verdict. Scientific verdict remains reserved for the lead.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    src_path = out_dir / "ring_current_sources.csv"
    target_t2_path = out_dir / "rediscover_ring_current_sources_target_local_T2.npy"
    df = pd.read_csv(src_path)
    require_columns(
        df,
        [
            "atom_index",
            "h5_row",
            "stratum",
            "dft_present",
            "dft_local_frame_valid",
            "dft_sigma_iso",
            "bare_T0",
            "is_self_or_bonded",
            "r",
            "cos_theta",
            "ring_in_plane_angle",
            "dipolar_3cos2m1_over_r3",
            "jb_kernel_present",
            "jb_unit_T0",
            *[f"jb_unit_T2_local_{i}" for i in range(5)],
            "jb_T0",
            *[f"jb_T2_local_{i}" for i in range(5)],
        ],
        "ring_current_sources",
    )
    target_t2_all = np.load(target_t2_path)
    if target_t2_all.shape != (len(df), 5):
        raise RuntimeError(f"target local T2 sidecar shape {target_t2_all.shape} != ({len(df)}, 5)")

    df = df.copy()
    df["gid"], _ = pd.factorize(df["atom_index"].astype(str) + ":" + df["h5_row"].astype(str), sort=True)
    n_groups = int(df["gid"].max()) + 1
    first = df.sort_values("gid").groupby("gid", sort=True).head(1).copy()
    first_rowid = first.index.to_numpy(int)
    atoms = first["atom_index"].to_numpy(int)
    frames = first["h5_row"].to_numpy(int)
    y0 = first["dft_sigma_iso"].to_numpy(float).reshape(-1, 1)
    y2 = target_t2_all[first_rowid]
    valid0 = (first["dft_present"].to_numpy(int) == 1) & finite_rows(y0)
    valid2 = (
        (first["dft_present"].to_numpy(int) == 1)
        & (first["dft_local_frame_valid"].to_numpy(int) == 1)
        & finite_rows(y2)
    )

    valid_source = (df["is_self_or_bonded"].to_numpy(int) == 0) & (
        df["jb_kernel_present"].to_numpy(int) == 1
    )
    cos_abs = np.abs(df["cos_theta"].to_numpy(float))
    r = df["r"].to_numpy(float)
    regions = [
        Region("all_valid", "producer-valid sources; self/bonded rings excluded", valid_source),
        Region(
            "in_plane_band",
            f"producer-valid and |cos_theta| <= {args.in_plane_abs_cos_max:g}",
            valid_source & (cos_abs <= args.in_plane_abs_cos_max),
        ),
        Region(
            "near_field_r_lt4",
            f"producer-valid and r < {args.far_r_min:g} A",
            valid_source & (r < args.far_r_min),
        ),
        Region(
            "axial_far",
            f"producer-valid, r >= {args.far_r_min:g} A, |cos_theta| >= {args.axial_abs_cos_min:g}",
            valid_source & (r >= args.far_r_min) & (cos_abs >= args.axial_abs_cos_min),
        ),
    ]
    agg = aggregate_region(df, regions, n_groups)

    rows: list[dict[str, object]] = []
    for region in regions:
        counts = agg[(region.name, "source_count", "count")].reshape(-1)
        summarize(
            rows,
            region,
            "point_dipole",
            "T0",
            agg[(region.name, "point_dipole", "T0")],
            y0,
            atoms,
            frames,
            valid0,
            counts,
        )
        summarize(
            rows,
            region,
            "johnson_bovey",
            "T0",
            agg[(region.name, "johnson_bovey", "T0")],
            y0,
            atoms,
            frames,
            valid0,
            counts,
        )
        summarize(
            rows,
            region,
            "johnson_bovey",
            "T2",
            agg[(region.name, "johnson_bovey", "T2")],
            y2,
            atoms,
            frames,
            valid2,
            counts,
        )

    unit_part = df.loc[valid_source, ["gid", "jb_unit_T0"]]
    all_jb0 = np.zeros(n_groups, dtype=float)
    if len(unit_part):
        unit_sum = unit_part.groupby("gid")["jb_unit_T0"].sum()
        all_jb0[unit_sum.index.to_numpy(int)] = unit_sum.to_numpy(float)
    bare0 = first["bare_T0"].to_numpy(float)
    ok = valid0 & finite_rows(all_jb0, bare0)
    diff = all_jb0[ok] - bare0[ok]
    sample_region = regions[1]
    inplane_rows = df.loc[sample_region.mask].copy()
    if inplane_rows.empty:
        manual = {
            "manual_atom_index": -1,
            "manual_h5_row": -1,
            "manual_source_rows": 0,
            "manual_jb_T0": math.nan,
            "grouped_jb_T0": math.nan,
            "manual_abs_diff": math.nan,
        }
    else:
        gid = int(inplane_rows.iloc[0]["gid"])
        chunk = inplane_rows[inplane_rows["gid"] == gid]
        manual_sum = float(chunk["jb_T0"].sum())
        grouped_sum = float(agg[("in_plane_band", "johnson_bovey", "T0")][gid, 0])
        manual = {
            "manual_atom_index": int(chunk.iloc[0]["atom_index"]),
            "manual_h5_row": int(chunk.iloc[0]["h5_row"]),
            "manual_source_rows": int(len(chunk)),
            "manual_jb_T0": manual_sum,
            "grouped_jb_T0": grouped_sum,
            "manual_abs_diff": abs(manual_sum - grouped_sum),
        }
    audit = {
        "jb_valid_vs_bare_T0_r": corrcoef(all_jb0[ok], bare0[ok]),
        "jb_valid_vs_bare_T0_rms": float(math.sqrt(float(np.mean(diff * diff)))) if diff.size else math.nan,
        "jb_valid_vs_bare_T0_max_abs": float(np.max(np.abs(diff))) if diff.size else math.nan,
        **manual,
    }

    artifact_dir = Path(args.artifact_dir)
    csv_path = artifact_dir / "johnson_bovey_region_recovery.csv"
    audit_path = artifact_dir / "johnson_bovey_region_recovery_audit.json"
    write_csv(rows, csv_path)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True), encoding="utf-8")
    write_report(rows, audit, Path(args.report_md), csv_path, out_dir)

    print("johnson-bovey region recovery")
    print(f"substrate={out_dir}")
    print(f"csv={csv_path}")
    print(f"report={args.report_md}")
    for row in rows:
        if row["target"] == "T0":
            metric = (
                f"r={fmt(row['pearson_r'])} R2={fmt(row['pearson_R2'])}"
                f"+/-{fmt(row['pearson_R2_jackknife_se'])}"
            )
        else:
            metric = (
                f"component_r={fmt(row['component_r'])}+/-{fmt(row['component_r_jackknife_se'])} "
                f"absT2_r={fmt(row['absT2_r'])}+/-{fmt(row['absT2_r_jackknife_se'])}"
            )
        print(
            f"{row['region']} {row['model']} {row['target']}: {metric} "
            f"gamma={fmt(row['gamma'])}+/-{fmt(row['gamma_jackknife_se'])} "
            f"atom_Neff={fmt(row['atom_signal_neff'])} sources={row['source_rows']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
