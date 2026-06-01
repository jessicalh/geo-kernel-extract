#!/usr/bin/env python3
"""Fixed-literature de-circularisation for ring-current and McConnell home strata.

This is a read-only Python consumer. It reads the already-emitted C++ substrate:
CSV scalar producer columns plus row-aligned fixed-literature T2 NPY sidecars.
It does not open H5, rebuild kernels, or re-emit any data.

Targets:
  ring_current: aromatic ring-facing H stratum
  mcconnell: backbone HN / bond stratum

The scalar and tensor comparisons use within-atom centering to remove the
atom-specific static DFT baseline. The fixed-kernel correlation is reported
without fitting; gamma is a separate no-intercept scale diagnostic.
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

import change_of_basis as cob


TARGETS = {
    "ring_current": {
        "stratum": "aromatic_H_ring_current",
        "csv": "ring_current_aggregated.csv",
        "t2_kernel": "rediscover_ring_current_aggregated_bare_kernel_T2.npy",
        "t2_target": "rediscover_ring_current_aggregated_target_local_T2.npy",
        "notes": "ring-current fixed Giessner-Prettre/Pople-style emitted kernel",
    },
    "mcconnell": {
        "stratum": "backbone_HN_bond_mcconnell",
        "csv": "mcconnell_aggregated.csv",
        "t2_kernel": "rediscover_mcconnell_aggregated_bare_kernel_T2.npy",
        "t2_target": "rediscover_mcconnell_aggregated_target_local_T2.npy",
        "notes": "McConnell fixed literature-chi emitted kernel",
    },
}


@dataclass(frozen=True)
class Split:
    train: np.ndarray
    test: np.ndarray
    train_frames: int
    test_frames: int


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Read already-emitted fixed-literature ring/mc kernels vs DFT."
    )
    ap.add_argument(
        "--out-dir",
        default="/tmp/rdc-ring-mc-capstone",
        help="directory containing ring_current/mcconnell CSV and NPY sidecars",
    )
    ap.add_argument(
        "--artifact-dir",
        default="/tmp/rediscover-literature-fixed-decirc",
        help="directory for CSV/JSON analysis artifacts",
    )
    ap.add_argument(
        "--report-md",
        default="src/rediscover/analysis/LITERATURE_FIXED_DECIRC.md",
        help="markdown report path",
    )
    ap.add_argument("--test-fraction", type=float, default=0.30)
    ap.add_argument("--seed", type=int, default=0)
    return ap.parse_args()


def require_columns(df: pd.DataFrame, cols: list[str], label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise RuntimeError(f"{label} missing required columns: {missing}")


def finite_rows(*arrays: np.ndarray) -> np.ndarray:
    mask = None
    for arr in arrays:
        a = np.asarray(arr)
        ok = np.isfinite(a).all(axis=1) if a.ndim > 1 else np.isfinite(a)
        mask = ok if mask is None else (mask & ok)
    if mask is None:
        raise RuntimeError("finite_rows needs at least one array")
    return mask


def make_frame_split(frames: np.ndarray, test_fraction: float, seed: int) -> Split:
    unique = np.sort(np.unique(frames))
    if len(unique) < 4:
        raise RuntimeError("need at least four frames for a frame split")
    rng = np.random.default_rng(seed)
    n_test = max(1, int(round(float(test_fraction) * len(unique))))
    n_test = min(n_test, len(unique) - 1)
    test_frames = set(rng.choice(unique, size=n_test, replace=False).tolist())
    test = np.array([int(f) in test_frames for f in frames], dtype=bool)
    train = ~test
    return Split(
        train=train,
        test=test,
        train_frames=int(train.any() and np.unique(frames[train]).size),
        test_frames=int(test.any() and np.unique(frames[test]).size),
    )


def centre_by_atom(
    x: np.ndarray,
    y: np.ndarray,
    atoms: np.ndarray,
    mean_mask: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    x = as_matrix(x)
    y = as_matrix(y)
    atoms = np.asarray(atoms, dtype=int)
    if mean_mask is None:
        mean_mask = np.ones(len(atoms), dtype=bool)
    xw = np.full_like(x, np.nan, dtype=float)
    yw = np.full_like(y, np.nan, dtype=float)
    for atom in np.sort(np.unique(atoms)):
        all_rows = atoms == atom
        mean_rows = all_rows & mean_mask
        if mean_rows.sum() == 0:
            continue
        xw[all_rows] = x[all_rows] - x[mean_rows].mean(axis=0, keepdims=True)
        yw[all_rows] = y[all_rows] - y[mean_rows].mean(axis=0, keepdims=True)
    return xw, yw


def as_matrix(a: np.ndarray) -> np.ndarray:
    arr = np.asarray(a, dtype=float)
    if arr.ndim == 1:
        return arr.reshape(-1, 1)
    return arr


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
        mask = mask & finite_rows(x, y)
    if mask.sum() < 3:
        return math.nan
    den = float((x[mask] * x[mask]).sum())
    if den <= 0.0:
        return math.nan
    return float((x[mask] * y[mask]).sum() / den)


def atom_jackknife_gamma_se(
    x: np.ndarray,
    y: np.ndarray,
    atoms: np.ndarray,
    base_mask: np.ndarray,
) -> float:
    vals = []
    labels = np.sort(np.unique(atoms[base_mask]))
    if labels.size <= 3:
        return math.nan
    for atom in labels:
        keep = base_mask & (atoms != atom)
        vals.append(gamma_fit(x, y, keep))
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size < 2:
        return math.nan
    mean = vals.mean()
    return float(math.sqrt((vals.size - 1) / vals.size * ((vals - mean) ** 2).sum()))


def tensor_metrics(x: np.ndarray, y: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    x = as_matrix(x)
    y = as_matrix(y)
    mask = mask & finite_rows(x, y)
    if mask.sum() < 3:
        return {
            "component_r": math.nan,
            "absT2_r": math.nan,
            "row_cosine_mean": math.nan,
        }
    xm = x[mask]
    ym = y[mask]
    dot = np.sum(xm * ym, axis=1)
    den = np.linalg.norm(xm, axis=1) * np.linalg.norm(ym, axis=1)
    ok = np.isfinite(dot) & np.isfinite(den) & (den > 0.0)
    row_cos = float(np.mean(dot[ok] / den[ok])) if ok.any() else math.nan
    return {
        "component_r": corrcoef(xm.reshape(-1), ym.reshape(-1)),
        "absT2_r": corrcoef(np.linalg.norm(xm, axis=1), np.linalg.norm(ym, axis=1)),
        "row_cosine_mean": row_cos,
    }


def scalar_metrics(x: np.ndarray, y: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    x = as_matrix(x)
    y = as_matrix(y)
    mask = mask & finite_rows(x, y)
    r = corrcoef(x[mask], y[mask]) if mask.sum() >= 3 else math.nan
    return {
        "pearson_r": r,
        "pearson_R2": float(r * r) if np.isfinite(r) else math.nan,
    }


def loao_scaled_prediction(
    x: np.ndarray,
    y: np.ndarray,
    atoms: np.ndarray,
    mask: np.ndarray,
) -> np.ndarray:
    x = as_matrix(x)
    y = as_matrix(y)
    pred = np.full_like(y, np.nan, dtype=float)
    labels = np.sort(np.unique(atoms[mask]))
    for atom in labels:
        score = mask & (atoms == atom)
        train = mask & (atoms != atom)
        gamma = gamma_fit(x, y, train)
        if np.isfinite(gamma):
            pred[score] = gamma * x[score]
    return pred


def signal_effective_n(x: np.ndarray, atoms: np.ndarray) -> tuple[float, int]:
    x = as_matrix(x)
    weights = []
    for atom in np.sort(np.unique(atoms)):
        m = atoms == atom
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


def frame_effective_n(y: np.ndarray, atoms: np.ndarray, frames: np.ndarray) -> tuple[float, float]:
    y = as_matrix(y)
    mag = y[:, 0] if y.shape[1] == 1 else np.linalg.norm(y, axis=1)
    total = 0.0
    rhos = []
    for atom in np.sort(np.unique(atoms)):
        m = atoms == atom
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


def verdict(gamma: float, se: float) -> str:
    if not np.isfinite(gamma):
        return "insufficient"
    if np.isfinite(se) and abs(gamma - 1.0) <= 2.0 * se:
        return "zero-circularity recovered law (gamma approx 1)"
    return "form-recovered-scale-fitted (gamma not approx 1)"


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


def manual_gamma_check(x: np.ndarray, y: np.ndarray, atoms: np.ndarray, mechanism: str, target: str) -> dict[str, float | str]:
    xw, yw = centre_by_atom(x, y, atoms)
    num = float((xw * yw).sum())
    den = float((xw * xw).sum())
    gamma = num / den
    via_func = gamma_fit(xw, yw)
    return {
        "mechanism": mechanism,
        "target": target,
        "numerator": num,
        "denominator": den,
        "manual_gamma": gamma,
        "function_gamma": via_func,
        "abs_diff": abs(gamma - via_func),
    }


def analyze_target(
    out_dir: Path,
    mechanism: str,
    spec: dict[str, str],
    split_fraction: float,
    seed: int,
    C: np.ndarray,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    csv_path = out_dir / spec["csv"]
    df = pd.read_csv(csv_path)
    require_columns(
        df,
        [
            "atom_index",
            "h5_row",
            "dft_present",
            "dft_local_frame_valid",
            "bare_kernel_present",
            "bare_T0",
            "dft_sigma_iso",
        ],
        mechanism,
    )
    atoms = df["atom_index"].to_numpy(int)
    frames = df["h5_row"].to_numpy(int)
    split = make_frame_split(frames, split_fraction, seed)
    rows: list[dict[str, object]] = []

    valid_scalar = (
        (df["dft_present"].to_numpy(int) == 1)
        & (df["bare_kernel_present"].to_numpy(int) == 1)
        & finite_rows(df["bare_T0"].to_numpy(float), df["dft_sigma_iso"].to_numpy(float))
    )
    x0 = df["bare_T0"].to_numpy(float).reshape(-1, 1)
    y0 = df["dft_sigma_iso"].to_numpy(float).reshape(-1, 1)
    rows.extend(
        summarize_protocol(
            mechanism,
            spec,
            "T0",
            x0,
            y0,
            atoms,
            frames,
            valid_scalar,
            split,
        )
    )

    kt2 = np.load(out_dir / spec["t2_kernel"]) @ C.T
    yt2 = np.load(out_dir / spec["t2_target"]) @ C.T
    if kt2.shape != (len(df), 5) or yt2.shape != (len(df), 5):
        raise RuntimeError(
            f"{mechanism} T2 shapes {kt2.shape}, {yt2.shape} do not match rows {len(df)}"
        )
    valid_t2 = (
        (df["dft_present"].to_numpy(int) == 1)
        & (df["dft_local_frame_valid"].to_numpy(int) == 1)
        & (df["bare_kernel_present"].to_numpy(int) == 1)
        & finite_rows(kt2, yt2)
    )
    rows.extend(
        summarize_protocol(
            mechanism,
            spec,
            "T2",
            kt2,
            yt2,
            atoms,
            frames,
            valid_t2,
            split,
        )
    )
    audit = manual_gamma_check(x0[valid_scalar], y0[valid_scalar], atoms[valid_scalar], mechanism, "T0")
    return rows, audit


def summarize_protocol(
    mechanism: str,
    spec: dict[str, str],
    target: str,
    x: np.ndarray,
    y: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    valid: np.ndarray,
    split: Split,
) -> list[dict[str, object]]:
    x = as_matrix(x)
    y = as_matrix(y)
    rows: list[dict[str, object]] = []
    atom_count = int(np.unique(atoms[valid]).size)
    signal_neff, signal_top90 = signal_effective_n(x[valid], atoms[valid])
    frame_neff, lag1 = frame_effective_n(y[valid], atoms[valid], frames[valid])

    x_frame, y_frame = centre_by_atom(x, y, atoms, mean_mask=valid & split.train)
    frame_score = valid & split.test & finite_rows(x_frame, y_frame)
    frame_train = valid & split.train & finite_rows(x_frame, y_frame)
    gamma_frame = gamma_fit(x_frame, y_frame, frame_train)
    se_frame = atom_jackknife_gamma_se(x_frame, y_frame, atoms, frame_train)
    row = base_result(
        mechanism,
        spec,
        target,
        "frame_split",
        x_frame,
        y_frame,
        frame_score,
        gamma_frame,
        se_frame,
        atom_count,
        signal_neff,
        signal_top90,
        frame_neff,
        lag1,
    )
    row["train_frames"] = split.train_frames
    row["test_frames"] = split.test_frames
    row["gamma_scaled_R2"] = r2_score(gamma_frame * x_frame[frame_score], y_frame[frame_score])
    rows.append(row)

    x_loao, y_loao = centre_by_atom(x, y, atoms, mean_mask=valid)
    loao_score = valid & finite_rows(x_loao, y_loao)
    gamma_all = gamma_fit(x_loao, y_loao, loao_score)
    se_all = atom_jackknife_gamma_se(x_loao, y_loao, atoms, loao_score)
    pred = loao_scaled_prediction(x_loao, y_loao, atoms, loao_score)
    row = base_result(
        mechanism,
        spec,
        target,
        "leave_atoms_out",
        x_loao,
        y_loao,
        loao_score,
        gamma_all,
        se_all,
        atom_count,
        signal_neff,
        signal_top90,
        frame_neff,
        lag1,
    )
    row["train_frames"] = int(np.unique(frames[valid]).size)
    row["test_frames"] = int(np.unique(frames[valid]).size)
    row["gamma_scaled_R2"] = r2_score(pred[loao_score], y_loao[loao_score])
    rows.append(row)
    return rows


def base_result(
    mechanism: str,
    spec: dict[str, str],
    target: str,
    protocol: str,
    x: np.ndarray,
    y: np.ndarray,
    score: np.ndarray,
    gamma: float,
    gamma_se: float,
    atom_count: int,
    signal_neff: float,
    signal_top90: int,
    frame_neff: float,
    lag1: float,
) -> dict[str, object]:
    if target == "T2":
        metrics = tensor_metrics(x, y, score)
        metrics["pearson_r"] = math.nan
        metrics["pearson_R2"] = math.nan
    else:
        metrics = scalar_metrics(x, y, score)
        metrics["component_r"] = math.nan
        metrics["absT2_r"] = math.nan
        metrics["row_cosine_mean"] = math.nan
    return {
        "mechanism": mechanism,
        "target": target,
        "stratum": spec["stratum"],
        "protocol": protocol,
        "rows_scored": int(score.sum()),
        "atoms_total": atom_count,
        "atom_signal_neff": signal_neff,
        "atom_signal_top90_count": signal_top90,
        "frame_effective_n_sum": frame_neff,
        "target_median_lag1": lag1,
        "gamma": gamma,
        "gamma_jackknife_se": gamma_se,
        "verdict_bucket": verdict(gamma, gamma_se),
        "unit_handling": (
            "emitted fixed kernel and DFT target are compared in matched ppm "
            "units after within-atom centering; gamma is dimensionless"
        ),
        "notes": spec["notes"],
        **metrics,
    }


def write_csv(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "mechanism",
        "target",
        "stratum",
        "protocol",
        "rows_scored",
        "atoms_total",
        "atom_signal_neff",
        "atom_signal_top90_count",
        "frame_effective_n_sum",
        "target_median_lag1",
        "pearson_r",
        "pearson_R2",
        "component_r",
        "absT2_r",
        "row_cosine_mean",
        "gamma",
        "gamma_jackknife_se",
        "gamma_scaled_R2",
        "train_frames",
        "test_frames",
        "verdict_bucket",
        "unit_handling",
        "notes",
    ]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def write_report(rows: list[dict[str, object]], audits: list[dict[str, object]], path: Path, csv_path: Path, out_dir: Path) -> None:
    lines = [
        "# Literature-Coefficient-Fixed De-Circularisation",
        "",
        "Read-only Python analysis of the already-emitted ring-current and McConnell substrate.",
        f"Substrate: `{out_dir}`",
        f"CSV artifact: `{csv_path}`",
        "",
        "Unit handling: emitted fixed kernels and DFT targets are compared in matched ppm units after within-atom centering. The reported `gamma` is therefore a dimensionless post-hoc multiplier. Correlation columns are un-fitted; `gamma_scaled_R2` is only the scale diagnostic.",
        "",
        "Thin-stratum handling: `atoms_total` is the emitted home-stratum atom count. `atom_signal_neff` is a participation-ratio effective atom count from emitted fixed-kernel modulation, and `atom_signal_top90_count` is the number of atoms carrying 90% of that modulation.",
        "",
        "Parameterized for the 750-frame substrate: rerun with `--out-dir` pointing at the later emitted ring/mc directory; no row count is hard-coded.",
        "",
        "## Results",
        "",
    ]
    table_cols = [
        "mechanism",
        "target",
        "protocol",
        "atoms_total",
        "atom_signal_neff",
        "atom_signal_top90_count",
        "frame_effective_n_sum",
        "target_median_lag1",
        "pearson_r",
        "pearson_R2",
        "component_r",
        "absT2_r",
        "gamma",
        "gamma_jackknife_se",
        "gamma_scaled_R2",
        "verdict_bucket",
    ]
    lines.append("| " + " | ".join(table_cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(table_cols)) + " |")
    for row in rows:
        vals = [fmt(row.get(c)) for c in table_cols]
        lines.append("| " + " | ".join(vals) + " |")
    lines.extend(
        [
            "",
            "## Self Audit",
            "",
            "- No fitting is used for `pearson_r`, `pearson_R2`, `component_r`, or `absT2_r`; those correlate the fixed emitted kernel directly against DFT on the scoring rows.",
            "- `gamma` is computed separately as `(K dot DFT) / (K dot K)` on centered emitted arrays. The leave-atoms-out `gamma_scaled_R2` refits gamma on all other atoms for each held-out atom.",
            "- T2 sidecars and targets are mapped with the frozen `change_of_basis.get_C()` matrix; no basis is re-derived in this script.",
        ]
    )
    for audit in audits:
        lines.append(
            "- Manual gamma check {mechanism} {target}: numerator={num}, denominator={den}, manual={manual}, function={func}, abs_diff={diff}.".format(
                mechanism=audit["mechanism"],
                target=audit["target"],
                num=fmt(audit["numerator"], 8),
                den=fmt(audit["denominator"], 8),
                manual=fmt(audit["manual_gamma"], 8),
                func=fmt(audit["function_gamma"], 8),
                diff=fmt(audit["abs_diff"], 8),
            )
        )
    lines.extend(
        [
            "",
            "## Verdict Buckets",
            "",
            "The bucket rule is mechanical and deliberately reserved: if `gamma = 1` lies within about two delete-atom jackknife SEs, the row is tagged `zero-circularity recovered law`; otherwise it is tagged `form-recovered-scale-fitted`. The scientific call remains for the lead.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    artifact_dir = Path(args.artifact_dir)
    report_path = Path(args.report_md)
    C = cob.get_C()
    rows: list[dict[str, object]] = []
    audits: list[dict[str, object]] = []
    for mechanism, spec in TARGETS.items():
        more_rows, audit = analyze_target(
            out_dir,
            mechanism,
            spec,
            args.test_fraction,
            args.seed,
            C,
        )
        rows.extend(more_rows)
        audits.append(audit)

    csv_path = artifact_dir / "literature_fixed_decirc.csv"
    audit_path = artifact_dir / "literature_fixed_decirc_audit.json"
    write_csv(rows, csv_path)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(audits, indent=2, sort_keys=True), encoding="utf-8")
    write_report(rows, audits, report_path, csv_path, out_dir)

    print("literature fixed de-circularisation")
    print(f"substrate={out_dir}")
    print(f"csv={csv_path}")
    print(f"report={report_path}")
    for row in rows:
        if row["target"] == "T2":
            metric = f"component_r={fmt(row['component_r'])} absT2_r={fmt(row['absT2_r'])}"
        else:
            metric = f"pearson_r={fmt(row['pearson_r'])} R2={fmt(row['pearson_R2'])}"
        print(
            f"{row['mechanism']} {row['target']} {row['protocol']}: "
            f"{metric} gamma={fmt(row['gamma'])} +/- {fmt(row['gamma_jackknife_se'])} "
            f"Neff_atoms={fmt(row['atom_signal_neff'])} top90={row['atom_signal_top90_count']} "
            f"{row['verdict_bucket']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
