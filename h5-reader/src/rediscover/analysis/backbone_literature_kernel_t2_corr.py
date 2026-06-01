#!/usr/bin/env python3
"""Read-only broad-backbone literature-kernel T2 correlation.

Front-load analysis/PATTERNS.md: Python reads the C++ emitted CSV/NPY substrate
only. It does not open H5, rebuild kernels, project tensors, or infer physics.

Inputs, row-aligned by broad_backbone_aggregated.csv:
  broad_backbone_aggregated_literature_kernel_T2.npy
  broad_backbone_aggregated_ring_literature_kernel_T2.npy
  broad_backbone_aggregated_bond_literature_kernel_T2.npy
  broad_backbone_aggregated_charge_literature_kernel_T2.npy
  broad_backbone_aggregated_target_local_T2.npy
"""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict

import numpy as np


FV_HN = {1, 2}
FV_N = {4, 5}
FV_CA = {6}
FV_C = {7}
FV_O = {8}
FV_HA = {9}

STRATA_ORDER = ["N", "CA", "C", "O", "HN", "HA", "HA2", "HA3"]

KERNELS = {
    "total": (
        "literature_kernel_present",
        "broad_backbone_aggregated_literature_kernel_T2.npy",
    ),
    "ring": (
        "ring_literature_kernel_present",
        "broad_backbone_aggregated_ring_literature_kernel_T2.npy",
    ),
    "bond": (
        "bond_literature_kernel_present",
        "broad_backbone_aggregated_bond_literature_kernel_T2.npy",
    ),
    "charge": (
        "charge_literature_kernel_present",
        "broad_backbone_aggregated_charge_literature_kernel_T2.npy",
    ),
}


def parse_args():
    ap = argparse.ArgumentParser(
        description="Unfitted fixed-literature broad-backbone T2 correlation."
    )
    ap.add_argument(
        "out_dir",
        nargs="?",
        default="/tmp/rdc-broad-backbone-axes",
        help="directory containing broad_backbone_aggregated.csv and sidecars",
    )
    ap.add_argument(
        "--kernel",
        choices=["all", *KERNELS.keys()],
        default="all",
        help="which emitted kernel sidecar to report",
    )
    ap.add_argument(
        "--report-csv",
        default=None,
        help="optional output CSV; defaults inside out_dir",
    )
    return ap.parse_args()


def stratum_of(frame_variant, atom_name):
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


def read_rows(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"no aggregated rows in {path}")
    return rows


def require_columns(rows, names):
    have = set(rows[0].keys())
    missing = [name for name in names if name not in have]
    if missing:
        raise RuntimeError("missing emitted columns: " + ", ".join(missing))


def as_int(row, name):
    return int(float(row[name]))


def finite_rows(kernel, target):
    return np.isfinite(kernel).all(axis=1) & np.isfinite(target).all(axis=1)


def pearson(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 3:
        return math.nan
    x = x - x.mean()
    y = y - y.mean()
    denom = math.sqrt(float(np.dot(x, x) * np.dot(y, y)))
    if denom == 0.0:
        return math.nan
    return float(np.dot(x, y) / denom)


def mean_row_cosine(kernel, target):
    dot = np.sum(kernel * target, axis=1)
    denom = np.linalg.norm(kernel, axis=1) * np.linalg.norm(target, axis=1)
    mask = np.isfinite(dot) & np.isfinite(denom) & (denom > 0.0)
    if not np.any(mask):
        return math.nan
    return float(np.mean(dot[mask] / denom[mask]))


def effective_n(rows, mask, target):
    groups = defaultdict(list)
    target_mag = np.linalg.norm(target, axis=1)
    for i, row in enumerate(rows):
        if not mask[i] or not np.isfinite(target_mag[i]):
            continue
        atom = as_int(row, "atom_index")
        h5_row = as_int(row, "h5_row")
        groups[atom].append((h5_row, float(target_mag[i])))

    total = 0.0
    rhos = []
    for values in groups.values():
        values = [v for _, v in sorted(values)]
        n = len(values)
        if n < 3:
            total += float(n)
            continue
        y = np.asarray(values, dtype=np.float64)
        yd = y - y.mean()
        denom = float(np.dot(yd, yd))
        if denom <= 0.0:
            rho = 0.0
            neff = float(n)
        else:
            rho = float(np.dot(yd[:-1], yd[1:]) / denom)
            rho = max(min(rho, 0.999), -0.999)
            neff = n * (1.0 - rho) / (1.0 + rho)
            neff = min(max(neff, 1.0), float(n))
        rhos.append(rho)
        total += neff
    median_rho = float(np.median(rhos)) if rhos else math.nan
    return total, len(groups), median_rho


def summarize(rows, kernel_name, present_col, kernel, target):
    require_columns(rows, [present_col, "dft_present", "dft_local_frame_valid"])
    strata = np.array(
        [stratum_of(row["frame_variant"], row["atom_name"]) for row in rows],
        dtype=object,
    )
    present = np.array([as_int(row, present_col) == 1 for row in rows])
    dft_ok = np.array(
        [
            as_int(row, "dft_present") == 1
            and as_int(row, "dft_local_frame_valid") == 1
            for row in rows
        ]
    )
    finite = finite_rows(kernel, target)
    base_mask = present & dft_ok & finite & np.array([s is not None for s in strata])

    out = []
    for stratum in ["ALL", *STRATA_ORDER]:
        if stratum == "ALL":
            mask = base_mask
        else:
            mask = base_mask & (strata == stratum)
        n_rows = int(mask.sum())
        if n_rows:
            component_r = pearson(kernel[mask].reshape(-1), target[mask].reshape(-1))
            magnitude_r = pearson(
                np.linalg.norm(kernel[mask], axis=1),
                np.linalg.norm(target[mask], axis=1),
            )
            cosine = mean_row_cosine(kernel[mask], target[mask])
        else:
            component_r = math.nan
            magnitude_r = math.nan
            cosine = math.nan
        neff, n_atoms, median_rho = effective_n(rows, mask, target)
        out.append(
            {
                "kernel": kernel_name,
                "stratum": stratum,
                "rows": n_rows,
                "atoms": n_atoms,
                "effective_n": neff,
                "median_lag1_rho": median_rho,
                "component_r": component_r,
                "magnitude_r": magnitude_r,
                "mean_row_cosine": cosine,
            }
        )
    return out


def fmt(x):
    if isinstance(x, int):
        return str(x)
    if not np.isfinite(x):
        return "nan"
    return f"{x:.4g}"


def main():
    args = parse_args()
    agg_path = os.path.join(args.out_dir, "broad_backbone_aggregated.csv")
    target_path = os.path.join(
        args.out_dir, "broad_backbone_aggregated_target_local_T2.npy"
    )
    rows = read_rows(agg_path)
    target = np.load(target_path)
    if target.shape != (len(rows), 5):
        raise RuntimeError(
            f"target shape {target.shape} does not match rows {(len(rows), 5)}"
        )

    names = list(KERNELS) if args.kernel == "all" else [args.kernel]
    all_results = []
    for name in names:
        present_col, file_name = KERNELS[name]
        kernel = np.load(os.path.join(args.out_dir, file_name))
        if kernel.shape != (len(rows), 5):
            raise RuntimeError(
                f"{file_name} shape {kernel.shape} does not match rows {(len(rows), 5)}"
            )
        all_results.extend(summarize(rows, name, present_col, kernel, target))

    report_csv = args.report_csv
    if report_csv is None:
        report_csv = os.path.join(
            args.out_dir, "broad_backbone_literature_kernel_t2_correlation.csv"
        )
    fields = [
        "kernel",
        "stratum",
        "rows",
        "atoms",
        "effective_n",
        "median_lag1_rho",
        "component_r",
        "magnitude_r",
        "mean_row_cosine",
    ]
    with open(report_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in all_results:
            w.writerow(row)

    print("unfitted broad-backbone literature-kernel T2 correlation")
    print(f"substrate={args.out_dir}")
    print(f"report_csv={report_csv}")
    print(
        "kernel,stratum,rows,atoms,effective_n,component_r,magnitude_r,mean_row_cosine"
    )
    for row in all_results:
        print(
            ",".join(
                [
                    row["kernel"],
                    row["stratum"],
                    str(row["rows"]),
                    str(row["atoms"]),
                    fmt(row["effective_n"]),
                    fmt(row["component_r"]),
                    fmt(row["magnitude_r"]),
                    fmt(row["mean_row_cosine"]),
                ]
            )
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
