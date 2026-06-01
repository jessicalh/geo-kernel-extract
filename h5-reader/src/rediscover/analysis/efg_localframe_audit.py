#!/usr/bin/env python3
"""Audit APBS EFG lab-frame versus local-frame T2 sidecars.

Reads only emitted EFG CSV/NPY payloads. No H5, no field or tensor
reconstruction, no projection. Tensor arrays are mapped through the frozen
library->e3nn 2e basis from change_of_basis.get_C().
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


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("efg_dir", nargs="?", default="/tmp/rdc-efg")
    ap.add_argument("--out-dir", default="/tmp/rdc-efg-localframe-audit")
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


def corrcoef(a, b):
    a = np.asarray(a, dtype=float).reshape(-1)
    b = np.asarray(b, dtype=float).reshape(-1)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return math.nan
    x = a[ok] - a[ok].mean()
    y = b[ok] - b[ok].mean()
    den = math.sqrt(float(np.dot(x, x) * np.dot(y, y)))
    if den <= 0.0:
        return math.nan
    return float(np.dot(x, y) / den)


def r2_score(pred, target):
    pred = np.asarray(pred, dtype=float)
    target = np.asarray(target, dtype=float)
    ok = np.isfinite(pred).all(axis=1) & np.isfinite(target).all(axis=1)
    if ok.sum() < 3:
        return math.nan
    p = pred[ok]
    y = target[ok]
    den = float(((y - y.mean(axis=0, keepdims=True)) ** 2).sum())
    if den <= 0.0:
        return math.nan
    return float(1.0 - ((y - p) ** 2).sum() / den)


def demean_by_atom(values, atoms):
    values = np.asarray(values, dtype=float)
    out = np.full_like(values, np.nan, dtype=float)
    for atom in np.unique(atoms):
        m = atoms == atom
        ok = m & np.isfinite(values).all(axis=1)
        if ok.sum() == 0:
            continue
        out[m] = values[m] - values[ok].mean(axis=0, keepdims=True)
    return out


def median_lag1_component_rho(values, atoms, frames):
    values = np.asarray(values, dtype=float)
    rhos = []
    for atom in np.unique(atoms):
        idx = np.flatnonzero(atoms == atom)
        idx = idx[np.argsort(frames[idx])]
        for c in range(values.shape[1]):
            y = values[idx, c]
            y = y[np.isfinite(y)]
            if len(y) < 3:
                continue
            yd = y - y.mean()
            den = float(np.dot(yd, yd))
            if den <= 0.0:
                rhos.append(0.0)
                continue
            rhos.append(max(min(float(np.dot(yd[:-1], yd[1:]) / den), 0.999), -0.999))
    return float(np.median(rhos)) if rhos else math.nan


def tensor_metrics(feature, target, atoms, frames):
    feature_dm = demean_by_atom(feature, atoms)
    target_dm = demean_by_atom(target, atoms)
    ok = np.isfinite(feature_dm).all(axis=1) & np.isfinite(target_dm).all(axis=1)
    if ok.sum() < 3:
        return {
            "feature_target_component_r": math.nan,
            "feature_target_absT2_r": math.nan,
            "constant_gamma_R2": math.nan,
            "target_median_lag1_rho": math.nan,
            "feature_median_lag1_rho": math.nan,
        }
    den = float((feature_dm[ok] * feature_dm[ok]).sum())
    gamma = 0.0 if den <= 0.0 else float((feature_dm[ok] * target_dm[ok]).sum() / den)
    return {
        "feature_target_component_r": corrcoef(feature_dm[ok].reshape(-1), target_dm[ok].reshape(-1)),
        "feature_target_absT2_r": corrcoef(
            np.linalg.norm(feature_dm[ok], axis=1),
            np.linalg.norm(target_dm[ok], axis=1),
        ),
        "constant_gamma_R2": r2_score(gamma * feature_dm[ok], target_dm[ok]),
        "target_median_lag1_rho": median_lag1_component_rho(target_dm[ok], atoms[ok], frames[ok]),
        "feature_median_lag1_rho": median_lag1_component_rho(feature_dm[ok], atoms[ok], frames[ok]),
    }


def load_pair(efg_dir, frame, C):
    if frame == "local":
        feature = np.load(efg_dir / "efg_feature_T2.npy")
        target = np.load(efg_dir / "efg_target_T2.npy")
    elif frame == "lab":
        feature = np.load(efg_dir / "efg_feature_lab_T2.npy")
        target = np.load(efg_dir / "efg_target_lab_T2.npy")
    else:
        raise ValueError(frame)
    return feature @ C.T, target @ C.T


def main():
    args = parse_args()
    efg_dir = Path(args.efg_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    C = cob.get_C()
    orth = float(np.abs(C.T @ C - np.eye(5)).max())
    df = pd.read_csv(efg_dir / "efg_aggregated.csv")
    df["stratum"] = [stratum_of(v) for v in df["frame_variant"]]

    rows = []
    for frame in ["lab", "local"]:
        feature, target = load_pair(efg_dir, frame, C)
        for stratum in STRATA_ORDER:
            base = (
                (df["stratum"].to_numpy() == stratum)
                & (df["dft_present"].to_numpy(int) == 1)
                & (df["apbs_efg_present"].to_numpy(int) == 1)
            )
            if frame == "local":
                base &= df["frame_valid"].to_numpy(int) == 1
            base &= np.isfinite(feature).all(axis=1) & np.isfinite(target).all(axis=1)
            idx = np.flatnonzero(base)
            if len(idx) == 0:
                continue
            atoms = df["atom_index"].to_numpy(int)[idx]
            frames = df["h5_row"].to_numpy(int)[idx]
            metrics = tensor_metrics(feature[idx], target[idx], atoms, frames)
            rows.append(
                {
                    "frame": frame,
                    "stratum": stratum,
                    "rows": int(len(idx)),
                    "atoms": int(len(np.unique(atoms))),
                    "frame_valid_fraction": float(
                        df["frame_valid"].to_numpy(int)[idx].mean()
                        if "frame_valid" in df
                        else math.nan
                    ),
                    **metrics,
                }
            )

    result = pd.DataFrame(rows).sort_values(["stratum", "frame"]).reset_index(drop=True)
    csv_path = out_dir / "efg_localframe_audit.csv"
    md_path = out_dir / "efg_localframe_audit.md"
    json_path = out_dir / "efg_localframe_audit_run.json"
    result.to_csv(csv_path, index=False)
    json_path.write_text(
        json.dumps(
            {
                "efg_dir": str(efg_dir),
                "change_of_basis_orthogonality_max": orth,
                "csv": str(csv_path),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    cols = [
        "frame", "stratum", "atoms", "rows", "frame_valid_fraction",
        "target_median_lag1_rho", "feature_median_lag1_rho",
        "feature_target_component_r", "feature_target_absT2_r",
        "constant_gamma_R2",
    ]
    lines = ["# EFG Local-Frame Audit", ""]
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
    for _, row in result[cols].iterrows():
        vals = []
        for col in cols:
            value = row[col]
            if isinstance(value, float):
                vals.append("nan" if not np.isfinite(value) else f"{value:.4g}")
            else:
                vals.append(str(value))
        lines.append("| " + " | ".join(vals) + " |")
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"change_of_basis.get_C() |C^T C - I|max={orth:.3e}")
    print(f"wrote {csv_path}")
    print(f"wrote {md_path}")
    print("frame,stratum,target_lag1,feature_lag1,component_r,absT2_r,gamma_R2")
    for _, row in result.iterrows():
        print(
            f"{row['frame']},{row['stratum']},"
            f"{row['target_median_lag1_rho']:.4g},"
            f"{row['feature_median_lag1_rho']:.4g},"
            f"{row['feature_target_component_r']:.4g},"
            f"{row['feature_target_absT2_r']:.4g},"
            f"{row['constant_gamma_R2']:.4g}"
        )


if __name__ == "__main__":
    main()
