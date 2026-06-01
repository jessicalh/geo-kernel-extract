#!/usr/bin/env python3
"""Scalar Buckingham APBS E-field -> DFT sigma_iso fitter.

Inputs are only the producer-emitted buckingham_efield substrate:
  buckingham_efield_aggregated.csv
  buckingham_efield_feature_field_local.npy

The C++ spine owns the local backbone frame and emits E_proj = E_local.z and
|E|. This fitter reads those scalars and fits, per backbone stratum:

  sigma_iso ~= atom_baseline + A * E_proj + B * |E|^2

The reported gate is frame-split within-atom sigma_iso R2. Target T1/T2 payloads
are deliberately ignored; T1 is convention-unverified.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


FV_HN = {1, 2}
FV_N = {4, 5}
FV_CA = {6}
FV_C = {7}
FV_O = {8}
FV_HA = {9}

STRATA_ORDER = ["N", "CA", "C", "O", "HN", "HA"]
TEST_FRAME_FRACTION = 0.20
FRAME_SPLIT_SEED = 0
MIN_FRAME_SPLIT_FRAMES = 5
THIN_ATOM_WARN = 10


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "out_dir",
        nargs="?",
        default=os.environ.get("REDISCOVER_OUT", "/tmp/rdc-buckingham"),
        help="directory containing buckingham_efield_aggregated.csv",
    )
    ap.add_argument("--alpha", type=float, default=0.0, help="ridge alpha; 0 = OLS")
    ap.add_argument("--seed", type=int, default=FRAME_SPLIT_SEED)
    return ap.parse_args()


def require(path, what):
    if not os.path.exists(path):
        sys.exit(
            f"FATAL: required emitted {what} not found:\n  {path}\n"
            "Run the extractor with --case buckingham_efield into this output directory."
        )


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


def r2(pred, target):
    pred = np.asarray(pred, dtype=float)
    target = np.asarray(target, dtype=float)
    ok = np.isfinite(pred) & np.isfinite(target)
    if ok.sum() < 2:
        return float("nan")
    target = target[ok]
    pred = pred[ok]
    denom = ((target - target.mean()) ** 2).sum()
    if denom <= 0:
        return float("nan")
    return float(1.0 - ((target - pred) ** 2).sum() / denom)


def make_frame_split(frames, seed):
    frames = np.sort(np.unique(frames))
    if len(frames) < MIN_FRAME_SPLIT_FRAMES:
        return set()
    n_test = max(1, int(TEST_FRAME_FRACTION * len(frames)))
    return set(np.random.default_rng(seed).choice(frames, n_test, replace=False))


def load(out_dir):
    agg_csv = os.path.join(out_dir, "buckingham_efield_aggregated.csv")
    field_npy = os.path.join(out_dir, "buckingham_efield_feature_field_local.npy")
    require(agg_csv, "aggregated CSV")
    require(field_npy, "local APBS E-field NPY")
    agg = pd.read_csv(agg_csv)
    field_local = np.load(field_npy)
    if field_local.ndim != 2 or field_local.shape[1] != 3:
        sys.exit(
            f"FATAL: buckingham_efield_feature_field_local.npy has shape "
            f"{field_local.shape}, expected (rows, 3)."
        )
    if len(agg) != field_local.shape[0]:
        sys.exit(
            f"FATAL: CSV rows {len(agg)} != local-field NPY rows "
            f"{field_local.shape[0]} (sidecars out of sync)."
        )
    return agg, field_local


def centred_by_train_atom(values, atoms, train_mask):
    values = np.asarray(values, dtype=float)
    out = np.full_like(values, np.nan, dtype=float)
    means = {}
    for atom in np.unique(atoms):
        train_atom = (atoms == atom) & train_mask
        if train_atom.sum() == 0:
            continue
        mean = values[train_atom].mean(axis=0)
        means[atom] = mean
        out[atoms == atom] = values[atoms == atom] - mean
    return out, means


def fit_beta(x_train, y_train, alpha):
    xtx = x_train.T @ x_train
    if alpha > 0:
        xtx = xtx + alpha * np.eye(xtx.shape[0])
    rhs = x_train.T @ y_train
    try:
        return np.linalg.solve(xtx, rhs)
    except np.linalg.LinAlgError:
        return np.linalg.lstsq(x_train, y_train, rcond=None)[0]


def fit_stratum(df, seed, alpha):
    present = (
        (df["dft_present"].to_numpy() == 1)
        & (df["apbs_efield_present"].to_numpy() == 1)
        & (df["frame_valid"].to_numpy() == 1)
        & np.isfinite(df["E_proj"].to_numpy(float))
        & np.isfinite(df["E_mag"].to_numpy(float))
        & np.isfinite(df["dft_sigma_iso"].to_numpy(float))
    )
    d = df.loc[present].copy().reset_index(drop=True)
    if len(d) < 4:
        return None

    frames = d["h5_row"].to_numpy()
    test_frames = make_frame_split(frames, seed)
    if not test_frames:
        return None
    test = d["h5_row"].isin(test_frames).to_numpy()
    train = ~test

    x = np.column_stack(
        [
            d["E_proj"].to_numpy(float),
            d["E_mag"].to_numpy(float) ** 2,
        ]
    )
    y = d["dft_sigma_iso"].to_numpy(float)
    atoms = d["atom_index"].to_numpy()

    xw, _ = centred_by_train_atom(x, atoms, train)
    yw, _ = centred_by_train_atom(y.reshape(-1, 1), atoms, train)
    yw = yw.reshape(-1)

    fit_mask = train & np.isfinite(xw).all(axis=1) & np.isfinite(yw)
    score_mask = test & np.isfinite(xw).all(axis=1) & np.isfinite(yw)
    if fit_mask.sum() < 3 or score_mask.sum() < 2:
        return None

    beta = fit_beta(xw[fit_mask], yw[fit_mask], alpha)
    pred = xw[score_mask] @ beta
    score = r2(pred, yw[score_mask])
    train_score = r2(xw[fit_mask] @ beta, yw[fit_mask])

    return {
        "r2": score,
        "train_r2": train_score,
        "A": float(beta[0]),
        "B": float(beta[1]),
        "effective_N": int(np.unique(atoms).size),
        "rows": int(len(d)),
        "train_rows": int(fit_mask.sum()),
        "test_rows": int(score_mask.sum()),
        "frames": int(np.unique(frames).size),
        "test_frames": int(len(test_frames)),
        "median_E_mag": float(np.median(d["E_mag"].to_numpy(float))),
    }


def fmt(x):
    return "nan" if not np.isfinite(x) else f"{x:+.4f}"


def coef(x):
    return "nan" if not np.isfinite(x) else f"{x:+.6g}"


def main():
    args = parse_args()
    agg, field_local = load(args.out_dir)
    del field_local  # contract/shape check only; scalar fit reads CSV scalars.

    agg["stratum"] = [stratum_of(v) for v in agg["frame_variant"]]
    print(
        "Buckingham APBS E-field scalar fit: sigma_iso modulation ~= "
        "A*E_proj + B*|E|^2"
    )
    print(
        "T0-only. T1 sidecar is unverified_emit_only and is not read by this fitter."
    )
    print(
        "PB-vs-Coulomb note: E_proj is APBS solvated-PB; broad_backbone "
        "charge_field_z is the FF14SB vacuum-Coulomb contrast."
    )
    print(f"rows={len(agg)}  strata_rows={agg['stratum'].notna().sum()}")

    results = []
    for name in STRATA_ORDER:
        d = agg[agg["stratum"] == name]
        res = fit_stratum(d, args.seed, args.alpha)
        if res is None:
            print(f"{name:2s} insufficient data for honest frame-split fit")
            continue
        thin = " THIN" if res["effective_N"] < THIN_ATOM_WARN else ""
        print(
            f"{name:2s} sigma_iso_R2={fmt(res['r2'])} "
            f"train_R2={fmt(res['train_r2'])} "
            f"A={coef(res['A'])} B={coef(res['B'])} "
            f"effective_N={res['effective_N']} rows={res['rows']} "
            f"test_rows={res['test_rows']} frames={res['frames']} "
            f"test_frames={res['test_frames']} median_|E|={res['median_E_mag']:.6g}"
            f"{thin}"
        )
        results.append((name, res))

    print("\nstratum,sigma_iso_R2_frame_split,A,B,effective_N,rows,test_rows,flag")
    for name, res in results:
        flag = "THIN" if res["effective_N"] < THIN_ATOM_WARN else ""
        print(
            f"{name},{fmt(res['r2'])},{coef(res['A'])},{coef(res['B'])},"
            f"{res['effective_N']},{res['rows']},{res['test_rows']},{flag}"
        )


if __name__ == "__main__":
    main()
