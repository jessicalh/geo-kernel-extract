#!/usr/bin/env python3
"""EFG scalar-gate Depth-A evidence collector.

Reads only the emitted /tmp/rdc-efg sidecars:
  efg_aggregated.csv, efg_feature_T2.npy, efg_target_T2.npy

The model path is the existing Schur-scalar APBS-EFG law:
  pred_T2 = g(|EFG|) * EFG_T2

This script distills whether the trained gate is approximately constant per
backbone stratum, compares the fitted constant-g Buckingham form against the
learned nonlinear gate, and reports effective N. It reuses the frozen
change_of_basis.get_C() path through equiv_t2_efg_e3nn.build_pack().
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import torch

import change_of_basis as cob
from equiv_t2_efg_e3nn import (
    STRATA_ORDER,
    build_pack,
    centred_by_train_atom,
    corrcoef,
    load,
    r2,
    stratum_of,
    train_model,
)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", nargs="?", default="/tmp/rdc-efg")
    ap.add_argument("--evidence-dir", default="/tmp/rdc-efg-arc-evidence")
    ap.add_argument("--epochs", type=int, default=4000)
    ap.add_argument("--lr", type=float, default=3e-3)
    ap.add_argument("--hidden", type=int, default=32)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--grid-size", type=int, default=101)
    ap.add_argument("--strata", nargs="*", default=STRATA_ORDER)
    ap.add_argument("--split", choices=["blocked", "random"], default="blocked")
    ap.add_argument("--purge-frames", type=int, default=1)
    return ap.parse_args()


def finite_fmt(x):
    return "nan" if not np.isfinite(x) else f"{x:.6g}"


def prediction_scores(pred, target, mask):
    if int(mask.sum()) < 2:
        return {
            "r2": math.nan,
            "component_r": math.nan,
            "magnitude_r": math.nan,
        }
    return {
        "r2": r2(pred[mask], target[mask]),
        "component_r": corrcoef(pred[mask].reshape(-1), target[mask].reshape(-1)),
        "magnitude_r": corrcoef(
            np.linalg.norm(pred[mask], axis=1),
            np.linalg.norm(target[mask], axis=1),
        ),
    }


def constant_gate_prediction(pack):
    feature = pack["feature"].detach().cpu().numpy()
    target = pack["target"].detach().cpu().numpy()
    group_atom = pack["group_atom_idx"]
    train = pack["g_tr"].detach().cpu().numpy().astype(bool)
    feature_dm = centred_by_train_atom(feature, group_atom, train)
    denom = float((feature_dm[train] ** 2).sum())
    if denom <= 0.0:
        gate = math.nan
        pred = np.full_like(target, np.nan)
    else:
        gate = float((target[train] * feature_dm[train]).sum() / denom)
        pred = gate * feature_dm
    return gate, pred


def model_prediction(model, pack):
    with torch.no_grad():
        return model(pack["feature"], pack["mag_norm"], pack["group_atom"],
                     center_mask=pack["g_tr"]).detach().cpu().numpy()


def gate_samples(model, pack, stratum, grid_size):
    mag = np.asarray(pack["mag"], dtype=float)
    lo = float(np.min(mag))
    hi = float(np.max(mag))
    q01, q05, q50, q95, q99 = np.quantile(mag, [0.01, 0.05, 0.50, 0.95, 0.99])
    grid = np.linspace(lo, hi, grid_size)
    norm = ((grid - pack["mag_mean"]) / pack["mag_std"]).reshape(-1, 1)
    dev = pack["feature"].device
    with torch.no_grad():
        gate = model.gate(torch.tensor(norm, dtype=torch.float32, device=dev))
    values = gate.detach().cpu().numpy().reshape(-1)
    rows = pd.DataFrame(
        {
            "stratum": stratum,
            "efg_T2_magnitude": grid,
            "gate": values,
            "inside_1_99pct": (grid >= q01) & (grid <= q99),
            "inside_5_95pct": (grid >= q05) & (grid <= q95),
        }
    )

    central = rows.loc[rows["inside_5_95pct"], "gate"].to_numpy(float)
    all_values = rows["gate"].to_numpy(float)
    mean = float(np.mean(central)) if len(central) else math.nan
    median = float(np.median(central)) if len(central) else math.nan
    std = float(np.std(central)) if len(central) else math.nan
    p05 = float(np.quantile(central, 0.05)) if len(central) else math.nan
    p95 = float(np.quantile(central, 0.95)) if len(central) else math.nan
    span = p95 - p05 if np.isfinite(p05) and np.isfinite(p95) else math.nan
    rel_span = abs(span) / max(abs(median), 1e-12) if np.isfinite(span) else math.nan
    return rows, {
        "efg_min": lo,
        "efg_p01": float(q01),
        "efg_p05": float(q05),
        "efg_median": float(q50),
        "efg_p95": float(q95),
        "efg_p99": float(q99),
        "efg_max": hi,
        "gate_mean_5_95": mean,
        "gate_median_5_95": median,
        "gate_std_5_95": std,
        "gate_p05_5_95": p05,
        "gate_p95_5_95": p95,
        "gate_span_5_95": span,
        "gate_rel_span_5_95": rel_span,
        "gate_min_full": float(np.min(all_values)),
        "gate_max_full": float(np.max(all_values)),
        "gate_vs_efg_r_full": corrcoef(grid, all_values),
    }


def effective_n_lag1(pack):
    target = pack["target"].detach().cpu().numpy()
    target_mag = np.linalg.norm(target, axis=1)
    h5_row = np.asarray(pack["h5_row"])
    group_atom = np.asarray(pack["group_atom_idx"])
    total = 0.0
    rhos = []
    for atom in np.unique(group_atom):
        idx = np.flatnonzero(group_atom == atom)
        idx = idx[np.argsort(h5_row[idx])]
        y = target_mag[idx]
        y = y[np.isfinite(y)]
        n = len(y)
        if n < 3:
            total += float(n)
            continue
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
        total += neff
        rhos.append(rho)
    return {
        "effective_n_lag1": float(total),
        "median_lag1_rho": float(np.median(rhos)) if rhos else math.nan,
        "atoms": int(len(np.unique(group_atom))),
    }


def gate_form_label(row):
    delta = row["nonlinear_test_r2"] - row["constant_test_r2"]
    rel_span = row["gate_rel_span_5_95"]
    if row["atoms"] < 10:
        thin = " THIN"
    else:
        thin = ""
    if not np.isfinite(row["nonlinear_test_r2"]) or not np.isfinite(row["constant_test_r2"]):
        return "insufficient" + thin
    if delta <= 0.02 and np.isfinite(rel_span) and rel_span <= 0.20:
        return "approximately constant" + thin
    if delta <= 0.02:
        return "linear-form predictive; learned gate drifts" + thin
    return "|EFG|-dependent" + thin


def main():
    args = parse_args()
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    evidence_dir = Path(args.evidence_dir)
    evidence_dir.mkdir(parents=True, exist_ok=True)

    C = cob.get_C()
    agg, feature_lib, target_lib = load(args.out_dir)
    agg = agg.copy()
    agg["__stratum"] = [stratum_of(fv, nm) for fv, nm in zip(agg["frame_variant"], agg["atom_name"])]

    print(
        f"device={dev} out_dir={args.out_dir} evidence_dir={evidence_dir} "
        f"epochs={args.epochs} hidden={args.hidden} split={args.split} "
        f"purge_frames={args.purge_frames}"
    )
    print(f"change_of_basis.get_C() |C^T C - I|max={np.abs(C.T @ C - np.eye(5)).max():.3e}")
    print(f"emitted rows={len(agg)} feature_shape={feature_lib.shape} target_shape={target_lib.shape}")

    summary_rows = []
    sample_frames = []
    for stratum in args.strata:
        idx = np.flatnonzero(agg["__stratum"].to_numpy() == stratum)
        if len(idx) == 0:
            continue
        pack = build_pack(
            agg.iloc[idx].reset_index(drop=True),
            feature_lib[idx],
            target_lib[idx],
            C,
            dev,
            split_strategy=args.split,
            split_seed=args.seed,
            purge_frames=args.purge_frames,
        )
        if pack is None or pack["n_atoms"] < 2:
            continue

        print(
            f"-- {stratum}: rows={pack['n_groups']} atoms={pack['n_atoms']} "
            f"frames={pack['n_frames']} train={int(pack['g_tr'].sum().item())} "
            f"test={int(pack['g_te'].sum().item())}"
        )
        model = train_model(pack, args, dev)
        target = pack["target"].detach().cpu().numpy()
        train = pack["g_tr"].detach().cpu().numpy().astype(bool)
        test = pack["g_te"].detach().cpu().numpy().astype(bool)

        constant_gate, constant_pred = constant_gate_prediction(pack)
        nonlinear_pred = model_prediction(model, pack)
        const_train = prediction_scores(constant_pred, target, train)
        const_test = prediction_scores(constant_pred, target, test)
        nl_train = prediction_scores(nonlinear_pred, target, train)
        nl_test = prediction_scores(nonlinear_pred, target, test)
        gate_df, gate_stats = gate_samples(model, pack, stratum, args.grid_size)
        neff = effective_n_lag1(pack)

        row = {
            "stratum": stratum,
            "rows": int(pack["n_groups"]),
            "frames": int(pack["n_frames"]),
            "test_rows": int(test.sum()),
            "train_rows": int(train.sum()),
            "split_strategy": pack["split_strategy"],
            "purged_train_frames": pack["purged_train_frames"],
            "cross_split_lag1_pairs": pack["cross_split_lag1_pairs"],
            "atoms": int(pack["n_atoms"]),
            "effective_n_lag1": neff["effective_n_lag1"],
            "median_lag1_rho": neff["median_lag1_rho"],
            "constant_gate": constant_gate,
            "constant_train_r2": const_train["r2"],
            "constant_train_component_r": const_train["component_r"],
            "constant_train_magnitude_r": const_train["magnitude_r"],
            "constant_test_r2": const_test["r2"],
            "constant_test_component_r": const_test["component_r"],
            "constant_test_magnitude_r": const_test["magnitude_r"],
            "nonlinear_train_r2": nl_train["r2"],
            "nonlinear_train_component_r": nl_train["component_r"],
            "nonlinear_train_magnitude_r": nl_train["magnitude_r"],
            "nonlinear_test_r2": nl_test["r2"],
            "nonlinear_test_component_r": nl_test["component_r"],
            "nonlinear_test_magnitude_r": nl_test["magnitude_r"],
            **gate_stats,
        }
        row["nonlinear_minus_constant_test_r2"] = row["nonlinear_test_r2"] - row["constant_test_r2"]
        row["gate_form"] = gate_form_label(row)
        summary_rows.append(row)
        sample_frames.append(gate_df)
        print(
            "   "
            f"constant R2={finite_fmt(row['constant_test_r2'])} "
            f"nonlinear R2={finite_fmt(row['nonlinear_test_r2'])} "
            f"gate={row['gate_form']} "
            f"N_eff_lag1={finite_fmt(row['effective_n_lag1'])}"
        )

    summary = pd.DataFrame(summary_rows)
    samples = pd.concat(sample_frames, ignore_index=True) if sample_frames else pd.DataFrame()
    summary_path = evidence_dir / "efg_distill_summary.csv"
    samples_path = evidence_dir / "efg_gate_samples.csv"
    json_path = evidence_dir / "efg_distill_run.json"
    summary.to_csv(summary_path, index=False)
    samples.to_csv(samples_path, index=False)
    json_path.write_text(
        json.dumps(
            {
                "out_dir": args.out_dir,
                "epochs": args.epochs,
                "hidden": args.hidden,
                "lr": args.lr,
                "seed": args.seed,
                "split": args.split,
                "purge_frames": args.purge_frames,
                "change_of_basis_orthogonality_max": float(np.abs(C.T @ C - np.eye(5)).max()),
                "summary_csv": str(summary_path),
                "gate_samples_csv": str(samples_path),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"\nwrote {summary_path}")
    print(f"wrote {samples_path}")
    print(f"wrote {json_path}")
    if not summary.empty:
        print("\nstratum,gate_form,atoms,effective_n_lag1,constant_test_R2,nonlinear_test_R2,gate_rel_span_5_95")
        for row in summary_rows:
            print(
                ",".join(
                    [
                        row["stratum"],
                        row["gate_form"],
                        str(row["atoms"]),
                        finite_fmt(row["effective_n_lag1"]),
                        finite_fmt(row["constant_test_r2"]),
                        finite_fmt(row["nonlinear_test_r2"]),
                        finite_fmt(row["gate_rel_span_5_95"]),
                    ]
                )
            )


if __name__ == "__main__":
    main()
