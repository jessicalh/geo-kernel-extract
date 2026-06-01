#!/usr/bin/env python3
"""backbone_pysr_distill -- symbolic fits for the broad-backbone readout.

Run this with the PySR venv:

    analysis/venv/bin/python backbone_pysr_distill.py /tmp/rdc-backbone-law-evidence

It consumes readout_samples.csv.gz produced by backbone_distill_evidence.py. The
target is the learned per-source radial gate selected by that script:
ring/bond use the axis gate when axes are active, charge uses the displacement
gate. These are learned-model readouts, not recomputed kernels.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from pysr import PySRRegressor


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("evidence_dir", nargs="?", default="/tmp/rdc-backbone-law-evidence")
    ap.add_argument("--niterations", type=int, default=45)
    ap.add_argument("--timeout", type=int, default=240)
    ap.add_argument("--sample", type=int, default=6000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--bond-per-category", action="store_true")
    ap.add_argument("--target-column", default="gate_chosen")
    return ap.parse_args()


def score(y, pred):
    ok = np.isfinite(y) & np.isfinite(pred)
    if ok.sum() < 3:
        return float("nan")
    yy = y[ok]
    pp = pred[ok]
    den = ((yy - yy.mean()) ** 2).sum()
    if den <= 0:
        return float("nan")
    return float(1.0 - ((yy - pp) ** 2).sum() / den)


def make_model(args, variable_names):
    return PySRRegressor(
        niterations=args.niterations,
        populations=20,
        population_size=36,
        maxsize=22,
        binary_operators=["+", "-", "*", "/"],
        unary_operators=["square"],
        elementwise_loss="L2DistLoss()",
        constraints={"/": (-1, 9)},
        timeout_in_seconds=args.timeout,
        progress=False,
        deterministic=True,
        parallelism="serial",
        random_state=args.seed,
        verbosity=0,
        temp_equation_file=True,
    )


def fit_one(args, evidence_dir, label, d, variables):
    d = d.copy()
    d = d[np.isfinite(d[args.target_column])]
    for col in variables:
        d = d[np.isfinite(d[col])]
    if len(d) < 200:
        return None
    rng = np.random.default_rng(args.seed)
    if len(d) > args.sample:
        d = d.iloc[rng.choice(len(d), size=args.sample, replace=False)].copy()
    X = d[variables].to_numpy(float)
    y = d[args.target_column].to_numpy(float)
    print(f"\n[{label}] PySR rows={len(d)} vars={variables} target={args.target_column}", flush=True)
    model = make_model(args, variables)
    model.fit(X, y, variable_names=variables)
    eq = model.equations_.copy()
    eq.insert(0, "label", label)
    eq.insert(1, "target", args.target_column)
    eq.to_csv(evidence_dir / f"pysr_equations_{args.target_column}_{label}.csv", index=False)
    best = model.get_best()
    pred = model.predict(X)
    r2 = score(y, pred)
    print(f"[{label}] best={best['equation']} R2={r2:+.4f}", flush=True)
    return {
        "label": label,
        "target": args.target_column,
        "rows": int(len(d)),
        "variables": ",".join(variables),
        "best_equation": str(best["equation"]),
        "best_complexity": int(best["complexity"]),
        "best_loss": float(best["loss"]),
        "r2": r2,
    }


def main():
    args = parse_args()
    evidence_dir = Path(args.evidence_dir)
    readout_path = evidence_dir / "readout_samples.csv.gz"
    readout = pd.read_csv(readout_path)
    out = []

    ring = readout[(readout.kind == "ring") & (readout.r >= 3.0) & (readout.ring_intensity != 0)]
    res = fit_one(args, evidence_dir, "ring", ring, ["r", "cos_theta", "ring_intensity"])
    if res:
        out.append(res)

    bond = readout[(readout.kind == "bond") & (readout.r >= 2.0) & (readout.bond_category >= 0)]
    res = fit_one(args, evidence_dir, "bond_pooled", bond, ["r", "cos_theta", "bond_category"])
    if res:
        out.append(res)
    if args.bond_per_category:
        for cat in sorted(bond.bond_category.dropna().unique()):
            d = bond[bond.bond_category == cat]
            res = fit_one(args, evidence_dir, f"bond_cat{int(cat)}", d, ["r", "cos_theta"])
            if res:
                out.append(res)

    charge = readout[(readout.kind == "charge") & (readout.r >= 1.5) & (np.abs(readout.source_q_e) > 1e-5)]
    res = fit_one(args, evidence_dir, "charge", charge, ["r", "source_q_e"])
    if res:
        out.append(res)

    summary = pd.DataFrame(out)
    summary_path = evidence_dir / f"pysr_summary_{args.target_column}.csv"
    summary.to_csv(summary_path, index=False)
    (evidence_dir / "pysr_manifest.json").write_text(
        json.dumps(
            {
                "niterations": args.niterations,
                "timeout": args.timeout,
                "sample": args.sample,
                "seed": args.seed,
                "bond_per_category": args.bond_per_category,
                "target_column": args.target_column,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"\nWrote {summary_path}", flush=True)


if __name__ == "__main__":
    main()
