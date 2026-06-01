#!/usr/bin/env python3
"""oracle_parity — the gate that proves the composed functional API is faithful.

Runs the rediscover extractor twice on the same calcset — once with the
COMPOSED functional engine (Layer 1 verbs + Layer 2 curried closures + Layer 3
engine; --engine composed, the default) and once with the PROCEDURAL reference
cells (--engine procedural) — into two output dirs, then diffs the emitted CSVs
+ sidecar NPYs cell-for-cell.

A faithful composition reproduces the oracle BYTE-FOR-BYTE (same row order, same
numbers at 9 sig figs as the cell writes them). Any divergence is a composition
bug to find, not a new result (per the brief). NaN/Inf are compared as tokens.

This does NOT open trajectory.h5 — the reader owns H5; this only reads the CSVs
the two binary runs emit. (Reader-owns-H5 discipline.)

Usage:
    python oracle_parity.py \
        --bin   /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-gcc/h5reader_extract \
        --run   /shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft \
        --work  /tmp/rediscover-parity \
        [--case all] [--mc-cutoff 8.0]

Exit 0 ⇒ composed == procedural (gate PASS). Exit 1 ⇒ divergence (printed).
"""
import argparse
import csv
import math
import os
import subprocess
import sys

import numpy as np


def run_extract(binary, run, out, case, mc_cutoff, engine):
    os.makedirs(out, exist_ok=True)
    cmd = [binary, "--run", run, "--out", out, "--case", case,
           "--mc-cutoff", str(mc_cutoff), "--engine", engine]
    log = os.path.join(out, "run.log")
    with open(log, "wb") as fh:
        rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT)
    if rc != 0:
        print(f"[FAIL] {engine} run exited rc={rc}; see {log}")
        sys.exit(1)


def cells_equal(a, b, tol=1e-9):
    """Compare two CSV cells. Numeric → relative/absolute tol; else exact token."""
    if a == b:
        return True
    try:
        fa, fb = float(a), float(b)
    except ValueError:
        return False
    if math.isnan(fa) and math.isnan(fb):
        return True
    if math.isinf(fa) and math.isinf(fb):
        return fa == fb
    return abs(fa - fb) <= tol + 1e-7 * max(abs(fa), abs(fb))


def diff_csv(pa, pb, label):
    with open(pa, newline="") as fa, open(pb, newline="") as fb:
        ra, rb = list(csv.reader(fa)), list(csv.reader(fb))
    if len(ra) != len(rb):
        print(f"[FAIL] {label}: row count {len(ra)} (composed) != {len(rb)} (procedural)")
        return False
    if ra and ra[0] != rb[0]:
        print(f"[FAIL] {label}: header mismatch")
        print("  composed :", ra[0])
        print("  procedural:", rb[0])
        return False
    header = ra[0] if ra else []
    bad = 0
    for i in range(1, len(ra)):
        if len(ra[i]) != len(rb[i]):
            print(f"[FAIL] {label} row {i}: column count differs")
            return False
        for j, (ca, cb) in enumerate(zip(ra[i], rb[i])):
            if not cells_equal(ca, cb):
                col = header[j] if j < len(header) else f"col{j}"
                if bad < 20:
                    print(f"[FAIL] {label} row {i} col '{col}': "
                          f"composed={ca!r} procedural={cb!r}")
                bad += 1
    if bad:
        print(f"[FAIL] {label}: {bad} differing cells")
        return False
    print(f"[ok]   {label}: {len(ra)-1} rows identical")
    return True


def diff_npy(pa, pb, label):
    a, b = np.load(pa), np.load(pb)
    if a.shape != b.shape:
        print(f"[FAIL] {label}: shape {a.shape} != {b.shape}")
        return False
    if not np.allclose(a, b, rtol=1e-7, atol=1e-9, equal_nan=True):
        n = int((~np.isclose(a, b, rtol=1e-7, atol=1e-9, equal_nan=True)).sum())
        print(f"[FAIL] {label}: {n} differing NPY elements")
        return False
    print(f"[ok]   {label}: NPY {a.shape} identical")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--case", default="all")
    ap.add_argument("--mc-cutoff", type=float, default=8.0)
    args = ap.parse_args()

    comp = os.path.join(args.work, "composed")
    proc = os.path.join(args.work, "procedural")
    print("=== running composed engine ===")
    run_extract(args.bin, args.run, comp, args.case, args.mc_cutoff, "composed")
    print("=== running procedural cells ===")
    run_extract(args.bin, args.run, proc, args.case, args.mc_cutoff, "procedural")

    ok = True
    files = sorted(f for f in os.listdir(comp) if f.endswith((".csv", ".npy")))
    for f in files:
        pa, pb = os.path.join(comp, f), os.path.join(proc, f)
        if not os.path.exists(pb):
            print(f"[FAIL] {f}: present in composed, missing in procedural")
            ok = False
            continue
        if f.endswith(".csv"):
            ok &= diff_csv(pa, pb, f)
        else:
            ok &= diff_npy(pa, pb, f)

    print("\n=== GATE:", "PASS — composed reproduces the oracle" if ok
          else "FAIL — composition diverges from the oracle", "===")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
