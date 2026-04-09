#!/usr/bin/env python3
"""
validate_smoke.py — Deep validation of smoke test output.

Reads the NPY files and log.jsonl produced by test_smoke.cpp and performs
semantic validation that C++ cannot easily do:
  - Load every .npy via numpy, check shapes, dtypes, NaN/Inf counts
  - Verify atom count consistency across all per-atom arrays
  - Check physical bounds (charges, positions, shielding magnitudes)
  - Parse log.jsonl for errors and timing
  - Binary-compare two directories (run vs blessed)

Usage:
    python validate_smoke.py <run_dir>                    # validate one run
    python validate_smoke.py <run_dir> --blessed <dir>    # validate + compare
    python validate_smoke.py --latest                     # validate most recent smoke run
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def load_npy_files(directory: Path) -> dict[str, np.ndarray]:
    """Load all .npy files in directory, return {name: array}."""
    arrays = {}
    for p in sorted(directory.glob("*.npy")):
        try:
            arrays[p.name] = np.load(p)
        except Exception as e:
            print(f"  FAIL: could not load {p.name}: {e}")
    return arrays


def validate_arrays(arrays: dict[str, np.ndarray], label: str) -> int:
    """Validate array shapes, dtypes, NaN/Inf. Returns error count."""
    errors = 0
    print(f"\n  [{label}] Loaded {len(arrays)} arrays")

    # Determine atom count from pos.npy
    if "pos.npy" not in arrays:
        print("  FAIL: pos.npy missing — cannot determine atom count")
        return 1

    pos = arrays["pos.npy"]
    N = pos.shape[0]
    print(f"  Atom count (from pos.npy): {N}")

    # Per-atom arrays must have shape[0] == N
    per_atom_files = [
        "pos.npy", "element.npy", "residue_index.npy", "residue_type.npy",
        "bs_shielding.npy", "hm_shielding.npy", "mc_shielding.npy",
        "pq_shielding.npy", "disp_shielding.npy", "ringchi_shielding.npy",
        "coulomb_shielding.npy", "hbond_shielding.npy",
        "coulomb_E.npy", "apbs_E.npy", "bs_total_B.npy",
        "mopac_charges.npy",
        "mopac_coulomb_shielding.npy", "mopac_mc_shielding.npy",
    ]
    for name in per_atom_files:
        if name not in arrays:
            continue
        arr = arrays[name]
        if arr.shape[0] != N:
            print(f"  FAIL: {name} shape[0]={arr.shape[0]} != atom count {N}")
            errors += 1

    # Check for NaN/Inf in all float arrays
    for name, arr in arrays.items():
        if arr.dtype.kind == "f":
            nan_count = np.isnan(arr).sum()
            inf_count = np.isinf(arr).sum()
            if nan_count > 0:
                print(f"  FAIL: {name} has {nan_count} NaN values")
                errors += 1
            if inf_count > 0:
                print(f"  FAIL: {name} has {inf_count} Inf values")
                errors += 1

    # Physical bounds checks
    if "mopac_charges.npy" in arrays:
        q = arrays["mopac_charges.npy"]
        if np.abs(q).max() > 5.0:
            print(f"  WARN: mopac_charges max |q| = {np.abs(q).max():.3f}")

    if "pos.npy" in arrays:
        r = np.linalg.norm(pos, axis=1)
        if r.max() > 500:
            print(f"  WARN: positions max |r| = {r.max():.1f} A — very large")

    # Shielding tensor arrays should be shape (N, 9)
    tensor_files = [
        "bs_shielding.npy", "hm_shielding.npy", "mc_shielding.npy",
        "pq_shielding.npy", "disp_shielding.npy", "ringchi_shielding.npy",
        "coulomb_shielding.npy", "hbond_shielding.npy",
    ]
    for name in tensor_files:
        if name not in arrays:
            continue
        arr = arrays[name]
        if arr.ndim != 2 or arr.shape[1] != 9:
            print(f"  FAIL: {name} shape={arr.shape}, expected (N, 9)")
            errors += 1

    # Ring contributions sparse array
    if "ring_contributions.npy" in arrays:
        rc = arrays["ring_contributions.npy"]
        print(f"  Ring contributions: {rc.shape[0]} (atom,ring) pairs, "
              f"{rc.shape[1]} columns")
        if rc.shape[1] != 59:
            print(f"  FAIL: ring_contributions columns={rc.shape[1]}, expected 59")
            errors += 1

    # Summary
    total_bytes = sum(arr.nbytes for arr in arrays.values())
    print(f"  Total data: {total_bytes / 1024:.1f} KB across {len(arrays)} arrays")

    if errors == 0:
        print(f"  [{label}] ALL CHECKS PASSED")
    else:
        print(f"  [{label}] {errors} ERRORS")

    return errors


def validate_sdk_load(directory: Path, label: str) -> int:
    """Load extraction via nmr_extract SDK and validate typed wrappers.

    This is the trust-boundary test: C++ wrote NPY, SDK reads it back
    through typed classes. If the column counts or semantics drift,
    this catches it.
    """
    try:
        from nmr_extract import load
    except ImportError:
        print(f"  [{label}] SKIP SDK validation — nmr_extract not installed")
        return 0

    errors = 0
    print(f"\n  [{label}] SDK load validation")

    try:
        p = load(str(directory))
    except Exception as e:
        print(f"  FAIL: SDK load raised {type(e).__name__}: {e}")
        return 1

    print(f"  SDK loaded: {p.n_atoms} atoms, id={p.protein_id}")

    # ── Ring contributions: cos_phi / sin_phi ──────────────────────
    rc = p.ring_contributions
    if rc is not None and rc.n_pairs > 0:
        cos_phi = rc.cos_phi
        sin_phi = rc.sin_phi
        if cos_phi.shape != (rc.n_pairs,):
            print(f"  FAIL: cos_phi shape {cos_phi.shape}, expected ({rc.n_pairs},)")
            errors += 1
        if sin_phi.shape != (rc.n_pairs,):
            print(f"  FAIL: sin_phi shape {sin_phi.shape}, expected ({rc.n_pairs},)")
            errors += 1

        # cos²φ + sin²φ ≈ 1 for atoms not on the ring axis
        rho = rc.rho
        off_axis = rho > 0.1
        if off_axis.sum() > 0:
            unit_check = cos_phi[off_axis]**2 + sin_phi[off_axis]**2
            max_err = np.abs(unit_check - 1.0).max()
            if max_err > 1e-10:
                print(f"  FAIL: cos²φ+sin²φ max error = {max_err:.2e}")
                errors += 1
            else:
                print(f"  cos²φ+sin²φ: max error {max_err:.2e} "
                      f"({off_axis.sum()}/{rc.n_pairs} off-axis pairs)")

        # On-axis atoms (rho ≈ 0) should have default (1, 0)
        on_axis = rho <= 0.1
        if on_axis.sum() > 0:
            if not np.allclose(cos_phi[on_axis], 1.0):
                print(f"  FAIL: on-axis cos_phi not 1.0")
                errors += 1
            if not np.allclose(sin_phi[on_axis], 0.0):
                print(f"  FAIL: on-axis sin_phi not 0.0")
                errors += 1
    else:
        print("  SKIP: no ring contributions (no aromatic residues)")

    # ── MOPAC: valency + dipole ────────────────────────────────────
    if p.mopac:
        # Valency
        val = p.mopac.core.scalars.valency
        if val.shape != (p.n_atoms,):
            print(f"  FAIL: valency shape {val.shape}, expected ({p.n_atoms},)")
            errors += 1
        elif val.max() < 0.1:
            print(f"  FAIL: valency max={val.max():.4f} — suspiciously small")
            errors += 1
        else:
            print(f"  Valency: range [{val.min():.2f}, {val.max():.2f}]")

        # Dipole
        dip = p.mopac.core.global_.dipole
        if dip.shape != (3,):
            print(f"  FAIL: dipole shape {dip.shape}, expected (3,)")
            errors += 1
        dip_mag = p.mopac.core.global_.dipole_magnitude
        if dip_mag < 1e-10:
            print(f"  WARN: dipole magnitude = {dip_mag:.4f} D — zero or near-zero")
        else:
            print(f"  Dipole: |d| = {dip_mag:.2f} D, "
                  f"components = ({dip[0]:.2f}, {dip[1]:.2f}, {dip[2]:.2f})")

        # Heat of formation (already existed, verify accessible)
        hof = p.mopac.core.global_.heat_of_formation
        print(f"  Heat of formation: {hof:.1f} kcal/mol")
    else:
        print("  SKIP: no MOPAC data")

    if errors == 0:
        print(f"  [{label}] SDK VALIDATION PASSED")
    else:
        print(f"  [{label}] SDK VALIDATION: {errors} ERRORS")

    return errors


def validate_log(log_path: Path) -> int:
    """Parse log.jsonl and check for errors. Returns error count."""
    errors = 0
    entries = []
    with open(log_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                entries.append(json.loads(line))
            except json.JSONDecodeError:
                print(f"  WARN: unparseable log line: {line[:80]}")

    print(f"\n  Log: {len(entries)} entries")

    error_entries = [e for e in entries if e.get("level") == "ERROR"]
    if error_entries:
        print(f"  FAIL: {len(error_entries)} ERROR entries in log:")
        for e in error_entries:
            print(f"    [{e.get('op')}] {e.get('detail')}")
        errors += len(error_entries)
    else:
        print("  Log: no errors")

    # Timing summary from [END] entries
    timings = {}
    for e in entries:
        op = e.get("op", "")
        detail = e.get("detail", "")
        if "[END]" in op and "elapsed=" in detail:
            name = op.replace(" [END]", "")
            ms = detail.split("elapsed=")[1].rstrip("ms")
            try:
                timings[name] = int(ms)
            except ValueError:
                pass

    if timings:
        print("  Timings:")
        for name, ms in sorted(timings.items(), key=lambda x: -x[1]):
            print(f"    {name}: {ms}ms")

    return errors


def binary_compare(run_dir: Path, blessed_dir: Path) -> int:
    """Compare all .npy files between two directories. Returns diff count."""
    run_files = {p.name for p in run_dir.glob("*.npy")}
    blessed_files = {p.name for p in blessed_dir.glob("*.npy")}

    errors = 0
    missing = blessed_files - run_files
    new = run_files - blessed_files

    if missing:
        print(f"  FAIL: missing from run: {missing}")
        errors += len(missing)
    if new:
        print(f"  NEW files (not in blessed): {new}")

    common = run_files & blessed_files
    identical = 0
    different = []

    for name in sorted(common):
        a = np.load(run_dir / name)
        b = np.load(blessed_dir / name)
        if a.shape != b.shape:
            different.append((name, f"shape {a.shape} vs {b.shape}"))
        elif a.dtype != b.dtype:
            different.append((name, f"dtype {a.dtype} vs {b.dtype}"))
        elif not np.array_equal(a, b):
            if a.dtype.kind == "f":
                max_diff = np.abs(a - b).max()
                different.append((name, f"max delta = {max_diff:.2e}"))
            else:
                n_diff = (a != b).sum()
                different.append((name, f"{n_diff} values differ"))
        else:
            identical += 1

    print(f"\n  Binary comparison: {identical} identical, "
          f"{len(different)} different")

    for name, reason in different:
        print(f"    DIFF: {name} — {reason}")
        errors += 1

    return errors


def find_latest_smoke(base: Path) -> Path | None:
    """Find the most recent timestamped smoke run."""
    smoke_dir = base / "tests" / "golden" / "smoke"
    if not smoke_dir.exists():
        return None
    dirs = sorted(smoke_dir.iterdir(), reverse=True)
    return dirs[0] if dirs else None


def main():
    parser = argparse.ArgumentParser(description="Validate smoke test output")
    parser.add_argument("run_dir", nargs="?", help="Path to smoke run directory")
    parser.add_argument("--blessed", help="Path to blessed baseline for comparison")
    parser.add_argument("--latest", action="store_true",
                        help="Validate most recent smoke run")
    args = parser.parse_args()

    if args.latest:
        project = Path(__file__).resolve().parent.parent
        latest = find_latest_smoke(project)
        if not latest:
            print("No smoke runs found")
            sys.exit(1)
        # Validate all subdirs (nodft, withdft) in the latest run
        total_errors = 0
        blessed_base = project / "tests" / "golden" / "blessed"
        for sub in sorted(latest.iterdir()):
            if not sub.is_dir():
                continue
            label = sub.name
            print(f"\n{'='*60}")
            print(f"  Validating: {sub}")
            print(f"{'='*60}")

            arrays = load_npy_files(sub)
            total_errors += validate_arrays(arrays, label)
            total_errors += validate_sdk_load(sub, label)

            log_path = sub / "log.jsonl"
            if log_path.exists():
                total_errors += validate_log(log_path)

            blessed = blessed_base / label
            if blessed.exists():
                total_errors += binary_compare(sub, blessed)

        sys.exit(1 if total_errors > 0 else 0)

    if not args.run_dir:
        parser.print_help()
        sys.exit(1)

    run_dir = Path(args.run_dir)
    if not run_dir.exists():
        print(f"Directory not found: {run_dir}")
        sys.exit(1)

    total_errors = 0

    arrays = load_npy_files(run_dir)
    total_errors += validate_arrays(arrays, run_dir.name)
    total_errors += validate_sdk_load(run_dir, run_dir.name)

    log_path = run_dir / "log.jsonl"
    if log_path.exists():
        total_errors += validate_log(log_path)

    if args.blessed:
        blessed = Path(args.blessed)
        if blessed.exists():
            total_errors += binary_compare(run_dir, blessed)
        else:
            print(f"Blessed directory not found: {blessed}")

    sys.exit(1 if total_errors > 0 else 0)


if __name__ == "__main__":
    main()
