#!/usr/bin/env python3
"""
Batch feature extraction for calibration proteins.

Runs nmr_extract --mutant on each WT+ALA pair from calibration/{ID}/,
producing NPY feature directories under calibration/features/{run}/{ID}/.

Each protein dir has real files (no symlinks):
    {ID}_WT.xyz, {ID}_WT.prmtop, {ID}_WT_nmr.out
    {ID}_ALA.xyz, {ID}_ALA.prmtop, {ID}_ALA_nmr.out

Usage:
    python learn/extract.py --run GatedCalibration
    python learn/extract.py --run GatedCalibration --protein P84477
    python learn/extract.py --run GatedCalibration --dry-run
    python learn/extract.py --run GatedCalibration --resume
    python learn/extract.py --run GatedCalibration --config params.toml
"""

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Paths
REPO = Path(__file__).resolve().parent.parent
CALIBRATION = REPO / "calibration"
FEATURES_BASE = CALIBRATION / "features"
NMR_EXTRACT = REPO / "build" / "nmr_extract"

KNOWN_RUNS = [
    "FirstExtraction",
    "CalibrationExtractionTest",
    "CalibrationExtractionKeep",
    "GatedCalibration",
    "SubmissionExtraction",
    "AzimuthalExtraction",
]


def flush_vcache():
    """Flush page cache to free unified memory for GPU on DGX Spark."""
    try:
        subprocess.run(["sudo", "sh", "-c", "sync; echo 3 > /proc/sys/vm/drop_caches"],
                       capture_output=True, timeout=10)
    except Exception:
        pass


def is_complete_pair(protein_dir: Path, protein_id: str) -> tuple[bool, str]:
    """Check that a calibration protein has the required extraction files."""
    for variant in ("WT", "ALA"):
        root = protein_dir / f"{protein_id}_{variant}"
        for ext in (".xyz", ".prmtop"):
            path = Path(str(root) + ext)
            if not path.exists():
                return False, f"missing {path.name}"
    return True, ""


def extract_one(protein_id: str, config_path: str | None,
                dry_run: bool = False) -> dict:
    """Extract features for one WT+ALA protein pair using --mutant mode.

    Returns dict with keys: protein_id, ok, error, elapsed, n_arrays.
    """
    result = {"protein_id": protein_id,
              "ok": False, "error": "", "elapsed": 0.0, "n_arrays": 0}

    protein_dir = CALIBRATION / protein_id
    ok, reason = is_complete_pair(protein_dir, protein_id)
    if not ok:
        result["error"] = reason
        return result

    wt_root  = str(protein_dir / f"{protein_id}_WT")
    ala_root = str(protein_dir / f"{protein_id}_ALA")
    output_dir = FEATURES_BASE / _current_run / protein_id

    if dry_run:
        result["ok"] = True
        result["error"] = f"[dry-run] --wt {wt_root}  --ala {ala_root}  --output {output_dir}"
        return result

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [str(NMR_EXTRACT), "--mutant",
           "--wt", wt_root, "--ala", ala_root,
           "--output", str(output_dir)]
    if config_path:
        cmd.extend(["--config", config_path])

    try:
        t0 = time.time()
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        result["elapsed"] = time.time() - t0

        if proc.returncode != 0:
            result["error"] = proc.stderr.strip()[-500:] if proc.stderr else f"exit code {proc.returncode}"
            return result

        npy_files = list(output_dir.glob("*.npy"))
        result["n_arrays"] = len(npy_files)
        result["ok"] = True

        for line in (proc.stderr or "").splitlines():
            if "Wrote" in line and "arrays" in line:
                result["error"] = line.strip()
                break

    except subprocess.TimeoutExpired:
        result["error"] = "timeout (600s)"
    except Exception as e:
        result["error"] = str(e)

    return result


def is_extracted(protein_id: str) -> bool:
    """Check if features already exist for this protein."""
    d = FEATURES_BASE / _current_run / protein_id
    return (d / "pos.npy").exists()

# Set by main() from --run argument
_current_run = ""


def main():
    parser = argparse.ArgumentParser(description="Batch feature extraction")
    parser.add_argument("--run", type=str, required=True,
                        help="extraction run name (e.g. FirstExtraction)")
    parser.add_argument("--workers", type=int, default=1,
                        help="parallel workers (default 1, xTB is single-process)")
    parser.add_argument("--protein", type=str, default=None,
                        help="extract single protein")
    parser.add_argument("--dry-run", action="store_true",
                        help="show what would run")
    parser.add_argument("--resume", action="store_true",
                        help="skip already-extracted")
    parser.add_argument("--config", type=str, default=None,
                        help="TOML config file for calculator parameter overrides")
    args = parser.parse_args()

    global _current_run
    _current_run = args.run
    if args.run not in KNOWN_RUNS:
        print(f"WARNING: '{args.run}' not in known runs: {KNOWN_RUNS}",
              file=sys.stderr)

    if not NMR_EXTRACT.exists():
        print(f"ERROR: {NMR_EXTRACT} not found. Build first.", file=sys.stderr)
        sys.exit(1)

    if not CALIBRATION.is_dir():
        print(f"ERROR: {CALIBRATION} not found.", file=sys.stderr)
        sys.exit(1)

    # Build job list from calibration directory
    if args.protein:
        proteins = [args.protein]
    else:
        proteins = sorted(d.name for d in CALIBRATION.iterdir()
                          if d.is_dir() and not d.name.startswith("."))

    jobs = []
    skipped = 0
    for pid in proteins:
        if args.resume and is_extracted(pid):
            skipped += 1
            continue
        jobs.append(pid)

    print(f"Proteins: {len(proteins)}  Jobs: {len(jobs)}  Skipped: {skipped}")

    if not jobs:
        print("Nothing to do.")
        return

    run_dir = FEATURES_BASE / _current_run
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "extract_log.jsonl"

    ok_count = 0
    fail_count = 0
    t_start = time.time()

    if args.workers <= 1 or args.dry_run:
        for i, pid in enumerate(jobs):
            flush_vcache()
            r = extract_one(pid, args.config, dry_run=args.dry_run)
            status = "OK" if r["ok"] else "FAIL"
            print(f"[{i+1}/{len(jobs)}] {pid}: {status}  "
                  f"{r['elapsed']:.1f}s  {r['n_arrays']} arrays"
                  f"{'  ' + r['error'] if r['error'] else ''}")
            if r["ok"]:
                ok_count += 1
            else:
                fail_count += 1
            if not args.dry_run:
                with open(log_path, "a") as f:
                    f.write(json.dumps(r) + "\n")
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(extract_one, pid, args.config): pid
                       for pid in jobs}
            for i, future in enumerate(as_completed(futures)):
                pid = futures[future]
                try:
                    r = future.result()
                except Exception as e:
                    r = {"protein_id": pid,
                         "ok": False, "error": str(e),
                         "elapsed": 0.0, "n_arrays": 0}

                status = "OK" if r["ok"] else "FAIL"
                print(f"[{i+1}/{len(jobs)}] {pid}: {status}  "
                      f"{r['elapsed']:.1f}s  {r['n_arrays']} arrays"
                      f"{'  ' + r['error'] if r['error'] else ''}")
                if r["ok"]:
                    ok_count += 1
                else:
                    fail_count += 1
                with open(log_path, "a") as f:
                    f.write(json.dumps(r) + "\n")

    elapsed = time.time() - t_start
    print(f"\nDone: {ok_count} OK, {fail_count} failed, {elapsed:.0f}s total")


if __name__ == "__main__":
    main()
