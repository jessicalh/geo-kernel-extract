#!/usr/bin/env python3
"""
Batch feature extraction for all 725 proteins.

Runs nmr_extract --orca independently on each WT and ALA conformation,
producing ~1450 feature directories under learn/features/{run}/{protein_id}/{wt,ala}/.

Extraction runs:
    FirstExtraction           — inaugural run, default parameters
    CalibrationExtractionTest — sweep experiments, exploratory
    CalibrationExtractionKeep — validated/optimised parameters
    SubmissionExtraction      — final thesis run

Usage:
    python learn/extract.py --run FirstExtraction
    python learn/extract.py --run FirstExtraction --workers 1
    python learn/extract.py --run FirstExtraction --protein P84477
    python learn/extract.py --run FirstExtraction --dry-run
    python learn/extract.py --run FirstExtraction --resume
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Paths
REPO = Path(__file__).resolve().parent.parent
CONSOLIDATED = Path("/shared/2026Thesis/consolidated")
FEATURES_BASE = REPO / "learn" / "features"
NMR_EXTRACT = REPO / "build" / "nmr_extract"

KNOWN_RUNS = [
    "FirstExtraction",
    "CalibrationExtractionTest",
    "CalibrationExtractionKeep",
    "SubmissionExtraction",
]


def find_files(protein_dir: Path, prefix: str) -> dict:
    """Find the WT or ALA files for a protein by prefix (e.g. 'P84477_WT')."""
    files = {}
    for f in protein_dir.iterdir():
        name = f.name
        if not name.startswith(prefix):
            continue
        if name.endswith(".pdb") and "_amber" not in name and "_water" not in name:
            files["pdb"] = f
        elif name.endswith(".xyz"):
            files["xyz"] = f
        elif name.endswith(".prmtop"):
            files["prmtop"] = f
        elif name.endswith("_nmr.out"):
            files["nmr_out"] = f
    return files


def make_symlink_dir(files: dict, tmpdir: str) -> str:
    """Create a temp dir with symlinks matching what FindOrcaFiles expects."""
    d = tempfile.mkdtemp(dir=tmpdir)
    for key, src in files.items():
        ext_map = {"pdb": ".pdb", "xyz": ".xyz", "prmtop": ".prmtop",
                    "nmr_out": "_nmr.out"}
        # Use a simple name so FindOrcaFiles picks it up unambiguously
        linkname = f"input{ext_map[key]}"
        os.symlink(src, os.path.join(d, linkname))
    return d


def flush_vcache():
    """Flush page cache to free unified memory for GPU on DGX Spark."""
    try:
        subprocess.run(["sudo", "sh", "-c", "sync; echo 3 > /proc/sys/vm/drop_caches"],
                       capture_output=True, timeout=10)
    except Exception:
        pass


def extract_one(protein_id: str, variant: str, dry_run: bool = False) -> dict:
    """Extract features for one protein variant (wt or ala).

    Returns dict with keys: protein_id, variant, ok, error, elapsed, n_arrays.
    """
    result = {"protein_id": protein_id, "variant": variant,
              "ok": False, "error": "", "elapsed": 0.0, "n_arrays": 0}

    protein_dir = CONSOLIDATED / protein_id
    prefix = f"{protein_id}_{'WT' if variant == 'wt' else 'ALA'}"
    files = find_files(protein_dir, prefix)

    # Need at minimum xyz and prmtop
    for required in ("xyz", "prmtop"):
        if required not in files:
            result["error"] = f"missing {required} for {prefix}"
            return result

    output_dir = FEATURES_BASE / _current_run / protein_id / variant

    if dry_run:
        result["ok"] = True
        result["error"] = f"[dry-run] would extract to {output_dir}"
        return result

    output_dir.mkdir(parents=True, exist_ok=True)

    # Create temp dir with unambiguous symlinks
    tmpdir = tempfile.mkdtemp(prefix=f"nmr_{protein_id}_{variant}_")
    try:
        link_dir = make_symlink_dir(files, tmpdir)

        cmd = [str(NMR_EXTRACT), "--orca", link_dir,
               "--output", str(output_dir)]

        t0 = time.time()
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        result["elapsed"] = time.time() - t0

        if proc.returncode != 0:
            result["error"] = proc.stderr.strip()[-500:] if proc.stderr else f"exit code {proc.returncode}"
            return result

        # Count output arrays
        npy_files = list(output_dir.glob("*.npy"))
        result["n_arrays"] = len(npy_files)
        result["ok"] = True

        # Parse array count from stderr if available
        for line in (proc.stderr or "").splitlines():
            if "Wrote" in line and "arrays" in line:
                result["error"] = line.strip()  # not really error, just info
                break

    except subprocess.TimeoutExpired:
        result["error"] = "timeout (600s)"
    except Exception as e:
        result["error"] = str(e)
    finally:
        # Clean up temp dir
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)

    return result


def is_extracted(protein_id: str, variant: str) -> bool:
    """Check if features already exist (at least pos.npy)."""
    return (FEATURES_BASE / _current_run / protein_id / variant / "pos.npy").exists()

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
    parser.add_argument("--variant", choices=["wt", "ala", "both"],
                        default="both", help="which variant(s)")
    args = parser.parse_args()

    global _current_run
    _current_run = args.run
    if args.run not in KNOWN_RUNS:
        print(f"WARNING: '{args.run}' not in known runs: {KNOWN_RUNS}",
              file=sys.stderr)

    if not NMR_EXTRACT.exists():
        print(f"ERROR: {NMR_EXTRACT} not found. Build first.", file=sys.stderr)
        sys.exit(1)

    # Build job list
    if args.protein:
        proteins = [args.protein]
    else:
        proteins = sorted(d.name for d in CONSOLIDATED.iterdir() if d.is_dir())

    variants = ["wt", "ala"] if args.variant == "both" else [args.variant]

    jobs = []
    skipped = 0
    for pid in proteins:
        for v in variants:
            if args.resume and is_extracted(pid, v):
                skipped += 1
                continue
            jobs.append((pid, v))

    print(f"Proteins: {len(proteins)}  Variants: {len(variants)}  "
          f"Jobs: {len(jobs)}  Skipped: {skipped}")

    if not jobs:
        print("Nothing to do.")
        return

    # Log file — per run
    run_dir = FEATURES_BASE / _current_run
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "extract_log.jsonl"

    ok_count = 0
    fail_count = 0
    t_start = time.time()

    if args.workers <= 1 or args.dry_run:
        # Sequential
        for i, (pid, v) in enumerate(jobs):
            flush_vcache()
            r = extract_one(pid, v, dry_run=args.dry_run)
            status = "OK" if r["ok"] else "FAIL"
            print(f"[{i+1}/{len(jobs)}] {pid}/{v}: {status}  "
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
        # Parallel
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(extract_one, pid, v): (pid, v)
                       for pid, v in jobs}
            for i, future in enumerate(as_completed(futures)):
                pid, v = futures[future]
                try:
                    r = future.result()
                except Exception as e:
                    r = {"protein_id": pid, "variant": v,
                         "ok": False, "error": str(e),
                         "elapsed": 0.0, "n_arrays": 0}

                status = "OK" if r["ok"] else "FAIL"
                print(f"[{i+1}/{len(jobs)}] {pid}/{v}: {status}  "
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
