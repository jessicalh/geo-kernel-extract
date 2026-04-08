#!/usr/bin/env python3
"""
Populate the calibration/ directory with symlinks from consolidated/.

Creates one subdirectory per complete WT+ALA protein pair.  The six files
that ExpandOrcaRoot needs sit at the protein root; everything else goes
into prep/ so protonation artifacts never pollute the extraction namespace.

Layout per protein:

    calibration/{ID}/
        {ID}_WT.xyz              extraction input (tleap coords)
        {ID}_WT.prmtop           extraction input (tleap topology)
        {ID}_WT_nmr.out          extraction input (Orca DFT tensors)
        {ID}_ALA.xyz
        {ID}_ALA.prmtop
        {ID}_ALA_nmr.out
        prep/                    protonation chain + Orca raw
            ... everything else
        metadata.json

Usage:
    python calibration/populate.py                 # populate (default)
    python calibration/populate.py --dry-run       # preview only
    python calibration/populate.py --verify        # check existing links
"""

import argparse
import os
import shutil
import sys
from pathlib import Path

CONSOLIDATED = Path("/shared/2026Thesis/consolidated")
CALIBRATION  = Path(__file__).resolve().parent

# The six files ExpandOrcaRoot constructs from a root name.
# These sit at the protein directory root.
EXTRACTION_SUFFIXES = {
    "_WT.xyz", "_WT.prmtop", "_WT_nmr.out",
    "_ALA.xyz", "_ALA.prmtop", "_ALA_nmr.out",
}

# Required for a pair to be "complete" (nmr.out is optional — DFT may be absent)
REQUIRED_SUFFIXES = {"_WT.xyz", "_WT.prmtop", "_ALA.xyz", "_ALA.prmtop"}


def classify_file(protein_id: str, filename: str) -> str:
    """Return 'extraction' if the file belongs at the protein root, else 'prep'."""
    for suffix in EXTRACTION_SUFFIXES:
        if filename == protein_id + suffix:
            return "extraction"
    return "prep"


def _find_timestamped_nmr(src_dir: Path, protein_id: str, variant: str) -> Path | None:
    """Find a timestamped _nmr.out file when the plain one is absent.

    Consolidated has both plain ({ID}_{VAR}_nmr.out) and timestamped
    ({ID}_{VAR}_YYYYMMDD_HHMMSS_nmr.out).  Many proteins only have
    the timestamped version.  Return the newest match, or None.
    """
    import re
    prefix = f"{protein_id}_{variant}_"
    pattern = re.compile(rf"^{re.escape(prefix)}\d{{8}}_\d{{6}}_nmr\.out$")
    matches = sorted(
        (f for f in src_dir.iterdir() if pattern.match(f.name)),
        key=lambda f: f.name, reverse=True,   # newest timestamp first
    )
    return matches[0] if matches else None


def is_complete(protein_dir: Path, protein_id: str) -> tuple[bool, list[str]]:
    """Check whether a protein has the minimum files for extraction."""
    missing = []
    for suffix in REQUIRED_SUFFIXES:
        if not (protein_dir / (protein_id + suffix)).exists():
            missing.append(protein_id + suffix)
    return len(missing) == 0, missing


def populate_one(protein_id: str, dry_run: bool = False) -> dict:
    """Create calibration/{ID}/ with symlinks for one protein."""
    src_dir = CONSOLIDATED / protein_id
    dst_dir = CALIBRATION / protein_id

    result = {"protein_id": protein_id, "ok": False,
              "extraction": 0, "prep": 0, "error": ""}

    complete, missing = is_complete(src_dir, protein_id)
    if not complete:
        result["error"] = f"incomplete: missing {', '.join(missing)}"
        return result

    if not dry_run:
        dst_dir.mkdir(exist_ok=True)
        (dst_dir / "prep").mkdir(exist_ok=True)

    # Track which extraction-level names we've placed
    placed_extraction = set()

    for src_file in sorted(src_dir.iterdir()):
        if src_file.is_dir():
            continue

        category = classify_file(protein_id, src_file.name)

        if category == "extraction":
            # Copy with canonical name — no symlink fragility
            dst = dst_dir / src_file.name
            placed_extraction.add(src_file.name)
            if not dry_run:
                if dst.is_symlink() or dst.exists():
                    dst.unlink()
                shutil.copy2(src_file, dst)
            result["extraction"] += 1
        else:
            # Prep files: symlink is fine for archival
            dst = dst_dir / "prep" / src_file.name
            if not dry_run:
                rel = os.path.relpath(src_file, dst.parent)
                if dst.is_symlink() or dst.exists():
                    dst.unlink()
                dst.symlink_to(rel)
            result["prep"] += 1

    # Fill in missing _nmr.out from timestamped versions — copy with
    # the canonical name so ExpandOrcaRoot finds it without guessing
    for variant in ("WT", "ALA"):
        canonical = f"{protein_id}_{variant}_nmr.out"
        if canonical not in placed_extraction:
            ts_file = _find_timestamped_nmr(src_dir, protein_id, variant)
            if ts_file is not None:
                dst = dst_dir / canonical
                if not dry_run:
                    if dst.is_symlink() or dst.exists():
                        dst.unlink()
                    shutil.copy2(ts_file, dst)
                result["extraction"] += 1

    result["ok"] = True
    return result


def verify_one(protein_id: str) -> list[str]:
    """Check that all symlinks in a calibration protein dir are valid."""
    dst_dir = CALIBRATION / protein_id
    broken = []
    if not dst_dir.is_dir():
        return [f"{protein_id}: directory missing"]

    for f in dst_dir.rglob("*"):
        if f.is_symlink() and not f.exists():
            broken.append(f"{f.relative_to(CALIBRATION)}: broken → {os.readlink(f)}")
    return broken


def main():
    parser = argparse.ArgumentParser(description="Populate calibration symlink tree")
    parser.add_argument("--dry-run", action="store_true", help="preview without creating")
    parser.add_argument("--verify", action="store_true", help="check existing links")
    parser.add_argument("--protein", type=str, help="single protein only")
    args = parser.parse_args()

    # Discover proteins
    if args.protein:
        proteins = [args.protein]
    else:
        proteins = sorted(d.name for d in CONSOLIDATED.iterdir()
                          if d.is_dir() and not d.name.startswith("."))

    if args.verify:
        total_broken = 0
        for pid in proteins:
            broken = verify_one(pid)
            for b in broken:
                print(f"  BROKEN: {b}")
            total_broken += len(broken)
        checked = len(proteins)
        print(f"\nVerified {checked} proteins: {total_broken} broken links")
        sys.exit(1 if total_broken > 0 else 0)

    # Populate
    ok_count = 0
    skip_count = 0
    extraction_total = 0
    prep_total = 0

    for pid in proteins:
        r = populate_one(pid, dry_run=args.dry_run)
        if r["ok"]:
            ok_count += 1
            extraction_total += r["extraction"]
            prep_total += r["prep"]
        else:
            skip_count += 1
            print(f"  SKIP: {pid} — {r['error']}")

    tag = "[dry-run] " if args.dry_run else ""
    print(f"\n{tag}Complete: {ok_count}  Skipped: {skip_count}")
    print(f"{tag}Symlinks: {extraction_total} extraction + {prep_total} prep")


if __name__ == "__main__":
    main()
