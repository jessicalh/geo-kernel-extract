#!/usr/bin/env python3
"""
molprobity_validate.py -- run MolProbity geometry validators on a tree of
trajectory-frame PDBs and aggregate the results into per-PDB JSON + a
corpus CSV.

Designed for FramePdbEmitter output (filename convention
{stem}{_decorator?}_f{NNNNNN}_t{ps:.1f}.pdb -- frame_idx and time_ps are
recovered from the filename for the per-PDB record).

Validators (geometry-only; legacy `probe` binary is not required):
  - molprobity.ramalyze            Ramachandran outliers
  - molprobity.rotalyze            Rotamer outliers
  - molprobity.cbetadev            C-beta deviations >= 0.25 A
  - molprobity.cablam              Backbone outliers / disfavored
  - molprobity.mp_validate_bonds   Bond + angle outliers (z-score thresholds)

Clash detection (probe2) is currently NOT wrapped:
  - mmtbx.probe2 default text output has no count surface
  - --json output is unimplemented ("get_results_as_JSON has not been defined")
  - molprobity.clashscore needs the legacy probe binary (not installed)
  - molprobity.clashscore2 + mmtbx.nonbonded_overlaps are broken in this env
  Listed as TODO; geometry validators alone surface most "topology sins".

Each validator's SUMMARY-prefixed lines are captured verbatim alongside the
parsed numerical fields so a reviewer can see what we extracted vs what we
left on the floor.

Usage:
  tools/molprobity_validate.py \
      --pdb-dir molprobity_runs/1P9J_5801/pdbs \
      --output-dir molprobity_runs/_validation/1P9J_5801 \
      [--n-frames 10] \
      [--workers 8]
"""

import argparse
import csv
import json
import multiprocessing as mp
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ENV_BIN = Path("/home/jessica/micromamba/envs/molprobity/bin")
RAMALYZE = str(ENV_BIN / "molprobity.ramalyze")
ROTALYZE = str(ENV_BIN / "molprobity.rotalyze")
CBETADEV = str(ENV_BIN / "molprobity.cbetadev")
CABLAM   = str(ENV_BIN / "molprobity.cablam")
MPBONDS  = str(ENV_BIN / "molprobity.mp_validate_bonds")

TIMEOUT_SECONDS = 60


# ---------------------------------------------------------------------------
# Subprocess driver (timeouts surface as captured records, not exceptions)
# ---------------------------------------------------------------------------

def run_validator(cmd):
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=TIMEOUT_SECONDS
        )
        return {
            "rc": proc.returncode,
            "stdout": proc.stdout,
            "stderr": proc.stderr,
        }
    except subprocess.TimeoutExpired:
        return {"rc": -1, "stdout": "", "stderr": "TIMEOUT"}
    except Exception as exc:
        return {"rc": -2, "stdout": "", "stderr": f"EXC: {exc}"}


def collect_summary_lines(stdout):
    return [ln for ln in stdout.splitlines() if ln.lstrip().startswith("SUMMARY:")]


# ---------------------------------------------------------------------------
# Parsers (verbatim SUMMARY lines preserved; parsed fields are best-effort)
# ---------------------------------------------------------------------------

def parse_ramalyze(stdout):
    fav = allowed = outlier = total = -1
    pct_out = pct_fav = float("nan")
    for line in stdout.splitlines():
        m = re.search(
            r"SUMMARY:\s+(\d+)\s+Favored,\s+(\d+)\s+Allowed,\s+(\d+)\s+Outlier"
            r"\s+out of\s+(\d+)",
            line,
        )
        if m:
            fav, allowed, outlier, total = (int(x) for x in m.groups())
        m = re.search(r"SUMMARY:\s+([\d.]+)%\s+outliers", line)
        if m:
            pct_out = float(m.group(1))
        m = re.search(r"SUMMARY:\s+([\d.]+)%\s+favored", line)
        if m:
            pct_fav = float(m.group(1))
    return {
        "n_favored":  fav,
        "n_allowed":  allowed,
        "n_outlier":  outlier,
        "n_total":    total,
        "pct_outlier": pct_out,
        "pct_favored": pct_fav,
        "summary_lines": collect_summary_lines(stdout),
    }


def parse_rotalyze(stdout):
    """rotalyze emits only one SUMMARY line (% outliers).
    Per-residue lines have variable column counts (1-4 chi angles); count
    them to recover the totals.
    """
    pct_out = float("nan")
    for line in stdout.splitlines():
        m = re.search(r"SUMMARY:\s+([\d.]+)%\s+outliers", line)
        if m:
            pct_out = float(m.group(1))
    n_outlier = n_allowed = n_favored = n_total = 0
    for line in stdout.splitlines():
        if not re.match(r"^\s*[A-Z\s]\s+\d+\s+[A-Z]{3}:", line):
            continue
        n_total += 1
        if "OUTLIER" in line:
            n_outlier += 1
        elif "Allowed" in line:
            n_allowed += 1
        elif "Favored" in line:
            n_favored += 1
    return {
        "n_favored":  n_favored,
        "n_allowed":  n_allowed,
        "n_outlier":  n_outlier,
        "n_total":    n_total,
        "pct_outlier": pct_out,
        "summary_lines": collect_summary_lines(stdout),
    }


def parse_cbetadev(stdout):
    n_dev = -1
    for line in stdout.splitlines():
        m = re.search(
            r"SUMMARY:\s+(\d+)\s+C-beta\s+deviations\s+>=\s+0\.25\s+Angstrom",
            line,
        )
        if m:
            n_dev = int(m.group(1))
    max_dev = 0.0
    for line in stdout.splitlines():
        m = re.match(
            r"^\S+\s*:\s*\S*\s*:\s*\S+\s*:\s*\S+\s*:\s*\d+\s*:\s*([\d.]+)\s*:",
            line,
        )
        if m:
            v = float(m.group(1))
            if v > max_dev:
                max_dev = v
    return {
        "n_deviations_ge_0p25A": n_dev,
        "max_deviation_a":       max_dev,
        "summary_lines":         collect_summary_lines(stdout),
    }


def parse_cablam(stdout):
    n_full = n_caonly = -1
    n_disfavored = pct_disfavored = float("nan")
    n_outlier = pct_outlier = float("nan")
    n_severe  = pct_severe   = float("nan")
    for line in stdout.splitlines():
        m = re.search(
            r"SUMMARY:\s+CaBLAM\s+found\s+(\d+)\s+full\s+protein\s+residues"
            r"\s+and\s+(\d+)\s+CA-only\s+residues",
            line,
        )
        if m:
            n_full, n_caonly = int(m.group(1)), int(m.group(2))
        m = re.search(
            r"SUMMARY:\s+(\d+)\s+residues\s+\(([\d.]+)%\)\s+have\s+disfavored",
            line,
        )
        if m:
            n_disfavored, pct_disfavored = int(m.group(1)), float(m.group(2))
        m = re.search(
            r"SUMMARY:\s+(\d+)\s+residues\s+\(([\d.]+)%\)\s+have\s+outlier",
            line,
        )
        if m:
            n_outlier, pct_outlier = int(m.group(1)), float(m.group(2))
        m = re.search(
            r"SUMMARY:\s+(\d+)\s+residues\s+\(([\d.]+)%\)\s+have\s+severe",
            line,
        )
        if m:
            n_severe, pct_severe = int(m.group(1)), float(m.group(2))
    return {
        "n_full_residues":   n_full,
        "n_ca_only":         n_caonly,
        "n_disfavored":      n_disfavored,
        "pct_disfavored":    pct_disfavored,
        "n_outlier":         n_outlier,
        "pct_outlier":       pct_outlier,
        "n_severe_geometry": n_severe,
        "pct_severe":        pct_severe,
        "summary_lines":     collect_summary_lines(stdout),
    }


def _bond_body_rows(stdout):
    """Yield bond outlier body rows (between 'Bond lengths' and 'Bond
    angles' headers).  Each row's tail is `... atom1_name atom2_name
    sigma`.  Returns whitespace-tokenised lines as a list of token
    lists."""
    in_bond_section = False
    rows = []
    for line in stdout.splitlines():
        if "----------Bond lengths----------" in line:
            in_bond_section = True
            continue
        if "----------Bond angles----------" in line:
            in_bond_section = False
            continue
        if not in_bond_section:
            continue
        toks = line.split()
        # Bond rows have at least 11 whitespace tokens:
        #   chain res# restype atom1   chain res# restype atom2   atom1 atom2 sigma
        # Header rows ("residue   atom 1  atom 2  sigmas") have <= 5.
        if len(toks) < 11:
            continue
        # Last token must be a numeric sigma to be a real outlier row.
        try:
            float(toks[-1])
        except ValueError:
            continue
        rows.append(toks)
    return rows


def parse_mp_bonds(stdout):
    n_bond_out = n_bond_total = -1
    n_ang_out  = n_ang_total  = -1
    for line in stdout.splitlines():
        m = re.search(r"^\s*(\d+)/(\d+)\s+bond\s+outliers\s+present", line)
        if m:
            n_bond_out, n_bond_total = int(m.group(1)), int(m.group(2))
        m = re.search(r"^\s*(\d+)/(\d+)\s+angle\s+outliers\s+present", line)
        if m:
            n_ang_out, n_ang_total = int(m.group(1)), int(m.group(2))

    # Heavy-vs-H split: AMBER ff14SB MD with LINCS-constrained H-X bonds
    # produces a uniform mean-offset sigma for H-X (~6 for C-H/O-H, ~7.5
    # for N-H) against MolProbity's X-ray riding-H reference.  That is a
    # known systematic, not a geometry sin.  The actionable signal is the
    # heavy-atom outlier count: bonds where neither endpoint is hydrogen.
    n_bond_out_heavy = 0
    n_bond_out_h     = 0
    max_sigma_heavy  = 0.0
    max_sigma_h      = 0.0
    for toks in _bond_body_rows(stdout):
        atom1 = toks[-3]
        atom2 = toks[-2]
        sigma = abs(float(toks[-1]))
        if atom1.startswith("H") or atom2.startswith("H"):
            n_bond_out_h += 1
            if sigma > max_sigma_h:
                max_sigma_h = sigma
        else:
            n_bond_out_heavy += 1
            if sigma > max_sigma_heavy:
                max_sigma_heavy = sigma

    return {
        "n_bond_outliers":       n_bond_out,
        "n_bond_total":          n_bond_total,
        "n_bond_outliers_heavy": n_bond_out_heavy,
        "n_bond_outliers_h":     n_bond_out_h,
        "max_sigma_heavy":       max_sigma_heavy,
        "max_sigma_h":           max_sigma_h,
        "n_angle_outliers":      n_ang_out,
        "n_angle_total":         n_ang_total,
        "summary_lines": [
            ln for ln in stdout.splitlines() if "outliers present" in ln
        ],
    }


# ---------------------------------------------------------------------------
# Per-PDB driver
# ---------------------------------------------------------------------------

FILENAME_RE = re.compile(r"^.*_f(\d+)_t([-\d.]+)\.pdb$")


def parse_pdb_filename(name):
    m = FILENAME_RE.match(name)
    if not m:
        return None, None
    return int(m.group(1)), float(m.group(2))


def validate_pdb(pdb_path):
    """Run all five validators on one PDB. Returns one dict per PDB."""
    frame_idx, time_ps = parse_pdb_filename(pdb_path.name)
    record = {
        "pdb": pdb_path.name,
        "frame_idx": frame_idx,
        "time_ps":   time_ps,
    }

    raw_rama = run_validator([RAMALYZE, str(pdb_path)])
    record["ramalyze"] = {"rc": raw_rama["rc"], **parse_ramalyze(raw_rama["stdout"])}

    raw_rota = run_validator([ROTALYZE, str(pdb_path)])
    record["rotalyze"] = {"rc": raw_rota["rc"], **parse_rotalyze(raw_rota["stdout"])}

    raw_cb = run_validator([CBETADEV, str(pdb_path)])
    record["cbetadev"] = {"rc": raw_cb["rc"], **parse_cbetadev(raw_cb["stdout"])}

    raw_cabl = run_validator([CABLAM, str(pdb_path)])
    record["cablam"] = {"rc": raw_cabl["rc"], **parse_cablam(raw_cabl["stdout"])}

    raw_bonds = run_validator([MPBONDS, str(pdb_path)])
    record["mp_bonds"] = {"rc": raw_bonds["rc"], **parse_mp_bonds(raw_bonds["stdout"])}

    return record


# ---------------------------------------------------------------------------
# Aggregate + write
# ---------------------------------------------------------------------------

CSV_COLUMNS = [
    "pdb",
    "frame_idx",
    "time_ps",
    "rama_pct_favored",
    "rama_pct_outlier",
    "rama_n_outlier",
    "rama_n_total",
    "rota_pct_outlier",
    "rota_n_outlier",
    "rota_n_total",
    "cbeta_n_dev_ge_0p25A",
    "cbeta_max_dev_a",
    "cablam_n_disfavored",
    "cablam_pct_disfavored",
    "cablam_n_outlier",
    "cablam_pct_outlier",
    "cablam_n_severe",
    "n_bond_outliers",
    "n_bond_outliers_heavy",
    "n_bond_outliers_h",
    "n_bond_total",
    "max_sigma_heavy",
    "max_sigma_h",
    "n_angle_outliers",
    "n_angle_total",
]


def record_to_row(rec):
    r = rec.get("ramalyze",  {})
    rt = rec.get("rotalyze",  {})
    cb = rec.get("cbetadev", {})
    cl = rec.get("cablam",   {})
    mb = rec.get("mp_bonds", {})
    return [
        rec.get("pdb", ""),
        rec.get("frame_idx", ""),
        rec.get("time_ps", ""),
        r.get("pct_favored", ""),
        r.get("pct_outlier", ""),
        r.get("n_outlier", ""),
        r.get("n_total", ""),
        rt.get("pct_outlier", ""),
        rt.get("n_outlier", ""),
        rt.get("n_total", ""),
        cb.get("n_deviations_ge_0p25A", ""),
        cb.get("max_deviation_a", ""),
        cl.get("n_disfavored", ""),
        cl.get("pct_disfavored", ""),
        cl.get("n_outlier", ""),
        cl.get("pct_outlier", ""),
        cl.get("n_severe_geometry", ""),
        mb.get("n_bond_outliers", ""),
        mb.get("n_bond_outliers_heavy", ""),
        mb.get("n_bond_outliers_h", ""),
        mb.get("n_bond_total", ""),
        mb.get("max_sigma_heavy", ""),
        mb.get("max_sigma_h", ""),
        mb.get("n_angle_outliers", ""),
        mb.get("n_angle_total", ""),
    ]


def write_outputs(records, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)

    # Per-PDB JSON
    per_dir = output_dir / "per_pdb"
    per_dir.mkdir(exist_ok=True)
    for rec in records:
        stem = Path(rec["pdb"]).stem
        with open(per_dir / f"{stem}.json", "w") as fh:
            json.dump(rec, fh, indent=2)

    # Aggregate JSON
    with open(output_dir / "validation_aggregate.json", "w") as fh:
        json.dump(records, fh, indent=2)

    # Aggregate CSV (frame-ordered when frame_idx is present)
    sorted_records = sorted(
        records,
        key=lambda r: (r.get("frame_idx") if isinstance(r.get("frame_idx"), int) else 0)
    )
    with open(output_dir / "validation_aggregate.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_COLUMNS)
        for rec in sorted_records:
            w.writerow(record_to_row(rec))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb-dir", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--n-frames", type=int, default=0,
                    help="Limit to first N PDBs by sorted filename (0 = all)")
    ap.add_argument("--workers", type=int, default=1,
                    help="Parallel worker count (1 = sequential, "
                         "useful for the smoke test)")
    args = ap.parse_args()

    if not args.pdb_dir.is_dir():
        print(f"ERROR: --pdb-dir {args.pdb_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    pdbs = sorted(args.pdb_dir.glob("*.pdb"))
    if args.n_frames > 0:
        pdbs = pdbs[:args.n_frames]
    if not pdbs:
        print(f"ERROR: no PDBs found under {args.pdb_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"[validate] {len(pdbs)} PDBs, workers={args.workers}", file=sys.stderr)

    records = []
    if args.workers <= 1:
        for i, pdb in enumerate(pdbs, 1):
            rec = validate_pdb(pdb)
            records.append(rec)
            print(
                f"  [{i:>4}/{len(pdbs)}] {pdb.name}: "
                f"rama_out={rec['ramalyze'].get('pct_outlier')}% "
                f"bonds_out={rec['mp_bonds'].get('n_bond_outliers')}/"
                f"{rec['mp_bonds'].get('n_bond_total')} "
                f"cbeta={rec['cbetadev'].get('n_deviations_ge_0p25A')}",
                file=sys.stderr,
            )
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as exe:
            future_to_pdb = {exe.submit(validate_pdb, p): p for p in pdbs}
            done = 0
            for fut in as_completed(future_to_pdb):
                rec = fut.result()
                records.append(rec)
                done += 1
                if done % 50 == 0 or done == len(pdbs):
                    print(f"  [{done}/{len(pdbs)}] done", file=sys.stderr)

    write_outputs(records, args.output_dir)
    print(
        f"[validate] wrote per-PDB JSON, "
        f"validation_aggregate.json, validation_aggregate.csv "
        f"to {args.output_dir}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
