#!/usr/bin/env python3
"""
Batch MOPAC PM7+MOZYME extraction for all proteins.

Runs MOPAC on each WT and ALA conformation, producing per-atom
charges, bond orders, and scalars as .npy files alongside the
geometric kernel features.

Output per conformation:
    mopac_charges.npy     — (N,) Mulliken charges
    mopac_bond_orders.npy — (B, 3) [atom_i, atom_j, bond_order]
    mopac_scalars.npy     — (N, 3) [charge, s_pop, p_pop]
    mopac_global.npy      — (4,) [heat_of_formation, homo_lumo_gap, 0, 0]

Usage:
    python learn/mopac_extract.py
    python learn/mopac_extract.py --protein P84477
    python learn/mopac_extract.py --resume
    python learn/mopac_extract.py --threads 4
"""

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
CONSOLIDATED = Path("/shared/2026Thesis/consolidated")
FEATURES_BASE = REPO / "learn" / "features" / "FirstExtraction"
MOPAC = Path("/home/jessica/micromamba/envs/mm/bin/mopac")


def xyz_to_mop(xyz_path: Path, charge: int, threads: int = 4) -> str:
    """Convert XYZ file to MOPAC input string."""
    lines = xyz_path.read_text().strip().splitlines()
    natoms = int(lines[0])

    mop = f"PM7 MOZYME 1SCF CHARGE={charge} BONDS MULLIK LET GEO-OK THREADS={threads}\n"
    mop += f"{xyz_path.stem} {natoms} atoms\n"
    mop += "\n"

    for line in lines[2:2+natoms]:
        parts = line.split()
        if len(parts) < 4:
            continue
        elem, x, y, z = parts[0], parts[1], parts[2], parts[3]
        mop += f"  {elem:2s}  {x:>14s} 0  {y:>14s} 0  {z:>14s} 0\n"

    return mop


def parse_mopac_output(out_path: Path, natoms: int):
    """Parse MOPAC .out file for charges, bond orders, scalars."""
    text = out_path.read_text()

    # Mulliken charges and populations
    charges = np.zeros(natoms)
    scalars = np.zeros((natoms, 3))  # charge, s_pop, p_pop

    mulliken = text.find("MULLIKEN POPULATIONS AND CHARGES")
    if mulliken >= 0:
        for line in text[mulliken:].splitlines()[3:3+natoms]:
            m = re.match(r'\s+(\d+)\s+\w+\s+([\d.]+)\s+([-\d.]+)', line)
            if m:
                idx = int(m.group(1)) - 1
                if idx < natoms:
                    charges[idx] = float(m.group(3))

    # Net atomic charges with orbital populations
    net_charges = text.find("NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS")
    if net_charges >= 0:
        for line in text[net_charges:].splitlines()[3:3+natoms+5]:
            parts = line.split()
            if len(parts) >= 6:
                try:
                    idx = int(parts[0]) - 1
                    if 0 <= idx < natoms:
                        scalars[idx, 0] = float(parts[2])  # charge
                        scalars[idx, 1] = float(parts[4])  # s-pop
                        scalars[idx, 2] = float(parts[5])  # p-pop
                except (ValueError, IndexError):
                    pass

    # Bond orders
    bond_orders = []
    bo_start = text.find("BOND ORDERS")
    if bo_start >= 0:
        bo_section = text[bo_start:]
        for line in bo_section.splitlines()[2:]:
            if not line.strip() or "BOND ORDERS" in line:
                continue
            # Format: "  4  C     (3.119)    17  C 1.014     6  C 0.968 ..."
            m = re.match(r'\s+(\d+)\s+\w+\s+\([\d.]+\)(.*)', line)
            if m:
                atom_i = int(m.group(1)) - 1
                rest = m.group(2)
                pairs = re.findall(r'(\d+)\s+\w+\s+([\d.]+)', rest)
                for atom_j_str, bo_str in pairs:
                    atom_j = int(atom_j_str) - 1
                    bo = float(bo_str)
                    if bo > 0.1:  # skip negligible
                        bond_orders.append([atom_i, atom_j, bo])

    if bond_orders:
        bond_orders = np.array(bond_orders)
    else:
        bond_orders = np.zeros((0, 3))

    # Global scalars
    global_scalars = np.zeros(4)
    hof = re.search(r'FINAL HEAT OF FORMATION\s*=\s*([-\d.]+)', text)
    if hof:
        global_scalars[0] = float(hof.group(1))

    gap = re.search(r'HOMO-LUMO GAP\s*=?\s*([-\d.]+)', text)
    if gap:
        global_scalars[1] = float(gap.group(1))

    return charges, scalars, bond_orders, global_scalars


def get_charge(protein_id: str, variant: str) -> int:
    """Get net charge from metadata.json."""
    meta = CONSOLIDATED / protein_id / "metadata.json"
    if meta.exists():
        m = json.loads(meta.read_text())
        key = "wt_charge" if variant == "wt" else "ala_charge"
        return m.get(key, 0)
    return 0


def extract_one(protein_id: str, variant: str, threads: int = 4) -> dict:
    """Run MOPAC on one protein variant."""
    result = {"protein_id": protein_id, "variant": variant,
              "ok": False, "error": "", "elapsed": 0.0}

    suffix = "WT" if variant == "wt" else "ALA"
    xyz_files = list((CONSOLIDATED / protein_id).glob(f"*_{suffix}.xyz"))
    if not xyz_files:
        result["error"] = f"no {suffix} xyz"
        return result

    xyz_path = xyz_files[0]
    natoms = int(xyz_path.read_text().splitlines()[0].strip())
    charge = get_charge(protein_id, variant)

    out_dir = FEATURES_BASE / protein_id / variant
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f"mopac_{protein_id}_{variant}_") as tmpdir:
        mop_path = Path(tmpdir) / "calc.mop"
        mop_path.write_text(xyz_to_mop(xyz_path, charge, threads))

        t0 = time.time()
        try:
            proc = subprocess.run(
                [str(MOPAC), str(mop_path)],
                capture_output=True, text=True, timeout=600
            )
            result["elapsed"] = time.time() - t0
        except subprocess.TimeoutExpired:
            result["error"] = "timeout 600s"
            result["elapsed"] = time.time() - t0
            return result

        out_file = Path(tmpdir) / "calc.out"
        if not out_file.exists():
            result["error"] = "no output file"
            return result

        # Check for normal termination
        out_text = out_file.read_text()
        if "JOB ENDED NORMALLY" not in out_text and "ended normally" not in out_text:
            result["error"] = "abnormal termination"
            return result

        # Check it actually computed (not just geometry error)
        if "HEAT OF FORMATION" not in out_text:
            result["error"] = "no SCF result (geometry rejected?)"
            return result

        charges, scalars, bond_orders, global_scalars = parse_mopac_output(out_file, natoms)

        np.save(out_dir / "mopac_charges.npy", charges)
        np.save(out_dir / "mopac_scalars.npy", scalars)
        np.save(out_dir / "mopac_bond_orders.npy", bond_orders)
        np.save(out_dir / "mopac_global.npy", global_scalars)

        result["ok"] = True
        result["n_bonds"] = len(bond_orders)

    return result


def is_done(protein_id: str, variant: str) -> bool:
    return (FEATURES_BASE / protein_id / variant / "mopac_charges.npy").exists()


def main():
    parser = argparse.ArgumentParser(description="Batch MOPAC extraction")
    parser.add_argument("--protein", type=str, default=None)
    parser.add_argument("--variant", choices=["wt", "ala", "both"], default="both")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if args.protein:
        proteins = [args.protein]
    else:
        proteins = sorted(d.name for d in CONSOLIDATED.iterdir() if d.is_dir())

    variants = ["wt", "ala"] if args.variant == "both" else [args.variant]

    jobs = []
    skipped = 0
    for pid in proteins:
        for v in variants:
            if args.resume and is_done(pid, v):
                skipped += 1
                continue
            jobs.append((pid, v))

    print(f"Proteins: {len(proteins)}  Variants: {len(variants)}  "
          f"Jobs: {len(jobs)}  Skipped: {skipped}  Threads: {args.threads}")

    if not jobs or args.dry_run:
        return

    log_path = FEATURES_BASE / "mopac_log.jsonl"
    ok = 0
    fail = 0
    t_start = time.time()

    for i, (pid, v) in enumerate(jobs):
        r = extract_one(pid, v, threads=args.threads)
        status = "OK" if r["ok"] else "FAIL"
        extra = f"  bonds={r.get('n_bonds', 0)}" if r["ok"] else f"  {r['error']}"
        print(f"[{i+1}/{len(jobs)}] {pid}/{v}: {status}  {r['elapsed']:.1f}s{extra}")
        if r["ok"]:
            ok += 1
        else:
            fail += 1
        with open(log_path, "a") as f:
            f.write(json.dumps(r) + "\n")

    print(f"\nDone: {ok} OK, {fail} failed, {time.time()-t_start:.0f}s total")


if __name__ == "__main__":
    main()
