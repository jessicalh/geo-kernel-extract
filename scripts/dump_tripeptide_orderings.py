#!/usr/bin/env python3
"""Dump one canonical DFT-output row per (tripeptide, frame_type).

Phase 1 step 1 of the typed-tripeptide-topology refactor
(spec/plan/bones/typed-tripeptide-topology-design-2026-05-10.md — RETIRED
2026-05-11 before any of its code landed; this script is an orphan from
that retired plan and is not called from any test or production path).
The output JSON was the source the user reviewed before hand-authoring
data/topology/tripeptide_orderings.toml; the refactor was replaced by
the perception-driven LarsenResidue model and this script never went
into production use.

DSN is read from ~/.nmr_tools.toml [databases].tensorcs15 to keep the
script in line with how Session loads its libpq connection at runtime.
A row's geometry JSONB is returned by psycopg2 as a Python list of
dicts; we forward the per-atom (atom_idx, element, x, y, z) shape and
nothing else, since the DB doesn't carry IUPAC/PDB atom names — those
are what the human is about to assign by inspecting this dump.

The dump is intentionally deterministic: DISTINCT ON (tripeptide,
frame_type) ORDER BY ..., calc_id picks the lowest calc_id for each
combination, so re-running the script against the same DB content
produces a byte-identical dump.
"""

from __future__ import annotations

import argparse
import json
import sys
import tomllib
from pathlib import Path

import psycopg2


CONFIG_PATH = Path.home() / ".nmr_tools.toml"


def load_dsn() -> str:
    with CONFIG_PATH.open("rb") as fh:
        cfg = tomllib.load(fh)
    dbs = cfg.get("databases")
    if not isinstance(dbs, dict) or "tensorcs15" not in dbs:
        raise SystemExit(
            f"{CONFIG_PATH} missing [databases].tensorcs15 DSN string"
        )
    return dbs["tensorcs15"]


def fetch_orderings(dsn: str) -> list[dict]:
    sql = """
        SELECT DISTINCT ON (tripeptide, frame_type)
            tripeptide,
            central_residue,
            frame_type,
            phi, psi,
            chi1, chi2, chi3, chi4,
            n_atoms,
            calc_id,
            geometry
        FROM raw_dft_calculations
        ORDER BY tripeptide, frame_type, calc_id
    """
    rows: list[dict] = []
    with psycopg2.connect(dsn) as conn:
        with conn.cursor() as cur:
            cur.execute(sql)
            for (
                tripeptide,
                central_residue,
                frame_type,
                phi,
                psi,
                c1,
                c2,
                c3,
                c4,
                n_atoms,
                calc_id,
                geometry,
            ) in cur.fetchall():
                if not isinstance(geometry, list):
                    raise SystemExit(
                        f"calc_id={calc_id} geometry JSONB is not a list "
                        f"(type={type(geometry).__name__})"
                    )
                atoms = []
                for a in geometry:
                    atoms.append(
                        {
                            "atom_idx": a["atom_idx"],
                            "element": a["element"],
                            "atomic_number": a.get("atomic_number"),
                            "x": a["x"],
                            "y": a["y"],
                            "z": a["z"],
                        }
                    )
                atoms.sort(key=lambda a: a["atom_idx"])
                rows.append(
                    {
                        "tripeptide": tripeptide,
                        "central_residue": central_residue,
                        "frame_type": frame_type,
                        "calc_id": calc_id,
                        "phi": phi,
                        "psi": psi,
                        "chi1": c1,
                        "chi2": c2,
                        "chi3": c3,
                        "chi4": c4,
                        "n_atoms": n_atoms,
                        "atoms": atoms,
                    }
                )
    return rows


def summarise(rows: list[dict]) -> str:
    lines = [
        f"{r['tripeptide']:>4}  {r['frame_type']:32}  "
        f"n_atoms={r['n_atoms']:>3}  calc_id={r['calc_id']}  "
        f"elements={''.join(a['element'] for a in r['atoms'])}"
        for r in rows
    ]
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--output",
        required=True,
        type=Path,
        help="JSON output path (typically data/topology/tripeptide_orderings_dump.json)",
    )
    ap.add_argument(
        "--summary",
        action="store_true",
        help="Print per-row element-pattern summary to stderr",
    )
    args = ap.parse_args()

    dsn = load_dsn()
    print(f"DSN: {dsn}", file=sys.stderr)
    rows = fetch_orderings(dsn)
    print(f"fetched {len(rows)} (tripeptide, frame_type) combinations", file=sys.stderr)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as fh:
        json.dump(rows, fh, indent=2)
    print(f"wrote {args.output}", file=sys.stderr)

    if args.summary:
        print(summarise(rows), file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
