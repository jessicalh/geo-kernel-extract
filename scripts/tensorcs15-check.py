#!/usr/bin/env python3
"""Read-only tensorcs15 restore/integrity checker.

This validates the live PostgreSQL database against the tracked frozen-corpus
manifest. It does not print the DSN and does not modify the database.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11
    print("FAIL Python 3.11+ tomllib is required", file=sys.stderr)
    sys.exit(2)

try:
    import psycopg2
    from psycopg2 import sql
except ModuleNotFoundError:
    print("FAIL psycopg2 is required for tensorcs15 checks", file=sys.stderr)
    sys.exit(2)


def load_toml(path: Path) -> dict:
    try:
        with path.open("rb") as fh:
            return tomllib.load(fh)
    except Exception as exc:
        raise SystemExit(f"FAIL cannot read TOML {path}: {exc}") from exc


def default_manifest_path() -> Path:
    script_dir = Path(__file__).resolve().parent
    root = script_dir.parent
    candidates = [
        root / "config" / "tensorcs15.manifest.toml",
        root / "share" / "nmr-shielding" / "tensorcs15.manifest.toml",
        Path("/etc/nmr-shielding/tensorcs15.manifest.toml"),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0]


def default_config_path() -> Path:
    env_path = os.environ.get("NMR_TOOLS_TOML")
    if env_path:
        return Path(env_path)
    return Path.home() / ".nmr_tools.toml"


def load_dsn(config_path: Path) -> str:
    cfg = load_toml(config_path) if config_path.is_file() else {}
    dbs = cfg.get("databases", {})
    dsn = dbs.get("tensorcs15", "") if isinstance(dbs, dict) else ""
    if dsn:
        return str(dsn)
    return os.environ.get("NMR_TENSORCS15_DSN", "")


class Reporter:
    def __init__(self, quiet: bool = False) -> None:
        self.quiet = quiet
        self.failures = 0

    def ok(self, label: str, value: object) -> None:
        if not self.quiet:
            print(f"OK    {label}: {value}")

    def fail(self, label: str, actual: object, expected: object) -> None:
        self.failures += 1
        print(
            f"FAIL  {label}: actual={actual!r} expected={expected!r}",
            file=sys.stderr,
        )

    def check_equal(self, label: str, actual: object, expected: object) -> None:
        if actual == expected:
            self.ok(label, actual)
        else:
            self.fail(label, actual, expected)


def fetch_one(cur, query, params=()):
    cur.execute(query, params)
    row = cur.fetchone()
    return row[0] if row else None


def fetch_dict(cur, query, params=()):
    cur.execute(query, params)
    return {str(k): int(v) for k, v in cur.fetchall()}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, default=default_config_path(),
                    help="Runtime TOML path containing [databases].tensorcs15")
    ap.add_argument("--manifest", type=Path, default=default_manifest_path(),
                    help="tensorcs15 manifest TOML")
    ap.add_argument("--quiet", action="store_true",
                    help="Print only failures and summary")
    args = ap.parse_args()

    manifest = load_toml(args.manifest)
    dsn = load_dsn(args.config)
    if not dsn:
        print(
            "FAIL tensorcs15 DSN not configured; set [databases].tensorcs15 "
            "or NMR_TENSORCS15_DSN",
            file=sys.stderr,
        )
        return 2

    table = str(manifest.get("table", "raw_dft_calculations"))
    table_id = sql.Identifier(table)
    rep = Reporter(quiet=args.quiet)

    if not args.quiet:
        print(f"INFO  manifest: {args.manifest}")
        print(f"INFO  manifest_id: {manifest.get('manifest_id', '<not set>')}")
        print("INFO  DSN configured (value not printed)")

    try:
        with psycopg2.connect(dsn) as conn:
            conn.set_session(readonly=True, autocommit=True)
            with conn.cursor() as cur:
                cur.execute(sql.SQL("SELECT to_regclass(%s)"), (table,))
                regclass = cur.fetchone()[0]
                rep.check_equal("table exists", regclass is not None, True)
                if regclass is None:
                    print(
                        f"SUMMARY failed={rep.failures}",
                        file=sys.stderr,
                    )
                    return 1

                cur.execute(
                    """
                    SELECT column_name
                      FROM information_schema.columns
                     WHERE table_name = %s
                    """,
                    (table,),
                )
                actual_columns = {row[0] for row in cur.fetchall()}
                expected_columns = set(manifest.get("required_columns", []))
                rep.check_equal(
                    "required columns",
                    sorted(actual_columns & expected_columns),
                    sorted(expected_columns),
                )

                rep.check_equal(
                    "total rows",
                    fetch_one(
                        cur,
                        sql.SQL("SELECT count(*) FROM {}").format(table_id),
                    ),
                    int(manifest["total_rows"]),
                )
                rep.check_equal(
                    "distinct tripeptides",
                    fetch_one(
                        cur,
                        sql.SQL(
                            "SELECT count(DISTINCT tripeptide) FROM {}"
                        ).format(table_id),
                    ),
                    int(manifest["distinct_tripeptides"]),
                )
                rep.check_equal(
                    "distinct frame types",
                    fetch_one(
                        cur,
                        sql.SQL(
                            "SELECT count(DISTINCT frame_type) FROM {}"
                        ).format(table_id),
                    ),
                    int(manifest["distinct_frame_types"]),
                )
                cur.execute(
                    sql.SQL("SELECT min(calc_id), max(calc_id) FROM {}").format(
                        table_id
                    )
                )
                min_calc_id, max_calc_id = cur.fetchone()
                rep.check_equal("min calc_id", min_calc_id, int(manifest["min_calc_id"]))
                rep.check_equal("max calc_id", max_calc_id, int(manifest["max_calc_id"]))
                rep.check_equal(
                    "null payloads",
                    fetch_one(
                        cur,
                        sql.SQL(
                            "SELECT count(*) FROM {} "
                            "WHERE geometry IS NULL OR tensors IS NULL"
                        ).format(table_id),
                    ),
                    int(manifest["null_payloads"]),
                )

                actual_frame_counts = fetch_dict(
                    cur,
                    sql.SQL(
                        "SELECT frame_type, count(*) FROM {} "
                        "GROUP BY frame_type ORDER BY frame_type"
                    ).format(table_id),
                )
                expected_frame_counts = {
                    str(k): int(v)
                    for k, v in manifest.get("frame_type_counts", {}).items()
                }
                rep.check_equal(
                    "frame type keys",
                    sorted(actual_frame_counts),
                    sorted(expected_frame_counts),
                )
                for key, expected in expected_frame_counts.items():
                    rep.check_equal(
                        f"frame type count {key}",
                        actual_frame_counts.get(key),
                        expected,
                    )

                actual_tripeptide_counts = fetch_dict(
                    cur,
                    sql.SQL(
                        "SELECT tripeptide, count(*) FROM {} "
                        "GROUP BY tripeptide ORDER BY tripeptide"
                    ).format(table_id),
                )
                expected_tripeptide_counts = {
                    str(k): int(v)
                    for k, v in manifest.get("tripeptide_counts", {}).items()
                }
                rep.check_equal(
                    "tripeptide keys",
                    sorted(actual_tripeptide_counts),
                    sorted(expected_tripeptide_counts),
                )
                for key, expected in expected_tripeptide_counts.items():
                    rep.check_equal(
                        f"tripeptide count {key}",
                        actual_tripeptide_counts.get(key),
                        expected,
                    )

                grid = manifest.get("grid_minus_180_minus_180", {})
                rep.check_equal(
                    "grid -180/-180 distinct tripeptides",
                    fetch_one(
                        cur,
                        sql.SQL(
                            "SELECT count(DISTINCT tripeptide) FROM {} "
                            "WHERE phi = -180 AND psi = -180"
                        ).format(table_id),
                    ),
                    int(grid["distinct_tripeptides"]),
                )
                grid_counts = fetch_dict(
                    cur,
                    sql.SQL(
                        "SELECT frame_type, count(*) FROM {} "
                        "WHERE phi = -180 AND psi = -180 "
                        "GROUP BY frame_type ORDER BY frame_type"
                    ).format(table_id),
                )
                for key, expected in grid.get("frame_type_counts", {}).items():
                    rep.check_equal(
                        f"grid -180/-180 frame rows {key}",
                        grid_counts.get(str(key)),
                        int(expected),
                    )
                grid_distinct = fetch_dict(
                    cur,
                    sql.SQL(
                        "SELECT frame_type, count(DISTINCT tripeptide) FROM {} "
                        "WHERE phi = -180 AND psi = -180 "
                        "GROUP BY frame_type ORDER BY frame_type"
                    ).format(table_id),
                )
                for key, expected in grid.get(
                    "frame_type_distinct_tripeptides", {}
                ).items():
                    rep.check_equal(
                        f"grid -180/-180 frame tripeptides {key}",
                        grid_distinct.get(str(key)),
                        int(expected),
                    )

                sentinel = manifest.get("sentinel", {})
                cur.execute(
                    sql.SQL(
                        "SELECT calc_id, tripeptide, phi, psi, n_atoms, "
                        "frame_type, jsonb_typeof(geometry), "
                        "jsonb_array_length(geometry), jsonb_typeof(tensors), "
                        "jsonb_array_length(tensors) "
                        "FROM {} "
                        "WHERE tripeptide = %s AND phi = %s AND psi = %s "
                        "ORDER BY calc_id ASC LIMIT 1"
                    ).format(table_id),
                    (
                        sentinel["tripeptide"],
                        int(sentinel["phi"]),
                        int(sentinel["psi"]),
                    ),
                )
                row = cur.fetchone()
                rep.check_equal("sentinel row exists", row is not None, True)
                if row is not None:
                    labels = [
                        "calc_id",
                        "tripeptide",
                        "phi",
                        "psi",
                        "n_atoms",
                        "frame_type",
                        "geometry_type",
                        "geometry_len",
                        "tensors_type",
                        "tensors_len",
                    ]
                    expected_values = [
                        int(sentinel["calc_id"]),
                        sentinel["tripeptide"],
                        int(sentinel["phi"]),
                        int(sentinel["psi"]),
                        int(sentinel["n_atoms"]),
                        sentinel["frame_type"],
                        sentinel["geometry_type"],
                        int(sentinel["geometry_len"]),
                        sentinel["tensors_type"],
                        int(sentinel["tensors_len"]),
                    ]
                    for label, actual, expected in zip(labels, row, expected_values):
                        rep.check_equal(f"sentinel {label}", actual, expected)
    except Exception as exc:
        print(f"FAIL tensorcs15 check error: {exc}", file=sys.stderr)
        return 1

    if rep.failures:
        print(f"SUMMARY failed={rep.failures}", file=sys.stderr)
        return 1
    print("SUMMARY failed=0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
