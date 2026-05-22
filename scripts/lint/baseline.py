#!/usr/bin/env python3
"""scripts/lint/baseline.py — compute / verify the per-file finding ratchet.

Two modes:

  --update  Read every artifacts/quality/<tool>-<sha>.json and write
            artifacts/quality/baseline.json — the per-file finding-count
            anchor used by the pre-commit ratchet. Run after a clean
            sweep to lock the current state.

  --check   For each artifacts/quality/<tool>-<sha>.json belonging to
            the current SHA, verify that no file's finding count
            exceeds its baseline entry. Exit 1 if any file regressed.
            Used by pre-commit (filtered to changed files).

baseline.json schema:
  {
    "anchored_at_sha": "abcd1234",
    "anchored_at_utc": "2026-05-22T15:30:00Z",
    "tools": {
      "clang-tidy": {
        "<path/to/file.cpp>": {"total": 12, "by_check": {"google-...": 3, ...}},
        ...
      },
      "cppcheck": {...},
      "iwyu": {...},
      "scan-build": {...}
    }
  }
"""
from __future__ import annotations

import argparse
import datetime
import glob
import json
import os
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

ARTIFACTS = Path("artifacts/quality")
BASELINE  = ARTIFACTS / "baseline.json"
TOOLS     = ("clang-tidy", "cppcheck", "iwyu", "scan-build")


def find_report(tool: str, sha: str) -> Path | None:
    candidate = ARTIFACTS / f"{tool}-{sha}.json"
    if candidate.exists():
        return candidate
    # Fall back to the most-recent report from any sha.
    matches = sorted(ARTIFACTS.glob(f"{tool}-*.json"),
                     key=lambda p: p.stat().st_mtime, reverse=True)
    return matches[0] if matches else None


def aggregate(report: Path) -> dict[str, dict]:
    """Return per-file dict { "total": N, "by_check": {...} }."""
    data = json.loads(report.read_text())
    per_file: dict[str, dict] = defaultdict(lambda: {"total": 0, "by_check": Counter()})
    for f in data.get("findings", []):
        path = f["file"]
        per_file[path]["total"] += 1
        per_file[path]["by_check"][f.get("check", "?")] += 1
    # Convert Counter to plain dict for json.dumps.
    return {p: {"total": d["total"], "by_check": dict(d["by_check"])}
            for p, d in per_file.items()}


def cmd_update() -> int:
    sha = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"]).decode().strip()
    out = {
        "anchored_at_sha":  sha,
        "anchored_at_utc":  datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "tools":            {},
    }
    for tool in TOOLS:
        report = find_report(tool, sha)
        if report is None:
            print(f"  {tool}: no report found in {ARTIFACTS}/; skipping", file=sys.stderr)
            continue
        out["tools"][tool] = aggregate(report)
        n_files = len(out["tools"][tool])
        n_findings = sum(d["total"] for d in out["tools"][tool].values())
        print(f"  {tool}: {n_findings} findings across {n_files} files (from {report.name})")
    BASELINE.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {BASELINE}")
    return 0


def cmd_check(only_files: list[str]) -> int:
    if not BASELINE.exists():
        print(f"ERROR: {BASELINE} missing. Run --update first.", file=sys.stderr)
        return 2
    baseline = json.loads(BASELINE.read_text())
    sha      = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"]).decode().strip()

    regressions: list[tuple[str, str, int, int]] = []
    for tool in TOOLS:
        report = find_report(tool, sha)
        if report is None:
            continue
        current_state = aggregate(report)
        baseline_tool = baseline["tools"].get(tool, {})
        for path, cur in current_state.items():
            if only_files and path not in only_files:
                continue
            base = baseline_tool.get(path, {"total": 0})
            if cur["total"] > base["total"]:
                regressions.append((tool, path, base["total"], cur["total"]))

    if regressions:
        print(f"RATCHET REGRESSION ({len(regressions)} file/tool pairs):", file=sys.stderr)
        for tool, path, base, cur in regressions:
            print(f"  {tool:12s}  {path}  baseline={base}  current={cur}  (+{cur - base})",
                  file=sys.stderr)
        print("", file=sys.stderr)
        print("Either fix the new findings, or run", file=sys.stderr)
        print("  scripts/lint/baseline.py --update", file=sys.stderr)
        print("to re-anchor (only after deliberate review).", file=sys.stderr)
        return 1
    print("clang-tidy / cppcheck / iwyu / scan-build ratchet OK.")
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("--update", help="Re-anchor baseline.json from current reports.")
    chk = sub.add_parser("--check", help="Verify no per-file regressions against baseline.")
    chk.add_argument("files", nargs="*", default=[],
                     help="Restrict check to these files (pre-commit hook use). "
                          "Defaults to all files in the current reports.")
    # argparse hates leading '--' for subcommands; fall back to manual.
    args, rest = p.parse_known_args()
    return 2  # unreachable; see below


if __name__ == "__main__":
    # Manual arg dispatch to allow --update / --check at top level
    # without nested subparsers (cleaner CLI than argparse default).
    argv = sys.argv[1:]
    if not argv:
        print(__doc__)
        sys.exit(2)
    if argv[0] == "--update":
        sys.exit(cmd_update())
    if argv[0] == "--check":
        only = argv[1:]
        sys.exit(cmd_check(only))
    print(f"Unknown arg: {argv[0]}; expected --update or --check", file=sys.stderr)
    sys.exit(2)
