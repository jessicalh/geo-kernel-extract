#!/usr/bin/env python3
"""Rediscover run-root manager.

This is a script-level wrapper around the existing rediscover extractor and
Python consumers. It standardizes run directories, records a small run index,
and drops only superseded emitted substrate files from managed run dirs.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


DEFAULT_ROOT = Path(os.environ.get("REDISCOVER_RUN_ROOT", "/tmp/rediscover-runs"))
DEFAULT_KEEP_SUBSTRATE = int(os.environ.get("REDISCOVER_KEEP_SUBSTRATE", "1"))
DEFAULT_ACTIVE_MINUTES = int(os.environ.get("REDISCOVER_ACTIVE_MINUTES", "720"))
DEFAULT_EXTRACT_BIN = os.environ.get(
    "REDISCOVER_EXTRACT_BIN", "build/linux-gcc/h5reader_extract"
)

RUN_META = ".rediscover-run.json"
ROOT_MANIFEST = "manifest.json"
ANALYSIS_DIR = "analysis"
RUN_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
SHARED_ROOT = Path("/shared").resolve()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def parse_time(value: object) -> float | None:
    if not isinstance(value, str) or not value:
        return None
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00")).timestamp()
    except ValueError:
        return None


def clean_script_args(args: list[str]) -> list[str]:
    return args[1:] if args and args[0] == "--" else args


def require_run_name(name: str) -> str:
    if not RUN_NAME_RE.match(name) or name in {".", ".."}:
        raise SystemExit(
            f"invalid run name {name!r}; use letters, digits, '.', '_' and '-' only"
        )
    return name


def resolve_root(path: Path | None) -> Path:
    return (path or DEFAULT_ROOT).expanduser().resolve()


def run_dir(root: Path, name: str) -> Path:
    return root / require_run_name(name)


def load_json(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)
    except FileNotFoundError:
        return {}
    except json.JSONDecodeError as exc:
        raise SystemExit(f"cannot parse {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise SystemExit(f"{path} does not contain a JSON object")
    return data


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write("\n")
    tmp.replace(path)


def load_meta(path: Path) -> dict[str, Any]:
    return load_json(path / RUN_META)


def save_meta(path: Path, meta: dict[str, Any]) -> None:
    meta["updated_at"] = utc_now()
    write_json(path / RUN_META, meta)


def relpath(path: Path, base: Path) -> str:
    return path.relative_to(base).as_posix()


def producer_manifest_substrate(run_path: Path) -> set[str]:
    manifest = load_json(run_path / "manifest.json")
    out: set[str] = set()
    relationships = manifest.get("relationships", [])
    if not isinstance(relationships, list):
        return out
    for relationship in relationships:
        if not isinstance(relationship, dict):
            continue
        outputs = relationship.get("outputs", {})
        if not isinstance(outputs, dict):
            continue
        for key in ("sources_csv", "aggregated_csv"):
            value = outputs.get(key)
            if isinstance(value, str) and value:
                out.add(value)
        npys = outputs.get("array_payload_npys", [])
        if isinstance(npys, list):
            for value in npys:
                if isinstance(value, str) and value:
                    out.add(value)
    return out


def recorded_substrate(run_path: Path) -> set[str]:
    meta = load_meta(run_path)
    files = meta.get("substrate_files", [])
    if not isinstance(files, list):
        return set()
    return {value for value in files if isinstance(value, str) and value}


def pattern_substrate(run_path: Path) -> set[str]:
    out: set[str] = set()
    for child in run_path.iterdir() if run_path.exists() else []:
        if not child.is_file():
            continue
        if child.name in {RUN_META, "manifest.json"}:
            continue
        if child.suffix in {".csv", ".npy", ".npz"} or child.name.endswith(".csv.gz"):
            out.add(child.name)
    return out


def discover_substrate(run_path: Path, *, include_patterns: bool = False) -> list[str]:
    names = recorded_substrate(run_path) | producer_manifest_substrate(run_path)
    if include_patterns and not names:
        names |= pattern_substrate(run_path)
    existing = [name for name in sorted(names) if (run_path / name).is_file()]
    return existing


def max_mtime(path: Path) -> float:
    if not path.exists():
        return 0.0
    newest = path.stat().st_mtime
    for root, dirs, files in os.walk(path):
        for name in dirs + files:
            try:
                newest = max(newest, (Path(root) / name).stat().st_mtime)
            except FileNotFoundError:
                pass
    return newest


def size_of(paths: list[Path]) -> int:
    total = 0
    for path in paths:
        try:
            total += path.stat().st_size
        except FileNotFoundError:
            pass
    return total


def dir_size(path: Path, *, exclude_prefix: Path | None = None) -> int:
    if not path.exists():
        return 0
    total = 0
    exclude_resolved = exclude_prefix.resolve() if exclude_prefix and exclude_prefix.exists() else None
    for root, dirs, files in os.walk(path):
        root_path = Path(root)
        if exclude_resolved is not None:
            dirs[:] = [
                d for d in dirs
                if not (root_path / d).resolve().is_relative_to(exclude_resolved)
            ]
        for name in files:
            try:
                total += (root_path / name).stat().st_size
            except FileNotFoundError:
                pass
    return total


def human_size(n: int) -> str:
    units = ("B", "KB", "MB", "GB", "TB")
    value = float(n)
    for unit in units:
        if value < 1024.0 or unit == units[-1]:
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024.0
    return f"{n} B"


def is_managed_run(path: Path) -> bool:
    return (path / RUN_META).is_file()


def has_nmr_extract_signature(path: Path) -> bool:
    if not path.is_dir():
        return False
    npys_dir = path / "npys"
    if not (
        (path / "trajectory.h5").is_file()
        and (path / "extraction_manifest.json").is_file()
        and npys_dir.is_dir()
        and (path / "pdbs").is_dir()
    ):
        return False
    try:
        return any(child.name.startswith("frame_") for child in npys_dir.iterdir())
    except OSError:
        return False


def nmr_extract_signature_ancestor(path: Path) -> Path | None:
    current = path if path.is_dir() else path.parent
    while True:
        if has_nmr_extract_signature(current):
            return current
        parent = current.parent
        if parent == current:
            return None
        current = parent


def delete_refusal_reason(path: Path, managed_root: Path) -> str | None:
    resolved_root = managed_root.resolve()
    resolved_path = path.resolve()
    if resolved_path.is_relative_to(SHARED_ROOT):
        return f"resolved path {resolved_path} is under /shared"
    if resolved_path == resolved_root or not resolved_path.is_relative_to(resolved_root):
        return f"resolved path {resolved_path} is outside managed run root {resolved_root}"
    signature_dir = nmr_extract_signature_ancestor(resolved_path)
    if signature_dir is not None:
        return (
            "nmr_extract extraction signature at "
            f"{signature_dir} (trajectory.h5 + extraction_manifest.json + "
            "npys/frame_* + pdbs); categorically off the cleanup table"
        )
    return None


def guarded_unlink(path: Path, managed_root: Path) -> tuple[bool, int]:
    reason = delete_refusal_reason(path, managed_root)
    if reason is not None:
        print(f"REFUSE delete {path}: {reason}", file=sys.stderr)
        return False, 0
    try:
        bytes_deleted = path.stat().st_size
    except FileNotFoundError:
        return False, 0
    try:
        path.unlink()
    except FileNotFoundError:
        return False, 0
    return True, bytes_deleted


def scan_run(root: Path, path: Path, *, include_patterns: bool = False) -> dict[str, Any]:
    meta = load_meta(path)
    substrate = discover_substrate(path, include_patterns=include_patterns)
    substrate_paths = [path / name for name in substrate]
    analysis_path = path / ANALYSIS_DIR
    substrate_bytes = size_of(substrate_paths)
    analysis_bytes = dir_size(analysis_path)
    total_bytes = dir_size(path)
    newest = max_mtime(path)
    created = meta.get("created_at") or datetime.fromtimestamp(
        path.stat().st_ctime, timezone.utc
    ).isoformat(timespec="seconds")
    completed = meta.get("completed_at")
    return {
        "name": path.name,
        "path": str(path),
        "managed": is_managed_run(path),
        "created_at": created,
        "completed_at": completed,
        "updated_at": meta.get("updated_at"),
        "status": meta.get("status", "unknown"),
        "mtime": datetime.fromtimestamp(newest, timezone.utc).isoformat(timespec="seconds"),
        "size_bytes": total_bytes,
        "substrate_bytes": substrate_bytes,
        "analysis_bytes": analysis_bytes,
        "substrate_present": bool(substrate),
        "substrate_file_count": len(substrate),
    }


def scan_runs(root: Path, *, include_patterns: bool = False) -> list[dict[str, Any]]:
    if not root.exists():
        return []
    runs: list[dict[str, Any]] = []
    for child in sorted(root.iterdir()):
        if child.is_dir() and is_managed_run(child):
            runs.append(scan_run(root, child, include_patterns=include_patterns))
    return runs


def write_root_manifest(root: Path, runs: list[dict[str, Any]]) -> None:
    data = {
        "framework": "rediscover-runs",
        "version": 1,
        "root": str(root),
        "generated_at": utc_now(),
        "keep_substrate_default": DEFAULT_KEEP_SUBSTRATE,
        "runs": runs,
    }
    write_json(root / ROOT_MANIFEST, data)


def print_manifest(runs: list[dict[str, Any]]) -> None:
    if not runs:
        print("no managed rediscover runs")
        return
    header = (
        f"{'name':30} {'status':12} {'substrate':10} "
        f"{'sub_size':>10} {'analysis':>10} {'updated'}"
    )
    print(header)
    print("-" * len(header))
    for run in runs:
        print(
            f"{run['name'][:30]:30} "
            f"{run['status'][:12]:12} "
            f"{str(run['substrate_present']):10} "
            f"{human_size(int(run['substrate_bytes'])):>10} "
            f"{human_size(int(run['analysis_bytes'])):>10} "
            f"{run.get('updated_at') or run.get('mtime') or ''}"
        )


def sort_key_for_retention(run: dict[str, Any]) -> float:
    for key in ("completed_at", "updated_at", "created_at"):
        parsed = parse_time(run.get(key))
        if parsed is not None:
            return parsed
    mtime = parse_time(run.get("mtime"))
    return mtime if mtime is not None else 0.0


def is_active(run: dict[str, Any], *, active_minutes: int) -> bool:
    if run.get("status") == "running":
        return True
    mtime = parse_time(run.get("mtime"))
    if mtime is None:
        return False
    return (time.time() - mtime) < (active_minutes * 60)


def clean_old_substrate(
    root: Path,
    *,
    keep: int,
    active_minutes: int,
    force: bool,
    include_patterns: bool = False,
) -> tuple[list[dict[str, Any]], int]:
    root.mkdir(parents=True, exist_ok=True)
    managed_root = root.resolve()
    runs = scan_runs(root, include_patterns=include_patterns)
    with_substrate = [r for r in runs if r["substrate_present"]]
    protected_active = {
        r["name"] for r in with_substrate if is_active(r, active_minutes=active_minutes)
    }
    eligible = [
        r for r in with_substrate
        if r["name"] not in protected_active and r.get("status") != "failed"
    ]
    keep_names = {
        r["name"]
        for r in sorted(eligible, key=sort_key_for_retention, reverse=True)[: max(keep, 0)]
    }
    candidates = [
        r for r in eligible
        if r["name"] not in keep_names
    ]

    deleted_bytes = 0
    for run in candidates:
        path = Path(str(run["path"]))
        substrate = discover_substrate(path, include_patterns=include_patterns)
        files = [path / name for name in substrate]
        bytes_for_run = size_of(files)
        print(
            f"{'DELETE' if force else 'DRY-RUN'} substrate "
            f"{run['name']} ({len(files)} files, {human_size(bytes_for_run)})"
        )
        for file_path in files:
            print(f"  {file_path}")
        if not force:
            continue
        deleted_names: list[str] = []
        deleted_bytes_for_run = 0
        for name, file_path in zip(substrate, files):
            was_deleted, bytes_deleted = guarded_unlink(file_path, managed_root)
            if not was_deleted:
                continue
            deleted_names.append(name)
            deleted_bytes_for_run += bytes_deleted
        if not deleted_names:
            continue
        meta = load_meta(path)
        drops = meta.setdefault("substrate_drops", [])
        if isinstance(drops, list):
            drops.append(
                {
                    "dropped_at": utc_now(),
                    "files": deleted_names,
                    "bytes": deleted_bytes_for_run,
                    "reason": f"superseded keep_substrate={keep}",
                }
            )
        meta["substrate_dropped_at"] = utc_now()
        save_meta(path, meta)
        deleted_bytes += deleted_bytes_for_run

    if protected_active:
        for name in sorted(protected_active):
            print(f"SKIP active/recent run {name}")
    if not candidates:
        print("no superseded substrate candidates")

    runs = scan_runs(root, include_patterns=include_patterns)
    write_root_manifest(root, runs)
    return candidates, deleted_bytes


def has_option(args: list[str], option: str) -> bool:
    return any(arg == option or arg.startswith(f"{option}=") for arg in args)


def script_supports_option(script: Path, option: str) -> bool:
    try:
        text = script.read_text(encoding="utf-8")
    except OSError:
        return False
    quoted_double = f'"{option}"'
    quoted_single = f"'{option}'"
    return quoted_double in text or quoted_single in text


def cmd_emit(args: argparse.Namespace) -> int:
    root = resolve_root(args.root)
    path = run_dir(root, args.name)
    root.mkdir(parents=True, exist_ok=True)

    existing_substrate = discover_substrate(path, include_patterns=True) if path.exists() else []
    if path.exists() and existing_substrate and not args.allow_existing:
        raise SystemExit(
            f"{path} already has substrate; choose a fresh run name or pass --allow-existing"
        )

    path.mkdir(parents=True, exist_ok=True)
    extractor_args = clean_script_args(args.extractor_args)
    bin_path = Path(args.bin).expanduser()
    command = [
        str(bin_path),
        "--run",
        str(args.run),
        args.extract_out_option,
        str(path),
        "--case",
        args.case,
        *extractor_args,
    ]
    meta = load_meta(path)
    meta.update(
        {
            "framework": "rediscover-runs",
            "framework_version": 1,
            "name": args.name,
            "created_at": meta.get("created_at") or utc_now(),
            "status": "running",
            "extractor": {
                "bin": str(bin_path),
                "run": str(args.run),
                "case": args.case,
                "out_option": args.extract_out_option,
                "extra_args": extractor_args,
            },
            "command": command,
        }
    )
    save_meta(path, meta)
    print(f"run_dir={path}")
    print("+ " + " ".join(command))
    rc = subprocess.run(command, cwd=args.cwd).returncode

    meta = load_meta(path)
    meta["exit_code"] = rc
    if rc == 0:
        substrate = discover_substrate(path, include_patterns=True)
        meta["status"] = "complete"
        meta["completed_at"] = utc_now()
        meta["substrate_files"] = substrate
    else:
        meta["status"] = "failed"
        meta["failed_at"] = utc_now()
    save_meta(path, meta)

    runs = scan_runs(root, include_patterns=True)
    write_root_manifest(root, runs)

    if rc == 0 and not args.no_auto_clean:
        clean_old_substrate(
            root,
            keep=args.keep_substrate,
            active_minutes=args.active_minutes,
            force=True,
            include_patterns=True,
        )
    return rc


def cmd_analyze(args: argparse.Namespace) -> int:
    root = resolve_root(args.root)
    path = run_dir(root, args.name)
    if not is_managed_run(path):
        raise SystemExit(f"{path} is not a managed rediscover run")
    if not discover_substrate(path, include_patterns=True) and not args.allow_missing_substrate:
        raise SystemExit(f"{path} has no substrate files; cannot run substrate consumer")

    script = Path(args.script).expanduser()
    if not script.is_file():
        raise SystemExit(f"analysis script not found: {script}")

    script_args = clean_script_args(args.script_args)
    artifact_name = args.artifact_name or script.stem
    artifact_dir = path / ANALYSIS_DIR / artifact_name
    artifact_dir.mkdir(parents=True, exist_ok=True)

    injected: list[str] = []
    if not args.no_out_dir and not has_option(script_args, "--out-dir"):
        injected.extend(["--out-dir", str(path)])
    if (
        not args.no_artifact_dir
        and script_supports_option(script, "--artifact-dir")
        and not has_option(script_args, "--artifact-dir")
    ):
        injected.extend(["--artifact-dir", str(artifact_dir)])
    if (
        not args.no_report_md
        and script_supports_option(script, "--report-md")
        and not has_option(script_args, "--report-md")
    ):
        injected.extend(["--report-md", str(artifact_dir / f"{artifact_name}.md")])

    command = [args.python, str(script), *injected, *script_args]
    print(f"run_dir={path}")
    print(f"artifact_dir={artifact_dir}")
    print("+ " + " ".join(command))
    rc = subprocess.run(command, cwd=args.cwd).returncode

    meta = load_meta(path)
    history = meta.setdefault("analysis_runs", [])
    if isinstance(history, list):
        history.append(
            {
                "script": str(script),
                "artifact_dir": str(artifact_dir),
                "command": command,
                "exit_code": rc,
                "ran_at": utc_now(),
            }
        )
    save_meta(path, meta)

    runs = scan_runs(root, include_patterns=True)
    write_root_manifest(root, runs)
    return rc


def cmd_clean(args: argparse.Namespace) -> int:
    root = resolve_root(args.root)
    _, deleted = clean_old_substrate(
        root,
        keep=args.keep_substrate,
        active_minutes=args.active_minutes,
        force=args.force,
        include_patterns=args.include_patterns,
    )
    if args.force:
        print(f"deleted {human_size(deleted)}")
    else:
        print("dry run only; pass --force to delete")
    return 0


def cmd_manifest(args: argparse.Namespace) -> int:
    root = resolve_root(args.root)
    root.mkdir(parents=True, exist_ok=True)
    runs = scan_runs(root, include_patterns=args.include_patterns)
    write_root_manifest(root, runs)
    if args.json:
        print(json.dumps({"root": str(root), "runs": runs}, indent=2, sort_keys=True))
    else:
        print_manifest(runs)
        print(f"manifest={root / ROOT_MANIFEST}")
    return 0


def cmd_path(args: argparse.Namespace) -> int:
    print(run_dir(resolve_root(args.root), args.name))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="rediscover run root (env: REDISCOVER_RUN_ROOT)",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    emit = sub.add_parser("emit", help="run the C++ extractor into <root>/<run-name>")
    emit.add_argument("name", help="run directory name")
    emit.add_argument("--run", required=True, type=Path, help="calcset directory or .LGS path")
    emit.add_argument("--bin", default=DEFAULT_EXTRACT_BIN, help="extractor binary path")
    emit.add_argument("--case", default="all", help="extractor --case value")
    emit.add_argument(
        "--extract-out-option",
        default="--out",
        choices=("--out", "--out-dir"),
        help="extractor output option; current C++ CLI uses --out",
    )
    emit.add_argument(
        "--keep-substrate",
        type=int,
        default=DEFAULT_KEEP_SUBSTRATE,
        help="number of completed, non-recent runs whose substrate is retained",
    )
    emit.add_argument(
        "--active-minutes",
        type=int,
        default=DEFAULT_ACTIVE_MINUTES,
        help="never clean runs touched within this many minutes",
    )
    emit.add_argument("--allow-existing", action="store_true")
    emit.add_argument("--no-auto-clean", action="store_true")
    emit.add_argument("--cwd", default=None, help="working directory for extractor")
    emit.add_argument("extractor_args", nargs=argparse.REMAINDER)
    emit.set_defaults(func=cmd_emit)

    analyze = sub.add_parser("analyze", help="run a Python analysis consumer for a run")
    analyze.add_argument("name", help="managed run name")
    analyze.add_argument("script", help="analysis script path")
    analyze.add_argument("--artifact-name", help="analysis/<artifact-name> directory")
    analyze.add_argument("--python", default=sys.executable)
    analyze.add_argument("--cwd", default=None, help="working directory for analysis")
    analyze.add_argument("--allow-missing-substrate", action="store_true")
    analyze.add_argument("--no-out-dir", action="store_true", help="do not inject --out-dir")
    analyze.add_argument("--no-artifact-dir", action="store_true")
    analyze.add_argument("--no-report-md", action="store_true")
    analyze.add_argument("script_args", nargs=argparse.REMAINDER)
    analyze.set_defaults(func=cmd_analyze)

    clean = sub.add_parser("clean", help="drop superseded substrate; dry-run by default")
    clean.add_argument(
        "--keep-substrate",
        type=int,
        default=DEFAULT_KEEP_SUBSTRATE,
        help="number of completed, non-recent runs whose substrate is retained",
    )
    clean.add_argument(
        "--active-minutes",
        type=int,
        default=DEFAULT_ACTIVE_MINUTES,
        help="never clean runs touched within this many minutes",
    )
    clean.add_argument("--force", action="store_true", help="actually delete substrate files")
    clean.add_argument(
        "--include-patterns",
        action="store_true",
        help="fallback to top-level *.csv/*.npy/*.npz in managed runs lacking manifests",
    )
    clean.set_defaults(func=cmd_clean)

    manifest = sub.add_parser("manifest", help="refresh and print the run index")
    manifest.add_argument("--json", action="store_true")
    manifest.add_argument(
        "--include-patterns",
        action="store_true",
        help="include top-level substrate-like files if a managed run lacks manifests",
    )
    manifest.set_defaults(func=cmd_manifest)

    path = sub.add_parser("path", help="print <root>/<run-name>")
    path.add_argument("name")
    path.set_defaults(func=cmd_path)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
