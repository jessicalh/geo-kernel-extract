#!/usr/bin/env python3
"""Exercise the actual portable runtime through Reader's REST and snapshot API.

Use staging fixtures only. Hash the fixture before and after to check that the
portable session preserves it. The optional ML check needs complete F006 inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import time
import urllib.request


def tree_digest(root: Path) -> dict[str, str]:
    values = {}
    for path in sorted(root.rglob("*")):
        if path.is_file():
            with path.open("rb") as stream:
                values[str(path.relative_to(root))] = hashlib.file_digest(stream, "sha256").hexdigest()
    return values


def reader_mount_options(parent_pid: int) -> str | None:
    """Inspect this launch's actual Reader process, not a separate probe container."""
    pending = [parent_pid]
    seen = set()
    while pending:
        pid = pending.pop()
        if pid in seen:
            continue
        seen.add(pid)
        try:
            arguments = Path(f"/proc/{pid}/cmdline").read_bytes().split(b"\0")
            if arguments and arguments[0].rsplit(b"/", 1)[-1] == b"h5reader":
                for line in Path(f"/proc/{pid}/mountinfo").read_text().splitlines():
                    fields = line.split()
                    if len(fields) > 6 and fields[4] == "/provenance":
                        return fields[5]
            pending.extend(int(child) for child in
                           Path(f"/proc/{pid}/task/{pid}/children").read_text().split())
        except (FileNotFoundError, ProcessLookupError):
            continue
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--launcher", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--fixture", required=True, help="Relative path below source root")
    parser.add_argument("--check-ml", action="store_true")
    parser.add_argument("--check-local-library-key", help="Open this entry directly from the full 176-entry catalog")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for path in (args.launcher, args.source_root, args.workspace, args.output):
        if str(path).startswith("/mnt/provenance"):
            parser.error("Acceptance uses staging only; do not inspect the Provenance mount")
    fixture = args.source_root / args.fixture
    if not fixture.is_file() or not fixture.resolve().is_relative_to(args.source_root.resolve()):
        parser.error("Fixture must be a file within the staging source root")
    hash_root = args.source_root / "01_Reader/datasets" if args.check_local_library_key else fixture.parent
    before = tree_digest(hash_root)
    args.output.mkdir(parents=True, exist_ok=True)
    command = ["/bin/sh", str(args.launcher), "--source-root", str(args.source_root),
               "--workspace", str(args.workspace), "--", "--rest", "0",
               "/provenance/" + args.fixture]
    if not os.environ.get("DISPLAY"):
        xvfb = shutil.which("xvfb-run")
        if not xvfb:
            parser.error("Use an X11 display or install Xvfb for the development test")
        command[:0] = [xvfb, "-a", "-s", "-screen 0 1600x1000x24 +extension GLX +extension RANDR +render -noreset"]
    log_path = args.output / "reader.log"
    report: dict = {}
    with log_path.open("wb") as log:
        environment = dict(os.environ)
        environment["H5READER_PORTABLE_KEEP_STDIO"] = "1"
        environment["H5READER_PORTABLE_NO_DIALOG"] = "1"
        environment["H5READER_PORTABLE_NONINTERACTIVE_SETUP"] = "1"
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                   env=environment, start_new_session=True)
        url = None
        try:
            deadline = time.monotonic() + 180
            while time.monotonic() < deadline:
                match = re.search(r"H5READER_REST_PORT=(\d+)", log_path.read_text(errors="replace"))
                if match:
                    url = "http://127.0.0.1:" + match.group(1)
                    break
                if process.poll() is not None:
                    raise RuntimeError(f"Reader exited {process.returncode}; see {log_path}")
                time.sleep(0.2)
            if url is None:
                raise RuntimeError(f"No REST handshake; see {log_path}")

            def request(path: str, body: dict | None = None):
                data = None if body is None else json.dumps(body).encode()
                with urllib.request.urlopen(urllib.request.Request(
                        url + path, data=data, headers={"Content-Type": "application/json"}), timeout=90) as response:
                    raw = response.read()
                    return raw if response.headers.get_content_type() == "image/png" else (json.loads(raw) if raw else None)

            report["health"] = request("/health")
            assert report["health"]["ok"]
            report["initial_state"] = request("/ui/state")
            report["source_mount_options"] = reader_mount_options(process.pid)
            for target in ("window", "scene"):
                png = request("/api/screenshot", {"target": target})
                assert png.startswith(b"\x89PNG\r\n\x1a\n") and len(png) > 1000
                (args.output / f"{target}.png").write_bytes(png)
                report[f"{target}_png_bytes"] = len(png)
            if args.check_ml:
                metric = request("/dashboard/metric", {
                    "descriptor_id": "ml:experimental_shielding_iso",
                    "anchor": {"atom": 16}, "modes": ["strip.scalar"]})
                deadline = time.monotonic() + 90
                while time.monotonic() < deadline:
                    tracks = request("/dashboard/display")["strip_tracks"]
                    track = next((row for row in tracks if row["signal_id"] == metric["id"]
                                  and row["descriptor_id"] == "ml:experimental_shielding_iso"), None)
                    if track and track.get("valid") and track["valid"][0] == 1:
                        report["cpu_ml_track"] = track
                        break
                    time.sleep(0.5)
                assert "cpu_ml_track" in report, "CPU inference did not produce a valid result"
                ml = request("/ui/state")["experimentalShieldingMl"]
                assert ml["activeDevice"] == "cpu"
                report["ml"] = ml
            if args.check_local_library_key:
                library = request("/api/trajectories")
                assert library["local_library"] is True
                assert len(library["entries"]) == 176
                assert library["source_root"] == "/provenance"
                report["library_entry_count"] = len(library["entries"])
                request("/api/trajectories/show", {})
                png = request("/diagnostics/screenshot", {"target": "widget", "object_name": "TrajectoryLibraryDialog"})
                (args.output / "local-library.png").write_bytes(png)
                key = args.check_local_library_key
                request("/api/trajectories/open", {"key": key})
                report["direct_run_state"] = request("/ui/state")
                assert report["direct_run_state"]["loaded"] is True
                assert report["direct_run_state"]["frames"] == 100
                request("/frame/set", {"frame": 99})
                assert request("/ui/state")["currentFrame"] == 99
                (args.output / "direct-trajectory.png").write_bytes(request("/api/screenshot", {"target": "scene"}))
                copies = args.workspace / "local-copies"
                report["dataset_copies_created"] = len(list(copies.iterdir())) if copies.exists() else 0
                assert report["dataset_copies_created"] == 0, "Direct viewing created a dataset copy"
            request("/shutdown", {})
            process.wait(timeout=30)
            assert process.returncode == 0, process.returncode
            after = tree_digest(hash_root)
            report["source_files_unchanged"] = before == after
            report["changed_source_paths"] = [name for name in sorted(before.keys() | after.keys())
                                              if before.get(name) != after.get(name)]
            assert report["source_files_unchanged"], "The fixture changed during the session"
            report["source_file_count"] = len(before)
            log_text = log_path.read_text(errors="replace")
            assert "llvmpipe" in log_text, "Software renderer was not reported"
            report["renderer"] = next(line for line in log_text.splitlines() if "OpenGL renderer string:" in line)
            report["backend"] = re.search(r"H5READER_PORTABLE_BACKEND=(\w+)", log_text).group(1)
            report["local_runtime_extracted"] = any((args.workspace / "runtime").iterdir())
            if report["backend"] == "apptainer":
                assert report["local_runtime_extracted"] is False
                assert report["source_mount_options"] and "ro" in report["source_mount_options"].split(",")
                report["source_bind_read_only"] = True
            report["ok"] = True
        finally:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
            (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key not in ("initial_state", "ml")}, indent=2))


if __name__ == "__main__":
    main()
