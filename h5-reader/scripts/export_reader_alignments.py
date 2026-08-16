"""Run the accepted179 Reader alignment export one member at a time."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
import uuid
from datetime import datetime, timezone
from pathlib import Path


ACCEPTED_ROSTER_SHA256 = (
    "60bc861ef913f10f13f1ab2c44ecca77d58842dc10b9192a800d5a5fe3cb2a5d"
)
ACCEPTED_MEMBER_COUNT = 179
PILOT_MEMBERS = ("bmr4095", "bmr5563", "bmr15784", "bmr4664")
PORT_PREFIX = b"H5READER_REST_PORT="
MEMBERS_TABLE_HEADER = (
    "member_id\trun_id\tframes\tatoms\tprimary_fit_atoms\tca_fit_atoms\t"
    "primary_iterations\tprimary_final_delta_A\tgenerated_at_utc"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while block := source.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def read_authoritative_roster(path: Path) -> list[str]:
    if not path.is_file():
        raise RuntimeError(f"accepted roster does not exist: {path}")
    digest = sha256_file(path)
    if digest != ACCEPTED_ROSTER_SHA256:
        raise RuntimeError(
            f"accepted roster SHA-256 mismatch: expected {ACCEPTED_ROSTER_SHA256}, "
            f"got {digest}"
        )
    members = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    if len(members) != ACCEPTED_MEMBER_COUNT or any(not member for member in members):
        raise RuntimeError(
            f"accepted roster must contain exactly {ACCEPTED_MEMBER_COUNT} nonempty rows"
        )
    if len(set(members)) != len(members):
        raise RuntimeError("accepted roster contains duplicate member IDs")
    return members


def copy_roster_once(source: Path, output_root: Path) -> Path:
    destination = output_root / "accepted179.txt"
    source_bytes = source.read_bytes()
    if destination.exists():
        if not destination.is_file() or destination.read_bytes() != source_bytes:
            raise RuntimeError(
                f"output roster exists but is not the authoritative byte copy: {destination}"
            )
        return destination

    temporary = output_root / f".accepted179.txt.tmp.{uuid.uuid4().hex}"
    temporary.write_bytes(source_bytes)
    os.replace(temporary, destination)
    return destination


def unique_lgs(source_root: Path, member: str) -> Path:
    member_directory = source_root / member
    if not member_directory.is_dir():
        raise RuntimeError(
            f"source member directory does not exist: {member_directory}"
        )
    matches = sorted(member_directory.glob("*.lgs"))
    if len(matches) != 1:
        raise RuntimeError(
            f"source member must contain exactly one LGS, found {len(matches)}: "
            f"{member_directory}"
        )
    return matches[0]


def request_json(
    base_url: str,
    method: str,
    route: str,
    body: dict[str, object] | None,
    timeout_seconds: float,
) -> tuple[int, dict[str, object]]:
    payload = None
    headers: dict[str, str] = {}
    if body is not None:
        payload = json.dumps(body).encode("utf-8")
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(
        f"{base_url}{route}", data=payload, headers=headers, method=method
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
            status = response.status
            data = response.read()
    except urllib.error.HTTPError as error:
        status = error.code
        data = error.read()
    try:
        decoded = json.loads(data.decode("utf-8")) if data else {}
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RuntimeError(
            f"Reader returned non-JSON data for {method} {route} (HTTP {status})"
        ) from error
    if not isinstance(decoded, dict):
        raise RuntimeError(
            f"Reader returned a non-object JSON value for {method} {route}"
        )
    return status, decoded


def wait_for_handshake(
    process: subprocess.Popen[bytes], log_path: Path, timeout_seconds: float
) -> int:
    deadline = time.monotonic() + timeout_seconds
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise RuntimeError(
                f"Reader exited before the REST handshake with code {process.returncode}"
            )
        try:
            log_bytes = log_path.read_bytes()
        except OSError:
            log_bytes = b""
        for line in log_bytes.splitlines():
            marker = line.find(PORT_PREFIX)
            if marker >= 0:
                value = line[marker + len(PORT_PREFIX) :].strip()
                if value.isdigit():
                    return int(value)
        time.sleep(0.1)
    raise RuntimeError(
        f"timed out after {timeout_seconds:g}s waiting for the Reader REST handshake"
    )


def wait_for_health(base_url: str, timeout_seconds: float = 15.0) -> None:
    deadline = time.monotonic() + timeout_seconds
    while time.monotonic() < deadline:
        try:
            status, body = request_json(base_url, "GET", "/health", None, 2.0)
            if status == 200 and body.get("ok") is True:
                return
        except (OSError, RuntimeError):
            pass
        time.sleep(0.1)
    raise RuntimeError("Reader announced a REST port but did not become healthy")


def log_tail(path: Path, byte_count: int = 8000) -> str:
    try:
        return path.read_bytes()[-byte_count:].decode("utf-8", errors="replace")
    except OSError as error:
        return f"<could not read log: {error}>"


def stop_reader(process: subprocess.Popen[bytes], base_url: str | None) -> None:
    if process.poll() is not None:
        return
    if base_url is not None:
        try:
            request_json(base_url, "POST", "/shutdown", {}, 10.0)
        except (OSError, RuntimeError):
            pass
    try:
        process.wait(timeout=20.0)
        return
    except subprocess.TimeoutExpired:
        process.terminate()
    try:
        process.wait(timeout=10.0)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait(timeout=10.0)


def run_member(
    reader: Path,
    lgs_path: Path,
    output_root: Path,
    apply_display: bool,
    startup_timeout_seconds: float,
    export_timeout_seconds: float,
) -> dict[str, object]:
    log_file = tempfile.NamedTemporaryFile(
        prefix=f"h5reader_alignment_{lgs_path.parent.name}_",
        suffix=".log",
        delete=False,
    )
    log_path = Path(log_file.name)
    log_file.close()
    log_handle = log_path.open("wb")
    process = subprocess.Popen(
        [str(reader), "--rest", "0", str(lgs_path)],
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        env={**os.environ},
    )
    base_url: str | None = None
    succeeded = False
    try:
        port = wait_for_handshake(process, log_path, startup_timeout_seconds)
        base_url = f"http://127.0.0.1:{port}"
        wait_for_health(base_url)
        status, response = request_json(
            base_url,
            "POST",
            "/api/alignment/export",
            {"output_root": str(output_root), "apply_display": apply_display},
            export_timeout_seconds,
        )
        if status != 200 or response.get("ok") is not True:
            raise RuntimeError(
                f"alignment export failed with HTTP {status}: "
                f"{json.dumps(response, sort_keys=True)}"
            )
        succeeded = True
        return response
    except Exception as error:
        raise RuntimeError(
            f"{error}\nReader log retained at {log_path}\n{log_tail(log_path)}"
        ) from error
    finally:
        stop_reader(process, base_url)
        log_handle.close()
        if succeeded and process.returncode == 0:
            try:
                log_path.unlink()
            except OSError:
                pass


def write_complete_record(
    output_root: Path, roster_copy: Path, member_count: int
) -> None:
    members_table = output_root / "members.tsv"
    if not members_table.is_file():
        raise RuntimeError(f"members.tsv is missing: {members_table}")
    record = {
        "schema_version": 1,
        "complete": True,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(
            timespec="milliseconds"
        ),
        "member_count": member_count,
        "accepted_roster": roster_copy.name,
        "accepted_roster_sha256": sha256_file(roster_copy),
        "members_tsv_sha256": sha256_file(members_table),
        "contains_targets": False,
        "protected_targets_opened": 0,
    }
    destination = output_root / "COMPLETE.json"
    temporary = output_root / f".COMPLETE.json.tmp.{uuid.uuid4().hex}"
    temporary.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, destination)


def read_member_table(members_table: Path) -> list[str]:
    if not members_table.is_file():
        raise RuntimeError(f"members.tsv is missing: {members_table}")
    lines = members_table.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != MEMBERS_TABLE_HEADER:
        raise RuntimeError(f"members.tsv has an unexpected header: {members_table}")
    table_members = [line.split("\t", 1)[0] for line in lines[1:] if line]
    if len(table_members) != len(set(table_members)):
        raise RuntimeError("members.tsv contains duplicate member IDs")
    return table_members


def validate_resume_prefix(
    output_root: Path, members_table: Path, completed_prefix: list[str]
) -> None:
    table_members = set(read_member_table(members_table))
    for member in completed_prefix:
        if member not in table_members:
            raise RuntimeError(f"resume prefix is absent from members.tsv: {member}")
        manifest_path = output_root / member / "export.json"
        if not manifest_path.is_file():
            raise RuntimeError(f"resume prefix has no export.json: {manifest_path}")
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
            raise RuntimeError(f"could not read completed manifest: {manifest_path}") from error
        identity = manifest.get("identity")
        if (
            manifest.get("complete") is not True
            or not isinstance(identity, dict)
            or identity.get("member_id") != member
            or manifest.get("contains_targets") is not False
            or manifest.get("protected_targets_opened") != 0
        ):
            raise RuntimeError(f"resume prefix manifest is not complete: {manifest_path}")


def validate_complete_member_table(
    members_table: Path, accepted_members: list[str]
) -> None:
    table_members = read_member_table(members_table)
    if set(table_members) != set(accepted_members):
        missing = sorted(set(accepted_members) - set(table_members))
        extra = sorted(set(table_members) - set(accepted_members))
        raise RuntimeError(
            f"members.tsv does not match accepted179; missing={missing}, extra={extra}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Export the four fixed pilots or the complete accepted179 roster through "
            "one Reader process per member."
        )
    )
    parser.add_argument("--reader", required=True, type=Path)
    parser.add_argument("--roster", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--mode", required=True, choices=("pilot", "full"))
    parser.add_argument(
        "--start-at",
        help=(
            "full-mode member at which to resume; every earlier roster member must "
            "already have a completed manifest and members.tsv row"
        ),
    )
    parser.add_argument(
        "--startup-timeout",
        type=float,
        default=6 * 60 * 60,
        help="catastrophe backstop in seconds while Reader loads a member",
    )
    parser.add_argument(
        "--export-timeout",
        type=float,
        default=6 * 60 * 60,
        help="catastrophe backstop in seconds while Reader exports a member",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    reader = args.reader.resolve()
    roster = args.roster.resolve()
    source_root = args.source_root.resolve()
    output_root = args.output_root.resolve()
    if not reader.is_file():
        raise RuntimeError(f"Reader binary does not exist: {reader}")
    if not source_root.is_dir():
        raise RuntimeError(f"source root does not exist: {source_root}")
    if args.startup_timeout <= 0.0 or args.export_timeout <= 0.0:
        raise RuntimeError("timeouts must be positive")

    accepted_members = read_authoritative_roster(roster)
    members = list(PILOT_MEMBERS) if args.mode == "pilot" else accepted_members
    missing_pilots = set(PILOT_MEMBERS) - set(accepted_members)
    if missing_pilots:
        raise RuntimeError(
            f"pilot members are absent from the accepted roster: {sorted(missing_pilots)}"
        )

    output_root.mkdir(parents=True, exist_ok=True)
    roster_copy = copy_roster_once(roster, output_root)

    start_position = 0
    if args.start_at is not None:
        if args.mode != "full":
            raise RuntimeError("--start-at is only valid with --mode full")
        try:
            start_position = accepted_members.index(args.start_at)
        except ValueError as error:
            raise RuntimeError(
                f"--start-at member is absent from accepted179: {args.start_at}"
            ) from error
        completed_prefix = accepted_members[:start_position]
        validate_resume_prefix(
            output_root, output_root / "members.tsv", completed_prefix
        )
        members = accepted_members[start_position:]

    print(
        f"Reader alignment {args.mode}: {len(members)} member(s) this invocation\n"
        f"  reader: {reader}\n"
        f"  source: {source_root}\n"
        f"  output: {output_root}\n"
        f"  roster position: {start_position + 1}/{len(accepted_members)}",
        flush=True,
    )
    completed_members = (
        accepted_members[:start_position] if args.mode == "full" else []
    )
    total_members = len(accepted_members) if args.mode == "full" else len(members)
    for position, member in enumerate(members, start=start_position + 1):
        lgs_path = unique_lgs(source_root, member)
        started = time.monotonic()
        print(f"[{position}/{total_members}] {member}: {lgs_path.name}", flush=True)
        response = run_member(
            reader,
            lgs_path,
            output_root,
            apply_display=args.mode == "pilot",
            startup_timeout_seconds=args.startup_timeout,
            export_timeout_seconds=args.export_timeout,
        )
        manifest = response.get("manifest")
        if not isinstance(manifest, dict):
            raise RuntimeError(
                f"Reader did not return the validated manifest for {member}"
            )
        identity = manifest.get("identity")
        if not isinstance(identity, dict) or identity.get("member_id") != member:
            raise RuntimeError(
                f"Reader returned the wrong completed member for {member}"
            )
        completed_members.append(member)
        primary = manifest["primary_fit"]
        dimensions = manifest["dimensions"]
        validation = manifest["validation"]
        elapsed = time.monotonic() - started
        print(
            f"  complete in {elapsed:.1f}s; resume={response.get('already_complete')}; "
            f"atoms={dimensions['atoms']}; fit={dimensions['primary_fit_atoms']}; "
            f"iterations={primary['iterations']}; delta={primary['final_delta_A']:.6g} A; "
            f"det=[{validation['determinant_min']:.12g},"
            f"{validation['determinant_max']:.12g}]",
            flush=True,
        )

    if args.mode == "full":
        if completed_members != accepted_members:
            raise RuntimeError(
                "Reader did not validate the complete accepted roster in order"
            )
        validate_complete_member_table(output_root / "members.tsv", accepted_members)
        write_complete_record(output_root, roster_copy, len(accepted_members))
        print(f"wrote {output_root / 'COMPLETE.json'}", flush=True)
    else:
        print("pilot complete; COMPLETE.json was not created", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        raise SystemExit(130)
    except Exception as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
