#!/usr/bin/env python3
"""Build an offline Reader catalog from the publisher's neutral SSD inventories.

This build-time utility opens only the explicitly supplied report files and,
optionally, staged LGS metadata under --staged-lgs-root. It never opens, stats,
or resolves source-drive paths. The output contains drive-relative paths, not
machine mount points or URLs, and is independent of the HTTPS catalog.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import stat
import sys


EXPECTED = {"datasets": 176, "MD": 173, "small_proteins": 3, "DFT_frame_pairs": 7253}
DATASET_COLUMNS = {
    "member", "dataset_id", "title", "frames", "dft_frames",
    "disk_member_directory", "disk_lgs", "disk_lgs_relative_to_drive", "staged_lgs",
}
DFT_COLUMNS = {
    "member", "frame_index", "disk_metadata", "disk_primary_output",
    "portable_metadata", "portable_primary_output",
}
SMALL_PROTEIN_TITLES = {
    "1p9j_T1E_1P9J_5801": "1P9J",
    "trpcage_TC5b_1L2Y_5292": "Trp-cage TC5b",
    "chignolin_CLN025_5AWL_0": "Chignolin CLN025",
}


class ContractError(ValueError):
    """The reports cannot safely describe a complete local library."""


def relative_path(value: str, base: str = "") -> str:
    """Normalize a lexical relative path, allowing .. only inside the drive.

    LGS DFT references deliberately walk from 01_Reader into 04_MD_records.
    Do not resolve symlinks or consult the mounted source drive here.
    """
    if not value or value.startswith("/") or "\\" in value or ":" in value:
        raise ContractError(f"Expected a relative POSIX path: {value!r}")
    if any(ord(char) < 32 for char in value):
        raise ContractError("Control character in a path")
    parts = []
    for component in (base + "/" + value).split("/"):
        if component in ("", "."):
            continue
        if component == "..":
            if not parts:
                raise ContractError(f"Path escapes the drive root: {value!r}")
            parts.pop()
        else:
            parts.append(component)
    if not parts:
        raise ContractError(f"Path refers to the drive root: {value!r}")
    return "/".join(parts)


def from_drive(value: str, source_prefix: str) -> str:
    prefix = source_prefix.rstrip("/") + "/"
    if not value.startswith(prefix):
        raise ContractError(f"Source path is outside {source_prefix}: {value!r}")
    return relative_path(value[len(prefix):])


def integer(value: str, label: str, minimum: int = 0) -> int:
    if not re.fullmatch(r"[0-9]+", value):
        raise ContractError(f"{label} must be an integer: {value!r}")
    number = int(value)
    if number < minimum:
        raise ContractError(f"{label} must be at least {minimum}")
    return number


def read_tsv(data: bytes, required: set[str], name: str) -> list[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")), delimiter="\t")
    fields = reader.fieldnames or []
    if len(fields) != len(set(fields)) or not required.issubset(fields):
        raise ContractError(f"Missing or duplicate columns in {name}")
    result = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in result):
        raise ContractError(f"Malformed TSV row in {name}")
    return result


def reject_symlinks(path: Path) -> None:
    """Inspect staging path components without following any symlink."""
    current = Path(path.absolute().anchor)
    for component in path.absolute().parts[1:]:
        current /= component
        if stat.S_ISLNK(current.lstat().st_mode):
            raise ContractError(f"Staged LGS path contains a symlink: {current}")


def build_catalog(
    inventory: bytes,
    dft_pairs: bytes,
    contract: bytes,
    source_prefix: str = "/mnt/provenance",
    staged_lgs_root: Path | None = None,
) -> tuple[dict, dict]:
    """Validate the full corpus and return catalog and optional DFT roster."""
    metadata = json.loads(contract)
    if metadata.get("records") != EXPECTED:
        raise ContractError(f"Contract records must equal {EXPECTED}")
    rows = read_tsv(inventory, DATASET_COLUMNS, "dataset inventory")
    companions = read_tsv(dft_pairs, DFT_COLUMNS, "DFT companion inventory")
    entries = {}
    paths = set()
    counts = {key: 0 for key in EXPECTED}
    for row in rows:
        member = row["member"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", member):
            raise ContractError(f"Invalid member identifier: {member!r}")
        if member in entries:
            raise ContractError(f"Duplicate member: {member}")
        source = relative_path(row["disk_lgs_relative_to_drive"])
        if source != from_drive(row["disk_lgs"], source_prefix):
            raise ContractError(f"Inconsistent LGS paths for {member}")
        directory = from_drive(row["disk_member_directory"], source_prefix)
        if str(PurePosixPath(source).parent) != directory:
            raise ContractError(f"LGS must be directly inside its member directory: {member}")
        parts = PurePosixPath(directory).parts
        if len(parts) != 4 or parts[:2] != ("01_Reader", "datasets") or parts[3] != member:
            raise ContractError(f"Unexpected member directory for {member}: {directory}")
        group = {"MD_173": "MD", "Small_proteins": "small_protein"}.get(parts[2])
        if group is None or PurePosixPath(source).suffix.lower() != ".lgs":
            raise ContractError(f"Unexpected dataset group or LGS name: {source}")
        if source in paths:
            raise ContractError(f"Duplicate LGS path: {source}")
        paths.add(source)
        if not row["dataset_id"].strip() or not row["title"].strip():
            raise ContractError(f"Missing dataset identity or title for {member}")
        frames = integer(row["frames"], f"{member} frames", 1)
        dft_frames = integer(row["dft_frames"], f"{member} DFT frames")
        if dft_frames > frames:
            raise ContractError(f"DFT frame count exceeds trajectory frame count: {member}")
        description = f"{member} · {frames:,} frames · " + ("MD trajectory" if group == "MD" else "small protein")
        if dft_frames:
            description += f" · {dft_frames:,} DFT frames"
        bmrb = re.fullmatch(r"bmr([0-9]+)", member, flags=re.IGNORECASE)
        title = f"BMRB {bmrb.group(1)}" if bmrb else SMALL_PROTEIN_TITLES.get(member, row["title"])
        if member in SMALL_PROTEIN_TITLES:
            description = f"{row['title']} · {frames:,} frames · {dft_frames:,} DFT frames"
        entries[member] = {
            "key": f"local-{member}", "member": member,
            "dataset_id": row["dataset_id"], "protein_id": member,
            "title": title, "source_title": row["title"], "description": description, "group": group,
            "frames": frames, "dft_frames": dft_frames, "source_lgs": source,
        }
        counts["datasets"] += 1
        counts["MD" if group == "MD" else "small_proteins"] += 1

    rosters = {member: {} for member in entries}
    portable_paths = set()
    normalized_pairs = []
    for row in companions:
        member = row["member"]
        if member not in entries:
            raise ContractError(f"DFT row names unknown member: {member}")
        frame = integer(row["frame_index"], f"{member} frame index")
        if frame in rosters[member]:
            raise ContractError(f"Duplicate DFT frame for {member}: {frame}")
        source_metadata = from_drive(row["disk_metadata"], source_prefix)
        source_output = from_drive(row["disk_primary_output"], source_prefix)
        if not source_metadata.startswith("04_MD_records/") or not source_output.startswith("04_MD_records/"):
            raise ContractError(f"DFT records must be beneath 04_MD_records: {member}/{frame}")
        destination_metadata = relative_path(row["portable_metadata"])
        destination_output = relative_path(row["portable_primary_output"])
        expected_parent = f"dft/frame_{frame:06d}"
        for source, destination in ((source_metadata, destination_metadata), (source_output, destination_output)):
            if str(PurePosixPath(destination).parent) != expected_parent:
                raise ContractError(f"Unexpected portable DFT directory: {destination}")
            if PurePosixPath(source).name != PurePosixPath(destination).name:
                raise ContractError(f"Portable DFT filename differs from source: {destination}")
            identity = (member, destination)
            if identity in portable_paths:
                raise ContractError(f"Duplicate portable DFT path for {member}: {destination}")
            portable_paths.add(identity)
        pair = {
            "member": member, "frame_index": frame,
            "source_metadata": source_metadata, "source_primary_output": source_output,
            "portable_metadata": destination_metadata, "portable_primary_output": destination_output,
        }
        rosters[member][frame] = pair
        normalized_pairs.append(pair)
        counts["DFT_frame_pairs"] += 1

    for member, entry in entries.items():
        if len(rosters[member]) != entry["dft_frames"]:
            raise ContractError(f"Incomplete DFT roster for {member}: expected {entry['dft_frames']}, got {len(rosters[member])}")
        if staged_lgs_root is not None:
            staged = staged_lgs_root / entry["source_lgs"]
            # Inspect SSD copies only. Never use disk_lgs or follow a staged symlink.
            reject_symlinks(staged)
            document = json.loads(staged.read_text(encoding="utf-8"))
            if document.get("kind") != "trajectory" or document.get("dataset_id") != entry["dataset_id"]:
                raise ContractError(f"Staged LGS identity differs for {member}")
            entry["protein_id"] = document.get("protein_id") or member
            parent = str(PurePosixPath(entry["source_lgs"]).parent)
            staged_roster = {}
            for frame in document.get("dft", {}).get("frames", []):
                index = frame.get("frame_index")
                if isinstance(index, bool) or not isinstance(index, int) or index < 0 or index in staged_roster:
                    raise ContractError(f"Invalid or duplicate staged DFT frame for {member}: {index!r}")
                staged_roster[index] = relative_path(frame["meta_json"], parent)
            expected_roster = {index: pair["source_metadata"] for index, pair in rosters[member].items()}
            if staged_roster != expected_roster:
                raise ContractError(f"Staged LGS DFT roster differs for {member}")
    if counts != EXPECTED:
        raise ContractError(f"Incomplete corpus: expected {EXPECTED}, got {counts}")
    total_frames = sum(entry["frames"] for entry in entries.values())
    verified_frames = metadata.get("data_state", {}).get("verified_counts", {}).get("frame_directories_checked")
    if verified_frames is not None and total_frames != verified_frames:
        raise ContractError(f"Frame total differs from verified contract: {total_frames} != {verified_frames}")
    fingerprints = {
        "dataset_inventory_sha256": hashlib.sha256(inventory).hexdigest(),
        "dft_inventory_sha256": hashlib.sha256(dft_pairs).hexdigest(),
        "contract_sha256": hashlib.sha256(contract).hexdigest(),
    }
    catalog = {
        "schema_version": 1, "kind": "h5reader-local-library",
        "records": {**counts, "total_frames": total_frames},
        "source_reports": fingerprints,
        "datasets": sorted(entries.values(), key=lambda entry: entry["key"]),
    }
    roster = {
        "schema_version": 1, "kind": "h5reader-local-dft-roster",
        "source_reports": fingerprints,
        "pairs": sorted(normalized_pairs, key=lambda pair: (pair["member"], pair["frame_index"])),
    }
    return catalog, roster


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--dft-pairs", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dft-output", type=Path)
    parser.add_argument("--source-prefix", default="/mnt/provenance", help="Lexical inventory mount prefix; never opened")
    parser.add_argument("--staged-lgs-root", type=Path, help="Optional SSD directory mirroring drive-relative LGS paths")
    args = parser.parse_args()
    try:
        catalog, roster = build_catalog(
            args.inventory.read_bytes(), args.dft_pairs.read_bytes(), args.contract.read_bytes(),
            args.source_prefix, args.staged_lgs_root,
        )
        for destination, document in ((args.output, catalog), (args.dft_output, roster)):
            if destination is not None:
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_text(json.dumps(document, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(json.dumps({"catalog": str(args.output), "records": catalog["records"]}))
        return 0
    except (ContractError, OSError, UnicodeError, json.JSONDecodeError, KeyError, TypeError) as error:
        print(f"Local catalog validation failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
