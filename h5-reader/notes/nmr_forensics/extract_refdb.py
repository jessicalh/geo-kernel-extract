#!/usr/bin/env python3
"""
extract_refdb.py — extract per-protein RefDB .str.corr files for the 10.

Source: `refdb_2016_snapshot/RefDB_files.tar.gz` (30 MB, SHA256 recorded in
provenance). Per protein, extracts
`data/www/lib/RefDB/RefDB/corrdata/bmr<id>.str.corr` into
`source_data/<pdb>_<bmrb>/raw/bmr<id>.str.corr`. Hash + extraction metadata
recorded in `raw_provenance.toml`.

Run from nmr_forensics/ in the venv:
    ./.venv/bin/python extract_refdb.py
"""

import hashlib
import sys
import tarfile
import tomllib
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).parent
SOURCE_DATA = HERE / "source_data"
TARBALL = HERE / "refdb_2016_snapshot" / "RefDB_files.tar.gz"

PROTEINS = [
    ("1DV0", 4757),
    ("1HS5", 4934),
    ("1HD6", 4820),
    ("1B1V", 4292),
    ("1HA9", 5028),
    ("1G26", 4656),
    ("1CBH", 192),
    ("1I2V", 4976),
    ("1I8X", 4351),
    ("1BWX", 3449),
]

TAR_PATH_TEMPLATE = "data/www/lib/RefDB/RefDB/corrdata/bmr{bmrb_id}.str.corr"


def tarball_sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def extract_one(tar, pdb_id, bmrb_id, tarball_hash, now_iso):
    tag = f"{pdb_id}_{bmrb_id}"
    protein_dir = SOURCE_DATA / tag
    raw_dir = protein_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    member_path = TAR_PATH_TEMPLATE.format(bmrb_id=bmrb_id)
    filename = f"bmr{bmrb_id}.str.corr"
    out_path = raw_dir / filename

    try:
        member = tar.getmember(member_path)
    except KeyError:
        raise RuntimeError(f"[{tag}] {member_path} not in RefDB tarball")

    fileobj = tar.extractfile(member)
    if fileobj is None:
        raise RuntimeError(f"[{tag}] {member_path} is not a regular file")
    content = fileobj.read()
    new_hash = hashlib.sha256(content).hexdigest()

    if out_path.exists():
        old_hash = hashlib.sha256(out_path.read_bytes()).hexdigest()
        if old_hash == new_hash:
            status = "unchanged"
        else:
            status = "content_changed"
    else:
        status = "new"
    out_path.write_bytes(content)

    entry = {
        "kind": "refdb_nmr_star_2_corr",
        "filename": filename,
        "url": "refdb_2016_snapshot/RefDB_files.tar.gz",
        "tar_member": member_path,
        "tarball_sha256": tarball_hash,
        "sha256": new_hash,
        "size_bytes": len(content),
        "extracted_at": now_iso,
        "status": status,
    }
    return tag, entry


def update_provenance(tag, new_refdb_entry):
    prov_path = SOURCE_DATA / tag / "raw_provenance.toml"

    existing_entries = []
    if prov_path.exists():
        existing = tomllib.loads(prov_path.read_text())
        for s in existing.get("sources", []):
            if s.get("kind") != "refdb_nmr_star_2_corr":
                existing_entries.append(s)

    all_entries = existing_entries + [new_refdb_entry]

    lines = [
        "# Raw-source provenance for per-protein experimental-shift inputs.",
        "# Auto-written by fetch_bmrb.py / extract_refdb.py. Do not hand-edit;",
        "# re-run the fetch scripts to refresh.",
        "",
        f'protein_tag = "{tag}"',
        "",
    ]
    for e in all_entries:
        lines.append("[[sources]]")
        for k in ("kind", "filename", "url", "tar_member"):
            if k in e:
                lines.append(f'{k} = "{e[k]}"')
        for k in ("tarball_sha256", "sha256"):
            if k in e:
                lines.append(f'{k} = "{e[k]}"')
        if "size_bytes" in e:
            lines.append(f'size_bytes = {e["size_bytes"]}')
        for k in ("fetched_at", "extracted_at"):
            if k in e:
                lines.append(f'{k} = "{e[k]}"')
        if "status" in e:
            lines.append(f'status = "{e["status"]}"')
        lines.append("")
    prov_path.write_text("\n".join(lines))


def main():
    if not TARBALL.is_file():
        print(f"RefDB tarball not found at {TARBALL}", file=sys.stderr)
        sys.exit(1)

    tarball_hash = tarball_sha256(TARBALL)
    print(f"RefDB tarball SHA256: {tarball_hash}", flush=True)
    now_iso = datetime.now(timezone.utc).isoformat()

    with tarfile.open(TARBALL, "r:gz") as tar:
        summary = []
        for pdb_id, bmrb_id in PROTEINS:
            tag, entry = extract_one(tar, pdb_id, bmrb_id, tarball_hash, now_iso)
            update_provenance(tag, entry)
            summary.append((tag, entry["status"], entry["size_bytes"]))
            print(f"[{tag}] {entry['filename']} {entry['status']} ({entry['size_bytes']} bytes)")

    print("\nSummary:")
    for tag, status, size in summary:
        print(f"  {tag:14s}  {status:16s}  {size:>8d} bytes")


if __name__ == "__main__":
    main()
