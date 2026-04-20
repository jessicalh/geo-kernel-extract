#!/usr/bin/env python3
"""
fetch_bmrb.py — fetch BMRB NMR-STAR files for the 10 calibration proteins.

Scope (2026-04-20): ONLY the 10 calibration proteins. Fleet is deferred.
See memory:feedback_nmr_scope_10_first.

Per protein, fetches `bmr<id>_3.str` (NMR-STAR 3.x format) from
`https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr<id>/` and stores under
`source_data/<pdb>_<bmrb>/raw/`. SHA256 and fetch metadata are recorded in
`source_data/<pdb>_<bmrb>/raw_provenance.toml`.

Run from nmr_forensics/ in the venv:
    ./.venv/bin/python fetch_bmrb.py

Idempotent: re-runs detect unchanged content by hash and leave it.
"""

import hashlib
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

HERE = Path(__file__).parent
SOURCE_DATA = HERE / "source_data"

# Scope: the 10 calibration proteins.
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

URL_TEMPLATE = "https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{bmrb_id}/bmr{bmrb_id}_3.str"


def fetch_one(pdb_id, bmrb_id):
    tag = f"{pdb_id}_{bmrb_id}"
    protein_dir = SOURCE_DATA / tag
    raw_dir = protein_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    filename = f"bmr{bmrb_id}_3.str"
    out_path = raw_dir / filename
    url = URL_TEMPLATE.format(bmrb_id=bmrb_id)

    print(f"[{tag}] fetching {url}", flush=True)
    t0 = time.time()
    resp = requests.get(url, timeout=60)
    elapsed_s = time.time() - t0
    resp.raise_for_status()

    content = resp.content
    new_hash = hashlib.sha256(content).hexdigest()

    if out_path.exists():
        old_hash = hashlib.sha256(out_path.read_bytes()).hexdigest()
        if old_hash == new_hash:
            print(f"[{tag}] unchanged (sha256 {new_hash[:16]}...)", flush=True)
            return tag, url, new_hash, len(content), elapsed_s, "unchanged"
        else:
            print(
                f"[{tag}] CONTENT CHANGED at BMRB "
                f"(old {old_hash[:16]}... → new {new_hash[:16]}...)",
                flush=True,
            )
            status = "content_changed"
    else:
        status = "new"

    out_path.write_bytes(content)
    return tag, url, new_hash, len(content), elapsed_s, status


def update_provenance(tag, entries):
    """
    entries: list of dicts with keys:
        kind, filename, url, sha256, size_bytes, fetched_at, status
    """
    prov_path = SOURCE_DATA / tag / "raw_provenance.toml"
    lines = [
        "# Raw-source provenance for per-protein experimental-shift inputs.",
        "# Auto-written by fetch_bmrb.py / extract_refdb.py. Do not hand-edit;",
        "# re-run the fetch scripts to refresh.",
        "",
        f'protein_tag = "{tag}"',
        "",
    ]
    for e in entries:
        lines.append(f'[[sources]]')
        lines.append(f'kind = "{e["kind"]}"')
        lines.append(f'filename = "{e["filename"]}"')
        lines.append(f'url = "{e["url"]}"')
        lines.append(f'sha256 = "{e["sha256"]}"')
        lines.append(f'size_bytes = {e["size_bytes"]}')
        lines.append(f'fetched_at = "{e["fetched_at"]}"')
        lines.append(f'status = "{e["status"]}"')
        lines.append("")
    prov_path.write_text("\n".join(lines))


def main():
    now_iso = datetime.now(timezone.utc).isoformat()
    summary = []
    for pdb_id, bmrb_id in PROTEINS:
        try:
            tag, url, sha, size, elapsed, status = fetch_one(pdb_id, bmrb_id)
        except Exception as exc:
            print(f"[{pdb_id}_{bmrb_id}] FETCH FAILED: {exc}", flush=True)
            sys.exit(1)
        entry = {
            "kind": "bmrb_nmr_star_3",
            "filename": f"bmr{bmrb_id}_3.str",
            "url": url,
            "sha256": sha,
            "size_bytes": size,
            "fetched_at": now_iso,
            "status": status,
        }
        # Preserve any existing RefDB provenance entries by reading/appending.
        prov_path = SOURCE_DATA / tag / "raw_provenance.toml"
        existing_refdb = []
        if prov_path.exists():
            # Crude parse: preserve any `[[sources]]` block whose kind != bmrb_nmr_star_3.
            text = prov_path.read_text()
            import tomllib
            existing = tomllib.loads(text)
            for s in existing.get("sources", []):
                if s.get("kind") != "bmrb_nmr_star_3":
                    existing_refdb.append(s)
        update_provenance(tag, [entry] + existing_refdb)
        summary.append((tag, status, size, elapsed))

    print("\nSummary:")
    for tag, status, size, elapsed in summary:
        print(f"  {tag:14s}  {status:16s}  {size:>8d} bytes  {elapsed:6.2f}s")


if __name__ == "__main__":
    main()
