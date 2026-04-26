#!/usr/bin/env bash
#
# fetch_ccd.sh — refresh the Chemical Component Dictionary used by
# CcdValidator. Downloads the gzipped CCD from RCSB and decompresses it
# in place. Output: data/ccd/components.cif (~376 MB raw).
#
# CcdValidator uses gemmi's header-only cif::read_file which reads
# plain CIF — that's why we keep the file uncompressed on disk. Disk
# is cheap, ABI coupling to libgemmi_cpp.so was not.
#
# The CCD is read-only at runtime, used at protein load time to enrich
# IUPAC topology diagnostics. Standard 20 amino-acid entries have not
# changed in decades; refresh primarily benefits non-standard residues.
#
# Source: PDB Chemical Component Dictionary at RCSB.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="$ROOT/data/ccd"
GZ_FILE="$OUT_DIR/components.cif.gz"
OUT_FILE="$OUT_DIR/components.cif"
URL="https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"

mkdir -p "$OUT_DIR"

echo "Fetching $URL"
echo "  -> $GZ_FILE"

if ! curl --fail --location --show-error --silent --output "$GZ_FILE.tmp" "$URL"; then
    echo "ERROR: download failed" >&2
    rm -f "$GZ_FILE.tmp"
    exit 1
fi

mv "$GZ_FILE.tmp" "$GZ_FILE"

echo "Decompressing -> $OUT_FILE"
gunzip --force --keep "$GZ_FILE"

GZ_SIZE=$(du -h "$GZ_FILE" | cut -f1)
RAW_SIZE=$(du -h "$OUT_FILE" | cut -f1)
echo "Done. components.cif.gz=$GZ_SIZE  components.cif=$RAW_SIZE"
