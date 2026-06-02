#!/usr/bin/env bash
# Export the local tensorcs15 PostgreSQL replica into the Docker vendor payload.
# The DSN is read from NMR_TENSORCS15_DSN or from the runtime TOML and is never
# printed.

set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)

OUT=${1:-"${ROOT_DIR}/docker/vendor/tensorcs15/tensorcs15.dump"}
CONFIG=${NMR_TOOLS_TOML:-"${HOME}/.nmr_tools.toml"}
TABLE=${NMR_TENSORCS15_TABLE:-raw_dft_calculations}

fail() {
    printf 'export-tensorcs15-dump: %s\n' "$*" >&2
    exit 2
}

need_cmd() {
    command -v "$1" >/dev/null 2>&1 || fail "$1 not found"
}

load_dsn() {
    if [ -n "${NMR_TENSORCS15_DSN:-}" ]; then
        printf '%s\n' "$NMR_TENSORCS15_DSN"
        return 0
    fi
    python3 - "$CONFIG" <<'PY_DSN'
import pathlib
import sys

try:
    import tomllib
except ModuleNotFoundError:
    sys.exit(1)

path = pathlib.Path(sys.argv[1])
if not path.is_file():
    sys.exit(1)

try:
    data = tomllib.loads(path.read_text())
    dsn = data.get("databases", {}).get("tensorcs15", "")
except Exception:
    sys.exit(1)

if dsn:
    print(dsn)
else:
    sys.exit(1)
PY_DSN
}

need_cmd python3
need_cmd pg_dump
need_cmd pg_restore

if [ ! -f "${ROOT_DIR}/config/tensorcs15.manifest.toml" ]; then
    fail "missing tensorcs15 manifest: ${ROOT_DIR}/config/tensorcs15.manifest.toml"
fi

dsn=$(load_dsn) || fail "tensorcs15 DSN not configured in NMR_TENSORCS15_DSN or $CONFIG"

mkdir -p "$(dirname -- "$OUT")"
tmp="${OUT}.tmp.$$"
trap 'rm -f "$tmp"' EXIT

printf 'export-tensorcs15-dump: checking live tensorcs15 manifest\n' >&2
"${ROOT_DIR}/scripts/tensorcs15-check.py" --quiet

printf 'export-tensorcs15-dump: writing %s\n' "$OUT" >&2
pg_dump \
    --format=custom \
    --compress=9 \
    --no-owner \
    --no-privileges \
    --table="$TABLE" \
    --file="$tmp" \
    "$dsn"

pg_restore --list "$tmp" >/dev/null
mv "$tmp" "$OUT"
trap - EXIT

printf 'export-tensorcs15-dump: wrote %s\n' "$OUT" >&2
sha256sum "$OUT"
