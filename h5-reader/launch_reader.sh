#!/usr/bin/env bash
# launch_reader.sh — build and run h5reader with sensible defaults.
#
# Usage:
#   ./launch_reader.sh                             # build + run, no H5
#   ./launch_reader.sh path/to/protein_analysis.h5 # build + run with H5

set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

PRESET="${H5READER_PRESET:-linux-gcc}"

# Configure if the build directory doesn't exist yet.
if [[ ! -d "build/${PRESET}" ]]; then
    cmake --preset "${PRESET}"
fi

# Build.
cmake --build --preset "${PRESET}" -j

BIN="build/${PRESET}/h5reader"
if [[ ! -x "${BIN}" ]]; then
    echo "h5reader binary missing after build: ${BIN}" >&2
    exit 1
fi

# Forward all arguments to the binary.
exec "${BIN}" "$@"
