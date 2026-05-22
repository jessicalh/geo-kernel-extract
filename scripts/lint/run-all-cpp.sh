#!/usr/bin/env bash
# scripts/lint/run-all-cpp.sh
#
# Run all four C++ static analysis tools in sequence:
#   1. clang-tidy
#   2. cppcheck
#   3. include-what-you-use
#   4. scan-build  (slowest; skip with --skip-scan-build)
#
# Each tool emits its own JSON to artifacts/quality/<tool>-<sha>.json.
# Use scripts/lint/baseline.py to derive per-file finding counts for
# the ratchet.
#
# Usage:
#   scripts/lint/run-all-cpp.sh
#   scripts/lint/run-all-cpp.sh --skip-scan-build
#
# Total elapsed: ~10 min for the first three on this codebase, +30 min
# for scan-build.

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

SKIP_SCAN=0
for arg in "$@"; do
    case "$arg" in
        --skip-scan-build) SKIP_SCAN=1 ;;
        *) echo "Unknown arg: $arg" >&2; exit 2 ;;
    esac
done

# Always re-emit compile_commands.json before sweeps (cheap if up-to-date).
( cd build && cmake . >/dev/null 2>&1 ) || true

echo "===== [1/4] clang-tidy ====="
./scripts/lint/run-clang-tidy.sh

echo ""
echo "===== [2/4] cppcheck ====="
./scripts/lint/run-cppcheck.sh

echo ""
echo "===== [3/4] include-what-you-use ====="
./scripts/lint/run-iwyu.sh

if [[ $SKIP_SCAN -eq 0 ]]; then
    echo ""
    echo "===== [4/4] scan-build ====="
    ./scripts/lint/run-scan-build.sh
else
    echo ""
    echo "===== [4/4] scan-build SKIPPED ====="
fi

echo ""
echo "All sweeps complete. Reports in artifacts/quality/"
ls -lh artifacts/quality/ | grep -E "$(git rev-parse --short HEAD)" || true
