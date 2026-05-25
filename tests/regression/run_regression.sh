#!/bin/bash
#
# Regression test: exercises real use cases and diffs output against baseline.
#
# Use case B: ORCA DFT (small protein, already protonated, with DFT tensors)
#
# Reads paths from testpaths.toml via the same convention as TestEnvironment.
# Requires nmr_extract to be built.
#
# Use Case D (--fleet against CHARMM/XTC fleet/) was retired 2026-05-04.
# The --fleet CLI flag was removed 2026-04-12; the fleet/ fixture and
# its loader (GromacsEnsembleLoader) moved to tests/bones/ on 2026-05-04.
#
# ⚠ MAINTENANCE NOTE (audit C10, 2026-05-25): the Use Case B invocation
# below uses the RETIRED `--orca DIR` shape. The canonical CLI is
# `--orca --root NAME` (NAME expands to {NAME}.xyz/.prmtop/_nmr.out), so
# nmr_extract rejects the positional-dir form and the test FAILs. Unlike
# batch_extract.sh's consolidated data, the orca_dir fixture DOES map to
# --root (e.g. {orca_dir}/A0A7C5FAR6_ALA), so this is fixable — but doing
# so needs a decision on which pose(s) to exercise (WT vs ALA), whether
# MOPAC runs (~10 min/protein), and a baseline_orca re-bless. Left for
# that decision; the read_toml path-resolution below is already correct.
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="$PROJECT_DIR/build"
TOML="$PROJECT_DIR/tests/testpaths.toml"
REGRESSION_DIR="$SCRIPT_DIR"

# Read a key from testpaths.toml, mirroring TestEnvironment's C10 path
# resolution: a value that does not start with '/' is repo-root-relative
# and resolved against PROJECT_DIR; absolute values are used verbatim.
read_toml() {
    local val
    val=$(grep "^$1" "$TOML" | sed 's/.*= *"\(.*\)"/\1/')
    if [ -n "$val" ] && [ "${val#/}" = "$val" ]; then
        val="$PROJECT_DIR/$val"
    fi
    printf '%s\n' "$val"
}

NMR_EXTRACT="$BUILD_DIR/nmr_extract"
if [ ! -x "$NMR_EXTRACT" ]; then
    echo "FAIL: nmr_extract not found at $NMR_EXTRACT — build first"
    exit 1
fi

ORCA_DIR=$(read_toml orca_dir)
PASS=0
FAIL=0

# ============================================================================
# Use case B: ORCA DFT
# ============================================================================

echo "=== Use case B: ORCA DFT ==="
ORCA_OUT="$REGRESSION_DIR/output_orca"
rm -rf "$ORCA_OUT"

if "$NMR_EXTRACT" --orca "$ORCA_DIR" --output "$ORCA_OUT" 2>/dev/null; then
    ORCA_COUNT=$(ls "$ORCA_OUT"/*.npy 2>/dev/null | wc -l)
    echo "  Produced $ORCA_COUNT arrays"
    if [ "$ORCA_COUNT" -ge 25 ]; then
        echo "  PASS: ORCA use case"
        PASS=$((PASS + 1))
    else
        echo "  FAIL: expected >= 25 arrays, got $ORCA_COUNT"
        FAIL=$((FAIL + 1))
    fi
else
    echo "  FAIL: nmr_extract --orca failed"
    FAIL=$((FAIL + 1))
fi

# If baseline exists, do byte-identical diff
ORCA_BASELINE="$REGRESSION_DIR/baseline_orca"
if [ -d "$ORCA_BASELINE" ]; then
    DIFFS=0
    for f in "$ORCA_BASELINE"/*.npy; do
        base=$(basename "$f")
        if [ -f "$ORCA_OUT/$base" ]; then
            if ! cmp -s "$f" "$ORCA_OUT/$base"; then
                echo "  DIFF: $base"
                DIFFS=$((DIFFS + 1))
            fi
        else
            echo "  MISSING: $base"
            DIFFS=$((DIFFS + 1))
        fi
    done
    if [ "$DIFFS" -eq 0 ]; then
        echo "  PASS: byte-identical to baseline"
        PASS=$((PASS + 1))
    else
        echo "  WARN: $DIFFS files differ from baseline"
    fi
else
    echo "  No baseline — run with --save-baseline to create one"
fi

# Use case D (CHARMM fleet) excerpted 2026-05-04 to
# tests/bones/run_regression_fleet_excerpt.sh.

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "=== Regression summary: $PASS passed, $FAIL failed ==="

# Save baseline if requested
if [ "${1:-}" = "--save-baseline" ]; then
    if [ -d "$ORCA_OUT" ]; then
        rm -rf "$ORCA_BASELINE"
        cp -r "$ORCA_OUT" "$ORCA_BASELINE"
        echo "Saved ORCA baseline to $ORCA_BASELINE"
    fi
fi

exit "$FAIL"
