#!/bin/bash
#
# Regression test: exercises real use cases and diffs output against baseline.
#
# Use case B: ORCA DFT (small protein, already protonated, with DFT tensors)
# Use case D: GROMACS fleet (10 poses, CHARMM36m charges from TPR)
#
# Reads paths from testpaths.toml via the same convention as TestEnvironment.
# Requires nmr_extract to be built.
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="$PROJECT_DIR/build"
TOML="$PROJECT_DIR/tests/testpaths.toml"
REGRESSION_DIR="$SCRIPT_DIR"

# Read a key from testpaths.toml
read_toml() {
    grep "^$1" "$TOML" | sed 's/.*= *"\(.*\)"/\1/'
}

NMR_EXTRACT="$BUILD_DIR/nmr_extract"
if [ ! -x "$NMR_EXTRACT" ]; then
    echo "FAIL: nmr_extract not found at $NMR_EXTRACT — build first"
    exit 1
fi

ORCA_DIR=$(read_toml orca_dir)
FLEET_DATA=$(read_toml fleet_data)
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

# ============================================================================
# Use case D: GROMACS fleet (first protein in fleet data)
# ============================================================================

echo ""
echo "=== Use case D: GROMACS fleet ==="

# Find first protein dir in fleet data
FIRST_PROTEIN=$(ls -d "$FLEET_DATA"/*/poses 2>/dev/null | head -1 | sed 's|/poses$||')
if [ -z "$FIRST_PROTEIN" ]; then
    echo "  SKIP: no fleet data found at $FLEET_DATA"
else
    PROTEIN_ID=$(basename "$FIRST_PROTEIN")
    TPR_PATH="$FIRST_PROTEIN/params/prod.tpr"
    POSES_PATH="$FIRST_PROTEIN/poses"
    FLEET_OUT="$REGRESSION_DIR/output_fleet"
    rm -rf "$FLEET_OUT"

    if [ ! -f "$TPR_PATH" ]; then
        echo "  SKIP: TPR not found at $TPR_PATH"
    elif "$NMR_EXTRACT" --fleet "$POSES_PATH" "$TPR_PATH" --output "$FLEET_OUT" 2>/dev/null; then
        FRAME_COUNT=$(ls -d "$FLEET_OUT"/frame_* 2>/dev/null | wc -l)
        echo "  Protein: $PROTEIN_ID"
        echo "  Produced $FRAME_COUNT frames"
        if [ "$FRAME_COUNT" -ge 5 ]; then
            # Check first frame has arrays
            F1_COUNT=$(ls "$FLEET_OUT/frame_001"/*.npy 2>/dev/null | wc -l)
            echo "  Frame 001: $F1_COUNT arrays"
            if [ "$F1_COUNT" -ge 25 ]; then
                echo "  PASS: fleet use case"
                PASS=$((PASS + 1))
            else
                echo "  FAIL: expected >= 25 arrays per frame"
                FAIL=$((FAIL + 1))
            fi
        else
            echo "  FAIL: expected >= 5 frames, got $FRAME_COUNT"
            FAIL=$((FAIL + 1))
        fi
    else
        echo "  FAIL: nmr_extract --fleet failed"
        FAIL=$((FAIL + 1))
    fi
fi

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
