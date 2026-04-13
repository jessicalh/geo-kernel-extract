#!/usr/bin/env bash
# Run all Stage 1 analyses.  All output goes to src/output/actual_physics/.
# Reproduces every number in the stage1-mutations/notes/ documents.
#
# Prerequisites:
#   - nmr_extract built (build/nmr_extract)
#   - Stage1Results extraction complete (723 proteins)
#   - calibration.toml pointing at Stage1Results
#
# Usage:
#   cd learn/src
#   bash ../stage1-mutations/analysis/run_all.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$(dirname "$SCRIPT_DIR")/../src"

echo "Working directory: $(pwd)"
echo "Config: calibration.toml"
echo ""

run_analysis() {
    local name="$1"
    local module="$2"
    echo "━━━ $name ━━━"
    python3 -u -c "
from mutation_set.config import load_config
import sys; cfg = load_config('calibration.toml')
sys.path.insert(0, str(cfg.paths.sdk))
from actual_physics.${module} import run
run(cfg)
"
    echo ""
}

# 1. Full kernel space: cosine matrices, per-group R², forward selection
#    Output: output/actual_physics/full_space/
run_analysis "Full space analysis" "full_space_analysis"

# 2. Geometry-only basis: what survives without MOPAC, ff14SB vs MOPAC EFG
#    Output: output/actual_physics/geometry_only_basis/
run_analysis "Geometry-only basis" "geometry_only_basis"

# 3. Dimensionality: PCA-ridge curves, distance bands, raw vs normalised
#    Output: output/actual_physics/dimensionality/
run_analysis "Dimensionality test" "dimensionality_test"

# 4. Accumulation basis: ridge coefficient structure, PCA projection
#    Output: output/actual_physics/accumulation_basis/
run_analysis "Accumulation basis" "accumulation_basis"

# 5. Dia/para decomposition: parse raw orca, compute delta dia/para
#    Output: output/actual_physics/dia_para/
run_analysis "Dia/para decomposition" "orca_dia_para"

# 6. Element physics: per-element kernel decomposition, forward selection
#    Output: output/secondary/element_physics/
run_analysis "Element physics" "element_physics"

# 7. Twenty realities: analytical physics tests
#    Output: output/secondary/twenty_realities/
run_analysis "Twenty realities" "twenty_realities"

# 8. Tensor realities: tensor-direct physics tests
#    Output: output/secondary/tensor_realities/
run_analysis "Tensor realities" "tensor_realities"

# 9. Export per-atom data for R figures
#    Output: output/actual_physics/atom_tensor_data.csv
run_analysis "Export for R" "export_for_r"

# 10. Physics calibration: progressive fair-scalar table
#     Output: output/actual_physics/physics_calibration/
run_analysis "Physics calibration" "physics_calibration"

# 11. Clean calibration: 55-kernel per-element ridge
#     Output: output/actual_physics/clean_calibration/
run_analysis "Clean calibration" "clean_calibration"

echo "━━━ All analyses complete ━━━"
echo "Output in: src/output/actual_physics/ and src/output/secondary/"
echo "JSON results can be read by R scripts in stage1-mutations/analysis/"
