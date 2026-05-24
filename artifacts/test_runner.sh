#!/bin/bash
# artifacts/test_runner.sh
#
# Comprehensive test: rebuild from current src/+tests/ state, run all
# ctests + Python SDK + h5-reader build + every CLI mode end-to-end.
# AIMNet2 is mandatory on every nmr_extract invocation. No symlinks.
#
# Logs to artifacts/test_run_<timestamp>/<phase>.log
# Continues on failure; final SUMMARY.log shows PASS/FAIL/elapsed per phase.

set -uo pipefail   # NOT -e — we want failures logged and continued

cd "$(git rev-parse --show-toplevel)"

TIMESTAMP=$(date +%Y%m%dT%H%M%S)
RUNDIR=artifacts/test_run_${TIMESTAMP}
mkdir -p "$RUNDIR"
SUMMARY=$RUNDIR/SUMMARY.log

AIMNET2_MODEL=/shared/2026Thesis/nmr-shielding/data/models/aimnet2_wb97m_0.jpt
CONFIG=data/calculator_params.toml
WRAPPER=scripts/run_with_cuda_env.sh
NMR_EXTRACT=./build/nmr_extract

# Common arg trailer: always include --aimnet2 and --config.
COMMON_ARGS=(--aimnet2 "$AIMNET2_MODEL" --config "$CONFIG")

# --- pre-flight: verify the things the run depends on ---
echo "=== test_runner.sh starting $(date -Is) ===" | tee "$SUMMARY"
echo "Run dir: $RUNDIR" | tee -a "$SUMMARY"
echo "" | tee -a "$SUMMARY"

PREFLIGHT_OK=1
for f in "$AIMNET2_MODEL" "$CONFIG" "$WRAPPER" \
         tests/data/illustrative_peptides/1l2y.pdb \
         tests/data/external/1UBQ.pdb \
         tests/data/1ubq_protonated.pdb \
         /mnt/expansion/larsen_archive/structures/1UBQ_pm6dh3plus.pdb \
         tests/data/orca/A0A7C5FAR6_WT.xyz \
         tests/data/orca/A0A7C5FAR6_WT.prmtop \
         tests/data/orca/A0A7C5FAR6_ALA.xyz \
         tests/data/orca/A0A7C5FAR6_ALA.prmtop \
         tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.tpr \
         tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.trr \
         tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.edr \
         tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/topol.top; do
  if [ ! -e "$f" ]; then
    echo "PREFLIGHT MISSING: $f" | tee -a "$SUMMARY"
    PREFLIGHT_OK=0
  fi
done
if [ $PREFLIGHT_OK -ne 1 ]; then
  echo "PREFLIGHT FAILED — aborting" | tee -a "$SUMMARY"
  exit 2
fi
python3 -c "import nvidia.cu13" 2>/dev/null || {
  echo "PREFLIGHT MISSING: python3 import nvidia.cu13 — AIMNet2 will fail at runtime" | tee -a "$SUMMARY"
  exit 2
}
echo "preflight OK ($(date -Is))" | tee -a "$SUMMARY"
echo "" | tee -a "$SUMMARY"

# --- run_phase: log, time, continue-on-failure ---
run_phase() {
  local name=$1; shift
  local log=$RUNDIR/${name}.log
  local t0; t0=$(date +%s)
  echo "===== [$name] starting $(date -Is)" | tee -a "$SUMMARY"
  ("$@") >"$log" 2>&1
  local rc=$?
  local t1; t1=$(date +%s)
  local elapsed=$((t1-t0))
  local hms; hms=$(printf "%02dh%02dm%02ds" $((elapsed/3600)) $((elapsed%3600/60)) $((elapsed%60)))
  if [ $rc -eq 0 ]; then
    echo "      [$name] PASS in $hms" | tee -a "$SUMMARY"
  else
    echo "      [$name] FAIL rc=$rc in $hms (see $log)" | tee -a "$SUMMARY"
  fi
  return 0
}

# --- phase 01: clean rebuild ---
run_phase 01_clean_build bash -c '
  set -euo pipefail
  rm -rf build
  cmake --preset default
  cmake --build build -j"$(nproc)"
'

# --- phase 02: ctest (711 tests, env via gtest_discover_tests PROPERTIES) ---
run_phase 02_ctest bash -c '
  cd build
  ctest --output-on-failure
'

# --- phase 03: Python SDK ---
run_phase 03_python_sdk bash -c '
  set -euo pipefail
  cd python
  python3 -m pytest tests/ -v
'

# --- phase 04: h5-reader build ---
run_phase 04_h5reader_build bash -c '
  set -euo pipefail
  if [ ! -d h5-reader/build/linux-release ]; then
    cmake -S h5-reader -B h5-reader/build/linux-release -DCMAKE_BUILD_TYPE=Release
  fi
  cmake --build h5-reader/build/linux-release -j"$(nproc)"
'

# --- phase 05: Mode 1 PDB ---
run_phase 05a_mode1_pdb_1l2y_mopac \
  "$WRAPPER" "$NMR_EXTRACT" \
    --pdb tests/data/illustrative_peptides/1l2y.pdb \
    "${COMMON_ARGS[@]}" \
    --output "$RUNDIR/mode1_1l2y_mopac"

run_phase 05b_mode1_pdb_1l2y_nomopac \
  "$WRAPPER" "$NMR_EXTRACT" \
    --pdb tests/data/illustrative_peptides/1l2y.pdb --no-mopac \
    "${COMMON_ARGS[@]}" \
    --output "$RUNDIR/mode1_1l2y_nomopac"

run_phase 05c_mode1_pdb_1ubq_mopac \
  "$WRAPPER" "$NMR_EXTRACT" \
    --pdb tests/data/external/1UBQ.pdb \
    "${COMMON_ARGS[@]}" \
    --output "$RUNDIR/mode1_1ubq_mopac"

# --- phase 06: Mode 2 Protonated PDB ---
run_phase 06a_mode2_prot_1ubq_mopac \
  "$WRAPPER" "$NMR_EXTRACT" \
    --protonated-pdb tests/data/1ubq_protonated.pdb \
    "${COMMON_ARGS[@]}" \
    --output "$RUNDIR/mode2_1ubq_mopac"

run_phase 06b_mode2_prot_1UBQ_pm6dh3plus_mopac \
  "$WRAPPER" "$NMR_EXTRACT" \
    --protonated-pdb /mnt/expansion/larsen_archive/structures/1UBQ_pm6dh3plus.pdb \
    "${COMMON_ARGS[@]}" \
    --output "$RUNDIR/mode2_pm6dh3plus_mopac"

# --- phase 07: Mode 3 ORCA single root, WT and ALA, with and without MOPAC ---
run_phase 07a_mode3_orca_WT_mopac \
  "$WRAPPER" "$NMR_EXTRACT" \
    --orca --root tests/data/orca/A0A7C5FAR6_WT \
    "${COMMON_ARGS[@]}" \
    --output "$RUNDIR/mode3_WT_mopac"

run_phase 07b_mode3_orca_WT_nomopac \
  "$WRAPPER" "$NMR_EXTRACT" \
    --orca --root tests/data/orca/A0A7C5FAR6_WT --no-mopac \
    "${COMMON_ARGS[@]}" \
    --output "$RUNDIR/mode3_WT_nomopac"

run_phase 07c_mode3_orca_ALA_mopac \
  "$WRAPPER" "$NMR_EXTRACT" \
    --orca --root tests/data/orca/A0A7C5FAR6_ALA \
    "${COMMON_ARGS[@]}" \
    --output "$RUNDIR/mode3_ALA_mopac"

run_phase 07d_mode3_orca_ALA_nomopac \
  "$WRAPPER" "$NMR_EXTRACT" \
    --orca --root tests/data/orca/A0A7C5FAR6_ALA --no-mopac \
    "${COMMON_ARGS[@]}" \
    --output "$RUNDIR/mode3_ALA_nomopac"

# --- phase 08: Mode 4 Mutant pair, with and without MOPAC ---
run_phase 08a_mode4_mutant_mopac \
  "$WRAPPER" "$NMR_EXTRACT" \
    --mutant \
      --wt tests/data/orca/A0A7C5FAR6_WT \
      --ala tests/data/orca/A0A7C5FAR6_ALA \
    "${COMMON_ARGS[@]}" \
    --output "$RUNDIR/mode4_mutant_mopac"

run_phase 08b_mode4_mutant_nomopac \
  "$WRAPPER" "$NMR_EXTRACT" \
    --mutant \
      --wt tests/data/orca/A0A7C5FAR6_WT \
      --ala tests/data/orca/A0A7C5FAR6_ALA \
      --no-mopac \
    "${COMMON_ARGS[@]}" \
    --output "$RUNDIR/mode4_mutant_nomopac"

# --- phase 09: Mode 5 Trajectory. Real production fixture (no fixture setup;
# the in-tree 1P9J batcave_local_15ns trajectory is point-at-directly,
# parent-dir topol.top is picked up by convention). MOPAC at ~10 min/frame
# is intractable on the full 750-frame trajectory (~125 h); MOPAC coverage
# stays in the single-conformation phases (05-08). ---
MODE5_TRAJ=tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z
MODE5_OUTPUT=$RUNDIR/mode5_output
MODE5_PDBS=$RUNDIR/mode5_frame_pdbs
MODE5_NPYS=$RUNDIR/mode5_frame_npys

run_phase 09a_mode5_trajectory_aimnet2 \
  "$WRAPPER" "$NMR_EXTRACT" \
    --trajectory "$MODE5_TRAJ" \
    "${COMMON_ARGS[@]}" \
    --emit-frame-pdbs "$MODE5_PDBS" --pdb-stride 100 \
    --output "$MODE5_OUTPUT"

# 09b — second run, NPY-stride emitter coverage. High stride keeps disk
# bounded (one frame_NNNNNN/ dir at ~6 MB + the per-protein sidecars).
run_phase 09b_mode5_trajectory_npy_stride \
  "$WRAPPER" "$NMR_EXTRACT" \
    --trajectory "$MODE5_TRAJ" \
    "${COMMON_ARGS[@]}" \
    --emit-frame-npys "$MODE5_NPYS" --npy-stride 9999 \
    --output "$RUNDIR/mode5_output_npy_run"

# --- phase 10: verify outputs ---
run_phase 10_verify_outputs bash -c "
  echo '=== NPY count per mode output dir ==='
  for d in '$RUNDIR'/mode*/; do
    [ -d \"\$d\" ] || continue
    npy_count=\$(find \"\$d\" -maxdepth 2 -name '*.npy' 2>/dev/null | wc -l)
    echo \"  \$d : \$npy_count npy files\"
  done
  echo ''
  echo '=== Trajectory output (mode 5) ==='
  ls -la '$MODE5_OUTPUT/' 2>/dev/null | head -30
  echo ''
  echo '=== Per-frame PDBs (mode 5) ==='
  pdb_count=\$(find '$MODE5_PDBS' -maxdepth 2 -name '*.pdb' 2>/dev/null | wc -l)
  echo \"  $MODE5_PDBS : \$pdb_count per-frame PDBs\"
  echo ''
  echo '=== Per-frame NPYs (mode 5) ==='
  ls -la '$MODE5_NPYS/' 2>/dev/null | head -15
  npy_dir_count=\$(find '$MODE5_NPYS' -maxdepth 1 -type d -name 'frame_*' 2>/dev/null | wc -l)
  echo \"  $MODE5_NPYS : \$npy_dir_count frame_NNNNNN/ dirs\"
  echo ''
  echo '=== trajectory.h5 size ==='
  ls -la '$MODE5_OUTPUT/trajectory.h5' 2>/dev/null || echo '  (no trajectory.h5)'
"

# --- final summary ---
echo "" | tee -a "$SUMMARY"
echo "===== FINAL SUMMARY at $(date -Is)" | tee -a "$SUMMARY"
echo "" | tee -a "$SUMMARY"
grep -E "PASS|FAIL|PREFLIGHT" "$SUMMARY" | grep -vE "^=====|^$"
echo "" | tee -a "$SUMMARY"
echo "All logs: $RUNDIR/" | tee -a "$SUMMARY"
echo "===== DONE $(date -Is) ====="
