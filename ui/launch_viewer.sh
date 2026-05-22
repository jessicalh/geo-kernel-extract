#!/bin/bash
# Launch the viewer.
#
# Usage:
#   ui/launch_viewer.sh                                        # 1ubq protonated
#   ui/launch_viewer.sh --pdb path/to/protein.pdb              # bare PDB
#   ui/launch_viewer.sh --protonated-pdb path/to/protein.pdb   # pre-protonated
#   ui/launch_viewer.sh --orca --root path/to/A0A7C5FAR6_WT    # ORCA DFT
#   ui/launch_viewer.sh --mutant --wt path/to/WT --ala path/to/ALA
#   ui/launch_viewer.sh --analysis-h5 path/to/{X}_analysis.h5   # ns0 pose auto-derived
#
# AIMNet2 is a stone requirement on every production path
# (feedback_aimnet2_required_no_weasel). Resolution order is
# CLI --aimnet2 PATH  >  calculator_params.toml [aimnet2_model_path]  >  error.
# The project default lives at data/calculator_params.toml and is
# loaded by main_viewer.cpp between ParseJobSpec and ValidateJobSpec,
# so no flag is required for the canonical model.
#
# CUDA: launches go through scripts/run_with_cuda_env.sh, which
# prepends PyTorch's bundled nvidia.cu13 lib to LD_LIBRARY_PATH so
# libnvrtc.so.13's dlopen of libnvrtc-builtins.so.13.0 succeeds at
# AIMNet2 JIT time.
#
# Viewer always skips MOPAC and Coulomb (batch/calibration paths, not interactive).
# APBS runs by default. Pass --no-apbs to skip.
#
# LOG TAB: The viewer binds UDP port 9998 for its own log tab. Do NOT run
# udp_listen.py alongside the viewer — both would bind 9998 and Linux unicast
# UDP delivers each datagram to exactly one socket (the last to bind). If
# udp_listen.py binds first, the viewer log tab will be silent.
#
# udp_listen.py is for batch/CLI sessions without an open viewer. It is a
# "safe" listener (SO_REUSEADDR + SO_REUSEPORT) and will not block the viewer
# from binding, but it will not receive viewer-session datagrams if the viewer
# bound 9998 after it.
#
# REST commands (port 9147):
#   echo '{"cmd":"status"}' | nc -q1 localhost 9147
#   echo '{"cmd":"screenshot","path":"/tmp/shot.png"}' | nc -q1 localhost 9147
#   echo '{"cmd":"export_features","path":"/tmp/features"}' | nc -q1 localhost 9147
#   echo '{"cmd":"reset_view"}' | nc -q1 localhost 9147

set -e
cd "$(dirname "$0")/.."

# Default: protonated 1ubq. Otherwise pass all args through to JobSpec.
if [ $# -eq 0 ]; then
    DISPLAY=:1 scripts/run_with_cuda_env.sh build/nmr-viewer \
        --protonated-pdb tests/data/1ubq_protonated.pdb 2>&1
else
    DISPLAY=:1 scripts/run_with_cuda_env.sh build/nmr-viewer "$@" 2>&1
fi
