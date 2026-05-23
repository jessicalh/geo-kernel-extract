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
# LOG TAB: the viewer binds UDP port 9998 for its own log tab. With the
# canonical multicast `udp_host = "239.255.0.1"` in [logging] of
# ~/.nmr_tools.toml, both the viewer and ui/udp_listen.py can co-listen —
# each socket joins the multicast group and the kernel delivers every
# datagram to both. Run udp_listen.py in another terminal during a
# viewer session for a terminal-side tail of the same stream.
#
# Unicast hosts (e.g. 127.0.0.1) fall back to the old "last binder wins"
# behaviour — only one of the viewer's Log dock and udp_listen.py will
# see datagrams. The canonical TOML is multicast for this reason.
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
