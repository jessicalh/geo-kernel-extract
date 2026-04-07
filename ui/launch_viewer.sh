#!/bin/bash
# Launch the viewer with UDP logging.
#
# Usage:
#   ui/launch_viewer.sh                                        # 1ubq protonated
#   ui/launch_viewer.sh --pdb path/to/protein.pdb              # bare PDB
#   ui/launch_viewer.sh --protonated-pdb path/to/protein.pdb   # pre-protonated
#   ui/launch_viewer.sh --orca --root path/to/A0A7C5FAR6_WT    # ORCA DFT
#   ui/launch_viewer.sh --mutant --wt path/to/WT --ala path/to/ALA
#   ui/launch_viewer.sh --fleet --tpr path/to/prod.tpr --poses path/to/poses
#
# REST commands (port 9147):
#   echo '{"cmd":"status"}' | nc -q1 localhost 9147
#   echo '{"cmd":"screenshot","path":"/tmp/shot.png"}' | nc -q1 localhost 9147
#   echo '{"cmd":"export_features","path":"/tmp/features"}' | nc -q1 localhost 9147
#   echo '{"cmd":"reset_view"}' | nc -q1 localhost 9147

set -e
cd "$(dirname "$0")/.."

# Start UDP listener in background
python3 ui/udp_listen.py &
UDP_PID=$!
trap "kill $UDP_PID 2>/dev/null" EXIT

# Default: protonated 1ubq. Otherwise pass all args through to JobSpec.
if [ $# -eq 0 ]; then
    DISPLAY=:1 build-ui/nmr-viewer --protonated-pdb tests/data/1ubq_protonated.pdb 2>&1
else
    DISPLAY=:1 build-ui/nmr-viewer "$@" 2>&1
fi
