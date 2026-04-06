#!/bin/bash
# Launch the viewer with UDP logging and a test protein.
#
# Usage:
#   ui/launch_viewer.sh                          # 1ubq_protonated.pdb
#   ui/launch_viewer.sh path/to/protein.pdb      # specific PDB
#   ui/launch_viewer.sh path/to/protein_dir/     # protein directory (ORCA comparison)
#
# REST commands (port 9147):
#   echo '{"cmd":"status"}' | nc -q1 localhost 9147
#   echo '{"cmd":"screenshot","path":"/tmp/shot.png"}' | nc -q1 localhost 9147
#   echo '{"cmd":"reset_view"}' | nc -q1 localhost 9147
#   echo '{"cmd":"load_pdb","path":"/path/to/file.pdb"}' | nc -q1 localhost 9147
#   echo '{"cmd":"set_overlay","mode":"classical"}' | nc -q1 localhost 9147
#   echo '{"cmd":"set_calculators","bs":true,"hm":false}' | nc -q1 localhost 9147

set -e
cd "$(dirname "$0")/.."

# Start UDP listener in background
python3 ui/udp_listen.py &
UDP_PID=$!
trap "kill $UDP_PID 2>/dev/null" EXIT

# Default test protein
INPUT="${1:-tests/data/1ubq_protonated.pdb}"

if [ -d "$INPUT" ]; then
    DISPLAY=:1 build-ui/nmr-viewer --protein "$INPUT" 2>&1
else
    DISPLAY=:1 build-ui/nmr-viewer --pdb "$INPUT" 2>&1
fi
