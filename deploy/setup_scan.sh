#!/bin/bash
#
# setup_scan.sh — Build nmr_extract on a scan machine (Ubuntu 22.04).
#
# Prerequisites already present:
#   - Eigen3 (apt)
#   - GROMACS 2026 at /opt/gromacs/2026.0/ + source at ~/build/gmx2026/
#   - libmopac at ~/micromamba/envs/mm/lib/
#   - cmake, g++, make, git
#
# This script installs missing apt packages, builds reduce from source,
# clones nmr-shielding, and builds nmr_extract.
#
# Usage:
#   ssh scan1
#   bash /path/to/setup_scan.sh
#
# After running, nmr_extract is at ~/nmr-shielding/build/nmr_extract
#

set -euo pipefail

echo "=== Installing apt packages ==="
sudo apt-get update -qq
sudo apt-get install -y -qq \
    libopenbabel-dev \
    libdssp-dev libcifpp-dev \
    libapbs-dev libfetk-dev

echo "=== Building reduce from source ==="
mkdir -p ~/builds
if [ ! -f ~/builds/reduce-src/build/reduce_src/libreducelib.a ]; then
    cd ~/builds
    if [ ! -d reduce-src ]; then
        git clone https://github.com/rlabduke/reduce.git reduce-src
    fi
    cd reduce-src
    mkdir -p build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
    make -j$(nproc) reducelib 2>&1 | tail -5
    echo "reduce built: $(ls -la reduce_src/libreducelib.a)"
else
    echo "reduce already built"
fi

echo "=== Cloning nmr-shielding ==="
cd ~
if [ ! -d nmr-shielding ]; then
    git clone /shared/2026Thesis/nmr-shielding nmr-shielding 2>/dev/null \
        || git clone git@github.com:USER/nmr-shielding.git nmr-shielding
fi
cd nmr-shielding
git pull --ff-only 2>/dev/null || true

echo "=== Configuring CMake ==="
mkdir -p build && cd build

# Paths for this machine
GROMACS_LIB="/opt/gromacs/2026.0/lib/libgromacs_mpi.so"
GROMACS_SRC="$HOME/build/gmx2026/gromacs-2026.0/src"
GROMACS_BUILD="$HOME/build/gmx2026/gromacs-2026-build"
MOPAC_INCLUDE="$HOME/micromamba/envs/mm/include"
MOPAC_LIB="$HOME/micromamba/envs/mm/lib/libmopac.so"
REDUCE_SRC="$HOME/builds/reduce-src"
REDUCE_LIB="$REDUCE_SRC/build/reduce_src/libreducelib.a"
REDUCE_LIBPDB="$REDUCE_SRC/build/libpdb/libpdb++.a"
REDUCE_TOOLCLASSES="$REDUCE_SRC/build/toolclasses/libtoolclasses.a"

# OpenBabel on Ubuntu 22.04
if [ -f /usr/lib/x86_64-linux-gnu/libopenbabel.so ]; then
    OB_LIB="/usr/lib/x86_64-linux-gnu/libopenbabel.so"
elif [ -f /usr/lib/x86_64-linux-gnu/libopenbabel.so.7 ]; then
    OB_LIB="/usr/lib/x86_64-linux-gnu/libopenbabel.so.7"
else
    echo "ERROR: libopenbabel not found"
    exit 1
fi

cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DGROMACS_LIB="$GROMACS_LIB" \
    -DGROMACS_SRC="$GROMACS_SRC" \
    -DGROMACS_BUILD="$GROMACS_BUILD" \
    -DMOPAC_INCLUDE="$MOPAC_INCLUDE" \
    -DMOPAC_LIB="$MOPAC_LIB" \
    -DOPENBABEL_LIB="$OB_LIB" \
    -DREDUCE_SRC="$REDUCE_SRC" \
    -DREDUCE_LIB="$REDUCE_LIB" \
    -DREDUCE_LIBPDB="$REDUCE_LIBPDB" \
    -DREDUCE_TOOLCLASSES="$REDUCE_TOOLCLASSES" \
    2>&1 | tail -20

echo "=== Building nmr_extract ==="
make -j$(nproc) nmr_extract 2>&1 | tail -10

echo ""
echo "=== Done ==="
echo "Binary: $(pwd)/nmr_extract"
echo ""
echo "Test with:"
echo "  LD_LIBRARY_PATH=$HOME/micromamba/envs/mm/lib:/opt/gromacs/2026.0/lib \\"
echo "    ./nmr_extract --help"
