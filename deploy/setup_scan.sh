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
# Ubuntu 22.04 differences from 24.04:
#   - No libdssp-dev (must build dssp+cifpp from source)
#   - No libfetk-dev (use libmaloc-dev for APBS link)
#
# Usage:
#   ssh scan1
#   bash ~/nmr-shielding/deploy/setup_scan.sh
#
# After running, nmr_extract is at ~/nmr-shielding/build/nmr_extract
#

set -euo pipefail

# Disable broken kitware repo if present
if [ -f /etc/apt/sources.list.d/kitware.list ]; then
    sudo mv /etc/apt/sources.list.d/kitware.list /etc/apt/sources.list.d/kitware.list.disabled
fi

echo "=== Installing apt packages ==="
export DEBIAN_FRONTEND=noninteractive
sudo apt-get update -qq
sudo apt-get install -y -qq \
    libopenbabel-dev \
    libcifpp-dev \
    libapbs-dev libmaloc-dev \
    libboost-dev

# ============================================================================
# Build dssp from source (no libdssp-dev on 22.04)
# ============================================================================
echo "=== Building dssp from source ==="
mkdir -p ~/builds
if [ ! -f /usr/local/lib/libdssp.a ] && [ ! -f /usr/local/lib/libdssp.so ]; then
    cd ~/builds
    if [ ! -d dssp-src ]; then
        git clone https://github.com/PDB-REDO/dssp.git dssp-src
    fi
    cd dssp-src
    mkdir -p build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local 2>&1 | tail -5
    make -j$(nproc) 2>&1 | tail -5
    sudo make install 2>&1 | tail -3
    echo "dssp installed"
else
    echo "dssp already installed"
fi

# ============================================================================
# Build reduce from source
# ============================================================================
echo "=== Building reduce from source ==="
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

# ============================================================================
# Configure and build nmr_extract
# ============================================================================
echo "=== Configuring CMake ==="
cd ~/nmr-shielding
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

echo "=== Building nmr_extract + fes_fleet_smoke_tests ==="
make -j$(nproc) nmr_extract fes_fleet_smoke_tests 2>&1 | tail -10

echo ""
echo "=== Done ==="
echo "Binary: $(pwd)/nmr_extract"
echo ""
echo "Test with:"
echo "  LD_LIBRARY_PATH=$HOME/micromamba/envs/mm/lib:/opt/gromacs/2026.0/lib \\"
echo "    ./fes_fleet_smoke_tests"
