# External Dependencies

Complete inventory of every external library, tool, and binary.
For setting up on a new machine.

Last updated: 2026-04-03.

---

## Linked C++ Libraries

### Eigen3 (linear algebra)
- **Version:** 3.4.0
- **Install:** `apt install libeigen3-dev`
- **CMake:** `find_package(Eigen3 REQUIRED)`
- **Used for:** Vec3, Mat3, SVD throughout all calculators

### DSSP + cifpp (secondary structure + PDB parsing)
- **Version:** dssp 4.2.2, cifpp 5.0.7
- **Install:** `apt install libdssp-dev libcifpp-dev`
- **CMake:** `find_package(dssp REQUIRED)` (brings cifpp)
- **Used for:** DsspResult (SS, phi/psi, SASA, H-bonds), PDB/mmCIF parsing

### OpenBabel3 (bond perception)
- **Version:** 3.1.1
- **Install:** `conda install -c conda-forge openbabel`
- **Paths:** lib `miniforge3/lib/libopenbabel.so`, headers `miniforge3/include/openbabel3/`
- **CMake:** Manual `OPENBABEL_INCLUDE`, `OPENBABEL_LIB`
- **Used for:** Bond detection and hybridisation in CovalentTopology

### GROMACS 2026.0 (TPR reading)
- **Install:** Build from source (cmake, make)
- **Paths:**
  - Library: `~/gromacs/lib/libgromacs.so`
  - Source headers: `builds/gromacs-2026.0/src/` (internal API)
  - Build headers: `builds/gromacs-2026.0/build/api/legacy/include/`
- **CMake:** `GROMACS_LIB`, `GROMACS_SRC`, `GROMACS_BUILD`
- **Used for:** Reading TPR binary topology in GromacsEnsembleLoader

### reduce (hydrogen placement)
- **Version:** 4.10 (Richardson lab, Duke University)
- **Source:** https://github.com/rlabduke/reduce
- **Build:** `git clone`, then `mkdir build && cd build && cmake .. && make reducelib`
- **Produces:** three static libraries:
  - `build/reduce_src/libreducelib.a`
  - `build/libpdb/libpdb++.a`
  - `build/toolclasses/libtoolclasses.a`
- **Data:** `reduce_wwPDB_het_dict.txt` (~60MB, ships with source)
- **CMake:** `REDUCE_SRC`, `REDUCE_LIB`, `REDUCE_LIBPDB`, `REDUCE_TOOLCLASSES`, `REDUCE_HET_DICT`
- **Compile defs needed:** `BOOLPREDEFINED CHARFUNCMACROS BRACKETOPERPARMS
  LISTFRIENDFIX INCTEMPLATEDEFNS LEFT_JUSTIFY_NUC_RES_OK
  AROMATICS_ACCEPT_HBONDS HET_DICTIONARY="<path>"`
- **Used for:** Adding H atoms to bare PDBs with PDB-standard naming
- **Note:** C++ globals for configuration. SetReduceBuildMode() resets before each call.

### APBS (Poisson-Boltzmann solvation)
- **Version:** 3.4.1
- **Install:** `apt install libapbs-dev libfetk-dev`
- **Bridge:** Custom C file `biot-savart/src/External/apbs_bridge.c`
- **Link:** `apbs_routines apbs_mg apbs_generic apbs_pmgc mc punc maloc`
- **Used for:** Solvated E-field and EFG (ApbsFieldResult)
- **Known issue:** C global state does not reinitialise. Needs fix for viewer.

---

## Header-Only Libraries (in extern/)

| Library | File | Used for |
|---------|------|----------|
| sphericart | sphericart.h, sphericart.hpp, templates*.hpp | Spherical tensor decomposition (T0, T1, T2) |
| nanoflann | nanoflann.hpp | KD-tree spatial index (15A neighbour lists) |

---

## External Binaries (called via std::system)

### xTB (semiempirical quantum chemistry)
- **Install:** `conda create -n xtb-env && conda install -c conda-forge xtb`
- **Path:** `miniforge3/envs/xtb-env/bin/xtb`
- **Invoked as:** `xtb input.xyz --chrg N --uhf 0 --gfn 2 --json`
- **Used for:** Mulliken charges, polarisability, bond orders
- **Size limit:** Core dumps >~1000 atoms (GFN2 memory scaling)

### PROPKA (pKa prediction)
- **Install:** `pip install propka`
- **Path:** `miniforge3/bin/propka3`
- **Invoked as:** `propka3 --quiet input.pdb`
- **Used for:** Predicting pKa of titratable residues

### gmx (GROMACS command-line)
- **Path:** `~/gromacs/bin/gmx`
- **Invoked as:** `gmx dump -s file.tpr`
- **Used for:** Legacy GmxTprChargeSource (fleet loader uses libgromacs directly)

### KaML (ML pKa) — optional
- **Path:** `~/KaML/KaML-CBTrees/KaML-CBtree.py`
- **~80% success rate on ARM64**

### tleap (AMBER topology) — offline only
- **Path:** `~/amber24/bin/tleap`
- **Not called at runtime** — used to prepare training data prmtop files

---

## Data Files

| File | Path | Used by |
|------|------|---------|
| ff14SB params | `biot-savart/data/ff14sb_params.dat` | ParamFileChargeSource (422 entries) |
| reduce het dict | `builds/reduce-src/reduce_wwPDB_het_dict.txt` | reduce protonation |
| 1UBQ protonated | `tests/data/1ubq_protonated.pdb` | All PDB-loading tests |

---

## Build Summary

**System packages (apt):**
```bash
apt install libeigen3-dev libdssp-dev libcifpp-dev libapbs-dev libfetk-dev
```

**Conda packages (miniforge3):**
```bash
conda install -c conda-forge openbabel
pip install propka
conda create -n xtb-env && conda install -n xtb-env -c conda-forge xtb
```

**Build from source:**
```bash
# GROMACS
cd builds && tar xf gromacs-2026.0.tar.gz && cd gromacs-2026.0
mkdir build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX=~/gromacs && make -j8 && make install

# reduce
cd builds && git clone https://github.com/rlabduke/reduce.git reduce-src
cd reduce-src && mkdir build && cd build && cmake .. && make -j8 reducelib
```

**Existing installations (if available):**
```
AMBER24  →  ~/amber24/   (for tleap, reduce binary fallback)
KaML     →  ~/KaML/      (optional ML pKa predictor)
```

**Test framework:** GTest 1.14.0 (fetched automatically by CMake FetchContent)

---

## Not Yet Used (future)

- **e3nn** (Python): equivariant neural network for ParameterCorrectionResult
- **LibTorch** (C++): TorchScript inference for C++ model integration
