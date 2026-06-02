# Build Environment

The committed build files do not encode machine-local user paths. Configure
dependency paths with CMake cache variables, usually through the gitignored
`CMakeUserPresets.json`, or seed those cache variables with `NMR_*`
environment variables on the first configure.

Resolution order is:

1. Explicit CMake cache variable (`-D...` or `CMakeUserPresets.json`)
2. Matching `NMR_*` environment variable
3. Standard system discovery and conservative compatibility probes

## Core configure

`NMR_PROFILE=local` is the default and allows machine-local dependency paths
for development. `NMR_PROFILE=producer-full` is the release/container gate: it
enables `NMR_REQUIRE_PORTABLE_DEPS`, uses the installed data directory, disables
`-march=native`, and fails configure if producer dependencies resolve through
personal or host-local paths such as `/home`, `/shared`, `/mnt`, `/opt/orca`,
or `/usr/local/cuda`. The producer Dockerfile may explicitly allow
`/usr/local/cuda` when it comes from the selected NVIDIA CUDA base image; all
project-specific vendored dependencies still belong under the producer prefix.

The intended release configure shape is:

```bash
cmake -S . -B build-producer \
  -DNMR_PROFILE=producer-full \
  -DCMAKE_INSTALL_PREFIX=/opt/nmr-shielding \
  -DNMR_GROMACS_ROOT=/opt/nmr-shielding/deps/gromacs \
  -DNMR_MOPAC_ROOT=/opt/nmr-shielding/deps/chem-env \
  -DNMR_TORCH_CMAKE_PREFIX_PATH=/opt/nmr-shielding/deps/torch/share/cmake
```

| Dependency | Cache variable | Environment variable | Notes |
| --- | --- | --- | --- |
| profile | `NMR_PROFILE` | none | `local` or `producer-full`. |
| portable dependency gate | `NMR_REQUIRE_PORTABLE_DEPS`, `NMR_PORTABLE_DEP_ALLOW_REGEXES` | none | Fails configure on machine-local producer dependency paths when enabled. |
| system library dir | `NMR_SYSTEM_LIBRARY_DIR` | `NMR_SYSTEM_LIBRARY_DIR` | Prepended to RUNPATH/`LD_LIBRARY_PATH` so `libpq`, `libssl`, and `libcrypto` resolve from one install. |
| libpq | `PostgreSQL_LIBRARY`, `PostgreSQL_INCLUDE_DIR` | `NMR_POSTGRESQL_LIBRARY`, `NMR_POSTGRESQL_INCLUDE_DIR` | Linked unconditionally for the local tensorcs15 table reader. |
| DSSP | `NMR_DSSP_ROOT`, `DSSP_LIB`, `DSSP_INCLUDE` | `NMR_DSSP_ROOT`, `NMR_DSSP_LIB`, `NMR_DSSP_INCLUDE` | `find_package(dssp)` is tried first; these variables drive the manual fallback. |
| OpenBabel | `NMR_OPENBABEL_ROOT`, `OPENBABEL_INCLUDE`, `OPENBABEL_LIB` | `NMR_OPENBABEL_ROOT`, `NMR_OPENBABEL_INCLUDE`, `NMR_OPENBABEL_LIB` | Used for bond perception and hybridisation. |
| GROMACS | `NMR_GROMACS_ROOT`, `GROMACS_LIB`, `GROMACS_SRC`, `GROMACS_BUILD` | `NMR_GROMACS_ROOT`, `NMR_GROMACS_LIB`, `NMR_GROMACS_SRC`, `NMR_GROMACS_BUILD` | Needs both the library and internal source/build headers. |
| reduce | `REDUCE_SRC`, `REDUCE_LIB`, `REDUCE_LIBPDB`, `REDUCE_TOOLCLASSES`, `REDUCE_HET_DICT` | `NMR_REDUCE_SRC`, `NMR_REDUCE_LIB`, `NMR_REDUCE_LIBPDB`, `NMR_REDUCE_TOOLCLASSES`, `NMR_REDUCE_HET_DICT` | Linked archives used by protonation. |
| MOPAC | `NMR_MOPAC_ROOT`, `MOPAC_INCLUDE`, `MOPAC_LIB` | `NMR_MOPAC_ROOT`, `NMR_MOPAC_INCLUDE`, `NMR_MOPAC_LIB` | Linked unconditionally; runtime `--no-mopac` only skips execution. |
| HDF5 | `NMR_HDF5_INCLUDE_DIR`, `NMR_HDF5_C_LIBRARY` | same names | Kept explicit to avoid header/library ABI mismatches. |
| HDF5 transitives | `NMR_CURL_LIBRARY`, `NMR_Z_LIBRARY` | same names | Direct link dependencies for the HDF5 C library. |
| APBS/MALOC | `NMR_APBS_INCLUDE_DIR`, `NMR_MALOC_INCLUDE_DIR`, `NMR_APBS_*_LIBRARY`, `NMR_MALOC_LIBRARY` | same names | APBS bridge libraries. FETK `mc`/`punc` remain optional when bundled into APBS. |
| Torch | `NMR_TORCH_CMAKE_PREFIX_PATH` or `Torch_DIR` | `NMR_TORCH_CMAKE_PREFIX_PATH` | If unset, CMake asks the configured Python interpreter where pip-installed torch lives. |
| cu13 runtime | `NMR_NVRTC_LIB_DIR` or `NMR_NVIDIA_CU13_DIR` | same names | Used by tests and `scripts/run_with_cuda_env.sh` for `libnvrtc-builtins.so.13.0`; this must match the producer Torch CUDA ABI (`2.11.0+cu130`, CUDA 13.0). |
| GTest fetch | `NMR_FETCH_GTEST` | none | `ON` preserves current behavior; `OFF` requires system GTest. |

Compiler behavior that can vary by target machine is also explicit:

| Setting | Default | Purpose |
| --- | --- | --- |
| `NMR_ENABLE_NATIVE_ARCH` | `ON` | Adds `-march=native`. Disable for portable binaries or Docker images intended to run on older CPUs. |
| `NMR_ENABLE_LTO` | `ON` | Adds `-flto`. Disable for toolchains or containers where LTO is unavailable. |
| `NMR_EXTRA_CXX_FLAGS` | empty | Semicolon-separated additional C++ flags. |

## Scripts

Runtime/build scripts use the same naming style:

| Script | Variables |
| --- | --- |
| `scripts/run_with_cuda_env.sh` | `NMR_NVRTC_LIB_DIR`, `NMR_NVIDIA_CU13_DIR`, `NMR_PYTHON` |
| `artifacts/test_runner.sh` | `NMR_AIMNET2_MODEL`, `NMR_CALCULATOR_CONFIG`, `NMR_CUDA_WRAPPER`, `NMR_BUILD_DIR`, `NMR_CMAKE_PRESET`, `NMR_PM6DH3PLUS_PDB`, `NMR_MODE5_TRAJ`, `NMR_REQUIRE_OPTIONAL_FIXTURES` |
| `deploy/setup_scan.sh` | `NMR_PROJECT_ROOT`, `NMR_BUILD_PARENT`, `NMR_BUILD_DIR`, `NMR_DSSP_SRC`, `NMR_REDUCE_SRC`, `NMR_DSSP_REPO`, `NMR_REDUCE_REPO`, `NMR_SCAN_INSTALL_APT`, `NMR_SCAN_BUILD_DSSP`, `NMR_SCAN_BUILD_REDUCE`, `NMR_SCAN_INSTALL_PREFIX`, `NMR_SCAN_APT_PACKAGES`, `NMR_GROMACS_ROOT`, `NMR_GROMACS_LIB`, `NMR_GROMACS_SRC`, `NMR_GROMACS_BUILD`, `NMR_MOPAC_ROOT`, `NMR_MOPAC_INCLUDE`, `NMR_MOPAC_LIB`, `NMR_OPENBABEL_LIB`, `NMR_SCAN_EXTRA_TARGETS` |
| `tests/regression/run_regression.sh` | `NMR_EXTRACT`, `NMR_BUILD_DIR`, `NMR_TESTPATHS_TOML`, `NMR_REGRESSION_DIR`, `NMR_REGRESSION_ORCA_ROOT`, `NMR_CALCULATOR_CONFIG`, `NMR_AIMNET2_MODEL` |
| `calibration/populate.py` | `NMR_CONSOLIDATED_DIR` |
| `scripts/batch_extract.sh` | `NMR_CONSOLIDATED_DIR`, `NMR_EXTRACT` |
| `scripts/larsen_hbond_grid_parse/parse_larsen_hbond_grids.py` | `NMR_LARSEN_ARCHIVE_TAR` |
| `tools/amber/generate_ff14sb_pb_table.py` | `NMR_TLEAP`, `AMBERHOME` |

`artifacts/test_runner.sh` skips optional external-heavy fixtures when they
are absent: the PM6-DH3+ protonated PDB (`NMR_PM6DH3PLUS_PDB`) and the
large Mode 5 trajectory (`NMR_MODE5_TRAJ`). Set
`NMR_REQUIRE_OPTIONAL_FIXTURES=1` in CI/full-validation environments to turn
those skips into preflight failures.

`deploy/setup_scan.sh` still source-builds DSSP and reduce by default. Set
`NMR_SCAN_BUILD_DSSP=0` or `NMR_SCAN_BUILD_REDUCE=0` only when those
dependencies are already provided by the image or machine.

The retired `ui/` viewer is ignored local archive, not part of the release repo
or active producer build. Desktop UI work lives in `h5-reader/`, which has its
own CMake presets and distribution process. `artifacts/test_runner.sh` does not
build h5-reader; this keeps
producer/Docker validation separate from desktop Qt/VTK dependency choices.

## Runtime Config

The default runtime TOML remains `~/.nmr_tools.toml`. For containers and
non-home deployments, set `NMR_TOOLS_TOML=/path/to/nmr_tools.toml`.

| Runtime item | TOML key | Environment variable |
| --- | --- | --- |
| MOPAC binary | `mopac` | `NMR_MOPAC` |
| AmberTools tleap | `tleap` | `NMR_TLEAP`, then `AMBERHOME/bin/tleap` |
| ff14SB parameter table | `ff14sb_params` | `NMR_FF14SB_PARAMS` |
| temp directory | `tmpdir` | `NMR_TMPDIR` |
| BMRB atom-name table (`data/bmrb_atom_nom.tbl`) | `bmrb_atom_nom` | `NMR_BMRB_ATOM_NOM` |
| tensorcs15 DSN | `[databases].tensorcs15` | `NMR_TENSORCS15_DSN` |
| Larsen H-bond grids | `larsen_hbond_grids` | `NMR_LARSEN_HBOND_GRIDS` |
| UDP log host/port | `[logging].udp_host`, `[logging].udp_port` | `NMR_LOG_UDP_HOST`, `NMR_LOG_UDP_PORT` |
| file log path | `[logging].file` | `NMR_LOG_FILE` |

`data/calculator_params.toml` is committed with a repo-relative AIMNet2 model
path. `nmr_extract` resolves relative `aimnet2_model_path` values from the
TOML file directory. Override with `--aimnet2` or `NMR_AIMNET2_MODEL` for
alternate model assets.

## Test Fixtures

`tests/testpaths.toml` keeps in-tree fixtures repo-relative. External/local
fixtures are supplied per machine either by conventional repo-relative local
paths, empty opt-in fields, or environment variables:

| Fixture | Environment variable |
| --- | --- |
| consolidated mutant-pair tree | `NMR_CONSOLIDATED_DIR` |
| ORCA test fixture directory | `NMR_ORCA_TEST_DIR` |
| AIMNet2 model | `NMR_AIMNET2_MODEL` |
| baseline feature output | `NMR_BASELINE_FEATURES` |
| AMBER trajectory fixture root | `NMR_FLEET_AMBER_DIR` |
| Larsen PM6-DH3+ 1UBQ PDB | `NMR_LARSEN_1UBQ_PM6_PDB` |

The consolidated mutant-pair tree and Larsen PM6-DH3+ PDB default to empty
because there is no portable in-repo location. The AMBER trajectory root,
AIMNet2 model, and baseline feature output point at conventional local paths;
tests still skip with explicit reasons if those paths are absent or stale.

CTest tier labels and large fixture/archive directory policy are documented in
`tests/TEST_HEALTH.md` and `doc/REPO_HOUSEKEEPING.md`.

For a local machine, copy `CMakeUserPresets.example.json` to
`CMakeUserPresets.json`, edit the paths, and run:

```bash
cmake --preset local
cmake --build --preset local
```
