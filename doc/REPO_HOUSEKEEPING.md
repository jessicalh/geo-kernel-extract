# Repository Housekeeping

This map separates build-critical producer material from local state, archives,
and independent subprojects. It is a packaging aid: if a comment here stops
matching the tree, fix the tree or fix the comment.

## Active Producer Surface

These paths are part of building, testing, or understanding `nmr_extract`:

| Path | Role |
| --- | --- |
| `CMakeLists.txt` | Top-level producer build and CTest catalog. |
| `src/` | Producer library and CLI implementation. |
| `tools/` | Source generators and maintenance tools used by the producer. |
| `scripts/` | Runtime/build helper scripts; environment variables are listed in `doc/BUILD_ENVIRONMENT.md`. |
| `data/` | Runtime parameter tables, models, and committed reference data. |
| `tests/` | Active test sources, small fixtures, local fixture config, and ignored large fixture roots. |
| `doc/`, `spec/`, `OBJECT_MODEL.md`, `PATTERNS.md` | Current architecture and design documentation. |
| `python/` | Read-only SDK/loader for producer output. |

## Separate Subprojects

| Path | Status |
| --- | --- |
| `h5-reader/` | Active desktop reader. It has its own CMake/dependency/distribution process, is being handled separately, and is not part of the producer Docker image. |
| `learn/` | Independent calibration and analysis workspace. It is ignored/untracked, kept on disk, and should not be treated as producer build input or staged with producer packaging work. |
| `analysis-speculative/` | Independent scratch/prototype workspace. Keep it on disk, but do not treat it as producer release surface. |
| `references/`, `references-meta/` | Literature corpus and committed metadata. Derived text/page renders are ignored. |

## Release Noise / Git Policy

These paths are currently present in the working tree and may contain tracked
files, but they are outside the producer release boundary. Untracking is an
explicit Git operation; do not move or delete them as part of producer cleanup.

| Path | Current decision |
| --- | --- |
| `site/` | Ignored generated/static output from an older documentation pass. It is not the outward-facing site and should not ship as producer release material. |
| `learn/` | Ignored/untracked independent analysis/calibration workspace, left on disk. |
| `analysis-speculative/` | Ignored/untracked independent scratch/prototype workspace, left on disk. |

## Retired Or Local State

| Path | Status |
| --- | --- |
| `ui/` | Retired legacy viewer. It is not built by the top-level producer CMake. |
| `tests/bones/` | Ignored archival test/code/data material. Not part of active CTest discovery. |
| `bones-make/` | Quarantined stray-root CMake archive/local state. Not an active build tree. |
| `bad-builds/` | Salvaged old build artifacts. Not source. |
| `calibration/` | Large local calibration outputs; committed source/config should live outside ignored run directories. |
| `build*/`, `Testing/`, `_deps/` | Generated build/test state. |

## Large Test Fixtures

| Path | Current size | Contract |
| --- | ---: | --- |
| `tests/data/fleet_amber/` | about 30G | AMBER trajectory fixtures. Gitignored and supplied per machine, usually via `NMR_FLEET_AMBER_DIR`. |
| `tests/golden/smoke/` | about 1.1G | Generated smoke-test output. Gitignored. |
| `tests/golden/blessed/` | about 22M | Existing tracked baselines remain tracked; new files require explicit `git add -f` because the directory is ignored. |
| `tests/bones/` | about 600M | Archive only. Do not depend on it from active code, tests, comments, or docs. |

## Test Tiers

CTest has 698 discovered producer tests in the current build. Primary labels:

| Label | Tests | Typical use |
| --- | ---: | --- |
| `fast` | 119 | Day-to-day logic and source-boundary checks. |
| `conformation` | 303 | Static structure/calculator coverage. |
| `trajectory` | 232 | Trajectory and H5 coverage. |
| `mopac` | 38 | Slow/optional MOPAC coverage. |
| `batch` | 4 | Consolidated mutant-pair sweeps; external fixture tree. |
| `smoke` | 2 | End-to-end producer output comparison. |

See `tests/TEST_HEALTH.md` for commands and refresh rules.

## Packaging Rules

Do not include local build trees, generated smoke output, `tests/bones/`,
`bones-make/`, `bad-builds/`, `site/`, `learn/`, `analysis-speculative/`, or
large local fixture trees in a Docker build context unless a specific image
needs them.

For Docker and other new machines, start from `doc/BUILD_ENVIRONMENT.md`,
`CMakeUserPresets.example.json`, `tests/testpaths.toml`, and the runtime
`NMR_*` variables rather than copying local absolute paths.
