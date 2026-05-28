# Test Health - 2026-05-28

This file records the live test topology and the last verification performed
during the packaging cleanup. Do not preserve stale pass/fail counts here:
rerun the relevant tier and update this file when the catalog or skip behavior
changes.

## C++ CTest Catalog

Current catalog after CMake discovery:

| Label | Tests | Command | Purpose |
| --- | ---: | --- | --- |
| `fast` | 119 | `ctest --test-dir build -L fast --output-on-failure` | Pure logic plus source-tree discipline checks. |
| `conformation` | 303 | `ctest --test-dir build -L conformation --output-on-failure` | Static-structure calculators, topology, charge setup, CLI parsing. |
| `trajectory` | 232 | `ctest --test-dir build -L trajectory --output-on-failure` | AMBER trajectory loading, per-frame trajectory results, H5 round trips. |
| `mopac` | 38 | `ctest --test-dir build -L mopac --output-on-failure` | MOPAC result family and MOPAC-dependent trajectory results. |
| `batch` | 4 | `ctest --test-dir build -L batch --output-on-failure` | Multi-protein consolidated mutant-pair sweeps. |
| `smoke` | 2 | `ctest --test-dir build -L smoke --output-on-failure` | End-to-end producer smoke tests against blessed output policy. |

CTest discovery reports 698 tests total. The labels are primary
suite selectors only. Fixture availability is still decided inside the tests
with `GTEST_SKIP()` and a reason.

## Test Binaries

The top-level CMake build defines these test executables:

| Binary | Label | Notes |
| --- | --- | --- |
| `unit_tests` | `fast` | Logic-only tests plus runtime/config sanity. |
| `string_barrier_tests` | `fast` | Source-tree dependency boundary checks; intentionally does not link `nmr_shielding`. |
| `conformation_scope_tests` | `conformation` | Per-conformation and per-calculator tests. |
| `trajectory_scope_tests` | `trajectory` | TrajectoryResult tests plus AMBER/GROMACS streaming coverage. |
| `mopac_tests` | `mopac` | Slow/optional MOPAC coverage isolated from day-to-day tiers. |
| `batch_tests` | `batch` | Requires `NMR_CONSOLIDATED_DIR` or a local `tests/testpaths.toml` override. |
| `smoke_tests` | `smoke` | Requires current in-tree smoke fixtures and the AIMNet2 model. |

## Fixture Contract

`tests/testpaths.toml` is the default local fixture map. In-tree paths are
repo-root-relative. External/local fixture roots are supplied per machine by
conventional repo-relative local paths, empty opt-in fields, environment
variables, or a local `NMR_TESTPATHS_TOML`.

Important fixture variables:

| Variable | Used for |
| --- | --- |
| `NMR_TESTPATHS_TOML` | Replace the whole fixture map. |
| `NMR_CONSOLIDATED_DIR` | Batch mutant-pair sweeps. |
| `NMR_ORCA_TEST_DIR` | ORCA single-pair fixture directory. |
| `NMR_AIMNET2_MODEL` | Smoke and AIMNet2 calculator tests. |
| `NMR_BASELINE_FEATURES` | Python SDK and write-feature baseline output. |
| `NMR_FLEET_AMBER_DIR` | Large AMBER trajectory fixtures. |
| `NMR_LARSEN_1UBQ_PM6_PDB` | External Larsen PM6-D3H+ 1UBQ geometry. |

`tests/bones/` is ignored archival material and not part of active test
discovery.
Do not move tests there as a soft-delete unless they are genuinely retired,
and do not restore ignored bones code into the active build without an
explicit user decision.

## Golden Output

`tests/golden/smoke/` is generated smoke output. `tests/golden/blessed/` is
ignored for new files, but currently has tracked blessed baselines. Existing
tracked files remain durable; adding a new baseline under that ignored tree
requires an explicit `git add -f`.

## Test Cruft / Watchlist

These are housekeeping observations from source inspection and the commands
listed below:

| Item | Status |
| --- | --- |
| `tests/regression/run_regression.sh` | Manual runner, not wired into CTest. Its active ORCA path now uses `--orca --root`, but the optional baseline comparison is still byte-identical rather than the smoke/bless tolerance policy. Decide whether to keep it manual, wire it into CTest, or fold it into the blessed-output mechanism before treating it as a health signal. |
| `tests/regression/baseline_orca/` | Local ignored `.npy` baseline directory. It is not a durable blessed baseline like tracked files under `tests/golden/blessed/`. |
| Consolidated fixture gating | These tests are runnable, not dead. The portable default leaves `consolidated` empty so CTest skips them; on this machine `/shared/2026Thesis/consolidated` exists and the consolidated-gated reruns passed with `NMR_CONSOLIDATED_DIR` set. Keep this as an explicit fixture dependency rather than baking the local absolute path into tracked config. |
| Larsen PM6 1UBQ fixture | These tests are runnable, not dead. The portable default leaves `larsen_1ubq_pm6_pdb` empty; on this machine the historical archive fixture exists at `/mnt/expansion/larsen_archive/structures/1UBQ_pm6dh3plus.pdb` and all eight gated tests passed with `NMR_LARSEN_1UBQ_PM6_PDB` set. |
| Python SDK fixture skips | `python3 -m pytest python/tests -q -rs` reports 141 passed, 84 skipped. The skips are stale/missing fixture skips: 10 missing `atoms_category_info.npy` for `sdk_geo_only`, 10 missing `/tmp/mutant_smoke`, 60 stale GEO-only extraction rows missing `residues.npy`, 2 stale P84477 baseline rows missing `aimnet2_charges.npy`, and 2 unavailable mutant extraction rows. |

Resolved during this cleanup: `AIMNet2ChargeResponseGradientTest.WriteFeaturesEmitsBothNpys`
was a stale test, not an AIMNet2 availability failure. AIMNet2 loaded,
inference ran, and the backward charge-response gradient completed; the test
still expected the pre-2026-05-20 `aimnet2_polarisability*.npy` names. The
test and nearby active comments were updated to the current
`aimnet2_charge_response_gradient*.npy` contract.

Also resolved during this cleanup: `mopac_tests` built with two format-string
warnings in `tests/test_full_pipeline.cpp` (`%zu` for `int` count accessors).
The format specifiers were corrected and `mopac_tests` rebuilt cleanly.

Also resolved during this cleanup: the two active no-`tleap` negative-path
tests were removed after user agreement. The production failure mode still
exists, and the remaining active Amber preparation tests cover positive `tleap`
use and invalid-input loud failures. The removed tests only exercised the
absence of `tleap` on a host where `tleap` is normally configured, making them
skip-only cruft rather than a useful health signal.

## Last Verification In This Cleanup

Run on 2026-05-27:

```bash
cmake --build build --target unit_tests string_barrier_tests -j$(nproc)
cmake --build build --target conformation_scope_tests -j8
cmake --build build --target trajectory_scope_tests -j8
cmake --build build --target mopac_tests -j8
cmake --build build --target smoke_tests -j8
cmake --build build --target batch_tests -j8
cmake -S . -B build
ctest --test-dir build -N
ctest --test-dir build -N -L fast
ctest --test-dir build -N -L conformation
ctest --test-dir build -N -L trajectory
ctest --test-dir build -N -L mopac
ctest --test-dir build -N -L batch
ctest --test-dir build -N -L smoke
ctest --test-dir build -L fast --output-on-failure
ctest --test-dir build -R 'AIMNet2ChargeResponseGradientTest\.' --output-on-failure
ctest --test-dir build -L conformation --output-on-failure
ctest --test-dir build -L trajectory --output-on-failure
ctest --test-dir build -L mopac --output-on-failure
ctest --test-dir build -L smoke --output-on-failure
ctest --test-dir build -L batch --output-on-failure
/usr/bin/env NMR_CONSOLIDATED_DIR=/shared/2026Thesis/consolidated \
  ctest --test-dir build -R '^(OperationRunnerTest\.(RunMutantComparison|AttachesAllResults|AtomFieldsPopulated|SkipsCoulombWithoutCharges|GeometryAccessible)|BatchHBond\.AllCleanPairs|SampleAtTest\.(BSMatchesAtomValues|BSGridAboveRing|BSButterflyField|AllCalculatorsSample)|WriteFeatures\.P84477Baseline|SmokeTest\.WithDft)$' --output-on-failure
/usr/bin/env NMR_CONSOLIDATED_DIR=/shared/2026Thesis/consolidated \
  ctest --test-dir build -L batch --output-on-failure
/usr/bin/env NMR_LARSEN_1UBQ_PM6_PDB=/mnt/expansion/larsen_archive/structures/1UBQ_pm6dh3plus.pdb \
  ctest --test-dir build -R 'LarsenHBondShieldingTest' --output-on-failure
/usr/bin/env NMR_LARSEN_1UBQ_PM6_PDB=/mnt/expansion/larsen_archive/structures/1UBQ_pm6dh3plus.pdb \
  ctest --test-dir build -R '^(TripeptideBackboneShieldingTest\.RunsOn1UbqPm6|TripeptideNeighborShieldingTest\.RunsOn1UbqPm6)$' --output-on-failure
python3 -m pytest python/tests -q
python3 -m pytest python/tests -q -rs
```

The label counts above are from those commands. The `fast` tier passed
119/119 in 41.23 seconds after the final CMake reconfigure and rebuild.

Additional targeted run on 2026-05-28:

```bash
ctest --test-dir build -N -R AIMNet2
ctest --test-dir build -R AIMNet2 --output-on-failure
```

The AIMNet2 selection listed 12 tests and passed 12/12 in 21.62 seconds
(2 `conformation`, 10 `trajectory`). This exercises model loading/callability,
including response-gradient and trajectory time-series integration tests.

The `conformation` tier initially exposed the stale AIMNet2 output-name test
described above. After the test fix, the final `conformation` run completed
without failures:

Default portable fixture map:

| Label | Pass | Fail | Skip | Wall |
| --- | ---: | ---: | ---: | ---: |
| `fast` | 119 | 0 | 0 | 41.23s |
| `conformation` | 284 | 0 | 19 | 288.57s |
| `trajectory` | 232 | 0 | 0 | 3402.36s |
| `mopac` | 38 | 0 | 0 | 1081.80s |
| `smoke` | 1 | 0 | 1 | 92.61s |
| `batch` | 0 | 0 | 4 | 1.08s |

Local external fixtures enabled:

| Run | Pass | Fail | Skip | Wall |
| --- | ---: | ---: | ---: | ---: |
| Consolidated-gated conformation/smoke subset | 12 | 0 | 0 | 134.18s |
| `batch` with `NMR_CONSOLIDATED_DIR=/shared/2026Thesis/consolidated` | 4 | 0 | 0 | 341.73s |
| Larsen H-bond 1UBQ PM6 group with `NMR_LARSEN_1UBQ_PM6_PDB=/mnt/expansion/larsen_archive/structures/1UBQ_pm6dh3plus.pdb` | 6 | 0 | 0 | 6.62s |
| Tripeptide 1UBQ PM6 pair with `NMR_LARSEN_1UBQ_PM6_PDB=/mnt/expansion/larsen_archive/structures/1UBQ_pm6dh3plus.pdb` | 2 | 0 | 0 | 4.04s |

Skipped conformation entries grouped by cause:

| Cause | Count | Tests |
| --- | ---: | --- |
| `NMR_CONSOLIDATED_DIR` / consolidated P84477 fixture not supplied in the default portable fixture map. These passed when rerun with `/shared/2026Thesis/consolidated`. | 11 | `OperationRunnerTest.RunMutantComparison`; `OperationRunnerTest.AttachesAllResults`; `OperationRunnerTest.AtomFieldsPopulated`; `OperationRunnerTest.SkipsCoulombWithoutCharges`; `OperationRunnerTest.GeometryAccessible`; `BatchHBond.AllCleanPairs`; `SampleAtTest.BSMatchesAtomValues`; `SampleAtTest.BSGridAboveRing`; `SampleAtTest.BSButterflyField`; `SampleAtTest.AllCalculatorsSample`; `WriteFeatures.P84477Baseline` |
| `NMR_LARSEN_1UBQ_PM6_PDB` not supplied in the default portable fixture map. These passed when rerun with `/mnt/expansion/larsen_archive/structures/1UBQ_pm6dh3plus.pdb`. | 8 | `LarsenHBondShieldingTest.*` six-test fixture group; `TripeptideBackboneShieldingTest.RunsOn1UbqPm6`; `TripeptideNeighborShieldingTest.RunsOn1UbqPm6` |

Skipped smoke entry:

| Cause | Count | Tests |
| --- | ---: | --- |
| `NMR_CONSOLIDATED_DIR` / consolidated P84477 fixture not supplied in the default portable fixture map. This passed when rerun with `/shared/2026Thesis/consolidated`. | 1 | `SmokeTest.WithDft` |

Trajectory runtime note: the tier is green but not short. Several AMBER
streaming and frame-PDB-emitter tests each run a full trajectory path and take
about 200-300 seconds apiece.

Skipped batch entries under the default portable fixture map:

| Cause | Count | Tests |
| --- | ---: | --- |
| `NMR_CONSOLIDATED_DIR` / consolidated P84477 fixture not supplied. These passed when rerun with `/shared/2026Thesis/consolidated`. | 4 | `BatchBiotSavartHaighMallion.AllCleanPairs`; `BatchCoulombRingChi.AllCleanPairs`; `BatchMcConnell.AllCleanPairs`; `BatchPiQuadDisp.AllCleanPairs` |

Python pytest result:

| Command | Pass | Fail | Skip | Wall |
| --- | ---: | ---: | ---: | ---: |
| `python3 -m pytest python/tests -q` | 141 | 0 | 84 | 2.09s |
| `python3 -m pytest python/tests -q -rs` | 141 | 0 | 84 | 1.92s |

Python skip groups from `-rs`:

| Cause | Count |
| --- | ---: |
| `tests/data/sdk_geo_only/1Q8K/.../atoms_category_info.npy` missing; set `NMR_CATEGORY_INFO_NPY` or generate the category-info extraction. | 10 |
| `/tmp/mutant_smoke` extraction missing; set `NMR_MUTANT_SMOKE_DIR` or run the mutant smoke extraction. | 10 |
| GEO-only extraction fixture is stale and missing `residues.npy`. | 60 |
| P84477 baseline fixture is stale and missing `aimnet2_charges.npy`. | 2 |
| Mutant extraction fixture is unavailable. | 2 |

## Refresh Procedure

1. Build the affected test binaries, or run `cmake --build build -j$(nproc)` for a full local build.
2. Check the catalog with `ctest --test-dir build -N` and label counts with `ctest --test-dir build -N -L <label>`.
3. Run the tier you are claiming as healthy.
4. Update this file only with counts and skip reasons from commands you actually ran.
