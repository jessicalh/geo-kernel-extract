# Test Inventory — 2026-05-03

Project-wide accounting of what tests we have, what fixtures they hit,
and disposition (live / dead-letter / fixture-free). Captured so the
triage decisions in `review_items_to_assess.md` T2 can land.

`HEAD = 776ca75` at time of inventory.

51 test source files (`tests/test_*.cpp` + `pass0_demo.cpp`) compiled
into **11 test executables** by `CMakeLists.txt:362-557`. Plus
`testpaths.toml` defines fixture path mappings for the
`TestEnvironment` helper.

---

## Fixtures: what each is, force-field family, status

| Fixture | Path | Contents | FF family | Status |
|---|---|---|---|---|
| `1ubq_protonated.pdb` | `tests/data/1ubq_protonated.pdb` | Single protonated PDB, ~100 KB | AMBER ff14SB protonated | **live** — basic structure fixture, not trajectory |
| `external/1UBQ.pdb` | `tests/data/external/1UBQ.pdb` | Unprotonated 1UBQ crystal | n/a (raw PDB) | **live** — used as `ubq_crystal` |
| `external/1OKH_4587_protonated.pdb` | `tests/data/external/1OKH_4587_protonated.pdb` | "GROMACS protonated" PDB | unclear; needs check | **TBD** — `gmx_protonated`; verify FF intent |
| `orca/A0A7C5FAR6_*` | `tests/data/orca/` | WT+ALA pair: PDBs, prmtops, ORCA NMR output, XYZ | AMBER ff14SB + ORCA reference | **live** — ORCA mutant test pair |
| `fleet/` (1A6J_5789, 1AEP_4814) | `tests/data/fleet/` | Older CHARMM/XTC fleet, two proteins. 1A6J trimmed from 10 to 1 pose at master commit `b69d55c` | **CHARMM36m / XTC** | **dead-letter** |
| `fleet_amber/` (1P9J_5801, 1Z9B_6577) | `tests/data/fleet_amber/` | Round-3 Option B 15 ns runs, prep+run from 2026-05-01. Each protein has `sources/`, `decisions.json`, `validator_report.json`, `prep_run_<TS>/...production/` (paths pinned in testpaths.toml). Plus `_backup_round1_*` and `_backup_round2_*` dirs (gitignored) | **AMBER ff14SB / TRR** | **live** — the canonical trajectory fixture |
| `fes_fleet/` (1HD6_4820, 1I8X_4351) | `tests/data/fes_fleet/` | 2 proteins for fes-sampler naming validation | unclear; needs check | **TBD** — fes-sampler is per memory "evolving into libfes for nmr_extract ensemble mode"; live status uncertain |
| `fleet_test_fullsys/` (1ZR7_6721) | `tests/data/fleet_test_fullsys/` | Full-system XTC (protein + water + ions) | **CHARMM36m / XTC** | **dead-letter** — used only by `test_water_field` via `xtc_reader.h` + `GromacsEnsembleLoader` (both quarantined-legacy) |
| `fleet_test_large/` | `tests/data/fleet_test_large/` (gitignored, ~295 MB) | PLUMED MD ensemble: walker_0..4, initial-samples, run_params, harvest_receipt | **CHARMM36m** | **dead-letter, also orphan** — not referenced by any test cpp file |
| `sdk_geo_only/` | `tests/data/sdk_geo_only/` (gitignored, ~137 MB) | Geometry-only SDK extraction (1Q8K) | derived from extraction | **TBD** — Python SDK fixture, not C++ test |
| Consolidated 723 pairs | `/shared/2026Thesis/consolidated/` | 723 WT+ALA mutant pairs (PDB-based) | static structures | **live** — Stage 1 calibration set |
| ff14sb params | `tests/data/../data/ff14sb_params.dat` | Flat ff14SB parameter table | AMBER ff14SB | **live** |
| baseline_features | `baseline_features/P84477/` | NPY baseline for regression | output | **live** (regression check) |
| Smoke run history | `tests/golden/smoke/2026-04-09..2026-04-29*/` | 28 historical smoke runs | mixed; from CHARMM-era + transition | **archive** — not a fixture but worth knowing they're there |
| Smoke blessed | `tests/golden/blessed/withdft/`, `tests/golden/blessed/fleet/` | Blessed baseline NPY for binary diff | per memory `project_smoke_test_bless_deferred_20260424`: 28 NPY binary-diff failures, bless deferred | **deferred** — expected drift vs current code, not yet re-blessed |

---

## Test executables (11) and what they're for

| Executable | Tests | Fixtures touched | Status |
|---|---|---|---|
| `unit_tests` | 7 cpp | none | **live, fixture-free** |
| `structure_tests` | 22 cpp | `1ubq_protonated.pdb`, `orca/`, `external/1OKH_*`, `consolidated/` (one) | **live, static PDB** |
| `trajectory_tests` | 2 cpp | `fleet/` (CHARMM, dead-letter) + `fleet_amber/` (live) | **mixed** — needs split |
| `mopac_tests` | 3 cpp | `1ubq_protonated.pdb`, `orca/` | **live, static PDB** |
| `batch_tests` | 4 cpp | `consolidated/` (the 723 pairs) | **live, static PDB sweep** |
| `smoke_tests` | 1 cpp (`test_smoke.cpp`) | `1ubq_protonated.pdb`, `consolidated/` | **live** — main smoke; produces blessed NPY |
| `fleet_smoke_tests` | 1 cpp (`test_smoke_fleet.cpp`) | `fleet/` CHARMM | **dead-letter** — 1 inert failure (`SmokeFleet.AllPoses` against trimmed fixture) |
| `job_spec_tests` | 1 cpp (`test_job_spec.cpp`) | `1ubq_protonated.pdb`, `external/1UBQ`, `fleet/1A6J_5789`, `orca/` | **mixed** — partly dead-letter via fleet refs |
| `gromacs_streaming_tests` | 1 cpp (`test_amber_streaming.cpp`) | `fleet_amber/` (primary), `fleet_test_fullsys/` (referenced in source) | **live** primary + needs check on fullsys ref |
| `water_field_tests` | 1 cpp (`test_water_field.cpp`) | `fleet_test_fullsys/` CHARMM via `xtc_reader.h` + `GromacsEnsembleLoader` | **dead-letter** — both fixture and loader code are dead-letter |
| `fes_fleet_smoke_tests` | 1 cpp (`test_smoke_fes_fleet.cpp`) | `fes_fleet/` (skipped if not present) | **TBD** — fes-sampler status |

Plus `nmr_pass0_demo` — example driver, not a test.

---

## Per-cpp disposition

### Pure unit tests (fixture-free, all live)

7 in `unit_tests` executable. Pure logic, no protein loading.

| File | Tests |
|---|---|
| `test_spherical_tensor.cpp` | SphericalTensor decomposition / reconstruct |
| `test_amino_acid.cpp` | `AminoAcidType` table, variant indices, residue properties |
| `test_ring_hierarchy.cpp` | Ring type class hierarchy + virtual interface |
| `test_object_model.cpp` | Protein / Residue / ConformationAtom contracts; the new `DisulfideAuthority*` regression tests added 2026-05-03 |
| `test_naming_registry.cpp` | NamingRegistry CHARMM↔Standard translations (still relevant for the ff14SB flat table path) |
| `test_runtime_environment.cpp` | RuntimeEnvironment Load + RequireLoaded |
| `test_calculator_config.cpp` | CalculatorConfig TOML overrides |

Disposition: **all live, keep**.

### Static-PDB tests (`structure_tests` executable)

22 cpp files. All hit `1ubq_protonated.pdb` and/or `orca/A0A7C5FAR6_*`.

**Loader / boundary tests:**
- `test_pdb_loading.cpp` — UbqProtonated. Static PDB load.
- `test_atom_flat.cpp` — UbqProtonated. Atom shape.
- `test_protonation_detection.cpp` — UbqProtonated + GmxProtonated (1OKH).
- `test_protonation_pipeline.cpp` — UbqProtonated + OrcaDir.
- `test_foundation_results.cpp` — UbqProtonated + GmxProtonated + Ff14sbParams.

**Charge / FF tests (recent AMBER slice):**
- `test_amber_charge_resolver.cpp` — UbqProtonated + OrcaDir + Ff14sbParams.
- `test_amber_leap_input.cpp` — UbqProtonated + Ff14sbParams.
- `test_amber_prepared_charge_source.cpp` — UbqProtonated + Ff14sbParams.

**APBS:**
- `test_apbs_field_result.cpp`, `test_apbs_wired.cpp`, `test_apbs_ff14sb.cpp` — all UbqProtonated + Ff14sbParams.

**Geometry / DSSP / Demo:**
- `test_geometry_result.cpp`, `test_dssp_result.cpp`, `test_demo_result.cpp` — UbqProtonated.
- `test_two_conformations.cpp`, `test_traversal_dump.cpp` — UbqProtonated.

**Calculator tests (per-class):**
- `test_biot_savart_result.cpp` — UbqProtonated + OrcaDir.
- `test_haigh_mallion_result.cpp` — UbqProtonated + OrcaDir.
- `test_mcconnell_result.cpp` — UbqProtonated + OrcaDir.
- `test_coulomb_result.cpp` — UbqProtonated + OrcaDir + Ff14sbParams.
- `test_ring_susceptibility_result.cpp` — UbqProtonated + OrcaDir.
- `test_pi_quadrupole_result.cpp` — UbqProtonated + OrcaDir.
- `test_dispersion_result.cpp` — UbqProtonated + OrcaDir.
- `test_hbond_result.cpp` — Consolidated + OrcaDir.

**Pipeline:**
- `test_calculation_runner.cpp` — UbqProtonated + Consolidated + OrcaDir.
- `test_pipeline_and_sample.cpp` — Consolidated.
- `test_write_features.cpp` — Consolidated + BaselineFeatures.

Disposition: **all live, keep**. Static PDB fixtures + ORCA reference pair are unchanged scope. Worth verifying that `gmx_protonated` (1OKH_4587) is still meaningful — its name suggests GROMACS prep, FF unclear.

### MOPAC tests (`mopac_tests` executable)

3 cpp files. ~1 hour runtime.

| File | Fixtures | Disposition |
|---|---|---|
| `test_mopac_result.cpp` | UbqProtonated | **live** |
| `test_full_pipeline.cpp` | UbqProtonated + Ff14sbParams | **live** |
| `test_mutation_delta.cpp` | OrcaDir | **live** |

Disposition: **all live, keep**.

### Batch tests (`batch_tests` executable)

4 cpp files. ~2 hours runtime. All sweep across `Consolidated()` (the 723 pairs).

| File | Disposition |
|---|---|
| `test_batch_mcconnell.cpp` | **live** |
| `test_batch_coulomb_ringchi.cpp` | **live** |
| `test_batch_biot_savart_haigh_mallion.cpp` | **live** |
| `test_batch_piquad_disp.cpp` | **live** |

Disposition: **all live, keep**.

### Trajectory tests (`trajectory_tests` executable)

2 cpp files — but they hit different fixtures and have different status.

| File | Fixtures | Disposition |
|---|---|---|
| `test_fleet_loader.cpp` | `fleet/` CHARMM, via `FleetData()` | **dead-letter** — contains the 4-of-the-5 inert FleetLoader failures (`HasTenFrames`, `PositionsDifferBetweenFrames`, `FullPipelineAllFrames`) per CLAUDE.md 2026-04-29 block |
| `test_amber_trajectory.cpp` | `fleet_amber/` AMBER/TRR | **live** — 6/6 passing per session-handoff-20260502-evening.md |

Disposition: **split required**. The two cpp files don't belong in the same executable: one is dead-letter, one is live. Either:
- Re-derive `test_fleet_loader.cpp` against `fleet_amber/` as `test_amber_fleet_loader.cpp`, then delete the original; or
- Delete `test_fleet_loader.cpp` outright if `test_amber_trajectory.cpp` already covers the same conceptual ground.

### Smoke tests

| Executable | File | Fixtures | Disposition |
|---|---|---|---|
| `smoke_tests` | `test_smoke.cpp` | UbqProtonated + Consolidated + Ff14sbParams | **live** — main smoke. Produces blessed NPY (28 currently in deferred-bless state per memory). |
| `fleet_smoke_tests` | `test_smoke_fleet.cpp` | `fleet/1A6J_5789` CHARMM | **dead-letter** — `SmokeFleet.AllPoses` is the fifth inert failure (1A6J trim). |
| `fes_fleet_smoke_tests` | `test_smoke_fes_fleet.cpp` | `fes_fleet/` | **TBD** — does fes-sampler still drive a live use case? Memory says "evolving into libfes for nmr_extract ensemble mode" but no recent activity logged. |
| `gromacs_streaming_tests` | `test_amber_streaming.cpp` | `fleet_amber/` primary; references `fleet_test_fullsys` | **mostly live** — primary path is AMBER/TRR. The fleet_test_fullsys reference needs eyes-on; might be vestigial from when this file was named `test_gromacs_streaming.cpp` (CHARMM/XTC). |

Disposition:
- `smoke_tests`: **keep**. Re-bless the 28 NPY when ready.
- `fleet_smoke_tests`: **delete or re-derive against fleet_amber**.
- `fes_fleet_smoke_tests`: **decision needed** on fes-sampler scope.
- `gromacs_streaming_tests`: **keep**, audit fullsys ref.

### JobSpec test

| File | Fixtures | Disposition |
|---|---|---|
| `test_job_spec.cpp` | `1ubq_protonated.pdb`, `external/1UBQ`, `fleet/1A6J_5789/params/prod`, `fleet/1A6J_5789/poses/`, `orca/A0A7C5FAR6_*` | **mixed** — `JobSpecE2E.FleetLibraryDirect` is the fifth inert failure (uses fleet/1A6J_5789 trimmed fixture). Other tests in this file likely live. |

Disposition: **partial fix**. The fleet-using tests within this file are dead-letter; rest are live. Either re-derive against fleet_amber, or split the test file.

### Water field test

| File | Fixtures | Disposition |
|---|---|---|
| `test_water_field.cpp` | `fleet_test_fullsys/1ZR7_6721` CHARMM XTC, via `xtc_reader.h` + `GromacsEnsembleLoader` (both quarantined-legacy) | **dead-letter on multiple counts** — fixture is CHARMM-era, loader code is quarantined-legacy, and the test exercises explicit-solvent calculators (`WaterFieldResult`, `HydrationShellResult`, `GromacsEnergyResult`) on that path |

Disposition: **re-derive against existing `fleet_amber/`, no new fixture needed**.

The `fleet_amber/` trajectories are already full-system. `production.trr` (1.8 GB) carries positions + velocities for the whole system; `production.xtc` (546 MB) is configured `compressed-x-grps = System` (`production.mdp:13`) so it also writes the full system. `GromacsFrameHandler.cpp:65` sets `natoms_ = topo.total_atoms` and `:178-179` calls `tp_.SysReader().ExtractFrame(fixed_xyz, protein_positions_, solvent_)` — the protein slice and solvent slice are both populated every frame, and `traj.env_.solvent` carries the per-frame solvent data the explicit-solvent calculators consume.

So the work is: rewrite `test_water_field.cpp` against `fleet_amber/` using the same `Trajectory::Run` + `GromacsFrameHandler` path that `test_amber_streaming.cpp` uses, drop the `xtc_reader.h` and `GromacsEnsembleLoader` includes (quarantined-legacy), exercise `WaterFieldResult` / `HydrationShellResult` / `GromacsEnergyResult` on 1P9J_5801 (or both fixtures). The calculators pick up their `Contract` typedef during the rewrite as part of the T1 cleanup pass.

The "we need to know AMBER works" piece is real: confirming the three solvent calculators produce sane fields on the new path is the actual verification gate. But it is not blocked on fixture preparation — the fixture is fine.

### pass0_demo

| File | Disposition |
|---|---|
| `pass0_demo.cpp` | Example driver, not a gtest. Compiled as `nmr_pass0_demo`. Likely outdated; verify or delete. |

---

## Summary by disposition

**Live, keep as-is:** all of `unit_tests` (7), all of `structure_tests` (22), all of `mopac_tests` (3), all of `batch_tests` (4), `test_smoke.cpp`, `test_amber_trajectory.cpp`, `test_amber_streaming.cpp` (with audit on the fleet_test_fullsys reference). Roughly **38 files**.

**Dead-letter, candidates for delete or re-derive against AMBER fixture:**
- `test_fleet_loader.cpp` — the four of five inert FleetLoader failures.
- `test_smoke_fleet.cpp` — the fifth inert failure.
- `test_water_field.cpp` — fixture + loader both dead-letter; calculator coverage question is real.
- The fleet-using tests within `test_job_spec.cpp` (specifically `JobSpecE2E.FleetLibraryDirect`).

**TBD, decision needed:**
- `test_smoke_fes_fleet.cpp` — fes-sampler scope.
- `external/1OKH_4587_protonated.pdb` (`gmx_protonated`) — verify FF intent.
- `fleet_test_large/` — orphan, no test references; probably delete.
- `pass0_demo.cpp` — verify or delete.

**Deferred (separate decision already made):**
- 28 NPY binary-diff bless deferrals in `tests/golden/blessed/withdft/`.

---

## Adjacent cleanup unlocked

When the dead-letter tests above are deleted or re-derived:

- `tests/data/fleet/` (1A6J_5789 + 1AEP_4814 CHARMM fixtures) become deletable.
- `tests/data/fleet_test_fullsys/1ZR7_6721` becomes deletable.
- `tests/data/fleet_test_large/` becomes deletable (already orphan).
- `src/xtc_reader.h`, `src/GromacsEnsembleLoader.{h,cpp}`, and the `GmxTprChargeSource::LoadCharges` subprocess path lose their last live consumer and become deletable (matches the CHARMM-dead-letter cleanup goal in `review_items_to_assess.md` T1).
- The 28 historical smoke run dirs in `tests/golden/smoke/` from before 2026-04-30 (dated `2026-04-09..2026-04-29*`) can be pruned — they pre-date the current test scope.

---

## What this doc isn't

It isn't the deletion / re-derivation work itself. It's the inventory
that makes that work decidable. Three decisions need conversation
before the cleanup runs:

1. **fes-sampler status** — keep the `test_smoke_fes_fleet.cpp` path
   alive, or retire it?
2. **Solvent calculator coverage** — clarified 2026-05-03: not
   blocked on fixture prep. The current `fleet_amber/` trajectories
   are full-system (TRR by default, XTC via
   `compressed-x-grps = System`), and `GromacsFrameHandler`
   populates the solvent slice every frame. `test_water_field.cpp`
   gets re-derived against `fleet_amber/` on the same
   `Trajectory::Run` path. The real verification — confirming
   `WaterFieldResult` / `HydrationShellResult` / `GromacsEnergyResult`
   produce sane fields on the AMBER path — is the "we need to know
   AMBER works for solvent" gate.
3. **Fleet test re-derivation vs deletion** — do we still want a
   `test_fleet_loader.cpp`-equivalent against the AMBER fleet, or
   does `test_amber_trajectory.cpp` cover the conceptual ground?

The remaining items are mechanical once those three are decided.
