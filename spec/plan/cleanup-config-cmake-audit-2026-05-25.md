# Cleanup / Config / CMake Audit - Working Inventory (2026-05-25)

Purpose: internal implementation inventory for the "proper program" cleanup
phase. Read after `OBJECT_MODEL.md` and `PATTERNS.md`. This is not an
architecture proposal; it is a path/line map of cleanup candidates already
found, split into three tracks:

1. General cleanup: local duplication, schema drift, helper gaps.
2. Runtime/config cleanup: literals, paths, DB/config leakage, model params.
3. CMake/configure cleanup: build portability and install cadence.

Target cadence for the CMake track:

```bash
cmake --preset default
cmake --build --preset default
ctest --preset default
cmake --install build
```

`cmake --install` does not exist as a real project surface yet; that is the
main missing piece.

## Progress (updated 2026-05-25)

Tracks 1–2 are essentially done; Track 3 (CMake) is untouched and is its own
careful phase (do it behind a backup — see the safety gate below). Each landed
item was chunked + dual/codex-reviewed + (where numeric) smoke-blessing-verified.

**Done + pushed:**
- G8 / C1 / C2 — unit literals → named PhysicalConstants (`73ed79d`).
- C3 / C5 — H5 metadata reads CalculatorConfig; conda fallbacks dropped,
  mopac/tleap now TOML→AMBERHOME→PATH (`d843d15`).
- G7 — h5-reader stale byte-count comment 84→83 (`6eb41aa`).
- G1 — TopologySidecar manifest → nlohmann/ordered_json (`be57b4d`).
- G2 — SphericalTensor::PackFull9/PackT2 own the [T0,T1,T2] order; all 40
  former hand-rolled packers migrated (`2c60636` + `8ee12f4`, total).
- C4 — APBS grid/dielectric/temperature/ionic params → 8 TOML keys (`f212747`).
- C6 — already resolved by cleanup chunk B (UniqueTempPath removed, `1979c52`).
- C7 / C8 / C9 — script hygiene: DSN redaction, machine-path de-hardcoding
  (`5f01828`). **C8 caveat:** batch_extract.sh paths fixed, but its retired
  per-file --mutant/--orca flags were NOT migrated — the canonical --root
  expansion can't map the timestamp-stemmed consolidated data without the
  forbidden discovery glob. Left with a dated in-file note: retire the script
  or re-stem the data (Jessica's decision).

**Done + pushed (2026-05-26 — closes Tracks 1–2 except G9):**
- C10 — testpaths relative resolution + G4 Larsen fixture + G5 temp helper,
  one TestEnvironment portability chunk (`9d288de`). Also fixed
  run_regression.sh's read_toml to mirror the same resolution; its stale
  `--orca DIR` CLI shape (pre-existing, separate) noted in-file for a re-bless.
- G3 — 54 test-local `LoadCalculatorConfig()` helpers → one
  `TestEnvironment::LoadCalculatorConfig()`, 234 calls requalified (`156c55e`).
  The old local path was non-existent → silent default; now loads the real
  TOML (== defaults, so inert).
- G6 — `TrajectoryInputFiles::FromProductionDir` owns production.{tpr,trr,edr}
  in Cli/ModeSpec.h; Validate + nmr_extract use it (`e58b1bb`).

**Remaining:**
- G9 — shared H5 frame-dataset name constants (codex: not first pass; deferred).
- Track 3 (CMake/install B1–B8) — **not started.** The whole track. Do behind
  a backup; verify the configure→build→test→install cadence before folding in.

**Audit-validation pass (2026-05-26):** an agent re-validated every item
against the tree (all accurate; no inaccurate claims) and swept for new
instances in the same categories. Landed (`d09ab02`): a `/mnt/expansion`
argparse default in parse_larsen_hbond_grids.py the C9 sweep missed (now
env-overridable), and a fail-loud log on PdbFileReader's non-numeric-seq_id
fallback. Investigated and REJECTED: a proposed "APBS reads
efield_magnitude_sanity_clamp like its siblings" change — codex caught that
the APBS clamp runs in pre-conversion kT/(e·A) units while the shared key is
V/A, so sharing it would mis-scale APBS under any non-default override. The
APBS clamp keeps APBS_SANITY_LIMIT; the distinction is now documented in-code
(ApbsFieldResult.cpp) so it isn't re-flagged. Whether APBS *should* clamp
post-conversion in V/A to truly match siblings is a physics/blessing question,
not config plumbing.

## Existing Homes To Use

Do not invent new global utility buckets unless one of these homes is wrong.

| Concern | Existing home | Anchors |
|---|---|---|
| Citable fixed constants and unit conversions | `src/PhysicalConstants.h` | Header rule at `src/PhysicalConstants.h:1-12`; unit section starts `src/PhysicalConstants.h:29` |
| Tunable calculator parameters | `src/CalculatorConfig.{h,cpp}` + `data/calculator_params.toml` | API in `src/CalculatorConfig.h:1-38`; defaults in `src/CalculatorConfig.cpp:51-91`; TOML keys in `data/calculator_params.toml:28-74` |
| Runtime tools, env vars, data dirs, DSNs | `src/RuntimeEnvironment.{h,cpp}` + `src/Session.cpp` | TOML key parsing at `src/RuntimeEnvironment.cpp:139-146`; path resolution at `src/RuntimeEnvironment.cpp:155-228`; DSN consumer at `src/Session.cpp:49-68` |
| Test fixture paths | `tests/TestEnvironment.{h,cpp}` + `tests/testpaths.toml` | Interface `tests/TestEnvironment.h:34-59`; loader `tests/TestEnvironment.cpp:34-91`; current path file `tests/testpaths.toml:1-74` |
| Naming/schema surfaces | `src/NamingRegistry.h`, `src/SemanticEnums.h`, `src/CategoryInfoProjection.cpp`, `src/TopologySidecar.cpp` | NPY atom category layout starts `src/CategoryInfoProjection.cpp:313`; topology manifest writer starts `src/TopologySidecar.cpp:468` |
| Root CMake presets | `CMakePresets.json` | configure presets `CMakePresets.json:4-23`; build presets `CMakePresets.json:30-32`; test presets `CMakePresets.json:34-47` |
| h5-reader build pattern | `h5-reader/CMakePresets.json` | preset binary dir `h5-reader/CMakePresets.json:10` |

## Track 1 - General Cleanup

### G1. TopologySidecar manifest should use structured JSON

Issue: `TopologySidecar::WriteManifest` manually streams JSON. This repeats
the risk just fixed elsewhere: escaping and key-order stability should come
from `nlohmann::ordered_json`, not ad hoc string assembly.

Anchors:

- Manual manifest writer: `src/TopologySidecar.cpp:468`, stream starts at
  `src/TopologySidecar.cpp:478`.
- Existing good pattern: `src/GromacsToAmberReadbackBlock.cpp:103-134`.
- Existing good pattern: `src/OperationLog.cpp:98-112`.
- Existing good pattern: `src/Trajectory.cpp:403-410`.

Likely fix:

- Add `#include <nlohmann/json.hpp>` to `TopologySidecar.cpp`.
- Build a `nlohmann::ordered_json` tree and dump with
  `ordered_json::error_handler_t::replace`, matching the existing emitters.
- Preserve exact keys, schema version, and insertion order.

Blast radius: light. Output file only: `extraction_manifest.json`. Tests should
compare presence/key values, not whitespace.

### G2. SphericalTensor should own its packing order

Issue: the same `[T0, T1[3], T2[5]]` packing is hand-coded in several places.
The object model already says `SphericalTensor` owns the decomposition; the
pack order is part of that representation and should be one named API.

Anchors:

- Type definition: `src/Types.h:268-299`.
- Decompose/reconstruct implementation: `src/Types.cpp:25-73`.
- Local packer: `src/ConformationResult.cpp:16-20`.
- Local packer: `src/BiotSavartResult.cpp:453-459`.
- Local packer: `src/CoulombResult.cpp:365-369`.
- Local packer: `src/LarsenHBondShieldingResult.cpp:924-934`.
- Direct T2 loops to consider after full-pack API exists:
  `src/CoulombResult.cpp:393-394`, `src/BiotSavartResult.cpp:486-493`.

Likely fix:

- Add methods such as `void PackFull9(double* out) const` and optionally
  `void PackT2(double* out) const` on `SphericalTensor`.
- Replace local packers only where they are literal full-9 or T2-5 packing.
- Do not centralize NPY file names or calculator-specific column layouts in
  `SphericalTensor`.

Blast radius: light to medium. Many touched files, but behavior should be
byte-identical if the order is preserved. Good verification target:
feature-writing tests and blessed NPY comparisons.

### G3. Repeated test `LoadCalculatorConfig` helpers

Issue: tests repeatedly define a local `LoadCalculatorConfig()` that rebuilds
`std::string(NMR_TEST_DATA_DIR) + "/../data/calculator_params.toml"`. This is
pure test glue and belongs in `TestEnvironment`.

Existing helper home:

- `tests/TestEnvironment.h:34-59`
- `tests/TestEnvironment.cpp:34-91`
- `tests/test_main.cpp:22` already loads calculator params globally with a
  different relative path: `NMR_TEST_DATA_DIR + "/../../data/calculator_params.toml"`.

Repeated local helper anchors:

```text
tests/test_larsen_hbond_water_term_time_series.cpp:77,80
tests/test_hydration_shell_time_series.cpp:56,58
tests/test_mopac_bond_order_welford.cpp:60,63
tests/test_eeq_welford.cpp:54,56
tests/test_frame_pdb_emitter.cpp:78,81
tests/test_apbs_efg_time_series.cpp:61,64
tests/test_larsen_hbond_count_time_series.cpp:78,81
tests/test_mopac_mc_shielding_time_series.cpp:62,65
tests/test_piquad_shielding_time_series.cpp:61,63
tests/test_mc_welford.cpp:64,67
tests/test_dihedral_bin_transition.cpp:57,60
tests/test_rmsd_tracking.cpp:51,54
tests/test_hydration_geometry_time_series.cpp:57,59
tests/test_mopac_vs_ff14sb_reconciliation.cpp:60,63
tests/test_larsen_hbond_1pHaB_shielding_time_series.cpp:66,69
tests/test_mopac_coulomb_shielding_time_series.cpp:60,63
tests/test_sasa_time_series.cpp:58,60
tests/test_gromacs_energy_time_series.cpp:55,57
tests/test_disp_shielding_time_series.cpp:61,63
tests/test_tripeptide_bb_method_tag_time_series.cpp:78,81
tests/test_tripeptide_bb_shielding_time_series.cpp:110,113
tests/test_mc_shielding_time_series.cpp:75,78
tests/test_ringchi_shielding_time_series.cpp:61,63
tests/test_tripeptide_neighbor_residual_vec_next_time_series.cpp:72,75
tests/test_dssp8_time_series.cpp:55,58
tests/test_larsen_hbond_1pHB_shielding_time_series.cpp:77,80
tests/test_j_coupling_time_series.cpp:63,66
tests/test_bs_anomalous_atom_marker.cpp:80,83
tests/test_ring_pucker_time_series.cpp:56,59
tests/test_larsen_hbond_2pHB_shielding_time_series.cpp:66,69
tests/test_water_field_time_series.cpp:57,59
tests/test_hydration_geometry_welford.cpp:70,72
tests/test_mopac_charge_welford.cpp:75,78
tests/test_apbs_efield_time_series.cpp:59,61
tests/test_hm_shielding_time_series.cpp:82,85
tests/test_aimnet2_embedding_charge_response_gradient_time_series.cpp:65,68
tests/test_tripeptide_neighbor_residual_vec_prev_time_series.cpp:86,89
tests/test_bonded_energy_time_series.cpp:56,58
tests/test_aimnet2_charge_time_series.cpp:61,63
tests/test_tripeptide_neighbor_shielding_time_series.cpp:103,106
tests/test_hydration_shell_welford.cpp:65,67
tests/test_tripeptide_bb_residual_vec_time_series.cpp:105,108
tests/test_bs_welford.cpp:72,75
tests/test_sasa_welford.cpp:52,54
tests/test_water_field_welford.cpp:59,61
tests/test_ring_neighbourhood_trajectory_stats.cpp:61,64
tests/test_amber_streaming.cpp:96,99
tests/test_hbond_shielding_time_series.cpp:63,65
tests/test_rmsd_spike_dft_coord.cpp:58,61
tests/test_dihedral_time_series.cpp:62,65
tests/test_hm_welford.cpp:65,68
tests/test_larsen_hbond_2pHaB_shielding_time_series.cpp:66,69
tests/test_dssp8_transition.cpp:53,56
tests/test_hbond_count_welford.cpp:53,55
tests/test_calculator_config.cpp:78
```

Likely fix:

- Add `TestEnvironment::CalculatorParams()` or
  `TestEnvironment::LoadCalculatorConfig()`.
- Remove local helper definitions from tests.
- Decide whether `tests/test_main.cpp:22` should be the only loader. If a
  test deliberately reloads config, keep an explicit call but use the helper.

Blast radius: light but broad. Test-only.

### G4. External Larsen fixture path belongs in TestEnvironment

Issue: several tests hardcode the same external Larsen archive PDB path. The
tests already know how to skip when a fixture is missing; the path should be
configured through `testpaths.toml` like the other external fixtures.

Anchors:

- `tests/test_tripeptide_backbone_shielding.cpp:52`
- `tests/test_tripeptide_neighbor_shielding.cpp:45`
- `tests/test_larsen_hbond_shielding.cpp:55`
- Related comment-only archive dependency:
  `tests/test_larsen_residue_against_source_log.cpp:12`

Likely fix:

- Add `larsen_1ubq_pm6_pdb` to `tests/testpaths.toml`.
- Add `TestEnvironment::Larsen1UbqPm6Pdb()`.
- Keep skip behavior when the configured path is empty/missing.

Blast radius: light. Test-only.

### G5. Test debug output paths should use a temp helper

Issue: some tests write debug artifacts directly to fixed `/tmp/...` paths.
This is harmless locally but not clean for repeatable test runs and parallel
CI.

Anchors:

- `tests/test_hm_welford.cpp:194-197` copies to
  `/tmp/hm_welford_inspect.h5`.
- `tests/test_tripeptide_neighbor_shielding.cpp:170` uses
  `/tmp/tripeptide_neighbor_smoke_out`.
- `tests/test_tripeptide_backbone_shielding.cpp:153` uses
  `/tmp/tripeptide_bb_smoke_out`.
- CLI parser tests use `/tmp/out`, `/tmp/pdbs`, `/tmp/npys` as literal command
  arguments at `tests/test_cli_parse.cpp:86`, `137`, `146`, `156`, `162`,
  `168-169`, `182`, `187`, `216`; those are contract-test literals and lower
  priority.

Likely fix:

- Add `TestEnvironment::TempPath(stem)` or local test utility using
  `std::filesystem::temp_directory_path()`.
- Keep intentional parser-literal tests as-is unless they write files.

Blast radius: light. Test-only.

### G6. Production trajectory filename convention is duplicated

Issue: `production.tpr`, `production.trr`, and `production.edr` are duplicated
in validation, parse errors, usage text, and execution. This is a small helper
candidate, not an architecture change.

Anchors:

- Usage text: `src/Cli/PrintUsage.cpp:39-40`
- Validation: `src/Cli/Validate.cpp:65-67`
- Parse diagnostic: `src/Cli/Parse.cpp:202`
- Execution path assembly: `src/nmr_extract.cpp:257-259`

Likely fix:

- Add a small typed helper near `Cli` or trajectory input code, e.g.
  `TrajectoryInputFiles::FromProductionDir(dir)`.
- Keep CLI semantics unchanged.

Blast radius: light. CLI trajectory mode only.

### G7. `atoms_category_info.npy` h5-reader comment is stale

Issue: h5-reader comment says 84 bytes, but the struct/static assert says 83.
The writer layout also sums to 83. This is documentation drift, not a runtime
bug.

Anchors:

- Stale comment: `h5-reader/src/io/QtNpyRecords.h:109`
- Static assert: `h5-reader/src/io/QtNpyRecords.h:168`
- Writer layout starts: `src/CategoryInfoProjection.cpp:313`

Likely fix:

- Change the comment to 83 bytes.
- Consider adding a one-line comment at the writer with the total size if it
  helps avoid another drift.

Blast radius: trivial. h5-reader comments only unless tests inspect comments.

### G8. Ring contribution exponential uses literal `4.0`

Issue: `ConformationResult` uses a hardcoded exponential length where
`PhysicalConstants.h` already defines `EXP_DECAY_LENGTH`.

Anchors:

- Literal use: `src/ConformationResult.cpp:86`
- Existing constant: `src/PhysicalConstants.h:89`
- Do not merge with graph-distance decay:
  `src/MolecularGraphResult.h:23-24` is bond-count BFS distance, not Angstroms.

Likely fix:

- Include `PhysicalConstants.h` in `ConformationResult.cpp`.
- Replace `/ 4.0` with `/ EXP_DECAY_LENGTH`.

Blast radius: trivial if value stays `4.0`.

### G9. Shared HDF5 frame dataset names repeat widely

Issue: `frame_indices`, `frame_times`, and `source_attached_per_frame` repeat
across many trajectory result writers and tests. This is a real schema surface,
but centralizing every H5 string may become a broad schema refactor.

Representative anchors:

- `src/BsShieldingTimeSeriesTrajectoryResult.cpp:170-171`
- `src/DihedralTimeSeriesTrajectoryResult.cpp:477-482`
- `src/BondedEnergyTimeSeriesTrajectoryResult.cpp:187-192`
- `src/AIMNet2EmbeddingTimeSeriesTrajectoryResult.cpp:142-152`
- Tests assert the same names, e.g. `tests/test_rmsd_tracking.cpp:176-178`,
  `tests/test_dihedral_time_series.cpp:223-225`,
  `tests/test_j_coupling_time_series.cpp:172`.

Likely fix:

- If touched, add small constants on `TrajectoryResult` such as
  `kFrameIndicesDataset`, `kFrameTimesDataset`, `kSourceAttachedDataset`.
- Do not centralize per-calculator group/dataset schemas in this cleanup pass.

Blast radius: medium if done repo-wide. Not first pass unless another H5 edit
is already underway.

## Track 2 - Runtime / Config Cleanup

### C1. AMBER charge conversion factor is duplicated

Issue: `18.2223` appears in C++ charge readers and the Python table generator.
This is a unit conversion constant. It should have one named home.

Anchors:

- `src/ChargeSource.h:151`, `src/ChargeSource.h:176-177`
- `src/AmberPreparedChargeSource.cpp:137`, comment at
  `src/AmberPreparedChargeSource.cpp:214`
- `tools/amber/generate_ff14sb_pb_table.py:21`,
  `tools/amber/generate_ff14sb_pb_table.py:250`

Likely fix:

- C++: move/alias to `PhysicalConstants.h`, e.g.
  `AMBER_PRMTOP_CHARGE_FACTOR`.
- Python: keep a local constant but update name/comment to match C++ and add a
  source comment. Cross-language sharing is not worth adding machinery for one
  scalar.

Blast radius: light. No numeric behavior change.

### C2. GROMACS nm/Angstrom conversions should be named

Issue: `0.1` and `10.0` encode Angstrom <-> nm conversion in several GROMACS
paths. `PhysicalConstants.h` already has a unit conversion section.

Anchors:

- `src/BondedEnergyResult.cpp:16`, use at `src/BondedEnergyResult.cpp:21`,
  `src/BondedEnergyResult.cpp:40`
- `src/GromacsFrameHandler.cpp:148`, `src/GromacsFrameHandler.cpp:196-198`
- `src/FullSystemReader.cpp:577-579`
- Existing unit conversion section: `src/PhysicalConstants.h:29-39`

Likely fix:

- Add `ANGSTROM_PER_NANOMETRE = 10.0` and
  `NANOMETRE_PER_ANGSTROM = 0.1` to `PhysicalConstants.h`.
- Replace local `A_TO_NM` and raw `* 10.0`.

Blast radius: light. No numeric behavior change.

### C3. Config-backed metadata writes hardcoded cutoff values

Issue: trajectory H5 metadata writes hardcoded values that already exist in
`CalculatorConfig`. If TOML changes, metadata lies.

Anchors:

- Water metadata literals:
  `src/WaterFieldTimeSeriesTrajectoryResult.cpp:154-158`
- H-bond count metadata literal:
  `src/HBondCountWelfordTrajectoryResult.cpp:229`, write at
  `src/HBondCountWelfordTrajectoryResult.cpp:248`
- Config defaults:
  `src/CalculatorConfig.cpp:56`,
  `src/CalculatorConfig.cpp:89-91`
- TOML keys:
  `data/calculator_params.toml:33`,
  `data/calculator_params.toml:72-74`

Likely fix:

- Write attributes from `CalculatorConfig::Get(...)`.
- Keep attribute names unchanged.

Blast radius: light. H5 metadata only unless tests assert exact attrs.

### C4. APBS model parameters are hardcoded calculator configuration

Issue: APBS grid and PB condition parameters are embedded in code. These are
model/condition choices, not universal constants. Existing config already has
the E-field clamp key, but APBS still uses `APBS_SANITY_LIMIT`.

Anchors:

- Grid dimension: `src/ApbsFieldResult.cpp:165-166`
- Fine/coarse padding: `src/ApbsFieldResult.cpp:163-164`,
  `src/ApbsFieldResult.cpp:179-180`
- Dielectrics/temperature/ionic strength:
  `src/ApbsFieldResult.cpp:184-187`
- Clamp uses physical constant:
  `src/ApbsFieldResult.cpp:258-259`
- Thermal conversion constant:
  `src/ApbsFieldResult.cpp:269-270`,
  `src/PhysicalConstants.h:62-64`
- Existing config clamp key:
  `src/CalculatorConfig.cpp:72`, `data/calculator_params.toml:49`
- Existing `APBS_SANITY_LIMIT` constant:
  `src/PhysicalConstants.h:95`

Likely fix:

- Add TOML/default keys:
  `apbs_grid_dim`, `apbs_fine_padding_A`, `apbs_fine_min_dim_A`,
  `apbs_coarse_padding_A`, `apbs_protein_dielectric`,
  `apbs_solvent_dielectric`, `apbs_temperature_K`,
  `apbs_ionic_strength_M`.
- Use `CalculatorConfig::Get("efield_magnitude_sanity_clamp")` for the clamp.
- Keep `KT_OVER_E_298K` in `PhysicalConstants.h` unless temperature becomes
  variable in the unit conversion. If APBS temperature is configurable, check
  whether the conversion should be `kT/e` from that configured temperature
  instead of fixed 298.15 K.

Blast radius: medium. Behavior unchanged with default TOML, but APBS numerical
output is sensitive. Needs APBS field/time-series tests.

### C5. Runtime binary path fallbacks still contain local machine paths

Issue: `RuntimeEnvironment` is mostly the right home, but it still falls back
to Jessica-specific micromamba paths for `mopac` and `tleap`. For sync to a
new machine, config/PATH should be the boundary.

Anchors:

- MOPAC fallback: `src/RuntimeEnvironment.cpp:155-160`
- `tleap` fallback: `src/RuntimeEnvironment.cpp:164-178`
- Config parsing for these keys: `src/RuntimeEnvironment.cpp:139-140`
- Tool accessors: `src/RuntimeEnvironment.cpp:283-284`

Likely fix:

- Remove hardcoded conda fallbacks, or move them to an untracked example
  config/user preset.
- Resolution order should be TOML -> env/PATH -> empty/loud skip, depending on
  tool criticality.

Blast radius: medium. May affect local convenience but improves portability.

### C6. `BuildFromPdb` bypasses RuntimeEnvironment temp-dir handling

Issue: `PdbFileReader.cpp` generates temp PDB paths under `/tmp` directly,
while `RuntimeEnvironment` owns `tmpdir` and `TempFilePath`.

Anchors:

- Direct `/tmp` path: `src/PdbFileReader.cpp:229`
- Runtime temp dir fallback: `src/RuntimeEnvironment.cpp:205`
- Runtime helper: `src/RuntimeEnvironment.cpp:276-282`

Likely fix:

- Replace local `UniqueTempPath()` or route it through
  `RuntimeEnvironment::TempFilePath()`.
- Check `BuildFromPdb` call paths: if `RuntimeEnvironment::Load()` is not
  guaranteed before `BuildFromPdb`, either load earlier or keep a small local
  fallback using `std::filesystem::temp_directory_path()`.

Blast radius: light to medium because `BuildFromPdb` is used in tests and
possibly before runtime setup.

### C7. Database DSN is redacted in C++, leaked in Python script

Issue: production C++ redacts libpq DSNs, but the orphan dump script prints
the raw DSN to stderr.

Anchors:

- C++ redaction helper: `src/TripeptideDftTable.cpp:257-303`
- C++ redacted log use: `src/TripeptideDftTable.cpp:312`
- Python config/DSN load: `scripts/dump_tripeptide_orderings.py:40-48`
- Python raw print: `scripts/dump_tripeptide_orderings.py:145-146`
- Script SQL: `scripts/dump_tripeptide_orderings.py:53-66`
- C++ SQL: `src/TripeptideDftTable.cpp:385-388`

Likely fix:

- Remove the Python DSN print or redact at least password/passfile/URI
  credentials.
- SQL duplication is not worth centralizing unless more DB scripts appear.

Blast radius: light. Script-only.

### C8. `batch_extract.sh` has hardcoded paths and stale CLI flags

Issue: script cannot sync to another machine and appears to call retired CLI
flags.

Anchors:

- Hardcoded consolidated root and binary:
  `scripts/batch_extract.sh:25-26`
- Current canonical CLI:
  `CLAUDE.md:41-44`,
  `src/Cli/PrintUsage.cpp:24-31`,
  `src/Cli/Parse.cpp:130-151`
- Stale mutant flags:
  `scripts/batch_extract.sh:134-140`
- Stale ORCA flags:
  `scripts/batch_extract.sh:153-156`

Likely fix:

- Make consolidated root and binary path CLI args or env/config.
- Update calls to the current `--orca --root NAME` and
  `--mutant --wt NAME --ala NAME` mode shape, or retire the script if the
  input layout no longer maps to root names.

Blast radius: light if script-only, but verify against real batch workflow.

### C9. Other script-local hardcoded paths

Issue: several helper scripts still contain local absolute paths. These are
not all production paths, but they are part of the "sync to new machine" sweep.

Anchors:

- `scripts/references/qwen_summary.py:20`
- `tools/molprobity_validate.py:47-52`
- `tools/amber/generate_ff14sb_pb_table.py:270`
- `scripts/larsen_hbond_grid_parse/pre_compute_dense_grids.py:66`

Likely fix:

- Use CLI args first, then env vars, then `shutil.which()`/relative repo paths.
- For `generate_ff14sb_pb_table.py`, default `--tleap` should come from
  `shutil.which("tleap")` or `$AMBERHOME/bin/tleap`, not a user path.

Blast radius: light. Script-only, but check generated table reproducibility.

### C10. `tests/testpaths.toml` is machine-local but includes in-tree files

Issue: `tests/testpaths.toml` is explicitly "for batcave" and contains many
absolute paths that could be relative to the repo. External data can remain
configured; in-tree fixtures should not need local absolute paths.

Anchors:

- Header says machine-local: `tests/testpaths.toml:1-2`
- In-tree fixtures as absolute paths:
  `tests/testpaths.toml:5`, `tests/testpaths.toml:8`,
  `tests/testpaths.toml:11`, `tests/testpaths.toml:15`,
  `tests/testpaths.toml:46`, `tests/testpaths.toml:60`,
  `tests/testpaths.toml:64`, `tests/testpaths.toml:67`
- External consolidated root: `tests/testpaths.toml:18`
- Loader currently treats values literally:
  `tests/TestEnvironment.cpp:34-91`

Likely fix:

- Teach `TestEnvironment` to resolve relative paths against repo/test root.
- Convert in-tree entries to relative paths.
- Keep external data paths configurable and skippable.

Blast radius: medium, test-only. Helps CMake/test portability.

## Track 3 - CMake / Build / Install Cleanup

Safety gate: do this track with unusual care, after the spinner backup. The
CMake work touches dependency discovery, runtime data lookup, RPATH, and install
semantics; small-looking edits can break local build/test launch behavior. Keep
changes incremental, preserve the existing local build until the replacement
path is verified, and test configure/build/test/install as a cadence before
folding in broader cleanup.

### B1. Dependency discovery uses source-level local defaults

Issue: root `CMakeLists.txt` has cache variables, but many defaults are local
absolute paths. New machines can override them, but a clean configure is not
self-explanatory and several missing optional deps still leak into includes or
links later.

Anchors:

- PostgreSQL explicit system paths:
  `CMakeLists.txt:55-71`
- OpenBabel explicit system paths:
  `CMakeLists.txt:90-97`
- GROMACS local paths:
  `CMakeLists.txt:102-111`
- Reduce local paths:
  `CMakeLists.txt:115-128`
- MOPAC local paths:
  `CMakeLists.txt:131-138`
- Torch Python path injected into `CMAKE_PREFIX_PATH`:
  `CMakeLists.txt:184-187`
- APBS/FETK include dirs:
  `CMakeLists.txt:201-204`
- HDF5/curl/z hardcoded link paths:
  `CMakeLists.txt:163-168`
- Target includes still consume those vars:
  `CMakeLists.txt:373-385`
- Target links still consume those vars:
  `CMakeLists.txt:407-412`

Likely fix:

- Split into `cmake/NmrDependencies.cmake` or small `FindNmr*.cmake` modules.
- Prefer `find_package`/`find_path`/`find_library` with `HINTS` from cache
  vars and environment, not local defaults in source.
- Keep cache variables so a local `CMakeUserPresets.json` can provide exact
  paths.

Blast radius: medium. Configure-only if done carefully; build impact if
optional deps are gated wrong.

### B2. Missing optional dependencies are warnings but still wired into build

Issue: CMake warns when GROMACS, Reduce, or MOPAC are missing, but the main
target still includes and links their paths. That gives a configure that
"succeeds" and a build that cannot.

Anchors:

- GROMACS warning: `CMakeLists.txt:108-111`; later include/link:
  `CMakeLists.txt:378-385`, `CMakeLists.txt:408`
- Reduce warning: `CMakeLists.txt:125-128`; later include/link:
  `CMakeLists.txt:375-377`, `CMakeLists.txt:409-411`
- MOPAC warning: `CMakeLists.txt:135-138`; later include/link:
  `CMakeLists.txt:374`, `CMakeLists.txt:412`

Likely fix:

- Add explicit options, e.g. `NMR_WITH_GROMACS`, `NMR_WITH_REDUCE`,
  `NMR_WITH_MOPAC`, `NMR_WITH_TORCH`, `NMR_WITH_APBS`.
- Either fail at configure when an enabled dependency is missing or compile out
  the dependent source files behind the option.
- Do not leave "warning but broken build" states.

Blast radius: medium to high. Requires knowing which modes/calculators are
allowed to be disabled. This is build architecture, but still within CMake
cleanup if options only gate existing surfaces.

### B3. Root presets exist but do not encode dependency profiles

Issue: root `CMakePresets.json` exists, but only has `default` and
`asan-ubsan`. It does not provide a clean local/user handoff for dependency
roots or install prefixes.

Anchors:

- Presets: `CMakePresets.json:4-47`
- `.gitignore` ignores build dirs but not `CMakeUserPresets.json`:
  `.gitignore:1-6`, `.gitignore:18-30`
- h5-reader already uses preset build dirs:
  `h5-reader/CMakePresets.json:10`

Likely fix:

- Keep committed presets generic.
- Add `.gitignore` entry for `CMakeUserPresets.json`.
- Add a committed `CMakeUserPresets.example.json` or docs showing local
  dependency roots.
- Add configure/build/test presets for normal release/debug if useful, but
  keep local paths out of committed JSON.

Blast radius: light. CMake UX only.

### B4. No install surface

Issue: no `install(...)` rules were found in the root CMake file. The project
can build from the tree, but it cannot complete a real configure/build/install
cadence.

Anchors:

- Library target starts: `CMakeLists.txt:223`
- CLI target: `CMakeLists.txt:438-439`
- Test infra starts: `CMakeLists.txt:470`
- `rg -n "install\\(" CMakeLists.txt` found no root install rules.

Likely fix:

- Include `GNUInstallDirs`.
- Install `nmr_extract` to `${CMAKE_INSTALL_BINDIR}`.
- If library installation is needed, install `nmr_shielding` archive/shared
  and public headers. If the program is the real artifact, start with CLI and
  runtime data only.
- Install runtime data needed by default config, especially
  `data/calculator_params.toml` and `data/ff14sb_params.dat`.
- Decide whether generated semantic tables are source-owned or install data
  (currently compiled into `src/generated/LegacyAmberSemanticTables.cpp`).

Blast radius: medium. Build/install only, but runtime data lookup must be
coherent.

### B5. Runtime data path is compiled to the source tree

Issue: `NMR_DATA_DIR` is defined as `${CMAKE_SOURCE_DIR}/data`. That works for
in-tree builds but not for installed binaries on a clean machine.

Anchors:

- Compile definition: `CMakeLists.txt:390-391`
- Runtime fallback uses the macro:
  `src/RuntimeEnvironment.cpp:183-193`
- Test data macro is similarly source-tree based:
  `CMakeLists.txt:551-552`

Likely fix:

- Add separate build-tree and install-tree defaults, e.g.
  `NMR_SOURCE_DATA_DIR` for tests/developer runs and `NMR_INSTALL_DATA_DIR` for
  installed binaries.
- Or remove compiled data path as an install dependency and require config/env
  for installed runs.
- Best practical path: compile default data dir to install prefix for install
  builds, and let `NMR_DATA_DIR` env/TOML override.

Blast radius: medium. Runtime behavior changes only outside the source tree if
defaults are preserved for build-tree runs.

### B6. RPATH and loader handling is local/system-specific

Issue: root CMake prepends `/usr/lib/x86_64-linux-gnu` to build and install
RPATH, turns off origin RPATH, and also injects `LD_LIBRARY_PATH` for tests.
This was a practical local fix for libpq/OpenSSL conflicts, but it is not a
portable install story.

Anchors:

- RPATH comments and values: `CMakeLists.txt:21-33`
- Test env `LD_LIBRARY_PATH`: `CMakeLists.txt:503-519`
- Torch path injection: `CMakeLists.txt:184-187`
- Compile flags include `-march=native -flto`:
  `CMakeLists.txt:389`

Likely fix:

- Keep local workaround available through cache/user presets.
- For install, prefer `$ORIGIN`/relative RPATH where dependencies are bundled
  or user-provided, and avoid hardcoding distro lib dirs in install RPATH.
- Move `-march=native` behind an option such as `NMR_NATIVE_OPTIMIZE`.

Blast radius: medium. Loader behavior is fragile; test with `ldd` and direct
binary launch after changes.

### B7. GTest fallback downloads from the network

Issue: if system GTest is missing, configure uses `FetchContent` and a GitHub
URL. That breaks offline/restricted syncs unless dependencies are already
vendored/cached.

Anchors:

- `find_package(GTest QUIET)`: `CMakeLists.txt:141`
- FetchContent block: `CMakeLists.txt:142-149`
- Test helper links GTest:
  `CMakeLists.txt:526-563`, `CMakeLists.txt:766-767`

Likely fix:

- Add an option such as `NMR_FETCH_GTEST` defaulting to `OFF` for reproducible
  configure, or document the network requirement in presets.
- Prefer system package or checked-out dependency path for sync-to-new-machine
  work.

Blast radius: light to medium. Test-only dependency, but configure behavior
changes.

### B8. Topology generator RDKit paths are local defaults

Issue: generator CMake uses a local micromamba RDKit root as the default. It is
opt-in, but should follow the same dependency-discovery/user-preset pattern as
the root build.

Anchors:

- Local default: `tools/topology/CMakeLists.txt:25-29`
- RDKit libraries assembled from root:
  `tools/topology/CMakeLists.txt:41-48`
- Runtime rpath to RDKit dir:
  `tools/topology/CMakeLists.txt:82-83`

Likely fix:

- Keep `RDKIT_ROOT` as a cache variable, but default from env or discovery
  rather than `/home/jessica/...`.
- Put local RDKit root in `CMakeUserPresets.json`.

Blast radius: light. Opt-in generator only.

## Suggested Implementation Order

1. CMake baseline:
   `CMakeUserPresets.json` ignore/example, dependency discovery cleanup,
   explicit optional-dependency gates.
2. Install/runtime data:
   install rules plus data-dir handoff; verify installed `nmr_extract` can find
   config/data or fails loudly with an actionable message.
3. Runtime/config literals:
   unit conversions, H5 metadata config reads, APBS config keys, runtime path
   fallbacks.
4. General light cleanup:
   JSON manifest, `SphericalTensor` packing, test config helpers, schema comment.
5. Script sweep:
   `batch_extract.sh`, DB script DSN print, helper script path defaults.

## Verification Targets

Minimum after CMake changes:

```bash
cmake --preset default
cmake --build --preset default
ctest --preset default
```

After install rules exist:

```bash
cmake --install build --prefix /tmp/nmr-install-smoke
/tmp/nmr-install-smoke/bin/nmr_extract --help
```

After runtime/config cleanup:

```bash
./build/unit_tests
./build/trajectory_tests
./build/water_field_tests
./build/hbond_count_welford_tests
./build/apbs_efield_time_series_tests
./build/apbs_efg_time_series_tests
```

Exact test binary names depend on the current CMake target set; prefer
`ctest --preset default` when possible.

## Do Not Chase In This Pass

- Generated topology tables unless a CMake/install step needs to move them.
- Test assertion tolerances and synthetic geometry numbers that encode the test
  itself.
- Tensor algebra dimensions (`3`, `5`, `9`, etc.) where the local math context
  is already explicit.
- Broad `TrajectoryResult` lifecycle abstraction. Prior read says this is not a
  light cleanup.
- TopologySidecar vs CategoryInfoProjection structured NPY duplication. Their
  fail-loud packed dtype mirrors are intentional unless a schema generator is
  designed separately.
