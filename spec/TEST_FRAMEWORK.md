# Test Framework

What we actually have, what each test does, and what to run when.


## Test Executables

### By cost tier

| Executable | Files | Runtime | What it covers |
|------------|-------|---------|----------------|
| `unit_tests` | 8 | <5s | Pure logic: spherical tensor, amino acid, ring hierarchy, object model, naming, config. No file I/O. |
| `structure_tests` | 24 | ~5 min | Static PDB path: 1UBQ + P84477. All classical calculators, APBS, protonation, DSSP, write features. No MOPAC. |
| `trajectory_tests` | 1 | ~5 min | GROMACS fleet loader: TPR reading, multi-pose PDB, charge extraction. |
| `smoke_tests` | 1 | ~80s | End-to-end pipeline: NoDft (1UBQ) + WithDft (P84477). NPY write, log validation, binary comparison. |
| `gromacs_streaming_tests` | 1 | ~5 min | Trajectory streaming: XTC reading, PBC fix, AIMNet2 + H5 write. Critical path for H5 and AIMNet2. |
| `water_field_tests` | 1 | ~5 min | Explicit solvent calculators on full-system trajectory. |
| `job_spec_tests` | 1 | ~10 min | JobSpec 6-mode parser and execution. |
| `fleet_smoke_tests` | 1 | ~40 min | Fleet extraction: 1A6J_5789, 10 GROMACS poses, CHARMM36m charges, full pipeline per pose. |
| `fes_fleet_smoke_tests` | 1 | ~15 min | FES-sampler naming: 1I8X_4351 + 1HD6_4820, 6 descriptive-named poses each. |
| `mopac_tests` | 3 | ~1 hr | MOPAC semiempirical: test_mopac_result, test_full_pipeline, test_mutation_delta. |
| `batch_tests` | 4 | ~2 hr | Multi-protein sweep across consolidated/: all classical calculators on many proteins. |

### By use case

**Static structure** (--pdb, --orca, --mutant):
`unit_tests`, `structure_tests`, `smoke_tests`, `mopac_tests`, `batch_tests`

**Trajectory analysis** (--trajectory):
`trajectory_tests`, `gromacs_streaming_tests`, `water_field_tests`,
`fleet_smoke_tests`, `fes_fleet_smoke_tests`, `job_spec_tests`

```
# Always run (sub-second)
cd build && ./unit_tests

# Static structure regression (~5 min)
./structure_tests

# Trajectory development (~5 min each)
./trajectory_tests
./gromacs_streaming_tests    # H5 + AIMNet2 critical path
./water_field_tests

# Pipeline confidence (~80s, includes MOPAC)
./smoke_tests

# Before production batch runs (~40 min)
./fleet_smoke_tests

# Run everything
ctest
```


## What to Run When

This is the decision tree. Before starting a session, identify the change type and
run the matching test tier.

### Tier 0: Unit tests (always, <5s)

Run `./unit_tests`. Fast sanity check on core logic — run before anything else.

### Tier 1: Smoke only (~80 seconds)

Run `./smoke_tests`. This covers:

- **Changed a parameter** in calibration.toml or CalculatorConfig
- **Changed an extractor's computation** (e.g. fixed a sign, added a filter)
- **Added a new field** to an existing extractor's WriteFeatures
- **Changed OperationRunner** sequencing or dependency logic

The binary comparison will catch any byte-level change in the 46-49 NPY output
files. If you intentionally changed output, re-bless the baseline (see below).

Binary comparison breaks when you add a new extractor or add new NPY output
files. That is expected. The comparison reports new files separately from
changed files, so you can distinguish "new output I intended" from "existing
output that changed when it should not have." Re-bless after confirming the
new output is correct.

### Tier 2: Structure tests (~5 min)

Run `./structure_tests`. This replaces the old `--gtest_filter` incantations.

- **Added a new extractor**: smoke verifies pipeline integration, structure_tests
  verifies analytical correctness for all classical calculators.
- **Changed PDB loading or topology**
- **Changed protonation, DSSP, HBond, charge assignment**
- **Changed any classical calculator**

### Tier 2.5: Trajectory tests (~5 min each)

Run `./trajectory_tests`, `./gromacs_streaming_tests`, `./water_field_tests`.

- **Changed GromacsProtein, GromacsFrameHandler, or accumulators**
- **Changed AIMNet2 or H5 output** — `gromacs_streaming_tests` is the critical path
- **Changed WaterField, HydrationShell, or HydrationGeometry** — `water_field_tests`
- **Changed fleet loader or TPR reading** — `trajectory_tests`

### Tier 2.5b: Fleet smokes (~40 min / ~15 min)

Run `./fleet_smoke_tests` and/or `./fes_fleet_smoke_tests`.

- **Before launching a production batch** (700 proteins x 10 poses)
- **Changed CHARMM36m charge handling** or force field paths
- **Changed pose PDB filename resolution** (fes_fleet)

### Tier 3: Full suite

Run `ctest` or run each target individually. Required for:

- **Changed fundamental data structures** (Atom, Ring, Protein, ProteinConformation)
- **Before merging to master** after any non-trivial change
- **When in doubt**

Individual slow targets:
- `./mopac_tests` (~1 hr) — MOPAC semiempirical, mutation delta, full pipeline
- `./batch_tests` (~2 hr) — multi-protein sweep across consolidated/


## Smoke Test Design

`test_smoke.cpp` runs two configurations through the identical validation pipeline:

### SmokeNoDft (1UBQ, ~70 seconds)
- Loads `tests/data/1ubq_protonated.pdb` (1231 atoms, 76 residues, 4 rings)
- Charges from ff14SB parameter file
- No ORCA DFT path
- 17 results attached, 46 NPY files written

### SmokeWithDft (P84477, ~7 seconds)
- Loads from `consolidated/P84477/` (351 atoms, ORCA DFT run)
- Charges from prmtop
- ORCA NMR shielding tensors loaded
- 18 results attached, 49 NPY files written

### SmokeFleet (1A6J_5789, ~40 min) — separate executable: `fleet_smoke_tests`
- Loads from `tests/data/fleet/1A6J_5789/` via BuildFromGromacs
- 2480 atoms, 10 GROMACS poses, CHARMM36m charges from TPR
- `RunAllFrames` runs the full pipeline on all 10 poses independently
- 13+ results per frame, 35+ NPY files per frame (10 frame directories)
- Binary comparison on frame_001 only (blessed at `tests/golden/blessed/fleet/frame_001/`)

### Validation pipeline (all configs)

1. **Pipeline**: `OperationRunner::Run` — all extractors in dependency order.
   Asserts minimum result count and no error return.

2. **Feature write**: `ConformationResult::WriteAllFeatures` — every attached
   result writes its NPY contribution. Asserts minimum file count.

3. **Log validation**: `OperationLog::ConfigureFile` captures JSON-lines to
   `log.jsonl` in the output directory. After the run, the test:
   - Asserts zero ERROR entries
   - Verifies every attached result appears in the log
   - Checks all `[BEGIN]` scopes have matching `[END]`

4. **NPY validation**: counts files, checks non-zero size, verifies NPY
   magic bytes in every file.

5. **Binary comparison**: if a blessed baseline exists at
   `tests/golden/blessed/{nodft,withdft}/`, byte-compares every NPY file.
   Reports identical/different counts and fails on any difference.

### Output directory structure

```
tests/golden/
  blessed/              # known-good baseline (re-bless after intentional changes)
    nodft/              # 46 .npy + log.jsonl
    withdft/            # 49 .npy + log.jsonl
  smoke/
    2026-04-09_082931/  # timestamped runs (accumulate, prune manually)
      nodft/
      withdft/
```

### Blessing a new baseline

After an intentional change to extractor output:

```bash
# Run the smoke test (writes to a new timestamped directory)
./smoke_tests

# Verify the Python validation agrees the output is correct
python3 ../tests/validate_smoke.py --latest

# Bless (replace the old baseline with the new run)
cp -r tests/golden/smoke/<timestamp>/nodft tests/golden/blessed/nodft
cp -r tests/golden/smoke/<timestamp>/withdft tests/golden/blessed/withdft
```

### Python validation

`tests/validate_smoke.py` provides deeper semantic validation:

```bash
python3 tests/validate_smoke.py --latest          # validate most recent run
python3 tests/validate_smoke.py <dir>              # validate specific directory
python3 tests/validate_smoke.py <dir> --blessed <dir>  # compare two directories
```

It loads every .npy via numpy and checks:
- Shape consistency (all per-atom arrays have shape[0] == atom count)
- No NaN or Inf values
- Physical bounds (charges, positions)
- Tensor shapes (shielding arrays are (N, 9))
- Ring contribution column count (59)
- SDK typed load (valency, dipole, cos_phi/sin_phi validation)
- Log timing breakdown per calculator


## Complete Test File Inventory

### Unit tests (no external data, fast)

| File | Tests |
|------|-------|
| test_spherical_tensor.cpp | SphericalTensor decomposition/reconstruction roundtrip |
| test_amino_acid.cpp | AminoAcidType enum roundtrip, aromaticity detection |
| test_ring_hierarchy.cpp | Ring type intensity, lobe offset, aromaticity properties |
| test_object_model.cpp | Protein/Residue/Atom hierarchy construction |
| test_naming_registry.cpp | Amino acid naming and protonation state resolution |
| test_iupac_atom_identity.cpp | IUPAC atom lookup by residue and atom name |
| test_runtime_environment.cpp | ff14sb_params, tmpdir, MOPAC library availability |
| test_calculator_config.cpp | Ring intensity and JB lobe offset config defaults |
| test_biot_savart_result.cpp | Biot-Savart field from wire segments (analytical, also 1UBQ integration) |
| test_coulomb_result.cpp | Coulomb E-field and EFG from point charges (analytical, also 1UBQ integration) |
| test_mcconnell_result.cpp | McConnell dipolar kernel and EFG (analytical, also 1UBQ integration) |
| test_ring_susceptibility_result.cpp | Ring susceptibility tensor Trace(M/r^3)/3 identity (analytical, also 1UBQ integration) |
| test_pi_quadrupole_result.cpp | Pi-quadrupole EFG tensor symmetry and tracelessness (analytical, also 1UBQ integration) |
| test_dispersion_result.cpp | Dispersion van der Waals tensor (analytical, also 1UBQ integration) |

### Integration tests — 1UBQ (load protein, run calculators)

| File | Tests |
|------|-------|
| test_pdb_loading.cpp | PDB parsing: residue count, atom count, ring count, bonds, elements |
| test_geometry_result.cpp | Geometry computation and ring coordinate frames |
| test_demo_result.cpp | DemoResult dependency on GeometryResult |
| test_traversal_dump.cpp | Full typed object model traversal (every atom, bond, ring, residue) |
| test_two_conformations.cpp | Multiple independent conformations on one protein |
| test_dssp_result.cpp | DSSP secondary structure assignment (helix in residues 23-34, beta strands) |
| test_protonation_detection.cpp | Titratable residue protonation state detection |
| test_atom_flat.cpp | Atom parent_atom_index hierarchy (also has unit test portion) |
| test_foundation_results.cpp | GeometryResult, ChargeAssignmentResult, EnrichmentResult integration |
| test_apbs_field_result.cpp | APBS electrostatic field computation |
| test_apbs_wired.cpp | APBS bridge with vacuum Coulomb fallback |
| test_apbs_ff14sb.cpp | APBS with ff14SB force field charges |
| test_mopac_result.cpp | MOPAC PM7+MOZYME semiempirical charges (1231 atoms, ~60s) |
| test_full_pipeline.cpp | All 9 foundation results attached in dependency order |
| test_haigh_mallion_result.cpp | Haigh-Mallion ring current contact term |
| test_calculation_runner.cpp | OperationRunner sequencing on 1UBQ (also uses P84477) |

### Integration tests — P84477 / ORCA (DFT pipeline)

| File | Tests |
|------|-------|
| test_protonation_pipeline.cpp | Full protonation pipeline with PROPKA and charge assignment |
| test_hbond_result.cpp | Hydrogen bond detection and shielding contribution |
| test_pipeline_and_sample.cpp | OperationRunner + viewer library contract (SampleAt, geometry access) |
| test_write_features.cpp | Full pipeline NPY output to baseline_features/P84477/ |
| test_mutation_delta.cpp | WT vs ALA mutation effect prediction |

### Batch tests — multiple proteins from consolidated

| File | Tests |
|------|-------|
| test_batch_mcconnell.cpp | McConnell across multiple consolidated proteins |
| test_batch_coulomb_ringchi.cpp | Coulomb and ring chi across multiple proteins |
| test_batch_piquad_disp.cpp | Pi-quadrupole and dispersion across multiple proteins |
| test_batch_biot_savart_haigh_mallion.cpp | BS and HM across multiple proteins |

### Special

| File | Tests |
|------|-------|
| test_fleet_loader.cpp | GROMACS ensemble loader: TPR reading, multi-pose PDB, 10-pose batch |
| test_job_spec.cpp | JobSpec 5-mode parser and execution (separate executable) |
| test_smoke.cpp | End-to-end smoke test: NoDft + WithDft pipeline, NPY write, log validation, binary comparison (separate executable: `smoke_tests`) |
| test_smoke_fleet.cpp | Fleet extraction smoke: 1A6J_5789 x 10 poses, full pipeline per pose (separate executable: `fleet_smoke_tests`, ~40 min) |


## Pipeline Result Order

This is the sequence that `OperationRunner::Run` executes. The smoke test
validates that all of these attach in order and produce the expected NPY output.

### Tier 0 — Foundation
1. GeometryResult (ring/bond geometry, coordinate frames)
2. SpatialIndexResult (nearest-neighbour lookup)
3. EnrichmentResult (atom roles, residue classification)
4. DsspResult (secondary structure, requires cifpp/dssp library)
5. ChargeAssignmentResult (if charge source provided)

### Tier 0.5 — External tools (require charges)
6. MopacResult (PM7+MOZYME semiempirical, ~60s on 1UBQ)
7. ApbsFieldResult (Poisson-Boltzmann solvated E-field, ~30s on 1UBQ)

### Tier 1 — Classical calculators
8. BiotSavartResult (ring current B-field and shielding)
9. HaighMallionResult (ring current contact terms)
10. McConnellResult (bond magnetic anisotropy)
11. RingSusceptibilityResult (diamagnetic susceptibility shielding)
12. PiQuadrupoleResult (pi-electron EFG)
13. DispersionResult (van der Waals shielding)
14. CoulombResult (if charges: electrostatic E-field and EFG)
15. MopacCoulombResult (if MOPAC: Coulomb from semiempirical charges)
16. MopacMcConnellResult (if MOPAC: McConnell from semiempirical bond orders)
17. HBondResult (if DSSP: hydrogen bond shielding)
18. SasaResult (Shrake-Rupley per-atom SASA)
19. AIMNet2Result (if AIMNet2 model: neural network charges + Coulomb EFG)

### Tier 2 — DFT comparison (optional)
20. OrcaShieldingResult (ORCA DFT shielding tensors, if path provided)


## NPY Output Files

The smoke test writes and validates these files. The number is the minimum
expected count for each configuration.

### Identity arrays (4 files, always written)
- `pos.npy` — (N, 3) float64, atom positions in angstroms
- `element.npy` — (N,) int32, atomic numbers
- `residue_index.npy` — (N,) int32, residue index per atom
- `residue_type.npy` — (N,) int32, AminoAcidType enum per atom

### Structural arrays (2 files)
- `ring_contributions.npy` — (P, 59) float64, per-(atom,ring) pair features
- `ring_geometry.npy` — (R, 10) float64, ring center/normal/radius

### Per-calculator arrays (varies)
Each calculator writes 1-5 arrays: a `_shielding.npy` (N, 9) spherical tensor,
plus calculator-specific scalars, per-type breakdowns, or E-field vectors.

### DFT arrays (3 files, WithDft only)
- `orca_diamagnetic.npy`, `orca_paramagnetic.npy`, `orca_total.npy`
