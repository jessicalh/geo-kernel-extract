# Next Session: MD Trajectory Loading + Loader Cleanup

Read memory and spec/INDEX.md first, as always. Then read this.

## What happened last session

See memory: `project_session_20260402_ui_and_builder.md`. Key outcomes:

- **CalculationRunner** (src/CalculationRunner.h): three tested cases
  (PDB only, single DFT, WT+ALA comparison). This is the "run stuff"
  entry point. The conformation accumulates results. No state on the
  runner.
- **Pipeline** (src/Pipeline.h): lower-level convenience, may fold
  into CalculationRunner.
- **Coulomb filter fix**: CoulombResult now uses KernelFilterSet like
  all other calculators. SelfSourceFilter, logged.
- **Protonation design history**: deep analysis of why reprotonation
  requires a new Protein (two indexing systems — symbolic vs geometric).
  Decision: protonate on load. See protonation/DESIGN_HISTORY.md.
- **Three architecture issues** identified. Two remain:
  1. ProtonationDetectionResult const_cast (fix by setting variant_index
     during construction)
  2. Everything loads as CrystalConformation (fix by adding the missing
     conformation types)

## What this session does

### 1. Add missing conformation types (fix architecture issue #2)

Add to ProteinConformation.h:
- **MDFrameConformation**: frame_index, time_picoseconds, walker_index,
  boltzmann_weight. Factory method: `protein.AddMDFrame(positions, metadata)`.
- **PredictionConformation**: tool (AlphaFold/OpenFold3/ESMFold),
  confidence_per_residue. Factory method: `protein.AddPrediction(positions, metadata)`.

These follow the CrystalConformation pattern exactly — subclass with
metadata fields, factory method on Protein that creates and owns it.
Ten lines each. The base class does all the work.

Fix OrcaRunLoader to use AddPrediction for AlphaFold structures
instead of AddCrystalConformation with NaN resolution.

### 2. Give the two loaders a parent

LoadProtein and LoadOrcaRun share a pattern: create Protein, add
atoms/residues, call FinalizeConstruction, set ProteinBuildContext,
create conformation. Abstract the shared steps so that adding a
third loader (GROMACS fleet) follows the same pattern.

This is NOT a class hierarchy. It's extracting the shared code so
it exists in one place. The loaders become thin wrappers that parse
their specific format and call the shared construction steps.

### 3. Load GROMACS fleet data

**Data location**: `/mnt/extfast/fleet_results/`
- 203 proteins, each with 5 walkers
- CHARMM36m force field, protonated structures
- Each protein directory has:
  - `best_pose.pdb` — protonated PDB of lowest free-energy frame
    (CHARMM naming: H1/H2/H3 N-term, HSD/HSE histidines)
  - `harvest_receipt.json` — metadata (protein_id, n_atoms, timing,
    fes_minimum_frame with walker/time_ps/rmsd/rg)
  - `run_params/prod.tpr` — GROMACS topology (charges via GmxTprChargeSource)
  - `walker_0/` through `walker_4/` — each with:
    - `md.xtc` — trajectory (XTC format, needs GROMACS or MDAnalysis to read)
    - `md.tpr` — per-walker topology
    - `COLVAR`, `HILLS` — PLUMED OPES metadynamics data
    - `md.gro` — final frame (GRO format)

- `sampled_poses/` subdirectory (at /mnt/extfast/fleet_results/sampled_poses/):
  - 182 proteins with 10 Boltzmann-reweighted poses each
  - `pose_001.pdb` through `pose_010.pdb` — protonated PDBs (same
    atom count as best_pose, CHARMM naming)
  - `ensemble.json` — metadata per pose: walker, time_ps, rmsd_nm,
    rg_nm, rbias_kJmol, weight (Boltzmann weight from OPES
    metadynamics reweighting)
  - Temperature: 308 K, method: boltzmann_reweighted_metad

**Loading strategy**:
- Start with best_pose.pdb (one frame per protein). This exercises
  the loader, NamingRegistry (CHARMM names), charges
  (GmxTprChargeSource from prod.tpr), and all calculators.
- Then: load the 10 sampled poses as MDFrameConformations on the
  SAME protein. Each pose is a PDB file (no XTC reading needed).
  ensemble.json provides frame metadata (walker, time, weight).
  One protein, 10 conformations, 10 Pipeline runs. This is the
  ensemble analysis path the refinement model needs.
- Later: XTC trajectory loading for full trajectory analysis.

**Naming**: CHARMM36m uses different names from PDB/AMBER:
- HSD/HSE/HSP (not HID/HIE/HIP)
- HN (not H) for backbone amide
- HB1/HB2 (not HB2/HB3) for beta methylenes
- NamingRegistry already handles all of these. GmxTprChargeSource
  already tested (test_protonation_pipeline.cpp GmxChargeTest).

### 4. Add a CalculationRunner case for fleet data

```cpp
// Case 4: GROMACS fleet best-pose analysis
static RunResult RunFleetBestPose(
    ProteinConformation& conf,
    ChargeSource& charges);
```

Same as RunSingleDft but no ORCA comparison (fleet data has no DFT).

### 5. Test on real fleet data

Load one protein (1A6J_5789), run all calculators, verify
results are nonzero and physically reasonable. Then batch test
on all 203 proteins.

## Constraints

- Do NOT touch calculator code. The calculators work. 277 tests.
- Do NOT change ConformationResult, ConformationAtom, or the
  singleton/dependency framework.
- Talk before changing Protein.h — it's the core object.
- NamingRegistry is read-only data. If CHARMM names are missing,
  add entries. Don't restructure the registry.
- XTC reading will need an external library (GROMACS libgmx or
  MDAnalysis via Python). For this session, start with best_pose.pdb
  (just a PDB file). XTC loading is a follow-up.

## What success looks like

- MDFrameConformation and PredictionConformation exist and are tested
- OrcaRunLoader uses AddPrediction (no more CrystalConformation lie)
- A fleet protein loads from best_pose.pdb with CHARMM charges
- All 8 calculators produce nonzero results on fleet data
- 203 fleet proteins batch-validated
- All existing 277 tests still pass
