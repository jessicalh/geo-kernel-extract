# Session 2026-04-03: Fleet Loader + Loading Architecture

## What was done

### Conformation type hierarchy
- PredictionConformation, MDFrameConformation added to ProteinConformation.h/.cpp
- Factory methods on Protein: AddPrediction(), AddMDFrame()
- Typed index vectors (prediction_indices_, md_frame_indices_)

### Orthogonal conformation access
- `Protein::Conformation()` — returns primary conformation as ProteinConformation&
- `Protein::AddConformation()` — base factory for unknown provenance
- Typed accessors (MDFrameAt, PredictionAt, CrystalConf) for metadata consumers
- All 198 test calls changed from CrystalConf() to Conformation()

### Loader honesty
- OrcaRunLoader: AddCrystalConformation → AddPrediction("AlphaFold+tleap")
- PdbFileReader: AddCrystalConformation(0,0,0) → AddConformation() — no fake crystal metadata

### Fleet loader (GromacsEnsembleLoader)
- Reads TPR via libgromacs (not text dump parsing)
- TPR is single authority for: atom names, residue names, charges, elements, residue boundaries
- NamingRegistry translates CHARMM→canonical with declared ToolContext::Charmm
- HIS variant detection: CHARMM name if explicit (HSP), else atom name inspection (HD1/HE2)
- CovalentTopology from first pose geometry
- 10 MDFrameConformations with walker, time_ps, weight, rmsd_nm, rg_nm
- PreloadedChargeSource wraps TPR charges for ChargeAssignmentResult
- Net charge computed from TPR charge sum
- RunAllFrames() calls Pipeline::RunClassicalCalculators per frame

### libgromacs integration
- Linked against GROMACS 2026.0 at /home/jessicalh/gromacs/
- Internal headers from source tree (fileio/tpxio.h, topology/topology.h, etc.)
- CMakeLists.txt updated with GROMACS_SRC, GROMACS_BUILD, GROMACS_LIB

### Test infrastructure
- In-repo test data: tests/data/fleet/{1A6J_5789,1AEP_4814}/
- 11 fleet loader tests including FullPipelineAllFrames
- Full pipeline verified: 14 results per frame × 10 frames = 290 NPY arrays
- APBS runs (1.3s per frame), all 8 calculators run, DSSP runs
- 283 tests pass, 3 skipped (KaML), baselines byte-identical

## CRITICAL: Silent tool failures

The Pipeline logs xTB/APBS failures and continues. No batch test verifies
these results exist. This means:

- xTB has never been validated on the 723-protein batch
- APBS has never been validated on the 723-protein batch
- The 4 batch test files each have their own manual AttachResult chain
  that skips xTB and APBS entirely
- Only test_full_pipeline.cpp checks HasResult<XtbChargeResult>() and
  HasResult<ApbsFieldResult>(), and it runs on one protein

The next session must:
1. Make the batch tests use Pipeline (one path, not four copies)
2. Add EXPECT_TRUE(HasResult<ApbsFieldResult>()) to batch validation
3. Decide whether xTB failure is skip-with-warning or hard fail
4. Investigate xTB on large proteins (2480 atoms core dumps, 335 works).
   128GB RAM is not nothing. Research GFN2-xTB memory scaling — it may
   be O(N^2) or O(N^3) in basis functions. 2480 atoms = 6167 basis
   functions; the SCF matrix is 6167^2 ≈ 300MB, but intermediates may
   be much larger. Check if xTB has a memory limit flag or if GFN-FF
   is an alternative for large systems.

## What is NOT done — next session MUST address

### 1. xTB size limit (real constraint)
xTB core dumps on 2480 atoms (1A6J). Works on 335 atoms (1B1V, 58s).
GFN2-xTB has a practical limit around 500-1000 atoms.

Options:
- Fragment-based xTB (compute per-residue or per-fragment)
- GFN-FF instead of GFN2 (faster, less accurate, may handle larger systems)
- Skip xTB for fleet proteins >1000 atoms (honest about what we can't do)
- The ORCA proteins are 300-600 atoms and work fine

### 2. Net charge plumbing
Net charge is needed by xTB (--chrg) and potentially APBS.
Currently:
- Fleet: computed from TPR charge sum in FleetLoadResult.net_charge
- ORCA: available from prmtop charges or ProtonationState::NetChargeForProtein
- PDB-only: no charges, no net charge (xTB will use 0, which may be wrong)

What's needed: PipelineOptions should carry net_charge. The caller computes it
from whatever charge source they have. This is done for xtb_net_charge but
should be verified for APBS.

### 3. CalculationRunner vs Pipeline unification
Two parallel paths exist:
- CalculationRunner::RunSingleDft — foundation + charges + 8 calculators (NO APBS, NO xTB)
- Pipeline::RunClassicalCalculators — foundation + charges + xTB + APBS + 8 calculators

The fleet loader correctly uses Pipeline. The ORCA batch tests still use
CalculationRunner. The batch runner (for 723 proteins) should use Pipeline too.
CalculationRunner can remain for specific test scenarios that don't need
APBS/xTB, but production code should use Pipeline.

### 4. Feature output naming for ensembles
Current: one directory per frame (frame_001/, frame_002/, ...) with ~29 NPY each.
Desired: single directory per protein with a manifest JSON mapping filenames to
frame metadata (walker, time_ps, weight). Consider stacked NPY arrays
(shape [n_frames, n_atoms, ...]) instead of per-frame files — 700 proteins ×
10 frames = 7000 directories is unacceptable.

### 5. xTB on ORCA path
The ORCA batch runner (723 proteins × 2 variants = 1446 runs) uses
CalculationRunner::RunSingleDft which skips xTB. If xTB is needed for the
model, the batch runner needs to switch to Pipeline. The ORCA proteins are
small enough (300-600 atoms) that xTB should work.

## Architecture notes for the next session

### The three loading pipelines share three layers:
```
Layer 1: Read external format → canonical atoms + residues + variant_index
Layer 2: CacheResidueBackboneIndices + DetectAromaticRings (symbolic topology)
Layer 3: CovalentTopology::Resolve (geometric topology)
```
Layers 2-3 are fused in FinalizeConstruction. Loaders differ only in Layer 1.

### ChargeSource hierarchy:
```
ChargeSource (abstract)
├── ParamFileChargeSource   — ff14SB from flat parameter file
├── PrmtopChargeSource      — ff14SB/ff19SB from AMBER prmtop
├── PreloadedChargeSource   — pre-extracted charges (e.g., from TPR via libgromacs)
├── GmxTprChargeSource      — CHARMM36m from GROMACS .tpr (calls gmx dump — DEPRECATED by PreloadedChargeSource)
└── StubChargeSource        — uniform test charges
```

### ForceField and ToolContext are DECLARED, not inferred:
The caller knows what data they're loading. The loader is TOLD "this is CHARMM36m"
or "this is ff14SB". It does not look at the charges and guess.
