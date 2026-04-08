# learn/ — Calibration Model for NMR Shielding Tensors

**NEVER USE SYMLINKS.** Copy data or point at the real path. Symlink trees
mix old and new extractions and break the SDK's required-file validation.

This is Goal 1: calibrate classical geometric kernels against DFT using
WT-ALA mutation deltas.  NOT a GNN.  NOT the upstream RefDB predictor.

## Current pipeline

```
calibration/populate.py     → symlink tree from consolidated/
learn/extract.py            → nmr_extract --mutant → calibration/features/{run}/{ID}/
learn/c_equivariant/        → training pipeline
```

### Data flow

1. **Extract**: `extract.py --run CalibrationExtractionTest` calls
   `nmr_extract --mutant --wt ROOT --ala ROOT` for each protein.
   Produces 51 NPY arrays per protein in a flat directory.  C++ computes
   the WT-ALA delta (delta_shielding, delta_scalars, delta_ring_proximity).

2. **Load**: `protein.load_protein(path)` returns a typed `ProteinFeatures`
   with physics-grouped access.  All 51 arrays are wrapped in typed classes
   from `features.py` (ShieldingTensor, EFGTensor, PerRingTypeT2, etc.).
   Shape-validated on construction.  Fails on unregistered files.

3. **Dataset**: `c_equivariant/dataset.py` builds `CalibrationDataset`:
   - 40 kernel T2 vectors (L=2): BS/HM/Disp per-ring-type + MC/MopacMC
     bond-categories + calculator totals.  PQ dropped (r < 0.03).
   - 78 scalar features (L=0): element, residue type, ring proximity,
     McConnell/Coulomb/HBond scalars, MOPAC charge, T1 magnitudes,
     per-ring cylindrical geometry, mutation identity, max bond order.
   - Per-protein kernel normalization (angular structure preserved).
   - Scalar z-score (train stats shared to val).

4. **Model**: `c_equivariant/model.py` — KernelMixingHead (scalar MLP →
   per-kernel weights × L=2 kernels) + optional EquivariantCorrectionHead
   (e3nn tensor products on 9 individual L=2 inputs).

5. **Train**: `c_equivariant/train.py` — per-component R² (m=-2..+2),
   distance-stratified R² (0-4/4-8/8-12/12+ Å), per-protein R²,
   naive physics baselines (ridge, unweighted sums), kernel weight logging.
   `--distance-weighted` flag for exp(-d/8Å) loss weighting.

6. **Analyze**: `c_equivariant/analyze.py` — cold analytical diagnostic.
   Per-kernel correlation, per-element breakdown, per-distance ridge R²,
   per-protein ridge R², residual structure.  No learning.

## Key physics constraints

- **T2 is sacred**: never collapse to scalar.  The angular residual IS the thesis.
- **Per-protein kernel normalization**: isolates angular structure, strips
  magnitude.  The calibration question is about direction, not scale.
- **Mechanical mutants**: backbone doesn't change.  No DSSP/dihedral features.
- **WT features only**: ring kernels ARE the delta (ALA has no rings).
  Non-ring kernels provide environment context for the MLP.
- **Equivariance**: scalar × L=2 = L=2.  T1 used as magnitude only (rotation-invariant).

## What NOT to do

- Do not use `load.py` — it's deprecated (old wt/ala subdirectory layout).
- Do not use `delta.py` for atom matching — C++ handles this now.
- Do not add graph structure or message passing — this is calibration,
  not the upstream GNN.  Per-atom independence is the assumption.
- Do not add DSSP/backbone features — see "mechanical mutants" above.
- Do not use global kernel normalization — per-protein is correct for
  the calibration angular-structure question.
- Do not treat this as a large-scale ML problem.  48 physics-derived
  kernels, ~150K atoms.  Every feature should encode physics.

## Files

### Current (use these)
- `features.py` — typed wrappers: SphericalTensor, ShieldingTensor, etc.
- `protein.py` — feature registry + load_protein() → ProteinFeatures
- `extract.py` — batch extraction to calibration/features/
- `c_equivariant/dataset.py` — CalibrationDataset, compute_r2, baselines
- `c_equivariant/model.py` — KernelMixingHead + EquivariantCorrectionHead
- `c_equivariant/train.py` — training loop with full diagnostics
- `c_equivariant/analyze.py` — cold analytical diagnostic

### Deprecated (do not use for new code)
- `load.py` — old loader, assumes wt/ala subdirs
- `delta.py` — old Python atom matching, replaced by C++ delta
- `status.py` — old diagnostic, uses deprecated loader
- `t2_residual.py` — old analysis, uses deprecated loader

### Reference
- `JOURNAL.md` — learning progress record
- `ARGUMENT.md` — thesis argument for validation approach
- `CALIBRATION_CHECKLIST.md` — parameter validation plan
- `bones/` — historical session notes and exploratory work
