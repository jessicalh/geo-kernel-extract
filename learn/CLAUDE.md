# learn/ — Calibration Model for NMR Shielding Tensors

**NEVER USE SYMLINKS.** Copy data or point at the real path.

This is Goal 1: calibrate classical geometric kernels against DFT using
WT-ALA mutation deltas.  NOT a GNN.  NOT the upstream RefDB predictor.

## Current pipeline

All configuration in `src/calibration.toml`.  All data access via `nmr_extract` SDK.

```
calibration/{ID}/               → real prmtop, xyz, nmr.out per protein
learn/extract.py                → nmr_extract --mutant → calibration/features/{run}/{ID}/
src/mutation_set/               → training pipeline (config-driven, SDK-based)
```

### Data flow

1. **Extract**: `extract.py` calls `nmr_extract --mutant` for each protein.
   Produces 53 NPY arrays per protein including `ring_contributions.npy`
   (per-ring T2 tensors) and `ring_geometry.npy`.

2. **Load**: `nmr_extract.load(path)` returns a typed `Protein` with
   every array wrapped in SDK types (ShieldingTensor, EFGTensor,
   RingContributions, etc.).  Named properties throughout — no index
   arithmetic.  e3nn Irreps on every tensor.

3. **Dataset**: `src/mutation_set/dataset.py` orchestrates:
   - `kernels.py` assembles 46 T2 kernels via SDK `.as_block()` accessors.
     Named KernelLayout registry — no hardcoded indices.
   - `scalars.py` assembles 190 scalar features as named ScalarBlocks.
     Categorical columns tracked by name, not position.
   - Per-protein kernel normalization + kernel scale factors as scalars.
   - Scalar z-score (train stats shared to val).

4. **Model**: `src/mutation_set/model.py` — gated KernelMixingHead
   (MLP → per-kernel weights, gated by kernel self-reported magnitude,
   scaled by sqrt(n_kernels)) + optional EquivariantCorrectionHead.

5. **Train**: `src/mutation_set/train.py` — GPU-resident loop (no
   DataLoader overhead), cosine annealing, per-component and distance-
   stratified R², per-protein R², correction_scale logging.
   Config-driven via calibration.toml with CLI overrides.

6. **Analyze**: `src/mutation_set/analyze.py` — cold analytical diagnostic.
   7 sections: kernel correlation, per-element breakdown, distance signal,
   per-protein ridge, greedy forward selection, leave-one-out ablation,
   residual subspace projection.

## Key physics constraints

- **T2 is sacred**: never collapse to scalar.  The angular residual IS the thesis.
- **Per-protein kernel normalization**: isolates angular structure, strips
  magnitude.  Kernel scale factors passed as scalars so the MLP knows
  what was stripped (bridges per-protein and global).
- **Kernel self-gating**: each kernel gates its own weight by its T2
  magnitude.  Distance dependence emerges from the kernels' own physics
  (1/r³ falloff, cutoffs) not from a training hyperparameter.
- **Mechanical mutants**: backbone doesn't change.  No DSSP/dihedral features.
- **WT features only**: ring kernels ARE the delta (ALA has no rings).
- **Equivariance**: scalar × L=2 = L=2.  T1 used as magnitude only.

## What NOT to do

- Do not use symlinks for data.
- Do not add graph structure or message passing — per-atom independence.
- Do not add DSSP/backbone features — mechanical mutants.
- Do not use global kernel normalization — per-protein is correct.
- Do not hardcode column indices — use SDK named accessors.
- Do not treat this as large-scale ML.  46 kernels, ~450K atoms.
- Do not vibe-code with hyperparameters — diagnose before sweeping.

## Files

### Active pipeline (use these)
- `src/calibration.toml` — single config: paths, features, training, analysis
- `src/mutation_set/config.py` — TOML → frozen Config dataclass
- `src/mutation_set/kernels.py` — KernelLayout registry + assembly
- `src/mutation_set/scalars.py` — ScalarLayout + named block assembly
- `src/mutation_set/dataset.py` — CalibrationDataset orchestrator
- `src/mutation_set/model.py` — ShieldingT2Model (gated mixing + correction)
- `src/mutation_set/train.py` — training entry point
- `src/mutation_set/evaluate.py` — R², baselines, per-protein metrics
- `src/mutation_set/analyze.py` — cold 7-section diagnostic
- `extract.py` — batch extraction driver

### Data
- `calibration/{ID}/` — real prmtop, xyz, nmr.out (723 proteins)
- `calibration/features/GatedCalibration/` — current extraction with ring_contributions
- `runs/` — training run outputs (header.json, epochs.jsonl, summary.json)

### Reference
- `JOURNAL.md` — learning progress record
- `EXPERIMENTS.md` — experiment log with results table
- `ARGUMENT.md` — thesis argument for validation approach
- `CALIBRATION_CHECKLIST.md` — parameter validation plan

### Historical (do not use for new code)
- `bones/pre_sdk/` — old c_equivariant pipeline, protein.py, features.py
- `bones/old_extractions/` — CalibrationExtractionTest, FirstExtraction
- `bones/` — session notes and exploratory work
