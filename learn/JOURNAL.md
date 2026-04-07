# Learning Journal

Living record of training rounds. Every round is an experimental
result with logged outputs — not a tuning exercise.

## What this system does

Eight classical calculators (+ 2 MOPAC-derived) compute geometric
kernels — the spatial shape of each shielding mechanism as rank-2
tensors. We validate these against DFT on ~720 wild-type / alanine
mutant pairs. The WT-ALA delta isolates the ring contribution:
everything not ring-related cancels. The R^2 measures how much of
the DFT mutation delta our geometric kernels explain.

R^2 = 0.7 means: 70% of what DFT says the aromatic ring does to
the angular shielding pattern is captured by our classical physics
calculators with learned scalar parameters. The residual is what
the classical picture misses — itself a thesis result.

## What flows from this

1. **Calibration**: the validated tensors let us sweep geometry
   parameters (cutoffs, filter extents, shaped boundaries) against
   the DFT delta and find empirically optimal values for every
   constant in the pipeline.
2. **MOPAC integration**: QM charges and bond orders inform geometry
   decisions — not as model features but as continuous descriptions
   of electronic structure for cutoff evaluation and inclusion.
3. **Upstream model**: the calibrated, validated tensors pass up to
   an e3nn NMR prediction model as physically grounded features.
   The upstream model doesn't rediscover ring currents — it gets
   tensors that already explain the angular physics, with the
   residual telling it where to look for what's missing.

## Round structure

Each round produces:
- `learn/runs/{run_name}/header.json` — hyperparameters, data
- `learn/runs/{run_name}/epochs.jsonl` — per-epoch train/val loss
- `learn/runs/{run_name}/summary.json` — final R^2, RMSE
- `learn/runs/{run_name}/best_model.pt` — best val checkpoint

**Target**: WT-ALA delta DFT T2 (the ring effect, ~2 ppm scale).
**Features**: WT-ALA delta kernel T2 (what each calculator says
changes when the ring is removed).
**Model**: environment-dependent mixing weights via MLP on scalar
features × L=2 kernel tensors. Equivariant by construction.

## Rounds

### Round 1: 8 classical calculators (K=40)

**Status:** Results lost. Run predates current git history.
**What we know:** Mixing-only model (no correction head) with
warmup and cosine schedule gave the best results. Correction head
oscillated and didn't help. Exact R^2 lost but was the benchmark.
**Data:** FirstExtraction_pre_mopac/ (36 arrays per conformation).
**CRITICAL FIX discovered in Round 2 setup:** the training code was
initially regressing on raw WT DFT T2 (~60 ppm scale) instead of
WT-ALA delta T2 (~2 ppm scale). This was corrected 2026-04-06.
Whether Round 1's good results used deltas or raw WT is uncertain —
the compiled bytecode doesn't preserve enough to tell. If Round 1
was on raw WT, Round 2 is the first real delta result.

### Round 2: 10 calculators with MOPAC (K=46) — CURRENT

**Status:** Extraction running (2026-04-06), ~480/1450 jobs done.
Ridge (Level A): **R^2 = 0.37** on delta T2 (~240 proteins).
Equivariant mixing: **Val R^2 = 0.60** (33 scalars with residue
type, 46 kernels, ~240 proteins, 200 epochs, best at epoch 171).
MopacCoulomb is the #1 Ridge contributor.

**What changed from Round 1:**
- MopacCoulombResult: Coulomb EFG from PM7 Mulliken charges (1 feature)
- MopacMcConnellResult: bond-order-weighted McConnell (5 features)
- 14 scalar features (was 12): |MopacCoulomb T2|, |MopacMC T2|
- MopacCoulomb T2 as 5th L=2 input to correction head (if used)
- All extraction from single C++ binary (MOPAC inside pipeline)
- Clean 46-array extraction, no mixed old data
- Training on WT-ALA delta T2 (confirmed correct)
**Data:** FirstExtraction/ (46 arrays per conformation)

### Architecture note: C++ owns deltas, Python owns learning

MutationDeltaResult in C++ already computes proper topology-aware
atom matching, delta DFT shielding, delta charges (ff14SB + MOPAC),
ring proximity with cylindrical coordinates. But it has no WriteFeatures.
Currently Python reimplements atom matching (position heuristic, 0.5A).

**Next step**: add WriteFeatures to MutationDeltaResult so the C++ writes
authoritative delta .npy files (delta_dft_T2, matched_indices,
ring_proximity, delta per-calculator). Python just loads and trains.
Clean boundary: C++ owns the physics object model, Python owns learning.
This matters for any academic use of this system.

### Round 3 (planned): GeometryChoice documentation

GeometryChoice objects on each calculator (see spec/GEOMETRY_CHOICE_BRIEF.md).
Round 2 R^2 is the "did we break anything" baseline.

---

## Key findings (append as discovered)

- **Target must be WT-ALA delta T2**, not raw WT T2. Raw WT has
  ~60 ppm variance dominated by local electronic structure that
  kernels cannot explain. Delta has ~2 ppm variance from ring
  effects that kernels are designed to capture.
- Correction head was not helpful in Round 1. Start mixing-only.
- 20-epoch linear warmup (lr/20 → lr) stabilises early training.
  Without it, full LR on random weights causes oscillation.
- CosineAnnealingWarmRestarts(T_0=50, T_mult=2) after warmup.
- lr=1e-4, weight_decay=1e-5, batch_size=2048, 500 epochs.
- Ridge (Level A) on deltas: R^2 = 0.35 at 81 proteins (Round 2).
  MopacCoulomb is the single largest contributor.
- Ridge on raw WT T2 gave R^2 ~ 0.20 — misleadingly low because
  the kernels explain the ring part, not the atom's intrinsic part.
- Ridge delta T2 R^2 = 0.356 at 105 proteins. MopacCoulomb is the
  #1 single-feature contributor. This is the linear floor.
- The pattern (good T2, bad T0, independent calculators) appeared
  early in previous rounds (~0.6 with the equivariant model, climbing
  to ~0.7 with more data and tuning).
- **BUG FOUND + FIXED (2026-04-06):** KernelMixingHead output scale
  mismatch. MLP outputs 46 weights with std~1, each multiplied by
  kernel with std~1, summed → output std~sqrt(46)~7. Target std=1.
  Model spent 500 epochs shrinking weights instead of learning.
  Fix: `weights / sqrt(n_kernels)` in forward(). One line.
  Result: R² went from -1.27 to **+0.48** (Val, 160 proteins).
  Ridge baseline: 0.37. Still improving at epoch 500.
- Adversarial agent (Opus) identified the scale issue. Also flagged
  dead scalar features and unbalanced scalar scales.
- **CLEAN BACKUP FOUND** at bones/nmr-shielding-clean-20260405/.
  The working code used: (a) per-protein kernel normalization
  (each protein's kernels divided by own std before stacking),
  (b) per-protein scalar normalization, (c) lr=1e-3, wd=1e-4,
  bs=512, CosineAnnealingLR, 50 epochs, no warmup. And it trained
  on **raw WT T2, not deltas**.
- The R²~0.7 from Round 1 was on raw WT T2 (all shielding anisotropy).
  Our delta formulation is a different, harder question (ring-specific
  effect only, ~2 ppm scale vs ~60 ppm). R²=0.55 on deltas with
  per-protein norm (170 proteins, 200 epochs) is the correct comparison.
- Per-protein norm is the key architectural choice. It makes the MLP
  learn relative kernel patterns within a protein rather than absolute
  magnitudes across proteins. This is physically right — the shape of
  the kernel combination matters, not the scale.
- **Round 1 target: R²=0.71 on T2 deltas** (per user). The clean
  backup trains on raw WT T2 and predates the delta switch. The
  0.71 came after switching to deltas in-session — that code did
  not survive. Our 0.60 on ~240 proteins is partial data; the full
  725-protein run will show whether we reach 0.71.
- **Residue type one-hot (20 AA)** moved Val R² from 0.55 to 0.60.
  Train/val gap narrowed from 0.06 to 0.03. The MLP needs local
  chemical context independent of the kernels it's weighting.
- **Val R² bug fixed**: val targets were unscaled with train_std
  instead of val_std. Sub-1% impact but now correct.
- Agent review (AGENT_REVIEW_20260406.md) confirms architecture is
  sound, equivariance is correct, per-protein normalization is right.
  The 40% residual is dominated by paramagnetic contributions
  (~15-20%) that no classical kernel models.
