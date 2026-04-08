# Agent Review: Round 2 Model Analysis (2026-04-06)

**Status:** User finds recommendations sound but wants discussion before
implementation. Do not implement without consulting user.

Two adversarial agents reviewed the model independently. The first
focused on feature improvements, the second on structural soundness.
Combined findings below.

---

## Structural Soundness (Agent 2)

### Verdict: model is architecturally sound

- **Equivariance** is correct. scalar(MLP) x L=2(kernel) = L=2. The
  self-gating (using |T2| magnitudes as scalar inputs) is rotationally
  invariant and standard, not circular.
- **Loss function** (equal MSE on 5 T2 components) is the unique
  rotationally-invariant choice. Equivalent to minimising tensor
  Frobenius norm error.
- **Per-protein kernel normalization on deltas** is mathematically
  sound. Far-atom noise is NOT amplified because std is computed over
  all M*5 elements. Dead kernels (11/46) are handled by the > 1e-10
  check.
- **Output scaling** is correct after the sqrt(n_kernels) fix.
  Dynamic range is adequate.

### One minor bug found

Val R^2 computation unscales val targets using train_ds.target_std,
but val targets were scaled by val_ds.target_std. Effect is sub-1%
(both stds are ~2.0 ppm from same population) but should be fixed
for correctness.

Fix: compute R^2 on raw unscaled targets, or pass the correct
target_std per split.

### Python atom matching

Greedy nearest-neighbor with 0.5A tolerance. The C++ MutationDeltaResult
uses KD-tree with bijection enforcement and second-chance passes.
Estimated impact: ~1% R^2. Not the gap source, but should move to C++
for correctness (MutationDeltaResult.WriteFeatures — next session).

### The residual

The residual is what we haven't explained yet. It contains a mix of
paramagnetic contributions, charge redistribution, geometry effects,
and possibly missing angular structure. The relative shares are
unknown and determining them is part of the ongoing work. Do not
treat any breakdown as a ceiling — it is an open question how much
of the residual is recoverable with better modeling.

### Useful diagnostics (when ready)

- Per-component R^2: which of the 5 T2 m-values is worst?
- Per-distance-zone R^2: near-ring (< 5A) vs far.
- These help understand WHERE to improve, not WHETHER to improve.

---

## Feature Improvements (Agent 1)

### Priority 1: Add residue_type as 20-dim one-hot

**Rationale:** The MLP knows element (C/N/O) and ring distance but
not what amino acid the atom belongs to. A backbone carbonyl in GLY
behaves differently from one in LEU when a nearby PHE is removed.
The MLP needs local chemical context independent of the kernels
to learn which kernel combination is appropriate.

The current 6 magnitude scalars are derived from the kernels themselves
(circular — the MLP needs context independent of what it's weighting).
Residue type encodes backbone phi/psi propensities, sidechain bulk,
polarity, and hydrogen bonding patterns.

**Data:** Already loaded (load.py:123, residue_type.npy), already
extracted. Change: ~10 lines in train.py, one arg in model.py.
13 -> 33 scalar features, ~1,280 extra parameters. Train/val gap
(0.06) gives headroom.

**Expected impact:** Val R^2 0.55 -> 0.60-0.65 (agent estimate).

### Priority 2: Add MOPAC delta charge as scalar

**Rationale:** The delta charge (WT mopac_charge - ALA mopac_charge)
directly measures electronic reorganisation from the mutation. Tells
the MLP "how much did the electron distribution change at this atom."
If the delta charge is large but the Coulomb kernel is small, the
MLP can upweight the correction.

**Data:** mopac_charges.npy already extracted. Two features: absolute
WT charge and delta charge. 33 -> 35 scalar features.

### Priority 3: Add T0 kernel magnitudes as scalars

**Rationale:** Current magnitude scalars use only T2 norms. T0
(isotropic part) carries independent information about total mechanism
strength. A large T0 with small T2 = strong isotropic effect with
weak angular structure, fundamentally different physics.

6 more scalars: |delta_T0| for MC, Coulomb, BS, HBond, MopacCoulomb,
MopacMC. All L=0 (scalar), equivariance preserved.

### What NOT to do

- **Don't increase MLP hidden size.** Not capacity-limited (7,936
  params on 48K atoms, 0.06 gap). Bottleneck is input information.
- **Don't remove per-protein normalization.** Confirmed by both
  agents as correct and important.
- **Don't change the training schedule.** lr=1e-3, wd=1e-4, bs=512,
  CosineAnnealingLR are well-tuned and match the working backup.
- **Don't add L=1 features.** E-field vectors (coulomb_E, bs_total_B)
  require architectural changes to maintain equivariance. Save for
  a future round.
- **Don't break the equivariant constraint** (one scalar weight per
  kernel across all 5 T2 components). Per-component weights (230
  weights) would overfit to frame-dependent correlations. The
  constraint is scientifically essential.

---

## User notes (to be added after discussion)

- User wants to understand the C++ / Python boundary implications
  of each recommendation before implementing
- Residue type comes from the C++ object model (Protein::ResidueAt)
  and is written by WriteFeatures — this is not Python reinventing
  atomic modeling, it is using what the C++ already provides
- MOPAC charges likewise come from MopacResult in C++
- Discussion needed: does adding residue type risk the model learning
  amino-acid-specific biases that don't generalise, rather than
  physics? (Counter: the mutation dataset has diverse residue
  environments, and the val split is by protein)
