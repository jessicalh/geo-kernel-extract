# Ensemble Accumulation: Session 1 (2026-04-11)

## Context

The calibration chapter is settled (R²=0.818, 28 realities, per-element
ridge on 55 kernels).  This session begins the ensemble/prediction
project: what to accumulate across 685 proteins × 5005 metadynamics
frames, to feed the GNN in nmr-training/ and validate the calibration
against experimental shifts (RefDB).

All analysis scripts are in `learn/src/actual_physics/`.  All use the
same infrastructure (iter_proteins, KernelLayout, assemble_kernels).

## Finding 1: The geometry-only kernel space

**Script:** `geometry_only_basis.py`

44 of 55 core kernels survive without MOPAC/APBS (the 11 dropped are
MopacMC categories, MopacCoulomb/MopacMC totals, MopacEFG, APBS,
DeltaAPBS).

Fair-set R² (per-protein normalization + kernel scales + mutation type):

| Element | Full (55) | Geo-only (44) | Gap    |
|---------|-----------|---------------|--------|
| H       | 0.949     | 0.942         | 0.007  |
| C       | 0.691     | 0.506         | 0.185  |
| N       | 0.704     | 0.662         | 0.042  |
| O       | 0.632     | 0.587         | 0.045  |
| Weighted| 0.818     | 0.750         | 0.067  |

**Carbon is the problem.**  The MOPAC advantage for C is signal quality,
not a missing dimension.  ff14SB-MOPAC EFG cosine = 0.611 for C (partially
correlated).  MOPAC EFG explains only 5.6% of the geo-only residual.

**Nitrogen doesn't need MOPAC.**  ff14SB and MOPAC EFGs are nearly
identical (cos=0.896, R² difference <0.001).

**Oxygen prefers ff14SB.**  ff14SB EFG R²=0.172 vs MOPAC R²=0.156.

Physics-group decomposition (geo-only, normalized):

| Element | Ring current | ff14SB EFG | Bond aniso | Quadrupole | Dispersion |
|---------|-------------|------------|------------|------------|------------|
| H       | 0.800       | 0.668      | 0.122      | 0.047      | 0.393      |
| C       | 0.132       | 0.228      | 0.035      | 0.020      | 0.122      |
| N       | 0.246       | 0.077      | 0.046      | 0.147      | 0.170      |
| O       | 0.231       | 0.165      | 0.048      | 0.078      | 0.272      |

Forward selection in geo-only: EFG_aro enters first for ALL elements.
H: EFG_aro(0.668) → RingSusc(0.732) → ...BS types.
C: EFG_aro(0.228) → Disp types.  N: EFG_aro(0.072) → PQ_total → Disp.
O: EFG_aro(0.163) → Disp types.

**Open question:** Can cheaper charge estimates (Gasteiger, EEM, MMFF94)
recover any of the carbon gap without the MOPAC cost?

## Finding 2: Variance and prediction are orthogonal

**Scripts:** `dimensionality_test.py`, `accumulation_basis.py`

Raw PCA eigenspectrum: 2 blocks of 5 eigenvalues (corresponding to
5 T2 components of 2 physical sources).  PC1-PC10 capture 99.3%+ of
variance.  But the ridge weight is **0.00% in PC1-PC10**.  90% of
ridge weight is in PC77-96.

The high-variance PCs are Coulomb_total (median 4.5-11.2) and EFG_bb
(median 1.1-10.4).  The prediction comes from BS/HM/PQ (median 0.007).
1000x magnitude gap.

PCA-then-ridge (predictive dimensionality):
- Raw kernels: 3 predictive PCA dims for all elements in all contexts
- Normalized: H≈46, C≈6, N≈3, O≈11 (but PCA misses most of the
  ridge signal — N gets val R²≈0 from PCA but R²=0.425 from ridge)

Ridge coefficient matrix SVD: effective rank ≈ 5 for all elements.
The ridge uses all 5 T2 components independently.

Predicted variance per kernel sums to 400% (H) through 1400% (N),
indicating massive cancellation.  McConnell dominates predicted
variance everywhere despite not being the most important coefficient.

Cross-kernel prediction correlation reveals two uncorrelated clusters:
McConnell (bond anisotropy) and ring current (BS, HM, RingSusc, EFG_aro).
Correlated within cluster (0.2-0.5), uncorrelated between (0.01-0.04).

**Consequence:** Accumulate per-kernel raw T2 moments independently.
Do NOT pre-project onto PCA or ridge basis.  PCA destroys signal.
Ridge projection amplifies cancellation noise.

## Finding 3: Ensemble dynamics regime

**Script:** `expected_ensemble_variance.py`

Using measured power-law exponents and dimensional analysis:

    δK/K = n × δr / r

At moderate dynamics (δr = 1.5 Å RMSF):

| Kernel | n     | Median δK/K | Scrambled | Critical dist |
|--------|-------|-------------|-----------|---------------|
| BS/HM  | 3.04  | 0.45        | 5%        | 4.6 Å         |
| Chi    | 3.00  | 0.44        | 5%        | 4.5 Å         |
| PQ     | 5.05  | 0.74        | 25%       | 7.6 Å         |
| Disp   | ~6    | (worse)     | >50%      | ~9 Å          |

**Nothing is stable.**  0% of atom-ring pairs have δK/K < 0.3.
Every pair is in the dynamic regime (variance 30-100% of mean)
or worse.

The 1/r³ kernels are the workhorse of the ensemble: reliable
Welford mean with physically meaningful variance.  PQ is half-
scrambled; only the nearest pairs carry signal.  Dispersion is
mostly noise at ensemble timescales.

The variance IS information: it tells the GNN which atoms are
in the dynamic regime and which kernel features are reliable.

Element-independent: the geometric kernels see the same distance/angle
relationships for all probe atoms.  Element dependence enters only
through the ridge coefficients.

**Intermediate assumptions that need later work:**
1. RMSF values assumed from MD literature, not measured from fleet XTCs
2. δK/K analysis uses distance only, ignores angular sensitivity
3. "Scrambled" means variance > mean, but the Boltzmann-weighted mean
   is still the correct expectation — variance tells you the spread
4. McConnell cancellation: physics or ridge artifact?

## Accumulation design (current understanding)

1. **Per-atom, per-kernel, full 9-component tensor Welford moments.**
   Raw magnitudes, no normalization, no projection.  Each kernel at
   its own scale.  The 1000x magnitude hierarchy IS information.

2. **Per-ring decomposed contributions** as well as per-atom totals.
   Ring identity is stable across frames (same topology).
   Memory: ~67 MB for largest protein, well within budget.

3. **Ring geometry moments** (center, normal, radius per ring).
   Root cause of kernel variation.

4. **Frame diversity filter** keeping 10-20 survivors.

5. **DSSP ensemble stats** and **Coulomb field moments** as context.

6. **Do NOT pre-project, pre-whiten, or gate during accumulation.**
   Downstream consumers handle scale via their own normalization.

## Finding 4: AIMNet2 recovers the carbon dimension

**RESOLVED:** The carbon gap was a lost EFG dimension (ff14SB EFG
points 45° wrong for carbon, cos=0.709 with MOPAC).  Tested AIMNet2
neural network potential (ωB97M-trained Hirshfeld charges) on 10
calibration proteins.

EFG T2 direction cosine with MOPAC:

| Element | ff14SB | AIMNet2 | Hybrid (H+C) |
|---------|--------|---------|---------------|
| H       | 0.609  | 0.932   | 0.909         |
| C       | 0.709  | 0.971   | 0.958         |
| N       | 0.974  | 0.990   | 0.986         |
| O       | 0.998  | 0.999   | 0.992         |

AIMNet2 timing on RTX 5090 (batcave):
- 2480 atoms: 0.04s warm
- 4876 atoms: 0.17s warm
- 5000 frames × 0.17s = 14 minutes per protein (every frame)

**MOPAC is eliminated from the ensemble pipeline.**

## Finding 5: Boltzmann N_eff ≈ 49

From COLVAR files of the fleet test protein (1Q8K_10023, 5 walkers,
5005 frames, bias factor 12, T=298K):

- N_eff = 49 (1% efficiency)
- 90% of weight in 18 frames
- Subsampling from 5005 to 500 barely matters (N_eff 49 → 19)

**But:** the NMR experiment is not a Boltzmann average.  Low-weight
frames may dominate through chemical exchange.  The GNN gets
Boltzmann weights as input features, not as a selection filter.
Accumulate across ALL frames; let the GNN learn weighting.

## The GNN training package

### Three protein populations

**A: Fleet + RefDB (up to 685).** Full 5005-frame ensemble extraction
+ AIMNet2 charges.  Experimental T0 target from RefDB.

**B: RefDB without fleet (~630).** Single PDB structure, one frame
of features.  Same T0 target.  No ensemble stats (variance=0).

**C: Calibration with DFT (110-723).** WT-ALA pairs, full T2 target
from ORCA.  Validates kernels but different target.

### Graph structures to export per protein

**1. Radius graph** — atom pairs within 8-10Å (from Boltzmann
minimum or PDB).  Edge features: distance, relative position
vector (L=1), bond type.  Standard for molecular GNNs.

**2. Bond graph** — covalent bonds from topology.  Currently
computed but not exported.  Trivial to add as `bonds.npy`.
Edge features: bond category, bond order.

**3. Atom-ring graph** — `ring_contributions.npy`, already exists.
Each (atom, ring) pair within 15Å with per-ring BS/HM/PQ/Chi
T2 tensors + cylindrical geometry.  **Uniquely ours.  No other
NMR predictor has per-ring tensor edge features.**

### Node features

Static (topology): element, residue type, atom role, hybridisation,
graph distance to ring/N/O, conjugation flags.

Per-frame accumulated (Welford, Boltzmann + unweighted):
- 8 calculator shielding tensors (full T0+T1+T2)
- Per-type T2 decomposition (32 ring-type kernels)
- McConnell category T2 (5 bond categories)
- ff14SB Coulomb EFG (bb, aro)
- AIMNet2 Hirshfeld charges
- AIMNet2 Coulomb EFG
- Backbone dihedrals (phi, psi)
- Chi angles (chi1-chi4)
- Per-atom SASA
- DSSP secondary structure

### Warm start from calibration (Option A)

Pre-multiply per-ring kernel T2 by calibrated per-element ridge
coefficients.  Edge features on the atom-ring graph become the
**calibrated contribution**, not the raw kernel.  The GNN starts
from the physics answer (R²=0.818) and learns corrections.

Calibration metadata shipped once (not per protein):
- Per-element ridge coefficients (4 elements × 55 kernels × 5×5)
- Per-kernel importance rankings
- Calibrated ring current intensities per ring type

### vs UCBShift/SPARTA+

They send per-residue scalars from one PDB.  We send per-atom
tensors from 5005 frames.  They predict T0 only from a ring
current scalar + backbone dihedrals + SASA.  We provide full
T0+T1+T2 from 8 physics calculators + geometry-dependent charges
+ ensemble dynamics + pre-calibrated warm start.

Features they have that we match or exceed: backbone dihedrals
(dssp_backbone.npy), SASA (adding per-atom), ring current
(full tensor not scalar), H-bonds (full tensor), E-field
(full vector + EFG tensor), secondary structure (DSSP).

## C++ changes required

**NOTE:** This section is superseded by
`learn/docs/cpp_marching_orders_2026-04-11.md` which has the
refined spec including design decisions from the follow-up session.
Read that document for implementation.  This section is kept for
historical context only.

### Summary of calculators identified

1. AIMNet2Result — charges, aim embedding, Coulomb EFG
2. Chi angles — sidechain dihedrals
3. Per-atom SASA — Shrake-Rupley
4. Bond graph export
5. Radius graph export

## Open questions (from session 1, some resolved in session 2)

1. Ridge cancellation structure: McConnell physics (why does it
   dominate predicted variance?)
2. Angular sensitivity (theta dependence, not just distance)
3. Actual RMSF from fleet XTC data (replace assumed scenarios)
4. ~~libtorch build integration~~ — resolved: CUDA mandatory, find_package(Torch)
5. ~~Ensemble observer architecture~~ — deferred: calculators first
6. What does the GNN see for Population B (single PDB) proteins
   that makes the ensemble features gracefully optional?
