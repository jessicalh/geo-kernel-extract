# Calibration Experiments Log

All experiments from the 2026-04-08 session.  Runs are in learn/runs/.
Commit after every experiment from now on.

| Run name | Proteins | Epochs | Key change | Val R² | Val 0-4Å | Notes |
|---|---|---|---|---|---|---|
| smoke_test_v2 | 8 | 20 | First CalibrationDataset (old 66 scalars, 48 kernels, global kernel norm) | ~0 | — | Only 1 protein had delta, rest from old extraction |
| early_signal_13 | 13 | 50 | Same as above, more proteins | 0.16 | — | First real signal |
| critique_fixes | 29 | 30 | Global kernel norm, H/S element, z-score scalars | -0.06 | -0.10 | Global norm killed it — model learning nuisance scale |
| perprot_restored | 50 | 50 | Reverted to per-protein kernel norm, kept z-score + H/S element | **0.44** | **0.56** | Best result. Per-protein norm isolates angular structure |
| diagnostics_test | 45 | 30 | Added distance-weighted loss (tau=8Å), naive baselines | 0.05 | 0.09 | Lower than perprot — fewer val proteins? distance weighting? |
| xkernel_test | 57 | 50 | Added 18 cross-kernel T2 dot products (BS·HM, BS·PQ, BS·MC/Coulomb) | -0.04 | -0.12 | Overfitting. 86 scalars too many for 57 proteins. Stripped. |
| analysis_driven | 77 | 50 | Mutation identity (4), z/rho cylindrical coords, PQ dropped (40 kernels, 78 scalars) | **0.42** | **0.51** | On track. 4-8Å improved to 0.26. 8-12Å improved from -1.0 to -0.28 |
| 94prot_500ep | 94 | 500 | Same model as analysis_driven, 500 epochs, distance-weighted | **0.58** | **0.70** | Near 0.60 target at 13% of data. 8-12Å=0.29 (positive!). Median per-protein=0.60. Train-val gap=0.04 |
| efg_kernels_104s | ~105 | 500 | +6 EFG T2 kernels (46 total), +24 per-type T0 scalars, +2 MOPAC s/p pop (104 scalars) | **0.625** | **0.72** | Past 0.60 barrier. EFG angular structure is real signal. 8-12Å train=0.52. Per-protein median=0.64 |

## Analytical diagnostics (no learning)

Run: `analyze.py --run CalibrationExtractionTest` on 66 proteins.

- Per-protein ridge R²: median=0.81, range [0.60, 0.95] — kernels have the signal
- Global ridge R²: 0.23 — weights don't transfer across proteins
- BS and HM dominate (r≈-0.23), both negatively correlated with target
- H atoms strongest: HM_TYR r=-0.51
- PQ near zero everywhere (r<0.03) — dropped
- 0-4Å ridge R²=0.30 (hard), 8-12Å ridge R²=0.87 (easy, small signal)

## Key lessons

1. Per-protein kernel normalization is correct for calibration (angular structure question, not magnitude)
2. Global kernel normalization forces model to learn nuisance scale variable — hurts with limited data
3. Cross-kernel dot products overfit at <100 proteins — revisit at 700
4. Mutation identity (which aromatic residue removed) is cheap and useful
5. Distance-weighted loss helps mid/far-field, doesn't hurt near-field
6. Naive baselines (unweighted sums) are negative because of sign structure, not broken physics
7. The 0.60 target is reachable: per-protein ceiling is 0.81, MLP at 0.44 with 77 proteins

## Open issues (2026-04-08 end of session)

1. ~~make_model() defaults stale~~ — FIXED, now imports from dataset.py
2. WT ring neighbourhood geometry collapsed to per_type_T2 sums — individual
   ring contributions with distances not passed. The "stacked vs isolated"
   problem. Needs C++ to write per-ring T2+distance, or creative use of
   existing delta_ring_proximity for WT rings (not just removed ones).
3. No equivariance test — the correction head irreps declarations are untested.
   e3nn guarantees equivariance IF irreps are correct, but a wrong string is
   a silent numerical bug.
4. No shuffled-target sanity test — would confirm the model learns physics,
   not memorization artifacts.
5. Correction head self-TP (tp(x,x)) can only learn quadratic combinations —
   cannot express a linear additive correction. Need to test whether
   correction_scale learns anything meaningful or stays near 0.1.

## 2026-04-09: Azimuthal angle + new features

Added cos(φ), sin(φ) per (atom, ring) pair to C++ extraction (59-column
ring_contributions), MOPAC valency (sum of Wiberg bond orders), and MOPAC
molecular dipole vector.  Also built φ-modulated L=2 kernels (cos φ × BS T2,
etc.) and dipole-aligned kernels ((dipole · n̂) × BS T2).

Secondary analysis (tools 1-5, 723 proteins on old extraction) completed
first.  Key findings: HIE self-fit R² = 0.062 (worst ring type), kernel
space effective dimensionality = 3, ring_proximity is the most impactful
scalar group (+0.095 interaction delta).

Then tested new features on first 50 proteins of AzimuthalExtraction:

| Run | Kernels | Scalars | Ridge val | MLP val | Δ vs baseline |
|---|---|---|---|---|---|
| baseline_91k_50 | 91 | 247 | 0.557 | **0.602** | — |
| azimuthal_scalars_only_50 | 91 | 261 | 0.557 | 0.597 | -0.005 |
| azimuthal_cossinbs_50 | 103 | 273 | 0.575 | 0.582 | -0.020 |
| azimuthal_first50 | 127 | 297 | 0.617 | 0.586 | -0.016 |
| azimuthal_dipfix_50 | 127 | 297 | 0.617 | 0.586 | -0.016 |

**Findings:**
- Ridge goes UP with more kernels (+0.023 to +0.060) — new kernels have linear signal
- MLP goes DOWN with more kernels (-0.016 to -0.020) — model struggles with larger weight space
- New scalars alone are neutral (-0.005) — φ/valency/dipole as scalar context doesn't help yet
- Dipole kernel normalization (unit direction vs raw Debye) makes no difference (both 0.586)
- DipoleBs kernels were 1000x larger than other kernels before normalization fix

**Interpretation:**
The φ-modulated kernels pre-compute a product (cos φ × T2) that the MLP
could learn internally from scalar cos_phi × existing kernel weights.
Adding explicit kernels expands the weight space from 91 to 127 without
adding truly independent angular directions — the eigenspectrum stays at
~3 effective dimensions.  The MLP pays for 36 more weights to learn but
gets no new angular information it couldn't already construct.

**What to try next:**
- Full 723 run with scalars-only (91 kernels + 261 scalars) — fair test at scale
- Forward selection on 127-kernel layout to identify if any new kernel is selected
- Ring-ring relative normals as scalar (already computed in C++, not written)
- Reconsider: is the MLP architecture (hidden_mix=64) too small for 127 kernels?

## R analytical workbench (2026-04-09, 99 proteins, 127 kernels)

First run of the R workbench on 99 AzimuthalExtraction proteins.
All results are ridge regression — no training, no GPU.

**Full ridge R²: 0.380**

### Kernel groups (single):
| Group | Kernels | R² |
|---|---|---|
| efg | 6 | 0.300 |
| per_ring | 72 | 0.208 |
| ring_type | 32 | 0.186 |
| total | 7 | 0.036 |
| bond_cat | 10 | 0.026 |

EFG (6 kernels) beats ring_type (32 kernels) and per_ring (72 kernels).
MopacEFG_aro alone is R²=0.285 — 75% of total ridge signal in one kernel.

### Forward selection — no new kernel in top 24:
Top 10: MopacEFG_aro, EFG_aro, RingSusc_total, HM_HIE, MC_aromatic_total,
DeltaAPBS_EFG, BS_TYR, HM_H_ring0, PQ_total, BS_TRP_perimeter.
First φ-modulated kernel: CosChi_ring0 at position 25 (+0.001).

**Conclusion: the 36 new φ-modulated and dipole-aligned L=2 kernels
add no linear signal.  The φ information is better exploited as scalar
context (ring_proximity +0.305 interaction delta).**

### Scalar interaction analysis (base = top-10, R²=0.334):
| Group | Width | Δ R² |
|---|---|---|
| ring_proximity | 61 | +0.305 |
| residue_type | 20 | +0.146 |
| element | 5 | +0.098 |
| mcconnell | 6 | +0.097 |
| mopac_mcconnell | 6 | +0.097 |
| bond_orders | 4 | +0.096 |
| delta | 4 | +0.094 |
| mopac_electronic | 3 | +0.068 |
| mopac_dipole | 1 | +0.009 |

ring_proximity is 2x the next group.  It contains per-ring (z, rho,
theta, mcconnell_factor, exp_decay, disp_scalar, disp_contacts,
cos_phi, sin_phi) × 6 rings + n_pairs count.

### What this means:
- The L=2 kernel space is dominated by EFG (electrostatic angular pattern)
- Ring current kernels refine but don't dominate
- The MLP's value (0.35→0.61) comes from using scalar context to weight
  kernels per-atom, not from having more L=2 directions
- New explicit L=2 kernels hurt the MLP (0.602→0.586) because they expand
  the weight space without adding independent angular information
- The φ angle matters as a scalar (ring_proximity +0.305) not as a kernel

## Hidden probe (2026-04-09, azimuthal_scalars_only_50 model)

Loaded the trained MLP (hidden_mix=64, 91 kernels, 261 scalars) and
ran 31K atoms through it, capturing hidden activations at each layer.

**Hidden layer dimensionality:**
| Layer | Total dims | Active (>1% max eig) | 90% var | 95% var |
|---|---|---|---|---|
| Hidden 1 | 64 | 3 | 1 | 2 |
| Hidden 2 | 64 | 2 | 1 | 1 |
| Pre-gate output | 91 | 1 | 1 | 1 |
| Gated weights | 91 | 9 | 5 | 9 |

**Key finding:** The MLP collapses 261 scalars into 1-2 directions.
It learns a near-constant relative weight vector.  The kernel self-
gating (magnitude / (magnitude + threshold)) creates the spatial
diversity — 9 effective dimensions in the gated weights, all from
the kernels' own magnitude variation.

**Implications:**
- hidden_mix=64 is ~30x oversized.  The model uses 2-3 units.
- The 0.35→0.61 gap is mostly gating, not learned scalar-dependence.
- Adding kernels hurts because more weights in one near-constant
  vector means more things to get slightly wrong.
- The scalar interaction analysis (ring_proximity +0.305) simulates
  what gating does — it's not what the MLP actually exploits.
- This explains why new features (φ, dipole, RBF) don't help: the
  MLP isn't using scalar context to differentiate kernel weights.

**Open question:** Is the MLP FAILING to learn scalar-dependent
weights (training issue — LR, capacity, epochs), or is a near-
constant weighting genuinely optimal for this problem?  If the
per-protein ridge ceiling (0.81) requires per-protein weights, and
the MLP learns one global weight vector plus gating, then the
0.61 ceiling IS the gating ceiling.

## The 0.58 ceiling question

Historical pattern: ~0.58-0.60 R² regardless of protein count (94 or 720).
More data doesn't help. The gap to per-protein ridge (0.81) is real per-protein
structure that varies across proteins in ways the scalar features can't predict.
NOT noise, but too specific to generalise.

Next: audit all 51 NPY arrays against what reaches the model. The answer may
be in features we're extracting but not using.

## 2026-04-10: Per-element physics and the end of the MLP

### The reframe

The MLP was an engineering detour.  The goal is not R² optimisation
but physics explanation: do the geometric kernels see what NMR theory
says they should see?  The answer is yes — and the MLP obscures it.

### Per-element kernel decomposition (element_physics.py)

Stratified 110 proteins (69,080 atoms) by element × kernel physics
group.  Literature prediction confirmed in the T2 channel:

| Element | n_atoms | Ring current R² | EFG R² | All R² |
|---------|---------|----------------|--------|--------|
| H       | 35,080  | **0.909**      | 0.770  | 0.935  |
| C       | 20,673  | 0.156          | **0.431** | 0.546 |
| N       | 6,100   | **0.325**      | 0.100  | 0.434  |
| O       | 6,927   | **0.317**      | 0.196  | 0.387  |
| pooled  | 69,080  | 0.213          | 0.318  | 0.385  |

Ring current dominates H (Boyd & Skrynnikov 2002).  EFG dominates C
(Sahakyan & Vendruscolo 2013).  N is multi-mechanism (Xu & Case 2002).
The pooled view hides all of this.

### 28 NMR realities verified

20 analytical + 8 tensor-direct tests against known NMR physics:
- BS distance slope = -3.04 (theory -3)
- PQ distance slope = -5.05 (theory -5)
- HM/BS magnitude ratio = 1.009 ± 0.037
- BS-HM far-field cos = 0.9999
- H-bond contribution = 0.002 R² (mechanical mutant: backbone unchanged)
- Element weighting gap = +0.333 R²

Full writeup: learn/docs/twenty_eight_realities_2026-04-10.md

### The MLP is worse than per-element ridge

| Approach | R² | Parameters |
|----------|-----|-----------|
| Pooled ridge (55 kernels) | 0.374 | 275 |
| Gated MLP (h=64, 261 scalars) | 0.610 | 26,624 |
| Per-element ridge (55 kernels) | 0.690 | 1,100 |
| + per-protein norm + scales + mut_type | **0.818** | ~2,200 |

### Progressive calibration (physics_calibration.py)

Tested each scalar interaction independently:
- Kernel scale factors: +0.06 per element (what normalization stripped)
- Mutation type (which residue removed): +0.03 for H, +0.22 for N
- MOPAC valence: +0.000
- Bond order: +0.000
- Molecular dipole: +0.000

The kernels already encode everything valence/bond/dipole know.
The fair physics model is: 55 kernels + per-protein norm + element +
mutation type.  R² = 0.818.

### The calibration answer

The ridge coefficient vector (55 × 4 elements × ~5 mutation types)
IS the table of calibrated physical constants — the tensor-valued
extension of Case's (1995) Table 3.  No MLP needed.

### Open for next session

- Rerun on full 723 proteins (AzimuthalExtraction: 613 remaining)
- Steal OLS-vs-ridge subspace test from analyze.py per element
- Steal leave-one-out ablation from analyze.py per element
- Pin down the thesis coefficient table with confidence intervals

### Scripts

All analysis in learn/src/actual_physics/ (see OVERVIEW.md there).
Figures in learn/R/twenty_eight_realities.R.
LaTeX in learn/docs/realities_latex.tex.
