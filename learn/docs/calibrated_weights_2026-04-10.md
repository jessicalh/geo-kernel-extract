# Per-Element Calibrated Weight Vectors

## What these are

The calibration pipeline fits per-kernel weights against DFT shielding
tensor deltas.  When trained separately per element (H, C, N, O), the
weight vectors reveal which physical mechanisms dominate the T2
shielding perturbation for each atom type.  These weights are the
tensor-valued extension of Case's (1995) ring current intensity
factors — calibrated physical constants, not model parameters.

110 proteins, 69,080 atoms.  Models: KernelMixingHead with hidden=4
(the minimum needed; pooled hidden=64 collapses to 2-3 effective
dimensions anyway).  Gating by kernel self-reported T2 magnitude.

## Performance

| | Ridge | Model (h=4) |
|---|---|---|
| H (35,080) | 0.905 | **0.906** |
| C (20,673) | 0.562 | **0.666** |
| N (6,100) | 0.596 | 0.222 |
| O (6,927) | 0.534 | 0.293 |
| Pooled | 0.374 | 0.501 |
| **Weighted per-element** | — | **0.711** |

Four tiny per-element models outperform the pooled model (0.711 vs
0.501).  H matches the ridge ceiling — gating adds nothing beyond
linear for hydrogen.  C exceeds ridge by 0.10 — gating helps carbon
by suppressing distant kernels.  N and O underperform ridge because
the training data is thin (6K atoms) relative to the 121-kernel
weight space.

## Dominant weights by element

The top 5 kernels by |weight| for each element, with physical
interpretation.  All weights are mean gated weights across atoms
(including the effect of gating by kernel magnitude).

### Hydrogen

| Kernel | Weight | Physical mechanism |
|--------|--------|-------------------|
| MopacEFG_aro | -0.036 | Buckingham EFG from aromatic ring (PM7 charges) |
| EFG_aro | -0.031 | Buckingham EFG from aromatic ring (ff14SB charges) |
| RbfBsMid_ring0 | -0.022 | Biot-Savart ring current at 4-8 A (nearest ring) |
| AngBsEquat_ring0 | -0.022 | Biot-Savart in-plane component (nearest ring) |
| RbfBsFar_ring0 | -0.017 | Biot-Savart ring current at 8-12 A |

**Interpretation:** H weight vector is dominated by EFG (far-field
angular pattern) plus distance-resolved BS ring current.  The
negative signs are correct for a deletion delta.  The EFG and BS
work together — EFG captures the overall angular direction, BS
refines the distance-dependent amplitude through gating.

### Carbon

| Kernel | Weight | Physical mechanism |
|--------|--------|-------------------|
| MC_CN_nearest | -0.056 | McConnell: nearest C-N peptide bond anisotropy |
| MopacEFG_aro | -0.035 | Buckingham EFG (PM7 charges) |
| MopacMC_CO_nearest | -0.034 | McConnell: nearest C=O (MOPAC bond order weighted) |
| MC_backbone_total | +0.029 | McConnell: all backbone bond anisotropy |
| MopacMC_total | +0.026 | McConnell: total bond anisotropy (MOPAC weighted) |

**Interpretation:** Carbon is bond-dominated.  The single strongest
weight (MC_CN_nearest = -0.056) is the peptide C-N bond, whose
magnetic anisotropy changes when the adjacent aromatic sidechain is
removed.  EFG is secondary.  This matches Sahakyan & Vendruscolo's
finding that heavy-atom shifts are electrically sensitive — but for
T2, the bond anisotropy mechanism (McConnell) is equally important
as the EFG.

### Nitrogen

| Kernel | Weight | Physical mechanism |
|--------|--------|-------------------|
| PQ_total | +0.029 | Pi-quadrupole EFG (all rings) |
| DeltaAPBS_EFG | -0.024 | Solvation-corrected EFG |
| MopacEFG_aro | -0.021 | Buckingham EFG (PM7 charges) |
| MopacMC_CO_nearest | +0.023 | Nearest C=O bond (MOPAC weighted) |
| MC_total | +0.017 | Total bond anisotropy |

**Interpretation:** No single mechanism dominates.  The largest weight
is PQ_total — the pi-electron quadrupole, an electrostatic effect
from the ring's charge distribution (1/r^5).  DeltaAPBS_EFG entering
at rank 2 is notable: the solvation correction matters for N but
not for H or C.  This suggests the amide nitrogen's lone pair is
sensitive to the dielectric environment around the aromatic ring.

### Oxygen

| Kernel | Weight | Physical mechanism |
|--------|--------|-------------------|
| MopacEFG_aro | -0.057 | Buckingham EFG — strongest of any element |
| EFG_aro | -0.031 | Classical EFG |
| DeltaAPBS_EFG | +0.027 | Solvation EFG |
| HM_H_ring0 | +0.025 | Haigh-Mallion surface integral (nearest ring) |
| RbfBsNear_ring0 | -0.023 | Near-field BS ring current |

**Interpretation:** Oxygen's strongest weight is MopacEFG_aro at
-0.057, the largest single weight in the entire table.  Carbonyl
oxygen responds more strongly to the aromatic EFG than any other
element.  The C=O bond itself doesn't change (McConnell absent),
so all signal comes from through-space fields.

## Cross-element patterns

### EFG weight varies by element

| Element | EFG_aro | MopacEFG_aro |
|---------|---------|-------------|
| H | -0.031 | -0.036 |
| C | -0.020 | -0.035 |
| N | **+0.001** | -0.021 |
| O | -0.031 | -0.057 |

Nitrogen's classical EFG_aro weight is essentially zero.  The MOPAC
EFG still contributes (-0.021) — PM7 charges capture electronic
structure that ff14SB point charges miss at the amide nitrogen.

### Bond anisotropy is element-specific

| Element | MC_CN_nearest | MC_backbone | MC_aromatic |
|---------|--------------|-------------|-------------|
| H | +0.001 | +0.001 | +0.005 |
| C | **-0.056** | +0.028 | -0.018 |
| N | +0.001 | +0.017 | -0.001 |
| O | -0.004 | -0.021 | -0.015 |

MC_CN_nearest is exclusively a carbon feature (weight 14x larger
than any other element).  This is physically correct: the C-N
peptide bond's magnetic anisotropy directly affects the shielding
tensor of the bonded carbon atoms.  Hydrogen doesn't feel bond
anisotropy (it's bonded to only one atom with short bond length).

### Ring current is universal in sign, element-specific in magnitude

| Element | RbfBsMid_ring0 | HM_H_ring0 | RingSusc |
|---------|---------------|------------|---------|
| H | -0.022 | +0.017 | -0.006 |
| C | -0.021 | +0.012 | +0.019 |
| N | -0.015 | -0.008 | +0.002 |
| O | -0.020 | +0.025 | -0.009 |

The mid-range BS kernel (RbfBsMid) has consistent negative weight
across elements — the ring current distance dependence is universal.
HM_H_ring0 (the raw surface integral, pure T2) reverses sign for
nitrogen, consistent with the paramagnetic-term dominance of 15N
shielding (Saito et al. 2010).

## Gating thresholds

Each kernel self-gates by its T2 magnitude.  The threshold is the
median magnitude in the training set — kernels below threshold
contribute proportionally less.  Low threshold = sensitive at low
magnitude = long-range.  High threshold = only active when strong
= short-range.

| Kernel family | Median threshold | Range character |
|--------------|-----------------|-----------------|
| BS per-type | 0.008-0.024 | Long-range (dipolar 1/r^3) |
| HM per-type | 0.005-0.024 | Long-range |
| EFG | 0.20-0.45 | Medium-range |
| McConnell | 0.04-0.10 | Medium-range |
| Dispersion | 0.0001-0.001 | Short-range (1/r^6 cutoff) |
| PQ per-type | 0.001-0.006 | Medium (1/r^5) |
| Per-ring BS | 0.002-0.015 | Varies by ring distance rank |

The BS and HM thresholds are the lowest (most sensitive), consistent
with their 1/r^3 reach.  Dispersion thresholds are extremely low
because the CHARMM switching function tapers them to zero beyond 5A
— they're either present or absent, with no intermediate regime.

## Comparison to Case (1995)

Case calibrated isotropic ring current intensity factors against
single-molecule DFT.  We calibrate the full T2 weight vector against
protein-embedded DFT.  Direct comparison of ring current weights:

| Ring type | Case I (HM) | Our |w_BS|/|w_BS_PHE| | Our |w_HM|/|w_HM_PHE| |
|-----------|------------|---------------------|----------------------|
| PHE | 1.46 | 1.00 | 1.00 |
| TYR | 1.24 | 1.92 | 0.81 |
| TRP-6 | 1.24 | 0.26 | 0.32 |
| HIE | 1.35 | 0.45 | 0.43 |

The HM ratios are in the right ballpark for TYR and HIE.  TRP
weights are lower than Case's because the per-ring kernels (BS_ring0
etc.) absorb the TRP-specific near-field signal that per-type sums
cannot capture.  Case's model had no per-ring decomposition.

## What this means

1. **The weight vector is a table of calibrated constants**, not a
   model artefact.  Each weight maps to a physical parameter (ring
   current intensity, Buckingham coefficient, McConnell anisotropy)
   scaled by normalization.

2. **The constants are element-dependent.**  This is predicted by the
   literature (Boyd & Skrynnikov 2002, Sahakyan & Vendruscolo 2013,
   Saito et al. 2010) but has not been measured in the T2 channel.

3. **Four tiny models (hidden=4) outperform one large model (hidden=64)**
   by 0.21 in weighted R^2.  The physics is element-dependent; a
   pooled model wastes capacity trying to learn what element
   stratification provides for free.

4. **The gating thresholds encode physical ranges**, matching the
   multipolar distance dependence verified in the tensor reality
   tests (BS slope -3.04, PQ slope -5.05).

## Reproduction

```bash
cd learn/src
python3 -c "
from mutation_set.config import load_config
import sys
cfg = load_config('calibration.toml')
sys.path.insert(0, str(cfg.paths.sdk))
from actual_physics.per_element_calibration import run
run(cfg)
"
```

Output: `learn/src/output/actual_physics/calibration/`
