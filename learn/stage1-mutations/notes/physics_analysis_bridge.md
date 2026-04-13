# Physics–Analysis Bridge

2026-04-13.  How the physics dimensions connect to the analysis
dimensions.  55 kernels → 8 physics groups → element-specific
predictive dimensionality.  400 proteins from Stage1Results.

Script: src/actual_physics/full_space_analysis.py
Data: output/actual_physics/full_space/

---

## The two trees

The **physics tree** describes the interactions: 55 kernels in 8
groups, with measured angular independence between them (cosine
matrix).  The tree is element-invariant — the same kernels, same
cosines, same independence structure for H, C, N, O.

The **analysis tree** describes what each element resolves: how
many independent predictive dimensions the ridge regression can
extract, which groups contribute, and how normalisation changes
the picture.  This IS element-dependent.

The bridge: the physics tree defines the available angular space.
The analysis tree shows how much of that space each element can
access.  H accesses 20 dimensions because 20 independent angular
views carry signal for hydrogen.  N accesses 3 because only 3
blurred mixtures carry signal for nitrogen.  The angular space is
the same; the projection onto the DFT target is element-specific.

---

## Why raw = 3 for every element

Raw kernel space has 3 predictive PCA dimensions universally.
This is the magnitude dimension: one L=2 source (5 eigenvalues at
16-24% each, capturing 85% of raw variance) encodes "how far am I
from a ring."  Every kernel knows this.  The 3 predictive raw PCs
are projections of this shared magnitude signal.

After normalisation strips magnitude, the element-specific angular
structure appears: H=20, C=6, N=3, O=12.

The raw→normalised transition is not a trick.  It separates the
question "how far?" (magnitude, answered the same way by every
kernel, protein-size-confounded) from "which direction?" (angular,
answered differently by each kernel family, element-specific).

---

## Per-element analysis

### Hydrogen: 20 normalised dimensions, R² = 0.856

**Which groups contribute:**
Ring current (17 kernels, R² = 0.782 normalised) dominates.
Forward selection: EFG first (0.70), then 5 ring current steps
(+0.07, +0.02, +0.02, +0.02, +0.01).

**Why 20 dimensions:**
Ring current provides the diversity.  17 kernels with within-group
cosine range [0.35, 1.00].  The BS↔HM pairs for the same ring type
are cos=1.00 (identical — 2 models, 1 view).  But BS_PHE ↔
BS_TRP_pyrrole = 0.36 (nearly independent — different ring
geometries produce different T2 directions).

So the 17 ring current kernels collapse to roughly 8-9 independent
views (8 ring types minus the HIS/HID pair that are inactive, plus
RingSusc which is partially independent at cos=0.55).  EFG adds 2-3
(ff14SB and MOPAC are cos=0.58, partially independent).  Dispersion
adds a few more (within-group range [0.34, 0.94]).  Total: ~20
independent angular views that hydrogen can resolve.

**What normalisation costs:**
R² drops from 0.921 to 0.856 (-0.06).  Ring current drops most
(0.885→0.782).  Magnitude is real signal for hydrogen — how far
from the ring IS how strong the ring current is.

### Carbon: 6 normalised dimensions, R² = 0.484

**Which groups contribute:**
EFG dominates (MOPAC EFG = 0.380 normalised, ff14SB = 0.225).
Forward selection: MOPAC EFG → ff14SB EFG → bond aniso → 
dispersion → dispersion → ring current.

**Why 6 dimensions:**
EFG provides 2-3 views (ff14SB/MOPAC are cos=0.61, partially
independent; backbone/aromatic decomposition).  Dispersion adds
1-2 (enters at rank 4-5).  Ring current adds 1-2.  Bond aniso
(MC_aromatic, cos→target=0.39) adds 1 marginal.  Total ~6.

Carbon doesn't resolve the ring current diversity because ring
current R² = 0.121 (weak) — the 17 ring current kernels carry
little signal for carbon, so their angular diversity is noise.

**The charge-polarisation gap:**
Full → geo-only gap = +0.197 (from the main calibration).
MOPAC EFG = 0.380 vs ff14SB = 0.225 — nearly all the gap is
in the EFG charge source.  The polarisation dimension is
specifically a carbon phenomenon in this data.

### Nitrogen: 3 normalised dimensions, R² = 0.267

**Which groups contribute:**
Everything, weakly.  Forward selection: EFG → quadrupole →
ring current → MOPAC bond → dispersion → MOPAC bond →
dispersion × 5.  Six different groups in the first 6 steps.

**Why only 3 dimensions:**
No single group exceeds R² = 0.124.  The largest per-kernel R² is
EFG_aro at 0.061.  Five kernel families each contribute 0.015-0.07.
The individual contributions are too weak to resolve the within-
group diversity.  The 3 PCA directions that survive are blurred
mixtures of all groups — the ridge takes what it can get from 3
noisy combinations, rather than resolving each group separately.

**The near-field exception:**
At 0-4 Å (280 atoms for N in 110-protein data), predictive dims
go from 3 to 4.  Close to the ring, the signal-to-noise per group
improves enough to resolve one more angular direction.

**Dispersion after normalisation:**
0.015 → 0.068 (+0.05).  Dispersion enters 5 of the last 10
forward selection steps.  After normalisation strips magnitude,
dispersion's angular structure for nitrogen emerges from the noise.
But it's still weak — each dispersion step adds 0.007-0.009.

### Oxygen: 12 normalised dimensions, R² = 0.304

**Which groups contribute:**
EFG (0.158 norm) + ring current (0.195) + dispersion (0.234).
Dispersion is the largest group after normalisation.  Forward
selection: EFG → ring current → dispersion × 4 → quadrupole × 2 →
dispersion → ring current.

**Why 12 dimensions:**
Dispersion provides the diversity.  8 kernels with within-group
cosine range [0.33, 0.97].  Disp_PHE ↔ Disp_TRP_perimeter = 0.33
(nearly independent).  The r^-6 cutoff makes dispersion highly
ring-type-specific — different ring sizes/shapes produce different
angular patterns at close range.  Ring current adds 4-5 views.
Quadrupole adds 2-3.  Total ~12.

**The normalisation revelation:**
Dispersion: 0.058 → 0.234 (+0.176).  4x increase.  The largest
single normalisation effect in the entire analysis.

In raw space, the dispersion kernel magnitudes are dominated by
the r^-6 distance falloff.  This magnitude variation is protein-
dependent (ring count, proximity distribution).  It dominates the
raw variance but carries no angular information.  Normalisation
strips it, leaving the angular direction of the dispersion T2 at
each atom — which direction is the nearest ring vertex, and what
angular pattern does the DispChi (dispersion × susceptibility)
tensor create?

For oxygen atoms near aromatic rings, this angular pattern
carries information about the ring's effect on the carbonyl that
raw magnitude hid.

---

## The bridging table

```
Group          |  H                |  C                |  N                |  O
               | R²_n  dims  path | R²_n  dims  path | R²_n  dims  path | R²_n  dims  path
───────────────┼───────────────────┼───────────────────┼───────────────────┼──────────────────
ring_current   | .782  8-9   ████ | .121  1-2   █    | .124  <1    █    | .195  4-5   ██
EFG (ff14+mop) | .701  2-3   ██  | .380  2-3   ███ | .064  <1    █    | .150  2     █
bond_aniso     | .107  1     █   | .025  <1    ·    | .029  <1    ·    | .041  <1    ·
quadrupole     | .030  <1    ·   | .013  <1    ·    | .063  <1    █    | .038  1-2   █
dispersion     | .387  2-3   ██  | .115  1-2   █   | .068  <1    █    | .234  4-5   ███
solvation      | .021  <1    ·   | .017  <1    ·    | .006  0     ·    | .018  <1    ·
───────────────┼───────────────────┼───────────────────┼───────────────────┼──────────────────
ALL norm       | .856  20         | .484  6          | .267  3          | .304  12
ALL raw        | .921  3          | .518  3          | .215  3          | .246  3
```

"dims" = approximate number of independent predictive directions
this group contributes for this element (estimated from within-group
cosine diversity and group R²).  "path" = density in the forward
selection.

The column sums of dims roughly match the PCA-ridge plateau:
H: 8+3+1+0+3+0 ≈ 15-17 (plateau at 20 — some cross-group
interactions add a few more).
C: 2+3+0+0+1+0 ≈ 6.
N: 1+1+0+1+1+0 ≈ 3 (blurred — each is partial, not clean).
O: 5+2+0+2+4+0 ≈ 13 (plateau at 12).

---

## What this means for the thesis

1. **The element-dependent dimensionality is explained by the
   bridging table.**  It's not arbitrary — it follows from which
   physics groups carry signal for each element and how much
   within-group diversity the element can resolve.

2. **Normalisation is element-specific.**  H loses signal (magnitude
   matters).  N and O gain (magnitude confounds).  This is a physics
   result about element electronic structure, not a preprocessing
   choice.

3. **The gap between 55 kernels and the effective dimensions is the
   redundancy structure.**  BS≡HM (cos=1.00), TRP_benzene≈perimeter
   (cos=0.97), backbone≈total (cos=0.81).  The redundancy is
   MEASURED, not assumed.

4. **Nitrogen's 3 blurred dimensions are a genuine limitation of
   geometric kernels for nitrogen.**  The tools see nitrogen poorly
   because nitrogen's electronic structure (paramagnetic-dominated,
   lone pair, peptide bond coupling) responds to multiple geometric
   perturbations simultaneously and none dominates.

5. **Oxygen's dispersion revelation is a finding.**  The r^-6
   angular structure of ring-vertex proximity carries real T2
   information for oxygen that only appears after normalisation.
   This was not expected from the isotropic shift literature, where
   dispersion is negligible.
