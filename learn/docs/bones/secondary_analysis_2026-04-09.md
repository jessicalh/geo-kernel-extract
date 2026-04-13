# Secondary Analysis Results — 2026-04-09

## What this is

Five analysis tools that interrogate the calibration pipeline beyond
the headline R² number.  The calibration pipeline fits 91 geometric
kernels to DFT WT-ALA shielding deltas across 723 proteins (447,802
matched atoms) and achieves ridge R² = 0.35, equivariant MLP R² =
0.61.  These tools ask: *why* 0.35, *why* 0.61, *where* does each
ring type sit, and *what* scalar context drives the gap.

All data: `learn/src/output/secondary/`.  All figures:
`output/secondary/figures/`.  To regenerate everything:

```bash
cd learn/src && python3 -m secondary --all --config calibration.toml
cd learn/R  && for f in *.R; do Rscript "$f" ../src/output/secondary; done
```

---

## Tool 1 — Classical vs MOPAC divergence

### What it answers

Where do ff14SB fixed charges and MOPAC PM7 Mulliken charges
produce different T2 angular patterns?  This comparison needs no
DFT ground truth — it maps polarisation effects.

### Figures

- **divergence_cosine_by_element.pdf** — Violin plots of 5D cosine
  similarity between classical and MOPAC T2, split by element (H, C,
  N, O, S) and pair (coulomb shielding, McConnell shielding, EFG
  backbone, EFG aromatic).  The median cosine for coulomb shielding
  is ~0.7 — the total E-field direction changes substantially with
  MOPAC charges.  EFG aromatic is ~0.99 — the local aromatic field
  gradient barely cares which charge model you use.

- **divergence_coulomb_vs_distance.pdf** — Divergence magnitude vs
  distance to the nearest removed ring, for the coulomb shielding
  pair.  The loess fit shows divergence is roughly constant with
  distance — it's a global charge redistribution effect, not a
  local ring effect.

- **divergence_by_ring_type.pdf** — Mean divergence magnitude
  grouped by nearest ring type.  All ring types show similar
  divergence for coulomb and backbone EFG.  No ring type is
  specially sensitive to the charge model.

- **divergence_cosine_by_ring_type.pdf** — Same grouping, cosine
  similarity.  McConnell cosine is lowest (~0.37) — the bond
  anisotropy tensor direction is most sensitive to which charges
  you use.  This makes physical sense: McConnell sums over all
  bonds weighted by their individual properties, so charge-dependent
  bond polarisation changes the angular pattern.

### Key numbers

| Pair | Mean |diff| | Mean cosine | What it means |
|------|-------------|-------------|---------------|
| Coulomb shielding | 4.7 ppm | 0.69 | Total E-field rotates substantially |
| McConnell shielding | 0.52 ppm | 0.37 | Bond anisotropy direction most sensitive |
| EFG backbone | 2.5 ppm | 0.57 | Backbone field gradient moderately different |
| EFG aromatic | 0.02 ppm | 0.99 | Aromatic EFG nearly identical either way |

### Physics

The aromatic EFG is determined by the local ring geometry (charge
positions on the ring itself), which is similar in both charge
models.  The total Coulomb field is determined by ALL charges in
the protein, and MOPAC redistributes them through polarisation.
The McConnell tensor direction depends on individual bond dipole
properties, making it the most charge-sensitive.

---

## Tool 2 — Ring-type conditional analysis

### What it answers

How does the kernel vector behave near different ring types?  How
many atoms see each type?  Does ridge R² vary by ring type
exposure?

### Figures

- **strata_ridge_r2.pdf** — Bar chart of ridge R² by stratum (atoms
  grouped by which ring types they see) and kernel family.  The
  headline: hie_only atoms have R² = 0.24, phe_only = 0.46,
  trp_only = 0.49.  HIE is hardest to model, TRP easiest.

- **strata_kernel_magnitude.pdf** — Top 15 kernels by T2 magnitude,
  faceted by stratum.  Shows which kernels are active in each
  stratum.  EFG aromatic kernels dominate near HIE; ring current
  kernels dominate near TRP.

- **strata_atom_counts.pdf** — Atom counts per stratum.  52K
  hie_only, 62K phe_only, 65K tyr_only, 18K trp_only, 77K no_hie.
  HIE has plenty of data — the low R² is not a sample size issue.

- **strata_protein_composition.pdf** — How many proteins contain
  each ring type.  PHE and TYR are nearly universal; TRP is in
  ~40% of proteins; HIE in ~30%.

### Key numbers

| Stratum | n_atoms | R²_all | R²_ring | R²_efg |
|---------|---------|--------|---------|--------|
| hie_only | 52,238 | 0.240 | 0.107 | 0.186 |
| phe_only | 61,554 | 0.459 | 0.176 | 0.367 |
| tyr_only | 64,625 | 0.394 | 0.205 | 0.349 |
| trp_only | 17,896 | 0.494 | 0.314 | 0.394 |
| no_hie | 77,346 | 0.374 | 0.227 | 0.321 |
| all | 447,802 | 0.355 | 0.184 | 0.308 |

### Physics

HIE's ring current kernels explain almost nothing (R²_ring = 0.107)
because the imidazole ring current is weak (I = -5.16 vs -12.0 for
PHE) and the 5-membered ring geometry produces a spatially compact
field that decays rapidly.  The EFG explains more (0.186) because
HIE→ALA removes two nitrogens and their lone pairs, creating a
large charge redistribution that the Coulomb EFG captures.

TRP is best-modelled because its large fused ring system produces
strong, spatially extended ring current fields that match the
circular-loop BS/HM models well.

---

## Tool 3 — Inter-kernel T2 structure

### What it answers

How many independent angular directions do the 91 kernels actually
provide?  Is the kernel space different for atoms near different
ring types?

### Figures

- **kernel_cosine_heatmap_all.pdf** — 91x91 pairwise cosine
  similarity heatmap across all atoms.  The per-ring-type blocks
  (BS_PHE through PQ_HIE, first 32 kernels) show high internal
  similarity — six-membered ring kernels are nearly parallel.

- **kernel_eigenspectrum.pdf** — Cumulative variance explained by
  eigenvalue rank, for each stratum.  First eigenvalue captures
  ~80%, first 3 capture ~98%.  All strata look similar — the
  low-rank structure is universal, not ring-type-specific.

- **kernel_hie_vs_phe.pdf** — Side-by-side heatmaps for HIE-only
  vs PHE-only atoms.  The block structure is similar but the
  EFG kernels (rows 49-54) show higher inter-correlation in the
  HIE stratum — consistent with EFG dominance there.

### Key number

**Effective dimensionality: 3 (at 98% variance) from 91 kernels.**

This does NOT mean the kernels are useless.  It means 91 kernels
span ~3 angular directions *with fixed weights*.  The MLP uses
scalar context to choose DIFFERENT linear combinations per atom,
which is why R² jumps from 0.35 (fixed) to 0.61 (environment-
dependent).  The value of having 91 kernels is that they give the
MLP a rich menu to select from — not that they provide 91
independent directions.

### Physics

The 3 effective directions correspond roughly to: (1) the dominant
ring current / susceptibility T2 pattern, (2) the EFG T2 pattern,
(3) the McConnell / bond anisotropy T2 pattern.  These are the
three distinct physical source geometries: rings (sparse, strong),
charges (dense, moderate), bonds (dense, weak per bond).

---

## Tool 4 — Per-ring-type R² ablation

### What it answers

Is HIE bad because its geometric kernels are wrong, or because its
DFT contrast is unreliable?  Separates kernel quality from DFT
quality by fitting each ring type's atoms with only its own 4
kernels (BS_RT, HM_RT, Disp_RT, PQ_RT).

### Figures

- **ablation_self_fit.pdf** — THE MONEY PLOT.  Self-fit R² per ring
  type.  HIE = 0.062, the worst.  PHE = 0.098, TYR = 0.106,
  TRP_perimeter = 0.125.  The 4 ring-type-specific kernels explain
  less of the DFT delta at HIE than at any other ring type.

- **ablation_self_vs_all.pdf** — Self-fit (4 kernels) vs full fit
  (all 91 kernels) per ring type.  The gap between the bars shows
  how much other kernel families contribute.  For HIE: self = 0.062,
  all = 0.308, so non-ring kernels (primarily EFG) provide 80% of
  the explanatory power.  For TRP: self = 0.114, all = 0.383, so
  ring kernels provide 30% — still the minority, but much larger.

- **ablation_r2_by_stratum.pdf** — Per-stratum R² broken by kernel
  family (ring, bond, total, EFG, per-ring).  Shows which kernel
  family matters where.

- **ablation_target_magnitude.pdf** — Mean DFT delta T2 magnitude
  per ring type.  TRP (1.93-1.98 ppm) has the strongest signal,
  HIE (1.49 ppm) the weakest.  The weak signal at HIE is
  consistent with a weak ring current perturbation.

### Key numbers — self-fit

| Ring type | n_atoms | R²_self (4 kernels) | R²_all (91) |
|-----------|---------|---------------------|-------------|
| PHE | 190,939 | 0.098 | 0.359 |
| TYR | 187,272 | 0.106 | 0.350 |
| TRP_benzene | 73,161 | 0.114 | 0.383 |
| TRP_pyrrole | 74,089 | 0.095 | 0.383 |
| TRP_perimeter | 73,664 | 0.125 | 0.384 |
| HIE | 142,225 | **0.062** | **0.308** |

### Physics

HIE's ring current model is the weakest — the circular loop
approximation (BS) and surface integral (HM) model a symmetric
current that doesn't match the asymmetric imidazole electron
distribution.  The 5-membered ring with two nitrogens at specific
positions creates a field that depends on azimuthal angle φ (where
the atom sits relative to the N atoms), which the current kernels
don't capture.

But self-fit R² is low for ALL ring types (0.06-0.13), not just
HIE.  The ring current kernels alone are never sufficient — they
always need the EFG and other kernel families.  This is a finding
about the nature of the WT-ALA shielding delta: removing an
aromatic residue changes EVERYTHING (charges, bonds, geometry), not
just the ring current.

---

## Tool 5 — Scalar group interaction analysis

### What it answers

Which scalar features drive the ridge R² = 0.32 → MLP R² = 0.61
gap?  For each of 20 named scalar groups, adds scalar × kernel
interaction terms to a base ridge (top 10 kernels) and measures
the R² improvement.  This shows what scalar context the MLP
exploits, in a fully linear framework.

### Figures

- **scalar_interaction_impact.pdf** — Horizontal bar chart of R²
  delta per scalar group, sorted by impact.  ring_proximity
  dominates (+0.095), followed by residue_type (+0.082),
  kernel_scales (+0.077), element (+0.076), bond_orders (+0.070).

- **scalar_cumulative_r2.pdf** — R² as groups are added cumulatively
  in order of impact.  Shows diminishing returns and where the
  linear ceiling is.

- **scalar_width_vs_impact.pdf** — Scatter of group width (number
  of scalars) vs impact.  element (5 scalars) has nearly the same
  impact as kernel_scales (91 scalars) — compact features can be
  as informative as large ones.

### Key numbers — interaction R² (base = 0.320)

| Rank | Group | Width | R² delta |
|------|-------|-------|----------|
| 1 | ring_proximity | 49 | +0.095 |
| 2 | residue_type | 20 | +0.082 |
| 3 | kernel_scales | 91 | +0.077 |
| 4 | element | 5 | +0.076 |
| 5 | bond_orders | 4 | +0.070 |
| 6 | delta | 4 | +0.057 |
| 7 | mopac_mcconnell | 6 | +0.051 |
| 8 | mcconnell | 6 | +0.050 |
| 9 | mopac_electronic | 2 | +0.050 |
| 10 | t0_magnitudes | 6 | +0.038 |
| ... | ... | ... | ... |
| 20 | hbond | 3 | +0.009 |

### Physics

**ring_proximity** is the top group because it gives the model
per-ring spatial context: (z, rho, theta, mcconnell_factor,
exp_decay, disp_scalar, disp_contacts) for the 6 nearest rings.
This is where the atom sits relative to each ring, encoded as 7
geometry features per ring.  When the model knows this, it can
weight the ring current kernels appropriately — closer rings get
more weight, rings above vs beside the atom get different weight.

This is also where the proposed azimuthal angle φ (cos φ, sin φ per
ring) would land.  The ring_proximity block is already the most
impactful scalar group; adding in-plane direction would extend it
from 7 to 9 features per ring.

**kernel_scales** confirming the earlier finding: the MLP needs to
know the magnitude that per-protein normalization stripped.  Without
it, the model can't distinguish a protein with one nearby ring from
a protein with ten.

**element** being as impactful as kernel_scales (with 5 vs 91
features) shows that knowing H vs C vs N vs O is fundamentally
important — different atom types respond to the same geometric
kernel differently (H atoms on N-H bonds vs C atoms on the backbone
see the ring current from different geometric positions and with
different electronic sensitivity).

---

## Azimuthal angle φ — the proposed next feature

See `spec/AZIMUTHAL_ANGLE.md` for the full specification.

The secondary analysis motivates adding cos(φ), sin(φ) per (atom,
ring) pair to ring_contributions.npy.  φ is the in-plane angle
relative to vertex 0 of the ring.  For symmetric rings (PHE, TYR)
it carries little information.  For asymmetric rings (HIE
imidazole, with nitrogen positions breaking rotational symmetry) it
tells the model where the atom sits relative to the nitrogen lone
pairs.

ring_proximity is already the most impactful scalar group (+0.095).
Adding 2 scalars per ring (cos φ, sin φ) extends it to 10 features
per ring.  The hypothesis: HIE-stratum R² improves because the
model can now distinguish the nitrogen side from the carbon side;
PHE/TYR barely change.

C++ change: 6 lines of Eigen in BiotSavartResult.cpp, two new
fields on RingNeighbourhood, two new columns in
ring_contributions.npy (57 → 59).  Full re-extraction required.

---

## How to reproduce

```bash
# Python analysis (10 min on 723 proteins)
cd learn/src
python3 -m secondary --all --config calibration.toml

# R figures (30 seconds)
cd learn/R
Rscript divergence.R ../src/output/secondary
Rscript ring_strata.R ../src/output/secondary
Rscript kernel_structure.R ../src/output/secondary
Rscript ablation.R ../src/output/secondary
Rscript scalar_ablation.R ../src/output/secondary

# Quick test (1 minute on 10 proteins)
cd learn/src
python3 -m secondary --all --config calibration.toml --max-proteins 10
```

## Files

```
learn/src/secondary/
    loader.py              — shared protein iteration + ring-type tagging
    divergence.py          — tool 1: classical vs MOPAC T2 divergence
    ring_strata.py         — tool 2: ring-type conditional analysis
    kernel_structure.py    — tool 3: inter-kernel T2 structure
    ablation.py            — tool 4: per-ring-type R² ablation
    scalar_ablation.py     — tool 5: scalar group interaction analysis
    __main__.py            — CLI

learn/R/
    common.R               — shared theme, palettes
    divergence.R           — tool 1 figures
    ring_strata.R          — tool 2 figures
    kernel_structure.R     — tool 3 figures
    ablation.R             — tool 4 figures
    scalar_ablation.R      — tool 5 figures

spec/
    AZIMUTHAL_ANGLE.md     — proposed next feature (cos φ, sin φ)
```
