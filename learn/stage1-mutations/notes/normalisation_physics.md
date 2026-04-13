# Normalisation Physics

2026-04-13.  What per-protein normalisation does and doesn't change.
400 proteins from Stage1Results.

---

## The cosine matrix is invariant to normalisation

Per-protein normalisation scales each kernel's magnitude within a
protein.  It does not rotate the T2 vectors.  Cosine similarity
divides out magnitude.  Therefore: the inter-kernel cosine matrix
is IDENTICAL between raw and normalised space (verified to 2 dp
across all 55x55 entries for all elements).

The angular independence structure reported in
physics_group_grounding.md is intrinsic to the kernel geometry.
It is not an artifact of magnitude.  The relationships between
physics groups are a property of the physics.

## What normalisation changes: per-group R²

Normalisation changes how much each group contributes to ridge
prediction, not the angular relationships between them.

### Per-group R² raw → normalised (delta)

| Group | H | C | N | O |
|-------|---|---|---|---|
| ring_current | 0.885→0.782 (-0.10) | 0.129→0.121 (-0.01) | 0.122→0.124 (+0.00) | 0.198→0.195 (-0.00) |
| ff14sb_efg | 0.718→0.660 (-0.06) | 0.223→0.225 (+0.00) | 0.073→0.064 (-0.01) | 0.168→0.158 (-0.01) |
| mopac_efg | 0.760→0.701 (-0.06) | 0.392→0.380 (-0.01) | 0.072→0.062 (-0.01) | 0.158→0.150 (-0.01) |
| bond_aniso | 0.116→0.107 (-0.01) | 0.024→0.025 (+0.00) | 0.029→0.029 (+0.00) | 0.041→0.041 (+0.00) |
| mopac_bond | 0.113→0.108 (-0.01) | 0.020→0.021 (+0.00) | 0.032→0.033 (+0.00) | 0.040→0.041 (+0.00) |
| quadrupole | 0.019→0.030 (+0.01) | 0.008→0.013 (+0.01) | 0.069→0.063 (-0.01) | 0.030→0.038 (+0.01) |
| dispersion | 0.397→0.387 (-0.01) | 0.100→0.115 (+0.02) | 0.015→0.068 **(+0.05)** | 0.058→0.234 **(+0.18)** |
| solvation | 0.012→0.021 (+0.01) | 0.006→0.017 (+0.01) | 0.002→0.006 (+0.00) | 0.005→0.018 (+0.01) |
| ALL | 0.921→0.856 (-0.06) | 0.518→0.484 (-0.03) | 0.215→0.267 **(+0.05)** | 0.246→0.304 **(+0.06)** |

### H: normalisation hurts (-0.06 ALL)

Magnitude is real signal for hydrogen.  How far an H atom is from
the nearest ring determines how strong the ring current is, and
that's useful information the ridge can exploit.  Normalisation
strips this, costing 0.10 in ring_current R² and 0.06 overall.

This is the trade-off of normalisation: you lose per-protein
magnitude for the sake of angular isolation.

### N and O: normalisation helps (+0.05, +0.06 ALL)

For heavy atoms, raw magnitude is confounded with protein size
and ring count.  Stripping this reveals angular signal that was
buried.  N goes from 0.215 to 0.267, O from 0.246 to 0.304.

### Dispersion explosion for O: 0.058 → 0.234

The headline.  A 4x increase.  After normalisation, dispersion
becomes the dominant group for oxygen (0.234 vs ring current 0.195
vs EFG 0.158).

In raw space, the dispersion kernel magnitudes are dominated by
the 1/r^6 distance falloff — atoms near rings have huge values,
atoms far away have zero.  This magnitude variation is protein-
dependent (how many rings, how close).  Normalisation strips it,
leaving the ANGULAR pattern of dispersion: which direction is the
nearest ring vertex, and how does the DispChi (dispersion ×
susceptibility) tensor point?

For oxygen atoms near aromatic rings, this angular pattern
carries information about the ring's effect on the carbonyl that
raw magnitude hides behind protein-to-protein scale variation.

### Dispersion for N: 0.015 → 0.068

Smaller but same physics.  Factor of 4.5x.  Dispersion becomes
meaningful for nitrogen after normalisation.

## Interpretation

Per-protein normalisation is not a preprocessing trick.  It is
a physics operation: it separates "how far" (magnitude, shared
between kernels, protein-dependent) from "which direction"
(angular, kernel-family-specific, element-dependent).

For hydrogen, "how far" is useful signal and normalisation costs.
For heavy atoms, "how far" is confounded with nuisance variables
and normalisation gains.  The optimal strategy is element-dependent
— which is why per-element analysis is necessary.

The dispersion finding is concrete: oxygen's angular sensitivity
to nearby ring vertices is a real dimension that only becomes
visible after normalising away the magnitude scale.  In the thesis,
this should be presented as a finding about what the tool CAN see
when you remove what it sees TOO MUCH of.

## Forward selection after normalisation (groups only)

Bug in kernel name printing (shows "EFG_aro" for all due to
variable reuse).  Groups are correct.  Fix in next run.

**H normalised:** MopacEFG → ring_current × 4 → solvation →
ring_current → dispersion × 2 → quadrupole → dispersion × 2 →
ring_current → quadrupole.  Dominated by ring current after EFG.

**C normalised:** MopacEFG → ff14sb_efg → bond_aniso →
dispersion × 2 → ring_current → dispersion → ring_current × 3 →
dispersion × 2 → bond_aniso.  More balanced, dispersion at rank 4.

**N normalised:** ff14sb_efg → quadrupole → ring_current →
mopac_bond × 2 → dispersion × 6 → solvation → ring_current × 2 →
bond_aniso.  Dispersion dominates the middle (6 of 15 slots).

**O normalised:** ff14sb_efg → ring_current → dispersion × 4 →
quadrupole × 2 → dispersion → ring_current × 2 → dispersion →
quadrupole → mopac_bond → mopac_efg.  Dispersion dominates (5 of
first 10 slots).
