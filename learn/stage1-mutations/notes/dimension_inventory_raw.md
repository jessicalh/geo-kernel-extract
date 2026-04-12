# Dimension Inventory — Raw Working Notes

2026-04-12.  First pass at the complete inventory of dimensions
examined in Stage 1, built from reading the actual results not just
the config.

---

## What the config says vs what the data shows

The extractor computes 10 physical mechanisms generating 55 core
kernels plus 66 per-ring kernels.  But the results reveal dimensions
beyond the kernel layout.

---

## Physical mechanism kernels

**Ring current** (two models of the same effect):
- Biot-Savart: wire-loop line integral. T0+T1+T2. r^-3. 8 ring-type kernels + per-ring BS/RBF-near/RBF-mid/RBF-far/angular-axial/angular-equatorial
- Haigh-Mallion: surface integral. T0+T1+T2. r^-3. 8 ring-type kernels + per-ring HM/HM_H

**Electric field gradient** (two charge sources):
- ff14SB Coulomb: fixed point charges. Traceless symmetric = pure T2. r^-3. Backbone + aromatic decomposition
- MOPAC Coulomb: PM7 Mulliken charges. Same math, different charges. Same decomposition
- APBS: Poisson-Boltzmann solvated field. And DeltaAPBS (solvation - vacuum)

**Bond magnetic anisotropy** (two weighting schemes):
- McConnell: categorical Delta_chi. Asymmetric, T0+T1+T2. r^-3. 5 bond categories
- MOPAC McConnell: Wiberg bond-order weighted. Same math, continuous weights. 5 categories

**Ring susceptibility:** Point magnetic dipole at ring center. r^-3. T0+T1+T2. 1 total kernel + per-ring Chi

**Pi-quadrupole:** Axial quadrupole EFG from ring. Pure T2, traceless. r^-5. 8 ring-type kernels + per-ring PQ

**Dispersion:** van der Waals from ring vertices. r^-6 isotropic, r^-8 anisotropic. 8 ring-type + per-ring DispChi

**H-bond:** Dipolar from H-bond partner. r^-3. 1 total kernel

---

## Dispersion is a bigger story than expected

In the geometry-only forward selection, dispersion kernels enter at
rank 2-5 for C, N, and O.  For nitrogen specifically, 5 of the top 7
kernels after EFG_aro are dispersion.  Dispersion group R-squared for
O is 0.272 — higher than ring current (0.231) or EFG (0.165) in the
geometry-only basis.  Listed originally as a minor r^-6 effect.  It
isn't.

---

## The normalization story is itself a dimension

Raw kernels: top 10 PCs capture 99.3% of variance, predictive plateau
at k=3.  Normalized kernels: top 10 PCs capture only 24%, plateau at
k=46 for hydrogen.  Per-protein normalization strips the magnitude
structure (which is dominated by distance-to-nearest-ring) and reveals
angular structure that lives in the tail of the eigenspectrum.  The
scale factors you add back as scalars carry +0.06 R-squared uniformly
— that's the magnitude information you intentionally stripped.

---

## The ff14SB vs MOPAC EFG comparison is a proper dimension

Cosine similarity = 0.50 for H (nearly orthogonal!), 0.61 for C,
0.90 for N.  Same mathematical kernel, different charges.  For carbon,
MOPAC EFG gets R-squared=0.391 vs ff14SB's 0.224 — a 75% improvement.
But for the geo-only residual, MOPAC EFG on what's left after geometry
explains only 0.056 for C and essentially zero for N and O.  So MOPAC
wins by angular direction, not by adding independent signal after
geometry.

---

## Ridge weight lives in fine angular structure, not dominant PCs

The accumulation_basis found that the ridge weight is spread across 96
PCA components for 90% coverage in H — the top 10 variance PCs carry
0% of the ridge weight.  The prediction doesn't live in the dominant
variance modes.  It lives in the fine angular structure that
normalization reveals.

---

## Dimensions tested and rejected (proven zeros)

- MOPAC valence (sum of Wiberg bond orders): adds +0.000 for every element
- MOPAC molecular dipole: adds +0.000 for every element
- Bond order as scalar: adds +0.000 for every element
- HBond_total: R-squared = 0.002 (correctly zero in mechanical mutants)
- DeltaAPBS_EFG (solvation correction): R-squared = 0.005 (2% of direct EFG)

These are negative results that matter: they show the kernels already
encode everything these scalars know.

---

## Non-kernel dimensions in the calibration

- **Per-protein normalization scale factors**: +0.06 R-squared uniform across elements.  Bridges per-protein and global — the magnitude that normalization stripped is real signal.
- **Mutation type identity** (which aromatic residue was removed): +0.027 for H, +0.121 for C, +0.220 for N, +0.195 for O.  This is Case (1995)'s ring-type-specific intensity factors extended to T2.
- **Spatial decompositions**: RBF distance shells (near/mid/far gaussian windows on BS), angular partitions (axial/equatorial split at magic angle 54.7 deg).  Not new physics — partitions of existing kernel signal by spatial region.

---

## Distance-dependent dimensionality

For all elements, raw kernel space has 3 predictive dimensions at
every distance band (0-4A, 4-8A, 8-12A, 12+A).  But R-squared varies:

| Element | 0-4A   | 4-8A   | 8-12A  | 12+A   |
|---------|--------|--------|--------|--------|
| H       | 0.919  | 0.959  | 0.969  | 0.908  |
| C       | 0.501  | 0.556  | 0.951  | 0.853  |
| N       | 0.661  | 0.224  | 0.906  | 0.818  |
| O       | NaN    | 0.195  | 0.745  | 0.748  |

The same 3 dimensions work everywhere for hydrogen.  Carbon gets
better at distance because the near-field is harder.  Nitrogen
near-field (0-4A, 280 atoms) has 8 predictive dims — the multi-
mechanism complexity shows up only close to the ring.

---

## Still not fully examined

- expected_ensemble_variance.py — connects Stage 1 dimensionality to Stage 2 pose variation (dimensional analysis + measured power laws)
- Per-ring spatial decompositions in detail (RBF, angular basis)
- Cross-kernel prediction correlations from accumulation_basis
- The eigenspectrum block structure (whether eigenvalues come in blocks of 5 matching the L=2 components)
- The element_physics forward selection with MOPAC kernels included (vs geometry-only forward selection above)
