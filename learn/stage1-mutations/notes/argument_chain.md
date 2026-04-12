# The Argument Chain — Why Angular Patterns Are Physics

2026-04-12.  The central thesis claim and its grounding.

---

## The chapter question

What does an always-T2 kernel feature extraction tell us about the
physics of shielding?  Aromatic mutant DFTs provide concrete
statistical assessment of same.

---

## The claim

Established shielding physics, expressed as T2 geometric kernels,
survives quantitative statistical assessment against DFT at protein
scale.  The literature established this for isotropic shifts on small
molecules.  We show it holds for the full anisotropic tensor across
69,000 atoms in 110 proteins.

The contribution is not new physics.  It is that the known physics,
when expressed correctly as rank-2 tensors and tested against DFT
shielding tensor differences, produces a quantitative decomposition
into physical sources — per element, per mechanism, with calibrated
coefficients.

---

## The argument chain

Every step is either "cited" or "shown."

### Step 1: Each mechanism has a known multipolar order (CITED)

- Biot-Savart ring current: magnetic dipole, r^-3.
  Pople 1956, Johnson & Bovey 1958.
- Haigh-Mallion ring current: surface integral of magnetic dipole,
  r^-3.  Haigh & Mallion 1979, Case 1995.
- McConnell bond anisotropy: magnetic dipole from bond susceptibility,
  r^-3.  McConnell 1957.
- Coulomb EFG: Buckingham electric field effect, r^-3 (point charge
  sum).  Buckingham 1960.
- Pi-quadrupole: electric quadrupole from ring, r^-5.
  Stone 2013 (T-tensor formalism).
- Ring susceptibility: magnetic dipole from ring bulk susceptibility,
  r^-3.  Same dipolar kernel as McConnell.
- Dispersion: van der Waals, r^-6 isotropic / r^-8 anisotropic.

### Step 2: Different multipolar orders produce different angular symmetries (CITED)

Different multipolar orders produce rank-2 tensors with different
angular dependence.  A magnetic dipole gives (1 - 3cos^2 theta).
An electric quadrupole gives terms involving the ring normal that
the dipole doesn't have.  This is a mathematical property of the
spherical harmonic expansion.  Stone 2013, Ch. 3.

Therefore the T2 direction at an atom is a fingerprint of which
mechanisms contribute.  This is the logical consequence of steps
1 and 2.

### Step 3: The kernels confirm the multipolar predictions (SHOWN)

Three-significant-figure agreement with theory:

| Kernel | Expected slope | Measured slope |
|--------|---------------|---------------|
| Biot-Savart | -3 | -3.04 |
| Haigh-Mallion | -3 | -3.06 |
| Pi-Quadrupole | -5 | -5.05 |

BS-HM magnitude ratio: 1.009 +/- 0.037 (same effect, different
math).  BS-HM T2 cosine: 0.9977 near field, 0.9999 far field
(converge to same dipole limit).

### Step 4: The kernel families are angularly independent (SHOWN)

|cos(BS, EFG)| = 0.684 across 58,193 atom-ring pairs.  Between
collinear (1.0) and random in 5D (0.45).  They share geometric
information (both depend on ring position) but encode different
angular patterns (magnetic dipole vs electric quadrupole).

BS additivity for fused rings: T2(TRP5) + T2(TRP6) = T2(TRP9) at
machine precision.

### Step 5: The DFT target T2 aligns with the correct mechanism per element (SHOWN + CITED)

Mean |cos(EFG, target)| by element:
- H: 0.930 (tight alignment)
- C: 0.863
- N: 0.790
- O: 0.737

Literature prediction: ring current dominates H (Boyd & Skrynnikov
2002 JACS 124:1832), EFG dominates heavy atoms (Sahakyan &
Vendruscolo 2013 JPC-B 117:1989).

Ridge R-squared by physics group confirms per element:
- H: ring current 0.909, EFG 0.770
- C: ring current 0.156, EFG 0.431
- N: multi-mechanism (largest single = 0.325)
- O: mixed (ring current 0.317, EFG 0.196)

The isotropic literature extends to T2 for H and C.  Breaks down
for N and O — those need all kernel families.

### Step 6: The calibrated coefficients are physical constants (SHOWN)

Per-element ridge regression with 55 kernels, per-protein
normalization, kernel scale factors, mutation type identity.
Weighted R-squared = 0.818.

The coefficients are the tensor-valued extension of Case (1995)'s
ring current intensity factors.  They are calibrated constants, not
model parameters.

---

## What is NOT claimed

- No new physics.  All mechanisms cited from literature.
- No prediction.  The output is calibrated constants, not a
  predictor.  R-squared is evidence for the decomposition, not a
  performance metric.
- No proof that the decomposition is unique.  Ridge regression finds
  one set of coefficients.  The claim is that this set is physically
  interpretable, not that it is the only one.

---

## Rigour status of specific claims

### Rigorous (computed, in the data)

- Ridge R-squared per element per physics group
- Forward selection marginal deltas
- Distance slopes (log-log regression on 58K pairs)
- BS-HM magnitude ratio and cosine convergence
- PCA-ridge curves with train/val split
- Eigenspectrum cumulative variance
- Proven-zero scalars (valence, bond order, dipole)

### Needs more careful treatment

- "3 effective dimensions" — the PCA-ridge plateau detector says k=3
  for raw kernels, but the H curve keeps climbing.  The claim is about
  where the rate changes, not a hard floor.  Needs honest framing.

- "Top variance PCs carry 0% ridge weight" — real number from
  accumulation_basis, but the L1-norm denominator matters.  Framing
  needs care.

- "3 dimensions correspond to ring current, EFG, bond anisotropy" —
  asserted from eigenspectrum block structure (3 blocks of 5).  Not
  verified by projecting PCs onto physics groups.  This is work that
  needs doing.

### Honest limitations

- All numbers on 110 proteins (AzimuthalExtraction).  Final 723 will
  sharpen but not change the physics.
- Per-protein ceiling: 0.81 within, 0.35 across.  The gap is real and
  unmeasured by the current analysis.
- Near-field accuracy: wire-loop and surface-integral models are
  approximate below 4 Angstroms.
- HIE: worst ring type (self-fit R-squared = 0.062).
