# The Tensor Decomposition of the Shielding Perturbation: 723 Mutants

## What this document is

The NMR shielding extractor computes 121 geometric kernels — rank-2
spherical tensors (L=2, 5 components) — at each atom in a protein,
decomposing the shielding perturbation from aromatic ring deletion
into its physical sources.  Before using these kernels on NMR
problems, we verify that they reproduce 28 known properties of the
shielding interaction, drawn from the NMR and physics literature.

The tests fall into three categories:

- **Analytical** (ridge regression on the 110-protein mutation set)
- **Tensor-direct** (geometric properties of the kernel T2 vectors,
  no fitting)
- **Structural** (inter-kernel relationships, symmetry properties)

All tests are on WT-ALA mechanical mutants: backbone fixed, aromatic
sidechain replaced with alanine.  The DFT shielding tensor delta is
the reference.  110 proteins, 69,080 matched atoms.

---

## Distance dependence

The geometric kernels encode specific multipolar interactions.  The
T2 magnitude falloff with distance from the nearest ring is measured
directly (log-log slope on 58,193 atom-ring pairs, r > 2 A).

| Kernel | Physical source | Expected slope | Measured slope |
|--------|----------------|----------------|---------------|
| Biot-Savart | Magnetic dipole (ring current) | -3 | **-3.04** |
| Haigh-Mallion | Surface integral (ring current) | -3 | **-3.06** |
| Pi-Quadrupole | Electric quadrupole | -5 | **-5.05** |

Three-significant-figure agreement with multipolar theory.  The
distance dependence is not fitted — it is encoded in the geometry
of the kernel computation (Biot-Savart law, Stone T-tensor
formalism).

*Johnson & Bovey 1958; Buckingham 1959; Stone 2013.*

---

## Angular dependence

The Biot-Savart kernel magnitude depends on polar angle theta
relative to the ring normal.  Axial atoms (above/below ring) see
stronger fields than equatorial atoms (in-plane).

| Angle range | Mean BS magnitude | Relative |
|-------------|------------------|----------|
| 0-17 deg    | 0.0385           | 1.00     |
| 17-34 deg   | 0.0331           | 0.86     |
| 34-52 deg   | 0.0237           | 0.62     |
| 52-69 deg   | 0.0172           | 0.45     |
| 69-90 deg   | 0.0139           | 0.36     |

Monotonic decrease, 2.8x axial-to-equatorial ratio.  Consistent
with the (1 - 3 cos^2 theta) angular pattern of a magnetic dipole.

*Pople 1956; Case 1995.*

---

## Element dependence

The T2 shielding perturbation from ring deletion depends on the
probe atom's element.  This is predicted by the literature:
ring current dominates proton shielding (Boyd & Skrynnikov 2002);
electric field gradient dominates heavy atom shielding (Sahakyan &
Vendruscolo 2013).  The prediction is for isotropic shifts; we
test it in the T2 (anisotropic) channel.

### Ridge R-squared by physics group

| Element | n_atoms | Ring current | EFG   | Bond aniso | All   |
|---------|---------|-------------|-------|-----------|-------|
| H       | 35,080  | **0.909**   | 0.770 | 0.139     | 0.935 |
| C       | 20,673  | 0.156       |**0.431**| 0.041   | 0.546 |
| N       | 6,100   | **0.325**   | 0.100 | 0.070     | 0.434 |
| O       | 6,927   | **0.317**   | 0.196 | 0.055     | 0.387 |
| pooled  | 69,080  | 0.213       | 0.318 | 0.027     | 0.385 |

The pooled view (bottom row) suggests EFG dominance.  The per-
element view shows this is a composition artifact: H atoms respond
predominantly to ring current, C atoms to EFG, and the pooled
number reflects the element mixture at each distance.

### Tensor alignment by element (no fitting)

Mean cosine similarity between EFG_aro kernel T2 and DFT target T2:

| Element | Mean cos | Std   |
|---------|----------|-------|
| H       | -0.930   | 0.106 |
| C       | -0.863   | 0.239 |
| N       | -0.790   | 0.323 |
| O       | -0.737   | 0.347 |

The negative sign is correct: the target is the WT-ALA delta
(effect of removing the ring), and the EFG kernel measures the
ring's contribution.  The alignment is tightest for H (std 0.11)
and loosens through C, N, O — reflecting the decreasing importance
of the electrostatic mechanism for each element.

Per-element R-squared with all kernels vs pooled:
pooled R-squared = 0.385, weighted per-element = 0.718,
gap = **0.333**.

*Boyd & Skrynnikov 2002; Sahakyan & Vendruscolo 2013;
Saito et al. 2010; Xu & Case 2002.*

---

## Nitrogen

Nitrogen does not fit the simple ring-current-vs-EFG dichotomy.
No single kernel family dominates (largest single-kernel R-squared
= 0.080).  Five mechanisms each contribute 0.015-0.080 in forward
selection: EFG, pi-quadrupole, Biot-Savart, McConnell bond
anisotropy, ring susceptibility.

This is consistent with:

- 15N CSA is large (~166 ppm) and site-variable (+/-20 ppm),
  dominated by backbone angles and H-bonding (Yao et al. 2010;
  Hall & Fushman 2006).  In mechanical mutants these are held
  constant; the residual comes from secondary mechanisms.

- 15N shielding is paramagnetic-term-dominated (Saito et al. 2010),
  making it sensitive to multiple perturbation channels.

- The isotropic ring current effect on 15N is negligible
  (<0.6% of variance, Han et al. 2011), but this does not
  preclude an anisotropic (T2) effect.  Hall & Fushman (2006)
  noted aromatic residues as 15N CSA outliers in protein GB3.

- MopacMC_total entering the N forward selection at rank 4
  reflects the electronic coupling through the peptide bond
  (partial double-bond character, Wiberg bond order ~1.3-1.5).

---

## Model agreement: Biot-Savart and Haigh-Mallion

The two ring current models compute the same physical effect with
different mathematics (wire-loop segments vs surface integral).
Three comparisons:

| Property | BS | HM | Agreement |
|----------|----|----|-----------|
| T2 distance slope | -3.04 | -3.06 | identical |
| R-squared for H atoms | 0.772 | 0.784 | within 2% |
| T2 magnitude ratio (per atom) | — | 1.009 +/- 0.037 | 1% |

The BS-HM T2 angular agreement (cosine similarity):
- Near-field (<5 A): mean |cos| = 0.9977
- Far-field (>8 A): mean |cos| = 0.9999

They converge to the same T2 direction in the far field (both
are magnetic dipoles at long range) with 0.2% disagreement in
the near field (where wire-loop vs surface geometry matters).

For the calibration pipeline, keeping both provides redundancy
but not independent angular information.

*Case 1995 found no significant difference between JB and HM for
isotropic shifts.  We extend this to the T2 channel.*

---

## Ring-type dependence

### Self-fit R-squared (4 ring-specific kernels per type)

| Ring type | R-squared | Target magnitude |
|-----------|----------|-----------------|
| TRP_perimeter | 0.125 | 1.93 ppm |
| TRP_benzene | 0.114 | 1.98 ppm |
| TYR | 0.106 | — |
| PHE | 0.098 | — |
| TRP_pyrrole | 0.095 | — |
| HIE | **0.062** | 1.49 ppm |

PHE and TYR ring currents are similar (R-squared ratio 0.84),
consistent with similar ring current intensities (Case 1995:
I_PHE = 1.46, I_TYR = 1.24 in the HM model).

HIE (imidazole) is worst: the 5-membered ring has weak, asymmetric
current that the circular-loop Biot-Savart model approximates
poorly.  Ring current self-fit for all types is low (0.06-0.13);
the ring current kernels are never self-sufficient — they always
need the EFG and other families.

---

## Mechanical mutant controls

Two negative controls verify that the kernels respond to the
mutation, not to the unchanged backbone:

| Test | Kernel | R-squared | Expected |
|------|--------|----------|----------|
| H-bond unchanged | HBond_total | 0.002 | ~0 |
| Solvation minor | DeltaAPBS_EFG | 0.005 | << EFG_aro (0.230) |

H-bond shielding (backbone N-H...O=C) does not change in a
mechanical mutant; the kernel correctly sees nothing.  The
solvation correction (APBS continuum minus vacuum Coulomb) is
2% of the direct aromatic EFG, confirming that the aromatic
shielding perturbation is a local, through-space effect.

---

## Inter-kernel geometry

### Three effective dimensions

The 91-kernel T2 space has 3 effective dimensions at 98% variance
(eigenspectrum analysis).  These correspond to three physical
source geometries:

1. Ring current magnetic dipole (BS, HM, ring susceptibility)
2. Electric field gradient (Coulomb EFG, MOPAC EFG)
3. Bond magnetic anisotropy (McConnell)

### BS-EFG angular independence

Mean |cos(BS, EFG)| = 0.684 across 58,193 atom-ring pairs.
Between collinear (1.0) and random in 5D (0.45).  The two
kernel families share geometric information (both depend on
the ring position) but encode different angular patterns
(magnetic dipole vs electric quadrupole).

### BS additivity for fused rings

For tryptophan's fused ring system:
T2(TRP5) + T2(TRP6) = T2(TRP9) at machine precision.
The wire-segment Biot-Savart law is exactly additive for
non-overlapping current loops.

---

## Progressive calibration (physics_calibration.py)

Per-protein normalization restores angular isolation.  Two fair
categorical scalars — kernel scale factors (what normalization
stripped) and mutation type (which residue was removed) — are tested
as ridge interaction terms.  Three continuous scalars (MOPAC valence,
bond order, molecular dipole) are tested and contribute zero.

| Element | Base (55 norm) | + scales | + mutation type | Fair set |
|---------|---------------|----------|----------------|----------|
| H       | 0.875         | +0.063   | +0.027         | **0.949** |
| C       | 0.520         | +0.069   | +0.121         | **0.691** |
| N       | 0.450         | +0.059   | +0.220         | **0.704** |
| O       | 0.407         | +0.048   | +0.195         | **0.632** |
| Weighted | 0.684        |          |                | **0.818** |

Kernel scales add ~0.06 uniformly — the magnitude that per-protein
normalization stripped is real signal.  Mutation type adds +0.027 for
H but +0.12-0.22 for heavy atoms — knowing which ring was removed
is a physical variable (Case 1995 calibrated separate intensity
factors per ring type).

Valence, bond order, and molecular dipole each add +0.000 for every
element.  The kernels already encode everything these scalars know.

## Limitations

1. **Per-protein ceiling.**  Per-protein ridge R-squared = 0.81;
   global (cross-protein) ridge = 0.35; MLP = 0.61.  The gap
   represents per-protein structural variation (bulk electrostatic
   environment, protein shape) that local geometric kernels do not
   and should not capture.

2. **Near-field accuracy.**  The BS/HM models are approximate
   at <4 A.  The wire-loop and surface-integral models diverge
   from the true pi-electron current distribution at close range.
   Per-ring kernels partially address this but the functional
   form is fixed.

3. **HIE (imidazole).**  The worst-modeled ring type (self-fit
   R-squared = 0.062).  The symmetric circular-loop approximation
   fails for the asymmetric 5-membered ring with two nitrogen
   positions.

4. **Isotropic gap.**  These kernels describe the T2 (anisotropic)
   shielding perturbation.  The T0 (isotropic shift) requires the
   Buckingham A coefficient and ring current intensity as additional
   parameters, fitted separately.

5. **No structural relaxation.**  Mechanical mutants hold the
   backbone fixed.  Real mutations cause structural rearrangement
   that changes the kernel geometry.  The calibrated kernels apply
   to a fixed structure, not to a mutation experiment.

---

## Reproduction

```bash
cd learn/src

# Per-element physics (Table 1, forward selection)
python3 -c "
from mutation_set.config import load_config
from secondary.loader import setup_sdk
from actual_physics.element_physics import run
cfg = load_config('calibration.toml'); setup_sdk(cfg); run(cfg)
"

# Analytical realities 13-20
python3 -c "
from mutation_set.config import load_config
from secondary.loader import setup_sdk
from actual_physics.twenty_realities import run
cfg = load_config('calibration.toml'); setup_sdk(cfg); run(cfg)
"

# Tensor-direct realities T21-T28
python3 -c "
from mutation_set.config import load_config
from secondary.loader import setup_sdk
from actual_physics.tensor_realities import run
cfg = load_config('calibration.toml'); setup_sdk(cfg); run(cfg)
"
```

---

## References

- Buckingham, A.D. (1960) Can. J. Chem. 38, 300
- Boyd, J. & Skrynnikov, N.R. (2002) JACS 124, 1832
- Case, D.A. (1995) J. Biomol. NMR 6, 341
- Giessner-Prettre, C. & Pullman, B. (1987) Q. Rev. Biophys. 20, 113
- Hall, J.B. & Fushman, D. (2006) JACS 128, 7855
- Han, B. et al. (2011) J. Biomol. NMR 50, 43
- Johnson, C.E. & Bovey, F.A. (1958) J. Chem. Phys. 29, 1012
- Pople, J.A. (1956) J. Chem. Phys. 24, 1111
- Sahakyan, A.B. & Vendruscolo, M. (2013) JPC-B 117, 1989
- Saito, H., Ando, I. & Ramamoorthy, A. (2010) Prog. NMR Spectrosc. 57, 181
- Xu, X.P. & Case, D.A. (2002) Biopolymers 65, 408
- Yao, L. et al. (2010) JACS 132, 10866
