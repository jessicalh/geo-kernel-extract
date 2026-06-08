# McConnell Structured Grounding

This is the build-and-defend grounding for the new McConnell feature. It builds on
`mcconnell_grounding_agent1.md` and the 2026-06-06 controlling spec. It does not
reopen the settled decisions.

The decision is: replace the current three-term `M` with a point-dipole
propagator acting on a source susceptibility tensor; emit the full even
rank-2 tensor decomposition, `0e + 1e + 2e`, for every source category and for
both source-strength channels; use MOPAC Wiberg bond order only as a calibrated
source-strength channel; and set aromatic McConnell contributions to zero while
the ring-current kernels are active.

## Why This Kernel Exists

A protein nucleus feels this effect because the spectrometer field induces
electronic currents and magnetic moments in nearby chemical groups. If the
group has an anisotropic magnetic susceptibility, the induced moment is not a
featureless scalar. Its dipole field reaches a nearby nucleus and changes the
local shielding tensor. Case states this biomolecular shielding picture
directly: a remote chemical group with susceptibility tensor `chi` gives a
McConnell shielding contribution, and isotropic susceptibility does not
contribute to the isotropic shift; the observed effect depends on susceptibility
anisotropy. Aromatic and unsaturated groups are large examples, but peptide
groups also matter in proteins. [Case 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3877577/)

This earns a place among the Three because it is a named through-space mechanism
with a low-order tensor law. It is not the same as the ring-current calculation:
ring current is the distributed aromatic-current path, while this McConnell
feature is the local bond/group susceptibility path. It is also not the charge
or EFG kernel: the source is an induced magnetic moment, not a static electric
charge. Worcester and Pauling both make the protein-scale point that peptide
planes have diamagnetic anisotropy, and Tjandra et al. show that a diamagnetic
protein susceptibility tensor can be explained from peptide-bond and aromatic
side-chain contributions. [Worcester 1978](../../../references-text/worcester-1978-structural-origins-diamagnetic-anisotropy-proteins-text-1.txt),
[Pauling 1979](../../../references-text/pauling-1979-diamagnetic-anisotropy-peptide-group-text-1.txt),
[Tjandra et al. 1996](../../../references-text/tjandra-1996-field-dependence-nh-j-splittings-ubiquitin-text-3.txt)

Plainly: a protein is full of peptide bonds, carbonyls, side-chain unsaturation,
and aromatic rings. Those groups respond anisotropically to the magnetic field.
The ring part is handled by Biot-Savart and Haigh-Mallion. The non-aromatic
bond/group anisotropy part belongs here.

## What We Compute

For each target nucleus `i` and each accepted source group or bond `s`, compute
a rank-2 shielding-response matrix

```text
A_is = D(r_is) Q_s
```

where:

```text
r_is = x_i - c_s
R    = |r_is|
n    = r_is / R
D(r) = (3 n n^T - I) / R^3
```

`D(r)` is the point-dipole propagator. `Q_s` is the traceless anisotropic part
of the source susceptibility tensor. The matrix convention is column-vector:
`m_s = Q_s B0`, then `B_ind = D(r) m_s`, so the response to `B0` is `D Q`.
This ordering matters because `D` and `Q` need not commute; the product can
carry a real antisymmetric `1e` part.

The physical sign and the `4 pi`, SI/cgs, and ppm conversion constants are not
hard-coded into the emitted unit shape. The scale is pinned from defended
literature/physics values, with source scatter disclosed. The emitted feature is
a geometry-and-source-shape basis in `Angstrom^-3`; the sample is too small to
learn free per-category scales.

The mathematical anchor is the same point-dipole susceptibility form used in
the McConnell/PCS literature. Suturina and Kuprov write the point shift tensor
as the dipolar matrix acting on `chi`, then show the isotropic PCS scalar as a
rank-2 spherical-harmonic contraction. Case writes the same biomolecular
remote-group shielding expression. [Suturina and Kuprov 2016](../../../references-text/suturina-2016-pseudocontact-shifts-from-mobile-spin-labels-text-2.txt),
[Case 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3877577/)

## Source Categories

Use the existing topology categories, with one explicit correction: aromatic
sources are kept as a schema category but set to zero when the ring-current
kernels are active.

Recommended category names for the new emit:

```text
peptide_co
peptide_cn
backbone_other
sidechain_co
sidechain_other
disulfide
aromatic_zeroed
```

The source center `c_s` is the bond midpoint unless a later, primary-verified
source table gives a different physical origin for that category. Do not invent
category-specific offsets. The source axis `u_s` is the local bond or group axis
used by the category. Rhombic source terms require a complete local orthonormal
frame `(u_s, v_s, w_s)`, not just a bond direction.

## Source Susceptibility Tensor

For an axial source category:

```text
Q_s = dchi_ax,c (u_s u_s^T - I/3)
```

For a source category with a primary-verified rhombic value and a defensible
local frame:

```text
Q_s = dchi_ax,c (u_s u_s^T - I/3)
    + (dchi_rh,c / 2) (v_s v_s^T - w_s w_s^T)
```

Here `dchi_ax,c` is the axial susceptibility anisotropy for category `c`, and
`dchi_rh,c` is the rhombicity. In principal-axis terms, this construction gives
`dchi_ax = chi_u - (chi_v + chi_w)/2` and `dchi_rh = chi_v - chi_w`.

For the emitted basis fields, use a unit-strength version of this tensor shape
unless a separate literature-scaled diagnostic is explicitly requested:

```text
Qhat_s = (u_s u_s^T - I/3)
       + (rho_rh,c / 2) (v_s v_s^T - w_s w_s^T)

rho_rh,c = dchi_rh,c / dchi_ax,c
```

If `dchi_rh,c` is not primary-verified, set the rhombic part absent, not guessed.
Most organic source tables are axial or under-specified. That absence should be
visible in metadata:

```text
rhombic_status = absent_no_primary_table
```

This preserves the field shape while preventing a fabricated rhombic tensor from
becoming a thesis claim. Suturina and Kuprov's PCS expression uses five
irreducible susceptibility components and states that those components may be
expressed as axiality, rhombicity, and Euler angles. That supports carrying
rhombicity where the data actually exist; it does not support inventing it.
[Suturina and Kuprov 2016](../../../references-text/suturina-2016-pseudocontact-shifts-from-mobile-spin-labels-text-2.txt)

## Pair Calculation

For each target atom `i`:

1. Initialize one `3 x 3` accumulator for each `(category, channel)` pair.
2. Visit source bonds/groups within the configured cutoff.
3. Apply the same geometric exclusion logic as the other through-space kernels:
   singularity guard, self-source exclusion, and near-field/source-extent filter.
   These filters protect the calculation, but they do not erase the danger-zone
   caveat below.
4. If the source category is aromatic and Biot-Savart or Haigh-Mallion is active,
   add nothing. The aromatic fields remain zero.
5. Compute `D(r_is)`.
6. Build `Qhat_s` from the category source frame.
7. Compute the raw pair response:

```text
A_is = D(r_is) Qhat_s
```

8. Accumulate the fixed-source channel:

```text
accum[c, fixed] += A_is
```

9. Accumulate the MOPAC bond-order channel:

```text
bo_s = MOPAC Wiberg bond order for source s
accum[c, bo] += bo_s A_is
```

`bo_s` is dimensionless. Missing MOPAC bond order means `bo_s = 0` for the
bond-order channel only; it must not remove the fixed-source contribution. If a
noise floor is used, it zeros or clamps the bond-order channel and records the
threshold in metadata.

The fitted contribution used by a model or regression is then:

```text
sigma_hat_i,c = beta_c,0 accum[c, fixed] + beta_c,1 accum[c, bo]
```

`beta_c,0` and `beta_c,1` absorb the susceptibility magnitude, global sign,
unit conversion, and source-strength calibration. This is the defensible MOPAC
claim: Wiberg bond order is a QM-derived, rotationally invariant source
modulator, not a first-principles susceptibility. OpenMOPAC defines its printed
bond order as rotationally invariant, equivalent to Wiberg bond order, and
roughly `1/2/3` for ethane/ethylene/acetylene; `ALLBONDS` lowers the printed
threshold and includes hydrogen-involving bonds. [OpenMOPAC BONDS](https://openmopac.net/Manual/bonds.html),
[OpenMOPAC ALLBONDS](https://openmopac.net/Manual/allbonds.html)

## Irrep Decomposition And Parity

After each category/channel accumulator is complete, decompose the accumulated
`3 x 3` matrix with the existing `SphericalTensor::Decompose` convention:

```text
T0 = trace(A) / 3
T1 = antisymmetric pseudovector:
     [ (A_yz - A_zy)/2,
       (A_zx - A_xz)/2,
       (A_xy - A_yx)/2 ]
T2 = symmetric traceless part in the real spherical basis:
     [ sqrt(2) S_xy,
       sqrt(2) S_yz,
       sqrt(3/2) S_zz,
       sqrt(2) S_xz,
       (S_xx - S_yy)/sqrt(2) ]
```

The emitted irreps are always:

```text
0e + 1e + 2e
```

Never emit `1o` for this tensor. A rank-2 shielding/susceptibility response is
even under inversion. In e3nn language, parity multiplies in tensor products:
two polar-vector factors are `1o x 1o`, which produces even branches. Ben
Mahmoud et al. state the NMR magnetic-shielding tensor decomposition as
`T(0) + T(1) + T(2)` with dimensions `1, 3, 5`, and explicitly state that all
three components have even parity. Geiger and Smidt give the e3nn tensor-product
path rule `p1 p2 = p3` and the two-vector decomposition into scalar,
pseudovector, and five-component even `l=2`. [Ben Mahmoud et al. 2024](../../../references-text/ben-mahmoud-2024-gnn-solid-state-nmr-spherical-tensors-text-1.txt),
[Geiger and Smidt 2022](../../../references-text/geiger-smidt-2022-e3nn-euclidean-neural-networks-text-4.txt)

The `0e` channel is a real scalar coupling. For a traceless source tensor:

```text
pcs0_raw = n^T Qhat_s n / R^3
         = trace(D(r_is) Qhat_s) / 3
```

up to the chosen global sign convention. This is the PCS/McConnell scalar. It is
not the trace of a pure `2e` tensor; a pure `2e` object is traceless by
definition. It is the scalar branch of the coupling between the dipole
propagator and the susceptibility anisotropy. Suturina and Kuprov give this
same scalar PCS branch as the isotropic average and rewrite it with second-rank
spherical harmonics. [Suturina and Kuprov 2016](../../../references-text/suturina-2016-pseudocontact-shifts-from-mobile-spin-labels-text-2.txt)

## Fields Produced

Let:

```text
N = number of target atoms
C = 7 source categories:
    peptide_co, peptide_cn, backbone_other, sidechain_co,
    sidechain_other, disulfide, aromatic_zeroed
```

The chart rows are the logical fields below. Each has one fixed-source channel
and one Wiberg-bond-order-weighted channel.

```text
mc_<cat>_0e_fixed       shape (N, 1)   units Angstrom^-3   parity even
mc_<cat>_1e_fixed       shape (N, 3)   units Angstrom^-3   parity even
mc_<cat>_2e_fixed       shape (N, 5)   units Angstrom^-3   parity even

mc_<cat>_0e_bo          shape (N, 1)   units Angstrom^-3   parity even
mc_<cat>_1e_bo          shape (N, 3)   units Angstrom^-3   parity even
mc_<cat>_2e_bo          shape (N, 5)   units Angstrom^-3   parity even
```

Expanded, the chart row names are:

```text
mc_peptide_co_0e_fixed
mc_peptide_co_1e_fixed
mc_peptide_co_2e_fixed
mc_peptide_co_0e_bo
mc_peptide_co_1e_bo
mc_peptide_co_2e_bo

mc_peptide_cn_0e_fixed
mc_peptide_cn_1e_fixed
mc_peptide_cn_2e_fixed
mc_peptide_cn_0e_bo
mc_peptide_cn_1e_bo
mc_peptide_cn_2e_bo

mc_backbone_other_0e_fixed
mc_backbone_other_1e_fixed
mc_backbone_other_2e_fixed
mc_backbone_other_0e_bo
mc_backbone_other_1e_bo
mc_backbone_other_2e_bo

mc_sidechain_co_0e_fixed
mc_sidechain_co_1e_fixed
mc_sidechain_co_2e_fixed
mc_sidechain_co_0e_bo
mc_sidechain_co_1e_bo
mc_sidechain_co_2e_bo

mc_sidechain_other_0e_fixed
mc_sidechain_other_1e_fixed
mc_sidechain_other_2e_fixed
mc_sidechain_other_0e_bo
mc_sidechain_other_1e_bo
mc_sidechain_other_2e_bo

mc_disulfide_0e_fixed
mc_disulfide_1e_fixed
mc_disulfide_2e_fixed
mc_disulfide_0e_bo
mc_disulfide_1e_bo
mc_disulfide_2e_bo

mc_aromatic_zeroed_0e_fixed
mc_aromatic_zeroed_1e_fixed
mc_aromatic_zeroed_2e_fixed
mc_aromatic_zeroed_0e_bo
mc_aromatic_zeroed_1e_bo
mc_aromatic_zeroed_2e_bo
```

`mc_<cat>_0e_*` is the PCS scalar coupling in the current `SphericalTensor::T0`
normalization. `mc_<cat>_1e_*` is the even pseudovector from the antisymmetric
part of `D Q`. `mc_<cat>_2e_*` is the symmetric traceless real-spherical
five-vector.

The preferred storage is one packed array per category and channel, because it
matches the existing `SphericalTensor::PackFull9` convention:

```text
mc_<cat>_fixed.npy      shape (N, 9)
mc_<cat>_bo.npy         shape (N, 9)
```

with component order:

```text
T0,
T1_x, T1_y, T1_z,
T2_m-2, T2_m-1, T2_m0, T2_m+1, T2_m+2
```

The schema must also carry:

```text
irrep_layout = "0e,1e_x,1e_y,1e_z,2e_m-2,2e_m-1,2e_m0,2e_m+1,2e_m+2"
units        = "Angstrom^-3"
source_model = "unit susceptibility shape; scale pinned"
aromatic_zeroed_when_ring_active = true
bo_source    = "MOPAC Wiberg bond order, dimensionless"
```

For `aromatic_zeroed`, all six logical irrep/channel fields are present but
zero-valued when ring-current kernels are active. Keeping the zero rows is
useful: it makes the double-counting decision auditable in the feature chart.

## MOPAC Bond Order: Exact Claim

MOPAC bond order enters as a second source-strength basis, not as a law of
susceptibility:

```text
fixed channel:  sum A_is
bo channel:     sum BO_s A_is
model:          beta_0 fixed + beta_1 bo
```

The honest thesis sentence is:

> MOPAC provides a rotationally invariant Wiberg bond order, a semiempirical
> QM measure of electron sharing between atom pairs. We use it as a calibrated
> source-strength modulator beside a fixed source channel. We do not claim that
> susceptibility anisotropy is first-principles linear in Wiberg bond order.

I found no held citation that proves `dchi proportional to Wiberg bond order`
for these protein source categories. Therefore the document must not make that
claim. The support is narrower and solid: OpenMOPAC defines the bond-order
quantity; the two-channel design is our calibration choice.
[OpenMOPAC BONDS](https://openmopac.net/Manual/bonds.html)

## Aromatic Zeroing

Aromatic source categories must be zero while Biot-Savart and Haigh-Mallion are
active. The reason is not software tidiness. It is physics partitioning.
Case describes aromatic and unsaturated susceptibility anisotropy effects as the
same remote induced-current picture often called ring-current contributions. If
the ring-current field is already emitted by the ring kernels, a nonzero
aromatic McConnell susceptibility term gives the model the same aromatic
magnetic response twice. [Case 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3877577/)

Implementation rule:

```text
if category == aromatic and ring_current_active:
    pair contribution = 0
```

Do not let the MOPAC bond-order channel reintroduce aromatic McConnell through a
separate path.

## Problematic Parts

These are not reasons to drop the feature. They are the limits that must be
named before the examiner names them.

**1. Susceptibility magnitudes are not settled constants.** The local repo note
on McConnell `dchi` values explicitly says the current numeric table is
secondary/web-cited and primary verification is still needed. It reports
carbonyl disagreement at roughly a multi-fold scale and blocks hard thesis
claims that rest on those values. Held local sources already show the broader
problem: Worcester used an ester proxy of `8.8 x 10^-6`, while Pauling argued
that the corrected peptide-group value is `-5.36 x 10^-6`; Tjandra et al. also
note relatively large uncertainty in reported aromatic-residue anisotropies.
Therefore the defended design is calibrated basis emission, not a fixed
universal `dchi` table. [TECH_DEBT 2026-06-02](../../../references/incoming/research-reports/TECH_DEBT_mcconnell_dchi_primary_refs_2026-06-02.md),
[Worcester 1978](../../../references-text/worcester-1978-structural-origins-diamagnetic-anisotropy-proteins-text-1.txt),
[Pauling 1979](../../../references-text/pauling-1979-diamagnetic-anisotropy-peptide-group-text-1.txt),
[Tjandra et al. 1996](../../../references-text/tjandra-1996-field-dependence-nh-j-splittings-ubiquitin-text-3.txt)

**2. Rhombicity is often unknown.** PCS theory supports axiality plus
rhombicity as a complete five-component susceptibility description, but most
organic-group tables available to this project do not give a complete rhombic
tensor and local frame for every source category. Missing rhombicity is an
absence of data. It is not permission to synthesize a plausible-looking
in-plane term. [Suturina and Kuprov 2016](../../../references-text/suturina-2016-pseudocontact-shifts-from-mobile-spin-labels-text-2.txt)

**3. The source is not really a point.** The point-dipole model is a far-field
approximation. Suturina and Kuprov show that a distributed paramagnetic source
adds higher multipoles; their anisotropic Gaussian example has corrections
larger than 1 ppm out to about 10 Angstrom, and the point model breaks down at
distances comparable to the source distribution. That result is PCS-specific,
but the mathematical warning is the same for any spatially extended
susceptibility source. [Suturina and Kuprov 2016](../../../references-text/suturina-2016-pseudocontact-shifts-from-mobile-spin-labels-text-4.txt)

**4. The clean cone story is not clean at orbital level.** GIAO/NICS and
orbital-decomposition work shows that conventional pi-electron shielding or
deshielding explanations can be wrong or incomplete. PubMed's abstract for
Baranac-Stojanovic, Koch, and Kleinpeter states that TSNMRS values decomposed
into localized and canonical molecular-orbital contributions show nearby C=C
and aromatic-ring proton shifts should not be explained by the conventionally
accepted pi-electron shielding/deshielding effects. Their held 2012
ChemPhysChem paper gives the same kind of decomposition warning: in borazine
and related rings, sigma/pi and tensor-component balances, not a single
aromaticity or pi-current label, determine the shielding surfaces. [Baranac-Stojanovic et al. 2012, PubMed](https://pubmed.ncbi.nlm.nih.gov/22135110/),
[Baranac-Stojanovic et al. 2012, held corpus](../../../references-text/baranac-stojanovic-2012-density-functional-calculations-anisotropic-effects-text-2.txt)

## The Danger Zone

Below about 3 Angstrom, the point-dipole McConnell picture can have the wrong
sign for nearby protons over a C=C bond. Martin and Brown compare the
McConnell shielding-cone prediction with GIAO-SCF calculations and experiment:
for one close alkene geometry, McConnell predicts shielding while GIAO and
experiment show strong deshielding. They explain the failure directly: the
McConnell model includes only magnetic anisotropy, while near-field shifts also
contain electric-field effects, orbital interactions, and dispersion. [Martin and Brown 2000](../../../references-text/martin-2000-graphical-model-proton-deshielding-text-1.txt),
[Martin and Brown 2000, discussion](../../../references-text/martin-2000-graphical-model-proton-deshielding-text-2.txt)

Protein contacts live in this regime. The build response is not to hide the
problem with a cleaner irrep decomposition. The build response is:

```text
emit the point-dipole basis;
record near-field exclusions and distances;
calibrate against DFT/experiment;
state that the model is a far-field/semiempirical descriptor below ~3 Angstrom.
```

The irrep correction makes the feature equivariant and defensible. It does not
make the point-dipole approximation exact.

## Current Code Context

The current `McConnellResult` computes the old mixed tensor

```text
M = 9 cos(theta) d_hat b_hat^T
  - 3 b_hat b_hat^T
  - (3 d_hat d_hat^T - I)
```

and stores `M / R^3`, the separate dipolar `K`, and scalar category sums. The
current `MopacMcConnellResult` repeats the same old tensor and multiplies it by
MOPAC bond order. Both paths still accumulate aromatic channels. That is the
implementation being replaced, not the design being defended. [McConnellResult.cpp](../../../src/McConnellResult.cpp),
[MopacMcConnellResult.cpp](../../../src/MopacMcConnellResult.cpp)

The new feature should be one McConnell source calculation with two channels,
not two different physics paths:

```text
channel fixed:  unit source strength
channel bo:     Wiberg-weighted source strength
```

## Implementation Checks To Add

I am adding these checks because they are where implementation mistakes will
look numerically plausible.

1. **Rotation equivariance.** Rotate a conformation by `R`. `0e` must be
   unchanged, `1e` must rotate as an even pseudovector, and `2e` must rotate by
   the `l=2` Wigner matrix in the existing real-spherical convention.
2. **Parity declaration.** Metadata and model input specs must say `1e`, not
   `1o`.
3. **PCS scalar check.** For each pair, verify `T0(D Q) = trace(D Q)/3` and
   compare with `n^T Q n / R^3` under the chosen sign convention.
4. **Aromatic zero check.** With ring current active, every
   `mc_aromatic_zeroed_*` field is exactly zero in both fixed and BO channels.
5. **MOPAC separation check.** Missing or below-threshold bond order zeros only
   the BO channel. It must not change the fixed channel.
6. **Near-field audit.** Count accepted and rejected source-target pairs below
   3 Angstrom and report them. Do not silently treat the absence of a pair as
   proof that the danger zone does not exist.

## Citation Notes

Load-bearing citations used here:

- **Remote-group McConnell tensor and structural-biology reason:** Case,
  "Chemical shifts in biomolecules." It gives the induced-current picture,
  susceptibility tensor, remote-group McConnell equation, isotropic
  susceptibility cancellation, and pseudo-contact analogy.
  [PMC3877577](https://pmc.ncbi.nlm.nih.gov/articles/PMC3877577/)
- **Point-dipole susceptibility tensor, PCS scalar, spherical-harmonic coupling,
  axiality/rhombicity, distributed-source limit:** Suturina and Kuprov 2016.
  [text-2](../../../references-text/suturina-2016-pseudocontact-shifts-from-mobile-spin-labels-text-2.txt),
  [text-4](../../../references-text/suturina-2016-pseudocontact-shifts-from-mobile-spin-labels-text-4.txt)
- **NMR tensor irreps and even parity:** Ben Mahmoud et al. 2024.
  [text-1](../../../references-text/ben-mahmoud-2024-gnn-solid-state-nmr-spherical-tensors-text-1.txt)
- **e3nn tensor-product parity rule:** Geiger and Smidt 2022.
  [text-4](../../../references-text/geiger-smidt-2022-e3nn-euclidean-neural-networks-text-4.txt)
- **MOPAC Wiberg bond order:** OpenMOPAC `BONDS` and `ALLBONDS` manual pages.
  [BONDS](https://openmopac.net/Manual/bonds.html),
  [ALLBONDS](https://openmopac.net/Manual/allbonds.html)
- **Protein peptide-group susceptibility context:** Worcester 1978, Pauling
  1979, and Tjandra et al. 1996.
  [Worcester](../../../references-text/worcester-1978-structural-origins-diamagnetic-anisotropy-proteins-text-1.txt),
  [Pauling](../../../references-text/pauling-1979-diamagnetic-anisotropy-peptide-group-text-1.txt),
  [Tjandra](../../../references-text/tjandra-1996-field-dependence-nh-j-splittings-ubiquitin-text-3.txt)
- **Wrong-sign near-field warning:** Martin and Brown 2000.
  [summary](../../../references-meta/martin-2000-graphical-model-proton-deshielding-summary.txt),
  [text-1](../../../references-text/martin-2000-graphical-model-proton-deshielding-text-1.txt),
  [text-2](../../../references-text/martin-2000-graphical-model-proton-deshielding-text-2.txt)
- **Orbital-decomposition warning for pi-anisotropy explanations:** PubMed
  abstract for Baranac-Stojanovic, Koch, and Kleinpeter 2012, plus the held
  Baranac-Stojanovic/Koch/Kleinpeter 2012 ChemPhysChem decomposition paper.
  [PubMed 22135110](https://pubmed.ncbi.nlm.nih.gov/22135110/),
  [held text-2](../../../references-text/baranac-stojanovic-2012-density-functional-calculations-anisotropic-effects-text-2.txt)

Known citation gap:

The exact organic-group `dchi` values for the intended category table are not
primary-verified in the held corpus. The project tech-debt note says this
explicitly. Until those primary papers are fetched and summarized, this document
defends calibrated source-shape channels, not hard-coded numerical
susceptibility constants.
