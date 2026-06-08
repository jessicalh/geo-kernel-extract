Historical -- not current truth; see doc/emerging/kernel_design/ring_spec.md and doc/emerging/CONTROLLING_SPEC.md.
# Biot-Savart + Haigh-Mallion Structured Grounding

This is the build-and-defend grounding for the ring-current pair:
`BiotSavartResult` and `HaighMallionResult`. It builds on
`bs_hm_grounding_agent1.md`, `ring.md`, the 2026-06-06 controlling spec, and the
current implementation. It does not reopen the settled decisions.

The decision is: keep both calculators. Biot-Savart is the Johnson-Bovey
distributed-source route, integrated over current loops and lifted to a full
shielding tensor by the Boyd-Skrynnikov form. Haigh-Mallion is the independent
geometric/surface route. The two are separate responsibility kernels because
their agreement, after named calibration and geometry stratification, is a
two-path validation of the same ring-current signal. Both emit the full even
rank-2 decomposition:

```text
0e + 1e + 2e
```

never `1o`.

## Why This Kernel Pair Exists

Aromatic rings respond strongly to the spectrometer field. The applied field
induces a circulating π-electron current; that current produces a secondary
magnetic field in the space around the ring; a nearby nucleus feels that
secondary field as a contribution to its shielding. This is a through-space
magnetic mechanism, not a bond term and not an electrostatic field-gradient
term. Gershoni-Poranne and Stanger give the aromaticity-scale statement:
aromatic molecules show an induced ring current and a resulting induced
magnetic field under an external magnetic field. Case's biomolecular
calibration puts the same effect into the structural-biology setting: methane
probes placed near protein and nucleic-acid aromatic rings receive secondary
shifts that are well described by ring-current geometric factors, with
electrostatic terms handled separately. [Gershoni-Poranne and Stanger 2015](../../../references-text/gershoni-poranne-stanger-2015-magnetic-criteria-aromaticity-text-1.txt),
[Case 1995](../../../references-text/case-1995-ring-current-calibration-text-1.txt)

Plainly: proteins and nucleic acids put protons, nitrogens, carbons, and other
nuclei near Phe, Tyr, His, Trp, and nucleic-base rings. Those rings are not just
nearby matter; under `B0` they are magnetic sources. The ring-current pair is
the classical, cheap, auditable way to give the model that source.

The reason this is a pair, not one calculator, is also structural. Biot-Savart
and Haigh-Mallion are different mathematical routes to the same ring-current
quantity. Case calibrated Johnson-Bovey and Haigh-Mallion scalar ring-current
models against DFT secondary shifts and found no significant fit-quality
difference for biomolecular use. Moyna et al. compared Haigh-Mallion,
Johnson-Bovey, and point-dipole methods on protein proton shifts; the
Haigh-Mallion and Johnson-Bovey RMS values are close, with the all-proton table
at 0.220 ppm vs 0.228 ppm. That is the reported two-path value we keep.
[Case 1995](../../../references-text/case-1995-ring-current-calibration-text-2.txt),
[Moyna et al. 1998](../../../references-text/moyna-zauhar-williams-1998-comparison-ring-current-methods-structure-refinement-text-2.txt)

## What We Compute

For each target atom `i` and aromatic ring `r`, compute a rank-2 shielding
response kernel:

```text
G_ab(i, r) = - V_a(i, r) n_b(r)
```

where `n` is the ring normal and `V` is the secondary-field-like vector from the
chosen ring-current route. The minus sign follows the shielding convention
`sigma_ab = -dB_a^sec / dB_{0,b}`. The `b` index carries the ring normal because
the induced current strength depends on the flux of the applied field through
the ring. The result can have a trace, an antisymmetric part, and a symmetric
traceless part, so the emitted object is the full packed tensor:

```text
[T0, T1_x, T1_y, T1_z, T2_m-2, T2_m-1, T2_m0, T2_m+1, T2_m+2]
```

The two calculators differ in how they obtain `V`.

```text
Biot-Savart / Johnson-Bovey:
    V = B_BS
    B_BS = magnetic field from two displaced current loops
    G_BS_ab = -B_BS,a n_b PPM_FACTOR

Haigh-Mallion surface tensor:
    H_ab = integral over ring surface of (3 rho_a rho_b/rho^5 - delta_ab/rho^3) dS
    V = H n
    G_HM_ab = -V_a n_b
```

The current Biot-Savart code already computes the raw geometry with
`current=1.0` nA and does not multiply the tensor sum by `ring.Intensity()`.
That is correct for the new contract if it is named honestly: BS emits a
unit-current geometric kernel. The physical ring-current intensity is a separate
calibration coefficient. [BiotSavartResult.cpp](../../../src/BiotSavartResult.cpp),
[BiotSavartResult.h](../../../src/BiotSavartResult.h)

The current Haigh-Mallion code already stores both the raw surface tensor `H`
and the shielding tensor `G`. Keep both. `hm_H` is the symmetric/traceless raw
surface integral; `hm_G` is the rank-2 shielding response that can carry
`T0/T1/T2`. [HaighMallionResult.cpp](../../../src/HaighMallionResult.cpp),
[HaighMallionResult.h](../../../src/HaighMallionResult.h),
[ConformationAtom.h](../../../src/ConformationAtom.h)

## Biot-Savart Calculation

The Biot-Savart path is the Johnson-Bovey distributed-source model in code form.
For every accepted atom-ring pair:

1. Build two polygonal current loops from the ring vertices, displaced by
   `+lobe_offset` and `-lobe_offset` along the ring normal.
2. Put half the unit current on each loop.
3. Integrate the magnetic field from all loop segments with the finite
   wire-segment Biot-Savart expression in SI units.
4. Build the rank-2 response matrix:

```text
G_ab = -n_b B_a PPM_FACTOR
```

5. Decompose the accumulated `G` with `SphericalTensor::Decompose`.

This is a distributed-source calculation, not the far-field point-dipole
shadow. It is still classical and empirical after calibration, but it keeps the
ring's spatial extent where protein contacts actually live.

Boyd and Skrynnikov are the tensor citation. Their paper states the classical
ring-current models give the ring-current shielding contribution and supplies
the off-diagonal/full-tensor construction; their eq. 3 assembles the rank-2
shielding tensor from the scalar Johnson-Bovey component and the off-diagonal
component. In this project, the Boyd-Skrynnikov expression is the tensor lift of
the Johnson-Bovey/Biot-Savart field. It is not a separate calculator.
[Boyd and Skrynnikov 2002](../../../references-text/boyd-skrynnikov-2002-ring-current-chemical-shielding-anisotropy-text-1.txt)

The labels must be fixed. Existing names such as `bs_shielding_contribution`,
`bs_shielding.npy`, the HDF5 `"ppm"` unit label, and GeometryChoice `"nA"` for
intensity read as final shielding. They are not final shielding unless a named
ring-current coefficient has been applied. The correct base-unit meaning is
unit-current geometry, practically `ppm_T_per_nA` or an equivalent explicit
name. The intensity unit is `nA/T`, not bare `nA`. [BsShieldingTimeSeriesTrajectoryResult.h](../../../src/BsShieldingTimeSeriesTrajectoryResult.h),
[BsShieldingTimeSeriesTrajectoryResult.cpp](../../../src/BsShieldingTimeSeriesTrajectoryResult.cpp),
[CalculatorConfig.cpp](../../../src/CalculatorConfig.cpp)

## Haigh-Mallion Calculation

The published Haigh-Mallion scalar formula is a signed-area bond sum. Case writes
the geometric factor as a sum over ring-atom pairs/bonds of signed triangle
areas times inverse-cube atom distances. Moyna et al. write the same scalar
shape as:

```text
G(r) = sum_ij S_ij (1/r_i^3 + 1/r_j^3)
```

Sahakyan and Vendruscolo also present the simplified Haigh-Mallion shielding
form in terms of signed triangle areas `Sij` and the two endpoint distances.
[Case 1995](../../../references-text/case-1995-ring-current-calibration-text-1.txt),
[Moyna et al. 1998](../../../references-text/moyna-zauhar-williams-1998-comparison-ring-current-methods-structure-refinement-text-1.txt),
[Sahakyan and Vendruscolo 2013](../../../references-text/sahakyan-vendruscolo-2013-ring-current-electric-field-contributions-text-1.txt)

Our `HaighMallionResult` is not literally that scalar bond-sum formula. It
computes a surface integral of the dipolar kernel over a fan-triangulated ring
area:

```text
H_ab = integral_S [3 rho_a rho_b/rho^5 - delta_ab/rho^3] dS
V    = H n
G_ab = -V_a n_b
```

That is a strong independent geometric route, but it is a tensor surface
variant. The honest naming question remains:

```text
Benchmark: current hm_G.T0 vs published HM scalar bond-sum
Population: populated biomolecular geometries, split by ring type and region
Pass condition: correlate up to scale/sign; do not require pointwise equality

If pass:
    call it "Haigh-Mallion extended to the full tensor"
If fail:
    call it "HM-style tensor surface kernel"
```

Until that benchmark is done, do not over-claim formula equivalence. The
published HM model is scalar. Our value is a defensible full-tensor extension
candidate, not yet a proved identity.

## Irrep Decomposition And Parity

After each atom-level or atom-type-level accumulation, decompose the `3 x 3`
matrix with the existing `SphericalTensor::Decompose` convention:

```text
T0 = trace(G) / 3

T1 = antisymmetric pseudovector:
     [ (G_yz - G_zy)/2,
       (G_zx - G_xz)/2,
       (G_xy - G_yx)/2 ]

T2 = symmetric traceless part in the real spherical basis:
     [ sqrt(2) S_xy,
       sqrt(2) S_yz,
       sqrt(3/2) S_zz,
       sqrt(2) S_xz,
       (S_xx - S_yy)/sqrt(2) ]
```

The current `Types.cpp` decomposition is the right mathematical split and the
packed order is already explicit. [Types.cpp](../../../src/Types.cpp)

The parity is the bug. A shielding response tensor is even under inversion.
Ben Mahmoud et al. decompose rank-2 NMR tensors into `T(0)`, `T(1)`, and
`T(2)` with dimensions `1`, `3`, and `5`, and state that all three components
have even parity; `T(1)` is a pseudovector, not an odd polar vector. Geiger and
Smidt give the same e3nn tensor-product rule: two vectors decompose into a
scalar, a pseudovector, and a five-component even `l=2`, with parity multiplied
along tensor-product paths. [Ben Mahmoud et al. 2024](../../../references-text/ben-mahmoud-2024-gnn-solid-state-nmr-spherical-tensors-text-1.txt),
[Geiger and Smidt 2022](../../../references-text/geiger-smidt-2022-e3nn-euclidean-neural-networks-text-4.txt)

Therefore the emitted layout is:

```text
0e + 1e + 2e
```

not:

```text
0e + 1o + 2e
```

The current trajectory metadata headers and HDF5 writers for BS, HM, and
ring-susceptibility declare `0e+1o+2e`. That is a metadata/schema bug, part of a
recurring trajectory-header family. It does not mean the numerical decomposition
is wrong. [BsShieldingTimeSeriesTrajectoryResult.h](../../../src/BsShieldingTimeSeriesTrajectoryResult.h),
[HmShieldingTimeSeriesTrajectoryResult.h](../../../src/HmShieldingTimeSeriesTrajectoryResult.h),
[RingSusceptibilityShieldingTimeSeriesTrajectoryResult.h](../../../src/RingSusceptibilityShieldingTimeSeriesTrajectoryResult.h)

## Calibration Layer

The raw fields are not final shielding claims. They are geometry kernels. The
calibration layer supplies the physical scale, sign, and ring-type intensity.

Required named variants:

```text
raw_unit_current_BS
raw_unit_geometric_HM
Case95-JB
Case95-HM
legacy_GP_Pullman
```

The current code defaults look like the older Giessner-Prettre/Pullman lineage
normalized around Phe approximately `-12 nA/T`, with His and Trp-5 much smaller
than the Case refit. Case reports old and new ring-current intensity factors
side by side; the new fitted values increase His and Trp-5 substantially and are
different for Johnson-Bovey and Haigh-Mallion. That is why calibration sets must
be named, not silently mixed. [CalculatorConfig.cpp](../../../src/CalculatorConfig.cpp),
[Case 1995](../../../references-text/case-1995-ring-current-calibration-text-2.txt)

Do not multiply raw BS by `ring.Intensity()` inside the raw tensor sum. The
current unit-current behavior is the right base contract. Apply intensity only
when emitting or consuming a named scaled variant, and record:

```text
calibration_set
ring_type
coefficient_value
coefficient_units
source_reference
```

The 5-member-ring lobe offset is an honest open point. Case's implementation
uses current loops 0.64 Å from the ring plane and radii of 1.39 Å and 1.182 Å
for six- and five-membered rings. Perkins likewise describes the Johnson-Bovey
application with loop separation `2q = 0.128 nm`, i.e. ±0.64 Å, while changing
the loop radius for five- vs six-membered rings. Current defaults use 0.50-0.52
Å for His and Trp-5. Now that Johnson-Bovey 1958 is held by the project, the
primary should be checked once it is chunked or otherwise directly readable in
the repo; during this pass I found Perkins and Case as held text bridges, not a
chunked Johnson-Bovey primary. [Case 1995](../../../references-text/case-1995-ring-current-calibration-text-1.txt),
[Perkins et al. 1977](../../../references-text/perkins-1977-ring-current-johnson-bovey-aromatic-amino-acids-text-1.txt),
[CalculatorConfig.cpp](../../../src/CalculatorConfig.cpp)

## Fields Produced

Let:

```text
N = number of target atoms
R = number of aromatic rings
K = 8 aromatic ring-type slots in the current code
```

The logical per-atom fields are:

```text
bs_raw_unit_current        shape (N, 9)   units explicit unit-current kernel   parity 0e+1e+2e
hm_H_raw_surface           shape (N, 9)   units Angstrom^-1 surface kernel     parity 0e+1e+2e, with T0/T1 near zero by construction
hm_G_raw_surface           shape (N, 9)   units Angstrom^-1 shielding kernel   parity 0e+1e+2e

bs_Case95_JB               shape (N, 9)   calibrated by Case95-JB              parity 0e+1e+2e
hm_G_Case95_HM             shape (N, 9)   calibrated by Case95-HM              parity 0e+1e+2e
bs_legacy_GP_Pullman       shape (N, 9)   legacy reproducibility               parity 0e+1e+2e
```

The storage order is the existing packed order:

```text
T0,
T1_x, T1_y, T1_z,
T2_m-2, T2_m-1, T2_m0, T2_m+1, T2_m+2
```

Per atom and ring type, emit the full tensor for every type:

```text
bs_raw_unit_current_by_type      shape (N, K, 9)
hm_G_raw_surface_by_type         shape (N, K, 9)
bs_Case95_JB_by_type             shape (N, K, 9)
hm_G_Case95_HM_by_type           shape (N, K, 9)
```

Flatten to `(N, 72)` if the file format needs a two-dimensional NPY, but the
metadata must state that this is `8 x 9`. The current per-type emit drops `T1`
and writes only `T0` plus `T2`; the revamped contract must not do that.
[BiotSavartResult.cpp](../../../src/BiotSavartResult.cpp),
[HaighMallionResult.cpp](../../../src/HaighMallionResult.cpp),
[ConformationAtom.h](../../../src/ConformationAtom.h)

The per atom-ring sparse table is the most useful audit surface. Required row
fields:

```text
atom_index
ring_index
ring_type
distance_to_center
rho
signed_z
theta
azimuth_cos
azimuth_sin
(3 cos^2 theta - 1) / r^3
lobe_offset
raw_intensity_or_unit_current
calibration_set_id
calibration_coefficient

bs_raw_full9
bs_scaled_full9
hm_H_raw_full9
hm_G_raw_full9
hm_G_scaled_full9

agreement_scaled_T0_residual
agreement_full9_cosine
agreement_T2_cosine
agreement_norm_ratio
agreement_sign_match
filter_flags
```

The existing `ring_contributions.npy` already carries atom/ring ids, basic
ring-frame geometry, BS full 9, HM `H` full 9, HM `G` full 9, and azimuthal
values. The new contract extends that sparse row with calibration identity,
lobe offset, scaled variants, agreement diagnostics, and explicit filter flags.
[ConformationResult.cpp](../../../src/ConformationResult.cpp)

Required schema metadata:

```text
irrep_layout = "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
parity       = "0e+1e+2e"
normalization = "isometric_real_sph"
raw_contract_bs = "unit-current Johnson-Bovey/Biot-Savart geometric kernel"
raw_contract_hm = "unit-geometric Haigh-Mallion surface tensor kernel"
calibration_set = one of raw, Case95-JB, Case95-HM, legacy_GP_Pullman
```

## Two-Path Agreement

Do not compare raw BS and raw HM as equal numbers. They have different raw
units and different built-in scale conventions. Compare only after applying a
named calibration set.

Recommended diagnostics:

```text
scaled T0 R^2
full-9 cosine similarity
T2 cosine similarity
norm ratio
sign agreement
signed residual
absolute residual
residual bins by rho, signed_z, theta, ring_type, distance shell
```

Stratify before summarizing. At minimum:

```text
ring_type
rho bin
signed_z bin
theta bin
near-vertex / not-near-vertex
in-plane / out-of-plane
above-pi-cloud / outside-rim
```

What agreement buys:

1. A sign sanity check. Flipped normals and sign convention errors show up fast.
2. A topology check. Bad ring walks, wrong atom ordering, and fused-ring mistakes
   change the lobe pattern.
3. A unit/calibration check. Accidentally applying intensity twice, or not at
   all, changes norm ratios.
4. A tensor-orientation check. A bogus `T2` lift can match scalar `T0` while
   failing the `T2` cosine.
5. A geometry-danger check. Residual bins identify where the classical surfaces
   disagree and where the law study should not pretend one scalar residual
   tells the whole story.

Legitimate divergence is not failure. Perkins preferred Johnson-Bovey over
Haigh-Mallion for protons above and below aromatic amino-acid rings because of
the region near the π cloud. Moyna reports the biggest method discrepancies for
protons close to aromatic rings and for directly above/below macrocyclic rings.
Agarwal's paracyclophane work is the local-anisotropy warning: after correcting
for local anisotropic contributions, the residual ring-current term correlated
well with both Johnson-Bovey and Haigh-Mallion; the first-stage grounding records
the size of that local anisotropy as about 35% of the observed shift in the
example. [Perkins et al. 1977](../../../references-text/perkins-1977-ring-current-johnson-bovey-aromatic-amino-acids-text-1.txt),
[Moyna et al. 1998](../../../references-text/moyna-zauhar-williams-1998-comparison-ring-current-methods-structure-refinement-text-1.txt),
[Agarwal et al. 1977](../../../references-text/agarwal-1977-ring-currents-local-anisotropy-paracyclophane-text-3.txt)

## Cross-Kernel Boundaries

Ring-current and EFG do not overlap mechanistically. Ring current is magnetic:
an induced current creates a secondary magnetic field. EFG is electric: static
charges create an electric field and field gradient. Do not zero either one
because the other exists.

They do correlate as geometric descriptors near aromatic groups. Ring `2e` and
the aromatic part of EFG `2e` are both low-order angular descriptions of the
same local neighbourhood. That is collinearity to measure, not double-counting
to erase. Sahakyan and Vendruscolo are a useful held bridge here because their
RNA-base analysis fits ring-current and electric-field effects separately and
studies their interdependence. [Sahakyan and Vendruscolo 2013](../../../references-text/sahakyan-vendruscolo-2013-ring-current-electric-field-contributions-text-1.txt)

Ring-current and aromatic McConnell do overlap in mechanism. Both are magnetic
aromatic induced-current/susceptibility pictures. Therefore the McConnell
kernel keeps an `aromatic_zeroed` category while BS/HM are active, and the ring
pair keeps its full aromatic emit. The zeroing happens in McConnell, not here.
[McConnell structured grounding](mcconnell_structured_grounding.md)

`RingSusceptibilityResult` is a useful point-dipole-family cross-check, not one
of the two counted ring-current responsibility kernels. It should not be used to
collapse the BS/HM pair or to replace the distributed-source checks.
[RingSusceptibilityResult.h](../../../src/RingSusceptibilityResult.h)

## What To Feed The Three New Steps

Part 1, the law/correlation study:

```text
raw BS full9
raw HM H full9
raw HM G full9
Case95-JB BS full9
Case95-HM HM G full9
legacy scaled variants
per-type full9 tensors
T0-only scalar views
|T2| and T2 orientation/cosines
ring-frame geometry: distance, rho, z, theta, azimuth
lobe offset and ring type
BS-HM agreement diagnostics
filter flags and near-ring strata
```

Run the scalar laws and tensor laws separately. Report whether `T2` adds signal
over `T0` only. Use BS-HM agreement as a validation and stratification variable,
not as a raw feature with no interpretation.

Part 2, the shielding-tensor predictor:

```text
bs: 1x0e + 1x1e + 1x2e
hm_G: 1x0e + 1x1e + 1x2e
hm_H: mainly 1x2e, emitted in the same packed 9 for audit consistency
```

All are even parity. If a downstream target deliberately wants the symmetric
observable only, expose a secondary `0e+2e` view. Do not replace the primary
full emit with that view.

Part 3, the experiment-facing shift predictor:

Feed the rich ablatable record: raw tensors, scaled tensors, per-type tensors,
sparse pair rows, counts/shells, calibration identities, filter flags,
near-ring indicators, and BS-HM residual/coherence metrics. Experimental shifts
mix ring current with electric field, local anisotropy, hydrogen bonding,
solvent, conformational averaging, referencing, and assignment noise. The
feature record should preserve that ambiguity so ablations can separate it.

## Problematic Parts

These are not reasons to drop the pair. They are the limits that make the pair
defensible instead of over-sold.

**1. Calibration constants are empirical and model-specific.** Case's fitted
intensity factors differ from older Giessner-Prettre/Pullman values, and the
Johnson-Bovey and Haigh-Mallion fitted factors are not identical. His and Trp-5
are the obvious current-code pressure points. The build response is named
calibration sets and raw geometry preservation, not one silent constant table.
[Case 1995](../../../references-text/case-1995-ring-current-calibration-text-2.txt)

**2. The HM tensor equivalence is not proven.** The published HM object is a
scalar signed-area bond sum. Our code computes a surface-integral tensor. The
right next step is the `hm_G.T0` vs published-bond-sum benchmark in populated
geometries. Until then, the name stays conditional.

**3. The 5-member lobe offset needs primary checking.** Case and Perkins bridge
to ±0.64 Å loop displacement, while current His and Trp-5 defaults are 0.50-0.52
Å. Johnson-Bovey 1958 is now held by the project but was not found as a chunked
text file in this pass; check the primary before finalizing this default.

**4. Near-ring regions are real danger zones.** The point-dipole shadow
`(3 cos^2 theta - 1)/r^3` is a far-field descriptor, not the production ring
kernel. Even distributed classical models are stressed above/below the π cloud,
near vertices, in-plane near the ring, and in hetero/fused rings. The response
is to emit filters, distances, sparse rows, and stratified residuals.

**5. Current labels blur raw and calibrated quantities.** `bs_shielding.npy`,
`hm_shielding.npy`, and the trajectory headers must not make raw kernels look
like final ppm shieldings. This is a doc/schema fix, not a physics rewrite.

## The Danger Zone

Treat these regions as audit strata:

```text
distance below the near-field threshold
inside or close to the ring radius
directly above/below the pi cloud
near a ring vertex or heteroatom
in-plane near the rim
fused-ring boundary
asymmetric hetero-ring azimuth
large BS-HM norm ratio
large T2 disagreement with matched T0 sign
```

For accepted pairs in these regions, still emit the row, but mark the condition.
For rejected pairs, emit counts and filter reasons. A silent missing pair is not
defensible in a thesis or in a debugging session.

## Current Code Context

Current strengths:

```text
BiotSavartResult computes a Johnson-Bovey-style two-loop field.
BiotSavartResult already computes raw unit-current geometry.
HaighMallionResult stores both hm_H and hm_G.
Both calculators decompose full tensors with SphericalTensor::Decompose.
ring_contributions.npy already stores sparse per atom-ring full9 BS/HM values.
```

Current fixes to make in the revamp:

```text
rename raw BS/HM fields so they do not imply final ppm shielding
correct trajectory parity metadata to 0e+1e+2e
emit per-type full 8x9, including T1
add named scaled variants instead of one ambiguous "shielding" value
add calibration identifiers and agreement diagnostics to sparse rows
benchmark hm_G.T0 against the published HM scalar bond-sum
refresh/default calibration toward Case95 while preserving legacy GP/Pullman
check 5-member lobe offsets against Johnson-Bovey primary once ingested
```

Do not claim the revamped emit will fit DFT or experiment better. That is the
law study's empirical question. The defended claim here is the form, the
contract, and the corroborated reason to keep both paths.

## Implementation Checks To Add

1. **Rotation equivariance.** Rotate a conformation by `R`. `T0` must be
   invariant, `T1` must rotate as an even pseudovector, and `T2` must rotate by
   the existing real-spherical `l=2` convention.
2. **Parity metadata.** Every BS/HM full tensor payload says `0e+1e+2e`, never
   `0e+1o+2e`.
3. **Raw-vs-scaled separation.** Raw BS with `current=1.0` is unchanged by
   ring intensity. Scaled variants change only through the named calibration
   layer.
4. **BS sign sanity.** A ring-normal flip flips the expected tensor sign pattern
   in a controlled way; applying intensity twice fails the norm-ratio check.
5. **HM scalar benchmark.** Compare `hm_G.T0` with the published HM signed-area
   bond sum on a fixed geometry set, reporting scale, sign, correlation, and
   residual strata.
6. **Per-type completeness.** Per-type arrays include all 9 components for all 8
   ring-type slots; `T1` is not silently dropped.
7. **Sparse-row audit.** Every accepted atom-ring pair records geometry,
   calibration identity, BS full9, HM `H` full9, HM `G` full9, and agreement
   diagnostics. Every rejected pair contributes to filter counts.
8. **Near-ring stress tests.** Include above-plane, in-plane, near-vertex,
   hetero-ring, and fused-ring fixtures. Do not use only far-field benzene.
9. **Calibration reproducibility.** `Case95-JB`, `Case95-HM`, and
   `legacy_GP_Pullman` can be regenerated from explicit coefficient tables and
   produce distinct metadata.
10. **Pack/unpack round trip.** `PackFull9` and `Reconstruct` preserve `T0/T1/T2`
    within floating tolerance for both symmetric and rank-1 asymmetric examples.

## Cites Used

- **Ring-current mechanism and magnetic context:** Gershoni-Poranne and Stanger
  2015; Plasser and Glöcklhofer 2021. [Gershoni-Poranne and Stanger](../../../references-text/gershoni-poranne-stanger-2015-magnetic-criteria-aromaticity-text-1.txt),
  [Plasser and Glöcklhofer](../../../references-text/plasser-2021-vist-shielding-tensors-aromaticity-text-1.txt)
- **Biomolecular scalar calibration and no significant JB-vs-HM fit
  difference:** Case 1995. [text-1](../../../references-text/case-1995-ring-current-calibration-text-1.txt),
  [text-2](../../../references-text/case-1995-ring-current-calibration-text-2.txt)
- **Protein method comparison and near-ring discrepancy warning:** Moyna et al.
  1998. [text-1](../../../references-text/moyna-zauhar-williams-1998-comparison-ring-current-methods-structure-refinement-text-1.txt),
  [text-2](../../../references-text/moyna-zauhar-williams-1998-comparison-ring-current-methods-structure-refinement-text-2.txt)
- **Johnson-Bovey amino-acid bridge and 0.64 Å loop displacement bridge:**
  Perkins et al. 1977. [text](../../../references-text/perkins-1977-ring-current-johnson-bovey-aromatic-amino-acids-text-1.txt)
- **Full shielding-tensor ring-current lift:** Boyd and Skrynnikov 2002.
  [text](../../../references-text/boyd-skrynnikov-2002-ring-current-chemical-shielding-anisotropy-text-1.txt)
- **Local anisotropy warning and JB/HM compatibility after correction:**
  Agarwal et al. 1977. [text-2](../../../references-text/agarwal-1977-ring-currents-local-anisotropy-paracyclophane-text-2.txt),
  [text-3](../../../references-text/agarwal-1977-ring-currents-local-anisotropy-paracyclophane-text-3.txt)
- **HM scalar formula and ring/electric-field separation:** Sahakyan and
  Vendruscolo 2013. [text-1](../../../references-text/sahakyan-vendruscolo-2013-ring-current-electric-field-contributions-text-1.txt),
  [text-2](../../../references-text/sahakyan-vendruscolo-2013-ring-current-electric-field-contributions-text-2.txt)
- **NMR tensor irreps and even parity:** Ben Mahmoud et al. 2024.
  [text](../../../references-text/ben-mahmoud-2024-gnn-solid-state-nmr-spherical-tensors-text-1.txt)
- **e3nn tensor-product parity rule:** Geiger and Smidt 2022.
  [text](../../../references-text/geiger-smidt-2022-e3nn-euclidean-neural-networks-text-4.txt)

Known citation gaps:

```text
Johnson-Bovey 1958: held by the project, but not found as chunked text in this pass.
Haigh-Mallion 1971/1980 originals: not held here; Case, Perkins, Moyna, and
Sahakyan are the held bridges for the scalar formula and use in biomolecular
calibration.
```

## Output Boundary

This pass writes exactly one new file:
`doc/emerging/kernel_design/bs_hm_structured_grounding.md`. It does not edit
source files, other documents, generated files, commits, or build outputs.
