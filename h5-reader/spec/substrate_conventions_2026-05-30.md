# Substrate conventions (Phase 1 gate — precedes all substrate code)

Date: 2026-05-30
Status: revision 2 (2026-05-30 Opus validation pass) — settled per
the substrate-tiers audit (Codex + Claude general-purpose) and the
follow-up Opus authoritative pass that validated each convention
against actual project source (`src/`, `python/nmr_extract/`,
`h5-reader/src/model/`).

This document records every convention decision the substrate makes.
Each section names: the question, the options considered, the chosen
convention, the reasoning, and what breaks if the convention is later
changed.

## Revision history

- **2026-05-30 (revision 2, Opus validation pass)**: corrected SH
  basis ordering anchor from `h5-reader/src/model/SphericalDecomposition.cpp`
  (which is acknowledged in its own header to NOT match the library's
  order) to `src/Types.cpp::SphericalTensor::Decompose`; added P₂ vs
  (3cos²θ−1) asymmetry note for iRED/TCF; added prochiral R/S enum
  alongside diIndex numeric for HA frame; added TRP bridgehead
  constituent-rings (NOT perimeter ring) policy; removed `ring_pucker_6`
  (YAGNI for protein NMR scope); changed `ring_pucker_5` return to
  (Q, θ) tuple; added explicit HBondDetection parameter for h-bond
  accessors (both Larsen-geometric and DSSP-Kabsch-Sander coexist per
  `feedback_methods_accumulate`); added CoulombResult-divergence
  acknowledgement; added Markley pseudoatom name exposure; added PBC
  verbatim-port rule; added substrate-version-in-Parquet-filename rule;
  added equivalence-class-typed-enum scope clarification.
- **2026-05-30 (revision 1)**: initial draft per the substrate-tiers
  audit (Codex + Claude general-purpose).

Companions:
- `reader_as_platform_2026-05-29.md` — the platform architecture that
  uses these conventions
- `stage2_pysr_campaign_2026-05-29.md` — the campaign that consumes
  the substrate

## Hard rules across all conventions

- **All right-handed coordinate frames.**
- **All angles in radians internally**; degrees only at user-facing
  display layer. Both forms documented per accessor where users might
  see angles.
- **All distances in Ångström.** Internal SI consistency check:
  multipoles in `e·Å`, `e·Å²`; fields in `e/Å²`; field gradients in
  `e/Å³`.
- **Charges in elementary charge units `e`**. Proton = +1e, electron
  = −1e (standard SI).
- **No default cutoffs anywhere.** Every spatial cutoff parameter is
  required at the call site. Bypassing this rule is the single largest
  silent-misleading risk both audit agents flagged. `target_class` and
  `cutoff_Å` are both required parameters; defaults are forbidden.
- **All tensors as `Mat3` (lab/local frame) plus `SphericalTensor`
  decomposition** when a rank-2 form is the natural output. Never
  collapse to scalar in the substrate; downstream transforms decide
  whether to project.
- **PBC awareness**: every distance/angle/vector accessor that crosses
  potential PBC boundaries uses minimum-image. If PBC absent or
  unbounded, accessor returns Euclidean unchanged. Substrate detects
  PBC at load from the trajectory metadata.

## Spherical harmonic convention

**Question**: which spherical harmonic basis to use for neighbor-direction
features.

**Choice**: Match the LIBRARY's `src/Types.cpp::SphericalTensor::Decompose`
basis (NOT the h5-reader's `SphericalDecomposition.cpp` — they differ;
see "Pre-existing inconsistency" below). The 5 real l=2 components,
ordered m = −2..+2, in isometric real-SH normalization:

```
T2[0] = √2     · Sxy           // m = −2
T2[1] = √2     · Syz           // m = −1
T2[2] = √(3/2) · Szz           // m =  0
T2[3] = √2     · Sxz           // m = +1
T2[4] = (1/√2)(Sxx − Syy)      // m = +2
```

The isometric normalization preserves Frobenius norm: Σ T2_m² = Σ S_ij²
(matches e3nn-style component normalization, NOT SciPy's raw
`sph_harm` output which carries the √(15/4π) prefactor).

**Why the library's order, not the h5-reader's**: T2 residuals come
from the library's H5 in the library's component order; Y_l^m emitted
by the substrate must match. Anchoring on h5-reader's `SphericalDecomposition.cpp`
would produce a basis-mismatch that PySR would discover as a spurious
coefficient — exactly the silent-misleading failure mode this
conventions doc exists to prevent. **Corrected 2026-05-30 after Opus
validation pass identified the error in the original draft.**

**Pre-existing inconsistency (do not fix during the convention pass)**:
the library (`src/Types.cpp:46-57`) and the h5-reader
(`h5-reader/src/model/SphericalDecomposition.cpp`) use different T2
component orderings. The h5-reader's own header (lines 22-24)
acknowledges the divergence and notes it was acceptable "since kernel
shielding is not compared against DFT for the trajectory." The
substrate breaks that assumption immediately. Resolution: the
substrate uses the library order. The h5-reader's display ordering
stays as-is (locally consistent within h5-reader, doesn't cross-check
against substrate Y_l^m). Reconciling the two is a future fileformat
schema bump, out of scope here.

**Excluded**:
- Complex harmonics (PySR cannot reason about complex numbers).
- Solid-state-chemistry real harmonic convention with different sign on
  m < 0.
- Raw SciPy normalization (the `√(15/4π)` factor) — substrate emits
  isometric normalization to match `Types.cpp`. Transforms wanting
  SciPy-comparable values multiply by the prefactor explicitly.
- `Y₀⁰` (constant by construction; emitting it wastes PySR search
  budget on rediscovering `Y₀⁰ × constant = constant`).

**What breaks if changed**: every transform that consumes `Y_l^m`
columns; every T2 residual fit; every cross-tool comparison (AI Feynman
+ PySR will disagree on sign of m<0 components).

**Test**: substrate has a fixture computing `Decompose` on a known Mat3
(e.g., the dipolar `(3ẑẑᵀ − I)`) and asserting it matches `Types.cpp`'s
output byte-for-byte. Separately, a fixture asserting that for a unit
vector along ẑ, the Y_l^m output (isometric normalization) returns
T2[2] = √(3/2) and T2[0,1,3,4] = 0. Fails if convention drifts from
`Types.cpp`.

## Sign convention for the dipolar angular form

**Question**: `(3cos²θ − 1)` vs `(1 − 3cos²θ)` vs `P₂(cos θ) = ½(3cos²θ − 1)`.

**Choice**: **`(3cos²θ − 1)`** throughout the substrate. Matches the
existing `GEOMETRIC_KERNEL_CATALOGUE.md` derivation. Does NOT match
the physics `P₂` (which has the ½ factor); transforms that want `P₂`
divide by 2 explicitly.

**Excluded**:
- `P₂` form (half-factor causes discovered coefficients to differ
  by 2 from kernel-catalogue derivations).
- `(1 − 3cos²θ)` (algebraically the negation; mixing this with the
  chosen form gives PySR two ways to express the same physics).

**What breaks if changed**: discovered coefficients flip sign or
factor by 2; comparison with `GEOMETRIC_KERNEL_CATALOGUE.md`
literature constants requires rescaling; cross-transform consistency
breaks.

**Substrate API**: `traj.dipolar_kernel(atom, source_pt, ref_dir,
frame_t)` returns `(3cos²θ − 1)/r³` as a scalar. The `ref_dir` is the
reference direction (ring normal, bond axis, etc.) per the transform's
choice. The function is a method, not a column, so the convention is
enforced in one C++ implementation rather than drifting across N
transform files.

**Asymmetry with iRED / TCF order parameters (added 2026-05-30 after
Opus validation pass)**: dynamics observables like the Lipari-Szabo S²
and the time correlation function C(τ) = ⟨P₂(û(t)·û(t+τ))⟩ are
defined in the literature using `P₂(x) = ½(3x² − 1)`, the physicist's
Legendre form (Lipari-Szabo 1982; Brüschweiler-Wright). The library
implements this correctly in `src/TrajectorySpectral.h:53` (`inline
double P2(double u) { return 0.5 * (3.0 * u * u - 1.0); }`) and uses
it in `IRedOrderParameterTrajectoryResult` and
`ReorientationalDynamicsTrajectoryResult`. **The substrate must NOT
replace `P₂` in iRED/TCF accessors with the dipolar-shielding form**;
keep `P₂` for those specific dynamics observables. The asymmetry is
intentional and matches literature conventions for each quantity:
shielding kernels use `(3cos²θ − 1)`, dynamics observables use `P₂`.

## Local frames per atom class

All local frames right-handed, all returned as `LocalFrame { Vec3 z;
Vec3 x; Vec3 y; bool is_valid; FrameVariant variant; }`. `y` is always
`z × x`. `is_valid` false means the frame couldn't be constructed for
this atom at this frame (edge case below); transforms must check.

### HN frame
- **z**: unit vector from N to H (the NH bond direction)
- **x**: unit vector in the amide plane perpendicular to z. Computed
  as the in-plane component of the (C_{prev} − N) vector, where
  C_{prev} is the previous residue's backbone C. Specifically: project
  `(C_prev − N)` onto the plane normal to z, then normalize.
- **y**: `z × x` (out of amide plane, right-hand rule)

**Edge cases**:
- N-terminus (no `C_prev`): fall back to using `Cα − N` for the in-plane
  reference. `variant = HN_NTerminus`.
- Pre-proline (ω of proline is often ≈ 0): the amide plane is
  geometrically defined regardless of ω; emit normally with
  `variant = HN_PrePro`.
- Cis peptide bond (|ω| < 90°): frame still well-defined; emit with
  `variant = HN_Cis` so transforms can filter.
- Proline (no HN exists): `is_valid = false`.

### HE/HD aromatic frame (PHE, TYR, TRP, HIS)
- **z**: unit ring normal. Direction chosen by right-hand rule from a
  **fixed ring traversal order** (the atoms in the ring's
  `member_atoms` list in `QtRing`, traversed in registration order;
  see `src/Ring.cpp:34-37` for the upstream "same direction every
  frame" guarantee).
- **x**: unit vector from ring centroid to a **chemistry-typed anchor
  atom** (NOT the protonated nitrogen for histidines, since HID and
  HIE protonate different N atoms — anchoring on the protonated N
  rotates the frame 60° between tautomers silently).
  - PHE: anchor on CG
  - TYR: anchor on CG
  - TRP pyrrole ring (HE1): anchor on CG
  - TRP benzene ring (HE3, HZ2, HZ3, HH2): anchor on CD2
  - HID/HIE/HIP imidazole: anchor on CG (regardless of tautomer)
- **y**: `z × x`

**Edge cases**:
- TRP bridgehead atoms (CD2, CE2): atom is in BOTH the pyrrole (TRP5)
  AND the benzene (TRP6) rings. **Emit BOTH constituent ring frames
  (TRP5 + TRP6); do NOT emit a TRP9 perimeter-ring frame on
  bridgeheads.** Per `GEOMETRIC_KERNEL_CATALOGUE.md:275-291`, the
  fused-ring T0 is computed as sum of TRP5+TRP6 contributions, so
  emitting a perimeter frame on the bridgehead atoms would
  double-count.
  `traj.atom_ring_frames(atom)` returns a list (length 2 for TRP
  bridgeheads; length 1 for non-bridgehead aromatic atoms).
  **Clarification added 2026-05-30 after Opus validation pass.**
- Ring atom positions partial in this frame (e.g., HIS with
  alternate protonation): use whatever atoms are present;
  `is_valid = true` provided the ring centroid is computable.

### HA / Cα chirality frame
- **z**: unit vector from Cα to HA (the chirality direction)
- **x**: unit vector along Cα to N
- **y**: `z × x`

**Edge cases**:
- Glycine (HA2 + HA3, prochiral pair): the substrate exposes BOTH
  `QtAtom::diIndex` (IUPAC numeric, 2/3 per Markley 1998 Table 1) AND
  `QtAtom::prochiral` (R/S enum). They are NOT equivalent — for most
  residues the convention is 3=R, 2=S, but **Gly inverts**: HA2=R,
  HA3=S (per `SemanticEnums.h:170-171`). Transforms fitting on
  prochiral identity want R/S; transforms fitting on IUPAC name want
  the numeric. Expose both. For HA2 and HA3: same construction; the
  two frames are mirror images; transforms combining them must
  symmetrize correctly. **Both fields exposed per Opus validation
  pass 2026-05-30; original draft only mentioned diIndex.**
- D-amino acids: frame inverts naturally (z flips); `variant =
  DAminoAcid` flag for transforms.

### Cα frame (for backbone-Cα queries)
- **z**: unit vector along bisector of (Cα−N) and (Cα−C), pointing
  away from the backbone
- **x**: unit vector along (Cα−N)
- **y**: `z × x`

### C=O carbonyl frame
- **z**: unit vector from C to O (the carbonyl bond direction)
- **x**: in-peptide-plane direction, computed as in-plane component
  of (C − N_{next}) projected onto plane normal to z
- **y**: `z × x`

**Critical**: this z-axis IS the McConnell kernel's reference
direction. Make sure McConnell-related transforms use this frame's z,
not redefine. Convention mismatch silently misaligns every
neighbor-in-local-frame component used in a McConnell fit.

### Methodological divergence from the library's CoulombResult (acknowledged, not reconciled)

`src/CoulombResult.cpp:91-97` uses "first bond direction (arbitrary
but consistent)" as the heavy-atom primary axis for Buckingham E-field
projection. This is a different choice from the substrate's typed
per-class frames above. The library author called this out as
arbitrary; the substrate's per-class frames are an improvement. **Do
NOT retro-fit the library** — that's an extractor change which is out
of scope per `feedback_extractor_untouchable`. Instead, the substrate's
per-class frames produce a numerically different `coulomb_E_bond_proj`
from the library's. This is two methodological coordinates, both
valid; transforms should document which they use in their schema.
**Divergence acknowledgement added 2026-05-30 after Opus validation.**

## Charge sources

**Question**: which atomic charges does the substrate expose, and how
are they distinguished.

**Choice**: three separate substrate methods, no fallback chains:
- `traj.charge_ff14sb(atom, frame_t)` — from prmtop (static, frame
  index ignored but kept for signature uniformity)
- `traj.charge_mopac(atom, frame_t)` — per-frame from MOPAC `.mop`
  output (returns `None` when MOPAC absent for this frame)
- `traj.charge_aimnet2(atom, frame_t)` — per-frame from AIMNet2 model
  output

Multipole methods take a `charge_source` enum parameter:
`traj.charge_within(atom, frame_t, cutoff_Å, charge_source)` where
`charge_source ∈ {ChargeSource.FF14SB, ChargeSource.MOPAC,
ChargeSource.AIMNet2}`. No default.

**Excluded**: a single `traj.charge(atom, frame_t)` that picks based
on availability — would hide which charge source is in use; results
across transforms become incomparable.

**What breaks if changed**: charge-source sensitivity study breaks
silently; PySR rediscovers different equations for different runs
that look identical except for charge source.

## Charge multipole conventions

**Question**: how to define the local charge dipole and quadrupole.

**Choice**:
- **Origin**: at the target atom (always `r_atom`)
- **Dipole**: `μ = Σ_{i ∈ cutoff} q_i (r_i − r_atom)`, returns `Vec3`
- **Quadrupole**: `Q_ab = Σ_{i ∈ cutoff} q_i [3(r_i − r_atom)_a
  (r_i − r_atom)_b − δ_ab |r_i − r_atom|²]`, **traceless**, returns
  `Mat3`
- Cutoff is required, no default
- Target atom itself is excluded from the sum
- The atom's own residue can be excluded via a separate parameter
  `exclude_residue: bool` (default `False`); transforms studying
  through-space interactions set `True`

**Excluded**:
- Non-traceless quadrupole convention
- Dipole/quadrupole computed in a fixed lab frame (always atom-local)

**Test**: fixture with a single +1e charge at distance r=1Å along z
returns dipole=(0,0,1) and quadrupole with Q_zz = +2, Q_xx=Q_yy = −1
(traceless: Q_xx + Q_yy + Q_zz = 0).

## Cremer-Pople ring puckering

**Critical**: per `feedback_adversarial_review_physics`, this codebase
had a CRITICAL Cremer-Pople θ inversion bug in `PlanarGeometryResult`
that all smoke tests passed. The fix is in the library at
`src/PlanarGeometryResult.cpp:116-122` (mean-plane normal built from
sin/cos-weighted displacement basis: `n = R'₁ × R'₂` where `R'₁ = Σⱼ
(rⱼ − G) sin(2π j / 5)` and `R'₂ = Σⱼ (rⱼ − G) cos(2π j / 5)`; the
alternative `Σ rⱼ × rⱼ₊₁` edge-cross construction gives a normal
anti-parallel to canonical and was the historical bug). The substrate
must port the corrected formula with its own regression test (the
reader does not link the library per `feedback_extractor_untouchable`
+ the platform doc's reader-untouchable analog).

**Choice (corrected 2026-05-30 after Opus validation pass)**:
- **5-rings only.** `ring_pucker_6` is REMOVED from substrate v0.1.0.
  All 6-rings in standard amino acids are aromatic (PHE, TYR, TRP6,
  HIS) and planar — no puckering DOF. The only 5-ring is PRO
  pyrrolidine. The library has NO `CremerPople6Ring`; the original
  doc proposed a 6-ring API the project will never use. YAGNI.
- 5-rings: emit BOTH amplitude Q (Å) AND phase angle θ (degrees,
  wrapped [0°, 360°)). Q is the substrate's "is this ring puckered
  or essentially flat" answer; θ is meaningful only when Q is
  non-degenerate. Original draft emitted θ alone; **Opus catch:
  discarding Q loses the magnitude that's needed to gate on
  degeneracy.**
- Ring atom traversal: registration order in `QtRing::member_atoms`.
- Method: `traj.ring_pucker_5(ring_id, frame_t) -> (Q_angstrom,
  theta_degrees)`. Returns `theta = NaN` when `Q < 1e-6 Å`
  (degenerate planar; library's existing sub-amplitude guard at
  `PlanarGeometryResult.cpp:158-168`); consumers MUST check Q before
  using θ.
- 6-ring API: not exposed in v0.1.0. Add if sugars / modified
  residues enter scope later.

**Phase convention**: `θ = atan2(−Σ zⱼ sin(4π j/5), Σ zⱼ cos(4π j/5))
wrapped [0°, 360°)`. The negative sign on the sine sum is part of the
Cremer-Pople 1975 phase choice; matches `PlanarGeometryResult.cpp:146-169`.

**Provenance sidecar entry**: `convention_ref:
"src/PlanarGeometryResult.cpp::CremerPople5Ring (Cremer & Pople 1975
JACS 97:1354, sin/cos-weighted mean-plane normal; phase = atan2(−Σ zⱼ
sin(4πj/5), Σ zⱼ cos(4πj/5)) wrapped [0°, 360°))"`.

**Test**: idealized envelope ring (one atom at +δ, the others
coplanar) gives θ at one of the 5 envelope-pucker labels (θ = 18° +
72°·k for canonical labeling). Match against `PlanarGeometryResult`
on the same five coordinates byte-for-byte.

## Residual baselines

**Question**: what does "kernel residual" mean.

**Choice**: never expose a method named `kernel_residual`. Instead:
- `traj.literature_residual(atom, frame_t, kernel_name, lit_constant_id) -> SphericalTensor`
  — subtracts the literature equation's prediction using a named
  literature constant
- `traj.stage1_ridge_residual(atom, frame_t, kernel_name) -> SphericalTensor`
  — subtracts what Stage 1's per-stratum ridge fit predicts (uses
  calibrated coefficients per stratum)

Both return `SphericalTensor` (T2-preserving per
`feedback_t2_sacred`), not scalar.

**Excluded**: a single `kernel_residual` method that's ambiguous
about which baseline. Kernels aren't shielding per
`feedback_kernel_not_shielding`; "residual" requires explicit
baseline naming.

**What breaks if changed**: the residual-structure diagnostic (per
the Stage 2 methodology section, `#8`) becomes ambiguous;
cross-transform residual comparisons become wrong.

## Symmetry equivalence

**Question**: how are NMR-equivalent atom classes determined.

**Choice**: explicit chemistry/topology rules using the existing
`QtAtom::equivalenceClass` field (currently `int8_t`, promote to
typed enum). NOT geometric symmetry detection.

- PHE HE1/HE2 (and HD1/HD2): equivalence class id assigned at
  topology load
- TYR HE1/HE2 (and HD1/HD2): same
- HIS HD2/HE1 (etc., depending on tautomer): explicit per-tautomer rules
- ALA HB1/HB2/HB3 (methyl): equivalent
- VAL HG11/HG12/HG13 + HG21/HG22/HG23: two separate methyl classes
- ARG HH11/HH12/HH21/HH22 + HE: explicit
- LYS HZ1/HZ2/HZ3 (terminal ammonium): equivalent

Methods:
- `traj.symmetry_class(atom) -> ClassId`
- `traj.equivalent_atoms(atom) -> list[Atom]` (returns siblings)
- NO `symmetry_averaged_kernel` method. Averaging is a transform
  decision, not a substrate operation (both audit agents flagged
  this).

**What breaks if changed**: silent introduction of geometry-based
equivalence creates IUPAC-style traps the project has explicitly
rejected (per `project_iupac_revert_2026-04-27` +
`feedback_identity_from_chemistry_not_position`).

### Symmetry averaging support (added 2026-05-30 after Opus reframe pass)

The substrate does NOT provide a `symmetry_averaged_kernel` method
(averaging stays a transform decision). But the substrate exposes the
building blocks + a documented SDK helper for transforms that want
symmetry-averaged targets explicitly:

```python
# Substrate building block
traj.equivalent_atoms(atom) -> list[Atom]

# SDK helper (in stage2_sdk.symmetry module)
from stage2_sdk.symmetry import symmetry_averaged_target_series
y_avg = symmetry_averaged_target_series(traj, seed_atom, kernel_name)
# returns numpy (n_frames,) averaged across the equivalence class
```

**Why this matters for canonical rediscovery**: bootstrap coefficient
CIs (per the Stage 2 methodology section) tighten when prochiral
siblings are averaged because the per-class geometric-but-not-chemical
variance is removed. For PHE HE1/HE2 fits, treating the two as
independent samples inflates the bootstrap CI for the discovered
Pople coefficient. Explicit symmetry averaging via the SDK helper
sharpens the methods-paper claim.

**Discipline preserved**: averaging never becomes the substrate
default. Transforms explicitly opt in via the helper; provenance
sidecar records `symmetry_averaging: enabled` per column when used.

## Prochiral markers

**Question**: how are prochiral pairs (e.g., HB2 vs HB3, HA2 vs HA3
for Gly) distinguished.

**Choice**: typed enum on `QtAtom`: `ProchiralityLabel ∈ {None, ProR,
ProS}`. Exposed as Python attribute `atom.prochirality_label`.
Existing `QtAtom::diIndex` already encodes the distinction; promote
to documented public predicate.

**What breaks if changed**: methylene-related PySR fits average
chirality silently and lose the prochiral signal that
`feedback_residual_as_ml_feature` says should be preserved.

## DSSP usage

**Question**: how is DSSP secondary-structure assignment used.

**Choice**: stratification metadata ONLY. Never as a direct PySR
feature column.

- `atom.residue_dssp(frame_t) -> DsspCode` — typed enum return
- Transforms can stratify on DSSP (helix vs sheet rediscovery
  comparison)
- Transforms must NOT pass DSSP code directly to PySR as a feature

**Reasoning**: DSSP is a derived label, not a physical mechanism. PySR
"discovering" σ depends on DSSP code amounts to memorizing that
helices have different chemistry than sheets; no rediscovery value.

**What breaks if changed**: PySR finds DSSP shortcuts instead of
mechanisms; rediscovery story degrades.

## H-bond detection

**Question**: how are H-bond donor/acceptor relationships determined.

**Choice (corrected 2026-05-30 after Opus validation pass)**: the
producer has TWO H-bond detection methods that BOTH stay per
`feedback_methods_accumulate`:
- `HBondResult.cpp` — DSSP energy-based (Kabsch-Sander criterion via
  `DsspResult`)
- `LarsenHBondShieldingResult.cpp` + `LarsenHBondGrid` — geometric
  (spatial cutoff + classification per acceptor type via
  `ClassifyAcceptor`)

These produce different H-bond sets at borderline distances. The
substrate must NOT pick one silently — it exposes both with an
explicit detection parameter:

```python
class HBondDetection(Enum):
    GeometricLarsen = "geometric_larsen"
    DSSP_KabschSander = "dssp_kabsch_sander"
```

- `traj.h_bond_partners(atom, frame_t, detection: HBondDetection) -> list[HBondInfo]`
  where `HBondInfo { Atom partner; HBondRole role; float distance_HA;
  float angle_DHA_rad; LarsenClass class; bool water_eligible; }`
- `role ∈ {Donor, Acceptor}`
- `class ∈ {Primary1, Primary2, Secondary1, Secondary2}` per Larsen
  2015 classification (populated for `GeometricLarsen` detection;
  may be `None` for DSSP detection where Larsen classes don't apply)
- `traj.h_bond_count(atom, frame_t, detection) -> int` — required parameter
- `atom.residue_is_h_bond_donor(frame_t, detection) -> bool` — required parameter

**Original draft** said "geometric (not DSSP)" — Opus catch: that's
wrong, both detection methods exist in the H5 and both produce
publishable signal. The substrate's job is to expose both, not
pick one for the user.

**Provenance**: every Parquet column derived from an H-bond accessor
records its `detection` choice in the column metadata.

## Default cutoffs: forbidden

**Rule**: no substrate method takes a cutoff parameter with a default
value. Every call site specifies the cutoff explicitly.

**Examples**:
- `traj.neighbor_count_by_class(atom, frame_t, target_class, cutoff_Å)` — all four required
- `traj.charge_within(atom, frame_t, cutoff_Å, charge_source)` — all four required
- `traj.fraction_frames_within(atom, target_class, cutoff_Å)` — all three required
- `traj.h_bond_partners(atom, frame_t)` — no cutoff because H-bond geometry has its own thresholds documented in `HBondInfo`

**Recommended cutoff values per use case** (transforms cite these but
don't get them by default):
- Ring current / aromatic neighborhood: 8.0 Å
- H-bond H...A distance: 3.5 Å (Larsen primary cutoff)
- Charged sidechain neighborhood: 6.0 Å
- Solvent water shell: 4.0 Å (first shell)
- General van der Waals contact: 5.0 Å
- Larsen "second-sphere" effects: 6.0 Å

These appear in transform code, not in substrate defaults.

## Target class identification

**Question**: how do transforms specify "neighbors of class C".

**Choice**: predicate dispatch, not string dispatch.

- `traj.fraction_frames_within(atom, predicate, cutoff_Å)` where
  `predicate: Atom -> bool` (Python callable)
- Common predicates exposed as named methods on the Atom binding:
  - `atom.is_water_oxygen()`
  - `atom.is_aromatic_ring_centroid()` (returns False for atoms; True
    only when atom IS a ring centroid pseudo-atom)
  - `atom.is_h_bond_donor()`
  - `atom.is_h_bond_acceptor()`
  - `atom.is_charged()`
  - etc.

**Excluded**: string-based target_class (`'aromatic'`, `'donor'`,
etc.). Same anti-pattern as `feedback_oo_modeling`.

## Provenance metadata sidecar

**Rule**: every walk's output directory contains a `manifest.json`
recording every convention in effect for that walk. JSON (not YAML)
per the 2026-05-30 Qt-primitives alignment pass — eliminates a YAML
parser dependency, matches the existing `StructuredLogger` format,
serializable via Qt's native `QJsonDocument`:

```json
{
  "substrate_version": "0.1.0",
  "run_id": "01234567-89ab-cdef-0123-456789abcdef",
  "conventions": {
    "spherical_harmonics": {
      "convention": "library_isometric_real_sh",
      "_note": "matches src/Types.cpp::Decompose order [xy,yz,zz,xz,xx-yy]; NOT scipy_real (corrected 2026-05-31 — the example contradicted this doc's SH section + SphericalBasis.*)",
      "condon_shortley_phase": true,
      "l_max": 2,
      "basis_ordering_ref": "src/Types.cpp::SphericalTensor::Decompose",
      "normalization": "isometric_real_sh"
    },
    "dipolar_form": {
      "convention": "3cos2_minus_1",
      "reference": "GEOMETRIC_KERNEL_CATALOGUE.md"
    },
    "local_frames": {
      "HN": "v1_amide_plane_via_C_prev",
      "HE_aromatic": "v1_ring_normal_x_anchored_CG_or_CD2",
      "HA": "v1_Cα_chirality_NCα_x_axis",
      "Cα": "v1_bisector_NCα_CCα_axis",
      "CO": "v1_C_to_O_x_in_peptide_plane"
    },
    "charge_source_per_column": {
      "charge_within_columns": {
        "ff14sb": "from_prmtop",
        "mopac": "from_per_frame_mop_output",
        "aimnet2": "from_aimnet2_model_v_X.Y"
      }
    },
    "cremer_pople": {
      "convention_ref": "PlanarGeometryResult.cpp_corrected_2026-XX-XX"
    },
    "multipoles": {
      "quadrupole": "traceless",
      "origin": "target_atom"
    },
    "cutoffs": {
      "policy": "required_no_defaults"
    },
    "charge_units": "elementary_charge_e",
    "distance_units": "angstrom",
    "angle_units_internal": "radians"
  },
  "build_metadata": {
    "git_commit": "<hash>",
    "build_preset": "linux-release",
    "pybind11_version": "...",
    "python_version": "...",
    "pysr_version": "..."
  }
}
```

Written via `QSaveFile` + `QJsonDocument::toJson(QJsonDocument::Indented)`
so the file is atomic-on-commit and human-inspectable. The
`pysr_version` field is null/absent when the manifest is written
before the SDK runs PySR; SDK scripts amend the manifest in place
when they run.

Every SDK script reads this sidecar and refuses to consume a Parquet
whose conventions don't match its own expectations. Prevents
cross-walk comparisons across convention changes from silently
producing wrong answers.

## Substrate version + migration

Substrate version starts at `0.1.0`. Increments per any convention
change in this doc. Convention changes are MAJOR version bumps for
this v0 series; once shipped, any change requires user opt-in (the
sidecar will fail validation against pinned-version SDK scripts).

**Substrate version appears in Parquet filenames** (added 2026-05-30
after Opus validation pass): `walk_<timestamp>_v0.1.0.parquet`. Makes
accidental cross-version SDK runs syntactically loud rather than
validation-time loud. The Parquet filename + the sidecar's
`substrate_version` field both encode the same version; mismatches
between the two are themselves a validation failure.

## PBC — the extractor makes the protein whole; the reader-substrate does NOT do minimum-image (CORRECTED 2026-05-31)

**This section previously said "the substrate's minimum-image-aware accessors
MUST reuse `pbc_whole.h`." That was a category error and it persistently misled
review agents into believing PBC is an open reader-side problem. Corrected at the
root here.**

The **extractor cannot compute anything** — ring current, McConnell,
neighborhoods, EFG — without first **putting the protein back together across the
periodic box**. Making the protein whole is a precondition of every calculation,
done **upstream** in the extractor/library via the ported `pbc_whole`
(`feedback_pbc_verbatim` — DONE there, tested). Therefore:

- The **H5 positions the reader consumes are ALREADY whole.**
- The rediscover substrate does **NOT** do minimum-image PBC for its spatial
  queries — every neighborhood we query is **protein-internal** (~≤10 Å); we
  never query solvent. **`PbcMode = None`** for the substrate. There is **no
  `PbcCellSeries`** and **no reader-side `pbc_whole.h`** — `pbc_whole` is the
  extractor's job, already done, NOT a reader/substrate accessor.
- The substrate's only PBC-related obligation is a **load-time sanity assert**:
  fail loud if a frame's protein arrives wrapped (e.g. a bonded pair separated by
  ~box/2). Whole is the precondition; the assert guards it; that is all.

If a future reader use ever touches solvent or cross-box neighbors, revisit —
but for the protein shielding substrate, PBC is upstream and done.

## Equivalence class typed enum scope (added 2026-05-30 after Opus validation pass)

The substrate's typed equivalence-class enum is NEW infrastructure,
NOT a retrofit of the topology sidecar's existing
`int8_t equivalenceClass` (`h5-reader/src/model/QtAtom.h:64`). The
existing field stays untyped at the producer level; the substrate
introduces its own typed enum and bridges at the substrate boundary
(loads the int8_t, validates against the typed enum, emits as the
typed value to Python). Avoids extractor churn while giving Python
transforms a typed surface.

## Markley pseudoatom name exposure (added 2026-05-30 after Opus validation pass)

The substrate exposes `atom.pseudoatom_kind` returning the Markley
1998 Table 1 pseudoatom name from the existing
`QtAtom::pseudoatomKind` field. Pseudoatom names (`MB` for methyl,
`QB` for methylene-equivalent, `QG` for two-equivalent-CH₃ on Val,
`QQD` for three-fold equivalent groups, etc.) are the
literature-grounded equivalence labels every NMR paper uses; surfacing
them as a Python attribute gives transforms a discoverable, familiar
identifier on top of the substrate's typed equivalence class id.

## What this document does NOT decide

- Specific values of the literature constants used in
  `literature_residual` — those live in the reference equation YAML
  per Stage 2 methodology section
- Specific PySR hyperparameters — SDK script concern
- Cremer-Pople sign convention details — refer to library code
- DSSP definition variant — uses whatever the producer's H5 emits

## Cross-references

- `feedback_t2_sacred` — T2 preservation discipline
- `feedback_kernel_not_shielding` — kernels aren't shielding (residual
  framing)
- `feedback_identity_from_chemistry_not_position` — explicit-chemistry
  rule for symmetry classes
- `feedback_residual_as_ml_feature` — prochiral preservation rule
- `feedback_oo_modeling` — predicate dispatch rule
- `feedback_no_simplification` — no symmetry averaging in substrate
- `feedback_adversarial_review_physics` — Cremer-Pople bug history
- `project_iupac_revert_2026-04-27` — IUPAC-positional trap to avoid
- `GEOMETRIC_KERNEL_CATALOGUE.md` — dipolar form derivation reference
- `SphericalDecomposition.cpp` — T2 basis ordering reference
- `PlanarGeometryResult.cpp` — corrected Cremer-Pople reference
