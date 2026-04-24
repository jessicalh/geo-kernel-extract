# Identity and Dynamics Rollup — Tentative Spec

**Status:** TENTATIVE. In active design conversation. Captured here so the
thinking is not lost across sessions. **Do not treat any enum, field list,
or scope boundary as settled** — this document exists to be iterated on.
Re-read before acting on any of it.

**Date pinned:** 2026-04-22.

**Companion document:** `spec/ChangesRequiredBeforeProductionH5Run.md` —
the single coordination point for pre-fleet library changes. This spec
proposes *expanding* the scope of the next rollup beyond the
NamingRegistry fixes already captured there.

**Intent:** treat this as the last serious opportunity to pay cost
up-front in the H5 surface before Stage 2 / Stage 3 analysis begins.
Once the 10-protein calibration H5s are regenerated and the 485-protein
fleet starts extracting natively, re-regenerating is real work and we
prefer not to do it again.

---

## 1 — Context and motivation

Discovery path:

1. The 2026-04-20 BMRB/RefDB forensics exposed a pattern — when the H5
   atom-name column loses library-level typing discipline, downstream
   Python work (stratification, shift binding, per-atom stats) reaches
   for strings, and careful granularity established in the C++
   extractor is eroded.
2. The `NamingRegistry` coverage gaps already queued in
   `ChangesRequiredBeforeProductionH5Run.md` are one symptom. The
   underlying cause is that the H5 does not yet carry a typed NMR-level
   atom classification. Downstream code has to invent its own, and
   does so with strings.
3. The observation generalises: we have rich per-frame physics in the
   H5 and scattered enrichment booleans on `/atoms/`, but no **unified
   typed per-atom NMR identity**, and no **typed trajectory-level
   continuous summary** of dynamics. Both gaps push work into Python,
   which is where granularity gets lost.
4. Because the 10-protein calibration DFTs are complete (260/260,
   2026-04-18) and tie back to atoms by element-verified ordering,
   a library + schema update **does not invalidate the DFTs**. It does
   require a one-time regeneration of the 10 calibration H5s, and lets
   the 485 pending fleet MDs extract natively into the new schema.

What this spec proposes is **one coherent rollup**: NamingRegistry
fixes, typed per-atom invariants, typed per-trajectory continuous
summaries, a hookable event-extractor pattern, and the H5 schema
changes that carry them. The rollup lands once. Regeneration happens
once. All consumers (Stage 1 ridge, Stage 2 time-series learner,
Stage 3 GNN, viewer, h5-reader, forensics) update against one schema.

---

## 2 — Scope

**In scope:**

- A `NmrAtomIdentity` concept on `Protein` and `QtProtein`: typed
  per-atom classification, topology-determined, invariant across
  conformations.
- A `TrajectoryResult` / `EnsembleResult` on the analysis trajectory:
  typed continuous summaries of dynamics. A new object-model category —
  there is currently no typed home for "summary across the trajectory."
- A `/derived_events/` namespace in the H5: threshold-parameterised
  event lists, shipped as a **menu**, never as primary identity.
- A **hookable extractor stream**: named, typed, parameterised
  operations invocable both by the C++ extractor (at H5 write) and by
  Python (at H5 read), with canonical forms shared across the
  boundary.
- NamingRegistry fixes from `ChangesRequiredBeforeProductionH5Run.md`:
  the eight coverage-gap categories plus the ALA β-methylene wildcard
  narrowing. Lands in the same commit as the rest of this rollup.
- H5 schema version bump.
- Python SDK (`nmr_extract`) catalog extension and typed wrappers for
  every new field.
- Both viewers (`ui/`, `h5-reader/`) updated to decode the new typed
  fields — no re-computation on the viewer side.
- Regeneration of the 10 calibration H5s once, in a single validated
  pass.

**Out of scope (explicit):**

- Re-running any DFT. The 260 stay untouched; alignment is by
  element-verified ordering. (Confirmed from source: `OrcaRunLoader`
  matches by element+position, not atom name.)
- MD trajectories. Positions are unchanged.
- `FindOrcaFiles` alphabetical-glob bug. Separate commit.
- Bond-category enum DRY consolidation. Separate commit.
- `GeometryChoice` build-out. Separate effort.
- Any format change that would break the 260 existing analysis H5s
  from being readable in "legacy mode" — we expect old H5s to be
  unusable after the rollup, which is acceptable because the 10 will
  be regenerated.

---

## 3 — Design principles (hashed out in conversation)

1. **Typed throughout. Strings only at boundaries.** The
   Constitution's anti-pattern list — no string parsing of atom names,
   residue names, element symbols in calculators, condition sets, or
   feature extractors. NmrAtomIdentity's primary tuple is
   `(AminoAcid residue_type, IupacAtomPosition)` where
   `IupacAtomPosition` is itself a per-residue enum. Zero runtime
   strings past the PDB/CHARMM→IUPAC boundary.

2. **Raw signal stays primary.** `/positions/xyz (T, N, 3)` and the
   per-frame physics streams are lossless. Every summary, every
   compression, every event list is **additional**, never a
   replacement. A consumer who doesn't trust a summary can always go
   to the raw.

3. **Summarise, never discretise.** Continuous typed fields
   (covariance tensors, autocorrelations, density grids, PCA modes)
   are acceptable compressions — they preserve the object type and
   invite Python to do further analysis on a typed tensor. Categorical
   bins (rotamer state per frame, Ramachandran bin per frame, ring
   flip events) introduce **competing identities** that invite
   downstream code to pool on the label instead of the continuous
   signal. The discipline: carry typed continuous compressions freely;
   carry categorical discretisations only as explicit
   threshold-parameterised menus under `/derived_events/`, never peer
   to the physics.

4. **Rich is fine. Unused is fine.** The user's explicit preference:
   ship more typed classifiers than strictly needed. A lens that
   doesn't get used costs one column per atom at write time and zero
   at read time. A lens that would have been needed but wasn't shipped
   costs a regeneration round.

5. **Cost discipline: test+apply overhead is per-regeneration, not
   per-library-commit.** Library implementation can chunk into
   multiple commits; what has to be a rollup is the regeneration pass.
   So v1 should include everything we have concrete evidence we'll
   want, plus speculative additions that are cheap (enum columns,
   already-accumulated Welford statistics).

6. **Library computes what would lose granularity in Python.**
   Stratifications, classifications, and expensive trajectory scans
   (PCA, autocorrelation, clustering) live in C++ where we can be
   careful. Cheap derivations (FFT of a per-atom time series, simple
   thresholds applied to continuous signals) stay Python-side.

7. **Hookable, not baked-in.** Event extractors and threshold choices
   are research questions. The library ships a menu and exposes the
   same named operations in the Python SDK so consumers can run fresh
   extractions at load time with different parameters. Canonical forms
   shared across C++/Python so the invocation looks the same on both
   sides.

8. **Object-model fidelity.** The existing split — `Protein` is
   identity and topology (invariant), `ProteinConformation` is one
   geometric instance — is preserved. NmrAtomIdentity lives on
   `Protein` because its content is topology-determined (stereocentre
   identity from bonded graph + L-chirality + IUPAC convention, not
   from geometry). The new `TrajectoryResult` / `EnsembleResult` is a
   third object-model category — not a single conformation, not the
   bare protein, but the typed summary of a whole run. This closes an
   object-model gap that was pushing dynamics analysis into untyped
   Python.

9. **pro-R / pro-S is topology-determined, not geometry-determined.**
   Standard MD preserves CA chirality and bonded topology. IUPAC
   naming encodes the stereochemistry (VAL CG1 is pro-R by
   convention; CG2 is pro-S). A correct NamingRegistry translation IS
   the stereochemistry assignment. No per-frame re-classification.
   Methyls flip about their local axis inside a geometric
   twin — the *labels* don't.

---

## 4 — Part I: per-atom invariants (`NmrAtomIdentity`)

**Attach site:** on `Protein` (library) and `QtProtein` (h5-reader,
parallel). Constant across all conformations. Static.

**Primary identity tuple** (this is what BMRB/RefDB binding keys on):

- `AminoAcid residue_type` — existing enum, unchanged.
- `IupacAtomPosition` — new per-residue enum. Every valid
  (residue_type, position) pair enumerated once from `AminoAcidType`.
  Zero runtime strings after load. Example: `IupacAtomPosition::ILE_CD1`,
  `IupacAtomPosition::ILE_HD11`, `IupacAtomPosition::ALA_HB3`.

**Semantic lenses** (all derived from topology + IUPAC convention):

- `NmrClass` — unified NMR category enum absorbing and extending the
  scattered existing booleans (`is_amide_H`, `is_alpha_H`, `is_methyl`,
  `is_aromatic_H`, etc.). Proposed values:
  - Backbone: `BackboneHN`, `BackboneHA`, `BackboneN`, `BackboneCA`,
    `BackboneCO_C`, `BackboneCO_O`
  - Sidechain nitrogen: `SidechainAmideNH`, `SidechainAmideN`,
    `GuanidiniumNH`, `GuanidiniumN`, `AmmoniumNH`, `AmmoniumN`,
    `RingNH`, `RingN`
  - Sidechain oxygen: `HydroxylOH`, `HydroxylO`, `CarboxylOH`,
    `CarboxylO`, `AmideCO_O` (sidechain)
  - Sidechain sulfur: `ThiolSH`, `ThiolS`, `ThioetherS`,
    `DisulfideS`
  - Sidechain carbon: `SidechainC_sp3`, `SidechainC_sp2_carbonyl`,
    `SidechainC_aromatic`, `SidechainC_sp2_guanidinium`
  - Sidechain hydrogen: `MethylCH3`, `MethyleneCH2`, `MethineCH`,
    `RingCH`
  - Other: `CappingN`, `CappingC`, `Other`
- `Locant` — sidechain depth letter: α / β / γ / δ / ε / ζ / η /
  `Backbone` / `NotApplicable`.
- `SidechainDepthBonds` — integer bond-count from nearest backbone
  atom (redundant with `Locant` but numerically useful).
- `MethylGroup` — `ALA_β`, `VAL_γ1`, `VAL_γ2`, `LEU_δ1`, `LEU_δ2`,
  `ILE_γ2`, `ILE_δ1`, `THR_γ2`, `MET_ε`, or `None`. pro-R / pro-S
  encoded by which enum slot the atom lands in.
- `methyl_partner_atom_index` — for the twin pairs (VAL γ1↔γ2,
  LEU δ1↔δ2), the partner atom's index into protein.atoms. `SIZE_MAX`
  otherwise.
- `MethyleneGroup` — analogous enum for methylenes.
- `methylene_partner_atom_index` — prochiral partner for HA2/HA3,
  HB2/HB3, HG2/HG3, HD2/HD3, HE2/HE3. Which slot is pro-R vs pro-S
  from IUPAC convention.
- `ChiParticipation` — bitfield. Bit k set if this atom is one of the
  four dihedral-defining atoms for `chi_k` (k ∈ {1,2,3,4}).
- `RingAtomRole` — ring-member role. For PHE/TYR: `Ipso`, `Ortho`,
  `Meta`, `Para`. For TRP: distinct labels per ring (TrpBenzene
  `C4/C5/C6/C7`, IndolePerimeter, TrpPyrrole positions including
  `Nε1`). For HIS tautomers: `Nδ1`, `Cε1`, `Nε2`, `Cδ2`, `Cγ`.
  `NotInRing` otherwise.
- `ring_membership` — bitfield of which ring indices (in
  `protein.Rings()`) this atom belongs to. Fused rings mean an atom
  can be in multiple.

**Residue-level invariants, broadcast per-atom for fast slicing**:

- `ResidueCategory` — `Hydrophobic`, `Polar`, `AcidicCharged`,
  `BasicCharged`, `Aromatic`, `GLY`, `PRO`, `CYS_disulfide_bonded`,
  `CYS_free`.
- `ResidueHydropathy` — Kyte-Doolittle scalar.
- `ResidueVolume` — Chothia scalar.
- `ResidueIsoelectricPoint` — scalar.
- `ResidueHbondDonorCount`, `ResidueHbondAcceptorCount` — ints.

**Symmetry / pseudo-atom info**:

- `SymmetryClassId` — atoms that are chemically equivalent under
  topology and local symmetry (methyl Hs, PHE Hε1/Hε2, Hδ1/Hδ2)
  share an ID. For downstream pooling when the analysis calls for it.
- `PseudoAtomName` — BMRB-compliant collapsed name (HB# for
  methylenes, HG# for methyls). At the edge. String in storage, not
  used in C++ logic.

**Orthogonal indices** (pre-built on `Protein`, O(1) slices):

- `Protein::AtomsByNmrClass(NmrClass) -> const vector<size_t>&`
- `Protein::AtomsByLocant(Locant)`
- `Protein::AtomsByMethylGroup(MethylGroup)`
- `Protein::AtomsByRingAtomRole(RingAtomRole)`
- `Protein::AtomsByResidueCategory(ResidueCategory)`
- `Protein::AtomsByResidueAndPosition(AminoAcid, IupacAtomPosition)`
- `Protein::MethylGroups() -> vector<MethylGroupInfo>` (three atoms
  each, plus twin-pair linkage)
- `Protein::MethyleneGroups()`

Parallel to the existing `atoms_by_role`.

**Implementation note:** this is largely an extension of
`EnrichmentResult`. The existing `AtomRole` stays; `NmrClass` is the
finer semantic overlay. Old booleans on `ConformationAtom`
(`is_amide_H` etc.) can either be kept as redundant convenience or
deprecated in favour of `NmrClass` predicates. Decide in implementation.

---

## 5 — Part II: per-trajectory continuous summaries (`TrajectoryResult` / `EnsembleResult`)

**Attach site:** new object-model category. Not a
`ConformationResult`. Computed during the two-pass trajectory run in
`nmr_extract --trajectory --analysis`. Serialised into H5 under
`/ensemble/` or `/trajectory/` (choose in implementation).

**Motivation:** `Protein` is invariant, `ProteinConformation` is one
frame. There is currently no typed home for "summary across the
trajectory." The Welford accumulators in `GromacsProtein` are half of
this — they already accumulate per-scalar mean and variance per atom.
Lift the discipline from scalars to tensors, to autocorrelations, to
density fields.

**Per-atom trajectory tensors** (`T → 1` reductions):

- `mean_position` (N, 3) — trajectory-mean atomic position.
- `position_covariance` (N, 3, 3) — full anisotropic displacement
  tensor per atom. RMSF is its trace/3; the tensor carries directional
  anisotropy.
- `position_covariance_spherical` (N, 9) — SphericalTensor companion.
  T0 (isotropic RMSF), T2 (anisotropy). Same dual-representation
  discipline as the per-frame physics.

**Per-observable mean + covariance tensors** (extend Welford from
scalar to tensor):

- For every per-frame scalar/tensor field already in the H5
  (`bs_shielding`, `mc_shielding`, `coulomb_EFG`, `total_B_field`,
  `aimnet2_embedding`, `apbs_efield`, `water_efield`, …): trajectory
  `mean_<field>` and `cov_<field>` as typed tensors alongside the raw
  per-frame stream.

**Per-bond-vector / per-ring-normal continuous dynamics** (NMR-classic
continuous signals, no thresholds):

- `bond_vector_autocorrelation` for NMR-critical bond sets
  (backbone N-H, CA-HA, C=O, aromatic C-H): `(n_bonds, n_lags)`. The
  continuous function S² and τ_e are *derived from*. Ship the
  function; ship the fit too.
- `order_parameter_S2` — `(n_bonds,)`. Derived from the above.
  Classical NMR-MD bridge. One number per bond.
- `lipari_szabo_tau_e` — effective correlation time per bond from the
  fit.
- `ring_normal_autocorrelation` — `(n_rings, n_lags)`. What a "ring
  flip rate" would compress to; carry the continuous function instead.
- `ring_normal_covariance` (n_rings, 3, 3) — anisotropic flip
  signature without the flip event.

**Per-amide cross-correlated relaxation — back-calculated NMR
observable** (added 2026-04-22 per Session 0 section-5 walkthrough):

The σ tensor output and per-frame bond vectors together permit
direct back-calculation of dipole-dipole / CSA cross-correlated
relaxation (DD/CSA CCR) rates at every amide N. This is the first
native NMR observable we predict from the classical calibrated
kernels rather than as a correlation to measured shifts, and
nothing else at ~4000-atom scale produces per-frame calibrated
σ(T0+T1+T2) to feed this computation. Recorded here as a design
attempt so the scoping and typing persist even if validation
fails; if it fails the autocorr and principal-axis fields remain
useful raw material, only the derived rate is fragile.

- `nh_dipolar_csa_angle_autocorr` (n_amide_N, n_lags) —
  autocorrelation function of the second Legendre polynomial
  P₂(cos θ_DD,CSA) along the trajectory, where θ_DD,CSA is the
  per-frame angle between the NH bond vector and the amide
  ¹⁵N σ tensor's principal axis.
- `nh_dipolar_csa_ccr_rate` (n_amide_N,) — derived CCR rate
  Γ_DD,CSA in s⁻¹, computed from the long-time limit of the
  autocorrelation above and the overall tumbling correlation time
  τ_c. Stored with assumed τ_c in H5 attributes.
- `nh_dipolar_csa_principal_axis` (n_amide_N, 3) — per-residue
  mean direction of the σ tensor's dominant T2 eigendirection
  over the trajectory. Diagnostic: the CSA principal axis should
  drift only slightly with dynamics; large drifts flag
  tensor-extraction pathology rather than real physics.

References: Tugarinov-Kay 2003 (bib M10), Ferrage-Piserchio-
Cowburn-Ghose 2008 (bib M11), Fushman-Tjandra-Cowburn 1998
(bib M12). Published Γ_DD,CSA values exist for GB3, ubiquitin,
and other well-characterised proteins — the external validation
target.

Caveats baked into the typing:
- τ_c for overall tumbling is not converged from a 50 ns window
  (protein τ_c is typically of order the window length or longer).
  Ship the autocorr function as primary; downstream recomputes
  Γ with a literature τ_c per protein and documents what it used.
- CSA principal-axis identification is operationally the largest-
  magnitude eigendirection of the per-frame T2 (symmetric
  traceless) part of σ. Stable for amide ¹⁵N where CSA is well
  defined; skip atoms where the T2 Frobenius norm is below a
  documented threshold (noise-dominated identification).
- The back-calculated rate is dominated by (σ_∥ − σ_⊥), i.e. the
  amide ¹⁵N CSA magnitude. Sanity-check against Fushman-Tjandra-
  Cowburn 1998 residue-specific CSA measurements before trusting
  Γ comparisons.

**Conformational-space density fields** (rotamer / Ramachandran as
continuous fields):

- `chi_density_grid` — per residue per chi angle, a coarse grid
  (e.g. 36 bins). For residues with multiple chis, joint density
  (`36 × 36` for chi1 × chi2).
- `phi_psi_density` — per residue, `36 × 36` grid.

(If storage is an issue on the fleet, an alternative is a compact
parametric fit — von-Mises mixture — per residue. Decide in
implementation after one 10-protein pass.)

**Collective-motion** (continuous):

- `residue_residue_covariance` (R, R, 3, 3) — full per-pair 3×3
  covariance tensor. **Typing amended 2026-04-22 per Session 0
  section-5 walkthrough:** the scalar-per-pair form originally
  drafted here (R, R) loses the directional coupling axis that
  dynamic-allostery analyses depend on (essential dynamics,
  Amadei-Linssen-Berendsen 1993 bib M19; entropic allostery,
  Cooper-Dryden 1984 bib M23; Frederick-Wand / Sharp-Wand
  conformational-entropy programme bib M8/M9). Storage trivial
  (100 residues → 720 KB per protein). Downstream analyses that
  only need the magnitude take the Frobenius norm of this tensor.
- `residue_residue_covariance_magnitude` (R, R) — derived scalar
  Frobenius-norm view of the above. Quick-look companion.
- `pca_eigenmodes` (top-k, 3N) — top-k principal components of motion.
- `pca_eigenvalues` (top-k,).
- Top-k default 20. Storage cheap.
- All-atom-pair covariance: skip unless explicitly requested. Size
  gets real on 4000-atom systems.

**Local-environment continuous summaries** (where serendipity is most
likely per the conversation):

- `coulomb_E_variance` per atom — variance of local Coulomb E-field
  over trajectory. Detects breathing of the charge environment.
- `nearest_ring_geom_variance` per atom — variance of `(z, rho,
  theta, cos_phi, sin_phi)` against the nearest ring over the
  trajectory. The atom's continuous view of ring motion.
- `shielding_contribution_variance` per calculator per atom — "this
  atom's physics environment is quiet vs. breathing."

**Explicit non-members:** no ring-flip events, no rotamer-transition
events, no H-bond formation/breaking events, no SS-transition events
as part of `TrajectoryResult`. Those belong in
`/derived_events/`, see Part III.

---

## 6 — Part III: derived events as a threshold-parameterised menu

**Motivation from conversation:** for a GNN on ~4000 nodes per
protein, discrete event messages may genuinely be the right input
shape. But we do not know which events — ring flip at 90° vs 120°,
occlusion within 5Å vs 7Å, rotamer transition at stringent vs
lenient threshold. The threshold choice IS the research question.
Picking one and baking it in closes the question before we have asked
it.

**Pattern:**

- Events live under `/derived_events/` in the H5. Namespaced to make
  the competing-identity hazard explicit — a consumer has to *reach
  for* them and has to *pick a threshold*.
- Multiple thresholds ship as peer datasets:
  `/derived_events/ring_flip_90deg`,
  `/derived_events/ring_flip_120deg`,
  `/derived_events/ring_occlusion_5A`,
  `/derived_events/ring_occlusion_7A`,
  `/derived_events/rotamer_transition_stringent`,
  `/derived_events/rotamer_transition_lenient`,
  `/derived_events/hbond_lifetime_gaps`,
  etc. Exact menu to be decided; err rich.
- Each event dataset is a list of records —
  `(frame_index, atom_or_residue_or_ring_index, additional_typed_fields)`.
  Variable-length per protein.
- Documentation for every event dataset: threshold value, definition,
  what continuous field it was derived from. The
  documentation is load-bearing — a GNN consumer that picks a
  threshold has to see the full specification.

**A GNN experiment can ablate across thresholds** by selecting
different datasets, with zero extractor changes.

---

## 7 — Part IV: the hookable stream — named extractors across library / Python

**New discipline from 2026-04-22 conversation:** event extractors
(and potentially other derivations) are named, typed, parameterised
operations. The same operation name is invocable two ways:

1. **At H5 write time in C++.** The extractor runs during
   `nmr_extract --trajectory --analysis` to populate the
   `/derived_events/` datasets. One run produces the menu of
   parameter variants we ship.
2. **At H5 read time in Python.** The SDK exposes the same named
   operation with the same parameter names. Python consumers
   (particularly GNN experiments) can invoke the extractor fresh on a
   loaded trajectory with different thresholds than we shipped.

**Shared canonical form across C++ and Python.** Each extractor
declares a name, a parameter schema, input-field dependencies, and an
output record shape. The C++ and Python implementations are both
"canonical Python forms" in the sense that the invocation and output
type are the same in both. The Python side may be a thin wrapper
around a C++ implementation or a native re-implementation — decide
case by case based on speed needed.

**Example extractor contract (sketch):**

```
name:          ring_flip
parameters:    { angle_threshold_deg: float, min_dwell_ps: float }
inputs:        ring_normal_timeseries (derived from /positions/xyz
               and /topology/ring_atom_indices, so self-contained on
               H5 read)
output shape:  list of { frame_index: int, ring_index: int,
                         pre_normal: Vec3, post_normal: Vec3,
                         flip_angle: float }
canonical invocation:
  C++:    DerivedEventExtractors::RingFlip(ring_normal_ts,
                                           angle_threshold_deg=90,
                                           min_dwell_ps=50)
  Python: derived_events.ring_flip(ring_normal_ts,
                                   angle_threshold_deg=90,
                                   min_dwell_ps=50)
```

**Why this matters:**

- The extractor is the authority. The shipped datasets are one choice
  of parameters against it; they are not themselves the spec.
- Stage 3 GNN experiments can rerun with new thresholds without a
  library touch.
- Stage 1 / Stage 2 stats can decide they want a different threshold
  than we shipped and compute it in Python against the raw signals,
  with the same output shape as a shipped dataset.
- The glossary can describe the extractor once, and every shipped
  dataset inherits the description.

**Implementation note:** start with the 5-8 most obvious extractors
(ring flip, ring occlusion, rotamer transition, H-bond lifetime, SS
transition, maybe one or two NMR-specific). Each extractor needs a
doctring and a unit test. The canonical-form split is a v1 discipline
— if something turns out to be hard to make canonical across
C++/Python, ship it only on one side and document which.

---

## 8 — Part V: H5 schema changes

**`/atoms/` extensions** (invariant, N-wide columns):

- `iupac_atom_position` — int, maps to `IupacAtomPosition` enum
- `nmr_class` — int, maps to `NmrClass`
- `locant` — int, maps to `Locant`
- `sidechain_depth_bonds` — int
- `methyl_group` — int, maps to `MethylGroup`
- `methyl_partner_atom_index` — int
- `methylene_group` — int
- `methylene_partner_atom_index` — int
- `chi_participation` — int, bitfield
- `ring_atom_role` — int, maps to `RingAtomRole`
- `ring_membership` — int, bitfield
- `residue_category` — int
- `residue_hydropathy` — float
- `residue_volume` — float
- `residue_isoelectric_point` — float
- `residue_hbond_donor_count` — int
- `residue_hbond_acceptor_count` — int
- `symmetry_class_id` — int
- `pseudo_atom_name` — string (at the edge)

**New `/ensemble/` group** (trajectory-level continuous summaries):

- `/ensemble/mean_position (N, 3)`
- `/ensemble/position_covariance (N, 3, 3)`
- `/ensemble/position_covariance_spherical (N, 9)`
- `/ensemble/mean_<field>` and `/ensemble/cov_<field>` for each
  per-frame observable
- `/ensemble/bond_vector_autocorr` for NH / HA / CO / aromatic-CH
- `/ensemble/order_parameter_S2`
- `/ensemble/lipari_szabo_tau_e`
- `/ensemble/ring_normal_autocorr`
- `/ensemble/ring_normal_covariance`
- `/ensemble/chi_density_grid`
- `/ensemble/phi_psi_density`
- `/ensemble/residue_residue_covariance`
- `/ensemble/pca_eigenmodes`
- `/ensemble/pca_eigenvalues`
- `/ensemble/coulomb_E_variance`
- `/ensemble/nearest_ring_geom_variance`
- `/ensemble/shielding_contribution_variance` (per calculator)

**New `/derived_events/` namespace** (threshold-parameterised menu —
see Part III):

- `/derived_events/<extractor_name>__<params>` — one dataset per
  shipped parameter variant. Variable-length records.
- Attributes on each dataset record the extractor name, parameters,
  and schema version.

**Schema version attribute** on the root of every H5:
`analysis_h5_schema_version = 2` (current is 1). Readers check this on
load and refuse to open unknown versions. Old H5s from v1 are
intentionally unreadable after the rollup; the 10 are regenerated.

**Catalog contract:** every new NPY/H5 dataset needs an `ArraySpec`
entry in `python/nmr_extract/_catalog.py` and a typed wrapper class.
Non-negotiable per `spec/EXTRACTION_SDK.md`.

---

## 9 — Part VI: NamingRegistry activation

This rollup *is* the activation of the NamingRegistry fixes from
`spec/ChangesRequiredBeforeProductionH5Run.md`. Concretely:

- Narrow the wildcard β-methylene rule to exclude ALA (ALA HB1/HB2/HB3
  passes through).
- Add the eight missing coverage categories (ILE γ-carbon, ILE
  δ-methyl, ILE γ1-methylene, ARG/LYS/PRO δ-methylene, LYS
  ε-methylene, GLY α-methylene).
- Re-run `_probe_naming_conventions.py` over the 10 calibration
  proteins to confirm uniformity. (Original ChangesRequired activation
  criterion was to re-probe all 685; with a reduced initial scope —
  10 calibration — we confirm there first and again on the fleet when
  PDBs land.)
- `pack_experimental_shifts.py` local constants become unnecessary
  after regeneration. Document their retirement but leave them in
  place during transition.

Post-rollup activation steps documented in the existing
ChangesRequired file get updated: "status: APPLIED 2026-MM-DD, commit
SHA …" once the rollup lands.

---

## 10 — What does NOT change

- MD trajectories (`md.xtc`, `md.tpr`, `md.edr`). Regenerable but
  untouched.
- 1 ns pose PDBs already on disk for the 10 calibration proteins. They
  stay. If regenerated (optional), positional content is identical and
  atom_name column becomes fully IUPAC — strictly better for tleap.
- 260 completed ORCA DFT outputs. Alignment is by element-verified
  ordering (`OrcaRunLoader` + `OrcaShieldingResult`), independent of
  atom-name translation. DFTs remain valid post-rollup.
- Stage 1 calibration results on the 720-protein mutation set. The new
  per-atom-type classification is designed to reproduce the current
  AMBER-atom-name strata exactly. An explicit mapping table
  `NmrClass × Locant × element → Stage1_stratum` is a required
  artefact of the rollup; Stage 1 results should regenerate
  identically on the 10-protein subset as a regression check.

---

## 11 — What regenerates

- The 10 calibration H5s under `fleet_calibration-{working, stats,
  backup}/`. One pass, two-parallel, ~2-4 hours wall. Validated
  individually before fleet extraction activates.
- 1 ns pose PDBs may be regenerated alongside if the user wants
  fully-IUPAC atom_name in the PDB column. Optional; positional
  content is identical to existing.
- The 485 pending fleet MDs extract natively into the new schema once
  the rollup is activated. Zero extra cost — they were going to be
  extracted anyway; they simply arrive with the new fields populated.

---

## 12 — Known unknowns and open questions

Captured 2026-04-22 from conversation. These are the places where we
explicitly do not yet know the answer and the spec may change.

1. **Exact `NmrClass` enum contents.** The proposed values above are a
   first pass biased rich. Specific values — especially around
   sidechain carbons (sp2 vs sp3 vs aromatic sub-distinctions) — need
   a hand-derived hold-out test on all 20 canonical amino acids
   before they lock. Advisor / NMR-experimentalist input desirable.

2. **Diastereotopic pairing convention.** The spec currently encodes
   pro-R / pro-S in which enum slot the atom lands in (VAL CG1 →
   `VAL_γ1` → pro-R by IUPAC convention). Alternative: keep the slot
   agnostic and emit a separate `prochirality` enum
   (`ProR` / `ProS` / `Ambiguous`). Second option may be cleaner for
   Python slicing. Decide in implementation.

3. **Ring-position role for histidine tautomers.** HIS, HID, HIE have
   different protonation and thus different ring-position
   identities. The `RingAtomRole` enum needs tautomer-aware values.
   Verify with an NMR-HIS-expert perspective before locking.

4. **NMR-specific continuous objects we may be missing.** Conversation
   flagged:
   - Residual dipolar couplings (RDCs) in aligned ensembles — carry
     phase-like orientation information, not an autocorrelation.
   - Cross-spectral density between bond vectors — correlated motion
     in a way autocorrelation misses.
   - Transferred NOE quantities.
   - Windowed spectrogram of per-atom positions (short-time FFT) for
     detecting transient coherent motions — decided to defer in v1,
     Python-side FFT is cheap enough that naive per-position FFT is
     fine for most cases.
   Need to research. This is explicitly flagged as an advisor /
   literature question.

5. **Extractor menu for `/derived_events/`.** Which events ship in
   v1? Conversation suggested 5–8 obvious ones:
   ring_flip(90°, 120°), ring_occlusion(5Å, 7Å),
   rotamer_transition(stringent, lenient), hbond_lifetime_gaps,
   ss_transition. Exact list deferred until the hookable-stream
   pattern is built and we can see what's cheap to add.

6. **`all_atom_pair` covariance.** On 4000-atom proteins this is a
   `(4000, 4000, 3, 3)` tensor per protein = ~2 GB. Default no;
   decide per-use-case.

7. **Stage 1 strata reproduction.** The exact mapping
   `NmrClass × Locant × element → Stage1_stratum` needs to be drafted
   and verified before the rollup is irreversible. A 10-minute
   exercise but must happen.

8. **Python prototype for dynamic annotations.** Rotamer bin
   boundaries, ring-flip angle threshold, H-bond lifetime threshold —
   algorithmic choices with research-question character. Conversation
   agreed these should be Python-prototyped against existing H5s
   before C++ finalisation. In-scope for the rollup
   development path; out-of-scope for the spec to dictate.

9. **Hookable stream: where do canonical extractor forms live in the
   codebase?** Candidates: a new `src/derived_events/` directory in
   the library, a `python/nmr_extract/derived_events/` module
   mirroring it, with a shared `SCHEMA.md` describing the contract.
   Decide in implementation.

10. **Viewer decode work.** Both `ui/` and `h5-reader/` need to
    decode the new `/atoms/*` typed fields and the `/ensemble/`
    group. Feature plan impact in `h5-reader/notes/FEATURE_PLAN.md`
    to be assessed.

---

## 13 — 2026-04-22 pm addendum: MD's distinctive payload and non-shielding NMR endpoints

After sections 1–12 were drafted, the conversation turned on two points
that shape the spec's motivation and, in one place, add a concrete
field. Captured here as discussion, not as settled doctrine, so a
future reader receives them forward in time rather than as a priori
principles.

### 13.1 — What MD is actually bringing (given that MD→shielding is modestly weak)

The MD→shielding link is modestly weak: ff14SB / CHARMM36m are not
tuned against DFT shielding, the 1 ns snapshot cadence assumes
ergodicity at that interval, and the advisor chose MD as the
conformational sampler. Taken as given. The question the conversation
sharpened: what is MD actually bringing that static QM on a crystal
structure would not? The `TrajectoryResult` fields in section 5 map
onto this directly:

- **Conformational populations, not single minima.** `chi_density_grid`,
  `phi_psi_density`.
- **Solvation dynamics.** `coulomb_E_variance` per atom captures
  breathing of the local charge environment; APBS-on-a-single-frame
  has no access to this.
- **Methyl rotation, ring flipping, rotamer-transition rates.**
  Bond-vector and ring-normal autocorrelations. Static QM is mute here.
- **Large-amplitude, anharmonic, 300 K fluctuations.** Position
  covariance tensors, PCA modes. Harmonic QM vibrational analysis
  gives modes at the minimum and is the wrong object.
- **H-bond lifetime distributions.** Derivable from the existing
  per-frame H-bond geometry and energy streams; a distribution
  object would be a cheap addition.
- **Collective motion.** Residue–residue covariance, PCA projections.

**Spec implication (deferred to implementation-pass rewrite, not
applied now):** section 5's opening should name these as the specific
payload MD brings, not just "useful summaries." Right now section 5
reads as a list of fields; the motivation is implicit and a reader
who hasn't seen this conversation does not recover it.

### 13.2 — NMR is a black art because it is not just shielding

User framing, recorded close to verbatim: "NMR is not just shielding
but what happens overall when you whack a molecule with massive EMF;
if it were just shielding, it would not be such a black art."

Our 260-DFT-calibrated kernel decomposition is a novel first-principles
story on *one* axis — the shielding tensor (T0 + T1 + T2).
Experimental NMR observables are superpositions. What each one needs
from the H5:

| Observable | Handles needed from our pipeline |
|---|---|
| J-couplings (scalar) | Dihedrals (already in `/dihedrals/`) |
| Dipolar couplings (ave. out in isotropic) | Bond orientations per frame |
| RDCs (aligned media) | Bond-orientation tensor over trajectory — **candidate addition, see 13.3** |
| NOE | `1/r⁶` distance distributions + correlation times (positions + new bond-vector autocorrelations) |
| R1 / R2 / R1ρ relaxation | Spectral densities at Larmor — bond-vector autocorrelations in section 5 |
| Chemical exchange (CEST, CPMG) | Conformational-state populations (`chi_density_grid`); exchange rates are typically slower than our 25 ns window, so populations are accessible but kinetics generally are not |
| Quadrupolar (²H, ¹⁴N, ¹⁷O; spin > 1/2) | EFG variance per atom over trajectory (follows from the general "mean + covariance for every per-frame observable" rule; call out explicitly — see 13.3) |

**Framing implication:** the rollup serves two purposes independently.
The Stage 1 / Stage 2 / Stage 3 learning systems on the shielding axis,
AND a typed handle-set for every non-shielding NMR observable an
analyst might want. Downstream work (advisor's group, future students,
collaborators) gets the MD-informed continuous summaries without
needing to re-extract. That independence deserves a sentence in
whatever the thesis scoping chapter becomes.

### 13.3 — Concrete spec additions from this discussion

Three items: one new field, two framing notes. None change the scope
boundary (section 2) or the activation criteria (section 9).

**(a) Add `bond_orientation_tensor` to section 5's NMR-continuous
block.** `(n_bonds, 3, 3)` — the full orientation-correlation tensor
over the trajectory for NH / HA / CO / aromatic-CH bonds. This is
what RDC analysis reads, and it is what a phase-sensitive
frequency-domain learner on bond dynamics would want. 9 doubles per
bond times a few thousand bonds: trivial storage. Flagged for the
section 5 field list at implementation time. **Typing confirmed
2026-04-22 per Session 0 section-5 walkthrough:** the 3×3 form is
exactly the object consumed by Meiler-Prompers-Peti-Griesinger-
Brüschweiler 2001 (bib M14) and the Lakomek et al. 2008 SCRM
RDC-model-free analysis (bib M15). Ubiquitin and GB3 RDC datasets
from those works are the natural external validation bench.

**(b) Call out EFG-covariance-per-atom explicitly in section 5.** The
quadrupolar-relaxation handle is already implicit in the general
"mean and covariance for every per-frame observable" rule, but
`cov_apbs_efg` and `cov_coulomb_EFG_total` should ship by named
default, not by rule-coverage. Makes the intended consumer obvious to
a reader who doesn't know what quadrupolar relaxation is.

**(c) Reframe section 5's opening at implementation time.** Currently
reads "here are some useful fields"; should read "here is the payload
MD brings that static QM cannot, and the handles that non-shielding
NMR observables need." This reframing is the thing that keeps the
motivation visible without a spec reviewer having to reconstruct it
from session transcripts.

These three are the only carry-forward actions from this addendum.
Deferred to the implementation-pass rewrite of section 5.

### 13.4 — Posture for the research phase

User framing at the end of the conversation, kept in their words:
"our actual power is in these 260 DFTs, and our novel application is
a pretty good proxy for shielding from maths. We can think quantum
mechanics and the things we grab in our kernels as much as what NMR
people are doing as such."

Direction for the literature pass that follows spec stabilisation:
read dynamics / NMR / QM literature with the project's novel angle in
mind — shielding-from-maths via DFT-calibrated geometric kernels plus
MD continuous summaries — not as NMR experimentalists reading about
spectrum assignment. The kernel side and the DFT side are where the
novelty lives; the NMR surface is where the output eventually has to
connect, but it is not the lens through which the research phase
should be conducted.

### 13.5 — Note on kernel-kernel covariance as signal

Tucking this in so implementation doesn't lose the thread. Kernel
time series and kernel-kernel cross-correlations over the trajectory
carry information about protein physics — local stiffness, coupled
sidechain modes, rotameric switching, allosteric connection,
solvent-shell reorganisation — that is related to but not itself
shielding. Geometry-sampled kernels inherit the correlation structure
the protein's constraints (chain, rings, disulfides, packing) impose
on atom motion. The section 5 ensemble fields are positioned to
surface this, but what it actually looks like is a signal-at-scale
question that only the 10-protein data can answer once regenerated.
The first move there is an observability-style pass — distributions
per kernel, kernel-kernel covariance matrix, spectral density per
calculator, outliers, bimodalities — before any modelling choice.

### 13.6 — Validation of structural assignments

Added late in the 2026-04-22 pass. The NmrAtomIdentity classification
needs external and internal validation before the rollup commits;
silent mis-assignment corrupts every downstream binding.

**External references (ground truth):**

- PDB Chemical Component Dictionary (CCD) atom name lists per residue
  — authoritative for (residue_type, atom_name) pairs.
- IUPAC-IUB nomenclature recommendations for amino acid atom naming,
  particularly for diastereotopic pro-R / pro-S assignments.
- BMRB atom name conventions (the experimental shift side of the
  boundary).

Validate by building a fixture table from CCD + IUPAC + BMRB and
asserting our `IupacAtomPosition` enum matches for all 20 canonical
amino acids plus every protonation variant we support.

**Internal consistency checks (at enrichment, fail-hard on mismatch):**

- Every atom gets exactly one `NmrClass`. Non-overlapping coverage.
- Every ring-member atom gets exactly one `RingAtomRole`.
- Methyl groups: exactly 3 H atoms bonded to one C (topology check).
- Methylene groups: exactly 2 H atoms bonded to one C.
- `ChiParticipation`: each chi declared in `AminoAcidType` must have
  exactly 4 atoms in the protein with that bit set.
- Bidirectional name translation identity: CHARMM → IUPAC → CHARMM
  round-trips. Verifies `NamingRegistry` completeness directly; would
  have caught the ALA β wildcard bug.
- Hand-derived fixtures for all 20 canonical AAs match the computed
  classification exactly. See section 12 item 7; this is the concrete
  test.

**pro-R / pro-S geometric cross-check:**

Static MD preserves chirality, but the *input structure* can have a
mis-assigned atom name (e.g., HB2 and HB3 swapped in an incoming PDB).
A one-time geometric CIP-priority calculation at enrichment, compared
against the IUPAC-encoded label, catches these. Fail hard on
mismatch — the input is wrong, not our classification.

**Placement:** validation belongs on the Layer 0 build path
(`FinalizeConstruction` or equivalent post-build step) and as a
dedicated test file `tests/test_nmr_identity_validation.cpp` in the
rollup's test surface. Ships as part of the rollup, not after.

### 13.7 — HDF5 as universal output format

Conversation direction: promote `fileformat/` to be the sole output
path for everything the library writes. The output surface today is
split — analysis H5 via `AnalysisWriter::WriteH5`, per-result NPY via
`ConformationResult::WriteFeatures` / `WriteAllFeatures`, calibration
pipeline outputs under `calibration/{ID}/`, mutation-delta outputs,
smoke-test fixtures. Unifying them:

- Every library output goes through `fileformat/`.
- NPY deprecated as a first-class output — can still exist as a helper
  `h5 → npy` export for legacy consumers, but is not the library's
  responsibility.
- Single schema, single reader, single validator.
- Formal test framework in `fileformat/test/`: roundtrip tests per
  `ArraySpec`, schema-version tests, boundary-case tests (NaN, empty,
  oversized, wrong-shape), cross-platform byte-exactness tests.
- SDK `_catalog.py` remains the format contract — one `ArraySpec` per
  dataset, no exceptions, generated from or synchronised with the
  library's field manifest (see 13.8).

**Scope note:** this is bigger than the NmrAtomIdentity + EnsembleResult
scope in sections 4–5. Likely its own coordinated pass, either
alongside this rollup or as its immediate follow-on. Flagged here so
the discussion preserves it, not decided.

### 13.8 — Field-inclusion mechanism as the project's terminal structure

User framing: "we should figure out how to include fields. This will
be the structure this whole project ends on."

Interpretation: once the classifications and ensemble summaries
stabilise, the project terminates on a **field manifest** — a
declarative, typed, versioned list of every field the H5 carries,
with provenance, units, decoder class, and inclusion policy. Every
downstream consumer (Stage 1, Stage 2, Stage 3, viewer, reader,
advisor's future students, collaborators) reads that manifest. The
manifest IS the interface.

**Shape (tentative):**

- Each `ConformationResult` declares its fields via a declarative
  manifest — not scattered `createDataSet` calls inside a `WriteH5`
  method.
- Field metadata per entry: name, units, H5 path, dtype, shape,
  scope (atom / residue / ring / frame / system), primary modality
  (per the h5-reader field glossary), provenance (calculator +
  file:line + commit SHA), description, inclusion policy.
- TOML-driven inclusion / exclusion: some fields ship always, some
  optional, some debug-only.
- The `h5-reader/notes/H5_FIELD_GLOSSARY.md` auto-generates from the
  manifest — glossary becomes a view over the manifest, not a
  parallel document that drifts.
- SDK `_catalog.py` also generates from the manifest — one source of
  truth for what the H5 carries.

**Relationship to 13.6 and 13.7:** if we adopt the field-manifest
pattern, it naturally absorbs the validation decorators from 13.6
(per-field constraints) and the format-test framework from 13.7
(per-field roundtrip and boundary tests). This is why the user
described it as "the structure this whole project ends on" — it's
not a schema change but a schema-description-language change that
makes schema + validation + tests + glossary + SDK all downstream
artefacts of one typed specification.

**Worth seriously considering as the rollup's eventual form.** The
rollup as currently scoped in sections 4–5 is a set of new fields;
13.8 is the declarative framework those fields would live inside.
The question is whether the rollup ships the fields and then 13.8
comes after, or whether the rollup is rescoped around 13.8 and
ships both the framework and the fields together. Open.

### 13.9 — Deep question: serialise the original model, or continue parallel `QtProtein` types

User's articulated question: does it make more sense to continue the
H5 independently in the `QtProtein` form, or should we think in terms
of making the original model serialisable? Explicitly unresolved;
captured here so the question stays visible until it is answered.

**Option A — continue parallel `QtProtein` types (status quo):**

- `src/` types (`Protein`, `ProteinConformation`, `Atom`, `Ring`,
  `Bond`) are library-internal, not directly serialisable.
- `fileformat/` is the bridge: writes H5 from library types, reads
  H5 into `QtProtein`/`QtConformation`/`QtFrame` (h5-reader's
  Qt-native hierarchy).
- Parallel type hierarchies, coupled by schema only.

Pros: viewers never link the library, no OpenBabel / cifpp /
libgromacs / APBS / MOPAC / libdssp dependencies on the viewer side;
cross-platform viewer build stays tractable; schema is the stable
contract independent of type churn; viewers cannot be tempted to
"just recompute that quickly."

Cons: every new field costs double — library type + H5 writer +
Qt type + Qt reader. This rollup doubles a large surface. Drift
risk between the two sides grows as the schema widens. The work
the user is paying up-front in the H5 gets partially re-paid on
the viewer-decode side.

**Option B — split-library, serialisable original model:**

- Refactor `src/` into `src/model/` (lean: `Protein`,
  `ProteinConformation`, `Atom`, `Bond`, `Ring` hierarchy,
  `SphericalTensor`, `ConformationAtom`, only Eigen / HighFive /
  sphericart dependencies) and `src/builders/` + `src/calculators/`
  (heavy stack: OpenBabel, cifpp, libgromacs, MOPAC subprocess,
  APBS binding, libdssp).
- `nmr_shielding_model` target: lean, serialisable.
- `nmr_shielding_extract` target: links everything, is the heavy
  extractor path.
- Viewers (`ui/`, `h5-reader/`) link `nmr_shielding_model` only.
  Get typed access to real model types. `QtProtein`/`QtConformation`/
  `QtFrame` disappear as a parallel hierarchy, possibly replaced by
  thin Q_OBJECT adapters over the real types for signal / slot
  wiring.

Pros: single type hierarchy; schema follows types automatically;
zero drift risk — any type change is compile-time visible across
library and viewers; every new field is declared exactly once.
Viewer feature work reads typed properties directly without a
schema-to-type mapping layer.

Cons: substantial refactor — moving files, splitting CMake
targets, verifying model / builder / calculator separations are
genuinely clean. Viewers take a compile-time dependency on model
headers plus Eigen / HighFive / sphericart; viewer builds get
heavier but not prohibitively so. If any sloppy include between
model and builders / calculators exists, the split blocks until
it's fixed; audit cost is real.

**Option C — code-generated bridge (middle path):**

- Single schema definition (the field manifest of 13.8, or
  equivalent).
- Codegen produces both a C++ serialiser for the library and a
  Qt-native reader for the viewer.
- Type hierarchies stay separate but the schema contract is
  guaranteed consistent by the generator.

Pros: eliminates drift without the refactor; viewers stay free of
library deps.

Cons: adds build-time tooling; the generator itself is a maintenance
item. Types still diverge — only the serialisation glue is unified.

**Drawability audit (completed 2026-04-22 pm):**

Run against the current `src/` tree.

- **All model headers are lean.** `Protein.h`, `ProteinConformation.h`,
  `ConformationAtom.h`, `Ring.h`, `Bond.h`, `Atom.h`, `Residue.h`,
  `Types.h` (which contains `SphericalTensor` inline),
  `ConformationResult.h`, `PhysicalConstants.h`, `AminoAcidType.h`,
  `IupacAtomIdentity.h`, `CovalentTopology.h` itself, `GeometryChoice.h`,
  `ProteinBuildContext.h`, `CalculatorConfig.h` — only Eigen + stdlib +
  local model headers. No OpenBabel, cifpp, libgromacs, libdssp,
  MOPAC, or APBS in any header form.
- **All model `.cpp` files are lean except one.** `Protein.cpp`,
  `ProteinConformation.cpp`, `Ring.cpp`, `Atom.cpp`, `AminoAcidType.cpp`,
  `Types.cpp`, `IupacAtomIdentity.cpp`, `NamingRegistry.cpp`,
  `ConformationResult.cpp` (base class) — all stdlib + local.
  `CovalentTopology.cpp` is the single heavy coupling: it pulls
  OpenBabel for bond detection. Since bond detection runs only during
  `FinalizeConstruction` at extract time and viewers read pre-built
  topology from H5, `CovalentTopology.cpp` stays in the extract target
  naturally.
- **Support files lean.** `OperationLog.h`, `NpyWriter.h` — stdlib
  only.
- **`fileformat/analysis_file.h` already exists as a third type
  hierarchy.** Plain-flat-vector HDF5 schema carrier, distinct from
  library model types and from `QtProtein`. Currently source-included
  by `h5-reader/`.
- **h5-reader currently links HDF5 + Eigen + Qt + VTK only**, no
  library. Adding a new lean `nmr_shielding_model` target is additive.

**Verdict: Option B is drawable. B2 is the chosen sub-option** under
the user's "force everything through again, no backward compatibility"
stance (2026-04-22 pm). With no old files to translate, there is no
reason to preserve `fileformat/` as an intermediate layer — model
types become directly HDF5-serialisable and the old `AnalysisFile`
plain-struct layer moves to `bones/` as history.

**Target structure under B2:**

- **`nmr_shielding_model`** (new lean target): every model header
  and its `.cpp` (except `CovalentTopology.cpp`), `ConformationResult`
  base class, `AminoAcidType`, `IupacAtomIdentity`, `NamingRegistry`,
  `PhysicalConstants`, `GeometryChoice`, `CalculatorConfig`,
  `OperationLog`, `NpyWriter`, plus the new direct-HDF5 read/write
  methods. Dependencies: Eigen, sphericart, HighFive. Nothing else.
- **`nmr_shielding_extract`** (heavy): `CovalentTopology.cpp`, all
  builders (`PdbFileReader`, `GromacsEnsembleLoader`, `OrcaRunLoader`,
  `FullSystemReader`), all concrete calculators (`MopacResult`,
  `ApbsFieldResult`, `DsspResult`, ring/bond calculators, AIMNet2,
  SASA, EEQ, …), `AnalysisWriter`. Links `nmr_shielding_model`.
  Dependencies: OpenBabel, cifpp, libgromacs, libdssp, APBS, MOPAC
  subprocess — the full stack.
- **Viewers** (`ui/`, `h5-reader/`): link `nmr_shielding_model` only.
  Get typed access to real `Protein`, `ProteinConformation`, all
  enums, all model properties. `QtProtein` / `QtConformation` /
  `QtFrame` collapse — replaced by thin `Q_OBJECT` adapters over the
  real types where signal / slot wiring is needed.

Sequencing of the refactor vs. the feature work is set out in 13.10.

### 13.10 — Two-milestone sequencing (decided 2026-04-22 pm)

The rollup splits into two ordered milestones rather than one bundled
commit window. Decided after the 13.9 audit landed "B2 is drawable"
and the user observed that (a) we know the new feature work is needed
regardless, (b) the assistant's time estimates are typically
conservative but the risk of a bundled refactor is real, and (c)
shipping features against the known-working architecture first and
refactoring second absorbs risk rather than compressing it.

**Milestone 1 — features in the current architecture.**

All new feature work lands in the current library + `fileformat/` +
`QtProtein` architecture. Concretely:

- Categorical updates: `NmrAtomIdentity` and its enum lenses per
  section 4, added to `EnrichmentResult` and to `/atoms/` in the
  existing H5 schema.
- `NamingRegistry` fixes per section 9 and
  `ChangesRequiredBeforeProductionH5Run.md` — activated.
- Continuous / difference / variance recording per section 5, folded
  into the existing `GromacsProtein` two-pass machinery. The Welford
  accumulator pattern already in place extends from scalars to
  tensors, autocorrelations, covariance matrices, density grids.
  No new object-model category introduced in M1 — `GromacsProtein`
  absorbs the trajectory-summary discipline.
- Event scheme per sections 6–7, also folded into `GromacsProtein`'s
  frame-processing loop. The threshold-parameterised menu ships as
  multiple typed event arrays under `/derived_events/` in the
  existing H5.
- The hookable extractor pattern (section 7) ships its C++ side in M1
  and its Python canonical-form counterpart in the SDK, both against
  the current `fileformat/` schema.
- Validation infrastructure per 13.6 — hand-derived fixtures for all
  20 canonical AAs, external reference comparison against PDB CCD /
  IUPAC / BMRB, internal consistency checks, pro-R / pro-S geometric
  cross-check at enrichment.
- Both viewers (`ui/`, `h5-reader/`) updated to decode the new fields
  through the existing Qt-native read path.

**M1 deliverable:** 10 calibration H5s regenerated with the complete
new field surface in the current three-hierarchy shape. All tests
pass. The H5 at this point carries the full data surface the
project's analysis phase needs — classification correct, ensemble
summaries populated, events extractable via the hookable menu,
viewers reading typed fields. **This is the checkpoint where the
*data* is right.**

**Milestone 2 — serialisation architecture refactor (B2).**

With M1's feature surface stable and tested, the architectural
refactor runs as a separate pass. No new feature work in M2.

- Split `src/` into `src/model/` (lean) and `src/extract/` (heavy)
  CMake targets per the drawability audit in 13.9.
- Collapse `fileformat/` into the model. Model types become directly
  HDF5-serialisable (`Protein::ReadH5` / `Protein::WriteH5` or
  equivalents). The existing `AnalysisFile` plain-struct layer moves
  to `bones/` as history.
- Collapse `QtProtein` / `QtConformation` / `QtFrame` — viewers link
  `nmr_shielding_model` and operate on real model types, with thin
  `Q_OBJECT` adapters where signal / slot wiring requires them.
- Field-manifest pattern per 13.8 becomes natural to implement at M2
  boundary (one declarative manifest drives the library serialiser,
  the SDK catalog, and the viewer glossary). Worth considering as
  part of M2's scope.

**M2 deliverable:** the same 10 calibration H5s regenerate from the
refactored library, bit-exact to M1 output. Both viewers decode the
refactored model types. `QtProtein` hierarchy deleted. The old
`fileformat/` module retired to `bones/`. Library is coherent.
**This is the checkpoint where the *architecture* is right.**

**Rationale for the split:**

- M1 validates feature correctness in a well-understood architecture.
  Any classification bug, ensemble-statistic bug, or event-extractor
  bug surfaces against the known `fileformat/` + `QtProtein` shape
  rather than tangled with an architectural change.
- M2 moves validated feature code to a new architecture. Only the
  structural changes are untested at M2 boundaries; the feature
  logic is already exercised against real calibration data.
- If M2 slips or surfaces surprises, M1 has already shipped value —
  classifications and summaries are in consumers' hands, and the
  architecture refactor lands when it's ready.
- If M2 succeeds bit-exactly, the test is clean — M1 is the oracle.

**What this changes elsewhere in the spec:**

- Scope boundary (section 2): reframe from "one coherent rollup" to
  "one coherent M1 rollup, followed by an M2 architectural pass."
- Regeneration (section 11): one M1 regeneration round (the 10
  calibration H5s). M2 produces a second regeneration pass that must
  reproduce M1's output bit-exactly; treated as a validation test,
  not a new schema.
- Activation criteria (section 9): NamingRegistry activation is part
  of M1, not gated on M2.
- Scope mapping: sections 4, 5, 6, 7, 9, 13.6 are M1 scope. Sections
  13.7, 13.8, 13.9 are M2 scope.

**Assistant-estimate caveat (user note, 2026-04-22):** the
assistant's time estimates in prior sections tend to overstate. Risk
is nonetheless real; sequencing absorbs that risk rather than
compressing it. Do not read the two-milestone framing as "twice as
slow"; read it as "each milestone has a clean acceptance test the
other does not depend on."

### 13.11 — Physics foundations document comes *before* M1 (2-3 sessions)

Decided 2026-04-22 pm, late-session. Inserted ahead of M1 in the
sequencing of section 13.10. The next 2–3 sessions' focus is not
implementation; it is building a proper referenced physics model,
grounded in citable reality, for each component of shielding and
molecular physics at the protein scale and below. The thesis cannot
stand on "the assistant's training plus a handful of references" —
the extractor has demonstrated the signal exists and relates to
r2SCAN shielding (Stage 1 R² = 0.818), and what has to happen before
analysis is that the physical meaning of every signal the H5 carries
is anchored in primary literature, not inferred.

**Target document:** `spec/PHYSICS_FOUNDATIONS.md` (to be created).
Tier-1 reading when stable. Sacred-doc adjacent — constitution
dictates what must be true; physics foundations explains *why* in
citable form.

**Shape (tentative, for session-1 opening to refine):**

1. **Shielding fundamentals.** σ tensor definition, irreducible
   decomposition (T0+T1+T2), sign conventions, chemical shift vs
   shielding, time-averaging and motional narrowing, DFT shielding
   via GIAO, r2SCAN functional, the WT-ALA delta as aromatic
   isolator. References: Ramsey, Facelli, Ditchfield, Furness et al.
   2020.
2. **Ring current physics.** Aromatic π-electron origin (Pauling,
   London), quantum treatment (Pople 1956), Johnson-Bovey classical
   loop, Haigh-Mallion surface integral, Case calibration. Modern
   reviews (Mulder, Lampert).
3. **Bond anisotropy (McConnell).** Magnetic susceptibility
   anisotropy, full asymmetric tensor, Δχ per bond category. McConnell
   1957, ApSimon-Craig 1967, Bothner-By.
4. **Electric field effects (Buckingham).** E and E² terms, A/B
   coefficients, Coulomb vacuum vs APBS solvated, PM7 semiempirical
   polarisation. Buckingham 1960, Raynes 1962, Baker 2001 (APBS),
   Stewart 2013 (PM7).
5. **Aromatic quadrupole.** Stone T-tensor, π-cloud axial quadrupole,
   EFG interpretation for quadrupolar nuclei (²H, ¹⁴N, ¹⁷O). Stone,
   Pyykkö.
6. **H-bond shielding.** Cornilescu-Bax formalism; relation to bond
   anisotropy (same tensor, different axis). Cornilescu-Bax 1999.
7. **Dispersion.** London 1/r⁶, CHARMM switching. London 1937,
   Brooks 1983.
8. **Time-resolved NMR — the MD-derived payload.** Lipari-Szabo S²,
   spectral densities and R1/R2/NOE, chemical exchange (CEST, CPMG),
   RDC in aligned media. Lipari-Szabo 1982, Palmer reviews,
   Tjandra-Bax. This chapter is the theoretical backing for the
   TrajectoryResult fields in section 5.
9. **Proxies vs first principles.** Where we use proxies honestly
   (classical Coulomb as E-field proxy, point-dipole approximations
   in ring susceptibility), what rigorous first-principles would
   require, and the explicit scope decision: understand classical
   contributors rigorously and honestly; first-principles QM beyond
   the classical kernels is outside scope (dft-ex1's territory).
10. **Connections to code and H5.** Each calculator mapped to its
    physics section; each TrajectoryResult field mapped to its
    time-resolved-physics section; each signal anchored to a
    reference.

**Relationship to existing documents:**

- `GEOMETRIC_KERNEL_CATALOGUE.md` is *what* we compute, equation as
  it appears in code. Physics foundations is the *why*.
- `MATHS_GOALS.md` is the validation plan. Physics foundations
  justifies the plan.
- `CALCULATOR_PARAMETER_API.md` lists 93 parameters with equations
  and sparse references. Physics foundations puts those references
  in their full theoretical context.
- `references/ANNOTATED_BIBLIOGRAPHY.md` is the paper-level index.
  Physics foundations cites into it; every citation in the physics
  doc has a PDF in `references/` and an entry in the bibliography.

**Session plan (4 × ~4-hour sessions) — broad-search pass then three
drafting sessions:**

- **Session 0 — Broad-search / landscape pass.** Added late
  2026-04-22 pm per user's research philosophy: "ask the biggest
  questions you are addressing and see what people have to say, even
  if it doesn't have immediate application to what you are doing.
  That way we have checked out the haystack for needles. Then we can
  come back and build an entire little mini-castle of needles and
  pins, and know we are not leaving out things that people will see
  and note." Our prior literature has been focused — which has been
  appropriate under directed work but means we haven't swept the
  broader territory. Before drafting, we sweep.
  Concrete queries (no NMR qualifier):
  - "Contributors to diamagnetic shielding" at molecular scale
  - "Diamagnetic effects in proteins" (bulk susceptibility, MRI
    contrast origins, anisotropic response)
  - "Paramagnetism in biomolecules" (metal centers, radicals)
  - "Classical predictors of magnetic properties at molecular
    scale" (any field, any era)
  - "Time-correlation functions in molecular simulation"
    (dielectric relaxation, dipole autocorrelation, Raman line
    shape, not just Lipari-Szabo)
  - "GIAO and alternatives" / "gauge-origin problem in magnetic
    property calculation" (CSGT, IGLO, LORG, recent comparisons)
  - "MD-derived molecular property calculation" (general
    compute-on-frames-average across chemistry)
  - "Semiempirical methods for magnetic properties" (PM7/MOZYME
    suitability beyond charges)
  - "Ensemble averaging of computed molecular observables"
    (stat-mech view)
  - "Chemical shielding anisotropy" as a field
  - "Paramagnetic relaxation enhancement" / "pseudocontact shifts"
    (orbit-same-physics reference)
  Output: annotated landscape, PDFs fetched to `references/`,
  `ANNOTATED_BIBLIOGRAPHY.md` updated. Not a draft — a scaffold
  that the three drafting sessions build on. Method rationale
  preserved in this section so future sessions know the method,
  not just the queries.
- **Session 1 — Fundamentals + ring currents (parts 1, 2).** The
  core aromatic signal. Much of this is already covered in existing
  references per user note 2026-04-22 pm — "work is required but we
  know what it is." Literature canonical list: Ramsey, Pople 1956,
  Johnson-Bovey 1958, Haigh-Mallion 1973/1979, Case 1995/1999,
  Mulder, Furness et al. 2020 (r2SCAN). Extend with relevant finds
  from Session 0. Update `ANNOTATED_BIBLIOGRAPHY.md`. Write parts
  1 and 2 of `PHYSICS_FOUNDATIONS.md`.
- **Session 2 — Electric, bond, H-bond, dispersion, quadrupole
  (parts 3–7).** Partially covered by existing references per user
  note 2026-04-22 pm. Canonical list: Buckingham 1960, Raynes,
  McConnell 1957, ApSimon-Craig, Cornilescu-Bax 1999, Stone
  T-tensor, London, Brooks. Extend with Session-0 finds. Write
  parts 3–7.
- **Session 3 — Time-resolved + proxies + connections (parts 8–10).**
  Heaviest fresh literature pass per user's "session 3 is what needs
  the work." Canonical list: Lipari-Szabo 1982, Palmer reviews,
  Tjandra-Bax, advisor-flagged RDC references. Extend substantially
  with Session-0 finds on time-correlation functions, MD-derived
  properties, ensemble averaging. Read `dft-ex1` project for
  scope-boundary accuracy in part 9. Write parts 8–10 and cross-map
  to every H5 field.

**Between sessions:**

- Advisor's preferred reviews — especially for the time-resolved
  chapter — flag for the user to surface before session 3.
- User's class notes that should be primary sources rather than the
  assistant reaching for standard textbook references.
- The `dft-ex1` project pointer for scope-boundary clarity in part 9.

**Revised sequencing, master order:**

1. **Physics foundations** (next 2–3 sessions). Target doc
   `spec/PHYSICS_FOUNDATIONS.md`.
2. **M1 — features in current architecture** (section 13.10).
3. **M2 — B2 serialisation refactor** (section 13.10).

The rollup implementation work does not begin until physics
foundations exists in referenced form. The reason: every enum name,
every ensemble field, every event extractor choice in M1 should have
a physics-foundations anchor by the time it lands in code. Building
M1 on "the assistant's informal reading" reintroduces the exact
problem the physics document is solving.

Planned sequence, not commitments:

1. **Literature research — joint.** User's framing: "our actual power
   is in these 260 DFTs, and our novel application is a pretty good
   proxy for shielding from maths." Research scope should think of
   the project as quantum-mechanics + kernel-geometry as much as (or
   more than) NMR-experimentalist tradition. Advisor's work is on
   real expressed proteins — one input, not the whole frame. Bring
   literature home for dynamics representations (spectral density,
   RDC, time-series ML for MD), for NMR invariant classification,
   and for what the 260-DFT novel-proxy angle brings that
   BMRB-calibrators don't have.

2. **Draft the exact `NmrClass` / `Locant` / `RingAtomRole` enum
   contents.** Hand-derive against all 20 canonical amino acids as
   fixtures. This is the cheapest deliverable that hard-tests the
   spec.

3. **Draft the Stage 1 strata reproduction mapping.** Verify the new
   classification can reproduce every AMBER-atom-name stratum
   `learn/src/mutation_set/` uses today.

4. **Python prototype of NmrAtomIdentity against the 10 existing
   H5s.** Before library commit. Confirm forensics binding,
   stratification, and anticipated Stage 2 analyses can be expressed
   in the enum set. No regeneration needed for this step — classify
   in Python off the existing H5s.

5. **Spec review pass with the advisor.** Especially sections 4 and
   12 — the invariant classifier contents and the NMR-specific
   continuous-dynamics unknowns.

6. **Schedule commit.** If the above converge, the rollup lands as
   one commit (or a few atomic commits landed together), the 10 H5s
   regenerate in one pass, and production fleet extraction activates.

Nothing in this section is an instruction to execute. The spec is
captured so the thinking persists across sessions; the execution
sequence is conversation-driven.

---

## 15 — This document's status

**Tentative. Living. Not gospel.**

- It captures 2026-04-22 design conversation between user and
  assistant, no advisor review yet, no implementation work yet.
- Sections 4 (enum contents), 5 (ensemble field list), 6 (event
  menu), 9 (NamingRegistry list), 12 (open questions), 13
  (pm addendum — including 13.6 validation, 13.7 universal format,
  13.8 field manifest, 13.9 serialisation question, 13.10
  two-milestone sequencing) are all expected to churn. Section 13
  in particular is conversation-in-progress and should be read as
  discussion, not as a priori principle.
- 13.9 was initially open; closed 2026-04-22 pm at "B2 is drawable
  and is the right target." 13.10 sequences feature work before
  refactor (M1: features on current architecture, M2: B2 refactor).
  Both decisions are subject to revision but now have a concrete
  landing place in the spec rather than an open question.
- 13.11 (added late 2026-04-22 pm) inserts a physics-foundations
  document ahead of M1. The next 2–3 sessions' focus is writing
  `spec/PHYSICS_FOUNDATIONS.md` — a properly-referenced theory
  chapter for every mechanism the extractor touches. Rollup
  implementation does not begin until that exists.
- The principles in section 3 are the most stable — they are the
  distillation of multiple conversation turns and explicit user
  corrections (specifically: "events become competing identities";
  "raw signal stays primary"; "Python loses granularity so pay cost
  in H5 up front"; "richer than we need is fine").
- Regeneration cost assumptions in section 11 are verified against
  source — positions come from MD, DFTs align by element — but the
  *exact* time estimates are wall-clock guesses.
- Any future session editing this spec should preserve the
  tentative-and-dated framing. Do not convert this into a
  sacred-doc-style specification without explicit user sign-off.
