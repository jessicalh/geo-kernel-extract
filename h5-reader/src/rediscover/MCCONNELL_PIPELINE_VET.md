# McConnell Pipeline Vet

Scope: read-only vet of calculator -> emit -> consumption. No code fixes,
no merge, no re-emit, no ORCA, and no `trajectory.h5` access were used for
this report.

## Executive Verdict

The leading explanation is a source-selection bug/mismatch in rediscover:
the broad-backbone McConnell path includes the target atom's own anisotropic
bonds, and also includes directly bonded/local bonds, while the producer-side
calculator excludes at least endpoint self-sources and applies a dipolar
near-field filter. For backbone N/C/O, those own amide-plane bonds are exactly
the bonds McConnell is modelling. At 0.5-0.7 A from a bond midpoint, the
far-field point-dipole tensor is not a valid approximation and can dominate the
sum. This cleanly matches the symptom: C becomes robust but sign-flipped
because C is an endpoint of both C=O and C-N sources.

This is not a repeat of the F1 false alarm. The project-canonical McConnell
tensor is sound: `FIXES_AUDIT_opus.md:22` and `FIXES_AUDIT_opus.md:83-91`
explicitly resolved that formula as the canonical `K*chi_aniso` tensor, not a
home-rolled tensor. The problem here is which source bonds rediscover feeds
into that tensor.

Second-tier issues are real but less explanatory: Stage 1 and rediscover test
different quantities, and rediscover's default cutoff is 8 A while the
calculator default is 10 A. Frame, T2 basis, current units/sign, and current
deweighting are clean.

## Calculator -> Emit Trace

Producer McConnell uses the canonical full tensor:
`../src/McConnellResult.cpp:54-98` computes the scalar `f`, the dipole kernel
`K`, and the full McConnell tensor `M_over_r3`. The MOPAC/bond-order weighted
variant mirrors this at `../src/MopacMcConnellResult.cpp:47-89` and applies
the Wiberg bond order at `../src/MopacMcConnellResult.cpp:184-187`.

Producer source filtering is not "all nearby bond midpoints":
`../src/McConnellResult.cpp:120-125` and
`../src/MopacMcConnellResult.cpp:115-118` build a filter set containing
`MinDistanceFilter`, `SelfSourceFilter`, and `DipolarNearFieldFilter`.
`SelfSourceFilter` rejects target atoms that are either source endpoint
(`../src/KernelEvaluationFilter.h:216-241`). `DipolarNearFieldFilter` rejects
sources unless distance exceeds `source_extent * near_field_exclusion_ratio`
(`../src/KernelEvaluationFilter.h:174-197`). The default ratio is 0.5
(`../src/CalculatorConfig.cpp:70`).

Producer cutoff is configuration-backed, default 10 A:
`../src/CalculatorConfig.cpp:50-51`. The header constants in
`../src/McConnellResult.h:36-39` and `../src/MopacMcConnellResult.h:35-38`
are vestigial comments pointing to the runtime config.

Producer category accumulation is per atom over filtered bonds:
`../src/McConnellResult.cpp:139-183` searches nearby bonds, fills filter
context with atom index and source endpoints, and rejects failed filters.
`../src/McConnellResult.cpp:196-241` accumulates by bond category. Category T2
outputs are symmetrized and trace-removed at `../src/McConnellResult.cpp:269-281`;
the total packed shielding contribution is `SphericalTensor::Decompose(M_total)`
at `../src/McConnellResult.cpp:283-284`. It writes `mc_shielding.npy`,
`mc_category_T2.npy`, and `mc_scalars.npy` at
`../src/McConnellResult.cpp:360-395`. The MOPAC path is analogous at
`../src/MopacMcConnellResult.cpp:130-228`, `../src/MopacMcConnellResult.cpp:257-272`,
and `../src/MopacMcConnellResult.cpp:300-337`.

Trajectory H5 emission carries the bare geometric kernel, not calibrated ppm:
`../src/McConnellShieldingTimeSeriesTrajectoryResult.h:12-15` and
`../src/MopacMcConnellShieldingTimeSeriesTrajectoryResult.h:3-15` document
Angstrom^-3, pre-Delta-chi, full tensor payloads.

## Rediscover -> Consumption Trace

Rediscover's broad-backbone path constructs local frames for N/CA/C/O/HN/HA
at `src/rediscover/BroadBackbone.cpp:159-255`. It selects every backbone atom
or alpha hydrogen at `src/rediscover/BroadBackbone.cpp:424-430` and nearby
anisotropic bond midpoint sources at `src/rediscover/BroadBackbone.cpp:437-440`.
Only PeptideCO, PeptideCN, SidechainCO, and Aromatic bonds are indexed into
the bond cloud (`src/rediscover/SpatialIndexSet.cpp:19-29`,
`src/rediscover/SpatialIndexSet.cpp:92-102`).

The broad bond attacher computes source rows at
`src/rediscover/BroadBackbone.cpp:310-359`. It uses target-to-midpoint
displacement and rejects only `r <= 1e-6` at
`src/rediscover/BroadBackbone.cpp:342-345`. That does not reject endpoint
self bonds, because an endpoint atom sits roughly half a bond length from the
bond midpoint. The reducer then sums every finite bond source:
`src/rediscover/BroadBackboneSink.h:85-130`,
`src/rediscover/BroadBackboneSink.cpp:52-68`, and
`src/rediscover/BroadBackboneSink.cpp:249-251`.

The broad fallback `bondKernelT2FromSources` path has the same all-source
semantics: it loops over `SourceKind::Bond`, reconstructs the project-canonical
McConnell tensor from each source, and sums into one local T2 with no
self/bonded/near-field filter (`src/rediscover/BroadBackbone.cpp:80-106`).
`MakeBroadBackboneRelationship` uses the H5 bare kernel if present, otherwise
falls back to this function (`src/rediscover/BroadBackbone.cpp:467-483`);
the source-level `mc_lit_T2_local_*` columns used by the current analyses are
still source-summed from the all-source CSV path.

The standalone/composed McConnell path has the same problem:
`src/rediscover/McConnellNeighborhood.cpp:163-180` and
`src/rediscover/ComposedRelationships.cpp:225-274` reject only degenerate
midpoint/axis cases, then `src/rediscover/McConnellNeighborhood.cpp:213-215`
and `src/rediscover/ComposedRelationships.cpp:276-290` explicitly treat all
bond sources as valid. The comments at
`src/rediscover/McConnellNeighborhood.cpp:102-104` and
`src/rediscover/ComposedRelationships.cpp:276-290` say there is no self/bonded
concept for bonds.

The contrast with ring is explicit: ring constructs own-ring/own-atom sets and
excludes `is_self_or_bonded` from valid sums at
`src/rediscover/ComposedRelationships.cpp:132-177`. Producer ring also has a
bonded-neighbour exclusion (`../src/KernelEvaluationFilter.h:318-347`,
`../src/KernelEvaluationFilter.cpp:22-40`).

Current consumption reads emitted CSV/NPY sidecars, not the H5 kernel directly.
`src/rediscover/analysis/mcconnell_literature_decirc.py:590-626` reads
`broad_backbone_aggregated.csv`, `mc_lit_T2_local_*`, and
`broad_backbone_aggregated_target_local_T2.npy`. The Delta-chi calibration
source-sums emitted source CSV rows by `row_id` and category at
`src/rediscover/analysis/mcconnell_dchi_calibration.py:200-245`, deweights by
the same scale at `src/rediscover/analysis/mcconnell_dchi_calibration.py:248-259`,
and fits flattened T2 components at
`src/rediscover/analysis/mcconnell_dchi_calibration.py:301-316`.

## Suspect Verdicts

### Suspect 1: Self/Bonded Exclusion

Verdict: BUG/MISMATCH, lead.

Producer excludes endpoint self-sources and near-field sources:
`../src/McConnellResult.cpp:120-125`,
`../src/McConnellResult.cpp:165-183`,
`../src/KernelEvaluationFilter.h:174-197`, and
`../src/KernelEvaluationFilter.h:216-241`. Rediscover includes all finite
bond midpoint sources:
`src/rediscover/BroadBackbone.cpp:342-345`,
`src/rediscover/BroadBackboneSink.h:106-110`,
`src/rediscover/BroadBackboneSink.cpp:52-68`,
`src/rediscover/McConnellNeighborhood.cpp:163-180`, and
`src/rediscover/ComposedRelationships.cpp:276-290`.

Backbone stratum source inclusion in the current broad path:

- N: includes the atom's own previous peptide C-N (PeptideCN) when present.
  It also includes local C=O on the same amide plane when within cutoff, even
  though that source is bonded through the carbonyl C to N. N-CA is not in the
  anisotropic bond cloud because BackboneOther is not indexed.
- C: includes the atom's own C=O (PeptideCO) and its own next peptide C-N
  (PeptideCN) when present. This is the most direct explanation for the robust
  sign-flipped C stratum.
- O: includes the atom's own C=O (PeptideCO), and can include the adjacent
  C-N through the carbonyl C.
- CA/HA: no own anisotropic CA bond is in the bond cloud, but adjacent C=O and
  peptide C-N bonds bonded through CA-C or CA-N are included.
- HN: N-H itself is not an anisotropic source, but the peptide C-N containing
  the bonded N is included. Whether this should be excluded depends on the
  intended "bonded" definition; it is at least local through-bond rather than
  a remote anisotropic source.

This means the rediscover geometric T2 is not the same quantity as the
producer's filtered McConnell contribution for backbone heavy atoms. It is an
all-local-plus-remote bond-midpoint sum.

### Suspect 2: Frame Consistency

Verdict: CLEAN.

Broad local frames are built per stratum at
`src/rediscover/BroadBackbone.cpp:159-255`. Source vectors are stored in that
same local frame, and `src/rediscover/McConnellLiteratureKernel.cpp:50-88`
builds the McConnell tensor from the local source axis/displacement. Target
DFT tensors are rotated with `LocalFrame::TensorToLocal` (`src/rediscover/LocalFrameBasis.h:65-85`)
and emitted as `broad_backbone_aggregated_target_local_T2.npy`
(`src/rediscover/BroadBackboneSink.cpp:274-289`). The analysis scripts consume
that local target sidecar (`src/rediscover/analysis/mcconnell_literature_decirc.py:590-626`).

The extractor also performs an ORCA/H5 frame-alignment check before building
rediscover rows (`src/rediscover/main_extract.cpp:262-270`). I found no
EFG-style lab-frame/tumbling mismatch in the current McConnell decirc scripts.

### Suspect 3: T2 Basis / Traceless Convention

Verdict: CLEAN.

Producer decomposition order is `[xy,yz,zz,xz,xx-yy]` in
`../src/Types.cpp:25-60`, with reconstruction at `../src/Types.cpp:63-94`.
Rediscover's `DecomposeLibrary` is the same library-order clone at
`src/rediscover/SphericalBasis.cpp:7-37`, and `src/rediscover/SphericalBasis.h:1-16`
documents why that exact basis is used.

The current McConnell literature tensor trace-removes the PCS tensor before
decomposition (`src/rediscover/McConnellLiteratureKernel.cpp:77-84`), and the
target local sidecar is decomposed with the same library basis
(`src/rediscover/BroadBackboneSink.cpp:274-289`). The current decirc report's
T0 audit shows structural zeros for McConnell T0
(`src/rediscover/analysis/MCCONNELL_LITERATURE_DECIRC.md:34-45`). No basis or
normalisation mismatch is evident.

### Suspect 4: Stage 1 vs Rediscover Same Quantity

Verdict: MISMATCH / DIFFERENT TEST, not the lead bug.

The exact doc path named by the repo, `learn/stage1-mutations/`, is absent in
this checkout. The nearest current docs still define Stage 1 unambiguously:
`/shared/2026Thesis/nmr-shielding/CLAUDE.md:70-74` says Stage 1 is a settled
per-element, per-atom-type ridge on 720 proteins / 446K atoms / 55 kernels,
R2 = 0.818. `/shared/2026Thesis/nmr-shielding/CLAUDE.md:42-44` defines mutant
mode as WT+ALA with WT-ALA delta tensors, and
`/shared/2026Thesis/nmr-shielding/README.md:16-17` says calibration is against
DFT WT-ALA deltas across 720 proteins. The emitted target is WT minus mutant
DFT total shielding (`../src/MutationDeltaResult.cpp:458-460`) and is written
as `delta_shielding.npy` (`../src/MutationDeltaResult.cpp:681-687`), with
`delta_scalars` carrying `delta_T0` (`../src/MutationDeltaResult.cpp:720-728`).

The SDK catalog marks the static McConnell arrays as model features:
`/shared/2026Thesis/nmr-shielding/python/nmr_extract/_catalog.py:188-194`.
Those are `mc_shielding` (full spherical tensor), `mc_category_T2`, and
`mc_scalars`, all mechanism `bond_anisotropy`. The same catalog marks
`delta_shielding` as WT-ALA shielding delta
(`/shared/2026Thesis/nmr-shielding/python/nmr_extract/_catalog.py:374-379`).

Rediscover is different: it tests one 1P9J trajectory's per-frame local T2
against per-frame DFT T2 rows
(`/shared/2026Thesis/nmr-shielding/CLAUDE.md:75-79`,
`src/rediscover/analysis/mcconnell_literature_decirc.py:590-626`). A mechanism
can carry Stage-1 scalar/tensor mutation-delta signal and still fail a
one-protein, per-frame, local-T2 Delta-chi fit. That said, the self/near-field
source mismatch above is still a concrete pipeline bug and should be fixed
before interpreting the rediscover non-convergence physically.

### Suspect 5: Bond Categorization / Cutoff / Double Count

Verdict: MIXED.

Categorization is mostly clean. Topology assigns backbone C-O as PeptideCO,
cross-residue C-N as PeptideCN, other backbone as BackboneOther, and non-BB
C-O below 1.35 A as SidechainCO (`../src/CovalentTopology.cpp:153-194`).
Aromatic ring bonds are retagged as Aromatic after ring construction
(`../src/CovalentTopology.cpp:247-268`; ring construction call at
`../src/Protein.cpp:512-518`). Rediscover indexes exactly PeptideCO,
PeptideCN, SidechainCO, and Aromatic (`src/rediscover/SpatialIndexSet.cpp:19-29`).

Aromatic is correctly de-circularised out of McConnell: rediscover assigns
q = 0.0 for Aromatic at `src/rediscover/McConnellLiteratureKernel.cpp:31-44`,
and calibration excludes category 4 at
`src/rediscover/analysis/mcconnell_dchi_calibration.py:28-43`. The current
calibration report confirms aromatic rows have zero absolute T2 contribution
(`src/rediscover/analysis/MCCONNELL_DCHI_CALIBRATION.md:78-83`).

Cutoff is a mismatch. Producer defaults to 10 A
(`../src/CalculatorConfig.cpp:50-51`), while rediscover CLI defaults to 8 A
for both standalone and broad-backbone McConnell
(`src/rediscover/main_extract.cpp:165-185`). This can drop remote sources and
should be aligned, but it does not explain the near-field sign-flip pattern as
well as own-bond inclusion does.

I found no independent double-count beyond the same self/bonded inclusion
problem. Bonds are represented once by topology/source rows; there is no
evidence of PeptideCO/CN category duplication.

### Suspect 6: Units / Sign / Deweight

Verdict: CLEAN for the current scripts.

The current literature kernel uses the same project tensor with the scalar
scale `-prefactor*q/3`, where `prefactor = 1e24 / Avogadro`
(`src/rediscover/McConnellLiteratureKernel.cpp:31-88`). The calibration script
uses the same constants (`src/rediscover/analysis/mcconnell_dchi_calibration.py:28-34`)
and deweights by the exact same scale
(`src/rediscover/analysis/mcconnell_dchi_calibration.py:248-259`). The current
calibration report says source-sum vs aggregate RMS is 1.8e-8, manual q
conversion matches, and aromatic excluded rows sum to zero
(`src/rediscover/analysis/MCCONNELL_DCHI_CALIBRATION.md:78-83`).

One stale side path remains: older `literature_fixed_decirc.py` labels
`rediscover_mcconnell_aggregated_bare_kernel_T2.npy` as "fixed literature"
(`src/rediscover/analysis/literature_fixed_decirc.py:32-46`), but that is a
bare Angstrom^-3 sidecar per the SDK catalog
(`/shared/2026Thesis/nmr-shielding/python/nmr_extract/_catalog.py:211-222`).
The current `mcconnell_literature_decirc.py` and
`mcconnell_dchi_calibration.py` do not use that stale path for the reported
Delta-chi failures.

## Leading Explanation

The Stage-1-real mechanism reads as non-convergent in rediscover because
rediscover's broad-backbone McConnell source sum is polluted by local/own
amide-plane bonds that the producer-side calculator would reject. For N/C/O,
those are not small perturbations: they are the nearest possible anisotropic
bonds, evaluated at bond-midpoint distances where the point-dipole far-field
model is least valid. The C stratum's robust sign flip is the tell: carbonyl C
is the endpoint of C=O and often C-N, so own-bond tensors dominate a local T2
fit and can invert the apparent Delta-chi.

Stage 1 remains real evidence for McConnell as a mechanism, but it is not a
direct validation of this rediscover test. Stage 1 used a 720-protein WT-ALA
mutation-delta ridge with static McConnell features; rediscover is a one-protein
per-frame local-T2 de-circularisation/calibration. That quantity mismatch is
important context, not a reason to accept the current rediscover result as
physics.

## Ranked Candidate Fixes

Do not apply these in this vet.

1. Add explicit bond-source validity in rediscover broad/standalone McConnell.
   At minimum mirror producer `SelfSourceFilter` for endpoint target atoms.
   Prefer also recording an `is_self_or_bonded` flag analogous to ring so
   local/bonded bonds can be excluded from "valid" sums while still being
   auditable in source CSVs.

2. Mirror producer near-field semantics in rediscover: use source bond length
   and `near_field_exclusion_ratio` instead of only `r > 1e-6`. Decide and
   document whether rediscover intentionally goes stricter than producer by
   excluding bonded-neighbour bonds, especially for HN/HA.

3. Apply source validity consistently to `mc_lit_T2_local_*`,
   `bond_sum_dipolar`, and any source-summed category tensors. Preserve
   all-source columns separately for diagnostics.

4. Align rediscover McConnell cutoffs with producer defaults or record the
   calculator config used for comparison. Current default mismatch is 8 A in
   rediscover vs 10 A in producer.

5. After a re-emit in a separate, non-read-only task, rerun the decirc and
   Delta-chi calibration with an explicit before/after C-stratum audit:
   own C=O/C-N removed, near-field count removed, q_CO sign, absT2 correlation,
   and category SEs.
