# Welford data-shape expansion — design and execution plan (2026-05-17)

**Status:** Framework recorded 2026-05-17 after Codex review of the
Welford-TR batch landed today (commits `bdeedf9`, `609a59b`, `a068818`)
surfaced data-shape issues that the science-focused adversarial prompt
was meant to catch and that Codex caught manually instead.

**Decision recorded:** the bet from
`spec/plan/test-suite-realignment-deferred-2026-05-17.md` revises. Test
cleanup remains deferred (~1 week, schedule for later). **Data-shape
cleanup is NOT deferred** — it's the load-bearing layer for the entire
TR queue and for downstream calibration. ~12-18 hours of work before
the queue resumes. The user explicitly endorsed taking the extra time;
the "extra" is not extra at all because doing the queue with wrong data
shape produces an H5 set the fleet can't use.

**Cross-context discipline:** the principle and the per-Welford plan
in this document must survive session boundaries. Memory entries
`feedback_export_everything_upstream` and
`feedback_data_shape_codex_review_2026-05-17` pin the principle and
the trigger for posterity; PATTERNS.md gains the named principle so it
applies to every future TR; this plan doc is the execution checklist.

---

## The principle: forward-additive emission at upstream extraction

The C++ trajectory extractor is the upstream-of-everything position
in this codebase. Downstream: calibration in `learn/`, ML feature
engineering, the SDK that the viewer reads from, future GNN-side
ingestion, thesis-figure generation, manual analysis in R / Python /
notebooks that don't exist yet. The extractor does not know which of
these will need which field, in which scope, in which form.

Three implications follow:

1. **The H5 / NPY format moves forward, not backward.** Adding a
   column, a dataset, or an attribute is non-breaking — old readers
   ignore new fields; old data files read by new readers see
   missing fields they can default or skip. Removing or renaming a
   column is forever and breaks every consumer that touched it.
   Schema evolves additively or it doesn't evolve.

2. **The consumer is implicit.** Decisions like "drop the |T2|
   rollup because no consumer reads it" or "skip the T1 channels
   because BS AnomalyMarker only needs T0" optimise for known
   consumers and silently constrain future consumers. The
   minimum-viable instinct that drove the original BsWelford
   pattern is the wrong instinct for upstream extraction.

3. **Removals are forever; additions are free.** Storage is cheap
   (TrajectoryAtom growing from 136 to ~500 fields ≈ 4 KB per atom
   × 1500 atoms × 685 proteins ≈ 4 GB of fleet Welford state, no
   meaningful impact). Compute is cheap (online Welford update is
   O(1) per frame per atom per channel). The cost is line count
   on the struct and on H5 emission. That cost is paid once;
   downstream signal-gain is paid forever.

This is documented as a named pattern in PATTERNS.md "Lesson 25 —
Export Everything Upstream" and codified in memory as
`feedback_export_everything_upstream`.

## What gets expanded

The expansion is across four axes, applied per-Welford as appropriate.
Concrete per-Welford specs in the next section; here's the matrix.

### Axis 1: delta semantics — drift / mean_abs / rms / cadence

The current `*_delta_*` Welford channel accumulates signed raw
differences `Δx_n = x_n - x_{n-1}`. Mathematically `mean(Δ) = (x_N -
x_0) / (N-1)`, which telescopes — the running mean is endpoint drift,
not fluctuation. Four physics signals collapse into one weak number.

The expansion emits all four:

- **`*_drift_*`** (renamed from `*_delta_*` in semantic, kept as a
  separate dataset for compatibility): signed mean(Δx), captures
  equilibration trend. Near-zero on equilibrated runs by construction;
  nonzero on equilibrating runs.
- **`*_mean_abs_delta_*`**: `mean(|Δx_n|)`, total path length per
  frame step. Doesn't telescope; captures oscillation amplitude.
  Outlier-sensitive.
- **`*_rms_delta_*`**: `sqrt(mean(Δx² ))`. The standard
  fluctuation measure. Doesn't telescope; Gaussian-friendly.
- **Cadence metadata**: `mean_dt_ps` H5 attribute on the Welford
  group, derived from `traj.frame_times_`. Allows downstream to
  convert delta-per-stride to dx/dt in physical units.

For the integration test (stride=300 on 1.5 ns trajectory =
~300 ps between captured frames), all four are meaningfully
distinct. Calibration weights them independently.

### Axis 2: tensor channel preservation

The current Welfords collapse tensor sources to T0 (isotropic
scalar) + |T2| (Frobenius amplitude). PATTERNS.md "T2 Completeness"
says T2 is the primary analytical result; collapsing tensors to
scalars is rejection-worthy. The strict reading would forbid the
existing rollup; the pragmatic reading defers full T2 to the
time-series TR.

Under the export-everything principle, we emit BOTH:

- **`*_t2_m{−2,−1,0,+1,+2}_*` per-component T2 Welfords**: five
  separate accumulators on each T2 spherical-harmonic component,
  with the standard 8-stat block (mean, m2, std, min, max,
  min_frame, max_frame, n_frames). Preserves orientation.
- **`*_t2magnitude_*` |T2| Frobenius amplitude rollup**: the
  existing scalar, kept for downstream consumers (BS
  AnomalyMarker) that read it.
- **`*_t1_m{−1,0,+1}_*` per-component T1 Welfords**: for rank-1
  sources only (BS, HM where `G = -n⊗{B,V}` has T1 ≠ 0). McConnell
  / RingChi / HBond use the three-term McConnell form per
  PATTERNS Lesson 19 with T1 = 0 by construction — these skip
  the T1 block.
- **`*_t0_*`**: existing T0 rollup unchanged.

Layout ordering follows the SphericalTensor real-spherical-harmonic
convention documented elsewhere: T1 = [m−1, m0, m+1], T2 = [m−2,
m−1, m0, m+1, m+2]. This is the same ordering used by the
time-series TR's H5 `xyz` dataset packing.

### Axis 3: schema provenance

The Welford accumulators track more state than they emit. The H5
emission boundary discards:

- **`*_min_frame`, `*_max_frame`** per channel: lets a downstream
  reader go back to the trajectory and inspect the frame that
  produced the extremum.
- **`*_m2`** (raw Welford sum-of-squared-deviations) per channel:
  the variance numerator. Emitting it lets downstream recompute
  std under a different convention (population n vs. unbiased
  n−1) if their math requires it.
- **`delta_n`** per channel: the count of delta samples (= n_frames
  − 1 for AV-pattern Welfords). Needed for unbiased delta std.
- **H5 group attributes**:
  - `ddof=1` — bias-convention disclosure for the std datasets
  - `mean_dt_ps` — cadence metadata (Axis 1)
  - `frame_index_range=[first, last]` — span of trajectory covered
  - `welford_n_frames` — total frame count (group-level)

All additive. No existing dataset changes name or shape.

### Axis 4: distribution-shape additions where the physics calls for it

For HBondCount specifically:

- **`hbond_count_occupancy_fraction_*`** Welford on the indicator
  `(count > 0 ? 1.0 : 0.0)`. Captures "what fraction of frames is
  this atom engaged in any H-bond?" — distinct from
  `hbond_count_mean` which is `⟨N⟩` expected count. Both signals
  are independent features for ML feature engineering.

For all Welfords, **skewness and excess kurtosis** are cheap
additional higher-moment accumulators (m3 and m4 online via standard
formulas) that capture non-Gaussian distribution shape. **Deferred**
to a later iteration unless the calibration pipeline shows
sensitivity to non-Gaussian behavior — emit if cheap, otherwise wait
for evidence of need. Decision: **emit** because consumer-implicit,
forward-additive; cost is two more accumulators per channel.

## Per-Welford concrete spec

The fields are per-atom on `TrajectoryAtom`; the H5 emission shape
mirrors the field layout.

### BsWelfordTrajectoryResult (BiotSavart, rank-1)

Source: `ConformationAtom::bs_shielding_contribution` (SphericalTensor,
units = ppm·T/nA per OBJECT_MODEL drift table).

Channels: T0 (scalar) + T1 (3 components) + T2 (5 components) +
|T2| Frobenius (scalar).

Welford block per channel: `mean`, `m2`, `std`, `min`, `max`,
`min_frame`, `max_frame`, `n_frames` (8 doubles + 2 size_t = 10
fields per channel).

10 channels × 10 fields = 100 fields per atom.

Plus delta variants on T0 only (drift, mean_abs_delta, rms_delta,
each with mean/m2/std/min/max/n = 6 fields): 18 fields.

Plus higher moments on T0 only (m3, m4 → skew, kurtosis on Finalize):
4 fields.

Total: ~122 fields per atom for BS Welford. Cadence metadata is
shared at group level.

### HmWelfordTrajectoryResult (HaighMallion, rank-1)

Source: `hm_shielding_contribution` (SphericalTensor, units = Å⁻¹).

Same shape as BS. ~122 fields per atom.

### McConnellWelfordTrajectoryResult (full McConnell-form)

Source: `mc_shielding_contribution` (SphericalTensor, units = Å⁻³,
asymmetric non-traceless three-term form per Lesson 19).

Channels: T0 + T2 (5) + |T2|. **T1 = 0 by construction; skipped.**

7 channels × 10 fields = 70 fields + delta+higher-moment on T0 = 22.
Total ~92 fields per atom.

### EeqWelfordTrajectoryResult (scalar source)

Source: `eeq_charge` (double, units = elementary_charge).

One channel: value. 10 fields + delta variants (18) + higher
moments (4). ~32 fields per atom.

### SasaWelfordTrajectoryResult (scalar source)

Source: `atom_sasa` (double, units = Å²).

Same shape as Eeq. ~32 fields per atom.

### HBondCountWelfordTrajectoryResult (integer scalar source)

Source: `hbond_count_within_3_5A` (int, promoted to double).

Two channels: count + occupancy indicator. 2 × 10 = 20 fields +
delta variants on count (18) + higher moments (4). ~42 fields per
atom.

### Total field-budget growth

Existing TrajectoryAtom fields (after today's batch): ~136
Welford-related.

Expanded: ~122 + 122 + 92 + 32 + 32 + 42 = ~442 Welford-related
fields. **Net add: ~306 fields per atom on TrajectoryAtom.**

Storage impact at fleet scale: ~306 × 8 bytes × 1500 atoms × 685
proteins = ~2.5 GB additional Welford state across the fleet.
Negligible.

## H5 schema emission

Each Welford emits a single H5 group `/trajectory/{name}_welford/`
with one dataset per channel-statistic combination + the group-level
attributes. Naming convention:

```
/trajectory/bs_welford/
  attributes:
    result_name = "BsWelfordTrajectoryResult"
    n_frames = 50
    finalized = true
    ddof = 1
    mean_dt_ps = 300.0
    frame_index_range = [0, 1500]
    welford_n_frames = 50
    units_t0 = "ppm_T_per_nA"
    units_t1 = "ppm_T_per_nA"
    units_t2 = "ppm_T_per_nA"
    units_t2magnitude = "ppm_T_per_nA"
    irrep_layout_t1 = "m-1, m0, m+1"
    irrep_layout_t2 = "m-2, m-1, m0, m+1, m+2"
  datasets (per atom, shape (N_atoms,)):
    t0_mean, t0_m2, t0_std, t0_min, t0_max, t0_min_frame, t0_max_frame, t0_n_frames
    t0_skew, t0_kurtosis  (Finalize-derived)
    t1_m-1_mean, ..., t1_m-1_max_frame  (8 per channel)
    t1_m0_*, t1_m+1_*  (same)
    t2_m-2_*, t2_m-1_*, t2_m0_*, t2_m+1_*, t2_m+2_*  (same)
    t2magnitude_*  (same 8)
    t0_drift_mean, t0_drift_m2, t0_drift_std, t0_drift_min, t0_drift_max, t0_drift_n
    t0_mean_abs_delta_mean, ..., t0_mean_abs_delta_n
    t0_rms_delta_mean, ..., t0_rms_delta_n
```

For McConnell: drop the `t1_*` datasets. For Eeq / Sasa: drop
`t1_*`, `t2_*`, `t2magnitude_*`, rename `t0_*` → `value_*`. For
HBondCount: scalar channel + occupancy channel.

Old datasets (`t0_delta_mean`, etc.) **stay emitted** for forward
compatibility, with a `deprecated_use="t0_drift_mean"` attribute
on each. Consumers that ported to the new name use it; old
consumers continue to work. No name reuse.

## Python SDK consistency

`python/nmr_extract/_catalog.py` is audited end-to-end against the
expanded H5 schema. Two passes:

1. **Drift fix**: every existing ArraySpec that labels a kernel
   as "ppm" gets corrected per the OBJECT_MODEL drift table
   (HM → Å⁻¹, McConnell → Å⁻³, PiQuad → Å⁻⁵, Dispersion → Å⁻⁶,
   HBond → Å⁻³). ~10-15 entries to update.

2. **Schema expansion**: new ArraySpec entries for the expanded
   Welford channels — per-component T1/T2, drift variants,
   higher moments, occupancy. ~80-100 new entries across the 6
   Welfords.

Deferred items in the SDK: `tripeptide_*` and `larsen_*` TR groups
already have ArraySpec entries from prior batches; verify they
match the new conventions but no expected drift there.

## Test discipline (under the deferred-cleanup bet)

Each Welford keeps its existing integration test
(`Integration1P9J`). Test additions:

- **Magnitude assertions for new channels**: confirm the new T1
  / T2 per-component datasets are populated (`populated > 0`)
  and have plausible magnitudes. Same shape as the existing T0
  magnitude assertion.
- **Drift vs. dynamics smoke check**: assert that
  `t0_rms_delta_mean > 0` on the trajectory (something is
  fluctuating; if not, calculator is broken). Assert that
  `t0_drift_mean ≠ t0_rms_delta_mean` (they're distinct
  signals).
- **Schema attribute checks**: H5RoundTrip test verifies
  `ddof=1` and `mean_dt_ps` attributes are emitted and readable.

The discipline trio (Frame0Semantics, FinalizeIdempotency,
H5RoundTrip) stays as-is. Don't expand them — that's
deferred-cleanup territory per the existing bet.

## Execution plan

### Phase 1 — framework establishment (this session)

1. ✓ Write this design doc.
2. Update `PATTERNS.md` to add "Lesson 25 — Export Everything
   Upstream" as a named principle.
3. Update `OBJECT_MODEL.md` to reference this design doc and
   note the Welford schema is expanding (pointer, not duplicate
   spec).
4. Update `spec/plan/test-suite-realignment-deferred-2026-05-17.md`
   to record the bet revision.
5. Update
   `spec/plan/adversarial-review-prompt-science-focus.md` with
   tighter |T2|-scalar framing and SDK-side reading list.
6. Create memory entries
   `feedback_export_everything_upstream` and
   `feedback_data_shape_codex_review_2026-05-17`. Update
   `MEMORY.md`.
7. Expand `TrajectoryAtom.h` with all new fields for all six
   Welfords (central, sequential — touching one shared file).
8. Expand `BsWelfordTrajectoryResult` (`.h`, `.cpp`, test) as
   the exemplar. Land + run to verify the pattern.
9. Commit Phase 1.

### Phase 2 — parallel expansion of remaining five Welfords

Once the framework is committed and BS is the validated exemplar:

- **Worktree A (one agent, mechanical)**: expand HmWelford +
  McConnellWelford (the two tensor-source siblings to BS). Same
  shape minus McConnell's T1 channels.
- **Worktree B (one agent, mechanical)**: expand EeqWelford +
  SasaWelford + HBondCountWelford (the three scalar-source
  siblings). HBondCount additionally gets the occupancy
  channel.

Two agents in parallel, per the discipline of two-parallel-max
for mechanical translation
(`feedback_agent_quality_variance_on_same_brief`). Each agent
gets the BS exemplar as the clone target plus this design doc as
the spec.

### Phase 3 — SDK catalog audit

Walk `python/nmr_extract/_catalog.py` for the two passes above
(drift fix + schema expansion). Can be a single focused session
or delegated to an agent with the OBJECT_MODEL drift table as
the spec.

### Phase 4 — TR queue resumes

Remaining ~12 TRs (ApbsEfg, Water + Hydration trio, DSSP family,
Bonded/Gromacs energy, AIMNet2 embedding opt-in) land with the
expanded shape from the start. Each new Welford-style TR
inherits the per-channel discipline; each new FO TR inherits the
schema-provenance discipline (min_frame, max_frame, m2, dt
cadence, ddof attribute).

The water TR shape (NPY vs. H5 group) is still an open
architectural decision — handled in conversation before that
batch.

## Open items captured

- **Water TR shape**: undecided. NPY emission vs. H5 group; the
  per-atom + per-frame + per-water-molecule data is a three-axis
  array that doesn't fit DenseBuffer cleanly. Discuss before
  that batch.
- **Higher moments emission for non-scalar channels**: this doc
  spec'd higher moments (skew, kurtosis) on T0 only for tensor
  Welfords. Should we also emit them per-T1-component and
  per-T2-component? Cost is small; emit if no objection. Default
  to "yes" under the export-everything principle.
- **Adversarial-prompt refinement**: the science-focused prompt
  caught most of Codex's findings in principle but the |T2|
  scalar question was framed as "T0 vs T2 mis-pick" rather than
  "scalar summary discards orientation." Tighten in the prompt
  doc.

## Cross-context survival hooks

If a future session lands without conversation context:

- `MEMORY.md` lists `feedback_export_everything_upstream` and
  `feedback_data_shape_codex_review_2026-05-17`. Loading these
  memories establishes the principle.
- `PATTERNS.md` Lesson 25 (added in this batch) makes the
  principle a named pattern that applies to every TR.
- This doc (`spec/plan/welford-data-shape-design-2026-05-17.md`)
  is the comprehensive execution checklist. Per-Welford specs
  are concrete enough that an agent can pick up Phase 2 work
  from a cold start.
- `spec/plan/test-suite-realignment-deferred-2026-05-17.md` is
  updated to reflect the bet revision so the prior plan doesn't
  contradict.

If a future session is just resuming the TR queue (Phase 4): the
remaining TR list is in
`spec/plan/test-suite-realignment-deferred-2026-05-17.md`
"Open architectural item: water" plus the OBJECT_MODEL
PerFrameExtractionSet catalog (⏳ rows still pending after this
expansion).

## On the bet revision

The previous bet
(`spec/plan/test-suite-realignment-deferred-2026-05-17.md`) said:
finish TR queue tomorrow night against the current data shape,
defer test cleanup. The bet's load-bearing mitigation was a
science-focused adversarial reviewer on every batch.

Codex review on commits `bdeedf9` + `609a59b` + `a068818` ran
the science role manually and found data-shape issues. The bet
as offered wasn't viable — the science-focused review wasn't a
mitigation, it was the gate that revealed the bet was wrong-shaped.
Doing the queue at speed against wrong data shape produces a
fleet H5 set that fails calibration silently.

The revised bet:

1. **Data-shape cleanup is NOT deferred.** Land now (Phase 1+2+3),
   ~12-18 hours of work.
2. **Test cleanup remains deferred** per the original plan doc.
   Test discipline does NOT block the data-shape work; each
   expanded Welford keeps the same integration test shape plus
   small additive assertions.
3. **Science-focused adversarial review still required** for the
   resumed TR queue (Phase 4). Now its job is verifying the
   expanded data shape carries through correctly to each new TR.
4. **Rollback point**:
   `git tag pre-test-config-consolidation-2026-05-17` is still
   valid on `origin`. If Phase 2 or 3 reveals a deeper design
   problem, that tag is the recovery anchor.

The user explicitly endorsed the extra time on 2026-05-17 with
the framing: "We are lucky to have caught this and squared away
on goals." That's the operative spirit — the catch is a win, the
expansion is forward progress, the schedule absorbs.
