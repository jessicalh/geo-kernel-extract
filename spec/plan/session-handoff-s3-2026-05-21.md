# Session 3 handoff — RingNeighbourhood + Rmsd trio + S5 closeout

**Status as of master `89318fd` (pushed origin/master 2026-05-21):**

S2 closed. TR4-TR9 landed + three fix bundles (sci/math/anti-pattern
adversarial reviews + codex post-bundle review + the use-case-doc-aware
trim). All tests green; anti-pattern grep clean across the bundle; the
`--mopac` CLI flag is wired to FullFatFrameExtraction for the 1P9J +
mutant use case.

This doc covers the remaining sessions of the 13-TR plan: S3
(RingNeighbourhood), S4 (Rmsd trio + 1P9J production validation), and
S5 (codex final annealing + 4 Python analysis scripts + RETIRED
markers). Together they reach **code complete** — every planned
calculator landed, every Python analysis script in place, every
retired calculator marked, internal interfaces stable. Cleanup phase
(static analysis, per-calc audits, test simplification, doc cleanup)
operates AFTER code complete on a stable surface.

The S2 handoff doc is the structural precedent for this file's shape.

---

## Code complete: definition

"Code complete" here is the McConnell milestone:

- All planned features coded — every TR in the 13-TR plan landed,
  every Python analysis script implemented, every retired calculator
  marked.
- Internal interfaces locked — SDK schemas, H5 attrs, OBJECT_MODEL
  rows stable. Cleanup phase can operate without schema churn.
- Known issues durably recorded — memory entries for anything
  deferred (APBS-vs-Coulomb EFG sign, parity 1o/1e convention,
  MopacCoulomb EFG clamp, RingNeighbourhood architectural debt
  — see below).
- Tests cover the canonical Layer 0 trio per TR (Frame0, Idempotency,
  Integration1P9J) + any genuinely-new code paths.

NOT minimum-viable. Comprehensive. The week-long cleanup phase
(static analysis tooling, per-calc + per-TR audits from Claude and
codex in parallel, test simplification, comment + doc cleanup) is
what gets us from code complete to RTM.

---

## S3 — `RingNeighbourhoodTrajectoryStats` (TR10)

### Substrate state — clean

Substrate side is post-Bundle-B/C clean:

- `src/Ring.h` — typed Ring class hierarchy (PheBenzene / TyrPhenol /
  TrpBenzene / TrpPyrrole / HisImidazole / HidImidazole /
  HieImidazole / ProPyrrolidine / IndolePerimeter). Each carries its
  physics through virtuals: `Intensity()`, `JBLobeOffset()`,
  `NitrogenCount()`, `Aromaticity()`, `RingSizeValue()`, `TypeName()`.
- `src/RingTopology.h` — substrate-driven storage on
  `LegacyAmberTopology`, constructed via `ConstructFromSubstrate`
  from the `AtomSemanticTable` substrate. Aromatic + saturated
  split. Canonical cyclic walk per `RingSystemKind`. Bundle C / Slice
  B (2026-05-07).
- API surface: `protein.RingCount()` + `protein.RingAt(i)` (aromatic
  only); saturated ring puckers emitted via `RingPuckerTimeSeries`
  separately.

No string identity. No cifpp re-traversal at runtime. No
pre-LegacyAmberTopology assumptions in the typed surface.

### Consumer-side state — DEBT (user-flagged 2026-05-21)

Five ring-using calculators each independently push entries to
`ConformationAtom::ring_neighbours`:

| Calc | File | Cutoff | Calc-specific fields written |
|---|---|---|---|
| BiotSavart | `src/BiotSavartResult.cpp:165, 263-273` | `ring_current_spatial_cutoff` (15 Å) | `G_tensor`, `G_spherical`, `B_field`, `B_cylindrical`, `cos_phi`, `sin_phi` |
| HaighMallion | `src/HaighMallionResult.cpp:230, 285-294` | same | `hm_H_tensor`, `hm_H_spherical`, `hm_B_field`, `hm_G_tensor`, `hm_G_spherical` |
| RingSusceptibility | `src/RingSusceptibilityResult.cpp:146` | same | `chi_tensor`, `chi_spherical`, `chi_scalar` |
| PiQuadrupole | `src/PiQuadrupoleResult.cpp:151` | same | `quad_tensor`, `quad_spherical`, `quad_scalar` |
| Dispersion | `src/DispersionResult.cpp:214, 302-327` | same | `disp_tensor`, `disp_spherical`, `disp_scalar`, `disp_contacts` |

The pattern (visible at `DispersionResult.cpp:302-327`):

```cpp
RingNeighbourhood* rn = nullptr;
for (auto& existing : ca.ring_neighbours)
    if (existing.ring_index == ri) { rn = &existing; break; }
if (!rn) { /* create with geometric fields, push */ }
rn->disp_tensor = ...;  // calc-specific fields
```

Fat-union struct on `ConformationAtom.h:31-61` carries fields from
ALL five calculators. Find-or-create-then-write pattern; whoever
creates the entry sets the geometric fields, subsequent calculators
fill in their own slots.

**Mitigations that hold today (verified):**
- All five calculators use **identical** `ring_current_spatial_cutoff`
  (`CalculatorConfig::Get(...)`, no per-calc overrides) — so the
  (atom, ring) pair set is identical across calculators. No
  "atom-ring entries that only some calcs see" leakage.
- The geometric fields (`distance_to_center`, `direction_to_center`,
  `rho`, `z`, `theta`) are calculator-invariant; whoever-runs-first
  sets them consistently.
- `cos_phi`/`sin_phi` are BS-specific (BS computes them at lines
  263-273) — if BS runs FIRST in RunConfiguration's attach order,
  these are populated. If BS runs LATER, they're zero-default until
  BS gets there. Need to verify attach order locks BS-first; this is
  a real ordering fragility.

**The debt:**
- **D1**: Fat-union struct on `ConformationAtom` carrying calc-specific
  fields. Hard to audit; hard to extend without touching the shared
  struct.
- **D2**: Five parallel writers with find-or-create + per-field write.
  Architecturally noisy; intent is implicit (each calc "claims" its
  field set), not declared.
- **D3**: Ordering fragility on `cos_phi`/`sin_phi` (BS-owned fields,
  set during BS's pass; other calcs use the entry but don't compute
  these).
- **D4**: `gaussian_density` field defined on the struct but not
  populated by any of the 5 calcs I read. Possibly dead.

**The S3 architectural choice** — three paths, codex Pattern A
prompt asks for the decision:

- **A1**: TR10 reads `ca.ring_neighbours` as-is. Emits per-(atom,
  ring) Pattern A with all calc-specific fields. Inherits the fat-
  union + parallel-write mess. Lowest cost in this session;
  architecturally inherits D1-D3.
- **A2**: Lift the per-(atom, ring) computation into a shared
  `RingNeighbourhoodResult` ConformationResult. Each ring-using
  calc READS from it (sets only its calc-specific fields locally,
  not on the shared struct). TR10 reads from it. Substrate cleanup
  pays off the debt. Higher cost — touches 5 calculators.
- **A3**: TR10 captures GEOMETRIC fields only (distance, rho, z,
  theta, azimuth, ring_type, ring_index). Calc-specific T2
  contributions stay on per-atom shielding fields and are emitted
  by the existing BS/HM/Mc/PiQuad/RingChi/Dispersion shielding TS
  TRs (already in the bundle). TR10 captures only the geometric
  neighborhood, not the per-calc physics. Lean; defers D1-D3 to
  cleanup phase.

Per the user's framing: storage isn't the constraint
("if scientifically worth it, worth the disk"). The choice is
architectural — what the right shape is going INTO cleanup phase.

### Codex Pattern A prompt (paste into codex)

```
The S3 TR `RingNeighbourhoodTrajectoryStats` will emit per-(atom, ring)
data as a trajectory time series. The 13-TR plan (spec/plan/13-tr-build-plan-2026-05-20.md
§ Pattern A strawman) outlines a flat-fixed-R-per-atom-max layout with
sentinel -1 ring_index for unfilled slots, plus per-frame metadata.

Before locking the layout, decide the architecture.

CONSUMER-SIDE STATE: five ring-using ConformationResult calculators
(BiotSavartResult, HaighMallionResult, RingSusceptibilityResult,
PiQuadrupoleResult, DispersionResult) independently push entries to
`ConformationAtom::ring_neighbours` (a vector of fat-union
`RingNeighbourhood` structs at `ConformationAtom.h:31-61`). All five
use identical `ring_current_spatial_cutoff` so the (atom, ring) pair
set is identical across calcs; find-or-create-then-write is the
pattern. The struct carries fields from all five calculators.

SUBSTRATE STATE: clean. RingTopology constructed substrate-driven from
AtomSemanticTable; typed Ring class hierarchy with physics on virtuals.

QUESTIONS:

1. Architectural choice for TR10:
   A1) Read ca.ring_neighbours as-is, emit per-(atom, ring) Pattern A
       with all calc-specific fields (BS G_spherical, HM hm_G_spherical,
       RingChi chi_spherical, PiQuad quad_spherical, Dispersion
       disp_spherical). Inherits the fat-union + parallel-write design.
   A2) Lift to a shared `RingNeighbourhoodResult` ConformationResult.
       Each ring-using calc READS from it; only its own calc-specific
       output stays local. TR10 reads from the canonical struct.
       Substrate cleanup; touches 5 calculators.
   A3) TR10 captures geometric fields only (distance, rho, z, theta,
       cos_phi, sin_phi, ring_type, ring_index, ring_membership).
       Calc-specific T2 contributions live on the already-landed
       per-calc shielding TS TRs (BsShieldingTS, HmShieldingTS, etc.).
       Geometric-only TR10; defers the fat-union cleanup to the
       cleanup phase.

   Which is the right architecture for code complete? Storage is not
   a constraint (HDF5 + gzip-9; 128GB host).

2. Layout (if A1 or A2 chosen — Pattern A flat-vs-ragged):
   - Flat fixed-R-per-atom-max with sentinel ring_index = -1 for
     unfilled slots? Determine `R_per_atom_max` from the trajectory's
     observed max-rings-within-15Å per atom.
   - Ragged variable-length per-atom lists?
   - Sentinel convention if flat (uint32 max? int32 -1? float NaN
     for unfilled-slot scalar fields?).

3. Ring membership timing:
   - `ring_membership_per_atom: (N, R_per_atom_max) int32` — static
     (frame-invariant since ring topology doesn't change post-Seed)
     or per-frame (allows future ring-perception evolution)?
   - Static recommended per the plan strawman; confirm.

4. ProPyrrolidine handling:
   - Saturated rings (Pro pyrrolidine) live on RingTopology.saturated_
     and emit via RingPuckerTimeSeries already. Should they ALSO
     appear in TR10's neighborhood emission, or are they categorically
     excluded (aromatic-only, matching `protein.RingAt(i)`)?

Return: (a) architectural choice (A1/A2/A3) with rationale, (b)
locked layout decisions, (c) ProPyrrolidine inclusion call. Ideally
one round; design questions are scoped.
```

### S3 implementation order (once codex returns)

1. If A2 chosen: refactor pass on 5 ring-using calcs to read from
   shared `RingNeighbourhoodResult`. Likely 1-2 sessions ahead of
   TR10 implementation. Bundle as a separate commit before S3 TR10
   work begins.
2. Implement TR10 RingNeighbourhoodTrajectoryStats:
   - Header + cpp following the AIMNet2 CRG TS / TripeptideBackbone
     TS template (FO with DenseBuffer or per-atom growing buffer
     depending on layout)
   - `RingNeighbourhoodTrajectoryStats::Compute` per frame: for each
     atom, read `ca.ring_neighbours` (A1) or shared result (A2) or
     compute geometry (A3); push to per-atom buffer
   - `Finalize`: transfer to DenseBuffer or finalize in-place
   - `WriteH5Group`: emit `/trajectory/ring_neighbourhood/` with
     locked layout
3. Register in RunConfiguration: which config? Both PerFrame +
   FullFat get ring calcs (PerFrame has them), so TR10 registration
   goes in PerFrameExtractionSet — unlike TR5-TR9 which were
   FullFat-only.
4. Python SDK: `RingNeighbourhoodGroup` + loader + TrajectoryData
   wire-in.
5. OBJECT_MODEL.md row.
6. Tests: Layer 0 quartet (Frame0Semantics, FinalizeIdempotency,
   Integration1P9J).
7. Literature probe test: 1P9J Phe33 ring → Hα of Cα-stem residue at
   ~5 Å in X-ray. Per the 13-TR plan.
8. Codex annealing pass at S3 close.

### S3 memory entries to land

- `project_ring_neighbourhood_debt_2026-05-21` — record the fat-union
  + parallel-writer pattern as cleanup-phase candidate (whichever
  path is chosen for S3).

---

## S4 — Rmsd trio + DftPoseCoordinator + 1P9J production validation

### Templates / precedents

- **TR11 `RmsdTrackingTrajectoryResult`** — new canonical
  (reference-frame scalar AV). Per plan §"Locked decisions":
  - Reference frame = trajectory frame 0
  - Atom selection = backbone heavy atoms (N, CA, C, O)
  - Output: (T,) double in nm or Å (clarify at implementation)
  - Lifecycle: AV (running update each frame); finalize trivial
  - No precedent in the codebase for "scalar AV (T,) double" — first
    of this shape. The simpler shape (per-frame value, no per-atom
    state) means likely no DenseBuffer; per-frame vector accumulation
    + Finalize emits as (T,) flat.

- **TR12 `RmsdSpikeSelectionTrajectoryResult`** — near-mechanical
  clone of `ChiRotamerSelectionTrajectoryResult`
  (`src/ChiRotamerSelectionTrajectoryResult.{h,cpp}`):
  - Internal state: per-frame previous RMSD + cooldown counter
  - Compute reads RmsdTracking's per-frame value, compares to
    threshold + cooldown, pushes SelectionRecord onto
    `traj.MutableSelections()` on spike
  - Threshold from plan: 2.0 Å vs baseline frame 0 + 0.5 Å vs
    rolling 100-frame mean; cooldown 100 frames (~2 ns at 20 ps
    stride)
  - Cross-TR read: needs TR11 RmsdTracking attached first → declare
    `Dependencies()` returning `typeid(RmsdTrackingTrajectoryResult)`.
    See PATTERNS §17 cross-result-read marker discipline.

- **TR13 `DftPoseCoordinatorTrajectoryResult`** — new canonical
  (cross-TR SelectionBag aggregator):
  - At Finalize, read `traj.Selections().ByKind<RmsdSpikeSelectionTrajectoryResult>()`
    and `traj.Selections().ByKind<ChiRotamerSelectionTrajectoryResult>()`
  - Dedupe key: `(residue_index, frame_idx // 50)` ns-bucket (~1 ns
    at 20 ps stride per plan §"Locked decisions")
  - Push reduced set back to bag under own kind
  - SelectionRecord shape: `src/SelectionRecord.h` — `kind` (typeid),
    `frame_idx`, `time_ps`, `reason` string, free-form metadata map
  - RecordBag<SelectionRecord> API on Trajectory:
    `traj.Selections()` + `traj.MutableSelections()`; iteration via
    `ByKind<T>()`

### 1P9J production validation (S4 closeout)

After TR10-TR13 land, run the full 13-TR bundle against the 1P9J
fleet_amber fixture:

```
nmr_extract --trajectory tests/data/fleet_amber/1P9J_5801/.../production \
            --analysis-h5 ... --mopac \
            --aimnet2 data/models/aimnet2_wb97m_0.jpt \
            --output runs/1p9j_full_13tr_validation
```

Expected output H5 should contain ALL 13 TR groups (plus the
pre-existing TR groups from S1's Tripeptide/Larsen bundles etc.).
Verify via Python SDK's `load_trajectory()` that each group reads
back into a populated dataclass.

If MOPAC stride is desired ≠ trajectory stride: per user 2026-05-21,
MOPAC stride matches the 750-DFT cadence or stride 1. Pass via
existing stride mechanism; no new flag needed (per user).

### S4 implementation order

1. TR11 RmsdTracking first (TR12/TR13 depend)
2. TR12 RmsdSpike (clones ChiRotamerSelection)
3. TR13 DftPoseCoordinator (reads bags from TR12 + ChiRotamer)
4. 1P9J production validation run
5. Codex annealing pass at S4 close

---

## S5 — Codex final annealing + 4 Python analyses + RETIRED markers

### Codex final annealing

Multi-round on the full 13-TR bundle. Per the plan: "Codex review at
this stage is an annealing set, not single-shot — multiple rounds
expected, each prompt narrower than the last. Final pass (Session 5)
expects 3-5 rounds to drive the bundle to commit-readiness."

Pattern from S1 and S2: I draft the round-1 prompt (broad audit),
user pastes into codex, codex returns findings, I bundle fixes, we
iterate. Each subsequent round narrows scope per the findings the
previous round surfaced.

### Four Python analysis scripts in `learn/`

Per the 13-TR plan §Markwick + §Session 5:

1. **`block_avg_convergence.py`** — Grossfield-Zuckerman 2009
   block-averaged SEM. Per-atom σ-channel mean convergence diagnostic.
   Gates every other Stage 2 claim ("did we sample enough frames for
   the mean to be reliable?"). Reads the shielding TS NPYs;
   block-average across increasing block sizes; SEM curve per channel.

2. **`sigma_lipari_szabo.py`** — Per-atom σ-S² order parameter +
   per-residue groupby + cross-report against BMRB 5801 NH-S². Reads
   shielding TS NPYs; computes second moment of σ across frames per
   atom; groups by residue_index_per_atom; compares against published
   BMRB 5801 entry's per-residue NH-S² values.

3. **`sigma_essential_dynamics.py`** — SVD on per-frame σ tensor
   stack `(T, N, K)`. Tells us if signal is low-rank. PCA / essential
   dynamics on the shielding-tensor manifold. Reports leading
   eigenvalues + leading eigenvectors per atom-class group.

4. **`markwick_overlay_1p9j.py`** — The 1P9J chapter-section
   deliverable. σ-residual ↔ BMRB 5801 Rex spatial overlap.
   - Inputs: 7+ shielding TS NPYs (BS, HM, Mc, PiQuad, RingSusc,
     Dispersion, HBond) + BMRB 5801 Rex per residue (pre-stage as
     per-protein CSV/NPY from `references/incoming/_backup_*` per
     the plan)
   - σ-residual per residue: groupby on `residue_index_per_atom`;
     standard deviation of σ across trajectory per residue
   - Overlay: scatter of σ-residual vs Rex per residue, color-coded
     by secondary structure
   - Output: `learn/markwick_overlay_1p9j/figure.pdf` + chapter prose
   - Plan estimate: ~80 LOC Python

Per the plan: "Zero protein-model rebuilding in Python (the
residue_index_per_atom broadcast handles binding)." All four scripts
use the Python SDK's existing TrajectoryData accessors.

### RETIRED markers in `spec/PLANNED_CALCULATORS_*`

Per the 13-TR plan's retirement table, add `**RETIRED 2026-05-XX**:
<reason>` markers (preserving original entries for archaeology) to:

**In `spec/PLANNED_CALCULATORS_2026-04-22.md`:**
| Section | Calculator | Reason |
|---|---|---|
| §1 | `GreenKuboSpectralDensityResult` | Spectral view redundant with `SigmaEssentialDynamics` PCA at our 15 ns timescale |
| §2 | `PseudocontactShiftResult` | No paramagnetic substrate in scope; 1P9J + fleet diamagnetic |
| §3 | Per-secondary-structure CSA stratification (as C++ TR) | Groupby on existing DSSP TS NPYs; Python at chapter draft |
| §4 | Lie-group GP regression on SE(3) | Stage 3 modeling work, not a calculator |

**In `spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md`:**
| Section | Calculator | Reason |
|---|---|---|
| §1 | `BlockAveragedConvergenceResult` (as C++ TR) | Ships as Python analysis script (script #1 above) |
| §2 | `SigmaLipariSzaboResult` (as C++ TR) | Ships as Python analysis script (script #2) |
| §3 | `SigmaEssentialDynamicsResult` (as C++ TR) | Ships as Python analysis script (script #3) |
| §4 | `SigmaTuckerDecompositionResult` | Higher-order ED; defer behind PCA. Conditional trigger: PCA shows multi-axis structure 2D SVD can't unmix |
| §5 | `CrossCorrelatedRelaxationResult` | Backbone CSA-dipolar cross-correlation needs ns-ms; 15 ns timescale-marginal |
| §6 | `SigmaMSMResult` | Markov state model needs ms-timescale sampling; 15 ns has ~3-4 rotamer transitions max |
| §7 | `SigmaMemoryKernelResult` | Mori-Zwanzig formulation needs long-memory regime; 15 ns insufficient |
| §8 | Aromatic-H geometry sanity-check extraction | Diagnostic only; mdtraj one-liner at write-up |

---

## Closing — code complete checkpoint signals what shifts

When S5 closes (TR10-13 landed + Python scripts + RETIRED markers +
codex annealing settles), the project is **code complete**. The
posture shift is:

- Stop adding features. Anything else surfacing in cleanup phase
  goes to memory entries or a follow-up backlog, NOT into S3-S5
  scope.
- Static analysis baseline goes in (user's call: clang-tidy,
  cppcheck, ASan, UBSan, valgrind — pick before per-calc audits
  start so tooling frames what gets reviewed).
- Per-calc + per-TR audits begin. Parallel agents (the same shape
  the S2 fix bundles used) on a canonicalized audit prompt.
- Test simplification: many Layer 0 patterns are cargo-cult clones
  (Frame0Semantics, FinalizeIdempotency); cleanup pass canonicalizes
  + dedupes.
- Comment + doc pass: trim provenance commentary that's settled,
  update OBJECT_MODEL where it's drifted, retire bones/.

These are the cleanup-phase decisions; not in scope for S3-S5. Flag
them here so the next session(s) know the bar.

---

## What's pushed and where

- Branch: `master` (origin/master at `89318fd`)
- Remote: `https://github.com/jessicalh/geo-kernel-extract.git`
- Tests green at HEAD; --mopac CLI wired; net_charge plumbing
  correct; sentinel-aware Welford on TR6; signed cos on TR9.
- Memory entries durable for cleanup-phase reference:
  - `project_s2_mopac_bundle_landed`
  - `feedback_apbs_efield_gate_retrofit`
  - (S3 will add) `project_ring_neighbourhood_debt_2026-05-21`
