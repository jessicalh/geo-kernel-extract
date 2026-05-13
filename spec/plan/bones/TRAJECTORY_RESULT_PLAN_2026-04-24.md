# TrajectoryResult Plan — 2026-04-24

**Status:** note, not spec. Consolidates the working plan for the next
round of `TrajectoryResult` additions and the configuration shape they
need. Composed during the first session asked for new trajectory
results since the 2026-04-24 exemplars landed. Reflects the strategic
shift to the 2-protein × 15 ns × 40 ps DFT-paired dense validation set
and the object-model discipline grounded in PATTERNS.md §§1-18 +
OBJECT_MODEL.md.

Companion docs (read together):

- `spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md` — the
  references-to-pipeline accretion file (ACF family discipline,
  per-paper candidate list).
- `spec/PLANNED_CALCULATORS_2026-04-22.md` — five calculator / model-
  architecture ideas from the earlier literature pass
  (`GreenKuboSpectralDensityResult`, `PseudocontactShiftResult`,
  per-SS CSA stratification, SE(3) GP, volumetric FNO).
- `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md` — consolidated
  library-wide pass of ~40 calculator / I/O / diagnostic ideas.
- `spec/pending_decisions_20260423.md` — `AllWelfords` revival,
  `ChiRotamer` attachment, `FullFat` MOPAC deps.
- `spec/TRAJECTORY_WRITE_SURFACE_2026-04-24.md` — compile-time
  enforcement of the buffer discipline this plan depends on
  (sibling doc; read before implementing).

---

## 1. Framing

The thesis threads two frames at once:

- **Advisor framing:** *NMR learning system.* The equivariant model
  (Stage 3, `learn/c_equivariant/`) is the ML deliverable and the
  driver of a lot of the extractor design.
- **Department-head framing:** *"I want physics, not a damn
  bioinformatics paper."* The kernel outputs must be physically
  justified and their observable predictions must survive
  experimental comparison.

The 2-protein × 15 ns × ~375-DFT-frames dense validation set is the
bridge:

- **Physics validation** — kernel-derived observables (S², T₁, T₂,
  CCR, NOE, J(ω)) are compared against published experimental
  relaxation data on the chosen pair. Independent ground truth.
- **ML training substrate** — ~750 paired (kernel features, DFT σ
  tensor) points at dense cadence, high-quality supervision for the
  equivariant model.
- **Exploration** — paired (kernel, DFT) per frame reveals structure
  in σ(t) autocorrelation, tensor modes, etc., that the 685-fleet
  cannot (the fleet has sparse DFT only).

The 685-protein fleet becomes **coverage** — Layer 1–3a outputs across
many proteins, answering how broadly kernels carry trajectory-scope
signal. The 2-protein dense set is the **validation bed** and primary
thesis claim carrier.

Plan assumes all six observables land on the pair. No hedges on partial
experimental coverage.

---

## 2. Object-model grounding (non-negotiable)

The discipline every TR in this plan honours (PATTERNS.md §§1-18,
OBJECT_MODEL.md, feedback memory
`feedback_object_model_scope_discipline`):

| Buffer | Role | What lives there |
|---|---|---|
| `Protein` | Invariant definition | Identity, topology, residues, bonds, rings, cached backbone + chi indices, `ProteinBuildContext`. Const after `FinalizeConstruction`. |
| `ProteinConformation` + `ConformationAtom` | **Instant-observation buffer** (universal across crystal, NMR, MD frame, prediction, minimised) | Positions (const) + per-atom fields written by `ConformationResult`s (one writer per field). Works identically for a static PDB or an MD frame. |
| `ConformationResult` | Per-instant calculation | Reads `ProteinConformation`, writes `ConformationAtom` fields (+ CR-internal state). Singleton per type per conformation. |
| `TrajectoryProtein` | Protein in trajectory context | Wraps invariant Protein, seats conf0, owns per-atom `TrajectoryAtom` vector and attached TRs. Holds trajectory-scope invariant topology derivatives (`ChargeSource`, `BondedParameters`, `FullSystemReader`). |
| `TrajectoryAtom` | **Time-spanning per-atom store** | Typed rollup fields (one writer each, AV pattern) + per-atom event bag. No accumulator state. |
| `TrajectoryResult` | Calculation over time | Reads per-frame `ProteinConformation` (const), accumulates internally, lands rollup on `TrajectoryAtom` (AV) or transfers `DenseBuffer<T>` (FO). |
| `Trajectory` | Run process + record | Source paths, preloaded EDR, frame times/indices, run-scope selection bag, per-frame env stash. |
| `RunConfiguration` | Per-run shape | Factory list, required CR set, stride, AIMNet2 requirement. |

**The rule, user's words (2026-04-24):** *"If it has to do with evaluating
a frame or PDB, it goes in an old-style calculator and ConformationAtom
or the ConformationResult holds it. Then the trajectory scoops that up
and puts it in TrajectoryAtom or TrajectoryResult or in an appropriate
place in TrajectoryProtein."*

**The static-PDB-preserving test.** Before adding any per-frame
computation, ask: *if a user loaded this protein as one crystal
structure with no trajectory, could they still get this quantity?* If
not, it's in the wrong scope.

CR work and TR work live in **separate sessions, separate tests**. The
duplication — CR writes a CA field, TR reads it and writes its own
trajectory-scope output — is the discipline, not waste.

---

## 3. Layer plan

Each layer unlocks the next. Within a layer, items can ship independently
in any order.

### Layer 0 — TR test foundation (blocks everything else)

**Status 2026-04-24 evening:** BondLengthStats discipline tests landed
inline in `tests/test_gromacs_streaming.cpp` (commit `e306f2a`):
`BondLengthStatsEndToEnd` (integration) + `Frame0Semantics` +
`FinalizeIdempotency` + `H5RoundTrip` (the three-assert discipline).
BsWelford discipline tests are still pending — ship alongside the
first Layer 1 clone (`HmWelfordTrajectoryResult`) so the pair
demonstrates the clone recipe against both a no-CR-dep TR and a
BS-dep TR.

Template asserts for every new TR test: frame-0 semantics (AV: valid
after one Compute, stride ≥ fixture length isolates; FO: internal
buffer non-empty), Finalize idempotency, H5 round-trip via temp file.

**File organisation.** Inline in `test_gromacs_streaming.cpp` for now;
each TR's discipline tests group under a shared name prefix
(`BondLengthStats*`, then `BsWelford*`, then `HmWelford*`, …).
Extraction to per-TR-class files (original plan) is deferred until
the inline count grows enough to justify splitting — probably around
Layer 1 halfway mark (~4 Welford TRs in the tree).

### Layer 1 — Per-kernel Welford shelf

Clone `BsWelfordTrajectoryResult` once per existing non-MOPAC
`_shielding_contribution` source. Each TR ~150 lines, one
`Dependencies()` line, one source-field line, one schema-attr block.
Adds ~18 typed rollup fields to `TrajectoryAtom` per kernel.

| New TR | Source field on `ConformationAtom` | Per-frame CR dep |
|---|---|---|
| `HmWelfordTrajectoryResult` | `hm_shielding_contribution` | `HaighMallionResult` |
| `McWelfordTrajectoryResult` | `mc_shielding_contribution` | `McConnellResult` |
| `CoulombWelfordTrajectoryResult` | `coulomb_shielding_contribution` | `ApbsFieldResult` |
| `HBondWelfordTrajectoryResult` | `hbond_shielding_contribution` | `HBondResult` |
| `PiQuadWelfordTrajectoryResult` | `piquad_shielding_contribution` | `PiQuadrupoleResult` |
| `RingSuscWelfordTrajectoryResult` | `ringchi_shielding_contribution` | `RingSusceptibilityResult` |
| `DispersionWelfordTrajectoryResult` | `disp_shielding_contribution` | `DispersionResult` |
| `TotalGWelfordTrajectoryResult` | `total_G_spherical` | sum chain |
| `AimnetPredictedWelfordTrajectoryResult` | `predicted_T0` + `predicted_T2[5]` | `AIMNet2Result` |

Each reads existing CA fields — **no CR changes, no CA changes**. MOPAC
variants (`mopac_coulomb_*`, `mopac_mc_*`) park in the catalog, not cut;
they ship after the D5 reconciliation (see §5).

Decision to land during Layer 1: `AllWelfords` revival
(`pending_decisions_20260423.md` item 1). Once 3–4 Welford TRs are in,
the per-TR hand-spelled field-name lists start to drift. A single
source-of-truth enumeration on `TrajectoryAtom` prevents it.

### Layer 2 — Convergence gate

**`BlockAveragedConvergenceResult`** (time-series doc idea 1,
Grossfield & Zuckerman 2009). Per-atom per-tensor-component block SEM
vs block length; "effective independent sample count" per atom. FO,
`DenseBuffer<double>`.

Gates:
- every Layer 1 Welford claim (is the mean converged?);
- every Layer 3 ACF (is correlation significance meaningful?);
- every Layer 4 observable fit (does the τ fit make sense?).

At ~375 frames on the 2-protein set, the convergence question has
actual teeth.

### Layer 3 — Autocorrelation family

ACF family discipline lives in
`spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md` §"Autocorrelation
family discipline (ultrathink)". Three families, distinct schema:

**3a. Linear scalar ACFs** (`*AutocorrelationTrajectoryResult`,
`estimator="biased"`, `mean_convention="full_range"`). One per scalar
source on CA — clone `BsT0AutocorrelationTrajectoryResult`:

- `HmT0Autocorrelation`, `McT0Autocorrelation`, `CoulombT0Autocorrelation`,
  `HBondT0Autocorrelation`, `PiQuadT0Autocorrelation`,
  `RingSuscT0Autocorrelation`, `DispersionT0Autocorrelation`.
- `*T2MagnitudeAutocorrelation` cousins for each of the above (source:
  `T2Magnitude()` on the same `_shielding_contribution`).

**Family-level tightenings to land BEFORE cloning starts** (single
small edit to `BsT0AutocorrelationTrajectoryResult` + a test file):

1. `max_valid_lag` H5 attribute (zeros beyond this are padding).
2. `constant_atom_count` H5 attribute (flat-signal atoms vs missing).
3. `tests/test_bs_t0_autocorrelation_trajectory_result.cpp` with
   synthetic white noise / sinusoid / constant inputs.
4. One-paragraph header note on the exemplar distinguishing the three
   families (scalar / orientational / cross) so the clone doesn't drift
   schema.

**3b. Orientational P2 family**
(`*OrientationalCorrelationTrajectoryResult`, `form="P2_orientational"`,
`normalization="pairwise"`). Genuinely new shape —
⟨P2(u(0)·u(t))⟩, not a scalar ACF.

First: **`PeptideNHOrientationalCorrelationTrajectoryResult`** — the
canonical Lipari-Szabo S² substrate. Iterates over residues with a
valid amide H (`Residue.H != NONE`), reads
`conf.GeometryData().bond_directions[nh_bond_index]`
(per-frame, already cached by `GeometryResult` — no new CR). Per-residue
history buffer → Finalize computes ⟨P2(u(0)·u(t))⟩ at lag k → DenseBuffer.

Cousins (clone only after first lands and Layer 2 confirms 375-frame
convergence is adequate): `BackboneCAOrientationalCorrelation`,
`RingNormalOrientationalCorrelation`.

**3c. Cross-correlation family** — deferred until a reference motivates
a specific pair.

### Layer 4 — Relaxation observables

Sit downstream of Layer 3, compare to experiment.

- **`GreenKuboSpectralDensityResult`** — per-atom J(ω) at Larmor
  frequencies, by FFT of Layer 3a ACFs. Cheap once ACFs exist. Feed
  R₁/R₂ predictions.
- **`SigmaLipariSzaboResult`** — per-atom (S², τₑ) from the σ ACFs.
  **Open question:** C++ TR emitting fit-parameter DenseBuffer, or pure
  Python post-pass consuming Layer 3a's H5 group? scipy's
  `curve_fit` is the obvious tool and doesn't exist in C++.
  Recommend Python unless there's a good reason.
- **`CrossCorrelatedRelaxationResult`** — CCR rates. Needs σ tensor
  time series (existing `BsShieldingTimeSeries` + Layer 1 clones) and
  bond-vector time series (new FO TR
  `BondVectorTimeSeriesTrajectoryResult` over
  `GeometryResult.bond_directions[·]` + Larmor frequency list).
- **`BenchmarkBackCalculationResult`** — per-bench per-protein
  prediction / experimental / residual. **Open question:** class
  granularity (one per bench, desiderata C7; or one per observable
  type). **Open question:** where experimental reference data lives —
  invariant across frames, so NOT on `ConformationAtom` or
  `TrajectoryAtom`; belongs on `Protein` via an `ExperimentalReference`
  attached by a loader (desiderata C2). See §5.

### Layer 5 — Non-kernel diagnostic TRs

No CR deps beyond positions + existing CRs. Cheap, broadly useful,
parallelisable with Layer 1–3:

| TR | Shape | Source |
|---|---|---|
| `RamachandranPhiPsiTimeSeriesTrajectoryResult` | FO, per-residue (φ,ψ) time series | positions + N/CA/C indices |
| `ChiTimeSeriesTrajectoryResult` | FO, per-residue χ[0..3] time series | `Residue.chi[k]` + positions |
| `DsspStateHistogramTrajectoryResult` | AV per-residue 8-bin histogram | `DsspResult` |
| `SasaWelfordTrajectoryResult` | AV per-atom SASA rollup | `SasaResult` |
| `HydrationGeometryWelfordTrajectoryResult` | AV per-atom rollup | `HydrationGeometryResult` |
| `HBondFluxTrajectoryResult` | AV count + make/break events | `HBondResult` |
| `BondedEnergyWelfordTrajectoryResult` | AV per-atom rollup | `BondedEnergyResult` |

Context channels the thesis uses to caption every shielding signal.

### Layer 6 — Conformation-scope calculators (not TRs)

Parallel pipe. Each, once landed as a CR, immediately unlocks Layer 1
+ Layer 3a clones for its new fields.

- **`CSAPrincipalAxisResult`** (desiderata A1) — diagonalise per-atom σ
  tensor, store principal components + PAS orientations. Blocks
  the D1 per-SS CSA split diagnostic (Yao-Bax 2010 α/β 11 ppm target).
- **`PseudocontactShiftResult`** (planned A5) — same `K_ab` kernel as
  `RingSusceptibilityResult`, lanthanide source + external χ tensor
  input. Validation bench: Tharayil 2021 ubiquitin / GB1.
- **`NICSProbeEvaluator`** (desiderata A4) — evaluate kernels at
  non-atom probe positions; unifies library with viewer's
  `QtBiotSavartCalc` / `QtHaighMallionCalc`.

### Layer 7 — Exploratory (park until earlier layers speak)

Time-series doc ideas 3–7 (PCA, Tucker, MSM, memory kernel, Mori-
Zwanzig). Python-side mostly. Wait for Layer 1–3 data.

---

## 4. MOPAC-family TR sequencing

Per `project_mopac_trajectory_sequencing_20260424.md`:

- **Keep existing MOPAC CRs live and working as-is.** No refactor.
- **MOPAC-family TRs on the list, not cut.** `MopacCoulombWelford`,
  `MopacMcConnellWelford`, `MopacMcT0Autocorrelation`, etc. follow the
  same clone recipes as Layer 1 / Layer 3a.
- **D5 reconciliation gates the build decision.** On the 2-protein
  dense set, DFT + ff14SB + MOPAC are all per-frame paired. If ff14SB
  carries the MOPAC signal adequately on that bed, MOPAC-trajectory
  stays parked. If MOPAC diverges in a structured way, that's a thesis
  finding and the MOPAC-trajectory TRs earn the build time.
- **Pending decision 3** (`FullFat` MOPAC CR deps) resolves
  concurrently with the first MOPAC-family TR landing, not before.

User framing (2026-04-24): *"10 MOPAC frames takes longer than 1 DFT."*
Every MOPAC-trajectory design question lives in that frame.

---

## 5. Pre-build decisions to land

Items that cheapen everything downstream and are painful to retrofit:

1. **ACF family tightening** (§3.3a above) — before Layer 3a clones.
2. **`AllWelfords` revival** (`pending_decisions_20260423.md` #1) —
   after Layer 1 reaches 3–4 Welford TRs.
3. **`FullFatFrameExtraction` MOPAC deps** (#3) — land with the first
   MOPAC-family TR.
4. **Experimental reference data home** — for Layer 4
   `BenchmarkBackCalculationResult`. Invariant across frames; object-
   model rule puts it on `Protein` via an
   `ExperimentalReferenceLoader` (desiderata C2). Choose between
   (a) typed table attached to `Protein`, (b) separate
   `ExperimentalReference` object parallel to `ProteinBuildContext`,
   (c) sidecar loaded by `Session`.
5. **Larmor frequency list home** — for J(ω) and CCR. Config, not
   per-atom. On `RunConfiguration` or TOML-loaded via `Session`.
6. **Irrep-tagged H5 metadata** (desiderata C5) — formalise the
   `irrep_layout` / `normalization` / `parity` attrs that
   `BsShieldingTimeSeries` carries, so each `*ShieldingTimeSeries`
   clone writes the same convention.
7. **ChiRotamerSelection attachment** (#2) — sit until μs-MD harvester
   design is au-courant.
8. **`DenseValidationSet` RunConfiguration** — new named factory
   parallel to the existing three, sized for the 2-protein 15 ns
   × 40 ps run with MOPAC on + full classical stack. Not a mutation
   of `PerFrameExtractionSet`.

---

## 6. Open questions for user

Decisions that affect sequencing but don't block other work:

- **Which 2 proteins for the dense validation set?** Candidates with
  published relaxation data: ubiquitin (Lipari-Szabo + Loth 2005 CCR),
  GB1, MLF peptide (de Gortari 2010), lysozyme.
- **Drawn from the existing 685-fleet or new builds?** Affects which
  extraction lineage they inherit.
- **Which experimental observables exist per protein?** S², T₁, T₂,
  CCR, NOE, J(ω)? List drives TR priority within Layer 4.
- **MD-DFT pairing mechanism:** trajectory-concurrent (emit PDB +
  submit DFT per extracted frame during run) or post-hoc (extract
  full trajectory, batch DFT on emitted PDBs)? Current pipeline
  emits PDB snapshots at 1 ns intervals; 40 ps cadence is 25× denser.
- **Do Layer 3b / Layer 4 observable TRs ship on the 685-fleet too,
  or only the 2-protein dense set?** Coverage-vs-depth call. ACF lag
  coverage is already adequate on fleet (625 frames × 40 ps = 25 ns
  window, `N_LAGS=120` reaches 4.8 ns).
- **`BenchmarkBackCalculationResult` class granularity:** one per
  bench (C7 `/benches/<bench>/`) or one per observable type?
- **Experimental-reference data home:** option (a), (b), or (c) in §5
  item 4?
- **Literature pass new calculators:** specifics pending; each one's
  CR extension decided when surfaced.

---

## 7. Sequencing (dependency order, not ranking)

Per the layer plan above, in the order each step unlocks the next
without forcing re-extraction:

1. **Write-surface enforcement** (sibling doc
   `spec/TRAJECTORY_WRITE_SURFACE_2026-04-24.md`) — before any new
   TR. Locks the discipline at compile time for every clone that
   follows.
2. **Layer 0 tests** — `BondLengthStats` then `BsWelford` test files.
3. **ACF family tightening** — four small edits on
   `BsT0AutocorrelationTrajectoryResult`.
4. **First clone at each shape** — `HmWelford` (Layer 1), one
   `HmT0Autocorrelation` (Layer 3a), `PeptideNHOrientationalCorrelation`
   (Layer 3b new shape), `BlockAveragedConvergence` (Layer 2 gate).
   Each lands with its test.
5. **Fan out Layer 1** — remaining 6 Welford TRs. `AllWelfords`
   revival around here.
6. **Fan out Layer 3a** — remaining T0 and |T2| ACFs per kernel.
7. **Layer 5 non-kernel TRs** in parallel with 5/6 as time permits.
8. **Layer 4 observable TRs** — J(ω), Lipari-Szabo (Python?), CCR,
   BenchmarkBackCalc. After Layer 3 lands and the 2-protein set is
   running.
9. **Layer 6 conformation-scope calculators** — `CSAPrincipalAxis`
   first; unlocks diagnostics and new Layer 1/3a clones.
10. **MOPAC-family TRs** — only after D5 on the 2-protein dense set.
11. **Layer 7 exploratory** — only if Layer 2 / Layer 3 show signal
    worth extending.

Each step is mechanically small (~150-line files), testable in
isolation, safe to abandon if data arbitrates away from it. Nothing
forces re-extraction of any in-flight fleet proteins (the 685 fleet
was stopped 2026-04-30 after evaluation showed bad chain extractions;
the 485 not-yet-completed proteins were cancelled, and OF3 is
generating fresh direct-from-sequence structures for the next dispatch).
Nothing is gated on thesis scope — the pipeline adds channels; data
decides which carry signal.
