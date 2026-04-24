# Planned Calculators — Time-Series Extensions (2026-04-24)

**Status:** Note, not spec. Companion document to `PLANNED_CALCULATORS_2026-04-22.md`, surfaced from the 2026-04-23/24 reference-summarization pass (10+ paper summaries; see `references-meta/`). Captured here so the ideas are preserved for extractor-session hand-off. Re-read against the current pipeline and bibliography before implementing. None are committed design.

**Provenance:** all candidates below came out of the per-frame-time-series reframing of Stage 2 (see `proposed_outline_draft.md` Section 4, reframed 2026-04-24) and the reference summaries produced 2026-04-23/24. Each entry names the summaries and `PENDING_ACQUISITIONS.md` entries the extractor session should read before implementing.

**Posture:** these are time-series extensions to the shielding-tensor extraction pipeline. Ideas 1–3 are mechanically implementable with existing data and trivial code; they should ship first. Ideas 4–5 are larger but literature-grounded. Ideas 6–8 are exploratory / high-risk and should only be attempted if 1–5 reveal signal worth pursuing.

**Companion docs:**
- `references/PENDING_ACQUISITIONS.md` Tier 6 (time-series methodology) — papers to acquire for grounding.
- `references-meta/*-summary.txt` — the specific summaries referenced per idea below.
- `references-meta/TICKLERS.md` — two 2026-04-24 entries on Stage 2 validation from Markwick 2010 that relate to these calculators.
- `spec/PLANNED_CALCULATORS_2026-04-22.md` — original five ideas; idea 1 there (`GreenKuboSpectralDensityResult`) is a close cousin of our idea 4 below and should share the `σ(t)` autocorrelation substrate.

---

## 1. `BlockAveragedConvergenceResult` (new, trajectory-level diagnostic)

### What
Per-atom per-tensor-component block-averaged standard error of the mean, as a function of block length. Directly answers the question "is the arithmetic-mean shielding tensor over our 600-frame trajectory converged?"

### Why
At 2 GB of per-frame shielding-tensor data per protein, whether the time series represents independent signal or correlated noise is a genuine empirical question, not a cosmetic one. Grossfield & Zuckerman 2009 (H7 [CORE]) is the standard protocol for answering it for any MD-derived observable. Without this diagnostic, none of our Stage 2 averaging claims can be defended against a careful reviewer. This is the gate on every subsequent time-series claim.

### Type
`TrajectoryResult` field (per-atom × per-tensor-component × per-block-length error estimate), plus a scalar "apparent independent sample count" per atom. Post-trajectory calculator operating on existing per-frame shielding tensors.

### Dependencies
- Per-frame shielding tensors already in the extraction output (exists).
- No new external libraries; `numpy` arithmetic.

### Origin
- **Grossfield & Zuckerman 2009** (H7 [CORE], ⧖) — the protocol. `PENDING_ACQUISITIONS.md` Tier 6.
- Reference-side grounding: Markwick 2010's ¹⁵N 2.89 → 1.84 ppm RMSD gain is only trustworthy if the averaging has converged on both sides of the comparison.

### Risk
Low. Standard technique. The real risk is finding that the 25 ns trajectory yields fewer independent samples than the 600-frame count suggests, which would be a substantive Stage 2 finding, not a bug.

### Status
Ship first. Mechanically trivial; gates every subsequent averaging claim.

---

## 2. `SigmaLipariSzaboResult` (new, trajectory-level)

### What
Per-atom model-free fit of the shielding-tensor autocorrelation function `C_σ(t)` to the standard Lipari-Szabo form `C(t) = S² + (1 − S²)·exp(−t/τ_e)`, yielding a (S²_σ, τ_e) pair per atom per tensor component (or per irrep component after T0/T1/T2 decomposition).

### Why
Lipari-Szabo is the canonical NMR time-series summary: S² measures amplitude of fast motion, τ_e measures its timescale. Applying it to shielding-tensor fluctuations rather than to atomic-coordinate motions (the usual use) gives a direct per-atom summary of *shielding* dynamics, which is what the thesis actually cares about. Cross-comparison against per-atom backbone NH S² from NMR relaxation experiments (where BMRB data exists) becomes a per-residue cross-validation of Stage 2 signal.

### Type
`TrajectoryResult` field (per-atom, per-irrep component). Two scalars per atom per component.

### Dependencies
- Per-atom σ(t) autocorrelation (prerequisite of `GreenKuboSpectralDensityResult` in the 2026-04-22 doc; should be computed once and shared).
- `scipy.optimize.curve_fit` or equivalent.

### Origin
- **Lipari & Szabo 1982** (M1 [CORE], ⧖) — the formalism. `PENDING_ACQUISITIONS.md` Tier 6.
- **Kasinath & Wand 2013** (M7, ⧖) — biological interpretation of S² as entropy meter.
- **Berne & Harp 1970** (M28 [CORE], ⧖) — foundational TCF in MD.

### Risk
The fit assumes a single fast-motion timescale plus a static component. For atoms whose shielding fluctuates on multiple timescales, the simple fit will be misleading. Extended models (Clore-Szabo 1990 extended Lipari-Szabo) exist but add complexity. Start with the two-parameter fit and flag atoms whose residuals exceed a threshold as multi-timescale.

### Status
Ship second, after idea 1 verifies the autocorrelation substrate is converged enough to fit.

---

## 3. `SigmaEssentialDynamicsResult` (new, trajectory-level)

### What
PCA / essential-dynamics decomposition of the per-atom shielding trajectory. Stack the centered (T, N·9) shielding-tensor matrix, SVD across frames, extract top-k modes and their temporal coefficients. Each mode is a spatial pattern (a tensor at each atom) plus a time-course; modes are ordered by variance explained.

### Type
`TrajectoryResult` field holding the top-k modes (k ≈ 20) as (N, 9) spatial patterns plus (T,) time-courses, plus cumulative explained-variance curve.

### Why
Amadei 1993 Essential Dynamics is the standard template for decomposing time-indexed protein data into collective modes. Applied to atomic coordinates it extracts collective conformational motions; applied to shielding, it would extract collective shielding fluctuation modes. If the first few PCs explain most variance, the trajectory data is low-rank, which directly informs whether MSM-style coarse-graining makes sense and whether the tensor-decomposition approach in idea 4 will give interpretable results. If variance is spread across many modes, the signal is high-dimensional and simple decomposition will not reveal structure.

### Dependencies
- Per-frame shielding tensors (exists).
- `numpy.linalg.svd` or `sklearn.decomposition.PCA`.

### Origin
- **Amadei, Linssen & Berendsen 1993** (M19 [CORE], ⧖) — the essential-dynamics formalism. `PENDING_ACQUISITIONS.md` Tier 6.
- Reference-side grounding: Kondor 2025 (✓) on irreducible decompositions — ED is a data-driven analog when the group structure is not known a priori.

### Risk
Interpretation. PCs are data-driven; they may not correspond to physical modes. The diagnostic value comes from comparing PC variance distributions across proteins and checking whether modes align with secondary structure, residue type, or kernel-family contributions. If modes are uninterpretable, report that and move on.

### Status
Ship third. Trivial code, instructive output regardless of interpretation outcome.

---

## 4. `SigmaTuckerDecompositionResult` (new, trajectory-level, exploratory)

### What
Tucker decomposition (higher-order SVD) of the (T, N, K, 9) shielding-tensor trajectory, where T = frames, N = atoms, K = kernels, 9 = Cartesian tensor components (or equivalent T0/T1/T2 irrep encoding). Produces a core tensor plus mode matrices along each axis, each mode admitting interpretation as a temporal, spatial, kernel-family, or tensor-component pattern.

### Why
This is the NYC-taxi-tensor analogy already in Jessica's Stage 2 stats memory: multi-way time-series data with rich mode structure, where Tucker / HOSVD decomposition reveals latent factors that mean-based analysis averages away. Applied to (T, N, K, 9), it asks whether there are low-dimensional spatial-temporal-kernel patterns — e.g. "mode 3 is ring-current-dominated aromatic shielding fluctuating with backbone φ/ψ flips on timescale τ_3." If the decomposition returns such structure, it is a direct physics-interpretable signal from the time-series data that averaging cannot produce. If it returns only statistical regularities without physical meaning, report that and do not over-interpret.

### Type
Post-trajectory Python analysis over combined `TrajectoryResult` outputs for a protein. Not a new C++ calculator. Lives with the Stage 2 analysis pipeline (`learn/` or equivalent).

### Dependencies
- Per-frame shielding tensors (exists).
- Per-frame per-kernel contributions (exists as `shielding_contribution` per kernel).
- `tensorly`, `rTensor`, or equivalent tensor-decomposition library.
- Idea 3 (Essential Dynamics) as a sanity pre-check — if 2D PCA already shows high-dim variance, Tucker will not recover structure.

### Origin
- **Kolda & Bader 2009** SIAM Review 51, 455 (⧖, noted in PENDING Tier 6 adjacent) — canonical tensor-decomposition review.
- Reference-side grounding: Jessica's Stage 2 stats challenge memory (NYC-taxi framing).

### Risk
High interpretation risk, moderate technical risk. Tucker will always return *some* decomposition; the question is whether modes admit physical reading. Must commit up front to a protocol for interpretation (e.g. "top mode along K axis must correlate with a named kernel family or be flagged as statistical") so that results are not cherry-picked.

### Status
Ship after ideas 1–3 if those show the data has structure to decompose. If idea 3 returns many small modes, Tucker is probably not going to reveal much; report negative result and move on.

---

## 5. `CrossCorrelatedRelaxationResult` (new, trajectory-level)

### What
Per-atom cross-correlation between the shielding-tensor CSA correlation and the backbone NH dipolar correlation, yielding the CSA–dipolar cross-correlated relaxation rate at named Larmor frequencies. Complements idea 1 in the 2026-04-22 doc (`GreenKuboSpectralDensityResult`) by extending to multi-observable correlations.

### Why
CCR rates are directly experimentally measurable NMR observables (Tugarinov-Kay TROSY methodology). If our Stage 2 predictions reproduce measured CCR rates for proteins with BMRB CCR entries, that is independent validation of Stage 2 time-series signal. This is the cleanest bridge from our internal (T, N, K, 9) data to an experimentally measured quantity other than shielding itself.

### Type
`TrajectoryResult` field (per-atom, per-Larmor-frequency, per-CCR-observable type).

### Dependencies
- Per-atom σ(t) autocorrelation (shared with idea 2 and with `GreenKuboSpectralDensityResult`).
- Per-atom NH bond vector time series (should already be available from MD trajectory).
- Set of Larmor frequencies (already TOML-configurable per the 2026-04-22 doc).

### Origin
- **Tugarinov & Kay 2003** (M10 [CORE], ⧖) — cross-correlated relaxation in methyl groups, the formalism. `PENDING_ACQUISITIONS.md` Tier 6.
- **Ferrage 2008** — CCR methodology refinement (annotated bib M11, not yet in PENDING).

### Risk
Computable but interpretation depends on MD-sampling completeness at the relevant timescales. Our 25 ns MD is below the window where some protein dynamics live, so CCR predictions will be most reliable for nuclei whose relaxation is dominated by ps–ns motions. Flag scope explicitly.

### Status
Ship after idea 2, which produces the autocorrelation substrate CCR reuses.

---

## 6. `SigmaMSMResult` (exploratory, trajectory-level)

### What
Markov State Model on the shielding trajectory: discretize per-atom σ(t) into kinetic states via k-means or similar, estimate transition matrix via maximum-likelihood on frame-to-frame state transitions, extract slow-mode eigenvectors and implied timescales.

### Why
If per-atom shielding has distinct kinetic substates (e.g. two rotamer-correlated shielding values with slow interconversion), MSM recovers them and estimates the interconversion rate. Complements ideas 2 and 3 at a different scale — Lipari-Szabo fits a single timescale; ED extracts continuous modes; MSM extracts discrete kinetic states. For proteins where rotameric substates dominate the shielding distribution (de Gortari's MLF showed this at 5 μs), MSM is the right framework.

### Type
Python analysis. `PyEMMA` or `deeptime` library dependency. Not a new C++ calculator.

### Dependencies
- Per-frame shielding tensors (exists).
- Commitment to a discretization protocol (k selection, state assignment).
- Idea 3 output as sanity pre-check.

### Origin
- **Prinz et al. 2011** *J. Chem. Phys.* 134, 174105 (⧖, Tier 6 adjacent) — canonical MSM reference.
- **Chodera & Noé 2014** (⧖) — accessible review.

### Risk
Substantial: at 25 ns we may not have enough kinetic transitions to estimate a transition matrix reliably. The MSM literature typically works with μs–ms trajectories. Likely reports "not enough sampling" unless atoms with fast kinetic transitions (sub-ns) dominate. Still worth running as a negative-result exercise.

### Status
Exploratory. After ideas 1–5 and only if they suggest kinetic substates exist.

---

## 7. `SigmaMemoryKernelResult` (exploratory, trajectory-level, high risk)

### What
Mori-Zwanzig projection decomposition of σ(t) into systematic drift, memory kernel K(t), and noise η(t), satisfying the generalized Langevin equation `dσ/dt = drift − ∫ K(t−t')σ(t') dt' + η(t)`.

### Why
If shielding fluctuations have long-range temporal correlations that neither exponential (Lipari-Szabo, idea 2) nor discrete-state (MSM, idea 6) models capture, memory-kernel decomposition is the rigorous non-Markovian treatment. Reveals whether the shielding time series is Markovian (K(t) ≈ δ(t)) or has memory (K(t) with finite width). Little prior work in applied NMR-MD has done this for computed observables — space to explore.

### Type
Python analysis, high-methodology-risk. Not a new C++ calculator.

### Dependencies
- Per-atom σ(t) autocorrelation and higher-order correlations.
- A projection operator choice — the methodological decision point.

### Origin
- **Zwanzig 1961** *Phys. Rev.* 124, 983 (⧖, Tier 6) — memory effects origin.
- **Berne & Harp 1970** (M28 [CORE], ⧖) — MD-context treatment.

### Risk
Very high. Memory-kernel extraction from finite-length time series is numerically delicate. Our 25 ns may be too short to resolve any non-trivial K(t) unambiguously. Commit only if there is strong evidence from ideas 2 and 6 that memory effects exist.

### Status
Exploratory. Probably aspirational for this MSc; flag for PhD-stage if it continues.

---

## 8. Aromatic-H geometry sanity-check extraction (diagnostic, zero new calculator)

### What
Python-side diagnostic: compute our BS / HM kernel output at a set of geometrically constrained aromatic-H probe positions reproducing Agarwal 1977 [10]-paracyclophane (or a benzene-methane probe matching Case 1995), compare against the published numerical values.

### Why
Independent sanity check on whether our `BiotSavartResult.cpp` and `HaighMallionResult.cpp` reproduce known benchmarks. Case 1995's 421 methane-probe shifts and Agarwal 1977's [10]-paracyclophane methylene shifts are public, concrete numerical anchors. If our kernels match these, we have confidence the Stage 2 time-series values are physically correct. If they don't, we have a specific target to debug.

### Type
Python harness, no new C++ calculator. Builds probe geometries in code, calls library kernels, compares to tabulated values.

### Dependencies
- Library kernel evaluators (exist).
- Hand-entered benchmark values from Case 1995 Tables 3 and 4 and Agarwal 1977 Table 1.

### Origin
- **Case 1995** (✓) — methane-probe DFT shift tables.
- **Agarwal 1977** (✓) — [10]-paracyclophane methylene tabulated values.
- See also `references-meta/TICKLERS.md` entry on the JB two-loop separation, which may surface here as a specific numerical test.

### Risk
Low. May find disagreement that traces to the ring-normal-vector convention or the loop-separation parameter (see TICKLERS.md), but that's a discovery, not a bug.

### Status
Ship as part of the idea-1 verification pass. Mechanically trivial once benchmark values are entered.

---

## How to use this document

- Idea 1 is a hard prerequisite for 2–7. Do it first. Idea 8 can run in parallel.
- Ideas 2 and 3 are cheap, instructive, and produce publishable figures regardless of outcome. Ship these before committing to 4.
- Idea 4 (Tucker) is the physics-interpretation ambition; it requires ideas 1 and 3 as prerequisites and real commitment to an interpretation protocol.
- Ideas 5, 6, 7 are increasingly exploratory. Commit only if earlier ideas reveal signal to extend.
- Each idea's `Origin` lines name the papers the extractor session should have read before implementing. Papers marked (⧖) are in `PENDING_ACQUISITIONS.md`. If a paper is not yet on disk when an extractor session reaches its idea, acquire-and-summarise first.

## Acquisition priority list (for Jessica to hustle this week)

**Must-have before implementation:**
1. Lipari & Szabo 1982 — idea 2.
2. Grossfield & Zuckerman 2009 — idea 1.
3. Amadei-Linssen-Berendsen 1993 — idea 3.

**Strongly recommended:**
4. Kolda & Bader 2009 SIAM Review — idea 4.
5. Tugarinov & Kay 2003 — idea 5.

**Helpful for exploratory arms:**
6. Prinz et al. 2011 or Chodera-Noé 2014 — idea 6.
7. Zwanzig 1961 + Berne-Harp 1970 — idea 7.

All are in `references/PENDING_ACQUISITIONS.md` Tier 6. Case 1995 and Agarwal 1977 for idea 8 are already on disk.

## Amendments

### 2026-04-24 — Reframe through kernel-truthfulness (correlate, do not match)

The primary Stage 2 claim is that our classical kernel extractors **correlate** with DFT-computed shielding across several dimensions — not that they **match** DFT point-for-point. The target is r² good enough to reflect that the kernel extractors saw a signal, not numerical agreement. Under this framing the calculator priorities reshuffle:

- **Idea 8 is higher priority than originally written.** Direct ring-current kernel correlation test against Case 1995 and Agarwal 1977 tabulated values. If BS / HM outputs do not correlate with those benchmarks at the published probe geometries, we have a specific debuggable gap before any Stage 2 claim is defensible.
- **Ideas 2, 3, 4 read as kernel-vs-DFT correlation comparisons.** Not "does our S²_σ predict experimental relaxation" but "does kernel-derived S²_σ correlate with DFT-derived S²_σ on protein X." Not "does ED mode 1 match a motion" but "does kernel-ED mode 1 correlate with DFT-ED mode 1." The r² is the validation currency.
- **Ideas 5, 6, 7 are characterisation, not kernel validation.** Useful if the pipeline has room; first to drop if scope tightens.
- **Idea 1 (convergence) is foundational, not validation.** Establishes the statistical basis on which correlation claims are meaningful. Use its output to partition atoms into "trust this kernel-DFT correlation to ε ppm" versus "this atom needs more sampling before a correlation claim makes sense."

Stage 2 claim shape this supports: *"Kernel extractors correlate with DFT-computed shielding in the mean, in the amplitude of dynamic fluctuations, in the structure of principal modes, and in tensor-factor decomposition, to r² = X on protein Y. The kernels saw the signal."*

### 2026-04-24 — Autocorrelation family discipline (ultrathink)

The `BsT0AutocorrelationTrajectoryResult` exemplar's job is to anchor
the autocorrelation family so future authors clone inside the sane
shape. That anchoring works at four levels that deserve explicit
statement before the next ACF TR lands.

#### Physics commitment and its scope

`BsT0Autocorrelation` commits to the **biased estimator with
full-range mean**:

```
μ = (1/T) Σ x_t
C(k) = (1/T) Σ_{t=0..T-k-1} (x_t − μ)(x_{t+k} − μ)
ρ(k) = C(k) / C(0)
```

Declared via H5 attrs `estimator="biased"`, `mean_convention="full_range"`.

The commitment is *physics*, not preference:

- **Wiener–Khinchin non-negative PSD** — the Fourier transform of a
  biased ACF is a non-negative power spectrum. That is what NMR
  relaxation theory (Lipari–Szabo, J(ω) spectral-density mapping)
  actually consumes. Unbiased (1/(T−k)) ACF does not have this
  guarantee.
- **|ρ(k)| ≤ 1 at all finite T** — holds for biased; can fail for
  unbiased at large k.
- Matches MDAnalysis, numpy.correlate (1/N normalisation), the MD-NMR
  literature.

**Scope of this commitment: linear ACF of scalar time series only.**
The moment a TR reaches for a nonlinear form — P2(u(0)·u(t)) for
Lipari–Szabo order parameters, say — it's in a different family with
different normalisation conventions. Different class, different
schema attrs, not a variant of the linear ACF exemplar.

#### The family tree, named

Three distinct families all get called "correlation" in casual usage.
They should be distinct TR classes with distinct naming and distinct
schema attrs:

| Family | Class-name form | Core quantity | Schema attrs |
|---|---|---|---|
| Linear scalar ACF | `*AutocorrelationTrajectoryResult` | ρ(k) = C(k)/C(0) | `estimator="biased"`, `mean_convention="full_range"` |
| Orientational (P2, Lipari–Szabo) | `*OrientationalCorrelationTrajectoryResult` | ⟨P2(u(0)·u(t))⟩ | `form="P2_orientational"`, `normalization="pairwise"` (1/(T−k)) |
| Cross-correlation (two sources) | `*CrossCorrelationTrajectoryResult` | ρ_xy(k) with lag sign conventions | `form="cross_correlation"`, source-pair attrs |

The exemplar currently only demonstrates family #1. When someone
clones `BsT0Autocorrelation` to build a P2 orientational Result, they
must rename and re-attr — not sneak nonlinear physics in under a
schema that claims biased linear.

#### Code pattern — clone-and-adapt, not template

Per PATTERNS §17 addition: duplication preferred over chaining. For
ACF: **one class per scalar source.** Do NOT template
`AutocorrelationTrajectoryResult<ScalarExtractor>`. Do NOT subclass
through a protected base.

What distinguishes one linear-ACF TR from the next is exactly three
things:

1. **`Dependencies()`** — which ConformationResult type provides the scalar.
2. **The one line in `Compute` that names the scalar:**
   ```cpp
   per_atom_history_[i].push_back(conf.AtomAt(i).<source_field>);
   ```
3. **The H5 group name and any source-specific attribute** (e.g.
   `units`, parity if relevant).

Everything else is boilerplate: per-atom history
`vector<vector<double>>`, frame-time recording, Finalize doing the
biased-ACF sum, DenseBuffer transfer, WriteH5Group schema attrs.
Repetition is the discipline — a class is a class, the reader finds
the whole thing in one file.

**Rule of thumb:** if a new ACF TR is significantly more than ~150
lines, either the source is exotic (good — it earns its length) or
it's drifting into abstraction (bad — stop).

#### Uniform handling every ACF TR must carry

Invariants the exemplar already gets right; future authors must not
drop them:

- **ρ(0) = 1** when C(0) > 0. First sanity check at Finalize.
- **|ρ(k)| ≤ 1 for all k.** If any lag violates, the estimator is wrong.
- **Constant signal** (C(0) < ε) → emit zero at all lags. Worth an
  explicit `constant_atom_count` attribute so downstream can
  distinguish "flat signal" from "missing data."
- **Short runs** (T < N_LAGS) → zero-pad lags beyond T. Worth a
  `max_valid_lag` attribute so downstream doesn't misread zeros as
  signal.
- **`sample_interval_ps`** derived as median of consecutive Δt from
  the Result's own `frame_times_` record — robust to stride, and
  robust to test callers that don't populate `Trajectory::FrameTimes()`.
- **Memory**: per-atom `vector<double>` history. O(N·T·8) bytes. At
  25 ns × 625 frames × 4000 atoms ≈ 20 MB — comfortable. For μs
  scaling a bounded circular-window + per-lag running left/right
  tail-sum refactor is a swap-in; documented in the exemplar header.

#### Test pattern the family needs

Testing an ACF member should assert:

1. **Synthetic white-noise input** → ρ(0)=1, ρ(k→0) rapidly for k>0.
2. **Synthetic sinusoid** at known frequency → ρ(k) reproduces expected cosine.
3. **Constant input** → ρ = 0 at all lags (covers the "can't divide by zero" path).
4. **|ρ(k)| ≤ 1 at every k for every atom.**
5. **H5 round-trip** — write, read, bit-identical.
6. **`sample_interval_ps_`** equals the known Δt for a synthetic `frame_times_`.

These assertions do not need a protein or ConformationResult chain —
they need a `TrajectoryProtein` fixture with a tiny synthetic
trajectory and a test-only mechanism to feed scalars. The first ACF
test should probably be a **synthetic-input variant** of the test
rather than "run BsT0 on 1UBQ and inspect output": the former tests
the algorithm, the latter just confirms the wiring works.

#### Candidates sitting on the shelf (no ranking; neighbourhood map)

Each is a ~150-line file. Which get built is data-arbitrated; this
is the neighbourhood:

- **Other-classical-kernel linear ACFs** (~one file each):
  `HmT0Autocorrelation`, `McConnellT0Autocorrelation`,
  `RingSusceptibilityT0Autocorrelation`, `PiQuadrupoleT0Autocorrelation`,
  `CoulombT0Autocorrelation`, `HBondT0Autocorrelation`,
  `DispersionT0Autocorrelation`, `ApbsT0Autocorrelation`. Plus the
  `*T2Magnitude*Autocorrelation` cousins. Each is clone + one
  source-field line.
- **Pair / bond-scope linear ACFs:** `BondLengthAutocorrelation`,
  `McConnellDipolarTensorT0Autocorrelation` (per-bond-neighbourhood).
- **Scalar-derived from vector/tensor:**
  `CoulombEfgT2MagnitudeAutocorrelation` → quadrupolar relaxation
  proxy; `ApbsEfgT2MagnitudeAutocorrelation`.
- **Orientational family (different schema!):**
  `PeptideNHVectorOrientationalCorrelation` → Lipari–Szabo S² proxy;
  `BackboneCAVectorOrientationalCorrelation`;
  `RingNormalOrientationalCorrelation` → ring-tumbling contribution
  to ring-current shielding.
- **Cross-correlation family (different schema!):** paired
  bond–shielding, paired atom–atom — only if a reference proposes a
  specific pairing strategy.

#### Small decisions to land before cloning starts

Durable family-level tightening the exemplar doesn't yet force:

- **`max_valid_lag` attribute** added to the emission — zeros beyond
  this are padding, not signal.
- **`constant_atom_count` attribute** — lets downstream quickly see
  how many atoms had flat signal.
- **A test file** (`test_bs_t0_autocorrelation_trajectory_result.cpp`)
  that exercises a synthetic signal path, so the next ACF clone has
  a template to adapt.
- **A one-paragraph header note on the exemplar** explicitly stating
  that P2 orientational and cross-correlation are *cousin families*,
  not variants, so the clone doesn't drift scope.

None of this requires redesign of the object model or buffer
ownership. The discipline is in the schema attrs, the
class-per-source rule, and the edge-case handling each clone must
carry.
