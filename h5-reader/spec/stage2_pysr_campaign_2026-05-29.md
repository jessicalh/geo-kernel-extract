# Stage 2 — comprehensive PySR rediscovery campaign (best-case trajectory examination)

Date: 2026-05-29
Status: Stage 2 plan, settled architecture after 2026-05-29 conversation.
Supersedes prior speculative sketches. Architecture is "one walk, many
models" through a unified h5-reader platform.

## Project stage framing

- **Stage 1 (settled, 2026-04-10)**: per-element, per-atom-type ridge
  on 720 proteins / 446K atoms / 55 kernels. R² = 0.818. The mutation
  set + the static-snapshot calibration.
- **Stage 2 (this doc)**: **best-case trajectory examination on 1P9J.**
  Comprehensive PySR symbolic-equation rediscovery campaign across
  8-10 kernels in three physics regimes, with charge-source sensitivity
  + polarisability sensitivity studies. Built on the deeply-instrumented
  single-protein trajectory (15 ns, 750 DFT frames, 10 ps NPY intervals,
  MOPAC + AIMNet2 + full kernel set).
- **Stage 3 (future)**: model evaluation against experimental BMRB
  shifts; generalisation across the fleet; expert-system / narrative-
  engine deferred items per the design family.

Stage 2 is the *intellectual core* of the thesis — physics-based ML that
rediscovers and quantifies the canonical NMR shielding mechanisms from
the best-case trajectory data available. Stage 3 generalises; Stage 2
demonstrates.

## Why comprehensive, not batched

Per the user's direction (2026-05-29): **don't pace artificially. Take
the full campaign. Disappointments are findings.** Compute is plentiful;
parallelism is free; per-kernel analysis is the bottleneck. Pacing risks
scope creep. "We rediscovered N canonical mechanisms across three physics
regimes" beats "we rediscovered ring current and ran out of time."

---

# The architecture: one walk, many models

The platform is the h5-reader app, extended with a typed-substrate batch
mode that runs Python data-transform scripts during a per-protein walk
and writes per-transform output files. A separate Python SDK consumes
those outputs after the walk to fit models or do analysis.

**Only two applications to maintain: extract and reader.** The reader
binary handles both interactive vetting (today's behavior) and batch
walks (the new capability). One source tree, one core library.

## The pattern

```
┌─────────────────────────────────────────────────────────────────┐
│  Per-protein walk (expensive; do it as rarely as possible)      │
│                                                                 │
│  ┌──────────────┐    ┌──────────────────────────────┐           │
│  │ load protein │ → │ C++ harness iterates atoms ×  │           │
│  │ (typed bufs) │    │ frames; for each (atom, t),  │           │
│  │              │    │ calls every registered       │           │
│  │ build per-   │    │ transform with the rich      │           │
│  │ frame index  │    │ substrate available          │           │
│  │ (KD-tree,    │    │                              │           │
│  │  ring centr, │    │   for each transform T:      │           │
│  │  H-bond cand)│    │     if T.filter(atom):       │           │
│  │              │    │       row = T.transform(...) │           │
│  │ compute win- │    │       T.sink.append(row)     │           │
│  │ dow deriv-   │    │                              │           │
│  │ atives (vel, │    │ at end-of-protein:           │           │
│  │ kernel rate) │    │   each T.sink → Parquet      │           │
│  └──────────────┘    └──────────────────────────────┘           │
│                                                                 │
│  Output: one .parquet per transform per protein                 │
└─────────────────────────────────────────────────────────────────┘

                            ↓ (separate process, days/weeks later)

┌─────────────────────────────────────────────────────────────────┐
│  Per-model SDK scripts (cheap; rerun as desired)                │
│                                                                 │
│  python sdk/pysr_bs_aromatic_he.py outputs/bs_aromatic_he.parquet │
│  python sdk/ridge_baseline.py     outputs/ridge_features.parquet │
│  python sdk/mace_train.py         outputs/mace_features.parquet  │
│  python sdk/symbolic_alt.py       outputs/bs_aromatic_he.parquet │
│                                                                 │
│  Output: model coefs / Pareto frontiers / R² tables / etc.      │
└─────────────────────────────────────────────────────────────────┘
```

## Python is a data transform, not an orchestrator

Inside the walk, Python transforms are **pure data shaping**: given
typed access to the loaded protein + topology + per-frame derived data,
emit one or more rows per (atom, frame) into a named output sink. No
PySR, no sklearn, no MACE, no Julia subprocesses during the walk.

The rich substrate is available — the transform doesn't have to compute
its own distances or build its own indices. The C++ side provides:

- `frame.topology.nearest_aromatic_ring(atom)` → `RingGeometry` struct
  (distance, normal angle, ring atom count, ring class)
- `frame.topology.nearest_h_bond_acceptor(atom)` → `HBondGeometry` struct
- `frame.topology.neighbors_within(atom, radius)` → filtered atom list
- `frame.topology.charged_neighbors(atom, radius)` → `[(atom, charge)]`
- `frame.kernel(name, atom)` → typed `KernelValue` (T0/T2 components,
  dia/para parts for DFT shielding)
- `window.atom_velocity(atom)` → centered finite difference Vec3
- `window.kernel_rate(name, atom)` → finite difference scalar
- `window.distance_rate(atom_a, atom_b)` → pair-distance derivative
- `window.ring_approach_rate(atom, ring)` → distance-to-ring derivative
- Atom predicates: `atom.is_aromatic_HE()`, `atom.is_backbone_HN()`,
  `atom.is_sidechain_carboxylate()`, etc.
- `atom.residue`, `atom.element`, `atom.atom_name`, `atom.id`, etc.

A transform looks like this:

```python
# transforms/bs_aromatic_he.py
class Transform:
    REQUIRED_KERNELS = ['BS']          # harness loads only this kernel TR
    REQUIRED_WINDOW = 1                # static features; no rates needed

    def filter(self, atom):
        return atom.is_aromatic_HE()

    def schema(self):
        return {
            'protein_id': 'str', 'frame': 'int32', 'atom_id': 'int32',
            'residue_name': 'str', 'residue_index': 'int32',
            'r_ring': 'float64', 'cos_theta': 'float64',
            'ring_atom_count': 'int32', 'ring_class': 'str',
            'target_BS_T0': 'float64',
        }

    def transform(self, frame, atom):
        ring = frame.topology.nearest_aromatic_ring(atom)
        if ring is None or ring.distance > 8.0:
            return None
        return {
            'protein_id': frame.protein_id,
            'frame': frame.index,
            'atom_id': atom.id,
            'residue_name': atom.residue.name,
            'residue_index': atom.residue.index,
            'r_ring': ring.distance,
            'cos_theta': math.cos(ring.normal_angle),
            'ring_atom_count': ring.ring_atom_count,
            'ring_class': ring.classification,
            'target_BS_T0': frame.kernel('BS', atom).T0,
        }
```

That's the entire transform. No model code. No accumulators. The harness
collects emitted rows; at end-of-protein, writes one Parquet per transform.

## The killer facility: one walk → many models

You can only afford to walk 40K atoms × 1500-20000 frames so many times.
Per-protein, the walk is the expensive thing — load + index + iterate.

The platform makes this cost amortizable: **register N transforms in one
walk, each emitting its own output**. Each transform can be a different
feature set for a different downstream model.

```bash
h5reader --batch \
  --protein /path/to/1p9j.h5 \
  --transforms transforms/bs_aromatic_he.py \
               transforms/mc_backbone_hn.py \
               transforms/apbs_efg_charged_sidechains.py \
               transforms/larsen_hbond_donors.py \
               transforms/ridge_baseline_features.py \
               transforms/mace_descriptor_features.py \
               transforms/aimnet2_embedding_capture.py \
               transforms/charge_source_comparison_ff14sb.py \
               transforms/charge_source_comparison_mopac.py \
               transforms/charge_source_comparison_aimnet2.py \
  --output-dir runs/1p9j_2026-06-XX/
```

One protein load. One traversal. Ten Parquet outputs. The walk pays
once and feeds N downstream models.

Then, after the walk, each model gets its own SDK script:

```bash
python sdk/pysr_run.py     runs/1p9j_2026-06-XX/bs_aromatic_he.parquet
python sdk/pysr_run.py     runs/1p9j_2026-06-XX/mc_backbone_hn.parquet
python sdk/pysr_run.py     runs/1p9j_2026-06-XX/apbs_efg_charged_sidechains.parquet
python sdk/pysr_run.py     runs/1p9j_2026-06-XX/larsen_hbond_donors.parquet
python sdk/ridge_baseline.py runs/1p9j_2026-06-XX/ridge_baseline_features.parquet
python sdk/mace_train.py    runs/1p9j_2026-06-XX/mace_descriptor_features.parquet
python sdk/symbolic_alt.py  runs/1p9j_2026-06-XX/bs_aromatic_he.parquet   # try a different SR alg
python sdk/charge_compare.py runs/1p9j_2026-06-XX/charge_source_comparison_*.parquet
```

Each SDK script is independent. Rerun with different hyperparameters
without touching the walk. Run a new SDK script against an old walk's
output months later. The walk-output Parquet IS the typed contract.

## Two stages, cleanly separated

| Stage | Where | What | When |
|-------|-------|------|------|
| **Walk** | C++ harness, Python transforms | Pure data shaping per (atom, frame); emit rows | Once per protein, runs ~minutes to ~hour |
| **Modeling** | Standalone Python scripts | PySR / sklearn / PyTorch / custom over the Parquets | Whenever you want a result, against any walk's output |

The separation buys:

- **No ML during the walk.** No PySR Julia subprocess to coordinate, no
  PyTorch GPU memory to manage, no GIL-juggling around model fit calls.
  The walk's job is to produce data; that's all it does.
- **Reproducibility.** The Parquet outputs are the canonical input for
  every model. Re-running PySR with a new hyperparameter doesn't require
  re-walking. The walk is the expensive step; protect it.
- **Many models from one walk.** As above.
- **The SDK scripts are tiny.** ~30-50 lines each: load Parquet, slice
  columns, call the model library, write results. No infrastructure to
  reinvent per model.
- **Independent failure modes.** A bad model script doesn't corrupt the
  walk's output. A walk that finishes successfully produces durable
  artifacts even if every downstream model script crashes.

## "It doesn't get finished until the run is done"

The walk runs to completion (or user interrupts it). Either way, what
gets written to disk is consistent — each transform's Parquet sink is
buffered and committed at protein boundaries (or atom-batch boundaries
if streaming is needed). No half-written state.

If the user presses Done mid-walk, the harness finishes the current
frame, closes each transform's sink at whatever rows it had, exits
cleanly. The partial Parquets are still consumable by the SDK; just
smaller. No corruption window, no need for crash-recovery infrastructure.

---

# The unified app: extract and reader

Two applications: extract (producer) and reader (consumer + analysis
platform). The reader binary handles both interactive vetting (today's
behavior) and batch walks (the new capability). One source tree.

## What's already in h5-reader

The substrate Stage 2 needs is largely already there:

- **`QtTrajectoryH5`** — eager bulk loader; reads positions + all
  present TR groups into typed flat buffers in one pass. **This is the
  walk's preloaded substrate.**
- **`TrajectoryConformation::frame(t)`** — returns `QtFrame`, a value-
  type view over the buffers. Maps onto the transform's per-frame
  context.
- **`QtProtein` / `QtTopology`** — typed identity model with ring class
  hierarchy + virtual methods (`IsAromatic()`, etc.).
- **`QtProteinLoader::LoadRunPath()`** — does the producer-run-directory
  loading; reusable as-is.
- **The `QtRing` hierarchy, residue typing, AMBER topology decoding**,
  etc. — every typed primitive the transforms need is there.

## What's added (in `src/analysis/`)

- `QtSpatialIndex` — per-frame KD-tree (vendored nanoflann), lazy +
  cached per frame
- `QtWindowDerivatives` — finite differences across the preloaded
  trajectory
- `QtTopologyQueries` — `nearest_aromatic_ring`, `nearest_h_bond_acceptor`,
  `neighbors_within`, etc., composed over the KD-tree + ring data
- `QtTransformHarness` — the walk loop; iterates protein × frame ×
  atom, calls each registered Python transform via pybind11, collects
  emitted rows into per-transform sinks
- `QtTransformLoader` — given a path to `transforms/foo.py`, instantiates
  the `Transform` class, validates its contract (`REQUIRED_KERNELS`,
  `schema()`, `filter()`, `transform()`)
- `QtParquetSink` — buffered Arrow C++ writer, one per transform
- `bindings.cpp` — pybind11 exposing the substrate to Python
- `main_batch.cpp` — the `--batch` entry point

## The CMake shape

One `h5reader-core` library, one `h5reader` binary:

```
h5reader-core (STATIC lib)
├── src/model/*               — QtProtein, QtTopology, QtFrame, etc.
├── src/io/*                  — QtTrajectoryH5, QtProteinLoader
├── src/diagnostics/*         — logging, crash capture
├── src/analysis/*            — new: substrate primitives + harness
└── links: Qt6::Core, Eigen3, HDF5, HighFive, Arrow C++

h5reader (executable)
├── src/main_reader.cpp       — default mode: GUI
├── src/main_batch.cpp        — --batch mode: harness
├── src/app/*                 — widgets, scenes, overlays, REST
├── links: h5reader-core, Qt6::Widgets, VTK, pybind11::embed,
│         Python3::Embed
```

One binary because:

- The same code paths read trajectories whether for display or for
  transform iteration
- Movie capture (future) will want VTK in batch mode anyway — no point
  splitting it out now
- "Two applications to maintain" means extract + reader; not extract +
  reader-gui + reader-batch
- Mode selection is a CLI flag, not a separate build target

VTK ships with the batch binary. Headless render on scan boxes via
xvfb (already used in CI) or VTK offscreen rendering. No new system
dep on the scan fleet.

## Threading during the walk

Simpler than the earlier framing because there's no ML inside the walk:

- **Main thread**: `QApplication` (interactive mode) or `QCoreApplication`
  (batch mode). Owns Qt event loop, widgets/VTK in interactive mode.
- **Walk worker thread** (`QThread`): runs the harness loop. Holds the
  GIL only around Python `transform()` calls; releases it for C++
  preprocessing (KD-tree builds, finite differences, kernel reads).
- **Shared data**: preloaded `QtTrajectoryH5` buffers + per-frame
  indices are read-only after load. Multiple readers, no locks.

The interactive case (user scrubbing while a walk runs) gets a free win
from preload: visualization reads frame M, walk reads frame N, both
read from the same immutable buffers. No coordination needed.

In batch-only mode the worker thread is the only thread doing work;
the main thread services `finished` signals and exits. Standard Qt
headless worker pattern.

---

# Optional: interactive integration

The user's framing: "we can try it if we want, and probably should, and
it probably works and maybe we keep the original vision, but either way
it is good."

Two patterns can be added on top of the basic walk-and-Parquet platform:

**Pattern A (deferred): live walk results in the GUI.** A dock panel
that lets the user pick a `transforms/` directory and run a walk against
the currently-loaded protein, displaying progress + emitted-row counts
in real time. Output Parquets land on disk; user can also choose
"in-memory" sinks to feed the strip stack as new dynamic signals. No
PySR integration — model fitting still happens via SDK scripts in a
separate process. The dock is a convenience for "run a transform against
this loaded data without leaving the app."

**Pattern B (further deferred): in-app modeling.** A second dock that
runs an SDK script against a just-finished transform's output and
overlays the result on the 3D scene (PySR equation predictions colored
on atoms, ridge residuals as a strip signal, MACE predictions as a
sequence bar). Genuinely powerful, genuinely orthogonal to the platform's
core function. Earns its weight only after batch walks are routine.

Neither pattern is on the Stage 2 critical path. Both are easy to add
once the batch walk is solid because they reuse the same Python
transform plumbing.

---

# Equation forms — instantaneous, lagged, averaged, differential, spectral, per-trajectory

The prechewed substrate exposes time as a navigable dimension: Python
can call `traj.time_series(atom, kernel)` and get the full per-atom
trajectory as a zero-copy numpy view (atom-major layout means time
series for one atom is contiguous in memory — fast iteration). This
expands the hypothesis space well beyond instantaneous fits.

**Note on what's physically primary.** Under the Born-Oppenheimer
approximation that DFT shielding uses, σ is an instantaneous functional
of geometry — `σ(t) = σ(geometry(t))`. The instantaneous form is the
physics-grounded form; differential and other forms are derived
parameterizations (e.g. `dσ/dt` is the chain-rule expansion of the
instantaneous form). They're valuable as cross-checks and as probes of
where BO might break down (memory terms appearing where they shouldn't),
not as more-fundamental alternatives.

For each kernel rediscovery, the campaign tests multiple equation
forms:

- **Instantaneous (primary)**: `σ(t) = f(r(t), θ(t))`. Pople's classic
  form. The physics-grounded target.
- **Lagged**: `σ(t) = f(r(t−τ), θ(t−τ))`. Probes whether shielding
  responds with delay. Under BO should reduce to instantaneous; a
  nonzero lag is diagnostic of DFT convergence artifacts or
  unaccounted slow modes.
- **Smoothed / dynamics-averaged**: `σ(t) = f(⟨r⟩_τ(t), ⟨θ⟩_τ(t))`.
  The Wishart-Wüthrich claim that experimental shifts reflect
  dynamics-averaged geometry — currently a literature citation; with
  trajectory access it's a PySR experiment.
- **Differential / SINDy-style cross-check**: `dσ/dt = f(σ, r, dr/dt)`.
  Should match the chain-rule expansion of the instantaneous form.
  Memory terms appearing would be a finding (BO breakdown evidence).
- **Per-trajectory descriptors**: each atom emits ONE row with
  `(⟨r⟩, σ(r), autocorr_time(r), ⟨σ⟩, σ(σ), peak_freq(σ), ...)` —
  PySR fits across atoms on trajectory-level features.
- **Spectral**: `σ̃(ω) = g(r̃(ω))`. Frequency-domain coupling. Could
  reveal which vibrational modes dominate shielding response.
- **Per-atom curve fitting**: each atom traces a curve in (r, σ)
  space across its 750 frames. PySR fits ONE universal curve across
  all aromatic-HE atoms (universality test) OR per-atom curves with
  shared form but per-atom coefficients (atom-class structure).

For each kernel, the same walk produces BOTH per-frame Parquets
(instantaneous fits) AND per-atom Parquets (trajectory fits). The
walk visits the same data; transforms differ only in what they
extract.

The thesis claim sharpens from "we rediscovered ring current" to "we
rediscovered ring current in its INSTANTANEOUS form (R² = X),
DYNAMICS-AVERAGED form (R² = Y), LAGGED form (R² = Z; lag ~τ ps
suggesting Q), and SPECTRAL form (frequency ω dominates).
Comparison across functional forms tells us which formulation best
matches DFT — informing which formulation experimentalists should use
when comparing to ssNMR or solution NMR data."

The 5 new TRs (iRED, KernelDynamics, ReorientationalDynamics,
DihedralAutocorrelation, KernelCoherence) are the trajectory-level
features that feed these richer fits; their utility is demonstrated
by the equations PySR discovers using them.

No substrate change required — it already preloads everything;
`traj.time_series(atom, kernel)` is a zero-copy slice. Python uses
numpy/scipy for autocorrelation, FFT, spline interpolation, lag
analysis.

---

# What gets walked (the transforms)

The 8-10 kernel rediscoveries become 8-10 transforms, each emitting its
feature columns + target column:

## Magnetic field regime

- `transforms/bs_aromatic_he.py` — Pople ring current target. Features:
  ring distance, cos(normal angle), ring atom count. Target: BS T0.
- `transforms/bs_aromatic_he_with_t2.py` — same features, full T2
  tensor target. Tests covariant fitting.
- `transforms/hm_aromatic_he.py` — Haigh-Mallion variant. Same features
  as BS; tests whether PySR distinguishes HM's semi-empirical refinement
  from underlying Pople.
- `transforms/mc_backbone_hn.py` — McConnell anisotropy. Target: MC T0.

## Electric field regime

- `transforms/apbs_efg_charged_sidechains.py` — Buckingham field, FF14SB
  charges. Features: E_parallel magnitude, E² magnitude, distance to
  nearest charged sidechain. Target: APBS EFG T0.
- `transforms/apbs_efg_charged_sidechains_mopac.py` — same target, MOPAC
  charges as the field source.
- `transforms/apbs_efg_charged_sidechains_aimnet2.py` — same target,
  AIMNet2 charges.
- `transforms/water_field.py` — solvent electric field. Features: water
  density, per-water charge sum, distance/angle to nearest water cluster.

## Polarisability sensitivity (Buckingham field + CRG)

- `transforms/apbs_efg_with_crg.py` — same as APBS EFG but includes
  AIMNet2 charge response gradient as a feature.
- `transforms/water_field_with_crg.py` — water field + CRG.

## Orbital / hydrogen-bond

- `transforms/larsen_hbond_donors.py` — Larsen 1pHB/1pHaB/2pHB/2pHaB
  family. Features: H-bond distance, angle, donor/acceptor classes.
  Target: Larsen H-bond shielding.

## Validation rediscoveries

- `transforms/mopac_coulomb.py` — by construction Coulomb on MOPAC
  charges. Tests implementation.
- `transforms/mopac_mc.py` — McConnell on MOPAC-derived susceptibilities.

## Added 2026-05-30 (canonical-rediscovery scope expansion — "science lotto tickets" framing)

- `transforms/vibrational_correction.py` — target form `σ_observed =
  σ_static + α · ⟨(ΔR)²⟩` for bond-stretch-modulated shielding
  [Sundholm-Gauss-Schäfer 1996; Jameson 1996 NMR Perspective].
  Features from `traj.bond_strain(atom, bond_idx, frame_t)` and
  `traj.bond_strain_along_bond(atom, frame_t) -> Vec3` (substrate
  Phase 3 additions). Per-stratum vibrational coefficient α
  extracted via PySR. Methods-paper claim: Stage 2 trajectory data
  yields direct rediscovery of the vibrational-averaging correction
  term whose ab initio computation has been documented for decades
  but whose closed-form kernel form has not been quantitatively
  extracted from trajectory data.

- `transforms/wishart_sequence_context.py` — target form `σ_HN(t) =
  f(ring_current, McConnell) + g(ψ_{i-1})` for backbone HN;
  analogous `σ_HA(t) = f(...) + g(φ_{i+1})` for HA [Wishart-Sykes-
  Richards 1991 chemical shift index]. Features from
  `atom.neighbor_residue_phi/psi/chi(frame_t, offset)` for offset ∈
  {−2, −1, +1, +2} (substrate Phase 3 additions). PySR fits the
  additive sequence-context term g separately per atom class.
  Methods-paper claim: Stage 2 trajectory data yields quantitative
  rediscovery of the Wishart chemical shift index's sequence-
  context dependence, extending the empirical correlation to
  per-kernel resolution.

## Stage 2 ridge baseline

- `transforms/ridge_baseline_features.py` — 55 Stage 1 kernels + Welford
  spreads + KernelDynamics scalars + Reorient/iRED dynamics + OpenFold
  embeddings. Target: full T2 tensor. Consumed by `sdk/ridge_baseline.py`
  to reproduce + extend the Stage 1 ridge fit.

## Stage 1 rebuild (toolchain validation)

- `transforms/stage1_rebuild_features.py` — emits exactly the Stage 1
  ridge feature set against the static-snapshot fleet. Consumed by
  `sdk/ridge_stage1.py`; should reproduce R²=0.818.

---

# Methodology rigor (per-kernel standard practice)

Settled per the 2026-05-30 methodology-additions agent. These are
SDK-side practices that apply uniformly across all rediscovery runs,
not platform-substrate work. They preempt the standard reviewer
objections to symbolic-regression methods papers.

## Repeated-seed stability runs per kernel

For each kernel × stratum × variant, run PySR N≥10 times with
different RNG seeds, identical data and hyperparameters. Report:
- The discovered equation(s), sorted by frequency across seeds
- The coefficient distribution for the top-frequency form (mean ± std
  across the seeds that found it)

A single-run result is methodologically weak ("did you get lucky?").
Reporting "form discovered in 9/10 seeds; coefficient 1.20 ± 0.04
(mean ± std across runs)" is the rigor bar. Single cheapest move
toward defensible rediscoveries.

PySR runs are stochastic; 10 seeds in serial is ~10× a single run.
Overnight slot per kernel. Outputs are tiny JSONs.

## Bootstrap coefficient confidence intervals

After PySR converges on a form, refit JUST that functional form (no
further search — fixed structure) to ~200 bootstrap resamples of the
(frame, atom) rows. Report bootstrap mean + 95% CI per coefficient.

Different uncertainty than seed-stability:
- Seed-stability captures search noise
- Bootstrap captures data noise

Together: "Form discovered in 9/10 seeds (search stable);
coefficient 1.20 [1.16, 1.24] (data CI)." Methods readers want both.

Implementation: `scipy.optimize.curve_fit` on the form's sympy
expression, 200 bootstrap resamples, ~30 lines per SDK script.

## Seeded vs cold runs side-by-side

For each kernel, run PySR twice:
- **Cold**: open-ended search with no seed
- **Seeded**: literature equation pre-loaded as a `populations[0]`
  seed individual

Same operator set, same hyperparameters, same data. Both Pareto
frontiers reported side by side. Three useful outcomes:
- Cold finds the literature form: strongest rediscovery claim
- Seeded survives Pareto but cold doesn't: data supports the form,
  search wasn't sufficient — points at data sufficiency separately
  from search adequacy
- Seeded gets killed by Pareto: the literature form doesn't fit the
  data; report what does

Currently a cold-search failure is ambiguous (form absent vs search
inadequate). Seeded runs resolve the ambiguity.

Implementation: `populations` parameter on PySR, ~20 lines per SDK
script for seed-injection + side-by-side report. Needs reliable
sympy→PySR tree conversion; one-time PySR API exploration.

## Residual-structure diagnostic per stratum

After a kernel's discovered equation is settled, compute the residual
`r = DFT_T0 − discovered_equation(features)` for every (frame, atom)
and emit a residual diagnostic:
- Residual vs each feature (scatter + LOWESS)
- Residual autocorrelation in trajectory time
- Residual binned by adjacent-residue identity
- Residual vs each *other* kernel's output

Auto-generated PNG grid per kernel rediscovery.

**Where the discovered equation fails IS the new physics.** A
ring-current rediscovery whose residual correlates with a different
kernel suggests the SR didn't extract that mechanism's contribution
from the candidate features. This is the "interesting" success
category from §Success-criteria, made systematic instead of
opportunistic.

Implementation: one SDK script with parameterised grid plot, applied
to every kernel rediscovery output.

## Reference equation library (YAML)

Curated per-publication YAML with citation, sympy form, variable
dictionary. Concretizes the "literature matcher" from the
campaign-operations section into a usable artifact BEFORE the campaign
starts:

```yaml
- id: pople_1956_ring_current
  citation: "Pople JA. J Chem Phys 24:1111 (1956)"
  sympy_form: "K * (3*cos_theta**2 - 1) / r**3"
  variables:
    r: "distance from atom to ring centroid (Å)"
    cos_theta: "cosine of angle between ring normal and atom-centroid vector"
  K_published: -27.5e-6  # ppm·Å³, ring-class-specific
  applicable_strata: ["aromatic_HE", "aromatic_HD", "aromatic_HZ"]

- id: mcconnell_1957_anisotropy
  citation: "McConnell HM. J Chem Phys 27:226 (1957)"
  sympy_form: "(delta_chi / (3 * r**3)) * (1 - 3*cos_theta**2)"
  # ... etc
```

The matcher reads this and reports, for each PySR output, the best-
matching reference with similarity score (sympy structural distance +
coefficient ratio). Output:
`"Best lit match: Pople 1956 Eq 3, structural match exact, coefficient ratio 1.06."`

Forces explicit hypothesis statement BEFORE fitting, which preempts
"we discovered X, X is sort of like Pople" comparison weakness. ~1
day of curation with the actual papers open + ~100 lines of sympy
matcher code.

## Pre-screening: per-stratum sample-size + signal-to-noise gate

Pre-flight diagnostic before SR runs. For each (kernel, atom-class)
stratum, report:
- N rows
- var(target)
- var(target | best linear predictor)
- R² of a 3-feature OLS
- Condition number of feature matrix

Strata that fail thresholds (e.g. N < 500, OLS R² < 0.05, condition
> 1e6) get flagged or skipped with a logged reason. Without this
gate, PySR will hallucinate a "discovered equation" on noise.

Reporting "we ran PySR on 23 strata; 18 had sufficient signal; here
are the 18 results + the 5 we excluded with reasons" is
methodologically honest. ~100 lines, runs against any transform
Parquet pre-PySR.

## Cross-form convergence as a confidence signal

The campaign already runs multiple equation forms (instantaneous,
lagged, dynamics-averaged, differential, per-trajectory descriptors,
spectral, per-atom curves) per the equation-forms section. SDK
aggregator compares discovered equations across forms:

- Instantaneous form `f(r, θ)` and dynamics-averaged form
  `f(⟨r⟩, ⟨θ⟩)` discovering the SAME functional form with
  coefficients within bootstrap CI → cross-form convergence
- Diverging forms → BO-breakdown signal or slow-mode coupling

Results-table layout: rows = kernels, columns = equation forms,
cells = discovered equation + coefficient. Convergence/divergence
pattern is the methods-paper money shot.

## Cross-tool agreement on flagship kernels

For the three highest-priority rediscoveries (ring current, McConnell
anisotropy, Buckingham field) — not all 8-10 — run AI Feynman as a
parallel SR tool. Also take PySR's discovered form and refit
coefficients via plain ridge/OLS on the same data; report coefficient
agreement.

If three independent approaches agree on form and coefficient, the
rediscovery is rock-solid. AI Feynman is install-painful; limit to
three flagship kernels for high-confidence cross-validation without
burning compute on all 10.

PySINDy stays in scope for the differential equation forms
specifically (chain-rule cross-check on the instantaneous form per
the BO physics framing in the equation-forms section).

## Unit consistency check (lightweight)

Each transform's schema column declares its unit (e.g. `r_ring:
Angstrom`, `cos_theta: dimensionless`, `target_BS_T0: ppm`). Post-walk,
a `pint`-based checker validates discovered equations: do the
dimensional combinations work out? Flag any equation whose unit
signature doesn't reduce to ppm.

Lighter than AI Feynman's full dimensional analysis; catches "PySR
put r⁻² instead of r⁻³ and the units are off by Å" silent bugs.

## Standardized SR-result JSON schema

Every PySR SDK script writes JSON following one fixed schema:

```json
{
  "kernel": "BS",
  "stratum": "aromatic_HE",
  "variant": "FF14SB",
  "n_rows": 12450,
  "top_form": {
    "sympy_str": "-27.5e-6 * (3*cos_theta**2 - 1) / r**3",
    "ascii_str": "...",
    "complexity": 8
  },
  "seed_stability": {
    "n_seeds": 10,
    "freq": 0.9,
    "mean": 1.20,
    "std": 0.04
  },
  "coefficient_bootstrap_ci": [1.16, 1.24],
  "lit_match": {
    "best": "pople_1956_ring_current",
    "score": 0.98,
    "citation": "Pople JA. J Chem Phys 24:1111 (1956)"
  },
  "oos_r2": 0.873,
  "residual_summary": {...},
  "pareto_frontier": [...]
}
```

Then a single aggregator script reads all of them across the campaign
and produces master comparison tables, methods-paper figure data, and
cross-walk DuckDB queries. New kernels added later automatically
appear in the aggregator's tables.

## Per-protein control on 1UBQ (3 flagship kernels only)

Run the full campaign on 1P9J as planned, plus a parallel cut-down
run on one second protein (1UBQ, already a Larsen test fixture). For
the 3 flagship kernels only. If equations and coefficients agree, we
have minimal generalisation evidence at near-zero scope expansion.

Heads off the obvious reviewer objection "you fit on one protein,
would it generalise?" without opening Stage 3. Reuses every SDK
script. Real cost is walk time, not engineering.

Requires 1UBQ trajectory at comparable resolution (15ns, 20ps
stride) with all kernels. If unavailable, defer to Stage 3.

## Reproducibility manifest

Already covered by the substrate's provenance sidecar
(`substrate_conventions_2026-05-30.md`). The walk's manifest.json
(JSON via `QJsonDocument`, not YAML, per the 2026-05-30 Qt-primitives
alignment pass) captures every convention. SDK scripts extend the
same manifest
with PySR version, Julia version, SymbolicRegression.jl version,
RNG seeds for each run, all hyperparameters, transform-file hashes,
Parquet hash, OS, CPU.

Pairs with the standardised JSON schema above so every reported
finding can be replicated bit-identically.

## DFT component-decomposed targets per kernel (added 2026-05-30)

Standard practice for every kernel rediscovery transform: emit
target columns for σ_total, σ_dia, σ_para separately (and the T2
analogs). SDK runs PySR three times per kernel and reports
component-decomposed Pareto frontiers.

**Why this serves canonical rediscovery**: per Ramsey 1950 theory,
σ_dia (sum over occupied orbital current contributions) and σ_para
(orbital-mixing with low-lying empty orbitals) have different
geometric dependencies that the geometric kernels resolve differently:

- Ring current contributes almost entirely to σ_para (orbital-mixing
  effect) — fitting BS/HM kernels against σ_para alone isolates the
  canonical mechanism cleanly; σ_total mixes signal with σ_dia
  background
- Buckingham polarisability splits: α·E term mostly contributes to
  σ_para; β·E² term contributes to σ_dia
- Bond anisotropy (McConnell) is mixed
- Dispersion (intermediate-range orbital effects) contributes
  predominantly to σ_dia

Component-decomposed rediscovery sharpens each kernel's literature-
form fit by removing the orthogonal component as noise. Cross-
component agreement on the discovered form (Pople's `(3cos²θ−1)/r³`
should appear with comparable coefficient in σ_total and σ_para
fits, near-zero in σ_dia fit) is itself a Ramsey-theoretic
confirmation that strengthens the rediscovery.

**Substrate-side**: `traj.dft_shielding_components(atom, frame_t) ->
{total, dia, para}` is already in the substrate's Phase 4 surface
(per `reader_as_platform_2026-05-29.md`). This methodology section
just commits the campaign to using the separation explicitly rather
than fitting only against σ_total.

**Implementation**: each transform emits target columns
`target_<kernel>_T0_total`, `target_<kernel>_T0_dia`,
`target_<kernel>_T0_para` (and T2 analogs for kernels that fit
tensor targets). SDK script wraps PySR in a triple-call loop:

```python
for component in ['total', 'dia', 'para']:
    target_col = f'target_{kernel}_T0_{component}'
    result = run_pysr(X, df[target_col], ...)
    results_per_component[component] = result
```

Reported in the standardised JSON schema with a `component_breakdown`
section: discovered forms + coefficients per component + agreement
score.

## Things deliberately excluded from methodology

- **Full AI-Feynman dimensional-analysis tagging.** Overkill; the
  unit consistency check captures the value cheaply.
- **LOPO across the fleet as Stage 2 default.** Wrong scope; Stage 2
  is single-protein deep study. 1UBQ control above is the
  right-scoped version.
- **Ensemble voting across SR tools.** When tools disagree, the right
  move is to investigate, not average. Cross-tool agreement check
  reports disagreements, doesn't hide them.
- **Trajectory-cumulative-mean as a fourth equation form.** Already
  covered by the existing instantaneous + dynamics-averaged + lagged
  forms; adding more dilutes the cross-form convergence test.
- **Pre-screening with model-based signal probes.** The
  sample-size + OLS-R² gate catches the worst failures cheaply;
  richer pre-screening is premature.
- **SR-of-the-residual chaining.** The residual-structure diagnostic
  does it informally; chaining searches compounds noise into
  "discoveries."
- **Custom equivariant SR for T2 tensors.** Literature equations are
  T0-scalar; comparison currency is T0. T2-covariant SR is Stage 3+
  work needing its own literature anchor. The substrate's T2 tensor
  surface (per the platform doc's dipolar tensor operators) is
  available for transforms that emit T2 components as separate
  scalar targets if useful.

---

# Operational shape

No phase ladder. The work is one continuous engineering effort:

- Build `h5reader-core` (CMake refactor — extract a library from the
  current monolithic target). Mechanical, ~½ day.
- Add `src/analysis/` substrate primitives (`QtSpatialIndex`,
  `QtWindowDerivatives`, `QtTopologyQueries`). C++ tests on a 1P9J
  fixture. Parallelisable across Codex sessions.
- Add `QtTransformHarness` + `QtTransformLoader` + `QtParquetSink` +
  `main_batch.cpp`. Tests against a stub Python transform that just
  counts calls.
- Add pybind11 bindings exposing the substrate to Python.
- Write the transforms. Each is ~30-100 lines. Parallelisable across
  Codex sessions (one transform per session).
- Write the SDK scripts. Each is ~30-50 lines.
- Run the walks against 1P9J. ~10 seconds load (positions + 13 kernel
  TRs, no per-frame NPYs unless a transform asks for them) + iteration
  time per protein.
- Run the SDK scripts against the walk outputs. As often as we want,
  against the same Parquets.
- Analyse, write up, iterate.

The Stage 1 rebuild lands the same way: a `stage1_rebuild_features`
transform + a `ridge_stage1` SDK script. If it reproduces R²=0.818, the
toolchain is validated. If it doesn't, debug. No special infrastructure
distinguishes Stage 1 from Stage 2; it's the same platform with a
different transform + SDK pair.

---

# Storage + aggregation

Filesystem layout per walk:

```
runs/
└── 1p9j_2026-06-XX/
    ├── manifest.json                 # protein, walk-config, transform list (JSON via QJsonDocument)
    ├── bs_aromatic_he.parquet
    ├── mc_backbone_hn.parquet
    ├── apbs_efg_charged_sidechains.parquet
    ├── apbs_efg_charged_sidechains_mopac.parquet
    ├── apbs_efg_charged_sidechains_aimnet2.parquet
    ├── larsen_hbond_donors.parquet
    ├── ridge_baseline_features.parquet
    ├── ...
    └── walk.log                      # progress, warnings, timing
```

SDK script outputs land alongside:

```
runs/1p9j_2026-06-XX/
├── ...
└── sdk_outputs/
    ├── pysr_bs_aromatic_he.json      # discovered Pareto frontier
    ├── pysr_mc_backbone_hn.json
    ├── ridge_baseline.json
    ├── mace_train.checkpoint
    └── ...
```

Cross-walk aggregation happens via DuckDB over the per-walk Parquets +
SDK output JSONs. No Postgres. No server. No schema migrations.

Per the prior architecture decision: Stage 2 does NOT use Postgres. The
manifest is a YAML file; status is per-job marker files; queries are
DuckDB SQL over the result tree.

---

# Cross-protein scaling

Same transforms, more proteins:

```bash
h5reader --batch \
  --proteins '/path/to/fleet/*/trajectory.h5' \
  --transforms transforms/bs_aromatic_he.py \
               transforms/mc_backbone_hn.py \
               ... \
  --output-dir runs/fleet_2026-XX-XX/
```

The harness iterates proteins, loads each in turn, walks, writes per-
protein output directories, unloads, advances. Per-box memory is bounded
by the largest protein in the slice. The 600-protein run is days of wall
clock, parallelisable across scan1/scan2/scan3/batcave by sharding the
`--proteins` glob.

Per-protein walk cost is dominated by the per-frame transform calls.
For ~40K atoms × 20K frames the walk is hours per protein at typical
transform rates. **This is why "one walk, many transforms" is the
killer facility** — you can't afford to walk three times for three
transforms; you walk once and emit three outputs.

---

# Success criteria

Per kernel rediscovery:
- **Strong**: discovered equation matches published form exactly,
  coefficient within 30% of published magnitude
- **Acceptable**: same functional form, coefficient differs (literature
  gives baseline; data refines)
- **Interesting**: discovered equation differs in chemistry-
  interpretable ways (new term appears, classical approximation breaks
  down in this regime)
- **Failure**: no simple equation at complexity ≤ 15 with test R² > 0.8.
  Itself a finding: kernel depends on collective effects or longer-
  range context.

Per charge-source comparison:
- All three sources converge to same equation: equation robust; pick
  cheapest source
- Sources diverge: methodologically interesting; thesis ranks them by
  R² + literature agreement

Per polarisability sensitivity:
- CRG term appears with significant coefficient: polarisability matters
- CRG term coefficient is ≈ 0: static charges suffice

---

# What is NOT in this doc

- **Movie capture per protein.** Future capability. Architecture
  accommodates by linking VTK in batch mode (already required for
  shared codebase with the GUI binary), but the rendering loop +
  ffmpeg encoder are not in Stage 2 scope.
- **Live in-app PySR overlays.** Deferred Pattern B above.
- **Live in-walk experiment dock.** Deferred Pattern A above.
- **Expert-system rule library / narrative engine.** Stage 3+ items.
- **No commitment to chapter structure / writing schedule.**

---

# Cross-references

- Supersedes the per-kernel discovery sketch in
  `per_kernel_pysr_discovery_2026-05-29.md` (kept as early-iteration
  history).
- Companion: `stage2_initial_model_2026-05-29.md` — the parallel ridge
  baseline workstream, now a `ridge_baseline_features` transform +
  `ridge_baseline.py` SDK script through the same platform.
- Companion: `molprobity_conformational_tool_2026-05-29.md` — per-residue
  conformational classification, available as a transform feature.
- Companion: `frame_partition_index_2026-05-29.md` — Postgres-backed
  spatial-query substrate; explicitly NOT a Stage 2 dependency.
- Deferred Stage 3+: `expert_system_chemistry_overlay_2026-05-29.md`,
  `per_atom_narrative_engine_2026-05-29.md`.
- `project_calibration_done` memory: Stage 1 settled at R²=0.818.
- `project_three_stages` memory: Stage 1 / Stage 2 / Stage 3 framing.
- `project_first_ml_combination` memory: MACE + OF3 + Larsen as the
  next-ML-stage after rediscovery.
- `project_compute_fleet_machines` memory: scan1/scan2/scan3 + batcave
  for parallel walks.
- `reference_qwen_systemd`, `project_qwen_summary_workflow` — Qwen for
  per-kernel-finding prose writeup once equations are discovered.
- `feedback_extractor_untouchable` — Stage 2 doesn't modify the producer.
  Analogous reader-discipline: Stage 2 additions live in
  `src/analysis/`; `src/app/` widgets/scenes/overlays/REST are NOT
  modified unless explicitly working on an interactive-payoff feature.
- dft-ex1 precedent at `/shared/dft-ex1/` — pybind11 embed pattern
  applied to a research C++ binary. Same pattern here, opposite control
  direction (pybind11 exposes C++ to Python rather than C++ calling
  Python callbacks, though both directions may end up coexisting in the
  unified binary).
