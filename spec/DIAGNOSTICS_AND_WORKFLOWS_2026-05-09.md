# Diagnostics and Workflows — 2026-05-09 consolidation

**Status:** PENDING-WORK consolidation doc. Captures diagnostic and
workflow items migrated from `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md`
§D and §E as part of the 2026-05-09 post-topology doc-cleanup pass.
These are tools that operate over existing calculator output — not
calculators themselves. DESIDERATA retires after this consolidation;
this doc becomes the home for diagnostic / workflow design intent.

The five items here all share one architectural shape: they consume
calculator output that already exists (per-atom T2 emissions,
`atoms_category_info.npy`, `OrcaShieldingResult` deltas) plus
substrate identity columns, and surface cross-protein, cross-calculator,
cross-frame, or cross-bench patterns that single-protein /
single-calculator output cannot. Items 1 and 2 are first-class
*diagnostic* outputs over the 260 DFT pose set that auto-validate
GEOMETRIC_KERNEL_CATALOGUE claims and expose T2 residual structure.
Item 3 is a fleet-scale event extractor that turns trajectory data
into named events (ring flips, rotamer transitions, H-bond
break/form). Item 4 is the architectural surface that unifies probe-
point kernel evaluation across calculators, viewers, and training-
data generators. Item 5 is a CLI workflow mode for per-bench pipeline
runs.

The diagnostics described here run AGAINST existing fleet output.
They are consumers of the 2026-05-08 substrate slice
(`atoms_category_info.npy` + per-calculator T2 emissions), not new
calculators that produce per-atom fields.

---

## 1. Per-calculator T2-residual map (cross-protein diagnostic)

*(from DESIDERATA §D.2)*

### What
Per-atom T2 residual map over the 260 DFT pose set, binned by atom
type, secondary structure, distance-to-ring, and the substrate
typed-identity columns from `atoms_category_info.npy`. Residual is
DFT delta T2 minus the classical kernel prediction T2; one map per
calculator.

First-class pre-computed output rather than a notebook artifact —
the diagnostic table is generated alongside the calculator outputs
themselves so the binning and stratification convention is shared
across thesis figures.

### Why
MATHS_GOALS pillar 2 made physical: where does each calculator's T2
residual concentrate? If the residual is structureless across a
binning axis, the calculator already explains that variation; if the
residual concentrates in one bin (e.g., backbone N in α-helix), that
bin is a thesis claim or a kernel gap. Thesis-grade diagnostic per
13.4 posture.

The substrate slice (`atoms_category_info.npy`) lands the typed
identity columns this diagnostic stratifies by — atom mechanical
identity, ring position, planar group kind, secondary structure,
protonation variant. Without the substrate, every calculator-side
notebook had to re-derive its own binning convention and the
diagnostics did not compose.

### Implementation sketch
Python-side post-pass over fleet output. Inputs per protein: 260
DFT pose set's `OrcaShieldingResult` T2 + each calculator's per-atom
T2 + `atoms_category_info.npy`. Output: one HDF5 group
`/diagnostics/t2_residual_map/<calculator_name>/` carrying:

- per-atom residual (N, 5) — Mat3 → SphericalTensor T2 component
- binning axes from `atoms_category_info.npy` columns
- per-bin mean / variance / count

Per-bench HDF5 pairs with item 5 — a per-bench mode could emit
only the residual-map slices that bench needs.

### Dependencies
- `OrcaShieldingResult` per-atom on the 260 DFT pose set (exists)
- Per-calculator T2 emission (exists for 8/10 classical calculators;
  MopacCoulomb / MopacMcConnell shielding contributions still TBD —
  DESIDERATA §B.6, captured in PLANNED_CALCULATORS amendment)
- `atoms_category_info.npy` (landed 2026-05-08, commit `8accdb6`)
- Fleet calibration workspace (`fleet_calibration-stats/`)

### Cost
Negligible compared to the calculator runs themselves. One pass
over the 260-pose fleet H5 + binning + summary statistics.
Recomputable on demand from existing fleet output.

### Origin
- DESIDERATA §D.2 (2026-04-22)
- MATHS_GOALS pillar 2 (T2 angular residual)
- 260-pose DFT calibration set (collected 2026-04-18)

### Status
**PENDING.** No code yet. Substrate slice is the unblocker — before
2026-05-08 the typed identity columns were not available outside
the C++ library. After landing, the binning is one boolean-mask
operation per axis. Implementation lives in `learn/` against the
fleet calibration workspace, not in the C++ library.

---

## 2. Per-kernel-pair T2 correlation map (kernel-independence diagnostic)

*(from DESIDERATA §D.3)*

### What
Per-atom `|cos|` between the T2 component of each calculator pair
over the 260 DFT pose set. One scalar per (atom, calculator-pair);
aggregable into per-pair distributions, per-element distributions,
or per-substrate-axis stratifications.

### Why
Auto-validates GEOMETRIC_KERNEL_CATALOGUE independence claims
(`PQ vs McConnell ~0.38`, `MopacCoulomb vs ff14SB Coulomb ~0.85`,
`PQ vs Coulomb ~0.40`, `PQ vs RingChi ~0.37`, `PQ vs BiotSavart
~0.57` — all well below 0.9). These claims are currently asserted
from physical reasoning and a small initial validation; an automated
per-pair map turns them into a living, fleet-recomputed thesis
figure that updates whenever calculators or fleet poses change.

If a previously-independent pair correlates highly on a substrate-
stratified slice (e.g., PQ vs McConnell on backbone-only atoms),
the correlation reveals a hidden coupling and a possible kernel
redundancy. If MopacCoulomb vs ff14SB Coulomb falls below 0.85 on
a substrate slice, the gap is the calibration-hope test of
DESIDERATA §D.5 — does the polarisation signal carry information
beyond fixed charges?

### Implementation sketch
Same Python post-pass infrastructure as item 1. Per-atom T2 per
calculator → pairwise `|cos|` matrix per atom → aggregate into
the binning axes from `atoms_category_info.npy`. Output groups
`/diagnostics/t2_kernel_pair_correlations/<calc_a>/<calc_b>/` with
per-atom + per-bin summary stats.

For 8 classical calculators that is C(8, 2) = 28 calculator pairs.
With MopacCoulomb + MopacMcConnell on the sparse pose set, 45
pairs total. Manageable surface.

### Dependencies
- Per-calculator per-atom T2 (same as item 1)
- `atoms_category_info.npy`
- 260-pose DFT calibration set

### Cost
Negligible. O(N_atoms × N_calc_pairs) per pose; the 260 DFT pose
set across all proteins is small relative to the per-frame
trajectory output.

### Origin
- DESIDERATA §D.3 (2026-04-22)
- GEOMETRIC_KERNEL_CATALOGUE independence claims (lines 492–493)

### Status
**PENDING.** No code yet. Same unblocker as item 1.

### Cross-link to D.5 reconciliation diagnostic
DESIDERATA §D.5 (MOPAC-vs-ff14SB kernel-signal reconciliation on
the 260 DFT poses) is a specialised case of this map applied to
the (MopacCoulomb − Coulomb, DFT delta − Coulomb) and analogous
McConnell pair. Implementation is the same Python infrastructure
applied to differences rather than direct values. **Captured here
as a sub-bullet rather than its own item** because the
infrastructure and substrate dependencies are identical; the user
direction is to keep diagnostics composable, not to fragment.

---

## 3. Event menu hookable extractors (workflow / fleet)

*(from DESIDERATA §E.2)*

### What
Threshold-parameterised menu of named-event extractors over
trajectory data:

- **Ring flip detection** — per-residue ring-normal autocorrelation
  + dihedral hop detection (M17 Akke-Weininger); flip events
  identified by χ₂ crossing on Phe / Tyr / His.
- **Rotamer transitions** — Ramachandran φ/ψ + sidechain χ₁ bin-
  crossings per residue, time-stamped.
- **H-bond break/form events** — per-bond `HBondResult` presence
  threshold-crossings.

Output: per-trajectory HDF5 group `/derived_events/<event_kind>/`
with start-frame, end-frame, residue indices, and per-event
metadata. Hookable in both C++ (as a TrajectoryResult subclass
that scans existing per-frame output) and Python (so a threshold
sweep against a frozen fleet H5 is one invocation, not a
regeneration pass).

### Why
Shielding fields don't sit still — backbone σ on a Phe ring
neighbour swings tens of ppm during a flip. A per-protein "list of
events" + their frame ranges turns "did this protein flip its
Phe28 ring during the trajectory?" from a manual viewer-eyeballing
question into a computed substrate column.

Pairs naturally with the substrate identity stratification of items
1–2 (which strata flip more, which fewer? does flip rate correlate
with predicted T2 residual?) and with item 4's probe-point
evaluation API (events as triggers for finer-grained kernel
evaluation around the flipping atom).

### Implementation sketch
Two hooks:

1. **C++ TrajectoryResult subclasses** — `RingFlipEventTrajectoryResult`,
   `RotamerTransitionEventTrajectoryResult`, `HBondEventTrajectoryResult`.
   AV (always-valid mid-stream) shape; consume per-frame output of
   already-attached calculators (`RingNormalResult`,
   `RamachandranResult`, `HBondResult`); emit event records to
   `TrajectoryResult` storage.
2. **Python threshold-sweep wrapper** — `nmr_extract` SDK consumer
   that reads frozen fleet H5, applies thresholds parametrically,
   emits events HDF5 group. Same event schema as C++ output;
   different threshold values without re-running extraction.

Threshold parameters TOML-configured per `data/calculator_params.toml`
convention; HDF5 schema documented in
`spec/I_O_AND_SCHEMA_2026-05-09.md`.

### Architectural-state note (2026-05-09)

`ChiRotamerSelectionTrajectoryResult` already lives in the tree as
the existing rotamer-selection TrajectoryResult. Before implementing
a separate rotamer-transition extractor here, check that TR for
overlap — duplicating may be unwarranted.

### Dependencies
- Per-frame `RingNormalResult`, sidechain χ angles (substrate),
  `HBondResult` (all exist)
- `TrajectoryResult` infrastructure (landed)
- TOML parameter conventions (landed)

### Cost
C++ extraction: negligible vs the per-frame work that already runs.
Python sweep: per-protein H5 read + threshold; seconds per protein.

### Origin
- DESIDERATA §E.2 (2026-04-22)
- IDENTITY_AND_DYNAMICS_ROLLUP §13.7 (event-menu posture)
- M17 Akke-Weininger 2023 (ring flips)

### Status
**PENDING.** Schema decisions defer until per-bench schema work
(item 5) settles — the event records share the per-bench HDF5 slice
discipline.

---

## 4. Probe-point evaluation API (architectural surface)

*(from DESIDERATA §E.3)*

### What
Unified C++ surface for evaluating the existing kernel set
(McConnell, RingSusceptibility, BS, HM, Coulomb, HBond) at arbitrary
non-atom probe positions. One non-atom-aware code path that backs:

- A.4 `NICSProbeEvaluator` — NICS(0) at ring centroids, NICS(1)zz
  at ±1 Å offsets (PLANNED_CALCULATORS Amendment 2026-05-09(d)).
- C.7 volumetric training-data generation — Cartesian grid sampling
  for FNO / equivariant-NN training inputs.
- Viewer butterfly rendering — `h5-reader/`'s `QtBiotSavartCalc` /
  `QtHaighMallionCalc` consume a unified library evaluator instead
  of maintaining a parallel volumetric path.
- User-defined lattice sampling — externally specified probe sets
  for custom investigations.

### Why
Today the library evaluates kernels at atoms only; the viewer and
the prospective NICS / volumetric / training-data work all need
the same kernel math at non-atom points. Each consumer that re-
implements is a chance for the volumetric path to drift from the
calibrated atom path, and the project's discipline of "one source
of truth for kernel evaluation" is not held.

A single probe-point API: input is a position (or batch of
positions) plus a kernel selection; output is the same Mat3 +
SphericalTensor used at atoms. Calibration parameters apply
identically — calibrated kernel math is calibrated kernel math
regardless of where the probe sits.

### Implementation sketch
Add `EvaluateAtPoint(Vec3 probe, KernelFilterSet selection,
ProbePointResult& out)` virtual method on each kernel calculator
that already supports atom evaluation. Default implementation:
factor existing per-atom evaluation through a shared helper that
takes an arbitrary `Vec3` source position. Calculators where
"source position" is well-defined (BS, HM, McConnell, Coulomb,
RingSusceptibility, HBond) implement directly; calculators where
the source is anchored to atoms (e.g., DispersionResult) document
what they do at non-atom probes (return zero or skip).

`ProbePointResult` carries the same Mat3 + SphericalTensor minimum
output as a per-atom result, plus the probe position itself for
caller bookkeeping.

C.7 volumetric NPY emission and viewer integration are downstream
consumers. The library does NOT take on grid-generation or render
responsibility — only point-evaluation.

### Dependencies
- Existing kernel calculators (refactor, not greenfield)
- `ConformationResult` / atom evaluation paths must not regress;
  the refactor extracts a shared helper, doesn't change the atom
  output

### Cost
Per-evaluation: same cost as per-atom evaluation. No new physics.
Refactor cost: per calculator, ~50–100 lines extracting the
position-agnostic core. Pre-substrate this would have been more
invasive; the substrate slice already factored a lot of "is this
atom in a ring?" / "is this atom a backbone N?" identity logic out
of the kernel inner loops, so the remaining math is more amenable
to point-extraction.

### Origin
- DESIDERATA §E.3 (2026-04-22)
- NICS literature (Schleyer 1996 B14, Gershoni-Poranne 2015 B15)
- E12 VN-EGNN virtual-node philosophy

### Status
**PENDING.** A.4 `NICSProbeEvaluator` is the most concrete
near-term consumer (PLANNED_CALCULATORS 2026-05-09(d) amendment
flagged as the trigger). Recommend co-designing the API with the
NICS calculator implementation rather than speccing the API in the
abstract.

---

## 5. Per-bench pipeline mode (workflow / CLI)

*(from DESIDERATA §E.6)*

### What
CLI entry point: `nmr_extract --bench loth-2005` (or
`--bench yao-bax-2010`, etc.) runs only what that named bench
requires:

- Subset of calculators (e.g., Loth-2005 ubiquitin CCR needs σ
  tensor + NH bond vectors + Larmor frequencies → A.10
  `CCRRateResult` + dependencies; not the full 25-calculator
  pipeline).
- Required `ExperimentalReferenceLoader` for that bench (DESIDERATA
  §C.2).
- A.11 `BenchmarkBackCalculationResult` for that bench
  (PLANNED_CALCULATORS Amendment 2026-05-09(i)).

Initial bench roster from PHYSICS_FOUNDATIONS §0.9:
Loth-2005-ubiquitin-CCR, Yao-Bax-2010-15N-CSA-split,
Babaei-2017-bulk-χ, Tharayil-2021-PCS, Wylie-2006-δ22-CO-HN.

### Why
Validation work today must run the full pipeline and slice out the
relevant fields post-hoc. Each bench has narrow input requirements
that the full pipeline doesn't share — running everything for one
bench is wasteful (DFT-equivalent runs per bench at fleet scale)
and obscures the question. A per-bench mode makes each validation
into a discrete, reproducible, bench-named CLI invocation; the
methods text and the figure caption point at the same command.

Pairs with item 1's per-bench T2 residual slice and with the C.7
bench HDF5 slices (`/benches/loth_2005/`, etc.) captured separately
in `spec/I_O_AND_SCHEMA_2026-05-09.md`.

### Implementation sketch
- New bench-registration mechanism in `nmr_extract` CLI:
  `BenchSpec` carrying name + calculator subset + experimental
  loader + back-calculation result type.
- `--bench <name>` flag dispatches through a `RunBench(BenchSpec)`
  path that schedules only the named calculators and emits the
  bench-scoped HDF5 group.
- Initial registry hard-codes the five PHYSICS_FOUNDATIONS §0.9
  benches; future benches register additively.

Small CLI surface; the heavy lifting is the existing calculator
machinery + the loader + the back-calculation result. The bench
entry is an orchestration layer, not new physics.

### Dependencies
- C.2 `ExperimentalReferenceLoader` (PENDING; DESIDERATA §C)
- A.11 `BenchmarkBackCalculationResult` (PLANNED_CALCULATORS
  Amendment 2026-05-09(i))
- Per-bench input data (BMRB / RefDB tables, χ-tensor values,
  experimental CCR rates) — sourcing per-bench
- C.7 bench HDF5 schema (`spec/I_O_AND_SCHEMA_2026-05-09.md`)

### Cost
Implementation: small CLI + dispatch layer. Per-bench data
gathering varies by bench (some are ~kB tables; some require
literature digitisation).

### Origin
- DESIDERATA §E.6 (2026-04-22)
- PHYSICS_FOUNDATIONS §0.9 validation benches table
- Workflow ergonomics for thesis-figure reproducibility

### Status
**PENDING.** Blocked behind C.2 `ExperimentalReferenceLoader` and
A.11 `BenchmarkBackCalculationResult` (both planned-calculator
amendments). The CLI surface itself is straightforward once those
two land.

---

## Cross-references

- **A.4 `NICSProbeEvaluator`** (PLANNED_CALCULATORS Amendment
  2026-05-09(d), PENDING) consumes item 4 probe-point API as its
  evaluation primitive; co-design recommended.
- **A.11 `BenchmarkBackCalculationResult`** (PLANNED_CALCULATORS
  Amendment 2026-05-09(i), PENDING) is the per-bench result type
  invoked by item 5 per-bench pipeline mode; output flows into
  C.7 bench HDF5 slices.
- **C.7 bench HDF5 slices** (planned `spec/I_O_AND_SCHEMA_2026-05-09.md`
  entry) pairs with item 5; one bench-named HDF5 group per
  validation bench.
- **DESIDERATA §D.5 MOPAC-vs-ff14SB reconciliation** is captured
  as a sub-case of item 2 (per-kernel-pair correlation applied to
  differences); same infrastructure, no separate item.
- **PATTERNS.md "Calibration-ready vs raw-geometric separation"**
  (DESIDERATA §E.5) is a discipline rule on calculator output
  shape — handled in PATTERNS.md, not here.
- **Substrate slice** (`atoms_category_info.npy`, commit
  `8accdb6`) is the load-bearing dependency for items 1, 2, 3 —
  the typed identity columns are what stratification axes mean.

## How this doc is used

Forward diagnostics/workflow work touches one section, picks up
its pending design intent, and either implements OR refines the
design here. Updates land as amendments at the bottom of the
relevant section, same convention as
`spec/PLANNED_CALCULATORS_2026-04-22.md`. Re-read the relevant
section against the current pipeline and substrate before
implementation; ideas were captured against the 2026-04-22 state
of the world, refined here against the 2026-05-09 state, and
neither is the pipeline as it will look when the diagnostic is
written.

DESIDERATA §D and the diagnostic / workflow portions of §E retire
to bones once the full doc-cleanup pass completes; this doc is the
canonical home for the design intent. Items in §A (calculators),
§B (calculator variations), and §C (input/output surfaces) of
DESIDERATA migrate elsewhere — see
`spec/plan/post-topology-doc-cleanup-2026-05-09.md` for the
overall migration map.

## Amendments

*(None yet. Append below, never rewrite above.)*
