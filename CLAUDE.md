# nmr-shielding — project root

Thesis research tree for NMR chemical shielding prediction via
geometric kernels on protein structures. The physics pipeline is
classical calculators (ring current, electric field gradient, bond
anisotropy, etc.) calibrated against DFT WT-ALA deltas. The
engineering layers are a C++ library for extraction, two Qt6/VTK
visualisation subprojects, a Python SDK for reading output, and a
Python/R analysis pipeline for calibration and model work.

This `CLAUDE.md` is the entry point for a Claude session opened
against the tree. Human-oriented project description lives in
`README.md`. Design-document reading order lives in
`spec/INDEX.md`.

## Working-directory convention — read this before anything else

Claude Code stores memory keyed on the session's starting working
directory. A session started at
`/shared/2026Thesis/nmr-shielding/` reads the memory store
populated by prior sessions (several dozen entries covering user
profile, feedback, project state, and subproject-specific lessons
from the library, viewer, and reader work). A session started at
a subproject path (e.g. `/shared/2026Thesis/nmr-shielding/h5-reader/`)
opens a different keyspace and begins without that history.

Start sessions at the project root and `cd` into subprojects as the
task requires. This trades some directory-listing noise for memory
continuity. If a subproject session becomes necessary (uncommon
case), expect to re-derive context that the root session would
have had for free.

## Subprojects

Each subproject that has its own `CLAUDE.md` owns the authoritative
rules for work inside it. The pointers below name those files; read
them before modifying anything in the named directory.

### Active C++

- **`src/`** — the `nmr_shielding` library. Physics kernels and the
  `Protein` / `ProteinConformation` / `ConformationResult` object
  model. Built by the top-level `CMakeLists.txt`. No local
  `CLAUDE.md`; rules come from `spec/CONSTITUTION.md`,
  `PATTERNS.md`, and `OBJECT_MODEL.md`. The library is consumed by
  `ui/` via direct linking and by the Python SDK via a separate
  JSON/NPY/H5 output surface.
- **`ui/`** — `nmr-viewer`, the single-conformation Qt6/VTK viewer.
  Links the library directly so the renderer consumes the library's
  typed objects without duplication. Owns `ui/CLAUDE.md` and
  `ui/UI_ROADMAP.md`. Reads the analysis H5 via `fileformat/` as a
  read-only companion (`--analysis-h5 PATH`, Session 2026-04-16 in
  `UI_ROADMAP.md`). Never writes H5; never triggers a new extraction.
- **`h5-reader/`** — standalone Qt6/VTK trajectory reader. Does NOT
  link the library; consumes `fileformat/analysis_file.cpp` via
  source include only. Owns `h5-reader/CLAUDE.md` plus
  `notes/SCOPE.md`, `notes/POLISH_BACKLOG.md`,
  `notes/TIME_SERIES_EXPANSION.md`, `notes/RESIDUAL_RENDER_DROP.md`.
  Target audience: advisers on Linux / macOS / Windows. Parallel
  type hierarchy (`QtProtein`, `QtConformation`, `QtFrame`) mirrors
  the library's but is Qt-native; the two are not shared types.

### Frozen C++

- **`fileformat/`** — canonical HDF5 (de)serialiser for the analysis
  trajectory output of `nmr_extract --trajectory --analysis`.
  Consumed by `ui/` and `h5-reader/` via source include. Schema
  changes happen here and propagate outward; no subproject modifies
  `fileformat/` during feature sessions. Round-trip tests
  (`roundtrip_write.cpp`, `roundtrip_test.cpp`) and
  `npy_cross_validate.py` live alongside.
- **`extern/HighFive/`** — vendored HighFive (header-only C++ wrapper
  around HDF5). Only `include/` is used; upstream `src/`, `tests/`,
  `deps/` are gitignored. Do not modify.

### Python

- **`python/`** — `nmr_extract` SDK. Read-only Python wrapper around
  the library's NPY and H5 output. 62 tests cover a catalog of
  60-plus arrays. The format contract is
  `python/nmr_extract/_catalog.py` — one `ArraySpec` per NPY file.
  Any new library output file needs a `_catalog.py` entry and a
  wrapper class here. See `python/API.md`.
- **`learn/`** — calibration and analysis scripts. Ridge regression
  on 720 proteins and 446K atoms, per-element and per-atom-type
  stratification. Owns `learn/CLAUDE.md`. Python scripts plus R for
  graphics; no notebooks (see cross-cutting rules). The current
  Stage 1 rule, set 2026-04-15, is: stratify by AMBER atom name
  within each element (CA, C=O, CB, C sidechain, N bb, N side, etc.).

### Workspaces and smaller directories

- **`spec/`** — design documents. Tier-based reading order in
  `spec/INDEX.md`. CONSTITUTION, MATHS_GOALS,
  GEOMETRIC_KERNEL_CATALOGUE, TRAJECTORY_EXTRACTION,
  ANALYSIS_TRAJECTORY, TEST_FRAMEWORK, and decision records live
  here. `spec/meta-docs-review/` holds the 2026-04-03 documentation
  audit.
- **`doc/`** — `ARCHITECTURE.md`, mermaid diagrams, doxygen output.
  ARCHITECTURE is Tier 1 reading per `spec/INDEX.md`.
- **`analysis-speculative/`** — investigation scratch area. Python
  prototype trajectory readers (v1 and v2 with FEEDBACK.md on each)
  and design essays (INVESTIGATION_LAYERS, SCOPE, SIGNALS,
  THESIS_ARC, COMPUTE_BUDGET, EFG_GEOMETRY_DOMINANCE). Not in the
  build. A space to work through scope decisions before they
  graduate into `h5-reader/`, `learn/`, or `spec/`.
- **`references/`** — fetched papers (PDF). Save every fetched
  reference here; see `ANNOTATED_BIBLIOGRAPHY.md`.
- **`data/`** — `calculator_params.toml` and `ff14sb_params.dat`
  plus the `models/` directory. Library reads from here at runtime.
- **`tests/`** — library and SDK test suites. `tests/golden/` holds
  blessed smoke-test baselines (machine-local, gitignored).
  `tests/data/fleet_test_large/` is a 295 MB MD trajectory fixture
  set; gitignored, regenerable.
- **`protonation/`** — design documents for the protonation
  pipeline (BUILDER_AGENT_PROMPT, BUILDER_ANALYSIS, DESIGN_HISTORY).
- **`fleet-manager-tasks/`** — coordination notes for fleet
  extraction runs (RERUN_FULL_SYSTEM, STATE).
- **`baseline_features/`**, **`calibration/`**, **`train/`** —
  per-protein output workspaces. Mostly machine-local; see
  `.gitignore`. Some data here is produced by batch runs of
  `nmr_extract` and consumed by `learn/`.
- **`scripts/`** — helper shell scripts (`batch_extract.sh` etc.).
- **`deploy/`** — deployment scaffolding (`setup_scan.sh`).
- **`site/`** — static documentation site (`build.sh`, HTML pages).
  Not referenced by any subproject's code.
- **`bad-builds/`** — salvaged artifact area from earlier builds.
  Gitignored entirely. Keeps ELF binaries and CMake cache around
  for reference without committing them.

## What's in the repo vs not

Gitignored, machine-local:

- `build/`, `build-ui/`, `build-check/`, `build-forensic/` — CMake
  build trees.
- `bad-builds/` — artifact salvage area (above).
- `tests/data/fleet_test_large/` — 295 MB MD trajectory fixtures.
- `fileformat/test/**/*.h5` — analysis H5 fixtures, individual
  files 1–3 GB, total roughly 42 GB.
- `.claude/` at any depth — Claude session shell history and
  worktree metadata.
- `learn/runs/`, `learn/bones/old_extractions/`, `calibration/*/`
  — regenerable intermediate data.
- `*.npy`, `*.mc`, `*.o`, `*.a`, `*.so` — build and solver output.

Committed: source, design docs, test code, configuration,
vendored HighFive headers, reference PDFs.

## Where to start, by task

Pick the row that matches the task and read in order.

- **Add or change a calculator.** `spec/INDEX.md` Tier 1
  (`doc/ARCHITECTURE.md`, `PATTERNS.md`) + an existing calculator
  in `src/` as a model. Touch `src/` only. If the object model
  changes, update `OBJECT_MODEL.md` and `spec/CONSTITUTION.md` in
  the same commit.
- **Trajectory-scope work** (TrajectoryProtein, TrajectoryAtom,
  TrajectoryResult, Trajectory::Run, RunConfiguration, RunContext,
  adding a `*TrajectoryResult` subclass, anything in `src/Trajectory*`
  or `src/Run*`). **MANDATORY** before touching code: read the memory
  entry `feedback_trajectory_scope_philosophy.md`, then
  `spec/TRAJECTORY_LANDING_STATE_2026-04-23.md`, then the relevant
  WIP spec sections. The trajectory object model is deliberate
  organising principle — buffers from ctor, named operations on
  entities are rooms not wrappers, the per-frame loop is four lines.
  Multiple sessions have re-derived this through painful tangles;
  the memory entry exists so you do not.
- **Single-conformation viewer feature.** `cd ui/` and read
  `ui/CLAUDE.md`. Do not modify the library in service of a viewer
  feature (see "extractor is not to be modified" below).
- **Trajectory reader feature.** `cd h5-reader/` and read
  `h5-reader/CLAUDE.md`, `notes/SCOPE.md`, `notes/POLISH_BACKLOG.md`.
  Do not link the library; per-atom fields come from the H5,
  volumetric fields from `QtBiotSavartCalc` / `QtHaighMallionCalc`.
- **Calibration or analysis work.** `cd learn/` and read
  `learn/CLAUDE.md`. Per-element, per-atom-type stratification is
  the current discipline. Ridge is the model; MLPs were tested and
  rejected (see `project_calibration_done` memory entry).
- **Python SDK work (`nmr_extract`).** `cd python/` and read
  `python/API.md` plus `_catalog.py`. Read-only reader; no changes
  that would let the SDK write H5.
- **Format change (`fileformat/`).** Not during feature sessions.
  The user schedules format changes explicitly and consumers
  (viewer, reader, SDK) update together.
- **Design discussion without immediate code.**
  `analysis-speculative/` is the place to write a scope note, try
  a Python prototype, or chew through a design question before
  committing to a real subproject.

## Cross-cutting rules

These are load-bearing across subprojects. Each subproject's
`CLAUDE.md` may add more; none relaxes these.

### Physics and model

- **T2 is preserved end-to-end.** The rank-2 tensor output (Mat3 +
  SphericalTensor with T0 / T1 / T2 decomposition) is the thesis
  argument. Do not collapse a tensor to a scalar in any
  calculator, serialiser, UI field, or plot, unless the user
  explicitly asks for that specific scalar.
- **Protein is identity and topology only.** Per-atom geometry
  lives on `ProteinConformation`; per-frame geometry on
  `ProteinConformationFrame` (library) or `QtFrame` (reader). Do
  not add geometric properties to `Protein`.
- **The system outputs kernels, not shielding.** Calibration against
  DFT turns kernels into shielding; without calibration, the output
  is geometric. Keep the distinction in labels, comments, and
  commit messages.
- **Objects answer questions about themselves.** Virtual methods on
  typed objects (ring subclasses, residue types). No string
  dispatch on identity. The library and `h5-reader/` follow the
  same rule with different type hierarchies.
- **No simplification bias.** Per-element and per-atom-type
  complexity is the story, not noise to be averaged. H has roughly
  20 effective dimensions, C has 6, N has 3, O has 12 — report all
  of them, not "the protein has 3 dimensions."

### Process and engineering

- **The extractor is not to be modified during viewer or reader
  feature work.** A 2000-protein extraction run is imminent and
  library changes in service of UI features risk that run. If a
  feature needs a library change, pause and reassess.
- **No file discovery.** No try-and-fail, no glob, no regex on
  paths. Documented conventions only (e.g. `--trajectory DIR`
  derives `md.tpr`, `md.xtc`, `md.edr` by name). If the expected
  file is missing, log it and stop.
- **No pluggable interfaces unless the user asks.** Factories and
  abstract-base-class indirection are off by default. Direct named
  code is the norm.
- **No notebooks.** Python scripts plus R for graphics plus LaTeX
  for write-up. The session transcript is the lab notebook.
- **UDP logging on port 9997 is the primary debug channel for Qt
  subprojects.** When something misbehaves in `ui/` or
  `h5-reader/`, tail the UDP stream before speculating about cause.
  Linux unicast UDP delivers a datagram to exactly one socket;
  the reader and `udp_listen.py` cannot both consume at the same
  time.
- **No symlinks as workarounds.** `compile_commands.json` is
  pointed at per editor via `--compile-commands-dir=build/<preset>`,
  not symlinked to the project root.
- **Surface complex data.** The viewer and reader are the
  editing / vetting surface for H5 fields before they become thesis
  claims. Per-field glossary, per-metric colouring, tensor glyphs
  on atoms — these are not polish, they are the point.

### AI / ML framing

- **The goal is physics explanation, not prediction.** R² is a
  diagnostic for whether the kernels carry the signal. It is not
  the metric the thesis is graded on. Do not optimise for R².
- **The model is ridge regression.** Per-element, per-atom-type
  strata on 55 kernels give R² = 0.818 (settled 2026-04-10). MLPs
  were tested and rejected; they add no signal.
- **Do not assert physical conclusions from model diagnostics.**
  Model fit is evidence that the kernel set is complete enough.
  Physical conclusions come from the kernels themselves.

### References

- **Save fetched papers to `references/`.** Always persist the PDF
  alongside a note in `ANNOTATED_BIBLIOGRAPHY.md`, not just the URL.

## Current state (as of 2026-04-17)

Thesis is structured in three stages (`project_three_stages`
memory):

- **Stage 1 — mutations.** Settled. Per-element and per-atom-type
  ridge regression on 720 proteins / 446K atoms. 55 kernels,
  R² = 0.818. Atom-type stratification (2026-04-15) showed that
  "nitrogen is hard" was an element-pooling artifact: backbone N
  is hard (R² = 0.387), sidechain N is second-best atom type
  (R² = 0.887). See `learn/stage1-mutations/`.
- **Stage 2 — trajectories.** In progress. MD protocol is a single
  non-PLUMED 25 ns run per protein; 600 of ~1200 frames adopted as
  the baseline. Analysis-trajectory mode emits the exhaustive
  per-frame H5 consumed by the viewer and the reader. DFT reference
  structures are 26 per protein at 1 ns intervals (ns0 … ns25).
  Calibration DFT batch complete 2026-04-18: 260/260 across 10 proteins,
  0 failures. Results collected into `/shared/2026Thesis/fleet_calibration_dft/`
  (authoritative staging with per-job PDB + XYZ + .out + meta.json, worker logs,
  collection manifest) and replicated as orthogonal `dft_output/` siblings into
  all three calibration copies (`fleet_calibration-{working,stats,backup}/`).
  Fleet state: 200 of 685 proteins H5-complete; remaining 485 MDs
  running, PDBs landing ~2026-04-24, all in the same tree. DFTs on
  the 485 have not been scheduled.
- **Stage 3 — model evaluation.** Upstream of Stage 2 results. Not
  yet active.

Viewer and reader status:

- `ui/` — single-conformation viewer with analysis H5 companion
  (Session 2026-04-16). Read-only H5 consumption, Time Series (H5)
  dock tab, atom identity check with specific-mismatch logging.
- `h5-reader/` — feature-complete at the scaffold level as of
  2026-04-17. Atom picker, selection overlay, inspector dock,
  atom time-series dock via Qt6 Charts, BS / HM butterfly
  isosurfaces, B-field streamlines, backbone ribbon, ring polygons.
  Next 2–3 sessions expand the time-series illustrator (per-atom
  glyphs, colour bubbles, field glossary, tensor display). See
  `h5-reader/notes/POLISH_BACKLOG.md`. Cross-platform build pass
  (macOS, Windows, DGX Spark ARM64) deferred until after those
  feature sessions; it runs from separate per-platform Claude
  sessions against the shared directory.

## Git

Remote is `origin` at
`https://github.com/jessicalh/geo-kernel-extract.git`. Commits are
atomic per subproject topic (one commit per directory or feature).
Push only when the user asks.
