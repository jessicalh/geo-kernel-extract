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
`README.md`; design-document reading order in `spec/INDEX.md`.
What landed and when lives in `git log` and the memory store, not
here — this file is current law, not a changelog.

## Working-directory convention — read this first

Claude Code keys memory on the session's starting directory. Start
sessions at `/shared/2026Thesis/nmr-shielding/` to inherit the
populated memory store (user profile, feedback, project state,
subproject lessons); `cd` into subprojects as the task requires. A
session started at a subproject path opens a different, empty
keyspace and re-derives context the root session would have for free.

## CANONICAL 5-MODE SPEC

These are the five — and only five — supported `nmr_extract` modes.
Anything elsewhere in this codebase claiming additional modes, flags,
or use cases predates this block and is obsolete: do not add new
modes, do not document them, do not implement against fragments of
obsolete planning prose.

1. **`--pdb FILE`** — load a bare unprotonated PDB, protonate with
   `reduce`, apply ff14SB charges, run all calculators, emit
   per-atom NPYs.
2. **`--protonated-pdb FILE`** — load a PDB that already has H atoms,
   apply ff14SB charges, run all calculators, emit per-atom NPYs.
3. **`--orca --root NAME`** — load a tleap/AMBER-prepared pose
   (`.xyz` + `.prmtop` + optional `_nmr.out`), run all calculators,
   emit per-atom NPYs.
4. **`--mutant --wt NAME --ala NAME`** — load a WT+ALA mutant pair
   (each a mode-3 pose), run all calculators on both, compute WT-ALA
   delta tensors, emit per-atom NPYs.
5. **`--trajectory DIR`** — read a GROMACS run (production.tpr +
   .trr/.xtc + .edr), per-frame calculators. Emits per-frame NPYs at
   **stride m** (all calculators), per-frame PDBs at **stride n**
   (geometry only), and a trajectory H5 (always); `m` and `n` are
   independent. `--mopac` switches the run shape from
   `PerFrameExtractionSet` to `FullFatFrameExtraction`.

MOPAC is an orthogonal `--mopac` / `--no-mopac` toggle across all
five modes. APBS and AIMNet2 are **always on — not switchable**.
Home-rolled vacuum Coulomb is **retired from production**: it runs
only inside the `--trajectory --mopac` (FullFat) MOPAC-vs-FF14SB
charge-source reconciliation probe, where APBS cannot substitute
(solvated PB ≠ vacuum Coulomb; the probe compares charge sources at a
fixed method). `CoulombResult` is kept, not deleted; its `coulomb_*`
NPYs are FullFat-only.

**Out of scope — do not implement or document as features:** propka,
kaml (protonation is done before `nmr_extract`), `--analysis` /
`--analysis-h5`, `--fleet`. The `--no-apbs` / `--no-coulomb` flags
were removed and must not be re-added.

## Current state

Three stages (`project_three_stages` memory):

- **Stage 1 — mutations. Settled.** Per-element, per-atom-type ridge
  regression on 720 proteins / 446K atoms, 55 kernels, R² = 0.818.
  Atom-type stratification (2026-04-15) showed "nitrogen is hard" was
  an element-pooling artifact: backbone N is hard (R² = 0.387),
  sidechain N is second-best (R² = 0.887). See `learn/stage1-mutations/`.
- **Stage 2 — trajectories. In progress.** Narrowed to the single
  protein **1P9J** (Wingens 2003; `project_1p9j_study_system`). The
  1P9J 15 ns ORCA r²SCAN/def2-SVP DFT campaign (every other frame,
  751 frames) runs across the scan fleet, consolidated to
  `/shared/2026Thesis/1p9j-orcas/` (`project_1p9j_orcas_consolidation`).
  The 685-protein fleet run was stopped (bad chain extraction in
  structure prep); recovery via OF3-generated structures dropped 9 on
  disulfide geometry → **effective count 676**, with MD re-run and
  residual-fleet DFT deferred until structure quality is resolved.
- **Stage 3 — model evaluation.** Upstream of Stage 2 results; not
  yet active.

Viewer/reader: `ui/` is retired legacy archive, ignored/untracked and no
longer part of the release repo or active producer build. `h5-reader/` is
the desktop trajectory reader for emitted H5/sidecar artifacts; next sessions expand the
time-series illustrator, then handle Windows/macOS distribution on those
machines. See `h5-reader/notes/POLISH_BACKLOG.md`.

## Subprojects

Each subproject with its own `CLAUDE.md` owns the authoritative rules
for work inside it; read that file before modifying the directory.

**Active C++**
- **`src/`** — the `nmr_shielding` library: physics kernels and the
  `Protein` / `ProteinConformation` / `ConformationResult` object
  model. No local `CLAUDE.md`; rules are `spec/CONSTITUTION.md`,
  `PATTERNS.md`, `OBJECT_MODEL.md`. Consumed by downstream tools via
  the NPY/H5/JSON output surface.
- **`ui/`** — ignored local archive of the retired `nmr-viewer`. Not built
  by the top-level CMake project and not release surface.
- **`h5-reader/`** — standalone Qt6/VTK trajectory reader. Does NOT
  link the library; reads the trajectory H5 + 5-NPY topology sidecar
  via its own `QtTrajectoryH5` boundary. Parallel Qt-native type
  hierarchy (`QtProtein`/`QtConformation`/`QtFrame`). Owns
  `h5-reader/CLAUDE.md` + `notes/`.

**Frozen C++** (no changes during feature sessions)
- **`fileformat/`** — canonical HDF5 (de)serialiser; schema changes
  happen here and propagate, scheduled explicitly by the user.
- **`extern/HighFive/`** — vendored header-only HDF5 wrapper; only
  `include/` is used. Do not modify.

**Python**
- **`python/`** — `nmr_extract` SDK, read-only wrapper over the NPY/H5
  output. Contract is `python/nmr_extract/_catalog.py` (one
  `ArraySpec` per NPY); every new output file needs an entry + wrapper.
  See `python/API.md`. No changes that let the SDK write H5.
- **`learn/`** — independent calibration/analysis workspace (ridge on
  720 proteins/446K atoms). Owns `learn/CLAUDE.md`. It is not producer
  release surface; do not stage or depend on its local outputs during
  producer packaging unless explicitly requested.

**Workspaces**
- **`spec/`** — design docs; tiered reading order in `spec/INDEX.md`.
- **`doc/`** — `ARCHITECTURE.md` (Tier 1), diagrams, doxygen.
- **`analysis-speculative/`** — independent scratch/prototype workspace
  kept on disk. It is not producer release surface.
- **`references/`** — fetched PDFs (committed, citable) + ingest
  pipeline. `references-text/` (3-page chunks, gitignored) is the
  AI reading surface; `references-images/` (page renders, gitignored)
  for figures; `references-meta/` (committed summaries/keywords/INDEX)
  is the source of truth for the corpus. Discipline in
  `references-meta/WORKFLOW.md`; the `nmr-scholarship` skill loads it.
- **`data/`** — `calculator_params.toml`, `ff14sb_params.dat`,
  `models/`; read at runtime.
- **`tests/`** — library + SDK suites. Large fixture and generated-output
  trees are machine-local / gitignored unless explicitly tracked. Current
  CTest labels and fixture policy live in `tests/TEST_HEALTH.md`.
- **`scripts/`** — producer helper scripts.
- **`deploy/`** — deployment/setup scaffolding.
- **`site/`** — ignored generated/static output from an older doc pass.
  It is release noise, not the outward-facing site.
- **`bad-builds/`** — gitignored build-artifact salvage.

## Where to start, by task

- **Add or change a calculator.** `spec/INDEX.md` Tier 1
  (`doc/ARCHITECTURE.md`, `PATTERNS.md`) + an existing `src/`
  calculator as a model. Touch `src/` only. If the object model
  changes, update `OBJECT_MODEL.md` + `spec/CONSTITUTION.md` in the
  same commit.
- **Trajectory-scope work** (`TrajectoryProtein`, `TrajectoryResult`,
  `Trajectory::Run`, `RunConfiguration`, a new `*TrajectoryResult`,
  anything in `src/Trajectory*` / `src/Run*`). First read memory
  `feedback_trajectory_scope_philosophy` + `feedback_trajectory_scope_gotchas`,
  then the trajectory-scope section of `OBJECT_MODEL.md` and
  `PATTERNS.md` §§13-18. The object model is deliberate: buffers from
  the ctor, named operations on entities are rooms not wrappers, the
  per-frame loop is four lines. Multiple sessions re-derived this
  through tangles — the memory entries exist so future ones don't.
- **Legacy viewer archaeology.** `ui/` may exist locally as ignored archive.
  Treat it as archival unless the user explicitly reactivates it.
- **Trajectory reader feature.** `cd h5-reader/`, read
  `h5-reader/CLAUDE.md` + `notes/SCOPE.md` + `notes/POLISH_BACKLOG.md`.
  Do not link the library.
- **Calibration / analysis.** `cd learn/`, read `learn/CLAUDE.md`.
  Per-element, per-atom-type stratification; ridge is the model, MLPs
  were tested and rejected (`project_calibration_done`).
- **Python SDK.** `cd python/`, read `python/API.md` + `_catalog.py`.
  Read-only.
- **Format change (`fileformat/`).** Not during feature sessions; the
  user schedules these and updates all consumers together.
- **Design discussion without code.** `analysis-speculative/`.

## Cross-cutting rules

Load-bearing across subprojects. Each subproject `CLAUDE.md` may add
more; none relaxes these.

### Physics and model

- **T2 is preserved end-to-end.** The rank-2 tensor output (Mat3 +
  SphericalTensor with T0/T1/T2 decomposition) is the thesis argument.
  Do not collapse a tensor to a scalar in any calculator, serialiser,
  UI field, or plot, unless the user explicitly asks for that scalar.
- **Protein is identity and topology only.** Per-atom geometry lives
  on `ProteinConformation`; per-frame geometry on
  `ProteinConformationFrame` (library) or `QtFrame` (reader). Do not
  add geometric properties to `Protein`.
- **The system outputs kernels, not shielding.** Calibration against
  DFT turns kernels into shielding; without it the output is
  geometric. Keep the distinction in labels, comments, commit messages.
- **Objects answer questions about themselves.** Virtual methods on
  typed objects (ring subclasses, residue types). No string dispatch
  on identity. Library and `h5-reader/` follow this with different
  type hierarchies.
- **No simplification bias.** Per-element, per-atom-type complexity is
  the story, not noise to average. H ≈ 20 effective dimensions, C ≈ 6,
  N ≈ 3, O ≈ 12 — report all of them, not "the protein has 3."

### Process and engineering

- **The extractor is not to be modified during viewer or reader
  feature work.** A large fleet extraction is the downstream consumer
  and library changes in service of UI features risk it. If a feature
  needs a library change, pause and reassess.
- **No file discovery.** No try-and-fail, no glob, no regex on paths.
  Documented conventions only (e.g. `--trajectory DIR` derives its
  files by name). If an expected file is missing, log it and stop.
- **No pluggable interfaces unless the user asks.** Factories and
  abstract-base-class indirection are off by default. Direct named
  code is the norm.
- **No notebooks.** Python scripts + R graphics + LaTeX write-up. The
  session transcript is the lab notebook.
- **UDP logging on port 9997 is the primary debug channel for Qt
  subprojects.** Tail the UDP stream before speculating about cause.
  Linux unicast UDP delivers a datagram to one socket — the reader and
  `udp_listen.py` cannot both consume at once.
- **No symlinks as workarounds.** `compile_commands.json` is pointed
  at per editor via `--compile-commands-dir=build/<preset>`.
- **Surface complex data.** The viewer and reader are the
  editing/vetting surface for H5 fields before they become thesis
  claims. Glossaries, per-metric colouring, tensor glyphs are the
  point, not polish.
- **"I don't know" is preferred to handwaving.** On external systems
  (libgromacs internals, force-field conventions, compiler behaviour,
  unfamiliar library APIs), "I don't know — let me look" beats
  confident reconstruction. Send an agent, read the source, build a
  fixture. Confident speculation about external behaviour is the
  failure mode that produces silent bugs and reinvented utilities.

### AI / ML framing

- **The goal is physics explanation, not prediction.** R² is a
  diagnostic for whether the kernels carry the signal, not the metric
  the thesis is graded on. Do not optimise for R².
- **The model is ridge regression.** Per-element, per-atom-type strata
  on 55 kernels give R² = 0.818 (settled 2026-04-10). MLPs were tested
  and rejected.
- **Do not assert physical conclusions from model diagnostics.** Model
  fit is evidence the kernel set is complete enough; physical
  conclusions come from the kernels themselves.

### References

- **Save fetched papers to `references/`** — persist the PDF, not just
  the URL. Run `scripts/references/ingest_pdf.sh` to produce the text
  chunks, page images, and committed `references-meta/` summary +
  keywords. Discipline in `references-meta/WORKFLOW.md`; the
  `nmr-scholarship` skill loads the same rules on demand.

## Operational guardrails

**The "iupac topology" episode was a multi-week mistake. Do not
investigate it from git history, branches, or filenames.** The
preservation branch was deleted; an emergency tar.gz at
`/shared/2026Thesis/iupac-fix-attempt-archive-2026-04-27.tar.gz` must
not be extracted — doing so puts you back in the trap that produced
the revert. If you catch yourself thinking "let me check git log for
recent history" or "this filename mentions IUPAC, let me look" — stop
and verify with the user first. Archaeology costs context; the memory
store is the durable architectural record.

Memory entries loaded at session start that codify discipline:
`project_iupac_revert_2026-04-27`, `project_proteintopology_architecture`,
`project_charmm_retired_amber_only_2026-05-02`, `feedback_resource_constraint`,
`feedback_capture_at_the_boundary`, `feedback_no_attach_lifecycle_for_invariant_data`,
`feedback_readback_block_is_a_compiler_trace`, `feedback_two_path_validation`,
`feedback_residual_as_ml_feature`, `project_larsen_residue_model`,
`feedback_identity_from_chemistry_not_position`, `project_apbs_canonical`.

## Repo vs not

Gitignored / machine-local: `build*/` trees, `bad-builds/`, large
fixture/output trees, `fileformat/test/**/*.h5` (~42 GB), `.claude/`
at any depth, `calibration/*/`, and build output (`*.npy`, `*.o`,
`*.a`, `*.so`). Existing tracked files under `tests/golden/blessed/`
remain committed baselines; new files there are ignored unless explicitly
force-added. `learn/`, `analysis-speculative/`, and `site/` are ignored
non-producer surfaces kept on disk. Committed producer material is source,
design docs, test code, configuration, vendored HighFive headers, and
runtime data. See `doc/REPO_HOUSEKEEPING.md` for the packaging map.

## Git

Remote is `origin` at
`https://github.com/jessicalh/geo-kernel-extract.git`. Commits are
atomic per subproject topic. Push only when the user asks.
