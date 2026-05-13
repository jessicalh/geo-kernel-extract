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

## Current state — read this BEFORE the rest of CLAUDE.md

### Landed (master at HEAD `d1ad904`)

The 2026-04-29 → 2026-05-08 work landed across `master`:

- **AMBER charge slice + AMBER-as-project-standard** (commits `4ba5491`
  through `91e4124` and follow-ups). Substrate-first sequencing,
  `LegacyAmberTopology` typed contract on `Protein`,
  `ForceFieldChargeTable` as a first-class object, AMBER terminal-
  residue + cap handling. CHARMM context retired.
- **GROMACS readback block** (commit family May 2). Selective-authority
  merge of GROMACS-pdb2gmx FF-port labels to canonical AMBER labels via
  the `topol.top` rtp comment line — compiler-trace shape: read at
  load, applied to typed slots (`Residue.type`,
  `Residue.protonation_variant_index`), discarded.
- **Bundle B + C — substrate generator + ring substrate** (commits
  through `6ec9bff`). Typed `AtomSemanticTable` records emitted by
  cifpp + RDKit build-time generator (NOT linked into the runtime
  library); per-residue ring inventory now reads from substrate
  (`SemanticAt(ai).ring_position`) instead of walking string tables.
- **NamingApplicator + load-time canonicalisation** (commit `85de965`,
  `4ac7d79`, through `67be414`). Typed rule-set application, per-atom
  transient maps, post-protonation re-canonicalisation. Closes the
  PdbFileReader bypass so every loader path canonicalises.
- **FramePdbEmitter** (commits `6639f88`, `85c1637`). Opt-in per-frame
  PDB writer for trajectory runs; production-validated on 1P9J.
- **Topology sidecar** (commits `f2781da`, `dc50917`, 2026-05-13).
  Per-protein invariant topology projections alongside every nmr_extract
  run: extends `atoms_category_info.npy` with 6 new fields (`chain_id`,
  `residue_number`, `insertion_code`, `parent_atom_index`,
  `ff_atom_type_string`, `equivalence_class`); new `TopologySidecar`
  class emits `residues.npy` + `bonds.npy` + `rings.npy` +
  `ring_membership.npy` + `extraction_manifest.json` (codex
  TOPOLOGY_SIDECAR_CONTRACT). `ArraySpec` gains `native_axis` / `irreps`
  / `units` / `sign_convention` / `tensor_rank` / `parity` / `mechanism`
  fields populated for all ~108 CATALOG entries (resolves OI-016).
  Python SDK enforces invariants in `load()`: required-file check +
  manifest-vs-actual axis sizes + bond endpoint range + ring membership
  refs + residues.atom_count.sum() == n_atoms. Malformed exports fail
  loud. `learn/extract.py` STAGE1_AUDIT_OUTPUTS extended; `learn/src/
  secondary/loader.py` no longer silently skips on FileNotFoundError /
  ValueError. R-side regex-mechanism refactor (OI-120) is the
  downstream consumer.

- **CategoryInfoProjection slice** (commits `8accdb6`, `e1b5bcc`,
  `d1ad904`, 2026-05-08). One structured NPY per protein
  (`atoms_category_info.npy`, ~31 fields) carrying the categorical
  identity record; AMBER → IUPAC / BMRB name projection at the NPY
  emission boundary; `MutationDeltaResult` typed-identity matchup with
  spatial-NN sanity check; six new dia/para shielding NPYs (WT side,
  mut side, deltas); Python SDK `CategoryInfo` wrapper + `atom_nom.tbl`
  pre-flight consistency tests.

- **Tripeptide BB + Neighbor calculators + DFT backend** (2026-05-10/11).
  `TripeptideDftTable` (libpq loader against local `tensorcs15` Postgres
  replica), `TripeptidePoseAssembler` (two-path validation: DFT canonical
  ordering + typed substrate role cross-check; Vec3 residual recording
  as ML feature), `TripeptideBackboneShieldingResult` (σ_BB^i per Larsen
  2015), `TripeptideNeighborShieldingResult` (Δσ_BB^{i±1} per Larsen Eq
  3 cap-side reading — paper-verified against the Cβ/Val worked example).
  Smoke tests green on `1UBQ_pm6dh3plus.pdb`. `RuntimeEnvironment` now
  has section-aware TOML parser exposing `[databases].tensorcs15` DSN;
  `Session::LoadTripeptideDftTable` opens the libpq connection once and
  holds for the session lifetime. CMakeLists adds `PostgreSQL::PostgreSQL`.
  Full record at `spec/plan/session_2026-05-10_part2_tripeptide_calculators.md`.

- **LarsenResidue perception-driven model** (2026-05-11; six rounds
  of adversarial review by codex xhigh + opus general-purpose all
  landed same day). Replaces the retired positional-ordering-table
  design with bond-graph perception at the `TripeptideDftTable` load
  boundary. `src/LarsenResidue.{h,cpp}` builds typed
  `AtomMechanicalIdentity` per DFT atom via **K=3 Weisfeiler-Lehman
  multi-round graph isomorphism** against canonical `AminoAcidType`
  chemistry (HID/HIE/HIP variants strictly enforced via protein-side
  `Residue::protonation_variant_index` hint; ACE/NME hand-coded).
  `LarsenTripeptide` cached on `TripeptideDftRecord`.
  `TripeptidePoseAssembler` central matching is **purely typed-identity**
  end-to-end. Per-perceived-atom dispatch on the round-4
  `canonical_assignment_ambiguous` flag on `LarsenResidue::PerAtom`:
  singleton WL class → STRICT identity match (BranchAddress +
  DiastereotopicIndex bind, per Markley 1998 CIP; protects
  chemistry-distinct branches like ILE CG1/CG2); multi-atom WL class
  → relaxed match dropping BranchAddress + DiastereotopicIndex with
  nearest-spatial tiebreak (resolves graph-automorphic pairs like
  PHE CD1↔CD2 that no K rounds of WL can split). Cap and central
  assembly **decline the residue** when perception fails (no heuristic
  fallback — perception or nothing). The chi1-grid sidechain re-
  rotation `Rsc` was **removed** so position and tensor share the same
  rotation (the BB Kabsch); chi-grid coarseness lives in the Vec3
  `residual_vec` ML feature.

  **Names die at canonical construction (round 4):**
  `MatchPiece` returns canonical NODE INDICES, not atom-name strings.
  `FinalizeAdjacency` translates the name-keyed bond list to index-
  keyed `adj_by_idx`; runtime WL signatures, `MatchPiece`, and
  `EmitPiece` all operate on canonical node indices. `EmitPiece`
  reads typed identity directly from `canon.atoms[canon_idx].identity`
  and dispatches the role-pinned slot cache (N_idx, CA_idx, ...)
  by `BackboneRole`/`Locant` enum values — no name comparisons in
  perception's hot path.

  **Canonical identities grounded in generated topology table
  (round 5):** Chain atoms (NCapAla/CCapAla/Central) route through
  `StampChainIdentitiesViaTable`, which mirrors
  `ComposeAtomSemantic` from `src/LegacyAmberTopology.cpp:233`:
  parse name → build identity → apply methyl-H pseudoatom collapse
  via canonical bond graph (parent has 3+ H neighbours → clear DI)
  → `topology_generated::LookupBy(aa, variant_idx, identity)`. The
  row's identity is the authoritative typed tuple. Lookup miss is
  FATAL with the same shape as the protein side's
  `FatalSubstrateMiss`. Both substrates speak identical typed
  vocabulary post-collapse. Cap atoms (ACE/NME) stay hand-coded
  since they are not in the standard-20 substrate table.

  **DB-side integrity:** `MergeGeometryAndTensors` keys by `atom_idx`
  (not parallel-vector position) with duplicate detection + same-element
  cross-check; silent corruption from a reordered JSONB is structurally
  impossible. DSN is redacted via `PQconninfoParse` (handles case-
  insensitive keys, quoted values with spaces, URI form). The
  `IdentifyCentralBackbone`/`IdentifyAlaCap`/`IdentifyCTermAlaCap` legacy
  public APIs and the deprecated `central_n/ca/c/o` record fields are
  **deleted** — runtime is perception-only.

  **Always-on in nmr_extract** when `[databases].tensorcs15` DSN is
  configured: `Session::LoadTripeptideDftTable` is called in main,
  `RunOptions::tripeptide_dft_table` is threaded through all 5
  single-frame run sites + the trajectory mode, `OperationRunner::Run`
  attaches both `TripeptideBackboneShieldingResult` and
  `TripeptideNeighborShieldingResult` when the table is loaded. The
  BB + Neighbor calculator chi-fallback loops break on
  `IsHit() && larsen.has_value()` so a perception failure at chi depth
  N retries at shallower chi rather than abandoning the residue. The
  AAA reference query in the Neighbor calculator fails loud if its
  larsen is absent — every per-residue assembly subtracts that
  reference, so a silent AAA-reference perception failure would zero
  every Δσ_BB^{i±1} contribution; that failure mode is now caught at
  module entry.

  **Python SDK** (round 4): `python/nmr_extract/_catalog.py` carries
  7 `tripeptide_*` ArraySpec entries; `python/nmr_extract/_protein.py`
  exposes `TripeptideGroup` attached on `Protein.tripeptide` when
  any tripeptide NPY is present; exported from
  `nmr_extract.__init__`.

  **Determinism (round 4 M3):** `TripeptideDftTable::QueryNearest`
  adds `ORDER BY calc_id ASC` so chi-fallback row selection is
  stable across sessions (previously planner-dependent).

  **Per-direction neighbour residual semantics (round 4 M4):**
  `TripeptideNeighborShieldingResult` NaN-initialises the per-
  direction `tripeptide_neighbor_residual_vec_{prev,next}` in the
  per-conf reset loop; the WriteFeatures emitter mirrors the BB
  NaN-fill pattern. Downstream ML reads NaN as "no contribution in
  this direction" — distinguished from a coincidentally-zero
  residual.

  **Round 6 follow-ups (post-R5 investigation):**
  - tensorcs15 DB carries HIP only for HIS (verified by query).
    `TripeptideBackboneShieldingResult::Compute` and
    `TripeptideNeighborShieldingResult::Compute` now emit a
    structured `OperationLog::Warn` when a HIS residue has
    `protonation_variant_index ∈ {0, 1}` (HID/HIE) so the
    silently-skipped residue is no longer silent.
  - `StampChainIdentitiesViaTable` FATAL message now anchors to
    `LarsenResiduePerceptionTest.AllCombinationsPerceiveCleanly`
    as positive coverage so future developers know where to
    verify if they touch one side of the chemistry sync.
  - `test_larsen_residue_wl_ambiguity.cpp` gains
    `ChiFallbackIsDeterministic` — verifies `ORDER BY calc_id
    ASC` produces the same row on repeated queries at a
    fallback depth (ARG at -180/-180, n_chi_axes=0).
  - HBondHα `project_hbond_halpha_design` memory entry now phrases
    the sidechain-O-acceptor exclusion as "Phase 1 scope, revisit
    trigger: post-calibration per-acceptor-class residual".
  - `larsen-residue-design-2026-05-11.md` Open Question §1
    (prochiral methylene H disambiguation) now phrases Phase-2
    as "stays unimplemented unless calibration shows systematic
    pro-R/pro-S bias on methylene Hα/HB/HG/HD/HE pairs" with
    explicit trigger condition.

  **Tests:** 8 dedicated structure tests registered in
  `structure_tests`:
  - `LarsenResiduePerceptionTest.AllCombinationsPerceiveCleanly` —
    parity across all 20 (residue, frame_type) DB combos.
  - `LarsenResidueAgainstSourceLogTest.AaaLogPerceivesCleanly` —
    parses AAA Gaussian log directly, perception runs independent of
    DB (drift detector).
  - `LarsenResidueSerSidechainTest.OgAndBackboneOHaveDistinctIdentities`
    — SER OG ↔ BB O identity-distinctness regression.
  - `LarsenResidueWlAmbiguityTest.IleCg1AndCg2AreChemistryDistinct`
    (round 4) — ILE CG1/CG2 perceive as singleton (chemistry-
    distinct via K=1 WL split); `canonical_assignment_ambiguous=false`.
  - `LarsenResidueWlAmbiguityTest.PheCdAndCeAreGraphAutomorphic`
    (round 4) — PHE CD1/CD2/CE1/CE2 perceive as ambiguous (graph-
    automorphic; no K splits); PHE CZ singleton.
  - `LarsenResidueWlAmbiguityTest.ChiFallbackIsDeterministic`
    (round 6) — `ORDER BY calc_id ASC` makes chi-fallback row
    selection deterministic; two consecutive queries at
    (ARG, -180, -180, n_chi=0) return the same calc_id.
  - `TripeptideBackboneShieldingTest.RunsOn1UbqPm6` and
    `TripeptideNeighborShieldingTest.RunsOn1UbqPm6` — end-to-end smoke
    on `1UBQ_pm6dh3plus.pdb`.

  Python SDK adds `python/tests/test_tripeptide_group.py` (7 tests:
  catalog registration, optional-spec discipline, group attachment,
  dtype, NaN-fill contract for per-direction neighbour residual,
  absent-NPYs → `tripeptide=None`).

  **DB-side faithfulness verified:** DB row 21755 byte-for-byte
  matches the original Gaussian log `AAA_4_54_nmr.log`.

  **Smoke on 1UBQ_pm6dh3plus.pdb (post-round-4, unchanged from rounds 1-3):**
  - BB **74/76 residues / 1205/1232 atoms assigned** (was 1033 with
    all-heuristic) — 97.8% atom coverage.
  - BB **mean RMSD 0.015 Å** (was 0.547, 36× better) / max 0.048 Å.
  - Neighbor 76/76 residues received contributions (was 75/76) /
    522 atom accumulations.
  - 217 non-tripeptide unit tests pass.
  - 7/7 tripeptide structure tests pass.
  - 95/95 Python SDK tests pass.

  Full design at `spec/plan/larsen-residue-design-2026-05-11.md`;
  retires `spec/plan/bones/typed-tripeptide-topology-design-2026-05-10.md`
  (no code landed from the old design). Adversarial review provenance
  in `project_larsen_residue_model` memory entry.

  **Known queued (not blocking):** trajectory mode wires the per-frame
  calculators but has no `TimeSeriesTrajectoryResult` emission surface
  yet — per-frame tensors compute correctly and land on each frame's
  `ConformationAtom` fields, but the H5/NPY surface that surfaces
  them for downstream is queued as
  `TripeptideShieldingTimeSeriesTrajectoryResult` (mirror of
  `BsShieldingTimeSeriesTrajectoryResult`).

### Pending forward work

- **Calculator walkthrough** — bring existing classical calculators
  (`BiotSavartResult`, `McConnellResult`, `CoulombResult`,
  `HaighMallionResult`, etc.) onto the typed substrate. Each calculator
  gets a before / after audit pass with codex 5.5 xhigh as reviewer.
  The substrate's typed identity fields (`SemanticAt(ai)`) replace
  string-discrimination inside calculator code.
- **Planned calculators** — captured in
  `spec/PLANNED_CALCULATORS_2026-04-22.md` (with 2026-05-08 amendment
  for `PlanarGeometryResult`),
  `spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md`, and
  `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md`.
- **`PlanarGeometryResult`** specifically — per-frame conformation
  companion to the substrate's `PlanarGroupKind` / `RingPosition`
  fields. Substrate side LANDED; calculator side PENDING. See the
  Amendment 2026-05-08 in `spec/PLANNED_CALCULATORS_2026-04-22.md`.

### Planning artefacts

- **Active pending design**:
  `spec/plan/openai-5.5-strong-architecture-layout.md` (architecture
  record), `spec/plan/planned-calculator-substrate-audit-2026-05-06.md`
  (substrate ↔ planned-calculator map),
  `spec/plan/comprehensive-calculator-inventory-2026-04-30.md`,
  `spec/plan/md-rerun-685-discussion-priors-2026-04-30.md` (fleet ops).
- **Retired to bones**: detailed planning that captured the AMBER /
  topology / substrate / projection slices' design decisions lives at
  `spec/plan/bones/`. The decisions themselves are in `master`; the
  prose history is preserved for archaeology only. Do not consult
  bones/ docs to drive new work — they describe pre-landing state.

### Operational rules that still apply

The "iupac topology" episode was a multi-week mistake. **Do not
investigate it from git history, branches, or filenames.** The
preservation branch was deleted; an emergency tar.gz sits at
`/shared/2026Thesis/iupac-fix-attempt-archive-2026-04-27.tar.gz` and
should not be extracted. Extracting it would put you back in the trap
that produced the revert in the first place.

If you find yourself thinking "let me check git log to understand the
recent history" or "let me look at this old branch" or "this filename
mentions IUPAC, let me investigate" — **stop and verify with the user
first.** Archaeology costs context. The memory entries plus this
section's landed-work bullet list are the durable architectural record.

Memory entries loaded automatically at session start that codify
operational discipline:

- `project_iupac_revert_2026-04-27`
- `project_proteintopology_architecture`
- `project_charmm_retired_amber_only_2026-05-02`
- `feedback_resource_constraint`
- `feedback_capture_at_the_boundary`
- `feedback_no_attach_lifecycle_for_invariant_data`
- `feedback_readback_block_is_a_compiler_trace`
- `feedback_two_path_validation` (2026-05-10 — cross-substrate matching)
- `feedback_residual_as_ml_feature` (2026-05-10 — Vec3 not scalar)
- `project_tripeptide_calculators_landed` (2026-05-10/11)
- `project_larsen_residue_model` (2026-05-11 — perception-driven; supersedes
  the retired `project_typed_tripeptide_topology_design`)
- `feedback_identity_from_chemistry_not_position` (2026-05-11 — perception
  over ordering tables)
- `project_larsen_neighbor_axa_reuse` (CORRECTED 2026-05-10 — read at
  cap atoms, not central)

<!-- Pre-2026-05-08 working-tree status retired here; the slices above
     have all landed in master. Detailed planning prose is in
     spec/plan/bones/. -->

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
  TrajectoryResult, Trajectory::Run, RunConfiguration, adding a
  `*TrajectoryResult` subclass, anything in `src/Trajectory*` or
  `src/Run*`). Before touching code, read the memory entries
  `feedback_trajectory_scope_philosophy.md` and
  `feedback_trajectory_scope_gotchas.md`, then the trajectory-scope
  section of `OBJECT_MODEL.md` and `PATTERNS.md` §§13-18. For
  working-notes and pending appendices (`NmrAtomIdentity` —
  superseded 2026-04-28 by the `LegacyAmberTopology` plan in
  `spec/plan/openai-5.5-strong-architecture-layout.md` and memory entry
  `project_proteintopology_architecture`; catalog Appendix F; H5
  metadata schema), see
  `spec/plan/bones/pending_include_trajectory_scope_2026-04-22.md` —
  not authoritative for anything landed. Historical landing records:
  `spec/plan/bones/TRAJECTORY_LANDING_STATE_2026-04-23.md` +
  `spec/plan/bones/TRAJECTORY_REFACTOR_GAPS_2026-04-23.md`. The trajectory object
  model is deliberate organising principle — buffers from ctor,
  named operations on entities are rooms not wrappers, the per-frame
  loop is four lines. Multiple sessions re-derived this through
  tangles; the memory entries exist so future sessions do not.
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
- **"I don't know" is preferred to handwaving.** When reasoning
  about external systems (libgromacs internals, force-field
  conventions, compiler behaviour, library APIs we don't routinely
  touch), the honest answer "I don't know — let me look" is
  preferred to confident reconstruction. Send an agent, read the
  source, build a fixture. Confident-sounding speculation about
  external behaviour is the failure mode that produces silent bugs
  and rolled-our-own reimplementations of canonical utilities.

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

## Current state

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
  Fleet state (as of 2026-04-30): the 685-protein fleet run was
  STOPPED after evaluation showed extractions were bad (chain-
  extraction issue in structure preparation). 200 proteins had
  completed Use Case B `PerFrameExtractionSet` extractions before
  the stop; those H5s exist on disk under
  `fleet_calibration-{working,stats,backup}/` but are flawed and
  should not be treated as authoritative. The other 485 were
  CANCELLED before completion.

  Recovery path: OF3 (OpenFold3) is generating fresh structures
  directly from sequence for the same 685 proteins, replacing the
  PDB-extraction step that produced the bad chains. Once OF3
  structures land, the 685-fleet re-runs MD on the new structures.
  These OF3 outputs have not been extracted yet because they are
  *direct* OF3 predictions — no extraction step has run on them.

  DFTs on the 485 have not been scheduled and are deferred until
  the structure-quality issue is resolved.
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
