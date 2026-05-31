# h5-reader chewer spec — read this list

Design layer for the h5-reader chewer (the bulk-eager substrate +
event-harness + pybind11 surface that runs Stage 2 PySR rediscovery
transforms). Settled across the 2026-05-29 / 2026-05-30 design
conversations + audit passes (Codex + Claude general-purpose + Opus
validation).

Before starting any chewer-development work — design discussion,
code, tests, transform authoring, agent invocation — read the
following in order. Skipping the qt6-cpp skill, or skipping the
conventions doc, has been the failure mode that produces silently
wrong substrate code.

## 0. Required pre-reading — the Qt6/C++ skill

The team's Qt6 idioms (CENSUS_REGISTER, ACONNECT, ASSERT_THREAD,
crash diagnosis, object census, structured logging, QPointer
discipline, signal/slot patterns, lifecycle) are codified in the
qt6-cpp skill. Read these before any chewer-related code or design
work:

- `/home/jessica/.claude/skills/qt6-cpp/SKILL.md` — top-level entry point
- `/home/jessica/.claude/skills/qt6-cpp/references/architecture.md` — signal/slot patterns, error bus, startup/shutdown
- `/home/jessica/.claude/skills/qt6-cpp/references/crash-diagnosis.md` — minidumps, object census, QPointer, sanitizers
- `/home/jessica/.claude/skills/qt6-cpp/references/qt68-lts.md` — if Qt6 version nuances matter to current work
- Other reference files (cmake, networking, model-view, 3d-vtk, windows-gotchas) as the work demands

The skill is the source of truth for "how Qt6 is done here." Substrate
calculators, the chewer's QObject + QThread layer, the cancel-token
plumbing, and the UDP progress logging all assume the skill's idioms.

## 1. h5-reader CLAUDE.md (the directory's rules)

Before reading the chewer-specific docs, refresh on the reader's own
discipline: what this directory owns vs doesn't, the
`feedback_qt_discipline` items, the cross-platform bar, the
hard-won viewer lessons, the build setup.

- `/shared/2026Thesis/nmr-shielding/h5-reader/CLAUDE.md`

Critical rules from there that bind the chewer:
- "**DO NOT modify any files outside of `h5-reader/`**"
- The reader never links the `nmr_shielding` library
- Every QObject ctor: `CENSUS_REGISTER(this)`
- Every signal/slot connection: `ACONNECT(...)`
- Every thread-sensitive method: `ASSERT_THREAD(this)`
- Full error handling at every external-library boundary

## 2. Reader-as-platform architecture

The big picture: what the chewer IS, where it sits, what it shares
with the viewer, what it adds. Substrate eager-precomputation tier
list (with audit-derived additions: per-pair distance tensor at 1P9J
scale, bond strain features, neighbor residue context, T2 tensor
operators, integrated dipolar sums, Markley pseudoatom exposure),
library/app duality, three triggers (CLI batch / REST / GUI toolbar
button), Qt-citizen chewer with the 3-level Python cancel pattern,
threading model.

- `reader_as_platform_2026-05-29.md`

## 3. Substrate conventions (Phase 1 GATE — read before any substrate code)

**Mandatory before writing any substrate accessor.** Both audit
agents independently demanded this doc exist first because three
high-impact features (spherical harmonics, local frames,
multipoles) silently disagree across transforms otherwise. Revision
2 reflects the Opus validation pass.

Every convention call: SH basis ordering (matches the LIBRARY's
`src/Types.cpp::SphericalTensor::Decompose`, NOT the h5-reader's
`SphericalDecomposition.cpp` — they differ; reasoning + line refs in
the doc), dipolar form `(3cos²θ−1)` for shielding kernels but `P₂`
for iRED/TCF per Lipari-Szabo, local frames per atom class with
edge-case policies (HID/HIE anchored on CG, TRP bridgeheads emit
constituent rings not perimeter, Gly via both diIndex + prochiral
R/S enums), three explicit charge sources no fallback chains,
residual baselines named (`literature_residual`, `stage1_ridge_residual`,
never `kernel_residual`), default-cutoff policy NONE (the "methods-paper
money item"), Cremer-Pople only for 5-rings (PRO), returns (Q,
theta), H-bond detection explicit (Larsen-geometric vs DSSP-Kabsch-Sander,
both coexist), symmetry-averaging via SDK helper not substrate
default, provenance sidecar schema, PBC verbatim port rule,
substrate version in Parquet filename, Markley pseudoatom names.

- `substrate_conventions_2026-05-30.md`

## 4. Substrate calculator framework

The pattern every substrate calculator follows: three-phase
lifecycle (`onLoad` / `onCompute` / per-calculator query surface),
cancel propagation via CancelToken, UDP structured progress logging
via the existing `StructuredLogger`, scheduler with topological
dependency ordering, file organization with one calculator per file,
readability discipline (no template metaprogramming, no registration
macros, no std::function callbacks, no central bindings.cpp).

Includes the DRY discipline (reuse h5-reader's existing typed
structures and loaders, never roll your own) and the h5-reader
foundational refactor required to enable that DRY reuse (Changes A
+ B — extracted `io::DftShieldingLoader` and `io::TrajectoryFrameMap`
static helpers; both landed in commit `518855e`).

- `substrate_calculator_framework_2026-05-30.md`

## 5. Stage 2 PySR campaign

The WHY for the substrate work. What gets walked, what equation
forms are tested (instantaneous, lagged, dynamics-averaged,
differential as cross-check, per-trajectory descriptors, spectral,
per-atom curves), kernel rediscovery scope (8-10 canonical kernels
in three regimes plus vibrational-correction and Wishart-Sykes-Richards
scope-gated additions), charge-source sensitivity, polarisability
sensitivity, MOPAC's two roles, AIMNet2's three roles, methodology
rigor (seed-stability runs, bootstrap CIs, seeded-vs-cold pairs,
residual-structure diagnostic, reference equation YAML library,
cross-form convergence, cross-tool agreement on flagship kernels,
unit consistency, standardised SR-result JSON schema, 1UBQ control,
reproducibility manifest, DFT component-decomposed targets per
Ramsey theory).

- `stage2_pysr_campaign_2026-05-29.md`

## Order rationale

- **Qt skill first** — the idioms bind everything; reading
  architecture docs without the idioms produces "looks plausible
  but breaks at compile time / discipline review" code.
- **h5-reader CLAUDE.md second** — the directory's own rules
  (`feedback_qt_discipline`, never link the library, never modify
  files outside h5-reader/) bind every chewer decision.
- **Reader-as-platform third** — the WHAT and the architectural HOW.
- **Conventions doc fourth** — the explicit decisions every
  substrate accessor depends on; reading them after the platform
  doc means the conventions land in already-existing architectural
  context.
- **Calculator framework fifth** — the pattern for adding any
  substrate feature; depends on the conventions to know what each
  calculator's contract is.
- **Stage 2 campaign last** — the WHY (the methods-paper goal these
  substrate features serve). Reading the campaign first risks
  reverse-engineering the architecture from the campaign's needs
  rather than reading the settled design.

## Memory entries to load alongside

- `project_reader_as_platform` — points at these docs; settled
  architectural decisions for the chewer.
- `project_three_stages` — Stage 1 / Stage 2 / Stage 3 framing; the
  chewer is Stage 2 infrastructure.
- `project_1p9j_study_system` — the specific 1P9J trajectory the
  chewer's first dogfood runs against.
- `project_calibration_done` — Stage 1 settled at R²=0.818; the
  chewer's Stage 1 rebuild script reproduces this through the new
  toolchain.
- `feedback_t2_sacred`, `feedback_kernel_not_shielding`,
  `feedback_identity_from_chemistry_not_position`,
  `feedback_no_simplification`, `feedback_not_predictor`,
  `feedback_model_not_physics`, `feedback_adversarial_review_physics`,
  `feedback_extractor_untouchable`, `feedback_methods_accumulate`,
  `feedback_residual_as_ml_feature`, `feedback_pbc_verbatim` — each
  one is referenced by the chewer docs as a binding constraint.

## What is NOT in this spec layer

- The viewer's design (the existing `Conformation` /
  `TrajectoryConformation` / `SingleConformation` /
  `DftShieldingStore` machinery). Those serve the viewer correctly;
  the chewer goes around them, not through them. See the calculator
  framework's "DRY discipline" section for the rationale.
- The producer's design (`nmr_extract` and `src/`). Read-only per
  `feedback_extractor_untouchable`.
- The strip UI internals (`AbstractStripPanel`, `StripStackWidget`,
  dashboard signal infrastructure). Off-limits per the
  reader-untouchable analog; the chewer does not modify these.
- Speculative future-work docs (MolProbity tool, frame-partition
  index, expert system, narrative engine, etc.) — those live in
  `analysis-speculative/` and are not part of the chewer's settled
  spec.

## Provenance

The four spec docs in this directory landed in
`analysis-speculative/` first (the project's gitignored scratch
workspace) during the 2026-05-29 / 2026-05-30 design conversations.
They moved here when the chewer architecture stabilized enough to
become the canonical reference for chewer development. The
`analysis-speculative/` versions are gone; this directory is the
single source of truth for chewer design.

Prior agent passes that shaped these docs (audit context, kept for
reproducibility):
- 2026-05-29 substrate-tiers audit (Codex + Claude general-purpose)
- 2026-05-29 methodology-additions agent (Claude general-purpose)
- 2026-05-30 reader robustness audit (Codex + Claude general-purpose)
  — landed as commits `f1acbd4`, `9a52473`, `b58b13a`, `8ed67e2`
- 2026-05-30 Opus authoritative pass on conventions — caught the
  `SphericalDecomposition.cpp` vs `Types.cpp` divergence, several
  other refinements
- 2026-05-30 Opus speculative-reach reframe (canonical-rediscovery
  lens correction) — kept most drops; surfaced symmetry-averaging
  SDK helper
- 2026-05-30 foundational h5-reader refactor — Changes A + B landed
  as commit `518855e` (DftShieldingLoader + TrajectoryFrameMap
  helpers extracted)
