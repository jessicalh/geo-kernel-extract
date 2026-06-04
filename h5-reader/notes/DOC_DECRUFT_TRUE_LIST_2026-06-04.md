# H5 Reader Docs Decruft And Make-True List - 2026-06-04

Scope: reader docs under `notes/` plus rediscover docs under `src/rediscover/*.md` and `src/rediscover/analysis/*.md`.

Constraints for this pass: docs-only review. No source/build/git operations were needed. This file is the only requested write target.

## Current Reality Anchors

Reader UI truth is anchored by `notes/UI_STATE_OVERVIEW_2026-06-04.md` and `notes/STABILISATION_FEATURE_EVAL_2026-06-04.md`.

- The reader is no longer a prototype: `.LGS` load path, `QtProtein`/`Conformation`, VTK molecule scene, playback, atom picking, ordered selections, camera composer modes, transforms, ribbon/ring/field overlays, inspector/dashboard docks, and loopback REST harness are present.
- Stabilisation is partly built: `TransformedConformation` provides display-space rigid transforms, and `CameraComposer`/typed camera modes provide camera lock/composition. The old "not implemented yet" statement is false.
- Local interaction display isolation is still not built: no pick atom plus radius plus hide nonlocal atoms/residues UI, no display-only local molecule actor, no local radius filter, and no neighborhood endpoint in `src/app`.
- Exact three-atom Kabsch local-frame design should not be treated as ready: the shared Kabsch path rejects rank-deficient/coplanar subsets, so exact three-atom local frames need a plane/wedge stabilizer instead.
- Linux source build works through presets, but installability/adviser-ready packaging is not present: no `install()`, package presets, CPack, AppImage, or linuxdeployqt bundle.

Rediscover truth is anchored by `src/rediscover/NOW.md`, the current top corrections in `src/rediscover/STATE.md`, and the 2026-06-04 maths audit docs.

- 1P9J is a within-instrument. WITHIN results stand: charge q/r3 within R2 about 0.28, ring current within R2 about 0.28, unified combine within R2 about 0.43.
- All prior 1P9J LOAO/between numbers are retracted. Old "charge 0.38", "ring 0.17", and "unified 0.26" values were within-modulation, not true leave-one-atom-out/between transfer.
- True LOAO on 1P9J between is approximately null or nonsensical for the current kernels. Between/transferability now belongs to the 720-WT static corpus.
- Field standalone is approximately null, but field remains a useful combine contributor. McConnell/hbond standalone are also approximately null.
- Maths audit code fixes are code-complete, while final verdict re-runs are parked. Do not quote old between/LOAO values as results.

## Priority Key

- P0 must-fix/must-cut: currently misleading or likely to steer the next session into wrong work.
- P1 make-true soon: stale enough to waste time or contradict current docs, but less likely to cause immediate bad implementation.
- P2 archive/cosmetic: historical cruft, stale session context, or harmless pointers.

Verdicts:

- KEEP: current enough; no urgent action.
- MAKE-TRUE: keep the doc, but fix named stale claims/status.
- CUT: remove as live documentation, or collapse to a one-line historical pointer if deletion is not desired.
- ARCHIVE: keep only as historical/session record; add or preserve a top banner saying it is not current truth.

## P0 Must-Fix And Must-Cut

| Priority | Doc or cluster | Verdict | Misleading or dangerous claim | Current reality | Specific action |
|---|---|---:|---|---|---|
| P0 | `notes/STABILIZED_LOCAL_VIEW_DESIGN_2026-05-29.md` | MAKE-TRUE | Says "not implemented yet" and presents exact three-atom Kabsch local view as viable next implementation. | Stabilisation infrastructure exists (`TransformedConformation`, camera composer/lock modes). Local radius/isolation is not built. Exact three-atom Kabsch is invalid with the current rank policy. | Rewrite the top status block: "partially implemented as transform/camera infrastructure; local interaction isolation remains open; exact three-atom Kabsch must be replaced by plane/wedge stabilizer." Link `STABILISATION_FEATURE_EVAL_2026-06-04.md`. |
| P0 | `notes/PLAN_LATER_ITEMS_2026-05-29.md` | CUT | Reads like the active next-session plan and lists LATER items as unbuilt. | Many listed items landed: `lookupBondVector`, Reorient tensor reveal, fixed-frequency panels, sequence bar, owned-panel status accounting. | Cut as active plan. If retained, collapse to archival pointer plus a landed/open matrix. Do not let future sessions execute this as current work. |
| P0 | `notes/NEXT_SESSION_KICKOFF_PROMPT.md` | CUT | Paste-ready prompt instructs execution of stale `PLAN_LATER_ITEMS` phases and old session assumptions. | Current work state is in June 4 state/eval docs, not this kickoff prompt. | Delete or replace with a one-line archived pointer to `UI_STATE_OVERVIEW_2026-06-04.md`, `STABILISATION_FEATURE_EVAL_2026-06-04.md`, and current issue/backlog docs. |
| P0 | `notes/VISION_AND_PROGRESS.md` | ARCHIVE | Calls itself a living tracker, predates current source, and has stale "progress so far" and warning sections. | Source now includes the main reader UI, loader, panels, camera/selection machinery, and REST loopback. Deferred feature notes are not desk-ready requirements. | Mark as historical progress ledger. Remove or banner stale warnings that suggest `SCOPE.md`/glossary are the only truth. Point to README, `SCOPE.md`, `UI_STATE_OVERVIEW_2026-06-04.md`, and `STABILISATION_FEATURE_EVAL_2026-06-04.md`. |
| P0 | `notes/RUN_DESCRIPTOR_TOML.md` | CUT | Body still describes a TOML descriptor/parser path and "not implemented" state. | The top says it is superseded by `spec/CALCSET_MANIFEST.md`; `.LGS` JSON landed and toml++ was removed. | Cut body as live schema. Retain only a historical note pointing to the calcset manifest and current README/load path. |
| P0 | `src/rediscover/STATE.md` | MAKE-TRUE | Top sections still state old Stage2 claims: field recovered-law, charge/ring/unified LOAO positives, old NEXT tasks. | `NOW.md` retracts all 1P9J LOAO/between numbers; field standalone is approximately null; audit code fixes are complete and re-runs are parked. | Split "current synthesis" from historical Stage2 notes. Banner all old LOAO/between/field-law claims as superseded. Update NEXT to 720-WT statics, Stage3, and parked clean re-runs only. |
| P0 | `src/rediscover/GUIDANCE.md` | MAKE-TRUE | Status says current fit result is Build3/live substrate Build1 and open work is between-calculator/equations grading. | Build4/Stage2/audit have moved the state forward. True LOAO retracted old between. 720-WT owns transfer. | Replace status with NOW-aligned state: within results stand; no 1P9J between claim; audit code-complete/re-runs held; next serious between work is 720-WT. |
| P0 | `src/rediscover/REDISCOVERY_MAP.md` | ARCHIVE | Old roadmap/status grid makes early arcs look current and orders work around pre-Stage2 assumptions. | Stage2/audit supersedes much of the roadmap; current control plane is `NOW.md`, corrected `STATE.md`, and active specs. | Archive as historical science roadmap. Add top banner "not current work catalog." If a map is still wanted, rebuild from NOW after the audit re-runs. |
| P0 | `src/rediscover/WORK_CATALOG.md` | MAKE-TRUE | Says earlier enabling conditions are met, references near-term 1P9J/750/engine ordering, and predates true-LOAO retraction. | Current next work is 720-WT static ingest/pilot, parked clean re-runs, Stage3, and ubiquitin-50fr. | Either update to a current backlog or archive. Remove any suggestion that old 1P9J between/LOAO numbers guide priority. |
| P0 | `src/rediscover/POSTMORTEM_STAGE2.md` and `src/rediscover/POSTMORTEM_STAGE2_3.md` | ARCHIVE | Immutable records contain now-dangerous field recovered-law and positive LOAO/unified transfer claims. | Later true-LOAO/maths audit retracted those claims. | Keep as run records only, but add a top supersession banner pointing to `NOW.md`, `POSTMORTEM_TRUE_LOAO_2026-06-04.md`, and `MATHS_AUDIT_CHECKLIST_2026-06-04.md`. |
| P0 | `src/rediscover/analysis/FINDINGS.md` | MAKE-TRUE | Has an older "partially superseded" banner but still presents old LOAO/leave-atoms-out evidence as if actionable. | 1P9J between/LOAO claims are retracted; within-only claims and qualitative pattern evidence can remain with caveats. | Add a June 4 supersession banner: "Do not quote LOAO/between results from this file; use NOW/STATE/audit docs." Keep only within/pattern material as current. |
| P0 | `src/rediscover/FRESH_LOOK_2026-06-03.md` and `src/rediscover/H5READER_COMPREHENSION_2026-06-03.md` | CUT | One-line pointer stubs remain after doc-truth cleanup. | `NOW.md` and `INDEX.md` already carry the pointer role. | Delete these stubs, or leave them only if external links need compatibility. They add navigation cruft. |

## P1 Make-True Soon

| Priority | Doc or cluster | Verdict | Stale claim or problem | Current reality | Specific action |
|---|---|---:|---|---|---|
| P1 | `notes/CLEANUP_AFTER_PHASE_AH_2026-05-29.md` | ARCHIVE | NOW/LATER cleanup list and uncommitted-state notes read like current instructions. | Several items were fixed or folded into later plans; working-tree/session state is stale by definition. | Mark as historical cleanup punch list. Add a short completion/disposition matrix if keeping it. |
| P1 | `notes/ROBUSTNESS_BACKLOG_2026-05-30.md` | MAKE-TRUE | Lists several robustness gaps as open: H5 reader exception boundaries, required positions validation, slider scrub reflection, settings persistence. | Reader-local exception helpers/per-reader try/catch exist; required positions validation hard-errors; scrub reflection and QSettings persistence are present. Some crash/render-probe/polish items remain open. | Add per-item status (`done`, `open`, `verify`) and trim fixed sections into a short changelog. Keep only true open robustness debt. |
| P1 | `notes/POLISH_BACKLOG.md` | MAKE-TRUE | Still useful, but some items are landed or partially handled. | UI overview says this is still the best polish backlog, with partial overlap against built camera/panel/measurement features. | Add status tags per item and move built items to a "landed" section. Keep as active polish debt. |
| P1 | `notes/VIEWPORT_OBSERVATIONS_2026-05-30.md` | MAKE-TRUE | Some camera-mode gaps are now obsolete. | Atom/bond/dihedral/subset camera composition and plane lock exist; local display isolation still does not. | Banner "camera sections partly superseded; local isolation observations still current." |
| P1 | `notes/SCOPE.md` | MAKE-TRUE | High-level scope holds, but loader/dependency notes still reflect early source layout. | `.LGS`/manifest path and current dependency shape changed. | Small update to data-loading/dependency paragraphs. Keep the boundary statement. |
| P1 | `notes/H5_READER_REWRITE_DESIGN_2026-05-23.md` | ARCHIVE | Status says draft/pending user review/no code lands. | The rewrite landed and evolved. | Add top banner "historical typed-mirror design; implemented/evolved." Point to current README/SCOPE/UI state for present architecture. |
| P1 | `notes/H5_FIELD_GLOSSARY.md` | MAKE-TRUE | Early "reader today" inventory is stale. | Panels, signals, overlays, and selection/camera affordances now exist beyond the April snapshot. | Keep as glossary, but update the status block and any "currently missing" claims before using it as an implementation reference. |
| P1 | `notes/FEATURE_PLAN.md` | ARCHIVE | Tentative April feature plan predates the actual implementation. | Current feature state is in UI overview/polish backlog. | Archive as idea bank; do not use as current roadmap. |
| P1 | `notes/BUILD_LAYOUT_PLAN_2026-05-23.md` | ARCHIVE | Draft build-layout plan predates current presets/modules. | Source build works; packaging/installability still missing. | Archive design history. Add note: current build docs live in README; packaging gap lives in UI overview. |
| P1 | `notes/SINGLE_POSE_AND_ORCA_DESIGN_2026-05-26.md` | MAKE-TRUE | Design status does not clearly say which pieces landed. | Single-conformation path, `SingleConformation`, and `LoadPose` concepts are implemented/evolved through `.LGS`. | Add implementation status and pointer to current manifest/load path. Keep design rationale. |
| P1 | `notes/SCOPE_NEW_TRS_2026-05-29.md` | MAKE-TRUE | Reads as forward design, though much landed. | Typed result-series/dashboard concepts are partly implemented. | Add landed/remaining status and point to UI overview/dashboard classes. |
| P1 | `notes/READER_INPUT_DISPLAY_INVENTORY_2026-05-28.md` | MAKE-TRUE | Inventory counts and drift notes may be stale. | Source and catalog evolved. | Add "snapshot only; verify counts before acting" banner. Keep as historical inventory. |
| P1 | `notes/QT_PRIMITIVES_ALIGNMENT_2026-05-30.md` | MAKE-TRUE | Refactor scope and "ALL committed work" language can read as current task state. | Some primitives/refactor intent may have landed or shifted. | Add status matrix and archive any stale session instructions. |
| P1 | `src/rediscover/INDEX.md` | MAKE-TRUE | Gate text says maths audit is "NOT done"; this is now too coarse. | Audit code fixes are complete; verdict re-runs are parked; between/LOAO remains not quotable. | Update gate language to "audit code-complete, verdict re-runs held; do not quote between/LOAO." Keep the navigation role. |
| P1 | `src/rediscover/MATHS_AUDIT_CHECKLIST_2026-06-04.md` | MAKE-TRUE | Open-agenda wording predates the final "code-complete/re-runs held" state. | NOW says fixes are code-complete and clean verdict re-runs are parked. | Add a final status block without rewriting the audit trail. |
| P1 | `src/rediscover/PARTITION_FILTER_DESIGN.md` | MAKE-TRUE | Status says live substrate Build1 and next step is to move dominance quantile bins into C++. | Build4/Stage2 dominance arc has moved forward. | Update result/status block and link Stage2/audit docs. |
| P1 | `src/rediscover/analysis/VARIANCE_DECOMPOSITION_METHOD.md` | MAKE-TRUE | Design-only between/variance language predates the within/between bug fix. | Old between decomposition was contaminated by within-modulation; true between needs 720-WT. | Banner as method sketch only. Add warning not to use it to justify old 1P9J LOAO/between claims. |
| P1 | `src/rediscover/APBS_VALUES_VERIFICATION.md` | MAKE-TRUE | Can be read as "APBS radii are placeholder/catastrophic" globally. | `APBS_STAGE1_RECONCILIATION.md` refines it: mutation/ORCA prmtop path used real radii; trajectory/TPR path used placeholders; effect remains path-specific. | Add supersession banner pointing to `APBS_STAGE1_RECONCILIATION.md`. Keep details as investigation record. |
| P1 | `src/rediscover/SPEC_INTERP_1P9J_SCOPING_2026-06-04.md` | MAKE-TRUE | Says `SPEC_LEARNING_MODEL_2026-06-04.md` is not present and cites pre-fix e3nn result. | The learning spec exists, and `BUILD_INTERP_RESULT_2026-06-04.md` reports the clean interpolation result. | Update dependency/status block; point to `BUILD_INTERP_RESULT_2026-06-04.md`. |
| P1 | `src/rediscover/SPEC_LEARNING_MODEL_2026-06-04.md` | MAKE-TRUE | Says companion 720/interp specs are not present. | Both companion specs exist, and the interpolation build result exists. | Fix companion-spec references and add "interp built" status. |
| P1 | `src/rediscover/PRECIS_ADVISOR_2026-06-04.md` | MAKE-TRUE | Draft status says Stage3/interp pending. | Interp result was built; Stage3/between remain pending. | Refresh draft status before sharing with advisor. Include the clean interp metrics and the 1P9J-within caveat. |
| P1 | `src/rediscover/ENGINE_TOTALITY_DESIGN.md` | ARCHIVE | Top says proposal/no implementation. | Functional API and later placement/filtering work were built/evolved. | Archive as design history and link `MODEL_PLACEMENT_PROPOSAL.md`, `PARTITION_FILTER_DESIGN.md`, `STATE.md`, and current code. |
| P1 | `src/rediscover/MODEL_PLACEMENT_PROPOSAL.md` | MAKE-TRUE | Proposal/no-code status predates later e3nn protocol fix and interpolation work. | Later docs/code changed placement/protocol details. | Add implemented/superseded status and link the e3nn protocol fix and interp result. |
| P1 | `src/rediscover/PER_ATOM_SUBSTRATE_SPEC_2026-06-02.md` | MAKE-TRUE | Status only says landed/extended through Build1. | Later Build4/Stage2/audit moved substrate interpretation. | Update status to reference Build4/Stage2 and current limitations. |
| P1 | `src/rediscover/ALLATOM_FIT_SPEC_2026-06-03.md` | MAKE-TRUE | Status stops at Piece1/Piece2/alpha. | Build3/Build4/Stage2 supersede parts of the fit story. | Add "landed and superseded by Build3/4/Stage2" banner. |

## P2 Archive, Keep, Or Cosmetic Cleanup

| Priority | Doc or cluster | Verdict | Why | Specific action |
|---|---|---:|---|---|
| P2 | `notes/UI_STATE_OVERVIEW_2026-06-04.md` | KEEP | Current reader baseline for this pass. | Keep as current until a newer UI state note replaces it. |
| P2 | `notes/STABILISATION_FEATURE_EVAL_2026-06-04.md` | KEEP | Current stabilisation/local-view baseline. | Keep as current; it should supersede the May 29 stabilized-local-view design. |
| P2 | `notes/RESIDUAL_RENDER_DROP.md` | KEEP | Still describes an open diagnostic/polish issue referenced by current state. | Keep, but verify probe line numbers before acting. |
| P2 | `notes/OBSERVABLE_MIGRATION_BRIEF_2026-05-28.md` | ARCHIVE | First safe slice for old dashboard/observable migration. | Archive as design history; current dashboard state lives in UI overview. |
| P2 | `notes/STRIP_CALCULATION_CLASS_TREE_DESIGN_2026-05-28.md` | ARCHIVE | Long design likely superseded by implemented dashboard/panel model. | Archive with banner saying implemented class tree may differ. |
| P2 | `notes/TIME_SERIES_EXPANSION.md` and `notes/PLANNED_ANALYSIS_METHODS.md` | ARCHIVE | Useful idea bank, but old "next" language should not drive current work. | Mark historical/idea-bank only. |
| P2 | `notes/REVIEW_CODEX_PHASE_ABC_2026-05-29.md` and `notes/REVIEW_AGENT_PHASE_ABC_2026-05-29.md` | ARCHIVE | Review records. Some NOWs were fixed or dispositioned. | Preserve as review history; add disposition pointer if missing. |
| P2 | `notes/REVIEW_BRIEF_BUFFERING_2026-05-29.md`, `notes/REVIEW_BRIEF_RESULT_GROUPS_2026-05-29.md`, `notes/REVIEW_BRIEF_SPINE_AND_CATALOG_2026-05-29.md` | ARCHIVE | Review briefs from an earlier implementation slice. | Archive; do not use as current task prompts. |
| P2 | `notes/SESSION_HANDOFF_2026-05-26.md`, `notes/SESSION_HANDOFF_2026-05-27.md`, `notes/SESSION_2026-04-19.md` | ARCHIVE | Session-history/paste-prompt cruft. | Keep only as historical handoff records; add banner "not current kickoff." |
| P2 | `notes/nmr_forensics/**` | KEEP | Audit corpus/data records, not reader UI design docs. | Leave alone. Optionally add an index note that this subtree is historical analysis evidence, not current reader docs. |
| P2 | `src/rediscover/NOW.md` | KEEP | Current live marker. | Keep current; this should remain the first rediscover truth source. |
| P2 | `src/rediscover/NEXT_SESSION_PROMPT.md` | KEEP | Current continuation prompt after maths audit. | Keep, but refresh after parked re-runs or 720-WT pilot. |
| P2 | `src/rediscover/PHYSICS_ARCHITECTURE_UNIFICATION_2026-06-04.md` | KEEP | Current physics architecture synthesis. | Keep. |
| P2 | `src/rediscover/SPEC_MATHS_FIX_AND_REVIEW_2026-06-04.md` and `src/rediscover/SPEC_720_STATIC_INGEST_2026-06-04.md` | KEEP | Current design/run specs. | Keep; mark completed once the corresponding work lands. |
| P2 | `src/rediscover/BUILD_INTERP_RESULT_2026-06-04.md`, `src/rediscover/NOTES_INTERP_LIVE_PREDICTION_2026-06-04.md`, `src/rediscover/LEARN_STAGE1_UNDERSTANDING_2026-06-04.md` | KEEP | Current June 4 result/understanding notes. | Keep as current result records. |
| P2 | `src/rediscover/DESIGN.md`, `src/rediscover/SURFACE_DESIGN.md`, `src/rediscover/STUB_LANGUAGE.md` | ARCHIVE | Already present as historical/partly-built design sketches. | Keep banners; do not promote to current implementation truth. |
| P2 | `src/rediscover/CODEX_BRIEF_BUILD*.md`, `src/rediscover/CODEX_BRIEF_PIECE*.md`, `src/rediscover/POSTMORTEM_BUILD*.md`, `src/rediscover/POSTMORTEM_PIECE*.md`, `src/rediscover/POSTMORTEM_STAGE2_1.md`, `src/rediscover/POSTMORTEM_STAGE2_2.md` | ARCHIVE | Run prompts and run records. Valuable provenance, not current requirements. | Archive cluster; add a top-level index note that current truth is NOW/STATE, not old postmortems. |
| P2 | `src/rediscover/CODEX_BRIEF_DOC_TRUTH_PASS_2026-06-03.md` and `src/rediscover/POSTMORTEM_DOCTRUTH.md` | ARCHIVE | Previous doc-truth pass records. | Keep as provenance; current decruft work supersedes the brief. |
| P2 | `src/rediscover/POSTMORTEM_TRUE_LOAO_2026-06-04.md`, `POSTMORTEM_DIAPARA_CHECK_2026-06-04.md`, `POSTMORTEM_E3NN_PROTOCOL_FIX_2026-06-04.md`, `POSTMORTEM_LOAO_FIX_2026-06-04.md`, `POSTMORTEM_MATHS_WALK_2026-06-04.md` | KEEP | Current audit postmortems. | Keep; these are the superseders for old Stage2 claims. |
| P2 | `src/rediscover/APBS_RADII_AB_WORKAROUND.md` and `src/rediscover/APBS_STAGE1_RECONCILIATION.md` | KEEP | Current APBS/radii reconciliation notes. | Keep; use reconciliation as the high-level source. |
| P2 | `src/rediscover/APPLIED_MATHS_AUDIT.md`, `src/rediscover/APPLIED_MATHS_AUDIT_codex.md`, `src/rediscover/FIXES_AUDIT_opus.md` | ARCHIVE | Earlier audit/review records superseded by June 4 maths walk/checklist/fix postmortems. | Archive with pointer to current audit docs. |
| P2 | `src/rediscover/analysis/PATTERNS.md` | KEEP | Still useful as pattern index. | Keep, but update any pointers that still rely on archived `REDISCOVERY_MAP.md`. |
| P2 | `src/rediscover/analysis/ENV.md` | KEEP | Environment note. | Keep, but verify paths/tool versions before reuse. |
| P2 | `src/rediscover/analysis/BACKBONE_LAW_EVIDENCE.md`, `EFG_ARC_EVIDENCE.md`, `JOHNSON_BOVEY_REGION_RECOVERY.md`, `MCCONNELL_DCHI_CALIBRATION.md`, `MCCONNELL_LITERATURE_DECIRC.md`, `RING_LITERATURE_DECIRC.md`, `STATIC_ENVIRONMENT_CALIBRATION.md`, `LITERATURE_FIXED_DECIRC.md`, `HANDOFF_BACKBONE_FIT.md` | ARCHIVE | Evidence and handoff records. Some claims may be old, but the docs are not the control plane. | Mark as evidence archives. Add "check NOW/STATE before quoting quantitative claims." |
| P2 | `src/rediscover/analysis/venv/**` markdown/docs | CUT | Third-party vendored environment docs pollute repo-doc surveys. | Exclude from reader/rediscover doc hygiene. If tracked in repo, move environment out of source tree or ignore generated environment docs. |

## Cross-Doc Cleanup Rules

Apply these rules when editing the docs above:

1. Any doc that says a feature is "not implemented" must be checked against `UI_STATE_OVERVIEW_2026-06-04.md` before it remains in live docs.
2. Any rediscover doc with positive 1P9J LOAO/between claims must be bannered as superseded unless it is one of the true-LOAO/audit correction docs.
3. Any paste-ready "next session" prompt older than 2026-06-04 should be archived or deleted as a prompt, not left as executable guidance.
4. Historical postmortems should not be rewritten as if they never happened; add supersession banners instead.
5. Design docs can stay if they answer "why", but their top block must say whether the design is unbuilt, partly built, built-and-evolved, or invalidated.

## Summary Counts

Counts are action rows/clusters in this file, not raw Markdown file count. Large families such as old Build/Piece postmortems, `notes/nmr_forensics/**`, and `analysis/venv/**` are intentionally grouped.

| Verdict | Count |
|---|---:|
| KEEP | 13 |
| MAKE-TRUE | 25 |
| CUT | 5 |
| ARCHIVE | 19 |
| Total rows/clusters | 62 |

Top must-fix:

1. `notes/STABILIZED_LOCAL_VIEW_DESIGN_2026-05-29.md`: false "not implemented yet" plus invalid exact-three-atom Kabsch path.
2. `src/rediscover/STATE.md`: old Stage2/LOAO/field claims mixed into the current control document.
3. `src/rediscover/GUIDANCE.md`, `REDISCOVERY_MAP.md`, and `WORK_CATALOG.md`: old roadmap/status documents can steer work back to retracted between/LOAO claims.
4. `src/rediscover/POSTMORTEM_STAGE2.md`, `POSTMORTEM_STAGE2_3.md`, and `analysis/FINDINGS.md`: quantitative claims need explicit supersession banners before anyone quotes them.
5. `notes/VISION_AND_PROGRESS.md` and `notes/RUN_DESCRIPTOR_TOML.md`: high-visibility docs still carry pre-source or superseded load-path truth.

Top must-cut:

1. `notes/NEXT_SESSION_KICKOFF_PROMPT.md`: stale executable prompt.
2. `notes/PLAN_LATER_ITEMS_2026-05-29.md`: stale active plan whose main items partly landed.
3. `notes/RUN_DESCRIPTOR_TOML.md`: superseded TOML schema/parser body.
4. `src/rediscover/FRESH_LOOK_2026-06-03.md` and `H5READER_COMPREHENSION_2026-06-03.md`: pointer stubs duplicated by `NOW.md`/`INDEX.md`.
5. `src/rediscover/analysis/venv/**`: third-party/generated docs should not be part of repo doc hygiene.

## ADDENDUM (Claude supplement, 2026-06-04) — for the executing pass

The list above is thorough and accurate. Supplements + corrections before execution:

**Already partly done (don't redo):**
- `notes/STABILIZED_LOCAL_VIEW_DESIGN_2026-05-29.md` — the false "not implemented" is **already trued** (the
  Status block now reads SUPERSEDED + links the eval). Execution should do **only** the remaining piece: note
  that the exact-three-atom Kabsch path is invalid and must be a plane/wedge stabiliser. Leave the status block.

**Docs the notes/+rediscover scope missed — record:**
- `doc/thesis-overview/**` (PRECIS_ADVISOR + learning-model spec + Stage-1 understanding + interp notes +
  capability graph + README) — **KEEP**, current; the advisor overview stashed today for the lead to finalise.
- `h5-reader/notes/SESSION_CHECKPOINT_2026-06-04.md` — **KEEP**, the cross-session state insurance written today.
- New session briefs/reviews: `src/rediscover/CODEX_BRIEF_SPEC_*_2026-06-04.md`, `CODEX_BRIEF_720_BUILD_2026-06-04.md`,
  `CODEX_BRIEF_ADVERSARIAL_720_2026-06-04.md`, `REVIEW_opus_SPEC_*_2026-06-04.md`,
  `REVIEW_codex_adversarial_SPEC_720_2026-06-04.md` — **ARCHIVE** (session provenance, like the older
  CODEX_BRIEF/POSTMORTEM cluster; current truth is the SPEC_*/NOW/STATE).

**OUT OF SCOPE — do NOT touch this pass:**
- `doc/calculators/**` (CONTINUE_PROMPT, HANDOFF, BACKGROUNDERS_SEED, OUTLINE) — a SEPARATE calculator-exposition
  docs effort (the lead's concurrent track), not reader/rediscover decruft. Leave untouched; flag for the lead.
- Code, `analysis/venv/**`.

**Refinement:**
- `src/rediscover/PRECIS_ADVISOR_2026-06-04.md` — its Stage-3 section was rewritten today (the un-phased
  two-model e3nn pipeline + clean interp metrics), so the make-true is mostly done; only refresh the remaining
  draft-status line. Keep it in sync with the `doc/thesis-overview/` copy (the lead's working copy).

**Execution discipline:** ADD banners / fix stale claims — do NOT erase history (postmortems keep their content
+ a supersession banner). CUTs collapse to a one-line pointer, not a hard delete. DOCS ONLY; NO git (the lead
reviews + commits the decruft as ONE change set).
