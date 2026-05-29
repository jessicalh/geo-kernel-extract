# Next-session kickoff prompt

Start a fresh Claude Code session at `/shared/2026Thesis/nmr-shielding/`
(NOT inside `h5-reader/` — the project memory is keyed on the
top-level directory). Paste the block below verbatim as the first
message.

---

```
Resuming the h5-reader 5-TR surfacing project. The prior session
landed Phases A-H end-to-end + the 4 NOW cleanup items from
adversarial review; this session does the 10 LATER items as one
final polish pass before we ship.

PHASE 0 — STARTUP (do BEFORE writing any code)

1. Invoke the qt6-cpp skill so its references are loaded. Every
   QObject we touch this session needs CENSUS_REGISTER, every
   signal/slot connection needs ACONNECT, every thread-sensitive
   method needs ASSERT_THREAD, and the UDP log on port 9997 is the
   primary debug channel. These are non-negotiable per
   feedback_qt_discipline.

2. Per feedback_session_depth, read OBJECT_MODEL.md and PATTERNS.md
   IN FULL (3500+ and 1400+ lines respectively, chunked Reads if
   needed). Do not skim. Past sessions that skimmed produced chaos
   and re-rediscovered documented infrastructure.

3. Read the four h5-reader/notes/ documents that frame this work:
   - h5-reader/notes/PLAN_LATER_ITEMS_2026-05-29.md
     (the implementation plan you are executing this session)
   - h5-reader/notes/CLEANUP_AFTER_PHASE_AH_2026-05-29.md
     (the punch list that produced the plan; cross-linked)
   - h5-reader/notes/SCOPE_NEW_TRS_2026-05-29.md
     (the original scope doc — what the 5 TRs are + the
     visualisation framing decisions that are settled)
   - h5-reader/notes/REVIEW_CODEX_PHASE_ABC_2026-05-29.md
   - h5-reader/notes/REVIEW_AGENT_PHASE_ABC_2026-05-29.md
     (the two adversarial reviews that produced the NOW/LATER
     punch list; reference for context on WHY each LATER item
     matters and how the NOWs were already absorbed)

4. Read the previous session's plan file at
   /home/jessica/.claude/plans/validated-questing-metcalfe.md
   (the Phase A-H plan that produced commit 7a9d457).

5. Pull these memory entries by name so they are in working memory:
   - feedback_qt_discipline       (Qt + UDP debug discipline)
   - feedback_session_depth       (full-read of OBJECT_MODEL + PATTERNS)
   - feedback_extractor_untouchable (no producer-side changes, ever)
   - feedback_correlate_not_match (claim-shape convention for the
     thesis — applies to descriptor labels + status text)
   - feedback_methods_accumulate  (T2 sacred; no kernel-form-vs-
     grid-form retirement)
   - feedback_t2_sacred           (never collapse a tensor to scalar)
   - project_streaming_observer_for_sigma_pred
     (the deliverable framing; the Reorient J(ω) + tensor work in
     L-3 directly serves this)
   - feedback_no_simplification   (H ~ 20 dims, C ~ 6, N ~ 3, O ~ 12;
     present full element-dependent complexity)
   - feedback_no_abstractions     (no factories / no pluggable
     interfaces unless the user asks)
   - feedback_codex_review_invocation
     (codex exec --dangerously-bypass-approvals-and-sandbox
      --skip-git-repo-check - < /tmp/prompt.txt is the working
      reviewer call; prompt over stdin to avoid quoting hell;
      output via a file path you write into the prompt; wait for
      completion before reading)

6. Read the producer-side TRs (READ ONLY — never modified) so you
   know the H5 layout you are consuming:
   - src/IRedOrderParameterTrajectoryResult.{h,cpp}
   - src/KernelDynamicsTrajectoryResult.{h,cpp}
   - src/KernelCoherenceTrajectoryResult.{h,cpp}
   - src/ReorientationalDynamicsTrajectoryResult.{h,cpp}
   - src/DihedralAutocorrelationTrajectoryResult.{h,cpp}
   - src/TrajectorySpectral.h (shared math helpers)

7. Read the chassis files you will be modifying:
   - h5-reader/src/app/AbstractStripPanel.{h,cpp}  (the contract +
     paint helpers landed in Phase B)
   - h5-reader/src/app/StripStackWidget.{h,cpp}    (the panel
     container; L-2b restructures this)
   - h5-reader/src/app/DashboardDisplayController.{h,cpp}  (the
     panel router + denseH5Plan; touched by every LATER)
   - h5-reader/src/app/SceneRevealOverlay.cpp      (L-1d + L-3a)
   - h5-reader/src/io/QtTrajectoryH5.{h,cpp}       (L-1e defensive
     reader wrapping)
   - h5-reader/src/model/QtBondVectorBuffers.h     (the BondVector
     identity + composite TR structs; L-3a/b touch this)
   - h5-reader/src/model/QtPerAtomChannelBuffers.h (per-atom×channel
     curves / matrices)

PHASE 1 — PLAN MODE

After Phase 0 reading, enter plan mode and propose a phased plan
against PLAN_LATER_ITEMS_2026-05-29.md. The plan already specifies
L-1 (housekeeping, 5 items), L-2 (data extension + chassis cleanup,
2 items), L-3 (new rendering, 2 items), L-4 (feature add, 1 item).
Validate the phase boundaries and surface any open questions to the
user before exiting plan mode. The user has been actively engaged
in design decisions throughout this project — bring real decision
points to them rather than auto-deciding.

PHASE 2 — EXECUTION

Each phase ends with `ctest --output-on-failure` green from
`h5-reader/build/linux-debug` before moving to the next. The starting
baseline is:
- h5reader_model_tests (~60 cases)
- h5reader_scene_math_tests (14 cases)
- h5reader_rest_smoke (8 tests)
All three must remain 100% green. Add new test cases per phase as
appropriate (catalog presence row for chi descriptors, REST round-
trip for the new BondVector anchorToJson arm, etc.).

EXECUTION RULES (carry-over from the original 8-phase plan)

  - One stack widget, mixed panel types — L-2b unifies the
    collection if you choose to do it (read the plan; the trade-off
    is described there).
  - Chord diagram, not heatmap, for KernelCoherence (already
    landed; just don't second-guess).
  - Continuous power-spectrum line plot, not EQ bars, for
    KernelDynamics (already landed).
  - Bond-vector axis is first-class (sibling to atom/residue/ring) —
    BondVectorAnchor is (residue, kind), NOT an opaque index.
  - Residue-grouping ergonomics: a Residue anchor satisfies a
    BondVector descriptor via AnchorMatchesAxis widening — already
    landed; L-1d helper extraction makes the lookup symmetric across
    iRED and Reorient identity tables.
  - Sonification synth deferred to v2.
  - Per feedback_extractor_untouchable: zero producer-side changes.
    The 5 TRs are producer-complete; this work is entirely under
    h5-reader/.
  - Per feedback_qt_discipline: CENSUS_REGISTER / ACONNECT /
    ASSERT_THREAD / UDP port 9997 for debug. No QTimer on interactive
    controls. Class-driven everywhere.
  - L-1e adds try/catch wrappers on every Read* function in
    QtTrajectoryH5.cpp. The pattern is: catch HighFive::Exception +
    std::exception, log via the existing WarnGroupAbsent /
    WarnShapeMismatch helpers, leave the buffer null. Mirror this
    pattern verbatim — don't get creative.
  - L-3a (Mat3 tensor-glyph) is new VTK work. Use the existing
    SceneRevealOverlay actor lifecycle pattern (vtkSmartPointer<>
    everywhere, addAtom-style refcounting). Eigendecomposition via
    Eigen — pull the math into a separately-testable header so a
    scene_math_tests entry can exercise it.

ADVERSARIAL REVIEW HABIT

After each PHASE green-bar, consider whether the work warrants an
adversarial review. L-3 (new rendering) likely does; L-1 / L-2 don't
unless something feels off. When you do invoke codex, use the working
command from feedback_codex_review_invocation:

  codex exec --dangerously-bypass-approvals-and-sandbox \
             --skip-git-repo-check - < /tmp/codex-review-prompt.txt \
             2>&1 > /tmp/codex-review-output.txt &

(Run in the background, wait for completion notification, then read
the output file.) The prior session has examples of well-shaped review
prompts — see how REVIEW_CODEX_PHASE_ABC_2026-05-29.md was
constructed by reading /tmp/codex-review-prompt.txt-style framings
in the previous-session transcript if useful.

COMMIT POLICY

The previous session committed the Phase A-H + NOW + notes work as
commit 7a9d457 and the review artifacts as 7de457b. Commit at
natural phase boundaries (L-1 done, L-2 done, etc.) so the history
tells a clean story; or commit once at the end if the user prefers.
Ask before committing — Jessica explicitly approves each commit.

KNOWN UNCOMMITTED STATE (not yours to touch)

The repo has dirty files outside h5-reader/ from other parallel work
(papers/, references-meta/, src/ producer-side TRs that landed
in 7a9d457's parent, python/, learn/, etc.). Do NOT git add anything
outside h5-reader/ unless the user explicitly asks. Stage only your
own h5-reader/ changes when committing.

VERIFICATION BEFORE CLAIMING THE SESSION DONE

  cd /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-debug
  cmake --build . -j
  ctest --output-on-failure
  # 3/3 must pass

Manual (if a display is available; otherwise document as untested):

  ./h5reader --rest 9988 \
    /shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft &
  python ../../udp_listen.py  # in another terminal

  # In Add Signal dialog, verify the new mode chips appear and
  # the L-3 panels render (tensor glyph attached to bond, J(ω)
  # markers along log-x).

Final summary back to the user: which LATERs landed, which (if any)
slipped, test state, recommend next steps.

Begin Phase 0 reading.
```

---

**Notes on what's in the prompt vs. what isn't:**

- The prompt is fully self-contained — pasting it into a fresh
  session sets up Phase 0 reading, plan mode, execution, review
  habit, commit policy, and verification.
- It does NOT pre-decide the phase order — the plan file
  (`PLAN_LATER_ITEMS_2026-05-29.md`) already specifies L-1 → L-2 →
  L-3 → L-4, but the new session enters plan mode and validates
  with you before executing.
- It does NOT pre-load the codex review findings inline — instead it
  points at the two review docs in `h5-reader/notes/` that are now
  committed (so a context-cleared session reads them directly).
- The `feedback_codex_review_invocation` memory entry covers the
  exact working invocation, so the new session can re-run codex if
  L-3's new rendering work warrants a second adversarial pass.

**File locations the next session needs (all committed in
`7a9d457` and `7de457b`):**

```
h5-reader/notes/
├── SCOPE_NEW_TRS_2026-05-29.md          # input scope (decisions settled)
├── CLEANUP_AFTER_PHASE_AH_2026-05-29.md # NOW/LATER punch list
├── PLAN_LATER_ITEMS_2026-05-29.md       # the 4-phase plan for this session
├── REVIEW_CODEX_PHASE_ABC_2026-05-29.md # codex review
├── REVIEW_AGENT_PHASE_ABC_2026-05-29.md # Anthropic-agent review
└── NEXT_SESSION_KICKOFF_PROMPT.md       # this file (the prompt to paste)
```
