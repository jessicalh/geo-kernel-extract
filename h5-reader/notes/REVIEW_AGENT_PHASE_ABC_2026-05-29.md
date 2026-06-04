# Anthropic Agent adversarial review of Phase A/B/C (2026-05-29)

> **Historical review record — not current truth (trued 2026-06-04).** Preserve
> as review provenance; do not use as a current task prompt. Check
> `UI_STATE_OVERVIEW_2026-06-04.md` and current backlog/issues first.

Second-pass independent review of the Phase A (vocabulary), Phase B
(chassis), and Phase C (iRED end-to-end) changes — run in parallel
with the Codex review (see `REVIEW_CODEX_PHASE_ABC_2026-05-29.md`).
Tagging convention: **[NOW]** = must fix before Phase D builds on
top, OR is a real safety/correctness issue, OR will compound across
Phases D-G; **[LATER]** = can wait until Phase H polish.

## NOW

**1. iRED H5 reader has no per-array length guard.**
`h5-reader/src/io/QtTrajectoryH5.cpp:701-732` sets `n_vectors=M` from
`s2_ired` then unconditionally reads `residue_index`, `n_atom`,
`h_atom` into vectors whose actual length is whatever HighFive
resizes them to. `QtBondVectorTable::rowFor`
(`QtBondVectorBuffers.h:44-52`) walks `i<n_vectors` indexing
`residue_index[i]` and `kind[i]`; `SceneRevealOverlay` indexes
`tail_atom[*row]`/`head_atom[*row]`. If any of those three datasets
is shorter than M (malformed producer, partial write), the lookup is
out-of-bounds UB. **Fix:** assert `dims[0]==M` for all four datasets
and bail with `WarnShapeMismatch`. Phases D-G clone this reader
pattern — fix the template now.

**2. `denseH5Plan` iRED branch is unreachable today but will
mis-route when ReorientationalDynamics lands.**
`DashboardDisplayController.cpp:875-890` dispatches purely by
`descriptor.storagePath`, with no check on `displayModeId`/`channel`.
Phase E adds a Reorient descriptor that also keys to
`BondVectorAnchor` but produces multiple per-vector scalars (S², τ_e,
R1, R2, NOE). The current "match path → return s2_ired" shape doesn't
generalise — Reorient will need channel-id dispatch within the
branch. **Fix:** add a comment noting the contract, or sketch the
channel-id switch now so Phase E doesn't ad-hoc bolt one on.

**3. `BondVectorAnchor` widening from `ResidueAnchor` passes
`canBind` but the sampler returns `AnchorUnavailable`.**
`DashboardSignal.cpp:236-240` accepts a `Residue` anchor against a
`BondVector` requirement; but `denseH5Plan` (line 883) uses
`std::get_if<BondVectorAnchor>(&anchor)` and immediately bails with
`AnchorUnavailable` when handed the bare residue. So users who bind
iRED to a residue see "anchor unavailable" gaps with no fallback.
**Fix:** in the iRED sampler branch, if `anchor` is `ResidueAnchor`,
call `identity.rowFor(residue, /*kind=*/0)`. Same pattern will
reappear for Reorient.

**4. SceneRevealOverlay BondVector handling is keyed to iRED only
and will fight Phase E's Reorient table.**
`SceneRevealOverlay.cpp:188-211` and `:261-287` walk
`h5->iredOrderParameters()` exclusively. When Phase E adds Reorient,
an anchor with `kind=2` (CaHa) will silently miss iRED's NH-only
table and the reveal will be empty rather than checking Reorient.
**Fix:** extract a `lookupBondVector(h5, residue, kind) → optional<{tail, head}>`
that walks all known BondVector tables in priority order; both
`atomsForBinding` and `reveal()` consume it. Land before Phase E to
avoid two-table-bolt-on regret.

**5. Multi-mode signals only render one panel.**
`DashboardDisplayController.cpp:1396-1405` walks the signal once and
emits a `SequenceBarPanel` only when
`signal.binding.displayModeId == "static.bar.sequence"`. If
`displayModeIds` contains both `static.bar.sequence` and
`static.table` but the active `binding.displayModeId` is
`"static.table"`, `hasPanelMode` is true but no panel is built.
Compare this to the temporal path, which loops over all modes in
`displayModeIds`. **Fix:** loop
`for (const QString& mode : signal.displayModeIds) if (isPanelMode(mode))`
and dispatch on `mode`, not on `binding.displayModeId`. The current
single-binding check will silently break the moment a user enables
two display modes on one iRED signal.

**6. `samplingStatus = SampleStatus::Valid` is misleading for a
static-only descriptor.** `TrajectorySignalCatalog.cpp:535` forces
Valid; but `canSample()` (line 908-912) requires
`descriptor->temporal && displayModeId.startsWith("strip.")` — both
false for iRED — so `canSample` returns false. Functionally safe, but
`filterDescriptors` uses `samplingStatus != Valid` as the "exclude
pending" filter (line 843), so iRED now slips through filters that
should exclude static-only descriptors. Either intentional (then
document) or wrong (then add a new `SampleStatus::StaticOnly` rather
than overloading `Valid`). Whichever way, decide before Phase D-G
clone the bypass four more times.

**7. `TemporalStripPanel` / `SpectrumStripPanel` hold dangling
`const Track&` if `QVector<Track>` reallocates mid-paint.**
`StripStackWidget.cpp:60,285` stores `const Track&`; the panels are
constructed on the stack inside
`paintEvent`/`revealAt`/`mousePressEvent` from `tracks_[i]` and live
only for that call — safe today. But the ref is captured in a
long-lived `const`-ref field, so if someone re-uses the panel across
a `setTracks()` boundary (which Phase D-G might to keep VTK state),
the ref dangles silently. **Fix:** before Phase E lands the first
stateful panel, switch `track_`/`SpectrumTrack` to value-copy
(cheap, two pointers and a QColor) so the lifetime invariant is
local.

## LATER

**8. `mutable std::vector<std::unique_ptr<AbstractStripPanel>>
ownedPanels_` smell.** `DashboardDisplayController.h:183` is
`mutable` but the only mutating user is `takeOwnedPanels()`, which is
non-const. The `mutable` is unused — drop it. (If `takeOwnedPanels()`
was meant to be const, that hides a real const-correctness violation
since it moves out of the member.)

**9. `BindingForRow` lambda copies `rows` per click.**
`DashboardDisplayController.cpp:1459-1470` captures `rows` (the whole
vector) by value into the lambda; `std::function` then stores it
once. Each invocation reads the captured copy — not a per-click clone
— so cost is fine. But the panel also stores `rows` in `rows_`
(move-constructed), so the row data lives twice (panel + lambda).
For ~50 residues this is trivial; for ReorientationalDynamics (~3×50
vectors per protein) still trivial. Skip unless profiling says
otherwise.

**10. `static.curve.lag.animated` mode string is asymmetric with
siblings.** `DashboardDisplayController.cpp:53` — three modes are
noun-form (`bar.sequence`, `spectrum.power`, `chord.coupling`); this
one trails `.animated`. Phase D will introduce it; either drop
`.animated` or document why this one carries an adverb.

**11. Render order `tracks → spectrum → owned` is hard-coded for
Phase E heterogeneous interleaving.** `StripStackWidget.cpp:564-581`
paints in three sequential blocks. Phase E's "ellipsoid glyph next to
its strip" composition will want a single ordered panel list. Phase
B's shortcut is fine for C; revisit before Phase E mixes panel
kinds.

## Disposition (post-session)

- NOWs #1, #3, #5 applied inline before Phase D as templates so the
  clones didn't carry the defects.
- NOW #4 was deferred to Phase E itself (added the second-table
  lookup directly in SceneRevealOverlay when Reorient landed).
- NOWs #2, #6, #7 + all LATERs moved to
  `PLAN_LATER_ITEMS_2026-05-29.md` as L-1 / L-2 cleanup work.
- NOW #6 partially: the `SampleStatus::Valid` bypass survived
  because the dialog filter was extended through the new mode-kind
  registration (L-1d / NOW-2 in the cleanup pass) rather than via a
  new SampleStatus.

The independent Codex review (sandbox-bypassed) caught a different
set of issues. See `REVIEW_CODEX_PHASE_ABC_2026-05-29.md`; the
critical drain-on-tick bug was Codex's #1.
