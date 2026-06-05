# Codex — IMPLEMENT Stage 1: the typed VisualizationDefinition registry

## (Preamble prepended. VERIFY FROM THE CODE BEFORE CHANGING ANYTHING, EVER — cite file:line. This IS a code change: build GREEN, do NOT run git.)

## What
Implement Stage 1 of the visualisation-registry rebuild, EXACTLY per the build-ready plan:
`notes/STAGE1_REGISTRY_PLAN_REVISED_2026-06-05.md`. That plan already resolved a code-grounded
critique (it opens with a CHANGES-VS-CRITIQUE table) — follow it precisely; it is the authority.
The registry is the offerability gate. Stage 1 = registry foundation + migrate STRIPS + make hollow
modes structurally un-offerable + the empty-reality-check (a runtime LOG, NEVER a startup fail).
Stages 2 (Atom Colour) and 3 (Tensor glyph) are LATER — leave their slots, do NOT build them.

## How — honour the plan's discipline
- Work the plan's ORDERED STEPS (wrap-then-peel). Build GREEN after EACH step:
  `cmake --build --preset linux-rwdi --target h5reader -j$(nproc)` plus the test targets the plan
  touches. Each step is designed to compile independently — keep it so.
- PRESERVE behaviour byte-for-byte where the plan requires: the `DashboardSelectionController`
  removal cascade (`signalsBeingRemoved_`); strip-history persistence (Stage 1 keeps refs keyed by
  the `displayModeId` STRING); the panel-ref tri-state so `static.tensor` / `static.spectrum.power`
  stay tracked-but-hidden (never picker-visible); the reorient composite-panel coordinator path; the
  `strip.*` prefix rule.
- Include-cycle fix exactly as specified: `DisplayModeCapabilityFor` becomes a DECLARATION; its body
  moves to a new `DisplayModeCapability.cpp` that includes the registry; the registry builds
  capabilities from an explicit table and never recurses through `DisplayModeCapabilityFor`. Add the
  new TU(s) to `h5reader_core` and the touched test targets.
- Startup validation: WARN-only over the UNFILTERED catalog (`allDescriptorList()`), hardened to a
  debug assert ONLY in the final step after hollow strings are peeled. The empty-reality-check logs
  to UDP / `smokeSummary`; locate-before-absent via the REAL availability API + the one minimal
  additive accessor the plan specifies (`recordForStoragePath`). Wire validation at
  `ReaderMainWindow.cpp` where the catalog is built + availability set (NOT `main_reader.cpp`).
- Qt idioms per the plan: the registry singleton is non-QObject/immutable/function-local-static —
  follow its lifetime/thread contract; `ACONNECT`/`ASSERT_THREAD`/`CENSUS_REGISTER` everywhere else
  they apply.
- Lead-preference defaults for the plan's two open choices: PER-TYPE definition files (extensible for
  the Stage 2/3 definitions) and add new sources PER-TARGET. Note these in your output.

## Output
- Build GREEN (`h5reader` + touched test targets). Run the model/app tests the plan touches if
  feasible; report pass/fail with the command. Do NOT launch. Do NOT run git.
- Append a `notes/BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md` line: what landed (registry, strips
  migrated, hollow modes un-offerable, empty-reality-check), which plan steps completed, the two
  lead-preference choices taken, and any deviation from the plan + why.
- If a step won't go green, or the plan contradicts the code, STOP and write the blocker — do not
  guess or sprawl. T2 stays sacred; the extractor/producer and H5 write path are untouched.
