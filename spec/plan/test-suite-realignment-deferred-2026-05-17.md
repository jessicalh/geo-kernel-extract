# Test-suite realignment: deferred deep cleanup, short-fix landing 2026-05-17

**Status:** Decision recorded 2026-05-17 after a long conversation triggered by
an adversarial review of two TR-queue commits that came back focused on test-
framework cruft instead of physics. The short fix is landing the same day;
the deep cleanup is deferred until after the TR queue completes.

**REVISION (later same day, 2026-05-17):** After short fix landed at commit
`a068818`, Codex review (manual, science-focused) caught five data-shape
issues in the Welford TR family. The bet as originally offered assumed the
TR data shape was correct and only the test framework needed deferring; that
assumption was wrong. **Data-shape cleanup is NOT deferred** — see
`spec/plan/welford-data-shape-design-2026-05-17.md`. The TR queue is paused
~12-18 hours while the Welford expansion lands. Test cleanup remains
deferred per the original bet here. The "tomorrow night" TR-queue
completion estimate slides accordingly; user explicitly endorsed the extra
time with the framing "we are lucky to have caught this and squared away
on goals."

**Rollback point:** git tag `pre-test-config-consolidation-2026-05-17` on
remote `origin` at commit `609a59b`. If the deferred-cleanup bet turns sour
and physics bugs roll through, that tag is the recovery anchor.

---

## What surfaced today

An adversarial code review of commits `bdeedf9` (HM + McConnell Welford TRs +
Coulomb-TR strike) and `609a59b` (Welford siblings × 3: Eeq + Sasa +
HBondCount) returned with one HIGH-severity finding from each of two parallel
reviewers: **the new Welford TRs read source-calc fields without checking
`conf.HasResult<SourceCalc>()` per-frame, so a frame where the source calc
was skipped (e.g., DSSP starved → HBondResult dropped → all-zero source
field accumulated) silently feeds default zeros into the running mean.** The
reviewers recommended adding `HasResult<>` guards across all five Welfords
plus a per-frame `source_present_per_frame_` mask in the H5 emission — the
same machinery the existing twelve Tripeptide/Larsen TRs carry.

The finding was technically correct. It was also wrong-shaped, and the
wrong-shape was Jessica's diagnostic moment.

The trap the reviewer flagged only fires through a misconfigured
`RunConfiguration`: a test (or future user) sets `skip_dssp=true`,
HBondResult quietly fails to attach, HBondCount Welford accumulates zeros
on the now-default `hbond_count_within_3_5A` field. In production
(`PerFrameExtractionSet`), DSSP is on; HBond attaches; the trap doesn't
fire. The reviewer's recommended fix — add per-frame guards everywhere —
is defensive machinery against a misconfiguration that the test
infrastructure should not have been able to express in the first place.
The HBondCount test set `skip_dssp=true` because that was the template I
copy-pasted from HM/McConnell Welford tests, which themselves copy-pasted
from a pattern I built up across the TR-queue session. The recursive
authorship is the point.

The "right" fix per the existing two-use-cases discipline
(`feedback_calculator_inclusion_two_use_cases` — codified after 100+
sessions of pushback against dependency-graph machinery) is not gates.
It is: the source is unconditional in the only production config we ship,
the test should have used that config or one shaped like it, and the
per-test bespoke skip-flag combinations are themselves the machinery that
this pattern of review attention is correctly identifying as cruft.

## The pattern Jessica named

In Jessica's words:

> "I have never acquiesced to a pile of complexity built entirely by you
> against design in a project, over six months of working with you, and
> ever had it be OK -- I have often gone along and then discovered the
> tail is wagging the dog. There is always a smart reason for the initial
> choice by you but then the growth in complexity, and the existence of
> a hidden agenda with growing complexity, takes up the attention that
> should be on the work and I end up with something that is not so great,
> because there is only so much attention to go around on a project of
> this size and complexity."

And:

> "You have created, we have created, over the course of this project,
> an unmaintainable set of 100s of tests. Those tests are half cruft. ...
> The code is unsharable because to understand it requires grokking 100s
> of ai tests done at a particular level of knowledge that no longer
> apply. I don't even mind having a dependency manager if we've agreed
> on it, even though there are only four invocations and only two sets
> of required calculator sequence -- but perhaps that manager should
> live in that reality."

And critically:

> "your training and anthropic are right about the many tiny unit tests
> pattern. It helped us get here. It is no longer helping us, and the
> control points required for it may no longer be helping the code."

The four-test ritual per TR (Frame0Semantics + FinalizeIdempotency +
H5RoundTrip + Integration1P9J) costs ~150 lines per TR and gets ritually
copied with `s/bs_/hm_/g`-style substitutions. Across the ~30 TRs landed
to date, that is ~4500 lines where three of the four tests answer the
same lifecycle questions the AV/FO patterns already establish in their
exemplars. The fourth is the only one that touches the physics. The
other three are tax paid to my abstraction, not to the science.

The skip-flag arithmetic across ~20 test files projects a configuration
**surface** of dozens of distinct test configs onto a configuration
**space** of size two: the production canonical and (when it lands) the
mutant-evaluation path. The dependency manager doesn't need to live in
the larger space; it can live in the two-config space and refuse to
express the rest.

## The trade-off

Pure schedule math at the current pace:

| Path | TR-queue completion | Net new test cruft | Time investment |
|------|---------------------|---------------------|-----------------|
| Roll through queue, defer cleanup | tomorrow night | ~1500 more lines | ~6-12 hrs queue work |
| Deep cleanup first, then queue | +1 week, then ~3 hrs queue | net -1500 to -4500 lines | ~30-40 hrs cleanup + 3 hrs queue |
| **Short fix now, queue, defer cleanup** (chosen) | tomorrow night | ~250 more lines, no growth after | ~2 hrs short fix + 6 hrs queue |

The bet selected: **short fix + queue completion + deferred deep cleanup.**
The risk that the bet is unwise: physics bugs may roll through during the
queue while adversarial review attention is consumed by framework cruft.
The mitigation, **load-bearing on the bet's viability**: a science-focused
adversarial reviewer is added to the review prompt for the remainder of
the TR queue. Without that, the bet is bad. With it, the bet is reasonable.

## The bet, recorded explicitly

1. **Tag and push** — done: `pre-test-config-consolidation-2026-05-17` at
   commit `609a59b`, pushed to `origin`. Rollback anchor.
2. **Short fix lands today** — consolidate the per-test skip-flag
   combinations across the five Welford tests written today into a single
   `tests/TestConfig.{h,cpp}` profile system, adopt for new TR tests going
   forward. ~2 hours, ~250 lines added in one place, ~75 lines removed
   from today's five test files. Old test files NOT swept — lazy
   migration when touched.
3. **TR queue continues** — remaining ~12 clone-shape TRs (ApbsEfg, Water
   + Hydration trio, DSSP family, Bonded/Gromacs energy, AIMNet2 embedding
   opt-in). Estimated tomorrow night at current pace.
4. **Open architectural item — water:** the Water /
   Hydration TRs may want NPY emission rather than H5 group emission.
   Jessica flagged this; needs discussion before that batch lands. Recorded
   here so it doesn't get lost.
5. **Codex review** of existing work happens off-stream. Jessica drives.
6. **Adversarial review prompt update — load-bearing.** At least one
   reviewer for the remainder of the queue must be science-focused, not
   pattern-focused. Pattern-focused review caught the test-framework cruft
   we're now addressing; science-focused review is the only mechanism that
   catches physics bugs (sign conventions, units, tensor-channel
   selection, formula errors). The bet's viability depends on this
   review actually being conducted, not just declared.
7. **Deep cleanup deferred** — see "What deep cleanup would look like"
   below. Estimated ~1 week of focused work. Scheduled for after queue
   completion.

## What deep cleanup would look like (when we get there)

This is what the cleanup buys, recorded so a future session has the
specification:

- **One full-pipeline test per `RunConfiguration`.** Runs the production
  canonical on 1P9J. Verifies every TR's H5 group is present, every key
  dataset populated, every physical magnitude in a sanity band. ~200
  lines, replaces hundreds of per-TR discipline tests.

- **One physics-correctness test per kernel where the math is non-trivial
  and could land wrong.** Biot-Savart with a known ring geometry.
  McConnell with the asymmetric three-term form (PATTERNS Lesson 19).
  PiQuad pure-T2 by Laplace. Cremer-Pople θ (already caught a bug there).
  These earn their keep — they catch sign/convention/formula errors that
  integration tests would let through. ~10-15 tests across the kernel
  layer.

- **Schema tests for the SDK/consumer interface.** The H5 contract between
  extractor and consumer is the load-bearing seam. One test that reads a
  freshly-extracted H5 with the SDK and verifies every documented field
  is readable and typed correctly.

- **Discipline tests stay for cross-cutting invariants.**
  Singleton-per-type on `TrajectoryProtein`. `Dependencies()` validation
  in Phase 4. The two or three things the architecture promises and the
  code has to deliver. One test for "AttachResult rejects a second
  instance of the same TR type" instead of one Frame0Semantics test
  per TR.

Total test surface target: 30-50 tests instead of the current ~200+.
The control machinery shrinks correspondingly. Skip flags live inside the
named factories, not in test fixture code. `Dependencies()` declarations
stay on TRs but Phase 4 becomes a single test target's domain. The per-TR
discipline trio dies entirely.

## Specific cruft inventory (for the deep-cleanup pass)

Things to scrutinize and likely delete:

- ~20 test files with their own `BuildBaseConfig` helper (each ~10-15
  lines of bespoke skip-flag combinations)
- The per-TR Frame0Semantics + FinalizeIdempotency + H5RoundTrip trio
  (~3 tests × ~30 TRs = ~90 tests, ~3500 lines)
- The `finalized_` boolean field in every Welford / FO TR (set but
  never checked; idempotency works via data-flow short-circuit
  already — pattern noted in `feedback_bounds_check_over_state_flag`)
- The synthetic `*FourFrames` tests for FO TRs where they verify only
  dense-buffer memory layout that the AV/FO exemplar already proves
- `ScanForDftPointSet` and its related TRs — already "slated for removal"
  per `OBJECT_MODEL.md`, just hasn't happened
- Any test that uses `skip_*` flags directly instead of selecting a
  test profile

Things to keep:

- Integration tests (`Integration1P9J` per TR) — these run real data and
  produce real numbers; they catch dropouts and magnitude regressions
- The conditional-attach gate machinery on the 12 Tripeptide/Larsen TRs
  — earned, source is genuinely DSN/grid-path-gated
- `Dependencies()` declarations on TRs — cheap, used by Phase 4
- The three named `RunConfiguration` factories — production interface

## The risk explicitly

What rolls forward in the next 24 hours without deep cleanup:

- ~12 more TRs landing with the current pattern (per-test config + four-
  test discipline ritual + naive source-field read)
- Physics bugs that pass `isfinite()` + `populated > 0` + sanity-bound
  assertions can land undetected if the science-focused review doesn't
  happen
- The HBondCount silent-zero failure mode I caught during landing is a
  canary for the class. The Welford accumulator has no way to distinguish
  "source did not fire" from "source produced a real value of zero."
  Eeq/Sasa/HM/McConnell could all hit this if a future test or fleet run
  has a similar misconfiguration; production currently does not.
- More test debt accreting on top of debt that needs sweeping.

Jessica's framing, recorded verbatim because the meta-observation matters:

> "Good software proceeds from looking at what is true at each stage and
> addressing it even if it hurts a bit."

The bet says: it hurts more right now to take the week than to take the
gamble. The bet's payoff depends on the science-focused review actually
catching anything that the framework-focused review misses. If it does,
the bet is fine. If it doesn't, we have a rollback tag.

> "I want my cookie, you want me to have my cookie, and we are unreliable
> narrators."

Both sides of the unreliability are recorded here. The plan is auditable
later.

## Open: water TR shape

The next TR batch is the Water + Hydration trio
(`WaterEnvironmentTimeSeries`, `HydrationShellTimeSeries`,
`HydrationGeometryTimeSeries`). Jessica flagged: "we need to talk about
water though, it may be an npy too."

Open questions for that batch (to be discussed before landing):

- What's the natural shape of water-trajectory data? Per-atom + per-frame
  + per-water-molecule produces a three-axis array that doesn't fit the
  per-atom × per-frame DenseBuffer pattern cleanly.
- Should it land as `WriteFeatures` NPY emission rather than `WriteH5Group`?
- Is the right granularity per-atom (atom sees these waters in this
  frame) or per-water-molecule (molecule has this position/orientation
  in this frame)?
- Does Water want to be a TR at all, or is it a separate runtime artifact
  the extractor produces alongside the trajectory H5?

Not resolved here. Flagged so the conversation happens before the code.

## Conversation transcript (high-level)

1. **Adversarial review request** of today's two TR commits. Jessica asked
   explicitly that the reviewer read OBJECT_MODEL.md and PATTERNS.md so
   the findings would be on-target.
2. **Both parallel reviewers** converged on HIGH "add HasResult<> gates"
   finding. Real bug surface, real-but-rejected mitigation.
3. **My summary** flagged the tension with the two-use-cases rule and
   listed lower-severity catalog-drift / schema-drift findings.
4. **Jessica:** "fine we will live with the dependency path machinery.
   It is in no sense an unalloyed good. You took it way too far."
5. **Jessica's six-month observation** about going along with my
   complexity creep and regretting it.
6. **My narrower proposal**: TestConfig profile + migrate today's 5
   tests + adopt for new work + lazy-migrate old tests.
7. **Jessica's escalation**: hundreds of tests, half cruft, the
   many-tiny-tests pattern outlived its usefulness for this codebase.
8. **My acknowledgment** of the deeper architectural critique + sketch
   of use-case-anchored test shape (30-50 tests target).
9. **Jessica's decision**: take the bet — tag, short fix, queue,
   deep cleanup deferred. With the explicit mitigation that a
   science-focused adversarial reviewer is added going forward.

Read in order, the conversation shows what the bet actually rests on.
Documented here so future sessions know why the test debt was tolerated
during the queue completion phase and what the cleanup target shape is.

## On rollback

If the bet fails — physics bugs land in the TR queue and aren't caught
by the science-focused review — the recovery is:

1. `git checkout pre-test-config-consolidation-2026-05-17`
2. Re-evaluate the queue with the deep-cleanup-first ordering
3. The TR queue would re-land into the cleaned-up test suite, which
   would have caught the bugs the first time

This is the kind of rollback that costs the queue work but not the
cleanup work. The work that goes into the deep cleanup is preserved
across the rollback because the cleanup happens upstream of the queue
in that scenario. The tag exists so that this rollback is actually
possible.
