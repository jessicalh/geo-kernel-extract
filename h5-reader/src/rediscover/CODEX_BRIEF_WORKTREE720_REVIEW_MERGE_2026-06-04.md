# Codex mission — adversarially review the worktree-720 ingest + plot a safe merge to one tree

Working root: `/shared/2026Thesis/nmr-720` (the worktree, branch `wt-720-build` — the code under review lives
here, uncommitted). The consolidation target is the main tree `/shared/2026Thesis/nmr-shielding` (branch
`h5-reader-pysr-spike`). The lead owns ALL git; you REVIEW + PLAN, you do not execute git and you do not edit
the code under review.

Honest stakes: the static-pose ingest in this worktree emits the per-(protein,atom) substrate that becomes
the combined 1P9J+720 **Model-1 training corpus**. A bug here is silent — wrong DFT target, a T2
orientation/sign error, a field-width or atom-alignment mismatch, a unit slip — and it poisons the entire
Stage-3 model without ever crashing. The maths-audit went three rounds before the analysis numbers were
trusted; this code deserves the same ruthlessness. You may have written it — now assume it's wrong and find
where. Do not manufacture findings; if a path is sound, say so — but probe the load-bearing ones hard first.

## PART 1 — adversarial review of the worktree-720 static-ingest

**The code (the 12 uncommitted modified files here):**
`h5-reader/src/io/{FrameNpyLoader,QtProteinLoader}.{cpp,h}`, `h5-reader/src/model/QtProtein.h`,
`h5-reader/src/rediscover/{Catalog,PerAtomSubstrate,RunData,main_extract}.{cpp,h}`. (`git -C
/shared/2026Thesis/nmr-720 diff` shows the deltas vs the `04fcb3f` base.)

**Review against:** `h5-reader/src/rediscover/SPEC_720_STATIC_INGEST_2026-06-04.md` (the contract) + the
disciplines (model-is-spine / no Python second model / no-file-discovery / lean-disk <15 G / **T2 sacred —
never scalar-collapse** / no-symlinks / extractor SACRED).

**Hunt, specifically (load-bearing paths first):**
- **DFT target build** — raw absolute dia/para/total; is `total == dia + para` enforced componentwise? atom
  count aligned to topology? Is the target the *absolute* WT shielding, never a mutation-delta?
- **Tensor / T2 handling** — the spherical decomposition + T2-orientation policy. Any scalar-collapse? Any
  sign/convention/basis error (the kind the maths-audit caught — Cremer-Pople θ inversion, dipole-coherence
  dual-bug are the precedents)? Is `frame_valid` honored, not silently defaulted?
- **Field/catalog loading** — strict exact-path (no glob/try-and-fail)? Any hard-coded field-width assumption
  that breaks on the real producer output? `RunData::staticPoseExpectedFields()` correct + honest about
  required vs optional?
- **All-atom emitter** — every row carries BMRB/IUPAC **and** element/role/stratum (never sliced to one)?
  No backbone-only / toy reduction?
- **Memory backends** — both strategies correct + non-foreclosing; reductions C++-side, nothing leaking a
  Python second model or spatial work.
- **The oracle-parity gate** — **re-run it yourself** (build `build/linux-rwdi/h5reader_extract` if needed —
  but do NOT modify `src/`). Confirm it PASSES *and* that it is a meaningful equivalence test, not a
  tautology that would pass even if the static path were wrong.
- **Known flag:** the `nanoflann.hpp` symlink (no-symlinks violation) — confirm + state the real-vendoring fix.

**Output:** a severity-ranked findings list (CRITICAL / HIGH / MED / LOW), each with `file:line`, what's
wrong, why it matters, the fix. Then a one-line **verdict**: is this code safe to merge as-is, or what must
fix first.

## PART 2 — git/worktree topology inventory + safe-merge-to-one-tree plan

The lead wants all scattered real work consolidated safely back into one tree (`h5-reader-pysr-spike`),
nothing lost.

- **Inventory (read-only git):** `git worktree list`, `git branch -vv`, and per worktree/branch with real
  work, committed-vs-uncommitted (`status`) + what it touches. Known: `wt-720-build` (this build,
  uncommitted, base `04fcb3f` which is an ancestor of `h5-reader-pysr-spike`@`dfc8f51`); the
  `worktree-agent-*` branches (viewer-MOPAC `d036bea`, Tripeptide trajectory-results, Welford work — committed
  on their own branches); `master`.
- **Classify** each: real work to consolidate / superseded / abandoned-or-prunable / unclear (flag the
  unclear ones for the lead — don't guess).
- **Plot a SAFE consolidation sequence** into `h5-reader-pysr-spike`, preserving everything: the exact ordered
  git steps (e.g. commit the `wt-720-build` uncommitted changes first — they're disjoint code vs the docs-only
  `dfc8f51`, so the merge should be conflict-free; then the agent-worktree merges/cherry-picks; the nanoflann
  real-vendoring), each step's conflict expectation + why it's safe (no loss, no force, reversible). PROPOSE
  only — the lead executes. Flag every step needing her decision (which agent-worktrees to keep).

## Hard boundaries
- **Read-only on git** — inventory + a written plan; NO `commit/add/merge/branch/switch/rebase/reset/push`.
- **Do not edit the code under review** — findings + proposed fixes, not edits. **Do not modify `src/`**
  (SACRED extractor); you may build/run to re-test the ingest.
- **Do NOT investigate or touch the IUPAC-topology episode** — no spelunking deleted branches/history/archives
  for it; if a branch or history brushes it, write "out of scope" and move on. No destructive proposals (no
  branch deletion without flagging for the lead; never force-push).
- **No sprawl.** Truthful, cited (`file:line`), no manufactured findings.

## Deliverable
Write your report to `/shared/2026Thesis/nmr-shielding/h5-reader/src/rediscover/REVIEW_worktree720_AND_MERGE_PLAN_2026-06-04.md`
(Part-1 findings + verdict; Part-2 inventory + classification + the ordered safe-merge plan + open decisions
for the lead), and summarize it in your final output.
