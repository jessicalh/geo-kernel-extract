# Codex brief — doc truth-pass (grinder): make every rediscover doc TRUE, cut cruft, concise

> **Historical brief — not current truth (trued 2026-06-04).** Previous
> doc-truth pass provenance only; this decruft pass supersedes it.

You own the grind. Walk every rediscover doc and make it **true** against reality, **cut the
cruft**, and **correct tersely**. This is an EDIT+REPORT pass; the lead owns all git operations.
Runs AFTER Build 3 has landed (so you fold its commit/result/postmortem). Branch `h5-reader-pysr-spike` — NEVER merge/switch/rebase/PR.

## NO GIT — the lead owns all of it

**Do NOT run any git command. No commit, add, reset, rebase, checkout, stash — nothing.**
The lead has already checked the docs into the branch as a restore point before you start,
and the lead commits your corrections after reviewing them. Your entire job is: EDIT files
and report. Do not touch git state in any way. (Given git tools an agent tends to "tidy" the
tree — here that risks the never-merge branch and a dirty working tree. So: hands off git.)

## The standard (in priority order)

1. **TRUTH.** Every factual claim must match reality — `git log --oneline` on the branch
   (the actual commits + HEAD), the emitted run dirs, and the code. Where a doc asserts
   something the commits/code/runs no longer match, **correct it to the true statement.**
2. **CUT CRUFT — do not preserve it with a label.** Remove stale/superseded/redundant
   content outright. Specifically:
   - `STATE.md` self-contradictions (stale block + corrected block side by side): keep the
     TRUE current one, **delete the superseded one.** Known cases: the old "MISSING
     EXTRACTIONS" list (those channels are now emitted), "Loop 4 next" (it's Build 3), and
     "ringχ in the slope≈1 validation" (corrected to opposite-convention, never naïve-slope).
   - Wholly-superseded session-start snapshots `FRESH_LOOK_2026-06-03.md` and
     `H5READER_COMPREHENSION_2026-06-03.md` (still say HEAD `00ec168` / "review the spec
     next"): replace each with a one-line pointer to `NOW.md`/`STATE.md`, or delete — don't
     leave the stale body.
3. **CONCISE.** A correction adds **no more than a sentence.** The pass should make the docs
   SHORTER and truer, not longer. No verbose change-logs, no "previously we thought…" padding.
4. **PRESERVE the darlings — do NOT cut these.** Anything a deep-dive reader could NOT
   reconstruct from the code: the decisions, the rationale/why, the disciplines, the
   physics judgment, the through-lines (laws-from-dominance-clean-exemplars; per-type-not-
   pooled; dominance-as-gate; dia due-diligence; categorical-engine parked; the gating
   cadence). These are not cruft; keep them, corrected to truth if needed.
5. **When unsure whether something is cruft or a darling: KEEP it and note it** in your
   summary — do not cut a decision/rationale on a guess.

## Specific work

- **Fold Build 3:** reflect its commit + `POSTMORTEM_BUILD3.md` + result in `NOW.md` (the
  live marker — keep it short, update the loop ledger + current step) and `STATE.md`.
- **Fill the design-doc gap:** the partition-filter architecture currently lives only in
  briefs/STATE. Write a durable `src/rediscover/PARTITION_FILTER_DESIGN.md` distilling it
  (the C++ filters over the in-memory indexes, isolation primitives, the CaseHunter, the
  C++/Python boundary, the law-vs-model layer split, dominance-as-gate). Concise — a design
  doc, not a transcript.
- **Truth-correct the rest:** `GUIDANCE.md`, `DESIGN.md`, `SURFACE_DESIGN.md`, the specs
  (`PER_ATOM_SUBSTRATE_SPEC`, `ALLATOM_FIT_SPEC` — un-mark "draft for review only" since it
  landed), the briefs/postmortems (replace placeholder "exact SHA printed by codex" with the
  real commit hashes from `git log`). One true sentence each; cut what's stale.
- Leave the early-piece postmortems alone (covered by STATE + commit messages — not worth
  back-filling).

## Output — edit + report only, NO commit

Edit the docs + write the new `PARTITION_FILTER_DESIGN.md` in place (the working tree; the
lead commits). Write `src/rediscover/POSTMORTEM_DOCTRUTH.md` (≤30 lines): a **diff-summary**
— what you CORRECTED (claim → truth, with the reality you checked it against), what you CUT
(and why it was cruft), what you ADDED, and anything you KEPT-because-unsure (for the lead to
eyeball before committing). Print only a ≤10-line summary + that path; DO NOT echo full diffs
to stdout. Doc edits only — no emit, no `/shared`, and (restating) **no git.**
