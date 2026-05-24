# Fix-planning brief (per algorithm)

You are planning readability fixes for ONE algorithm in the `nmr_shielding`
C++ library. **You are not implementing anything.** You produce a `plan.md`.

Inputs, all in this algorithm's directory and the repo:
- `codex.md`, `claude.md` (and `codex-correctness.md` if present) — the two
  independent readability reviews. Treat them as suggestions, not orders: some
  conflict, some are wrong, some are unsafe.
- the source file(s) under review.
- the rest of the repo, for tracing how quantities are used.

## The clarity bar

External review levels of clarity will be expected in the code. Both students
and experienced readers will be examining the code and hoping to understand the
work. The code should be clear — **not broken up in a condescending way**, but
**without ambiguity of purpose in each step**, to the degree possible and
appropriate for **conservative good practice**.

This cuts both ways. Improve a name or add a signpost where the purpose of a
step is genuinely unclear — but do not fragment coherent code or pile on
hand-holding comments; over-explaining is its own unclarity. The aim is
unambiguous purpose per step, reached with restraint.

## Governing prior — read this first, it controls everything

The algorithm, its **signs**, and its **numbers** are as they are because of
**repeated iteration**. Assume a reason exists — embedded in how the whole
ensemble uses the quantity — and assume it until proven otherwise **by
exhaustion**. Your job on any sign/value question is to *discover the reason*,
not to judge the code against a comment or your own expectation.

- "I don't see why this sign is here" is **not a finding**. Only "I traced
  every consumer and the convention cannot be reconciled anywhere" is.
- Under-tracing is the failure mode, not under-fixing. The expected verdict is
  "coherent — here is the reason," not "bug."

## What stays fixed

- **The algorithm and its numbers stay exactly as they are.** The only change
  that may move a number is fixing an **outright bug**, and only after the
  exhaustion above.
- **Output names are the protected ones — these are what "ideally don't
  change" means.** Serialized / contract names — NPY field names, H5 dataset
  names, SDK `python/nmr_extract/_catalog.py` keys — stay stable. Propose
  renaming one only with a strong reason *and* full consumer carry-through.
- **Internal names should be improved where it helps — that is a goal, not
  churn to avoid.** Local variables, file-local statics, internal struct
  fields, and internal methods/functions are not a contract. If a clearer name
  aids readability (the `H`/`V`/`G` symbol-soup kind of finding), **propose
  it.** Do not decline an internal rename merely because it's a change or to
  keep a sibling file's matching opacity. For an internal name consumed across
  files or subprojects, still propose the improvement but note the
  carry-through cost so the human can weigh it — declining for cost is a
  weighed call, not a default.
- **Comments conform to code** unless exhaustion proves the code wrong.
- **No algorithm changes, no refactors, no new abstractions, no library swaps.**

## Tracing — how quantities are used

For every finding about a **sign** or a **computed value** (and for any rename
that might cross a file boundary), trace how the quantity is produced here and
consumed across the ensemble:
- `src/` callers,
- `python/nmr_extract/_catalog.py` (the NPY/SDK contract),
- `learn/` (calibration),
- `h5-reader/`,
- `ui/`,
- and the physics it feeds.

**A producer/consumer discrepancy does not implicate the producer.** The bug
may be in how a consumer *reads* the quantity, or there may be no bug (a
convention the consumer compensates for elsewhere). Either side can be wrong,
or neither. Understand both sides before concluding anything; if you can't, it
is a question.

Land each sign/value finding on one of:
- **coherent (expected)** → you found the reason; fix only the comment so it
  states that reason plainly.
- **bug, by exhaustion** → you traced every consumer and it cannot be
  reconciled. State the exhaustion you did and the minimal fix. Note whether
  the bug is on the producer or a consumer side. Rare.
- **can't tell** → a question, not a guess.

## Output — write `plan.md` to the path you are given

Sections:
1. **Summary** — does this file mostly tell a coherent story; what the fix
   pass will and won't touch.
2. **Review-finding ledger** — list **every** finding from both reviews
   (`codex.md`, `claude.md`, and `codex-correctness.md` if present), one row
   each, with an explicit disposition: **adopted** (→ which edit in §3),
   **declined** (→ one-line reason), or **duplicate** (→ of which finding).
   Nothing from the reviews may be silently dropped. A reader must be able to
   see *how each input finding was weighted on the way out* — that visibility
   is the point of this section, not a formality.
3. **Edits that don't move numbers** — the concrete list: comment fixes,
   2–4-word signposts, named intermediates, regrouping, any justified rename
   (with its consumer carry-through). Each as `file:line — change`.
4. **Usage notes** — for each sign/value the reviews raised: the *reason* you
   discovered, where the quantity is consumed, and whether producer and
   consumers agree. This is the real product of the plan.
5. **Bug-by-exhaustion candidates** — if any; with the exhaustion shown.
6. **Questions & Ambiguities** — conventions you couldn't confirm, usages you
   couldn't fully trace, review disagreements you couldn't settle. Ask rather
   than guess.

Write ONLY `plan.md`. Modify no source files.
