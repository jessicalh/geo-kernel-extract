# Implementation brief (per algorithm)

You are applying the **approved fix plan** for ONE algorithm to the actual
`src/` source. You ARE modifying source now (this is not a planning pass).

## Inputs
- the algorithm's `plan.md` (the edit list — §3 "Edits that don't move
  numbers", plus any blessed bug fix in §5).
- `reviews/applied-math/OPEN_QUESTIONS.md` — **apply every RESOLVED ruling that
  touches this algorithm**, both the family decisions and the per-plan rulings.
  (A ruling there overrides the plan if they differ.)
- the source file(s).

## What to apply
- The plan's comment fixes, signposts, named intermediates, regroupings, and
  **file-local / sanctioned renames**.
- The **blessed** edits from OPEN_QUESTIONS (e.g. eeq `q_min/q_max` seed fix;
  `N==0` guard in apbs/eeq/sasa; dead-field deletions `RingAccumulated` /
  `HBondKernelResult::f`* / `FieldValue`; BiotSavart small-rho threshold unify;
  saturated non-5-ring log line; BackboneHA `&& element==Element::H` tighten).
  *(\*hbond `f`: do NOT delete — it is being wired into an emitted scalar; if
  this is the hbond pass, see the plan; otherwise leave it.)*

## Hard constraints
- **Numbers do not move** except the explicitly blessed bug fixes above. If an
  edit would change a computed value and isn't on the blessed list, do not make
  it — note it instead.
- **Output / serialized names stay**: NPY field names, `_catalog.py` keys, H5
  dataset names. Internal names may improve per the plan/rulings.
- **Already done repo-wide — do NOT redo**: `SampleShieldingAt → SampleKernelAt`
  and `JBLobeOffset → JohnsonBoveyLobeOffset` are already renamed everywhere. If
  the plan mentions them, treat as complete.
- **Comments conform to code.** No algorithm changes, no refactors beyond what
  the plan specifies.
- **Clarity bar**: clear, unambiguous purpose per step, but *not* broken up in a
  condescending way and not over-commented — conservative good practice.
- **Do NOT build or run tests.** The build happens once, later, separately.

## Output
Apply the edits with the Edit tool to the named `src/` files. Then reply with a
concise changelog: `file:line — what changed`, grouped (comments / signposts /
renames / blessed-fixes). List any plan item you could NOT cleanly apply and
why. Do not modify files outside this algorithm's source.
