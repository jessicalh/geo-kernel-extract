# Fix plan — ring-svd-plane-fit (`src/Ring.{h,cpp}`)

## 1. Summary

Both files tell a coherent story. `Ring.h` is a flat, table-like type
hierarchy (base identity → geometry/accumulated structs → concrete ring
classes grouped by size category); `Ring.cpp` is a textbook
`ComputeGeometry` (collect vertices → centroid → SVD plane normal →
orientation guard → mean radius) plus a fail-loud `CreateRing` factory.
A chemist or external reviewer can read both top-to-bottom.

The friction the reviews found is real but small and almost entirely
**comment/signpost** work, plus a few **internal-name** improvements:

- One genuinely stale comment: `Ring.h:9` says "8 types" but the file
  defines 9 (it omits `ProPyrrolidineRing`). The sibling enum in
  `Types.h:180` already says "9 ring types" — Ring.h is the outlier.
- One stale provenance comment: `Ring.h:25` attributes geometry to
  `GeometryResult`; the method is `Ring::ComputeGeometry`.
- A dead struct (`RingAccumulated` + the `accumulated` member) with
  **zero consumers** anywhere in the tree — a real readability tax
  because it invents a vocabulary (`total_B_at_center`,
  `total_G_T0_diagnostic`, `mutual_B_from`) that the reviews then
  (correctly) flagged as opaque. Recommended for **removal**, which
  dissolves three of the naming findings at once.
- Two verbose prose blocks (`Ring.cpp:65-70` fail-loud rationale;
  `Ring.h:178-185` ProPyrrolidine docblock) — trim/keep per below.
- Internal names worth improving: `RingSizeValue → RingAtomCount`,
  `JBLobeOffset → JohnsonBoveyLobeOffset` (both cross-subproject — cost
  noted), `edge01/edge02 → first_edge/second_edge` (file-local, free).
  `TypeName → TypeLabel` proposed but **declined on weighed cost**.

**Will touch:** comments, one stale count, the dead struct, two
file-local variable names, and (proposed, cost-noted) two
cross-subproject method renames. **Will not touch:** the algorithm,
any number, the SVD column index, the literature intensity values, the
`RingTypeIndex` enum, serialized field names (`rings.npy`:
`atom_count` / `ring_type_index`).

---

## 2. Review-finding ledger

Every finding from `codex.md` and `claude.md`, one row each. (No
`codex-correctness.md` present in this directory.)

### codex.md

| # | Finding | Disposition |
|---|---------|-------------|
| C1 | `Ring.h:9` "8 types" stale (9 defined incl. ProPyrrolidine) | **adopted** → §3 E1 |
| C2 | `Ring.h:25` "computed by GeometryResult" conflicts with `Ring::ComputeGeometry` | **adopted** → §3 E2 |
| C3 | `Ring.h:39-44` `RingAccumulated` is an unexplained bag of fields; label or rename | **adopted (stronger)** → §3 E3 (remove the dead struct entirely; no consumers) |
| C4 | `Ring.h:178-185` ProPyrrolidine docblock too long; compress to signpost + citation | **declined** → see §3 E8 reasoning (claude C-equiv judged it "earns it"; it justifies a literal `0.0` as physics) |
| C5 | `Ring.cpp:65-70` fatal-path comment is process narration; trim to `// invalid enum guard` | **adopted (middle ground)** → §3 E7 |
| C6 | `Ring.h:40` `total_B_at_center` lacks units/frame; `total_field_at_center` | **adopted-by-removal** → dissolved by E3 (dead struct removed) |
| C7 | `Ring.h:42` `total_G_T0_diagnostic` is symbol soup | **adopted-by-removal** → dissolved by E3 |
| C8 | `Ring.h:43` `mutual_B_from` omits key/value meaning | **adopted-by-removal** → dissolved by E3 |
| C9 | `Ring.h:64` `JBLobeOffset()` — `JB` unexplained | **adopted** → §3 E5 (rename, cost noted) + E6 (signpost regardless) |
| C10 | `Ring.h:68` `TypeName()` returns labels like "TRP6", not full names; `TypeLabel()`/`RingLabel()` | **declined** → cross-subproject + Qt string-policy contract; see §4 |
| C11 | `Ring.h:67` `RingSizeValue()` awkward; `RingAtomCount()` carries meaning | **adopted** → §3 E4 (rename, cost noted) |
| C12 | `Ring.h:62-64` `Intensity()`/`LiteratureIntensity()`/`JBLobeOffset()` lack units/convention; add comment block | **adopted (partial)** → §3 E6 (units signpost on virtual block) |
| C13 | `Ring.cpp:34-39` `edge01`/`edge02` force light decoding; `first_edge`/`second_edge` | **adopted** → §3 E9 |
| C14 | `Ring.h:61` "Virtual const properties" abstract; `// ring-type properties` | **adopted** → §3 E10 |
| C15 | `Ring.h:70` "Non-virtual queries" clear enough | **declined** → no change needed (reviewer's own verdict) |
| C16 | `Ring.cpp:23-25` SVD signpost good | **declined** → keep as-is (already good) |
| C17 | `Ring.h:9-12` stale inventory is a doc-correctness issue | **duplicate** of C1 |
| C18 | `Ring.cpp` `positions[idx]` bounds depend on caller contract; not flagged as bug | **declined** → no change (matches project fail-loud-on-programmer-error stance; same as claude Cc4) |

### claude.md

| # | Finding | Disposition |
|---|---------|-------------|
| L1 | `Ring.h:64` `JBLobeOffset()` — "JB" unexpanded, add `// Johnson-Bovey lobe z-offset, Å` | **duplicate** of C9 → §3 E5/E6 |
| L2 | `Ring.h:42` `total_G_T0_diagnostic` — `G`/`T0` unglossed; trailing note | **duplicate** of C7 → dissolved by E3 |
| L3 | `Ring.h:43` `mutual_B_from` — key identity unstated; `// keyed by source ring index` | **duplicate** of C8 → dissolved by E3 |
| L4 | `Ring.cpp:65-70` fail-loud is 6 lines of prose; trim to one-line signpost | **duplicate** of C5 → §3 E7 |
| L5 | `Ring.h:178-185` ProPyrrolidine block long but earns it; optionally shed citation sentence | **declined** → keep (see E8); explicitly "acceptable as-is" per reviewer |
| L6 | `Ring.cpp:23-25` SVD-normal comment correct/well-placed | **duplicate** of C16 → no change |
| L7 | `Ring.cpp:34-35` right-hand-rule note is the right signpost | **declined** → keep as-is |
| L8 | `Ring.h:75 / Ring.cpp:8` `ComputeGeometry` name accurate | **declined** → no change (clean) |
| L9 | `Ring.h:63,95` `LiteratureIntensity()` vs `Intensity()` distinction clear | **declined** → no change (clean); units note still added via E6 |
| L10 | `Ring.cpp:32` `svd.matrixV().col(2)` hardcodes col 2; add `// col 2 = smallest σ (V is 3×3)` | **adopted** → §3 E11 |
| L11 | `Ring.cpp:38` orientation guard undefined if first 3 atoms collinear; check callers never pass degenerate set | **declined as code change; → §6 Question** (real aromatic/saturated rings can't be collinear; documented as a question to confirm caller contract, no guard added) |
| L12 | `Ring.cpp:21,46` divide by `vertices.size()` no zero guard; covered by `empty()` early return | **declined** → no change (reviewer confirms it's safe) |
| L13 | `Ring.cpp:15` `positions[idx]` unchecked; consistent with fail-loud stance | **duplicate** of C18 → no change |

---

## 3. Edits that don't move numbers

**E1** — `Ring.h:9` — change `// 8 types in 3 size categories:` to
`// 9 types in 3 size categories:` and add `ProPyrrolidineRing` to the
`FiveMemberedRing:` line of the inventory (currently it lists only the
four π-systems; add the saturated Pro 5-ring). Brings the docblock in
line with `Types.h:180` ("9 ring types") and the actual class list.

**E2** — `Ring.h:25` — change
`// Ring::Geometry -- conformation-dependent, computed by GeometryResult`
to `// RingGeometry -- conformation-dependent, computed by Ring::ComputeGeometry`.
(Both halves stale: the struct is `RingGeometry`, the method is
`Ring::ComputeGeometry`. `GeometryResult` is not in this pipeline — no
such symbol is referenced by Ring.{h,cpp}.)

**E3** — `Ring.h:35-44` and `Ring.h:78` — **remove** the
`RingAccumulated` struct and the `RingAccumulated accumulated;` member
on `Ring`. Exhaustive trace (grep across `src/`, `ui/`, `h5-reader/`,
`python/`, `learn/`) finds **no reader or writer** of `.accumulated`,
`total_B_at_center`, `intensity_used`, `total_G_T0_diagnostic`, or
`mutual_B_from` anywhere in the tree. It is dead scaffolding; removing
it deletes the three "symbol soup" naming findings (C6/C7/C8, L2/L3)
rather than papering over them with comments. *Carry-through: none —
no consumer.* (If the human prefers to keep it pending a planned
post-pass, the fallback is to retitle the fields with the trailing
notes the reviews suggested; but the no-consumer trace argues for
removal.)

**E4** — `Ring.h:67` (+ overrides at 88, 131, 212) — rename method
`RingSizeValue()` → `RingAtomCount()`. Returns the atom count of the
ring (6 / 5 / 9). *Carry-through (cross-subproject):* the parallel
`h5-reader/src/model/QtRing.h` mirrors this method at lines
83/99/104/204 — those four lines must rename in lockstep to preserve
the deliberate library↔Qt mirror (per CLAUDE.md, the two hierarchies
mirror by design). **No serialized impact:** `rings.npy` stores
`atom_count` computed from `atom_indices.size()` in
`TopologySidecar.cpp:217`, *not* from `RingSizeValue()` — the method
has no NPY/H5 consumer. Total touch: 4 lines in Ring.h + 4 in QtRing.h.

**E5** — `Ring.h:64` (+ overrides at 96,107,118,139,150,161,172,191,209)
— rename method `JBLobeOffset()` → `JohnsonBoveyLobeOffset()`. "JB" is
Johnson-Bovey, confirmed by the consumer: `BiotSavartResult.cpp` feeds
it into `JohnsonBoveyField(...)` (lines 188, 222, 380, 406).
*Carry-through (cross-subproject, larger):* consumed in
`src/BiotSavartResult.cpp` (5 call sites: 181, 186, 223, 380, 406) and
in `h5-reader/` — `QtRing.h` (10 override sites) plus two overlay call
sites (`QtBFieldStreamOverlay.cpp:194`, `QtFieldGridOverlay.cpp:183`).
~25 lines across three subprojects. **Weighed recommendation: do the
rename** — the expansion is mechanical and the clarity gain is real
for external readers — but flag the cost so the human can choose E5
(rename) or E6-only (signpost). If declined for cost, E6 alone still
resolves the finding.

**E6** — `Ring.h:61-68` — above the virtual property block, replace the
abstract `// Virtual const properties (overridden by each type class)`
with a units/convention signpost, e.g.:
```
// Ring-type physical properties (one override per type class):
//   Intensity / LiteratureIntensity : ring-current strength, nA/T
//                                      (negative = diamagnetic)
//   JohnsonBoveyLobeOffset           : JB two-loop z-offset, Å
//   NitrogenCount                    : ring N atoms
```
Covers C12/C14 and the units half of C9/L1/L9 in one place. (Use
`JBLobeOffset` here if E5 is declined.)

**E7** — `Ring.cpp:65-70` — trim the 6-line fail-loud rationale to a
2-line signpost that keeps the *why* without the project-process
narration, e.g.:
```
// Fail loud: an out-of-range enum cast is a programmer error.
// (Silent fallthrough to PheBenzene would forge a chemically real PHE.)
```
Middle ground between codex's `// invalid enum guard` (loses the
forged-PHE rationale, which is the actual reason) and claude's
single-line version. Keep `abort()`.

**E8** — `Ring.h:178-185` — **keep the ProPyrrolidine docblock.** Both
reviews split (codex C4 wants it trimmed; claude L5 says it "earns it,
acceptable as-is"). It justifies the literal `0.0` returns as *physics,
not calibration* — exactly the kind of non-obvious decision a docblock
should carry, and the only place that fact is stated. No edit.

**E9** — `Ring.cpp:36-37` — rename file-local `edge01` → `first_edge`,
`edge02` → `second_edge` (and the two uses at line 38). File-local,
zero carry-through, pure clarity win.

**E10** — folded into **E6** (the "Virtual const properties" header is
replaced there).

**E11** — `Ring.cpp:32` — append documentation of the hardcoded column:
`geo.normal = svd.matrixV().col(2);  // col 2 = smallest singular value (coords is N×3 ⇒ V is 3×3)`.
Documents the assumption claude L10 raised without changing behavior.

---

## 4. Usage notes (sign/value/name questions raised by the reviews)

**`JBLobeOffset` / Johnson-Bovey lobe offset (C9, L1).**
*Reason found:* it is the z-offset of the two current loops in the
Johnson-Bovey two-loop ring-current model. *Consumed by:*
`src/BiotSavartResult.cpp` passes it as the third argument to
`JohnsonBoveyField(vertices, normal, lobeOffset, current, point)` (lines
188/222/380/406) and records it as `"lobe_offset"` in Å (line 186).
`h5-reader/` recomputes the JB field locally in two overlays using the
same offset. *Producer/consumer agree.* The library reads the value
from `CalculatorConfig` (calibrated); `QtRing.h` hard-codes literals
(0.64/0.52/0.50/0.60/0.0) — see §6 for the value-drift note.

**Literature intensity sign convention (claude verdict, C12).**
*Reason found:* `LiteratureIntensity()` is negative for diamagnetic
rings (-12.0 PHE, -19.2 indole, etc.). The minus is consumed correctly:
`BiotSavartResult.cpp:225-228` documents `G_ab = -n_b·B_a·PPM_FACTOR`
with the comment *"sigma = I·G gives the correct sign using literature I
values (I < 0 for diamagnetic rings)."* *Producer and consumer agree;*
the sign is intentional and reconciled at the kernel. Coherent — the E6
signpost just states it plainly. **No number changes.**

**`RingSizeValue` → atom count (C11).** *Reason found:* returns ring
atom count (6/5/9). *Consumed by:* nothing serialized — `rings.npy`'s
`atom_count` derives from `atom_indices.size()`
(`TopologySidecar.cpp:217`), not this method. The method's only mirror
is `QtRing.h`. Safe to rename to `RingAtomCount` (E4); producer/consumer
(library/Qt mirror) rename together.

**`TypeName` (C10) — why declined.** *Reason found:* returns short
residue/ring labels ("PHE", "HIE", "TRP6", "TRP9", "PRO"). *Consumed
by:* `ui/` (MainWindow, ComputeWorker, RestServer — 6 sites, several
emitted into REST JSON as `ring_type`), `h5-reader/QtRing.h` (10
override sites) plus a documented **string-policy contract** in
`QtRing.h:14-17` ("the only string surface is `TypeName()` … matches the
library's `Ring::TypeName()`"). A rename would break that explicit
named contract and touch three subprojects including a REST surface. The
clarity gain ("Label" vs "Name" for a string-returning method) is
marginal; the carry-through is large and crosses a documented contract.
**Weighed decline.** (`TypeLabel` is also arguably *less* accurate than
"Name" — these are canonical short names, used as names elsewhere via
`RingTypeName` in `Types.h:243`.)

**"8 types" (C1).** Stale doc, not a value. Ring.h:9 predates the
ProPyrrolidine addition; `Types.h:180` was already updated to "9". E1
makes Ring.h conform to the code. No algorithm/number impact.

**`GeometryResult` (C2).** No such symbol exists in the Ring pipeline;
the producer is `Ring::ComputeGeometry`. Stale comment, E2 fixes it.

**`RingAccumulated` (C3/C6/C7/C8, L2/L3).** Traced every field across
all subprojects: no reader, no writer. Dead. E3 removes it.

---

## 5. Bug-by-exhaustion candidates

**None.** No reviewer raised an algorithm/sign/value bug, and tracing
confirms none. The SVD column index, orientation guard, centroid/radius
divisions, and the literature-intensity sign are all coherent with
their consumers (see §4). `Ring.cpp` `positions[idx]` and
`vertices.size()` divisions are guarded by the caller contract and the
`atom_indices.empty()` early return respectively — both reviews agree
they are safe, not bugs.

---

## 6. Questions & Ambiguities

1. **Degenerate-triangle orientation guard (claude L11).** `Ring.cpp:38`
   uses `edge01.cross(edge02)` (vertices 0,1,2) to fix normal sign. If
   the first three ring atoms were collinear the cross product → 0 and
   the sign flip is undefined. For real aromatic/saturated rings this
   cannot occur, so **no code guard is proposed.** Question for the
   human: is the caller contract for `atom_indices` (ordering /
   non-collinearity of the first three) documented anywhere, or should a
   one-line comment at `Ring.cpp:36` state "first three ring atoms are
   never collinear for real rings"? Flagging rather than guarding.

2. **`RingAccumulated` removal vs retention (E3).** It has zero
   consumers today. Was it scaffolding for a *planned* post-pass
   (mutual ring-field accumulation) that hasn't landed? If a future
   slice intends to populate `mutual_B_from` / `total_B_at_center`,
   the human may prefer to keep it with the suggested clarifying field
   comments instead of removing it. Default recommendation is removal;
   confirm there is no pending consumer.
   **RESOLVED 2026-05-24 (Jessica): dead to us right now — remove it.**

3. **`QtRing.h` lobe-offset value drift (out of scope, noted).** The
   library reads JB lobe offsets from `CalculatorConfig`
   (`phe_benzene_jb_lobe_offset`, etc.), while `h5-reader/QtRing.h`
   hard-codes literals (0.64/0.52/0.50/0.60). If calibration changes a
   config value, the Qt reader silently diverges. Not a Ring.{h,cpp}
   readability issue and not in this file's scope — flagging so it
   isn't lost.

4. **E5 (`JohnsonBoveyLobeOffset`) rename go/no-go.** Recommended
   (mechanical, ~25 lines, three subprojects). If the human judges the
   cross-subproject churn not worth it, E6's signpost alone resolves
   the finding. Needs a decision since it crosses subproject boundaries.
