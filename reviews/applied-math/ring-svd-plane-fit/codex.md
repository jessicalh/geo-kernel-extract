# ring-svd-plane-fit — codex review (readability focus)

- **Targets:** src/Ring.{h,cpp}
- **Model:** codex `gpt-5.5`, reasoning `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 8,127
- **Brief:** `../_brief.md` (readability-first)
- **See also:** `codex-correctness.md` (earlier correctness-first pass)

---

**Verdict**

`src/Ring.h` mostly tells a coherent structural story: base ring identity, geometry/accumulated state, then concrete ring types grouped by chemistry. The main breaks are stale inventory, unexplained symbols/units, and one long prose block that interrupts the otherwise table-like property listing.

`src/Ring.cpp` is readable and sequenced well: collect vertices → centroid → SVD plane normal → orientation guard → radius → factory. The only major clarity issue is the fatal error prose at the end, which is more project-process narration than code signpost.

**1. Coherent Story / Readability**

`src/Ring.h:9` — header says “8 types” but the file defines 9 including `ProPyrrolidineRing` — update the inventory so the first story matches the implementation.

`src/Ring.h:25` — “computed by GeometryResult” conflicts with `Ring::ComputeGeometry` in `Ring.cpp` — rename/comment as “computed by Ring::ComputeGeometry” unless `GeometryResult` is still part of the pipeline.

`src/Ring.h:39-44` — `RingAccumulated` is a sudden bag of post-pass fields with no readable through-line — add terse labels or rename fields so a reader knows which pass/quantity each value belongs to.

`src/Ring.h:178-185` — Pro pyrrolidine explanation is scientifically useful but too long relative to the otherwise compact type table — compress to a signpost plus citation, e.g. `// saturated: no pi ring current`, then keep details in docs or a nearby citation note.

`src/Ring.cpp:65-70` — fatal-path comment explains project discipline/history more than computation — replace with a short signpost like `// invalid enum guard`.

**2. Naming Carries Meaning**

`src/Ring.h:40` — `total_B_at_center` uses `B` without units/frame — consider `total_field_at_center` or `total_B_field_at_center`.

`src/Ring.h:42` — `total_G_T0_diagnostic` is symbol soup — a reader must already know what `G`, `T0`, and “diagnostic” mean; spell the quantity out, e.g. `total_traceless_gradient_T0` if that is accurate.

`src/Ring.h:43` — `mutual_B_from` omits key/value meaning — consider `mutual_field_from_ring` or `field_at_center_from_ring`.

`src/Ring.h:64` — `JBLobeOffset()` depends on unexplained `JB` — consider `JohnsonBoveyLobeOffset()` or add one signpost comment defining JB.

`src/Ring.h:68` — `TypeName()` returns compact residue/ring labels like `"TRP6"` and `"TRP9"`, not full type names — consider `TypeLabel()` or `RingLabel()`.

`src/Ring.h:67` — `RingSizeValue()` is awkward and underspecified — `RingAtomCount()` would carry the physical meaning better.

**3. Visible Math Structure / Grouping**

`src/Ring.h` — clean for this axis; it is mainly declarations and table-like physical constants.

`src/Ring.cpp:12-46` — the computation is visibly staged and easy to follow.

`src/Ring.cpp:34-39` — orientation guard is clear enough, but `edge01` / `edge02` force light decoding — `first_edge` and `second_edge` would read more naturally.

**4. Function / Method Naming**

`src/Ring.h:62-64` — `Intensity()`, `LiteratureIntensity()`, and `JBLobeOffset()` do not expose units or convention — if API churn is acceptable someday, names should include ring-current / Johnson-Bovey meaning; otherwise add a compact comment above the virtual block.

`src/Ring.cpp:52` — `CreateRing` is clear and returns the expected subclass.

`src/Ring.cpp:8` — `ComputeGeometry` is clear for center/normal/radius/vertices.

**5. Comments As Signposts**

`src/Ring.h:3-12` — useful orientation comment, but stale due to missing proline ring — update, don’t remove.

`src/Ring.h:61` — “Virtual const properties” is slightly abstract — `// ring-type properties` would be more direct.

`src/Ring.h:70` — “Non-virtual queries” is clear enough.

`src/Ring.h:178-185` — too verbose for an inline class table — terse replacement: `// saturated: no pi current`.

`src/Ring.cpp:23-25` — good signpost; it explains the one non-obvious SVD step.

`src/Ring.cpp:65-70` — too verbose/process-heavy — terse replacement: `// invalid enum guard`.

**6. Correctness**

`src/Ring.h:9-12` — stale type count/inventory is a documentation correctness issue, not an algorithm bug.

`src/Ring.cpp` — I do not see an actual correctness bug from the inlined text. Potential bounds checking on `positions[idx]` depends on caller contract; I would not flag it as a bug without seeing that contract.
