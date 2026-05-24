# coulomb — claude review (readability focus)

- **Targets:** src/CoulombResult.{h,cpp} + src/MopacCoulombResult.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict (central question, per file)

**CoulombResult.h** — Reads cleanly. The header comment lays out the kernel,
units, and source decomposition before any code, so the through-line is clear
to a physicist on first pass. Nothing to decode.

**CoulombResult.cpp** — Mostly a coherent story: classify atoms → bond
directions → N²-ish sum → unit conversion → de-trace → sanitise → store. The
per-charge `E_j`/`V_j` math is named and commented inline, which is the right
call. Where it strains the reader is the long undifferentiated tail of
post-loop bookkeeping (sanitise / clamp / GeometryChoice / many derived
scalars) that buries the two physical outputs, plus a self-contradiction
between the header's "spatial cutoff" claim and the file's repeated "no
cutoff" comments.

**MopacCoulombResult.h / .cpp** — The same story, told the same way; the
header comment correctly frames it as "same kernel, different charge source,"
which is exactly what a reader needs. It is a near-verbatim clone of the
non-MOPAC variant, so it inherits the same strengths and the same tail-of-
bookkeeping weakness. Clean enough to follow once you've read its sibling.

---

## 1. Coherent story / readability (primary)

- CoulombResult.cpp:34-35 / .h:34-35 vs .cpp:120,141 — header + inline comments say "spatial cutoff" / "sum over ALL atoms (no spatial cutoff)" but the loop actually iterates `AtomsWithinRadius(pos_i, coulomb_cutoff)` — contradictory signposts; reader cannot tell whether physics is truly long-range — reconcile the comments to whichever is true.
- CoulombResult.cpp:108 — section banner "Coulomb sum with spatial cutoff" contradicts the line-34 "no spatial cutoff" comment 70 lines up — pick one phrasing.
- CoulombResult.cpp:198-326 — ~130-line post-loop tail (de-trace, sanitise, clamp, two GeometryChoice blocks, ~10 derived scalars, solvent, shielding) buries the two physical results among bookkeeping — add a one-line banner separating "physics outputs" from "diagnostics/derived scalars".
- CoulombResult.cpp:307-325 — 18-line essay comment on why there is no unified tensor is valuable physics but stops the read cold at the assignment site — keep it but move the bulk to the header rationale, leave a 2-4 word signpost here.
- MopacCoulombResult.cpp:38-264 — entire Compute is a near-verbatim copy of CoulombResult::Compute; a reader must diff the two by eye to find the only real differences (charge source, no solvent, no aromatic-source count) — a 3-line "differs from CoulombResult only in: …" bander atop Compute would save that diffing.
- CoulombResult.cpp:122-124,139,178,295-297 — `aromatic_source_count` and `n_sidechain_aromatic_sources` are computed and threaded through but only feed a log line / one scalar; their purpose is unobvious mid-read — one-line "diagnostic only" note.

## 2. Naming carries meaning

- CoulombResult.cpp:159-160 — `r3`, `r5` are fine, but `r` (line 156) is the displacement *vector* while `r_mag` is its norm; calling the vector `r` invites confusion with scalar distance — `r_vec` would disambiguate.
- CoulombResult.cpp:289 / :282 — `coulomb_E_backbone_frac` is a projection (dot with unit vector, units V/A), not a fraction; the comment at :284-286 even says so — name reads as dimensionless ratio; consider documenting "(signed projection, V/A — not a ratio)" at the field or store site.
- CoulombResult.h:53 — `SampleEFieldAt` is clear; no change.
- CoulombResult.cpp:381 / Mopac:296 — `PackST_C` / `PackST_MCC` are cryptic abbreviations — `PackSphericalTensor` (file-static, collision-free) would read at the call site.
- Both files — `E_j`, `V_j`, `q_j` carry units only via comment; acceptable given the header units block, no change.

## 3. Visible math structure (grouping)

- CoulombResult.cpp:162-167 / Mopac:141-146 — the E and V per-charge build is well grouped and commented; this is the clearest block in the file — no change.
- CoulombResult.cpp:188-207 — unit-convert-then-detrace is correctly sequenced and labelled — no change.
- CoulombResult.cpp:278-297 — ten derived-scalar assignments run together with no grouping; the bond-projection, backbone-fraction, and aromatic scalars are three distinct ideas fused into one block — blank-line + 2-word labels per sub-group.
- CoulombResult.cpp:209-226 / Mopac:183-200 — sanitise-vec and sanitise-mat lambdas are defined inline inside the per-atom loop on every iteration; harmless but visually they interrupt the math flow — fine to leave, but a reader expects them hoisted.

## 4. Function / method naming

- CoulombResult.h:48-50 — `EFieldAt` / `EFGAt` / `EFGSphericalAt` say quantity and that they're per-atom; clear — no change.
- CoulombResult.cpp:381 / Mopac:296 — `PackST_C` / `PackST_MCC` don't say what they pack or the layout (T0,T1×3,T2×5 → 9) — rename + one-line layout comment.
- Both — `Compute` as a static factory is the codebase convention; clear in context — no change.

## 5. Comments as signposts

- CoulombResult.cpp:146 — "Self-exclusion via filter framework (not inline check)" restates the code and editorialises the design choice — drop or shorten to "// self-exclusion".
- CoulombResult.cpp:176-177 — "Count aromatic sidechain sources near this atom (atoms on aromatic residues that are sidechain)" is verbose for what `!is_backbone[j]` says — trim to "// aromatic sidechain source".
- CoulombResult.cpp:307-324 — long multi-paragraph WARNING/history comment mixes a real caveat (T0 not captured) with design narrative and a TODO — keep the caveat as "// T2 only; Buckingham T0 not included here", relocate the rest.
- CoulombResult.cpp:393-394 / Mopac:308-309 — "EFG schema rev 2026-05-18" date-stamped comment is process history embedded in code — keep the physics ("symmetric per charge → T0/T1 zero"), drop the rev date.
- MopacCoulombResult.cpp:251-252 — "Pure T2 (EFG is traceless). gamma converts this to shielding." is a good terse signpost — no change.
- CoulombResult.cpp:284-286 — backbone-projection comment is accurate and genuinely helpful (explains the stability choice) — no change.

## 6. Correctness (secondary)

- CoulombResult.cpp:139,178,297 — `n_sidechain_aromatic_sources` counts aromatic-and-not-backbone *source* atoms globally (independent of target i, since no distance gate on the count) yet is stored per-target as `aromatic_n_sidechain_atoms`; comment at :176 says "near this atom" but there is no proximity test — check whether a per-target/proximity count was intended.
- CoulombResult.cpp:114 / Mopac:104 — `MinDistanceFilter` plus `SelfSourceFilter` are added but the per-charge math has no explicit `r_mag` singularity guard (unlike `SampleEFieldAt`:366 which guards explicitly); correctness then depends on `MinDistanceFilter` rejecting `r < min` — check that the filter's threshold actually prevents `r3`/`r5` underflow.
- CoulombResult.cpp:202 / Mopac:175 — traceless projection is applied to `EFG_sidechain` too, but sidechain EFG is never stored or written (only total/backbone/aromatic are) — dead computation, not a bug; harmless.
- MopacCoulombResult.cpp:124 — uses raw `for j in 0..n_atoms` (true N²) while CoulombResult uses `AtomsWithinRadius`; the two variants sum over different source sets, so MOPAC and ff14SB EFG are not strictly comparable term-for-term — check this divergence is intended (header claims "same kernel … same decomposition").
- CoulombResult.cpp:235 / Mopac:209 — clamp scale lambda captures `E_mag` by value before the `*= scale` mutation; correct, captured value is pre-scale — no bug, noting for the reader.
- Both — de-trace happens *after* `COULOMB_KE` scaling and after accumulation, matching the header claim that it corrects FP drift on an already-traceless sum — consistent; no issue.
