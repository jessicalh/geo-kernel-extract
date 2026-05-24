# larsen-hbond-grid — claude review (readability focus)

- **Targets:** src/LarsenHBondGrid.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict

**LarsenHBondGrid.h** — Reads cleanly. A chemist can follow the whole story top to bottom: geometry contract (what r/θ/ρ mean and over which atoms), the per-class atom-role table, the canonical donor-frame definition, and the tensor-rotation contract are all laid out as prose before any code, and the struct/field names carry units and frame. The header is the rare case where the documentation genuinely substitutes for reverse-engineering. Minor: a couple of duplicated/long comment blocks and one terminology mismatch (NMA-vs-AmideHydrogen) that could trip a first reader.

**LarsenHBondGrid.cpp** — Mostly clear; the load/validate/read path is linear and well-signposted. The two interpolation routines (`TrilinearMat3`, `LookupAxis`) are the only passages requiring real study — they are dense index arithmetic with terse one-letter locals (`c000`…`c111`, `fr/fth/frho`, `i/f/v`) and a periodic branch whose wrap math is correct but not signposted step-by-step. The dihedral block in `ComputeLarsenHBondGeometry` is correct but compressed into bare `b1/b2/b3/n1/n2/m1` with no naming bridge to the chemistry. None of this is broken; it is the standard "trilinear + dihedral" symbol soup that a reader must decode rather than read.

---

## 1. Coherent story / readability

- LarsenHBondGrid.cpp:215-253 — `TrilinearMat3` fuses the full 3-axis reduction in one nested block; reader must mentally label the three lerp stages — suggestion: one signpost comment per stage (`// lerp along ρ`, `// then θ`, `// then r`).
- LarsenHBondGrid.cpp:380-388 — dihedral computed from `b1/b2/b3/n1/n2/m1` with no tie back to H–O–C–third; reader must reconstruct which bond vector is which — suggestion: inline trailing comments (`// b1 = H→O`, `// b2 = O→C`, `// b3 = C→third`).
- LarsenHBondGrid.cpp:296-323 — the periodic branch reuses `step`/`period`/`v`/`f`/`i` from the non-periodic branch but computes them differently (period via `axis.back()-axis[0]+step`); the two branches read as near-duplicates that diverge subtly — suggestion: a one-line signpost (`// periodic axis: wrap then locate cell, last cell wraps to 0`).
- LarsenHBondGrid.cpp:64 — `dims[0] * dims[1] * dims[2] * 9` mixes the validated dim product with the magic `9` (3×3) inline — suggestion: trailing comment `// *9 = flattened 3x3`.
- LarsenHBondGrid.cpp:268-276 — triple `for (a)(b)(c)` over `rs/ths/rhos` is clear but the 2-element arrays are unlabeled — minor; acceptable as-is.

## 2. Naming carries meaning

- LarsenHBondGrid.cpp:222-223 — `Nth`, `ir/ith/irho`, `fr/fth/frho` are terse but locally unambiguous given the doc; acceptable, no rename.
- LarsenHBondGrid.cpp:232-239 — `c000`…`c111` are standard trilinear-corner names; fine as-is.
- LarsenHBondGrid.cpp:385 — `m1` is the only frame-cross name with no mnemonic (it is n1×b̂2) — suggestion: rename `m1` → `n1_perp` or comment `// n1 rotated 90° about b2`.
- LarsenHBondGrid.cpp:310-315 — `f` (fractional axis coordinate) and `i` (cell index) are heavily reused across both branches; clear within scope but a reader scanning fast may conflate them with the struct's `frac`/`idx` — minor.
- LarsenHBondGrid.h:359 — method named `QueryNearest` but the body does trilinear *interpolation*, not nearest-neighbour lookup; the name actively misleads — suggestion: the header comment already says "trilinear", but consider noting in the doc that "Nearest" is historical (or check whether a rename is in scope).

## 3. Visible math structure (grouping)

- LarsenHBondGrid.cpp:241-249 — the collapse from 8 corners → 4 → 2 → 1 is the visible-structure win of the file; well grouped, no change.
- LarsenHBondGrid.cpp:368-389 — three geometry outputs (r, θ, ρ) are cleanly sequenced and separated by blank lines + a label comment each; good.
- LarsenHBondGrid.cpp:419 — Gram-Schmidt step `v3 - (v3.dot(z))*z` is unnamed inline; it is the one real math step in the frame build — suggestion: name it (`Vec3 x_raw = v3 - v3.dot(z) * z; // remove z-component`) — the comment on line 410 already covers it, so optional.
- LarsenHBondGrid.cpp:529-534 — three `LookupAxis` calls with three distinct miss-reasons are nicely grouped with per-line miss comments; good.

## 4. Function / method naming

- LarsenHBondGrid.cpp:43 / .h — `ReadFlatTensorOptional` says what it does and that it may return empty; good.
- LarsenHBondGrid.cpp:259 — `AnyCornerImputed` says quantity and returns bool; good.
- LarsenHBondGrid.cpp:290 — `LookupAxis` is clear; the `AxisLookup` return struct documents `idx/idx_next/frac` well.
- LarsenHBondGrid.h:359 — `QueryNearest` misnaming repeated here (see axis 2); only flag.
- LarsenHBondGrid.cpp:438 / .h:234 — `RotateTensorToProteinLabFrame` names quantity, frame, and direction; exemplary.

## 5. Comments as signposts

- LarsenHBondGrid.h:127-129 vs :282-285 — the `any_corner_imputed` semantics are explained twice (header narrative + struct field) in near-identical prose; one can be trimmed to a back-reference.
- LarsenHBondGrid.cpp:114-119 — `ValidateSchema` preamble carries process/history prose ("Optional reads previously allowed silent zeroing if a mandatory dataset was missing (codex M3 finding)") — suggestion: drop the provenance, keep the rule.
- LarsenHBondGrid.cpp:562-566 — comment narrates removed earlier code ("Earlier code recomputed r/θ from axis*frac… equivalent but confusing") — history prose; suggestion: replace with terse `// r/θ are the clamped query values; ρ is wrapped (set above)`.
- LarsenHBondGrid.cpp:215-218 — good terse signpost on `TrilinearMat3`; keep.
- LarsenHBondGrid.cpp:301-304 — "If a future grid breaks this invariant the wrap is wrong; bail." is a good grounded guard comment; keep.
- LarsenHBondGrid.cpp:428-429 — `// R @ v_log = v_canonical, so R has rows [x; y; z]` is exactly the right 1-line signpost; keep.

## 6. Correctness (secondary)

- LarsenHBondGrid.cpp:316-318 — periodic `i_next` wraps to 0, but `frac` is left as-is; correct only if `WrapRho` guarantees `v < axis.back()+step` so the last-cell `frac∈[0,1)` interpolates between `axis[n-1]` and `axis[0]`. Looks consistent with a 360° period — check the period/step bookkeeping against an actual ρ axis that does *not* include +180° as a sample.
- LarsenHBondGrid.cpp:307-308 — `std::fmod` can return exactly `0` then `v=axis[0]`; with `f=0,i=0` that is fine, but if `v` lands a hair below `axis[0]` after FP, the `if (v<axis[0]) v+=period` guard catches it — looks safe; no action.
- LarsenHBondGrid.cpp:529-534 — r/θ use `geom` values but ρ uses `rec.rho_deg` (the wrapped copy); intentional and correct — the comment on :504 documents it. Good.
- LarsenHBondGrid.cpp:441 — `R_protein.transpose() * sigma * R_protein` matches the header contract `σ_lab = R_proteinᵀ σ_canonical R_protein`; consistent. Cannot confirm the stored-tensor convention from these files — header asserts it; check against the parser.
- LarsenHBondGrid.h:171 / :111 — comment uses "ASN ODE1" (header) vs "ASN OD1" (line 68 table) for the same atom; one is a typo — check the intended Asn carbonyl-O name.
- LarsenHBondGrid.h:248-249 / .cpp:208 — `has_acceptor_readouts` is derived solely from `acceptor_N` non-empty; `ValidateSchema` enforces all-or-nothing so this is safe, but the single-field proxy is implicit — minor, acceptable given the schema check.
- LarsenHBondGrid.cpp:335 — non-periodic clamp `i = min(i, n-2)` correctly reserves `i+1`; the subsequent `frac` clamp to [0,1] is belt-and-suspenders after the bound check — fine.
