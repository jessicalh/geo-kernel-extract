# Fix-plan — `src/LarsenHBondGrid.{h,cpp}`

Algorithm under review: trilinear interpolation over a 3-D H-bond geometry
grid (r, θ, ρ; ρ periodic) plus Gram-Schmidt construction of the canonical
donor frame and the canonical→lab tensor rotation.

---

## 1. Summary

The file mostly tells a coherent story. The header is an unusually complete
contract: geometry definitions (r/θ/ρ over four named atoms), the per-class
atom-role table, the canonical donor-frame definition, the validity-mask
semantics, and the tensor-rotation contract are all stated as prose before any
code. The `.cpp` has a clean broad order — load/validate archives, compute
geometry + frame, query the grid — and the load/validate path is linear and
well signposted.

The genuinely dense passages are the two interpolation/indexing routines
(`TrilinearMat3`, `LookupAxis`) and the dihedral block in
`ComputeLarsenHBondGeometry`. None of them is broken; they are standard
"trilinear + dihedral" index/symbol soup that a once-through reader must decode
rather than read. The fix pass adds a small number of signpost comments and a
few named intermediates so the *purpose of each step* is visible, trims three
process/history comments to current-state statements, and corrects one
documentation typo.

**The fix pass will NOT** touch the algorithm, any sign, any number, the
rotation direction, the periodic-wrap math, or any output/contract name
(NPY field names, H5 dataset names, SDK `_catalog.py` keys). It is comments,
internal-local naming, and one prose-typo fix only.

**Verified before planning (per the governing prior):**

- **Rotation direction `Rᵀ σ R` vs the Kabsch path's `R σ Rᵀ` is earned, not a
  bug.** `ComputeLarsenDonorFrame` builds `R` with the canonical basis vectors
  as *rows* (`LarsenHBondGrid.cpp:431-433`), so by construction
  `R · v_log = v_canonical` (header :202-204, :428-429). Grid tensors are
  stored as `σ_canonical = R_log σ_log R_logᵀ` (header :135, :229), therefore
  recovering the lab tensor requires `σ_lab = R_proteinᵀ σ_canonical R_protein`
  — which is exactly `RotateTensorToProteinLabFrame` (`:441`). The consumer
  `LarsenHBondShieldingResult.cpp` uses the pair as documented:
  `ComputeLarsenDonorFrame(donor_H, donor_anchor, donor_third)` at :576-577,
  then `RotateTensorToProteinLabFrame(canonical_tensor, R_protein)` at
  :611-612, :631-632, :672-673. The Kabsch path
  (`TripeptidePoseAssembler.cpp:73`, `R σ Rᵀ`) builds `R` as a
  point-set-aligning rotation (`V D Uᵀ`, :56) whose *columns* are the mapped
  basis — opposite storage convention, so the opposite multiplication order is
  correct. **Coherent (expected); do not flag.**

- **ρ axis is periodic (±180 wrap).** `LookupAxis(..., periodic=true)` at
  `:296-323` and `WrapRho` at `:347-351` implement this; `QueryNearest` calls
  the ρ lookup with `periodic=true` (`:533`) on the already-wrapped
  `rec.rho_deg`. Consistent end-to-end.

---

## 2. Review-finding ledger

Every finding from `codex.md` and `claude.md` (`codex-correctness.md` not
present), one row each.

### codex.md

| # | Finding (loc) | Disposition |
|---|---|---|
| C1 | h:3 — header opens as a design memo; split mentally into signposted blocks | **Declined** — claude review judged the header "the rare case where documentation substitutes for reverse-engineering." The blocks are already `─`-ruled and titled. Splitting/trimming risks losing earned contract detail (rotation, validity-mask, sign conventions). Restraint per the clarity bar: do not fragment coherent prose. |
| C2 | h:95 — "parser and the future calculator" reads as WIP history | **Adopted** → §3 E1 (the calculator has landed; drop "future"). |
| C3 | h:107 — "future session" breaks production narrative | **Adopted** → §3 E2. |
| C4 | h:132 — rotation contract appears far from canonical-frame def; move or add signpost | **Adopted (signpost only)** → §3 E3. Decline the *move* (reordering the header is churn); add a short back-reference. |
| C5 | cpp:215 — `TrilinearMat3` reads as index-symbol soup; add block label + axis-order note | **Adopted** → §3 E5 (per-stage signpost) + E6 (axis-order note). |
| C6 | cpp:280 — `LookupAxis` mixes tolerance / periodic wrap / lower-cell; add signposts | **Adopted** → §3 E8. |
| C7 | cpp:224 — rename `idx` → `tensor_offset` | **Adopted** → §3 E7. |
| C8 | cpp:232 — `c000`…`c111` cryptic without axis order | **Adopted as comment, not rename** → §3 E6. The `cNNN` names are standard trilinear corners (claude CL-corners agrees); a one-line axis-order comment is the lighter fix than 8 long renames. |
| C9 | cpp:290 — rename return fields `idx/idx_next/frac` → `lower_index/upper_index/cell_fraction` | **Declined** — `AxisLookup` is a file-local struct, but its three fields are read at three call sites in `QueryNearest` (`:539-541`). The current names are unambiguous in context and the doc comment (`:280-283`) defines them. Rename is defensible churn but not a clarity win that earns the carry-through; left to human discretion (noted §6 Q3). |
| C10 | cpp:305 — rename `v` → `wrapped_value` | **Adopted** → §3 E9. |
| C11 | cpp:380 — rename `b1/b2/b3/n1/n2/m1/x/y` to chemistry-bearing names | **Partially adopted** → §3 E10. Rename `b1/b2/b3` to bond-vector names and add a trailing comment on the atan2 terms; keep `n1/n2/m1/x/y` as conventional dihedral symbols (renaming `x/y` to `dihedral_*_term` is over-naming the atan2 args). See claude CL-m1 for the `m1` mnemonic. |
| C12 | cpp:413 — rename `v3` → `H_to_third`/`third_to_H` | **Adopted** → §3 E11 (`v3` = `donor_H − donor_third`, so `H_to_third`). |
| C13 | cpp:361 — dihedral block reads better with named normals/atan2 terms | **Duplicate** of C11. |
| C14 | cpp:394 — `ComputeLarsenDonorFrame` is well grouped (clean) | **Declined (no-op)** — affirmation, no change requested. |
| C15 | cpp:496 — `QueryNearest` query path clean | **Declined (no-op)** — affirmation. |
| C16 | cpp:562 — comment narrates removed code; replace with current-state label | **Adopted** → §3 E13 (duplicate target of CL-hist2). |
| C17 | h:359 / cpp — rename method `QueryNearest` → `QueryInterpolated`/`Interpolate` (does trilinear, not nearest) | **Declined** — `QueryNearest` is a *public method* of class `LarsenHBondGrid`, called from `LarsenHBondShieldingResult.cpp:558`. It is also the name of the analogous `TripeptideDftTable::QueryNearest` (sibling API). The mismatch is real but the carry-through crosses files and matches a sibling convention; the header doc already says "trilinear". **Add a one-line doc note that "Nearest" is historical** rather than rename → §3 E4. Noted as a weighed cost call (§6 Q1). |
| C18 | cpp:219 / :259 / :347 — `TrilinearMat3`, `AnyCornerImputed`, `WrapRho` clean | **Declined (no-op)** — affirmations. |
| C19 | cpp:290 — `LookupAxis` → `LookupAxisCell` | **Declined** — minor; conflicts with C9's framing and the doc comment already says it returns a cell + fraction. Not worth the file-local churn. |
| C20 | h:2 — opening comment too expansive; terse section labels | **Duplicate** of C1 → declined. |
| C21 | cpp:114 — "codex M3 finding" is process prose; → "required datasets" | **Adopted** → §3 E12 (duplicate target of CL-hist1). |
| C22 | cpp:215 — interpolation comment too long, axis order hard to see | **Duplicate** of C5/C8-axis-order → §3 E5/E6. |
| C23 | cpp:256 — comment restates 8-corner enumeration; → "imputed-corner check" | **Adopted** → §3 E14. |
| C24 | cpp:371 / :428 — good signposts (clean) | **Declined (no-op)** — affirmations; keep. |
| C25 | cpp:378 / :401 / :410 — useful but slightly long signposts; shorten | **Declined** — these are already 3–7 words and carry the chemistry (`// theta at acceptor O`, `// donor z-axis`). Shortening further trades clarity for terseness; the clarity bar warns against both over- and under-explaining, and these sit at the right level. |
| C26 | h / cpp — no correctness bug found | **Declined (no-op)** — agrees with this plan. |

### claude.md

| # | Finding (loc) | Disposition |
|---|---|---|
| CL1 | cpp:215-253 — one signpost per lerp stage (`// lerp along ρ`, `// θ`, `// r`) | **Adopted** → §3 E5. Duplicate target of C5. |
| CL2 | cpp:380-388 — inline trailing comments tying `b1/b2/b3` to H→O / O→C / C→third | **Adopted** → §3 E10. Duplicate of C11. |
| CL3 | cpp:296-323 — one-line signpost on the periodic branch (wrap then locate cell, last cell wraps to 0) | **Adopted** → §3 E8. Duplicate target of C6. |
| CL4 | cpp:64 — trailing comment on `dims[...] * 9` (`// *9 = flattened 3x3`) | **Adopted** → §3 E15. |
| CL5 | cpp:268-276 — 2-element corner arrays unlabeled; acceptable as-is | **Declined (no-op)** — reviewer self-rated acceptable. |
| CL6 | cpp:222-223 — `Nth/ir/fr…` terse but OK given doc; no rename | **Declined (no-op)** — affirmation. |
| CL7 | cpp:232-239 — `c000…c111` standard; fine as-is | **Declined (no-op)** — affirmation; the axis-order comment (E6) addresses the only gap. |
| CL8 | cpp:385 — `m1` has no mnemonic; rename `n1_perp` or comment | **Adopted (comment)** → §3 E10 includes a `// n1 × b2̂` note on `m1`. |
| CL9 | cpp:310-315 — `f`/`i` reused across branches, may conflate with struct `frac`/`idx`; minor | **Declined** — within-scope locals; the signpost E8 demarcates the two branches, which is the real source of the conflation worry. |
| CL10 | h:359 — `QueryNearest` misleads; note "Nearest" is historical or check rename scope | **Adopted (doc note)** → §3 E4. Duplicate of C17. |
| CL11 | cpp:241-249 — corner collapse well grouped; no change | **Declined (no-op)** — affirmation. |
| CL12 | cpp:368-389 — r/θ/ρ cleanly sequenced; good | **Declined (no-op)** — affirmation. |
| CL13 | cpp:419 — name the Gram-Schmidt step `x_raw`; comment already covers it, optional | **Declined** — the code *already* names it `x_raw` (`:419`) and comments it (`:410-412`). The review's suggestion is already in the source. |
| CL14 | cpp:529-534 — three `LookupAxis` calls grouped with per-line miss comments; good | **Declined (no-op)** — affirmation. |
| CL15 | cpp:43 — `ReadFlatTensorOptional` good name | **Declined (no-op)** — affirmation. |
| CL16 | cpp:259 — `AnyCornerImputed` good | **Declined (no-op)** — affirmation. |
| CL17 | cpp:290 — `LookupAxis` clear; struct documents fields | **Declined (no-op)** — affirmation (note: tension with C9/C19; resolved in favor of no rename, §6 Q3). |
| CL18 | cpp:438 / h:234 — `RotateTensorToProteinLabFrame` exemplary | **Declined (no-op)** — affirmation. |
| CL19 | h:127-129 vs :282-285 — `any_corner_imputed` explained twice; trim one to a back-reference | **Adopted** → §3 E16. |
| CL20 | cpp:114-119 — drop the codex-M3 provenance, keep the rule | **Adopted** → §3 E12. Duplicate of C21. |
| CL21 | cpp:562-566 — comment narrates removed code; replace with terse current-state | **Adopted** → §3 E13. Duplicate of C16. |
| CL22 | cpp:215-218 / :301-304 / :428-429 — good signposts; keep | **Declined (no-op)** — affirmations; E5 *adds* per-stage labels without removing the existing block comment, and the guard/rotation comments are kept verbatim. |
| CL23 | cpp:316-318 — periodic `i_next` wraps to 0 but `frac` left as-is; correct only if `WrapRho` keeps last-cell `frac∈[0,1)`; check against a ρ axis that omits +180 | **Adopted as verification + comment** → §6 Q2 (verification) and §3 E8 (the signpost states the wrap intent). See §4 for the trace; this is *coherent*, the comment makes the invariant explicit. |
| CL24 | cpp:307-308 — `fmod` edge with `v=axis[0]`; guard catches it; safe | **Declined (no-op)** — reviewer concluded safe; matches this plan. |
| CL25 | cpp:529-534 — r/θ use `geom`, ρ uses `rec.rho_deg` (wrapped); intentional + documented at :504 | **Declined (no-op)** — affirmation; coherent (see §4). |
| CL26 | cpp:441 — `Rᵀ σ R` matches header; can't confirm stored-tensor convention from these files | **Resolved → coherent** (see §1 and §4). Traced both producer and consumer; no change. |
| CL27 | h:171 / :111 — "ASN ODE1" (header narrative) vs "ASN OD1" (role table :68); one is a typo | **Adopted** → §3 E17. `ODE1` is not a PDB atom name; `OD1` (role table + consumer `LarsenHBondShieldingResult.cpp:96`) is correct. Comment-only fix. |
| CL28 | h:248-249 / cpp:208 — `has_acceptor_readouts` derived from `acceptor_N` non-empty; safe via schema check; implicit | **Declined** — `ValidateSchema` (`:156-173`) enforces all-or-nothing on the acceptor set, so the single-field proxy is sound. Adding a comment here would be hand-holding; the schema function is the documented guarantor. (Optional micro-note noted in §6 Q4.) |
| CL29 | cpp:335 — non-periodic clamp + belt-and-suspenders `frac` clamp; fine | **Declined (no-op)** — affirmation. |

---

## 3. Edits that don't move numbers

Comment fixes, short signposts, and internal-local renames only. No output
name, no sign, no number, no algorithm change. Line numbers are pre-edit.

**Header (`LarsenHBondGrid.h`)**

- **E1** — `h:95` — replace "the parser and the future calculator both use it"
  → "the parser and the calculator both use it" (the calculator has landed).
- **E2** — `h:107-109` — drop "(`LarsenHBondShieldingResult`, future session)"
  process aside; keep "done by the calculator, not here — this layer returns
  whichever readouts the archive has."
- **E3** — `h:132` — prepend a 1-line back-reference to the rotation-contract
  block: `// (canonical donor frame defined above; R = ComputeLarsenDonorFrame)`
  so the reader is not asked to recall the frame def from ~50 lines up. Do NOT
  reorder the header.
- **E4** — `h:359` (doc above `QueryNearest`) — add one sentence: the method
  performs trilinear interpolation; the name "Nearest" is historical and
  retained for parity with `TripeptideDftTable::QueryNearest`. (No rename — see
  §2 C17, §6 Q1.)
- **E16** — `h:282-285` — trim the `any_corner_imputed` struct-field comment to
  a back-reference to the header narrative (`h:119-129`), removing the
  near-verbatim duplication (CL19).
- **E17** — `h:111` and `h:171` — fix "ASN ODE1" → "ASN OD1" (CL27); aligns the
  narrative + enum comment with the role table (`h:68`) and the consumer
  (`LarsenHBondShieldingResult.cpp:96`).

**Implementation (`LarsenHBondGrid.cpp`)**

- **E5** — `cpp:241-249` — add three terse stage signposts inside
  `TrilinearMat3`: `// 1) lerp along ρ` (before `c00…c11`), `// 2) along θ`
  (before `c0/c1`), `// 3) along r` (before `out(row,col)`). Keep the existing
  block comment at `:215`.
- **E6** — `cpp:232` — add one comment above the corner reads:
  `// corner index order: (r, θ, ρ) → c{r}{θ}{ρ}`. Resolves C8/CL7 without
  renaming the standard `cNNN` corners.
- **E7** — `cpp:224` — rename local lambda `idx` → `tensor_offset` inside
  `TrilinearMat3` (file-local; the `AnyCornerImputed` lambda at `:263` is a
  separate scope and may keep its own `idx` or be renamed `mask_offset` for
  symmetry — optional, noted §6 Q3).
- **E8** — `cpp:296` / `:325` — add two branch signposts in `LookupAxis`:
  `// periodic axis: wrap into [axis[0], axis[0]+period), then locate cell;`
  `// the last cell's upper corner wraps to index 0` (before the periodic
  block), and `// non-periodic axis: tolerance-clamp to bounds, then lower
  cell` (before `:325`). Makes the two near-duplicate branches read as
  deliberately distinct (C6/CL3).
- **E9** — `cpp:305` — rename local `v` → `wrapped_value` in the periodic
  branch (file-local).
- **E10** — `cpp:380-388` — in `ComputeLarsenHBondGeometry`'s dihedral block:
  rename `b1/b2/b3` → `H_to_O` / `O_to_C` / `C_to_third` (matches the
  bond-vector chemistry and the existing `H_to_O`/`O_to_H`/`O_to_C` naming used
  in the θ block above); add trailing comments `// n1 ⟂ plane(H,O,C)`,
  `// n2 ⟂ plane(O,C,third)`, `// m1 = n1 × b2̂` on `n1`/`n2`/`m1`. Keep `n1`,
  `n2`, `m1`, `x`, `y` as conventional dihedral symbols (renaming the atan2 args
  is over-naming). Note: `b2` is referenced in `m1 = n1.cross(b2.normalized())`
  — after renaming, that becomes `O_to_C.normalized()`.
- **E11** — `cpp:413` — rename local `v3` → `H_to_third` in
  `ComputeLarsenDonorFrame` (it is `donor_H_pos − donor_third_pos`); update the
  two `v3` uses at `:414` and `:419`. File-local.
- **E12** — `cpp:114-119` — trim the `ValidateSchema` preamble: drop "Optional
  reads previously allowed silent zeroing … (codex M3 finding)"; keep the rule
  ("each archive has a mandatory readout set; throw on any missing/extra").
  (C21/CL20.)
- **E13** — `cpp:562-566` — replace the removed-code narration with a terse
  current-state comment: `// r/θ are the clamped in-range query values; ρ is`
  `// the wrapped canonical form set at the top of this function.` (C16/CL21.)
- **E14** — `cpp:256` — shorten the corner-enumeration restatement to
  `// imputed-corner check: true if any of the 8 trilinear corners is masked 0`.
  (C23.)
- **E15** — `cpp:64` — trailing comment on the size product:
  `// × 9 = flattened 3×3 tensor per grid cell`. (CL4.)

---

## 4. Usage notes (the sign/value reasons discovered)

These are the conventions the reviews asked about, with the reason and the
producer/consumer agreement.

1. **Tensor rotation direction `σ_lab = Rᵀ σ_canonical R` (`cpp:441`).**
   *Reason:* `ComputeLarsenDonorFrame` returns `R` with the canonical basis
   vectors as its **rows** (`cpp:431-433`), so `R · v_log = v_canonical`
   (header :202-204, :428-429). Grids are stored as
   `σ_canonical = R_log σ_log R_logᵀ`. Inverting at runtime with the
   protein-side `R_protein` (same constructor) gives
   `σ_lab = R_proteinᵀ σ_canonical R_protein`.
   *Consumed at:* `LarsenHBondShieldingResult.cpp:611-612, 631-632, 672-673`,
   via the `RotateTensorToProteinLabFrame` helper — the call site never
   hand-rolls the order.
   *Producer/consumer:* **agree.** The contrast with the Kabsch path's
   `R σ Rᵀ` (`TripeptidePoseAssembler.cpp:73`) is a *different storage
   convention* (Kabsch `R = V D Uᵀ` has the mapped basis as columns), not an
   inconsistency. **Coherent — do not change.** (Per instruction.)

2. **ρ uses IUPAC dihedral sign, NOT the Larsen-filename sign.**
   *Reason:* the parser keys grid bins on IUPAC-sign ρ computed from the actual
   Gaussian-log atomic positions (header :50-54). `ComputeLarsenHBondGeometry`
   computes the same IUPAC-sign ρ via `atan2(m1·n2, n1·n2)` (`cpp:386-388`), so
   producer (parser) and runtime geometry agree by construction.
   *Consumed at:* `LarsenHBondShieldingResult.cpp:543` (geometry) →
   `QueryNearest` → grid. No consumer negates ρ. **Coherent.**

3. **ρ periodic wrap and the last-cell `frac` (CL23).**
   *Reason:* in the periodic branch (`cpp:296-322`), `period = back − front +
   step` (full 360°, asserted at `:302`), `frac = f − floor(f)` stays in
   `[0,1)`, and when the lower cell is the last sample `i_next` wraps to `0`
   (`:316-318`). This interpolates between `axis[n-1]` and `axis[0]` across the
   ±180 seam, which is correct **iff the ρ axis does not also store +180 as a
   duplicate of −180**. The `period == 360` assertion is the guard that catches
   a malformed axis (bails to `idx=-1` → non-hit). *Producer:* the dense-grid
   parser sets the ρ axis; this loader cannot see it from these two files.
   **Coherent given the documented invariant; the one open item is to confirm
   the axis omits +180 — see §6 Q2.** No code change; E8 makes the intent
   explicit.

4. **Record carries wrapped ρ but clamped r/θ (CL25).**
   `QueryNearest` sets `rec.rho_deg = WrapRho(geom.rho_deg)` (`:504`) and leaves
   `rec.r_angstrom`/`rec.theta_deg` as the in-range query values, then the ρ
   lookup uses `rec.rho_deg` while r/θ use `geom` (`:529-533`). This is
   intentional (documented at header :274-280, :348-350) so diagnostics record
   the canonical grid coordinate. **Coherent.**

5. **`any_corner_imputed` derived purely from the validity mask.**
   `AnyCornerImputed` (`:259-277`) checks the 8 trilinear corners against the
   uint8 mask; `false` when no mask is present. The consumer may downweight/skip
   (header :282-285). No sign/value concern; flag-only output.

---

## 5. Bug-by-exhaustion candidates

**None.** Both reviews reported no concrete correctness bug, and the two
sign/convention questions they raised (rotation direction, ρ sign/wrap) trace
to coherent, documented conventions with agreeing consumers (§4). The only
factual error found is a documentation typo (`ODE1` → `OD1`, E17), which moves
no number.

---

## 6. Questions & Ambiguities

- **Q1 — Rename `QueryNearest`?** The name is a misnomer (it interpolates).
  Renaming touches the public method and breaks parity with the sibling
  `TripeptideDftTable::QueryNearest`. The plan keeps the name and adds a doc
  note (E4). *Decision for the human:* accept the doc note, or schedule a
  coordinated rename across `LarsenHBondGrid`, `LarsenHBondShieldingResult.cpp`,
  and (for consistency) `TripeptideDftTable`. Recommended: keep + doc note.

- **Q2 — Does the ρ dense-grid axis omit +180°?** The periodic wrap (§4.3) is
  correct only if +180 is *not* stored as a separate sample from −180. This
  cannot be confirmed from the two source files; it depends on
  `scripts/larsen_hbond_grid_parse/pre_compute_dense_grids.py` /
  `convert_dense_to_h5.py` and the actual `*_dense.h5` `rho_axis`. The
  `period == 360` assertion (`:302`) catches the gross failure but not a stored
  duplicate-endpoint. **Verify the ρ axis spacing once** (read one archive's
  `rho_axis`) and, if confirmed, the E8 signpost stands; if +180 *is* stored,
  this is a real interpolation question to escalate (not changed here).

- **Q3 — Internal-local rename breadth.** E7/E9/E11 rename file-local lambda
  and locals (`idx`→`tensor_offset`, `v`→`wrapped_value`, `v3`→`H_to_third`):
  zero cross-file cost, adopted. The `AxisLookup` struct field rename (C9:
  `idx/idx_next/frac` → `lower_index/upper_index/cell_fraction`) is file-local
  but read at three call sites in `QueryNearest`; the reviews disagreed (codex
  pro-rename, claude pro-keep). Plan **declines** it as churn without a clear
  clarity win, but flags it for the human as a weighed call.

- **Q4 — `has_acceptor_readouts` single-field proxy (CL28).** Sound because
  `ValidateSchema` enforces all-or-nothing on the acceptor set. Left without a
  comment to avoid hand-holding; mention here in case the human prefers a
  one-line `// safe: ValidateSchema enforces acceptor all-or-nothing` at
  `cpp:208`.
