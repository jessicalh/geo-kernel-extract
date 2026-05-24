# Fix plan — PiQuadrupoleResult.{h,cpp}

## 1. Summary

The file tells a coherent story. The header derivation (`.h:3-41`) is
exemplary — it names the model (axial point quadrupole on the ring normal),
gives the Stone T-tensor route to the closed-form `G_ab`, names every symbol,
states the three structural properties (symmetric / traceless / pure-T2), and
flags units. Both reviews agree it should stay. The `.cpp` follows that story:
the power ladder (`r2/r5/r7/r9`) is a deliberate precompute, and the fused
five-term tensor assignment (`.cpp:90-99`) sits directly under a comment block
(`.cpp:81-87`) that maps each displayed term to the code.

The two genuine soft spots both reviews converge on are local and cheap:
1. the `RingNeighbourhood` find-or-create block (`.cpp:181-205`) interleaves
   cross-calculator metadata bookkeeping and an *unlabelled* cylindrical-frame
   projection (`z`/`rho`/`theta`) into the middle of the quadrupole flow; and
2. a cluster of bare literals (`8`, `5`, `40`, `9`) carry array-layout/boundary
   meaning the reader must infer — most importantly the `8` at `.cpp:214`,
   which is the documented `kAromaticRingTypeCount` aromatic/saturated boundary.

This fix pass will: add 2–4-word signposts to the cylindrical-frame and
metadata blocks; replace the bare `8` ring-type ceiling with the existing
named constant `kAromaticRingTypeCount` (declined as risky elsewhere — see §3);
optionally introduce named intermediates for the five tensor terms (a weighed
call — see ledger); and tighten two near-duplicate comments. It will **not**
touch the algorithm, the numbers, the protected output names
(`pq_shielding` / `pq_per_type_T0` / `pq_per_type_T2`), or the cross-family
method name `SampleShieldingAt`.

The two correctness questions both reviews raise (filter divergence between
`Compute` and `SampleShieldingAt`; `ti < 8` dropping Pro from per-type sums but
not from the total) are both **coherent by design** — traced and explained in
§4. No bug-by-exhaustion candidates.

---

## 2. Review-finding ledger

### codex.md

| # | Finding | Disposition |
|---|---------|-------------|
| C1 | `.h:64` `SampleShieldingAt` comment "evaluate PQ EFG kernel" shifts story EFG→shielding; rename/recomment | **Adopted (comment only)** → §3 E1. Decline the rename (method name is a cross-family contract, §3/§4). |
| C2 | `.cpp:88-97` five-term formula is one fused assignment; split into named intermediates (`axial_dd_term`, `cross_nd_term`, …) | **Adopted, conditional** → §3 E2. Offered as a weighed option; claude (CL-mathstruct) judged the adjacent comment block already makes the fusion legible. Recommend the lighter variant. |
| C3 | `.cpp:130-134` filter comment too long/process-like; shorten to `// near-field exclusions` | **Declined.** The length is earned: it states *why* the topology check matters more for a quadrupole than a dipole (higher multipole → larger convergence radius). claude (CL-comments) explicitly calls this a "good, grounded signpost." Brief §13 protects earned rationale. |
| C4 | `.cpp:181-205` find/create `RingNeighbourhood` interrupts physics; add `// ring-neighbour metadata` | **Adopted** → §3 E3. Duplicate of CL1. |
| C5 | `.cpp:196-201` cylindrical coords computed without naming the frame; add `// ring-frame coordinates` | **Adopted** → §3 E4. Duplicate of CL1/CL-mathstruct. |
| C6 | `.cpp:212-218` per-type accumulation needs `// per-ring-type features` signpost | **Adopted** → §3 E5. (Existing comment is `// Per-type accumulation`; tighten wording.) |
| C7 | `.cpp:239-257` `SampleShieldingAt` uses a different filter set than `Compute`; add signpost/rationale | **Adopted (comment)** → §3 E6. Resolved as coherent (§4); add the one-line rationale the sibling HaighMallion already carries. Duplicate of C18/CL6a. |
| C8 | `.cpp:47` rename `G` → `efg_kernel` | **Declined.** `G` is the exact symbol the header derivation and the `.cpp:81` comment use for `G_ab`; the field comment already says "EFG geometric kernel." Renaming the struct field would *break* the visible math↔code correspondence the brief wants to preserve. |
| C9 | `.cpp:48` rename `scalar` → `buckingham_a_kernel` / `quad_field_scalar` | **Adopted (partial)** → §3 E7. Keep field name `scalar` (header calls it "the scalar"; low value to churn the struct field) but the *comment* should name it the Buckingham A-term kernel. The output it feeds is `pq_per_type_T0` (§4). |
| C10 | `.cpp:61` rename `d` → `atom_from_ring` | **Declined.** `d = r_atom - r_center` is defined verbatim in the header and `.cpp` comments and is the symbol in every displayed formula term. Renaming severs the formula↔code link. claude did not flag it. |
| C11 | `.cpp:74` rename `dn` → `height_along_normal` | **Adopted (partial)** → §3 E8. Rename to `d_dot_n` (claude CL2's suggestion). `dn` appears literally in the header formula (`dn^2`, `dn`), so keep the mapping obvious; `d_dot_n` reads as "d·n" and the existing inline comment "height above ring plane" stays. Prefer this over codex's `height_along_normal` to preserve the formula symbol. |
| C12 | `.cpp:88` rename `diag_term` → `delta_trace_term` / `kronecker_delta_term` | **Adopted** → §3 E9. `delta_term` (it's the `delta_ab(...)` piece in the formula). claude (CL-naming) judged `diag_term` "fine"; codex's point that it's a matrix-location name not a physical one is the stronger read. Cheap, local. |
| C13 | `.cpp:145` rename `total_pairs` → `accepted_pairs` | **Adopted** → §3 E10. Correct: it increments only after `AcceptAll` and `continue`-on-reject, so it counts accepted evaluations. The log string at `.cpp:230` already labels it `atom_ring_pairs=` (post-acceptance), so the rename aligns code with its own log. |
| C14 | `.cpp:148` rename `ca` → `atom` / `conf_atom` | **Adopted** → §3 E11. `conf_atom` (it's the `MutableAtomAt(ai)` conformation atom). `atom` collides conceptually with `atom_pos`; `conf_atom` disambiguates. Local mutable; no carry-through. |
| C15 | `.cpp:213` rename `ti` → `ring_type_index` | **Adopted** → §3 E12. Local index; clearer. Pairs with E13 (the `8` boundary). |
| C16 | `.cpp:25-44` `.cpp` comment block duplicates header almost verbatim; replace with short pointer | **Declined (kept), trim only** → §3 E13a. claude (CL-comments) wants a trim to avoid two copies drifting; codex wants near-removal. A reader of the `.cpp` benefits from the formula being adjacent to the implementation. Keep it but shorten the redundant *properties/derivation prose* to a one-line "see header for derivation," retaining the displayed `G_ab` formula. Weighed: legibility of the dense assignment wins over DRY. |
| C17 | `.cpp:78` shorten scalar comment to `// Buckingham A kernel` | **Adopted** → §3 E7 (folded). |
| C18 | `.cpp:169` `// ---- GeometryChoice: filter exclusion ----` is visually noisy; → `// record exclusion` | **Adopted** → §3 E14. The banner styling is noise; the GeometryChoice machinery is self-evident from the `choices.Record(...)` call. |
| C19 | `.cpp:181` `// Find or create RingNeighbourhood` mechanical; → `// ring-neighbour metadata` | **Duplicate** of C4. |
| C20 | `.cpp:81-87`,`:162`,`:207`,`:212`,`:220` comments noted good / good-if-split | **Acknowledged, no change** (`:162`,`:207`,`:220` are already good; `:212` covered by C6; `:81-87` covered by C2). |
| C21 | `.h:65`/`.cpp:239` rename `SampleShieldingAt` → `SamplePqEfgAt` etc. | **Declined.** This is a method name shared verbatim by 7 ring/field calculators (BiotSavart, HaighMallion, McConnell, RingSusceptibility, Dispersion, HBond, PiQuadrupole) and called generically by `ui/` and `h5-reader/`. Renaming one breaks the family symmetry the consumers depend on. The story-shift concern is addressed by the comment fix (E1). See §4 + §6 Q1. |
| C22 | `.cpp:54` rename `ComputePiQuadKernel` → `ComputePiQuadEfgKernel` | **Declined (low value).** The function returns *both* the EFG `G` and the Buckingham A `scalar`; "EFG" in the name would actually *narrow* it misleadingly. Current name is accurate-enough; claude calls it "clean." |
| C23 | `.cpp:261` rename `PackST_PQ` → `PackPiQuadSphericalTensor` | **Declined.** claude (CL4) notes the packing is *generic* (T0,T1,T2 → 9 doubles), so the `_PQ` suffix slightly over-claims but renaming to `PiQuad`-specific over-claims the other way. It's a file-local static used once. Leave it; net-zero clarity. |
| C24 | `.cpp:225-226` "Store accumulated shielding contribution (pure T2)" hides the decompose stage; → `// decompose total EFG` | **Adopted** → §3 E15. The line *is* a `SphericalTensor::Decompose(G_total)`; naming the stage helps. |
| C25 | `.cpp:130-138` vs `239-257` correctness: grid sampling ≠ atom filter set | **Duplicate** of C7 / CL6a. Resolved coherent in §4. |

### claude.md

| # | Finding | Disposition |
|---|---------|-------------|
| CL1 | `.cpp:181-205` find/create block interrupts story with unlabelled cylindrical setup; add `// shared ring-neighbour geometry (used by other calcs)` | **Adopted** → §3 E3+E4. The "used by other calcs" qualifier is the stronger framing (the `z`/`rho`/`theta` feed `RingNeighbourhood`, consumed by other ring calculators, not by PiQuadrupole's own tensor). Duplicate of C4/C5 but with better wording. |
| CL2 | `.cpp:74` `dn` → `d_dot_n` or `height` | **Adopted** → §3 E8 (`d_dot_n`). Duplicate of C11. |
| CL3 | `.cpp:88` `diag_term` acceptable as-is | **Noted.** Disagrees with C12; I side with C12 (adopt `delta_term`) — see ledger C12. The disagreement is cosmetic and the change is harmless. |
| CL4 | `.cpp:182,189-204` `rn`/`new_rn` terse but type is right there; acceptable | **Declined (no change).** Agreed acceptable. |
| CL5 | `.cpp:269` `N` for atom count fine | **No change.** Agreed. |
| CL6 | `.cpp:90-99` fused 5-term expr legible thanks to adjacent comment; clean | **Noted** — this is the counterweight to C2; informs the conditional disposition of E2 (recommend the light variant). |
| CL7 | `.cpp:69-72` power ladder grouped/ordered; clean | **No change.** Agreed. |
| CL8 | `.cpp:196-201` cylindrical projection is an unnamed sub-step; label it | **Adopted** → §3 E4. Duplicate of C5/CL1. |
| CL9 | `.cpp:201` `theta = atan2(d_plane.norm(), abs(z))` folds both half-spaces via `abs(z)`; intent (axial symmetry) invisible; note `// fold to [0,pi/2] (axial symmetry)` | **Adopted** → §3 E16. Good catch — the `abs(z)` is deliberate (the quadrupole is symmetric across the ring plane) but silent. One-line comment. |
| CL10 | `.cpp:214` `ti < 8` and `.cpp:216` `c < 5` are magic; name/comment them | **Adopted** → §3 E13 (use `kAromaticRingTypeCount` for the 8) and E17 (comment the 5 as T2-component count). The `8` is the documented aromatic/saturated boundary (§4). Duplicate-adjacent to CL12. |
| CL11 | `.cpp:261` `PackST_PQ` `_PQ` over-claims; minor, leave it | **Declined (leave).** Duplicate of C23; agreed leave. |
| CL12 | `.cpp:25-44` duplicated kernel/properties block; trim to short pointer | **Duplicate** of C16 → §3 E13a (trim, keep formula). |
| CL13 | `.cpp:130-134` good grounded signpost | **No change.** Agreed; counters C3 (decline). |
| CL14 | `.cpp:42` `// Traceless: yes (Laplace; verified numerically in tests)` good grounded signpost | **No change.** Agreed; keep. |
| CL15 (corr) | `.cpp:248-251` `SampleShieldingAt` guards (`distance < radius`) differ from `Compute`'s full filter set; confirm intended | **Adopted as resolved-coherent** → §4. Same as C7/C25. Add rationale comment (E6). |
| CL16 (corr) | `.cpp:64` singularity guard returns zeroed `distance`; caller `:164` reads `kernel.distance`(0.0) into filter ctx — check it can't pass `DipolarNearFieldFilter` and produce spurious zero accumulation | **Resolved-benign** → §4. Traced: `MinDistanceFilter` (first in set) rejects `distance=0.0` before `DipolarNearField` is reached; and `G` is zero regardless. No edit needed; documented in §4. |
| CL17 (corr) | `.cpp:214` `ti < 8` drops ring types ≥8 from per-type sums while `:221` adds them to `G_total`; if any type index >7 the per-type NPYs and total disagree — confirm 8 is the true ceiling | **Resolved-coherent** → §4. The disagreement is *intended*: index 8 is `ProPyrrolidine` (saturated, aromaticity=None), excluded from per-*aromatic*-type accumulation by the documented `kAromaticRingTypeCount` convention, but Pro still contributes to the physical total EFG. Edit E13 makes this self-evident. |
| CL18 (corr) | `.cpp:208-209` `rn->quad_tensor`/`quad_spherical` overwritten not accumulated; correct because one `rn` per ring index; note the find-or-create relies on uniqueness | **Resolved-coherent, no edit.** Confirmed: the find-or-create loop guarantees one `RingNeighbourhood` per `ring_index` per atom, so overwrite is the right semantics. No change. |

---

## 3. Edits that don't move numbers

All edits are comment/signpost/local-name changes. None alters control flow,
arithmetic, output names, or the public method signature.

- **E1** — `PiQuadrupoleResult.h:64` — change comment `// Grid sampling: evaluate PQ EFG kernel at arbitrary 3D point.` → `// Grid sampling: accumulated PQ EFG tensor at an arbitrary 3D point` (states it returns the *summed, decomposed* geometric tensor, not a final shielding). Keep the method name `SampleShieldingAt` (family contract, §6 Q1).
- **E2** *(conditional — weighed call; recommend light form)* — `PiQuadrupoleResult.cpp:90-99` — keep the single fused assignment (CL6 confirms the adjacent comment makes it legible) but renumber/align the comment block at `.cpp:81-87` so each displayed term lines up one-to-one with its code line, OR, if Jessica prefers codex's split, extract five named scalars per `(a,b)` (`axial_dd_term`, `cross_nd_term`, `dipole_dd_term`, `normal_normal_term`, `delta_term`) and sum them. Recommendation: **light form** (comment alignment only) — the full split adds 5 temporaries inside the 3×3 loop for no numeric or correctness gain and risks reviewer reading the split as a physics change.
- **E3** — `PiQuadrupoleResult.cpp:181` — replace `// Find or create RingNeighbourhood` → `// shared ring-neighbour geometry (consumed by other ring calculators)`.
- **E4** — `PiQuadrupoleResult.cpp:196` — add above the `d_vec`/`z`/`d_plane` lines: `// ring-frame cylindrical coordinates (z along normal, rho in plane)`.
- **E5** — `PiQuadrupoleResult.cpp:212` — change `// Per-type accumulation` → `// per-aromatic-ring-type feature bookkeeping (not part of the EFG sum)`.
- **E6** — `PiQuadrupoleResult.cpp:239` (above the loop) — add the sibling's rationale, matching `HaighMallionResult.cpp:362`: `// Same physics as Compute(); no atom-specific filters apply to a bare grid point — only the geometric guards (singularity / inside-radius / cutoff).`
- **E7** — `PiQuadrupoleResult.cpp:78` — change `// Scalar: (3 cos^2 theta - 1) / r^4` → `// Buckingham A-term kernel: (3 cos^2 theta - 1) / r^4  (feeds pq_per_type_T0)`. Keep struct field name `scalar`.
- **E8** — `PiQuadrupoleResult.cpp:74` — rename local `dn` → `d_dot_n` (and `dn2` → `d_dot_n_sq`) throughout `ComputePiQuadKernel`; keep the inline comment "height above ring plane." Local to the function; no carry-through.
- **E9** — `PiQuadrupoleResult.cpp:88` — rename local `diag_term` → `delta_term` (it is the `delta_ab(...)` piece). Local; update the one use at `.cpp:97`.
- **E10** — `PiQuadrupoleResult.cpp:145,222` — rename local `total_pairs` → `accepted_pairs`; update the log string at `.cpp:230` is unaffected (it already reads `atom_ring_pairs=`). Local.
- **E11** — `PiQuadrupoleResult.cpp:148` — rename local `ca` → `conf_atom` throughout `Compute` (uses at `.cpp:148,173,183,203,204,215,217`). Local.
- **E12** — `PiQuadrupoleResult.cpp:213` — rename local `ti` → `ring_type_index`. Local.
- **E13** — `PiQuadrupoleResult.cpp:214` — replace literal `8` with the existing named constant `kAromaticRingTypeCount` (from `Types.h`, already in scope via the include graph): `if (ring_type_index >= 0 && ring_type_index < kAromaticRingTypeCount)`. This is the documented convention and `Types.h:212-222` *explicitly* lists `PiQuadrupoleResult.cpp` as a deferred adoption site for exactly this swap — the constant exists and is static-asserted against the enum. **Carry-through note:** the matching `8`/`40` array dimensions are the NPY ABI shape (see E18); leave those literals or comment them — do not couple them to the constant without a schema migration (Types.h:217-222 warns of this). This edit changes no number: `kAromaticRingTypeCount == 8`.
- **E13a** — `PiQuadrupoleResult.cpp:25-44` — trim the prose duplicating the header (properties list, derivation sentence) to a single line `// see header for the Stone T-tensor derivation and properties`, but **retain** the displayed `G_ab` five-term formula block (`.cpp:32-36`) and the scalar line (`.cpp:38`) so the formula stays adjacent to the implementation.
- **E14** — `PiQuadrupoleResult.cpp:169` — replace `// ---- GeometryChoice: filter exclusion ----` → `// record the exclusion`.
- **E15** — `PiQuadrupoleResult.cpp:225` — change `// Store accumulated shielding contribution (pure T2)` → `// decompose the total EFG tensor (T0=T1=0; pure T2)`.
- **E16** — `PiQuadrupoleResult.cpp:201` — add inline: `// theta folded to [0, pi/2] via |z| — quadrupole is symmetric across the ring plane`.
- **E17** — `PiQuadrupoleResult.cpp:216` — comment the `5` once at the loop: `for (int c = 0; c < 5; ++c)  // 5 T2 components` (a named constant for T2 cardinality is not introduced; the `5` recurs in `PackST_PQ` and `WriteFeatures` as the SphericalTensor T2 layout — leave those as the established convention).
- **E18** *(documentation comment only, no rename)* — `PiQuadrupoleResult.cpp:271-273` (`WriteFeatures`) — add a one-line comment that `9` = packed SphericalTensor (T0,T1×3,T2×5), `8` = aromatic ring types, `40` = 8×5 (per-type × T2). This is the NPY ABI shape; the literals stay (changing them is a schema migration, Types.h:217-222).

No renames cross a file boundary. The protected serialized names
(`pq_shielding`, `pq_per_type_T0`, `pq_per_type_T2`) and the family method name
`SampleShieldingAt` are untouched.

---

## 4. Usage notes (the reasons behind the sign/value questions)

**(a) `ti < 8` drops Pro from per-type sums but not from the total EFG — coherent, by design.**
`ti = ring.TypeIndexAsInt()`. The `RingTypeIndex` enum (`Types.h:189-205`) is
nine values: indices 0–7 are aromatic chemistries (PHE, TYR, TRP6, TRP5, TRP9,
HIS, HID, HIE) and index 8 is `ProPyrrolidine`, a **saturated** ring with
aromaticity = None. `Types.h:183-241` documents the convention: per-aromatic-type
accumulator arrays are sized `kAromaticRingTypeCount` (= 8) and "gate Pro out via
the conventional `if (ti < kAromaticRingTypeCount)` guard," with static_asserts
pinning the boundary to the enum. So the guard *intentionally* excludes Pro from
`per_type_pq_scalar_sum`/`per_type_pq_T2_sum` (the per-aromatic-type features),
while line `.cpp:221` `G_total += kernel.G` still includes Pro's physical EFG
contribution to the atom's total `pq_shielding`. **Consumers agree:** the SDK
catalog sizes `pq_per_type_T0` at 8 and `pq_per_type_T2` at 40 (8×5)
(`_catalog.py:170-172`), and the per-atom storage `ConformationAtom.h:243-244`
is `std::array<…,8>`. The producer and every consumer share the 8-wide
per-aromatic shape; the total tensor is separately 9-wide (T0,T1,T2). Not a bug.
Edit E13 makes the boundary self-documenting by using the named constant.

**(b) `Compute` and `SampleShieldingAt` apply different near-field criteria — coherent, family convention.**
`Compute` runs the full atom-specific filter set (MinDistance + DipolarNearField
+ RingBondedExclusion) because it has an atom *index* and topology to exclude
bonded rings. `SampleShieldingAt` evaluates at a bare 3D grid point that has no
atom index and no bond graph, so atom-specific filters are inapplicable; it
keeps only the three geometric guards (singularity, inside-ring-radius, cutoff).
This is the **identical pattern** used by the sibling `HaighMallionResult.cpp:362`,
which carries the explicit comment "Same physics as Compute(). No atom-specific
filters (grid points)." `SampleShieldingAt` is consumed by `ui/` (ComputeWorker
volumetric field rendering) and conceptually mirrored in `h5-reader/`'s
`QtBiotSavartCalc`/`QtHaighMallionCalc`; all expect the grid-sampling path to
omit atom filters. Producer and consumers agree. Edit E6 adds the missing
rationale comment that the sibling already has. Not a bug.

**(c) `scalar` (A^-4) vs the `pq_per_type_T0` NPY (units `Angstrom^-5` in catalog).**
`scalar = (3cos²θ−1)/r⁴` is the Buckingham A-term kernel (header `.h:37-41`,
units A^-4). It accumulates into `per_type_pq_scalar_sum` (`ConformationAtom.h:243`,
comment "(3cos²θ-1)/r⁴ per ring type"), which `WriteFeatures` emits as
`pq_per_type_T0.npy`. The `_catalog.py:170` entry tags this array
`units="Angstrom^-5"`, but the underlying scalar is A^-4. This is a likely
**catalog metadata typo** (the `_T0` slot here is the Buckingham A scalar, not a
T0 component of the A^-5 EFG tensor — they share the `_T0` *name* across the
ring-kernel family but PQ's is a genuinely different physical quantity with
different units). **However:** this is a metadata-only field in the SDK, the
brief protects output *names* not metadata, and `learn/` calibration applies a
free per-feature coefficient that absorbs any unit scale — so no number moves
either way. I am flagging it as a §6 question rather than editing the catalog,
because changing SDK contract metadata is outside this file's fix scope and the
units string has no numeric consumer I could find. **Producer (A^-4) and the
calibration consumer agree on the values; only the SDK's descriptive units
string disagrees.**

**(d) Singularity-guard return with `distance = 0.0` — benign.**
If `ComputePiQuadKernel` hits `r < singularity_guard_distance` it returns the
default `PiQuadKernelResult` (G = 0, distance = 0). The caller writes
`ctx.distance = 0.0`. `MinDistanceFilter` is **first** in the filter set
(`.cpp:136`) and rejects on `0.0 >= singularity_guard_distance` → false, so the
pair is excluded before `DipolarNearFieldFilter` is evaluated; and `G` is zero
regardless, so even an accumulation would add nothing. No spurious contribution.
No edit needed.

**(e) `accepted_pairs` (was `total_pairs`).**
The counter increments only after `filters.AcceptAll(ctx)` passes (the reject
branch `continue`s first), so it counts accepted atom–ring evaluations, matching
the existing log label `atom_ring_pairs=` at `.cpp:230`. The rename (E10) aligns
the variable with what it already reports.

---

## 5. Bug-by-exhaustion candidates

None. Every sign/value/correctness item raised by the reviews resolved to
**coherent (expected)** under tracing (§4 a–e). The closest thing to a real
defect is the §4(c) catalog units string (A^-5 vs A^-4 for `pq_per_type_T0`),
which is (i) outside this file, (ii) metadata-only with no numeric consumer, and
(iii) absorbed by calibration regardless — so it is logged as a question (§6 Q2),
not a bug.

---

## 6. Questions & Ambiguities

- **Q1 — `SampleShieldingAt` family rename.** Both reviews want PQ's method
  renamed to signal "geometric EFG tensor, not final shielding"
  (codex `SamplePqEfgAt`, claude leans leave-it). The name is shared verbatim by
  7 ring/field calculators and called generically by `ui/ComputeWorker.cpp:427`
  and mirrored in `h5-reader/`. A rename is only sensible as a *coordinated
  family-wide* rename across all 7 + both consuming subprojects — out of scope
  for a single-file fix pass, and arguably out of scope for the whole readability
  effort given the extractor-untouchable posture. **Recommend: keep the name,
  fix the comment (E1). Confirm you don't want the family-wide rename queued
  separately.**

- **Q2 — `pq_per_type_T0` units string.** `_catalog.py:170` says
  `units="Angstrom^-5"` but the Buckingham A scalar it carries is A^-4
  (header `.h:41`, `ConformationAtom.h:243`). Sibling ring calcs reuse the `_T0`
  *name* for genuinely different per-type quantities, so this isn't a
  copy-paste-from-EFG so much as a units field that was set to match the EFG
  tensor's A^-5. Is the A^-5 a deliberate "report everything in the calculator's
  nominal unit" convention, or a metadata typo? Either way no number moves; I
  left the catalog untouched. Flag for a one-line metadata fix if it's a typo.

- **Q3 — E2 fused-term split.** codex wants the five tensor terms split into
  named intermediates; claude judges the adjacent comment already makes the
  fusion legible. I recommend the light variant (comment alignment, no split) to
  avoid adding loop-local temporaries that a reviewer could misread as a physics
  change — but this is a taste call. Confirm preference before the fix session.
