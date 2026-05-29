# Fix plan — types-spherical-tensor (`src/Types.{h,cpp}`)

## 1. Summary

`Types.h` is a type-catalog header: enums (Element, Hybridisation, AtomRole,
BondOrder/Category, Ring* family, AminoAcid), the `SphericalTensor` math
struct, and `FieldValue`. `Types.cpp` carries the `SphericalTensor`
Decompose/Reconstruct/T2 math plus three amino-acid lookup helpers. Both
files mostly tell a coherent story. The math in `Types.cpp` is unusually
well-staged for tensor algebra and is *pinned by dedicated tests*
(`tests/test_spherical_tensor.cpp`) that verify the round-trip, the
isometric norm, the real-spherical-harmonic basis convention, and basis
orthogonality. Nothing in the math is in doubt.

The only real friction is **comment volume and one set of misleading inline
annotations**, both in two spots:

- `Types.h` `kAromaticRingTypeCount` / `ProPyrrolidine` — design-seam and
  deferred-adoption prose buries a one-line constant and its static_assert
  intent.
- `Types.cpp` `Reconstruct` lines 86–91 — inline `// S_xy + A_xy` style
  annotations that don't visually line up with the `+`/`-` signs in the code
  (they *do* reconcile once you fold in `A_xz = -T1[1]`, but the reader is
  forced to do that algebra).

**This pass will touch:** comment wording, two 2–4-word signposts, two local
variable names in `Types.cpp`, and (proposed, with carry-through) one
internal function rename. **It will not touch:** any number, sign, basis
coefficient, the T2 `m=-2..+2` layout, serialized field names, or the
`AminoAcid`/`RingTypeIndex`/`Element` enum values. No algorithm change. The
governing prior held throughout: every sign/value the reviews flagged traced
out as **coherent**; there are **no bug-by-exhaustion candidates**.

---

## 2. Review-finding ledger

Every finding from `codex.md` and `claude.md`, one row each. (No
`codex-correctness.md` present.)

### codex.md

| # | Finding (loc) | Disposition |
|---|---|---|
| C1 | `Types.h:207` `kAromaticRingTypeCount` block reads as migration prose | **adopted** → E1 |
| C2 | `Types.h:198` `ProPyrrolidine` inline comment too long | **adopted** → E2 |
| C3 | `Types.h:224` static_assert explanation forces future-scenario parse before the rule | **adopted** → E3 |
| C4 | `Types.h:304` `T2InnerProduct` comment repeats the preamble | **adopted (trim, keep the reason)** → E4 |
| C5 | `Types.cpp:46` add a one-line T2 basis-order signpost | **adopted** → E5 |
| C6 | `Types.h:243` `RingTypeName` returns codes not names; rename to `RingTypeCode`/`RingTypeLabel` | **adopted (proposed with carry-through)** → E6 |
| C7 | `Types.h:343` `source_index` doesn't say what it indexes | **declined** → see note below |
| C8 | `Types.cpp:25` param `s` too compressed; use `tensor`/`sigma` | **adopted** → E7 |
| C9 | `Types.cpp:39` `Sxx`/`Sxy` rely on the comment; consider `traceless_xx` | **declined** (this file *is* math-facing; the `Sxx` family is the file's settled vocabulary and is what the tests reference by eye — renaming adds verbosity without clarity; claude explicitly endorses keeping them) |
| C10 | `Types.cpp:25` Decompose is well-staged | not a finding (praise); no action |
| C11 | `Types.cpp:62` Reconstruct is coherent | not a finding (praise); no action |
| C12 | `Types.h:207` group ring boundary as constant + asserts + one ABI note | **duplicate** of C1 → E1 |
| C13 | `Types.h:291` Isotropic/Antisymmetric/TracelessSymmetric clear | not a finding; no action |
| C14 | `Types.h:296` `Decompose` could be `DecomposeTensor`, current acceptable | **declined** (author calls it acceptable; static method on `SphericalTensor` already disambiguates; rename is churn) |
| C15 | `Types.h:299` `Reconstruct` clear enough | not a finding; no action |
| C16 | `Types.h:322` `T2CosineWith` clean | not a finding; no action |
| C17 | `Types.h:243` `RingTypeName` is misleading | **duplicate** of C6 |
| C18 | `Types.h:212` verbose process comment → terse | **duplicate** of C1 → E1 |
| C19 | `Types.h:220` ABI note buried | **adopted (fold into E1 as the one retained ABI line)** → E1 |
| C20 | `Types.cpp:28` `// T0: isotropic` good signpost | keep; no action |
| C21 | `Types.cpp:31` good, slightly verbose → `// antisymmetric dual` | **declined** (the existing two lines give the explicit dual-vector formula `v_x=A_yz,…`, which claude C-K1 explicitly wants kept; shortening to "antisymmetric dual" drops the formula the reader needs) |
| C22 | `Types.cpp:38` good signpost → `// traceless symmetric part` | **declined** (already essentially this; current comment also states the formula, which is worth keeping; no edit) |
| C23 | `Types.cpp:46` good signpost → `// T2 projection` | **partially adopted** — folded into E5 (we add a basis-order signpost rather than shorten; the existing comment already says more than "T2 projection") |
| C24 | `Types.cpp:80` good signpost → `// Cartesian reconstruction` | **declined** (current comment is fine and slightly more informative; no edit) |
| C25 | `Types.cpp:104` repeats preamble → `// isometric T2 dot product` | **duplicate** of C4 → E4 |
| C26 | `Types.cpp:31`/`86` Levi-Civita sign convention — verify external consumers, not a claimed bug | **adopted as a Usage-note investigation** → §4 (resolved: coherent) |

Note on **C7 (`source_index`)**: traced — `FieldValue::source_index` has
**no consumers anywhere** in `src/`, `python/`, `learn/`, `h5-reader/`, or
`ui/` (only the struct definition exists; `ui/FieldOverlay.cpp`'s
`"FieldValue"` is an unrelated VTK array name). It is a currently-unused
field. Renaming or commenting an unused field that no one reads is churn
that adds nothing for a reader of the math; left as a Question (§6) for the
human, since the right move may be to remove rather than rename it.

### claude.md

| # | Finding (loc) | Disposition |
|---|---|---|
| K1 | `Types.h:206-241` `kAromaticRingTypeCount` block ~25 lines around a 1-line constant; compress | **duplicate** of C1 → E1 |
| K2 | `Types.h:310-321` `T2CosineWith` doc is an essay; keep the NaN sentence, move cos²/polarisation rationale out | **adopted** → E8 |
| K3 | `Types.cpp:25` Decompose clean | not a finding; no action |
| K4 | `Types.cpp:62-93` Reconstruct followable but inline annotations inconsistent with arithmetic | **adopted** → E9 |
| K5 | `Types.cpp:25` param `s` → `sigma`/`tensor` | **duplicate** of C8 → E7 |
| K6 | `Types.cpp:72,77,78` `Sxx_m_Syy` reads as a typo'd product; → `Sxx_minus_Syy` | **adopted** → E10 |
| K7 | other names clean (`T0/T1/T2`, Isotropic, etc.) | not a finding; no action |
| K8 | `Types.cpp:48-56` T2 assignments aligned w/ `// m=…` labels, clean | not a finding (overlaps E5 signpost); no action |
| K9 | `Types.cpp:74-78` tracelessness inversion explicitly named, clean | not a finding; no action |
| K10 | `Types.h:286-324` `SphericalTensor` grouping sensible | not a finding; no action |
| K11 | function/method names clean across both files | not a finding; no action |
| K12 | `Types.cpp:31-33` good Levi-Civita signpost, keep | keep; no action (and see C21 decline) |
| K13 | `Types.cpp:86-91` inline comments restate-and-contradict the code | **duplicate** of K4 → E9 |
| K14 | `Types.h:198-203` `ProPyrrolidine` 6-line trailing doc, trim | **duplicate** of C2 → E2 |
| K15 | `Types.h:168-176` `RingAromaticity` block grounded (Joule & Mills), keep | keep; no action |
| K16 | `Types.cpp:129-130` "why implemented here not header" note, keep | keep; no action |
| K17 | `Types.cpp:86-91` correctness: comments disagree with assignments, but math round-trips → **comment bug, not math bug** | **adopted** → E9 (same edit; see §4 for the trace) |
| K18 | `Types.cpp:145,149-157` `MSE`→`MET` — confirm intended (Se≠S) | **declined as an edit; resolved in §4** (intended project convention) |
| K19 | `Types.cpp:163,169` bound check `i<=20` admits the `Unknown`/`UNK` slot; asymmetric vs the `i<20` loop; harmless, worth a one-line note | **adopted** → E11 |
| K20 | `Types.cpp:118` `T2CosineWith` result can exceed [-1,1] by FP rounding; callers wanting strict cosine should clamp | **declined as a code change; recorded in §4 + §6** (no consumer asserts the range; not a readability edit) |
| K21 | `Types.h:30` vs comment line 27 — "5 elements" vs 6 enumerators incl. `Unknown`; no bug | **declined** (claude itself concludes "no bug"; comment is accurate for chemistry; no edit) |

---

## 3. Edits that don't move numbers

All edits are comment/signpost/local-name only, except **E6** (proposed
internal rename, with carry-through listed).

- **E1** — `Types.h:207-241` (the `kAromaticRingTypeCount` doc block):
  compress the ~25-line block to its invariant + the single ABI line + the
  static_assert intent. Keep: (a) the boundary meaning ("indices
  `0..kAromaticRingTypeCount-1` aromatic, the rest saturated; the
  `if (ti < kAromaticRingTypeCount)` guard gates Pro out of per-aromatic-type
  accumulation"); (b) one ABI sentence ("`ConformationAtom.h`'s
  `std::array<double,8>` and the NPY per-aromatic-type shape both bake in
  `8`; changing the boundary is a coordinated NPY schema migration"); (c) the
  static_assert intent ("the asserts pin the constant to the first saturated
  index so a re-order surfaces at the build line"). Drop the per-calculator
  deferred-adoption inventory (BiotSavart/HaighMallion/PiQuadrupole/Dispersion
  file list) — that is migration history, not type-model fact.

- **E2** — `Types.h:198-203` (`ProPyrrolidine` trailing doc): trim to one
  line, e.g. `///< Saturated 5-ring (Pro); Aromaticity=None, outside`
  `///< kAromaticRingTypeCount.` Drop the duplicated guard restatement (now
  in E1) and the inline Joule & Mills cite (already cited on
  `RingAromaticity`).

- **E3** — `Types.h:233-238` (first static_assert message): keep the message
  but lead with the rule, e.g. first sentence
  `"kAromaticRingTypeCount IS the index of the first saturated ring type."`
  then the one-line "re-order ⇒ update in lockstep" consequence. No logic
  change; message text only.

- **E4** — `Types.h:304-308` + `Types.cpp:103-107` (`T2InnerProduct`):
  shorten both to state the reason once: "T2-space dot product. Equals the
  Frobenius inner product of the symmetric-traceless matrices because the
  isometric normalization preserves the Frobenius norm." Keep the one-line
  reason; drop the cross-reference verbosity.

- **E5** — `Types.cpp:46` (before the five `T2[m]` assignments): add one
  signpost line naming the basis order, e.g.
  `// T2 basis, m = -2..+2 (xy, yz, zz, xz, xx-yy); isometric real-SH norm.`
  Keep the existing per-line `// m = …` comments — they pin each index.

- **E6 (proposed internal rename, with carry-through)** — `Types.h:243`
  `RingTypeName` → **`RingTypeCode`**. The returned strings are residue/ring
  *codes* (`"TRP6"`, `"TRP9"`, `"HID"`), not names; the current name misleads.
  **Carry-through:** internal-only, 3 call sites, all test/debug printf
  labels — no serialized contract:
  - `tests/test_batch_biot_savart_haigh_mallion.cpp:511`
  - `tests/test_traversal_dump.cpp:240`
  - `tests/test_batch_coulomb_ringchi.cpp:504`
  Cost is low (3 mechanical edits, all in `tests/`); recommend adopting.
  `RingTypeLabel` is an acceptable alternative if "code" reads too strongly
  as numeric — "code" is preferred since these are the same short tokens used
  as residue codes.

- **E7** — `Types.cpp:25` `Decompose(const Mat3& s)` → rename the parameter
  to `sigma` (matches the header preamble's `sigma` and the `s.trace()`/
  `s(i,j)` body, which then read as `sigma.trace()` / `sigma(i,j)`).
  Local-only; update all `s(…)`/`s.trace()` uses in the function body.

- **E8** — `Types.h:310-321` (`T2CosineWith` doc): keep the
  "returns NaN below `magnitude_threshold`" sentence and the one-line "signed,
  in [-1,1]" statement; condense the cos²/polarisation/"physically different
  tensor" rationale to a single clause or a pointer to the design doc. Do not
  drop the SIGNED note (it's load-bearing — see §4).

- **E9** — `Types.cpp:86-91` (`Reconstruct` inline annotations): replace the
  six misleading `// S_xy + A_xy, where A_xy = +v_z` style comments with
  annotations that match the code's literal signs. The clean form is one
  upstream line establishing the dual mapping
  (`// A_ij from T1: A_xy=T1[2], A_xz=-T1[1], A_yz=T1[0]; result(i,j)=S_ij+A_ij`)
  then either bare per-line signs or no per-line comment. The arithmetic is
  correct (see §4); this is a comment-conforms-to-code fix only.

- **E10** — `Types.cpp:72,77,78` (and the line-72 declaration): rename local
  `Sxx_m_Syy` → `Sxx_minus_Syy`. Removes the "m = product" misread. Local-only.

- **E11** — `Types.cpp:161-171` (the `i >= 0 && i <= 20` bound in
  `ThreeLetterCodeForAminoAcid`/`IsAromaticAminoAcid`): add a one-line comment
  that index 20 is the deliberate `Unknown`/`UNK` sentinel slot (rows
  `AA_CODES[20]="UNK"`, `AA_AROMATIC[20]=false` exist), so the inclusive
  `<= 20` here vs the `< 20` lookup loop is intentional, not an off-by-one.
  Comment only; no bound change.

---

## 4. Usage notes (the sign/value reasons discovered)

**(a) T2 layout / basis order — coherent.** `T2[0..4]` is the real spherical
harmonic order `m = -2, -1, 0, +1, +2`, i.e. `(√2·Sxy, √2·Syz, √(3/2)·Szz,
√2·Sxz, (Sxx−Syy)/√2)`. This is not folklore — `tests/test_spherical_tensor.cpp`
pins it four ways: `T2DiagonalMatrix`, `IsometricNormPreservation`,
`T2BasisMatchesY2m` (ratio to closed-form `Y²_m` constant across random unit
vectors), and `T2BasisOrthogonality` (each basis matrix → unit on exactly one
m). **Consumers all read it in the same `[0..4]` order**: the serialized
producers (`BsShieldingTimeSeriesTrajectoryResult.cpp:157-161`,
`HmShieldingTimeSeriesTrajectoryResult.cpp:121-125`,
`ApbsEfgTimeSeriesTrajectoryResult.cpp:185-189`, `ApbsFieldResult.cpp:355`,
`OrcaShieldingResult.cpp:238`, `TripeptideNeighborShieldingResult.cpp:394`).
The SDK contract
(`python/nmr_extract/_catalog.py`, `EFGTensor`/`PerRingTypeT2`/`*_T2` specs,
irreps `1x2e` / `0e+1e+2e`) carries the same 5-component T2. Producer and
consumers **agree**; the layout is a stable contract and must not move. The
E5 signpost only documents it.

**(b) Reconstruct antisymmetric signs — coherent (comment bug only).**
`Decompose` produces `T1[0]=A_yz`, `T1[1]=A_zx=−A_xz`, `T1[2]=A_xy` with
`A_ij=(s_ij−s_ji)/2`. `Reconstruct` rebuilds `result(i,j)=S_ij+A_ij`:
`result(0,1)=Sxy+T1[2]`, `result(0,2)=Sxz+A_xz=Sxz−T1[1]`,
`result(2,0)=Sxz+A_zx=Sxz+T1[1]`, etc. I worked all six off-diagonals by hand
and they are correct; `tests/test_spherical_tensor.cpp` `RoundtripAsymmetric`
(a fully non-symmetric matrix) confirms the round-trip to 1e-12.
**Consumers:** `Reconstruct` is called only by round-trip *tests*
(`test_spherical_tensor`, `test_apbs_*`) and `DemoResult.cpp:57` — there is no
production consumer that depends on the per-element sign beyond the
self-inverse property, which holds. So claude's K4/K17 is real but is purely
a **misleading-comment** issue: the inline `// S_xz + A_xz` text doesn't track
the literal `−`/`+` because `A_xz=−T1[1]`. E9 fixes the annotations to match.
The Levi-Civita dual convention (codex C26) is the standard
`v_x=A_yz, v_y=A_zx, v_z=A_xy`; it is internally consistent and no external
consumer reads the raw `T1[k]` with an opposite-handedness assumption (the
`T1` triple is serialized as `1e` and consumed as an opaque vector).

**(c) `T2CosineWith` SIGNED + NaN behavior — coherent.** The signed-cosine
choice is deliberate and consumed as such: `MopacVsFf14SbReconciliationTrajectoryResult.cpp:71`
takes the per-atom signed cosine via `T2CosineWith` (its own comment at
:144 names "Frobenius inner product / …"). A sign flip is a physically
different (opposite-polarisation) tensor, so signed is the meaningful
invariant — the doc's claim is right and load-bearing; E8 preserves it. The
NaN-below-threshold guard (`Types.cpp:118`) is a real contract; keep its
sentence. **FP overshoot of [-1,1] (K20):** possible for near-equal tensors,
but no consumer asserts a strict range (`MopacVsFf14Sb…` just stores it), so
this is not a correctness defect and not a readability edit — recorded in §6
in case a future consumer wants a clamp.

**(d) `MSE`→`MET` (K18) — coherent, project convention.**
`AminoAcidFromThreeLetterCode` maps selenomethionine `MSE`→`MET`. This is the
**same mapping the project's `NamingRegistry` makes** (`src/NamingRegistry.cpp:79`
`to_canonical_["MSE"]="MET"`, documented at `NamingRegistry.h:387`). The
identity/topology path treats MSE as MET deliberately (Se is handled
elsewhere as needed); the convention is consistent across the codebase. No
edit; not a silent substitution.

**(e) `i <= 20` inclusivity (K19) — coherent.** `AA_CODES`/`AA_AROMATIC` each
have 21 rows (index 20 = `UNK`/`false`), and `AminoAcid::Unknown` casts to 20.
So `ThreeLetterCodeForAminoAcid`/`IsAromaticAminoAcid` accepting `i==20`
returns the sentinel row by design; the lookup loop using `i < 20` is also
correct because it must not match the `UNK` row as a real amino acid. The
asymmetry is intentional. E11 adds the one-line note so a reader doesn't read
it as an off-by-one.

---

## 5. Bug-by-exhaustion candidates

**None.** Every sign/value the reviews raised traced to **coherent**:
- T2 basis order and signs — pinned by four tests, consumed uniformly,
  contract-stable.
- Reconstruct antisymmetric signs — hand-verified for all six off-diagonals
  and confirmed by `RoundtripAsymmetric`; the only defect is comment wording.
- Levi-Civita dual convention — standard and internally consistent; no
  consumer assumes the opposite handedness.
- `MSE`→`MET` — matches `NamingRegistry`, project-wide convention.
- `i<=20` inclusivity — backed by the 21-row sentinel tables.
- `T2CosineWith` FP overshoot — possible but un-asserted by any consumer; not
  a defect.

---

## 6. Questions & Ambiguities

1. **`FieldValue::source_index` (codex C7) is unused.** No reader exists in
   `src/`, `python/`, `learn/`, `h5-reader/`, or `ui/`. The rename codex
   suggests is moot until it has a consumer. Question for the human: is this
   a planned-but-unwired field (keep + comment its intended source) or dead
   (remove)? I did not edit it either way. (`FieldValue` as a whole — `tensor`
   + `spherical` + `source_calculator` + `source_index` — appears unconsumed;
   worth confirming whether the struct is live at all.)

2. **`T2CosineWith` range clamp (claude K20).** Should the return be clamped
   to [-1,1] to guard FP overshoot for near-equal tensors? No current consumer
   asserts the range, so I left it. If any downstream code (now or in the
   planned trajectory-scope cosine TR types) will assert strict cosine bounds,
   a clamp belongs at that consumer or here — a decision, not a readability
   fix.

3. **E6 rename target.** `RingTypeCode` vs `RingTypeLabel` — I recommend
   `RingTypeCode` (the tokens double as residue codes), but the three call
   sites are debug labels, so either reads fine. Confirm the preferred token
   before the rename lands.
