# types-spherical-tensor — claude review (readability focus)

- **Targets:** src/Types.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict

**Types.h** — Tells a coherent story. It is a declarations file, so the "through-line" is a catalogue of enums and one math struct; both read cleanly. A chemist can follow every enum and the `SphericalTensor` interface without decoding. The only friction is comment volume: several doc-comments (the `kAromaticRingTypeCount` block, the `T2CosineWith` block) are essay-length where a label would do, which buries the simple declarations underneath.

**Types.cpp** — Tells a coherent story, and unusually well for tensor algebra. The `Decompose`/`Reconstruct` pair reads top-to-bottom as build → trace → antisymmetric → symmetric → spherical-harmonic, with named Cartesian intermediates (`Sxx`, `Sxy`, …) instead of raw index soup. A physicist who knows irreducible-tensor decomposition could follow it once through. The one genuine readability snag is the inline `Reconstruct` comments that contradict the code beside them (see correctness).

---

### 1. Coherent story / readability

- Types.h:206-241 — the `kAromaticRingTypeCount` doc block is ~25 lines of design-seam prose around a one-line constant; a reader hunting the value must wade through deferred-adoption history — compress to 3–4 lines (what the boundary means + the static_assert intent), drop the per-file adoption inventory.
- Types.h:310-321 — the `T2CosineWith` doc is an essay on sign semantics and threshold calibration above a 4-line function; keep the "returns NaN below threshold" sentence, move the cos²/polarisation rationale to a design doc.
- Types.cpp:25 — clean; `Decompose` sequences trace → T1 → symmetric → T2 in the natural order, each block labeled.
- Types.cpp:62-93 — `Reconstruct` is followable but the reader must trust the inline `// S_xy + A_xy` annotations, which are inconsistent with the arithmetic (see Axis 6); as written, a reader cross-checking the comments against the code will be misled rather than helped.

### 2. Naming carries meaning

- Types.cpp:25 — parameter `s` for the input Mat3 is terser than the rest of the file's `Sxx`/`Szz` vocabulary; `sigma` (matching the header preamble's `sigma`) or `tensor` would read clearer.
- Types.cpp:72,77,78 — `Sxx_m_Syy` reads as a typo'd product on first glance; `Sxx_minus_Syy` removes the ambiguity at no cost.
- Otherwise clean: `T0/T1/T2`, `Isotropic`, `Antisymmetric`, `TracelessSymmetric`, the Cartesian `Sxy` family, and the element/ring/AA enum names all carry their quantity unambiguously to a domain reader.

### 3. Visible math structure (grouping)

- Types.cpp:48-56 — clean: constants named, then the five `T2[m]` assignments aligned with `// m = …` labels — the spherical-harmonic basis is visible at a glance.
- Types.cpp:74-78 — clean: the "recover Sxx, Syy from Szz and (Sxx−Syy) via tracelessness" step is the one non-obvious inversion and it is explicitly named.
- Types.h:286-324 — `SphericalTensor` groups data → accessors → factory/ops sensibly; no tangling.

### 4. Function / method naming

- Clean across both files. `Decompose`, `Reconstruct`, `T2Magnitude`, `T2InnerProduct`, `T2CosineWith`, `ElementFromSymbol`, `AtomicNumberForElement`, `IsAromaticAminoAcid` each say what they return. No vague or misleading method names.

### 5. Comments as signposts

- Types.cpp:31-33 — good signpost ("antisymmetric pseudovector via Levi-Civita mapping") with the dual-vector formula; keep.
- Types.cpp:86-91 — the inline `// S_xy + A_xy` style comments restate-and-contradict the code (lines 86/88/90 say `+`/`+`/`+` in prose while two of the three subtract or add opposite to the annotation); replace each with the bare offset sign actually used, or delete — a wrong signpost is worse than none.
- Types.h:198-203 — the `ProPyrrolidine` enumerator carries a 6-line trailing doc; trim to "Saturated 5-ring; Aromaticity=None, outside kAromaticRingTypeCount."
- Types.h:168-176 — `RingAromaticity` block is appropriately grounded (Joule & Mills citation) and not overlong; keep.
- Types.cpp:129-130 — useful "why implemented here not in header" note; keep.

### 6. Correctness (secondary)

- Types.cpp:86-91 — the inline comments disagree with the assignments. `result(0,1)=Sxy+T1[2]` is labeled `A_xy=+v_z` but `result(0,2)=Sxz-T1[1]` is labeled "S_xz + A_xz" while subtracting, and `result(2,0)=Sxz+T1[1]` is labeled "S_xz − A_xz" while adding. The arithmetic itself round-trips correctly against `Decompose` (T1[0]=A_yz, T1[1]=A_zx=−A_xz, T1[2]=A_xy), so this is a comment bug, not a math bug — fix the annotations to match.
- Types.cpp:145,149-157 — `AminoAcidFromThreeLetterCode` maps `MSE` (selenomethionine) → `MET`; confirm that is intended for the charge/topology path (Se≠S) rather than a silent substitution — check convention.
- Types.cpp:163,169 — bound check `i >= 0 && i <= 20` admits index 20, the `Unknown`/`UNK` slot; `AA_CODES[20]="UNK"` and `AA_AROMATIC[20]=false` exist so it is in-bounds, but the loop in `AminoAcidFromThreeLetterCode` only iterates `i < 20`, so the off-by-one inclusivity is asymmetric across the three functions — harmless given the sentinel rows, worth a one-line note that index 20 is the deliberate Unknown slot.
- Types.cpp:118 — `T2CosineWith` guards the divisor with `magnitude_threshold` but the result can still exceed [-1,1] by FP rounding for near-equal tensors; downstream callers treating it as a strict cosine should clamp — check whether any consumer asserts the range.
- Types.h:30 vs comment line 27 — enum has 5 elements plus `Unknown` (6 enumerators); comment says "the 5 elements" — accurate for chemistry, but `SymbolForElement`/`AtomicNumberForElement` fold `Unknown` into the `default` returning `"?"`/`0`, fine; no bug.
