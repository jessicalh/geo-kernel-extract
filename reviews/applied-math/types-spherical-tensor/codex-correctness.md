# types-spherical-tensor — codex review

- **Targets:** src/Types.{h,cpp}
- **Model:** codex `gpt-5.5`, reasoning `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 15355
- **Brief:** `../_brief.md`

---

**1. Algorithm Correctness**
- `src/Types.cpp:118` — degeneracy guard assumes callers pass a strictly positive threshold; `0.0` lets zero T2 vectors reach `0/0` at line 121 — check call sites or guard with `<=` / validate threshold.
- `src/Types.cpp:121` — cosine can round slightly outside documented `[-1, 1]` — clamp final ratio or document as unclamped.
- `src/Types.h` — clean on concrete algorithm correctness findings.

**2. Variable Naming Clarity**
- `src/Types.cpp:25` — `s` hides that this is the input Cartesian tensor — rename to `sigma` or `tensor`.
- `src/Types.cpp:49` — `SQRT3_2` is ambiguous: `sqrt(3)/2` vs `sqrt(3/2)` — rename to `SQRT_THREE_HALVES`.
- `src/Types.cpp:65` — `INV_SQRT3_2` repeats that ambiguity — rename to `INV_SQRT_THREE_HALVES`.
- `src/Types.cpp:72` — `Sxx_m_Syy` uses cryptic `m` for minus — rename to `Sxx_minus_Syy`.
- `src/Types.h` — clean on variable naming clarity.

**3. Grouping Of Math Operations**
- `src/Types.cpp` — clean; decomposition/reconstruction steps are visibly grouped.
- `src/Types.h` — clean; no tangled math blocks.

**4. Method / Function Naming**
- `src/Types.h:243` — `RingTypeName` returns compact residue/ring codes, not full names — rename to `RingTypeCode`.
- `src/Types.cpp` — clean on method/function naming.

**5. Comments**
- `src/Types.h:207` — long process/ABI note is too verbose for core type comments — replace with `// NPY shape boundary`.
- `src/Types.h:224` — static-assert explanation is overlong — replace with `// boundary invariant`.
- `src/Types.cpp:31` — Levi-Civita comment is verbose but the block needs a label — replace with `// antisymmetric dual`.
- `src/Types.cpp:38` — formula comment restates the following assignments — replace with `// traceless projection`.
- `src/Types.cpp:52` — inline T2 comments repeat assignment formulas — replace with bare labels like `// m = -2`.
- `src/Types.cpp:74` — two-line diagonal derivation is longer than needed — replace with `// recover diagonals`.
- `src/Types.cpp:103` — T2 inner-product comment is explanatory/verbose — replace with `// T2 Frobenius dot`.
