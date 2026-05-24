# types-spherical-tensor — codex review (readability focus)

- **Targets:** src/Types.{h,cpp}
- **Model:** codex `gpt-5.5`, reasoning `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 11,088
- **Brief:** `../_brief.md` (readability-first)
- **See also:** `codex-correctness.md` (earlier correctness-first pass)

---

**Verdict**

`src/Types.h` mostly tells a coherent type-catalog story until the ring-type boundary section. That passage becomes process/history-heavy and forces a reader to separate the actual invariant from deferred migration notes and ABI concerns. `SphericalTensor` is readable, though `T0/T1/T2` rely on the preamble to carry meaning.

`src/Types.cpp` tells a much clearer story: decompose tensor, reconstruct tensor, compare T2, then amino-acid lookup helpers. The math blocks are mostly ordered well, but several comments are longer than needed and a few names make the reader remember tensor conventions instead of reading them directly.

**1. Coherent Story / Readability**

`src/Types.h:207` — the `kAromaticRingTypeCount` block stops reading like a type definition and turns into migration/process prose: “Adoption in calculator code is deferred… Bundle C / Slice A…” — suggestion: reduce this to the invariant and the ABI warning, e.g. `// aromatic/saturated boundary` plus one short note that calculator arrays currently assume eight aromatic ring types.

`src/Types.h:198` — the `ProPyrrolidine` inline comment is much longer than the enum value it explains — suggestion: keep only the facts a reader needs at the enum site: `///< Saturated; excluded from aromatic accumulators.`

`src/Types.h:224` — the static-assert explanation requires the reader to parse future-change scenarios before seeing the simple rule — suggestion: state the rule once: `// Pro is the first saturated ring type; this pins the aromatic boundary.`

`src/Types.h:304` — `T2InnerProduct` comment explains the same norm-preservation point already introduced in the `SphericalTensor` preamble — suggestion: shorten to `// T2-space dot product; isometric normalization preserves Frobenius inner product.`

`src/Types.cpp:46` — the T2 basis mapping is readable, but the basis order only appears inside per-line comments — suggestion: add one signpost before assignments: `// T2 basis order: xy, yz, zz, xz, xx-yy`.

As written, a reader must distinguish physical invariants from release/process notes in the ring-boundary section.

**2. Naming Carries Meaning**

`src/Types.h:243` — `RingTypeName` returns compact labels like `"TRP6"` and `"TRP9"`, not full names — suggestion: rename to `RingTypeCode` or `RingTypeLabel` if call sites allow.

`src/Types.h:343` — `source_index` does not say what collection it indexes — suggestion: `source_record_index`, `calculator_source_index`, or a nearby comment naming the indexed source.

`src/Types.cpp:25` — parameter `s` is mathematically conventional but too compressed for a public decomposition routine — suggestion: `tensor` or `sigma`.

`src/Types.cpp:39` — `Sxx`, `Sxy`, etc. are acceptable tensor notation, but they rely on the preceding comment for meaning — suggestion: leave them if this file is math-facing; otherwise use `traceless_xx`, `traceless_xy`, etc.

**3. Visible Math Structure**

`src/Types.cpp:25` — `Decompose` is well staged: isotropic, antisymmetric, traceless symmetric, T2 projection. Clean.

`src/Types.cpp:62` — `Reconstruct` is also coherent: invert T2, recover diagonal terms, reassemble matrix. Clean.

`src/Types.h:207` — the ring-boundary structure hides the actual data model under too much prose — suggestion: group as: boundary constant, invariant asserts, one ABI note.

**4. Function / Method Naming**

`src/Types.h:291` — `Isotropic`, `Antisymmetric`, and `TracelessSymmetric` are clear.

`src/Types.h:296` — `Decompose` is clear in context, but `DecomposeTensor` would be clearer outside `SphericalTensor`; current name is acceptable.

`src/Types.h:299` — `Reconstruct` is clear enough because the receiver is `SphericalTensor`.

`src/Types.h:322` — `T2CosineWith` clearly says what it computes. Clean.

`src/Types.h:243` — `RingTypeName` is the only misleading function name because returned values are residue/ring codes, not names.

**5. Comments As Signposts**

`src/Types.h:212` — verbose/process comment — terse replacement: `// calculator adoption pending`.

`src/Types.h:220` — useful but buried ABI note — terse replacement: `// NPY ABI: eight aromatic entries`.

`src/Types.cpp:28` — `// T0: isotropic = trace / 3` is a good signpost.

`src/Types.cpp:31` — good signpost, slightly verbose — replacement: `// antisymmetric dual`.

`src/Types.cpp:38` — good signpost — replacement: `// traceless symmetric part`.

`src/Types.cpp:46` — good signpost — replacement: `// T2 projection`.

`src/Types.cpp:80` — good signpost — replacement: `// Cartesian reconstruction`.

`src/Types.cpp:104` — longer than needed and repeats the preamble — replacement: `// isometric T2 dot product`.

**6. Correctness**

No definite correctness bug found from the inlined text.

One thing to verify, not a claimed bug: `src/Types.cpp:31` and `src/Types.cpp:86` use a specific Levi-Civita sign convention for mapping antisymmetric matrices to `T1`. The decomposition and reconstruction are internally consistent; only check this if external consumers expect the opposite dual-vector convention.
