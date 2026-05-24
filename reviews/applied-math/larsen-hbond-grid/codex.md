# larsen-hbond-grid — codex review (readability focus)

- **Targets:** src/LarsenHBondGrid.{h,cpp}
- **Model:** codex `gpt-5.5` `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 18,869
- **Brief:** `../_brief.md`

---

**Verdict**
`src/LarsenHBondGrid.h`: It tells a mostly coherent story, but the story is buried in a header-length design memo. A domain reader can eventually follow the contracts, frames, and archive mapping, but not in one pass without wading through history, pipeline notes, and future-facing comments.

`src/LarsenHBondGrid.cpp`: The implementation has a clear broad order: load archives, compute geometry/frame, query grid. The places that break readability are the compressed interpolation/indexing blocks and a few comments that preserve review/process history instead of explaining the physics or data contract.

**1. Coherent Story / Readability**
`src/LarsenHBondGrid.h:3` — the file opens with a large design document rather than a compact API contract. As written, a reader must separate runtime invariants from parser history, limitations, and future calculator notes. Suggest splitting mentally into shorter signposted blocks: `grid lookup contract`, `geometry contract`, `donor frame`, `readout mapping`, `validity mask`.

`src/LarsenHBondGrid.h:95` — “parser and the future calculator” makes the API read like work-in-progress history. Suggest: “parser and calculator use the same donor-frame convention.”

`src/LarsenHBondGrid.h:107` — “future session” breaks production-code narrative. Suggest removing the process note and saying only where dispatch occurs.

`src/LarsenHBondGrid.h:132` — tensor rotation contract is important and mostly clear, but it appears after a long validity-mask section. As written, a reader must remember the canonical-frame definition from 50 lines earlier. Suggest moving the rotation contract directly after canonical donor frame or adding a short `// tensor rotation` signpost near the helper declaration.

`src/LarsenHBondGrid.cpp:215` — interpolation is correct-looking but reads as index-symbol soup. As written, a reader must decode `c000..c111`, then infer interpolation order from variable suffixes. Suggest a terse block label plus names like `rho00`, `theta0`, `r_interp`, or a comment `// rho -> theta -> r`.

`src/LarsenHBondGrid.cpp:280` — `LookupAxis` combines three ideas: bounds tolerance, periodic wrapping, and lower-cell lookup. As written, a reader must hold the periodic and non-periodic contracts together. Suggest signposts: `// periodic wrap`, `// bounds clamp`, `// lower grid cell`.

**2. Naming Carries Meaning**
`src/LarsenHBondGrid.cpp:224` — `idx` is generic inside tensor interpolation. Suggest `tensor_offset`.

`src/LarsenHBondGrid.cpp:232` — `c000` through `c111` are standard but cryptic without axis order. Suggest either a comment `// corners: r, theta, rho` or names like `corner_r0_th0_rho0`.

`src/LarsenHBondGrid.cpp:290` — `LookupAxis` is clear enough for grid code, but the return fields `idx`, `idx_next`, `frac` lose quantity/frame. Suggest `lower_index`, `upper_index`, `cell_fraction`.

`src/LarsenHBondGrid.cpp:305` — `v` is only locally acceptable, but this is a public-domain math contract. Suggest `wrapped_value`.

`src/LarsenHBondGrid.cpp:380` — `b1`, `b2`, `b3`, `n1`, `n2`, `m1`, `x`, `y` are conventional but not self-carrying for a once-through reader. Suggest names like `H_to_O`, `O_to_C`, `C_to_third`, `plane_HOC`, `plane_OCX`, `dihedral_sin_term`, `dihedral_cos_term`.

`src/LarsenHBondGrid.cpp:413` — `v3` hides that it is the donor-third direction. Suggest `H_to_third` or `third_to_H` depending on intended direction.

**3. Visible Math Structure / Grouping**
`src/LarsenHBondGrid.cpp:361` — `ComputeLarsenHBondGeometry` has good coarse order: distance, angle, dihedral. The dihedral block would read more coherently with named plane normals and atan2 terms.

`src/LarsenHBondGrid.cpp:394` — `ComputeLarsenDonorFrame` is well grouped: build z, project x, build y, assemble R. This is one of the clearer mathematical narratives in the file.

`src/LarsenHBondGrid.cpp:496` — `QueryNearest` has a coherent query path: normalize diagnostics, reject non-finite input, choose archive, lookup axes, interpolate readouts, set flags. Clean on structure.

`src/LarsenHBondGrid.cpp:562` — the comment explains old behavior and why it changed. That interrupts the current computation story. Suggest replacing with a current-state label such as `// diagnostic coordinates`.

**4. Function / Method Naming**
`src/LarsenHBondGrid.h:359` — `QueryNearest` is misleading because the implementation does trilinear interpolation, not nearest-neighbor lookup. Suggest `QueryInterpolated` or `Interpolate`.

`src/LarsenHBondGrid.cpp:219` — `TrilinearMat3` says the operation and return type clearly. Clean.

`src/LarsenHBondGrid.cpp:259` — `AnyCornerImputed` is clear and domain-specific. Clean.

`src/LarsenHBondGrid.cpp:290` — `LookupAxis` is acceptable but slightly vague about returning a cell and fraction. Suggest `LookupAxisCell`.

`src/LarsenHBondGrid.cpp:347` — `WrapRho` is clear. Clean.

**5. Comments As Signposts**
`src/LarsenHBondGrid.h:2` — the opening comment is useful but too expansive for a header. Suggest converting long explanatory paragraphs into terse section labels plus keeping only contracts callers need.

`src/LarsenHBondGrid.cpp:114` — “codex M3 finding” is process/history prose and will age badly. Suggest `// required datasets`.

`src/LarsenHBondGrid.cpp:215` — the interpolation comment is too long and still leaves the axis order hard to see. Suggest `// trilinear tensor interpolation`.

`src/LarsenHBondGrid.cpp:256` — comment restates the 8-corner enumeration. Suggest `// imputed-corner check`.

`src/LarsenHBondGrid.cpp:371` — good signpost: `// theta at acceptor O`. Clean.

`src/LarsenHBondGrid.cpp:378` — useful signpost, but could be shorter. Suggest `// rho dihedral`.

`src/LarsenHBondGrid.cpp:401` — useful signpost but over-explains. Suggest `// donor z-axis`.

`src/LarsenHBondGrid.cpp:410` — useful, but long. Suggest `// donor x-axis projection`.

`src/LarsenHBondGrid.cpp:428` — good signpost because row/transpose direction is easy to get wrong. Clean.

`src/LarsenHBondGrid.cpp:562` — historical explanation. Suggest `// diagnostic coordinates`.

**6. Correctness**
`src/LarsenHBondGrid.h` — no concrete correctness bug found from the inlined text.

`src/LarsenHBondGrid.cpp` — no concrete correctness bug found from the inlined text. I would not claim a sign/convention bug in the dihedral or tensor rotation without checking against the parser fixtures and Larsen reference cases.
