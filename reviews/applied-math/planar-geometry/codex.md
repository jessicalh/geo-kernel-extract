# planar-geometry — codex review (readability focus)

- **Targets:** src/PlanarGeometryResult.{h,cpp}
- **Model:** codex `gpt-5.5`, reasoning `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 18,867
- **Brief:** `../_brief.md` (readability-first)
- **See also:** `codex-correctness.md` (earlier correctness-first pass)

---

**Verdict**

`src/PlanarGeometryResult.h`: The quantities are scientifically coherent, but the header does not tell a clean first-read story. It opens with cross-result wiring, dates, amendment history, and process notes before the reader gets the simple contract: this result stores per-frame planar-geometry observables.

`src/PlanarGeometryResult.cpp`: The main `Compute()` function has a good numbered through-line: pyramidalization → omega → aromatic chi2 → pucker → write/log. The story breaks mostly in the helper math and comments, where internal history and compressed symbols make the reader separate the computation from project archaeology.

**1. Coherent Story / Readability**

`src/PlanarGeometryResult.h:6` — `"CROSS-RESULT READ (writer side, 2026-05-19, PATTERNS §17)"` makes the reader start with dependency-policy internals instead of the physical result — replace with a short “optional downstream readers” note or move this below the scientific contract.

`src/PlanarGeometryResult.h:16` — `"landed in the 2026-05-05 → 2026-05-08 topology slice"` interrupts the useful classification/deviation distinction — keep only: substrate stores classification; this result stores per-frame deviations/observables.

`src/PlanarGeometryResult.h:27` — the omega convention is useful, but the old-doc correction at lines 36–39 forces the reader to distinguish current behavior from historical behavior — keep the current convention and X-Pro mask behavior only.

`src/PlanarGeometryResult.cpp:107` — the Cremer-Pople block mixes derivation, orientation warning, and review history — keep the equations and normal-orientation convention, but move `"caught in adversarial review 2026-05-09"` out of the code comment.

`src/PlanarGeometryResult.cpp:153` — the degeneracy guard explanation is scientifically useful but too long for the local flow — replace with a short label like `// phase degeneracy guard` plus one sentence: theta is undefined when Q is effectively zero.

`src/PlanarGeometryResult.cpp:246` — `Compute()` is well ordered, but the section comments are much longer than the code path they introduce, especially omega lines 278–309 — reduce to signposts: `// peptide-bond omega`, `// bond-graph successor`, `// missing atom guard`.

**2. Naming Carries Meaning**

`src/PlanarGeometryResult.cpp:70` — `A, B, C, D, G, n` hide central atom, neighbours, centroid, and plane normal — use names like `central_pos`, `nbr0_pos`, `neighbour_centroid`, `plane_normal`.

`src/PlanarGeometryResult.cpp:121` — `R1` and `R2` require the reader to remember the equations above — use `sin_weighted_basis` and `cos_weighted_basis`.

`src/PlanarGeometryResult.cpp:142` — `cs` and `sn` are the worst symbol compression in the file; the negative sine convention is hidden in `sn += -z_j * ...` — use `cos_projection_sum` and `neg_sin_projection_sum`.

`src/PlanarGeometryResult.cpp:150` — `Qcos` / `Qsin` are understandable only after reading the equations — `q2_cos_theta` and `q2_sin_theta` would carry the physical meaning.

`src/PlanarGeometryResult.cpp:348` — `ri` means aromatic ring index here, after meaning residue index in the omega loop — rename to `arom_ring_i` or `ring_i`.

`src/PlanarGeometryResult.cpp:377` — `ri` again means saturated ring index — rename to `sat_ring_i`.

`src/PlanarGeometryResult.cpp:264` — `nb` is too compressed for a sign-sensitive neighbour ordering — use `neighbour_atoms`.

`src/PlanarGeometryResult.cpp:272` — `pyr` hides that this is a signed displacement in Å — use `signed_pyramidalization_A` or `out_of_plane_A`.

**3. Visible Math Structure**

`src/PlanarGeometryResult.h` — clean on this axis; it declares stored quantities rather than implementing math.

`src/PlanarGeometryResult.cpp:35` — the dihedral helper fuses bond vectors, plane normals, and signed `atan2` terms into one compact block — split with named intermediates such as `normal_123`, `normal_234`, `cos_term`, `sin_term`.

`src/PlanarGeometryResult.cpp:134` — the pucker projection is conceptually clear, but the actual accumulation at lines 142–151 collapses “z displacement”, “m=2 phase”, and “signed sine convention” — introduce named projection sums before scaling.

`src/PlanarGeometryResult.cpp:330` — omega and omega deviation are computed correctly in sequence, but `WrapPi(omega - M_PI)` would read better as a named intermediate like `deviation_from_trans`.

`src/PlanarGeometryResult.cpp:246` — the main computation stages are visibly grouped and sequenced well.

**4. Function / Method Naming**

`src/PlanarGeometryResult.cpp:35` — `Dihedral` does not say signed/radians/range — suggest `SignedDihedralRadians`.

`src/PlanarGeometryResult.cpp:70` — `Pyramidalization` sounds like a general concept, but the return is signed out-of-plane displacement in Å — suggest `SignedOutOfPlaneDisplacement`.

`src/PlanarGeometryResult.cpp:100` — `CremerPople5Ring` is acceptable, but `ComputeCremerPople5RingPucker` would better say it returns Q/theta.

`src/PlanarGeometryResult.cpp:177` — `ThreeBondedNeighbours` hides the “exactly three” and “sorted atom indices” behavior — suggest `SortedThreeBondedNeighbourAtoms`.

`src/PlanarGeometryResult.h:103` — getter names are mostly clear; the private comments carry units/frame well enough.

**5. Comments As Signposts**

`src/PlanarGeometryResult.h:2` — the file header is more design history than signpost — replace with short labels: `// emitted quantities`, `// omega convention`, `// pyramidalization`, `// aromatic chi2`, `// Cremer-Pople pucker`.

`src/PlanarGeometryResult.cpp:49` — WrapPi’s comment is mostly history and cross-file process — replace with `// wrap to [-pi, pi]`.

`src/PlanarGeometryResult.cpp:66` — `"For a perfectly planar sp2 site, A == G"` is not a reliable explanation; planar means A lies in the neighbour plane, not necessarily at the neighbour centroid — replace with `// central-atom displacement`.

`src/PlanarGeometryResult.cpp:187` — the neighbour-sort comment is useful but too process-heavy with `"720-protein fleet"` — replace with `// stable sign order` plus one direct sentence.

`src/PlanarGeometryResult.cpp:287` — references to `"feedback_huxley_data_discipline"` and `PATTERNS.md` break the scientific through-line — keep the observable rule, move policy references to docs.

`src/PlanarGeometryResult.cpp:465` — the X-Pro write comment is longer than needed after the compute block already explained it — replace with `// X-Pro omega mask`.

**6. Correctness**

No definite algorithmic correctness bug is visible from the inlined code alone.

`src/PlanarGeometryResult.cpp:67` — comment correctness issue: `A == G` for planar sp2 is not generally true — say A lies in the neighbour plane, so `(A - G).dot(n_hat)` is zero.

`src/PlanarGeometryResult.cpp:177` — this file uses `std::array` but the inlined includes do not include `<array>` — add the direct include if it is currently coming transitively.

`src/PlanarGeometryResult.cpp:193` — this file uses `std::sort` but the inlined includes do not include `<algorithm>` — add the direct include for self-contained compilation.

`src/PlanarGeometryResult.cpp:390` — `pucker_valid` counts finite `Q` even when `theta` is NaN under the sub-amplitude guard — if the log is meant to mean fully valid pucker, count finite `Q` and finite `theta`, or rename the diagnostic to finite-Q count.
