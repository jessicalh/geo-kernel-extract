# hbond — claude review (readability focus)

- **Targets:** src/HBondResult.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict

**HBondResult.h** — Clear. The header tells the whole story: what the kernel is, where the H-bond comes from (DSSP/Kabsch-Sander), the closed-form tensor with its symmetry properties, and the d/h_hat conventions. A chemist could read this and know what to expect.

**HBondResult.cpp** — Mostly readable, with a strong scaffold of banner comments that signpost the four stages (resolve → filter → accumulate → store). The math core (`ComputeHBondKernel`) is clean. The damage is concentrated in `Compute`'s Step-1 resolution loop, which is two near-identical ~50-line blocks (donor side / acceptor side) that a reader must diff in their head to confirm only the N/O roles swap. The story holds, but the duplication and the inline GeometryChoice-recording lambdas bury the actual resolution logic.

---

## 1. Coherent story / readability (primary)

- HBondResult.cpp:148-259 — The donor-side loop (153-204) and acceptor-side loop (207-258) are ~50 lines of near-identical code differing only in which residue plays N vs O; as written, a reader must diff the two blocks line-by-line to confirm the only difference is the donor/acceptor role swap — add a one-line signpost at each loop head: `// donor side: this residue's N → acceptor O` / `// acceptor side: donor N → this residue's O` so the symmetry is asserted, not inferred.
- HBondResult.cpp:160-171, 182-192, 214-225, 236-246 — The skip-reason GeometryChoice `Record` lambda (4-8 lines each, four times) interleaves diagnostics with the rejection control flow, so the actual "why we skip" (seq_sep too small / distance out of range) is visually subordinate to telemetry — keep the `if (reject) { ...; continue; }` skeleton legible; the lambda body is fine but reads as the point of the block when it is the footnote.
- HBondResult.cpp:89 vs 93-97 — `result.f` (the `(3cos²θ−1)/r³` axial scalar) and `result.M_over_r3` (the full tensor) are both computed but `f` is never consumed anywhere in this file; a reader wonders whether they are meant to be consistent — confirm `f` is read by a sibling/consumer, else flag as dead.
- HBondResult.cpp:374-387 — `nearest_kernel` recomputes `ComputeHBondKernel(atom_pos, nearest_hb.midpoint, nearest_hb.h_hat)` that was already computed inside the loop at 316; harmless but a reader pauses to check it's the same call — one-line comment `// recompute kernel for the nearest H-bond only` would close the question.

## 2. Naming carries meaning

- HBondResult.cpp:67 — `double f` on `HBondKernelResult` is opaque; it is the axial dipolar scalar `(3cos²θ−1)/r³` — rename `axial_scalar` or comment inline.
- HBondResult.cpp:311, 365, 371 — `count_3_5` / `hbond_count_within_3_5A` hardcode "3.5" in the name while the threshold is config-driven (`hbond_counting_radius`, line 352); if the config value ever differs from 3.5 the name lies — comment the name's "3.5" as the default value, or rename `count_within_radius`.
- HBondResult.cpp:146, 163 — `resolution_key` is a monotonic counter used only as a GeometryChoice slot id; the name suggests a lookup key — `choice_record_seq` or a comment.
- HBondResult.cpp:428, 443 — `PackST_HB` / the `_HB` suffix is a file-local de-collision hack; fine as a static, but the name says nothing about layout (T0|T1×3|T2×5) — one-word comment `// layout: T0, T1[3], T2[5]` (the magic offsets 1+i and 4+i already imply it).
- HBondResult.h:46 / .cpp:408 — `SampleShieldingAt` is clear; `point` is fine.

## 3. Visible math structure (grouping)

- HBondResult.cpp:91-98 — The full `M_ab` assembly is one fused expression spanning the three physical terms (9cosθ d̂h, −3hh, −(3d̂d̂−δ)) across the 3×3 loop; it matches the header formula exactly so it is decodable, but the three terms read as one wall — acceptable as-is given the header banner at 55-63 names them; no change needed.
- HBondResult.cpp:308, 361, 395 — `M_total` accumulate → `Decompose` is the clean payoff; the build→accumulate→decompose shape is visible. Good.
- HBondResult.cpp:194-203, 248-257 — The two `ResolvedHBond hb; hb.x = ...; hbonds.push_back(hb)` blocks are clean and identical-shaped; the field-by-field fill reads fine.

## 4. Function / method naming

- HBondResult.cpp:72 `ComputeHBondKernel` — Says what it returns (the kernel struct), clear.
- HBondResult.cpp:109 `HBondResult::Compute` — Standard factory name, fine.
- HBondResult.h:49 `HBondCount` — Returns `hbond_midpoints_.size()`; clear and matches the "for diagnostics" comment.
- Clean on this axis.

## 5. Comments as signposts

- HBondResult.cpp:26-41, 55-63, 105-107, 128-140, 276-296 — The banner/step comments are exactly the right kind: they name the stage and state the convention (donor N → acceptor O, the filter rationale). Keep.
- HBondResult.cpp:117-119 — Good: explains *why* SpatialIndexResult is a dependency but `(void)`-discarded. Exactly the non-obvious thing a reader needs.
- HBondResult.cpp:162, 183, 216, 237 — `// ---- GeometryChoice: ... ----` labels are repeated four times verbatim for the rejection records; fine as signposts but consider one says "seq_sep reject" and the other "distance reject" to distinguish the two rejection causes at a glance.
- HBondResult.cpp:351 — `// Count H-bonds within 3.5A of this atom` restates code but with a stale literal (threshold is config `hbond_counting_radius`, not hardwired 3.5) — change to `// count H-bonds within counting radius`.
- HBondResult.cpp:360 — `// Accumulate tensor from all H-bonds (1/r³ decay handles range)` — good, justifies the no-cutoff accumulation.
- HBondResult.cpp:378 — `// all DSSP H-bonds are backbone` — good, justifies the unconditional `true`.

## 6. Correctness (secondary)

- HBondResult.cpp:355 — `nearest_dist` is initialized to `NO_DATA_SENTINEL` (line 309); the `kernel.distance < nearest_dist` comparison only works if `NO_DATA_SENTINEL` is a large positive number — check `NO_DATA_SENTINEL` is `+inf`/large-positive, not 0/−1, or no nearest is ever recorded.
- HBondResult.cpp:144-145, 173-175, 227-229 — Dedup `seen` keys on `(donor_N, acceptor_O)` atom pairs; the donor-side and acceptor-side loops can both resolve the same physical bond from opposite ends, and the shared `seen` set correctly suppresses the duplicate — looks right; the asymmetric key ordering (N first, O second) is consistent across both loops so no aliasing. Good.
- HBondResult.cpp:419 vs 334/288 — In `Compute`, the near-field exclusion runs through `DipolarNearFieldFilter` (line 288); in `SampleShieldingAt` it is reimplemented inline as `distance < near_field_exclusion_ratio * hbond_distances_[hi]` (line 419) — confirm the inline ratio test matches the filter's criterion, else grid samples and atom evaluations use different near-field guards.
- HBondResult.cpp:177-180, 199-200 — `d = O_pos - N_pos`, `h_hat = d/dist`, midpoint = ½(N+O); consistent with header convention "donor N → acceptor O". Sign/direction consistent. Good.
- HBondResult.cpp:386 — `hbond_inv_d3 = 1/(nearest_dist³)` uses `nearest_dist` (the kernel distance from atom to midpoint), while the header doc and `hbond_distances_` track the N…O distance; these are different lengths — confirm `inv_d3` is intended as atom-to-midpoint⁻³, not N…O⁻³.
