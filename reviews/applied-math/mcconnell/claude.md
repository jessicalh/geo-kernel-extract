# mcconnell — claude review (readability focus)

- **Targets:** src/McConnellResult.{h,cpp} + src/MopacMcConnellResult.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict

**McConnellResult.h** — Tells a coherent story. The header comment lays out the full tensor formula, the three-term decomposition (T1/T0/T2 attribution), the K kernel, and the scalar f, and explicitly warns that f is not the trace. A chemist reads this once and knows what the file computes and why each output exists. Clean.

**McConnellResult.cpp** — Mostly readable; the kernel function reads top-to-bottom as a clean transcription of the documented formula. The `Compute` loop is long but follows an honest setup → per-bond accumulate → per-atom finalize arc. The main friction is the per-category T2 finalization, where the same symmetrize/de-trace idiom is open-coded three times alongside an unused lambda, and a few one-letter spatial names (`d`, `r`) that the comment header carries but the code does not.

**MopacMcConnellResult.h** — Coherent and self-explaining. The header states the bond-order weighting, that it is a measured QM quantity not a parameter, and the physical reading ("a C=O with order 1.8 contributes more"). A reader understands the delta from the plain McConnell variant immediately.

**MopacMcConnellResult.cpp** — Reads as a faithful clone of McConnellResult.cpp with the bond-order weight threaded in; the `weighted_M` / `weighted_f` named intermediates make the modulation visible at the accumulation site. Same minor frictions as the plain variant (open-coded T2 idiom, one-letter spatial names). The near-zero-bo skip is clear.

---

## 1. Coherent story / readability (primary)

- McConnellResult.cpp:269-285 — `project_traceless` lambda is defined then never called; the three category blocks open-code the same symmetrize+de-trace inline instead, so a reader trusts a helper that does nothing. — Delete the dead lambda (or call it); reader must currently verify the inline math matches the unused helper.
- McConnellResult.cpp:275-285 — The symmetrize-then-de-trace idiom (`0.5*(M+Mᵀ)`, subtract `trace/3·I`) is repeated verbatim three times with only the variable name changing; the shared "extract symmetric traceless part" step is visible only by reading all three. — Keep a single `// symmetrize + de-trace` signpost above the block (the comment at :274 says only "Extract symmetric part," omitting the de-trace).
- McConnellResult.cpp:257-265 — Nearest-CO finalization mixes a direction-renormalization guard, a T2 decompose, and the CN T2 in three asymmetric branches (CO renormalizes direction, CN does not); the asymmetry is real but unsignposted. — One label `// nearest CO/CN: direction + T2`.
- MopacMcConnellResult.cpp:153 — Three statements on one line (`{ zero_bo_skipped++; zero_bo_this_atom++; continue; }`) compress the noise-floor skip; reader must parse a packed line to see it is just a guard. — Split onto separate lines.
- McConnellResult.cpp:64, MopacMcConnellResult.cpp:57 — Early `return result;` returns a zero-initialized kernel on the singularity guard, so `distance` stays 0.0 and the bond silently contributes nothing downstream; correct, but the reader only learns the sentinel meaning by tracing the struct defaults. — `// too close: return zero kernel (distance stays 0)`.

## 2. Naming carries meaning

- McConnellResult.cpp:61-62, MopacMcConnellResult.cpp:54-55 — `d` and `r` are the midpoint→atom displacement and its length; the header comment names them but the code does not. — `displacement` / `r` is borderline-fine, but `d` alone forces the reader to hold the header in mind. Consider `disp`/`dist`.
- McConnellResult.cpp:75, .cpp:37 — `f` (McConnell scalar) is a true one-letter name; the field comment carries it but call sites (`co_sum += kernel.f`) read opaquely. — `f` is the literature symbol so defensible; leave, but the struct field could be `f_scalar`.
- McConnellResult.cpp:185 — `bn` and `BondNeighbourhood` are clear; no action.
- MopacMcConnellResult.cpp:152 — `bo` for bond order is terse but the comment on the line names it; acceptable given the local scope.

## 3. Visible math structure (grouping)

- McConnellResult.cpp:78-94 — The K kernel and the full M tensor are built in two adjacent nested loops; the shape is clear because each is preceded by its formula comment. Good — this is the model the rest of the file should match.
- McConnellResult.cpp:275-288 — Per-category T2 totals and the full-M decompose are four sequential blocks with no grouping header separating "symmetric-traceless category features" from "full shielding contribution." — A blank-line + one signpost would make the two output families distinct.
- MopacMcConnellResult.cpp:180-183 — `weighted_M` / `weighted_f` named once then accumulated; the weighting is visible at one site. Clean grouping of the Mopac delta.

## 4. Function / method naming

- McConnellResult.cpp:54, MopacMcConnellResult.cpp:47 — Both files define a file-local `ComputeBondKernel` with identical signature but different return type (`BondKernelResult` vs `MopacBondKernelResult`); fine since both are `static` (internal linkage), but a reader scanning both files sees two same-named functions. — Acceptable; no rename needed.
- McConnellResult.cpp:362 `PackST_MC` / MopacMcConnellResult.cpp:293 `PackST_MMC` — Names say "pack spherical tensor"; the `_MC`/`_MMC` suffixes disambiguate two identical static helpers. Clear enough.
- McConnellResult.h:55 `SampleShieldingAt` — Says what it returns (SphericalTensor at a point); good. Header comment at .cpp:330-331 correctly notes no atom-specific filters apply.

## 5. Comments as signposts

- McConnellResult.h:1-24 and MopacMcConnellResult.h:1-22 — Header docblocks are exemplary: formula, term-by-term irrep attribution, the f-is-not-trace caveat. This is the standard the .cpp blocks fall short of.
- McConnellResult.cpp:268 — "Apply traceless projection to fix floating-point drift" describes the *unused* lambda, not the inline code that actually does it; stale/misleading. — Move the comment to the inline blocks or delete with the lambda.
- McConnellResult.cpp:274 — "Extract symmetric part for T2 features" omits that the next line de-traces; partial. — `// symmetric traceless part for T2`.
- McConnellResult.cpp:119-120 — Filter-set comment is accurate and useful. Good.
- MopacMcConnellResult.cpp:94-97 — Compute docblock states the bond-order weighting and the zero-order skip rationale clearly. Good.
- McConnellResult.cpp:184 "Store in BondNeighbourhood" — restates the next five lines; low value but harmless.

## 6. Correctness (secondary)

- McConnellResult.h:37 `MCCONNELL_CUTOFF_A = 10.0` and MopacMcConnellResult.h:35 are unused constants; the actual cutoff comes from `CalculatorConfig::Get("...cutoff")` (cpp:139, :344, :129). The named constant and the config value can silently diverge. — Not a bug today; flag that the constexpr is dead and could mislead.
- McConnellResult.cpp:202 vs 214 — CO branch tracks `best_co_f`, midpoint, direction, and kernel; CN branch tracks only dist + kernel (no f/midpoint/direction). Intentional per the output set (no `nearest_CN_f`), but verify the asymmetry is by design, not omission.
- McConnellResult.cpp:257 `if (best_co_dist < NO_DATA_SENTINEL)` — guards on sentinel; the `dir_norm` near-zero guard at :259 is a sound degeneracy check. OK.
- MopacMcConnellResult.cpp:242-248 — Nearest-CO/CN T2 decompose `best_*_bo * best_*_kernel.K`; the weight is reapplied at decompose time rather than stored pre-weighted. Consistent with `best_co_f_weighted` already holding the weighted scalar — check that double-weighting is not happening (it is not: `best_co_kernel.K` is the unweighted kernel, weight applied once here). OK on inspection.
- Both files — `M_total` accumulates the full asymmetric tensor and is decomposed directly (cpp:288, :265); the category T2 paths symmetrize+de-trace but the full-shielding path does not. This matches the header's claim that the full M carries T0+T1+T2 (trace and antisymmetric parts are physical). Consistent — no de-trace should be applied there, and none is. Good.
- McConnellResult.cpp:343-348 (`SampleShieldingAt`) — re-implements the singularity guard, cutoff, and near-field exclusion inline rather than via the filter set used in `Compute`; the `near_field_exclusion_ratio * bond_len` test must stay in sync with `DipolarNearFieldFilter`'s definition. — check these two encode the same exclusion; divergence would make grid samples disagree with per-atom values.
