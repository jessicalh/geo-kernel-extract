# tripeptide-pose-assembler — claude review (readability focus)

- **Targets:** src/TripeptidePoseAssembler.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict

**`.h`** — Reads cleanly. The header is unusually well-signposted: the top-of-file block lays out the two-substrate matching discipline before any type, each struct field carries a units/frame comment (`R*(p_dft - cs) + ct`, `R σ R^T, ppm`), and `TripeptidePoseSide` is documented per-enumerator with the physical σ_BB term it serves. A chemist could follow the data contract without reading the .cpp. Minor: the prose-to-code ratio is high and some comments restate naming-history that a reader doesn't need.

**`.cpp`** — The math core (Kabsch, tensor rotation, identity matching) tells a coherent story; the SVD reflection guard and the `R σ R^T` rotation are each one labeled line. The story does break in two places: the central-assembly comment block (lines 315–333) carries a **stale description** that contradicts the code beneath it, and the per-atom matching loop fuses identity-dispatch + candidate-gather + nearest-spatial + emit into one ~75-line body where the actual computation is buried under a 25-line essay. Naming is strong throughout. Correctness is sound modulo two conventions worth confirming.

---

## Findings

### 1. Coherent story / readability (primary)

- `.cpp:315–322` — Comment block describes `AssembleCentral` as "element+nearest-distance (with the gotham sidechain re-rotation around CA-CB)" but the code below does typed-identity matching with a single BB Kabsch and **no** sidechain re-rotation (lines 414–425 explicitly removed `Rsc`) — stale comment directly contradicts the code; as written, a reader must read both blocks and the code to discover which is true. Suggestion: delete this block (lines 315–322); the 324–333 block beneath it is the accurate one.
- `.cpp:450–526` — The per-atom match loop interleaves a 25-line history essay (453–478), the ambiguity dispatch, candidate gather, nearest-spatial scan, outlier stat, and emit — the actual algorithm (relaxed = ambiguity flag → find candidates → pick nearest → emit) is ~15 lines hidden inside 75; as written, a reader must skip three comment blocks to trace the four-step computation. Suggestion: move the lines 453–478 ILE/PHE rationale to a 2-line signpost (`// strict for chemistry-distinct branches, relaxed for graph-automorphic pairs`) and let the loop body read as steps.
- `.h:44–48`, `.cpp:217–221`, `.cpp:233–235`, `.cpp:473–478` — Several comments narrate retired designs ("the previous element-pattern + 5 Å radius", "round-3 fix", "holdover from the validation-threshold-rejection design") — history prose that costs a first-time reader attention without helping them read the current code. Suggestion: trim to the durable fact; the spec doc already holds the history.
- `.cpp:67–69, 72–74` — `ApplyKabsch` and `RotateTensor` are crisp, named, single-line — genuinely clear, no change.

### 2. Naming carries meaning

- `.cpp:40` — `P, Q` for the centered source/target coordinate matrices is terse SVD-textbook notation but opaque to a chemist reading once; the `H = P * Qᵀ` covariance follows from it. Suggestion: `src_centered`, `dst_centered` (the `H` covariance name is fine as-is given the comment context).
- `.cpp:139` — local `rec` (an `AlignedDftAtom`) shadows the meaning of `rec` used everywhere else in the file as the `TripeptideDftRecord` parameter — same name, two different types in scope-adjacent code. Suggestion: rename the local to `out_atom` or `aligned`.
- `.cpp:451, 507` — `perc` and `adat` are cryptic; `perc` is the perceived atom, `adat` the aligned-DFT-atom being built. Suggestion: `perceived` and `aligned_atom`.
- `.cpp:57` — `sumSq` is clear enough given the rmsd context; no change.
- `.h`/`.cpp` — `K`, `R`, `sigma`, `H`, `U`, `V`, `D` in the Kabsch/rotation core are conventional and the comments anchor them; acceptable as math-domain names.

### 3. Visible math structure (grouping)

- `.cpp:36–64` — `KabschAlign` is well-staged: centroids → centered cols → covariance `H` → SVD → reflection guard → rotation → rmsd. The shape is obvious. Clean.
- `.cpp:493–498` — the nearest-spatial scan is a clean named-intermediate block (`best_atom`/`best_dist`). Clean.
- `.cpp:507–525` — the `AlignedDftAtom` construction (8 field assignments + tensor rotate + decompose) duplicates the same build that `EmitAlignedAtom` (139–175) does for the cap path, but inline rather than via the helper — the two emit-paths diverge silently. Not a refactor request, but flag: a reader comparing cap vs central must verify the two builds agree field-for-field. Suggestion: a one-line signpost `// emit (mirrors EmitAlignedAtom; central path keeps outliers)`.

### 4. Function / method naming

- `.cpp:185` — `LarsenLocalToRecIdx` clearly says what it maps (Larsen local index → rec.atoms index); good.
- `.cpp:335` — `IdentityCompatible` reads well as a predicate; the `relaxed` bool param is documented at every call site.
- `.cpp:67` — `ApplyKabsch` says what it does; returns a transformed point. Good.
- Function names are uniformly accurate (`AssembleAlaCap`, `AssembleCentralTyped`, `SubstrateRoleMatches`, `EmitAlignedAtom`). Clean on this axis.

### 5. Comments as signposts

- `.cpp:50–52` — `// Canonical Kabsch: use sign(det(V*Uᵀ)) ...` is a good signpost but tacks on a cross-reference to RmsdTrackingTrajectoryResult/a dated maths review — the date/provenance is noise at the read site. Suggestion: keep `// reflection guard: sign(det(VUᵀ))`, drop the provenance line.
- `.cpp:414–425` — 12-line comment justifying "single rotation, no Rsc" — accurate and load-bearing but verbose. Suggestion: compress to `// single BB Kabsch, no sidechain re-rotation: χ1-grid coarseness stays in residual_vec as ML feature`.
- `.cpp:528–534` — the CYX/HG explanation is genuinely useful (explains the dominant production warning) but is 7 lines at a non-hot site. Acceptable; it earns its length.
- `.cpp:297–302` — defensive-guard comment is honest about unreachability; fine.
- `.cpp:77–83, 123–125, 179–181` — the box-drawn section banners are good navigation aids; keep.

### 6. Correctness (secondary)

- `.cpp:53–56` — reflection guard uses `det(V·Uᵀ)` and places the correction in the diagonal `D = diag(1,1,sign)` with `R = V·D·Uᵀ` — this is the canonical Kabsch form and is correct. No bug.
- `.cpp:72–74` — `RotateTensor` computes `R σ Rᵀ`; correct for rotating a rank-2 Cartesian tensor from source into target frame, consistent with the `K.rotation` (source→target) produced by `KabschAlign`. Frame direction is self-consistent. Check: that `src` is always the DFT/Larsen side and `dst` the protein side at every `KabschAlign` call (lines 266–276, 398–408) — confirmed at both sites, so the tensor lands in protein frame as intended.
- `.cpp:113–115` — `BackboneHA` slot accepts `AlphaHydrogen` role **OR** `locant == Alpha`; the `|| locant==Alpha` widens the match to any α-locant atom (e.g. could it admit CA itself?). Check that no α-carbon reaches this slot — the slot is only dispatched for `res.HA` (line 292), so in practice safe, but the predicate alone is looser than its name. Not a bug given call constraints; flag as a latent looseness.
- `.cpp:503–505` — central path increments `n_above_threshold` but does **not** reject, matching the documented contract (residual is an ML feature). Correct and intentional.
- `.cpp:614` — `out.ok = ok && !out.aligned_atoms.empty()` — a residue where every atom failed identity match returns `ok=false` even though assembly ran; consistent with "hard structural failure" semantics. No bug.
- `.cpp:386` — `rec.larsen->central` is dereferenced after `AssembleCentral` (line 562) checked `rec.larsen.has_value()`; safe. `AssembleAlaCap` guards its own `rec.larsen` at line 236 before deref at 244. Both optional-derefs are guarded.
