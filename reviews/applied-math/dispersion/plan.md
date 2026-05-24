# Fix plan — DispersionResult.{h,cpp}

## 1. Summary

The file tells a coherent story. The physics is correct and matches the
header: a per-vertex London dispersion kernel `K_ab = S(r)·(3 d_a d_b/r^8 −
δ_ab/r^6)` summed over ring vertices, traceless per vertex, with a CHARMM C¹
switching taper, a ring-level near-field filter, and through-bond exclusion.
Nothing in the algorithm, its signs, or its numbers needs to change, and this
pass proposes no such change.

The genuine readability costs both reviews converge on are real and worth
fixing:

1. **The header does too much.** It carries the full algorithm note (kernel,
   tracelessness proof, filter policy, switching rationale, cutoff
   justification) ahead of a small public interface. A reader must separate
   contract from implementation policy.
2. **Two back-to-back block comments (`.cpp:26-49` and `:51-61`) restate the
   same taper rationale.** One should absorb the other.
3. **The switching-function distances are read from `CalculatorConfig` four
   times to form two squares (`.cpp:66-67`),** so a single `Rc²`/`Rs²`
   squaring reads as an eight-call expression.
4. **Inside `Compute`, the four lines of actual physics are buried** under
   three verbose `GeometryChoice::Record` provenance blocks and a hand-rolled
   find-or-create. Signpost comments let the eye skip the bookkeeping.
5. **Bare `8`, `5`, `40`** appear where a named constant already exists in the
   codebase (`kAromaticRingTypeCount`, Types.h:232).

What the pass will **not** touch: the algorithm, the kernel expression, the
signs, the cutoffs/onsets, the `ti < 8` accumulation boundary (it is a
documented project-wide convention, see §4), the public method name
`SampleShieldingAt` (shared across all seven ring/field calculators — see
ledger), and every serialized output name (`disp_shielding`,
`disp_per_type_T0`, `disp_per_type_T2`, and the `RingNeighbourhood` /
`ConformationAtom` field names, all consumed by the SDK and Qt readers).

---

## 2. Review-finding ledger

Legend: **adopted** (→ §3 edit), **declined** (reason), **duplicate** (of).

### codex.md

| # | Finding | Disposition |
|---|---------|-------------|
| C1 | `.h:3` header is a full algorithm note, not an interface signpost; trim to purpose/units/stored quantities | **adopted** → E1 |
| C2 | `.cpp:26` switching comment overlong and partly repeated by 51-61 | **adopted** → E2 |
| C3 | `.cpp:184` ring-filter comment interrupts the narrative; shorten signpost | **adopted** → E7 (kept brief, not stripped) |
| C4 | `.cpp:226-241` filtering + GeometryChoice fused; add signpost + shorter label | **adopted** → E8 |
| C5 | `.cpp:244-258` through-bond exclusion: logging dominates; `// reject bonded atom` | **adopted** → E8 |
| C6 | `.cpp:270-298` core vertex-sum easy to miss; stronger signpost | **adopted** → E9 |
| C7 | `.cpp:302-327` find-or-create also computes ring-frame geometry; split with two comments | **adopted** → E10 |
| C8 | `.cpp:329-343` storage + decomposition + per-type + total packed; label them | **adopted** → E11 |
| C9 | `.cpp:363-388` `SampleShieldingAt` uses a different-looking ring filter; signpost + confirm intent | **adopted (signpost)** → E12; intent settled in §4 / §6 |
| C10 | `.cpp:62` `DispSwitchingFunction` → `DispersionSwitchingFunction` | **adopted** → E3 |
| C11 | `.cpp:66-70` `rc2/rs2/num/den` → `cutoff_sq/switch_sq/numerator/denominator` | **adopted** → E4 (with distance hoist) |
| C12 | `.cpp:94-98` `DispVertexResult::K` → `kernel_tensor` | **declined** — `K` matches the header's `K_ab` formula symbol and the field is read by name in `ComputeDispVertex`; the struct is file-local but the one-letter name is formula-faithful and the surrounding comment ties it to `K_ab`. Low value, mild churn. |
| C13 | `.cpp:110` `S` → `switch_weight` | **declined** — header/comment establish `S(r)` as the switching function; `S` is formula-faithful. (claude N1 agrees.) |
| C14 | `.cpp:113` `d` → `atom_minus_vertex` | **adopted (lighter)** → E5: `d_av` with a `// d = r_atom − r_vertex` note, matching the header symbol while disambiguating from the ring-frame `d` reused at :315. |
| C15 | `.cpp:271-273` `K_ring/s_ring` → `ring_kernel_tensor/ring_scalar_sum` | **adopted (lighter)** → E9: keep `K_ring`, rename `s_ring`→`scalar_ring` is not worth it; add signpost instead. See note. |
| C16 | `.cpp:279` `vr` → `vertex_result` | **adopted** → E9 |
| C17 | `.cpp:303` `rn` → `ring_neighbour` | **adopted** → E10 |
| C18 | `.cpp:315-323` `d/z/d_plane` → `center_to_atom/normal_offset/in_plane_offset` | **adopted (lighter)** → E10: rename `d`→`center_to_atom`; keep `z`/`d_plane` but add the `// ring-frame coordinates` signpost. `z` maps directly to the stored field `rn->z`; renaming the local away from the field name reduces, not improves, traceability. |
| C19 | `.cpp:392` `PackST_D` → `PackSphericalTensorFloat64` | **adopted** → E13 |
| C20 | `.cpp:118-124` group inverse powers `inv_r6/inv_r8` | **declined** — `r6`/`r8` are already named locals built incrementally (`r2`,`r6`,`r8`); the expression uses division by them, and adding reciprocal intermediates changes the arithmetic shape (could perturb rounding). Per the brief, no edit that risks moving a number. The existing names are clear. |
| C21 | `.cpp:120-124` split tensor into `dyadic_term/isotropic_term` | **declined** — the two-line loop with the inline `K_ab` formula comment directly above is already the clearest form; introducing per-element named terms inside the `a,b` loop would fragment the one place the tensor is built (claude §3 calls this block "clean"). |
| C22 | `.cpp:319-323` `// ring-frame coordinates` | **adopted** → E10 (duplicate of C7/C18) |
| C23 | `.cpp:331` `// spherical decomposition` | **adopted** → E11 |
| C24 | `.h:55` `Compute` acceptable | **noted, no change** |
| C25 | `.h:59` `SampleShieldingAt` → `SampleDispersionAt` (returns dispersion tensor, not total shielding) | **declined** — see ledger note below and §6. The name is a shared interface convention across all seven ring/field calculators. |
| C26 | `.cpp:101` `ComputeDispVertex` → `ComputeDispersionVertexKernel` | **declined** — file-local static; `DispVertex` reads cleanly in context and the longer name is verbose without adding meaning given the file is `DispersionResult.cpp`. Marginal; left to keep the diff focused on the high-value edits. |
| C27 | `.cpp:136` `BondedToVertices` → `AtomsBondedToRingVertices` | **declined** — file-local; current name + its header comment are self-documenting (claude §4 calls it "clear and self-documenting"). |
| C28 | `.cpp:392` `PackST_D` → `PackSphericalTensorFloat64` | **duplicate** of C19 |
| C29 | `.cpp:26-49` replace local label with `// switching function` | **adopted** → E2 |
| C30 | `.cpp:75-92` block repeats the header; replace with `// per-vertex kernel` | **adopted (lighter)** → E6: keep a trimmed block (the per-vertex tracelessness identity is the one thing not obvious from the header math), drop the duplicated exclusion paragraph. |
| C31 | `.cpp:131-134` shorten to `// bonded vertex set` | **declined** — the existing 2-line comment explains *why* (through-bond vs through-space); shortening to a label drops the rationale that justifies the exclusion. Keep. |
| C32 | `.cpp:152-161` keep as `// atom-ring accumulation` | **adopted (lighter)** → E2-adjacent: trim the `Compute` block comment to remove sentences duplicated by the header, keep the dataflow overview. |
| C33 | `.cpp:260` "dispersion taper" is recording metadata; `// record taper parameters` | **adopted** → E8 |
| C34 | `.cpp:281-282` "Only fires when…" too process-specific; `// switching floor` | **adopted** → E9 |
| C35 | `.cpp:302` procedural comment → `// neighbour record` | **adopted** → E10 |
| C36 | `.cpp:335` `// Per-type accumulation` good terse signpost | **noted, keep** |
| C37 | Correctness check: `.cpp:372-376` grid guard simpler than `Compute`'s filter; verify intentional | **investigated → coherent**, see §6 |

**Ledger note on C25/C14:** `SampleShieldingAt` is the identical method name on
`BiotSavartResult`, `HaighMallionResult`, `McConnellResult`, `HBondResult`,
`PiQuadrupoleResult`, `RingSusceptibilityResult`, and `DispersionResult`
(verified by grep across `src/`). It is a deliberate uniform interface — the
viewer's grid sampler calls it polymorphically by convention. Renaming it on
one calculator would break the symmetry that makes the family legible and
would not carry through to a contract boundary cleanly. The "kernel not
shielding" point is real but applies to all seven and is a project-wide naming
decision, not a dispersion-specific fix. Out of scope here; raised as a
question in §6.

### claude.md

| # | Finding | Disposition |
|---|---------|-------------|
| L1 | `.cpp:66-67` `rc2`/`rs2` built from 4 inline `Get(...)` calls; hoist distances to `r_cut`/`r_switch` | **adopted** → E4 |
| L2 | `.cpp:108` guard fuses two config lookups + two comparisons + early return on one line; split | **adopted** → E5 |
| L3 | `.cpp:218-346` physics buried under three `Record(...)` blocks; one-line signposts per physics step | **adopted** → E8/E9 |
| L4 | `.cpp:302-327` find-or-create + z/rho/theta is a 25-line digression; signpost as bookkeeping | **adopted** → E10 (duplicate of C7) |
| L5 | `.cpp:337,339` bare `8` and `5`; name `kRingTypeCount`/`kT2Components` | **adopted** → E11, using the **existing** `kAromaticRingTypeCount` (Types.h:232) and a new file-local `kT2Components = 5` |
| L6 | `.cpp:404,418` `N*40` and `t*5+c` encode 8×5 with no reminder; add `// 8 ring types × 5 T2` | **adopted** → E14 |
| L7 | `.cpp:110` `double S` acceptable, no change | **noted, agrees with declined C13** |
| L8 | `.cpp:272,296,332` `s_ring`/`K_ring`/`K_total` acceptable given header symbols | **noted** — supports lighter C15 |
| L9 | `.cpp:392` `PackST_D` cryptic; rename or comment layout | **duplicate** of C19 → E13 |
| L10 | `.cpp:315,320` `d`/`d_plane` shadow the kernel's `d`; one-word note `// d: atom → ring center` | **adopted** → E10 (rename local `d`→`center_to_atom`, removes the shadow entirely) |
| L11 | `.cpp:121-124` `K_ab` loop clean | **noted, keep** (supports declined C21) |
| L12 | `.cpp:270-298` vertex-sum well isolated | **noted** |
| L13 | `.cpp:329-341` store + per-type run back-to-back with no separation; blank line + comment | **adopted** → E11 |
| L14 | `.h:59` `SampleShieldingAt` returns kernel not shielding; consider `SampleKernelAt` | **declined** — same reason as C25 (shared interface name); raised in §6 |
| L15 | `.cpp:26-49 / 51-61` two long comments restate the same taper rationale; collapse to one | **adopted** → E2 |
| L16 | `.cpp:45-48` comment cites `4.3/5.0 A` as literals but code reads `CalculatorConfig`; phrase as "default … see config" | **adopted** → E2 (phrase as defaults) |
| L17 | `.cpp:111` `// below switching threshold` good; keep | **adopted (refine)** → E5: change to `// below noise floor` to match the config key `dispersion_switching_noise_floor` (the threshold IS the noise floor, not the switch onset). |
| L18 | `.h:34` `Unit C6 = 1 … learnable` good one-liner | **noted, keep** |
| L19 | `.cpp:108 vs 381` both guard `r > cutoff`; consistent | **investigated → coherent**, §6 |
| L20 | `.cpp:283` noise-floor provenance branch consistent | **investigated → coherent**, §6 |
| L21 | `.cpp:337` `ti < 8` silently drops ring types ≥ 8 while `disp_total` includes them; check if a ring type can have index ≥ 8 | **investigated → coherent (intentional, documented)**, §4 + §6 |
| L22 | `.cpp:247` through-bond skips entire ring if bonded to any vertex; documented intent | **investigated → coherent**, §6 |
| L23 | `.cpp:374-375` grid path omits `DipolarNearFieldFilter` unlike `Compute`; confirm intended | **investigated → coherent (deliberate), but worth a confirming question**, §6 |

Nothing from either review is dropped silently.

---

## 3. Edits that don't move numbers

All edits are comment text, signpost insertion, named locals, regrouping, and
one file-local static rename + one file-local constant. No output name, no
serialized field, no public method name changes. No arithmetic shape changes.

**E1 — `.h:3-36` trim the header to an interface signpost.**
Keep: one-line purpose, the kernel/scalar formulas, the units line, the
`C6 = 1, learnable per ring type` line. Move the tracelessness *proof*, the
filter *policy* paragraphs, and the switching-function *derivation/cutoff
rationale* down into the `.cpp` block comments where the code lives (they
already partly duplicate `.cpp:26-92`). Net: the header states the contract;
the `.cpp` owns the implementation rationale. (codex C1.)

**E2 — `.cpp:26-61` collapse the two taper block comments into one.**
Merge `:26-49` and `:51-61` into a single block: the C¹ formula, the
boundary properties, the Brooks 1983 reference, and one paragraph of cutoff
rationale. Phrase the distances as **defaults read from config**:
`// onset (default 4.3 A) and cutoff (default 5.0 A) are read from
CalculatorConfig; values below are illustrative.` (codex C2/C29, claude
L15/L16.)

**E3 — `.cpp:62` rename `DispSwitchingFunction` → `DispersionSwitchingFunction`.**
File-local static; no cross-file consumer. (codex C10.)

**E4 — `.cpp:63-70` hoist the switching distances; rename math locals.**
```
const double r_switch = CalculatorConfig::Get("dispersion_switching_onset_distance");
const double r_cut    = CalculatorConfig::Get("dispersion_vertex_distance_cutoff");
if (r <= r_switch) return 1.0;
if (r >= r_cut)    return 0.0;
const double switch_sq = r_switch * r_switch;
const double cutoff_sq = r_cut * r_cut;
const double r2 = r * r;
const double numerator   = (cutoff_sq - r2) * (cutoff_sq - r2) * (cutoff_sq + 2.0*r2 - 3.0*switch_sq);
const double denominator = (cutoff_sq - switch_sq) * (cutoff_sq - switch_sq) * (cutoff_sq - switch_sq);
return numerator / denominator;
```
Identical arithmetic; only de-duplicates the config reads and names the parts.
(codex C11, claude L1.)

**E5 — `.cpp:108-113` split the `ComputeDispVertex` guard; name the kernel
offset.** Two guard lines (singularity guard; cutoff guard) instead of one
fused line; refine the floor comment; rename `d`→`d_av` with a one-line note:
```
// singularity guard: skip coincident / near-coincident points
if (r < CalculatorConfig::Get("singularity_guard_distance")) return result;
// outer cutoff: kernel is zero beyond the switching range
if (r > CalculatorConfig::Get("dispersion_vertex_distance_cutoff")) return result;

double S = DispersionSwitchingFunction(r);
if (S < CalculatorConfig::Get("dispersion_switching_noise_floor")) return result;  // below noise floor

Vec3 d_av = atom_pos - vertex_pos;   // d = r_atom − r_vertex (header symbol)
```
Same logic, same numbers. (claude L2/L17, codex C14.)

**E6 — `.cpp:75-92` trim the per-vertex block comment.** Keep the formulas
and the per-vertex tracelessness identity (not obvious from the header alone);
drop the through-bond-exclusion paragraph (it now lives at `BondedToVertices`
and in the header's filter section). Add the terse `// per-vertex kernel`
label. (codex C30.)

**E7 — `.cpp:184-187` shorten the ring-filter block to a signpost** while
keeping the one-clause "why" (breaks down inside the ring). e.g.
`// Ring near-field filter (DipolarNearFieldFilter, source_extent = ring
diameter): the discrete vertex sum is invalid when the field point is inside
the ring.` (codex C3.)

**E8 — `.cpp:226-268` signpost the two exclusion branches; relabel the taper
record.** Add `// reject near-field atom` before the near-field branch
(`:232`), `// reject bonded atom` before the through-bond branch (`:248`), and
change the `dispersion taper` Record label comment at `:260` to
`// record taper parameters` (the branch records metadata, it does not apply
the taper). The `Record(...)` lambda bodies are unchanged. (codex C4/C5/C33,
claude L3.)

**E9 — `.cpp:270-298` signpost the core vertex sum; rename `vr`; refine the
floor comment.** Insert `// --- vertex kernel sum ---` before the loop;
rename `vr` → `vertex_result`; change the `:282` comment to `// switching
floor: vertex tapered below noise floor`. Keep `K_ring`/`s_ring` (formula-
faithful, header-anchored). (codex C6/C16/C34, claude L3/L8.)

**E10 — `.cpp:302-327` signpost + de-shadow the find-or-create block.**
Insert `// locate or create the ring-neighbourhood record` before the search;
rename loop local `rn` → `ring_neighbour`; rename the ring-frame local `d`
(`:315`) → `center_to_atom` (removes the shadow of the kernel's `d`); insert
`// ring-frame coordinates (z along normal, rho in plane, theta polar)` before
`:319`. Keep `z`/`d_plane` (they map to stored fields `rn->z`, and the
in-plane vector is transient). (codex C7/C17/C18/C22/C35, claude L4/L10.)

**E11 — `.cpp:329-341` separate the storage / per-type / total stages.**
Blank line + `// store kernel on the ring-neighbourhood record` before
`:330`; `// spherical decomposition` understanding folded into that block;
keep the existing `// Per-type accumulation` (`:335`). Replace bare `8`/`5`:
```
int ti = ring.TypeIndexAsInt();
if (ti >= 0 && ti < kAromaticRingTypeCount) {           // Pro pyrrolidine (8) is non-aromatic, excluded here
    ca.per_type_disp_scalar_sum[ti] += s_ring;
    for (int c = 0; c < kT2Components; ++c)
        ca.per_type_disp_T2_sum[ti][c] += ring_neighbour->disp_spherical.T2[c];
}
```
where `kAromaticRingTypeCount` is the existing `Types.h:232` constant and
`kT2Components` is a new file-local `constexpr int kT2Components = 5;` near the
top of the `.cpp`. (codex C8/C23/C36, claude L5/L13/L21.) See §4 for why the
total still adds every ring.

**E12 — `.cpp:372-377` signpost the grid-sampling guard.** Add
`// grid-sampling guard: skip inside-ring and out-of-range points (raw field
sample, no provenance recording)` so the reader sees this path is intentionally
simpler than `Compute`. No logic change. (codex C9, claude L23.)

**E13 — `.cpp:392` rename `PackST_D` → `PackSphericalTensor9` and add a layout
comment** `// flatten SphericalTensor to 9 doubles: [T0 | T1(3) | T2(5)]`.
File-local static, single caller in `WriteFeatures`. (codex C19/C28, claude
L9.)

**E14 — `.cpp:404,418` annotate the 8×5 flattening.** Add
`// per_type_T2 layout: 8 ring types × 5 T2 components = 40` near the `N*40`
allocation and the `i*40 + t*5 + c` index. Optionally use `kAromaticRingTypeCount`
and `kT2Components` in the loop bounds (`t < kAromaticRingTypeCount`,
`c < kT2Components`) and the allocation (`N * kAromaticRingTypeCount`,
`N * kAromaticRingTypeCount * kT2Components`) for self-documenting sizes. The
output array widths (8, 40) stay numerically identical and match the catalog.
(claude L6.)

---

## 4. Usage notes (the sign/value questions, traced)

### The `ti < 8` per-type accumulation boundary — COHERENT (intentional, documented)

Both reviews flagged the bare `8` and claude (L21) asked whether a ring type
can have index ≥ 8 such that the per-type NPY and `disp_total` disagree.

**Reason discovered.** `RingTypeIndex` (Types.h) enumerates 9 chemistries with
`Count = 9`. Indices 0–7 are aromatic (`PheBenzene` … `HieImidazole`); index 8
is `ProPyrrolidine`, the saturated proline 5-ring, explicitly annotated
*"Aromaticity = None; Intensity = 0 … Falls outside the
`< kAromaticRingTypeCount` boundary for per-aromatic-type calculator
accumulation."* Types.h:232 defines `inline constexpr int
kAromaticRingTypeCount = 8;` with a guarding `static_assert` that
`ProPyrrolidine == kAromaticRingTypeCount`.

So the `ti < 8` guard is **the documented project convention**, not a silent
drop. It is **identical across `BiotSavartResult.cpp:303`,
`HaighMallionResult.cpp:337`, `PiQuadrupoleResult.cpp:214`, and
`DispersionResult.cpp:337`** — verified by grep. The per-type arrays are sized
8 on `ConformationAtom` (`per_type_disp_scalar_sum[8]`,
`per_type_disp_T2_sum[8][5]`), the catalog declares width 8 / 40, and the SDK
reads 8 / 40. Producer and all consumers agree.

The apparent asymmetry claude noticed is real and correct: `disp_total`
(`:343`) accumulates **every** passing ring including ProPyrrolidine, so the
full per-atom `disp_shielding` tensor includes the saturated ring's
contribution, while the **per-aromatic-type breakdown** deliberately omits it
(a saturated ring has no aromatic-type bucket to land in). This is the
intended physics, not a bug. The only fix is readability: replace the bare `8`
with the existing named constant (E11) so the boundary's meaning is visible at
the call site.

### Units / sign of the dispersion tensor — COHERENT

Catalog: `disp_shielding` carries `units="Angstrom^-6"`,
`sign_convention=_SHIELD_SIGN`, `irreps=_SHIELD_IRREPS`, `tensor_rank=2`,
`mechanism="ring_dispersion"`; `disp_per_type_T0` is `0e`, `disp_per_type_T2`
is `2e`, both `Angstrom^-6`. These match the header's stated units and the
`Decompose`d tensor packed `[T0|T1|T2]` by `PackST_D`. No sign manipulation
happens in this file beyond the kernel's own `−δ_ab/r^6` trace-removal term,
which the header proves traceless. No consumer (SDK `_protein.py`, h5-reader,
ui) re-signs the field. Nothing to reconcile — coherent.

### `SampleShieldingAt` returning a kernel, not shielding — naming, project-wide

Both reviews note the method returns the dispersion *kernel* tensor, and the
project's stated discipline is "the system outputs kernels, not shielding."
True, but the name is the **shared interface name** on all seven ring/field
calculators and is called polymorphically by the viewer grid sampler
(`ui/src/ComputeWorker.cpp:427` calls `bs.SampleShieldingAt`). Changing it on
dispersion alone breaks the family symmetry; changing it on all seven is a
project-wide naming decision outside this single-algorithm pass. Left as a
question (§6).

---

## 5. Bug-by-exhaustion candidates

**None.** No edit in §3 moves a number. The two items most likely to be a bug
on first read — the `ti < 8` clamp vs `disp_total` (claude L21) and the
`SampleShieldingAt` filter difference (codex C9 / claude L23) — both resolve to
documented, intentional, codebase-wide conventions on tracing (see §4 and §6).
The exhaustion for the clamp: traced the enum (Types.h), the named constant +
static_assert (Types.h:232-233), the on-atom array widths
(ConformationAtom.h:247-248), the three sibling calculators using the identical
guard, the catalog widths (8/40), and the SDK reader — every side agrees and
the total-vs-per-type asymmetry is the intended physics.

---

## 6. Questions & Ambiguities

1. **`SampleShieldingAt` name (codex C25 / claude L14).** The reviews are
   right that it returns a kernel, not shielding, and the name overstates. But
   it is a shared method name across seven calculators and is invoked
   polymorphically by the viewer. **Question for the human:** rename the whole
   family (e.g. `SampleKernelAt`) in a dedicated cross-calculator pass, or
   leave the convention? This pass declines to touch it. Not a dispersion-
   specific fix.

2. **Grid path omits `DipolarNearFieldFilter` (codex C9 / claude L23).**
   `SampleShieldingAt` (`:372-376`) uses a plain distance triple
   (`r < singularity_guard`, `r < radius`, `r > ring_current_spatial_cutoff`)
   where `Compute` runs the full `DipolarNearFieldFilter` with
   `source_extent = 2·radius`. The header says the filter is "same physics as
   all other ring calculators." Tracing the siblings would confirm whether the
   grid path *intentionally* samples the raw field (no near-field rejection,
   appropriate for a volumetric grid where you want the field everywhere) or
   should mirror `Compute`. The grid sampler is a viewer-side diagnostic, not a
   feature-producing path, so a difference is plausibly deliberate — but I
   could not find a comment stating intent. **Recommend** the E12 signpost
   ("raw field sample, no provenance recording") *and* a one-line confirmation
   from the human that the grid path is meant to skip the near-field filter. If
   confirmed deliberate, E12's comment is the whole fix; if not, it is a
   consumer-side (viewer-path) discrepancy to revisit — still not a change to
   the production `Compute` numbers.

3. **`kT2Components` constant placement.** E11 proposes a new file-local
   `constexpr int kT2Components = 5`. If a shared `kT2Components` already
   belongs in Types.h alongside `kAromaticRingTypeCount` (the `5` appears in
   every per-type calculator and in `SphericalTensor`), that would be the
   better home — but adding it to Types.h is a cross-file change touching the
   sibling calculators' readability too, beyond this single-algorithm pass.
   **Question:** keep it file-local here, or open a small Types.h constant in a
   separate pass that all per-type calculators adopt together?
