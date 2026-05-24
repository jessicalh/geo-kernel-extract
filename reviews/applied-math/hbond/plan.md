# hbond — fix plan (readability, no number changes)

Targets: `src/HBondResult.h`, `src/HBondResult.cpp`
Algorithm: McConnell dipolar kernel applied to DSSP-resolved backbone
N···O H-bonds (kernel-form, amide-H only). Source at the N···O midpoint,
field direction h_hat = donor N → acceptor O, dedup on `(donor_N, acceptor_O)`.

## 1. Summary

This file tells a coherent story. The header states the physical object,
the closed-form tensor, and the d / h_hat conventions cleanly. The `.cpp`
has a clear four-stage scaffold (resolve DSSP partners → build filter set →
accumulate per-atom tensors → write features), already signposted with
banner comments. `ComputeHBondKernel` matches the header formula term for
term. No sign or value problem survives tracing (see §4 and §5).

The fix pass touches **comments and internal local names only**. It will:
- name the one genuinely opaque internal field (`HBondKernelResult::f`),
  which tracing shows is **dead** — recommend deletion, with the fallback of
  a name + comment if the human prefers to keep it as a diagnostic;
- fix the two comments that hardcode a literal "3.5 Å" while the threshold
  is config-driven (`hbond_counting_radius`);
- add a handful of 2–4-word signposts inside the kernel loop and at the two
  symmetric resolution loops;
- distinguish the two repeated GeometryChoice rejection labels by cause.

It will **not** change the algorithm, any output number, any serialized name
(`hbond_shielding`, `hbond_scalars`, the `ConformationResult` atom fields),
or the tensor formula. It will not split the resolution loops or the fused
`M_over_r3` assignment (both decline below with reasons).

Protected output names confirmed downstream: `_catalog.py` keys
`hbond_shielding` / `hbond_scalars` (mechanism `hbond_kernel`); the atom
fields `hbond_shielding_contribution`, `hbond_nearest_dist`, `hbond_inv_d3`,
`hbond_count_within_3_5A`, `hbond_nearest_spherical`, `hbond_nearest_dir`,
`hbond_is_donor/acceptor/backbone` consumed by `learn/bones/AnalysisWriter.cpp`
and `GromacsProtein*`. None of these is proposed for renaming.

## 2. Review-finding ledger

### codex.md

| # | Finding | Disposition |
|---|---------|-------------|
| C1 | `:72` kernel loop jumps from geometry to compact tensor; add block labels (`// field direction`, `// angular factor`, `// tensor kernel`) | **adopted** → §3 E1 (signposts, not a per-line split) |
| C2 | `:89` `result.f` computed, never used here; name it or remove | **adopted** → §3 E2 (traced dead across whole repo; recommend removal) |
| C3 | `:91` split fused tensor expr into `alignment_term` / `hbond_axis_term` / `radial_term` | **declined** — expression matches the header formula one-to-one under the named banner at `:55`; splitting a single decodable assignment into three temporaries inside a 3×3 loop adds noise without adding clarity. claude L9 reaches the same conclusion. |
| C4 | `:148` two DSSP loops symmetric and long; add signposts + explicit endpoint names | **adopted (signposts only)** → §3 E3; loop-body split declined (see §3 note) |
| C5 | `:160` seq-separation filter buried in bookkeeping; add `// sequence exclusion` | **adopted** → §3 E5 |
| C6 | `:182` distance gate is a compound condition (singular + max dist); label or split into named bools | **adopted (label)** → §3 E6; bool-split declined as over-fragmentation |
| C7 | `:316` kernel computed before filter acceptance; add `// geometry for filters` | **adopted** → §3 E7 |
| C8 | `:389` donor/acceptor flagging separated; add `// endpoint flags` or move | **adopted (label only)** → §3 E8; no move (would change behavior ordering risk) |
| C9 | `:428` `PackST_HB` appears abruptly; add `// feature packing` | **adopted** → §3 E9 |
| C10 | `:17` rename `h_b`/`h_a` → `h_hat_a`/`h_hat_b` in header formula | **declined** — header line 20 already defines `h_hat = direction donor→acceptor`; the formula uses `d̂` and `h` as the tensor indices `a,b`. Renaming `h_a`→`h_hat_a` collides visually with the `d̂_a` notation and over-decorates a math line. The convention is stated one line below. |
| C11 | `:44` `donor_N`/`acceptor_O` hold indices not atoms → `_atom` suffix | **declined** — they hold atom *indices* (`size_t`), and `_atom` would imply an atom object/reference, which is more misleading, not less. The struct doc at `:34-40` already says "resolve to atoms: donor_N: the backbone N". Leave as-is. |
| C12 | `:49` `h_hat` acceptable | **acknowledged, no action** |
| C13 | `:66` `M_over_r3` → `kernel_tensor_Ainv3` / `shielding_kernel_Ainv3` | **declined** — `M_over_r3` ties the field name directly to the header formula symbol `M_ab / r³`; that linkage is the readability win. The units are stated in the struct banner ("Angstrom⁻³"). Renaming away from the formula symbol loses the trace. |
| C14 | `:67` `f` → `angular_factor` if kept | **superseded by C2** — recommend removal; if kept, adopt the rename (`axial_scalar`, see claude N1). |
| C15 | `:79` `d` → `midpoint_to_atom` | **declined** — `d` is the formula's `d = r_atom − r_partner` (header `:20`); kept to match the symbol. The comment signpost in E1 covers the frame. |
| C16 | `:80` `r` → `atom_distance` | **declined** — `r` is the formula radial variable; the `1/r³` structure is most legible with `r`. |
| C17 | `:148` `ri`/`bi`/`acc_ri`/`don_ri` → spelled-out residue names | **partially adopted** → §3 E4 renames `bi`→`partner_slot` and `acc_ri`/`don_ri`→`acceptor_residue_idx`/`donor_residue_idx`; `ri` kept (tight loop index, ubiquitous idiom here). |
| C18 | `:311` `count_3_5` hides the configured radius → `nearby_hbond_count` | **adopted** → §3 E10 (rename local `count_3_5` → `nearby_hbond_count`; the *output* field `hbond_count_within_3_5A` is protected and unchanged) |
| C19 | `:428` `PackST_HB` terse → `PackHBondSphericalTensor` | **adopted** → §3 E11 (file-local static; safe rename) |
| C20 | `:84` group the kernel stages with labels | **duplicate** of C1 → §3 E1 |
| C21 | `:91` tensor construction worst hotspot, expose three terms | **duplicate** of C3 → declined |
| C22 | `:327` label seq separation `// endpoint sequence gap` | **adopted** → §3 E12 |
| C23 | `:380` add `// nearest contribution` | **adopted** → §3 E13 (combined with claude L4: also note "recompute for nearest only") |
| C24 | `:395` add `// decompose total tensor` | **adopted** → §3 E14 |
| C25 | `:41` `Compute` acceptable | **acknowledged, no action** |
| C26 | `:46` `SampleShieldingAt` clear | **acknowledged, no action** |
| C27 | `:72` `ComputeHBondKernel` → `ComputeHBondTensorKernel` | **declined** — the return struct holds `M_over_r3` (a tensor) plus `distance`; "Kernel" already names the dipolar kernel object. Adding "Tensor" is redundant with the struct field. Marginal; left to avoid churn on a clear name. |
| C28 | `:26` banner verbose → `// resolved H-bond` | **declined** — the banner states the Kabsch-Sander provenance and the donor/acceptor→atom resolution that a chemist needs; claude L48 explicitly says keep. Trimming to three words loses the provenance. |
| C29 | `:55` banner heavy → `// dipolar tensor kernel` + formula | **declined** — same rationale; the "Same derivation as McConnell with b_hat → h_hat" line is the load-bearing cross-reference to the kernel catalogue. Keep. |
| C30 | `:117` keep the SpatialIndex `(void)` comment | **acknowledged, no action** (keep) |
| C31 | `:128` Step 1 comment long → shorten | **declined** — the comment enumerates the four skip conditions, which is exactly the contract of the loop; claude L48 says keep. Shortening to "resolve DSSP partners" drops the skip-condition list that documents the gates. |
| C32 | `:162`/`:183` GeometryChoice labels process-flavored/noisy → `// record rejected H-bond` | **adopted (refined)** → §3 E15: distinguish by cause (`// record: seq-sep reject` vs `// record: distance reject`), per claude L50/L5. |
| C33 | `:276` Step 2 comment too explanatory → `// evaluation filters` | **declined** — the comment explains *what each filter means physically* (SelfSource, DipolarNearField), which is the non-obvious part; claude L48 says keep. |
| C34 | `:293` Step 3 comment → `// accumulate atom tensors` | **declined (keep)** — current wording adds "1/r³ decay handles range naturally", which justifies the no-cutoff accumulation (claude L52 flags that justification as good). |
| C35 | `:351` "3.5A" comment conflicts with configurable radius → `// nearby H-bond count` | **adopted** → §3 E10 (comment fixed alongside the rename) |
| C36 | `:418` `// near-field cutoff` keep/shorten | **acknowledged** — keep; minor, folded into E16 note. |
| C37 (Correctness) | `:148` bounds: does `dssp.AllResidues()` size match `ResidueCount()`? | **resolved coherent** → §4.A. Traced: `DsspResult.cpp:84` `residues_.resize(protein.ResidueCount())`. Invariant holds. No code change. |
| C38 (Correctness) | `:416` redundant singularity guard after kernel already returns zero at distance==0 | **resolved coherent** → §4.D. The guard reads `kernel.distance` which is left at its default `0.0` when the kernel bails on the singularity; `0.0 < singularity_guard_distance` is true, so the guard correctly skips the degenerate sample. Not redundant *in effect*, but a reader can misread it — add a one-line note (§3 E16). |

### claude.md

| # | Finding | Disposition |
|---|---------|-------------|
| L1 | `:148-259` two ~50-line near-identical loops; add a one-line role signpost at each head | **adopted** → §3 E3 (duplicate of C4) |
| L2 | `:160-246` GeometryChoice `Record` lambdas bury the rejection reason | **adopted (labels)** → §3 E15; declines moving/extracting the lambda (no refactor per brief) |
| L3 | `:89 vs 93-97` `result.f` never consumed; confirm or flag dead | **adopted** → §3 E2 (duplicate of C2; traced dead — see §4.B) |
| L4 | `:374-387` `nearest_kernel` recomputes the kernel already computed at `:316` | **adopted (comment)** → §3 E13: `// recompute kernel for the nearest H-bond only`. Not a bug — the inner loop did not retain per-hbond kernels. |
| L5 | `:67` `f` opaque, it is `(3cos²θ−1)/r³` axial scalar → `axial_scalar` or comment | **adopted** → §3 E2 (duplicate of C2/C14) |
| L6 | `:311/365/371` `count_3_5`/`hbond_count_within_3_5A` hardcode 3.5 while threshold is config | **adopted for the local; declined for the output field** → §3 E10. Local `count_3_5` → `nearby_hbond_count`. The serialized atom field `hbond_count_within_3_5A` is a protected contract name (read in `AnalysisWriter.cpp`, `GromacsProtein*`); renaming it is out of scope — note added that "3.5A" is the *default* value of `hbond_counting_radius`. |
| L7 | `:146,163` `resolution_key` name suggests a lookup key; it's a record counter → `choice_record_seq` or comment | **adopted** → §3 E17 (rename local `resolution_key` → `choice_record_seq`) |
| L8 | `:428,443` `PackST_HB` / `_HB` suffix says nothing about layout → comment `// layout: T0, T1[3], T2[5]` | **adopted** → §3 E11 (rename + layout comment) |
| L9 | `:91-98` fused `M_ab` is decodable given the header banner; no change needed | **adopted (no-op)** — confirms declining C3 |
| L10 | `:308,361,395` build→accumulate→decompose is visible; good | **acknowledged, no action** |
| L11 | `:194-203,248-257` field-by-field `ResolvedHBond` fill reads fine | **acknowledged, no action** |
| L12 | `:72,109,h.49` function names clean | **acknowledged, no action** |
| L13 | banner/step comments are the right kind; keep | **acknowledged (keep)** — reinforces declining C28/C29/C31/C33/C34 |
| L14 | `:117-119` SpatialIndex `(void)` comment good | **acknowledged (keep)** (duplicate of C30) |
| L15 | `:162,183,216,237` four verbatim GeometryChoice labels → distinguish seq-sep vs distance | **adopted** → §3 E15 (duplicate of C32) |
| L16 | `:351` stale "3.5A" literal in comment → "within counting radius" | **adopted** → §3 E10 (duplicate of C35) |
| L17 | `:360` "1/r³ decay handles range" comment good | **acknowledged (keep)** |
| L18 | `:378` "all DSSP H-bonds are backbone" comment good | **acknowledged (keep)** |
| L19 (Correctness) | `:355` `nearest_dist < nearest_dist` only works if `NO_DATA_SENTINEL` is large-positive | **resolved coherent** → §4.C. `NO_DATA_SENTINEL = 99.0` (`PhysicalConstants.h:77`), positive and above `hbond_dipolar_max_distance` use for atom-to-midpoint distances of accepted bonds. Works as intended. No code change. |
| L20 (Correctness) | `:144-145` dedup key `(donor_N, acceptor_O)` correctly suppresses the same bond resolved from both ends; ordering consistent | **resolved coherent** → §4.E. Confirmed: both loops key N-first/O-second. No change. |
| L21 (Correctness) | `:419 vs 334/288` inline near-field test in `SampleShieldingAt` vs `DipolarNearFieldFilter`; confirm they match | **resolved coherent** → §4.F. Traced `KernelEvaluationFilter.h:177`: filter accepts when `distance > source_extent * near_field_exclusion_ratio`, with `source_extent = hb.distance` (N···O). The inline `SampleShieldingAt` skips when `distance < near_field_exclusion_ratio * hbond_distances_[hi]` (same N···O length). Identical criterion. No change. |
| L22 (Correctness) | `:177-180,199-200` sign/direction `d=O−N`, `h_hat=d/dist`, midpoint=½(N+O) consistent with header | **resolved coherent** → §4.G. Confirmed against header `:14,:20`. No change. |
| L23 (Correctness) | `:386` `hbond_inv_d3 = 1/nearest_dist³` uses atom-to-midpoint distance, while `hbond_distances_` is N···O; confirm intent | **resolved → question** → §4.H and §6 Q1. The two lengths are deliberately different (`nearest_dist` is atom→midpoint, the kernel's `r`; `hbond_distances_` is the source extent N···O). The `inv_d3` feature is the atom's `1/r³` falloff to the nearest source, which is the physically right scalar feature. Reads coherent; flagged as a confirm-with-author because no consumer comment states the intent explicitly. Recommend a one-line comment (§3 E18) rather than a code change. |

## 3. Edits that don't move numbers

All edits are comment text or internal local-name changes. No output name,
no formula, no control flow changes.

- **E1** — `HBondResult.cpp:84-97` — add three terse signposts inside
  `ComputeHBondKernel`: `// field direction (atom relative to midpoint)`
  before `Vec3 d = ...`; `// angular factor` before the `cos_theta` /
  `result.f` lines; `// dipolar tensor kernel, M_ab / r³` before the 3×3
  loop. (covers C1, C20)
- **E2** — `HBondResult.cpp:67,89` — **remove** `double f` from
  `HBondKernelResult` and the `result.f = (3cos²θ−1)/r³` assignment at
  `:89`. Traced dead repo-wide (§4.B). *Fallback if author wants to keep it:*
  rename `f` → `axial_scalar` and add `// (3cos²θ−1)/r³ axial dipolar scalar;
  diagnostic, not consumed`. (covers C2, C14, L3, L5)
- **E3** — `HBondResult.cpp:152, 206` — add a one-line role signpost at each
  resolution loop head:
  `// donor side: this residue's backbone N → partner's acceptor O` and
  `// acceptor side: partner's donor N → this residue's backbone O`.
  (covers C4, L1)
- **E4** — `HBondResult.cpp:153,154,207,208` — rename loop/locals:
  `bi` → `partner_slot`, `acc_ri` → `acceptor_residue_idx`,
  `don_ri` → `donor_residue_idx`. Keep `ri`. (covers C17 partial)
- **E5** — `HBondResult.cpp:160,214` — add `// sequence exclusion` before each
  `seq_sep <` gate. (covers C5)
- **E6** — `HBondResult.cpp:182,236` — add `// source distance gate
  (singular or beyond max N···O)` before each compound `dist <` / `dist >`
  condition. (covers C6)
- **E7** — `HBondResult.cpp:316` — add `// kernel geometry (computed first so
  filters can read distance/extent)`. (covers C7)
- **E8** — `HBondResult.cpp:389` — add `// endpoint flags: is this atom a
  donor N / acceptor O of any resolved bond`. (covers C8)
- **E9** — `HBondResult.cpp:428` — add `// feature packing: T0 | T1[3] | T2[5]`
  above the static (combined with E11). (covers C9, L8)
- **E10** — `HBondResult.cpp:311,351,352` — rename local `count_3_5` →
  `nearby_hbond_count`; change the `:351` comment to
  `// count H-bonds within hbond_counting_radius (default 3.5 Å) of this atom`.
  Output field `hbond_count_within_3_5A` unchanged (protected). (covers C18,
  C35, L6, L16)
- **E11** — `HBondResult.cpp:428,443` — rename file-local static
  `PackST_HB` → `PackHBondSphericalTensor` (two call sites, both in this file)
  and add the layout comment from E9. (covers C19, L8)
- **E12** — `HBondResult.cpp:327` — add `// endpoint sequence gap (min to
  either donor/acceptor residue)`. (covers C22)
- **E13** — `HBondResult.cpp:380` — add `// recompute kernel for the nearest
  H-bond only (inner loop did not retain per-bond kernels)`. (covers C23, L4)
- **E14** — `HBondResult.cpp:395` — add `// decompose accumulated total
  tensor`. (covers C24)
- **E15** — `HBondResult.cpp:162,183,216,237` — replace the four verbatim
  `// ---- GeometryChoice: hbond resolution ----` banners with cause-specific
  labels: `// record: seq-sep reject` (at :162, :216) and
  `// record: distance reject` (at :183, :237). (covers C32, L2, L15)
- **E16** — `HBondResult.cpp:416` — add `// distance==0 here means the kernel
  hit the singularity guard and returned zero; skip it` so a reader does not
  read this as a no-op duplicate of the in-kernel guard. (covers C38, C36)
- **E17** — `HBondResult.cpp:146` + use sites — rename local
  `resolution_key` → `choice_record_seq`. (covers L7)
- **E18** — `HBondResult.cpp:386` — add `// 1/r³ to the NEAREST bond's
  midpoint (atom-to-midpoint r, distinct from the N···O source extent in
  hbond_distances_)` to document the two-length distinction. (covers L23)

**Declined structural change (noted once):** neither resolution loop is
split or de-duplicated into a helper, and the fused `M_over_r3` assignment is
not broken into temporaries. The brief forbids refactors/new abstractions,
and both reviews' own correctness passes (L9, L11) judge the fused math and
the field-fill as decodable. Signposts (E3) assert the symmetry without
moving code.

## 4. Usage notes (the sign/value reasons discovered)

**A. DSSP residue indexing invariant (C37).** `Compute` indexes
`dssp.AllResidues()[ri]` for `ri ∈ [0, protein.ResidueCount())`.
`DsspResult.cpp:84` resizes `residues_` to exactly `protein.ResidueCount()`.
The two index spaces are the same. Coherent; comment unnecessary.

**B. `HBondKernelResult::f` is dead (C2/L3).** `result.f` is written at
`:89` and read **nowhere** — not in `HBondResult.cpp`, not in any other `src/`
file (`HBondKernelResult` and `ComputeHBondKernel` are file-local; grep across
`src/` returns only the producer line). It does not reach any NPY, H5 field,
or atom field. The only persisted scalar features are `hbond_nearest_dist`,
`hbond_inv_d3`, `hbond_count_within_3_5A` (`WriteFeatures` + `AnalysisWriter`).
Producer/consumer: no consumer exists. Recommend removal (E2).

**C. `nearest_dist` sentinel comparison (L19).** Initialized to
`NO_DATA_SENTINEL = 99.0` (`PhysicalConstants.h:77`), a large positive number.
The `kernel.distance < nearest_dist` test therefore admits the first accepted
bond (its distance ≪ 99) and tracks the minimum thereafter. `nearest_hb_idx`
stays `SIZE_MAX` when no bond is accepted, gating the `:374` block. Coherent.

**D. The `:416` singularity guard in `SampleShieldingAt` (C38).** When
`ComputeHBondKernel` bails at the singularity (`r < singularity_guard_distance`),
it returns with `result.distance` at its default `0.0`. The caller's guard
`kernel.distance < singularity_guard_distance` is then true and skips the
sample. So the second guard is the mechanism by which the early-return is
*observed* by the caller — not a redundant duplicate. Comment added (E16) to
prevent the misreading the reviews flagged.

**E. Dedup key direction (L20).** Both the donor-side and acceptor-side
loops build the dedup key as `(N atom index, O atom index)` — N first, O
second — so a bond resolved from either end produces the same key and the
shared `seen` set suppresses the duplicate. No aliasing. Coherent.

**F. Near-field criterion parity, `Compute` vs `SampleShieldingAt`
(L21).** `DipolarNearFieldFilter` (`KernelEvaluationFilter.h:176-177`)
accepts when `ctx.distance > ctx.source_extent * near_field_exclusion_ratio`,
and `Compute` sets `ctx.source_extent = hb.distance` (the N···O length,
`:322`). `SampleShieldingAt` skips when `kernel.distance <
near_field_exclusion_ratio * hbond_distances_[hi]`, and `hbond_distances_`
stores the same N···O length (`:269`). Same threshold, same length, opposite
sense of the inequality (accept vs skip) — consistent. The grid path also
keeps the singularity guard inline. Producer and both consumers agree.

**G. Direction / sign convention (L22).** `d = O_pos − N_pos`,
`h_hat = d / |d|`, `midpoint = ½(N+O)`. Matches the header `:14,:20`
("donor N → acceptor O"). The kernel's own `d = atom_pos − midpoint` is the
field-point vector (a different `d` from the resolution-stage `d`; both are
local and each matches its own comment). Coherent. (The two `d`s share a
letter but never coexist in one scope; not worth a rename.)

**H. `hbond_inv_d3` uses atom-to-midpoint r, not N···O (L23).** `:386`
computes `1/nearest_dist³` where `nearest_dist` is the kernel `r` (atom to
the nearest bond's midpoint). `hbond_distances_` / `hbond.distance` is the
N···O source extent. These are intentionally different lengths: `inv_d3` is
the per-atom `1/r³` falloff feature to the nearest source; the N···O extent
is the near-field-filter scale. Consumers (`AnalysisWriter.cpp:198`,
`GromacsProtein.cpp:213`, `_catalog.py` `hbond_scalars` slot 1) read it as the
atom's distance feature, consistent with atom-to-midpoint. Reads coherent;
no explicit comment states the intent today, so this is also Q1 below — the
edit (E18) is a comment, no number moves regardless of the answer.

## 5. Bug-by-exhaustion candidates

**None.** Every sign/value question raised by the two reviews resolved to
**coherent** under tracing (§4.A–H). The only "dead code" finding (`f`,
§4.B) is exhaustively traced to no consumer and is a removal, not a number
change.

## 6. Questions & Ambiguities

- **Q1 (low stakes, confirm intent).** `hbond_inv_d3` is computed from the
  atom-to-midpoint distance (`nearest_dist`), while `hbond_distances_` tracks
  the N···O source extent (§4.H). Both readings are internally consistent and
  the feature semantics ("atom's 1/r³ to nearest source") fit the
  atom-to-midpoint length. Confirm this is the intended length for the
  `hbond_inv_d3` ML feature before E18's comment is finalized. If the author
  actually intends N···O⁻³, that *would* be a number change and falls outside
  this readability pass — raise separately.
- **Q2 (keep-or-remove `f`).** §4.B shows `HBondKernelResult::f` is dead.
  Default recommendation is removal (E2). If there is an out-of-tree or
  planned consumer (e.g. a trajectory-scope `*TrajectoryResult` that has not
  landed), keep it and take the rename+comment fallback instead. Author's
  call; flagged because removal touches a (file-local) struct definition.
