# Fix plan — BiotSavartResult.{h,cpp}

## 1. Summary

`BiotSavartResult.h` and `.cpp` tell a coherent physics story end to end: SI
wire-segment field → Johnson-Bovey double-loop sum → unit-current geometric
kernel `G_ab = -n_b B_a · PPM_FACTOR` → spherical decomposition → per-atom
accumulation + per-type sums → grid samplers → NPY export. The header narrative
is genuinely strong (model, citation, sign derivation with a worked numeric
check, units). Two things make a reader re-read: (a) one header line states the
kernel sign **without** the minus that the rest of the file (and the code) uses,
and (b) `Compute()`'s body interleaves the physics thread with three
`GeometryChoice` provenance side-quests, a find-or-create scan, and a
twice-computed ring-frame coordinate block.

The fix pass will: correct the one stale sign line, add short signposts that let
the eye skip provenance/plumbing back to the physics, give a few opaque locals
physical names, and de-duplicate two near-identical comments. It will **not**
touch any number, sign, threshold, output name, accumulation order, or filter
logic. The two "correctness" flags both reviews raised (`SampleShieldingAt`
returns an un-intensity-scaled kernel; ring proximity counters use `+=`) are
traced below and land as **coherent** — comment clarifications only, no code
change.

Output/contract names protected and unchanged: `bs_shielding`,
`bs_per_type_T0`, `bs_per_type_T2`, `bs_total_B`, `bs_ring_counts` (NPY +
`_catalog.py` keys), and the public methods `Compute`, `SampleShieldingAt`,
`SampleBFieldAt`, `WriteFeatures`.

---

## 2. Review-finding ledger

### codex.md

| # | Finding (loc) | Disposition |
|---|---|---|
| C1 | h:10 vs h:23 — `G_ab` stated without minus, then with minus | **adopted** → §3 E1 |
| C2 | h:54 `SampleShieldingAt` "Returns the shielding tensor" but returns unit-current kernel | **adopted** (comment only) → §3 E2; reason in §4.2 |
| C3 | cpp:145 filter-setup comment long before physics loop | **adopted** → §3 E3 |
| C4 | cpp:176 GeometryChoice block interrupts physics | **adopted** → §3 E4 |
| C5 | cpp:219 physics resumes, buried; wants `// unit-current field` | **adopted** → §3 E5 |
| C6 | cpp:237 "Find or create RingNeighbourhood" wants `// ring-neighbour features` | **adopted** → §3 E6 |
| C7 | cpp:254 cylindrical block wants `// target in ring frame` | **adopted** → §3 E7 |
| C8 | cpp:280 storing G/decomp/B is clear | **declined** (no action requested) |
| C9 | cpp:297 accumulation clear | **declined** (no action) |
| C10 | cpp:319 proximity counts may count over pre-existing neighbours | **adopted** (comment only) → §3 E8; reason in §4.3 |
| C11 | cpp:47 `dl_m` acceptable | **declined** (no action) |
| C12 | cpp:48 `dA_m`/`dB_m` confusable with B; suggest `start_to_point_m`/`end_to_point_m` | **declined** — header defines them, both reviews call it minor; rename touches `WireSegmentField` internals only but adds little; weighed call to keep |
| C13 | cpp:59 `factor` → `biot_savart_scale` | **adopted** → §3 E9 |
| C14 | cpp:60 `dotTerm` → `endpoint_projection_term` | **adopted** → §3 E10 |
| C15 | cpp:92 `B` fine | **declined** (no action) |
| C16 | cpp:167 `G_total` clear | **declined** (no action) |
| C17 | cpp:237 `rn` acceptable, `ring_neighbour` cheaper | **declined** — `rn` is local, consistent with the 4 sibling ring calcs' identical block; comment signpost (E6) covers it |
| C18 | cpp:251 `d` too generic → `center_to_atom` | **adopted** → §3 E11 |
| C19 | cpp:256 `d_plane` → `in_plane_offset` | **declined** — `d_plane` is clear once `center_to_atom` + signpost land; renaming one of a `z`/`rho`/`theta` set adds churn without clarity gain |
| C20 | cpp:433 `PackST` → `PackSphericalTensor` | **adopted** → §3 E12 |
| C21 | cpp:87 unit conversion well grouped | **declined** (no action) |
| C22 | cpp:96 loop construction clear | **declined** (no action) |
| C23 | cpp:225 tensor construction well explained | **declined** (no action) |
| C24 | cpp:251 ring-frame block: split with `// axial/radial` + `// azimuth reference` labels | **adopted (partial)** → §3 E7/E13; one signpost on the coordinate block + keep the existing azimuth comment, not two extra sub-labels (would over-fragment) |
| C25 | cpp:285 B-field cylindrical repeats 251–258 math; wants `// project B into ring frame` | **adopted** → §3 E14; reason in §4.4 |
| C26 | cpp:312 `total_G_spherical` vs `bs_shielding_contribution` — why both? | **adopted** (comment only) → §3 E15; reason in §4.5 |
| C27 | cpp:356/386 samplers repeat ring-acceptance guards; wants matching comments | **adopted** → §3 E16 |
| C28 | h:57 `SampleShieldingAt` likely misleading; suggest rename `SampleShieldingKernelAt` | **declined (rename)** — public method consumed in `ui/ComputeWorker.cpp`, `h5-reader/`, doc, `ui/CLAUDE.md`; renaming all 8 sibling `SampleShieldingAt` for naming consistency is a cross-subproject contract change out of scope for a readability pass. Document the unit-current return instead → E2 (adopts the comment half of C28) |
| C29 | cpp:43 `WireSegmentField` clear | **declined** (no action) |
| C30 | cpp:77 `JohnsonBoveyField` → `JohnsonBoveyBField` | **declined** — return type `Vec3` + `// Tesla` already name the quantity; both reviews call the current name clear; cross-file rename (sampler lambda + two call sites) not worth it |
| C31 | cpp:433 `PackST` too compressed | **duplicate** of C20 |
| C32 | h:3 file signpost useful | **declined** (no action) |
| C33 | h:18 citation dense but acceptable | **declined** (no action) |
| C34 | cpp:26 section header clear | **declined** (no action) |
| C35 | cpp:87 good signpost keep | **declined** (no action) |
| C36 | cpp:96 Upper/Lower loop good | **declined** (no action; but see CL10) |
| C37 | cpp:176 `// GeometryChoice: ring current` → `// provenance sampler` | **adopted** → §3 E4 |
| C38 | cpp:200 good signpost | **declined** (no action) |
| C39 | cpp:207 `// GeometryChoice: near-field exclusion` → `// record exclusion` | **adopted** → §3 E17 |
| C40 | cpp:225 sign comment too verbose through 231; trim | **adopted** → §3 E18; reason in §4.1 |
| C41 | cpp:237 → `// ring-neighbour features` | **duplicate** of C6 |
| C42 | cpp:254 → `// target in ring frame` | **duplicate** of C7 |
| C43 | cpp:263 azimuth comment shorten to `// asymmetric ring side` | **declined** — both reviews elsewhere (CL: "genuinely good signpost… explains why"); the *why* (N-side vs C-side on HIE/TRP) is the load-bearing part, keep |
| C44 | cpp:327 `// GeometryChoice: ring shells` → `// record ring shells` | **adopted** → §3 E19 |
| C45 | cpp:348 "No filter set"/"still applied" reads contradictory | **adopted** → §3 E20 |
| C46 | h:10/h:23 sign documentation inconsistency (Correctness) | **duplicate** of C1 |
| C47 | h:54/cpp:386 `SampleShieldingAt` kernel-not-shielding (Correctness) | **duplicate** of C2 |
| C48 | cpp:321 proximity counters double-count if Compute re-runs | **duplicate** of C10 |
| C49 | cpp:258 `theta = atan2(.,abs(z))` folds sign of z | **adopted** (comment only) → §3 E21; reason in §4.6 |

### claude.md

| # | Finding (loc) | Disposition |
|---|---|---|
| CL1 | cpp:53,57,165,291,321-324,370,376,400-402 — inlined `CalculatorConfig::Get("literal")` in conditionals; hoist to named `const double` | **adopted (partial)** → §3 E22; hoist the per-iteration / per-conditional lookups to named consts at scope top. Not a refactor — same values, named once |
| CL2 | cpp:251-274 vs 286-295 cylindrical frame computed twice (`z`/`rho` then `z_coord`/`rho_mag`); signpost or compute once | **adopted (signpost)** → §3 E14; reason in §4.4. No "compute once" (would change structure) |
| CL3 | cpp:176-196 GeometryChoice block dwarfs the computation it precedes; end marker | **adopted** → §3 E4 |
| CL4 | cpp:258 `theta` recomputes `d_plane.norm()` (already `rho`), folds sign | **adopted** (comment + reuse `rho`) → §3 E21 / E23; reason in §4.6 |
| CL5 | cpp:238-244 find-or-create plumbing; `// find-or-create` signpost | **adopted** → §3 E6 |
| CL6 | cpp:44-45 `dA`/`dB` collide with B but acceptable given comment | **declined** — matches C12 weighed call |
| CL7 | cpp:90 `halfI_A`/`I_A` good | **declined** (no action) |
| CL8 | cpp:272-273 `cos_phi`/`sin_phi` etc clean | **declined** (no action) |
| CL9 | h:64/cpp:433 `PackST`/`out` terse but acceptable | **partial** — `out` kept; `PackST`→`PackSphericalTensor` via E12 (codex C20) |
| CL10 | cpp:302 `ti`, :170 `ri`, :161 `ai` single letters, consistent | **declined** (no action) |
| CL11 | cpp:232-235/:409-412/:189-192 — kernel build appears 3×, only 1 commented; comment the bare two | **adopted** → §3 E24 |
| CL12 | cpp:312-317 `total_G_spherical` vs `bs_shielding_contribution` verify distinct | **duplicate** of C26 → §3 E15; reason §4.5 |
| CL13 | cpp:43-111 `WireSegmentField`→`JohnsonBoveyField` staging exemplary | **declined** (no action) |
| CL14 | cpp:386-418 `SampleShieldingAt` duplicates `SampleBFieldAt` guards; fine | **declined** (no action; matching comments via E16) |
| CL15 | Function naming axis clean | **declined** (no action) |
| CL16 | cpp:117 header comment hardcodes "(15A)" — stale-comment hazard | **adopted** → §3 E25 |
| CL17 | cpp:225-231 sign comment duplicates header verbatim; trim | **duplicate** of C40 → §3 E18 |
| CL18 | cpp:263-265 azimuth comment good, keep | **declined** (no action; supports declining C43) |
| CL19 | cpp:96,100 Upper/Lower loop use `z`/`d` not locally defined; `// upper loop (+lobe_offset)` | **adopted** → §3 E26 |
| CL20 | cpp:294 `B_phi (zero by axial symmetry)` precise, keep | **declined** (no action) |
| CL21 | h:1-32 exemplary header, no change | **declined** (no action, modulo E1) |
| CL22 | h:10 vs h:23 sign (Correctness) | **duplicate** of C1 → §3 E1 |
| CL23 | cpp:291 vs :269 two different small-rho thresholds (`near_zero_vector_norm_threshold` vs hardcoded `1e-10`) | **can't tell → §6 Q1**; comment-flag only, no value change |
| CL24 | cpp:370-376 samplers skip `MinDistanceFilter`/`RingBondedExclusionFilter` that Compute applies; comment says intentional | **coherent → §4.7**; clarify comment via E20 |
| CL25 | cpp:303,457,467 magic `8` ring-type bound; if enum grows, silent truncation | **adopted (comment) → §3 E27**; reason §4.8 / §6 Q2 |
| CL26 | cpp:315-317 fields diverge across repeated Compute if `total_G_tensor` carries prior accumulation | **coherent → §4.5** (Compute runs once per fresh conf); clarified by E15 |

Nothing from either review is dropped without a row above.

---

## 3. Edits that don't move numbers

Comment fixes, signposts, named intermediates, and justified internal renames.
None change a value, sign, threshold, or output name.

- **E1** `h:10` — change `G_ab = n_b * B_a * PPM_FACTOR` to `G_ab = -n_b * B_a *
  PPM_FACTOR` so it matches h:23 and the code (cpp:192,235,412). Comment conforms
  to code (the minus is real, see §4.1).
- **E2** `h:55` — change "Returns the shielding tensor (sum over all rings within
  cutoff)" to "Returns the spherical decomposition of the summed **unit-current**
  G kernel (sum over rings within cutoff); multiply by ring intensity to get
  shielding." (See §4.2.) Method name unchanged (C28 declined).
- **E3** `cpp:145-147` — replace the long filter-setup prose with a short
  signpost: `// Atom-ring filters: near-field validity + ring/bonded exclusion.`
  (Keep the three `filters.Add(...)` lines self-documenting.)
- **E4** `cpp:176` — change `// ---- GeometryChoice: ring current ----` to
  `// Provenance: record one sampler per ring (physics resumes below).`
- **E5** `cpp:219` — add block label `// Unit-current B-field (I = 1 nA).` above
  the `JohnsonBoveyField` call (it already has the intensity-independence note).
- **E6** `cpp:237` — change `// Find or create RingNeighbourhood for this ring`
  to `// Ring-neighbour features: find-or-create per-(atom,ring) record.`
- **E7** `cpp:254` — change `// Cylindrical coordinates in ring frame` to
  `// Target position in ring frame (axial z, radial rho, polar theta).`
- **E8** `cpp:319` — change `// Ring proximity counts` to `// Ring proximity
  counts (each ring appears once in ring_neighbours; see find-or-create above).`
  (See §4.3.)
- **E9** `cpp:59` — rename local `factor` → `biot_savart_scale`
  (file-local, `WireSegmentField` only).
- **E10** `cpp:60` — rename local `dotTerm` → `endpoint_projection_term`
  (file-local).
- **E11** `cpp:251` — rename local `d` → `center_to_atom` within the
  find-or-create block (lines 251-261; references at 252,255,256). File-local.
- **E12** `cpp:433` + `h:64` — rename file-local static `PackST` →
  `PackSphericalTensor`. Internal to the .cpp (declared static); the `.h:62-64`
  comment references "Pack order for SphericalTensor" — no external caller.
  Carry-through: the one call site at `cpp:448`.
- **E13** `cpp:254-261` — keep as one labelled block (E7); do **not** add the two
  extra sub-labels codex C24 suggested (over-fragmentation).
- **E14** `cpp:285` — add `// Project B into the same ring frame (z/rho re-derived
  from the geometry above).` above the second cylindrical block. (See §4.4.)
- **E15** `cpp:315-317` — add one comment distinguishing the two decompositions:
  `// total_G_spherical: running sum across calculators on this atom.`
  `// bs_shielding_contribution: this calculator's per-call BS sum only.`
  (See §4.5.)
- **E16** `cpp:356` and `cpp:386` — add matching one-liners on both samplers'
  guard cascade: `// Grid acceptance: singularity guard + inside-source guard +
  distance cutoff (no per-atom/topology filters; grid points are not atoms).`
- **E17** `cpp:207` — change `// ---- GeometryChoice: near-field exclusion ----`
  to `// Provenance: record this atom-ring pair as filter-excluded.`
- **E18** `cpp:225-231` — trim the in-loop sign block to one line plus the
  formula, deferring the worked example to the header:
  `// G_ab = -n_b B_a PPM_FACTOR (shielding-sign convention; derivation + worked`
  `// example in the header).` (See §4.1 — the minus stays; only the duplicated
  prose is trimmed.)
- **E19** `cpp:327` — change `// ---- GeometryChoice: ring shells ----` to
  `// Provenance: record this atom's ring-shell counts.`
- **E20** `cpp:348-354` — rewrite the sampler header comment from
  "No filter set ... DipolarNearFieldFilter still applied" to:
  `// Same physics as Compute(), single point. Grid points are not atoms, so the`
  `// per-atom/topology filters (MinDistance, RingBondedExclusion) do not apply;`
  `// the inline distance/singularity/inside-source guards keep the multipole`
  `// valid.` (Matches what the code actually does — see §4.7.)
- **E21** `cpp:258` — add inline comment on `theta`:
  `// polar angle from the ring axis, folded to [0, pi/2] via abs(z)` (See §4.6.)
- **E22** `cpp:53,57,165,291,321-324,370,376,400-402` — hoist the
  `CalculatorConfig::Get("…")` lookups used inside conditionals/loops to named
  `const double` at the top of their enclosing scope (e.g.
  `const double endpoint_guard = CalculatorConfig::Get("biot_savart_wire_endpoint_guard");`
  in `WireSegmentField`; `const double ring_cutoff = CalculatorConfig::Get("ring_current_spatial_cutoff");`
  and the four `ring_proximity_shell_*` in `Compute`; the singularity / inside-source /
  cutoff trio in each sampler). Same values; readability only.
  **Carry-through note:** purely local; no signature or contract change. The four
  shell consts must keep the exact `ring_proximity_shell_{1..4}` keys (do not
  rename config keys — those are an external contract in `calculator_params.toml`).
- **E23** `cpp:258` — optional: reuse the already-computed `rho` instead of
  recomputing `d_plane.norm()` inside the `atan2` (CL4). Value-identical
  (`rho == d_plane.norm()` one line up). Flag as a judgment call: it is a
  micro-edit that removes a duplicate norm; safe because `rho` is in scope and
  equal. If preferred to stay literal, decline and keep E21's comment only.
- **E24** `cpp:189-192` and `cpp:409-412` — add the one-line kernel comment that
  cpp:225-235 carries, to the two bare copies:
  `// G_ab = -n_b B_a PPM_FACTOR (see header sign convention).`
- **E25** `cpp:117` — drop the hardcoded "(15A)" from the header comment or mark
  it "(default 15 A; value from config)" so it cannot drift from the config key.
- **E26** `cpp:96,100` — change `// Upper loop (z = +d)` / `// Lower loop (z =
  -d)` to `// Upper loop (+lobe_offset along normal)` / `// Lower loop
  (-lobe_offset along normal)` so the symbols match the local name
  `lobe_offset_ang`.
- **E27** `cpp:303` (and mirror at WriteFeatures 457/467) — add comment on the
  `ti < 8` guard: `// 8 = RingTypeIndex count; additions beyond 8 are dropped —
  keep in sync with the enum.` (Comment only; see §4.8 / Q2 for the structural
  question.)

---

## 4. Usage notes (the discovered reasons)

### 4.1 The kernel sign `G_ab = -n_b B_a · PPM_FACTOR` — coherent

The minus is **real and load-bearing**, applied identically at three code sites
(cpp:192, 235, 412) and in the sibling `HaighMallionResult` (cpp:387,
`G_ab = -n_b V_a`). It comes from the shielding-tensor definition
`σ_ab = -∂B_a^sec/∂B_{0,b}`, and the header carries a worked check
(`I = -12, G_T0 = -0.116 → σ = +1.40 ppm`, shielded). The only defect is
documentation: `h:10` omits the minus that `h:23` and all code apply — a stale
line, not a sign bug. **Producer/consumer agree:** `_catalog.py:123,149`
stamps `bs_shielding` with `sign_convention = "σ_ab = -dB_a^sec/dB_{0,b}"`, the
same string shared family-wide (`_SHIELD_SIGN`). Disposition: **coherent**; fix
the comment (E1, E18). Same family as the Haigh-Mallion header self-contradiction
(HM header h:15 also drops the minus its code applies) — flag noted, but HM is
out of scope for this file's plan; the BS fix is the template.

### 4.2 `SampleShieldingAt` returns a unit-current kernel, not intensity-scaled shielding — coherent

This is the **same convention as the per-atom path**, on purpose. `Compute()`
calls `JohnsonBoveyField(..., 1.0, ...)` (cpp:221, current = 1 nA) and stores the
unit-current kernel; the comment at cpp:219-220 states "The geometric kernel is
independent of intensity." The thesis-wide rule is **"the system outputs kernels,
not shielding"** — calibration multiplies kernels by intensity / regression
weights downstream. The catalog encodes this in the **units**:
`bs_shielding` is `"ppm_T_per_nA"`, `bs_per_type_T0/T2` are `"ppm_T_per_nA"` —
i.e. per-nanoampere, *not* ppm. Calibration consumes these as features
fed into ridge, which supplies the intensity scaling.

`SampleShieldingAt` does the identical thing for an arbitrary grid point, and its
**only consumers are viewers**: `ui/ComputeWorker.cpp:427` reads `stResult.T0`
purely to colour an isosurface grid (`grid.bsT0[idx]`), and `h5-reader` /
`QtBiotSavartCalc` mirror "the same sign convention" (their headers say so
explicitly). No consumer multiplies by `ring.Intensity()` and none expects ppm
from this method. **Producer and all consumers agree** that this is a kernel.
Disposition: **coherent**; the method name is misleading but is a public,
cross-subproject contract shared by 8 sibling calculators — fix the header
comment (E2), decline the rename (C28).

### 4.3 Ring proximity counters (`n_rings_within_*`) — coherent, no double-count

The `+=` at cpp:321-324 looks like an accumulation hazard, but three facts make
it safe:
1. The counters default-initialise to 0 on `ConformationAtom` (ConformationAtom.h:135-138).
2. **Only `BiotSavartResult` writes them** — grep across `src/*.cpp` shows no
   other calculator touches `n_rings_within_*`.
3. The loop iterates `ca.ring_neighbours`, the fat-union shared by all 5 ring
   calculators — but every calculator uses the **same find-or-create-by-ring_index**
   pattern (BS cpp:239-244; HM 299; McConnell/PiQuad/RingSusc/Dispersion at the
   greps in §trace), so a given (atom, ring) pair has **exactly one**
   `RingNeighbourhood` entry no matter how many ring calculators run. Counting
   over it counts each ring once.
4. `Compute()` is invoked once per conformation: `OperationRunner::Attach`
   (OperationRunner.cpp:53) attaches one result instance per type and refuses a
   second; in trajectory mode each frame gets a **fresh** `conf` from
   `tp.TickConformation(...)` (Trajectory.cpp:317), so counters start at 0 each
   frame. The `+=` therefore behaves as a single-pass `=`.

The one residual subtlety (worth the E8 comment, not a fix): BS runs **before**
the other 4 ring calculators (OperationRunner.cpp:165 vs 166+), so at BS's count
time `ca.ring_neighbours` holds only BS's own filtered entries — exactly the set
BS intends to count. Disposition: **coherent**; clarify with the E8 comment.

### 4.4 Twice-computed ring-frame coordinates (cpp:251-261 vs 286-295) — coherent (intentional split)

The first block (find-or-create, runs once per (atom,ring)) computes the
**static** geometry stored on the neighbourhood: `z`, `rho`, `theta`, azimuth.
The second block (runs every accepted pair) computes the **B-field projection**
into the same frame (`B_rho`, `B_z`). They share the `z`/`d_plane` derivation
because they describe the same ring frame, but they live in different scopes
(one guarded by find-or-create, one unconditional) and feed different fields.
This is not accidental duplication. Disposition: **coherent**; add the E14
signpost so a reader knows the repeat is deliberate. No "compute once" merge —
that would move structure.

### 4.5 `total_G_spherical` vs `bs_shielding_contribution` (cpp:315-317) — coherent, distinct fields

Both decompose a G sum, but **different** sums:
- `ca.total_G_spherical` decomposes `ca.total_G_tensor`, which is `+=`-accumulated
  **across calculators** on the atom (cpp:314 adds BS's `G_total` into a running
  cross-calculator total). This is the atom's aggregate geometric tensor.
- `ca.bs_shielding_contribution` decomposes `G_total` — **this calculator's
  per-call BS contribution only** — and is what `WriteFeatures` exports as
  `bs_shielding.npy` (cpp:448).

CL26's divergence worry ("if Compute re-runs, the `+=` total drifts from the
per-call field") does not arise: Compute runs once per fresh conf (§4.3). The two
fields are *intended* to differ — one is the BS-only export feature, the other is
the cross-calculator running total. Disposition: **coherent**; clarify with E15.

### 4.6 `theta = atan2(d_plane.norm(), abs(z))` (cpp:258) — coherent (deliberate fold)

`abs(z)` folds the polar angle to `[0, π/2]`: a point at axial offset `+z` and
one at `-z` get the same `theta`. This is consistent with the diamagnetic ring
current being symmetric about the ring plane (the Johnson-Bovey field is
mirror-symmetric across the plane). `theta` is stored on `RingNeighbourhood` as a
geometric descriptor; the sign of `z` is preserved separately in `rn->z`
(cpp:259), so no information is lost — the azimuth (cos_phi/sin_phi) plus signed
`z` fully locate the atom. Disposition: **coherent**; the name `theta` alone
doesn't reveal the fold, so add the E21 comment. (E23 optionally reuses `rho`
for the redundant `d_plane.norm()` — value-identical, judgment call.)

### 4.7 Samplers skip the per-atom filters (cpp:370-376, 400-402) — coherent

`Compute` applies `MinDistanceFilter` + `DipolarNearFieldFilter` +
`RingBondedExclusionFilter`; the samplers apply only inline distance /
singularity / inside-source guards. The difference is correct: the per-atom
filters encode topology (a ring vertex must not feel its own ring; bonded
neighbours are excluded), which is meaningless for an arbitrary grid point. The
physical near-field guard (`distance < geom.radius` → inside source) **is** kept
in the samplers. Disposition: **coherent**; the existing comment is just worded
confusingly — E20 rewrites it.

### 4.8 Magic `8` ring-type bound (cpp:303, 457, 467) — coherent today, fragile

`8` is the count of `RingTypeIndex` values; `per_type_G_T0_sum[8]` /
`per_type_G_T2_sum[8][5]` are fixed-size, and the catalog hard-codes
`bs_per_type_T0` = 8, `bs_per_type_T2` = 40. The guard `ti < 8` silently drops a
ring whose type index ≥ 8. Today there are exactly 8 types so nothing is dropped
— **coherent**. The fragility (a future enum addition truncates silently) is real
but is a structural/contract concern, not a current bug; E27 documents it and Q2
raises whether `RingTypeIndex::Count` should back the literal. No value change.

---

## 5. Bug-by-exhaustion candidates

**None.** Every sign/value flag from both reviews traced to a coherent reason:

- G-sign (C1/C46/CL22): documentation-only staleness at h:10; code is consistent
  at all 3 sites and matches `_catalog.py` and the HM sibling. Exhausted across
  producer + catalog + viewer consumers.
- `SampleShieldingAt` kernel-vs-shielding (C2/C47): matches the per-atom path,
  the catalog units (`ppm_T_per_nA`), `learn/` feature use, and both viewer
  consumers; no consumer expects intensity-scaled ppm. Exhausted.
- Proximity counters (C10/C48): single writer, single Compute per conf, one
  neighbourhood entry per ring via shared find-or-create. Exhausted across all 5
  ring calculators + OperationRunner + Trajectory frame loop.
- theta fold (C49/CL4): symmetric by ring-plane mirror symmetry; signed `z`
  retained separately. Coherent.

No edit in §3 changes a number.

---

## 6. Questions & Ambiguities

- **Q1 (CL23) — two small-rho thresholds.** The in-plane azimuth block (cpp:269)
  guards with a hardcoded `1e-10`, while the B-field projection block (cpp:291)
  uses `CalculatorConfig::Get("near_zero_vector_norm_threshold")`. Both guard the
  same `rho → 0` degeneracy. I could not confirm whether the divergence is
  intentional (two different physical tolerances) or an oversight (one was
  config-ised, the other missed). Both produce a near-zero-rho fallback that
  zeroes the affected quantity, so behaviour is benign at the limit — but the
  inconsistency is a readability/maintenance smell. **Do not change the value
  without a decision.** If they should agree, the fix is to route cpp:269 through
  the same config key; that is a value-touching change and out of scope here.
- **Q2 (CL25) — `8` vs `RingTypeIndex::Count`.** Is the literal `8` (cpp:303,
  457, 467; array sizes on ConformationAtom; catalog widths) meant to be backed
  by a `RingTypeIndex::Count` constant, so that adding a ring type is caught at
  compile time / fails loud instead of silently truncating at `ti < 8`? E27 adds
  a warning comment, but whether to introduce the named constant (and widen the
  arrays + catalog together) is a contract decision for the author. Left as a
  question, not an edit.
- **Q3 — HM header twin.** `HaighMallionResult.h:15` has the identical
  drop-the-minus inconsistency (`G_ab = n_b (H.n)_a` in prose vs `-n_b V_a` in
  code). It is the same family as C1/E1 and should get the same one-line fix, but
  HM is a separate file outside this plan's targets. Flagging so it is not lost.
