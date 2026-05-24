# Fix plan — HaighMallionResult.{h,cpp}

Inputs weighed: `codex.md` (readability pass), `codex-correctness.md`
(earlier correctness pass), `claude.md` (readability pass). Source:
`src/HaighMallionResult.h`, `src/HaighMallionResult.cpp`. Traced
consumers: `src/BiotSavartResult.{h,cpp}`, `src/PiQuadrupoleResult.cpp`,
`src/DispersionResult.cpp`, `src/RingSusceptibilityResult.cpp`,
`src/ConformationAtom.h`, `src/Types.h`, `python/nmr_extract/_catalog.py`,
`h5-reader/`, `ui/src/RestServer.cpp`.

## 1. Summary

The file tells a coherent story. The header opens with the two-step
construction (raw surface integral `H` → shielding kernel `G`), units,
symmetry/trace properties, the quadrature rule with citation, and the
BiotSavart parallel. The `.cpp` reads top-down: quadrature rule →
triangle integral → adaptive subdivision → ring surface integral →
per-atom/per-ring accumulation → grid sampler → feature output. The
math and its signs are correct and **match the BiotSavart sibling
verbatim** where the two share a convention.

Two things genuinely impede a first-pass reader, and both reviews flag
them:

1. **The G-sign documentation discrepancy.** Header line 15 writes
   `G_ab = n_b * (H . n)_a` (no minus); the `.cpp` constructs
   `G(a,b) = -n_b * V_a`. This is a *stale/incomplete comment*, not a
   code bug — the minus is correct and shared with BiotSavart. Fix the
   comment to match the code and carry the BS-style justification.
2. **Symbol density at the physics climax** (`H`, `V`, `G`) and a few
   terse geometry locals. These are *internal* names; under the
   corrected naming rule they should be improved where it aids clarity.

The fix pass will: correct/align comments to the code, add a small
number of 2–4-word signposts, and rename internal locals/helpers where
the rename earns its keep. It will **not** touch the algorithm, the
numbers, the signs, or any serialized name (`hm_shielding`,
`hm_per_type_T0`, `hm_per_type_T2`, and the `RingNeighbourhood` field
names `direction_to_center`, `rho`, `z`, `theta`, `hm_*`). It will not
refactor `SampleShieldingAt` to dedupe the math (the brief forbids
refactors; the duplication is small and intentional).

The strongest single piece of evidence: the BS header is internally
inconsistent in the *exact same way* (BS:10 `G_ab = n_b * B_a` without
minus, BS:23 `G_ab = -n_b * B_a * PPM_FACTOR` with minus + full
justification). BS resolves it by *also* carrying the signed
justification paragraph; HM omits that paragraph from the header and
only carries it inline at cpp:287–290. The fix brings HM's header up to
BS's standard.

## 2. Review-finding ledger

Every finding from all three reviews, with disposition. C = `codex.md`,
CC = `codex-correctness.md`, CL = `claude.md`.

| # | Src | Finding | Disposition |
|---|-----|---------|-------------|
| 1 | C h:15 | header `G_ab = n_b*(H.n)_a` omits minus used in cpp | **adopted** → E1 (comment-to-code) |
| 2 | C h:25 | `G = n (x) V` reads transposed vs `-V_a n_b` | **adopted** → E1 |
| 3 | C cpp:281 | rename/annotate `H`,`V`,`G` at the physics transition | **adopted** (rename `V`→`effective_field`; keep `H`,`G` + add signpost) → E4, E6 |
| 4 | C cpp:296 | neighbourhood bookkeeping interrupts G→store narrative; add label | **adopted** → E7 |
| 5 | C cpp:310 | geometry block `d,z,rho,theta` needs a signpost | **adopted** → E8 |
| 6 | C cpp:360 | `SampleShieldingAt` "Same physics" hides different filters; label exclusions | **adopted** → E9 |
| 7 | C cpp:35 | `TriQuadPoint`,`lambda`,`weight` good | no-op (acknowledged) |
| 8 | C cpp:83 | `H` underspecified; `surface_integral_H` | **declined** → E5 note: keep `H`, add signpost; rename to a long name hurts the dense `H(a,b)+=` kernel loop readability |
| 9 | C cpp:90 | `rS` → `surface_point` | **adopted** → E2 |
| 10 | C cpp:91 | `rho` → `field_minus_surface` | **partially adopted** → E3 (rename the *vector* `rho`→`sep`; keep wording short) |
| 11 | C cpp:285 | `V` → `effective_field` | **adopted** → E6 |
| 12 | C cpp:291 | `G` → `shielding_kernel` | **declined** → keep `G` (matches header, BS, and the `G(a,b)` index loop); add signpost instead |
| 13 | C cpp:310 | `d` → `center_to_atom`/`atom_from_center` | **adopted** → E10 (`d`→`atom_from_center`) |
| 14 | C cpp:315 | scalar `rho` reuses vector name; `radial_distance` | **adopted** → E11 (`rho` local→`in_plane_radius`; serialized field `rn->rho` unchanged) |
| 15 | C cpp:401 | `PackST_HM` → `PackHaighMallionTensor` | **adopted** → E13 |
| 16 | C cpp:98 | split fused kernel into `kernel_ab`+`quadrature_weight` "if touched" | **declined** → not touching; the inline `// K_ab=...` comment already names it; splitting a hot 3×3 loop is churn |
| 17 | C cpp:131 | name the two thresholds as locals before the branch | **declined** → the `if(level==0)…else if(level==1)` already reads clearly; named locals would not reduce ambiguity |
| 18 | C cpp:281 | Step 1/2/3 grouping is the strongest part | no-op (acknowledged) |
| 19 | C cpp:325 | storage/decomposition grouping is readable | no-op (acknowledged) |
| 20 | C cpp:380 | sampler duplicates `H→V→G`; same signpost wording | **adopted** → E9 (mirror the signpost) |
| 21 | C h:50 | `Compute` acceptable | no-op |
| 22 | C h:54 | `SampleShieldingAt` → `SampleShieldingTensorAt` | **declined** → cross-file public method; consumers in `SampleShieldingAt` callers (see §4); rename cost not justified by clarity gain — name already says it samples shielding |
| 23 | C cpp:44 | `Gauss7` → `TriangleGauss7Rule` | **adopted** → E14 (file-local static, zero carry-through) |
| 24 | C cpp:79 | `AccumulateTensor` → `AccumulateTriangleSurfaceIntegral` | **adopted** (as `AccumulateTriangleIntegral`, see #38) → E15 |
| 25 | C cpp:125 | `AccumulateAdaptive` → `AccumulateAdaptiveTriangleIntegral` | **adopted** → E16 |
| 26 | C cpp:158 | `SurfaceIntegral` → `ComputeSurfaceIntegralH` | **adopted** → E17 |
| 27 | C h:3 | header block long; move empirical/process prose out | **declined** → the header is the documented "model" opener (claude verdict); trimming risks losing the citation/convention context the brief wants kept |
| 28 | C h:23 | shorten "This parallels BiotSavartResult…" to `// Biot-Savart analogue` | **declined** → the parallel + "whether they agree on T2 is the empirical question" is load-bearing thesis context, not noise |
| 29 | C cpp:27 | banner → terse `// triangle quadrature` | **declined** → banners carry the rule citation (Stroud/Dunavant) and barycentric/orbit structure; that is signal, not decoration |
| 30 | C cpp:68 | banner → `// dipolar surface integral` | **declined** → same; banner carries the `H_ab` formula and `rho` definition |
| 31 | C cpp:108 | banner → `// adaptive subdivision` | **declined** → banner carries the level thresholds + point-count growth; useful |
| 32 | C cpp:238 | `// ---- GeometryChoice: surface integral ----` opaque → `// record source ring` | **adopted** → E18 (clarify wording) |
| 33 | C cpp:250 | `// Apply filter` restates code → `// pair filters` | **adopted** → E19 |
| 34 | C cpp:269 | banner → `// adaptive refinement` | **declined** → already a short signpost; fine as is |
| 35 | C cpp:287 | important sign comment, keep; align with header | **adopted** → E1 (header alignment) |
| 36 | C cpp:396 | feature layout comment clear | no-op |
| 37 | C h:15 / cpp:287 | doc/convention minus-sign mismatch; check BS + Boyd-Skrynnikov | **adopted** → E1; reason in §4 (coherent, not a bug) |
| 38 | C h:25 / cpp:294 | outer-product wording `n(x)V` vs `-V_a n_b` | **adopted** → E1 |
| 39 | C cpp:316 | `theta=atan2(.,abs(z))` discards side-of-plane sign; check consumers | **disposition: coherent** → §4; shared verbatim with BS; comment signpost via E8 |
| 40 | CC h:15 | G minus-sign mismatch; check BS/literature | **duplicate** of #1/#37 |
| 41 | CC h:25 | `G=n(x)V` transposed → `G_ab=-V_a n_b` | **duplicate** of #2/#38 |
| 42 | CC cpp:117 | subdivision tests only vertex distances; interior approaches under-refined | **disposition: coherent/by-design** → §4; documented behaviour (claude #54 agrees mitigation is the L1/L2 ladder). Add no code; §6 question if calibration ever shows it |
| 43 | CC cpp:310 | `d` is center→atom but stored as `direction_to_center`; sign reversed? | **disposition: coherent** → §4 (shared convention across all 4 ring calcs). Internal local renamed via E10; serialized field name unchanged |
| 44 | CC cpp:316 | `abs(z)` folds above/below ring into one theta | **duplicate** of #39 (coherent) |
| 45 | CC cpp:375 | sampler near-field guard uses center distance/radius not surface distance | **disposition: coherent** → §4; explained by sampler-vs-Compute filter difference. Clarified via E9 |
| 46 | CC cpp:81 | `r` hides field point → `field_point` | **adopted** → E12 |
| 47 | CC cpp:90 | `rS` → `surface_point` | **duplicate** of #9 → E2 |
| 48 | CC cpp:285 | `V` → `effective_field` | **duplicate** of #11 → E6 |
| 49 | CC cpp:310 | `d` → `center_to_atom` | **duplicate** of #13 → E10 |
| 50 | CC cpp:315 | scalar `rho` → `in_plane_radius` | **duplicate** of #14 → E11 |
| 51 | CC cpp:316 | `theta` → `unsigned_theta` | **declined** → serialized field is `rn->theta` (matches BS); a renamed *local* that feeds the unchanged field would mislead. The unsigned convention is documented via E8 signpost instead |
| 52 | CC cpp:291 | label outer-product construction | **adopted** → E6 signpost |
| 53 | CC cpp:313 | label ring-frame coordinate derivation | **duplicate** of #5 → E8 |
| 54 | CC cpp:325 | split storage block labels (`raw integral`, `shielding kernel`) | **adopted** → E20 |
| 55 | CC cpp:44 | `Gauss7` → `TriangleGauss7Rule` | **duplicate** of #23 → E14 |
| 56 | CC cpp:79 | `AccumulateTensor` → `AccumulateTriangleH` | **duplicate** of #24 → E15 |
| 57 | CC cpp:125 | `AccumulateAdaptive` → `AccumulateAdaptiveTriangleH` | **duplicate** of #25 → E16 |
| 58 | CC cpp:158 | `SurfaceIntegral` → `ComputeSurfaceIntegralH` | **duplicate** of #26 → E17 |
| 59 | CC cpp:401 | `PackST_HM` → `PackHMSphericalTensor` | **duplicate** of #15 → E13 |
| 60 | CC cpp:32 | "normalised to ref tri area=1/2" wording misleading; weights sum to 1 | **adopted** → E21 (comment fix) |
| 61 | CC cpp:108 | "close to a triangle" overstates the vertex test → `vertex-distance refinement` | **adopted** → E22 |
| 62 | CC cpp:244 | TODO stale/unclear → remove or `source ring record` | **adopted** → E23 (remove TODO) |
| 63 | CC cpp:281 | step comments restate code → terse | **declined** → the Step 1/2/3 labels are the clearest part (claude #18/#32 agree); terse-ing them removes the H/V/G narrative |
| 64 | CC cpp:284 | verbose contraction comment → `normal contraction` | **declined** → minor; the existing comment names the physics (`effective B-field`), which is the point |
| 65 | CC cpp:287 | sign note good but 3 lines → `shielding sign` | **declined** → contradicts brief: this is a grounded sign-convention note that should *stay* (claude #45 agrees "keep"). Compressing it loses the BS-parallel justification |
| 66 | CC cpp:326 | inline `T0=0,T1=0` overclaims exact zeros after FP quadrature | **adopted** → E24 (soften to "≈0 by construction") |
| 67 | CL verdict | header is the model opener; keep | no-op (acknowledged) |
| 68 | CL 234-344 | physics buried among 3 `choices.Record` blocks; add `// --- physics ---` at 281 | **adopted** → E25 |
| 69 | CL 296-323 | find-or-create + polar setup reads as housekeeping; signpost at 296 | **duplicate** of #4 → E7 |
| 70 | CL 312-319 | cylindrical-frame convention only inferable; add signpost | **duplicate** of #5 → E8 |
| 71 | CL 90-96 | `rho`,`rho3`,`rho5` clear from banner; `r` terse but ok | **partially adopted** → `r`→`field_point` via E12; `rho`(vector)→`sep` via E3 |
| 72 | CL 310 | `Vec3 d` bare → `disp`/`to_center` to match `direction_to_center` | **adopted** → E10 (`atom_from_center`) |
| 73 | CL 401 | `PackST_HM`/`out` fine; `_HM` redundant | **noted**; E13 still renames for clarity (codex twice asked) — weighed call |
| 74 | CL 281-294 | Step 1/2/3 grouping perfect, no change | no-op (acknowledged) |
| 75 | CL 99-103 | `K_ab` fused loop acceptable | no-op (acknowledged) |
| 76 | CL 336-341 | per-type accumulation fine | no-op (acknowledged) |
| 77 | CL 79 | `AccumulateTensor` → `AccumulateTriangleIntegral` | **adopted** → E15 |
| 78 | CL 158 | `SurfaceIntegral` good | (codex asks for `H` suffix; adopted via E17 — weighed toward explicitness) |
| 79 | CL 401 | `_HM` suffix noise, "leave it" | see #73 — overridden by codex twice + cross-calc clarity; **adopted** E13 |
| 80 | CL 244 | TODO is process prose; make actionable or drop | **duplicate** of #62 → E23 |
| 81 | CL 287-290 | minus-sign justification good; keep | **adopted** → E1 (and keep cpp note) |
| 82 | CL 326-330 | trailing `Decompose` signposts excellent; keep | **adopted** → E24 only softens "=0"; keeps the signpost |
| 83 | CL h:18-19 | header omits leading minus on G; update to show minus + cite cpp:287 note | **duplicate** of #1 → E1 |
| 84 | CL 316 | `theta` folds both faces onto [0,π/2]; check convention | **duplicate** of #39 (coherent) |
| 85 | CL 337 | `ti>=0 && ti<8` drops ring types ≥8 from per-type; check guaranteed | **disposition: coherent/by-design** → §4 (index 8 = saturated ProPyrrolidine, excluded by `kAromaticRingTypeCount`). Add E26 signpost; no code change |
| 86 | CL 377 vs 230 | sampler `distance<radius` guard not in Compute; confirm filters agree | **disposition: coherent** → §4 (filter set covers it). Clarified via E9 |
| 87 | CL 93 | singularity-guard `continue` drops 1 point, biases H near surface | **disposition: coherent/by-design** → §4 (L1/L2 subdivision is the mitigation, matches CC#42). §6 question only |
| 88 | CL 171 | fan triangulation winding can't flip sign (uses `.norm()`); robust | no-op (acknowledged — confirms no bug) |

## 3. Edits that don't move numbers

Comment, signpost, and internal-name edits only. Line numbers are
pre-edit. No serialized name and no number changes.

- **E1** `HaighMallionResult.h:15-27` — align the G doc-comment with the
  code and BiotSavart's header. Change `G_ab = n_b * (H . n)_a` to
  `G_ab = -n_b * (H . n)_a` (and `G = n (x) V` to `G = -V (x) n`, i.e.
  `G_ab = -V_a n_b`), then add a one-sentence sign note mirroring
  BiotSavartResult.h:23-27: the minus is from `σ_ab = -dB_a^sec/dB_{0,b}`,
  same convention as BS, where `σ = I·G` gives the correct sign with
  literature ring-current intensities. (Resolves #1,2,37,38,40,41,83.)
- **E2** `HaighMallionResult.cpp:90` — `Vec3 rS` → `surface_point`.
- **E3** `HaighMallionResult.cpp:91-103` — vector `Vec3 rho` →
  `sep` (the field-minus-surface separation vector), with the banner's
  `rho = r - r_s` wording kept. Avoids the scalar/vector `rho` clash
  with the `rho` *field* downstream. Update `rho3`,`rho5` derived names
  to `sep3`,`sep5`. File-local; no carry-through.
- **E4** `HaighMallionResult.cpp:83,101` — keep `Mat3& H` (matches header
  `H_ab`, banner, and the dense `H(a,b)+=` kernel loop) — declining #8;
  clarity comes from E5 signpost, not a long name in the hot loop.
- **E5** `HaighMallionResult.cpp:98` — keep the `// K_ab = ...` comment; no
  change (acknowledged clean by both reviews).
- **E6** `HaighMallionResult.cpp:285,381` — `Vec3 V` → `effective_field`
  in both `Compute` and `SampleShieldingAt`; update the `G(a,b) = -... * V(a)`
  references. Add inline `// outer product → rank-1 shielding kernel`
  before the `G(a,b)` loop. (Resolves #3,11,48,52.)
- **E7** `HaighMallionResult.cpp:296` — signpost
  `// locate or create the per-ring record for this atom` above the
  find-or-create loop. (Resolves #4,69.)
- **E8** `HaighMallionResult.cpp:310-319` — signpost
  `// ring-frame coordinates: z along normal, theta unsigned (folds both faces)`
  above the geometry block. This states the unsigned-theta convention in
  words so no rename is needed. (Resolves #5,39,44,53,70,84,51.)
- **E9** `HaighMallionResult.cpp:362,375-378` — change the
  `// Same physics as Compute()` banner to note the grid-point exclusions
  differ: `// Same kernel as Compute(); grid-point exclusions (no bonded/`
  `// atom filters; skips inside-ring-footprint via distance < radius)`.
  (Resolves #6,20,45,86.)
- **E10** `HaighMallionResult.cpp:310` — `Vec3 d` → `atom_from_center`
  (it is `atom_pos - geom.center`). Update `d.normalized()`, `d.dot(...)`,
  `d - z*normal` references. Note: the *field* it feeds,
  `direction_to_center`, is **not** renamed (see §4). (Resolves #13,43,49,72.)
- **E11** `HaighMallionResult.cpp:315` — scalar `double rho` →
  `in_plane_radius`; assign `rn->rho = in_plane_radius` (field unchanged).
  (Resolves #14,50.)
- **E12** `HaighMallionResult.cpp:81` — param `const Vec3& r` →
  `field_point` in `AccumulateTensor` (and the `r` uses at 91, 120-122 in
  `NeedsSubdivision`/`AccumulateAdaptive` for consistency). File-local
  statics; no carry-through. (Resolves #46,71.)
- **E13** `HaighMallionResult.cpp:401,418` — `PackST_HM` →
  `PackHaighMallionTensor`. File-local static; no carry-through. (Weighed
  call: codex asked twice (#15,59); claude said "leave it" (#73,79) — the
  cross-calculator clarity wins for a zero-cost local rename.)
- **E14** `HaighMallionResult.cpp:44,166` — `Gauss7` →
  `TriangleGauss7Rule`. File-local; update the one call site. (Resolves #23,55.)
- **E15** `HaighMallionResult.cpp:79,146` — `AccumulateTensor` →
  `AccumulateTriangleIntegral`. File-local static; update calls.
  (Resolves #24,56,77. Note: chose `...Integral` over codex's `...H` —
  reads cleaner and the banner already says it accumulates into `H`.)
- **E16** `HaighMallionResult.cpp:125,141-144,171` →
  `AccumulateAdaptiveTriangleIntegral`. File-local; update recursion +
  the call in `SurfaceIntegral`. (Resolves #25,57.)
- **E17** `HaighMallionResult.cpp:158,282,380` — `SurfaceIntegral` →
  `ComputeSurfaceIntegralH`. File-local static; two call sites.
  (Resolves #26,58,78.)
- **E18** `HaighMallionResult.cpp:238` — `// ---- GeometryChoice: surface integral ----`
  → `// record this ring as a source (once)`. (Resolves #32.)
- **E19** `HaighMallionResult.cpp:250` — `// Apply filter` →
  `// near-field / bonded pair filters`. (Resolves #33.)
- **E20** `HaighMallionResult.cpp:325` — add two short labels in the
  storage block: `// raw surface integral H` over the `hm_H_*` stores and
  `// rank-1 shielding kernel G` over the `hm_G_*` stores. (Resolves #54.)
- **E21** `HaighMallionResult.cpp:32` — comment "Weights sum to 1.0
  (normalised to reference triangle area = 1/2)" is confusing; reword to
  state plainly: weights sum to 1; the reference-triangle area factor is
  applied separately (triangle area enters via `triArea` in
  `AccumulateTensor`). (Resolves #60.)
- **E22** `HaighMallionResult.cpp:108` — banner "when the field point is
  close to a triangle" → "when the field point is close to a triangle
  *vertex*" to match the vertex-only test. (Resolves #61.)
- **E23** `HaighMallionResult.cpp:244` — remove the stale
  `// TODO: HM sampler needs AccumulateTensor refactor` (it references a
  refactor that doesn't apply to this record block). (Resolves #62,80.)
- **E24** `HaighMallionResult.cpp:326-330` — keep the excellent
  `Decompose(H): T0=0,T1=0` / `Decompose(G): T0,T1,T2` signposts; soften
  the exact-zero claim to `T0≈0, T1≈0 by construction (FP quadrature)`.
  (Resolves #66; preserves #82.)
- **E25** `HaighMallionResult.cpp:281` — add `// --- physics ---` (or
  `// shielding tensor construction`) signpost where the three Record
  blocks end and the H→V→G computation resumes. (Resolves #68.)
- **E26** `HaighMallionResult.cpp:337` — signpost on the `ti < 8` guard:
  `// aromatic ring types only (index 8 = saturated Pro; see kAromaticRingTypeCount)`.
  No code change — `8` stays (it is the NPY ABI shape and matches
  `kAromaticRingTypeCount`). (Resolves #85.)

Declined, for the record: E-less items #16,17,22,27,28,29,30,31,34,63,64,65
(see ledger reasons — mostly "banner carries real signal" or
"contradicts the brief's keep-the-sign-note rule").

## 4. Usage notes (the real product)

### G sign — `G_ab = -n_b * V_a` (the minus is correct)

**Reason found.** The minus is the shielding-tensor definition
`σ_ab = -dB_a^sec / dB_{0,b}`. It is **identical to the BiotSavart
sibling**, which constructs `G(a,b) = -normal(b) * B(a) * PPM_FACTOR`
(`src/BiotSavartResult.cpp:192,235,412`) and documents the rationale at
`BiotSavartResult.h:23-27`, including a numeric verification
(I=-12, 3 Å above PHE → σ = +1.40 ppm shielded). HM's `V = H·n` plays
the role of BS's `B`. So both ring-current shielding kernels share one
sign convention; the inline note at `HaighMallionResult.cpp:287-290`
already says this. The only defect is documentation: the **header
omits the minus** (h:15, h:25), and unlike BS it does not carry the
signed justification paragraph in the header.

**Consumers, and do they agree.**
- `python/nmr_extract/_catalog.py:160-161` tags `hm_shielding` with
  `sign_convention = "σ_ab = -dB_a^sec/dB_{0,b}"` (`_SHIELD_SIGN`,
  line 123) — the SDK contract names exactly this convention. Agrees.
- The same `_SHIELD_SIGN` is shared by `bs_shielding`, `pq_shielding`,
  `disp_shielding`, etc. (lines 149,169,177). Consistent across the
  ring-current/EFG family. Agrees.
- `learn/` consumes `hm_shielding`/`hm_per_type_*` as ridge-regression
  features; calibration fits a coefficient against DFT and is sign-
  agnostic (a global sign flip would be absorbed into the fitted
  coefficient). No consumer asserts a specific sign of G.
- `ui/src/RestServer.cpp:1185` exposes `direction_to_center` but not G's
  raw sign in a way that re-derives shielding; it displays the stored
  tensors.

**Verdict: coherent (expected).** Fix = E1 (header comment only). No
number moves.

### `direction_to_center` field — stores `atom_pos - center` (center→atom)

**Reason found.** Despite the literal name, the field is assigned the
center→atom unit vector in **every** ring calculator that writes it:
`BiotSavartResult.cpp:251-252` (`d = atom_pos - geom.center;
new_rn.direction_to_center = d.normalized();`), `HaighMallionResult.cpp:310-311`,
`DispersionResult.cpp:317`, and via `kernel.direction` in
`PiQuadrupoleResult.cpp:194` / `RingSusceptibilityResult.cpp:191`. It is a
shared `RingNeighbourhood` field (`src/ConformationAtom.h:35`) and is a
serialized/consumed quantity (`ui/src/RestServer.cpp:1185`).

**Verdict: coherent (shared convention).** The internal local `d` is
renamed to `atom_from_center` (E10) to make the *local* unambiguous, but
the **field name is left as-is** — it is a contract shared by four
calculators and a UI consumer; renaming it would be a cross-file/cross-
subproject change for no physics gain, and HM is not the place to do it
unilaterally. Noted in §6 as a possible repo-wide naming question.

### `theta = atan2(d_plane.norm(), abs(z))` — intentionally unsigned

**Reason found.** Identical to BiotSavart (`BiotSavartResult.cpp:258`).
`theta` is the polar angle off the ring axis folded to [0, π/2] by
`abs(z)`; the signed face information is retained in `z` (stored
separately as `rn->z`). The two together (`z` signed, `theta` unsigned)
fully describe the position; nothing downstream needs a signed theta
that isn't already covered by `z`. The `rho`/`theta`/`z` triple is the
shared ring-frame parameterization.

**Verdict: coherent (expected).** Fix = E8 signpost stating the
convention in words; no rename of the field (declining CC#51's
`unsigned_theta` because the field is `rn->theta`, shared with BS).

### `ti < 8` per-type guard — drops only saturated rings, by design

**Reason found.** `RingTypeIndex` (`src/Types.h:189-206`) has indices
0–7 aromatic and index 8 = `ProPyrrolidine`, a *saturated, non-aromatic*
5-ring (aromaticity None, ring-current intensity 0). `kAromaticRingTypeCount
= 8` (`Types.h:233`) is exactly that boundary, with static_asserts pinning
it and a long comment (Types.h:214-238) explaining that calculators still
use the literal `8` because it is also the NPY ABI shape (`std::array<double,8>`
at `ConformationAtom.h:132-133`). Pro rings are still added to `G_total`
(line 333) — they're only excluded from the *per-aromatic-type* breakdown,
which is correct: a saturated ring has no aromatic ring current.

**Verdict: coherent/by-design.** Fix = E26 signpost; `8` stays (changing
it is a coordinated NPY schema migration, explicitly out of scope).

### Sampler `distance < radius` guard vs Compute filters

**Reason found.** `SampleShieldingAt` (grid points, no atoms) cannot use
the atom-keyed `KernelFilterSet` that `Compute` uses
(`DipolarNearFieldFilter`, `RingBondedExclusionFilter`,
`MinDistanceFilter`), so it substitutes inline guards: singularity
distance, inside-ring-footprint (`distance < radius`), and the spatial
cutoff. `Compute`'s `DipolarNearFieldFilter` with `source_extent =
2*radius` covers the near-field/inside-ring case for real atoms. The two
paths target the same physics with path-appropriate exclusions.

**Verdict: coherent.** Fix = E9 (clarify the banner); no code change.

### Vertex-only subdivision test / dropped-quadrature-point bias

**Reason (partial).** `NeedsSubdivision` tests only the three triangle
vertices (cpp:120-122), and `AccumulateTensor` `continue`s past a single
near-coincident quadrature point (cpp:93) while keeping the other six.
Both reviews (CC#42, CL#54/#87) note the L1/L2 adaptive ladder
(2.0 Å → 1.0 Å, up to 112 points/fan triangle) is the intended
mitigation, and CL#88 confirms winding cannot flip the sign. This is
documented behaviour, not a discrepancy with any consumer.

**Verdict: coherent/by-design** for this readability pass — no code
change. Listed in §6 as a numerical-accuracy question for the standing
"PATTERNS.md numerical-stability hygiene pass," not a bug.

## 5. Bug-by-exhaustion candidates

**None.** No finding survived tracing as a producer- or consumer-side
bug. The two sign/convention items both reconcile against the
BiotSavart sibling and the SDK `_SHIELD_SIGN` contract; the per-type
guard, unsigned theta, and `direction_to_center` name are all shared
conventions with documented rationale. The numerical items (vertex-only
subdivision, dropped quadrature point) are documented design with a
named mitigation and no consumer that depends on the unrefined value.

## 6. Questions & ambiguities

1. **`direction_to_center` field name (repo-wide).** The field literally
   reads "direction to center" but stores center→atom in all five ring
   calculators plus the UI. This is internally consistent but the *name
   inverts the stored sense*. Worth a one-line confirmation that the name
   is intended as "the center-to-atom axis of this neighbourhood" (i.e.
   read it as the ring-frame axis label, not a from→to direction). If a
   repo-wide rename is ever wanted (`axis_from_center`?), it must land
   across `ConformationAtom.h` + 5 calculators + `RestServer.cpp` + any
   H5/SDK consumer together — out of scope for this single-file pass.
   **Not changing it here.**

2. **Vertex-only subdivision adequacy (numerical, not readability).**
   `NeedsSubdivision` keys on vertex distance, so an atom close to a
   triangle *edge or interior* but far from all three vertices could be
   under-refined. Reviews agree the L1/L2 ladder is the intended
   mitigation. Confirmation that fan-triangle geometry (centroid +
   consecutive vertices, max ~ring-radius edges) makes the edge/interior
   case unreachable in practice would close CC#42 / CL#87 fully. This is
   a calibration/accuracy question, deferred to the numerical-hygiene
   pass — no code change proposed now.

3. **`SampleShieldingTensorAt` rename (codex h:54).** Declined here as a
   cross-file public-method rename whose clarity gain is marginal
   (the name already says it samples shielding). Flagging in case the
   human wants the explicit `...TensorAt` form; if so it needs its
   `SampleShieldingAt` call sites updated in lockstep (grid/overlay
   consumers).
