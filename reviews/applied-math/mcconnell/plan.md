# McConnell readability fix plan

Targets: `src/McConnellResult.{h,cpp}`, `src/MopacMcConnellResult.{h,cpp}`.
The Mopac variant is the same kernel with each bond's contribution scaled by
its MOPAC Wiberg bond order; findings that apply to one apply to both unless
noted.

---

## 1. Summary

Both files tell a coherent story. The headers are exemplary: they state the
full McConnell tensor `M_ab`, the three-term irrep attribution (T1/T0/T2), the
symmetric-traceless dipolar kernel `K`, and the McConnell scalar `f` with the
explicit "f is not the trace" caveat. The `Compute` arc is honest:
setup → per-bond kernel → filter → category accumulate → per-atom finalize →
spherical features.

The one genuine story-break is shared by both files and is a **comment/naming
vs code** mismatch, not a math error: the per-category `T2_*_total` outputs are
built by symmetrizing + de-tracing the accumulated **full `M`** category sums,
but the comment calls them "from symmetric dipolar kernel sums" and the local
matrices are named `K_backbone` / `K_sidechain` / `K_aromatic`, colliding with
`K` (the dipolar kernel) used elsewhere in the file. The code is right (the
catalog declares these outputs `irreps="2e"`, i.e. symmetric-traceless rank-2,
and the symmetric-traceless part of full `M` is exactly that); the comment and
the local names are the inaccuracy.

The pass will: fix comments to match code, add a handful of 2–4-word
signposts, rename the colliding local `K_*` category matrices, remove a dead
lambda, and name a couple of schema magic numbers. It will **not** touch any
output name (NPY field, struct field consumed by ui/h5-reader, catalog key),
will **not** change any number, and will **not** restructure the loop.

Two sign/value items the reviews flagged both land on **coherent (expected)**
after tracing — see §4. No bug-by-exhaustion candidate survived (see §5).

---

## 2. Review-finding ledger

### codex.md

| # | Finding | Disposition |
|---|---------|-------------|
| C1 | h: public query names less clear than the physics comment (scalar-vs-tensor) | adopted → §3 E12, E13 (header method renames) |
| C2 | cpp:267-285: comment "symmetric dipolar kernel sums" vs code de-traces full `M` | adopted → §3 E1; usage note §4.1 |
| C3 | Mopac h: slightly over-explanatory / overconfident about zero bond order | adopted (partial) → §3 E2 (tighten zero-bo comment); decline trimming the whole docblock — see note |
| C4 | Mopac cpp:151-153: bond-order gate scopes "nearest"; not obvious | adopted → §3 E3 (`// bond-order gate` signpost) |
| C5 | cpp:86-92: dense fused 3-term `M` expression; split into named terms or label | adopted (label form) → §3 E4; decline splitting into locals — see note |
| C6 | cpp:198-240: one switch mixes scalar sums / tensor totals / nearest tracking | declined → the switch is one cohesive "accumulate this bond into its category" step; signposting each `case` interior fragments coherent code. The per-line role is already legible (`co_sum += …`, `M_backbone_total += …`, `if (dist < best…)`). |
| C7 | cpp:45-50: `M_over_r3`, `K`, `f` rescued only by comments → spell out | adopted (partial): `f` → see §4.3 (literature symbol, keep, field-comment already carries it); `M_over_r3` keep (accurate, and the header defines it); `K` keep at struct level (it IS the dipolar kernel and the header names it). The real `K` collision is the category locals — fixed in E5. |
| C8 | cpp:61: `d` hides the frame → `midpoint_to_atom` | adopted (lighter) → §3 E6 (`disp` + one-line comment); full `midpoint_to_atom` declined as verbose for a 3-line scope |
| C9 | cpp:142: `sidechain_sum` only accumulates SidechainCO → `sidechain_co_sum` | declined → it pairs with output field `mcconnell_sidechain_sum` and catalog `sidechain_sum`; renaming the local away from its destination field reduces, not improves, traceability. Noted in §4.4. |
| C10 | cpp:143-146: category totals already `/r^3` → `M_over_r3_total` | declined → churn; the `M_*_total` names are accurate (they hold the full tensor sum) and the `/r^3` is carried by the kernel struct field name `M_over_r3`. |
| C11 | cpp:153,206: `best_co_direction` stores `kernel.direction` = midpoint→atom, not bond dir → `best_co_midpoint_to_atom` | adopted as **comment**, not rename → §3 E7; the sign/sense question is the real issue, resolved in §4.2 (coherent). |
| C12 | cpp:275/279/283: `K_backbone` etc. collide with `K` = dipolar kernel | adopted → §3 E5 (rename to `backbone_T2_source` etc.) |
| C13 | Mopac cpp:152-153: `bo`/`zero_bo_skipped` terse, threshold-based not strictly zero | adopted (partial) → §3 E2/E3 (comment says "floor", split the packed line); `bo`→keep (local 3-use scope, comment names it); `zero_bo_skipped`→keep (matches the GeometryChoice/log key emitted, see §4.5) |
| C14 | cpp:74-92: full-tensor block is the only one not decomposed into named terms | duplicate of C5 |
| C15 | cpp:269-272: `project_traceless` defined but unused, op repeated inline | adopted → §3 E8 (delete dead lambda) |
| C16 | cpp:372-398: feature widths 9/25/6 are schema magic numbers | adopted → §3 E9 (named constants) |
| C17 | Mopac cpp:303-333: same schema magic numbers | adopted → §3 E9 (both files) |
| C18 | h:51 `CategorySum` doesn't say it returns the scalar sum → `CategoryScalarSum` | adopted → §3 E12 |
| C19 | h:52 `NearestCOContribution` hides scalar `f` → `NearestCOScalarContribution` | adopted → §3 E13 |
| C20 | h:55 `SampleShieldingAt` clearer as `SampleMcConnellShieldingAt` | declined → in-class the type already qualifies it (`McConnellResult::SampleShieldingAt`); sibling calcs use the same `SampleShieldingAt` name (BiotSavart), so the bare name is the cross-calc convention. Noted §6 Q3. |
| C21 | cpp:54 / Mopac cpp:47: `ComputeBondKernel` generic → `ComputeMcConnellBondKernel` | declined → both are `static` (internal linkage); renaming is harmless but the duplicate-name-across-files concern is exactly what claude C-J4 calls acceptable. Low value. |
| C22 | cpp:362 / Mopac cpp:293: `PackST_MC`/`PackST_MMC` cryptic → `PackSphericalTensor9` | adopted → §3 E10 |
| C23 | cpp:119-124: filter comment omits `MinDistanceFilter` | adopted → §3 E11 |
| C24 | cpp:172: `// ---- GeometryChoice: filter exclusion ----` process-heavy | declined → the `GeometryChoice:` prefix is a deliberate house tag marking provenance-recording blocks across all calculators (grep-anchor); trimming it de-conventionalizes one site. |
| C25 | cpp:267-268: comment potentially misleading → `// T2 projection of M` | duplicate of C2 |
| C26 | cpp:330-331: comment says `MCCONNELL_CUTOFF_A`, code uses config value | adopted → §3 E15 (and ties to dead-constant note §4.6) |
| C27 | cpp:346: "skip if inside the bond" less exact → `// near-field cutoff` | adopted → §3 E16 |
| C28 | Mopac h:13-16: explanatory prose longer than needed | declined → the docblock states the bond-order physics ("measured QM quantity, not a parameter; order 1.8 contributes more") which is the whole point of the variant; both reviews elsewhere praise these headers as exemplary. Trimming loses the delta-explanation. |
| C29 | Mopac cpp:94-96: "electronically insignificant" is interpretation; code uses a noise floor | adopted → §3 E2 |
| C30 | **Correctness** cpp:189: `direction_to_midpoint` receives `kernel.direction` (midpoint→atom); if field means atom→midpoint this is a sign error | adopted as investigation → §4.2: **coherent**, display-only consumer, field name is the inaccuracy → comment fix E7/E17, not a sign fix |
| C31 | **Correctness** cpp:267-285 / Mopac cpp:251-262: comments imply dipolar `K` sums, code de-traces full `M` totals | adopted as investigation → §4.1: **coherent**, catalog confirms `2e` intent → comment fix E1 |

### claude.md

| # | Finding | Disposition |
|---|---------|-------------|
| J1 | cpp:269-285: `project_traceless` lambda defined then never called | duplicate of C15 → §3 E8 |
| J2 | cpp:275-285: symmetrize+de-trace idiom repeated 3× verbatim; shared step visible only by reading all three | adopted (signpost only) → §3 E1; decline extracting a helper (no new abstraction per brief; the inline form is the model the file already uses elsewhere) |
| J3 | cpp:257-265: nearest-CO/CN finalize mixes dir-renorm guard, CO T2, CN T2 in asymmetric branches, unsignposted | adopted → §3 E18 (`// nearest CO/CN: direction + T2`) |
| J4 | Mopac cpp:153: three statements on one line compress the skip | adopted → §3 E3 (split lines) |
| J5 | cpp:64 / Mopac cpp:57: early `return result` returns zero kernel; sentinel meaning only via struct defaults | adopted → §3 E19 (`// too close: zero kernel, distance stays 0`) |
| J6 | cpp:61-62 / Mopac cpp:54-55: `d`/`r` one-letter; header names them, code doesn't | adopted (partial) → §3 E6 (`d`→`disp`); `r` kept (universal radial-distance symbol, immediately followed by `r3`, `d_hat`) |
| J7 | cpp:75/37: `f` true one-letter; defensible literature symbol; struct field could be `f_scalar` | declined the field rename → see §4.3; `f` is the catalogue symbol and the struct comment carries it. Decline gently. |
| J8 | cpp:185 `bn`/`BondNeighbourhood` clear; no action | no action (agreed) |
| J9 | Mopac cpp:152 `bo` terse but comment names it; acceptable | no action (agreed) — see C13 |
| J10 | cpp:78-94: K and M built in two adjacent loops, each preceded by its formula comment — good, the model for the rest | no action (praise) |
| J11 | cpp:275-288: category T2 + full-M decompose are 4 blocks with no grouping header | adopted → §3 E20 (blank line + one signpost separating "category T2" from "full shielding") |
| J12 | Mopac cpp:180-183: `weighted_M`/`weighted_f` make the weighting visible — clean | no action (praise) |
| J13 | cpp:54 / Mopac cpp:47: two same-named static `ComputeBondKernel`; acceptable, no rename | no action (agreed) — see C21 |
| J14 | cpp:362/293: `_MC`/`_MMC` suffixes disambiguate; clear enough | partial conflict with C22 → adopt the C22 rename `PackSphericalTensor9` (both files; the suffix disambiguation is unnecessary since they're file-local static) → §3 E10 |
| J15 | h:55 `SampleShieldingAt` says what it returns; good | conflicts with C20 → side with J15, decline rename (see C20) |
| J16 | headers exemplary: formula, irrep attribution, f-is-not-trace caveat | no action (praise) — these set the standard the .cpp comments should meet |
| J17 | cpp:268: comment describes the *unused* lambda, not the inline code; stale | duplicate of C2/J1 → §3 E1/E8 |
| J18 | cpp:274: "Extract symmetric part" omits the de-trace; partial | adopted → §3 E1 (`// symmetric traceless part for T2`) |
| J19 | cpp:119-120: filter-set comment accurate and useful | conflicts with C23 — both can hold: comment is useful but omits `MinDistanceFilter`; adopt the C23 addition → §3 E11 |
| J20 | cpp:184 "Store in BondNeighbourhood" restates next 5 lines; low value, harmless | declined → harmless; leave (removing it is churn, keeping it is fine). Noted. |
| J21 | Mopac cpp:94-97: zero-order skip docblock clear; good | partial conflict with C29 — the docblock is clear but "electronically insignificant" overstates a noise-floor test; adopt the C29 wording fix → §3 E2 |
| J22 | **Correctness** h:37 `MCCONNELL_CUTOFF_A=10.0` / Mopac h:35 unused; actual cutoff from config; can silently diverge | adopted → §4.6 + §3 E14 (mark dead / comment), §3 E15 (fix the comment that names the dead constant) |
| J23 | **Correctness** cpp:202 vs 214: CO branch tracks f/midpoint/direction, CN tracks only dist+kernel; verify asymmetry is by design | adopted as investigation → §4.7: **coherent** (output set has no `nearest_CN_f`/midpoint/dir); comment signpost E18 makes it visible |
| J24 | **Correctness** cpp:257 sentinel guard + :259 near-zero dir guard sound | no action (agreed) |
| J25 | **Correctness** Mopac cpp:242-248: weight reapplied at decompose, not double-applied; OK on inspection | adopted as investigation → §4.8: **coherent** (confirmed single application) |
| J26 | **Correctness** both: `M_total` decomposed directly (keeps T0+T1+T2); category paths de-trace; consistent with header | adopted as investigation → §4.1 (this is the key that resolves C2/C31): **coherent by design** |
| J27 | **Correctness** cpp:343-348 `SampleShieldingAt` re-implements singularity/cutoff/near-field inline rather than via filter set; must stay in sync with `DipolarNearFieldFilter` | adopted as investigation → §4.9: **can't fully confirm filter equivalence** → §6 Q1 + §3 E16 cross-ref comment |

---

## 3. Edits that don't move numbers

Comment fixes, signposts, internal renames, dead-code removal. No output name
changes. Line numbers are at current HEAD.

**McConnellResult.cpp**

- **E1** `McConnellResult.cpp:267-268, 274` — replace the two stale/partial
  comments. `:267` "Category T2 totals (from symmetric dipolar kernel sums)"
  → `// Category T2: symmetric-traceless part of each category's full-M sum`.
  `:268` "Apply traceless projection to fix floating-point drift" → delete
  (it describes the dead lambda; see E8). `:274` "Extract symmetric part for
  T2 features" → `// symmetrize, then subtract trace/3 → symmetric-traceless`.
- **E4** `McConnellResult.cpp:84-92` — keep the fused expression (it is a
  direct transcription of the header formula) but add inline term labels on
  the three summed lines: `// 9 cosθ d̂⊗b̂ (T1)`, `// -3 b̂⊗b̂ (T0)`,
  `// -(3 d̂⊗d̂ - I) (T2)`, mirroring the header's term table.
- **E5** `McConnellResult.cpp:275,279,283` — rename locals `K_backbone`,
  `K_sidechain`, `K_aromatic` → `backbone_T2_source`, `sidechain_T2_source`,
  `aromatic_T2_source` (or `*_symtraceless`). Internal, file-local; no
  consumer. Resolves the `K`-means-dipolar-kernel collision.
- **E6** `McConnellResult.cpp:61` (and Mopac `:54`) — rename `d` → `disp`
  with a trailing `// midpoint→atom displacement`. Keep `r`, `r3`, `d_hat`.
- **E7** `McConnellResult.cpp:206` — add `// kernel.direction is d̂ =
  midpoint→atom` next to `best_co_direction = kernel.direction;` (the
  field/name says "to midpoint"; the value is the opposite sense — see §4.2).
- **E8** `McConnellResult.cpp:269-272` — delete the unused `project_traceless`
  lambda. The inline blocks at :275-285 already do the operation.
- **E9** `McConnellResult.cpp:372-398` (and Mopac `:303-333`) — name the
  feature widths: `constexpr int kShieldingCols = 9; kCategoryT2Cols = 25; //
  5 categories × 5 T2; kScalarCols = 6;` used in the `vector` sizing and the
  `NpyWriter` calls. Pure readability; values unchanged.
- **E10** `McConnellResult.cpp:362` (and Mopac `:293`) — rename `PackST_MC` /
  `PackST_MMC` → `PackSphericalTensor9` (file-local static, no external
  consumer). Resolves both C22 and J14.
- **E11** `McConnellResult.cpp:119-120` — extend the filter-set comment to name
  all three: `// filters: MinDistance + SelfSource (atom is bond endpoint) +
  DipolarNearField (source extent = bond length)`.
- **E12** `McConnellResult.h:51` — rename query method `CategorySum` →
  `CategoryScalarSum` (returns the McConnell scalar `f` sum, not a tensor).
  **Carry-through:** grep shows no caller outside the class except tests;
  confirm and update any `.CategorySum(` call sites in `src/`/`tests/` in the
  same edit (see §6 Q2 for the carry-through scope I could not fully close).
- **E13** `McConnellResult.h:52` — rename `NearestCOContribution` →
  `NearestCOScalarContribution`. Same carry-through caveat as E12.
- **E15** `McConnellResult.cpp:330-331` — comment says "within
  MCCONNELL_CUTOFF_A"; code uses `CalculatorConfig::Get("mcconnell_bond_
  anisotropy_cutoff")`. Replace with `// within configured bond cutoff`.
- **E16** `McConnellResult.cpp:346` — replace "skip if inside the bond" with
  `// DipolarNearField: skip if closer than near_field_exclusion_ratio ·
  bond_len` and add `// keep in sync with DipolarNearFieldFilter` (see §6 Q1).
- **E18** `McConnellResult.cpp:257` — add signpost `// nearest CO/CN: renorm
  direction (CO only) + T2 decompose` above the two finalize branches (makes
  the documented asymmetry — §4.7 — visible).
- **E19** `McConnellResult.cpp:64` (and Mopac `:57`) — append to the
  singularity-guard return: `// too close: return zero kernel (distance stays
  0, bond contributes nothing)`.
- **E20** `McConnellResult.cpp:286-288` — add a blank line + signpost
  separating the category-T2 outputs from `// Full McConnell shielding
  contribution (full M: T0+T1+T2)` so the two output families are distinct.

**MopacMcConnellResult.cpp**

- **E2** `MopacMcConnellResult.cpp:94-96` — replace "they are electronically
  insignificant" with `// below the bond-order noise floor → skip`. Tightens
  C29/J21.
- **E3** `MopacMcConnellResult.cpp:151-153` — add `// bond-order gate: scopes
  every "nearest"/sum below to bonds above the noise floor` above the
  `bo`-fetch, and split the packed skip line onto three lines (J4).
- Mirror **E1, E4, E5, E6, E7, E9, E10, E18, E19, E20** in the Mopac file at
  its corresponding lines (`:251-262` category comment; `:79-83` term labels;
  `:252/256/260` `K_*` renames; `:54` `d`; `:189-192` direction comment;
  `:303-333` schema constants; `:293` Pack rename; `:242` nearest signpost;
  `:57` zero-kernel comment; `:264` full-shielding signpost).

**MopacMcConnellResult.h**

- **E2h** `MopacMcConnellResult.h` — leave the docblock (C28 declined); no
  change beyond keeping it consistent with E2's "noise floor" wording if the
  header restates it (it does not — header says "0.0", which is the sentinel,
  fine).

**Headers — dead constant**

- **E14** `McConnellResult.h:37` `MCCONNELL_CUTOFF_A` and
  `MopacMcConnellResult.h:35` `MOPAC_MCCONNELL_CUTOFF_A` — these constants are
  **unused by the runtime path** (the cutoff comes from config). Add
  `// NOTE: runtime uses CalculatorConfig "..._bond_anisotropy_cutoff"; this
  constant is referenced only in comments/docs.` Do **not** delete without
  confirming no other TU references them (grep showed none outside these files
  and their comments — see §4.6). Flagged for the human to decide deletion.

---

## 4. Usage notes (the real product)

### 4.1 Category `T2_*_total` are symmetric-traceless projections of full `M`, not sums of `K` — **coherent (expected)**

**Reason found.** Two distinct tensors flow through `Compute`:
- `kernel.M_over_r3` — the *full* McConnell tensor `M_ab/r³` (asymmetric,
  non-traceless; carries T0+T1+T2 per the header's term table).
- `kernel.K` — the symmetric-traceless dipolar kernel
  `(3 d̂d̂ − I)/r³`, whose information lives entirely in T2.

The per-category accumulators `M_backbone_total` / `M_sidechain_total` /
`M_aromatic_total` sum the **full `M`**. The finalize step (`:275-285`)
symmetrizes (`0.5(M+Mᵀ)`) and de-traces (`− trace/3·I`), yielding the
symmetric-traceless part of those full-M sums, then `SphericalTensor::Decompose`
extracts the 5 T2 components emitted to `mc_category_T2`.

**Why this is correct, not a `K`-sum.** The catalog declares both
`mc_category_T2` and `mopac_mc_category_T2` with `irreps="2e"`
(`python/nmr_extract/_catalog.py:190,288`) — a pure symmetric-traceless rank-2
object. The SDK wrapper `PerBondCategoryT2`
(`python/nmr_extract/_tensors.py:316`) reads exactly 5 components per category.
The full M, restricted to its symmetric-traceless part, *is* a "2e" object, and
it is the physically right thing: the category total is the angular shielding
pattern of all bonds in that category, with the (physical) trace and
antisymmetric pieces dropped because the category feature is defined as T2-only.
Summing per-bond `K` instead would give a numerically different result
(K omits the `9cosθ d̂⊗b̂` and `−3 b̂⊗b̂` terms before projection); the chosen
path keeps those terms' T2 content. claude J26 independently confirms the
full-M path (`mc_shielding_contribution`, decomposed directly at `:288`) keeps
T0+T1+T2, consistent with the header.

**Producer/consumers agree.** No consumer treats `mc_category_T2` as a `K`-sum;
all read it as a 5-vector T2 feature (`learn/bones/.../dataset.py:55-57`,
`features.py:219`). The mismatch is purely in the .cpp comment wording and the
local `K_*` names — fixed by E1 + E5.

### 4.2 `direction_to_midpoint` / `best_co_direction` store midpoint→atom (`d̂`), not atom→midpoint — **coherent (expected); field name is the inaccuracy**

**Reason found.** `kernel.direction = d̂ = (atom_pos − bond_midpoint)/r`
(`McConnellResult.cpp:61,70`). This points **from the midpoint toward the
atom**. It is assigned to `BondNeighbourhood::direction_to_midpoint`
(`:189`) and to `best_co_direction`→`dir_nearest_CO` (`:206,260`), both whose
names suggest the opposite sense ("to midpoint", i.e. atom→midpoint).

**Why it is not a sign bug.** The dipolar physics that uses direction —
`f = (3cos²θ − 1)/r³` and `K = (3 d̂d̂ − I)/r³` — is **invariant under
d̂ → −d̂** (both depend only on `d̂d̂`/`cos²θ`). So the stored sense cannot
change any computed shielding number regardless of which way it points.

**Consumers are display-only.** Traced every consumer of these fields:
- `direction_to_midpoint`: only `ui/src/RestServer.cpp:1260` (JSON key
  `"direction_to_midpoint"` for the atom dump) and the `MainWindow` inspector
  tree. Display only.
- `dir_nearest_CO` / `nearest_CO_midpoint`: `h5-reader/notes/
  H5_FIELD_GLOSSARY.md:570,641,1309` render them as a *line segment* on
  selection (not a magnitude, not fed to physics).

No consumer relies on the directional sense for math. **Disposition:** the
*name* is misleading but the *value* is coherent and harmless. Fix with a
comment (E7/E17) stating the sense; do **not** flip the vector (would silently
change the ui/h5-reader glyph orientation for no physical reason) and do **not**
rename `direction_to_midpoint` (it is a struct field consumed across ui +
h5-reader glossary; a rename is a cross-subproject contract change with no
physics payoff — see §6 Q4 for the cost note).

### 4.3 `f` (McConnell scalar) — keep the one-letter name — **coherent**

`f = (3cos²θ − 1)/r³` is the literature symbol, defined verbatim in the header
(`McConnellResult.h:21-23`) and in `GEOMETRIC_KERNEL_CATALOGUE.md:77-85`, with
the explicit "NOT the trace" caveat. The struct field comment carries it. Both
reviews call the rename optional/defensible. Renaming to `f_scalar` adds nothing
a chemist needs and de-aligns from the catalogue. Declined.

### 4.4 `sidechain_sum` / `*_sum` locals — keep — **coherent**

`sidechain_sum` accumulates only `SidechainCO` `f` values, then lands in
`ca.mcconnell_sidechain_sum` (output field) and `mc_scalars[2]` (catalog
`sidechain_sum`). The local name matches its destination field; renaming it to
`sidechain_co_sum` (codex C9) would make the local diverge from the output it
fills, hurting traceability. The "CO-only" scope is already visible from the
single `case BondCategory::SidechainCO:` that feeds it.

### 4.5 `zero_bo_skipped` — keep the name — **coherent**

`zero_bo_skipped` is logged (`MopacMcConnellResult.cpp:279`) and the per-atom
count is emitted as a GeometryChoice number keyed `"zero_bo_skipped"` (`:273`).
Renaming the local to `below_floor_pairs_skipped` (codex C13) would diverge from
the emitted key. The "floor not strictly zero" nuance is handled by fixing the
*comment* (E2) to say "noise floor", which is the accurate part.

### 4.6 `MCCONNELL_CUTOFF_A = 10.0` is dead — **coherent (vestigial), flag only**

The runtime cutoff is `CalculatorConfig::Get("mcconnell_bond_anisotropy_
cutoff")` (`:139`) and the Mopac analogue (`:129`). The named constants are
referenced only in their own headers and in two .cpp comments (one of which,
`:331`, misstates that they are used — fixed by E15). They can silently diverge
from the config value. Not a bug today (nothing reads them in the math). E14
marks them; deletion is the human's call (grep found no other TU referencing
either symbol).

### 4.7 CO vs CN nearest-tracking asymmetry — **coherent (by design)**

The CO branch records `best_co_f`, `best_co_midpoint`, `best_co_direction`,
`best_co_kernel`; the CN branch records only `best_cn_dist` + `best_cn_kernel`.
This matches the output set: `mc_scalars` has `nearest_CO_dist` *and*
`nearest_CN_dist` but only a CO scalar (`mcconnell_co_nearest`), CO midpoint
(`nearest_CO_midpoint`), and CO direction (`dir_nearest_CO`) — there are no
`nearest_CN_f` / midpoint / direction outputs. Both nearest CO and CN emit a
T2 (`T2_CO_nearest`, `T2_CN_nearest`). So the branch asymmetry exactly mirrors
the output asymmetry; it is omission-by-design, not a missed field. E18
signposts it.

### 4.8 Mopac nearest-CO/CN weight applied once — **coherent**

`best_co_kernel` stores the *unweighted* `kernel` (`:192`); the bond order is
reapplied at decompose time as `best_co_bo * best_co_kernel.K` (`:244`).
`best_co_f_weighted` separately holds the already-weighted scalar (`:191`).
Traced: the weight is applied exactly once to each emitted quantity (once to
the scalar at accumulate, once to the T2 at decompose). No double-weighting.
claude J25 reached the same conclusion.

### 4.9 `SampleShieldingAt` inline filters vs `DipolarNearFieldFilter` — see §6 Q1

`SampleShieldingAt` (`:334-353`) re-implements the singularity guard, the
cutoff, and the near-field exclusion (`near_field_exclusion_ratio · bond_len`)
inline, where `Compute` uses the `KernelFilterSet`. The inline near-field test
*looks* equivalent to `DipolarNearFieldFilter`, but I did not read that
filter's source to confirm byte-for-byte equivalence. Flagged Q1; E16 adds a
"keep in sync" comment rather than asserting equivalence.

---

## 5. Bug-by-exhaustion candidates

**None.** Both correctness-flavored review items (C30/C31, J22-J27) resolved to
**coherent** after tracing the catalog, SDK wrappers, ui consumers, and
h5-reader glossary (§4.1, §4.2, §4.6, §4.7, §4.8). The remaining uncertainty
(§4.9 filter equivalence) is a **question**, not a bug — the two paths may
agree; I have not exhausted `DipolarNearFieldFilter`'s definition.

---

## 6. Questions & Ambiguities

- **Q1 — `SampleShieldingAt` vs `DipolarNearFieldFilter` equivalence
  (claude J27).** The grid-sampling path open-codes the near-field exclusion
  as `distance < near_field_exclusion_ratio · bond_len`. Does
  `DipolarNearFieldFilter::Accept` encode the *same* test (same ratio, same
  `<` vs `<=`, same source_extent = bond_len)? If they diverge, grid samples
  and per-atom values disagree near bonds. I did not read the filter source.
  Confirm before relying on grid/per-atom agreement; until then E16 adds a
  "keep in sync" comment only.

- **Q2 — carry-through of the `CategorySum`/`NearestCOContribution` renames
  (E12/E13).** These are public methods on `McConnellResult`. A repo-wide grep
  for call sites was not exhaustively run for this plan; before renaming,
  confirm the only callers are within `src/` + `tests/` (no ui/h5-reader/SDK
  caller, since those read the NPY/struct fields, not these query methods). If
  a caller exists outside `src/`, weigh the rename against the carry-through.

- **Q3 — `SampleShieldingAt` naming convention.** I declined the
  `SampleMcConnellShieldingAt` rename (C20) because `BiotSavartResult` uses the
  same bare `SampleShieldingAt`. Confirm this is the intended cross-calculator
  convention (each calc exposes `SampleShieldingAt`) so the bare name is
  correct and not just inertia.

- **Q4 — `direction_to_midpoint` field rename (deferred, not proposed).**
  The field name is misleading (stores midpoint→atom; §4.2). A rename to
  `midpoint_to_atom_dir` would be clearer, but the name is a contract across
  `ui/src/RestServer.cpp` (JSON key) and the `h5-reader` field glossary. Since
  the value is physics-irrelevant (direction-sign-invariant kernel) and
  display-only, the cost of the cross-subproject carry-through outweighs the
  benefit. Proposed as a comment fix (E7) instead. Flag for the human in case
  the JSON key churn is acceptable.
