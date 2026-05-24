# Fix plan — RingSusceptibilityResult.{h,cpp}

## 1. Summary

The file tells a coherent story. The header states the physical model
(McConnell with `b̂ → n̂` ring normal), the tensor form, units, and the
asymmetric / non-traceless warning. `ComputeRingChiKernel` is scientifically
legible: scalar `f`, dipolar kernel `K`, full tensor `M_over_r3`, each grouped
and labelled to match a banner derivation that matches
`GEOMETRIC_KERNEL_CATALOGUE.md:298-343` (the canonical McConnell derivation,
re-used here with `n̂` in place of `b̂`). `Compute` walks atom → nearby rings
→ filter → store-on-neighbourhood → accumulate in a sensible order.

The friction both reviews land on is real but **local and non-numeric**:
- `theta` (`atan2(rho, |z|)`) reads like a cylindrical azimuth but is the
  polar angle measured from the ring normal, hemisphere-folded. Verified
  coherent against every consumer (see §4). Comment/label fix only.
- A handful of cryptic internal names (`rn`, `d` reused in `Compute`,
  `PackST_RS`, the bare `M_over_r3`/`K`/`f` triple) — internal, improvable.
- The middle of the atom loop (find-or-create `RingNeighbourhood`) is data
  bookkeeping interleaved with the physics; signposts help.

**This pass will touch:** comments (precision, not addition), 2–4-word block
signposts inside the loop, named Kronecker-δ, one file-static rename
(`PackST_RS`), and several local renames. **It will not touch:** any number,
any sign, the algorithm, the filter set, output NPY names/layout, or any
shared struct field. One cross-file rename (`direction_to_center`) and one
cross-file name-soup fix are *proposed with carry-through cost noted* and
left for the human to weigh — not adopted here, because they cross five
producers plus a consumer.

---

## 2. Review-finding ledger

Every finding from `codex.md` and `claude.md`, one row each. (No
`codex-correctness.md` present.)

### codex.md

| # | Finding (loc) | Disposition |
|---|---|---|
| C1 | h:13 — API names don't say full-tensor vs decomposed feature | **Declined** — `SampleShieldingAt` returns a `SphericalTensor` (its type says decomposed); `Compute`/`WriteFeatures` follow the project's uniform `*Result` convention. Renaming `SampleShieldingAt`→`SampleShieldingTensorAt` is a public method change for marginal gain; see C25/C27. |
| C2 | cpp:45 — `RingChiKernelResult` fields don't map to derivation terms | **Adopted (partial)** → §3 E1 (rename fields to `full_tensor_over_r3`, `dipolar_kernel`, `scalar_kernel`); internal struct, no contract. |
| C3 | cpp:61-75 — label geometric quantities (`// displacement geometry`, `// scalar kernel`) | **Adopted** → §3 E2. |
| C4 | cpp:77-91 — split tensor builds with `// dipolar kernel` / `// full susceptibility tensor`; named δ | **Adopted** → §3 E3 (labels already present at :74/:77/:83; tighten + add δ name per C9/L2). |
| C5 | cpp:124-128 — filter comment is process prose, shorten to `// exclusion filters` | **Declined** — claude (L9) judges this comment load-bearing (explains why topology check supplements distance at the boundary). Keep; this is the documented reason a reader needs. Conflicts with C5; L9 wins under "comments conform to code / keep the reason." |
| C6 | cpp:141-219 — loop mixes six concerns; add block labels | **Adopted (subset)** → §3 E4 (`// nearby rings`, `// filter pair`, `// ring-neighbour record`, `// store ring kernel`, `// accumulate atom tensor`). |
| C7 | cpp:176-205 — signpost that neighbourhood block is feature-attach, not kernel math | **Adopted** → §3 E5. |
| C8 | cpp:193-201 — `theta` underspecified; rename/comment as folded polar | **Adopted** → §3 E6 (comment + reuse `rho`). |
| C9 | cpp:231-249 — `SampleShieldingAt` relationship to `Compute` not stated; label cutoffs | **Adopted** → §3 E7. |
| C10 | cpp:253-256 — `PackST_RS` opaque | **Adopted** → §3 E8 (rename to `PackSphericalTensor`, file-static). |
| C11 | h:3-20 — clean header | **No action** (acknowledged). |
| C12 | cpp:46 — `M_over_r3` preserves notation not meaning | **Adopted** → §3 E1. |
| C13 | cpp:47 — `K` too compressed | **Adopted** → §3 E1. |
| C14 | cpp:48 — `f` too compressed | **Adopted** → §3 E1. |
| C15 | cpp:61 — `d` → `ring_to_atom` | **Adopted (kernel scope)** → §3 E9. Note: catalogue derivation uses `d`; keep the banner's `d`/`d̂` math notation, rename only the C++ local. |
| C16 | cpp:69 — `d_hat` → `ring_to_atom_hat` | **Declined** — `d_hat` is the catalogue symbol `d̂`; keeping it ties the code to the derivation banner directly above it. Renaming the unit vector but not the math notation would *increase* the gap. Weighed call. |
| C17 | cpp:141-143 — `ai`,`ca` → `atom_index`,`atom` | **Adopted (partial)** → §3 E10 (`ca`→`atom`; keep `ai` as the conventional index used across the sibling ring calculators). |
| C18 | cpp:150-152 — `ri`,`geom` → `ring_index`,`ring_geom` | **Declined** — `ri`/`geom` are the established loop idiom shared verbatim across all five ring calculators (`BiotSavart`, `HaighMallion`, `PiQuadrupole`, `Dispersion`); renaming one file's loop variables breaks cross-file pattern-matching for the reader. Weighed call; cost > benefit. |
| C19 | cpp:179 — `rn` cryptic across 27 lines | **Adopted** → §3 E11 (`rn`→`neighbour`). Note: `rn` is shared idiom across the five ring calculators; rename is file-local only, flagged in §3. |
| C20 | cpp:191 — `direction_to_center = kernel.direction` looks sign-reversed | **Adopted as Usage note + proposed rename** → §4 (item A) and §3 E12. The vector points *center→atom* (away from center); the field name says the opposite. Coherent across all consumers; field-name fix is cross-file, carry-through noted, left to human. |
| C21 | cpp:198 — `theta` is folded inclination, not azimuth | **Adopted** → §3 E6 (duplicate of C8). |
| C22 | cpp:253 — `PackST_RS` shorthand | **Duplicate of C10.** |
| C23 | cpp:74-91 — keep order, improve labels/field names | **Adopted** → §3 E1+E3. |
| C24 | cpp:87-91 — name local terms (`asymmetric_term`, `normal_term`, `dipolar_term`) | **Declined** — the three-term build is already labelled by the banner at :39-42 and the inline comment at :83-84; introducing three named temporaries inside the 3×3 double loop fragments a five-line expression that reads cleanly against the formula. Over-explaining per the clarity bar. |
| C25 | cpp:54 — `ComputeRingChiKernel` "Chi" insider | **Declined** — "Chi"/`chi` is the project-wide token for this calculator (`chi_tensor`, `chi_scalar`, `ringchi_shielding`, `CalculatorId::RingSusceptibility`); the file-static helper name matches the family. Renaming it alone desyncs from the field names it writes. |
| C26 | cpp:253 — `PackST_RS` → `PackSphericalTensor9` | **Duplicate of C10** (name chosen in E8). |
| C27 | h:43 — `SampleShieldingAt` → `SampleShieldingTensorAt` | **Declined** — duplicate concern with C1; the method is named identically across all six `*Result` siblings (`HaighMallion`, `McConnell`, `BiotSavart`, …). Cross-file public-API rename for marginal gain; decline. |
| C28 | h:39-41, h:46 — `Compute`/`WriteFeatures` clean | **No action.** |
| C29 | cpp:25-43 — derivation banner heavy, partly duplicates header | **Declined** — the cpp banner names each term's irrep (T0/T1/T2) and the auxiliary `K`/`f` outputs; the header banner is the physical narrative. They serve different readers (API vs implementer) and the duplication is shallow. claude (L13) calls the cpp banner "a model signpost." Keep both. |
| C30 | cpp:74 — `// scalar kernel` | **Adopted** → §3 E3. |
| C31 | cpp:77 — `// dipolar kernel` | **Adopted** → §3 E3. |
| C32 | cpp:83-84 — `// full susceptibility tensor` | **Adopted** → §3 E3. |
| C33 | cpp:124-128 — `// exclusion filters` | **Declined** — duplicate of C5; L9 keeps the rationale. |
| C34 | cpp:145 — `// nearby rings` | **Adopted** → §3 E4. |
| C35 | cpp:157 — `// filter pair` | **Adopted** → §3 E4. |
| C36 | cpp:164 — `// record exclusion` | **Adopted** → §3 E4. |
| C37 | cpp:176-178 — shorten to `// ring-neighbour record` | **Partially declined** — keep the BiotSavart cross-calculator note (L10 calls it useful) but add the C7 signpost in front. → §3 E5. |
| C38 | cpp:193 — `// folded polar coordinates` | **Adopted** → §3 E6. |
| C39 | cpp:207 — `// store ring kernel` | **Adopted** → §3 E4. |
| C40 | cpp:212 — `// accumulate atom tensor` | **Adopted** → §3 E4. |
| C41 | cpp:217 — `// decompose atom total` | **Adopted** → §3 E4. |
| C42 | cpp:191 (correctness) — verify `RingNeighbourhood` direction contract before trusting sign | **Adopted as Usage note** → §4 item A. Traced: coherent, name misleading. |
| C43 | cpp:198 (correctness) — verify `theta` downstream expectation | **Adopted as Usage note** → §4 item B. Traced: coherent. |

### claude.md

| # | Finding (loc) | Disposition |
|---|---|---|
| L1 | cpp:198 — `theta` recomputes `d_plane.norm()` (already `rho`); `abs(z)` folds hemispheres unexplained; reuse `rho` + comment | **Adopted** → §3 E6. (Reuse `rho` is a free clarity win; does not move the number — `atan2(rho, |z|)` == `atan2(d_plane.norm(), |z|)`.) |
| L2 | cpp:90/:81 — inline Kronecker δ ternary; add `// δ_ab` | **Adopted** → §3 E3 (trailing `// δ_ab` on both, or named `delta_ab`; chose trailing comment to avoid fragmenting the expression). |
| L3 | cpp:46 `M_over_r3` good; :48 `f` acceptable; :47 `K` fine given comments | **Noted** — claude finds these acceptable; codex (C2/C12-C14) wants renames. Adopted the renames (E1) because the field names also surface in the consumer-facing struct and the banner refers to `M`/`K`/`f` only in the math, not the code identifiers. Weighed toward codex. |
| L4 | cpp:179 `rn` mildly taxing → `neighbourhood`/`ring_nb` | **Adopted** → §3 E11 (chose `neighbour`). Duplicate of C19. |
| L5 | cpp:61,194 `d` reused as two displacements; in `Compute` name `to_atom` | **Adopted** → §3 E9 (kernel `d`→`ring_to_atom`) + E13 (Compute `d`→`atom_offset`). |
| L6 | grouping clean (f→K→M_over_r3) | **No action** (acknowledged). |
| L7 | cpp:193-201 — `theta` the one unnamed-intermediate snag | **Adopted** → §3 E6. Duplicate of L1/C8. |
| L8 | cpp:253 `PackST_RS` → `PackSphericalTensor` | **Adopted** → §3 E8. Duplicate of C10. |
| L9 | cpp:124-128 filter comment good/load-bearing; keep | **Adopted** → keep; overrides C5/C33. |
| L10 | cpp:176-178 cross-calculator (BiotSavart) signpost useful; keep | **Adopted** → §3 E5 keeps it. |
| L11 | no verbose/stale comments | **No action.** |
| L12 | cpp:241-243 — `SampleShieldingAt` reimplements accept logic as 3 inline checks vs `KernelFilterSet`; `< radius` reject differs from `DipolarNearFieldFilter` boundary | **Adopted as Question** → §6 Q1. Not a readability edit; a real convention question worth surfacing. |
| L13 | cpp:208-210 — per-ring `chi_*` overwritten (not accumulated) while `M_total` accumulates; confirm intent | **Adopted as Usage note** → §4 item C. Traced: intentional (per-neighbour record vs per-atom sum). Add one-line note → §3 E14. |
| L14 | cpp:72 — `cos_theta` assumes unit `ring_normal`; check normalization | **Adopted as Usage note** → §4 item D. Traced: `RingGeometry::normal` is SVD `matrixV().col(2)` (`Ring.cpp:32`), unit by construction; sign-fixed at `:38-39`. No guard needed. Coherent. |
| L15 | cpp:64 singularity guard precedes division; `r==0` early-return zeroes cleanly | **No action** (confirms correctness). |

---

## 3. Edits that don't move numbers

`file:line — change`. Lines are pre-edit; numeric output is unchanged by all.

- **E1** `RingSusceptibilityResult.cpp:45-48` — rename struct fields:
  `M_over_r3 → full_tensor_over_r3`, `K → dipolar_kernel`,
  `f → scalar_kernel`. Internal struct (`RingChiKernelResult`); used only at
  :154/:208-210/:213/:246. Carry the renames to those use sites in the same
  file. No contract crossing — `chi_tensor`/`chi_scalar`/`ringchi_shielding`
  output names are untouched.
- **E2** `cpp:61-75` — add terse block labels: `// ring→atom displacement`
  above :61; the existing `// Ring susceptibility scalar…` at :74 stays.
- **E3** `cpp:77-91` — keep the existing `// Symmetric traceless dipolar
  kernel K_ab` (:77) and `// Full tensor M_ab / r³` (:83) labels; add a
  trailing `// δ_ab` on the two `(a == b ? 1.0 : 0.0)` ternaries (:81, :90).
  Do **not** introduce named term temporaries (declined C24).
- **E4** `cpp:145,157,164,207,212,217` — tighten the in-loop signposts to
  2–4 words: `// nearby rings` (:145), `// filter pair` (:157),
  `// record exclusion` (:164), `// store ring kernel` (:207),
  `// accumulate atom tensor` (:212), `// decompose atom total` (:217).
  These mostly match existing comments; verify each against current text and
  only shorten where wordy.
- **E5** `cpp:176-178` — prepend one signpost line: `// Attach per-ring
  feature record to this atom (not kernel math).` then keep the existing
  BiotSavart cross-calculator note (L10).
- **E6** `cpp:193,198` — replace `// Cylindrical coordinates in ring frame`
  with `// Ring-frame coordinates: z along normal, rho in-plane, theta = polar
  angle from normal (hemisphere-folded)`. On :198 reuse `rho` instead of
  recomputing: `double theta = std::atan2(rho, std::abs(z));`. Arithmetically
  identical (`rho == d_plane.norm()`).
- **E7** `cpp:231-249` — add `// Grid-point sampling: same kernel as Compute,
  evaluated at an arbitrary point (no per-atom filter set — see plan Q1)` above
  the loop; label the three guards `// singularity / inside-ring / cutoff`.
- **E8** `cpp:253,264` — rename file-static `PackST_RS → PackSphericalTensor`.
  File-static, single caller at :264. No collision (it's `static`, file scope).
- **E9** `cpp:61` — kernel local `d → ring_to_atom`; update `d.norm()` (:62),
  `d / r` (:69). Keeps `d_hat` (declined C16) — i.e. `d_hat = ring_to_atom / r`.
- **E10** `cpp:142` — `ca → atom` (`auto& atom = conf.MutableAtomAt(ai);`),
  carry to :168, :180, :203-204, :208-210, :218. Keep `ai` (idiom).
- **E11** `cpp:179` — `rn → neighbour`; carry to :180-191, :203-204, :208-210.
  File-local; note the `rn` idiom is shared across sibling ring calculators
  (cosmetic divergence, no functional coupling).
- **E12** `cpp:191` — add trailing comment:
  `// contract: unit vector pointing center→atom (see RingChiKernelResult::direction)`.
  Documents the existing convention without renaming the shared field. (Rename
  proposal is §3 E15 below, left for human.)
- **E13** `cpp:194` — `Compute`'s local `d → atom_offset`
  (`Vec3 atom_offset = atom_pos - geom.center;`); carry to :195-196.
  Distinguishes it from the kernel's `ring_to_atom` (L5).
- **E14** `cpp:207-210` — add one line above `rn->chi_tensor = …`:
  `// per-neighbour record (overwritten per ring); per-atom sum is M_total below`.
  Closes L13.
- **E15 (PROPOSED, not adopted — cross-file)** `ConformationAtom.h:35`
  `direction_to_center` is misnamed: every producer (`BiotSavart.cpp:252`,
  `RingSusceptibility.cpp:191`, `PiQuadrupole.cpp:194`, `HaighMallion.cpp:311`,
  `Dispersion.cpp:317`) stores `(atom − center)` normalized, i.e. center→atom.
  A faithful name is `unit_center_to_atom` (or `direction_from_center`).
  **Carry-through cost:** 5 producer files + 1 consumer (`ui/RestServer.cpp:1185`
  JSON key `"direction_to_center"` — a UI contract string) + the struct decl.
  Renaming the JSON key would be a UI-facing change. **Recommendation:** keep
  the field name, adopt E12's contract comment instead; surface this to the
  human as a weighed decline. Listed here so the rename is visible, not hidden.

---

## 4. Usage notes (the real product)

### A. `direction` / `direction_to_center` sign (C20, C42)
- **Producer:** `ComputeRingChiKernel` sets `result.direction = d_hat` where
  `d = atom_pos − ring_center` (`cpp:61,70`) — points **center → atom**. The
  struct comment at `:50` says exactly this ("unit vector from ring center to
  atom"). `Compute:191` copies it into `RingNeighbourhood.direction_to_center`.
- **Reason it's not used in the kernel sign:** the tensor math uses
  `cos_theta = d_hat·n̂` (`:72`) and `d̂⊗d̂` / `d̂⊗n̂` products. These are
  even or sign-tracking in `d̂` per the catalogue derivation
  (`GEOMETRIC_KERNEL_CATALOGUE.md:316-318`); the convention is fixed by that
  derivation, not by the field name.
- **Consumers:** the only reader of `direction_to_center` is
  `ui/RestServer.cpp:1185`, which echoes it verbatim as JSON for display. No
  consumer does sign-sensitive arithmetic on it. **Producer and consumers
  agree; the field NAME is the only thing that misstates the direction.**
- **Verdict: coherent (expected).** Fix = E12 contract comment; optional
  cross-file rename E15 left to the human.

### B. `theta` = polar angle from ring normal, hemisphere-folded (C8/C21/C43, L1/L7)
- **Producer:** `theta = atan2(rho, |z|)` with `z = d·n̂` (axial) and `rho`
  the in-plane component norm (`cpp:194-198`). So `tan(theta) = rho/|z|` ⇒
  `theta` is the angle **from the ring normal/axis** (0 on-axis, π/2 in-plane),
  folded to one hemisphere by `|z|`. It is **not** the cylindrical azimuth.
- **Internal consistency check:** the per-pair McConnell scalar emitted at
  `ConformationResult.cpp:84` uses `cos_th = z / distance = cos(theta)`, and
  computes `(3cos²θ − 1)/r³`. This matches `f` in `ComputeRingChiKernel:75`
  (`cos_theta = d̂·n̂ = z/r`). So col[6]=`theta` and col[7]=`mcconnell_factor`
  use the *same* θ. Coherent.
- **Consumers:**
  - `ConformationResult.cpp:78` packs raw `theta` into `ring_contributions.npy`
    col [6]; SDK `python/nmr_extract/_ring.py:10` documents col[6] as
    `theta (rad)` and col[7] as `(3cos²θ−1)/r³`. **Producer/SDK agree.**
  - The azimuth (in-plane angle) is a *separate* quantity stored as
    `cos_phi`/`sin_phi` (cols 56/57, `BiotSavart.cpp:266-273`) — so there is no
    collision: `theta`=polar-from-normal, `phi`=azimuth.
  - **Trajectory path is different and must not be conflated:**
    `RingNeighbourhoodTrajectoryStats.cpp:159-162` channel [3] is
    `in_plane_angle` (the *azimuth* `atan2(sin_phi,cos_phi)` in [0,2π)),
    **not** the per-conformation `theta`. h5-reader `QtSpecialBuffers.h:34`
    and SDK `_trajectory.py:1932` correctly label channel [3]
    `in_plane_angle`. So the TR pipeline emits azimuth; the per-conformation
    pipeline emits polar `theta`. Both are correct for their layouts — a
    reader skimming both files could mistake one for the other, which is
    exactly why E6's comment should say "polar angle from normal" explicitly.
- **Verdict: coherent (expected).** Fix = E6 comment + `rho` reuse.

### C. Per-ring `chi_*` overwrite vs `M_total` accumulate (L13)
- **Producer:** `cpp:208-210` writes `rn->chi_tensor / chi_spherical /
  chi_scalar` for the *current* ring (overwrites any prior value on a
  re-found neighbour); `cpp:213` accumulates `M_total += kernel.M_over_r3`
  across rings, decomposed once at `:218` into
  `ringchi_shielding_contribution`.
- **Reason:** `RingNeighbourhood` is a *per-(atom,ring)* record — one entry per
  ring — so storing this ring's kernel on this ring's record is correct; there
  is nothing to sum within a single record. The per-atom sum lives separately
  on `ca.ringchi_shielding_contribution`. The find-or-create at :180-205 keys
  by `ring_index`, so re-find only happens when another calculator (BiotSavart)
  created the record first; the chi fields are still this ring's.
- **Consumers:** `ring_contributions.npy` col [45:54] = per-pair
  `chi_spherical` (`ConformationResult.cpp:92`); `ringchi_shielding.npy` =
  per-atom sum (`WriteFeatures:265`). Two distinct outputs, both intended.
- **Verdict: coherent (expected).** Fix = E14 one-line note.

### D. `ring_normal` unit-length assumption (L14)
- **Producer of normal:** `Ring::ComputeGeometry` sets
  `geo.normal = svd.matrixV().col(2)` (`Ring.cpp:32`) — a column of an
  orthonormal SVD `V`, hence unit-length — then sign-orients it (`:38-39`)
  without rescaling. `cos_theta = d̂·n̂` at `cpp:72` is therefore a true cosine.
- **Verdict: coherent (expected).** No guard needed; no edit.

---

## 5. Bug-by-exhaustion candidates

None. Every sign/value the reviews flagged (`direction_to_center`, `theta`,
the chi overwrite, the normal assumption) traced to a coherent reason with
producer and consumers in agreement (§4 A–D). The only outright defects found
would be naming/comment imprecision, which move no numbers.

---

## 6. Questions & Ambiguities

- **Q1 (from L12) — `SampleShieldingAt` filter divergence.** `Compute` uses
  `KernelFilterSet{ MinDistanceFilter, DipolarNearFieldFilter,
  RingBondedExclusionFilter }` (`cpp:129-132`). `SampleShieldingAt` instead
  applies three inline guards: singularity, `distance < geom.radius`, and
  `distance > ring_current_spatial_cutoff` (`cpp:241-243`). The
  `RingBondedExclusionFilter` is meaningfully absent (topological — no
  bonded-atom notion at a free grid point, so reasonable to drop), but the
  `< geom.radius` reject is a *different* near-field boundary than
  `DipolarNearFieldFilter` (which uses `source_extent = 2·radius` as its scale,
  `cpp:160`). The grid sampler is a visualization/diagnostic path
  (`ui/MainWindow.cpp:1844` registers the calculator label; the overlay seeds
  grid points), not a feature-emitting path, so a looser/different boundary may
  be intended. **Question for the human:** is the `< radius` cutoff in
  `SampleShieldingAt` the intended visualization near-field boundary, or should
  it mirror `DipolarNearFieldFilter`'s `2·radius` scale? Not touched in this
  pass either way (would move grid values).

- **Q2 — `direction_to_center` rename (E15).** Adopt the faithful name
  (`unit_center_to_atom`) across 5 producers + the `ui/RestServer.cpp` JSON key,
  or keep the field and settle for E12's contract comment? Recommendation is
  the comment (declines the cross-file/UI-contract churn), but the call is the
  human's because the JSON key is a small UI-facing contract.
