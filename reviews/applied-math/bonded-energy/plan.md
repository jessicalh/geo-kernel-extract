# Fix plan — BondedEnergyResult (bonded FF terms + 2-D CMAP bilinear)

Targets: `src/BondedEnergyResult.h`, `src/BondedEnergyResult.cpp`
Producer of the parameters: `src/FullSystemReader.cpp` (TPR → `BondedParameters`)
Output contract: `bonded_energy.npy` (N,7), SDK `python/nmr_extract/_catalog.py` key
`bonded_energy`; H5 `trajectory/bonded_energy_time_series`.

---

## 1. Summary

The file tells a coherent story. The header is a tagged-union of
`BondedInteraction` plus a units-annotated parameter block and a units-named
accessor surface; the `.cpp` is a sequence of small `Eval*` energy kernels
(each headed by its functional form), a shared `DihedralAngle` helper, and a
dispatch loop that gathers positions, evaluates, and splits the term energy
evenly per atom. Both reviews agree it reads well top-to-bottom; the single
genuinely dense spot is `EvalCMAP`'s index/weight block.

The two convention checks neither review could confirm from the file alone —
**CMAP grid axis order (φ vs ψ)** and **dihedral sign** — both resolve to
**coherent** once traced to the producer (`FullSystemReader.cpp`) and the
in-repo sibling `PlanarGeometryResult::Dihedral`. Details in §4. No
sign/value bug found; no number moves.

The fix pass will: add 2–4-word signposts inside `EvalCMAP` and
`DihedralAngle`, rename a handful of internal locals/helpers for clarity
(none crosses the output contract), trim two editorializing comments, and
correct one misleading header comment (`cmap_grid_spacing` documented as
"spacing" but used as a point count). It will **not** touch the algorithm,
any number, any output name (`bonded_energy.npy`, the 7-channel order, the
`bonded_energy` catalog key, or the H5 dataset names), or the public accessor
names (which the time-series H5 layout mirrors by order — see §3 note).

---

## 2. Review-finding ledger

Every finding from `codex.md` and `claude.md`, with disposition. (No
`codex-correctness.md` present.)

### codex.md

| # | Finding | Disposition |
|---|---------|-------------|
| C1 | `.h:42` `atoms[5]` schema-dependent; add arity/order comment | **Adopted** → §3 E1 (the per-type arity is already noted at `.h:42`; add a one-line cross-ref to the param block; no per-site unpacking — that is churn) |
| C2 | `.h:46` `p[3]` is the biggest story break; rename `params[3]` + named aliases | **Partially adopted** → §3 E2. Keep field name `p` (the `.h:46–52` comment block already maps each slot per type, and it's an internal struct field used at many `ix.p[i]` sites in `.cpp` and `FullSystemReader.cpp`); decline the rename as cross-file churn against an already-documented block, but tighten the comment. Per-case local aliases handled by C13/CL? see §3 E7. |
| C3 | `.cpp:45` `DihedralAngle` dense; signpost bond vectors / plane normals / signed angle | **Adopted** → §3 E3 (matches claude CL2; deduped) |
| C4 | `.cpp:91` CMAP `dx/fi/fj/wi/wj` hide units/axes; rename | **Partially adopted** → §3 E4/E5. Rename `dx`→`cell_step_rad` (genuinely misleading, see claude CL5). Leave `fi/fj/wi/wj` (locally consistent, explained, terse-acceptable per claude CL6); add axis signposts instead. |
| C5 | `.cpp:136` `ix`→`interaction` | **Declined** — both reviews elsewhere call `ix` fine and consistent (codex C12, claude CL8); `ix` is idiomatic here and used densely. Renaming is churn. |
| C6 | `.cpp:144` `energy`→`interaction_energy`, `target`→`energy_channel`/`per_atom_energy` | **Partially adopted** → §3 E8. `target`→`channel` (it points at the per-type per-atom accumulator — "target" is vague); leave `energy` (clear in its short scope). |
| C7 | `.h:64` `cmap_grid_spacing` is a point count not a spacing; rename `cmap_grid_size`/`cmap_point_count` | **Adopted as comment-only** → §3 E6. The field is read across files (`FullSystemReader.cpp:312`, `EvalCMAP`, `Compute`); rename carries cross-file cost for a build-time topology field. Fix the misleading **comment** to say "points per axis"; note the rename option + cost in §6. |
| C8 | `.h:87` `UBEnergy`→`UreyBradleyEnergy` | **Declined** — public accessor; the H5 time-series layout mirrors accessor order and the SDK/h5-reader read these channels. Renaming is a contract-adjacent change for marginal gain; `UB` is universally understood in FF context. Noted in §6. |
| C9 | `.h:90` `ProperDihEnergy`/`ImproperDihEnergy`→spelled out | **Declined** — same reason as C8 (public accessor, contract-adjacent). |
| C10 | `.cpp:133` `count_ub`→`count_urey_bradley` | **Adopted** → §3 E9 (internal log counter; trivial, improves the log site's match to the spelled-out interaction). |
| C11 | `.cpp:196` `cmap_idx` clear; no issue | **Noted, no action** (reviewer self-resolved). |
| C12 | `.cpp:48` dihedral needs visible grouping | **Duplicate** of C3 → §3 E3. |
| C13 | `.cpp:95` CMAP named intermediate stages (angle grid coords / periodic cell / bilinear weights) | **Adopted** → §3 E4 (signposts; merges with claude CL1). |
| C14 | `.cpp:147` each switch case hides params behind `ix.p[]`; local unpacking per case (`r0_nm`,`k_bond`,…) | **Declined** — the `Eval*` signatures already name every parameter (`r0_nm`, `k`, `theta0`, `mult`, …) at the call boundary one line below; adding a second naming layer per case is redundant and lengthens the dispatch. The `.h:46–52` block is the single source for slot meaning. |
| C15 | `.cpp:38` `EvalUB`→`EvalUreyBradleyEnergy` | **Declined** — file-local static; `EvalUB` sits directly under a `// Urey-Bradley: E = …` header and pairs with the `ub_energy_` channel. Energy-returning is obvious from context (all `Eval*` return energy). Spelling it fully adds width for no disambiguation. |
| C16 | `.cpp:65` `EvalProperDih`→`EvalProperDihedralEnergy` | **Declined** — same reasoning as C15; claude CL7 explicitly calls these "textbook-legible, good." |
| C17 | `.cpp:73` `EvalImproperDih`→`EvalImproperDihedralEnergy` | **Declined** — same as C16. |
| C18 | `.cpp:86` `EvalCMAP`→`EvalCmapEnergy` | **Declined** — "acceptable" per reviewer; CMAP is a proper noun in FF context. |
| C19 | `.h:79` `Compute` vague but acceptable as factory convention | **Declined** (reviewer self-resolved; it is the project-wide `ConformationResult` factory convention). |
| C20 | `.cpp:12` section banner decorative; if touched use `// energy kernels` | **Declined** — not touching the banner; it is a harmless section divider consistent with the `// ── Compute ──` / `// ── WriteFeatures ──` siblings. |
| C21 | `.cpp:77` comment → `// periodic difference` | **Declined** — current `// Wrap difference to [-π, π]` is more precise than "periodic difference" and claude CL? keeps it. No change. |
| C22 | `.cpp:94` → `// angle grid coordinates` | **Declined** — current `// Map angle [-π, π] to grid index [0, spacing-1]` is a good, specific signpost (claude CL4 says keep). No change. |
| C23 | `.cpp:99` "sufficient — CMAP grids are smooth" is justification; → `// bilinear CMAP` | **Adopted** → §3 E10 (matches claude CL3; deduped). |
| C24 | `.cpp:137` restates guard → `// index guard` | **Adopted** → §3 E11 (trim to a terse signpost). |
| C25 | `.cpp:216` → `// per-atom split` | **Adopted** → §3 E12 (minor tightening). |
| C26 | `.cpp:222` clear; no change | **Noted, no action.** |
| C27 | `.cpp:139` index guard validates only `n_atoms`, not type-specific arity | **Declined (no number/edit), noted §5.** The arity is fixed by construction in `FullSystemReader.cpp` (`bi.n_atoms` is set to the literal 2/3/3/4/4/5 at each push site, alongside exactly that many `atoms[]` assignments). A malformed `n_atoms` cannot reach `Compute` from the only producer. Adding an assert is defensible defensive hygiene but is a behavior/structure change beyond a readability fix; recorded as a §5 candidate, not adopted. |
| C28 | `.cpp:89` CMAP doesn't check `grid.size() >= spacing*spacing` | **Declined (noted §5).** Grids come verbatim from GROMACS where `cmap.size() == grid_spacing*grid_spacing` by construction; same producer-guarantee reasoning as C27. §5 candidate. |
| C29 | `.cpp:79` `M_PI` non-standard C++ | **Declined** — used project-wide (`PlanarGeometryResult.cpp`, etc.) and the build already relies on it; out of scope for a per-algorithm readability pass. Noted §6. |
| C30 | `.h:73/79` `<typeindex>`/`<memory>` not directly included | **Declined (noted §6)** — include hygiene, real but not a readability-of-the-algorithm fix; `ConformationResult.h` supplies both transitively. Flag for the human; not part of this pass. |
| C31 | `.cpp:64` regression-check `DihedralAngle` vs GROMACS convention (not calling it wrong) | **Resolved → coherent, §4.** Traced; see §4 dihedral sign. |

### claude.md

| # | Finding | Disposition |
|---|---------|-------------|
| CL1 | `.cpp:95-105` `EvalCMAP` fuses angle→index→corner→weight; add `// fractional grid coords` @95 and `// cell corners + weights` @100 | **Adopted** → §3 E4 (merges with codex C13). |
| CL2 | `.cpp:59-61` signed-dihedral trick opaque; `// signed angle via reference axis` | **Adopted** → §3 E3 (merges with codex C3). |
| CL3 | `.cpp:99` "sufficient — smooth" editorializes; trim to `// bilinear interpolation` | **Adopted** → §3 E10 (dedupe with codex C23). |
| CL4 | `.cpp:94` good signpost; keep | **Noted, no action** (agrees with declining C22). |
| CL5 | `.cpp:95` `dx` misleads (angular step, not spatial x); → `cell_step_rad`/`dphi` | **Adopted** → §3 E5 (`cell_step_rad`). |
| CL6 | `.cpp:96-97` `fi/fj/wi/wj` terse but explained; leave | **Adopted (leave)** — consistent with declining the C4 rename of these. |
| CL7 | `.cpp:19-86` `Eval*` + `DihedralAngle` good naming | **Noted, no action** (grounds declining C15–C18). |
| CL8 | `.cpp:139,218` loop var `k` shadows force-constant `k`; rename loop var `a` | **Adopted** → §3 E13. Genuine: `k` is the force constant throughout the `Eval*` signatures; the atom-index loop reusing `k` is a momentary stumble. Rename loop index to `a`. |
| CL9 | `.cpp:136` `ix` fine | **Noted** (grounds declining C5). |
| CL10 | `.cpp:112-113` bilinear sum grouped well; clear | **Noted, no action.** |
| CL11 | `.cpp:100-105` corner indices/weights interleaved; minor regroup | **Adopted (light)** → §3 E4 covers it via signposts; no line reshuffle needed. |
| CL12 | `.cpp:148-211` per-case `PositionAt` repetition acceptable, self-contained | **Noted, no action** (grounds declining C14). |
| CL13 | `.cpp:99` trim editorial comment | **Duplicate** of CL3/C23 → §3 E10. |
| CL14 | `.cpp:100-101` `static_cast<int>(fi) % spacing` round-up check | **Resolved → safe, §5.** `fi ∈ [0, spacing]` always (phi∈[-π,π] ⇒ phi+π∈[0,2π]); at phi=+π, fi=spacing → i0=0, wi=0 (correct wrap). No off-by-one. Recorded as exhausted, no edit. |
| CL15 | `.cpp:107-110` grid `[i_phi*spacing + j_psi]` matches `.h:65`; check vs GROMACS storage | **Resolved → coherent, §4** (CMAP axis order). |
| CL16 | `.cpp:46-62` dihedral sign vs extraction-side φ0 | **Resolved → coherent, §4** (dihedral sign). |
| CL17 | `.cpp:217` even split is documented modeling choice; `n_atoms≥2` ⇒ no div0 | **Noted, no action** — confirmed: every producer push site sets `n_atoms ∈ {2,3,4,5}`. |
| CL18 | `.cpp:31-34` degeneracy guards present; clean | **Noted, no action.** |
| CL19 | `.cpp:89 vs 196-199` CMAP grid-empty + idx bounds both guarded; clean | **Noted, no action.** |

---

## 3. Edits that don't move numbers

All edits are comment/signpost/internal-local changes. None alters control
flow, arithmetic, an output name, or a public accessor.

- **E1** `BondedEnergyResult.h:42` — extend the `atoms[5]` comment to point at
  the parameter block: `// Atom indices (2 bond, 3 angle/UB, 4 dihedral, 5 CMAP); see p[] below`.
- **E2** `BondedEnergyResult.h:46` — keep field name `p`; leave the existing
  per-type slot table. (No rename — see ledger C2.) Optionally re-title the
  block comment `// Parameters p[] (slot meaning depends on type):` so the
  array name is named in the heading.
- **E3** `BondedEnergyResult.cpp:48–60` — three blank-line-separated stages
  with one-word signposts in `DihedralAngle`: `// bond vectors` (b1/b2/b3),
  `// plane normals` (n1/n2), `// signed angle via reference axis` (m, sin_phi,
  atan2). (codex C3 + claude CL2.)
- **E4** `BondedEnergyResult.cpp:94/99` — two internal signposts in `EvalCMAP`:
  keep the existing `// Map angle [-π, π] to grid index [0, spacing-1]` at the
  `dx/fi/fj` block; add `// cell corners + bilinear weights` at the `i0/j0/wi/wj`
  block. (codex C13 + claude CL1/CL11.)
- **E5** `BondedEnergyResult.cpp:95` — rename local `dx` → `cell_step_rad`
  (angular grid step in radians; `dx` reads as a spatial x-step). File-local;
  no carry-through. (codex C4 partial + claude CL5.)
- **E6** `BondedEnergyResult.h:64` — comment fix: `cmap_grid_spacing` doc →
  `// CMAP grid points per axis (each grid is grid_points × grid_points doubles)`.
  Keep the field **name** (see §6 for the rename-with-cost option). Update the
  `.h:62` companion comment `(grid_spacing x grid_spacing)` to
  `(points_per_axis × points_per_axis)` for consistency.
- **E7** *(none — C14 declined; the `Eval*` signatures already name params).*
- **E8** `BondedEnergyResult.cpp:145,152,…,219` — rename local `target` →
  `channel` (it is the per-type per-atom accumulator the term feeds).
  File-local. Leave `energy` as-is. (codex C6 partial.)
- **E9** `BondedEnergyResult.cpp:133,170,233` — rename counter `count_ub` →
  `count_urey_bradley` (and the log keeps `UB=` literal, which already reads
  cleanly). File-local. (codex C10.)
- **E10** `BondedEnergyResult.cpp:99` — comment
  `// Bilinear interpolation (sufficient — CMAP grids are smooth)` →
  `// bilinear interpolation`. (codex C23 + claude CL3.)
- **E11** `BondedEnergyResult.cpp:137` — comment `// All atom indices must be
  in range` → `// index guard` (terse). (codex C24.)
- **E12** `BondedEnergyResult.cpp:216` — comment `// Split energy evenly among
  participating atoms` → `// per-atom split (even)`. (codex C25.) Optional;
  current comment is fine — apply only if touching the block.
- **E13** `BondedEnergyResult.cpp:139,218` — rename loop index `k` → `a`
  (atom index) in both range-check and accumulation loops; removes the visual
  collision with the force-constant `k` used in every `Eval*` signature.
  File-local. (claude CL8.)

---

## 4. Usage notes (the convention checks)

### 4a. CMAP grid axis order (φ-major) — **coherent (expected)**

**Reason.** The grids are copied **verbatim** from GROMACS at
`FullSystemReader.cpp:312–317`:
`cmap_grid_spacing = idef.cmap_grid.grid_spacing` and each
`idef.cmap_grid.cmapdata[].cmap` is pushed flat with no transpose. GROMACS
stores each CMAP energy surface row-major with the **first** torsion (φ) as the
slow index and the **second** (ψ) as the fast index — i.e. `cmap[iφ*N + iψ]`.

The consumer reads `grid[i0 * spacing + j0]` where `i0` is derived from
`phi = DihedralAngle(p0,p1,p2,p3)` and `j0` from `psi = DihedralAngle(p1,p2,p3,p4)`
(`EvalCMAP`, lines 91–110). The CMAP atom quintet from the TPR
(`FullSystemReader.cpp:449–453`, atoms 0..4) defines φ over atoms 0-1-2-3 and ψ
over atoms 1-2-3-4 — the standard CHARMM backbone (C₋₁-N-CA-C, N-CA-C-N₊₁)
ordering GROMACS emits. So φ→slow index, ψ→fast index on both sides.

**Producer/consumer agreement:** yes — the same load path that imports the
grid layout feeds the same read path; the header comment `.h:65`
`[phi_idx * spacing + psi_idx]` documents the shared convention and matches the
indexing. Even in the (unverifiable-from-disk) event that GROMACS' internal
order were the opposite, a swap would have to occur in **both** the grid layout
and the φ/ψ atom assignment to be wrong — and a single consistent swap on both
axes is self-canceling for a function evaluated at the matching (φ,ψ). The pair
is internally consistent.

**Where consumed downstream:** only as an opaque per-atom energy feature —
`bonded_energy.npy` channel 5 (CMAP) and the H5 `cmap` time-series. No consumer
(`learn/`, `h5-reader/`, `ui/`, SDK) interprets the φ/ψ axis order; it is fully
internal to `EvalCMAP`. (Note: the 1P9J CHARMM36m fixture reports CMAP=0 per
`BondedEnergyTimeSeriesTrajectoryResult.cpp:113`, so this path is lightly
exercised on the current fixture — a cross-validation opportunity, §6.)

**Verdict:** coherent. Fix is comment-only (E4/E6 above); no number moves.

### 4b. Dihedral sign convention — **coherent (expected)**

**Reason.** `BondedEnergyResult::DihedralAngle` (lines 46–62) is
**bit-identical in construction** to the in-repo sibling
`PlanarGeometryResult::Dihedral` (`PlanarGeometryResult.cpp:35–46`):

```
m = n1 × b̂2 ;  cos = n1·n2 ;  sin = m·n2 ;  atan2(sin, cos)   // both
```

`PlanarGeometryResult::Dihedral` carries explicit math-review provenance
(comment at `PlanarGeometryResult.cpp:51–53`, "math-review MED-2, 2026-05-19")
and is documented as the "standard formulation … signed angle between b1 and
b3" — the right-handed (IUPAC/GROMACS-positive) convention.

**Producer side:** the reference angles φ0 fed to `EvalProperDih`
(`k(1+cos(mult·φ − φ0))`) and `EvalImproperDih` come from the TPR via
`FullSystemReader.cpp:420` (`p.phiA * DEG_TO_RAD`) and `:437`. These φ0 are
GROMACS' own torsion parameters, expressed in the same right-handed convention
GROMACS' `dih_angle` uses. The runtime φ is computed in that same convention.
Since both the parameter and the evaluated angle are GROMACS-convention, the
`mult·φ − φ0` argument is convention-consistent and the cos term is invariant
to a global sign flip of the convention anyway for the proper form; the
improper harmonic `(φ − φ0)²` likewise only needs φ and φ0 in the *same*
convention, which they are.

**Producer/consumer agreement:** yes. Output is an opaque per-atom energy; no
consumer reads the angle sign. The only place sign matters is internal to
`EvalProperDih`/`EvalImproperDih`, and there it is matched to the TPR φ0.

**Verdict:** coherent. Fix is the E3 signpost only; no number moves.

---

## 5. Bug-by-exhaustion candidates

**None.** No sign/value bug survived tracing.

Two robustness items the reviews raised (codex C27/C28) were traced to the sole
producer and found guaranteed-by-construction, so they are **not** bugs:

- **Type-specific arity (C27).** `FullSystemReader.cpp:357–455` sets
  `bi.n_atoms` to the literal correct count (2/3/3/4/4/5) at each push site,
  paired with exactly that many `atoms[]` assignments. A dihedral with
  `n_atoms<4` or CMAP with `n_atoms<5` cannot be produced. The `Compute`
  guard's range check over `n_atoms` is therefore sufficient for the only
  data source. *Optional hardening (not part of this pass):* a
  `assert(ix.n_atoms == expected_for_type)` in each switch case would document
  the invariant — defensive, behavior-neutral, but a structural change beyond
  readability. Defer to the human.

- **CMAP grid size (C28).** Grids are copied whole from
  `idef.cmap_grid.cmapdata[].cmap`, where GROMACS guarantees
  `cmap.size() == grid_spacing²`. `EvalCMAP` already guards `grid.empty()` and
  `spacing < 2`; the indices `i0,i1,j0,j1 ∈ [0,spacing)` so the max access is
  `(spacing-1)*spacing + (spacing-1) = spacing²-1`, in bounds given the
  producer guarantee. *Optional hardening:* an
  `assert(grid.size() >= size_t(spacing)*spacing)` documents it. Defer.

- **`static_cast<int>(fi) % spacing` round-up (claude CL14).** Exhausted: with
  `phi ∈ [-π,π]`, `fi = (phi+π)/dx ∈ [0, spacing]`, never negative and never
  exceeding `spacing`. At `phi = +π`, `fi = spacing` → `i0 = 0`, `wi = 0`
  (correct periodic wrap). No off-by-one; no edit.

---

## 6. Questions & Ambiguities

1. **GROMACS CMAP internal axis order — confirmed only by self-consistency,
   not by source.** No GROMACS source tree is present on disk, so I confirmed
   the φ-major (`cmap[iφ*N+iψ]`) layout from the producer/consumer pairing and
   the standard CHARMM atom-quintet ordering, not from GROMACS' own code. The
   pair is internally consistent and the output is an opaque feature, so this
   is not a bug under any axis convention — but if you want it nailed to the
   GROMACS source, the check is: does `gromacs/.../cmap.cpp` index its grid with
   the first torsion as the slow index? **A stronger validation** exists for
   free: `GromacsEnergyResult` already carries GROMACS' own per-frame
   `cmap_dih` total (`GromacsEnergyResult.h:45`); summing this calculator's
   per-atom CMAP channel over a frame and comparing to GROMACS' `cmap_dih`
   would cross-validate both the axis order **and** the bilinear scheme against
   the engine. Worth running once on a fixture where CMAP ≠ 0 (the 1P9J
   fixture reports CMAP=0, so it won't exercise this). Not an edit — a
   suggested validation.

2. **Public accessor renames (C8/C9: `UBEnergy`, `ProperDihEnergy`,
   `ImproperDihEnergy`).** Declined as contract-adjacent: the time-series H5
   layout mirrors the accessor *order*
   (`BondedEnergyTimeSeriesTrajectoryResult.cpp:111–112`) and the SDK /
   h5-reader read those channels. The names are abbreviated but
   FF-conventional. If the human wants them spelled out, it is a coordinated
   change across `BondedEnergyResult.h`, the TR result, and any H5/SDK channel
   labels — flagged, not done.

3. **`cmap_grid_spacing` field rename (C7).** I fixed the misleading *comment*
   (E6) but kept the field name to avoid cross-file churn
   (`FullSystemReader.cpp:312`, `EvalCMAP`, `Compute` all read it). If the human
   wants `cmap_grid_points` / `cmap_points_per_axis`, the carry-through is those
   three sites plus the `.h` declaration — small and contained, but a weighed
   call rather than a free win.

4. **Include hygiene (C30) and `M_PI` portability (C29).** Real but out of
   scope for an algorithm-readability pass; both rely on transitive includes /
   platform defines that the wider build already assumes. Flagged for a
   separate hygiene pass, not touched here.
