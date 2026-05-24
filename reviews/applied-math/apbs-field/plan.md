# Fix plan — ApbsFieldResult (trilinear interp + central-difference E / EFG)

Targets: `src/ApbsFieldResult.h`, `src/ApbsFieldResult.cpp`
Scope: readability only. No algorithm change, no number moves. Output
names (`apbs_E`, `apbs_efg`, `apbs_efield`, `apbs_efg`, `apbs_efg_spherical`)
are protected and untouched.

---

## 1. Summary

The file tells a mostly coherent three-act story: grid interpolation
utilities → APBS solve + per-atom extraction → query/write. The math is
sound and every sign/unit traced below is **coherent with its consumers**.

Two real readability drags, both flagged by both reviewers and confirmed
here:

1. **Stale "fall back to vacuum Coulomb" banner** at `.cpp:284-286`. The
   code does the exact opposite — `Compute` returns `nullptr` with an
   explicit "No fallback" log, and the header (`.h:15-17`) states the
   no-fallback policy. The banner is the single most actively-misleading
   line in the file. Comment-only fix.

2. **Unsignposted sanitise/convert/store loop** at `.cpp:242-273`. One
   loop body fuses four distinct steps (finite-value guard, field-magnitude
   cap, unit conversion, store + decompose) under a bare `// Sanitise`.
   Signpost-only fix.

Plus a cluster of **stale line-number cross-references and review-provenance
prose** in the long comments (`.cpp:93-104`, `107-114`, `327-331`,
`347-350`), and a few **internal-name** improvements (`xArr/yArr/zArr`,
`lo/hi`, `gridResult`). The trajectory sibling
`ApbsEfgTimeSeriesTrajectoryResult.cpp` carries the *same* stale line-number
references (`ApbsFieldResult.cpp:262-272`, "line 102"); noted as carry-through,
not edited under this algorithm's scope unless the human opts in.

The fix pass will **not** touch: the algorithm, any stored/serialized number,
the grid parameters (161³, pdie 4.0, sdie 78.54, 0.15 M, ±20 Å padding), the
symmetrize-then-detrace sequence, or the catalog/SDK contract.

---

## 2. Review-finding ledger

Every finding from `codex.md` and `claude.md`, one row each.
(`codex-correctness.md` is not present.)

### codex.md

| # | Finding (loc) | Disposition |
|---|---|---|
| C1 | `.h:3` file-level story clear — no change | declined (no-op; already good) |
| C2 | `.cpp:21` "from old ApbsSolver" is process history → `// Grid interpolation` | adopted → E1 |
| C3 | `.cpp:68` add signpost `// central difference gradient` | adopted (modified) → E2 |
| C4 | `.cpp:81` group `// finite differences / symmetric projection / traceless projection` | partially adopted → E3 (blocks already labelled; only tighten) |
| C5 | `.cpp:93` shorten symmetrization comment, drop provenance | adopted → E4 |
| C6 | `.cpp:107` shorten traceless comment | declined (keep physics; see §4 — this comment is load-bearing and accurate) |
| C7 | `.cpp:139` "Separate x,y,z for C bridge" clear — no change | declined (no-op) |
| C8 | `.cpp:160` grid sizing clear — no change | declined (no-op) — but see claude CL3 for the orphaned `grid_dim` |
| C9 | `.cpp:242` split into `// finite-value guard` + `// field magnitude cap` | adopted → E5 |
| C10 | `.cpp:284` stale fallback banner → `// Main compute: APBS solve only` | adopted → E6 |
| C11 | `.cpp:327` WriteFeatures section comment is schema archaeology → `// Feature export` | partially adopted → E7 (trim provenance, keep the T2-only fact) |
| C12 | `.cpp:69` `E` terse; `electric_field` reads better | declined (local conventional symbol in a 10-line fn; `E = -grad φ` is standard) |
| C13 | `.cpp:82` `EFG` → `field_gradient`/`efg_tensor` | declined (domain-standard; `EFG` is the project-wide name, matches output) |
| C14 | `.cpp:83/89` `j`/`i` → `coord_axis`/`field_axis` | declined as rename; adopted as a one-line index-convention comment → E3 (claude CL5 duplicate) |
| C15 | `.cpp:130` `non_authoritative_radii` clear — no change | declined (no-op) |
| C16 | `.cpp:140` `xArr/yArr/zArr` → `x_coords/y_coords/z_coords` | adopted → E8 |
| C17 | `.cpp:164` `lo/hi` → `min_corner/max_corner` | adopted (as `bbox_min/bbox_max`) → E9 (merges claude CL4) |
| C18 | `.cpp:179` `grid_dim` → `grid_points_per_axis` | declined (rename); see E10 for the placement fix instead |
| C19 | `.cpp:197` `gridResult` → `apbs_grid` | adopted → E11 |
| C20 | `.cpp:335` `N` → `atom_count` | declined (conventional loop-count in a 25-line fn; low value) |
| C21 | `.cpp:57` trilinear sum — add `// trilinear interpolation` | adopted → E12 |
| C22 | `.cpp:81` fn name announces only first stage; rename or block-label | partially adopted → E3 (block labels; rename declined, see C26) |
| C23 | `.cpp:93-116` foreground the two projections | duplicate of C4/C5 |
| C24 | `.cpp:251-254` magic threshold needs local meaning / doc `APBS_SANITY_LIMIT` | adopted → E5 (signpost) + §6 note (constant is documented at definition, `PhysicalConstants.h:95`) |
| C25 | `.cpp:262-266` unit conversion clear — no change | declined (no-op) but see claude CL9 to trim the duplicate |
| C26 | `.cpp:81` `FieldGradientFromGrid` → `ProjectedFieldGradientFromGrid`/`TracelessEfgFromGrid` | declined (file-local static; rename gives marginal clarity over a block label, and a longer name fights the call site readability — block labels in E3 carry the same information) |
| C27 | `.cpp:126` `ComputeViaApbs` clear — no change | declined (no-op) |
| C28 | `.cpp:333` `WriteFeatures` framework-generic, acceptable | declined (no-op) |
| C29 | `.h:42` "Factory" → `// compute APBS fields` | adopted → E13 |
| C30 | `.cpp:8` C-bridge signpost good — no change | declined (no-op) |
| C31 | `.cpp:74` terse math comment good — no change | declined (no-op) |
| C32 | `.cpp:196` "Call the C bridge" good — no change | declined (no-op) |
| C33 | `.cpp:268` "both Mat3 AND SphericalTensor" → `// store tensors` | adopted → E5 (as `// store Mat3 + T2`) |
| C34 | `.cpp:347-350` stale line-number + schema history → `// T2 feature export`, drop line ref | adopted → E7 |
| C35 | `.h:46` `ElectricFieldAt` clear (units in header) | declined (no-op) |
| C36 | `.h:47` `FieldGradientAt` → `ElectricFieldGradientAt`/`EfgTensorAt` | declined (public API on header; rename crosses RestServer/inspector consumers — carry-through cost not worth it; see §6) |
| C37 | `.h:48` `FieldGradientSphericalAt` → `ElectricFieldGradientSphericalAt` | declined (same as C36) |
| C38 | `.cpp:164` correctness: `xArr[0]` read before empty-conf guard | declined as bug; recorded as §6 question (framework forbids empty conformations — confirm) |
| C39 | `.cpp:291-293` correctness: returns result with `conf_==nullptr` when ChargeAssignmentResult missing | declined as bug; recorded as §6 question (scheduler guarantees the dependency via `Dependencies()`; see analysis) |
| C40 | `.cpp:284` stale "fall back" is correctness-adjacent | duplicate of C10 |

### claude.md

| # | Finding (loc) | Disposition |
|---|---|---|
| CL1 | `.h` clean — no change | declined (no-op) |
| CL2 | `.cpp:284-286` fallback banner contradicts code → "try APBS; nullptr on failure (no fallback)" | duplicate of C10 → E6 |
| CL3 | `.cpp:242-273` fused loop needs signposts (`// clamp runaway field`, `// kT/e → V/A`, `// store Mat3 + T2`) | duplicate of C9/C33 → E5 |
| CL4 | `.cpp:160-179` grid-sizing comment lists 4 facts, code computes 3; `grid_dim=161` orphaned below its block | adopted → E10 (move `grid_dim` decl up; align comment to what's computed) |
| CL5 | `.cpp:243-260` E gets clamp+zero, EFG gets only zero — unexplained asymmetry; one-line note | adopted → E14 (one-line note; see §4 for the reason) |
| CL6 | `.cpp:140,164-167` `xArr/yArr/zArr` → `x_coords/y_coords/z_coords` | duplicate of C16 → E8 |
| CL7 | `.cpp:164-171` `lo/hi` → `bbox_min/bbox_max` | duplicate of C17 → E9 |
| CL8 | `.cpp:81-90` add `// EFG(i,j) = dE_i/dr_j` at line 82 | duplicate of C14 → E3 |
| CL9 | `.cpp:35` `frac/fx/fy/fz/ix/iy/iz` clear — no change | declined (no-op) |
| CL10 | `.cpp:57-64` trilinear sum legible — no change | conflicts with C21; adopt the lighter C21 (one `// trilinear interpolation` header) → E12 |
| CL11 | `.cpp:81-119` build→symmetrize→detrace is the model; no change | declined (no-op; confirms C4/C5 should only *tighten*, not restructure) |
| CL12 | `.cpp:236-273` extraction loop stages not grouped | duplicate of CL3/C9 → E5 |
| CL13 | `.cpp:68,81` fn names say what they return — clean | conflicts with C26; supports declining C26 |
| CL14 | `.cpp:126` `ComputeViaApbs` clear — clean | declined (no-op) |
| CL15 | `.h:46-48` query names unambiguous — clean | conflicts with C36/C37; supports declining them |
| CL16 | `.cpp:93-104` symmetrization comment: keep physics, cut provenance line | duplicate of C5 → E4 |
| CL17 | `.cpp:39-40` `floor()` gotcha comment — keep | declined (no-op; explicitly keep) |
| CL18 | `.cpp:107-114` traceless comment clear, keep | duplicate of C6 → keep (declined to shorten) |
| CL19 | `.cpp:242` `// Sanitise` too terse | duplicate of C9 → E5 |
| CL20 | `.cpp:262-264` unit-conversion comment duplicates header; trim inline | adopted → E15 |
| CL21 | `.cpp:347-350` WriteFeatures T2 comment good but "line 102" stale | duplicate of C34 → E7 |
| CL22 | `.cpp:179,182-185` magic numbers (161, pdie 4.0, sdie 78.54, 0.15 M) commented but not parameterised; check 78.54 and 0.15 units | declined as edit; recorded as §6 question (values are standard APBS mg-auto; verify against bridge) |
| CL23 | `.cpp:45-47,71-73,84-86` interp returns 0 OOB → edge atom gets one-sided derivative vs zero; guard-comment that padding guarantees interior stencils | adopted → E16 (one-line note; see §4) |
| CL24 | `.cpp:251-254` clamp rescales E but leaves EFG unclamped — confirm intended | duplicate of CL5 → E14 (see §4 for reason) |
| CL25 | `.cpp:288-311` two "no data" signals (default result vs nullptr) — check consumers distinguish | declined as bug; recorded as §6 question (same as C39) |
| CL26 | `.h:43 / .cpp:43` `Interpolate` divides by `spacing(d)` no zero-guard | declined (bridge sets spacing > 0; low risk; no action) |

---

## 3. Edits that don't move numbers

All comment / signpost / internal-rename. No serialized name, no number.

- **E1** — `.cpp:21-23` — replace decorative banner + "from old ApbsSolver"
  history with `// Grid interpolation utilities (no APBS dependency).`
- **E2** — `.cpp:68` — above `ElectricFieldFromGrid`, add `// E = -grad(phi)
  by central difference.` (the `E = -grad(phi)` line at 74 stays).
- **E3** — `.cpp:82` — add one index-convention line:
  `// EFG(i,j) = dE_i/dr_j  (j = differentiation axis, i = field component).`
  Keep the existing block sequence; the symmetrize/detrace blocks already
  carry headers (do not restructure — CL11).
- **E4** — `.cpp:93-104` — shorten the symmetrization comment: keep the
  load-bearing physics ("FD construction can leave an antisymmetric residue;
  explicit symmetrization pins T1 = 0 by construction so the T2-only emit
  doesn't mistake interpolation noise for a physical pseudovector"); drop the
  "Per R4 codex review 2026-05-18" provenance line and the schema-rev date.
- **E5** — `.cpp:242-273` — replace `// Sanitise` with step signposts:
  `// finite-value guard (zero non-finite E + EFG)` before the NaN/inf checks;
  `// field magnitude cap` at the `E_mag > APBS_SANITY_LIMIT` branch;
  `// kT/(e·A) -> V/A ; kT/(e·A^2) -> V/A^2` at the conversion (merge with E15);
  `// store Mat3 + T2` at the `ca.apbs_*` assignments (replaces the
  "both Mat3 AND SphericalTensor" comment).
- **E6** — `.cpp:284-286` — replace the stale banner with
  `// Main compute: APBS solve only. No fallback — nullptr on failure.`
- **E7** — `.cpp:327-331` and `.cpp:347-350` — WriteFeatures comments:
  keep the accurate "T2 only, T0+T1 structural zeros" statement; drop the
  "PackST removed" archaeology and the stale "line 102" cross-reference;
  replace with "symmetrized + traceless-projected above".
- **E8** — `.cpp:140,143-147,164-170` — rename `xArr/yArr/zArr` →
  `x_coords/y_coords/z_coords` (file-local; all uses in `ComputeViaApbs`).
- **E9** — `.cpp:164-171` — rename `lo/hi` → `bbox_min/bbox_max`
  (file-local; the `extent = hi - lo` line follows).
- **E10** — `.cpp:160-179` — move the `int grid_dim = 161;` declaration up
  adjacent to the grid-sizing comment, and align the comment to the three
  quantities the code actually computes (`fine_dims`, `coarse_dims`,
  `grid_dim`); the "161 per axis" line and its variable become adjacent.
- **E11** — `.cpp:197` (+uses through 233) — rename `gridResult` →
  `apbs_grid` (file-local; passed to `apbs_solve`/`apbs_free_grid` by address).
- **E12** — `.cpp:57` — add `// trilinear interpolation (8 corner weights)`
  above the eight-term return.
- **E13** — `.h:42` — replace `// Factory:` lead with `// Compute APBS
  E-field and EFG at each atom from charges.`
- **E14** — `.cpp:250-255` — one-line note on the E-only cap:
  `// Cap is a finite-value guard on E only; EFG is left to the traceless
  projection + NaN guard (a clamped-E atom is already flagged anomalous).`
  (See §4 — confirm wording with the reason discovered.)
- **E15** — `.cpp:262-264` — trim the unit-conversion comment to a single
  inline line (it currently restates the header verbatim); folded into E5.
- **E16** — `.cpp:45-47` — one-line note at the OOB return:
  `// OOB -> 0; the >=20 A grid padding keeps every atom's +/-1-cell
  central-difference stencil interior, so atoms never see this zero.`

Carry-through note (not edited here): `src/ApbsEfgTimeSeriesTrajectoryResult.cpp`
lines 124-125 and 164-165 reference `ApbsFieldResult.cpp:262-272` and the H5
`source` attr references the same; these line numbers drift with E1-E16.
Flag for the human: either update them in the same pass or replace with a
non-numeric phrase ("traceless projection applied at source"). Out of this
algorithm's file scope; listed so it isn't silently dropped.

---

## 4. Usage notes (the sign/value reasons, traced)

### E-field sign: `E = -grad(phi)`, points away from positive charge
- **Producer** (`.cpp:75`): `E(d) = -(phi(+) - phi(-)) / 2h` = `-grad(phi)`.
- **Sibling producer** `CoulombResult.cpp:163`: `E_j = q_j * r / r^3` with
  `r = pos_i - pos_j` → field at the observation point points away from a
  positive source. This is `-grad(q/r)`. **Same physical convention** as APBS.
- **Consumer** `CoulombResult.cpp:303`: `coulomb_E_solvent = apbs_efield -
  E_total`. Solvated-minus-vacuum is only meaningful if both sides share the
  sign — they do. **Producer and consumer agree.**
- **Consumers** `python/_catalog.py:295` (`apbs_E`, irreps `1e`, parity `odd`,
  V/A), `ui/TrajectoryH5.cpp:314/407` (Vec3 frame-0 read), h5-reader inspector
  (`apbs_efield`, V/Å). All read it as a physical V/Å vector; none re-signs.
  **Coherent (expected).**

### EFG sign / construction: `EFG(i,j) = dE_i/dr_j`, symmetrized, traceless
- **Producer** (`.cpp:90`): `EFG(i,j) = (E+_i - E-_i)/2h` = `∂E_i/∂r_j`.
  Then symmetrized (`.cpp:105`), then trace removed (`.cpp:116`).
- **Sibling** `CoulombResult.cpp:166`: `V_ab = q(3 r_a r_b/r^5 - δ_ab/r^3)`
  = `∂E_a/∂r_b` analytically, then `project_traceless`. **Same convention,
  same projection.** The two EFGs are directly differenced in
  `CoulombResult.cpp:304` (`coulomb_EFG_solvent = apbs_efg - EFG_total`) —
  agreement is required and holds.
- **Symmetrize-before-detrace order**: symmetrization kills the antisymmetric
  FD-interpolation residue that would otherwise be decomposed as a spurious
  T1 pseudovector; the subsequent T2-only emit then carries no fake T1. This
  is the reason the long `.cpp:93` comment exists — keep the physics (E4),
  cut only the provenance.
- **Traceless projection** (`.cpp:107` comment): removes the discretized
  self-potential delta-function Laplacian so the emitted EFG is the external
  (other-atoms + reaction-field) gradient, which is genuinely traceless by
  Laplace. Comment is **accurate and load-bearing — keep it** (decline C6).
- **Consumers**: `python/_catalog.py:296` `apbs_efg` EFGTensor, 5-component,
  irreps `1x2e`, units `V/A^2`, **no sign_convention / no parity** — *identical
  to every sibling EFG spec* (`coulomb_efg_*`, `water_efg*`, `mopac_coulomb_efg_*`,
  `aimnet2_efg*`). The contract is uniform across all EFG producers.
  `ApbsEfgTimeSeriesTrajectoryResult` emits the 5 T2 components with H5 attrs
  `irrep_layout=T2_m-2..+2`, `parity=2e`, `units=V/Å^2`. h5-reader reads via
  `ReadT2TimeSeries` (5-component) and displays `V/Å²`. **All consumers agree;
  coherent (expected).**

### Catalog `parity` for EFG
- `apbs_E` carries `parity="odd"`; `apbs_efg` carries none. This is **not a
  gap** — every EFG spec in the catalog omits parity, while the TR H5 path
  tags `parity="2e"` (e3nn even, correct for the Hessian of a parity-even
  scalar). The two surfaces are different (SDK ArraySpec vs H5 attr) and
  neither contradicts the physics. **Coherent.**

### Unit conversion `KT_OVER_E_298K = 0.025693`
- Producer multiplies both E and EFG (`.cpp:265-266`) by `kT/e` to go from
  APBS-native kT/(e·Å^n) to V/Å^n. Catalog and h5-reader/ui both label V/Å and
  V/Å². Defined once in `PhysicalConstants.h:64`. **Coherent.**

### E-cap-but-not-EFG asymmetry (CL5 / CL24)
- The cap (`.cpp:252`) rescales runaway E to `APBS_SANITY_LIMIT = 100 V/Å`
  (`PhysicalConstants.h:95`). EFG gets a NaN/inf guard but **no magnitude cap**.
  Discovered reason: this is a *guard*, not a coupled physical rescale —
  it caps the stored E vector to a sane bound while the EFG is independently
  protected by (a) the traceless projection and (b) its own NaN/inf zeroing.
  The two are not rescaled as a pair on purpose (a clamped-E atom is already
  an anomaly the calibration absorbs; coupling the rescale would invent a
  physical relationship that doesn't exist). The reviews flagged this as
  "unexplained"; the fix is the **one-line note (E14)**, not a code change.
  I did not find a consumer that assumes E and EFG were jointly rescaled, so
  there is no reconciliation failure — see §6 for the one residual question.

---

## 5. Bug-by-exhaustion candidates

**None.** No sign/value finding survived tracing. Every E and EFG sign,
the unit conversion, the symmetrize/detrace order, and the catalog/H5
conventions reconcile across `src/` (CoulombResult difference),
`python/nmr_extract/_catalog.py`, `ui/`, and `h5-reader/`. The only
"wrong" artifacts are the stale comment banner (E6) and stale line-number
cross-references (E7) — comments wrong, code right; fixed by conforming
the comments to the code.

The two correctness items the reviews raised (empty-conformation guard,
missing-dependency return shape) are **defensive-shape questions, not
demonstrated bugs**, and are recorded in §6 rather than driving a code
change — the governing prior (exhaustion before bug) is not met: both
depend on framework guarantees I could not fully confirm from this file.

---

## 6. Questions & Ambiguities

1. **Empty conformation (C38).** `.cpp:164` seeds `bbox_min/max` from
   `xArr[0]` before any `n_atoms == 0` check. Does the result scheduler /
   `ProteinConformation` guarantee `AtomCount() >= 1` before a
   `ConformationResult::Compute` runs? If yes, no change. If a direct
   `Compute` on an empty conformation is reachable, an early
   `if (n_atoms == 0) return false;` is warranted. **Need confirmation of
   the framework invariant** before touching it (it would be a guard, not a
   number move).

2. **Missing-dependency return shape (C39 / CL25).** `.cpp:291-293` returns a
   default-constructed `ApbsFieldResult` (with `conf_ == nullptr`) when
   `ChargeAssignmentResult` is absent, but returns `nullptr` when APBS fails.
   `Dependencies()` declares `ChargeAssignmentResult`, so a correct scheduler
   never invokes `Compute` without it — making the `!HasResult` branch
   unreachable in production. Two questions: (a) is `Compute` ever called
   outside the dependency-respecting scheduler (tests? direct calls)? (b) if
   the branch is truly dead, should it return `nullptr` for consistency with
   the APBS-failure path, or is the default-result a deliberate "dependency
   absent is not an error" signal that some consumer distinguishes? I found
   no consumer that branches on the two cases differently, but could not
   prove the branch is dead. **Question, not a fix.**

3. **EFG-uncapped intent (CL5 / CL24).** §4 gives the discovered reason
   (independent guards, not a coupled rescale). One residual: confirm with
   the author that the E-cap is purely a storage sanity bound and that no
   downstream physics step expects E and its gradient to have been clamped
   consistently. If confirmed, E14's note stands; if not, this escalates.

4. **Grid magic numbers (CL22).** `sdie = 78.54` (vs the often-quoted 78.4
   for water at 298 K) and `ionic_strength = 0.15` molar — confirm these
   match what the `apbs_bridge` expects (units of ionic strength in
   particular: molar vs mol/L vs the bridge's internal convention). Values
   look standard for APBS mg-auto; this is a documentation/verification
   question, not an edit — the numbers stay as-is regardless.

5. **TR sibling line-number carry-through.** `ApbsEfgTimeSeriesTrajectoryResult.cpp`
   hard-references `ApbsFieldResult.cpp:262-272` and "line 102" (already stale
   today). Does the human want those updated in this pass (cross-file edit) or
   replaced with non-numeric phrasing? Listed in §3 carry-through; needs a
   scope decision.
