# Open questions — flat index across all 19 plans

Pointer only: every `Questions & Ambiguities` entry, one line each, tagged to
its plan. No re-judgment — read the plan's own section for the full reasoning.
Each plan's section lives at `reviews/applied-math/<algo>/plan.md`.

## Decisions — families (2026-05-24, Jessica)

1. **`SampleShieldingAt` → `SampleKernelAt`: DO IT** (it contradicts the
   kernel-not-shielding principle). No special dedicated commit — fold into the
   normal per-code fix work; rename across the ~7 calculators + viewer /
   h5-reader call sites.
2. **Grid-path filter divergence: FINDING — partly lazy.** Dropping
   `RingBondedExclusionFilter`/`SelfSourceFilter` is principled (no atoms at a
   free grid point). But the near-field cut differs: `Compute` rejects inside
   `2·radius` (`DipolarNearFieldFilter`), the grid path only inside `radius` —
   so the overlay renders the `radius…2·radius` shell the feature path treats
   as multipole-invalid, with no comment of intent. Looks like an incomplete
   port. Viewer-path behavior fix (moves grid values) → separate viz pass,
   flagged not blessed.
3. **`direction_to_center` / `_to_midpoint`: KEEP + clarifying comment.** No
   rename (would churn the UI graphics).
4. **`N==0` guard: FINDING — no guard, unenforced invariant.** `OperationRunner`
   runs every `Compute(conf)` with no atom-count gate; apbs/eeq/sasa index
   `q(0)`/÷N → UB on a 0-atom conf. Load paths don't produce empty proteins but
   nothing asserts it. **RESOLVED (Jessica): add the `if(N==0) return` guard** in
   apbs/eeq/sasa (number-neutral, conservative good practice).
5. **Dead fields: DELETE** `HBondKernelResult::f`, `FieldValue` (`source_index`
   + whole struct), `RingAccumulated`. **Why (blow-offs found):**
   - `HBondKernelResult::f` is the **McConnell angular scalar** `(3cos²θ−1)/r³`.
     McConnell accumulates its `f` into emitted per-category sums; HBond
     computes `f` and drops it (`hbond_scalars` = nearest_dist/inv_d3/count
     only). Clone of McConnell's kernel that never wired the
     accumulate-and-emit. **RESOLVED (Jessica): WIRE IT UP** — HBond emits the
     angular scalar like McConnell (`f`→per-atom sum→`hbond_scalars`). This is a
     **feature add** (new emitted column → SDK `_catalog.py` + wrapper), NOT a
     readability fix — separate small task; needs a one-line design of the
     HBond angular-scalar channel(s). `f` is therefore kept and wired, not deleted.
   - `FieldValue` (per-source field attribution; h5-reader design mirrors
     `QtFieldValue`) is **never constructed/emitted in src/** — the same
     richer-output idea already retracted. `CalculatorId` it uses stays (live
     via GeometryChoice). Safe to delete.

## Per-plan rulings (2026-05-24, Jessica)

**Group A — behavior/semantics**
1. HBond `hbond_inv_d3`: **keep atom→midpoint** length; comment the length.
2. Tripeptide `BackboneHA` predicate: **tighten to `|| (locant==Alpha &&
   element==Element::H)`** + comment. The locant clause is load-bearing for
   GLY HA2/HA3 (stamped `Locant::Alpha, BackboneRole::None` per Markley,
   `LarsenHBondShieldingResult.cpp:355`); the element guard closes the CA
   over-match (unreachable today; predicate validates the assigned `res.HA`).
   Behavior-neutral.
3. Larsen ACE `HasAllRequiredSlots`: **leave** (`size()==6`); fix the comment
   to match the actual check.
4. Coulomb MOPAC cutoff: **keep uncapped (full N²)** + comment. NOT a speed
   issue — the ~10-min MOPAC SCF dominates; full N² avoids the ff14SB cutoff's
   truncation error. No unify.
5. BiotSavart small-rho thresholds: **unify** the hardcoded `1e-10` (cpp:269)
   to the `near_zero_vector_norm_threshold` config key (cpp:291).
6. Saturated non-5-ring skip: **add a one-line `OperationLog` skip note**.

**Group B — cross-subproject rename go/no-gos**
7. `JBLobeOffset` → **`JohnsonBoveyLobeOffset`: do it** (src + h5-reader + learn).
8. Bonded accessors (`UBEnergy`/`ProperDihEnergy`/`ImproperDihEnergy`):
   **keep** (mirror the TR H5 channel order — contract-adjacent).
9. Coulomb `EFieldAt` prefix + `E_backbone_frac` SDK property: **keep**.

## Recurring (same question across several plans — resolve once)

- **`SampleShieldingAt` family name** (returns a kernel, not shielding; shared
  verbatim across ~7 calculators, called polymorphically by the viewer):
  dispersion Q1, mcconnell Q3, pi-quadrupole Q1, haigh-mallion Q3,
  larsen-hbond-grid Q1 (`QueryNearest` twin). → keep vs coordinated family rename.
- **Grid-path (`SampleShieldingAt`) omits `DipolarNearFieldFilter`** that
  `Compute` runs — deliberate (raw volumetric field) or should mirror Compute:
  dispersion Q2, mcconnell Q1, ring-susceptibility Q1.
- **`direction_to_center` / `direction_to_midpoint` field name inverts stored
  sense** (stores center/midpoint→atom): haigh-mallion Q1, ring-susceptibility
  Q2, mcconnell Q4 — repo-wide rename crosses 5 producers + UI JSON + h5-reader.
- **Empty-conformation `N==0` guard** (is non-empty an upstream invariant?):
  apbs-field Q1, eeq-result Q1, sasa Q1.
- **Dead fields** (keep-as-planned vs remove): ring-svd `RingAccumulated`
  (**RESOLVED → remove**), hbond `HBondKernelResult::f`, types
  `FieldValue::source_index` (whole `FieldValue` maybe unconsumed).

## Flagged earlier as "seemed real"

- **pi-quadrupole Q2** — `pq_per_type_T0` catalog units `Angstrom^-5`
  (`_catalog.py:170`) vs the Buckingham A scalar's `A^-4` (`.h:41`): deliberate
  convention or metadata typo. (output-metadata)
- **larsen-hbond-grid Q2** — does the ρ dense-grid axis store +180 separately
  from −180? periodic wrap is correct only if not; verify one archive's
  `rho_axis`.

## Per plan

**aimnet2-charge-response** — Q1: do model params still carry `requires_grad`
(backward accumulates unused param grads; harmless if so)? Q2: class rename
`…ChargeResponseGradientResult`→`…ChargeL2CoordinateGradientResult` — contract
surface (TR provenance strings) or internal? Q3: `ui/` "polarisability" labels
(MainWindow.cpp:1487/1490) — follow-up ui-scoped edit?

**apbs-field** — Q1: empty-conf guard [recurring]. Q2: missing-dependency
return shape (default-result vs `nullptr`) — branch dead / make consistent? Q3:
confirm EFG left uncapped is intended (E-cap is storage sanity only, nothing
downstream expects E+gradient clamped together). Q4: grid constants
`sdie=78.54`, `ionic_strength=0.15` — match `apbs_bridge` units? Q5: TR sibling
stale line-number refs — update in this pass or de-numericize?

**biot-savart** — Q1: two small-rho thresholds (cpp:269 hardcoded `1e-10` vs
cpp:291 config key) — intentional or oversight? Q2: literal `8` vs a
`RingTypeIndex::Count` named constant (+widen arrays/catalog)? Q3:
`HaighMallionResult.h:15` has the identical G-sign drop-the-minus (cross-file).

**bonded-energy** — Q1: CMAP φ-major axis order confirmed only by
self-consistency — optional cross-validate vs `GromacsEnergyResult.cmap_dih` on
a CMAP≠0 fixture. Q2: public accessor renames `UBEnergy`/`ProperDihEnergy`/
`ImproperDihEnergy` — contract-adjacent (TR H5 channel order). Q3:
`cmap_grid_spacing` field rename — weighed call (3 sites + `.h`). Q4: include
hygiene + `M_PI` portability — separate pass.

**coulomb** — Q1: `sources_beyond_cutoff` GeometryChoice *key* — downstream
contract (UI surfaces it) or free to rename? Q2: MOPAC variant omits the 20 Å
cutoff — deliberate or clone oversight? Q3: SDK property
`CoulombScalars.E_backbone_frac` clearer name (soft API surface)? Q4: disagreement
on `EFieldAt` etc. — `Total…` prefix or keep. Q5: relocate the two-kernels
header essay into `CoulombResult.h` — right home?

**dispersion** — Q1: `SampleShieldingAt` family name [recurring]. Q2: grid path
omits `DipolarNearFieldFilter` [recurring]. Q3: `kT2Components=5` — file-local
or shared in `Types.h`?

**eeq-result** — Q1: empty-conf guard [recurring]. Q2: verify "Caldeweyher 2019
Eq. 8" citation number. Q3: sign of `D4EeqParams.chi` — confirm via
`PhysicalConstants.h` (frozen). Q4: banner-vs-code naming preference (u/v/q in
banner, descriptive in code). Q5: `codex-correctness.md` not ledgered — **note:
it does exist in the dir; the agent missed it**.

**haigh-mallion** — Q1: `direction_to_center` name [recurring]. Q2: vertex-only
subdivision adequacy (numerical/accuracy — deferred to numerical-hygiene pass).
Q3: `SampleShieldingTensorAt` rename [recurring].

**hbond** — Q1: `hbond_inv_d3` uses atom-to-midpoint distance vs N···O extent —
confirm intended length for the ML feature (N···O⁻³ would move a number). Q2:
remove dead `HBondKernelResult::f` [recurring dead field].

**larsen-hbond-grid** — Q1: rename `QueryNearest` (it interpolates) [recurring].
Q2: ρ axis +180 storage [flagged real]. Q3: `AxisLookup` field rename
(`idx/idx_next/frac`) — file-local, reviews disagree. Q4: add one-line safety
comment on `has_acceptor_readouts` proxy?

**larsen-residue** — Q1: amide-detection windows (C=O<1.32 Å, C-N 1.20–1.50 Å)
hold for live tensorcs15 DFT geoms? (run `AllCombinationsPerceiveCleanly`). Q2:
ACE `HasAllRequiredSlots` also check `C_idx≥0 && O_idx≥0` (behavior change)? Q3:
disagreement on trimming the two big comment blocks (WL rationale, ambiguity).

**mcconnell** — Q1: `SampleShieldingAt` vs `DipolarNearFieldFilter` equivalence
[recurring]. Q2: carry-through of `CategorySum`/`NearestCOContribution`
public-method renames — confirm callers only src/+tests. Q3: `SampleShieldingAt`
convention [recurring]. Q4: `direction_to_midpoint` rename [recurring].

**pi-quadrupole** — Q1: `SampleShieldingAt` family rename [recurring]. Q2:
`pq_per_type_T0` units A^-5 vs A^-4 [flagged real]. Q3: E2 fused-term split —
taste call.

**planar-geometry** — Q1: `omega_deviation` range `[-π,π]` vs catalog `(-π,π]`
— normalize the `_catalog.py` string? Q2: h5-reader `pucker_theta` "rad" comment
(degrees on producer) — consumer-side bug, file with h5-reader owner? Q3:
saturated non-5-ring silent skip — add a log line (behavior)? Q4: `res_i` rename
taste (`res_ip`→`res_next` only, or also `res_i`→`res_curr`).

**ring-susceptibility** — Q1: `SampleShieldingAt` filter divergence (`<radius`
vs `2·radius`) [recurring]. Q2: `direction_to_center` rename [recurring].

**ring-svd-plane-fit** — Q1: degenerate-triangle orientation guard — document
caller contract (first-3 non-collinear)? Q2: `RingAccumulated` **RESOLVED →
remove**. Q3: `QtRing.h` lobe-offset literals drift vs config (cross-subproject
note). Q4: `JohnsonBoveyLobeOffset` rename go/no-go (~25 lines, 3 subprojects).

**sasa** — Q1: central validation for `sasa_n_points>0` or add a one-line
guard? Q2: optional perf hoist of `BondiRadius(S)` out of the loop — include or
leave (edges toward refactor)?

**tripeptide-pose-assembler** — Q1: disagreement on `H`→`cross_covariance` vs
keep `H` (sibling consistency with `RmsdTrackingTrajectoryResult`). Q2: `D`→
`reflection_fix` declined — confirm. Q3: `BackboneHA` slot predicate looser than
its name (`|| locant==Alpha`) — intentional tolerance or tighten (behavior)? Q4:
`EmitAlignedAtom`→`ValidateAndEmitCapAtom` — file-local, optional.

**types-spherical-tensor** — Q1: `FieldValue::source_index` (and maybe all of
`FieldValue`) unused — planned or dead [recurring dead field]? Q2: `T2CosineWith`
range clamp — no consumer asserts; decision. Q3: `RingTypeCode` vs
`RingTypeLabel` — confirm token.
