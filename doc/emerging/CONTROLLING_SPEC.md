# The Controlling Spec — the forward build

*The spine of `doc/emerging/`. Per feature: the **T0/T1/T2 we want** + the **calculation
method, practically**, down to drill-implementable (hand to the standard "go build this
calculator" Result+NPYs drill). T2-preservation is the unit of the spec.*

Settled 2026-06-06. The three parts below **supersede the old three stages** (similar, but
re-framed as three complete deliverables; some run concurrently).

---

## Scope of the spec (pinned)

This document governs ONLY the **forward build** — the **new calculators, the revamped
(equivariant) model, and the revamped stats**. It is NOT laid back over existing settled work.

- **IN (the spec governs):** the MOPAC-extended capture; the revamp of the kept kernels
  (ring / McConnell / charge-EFG) into clean equivariant T0/T1/T2 emit; any genuinely new
  calculator (π-quadrupole if it earns a slot); the e3nn / SEGNN model; the
  rediscover → PySR → equivariance → R² stats sweep and its guardrail.
- **OUT (untouched, not re-governed):** the existing producer extractor as-is, the settled
  Stage-1 55-kernel ridge (weighted R² = 0.718), the current NPY/H5 emit surface, the blessed
  tests.
- **The mechanism keeping them apart is deprecate-and-add** — the old keeps running for the
  existing pipeline; the spec describes the *new* equivariant surface built **alongside** it.
  (MOPAC the one exception: run-once → consolidate in place.)

"Controlling document" = the spec for what we **newly build and emit**, never a retroactive
audit of what's done.

## Our geometry-encoding kernels are being pulled from the defended thesis (2026-06-08)

**What these kernels are (corrected; `GEOMETRIC_KERNEL_ACCOUNT.md` §2–3, §5).** Our geometric kernels
— ring, McConnell, charge/EFG, and the others — are all the **same kind of object**: a low-order
(`l≤2`) **multipole moment of the weighted neighbourhood**, `K_i = Σ_j w_j (3 ρ̂ρ̂ − δ)/ρ³`. The
**geometry** is the substance; each source's chemistry enters **only** as the scalar weight `w_j` (a
charge for EFG, a ring-current intensity for ring, a `Δχ` for McConnell) — "thrown in." They are **not**
"physical-source tensors" separable from geometry (an earlier note this session said that — it was
wrong). They **attempt to encode geometry**, and there are many of them, not just McConnell.

**Disposition:** these geometry-encoding kernels are being **pulled from the defended thesis** —
defending the maths of our novel encoding is too much work for this project's scope, whatever its real
intellectual interest. They stay **built and emitted** (the producer surface is unchanged); they are
just not what the thesis defends.

**What stands in their place for the models:** raw geometry, which **e3nn encodes itself** (its own
`l≤2` machinery — feeding our `l≤2` kernels into it is circular anyway, §7), against the **DFT shielding
tensors** as the target, with **AIMNet2** and other genuinely-physical tensors available as features —
not our multipole kernels.

*OPEN (Jessica): the exact disposition scope — does this pull the kernels from the **law study** too
(and what does Part 1 then characterise), or only from the predictive models? Pending your call. Do not
propagate the older "fido-bits for the models" / "Inputs: the four kernels" language as settled
meanwhile.*

## The bar — defensibility, not agent-confidence (pinned 2026-06-06)

A feature is **ADDED to the build only when the choice is defensible** — grounded in the literature,
the e3nn/equivariant field-standard, and our own code — *not* because an agent produced a confident
answer. The failure mode we will **not** repeat is **the back of the napkin**: choices that felt
right but couldn't be defended. A grounding agent's report is itself judged against this bar — a
thin, hand-wavy, or un-citable report is a *not-yet*, not a go. We must be able to **defend every one
of these choices**.

---

## The three parts

Three **complete, independent** deliverables. The spec serves each independently and to a
**complete version** — no half-builds, no "fill in later." When this pass is done, writing
begins (modulo bug-fix runs). The three must not be conflated; they share machinery (the typed
model, the reader-side math, e3nn), not goals.

**Where equivariance lives (settled 2026-06-07, dragon-informed).** Equivariance is **not mandated
as a predictor.** It is **Part 1's analysis lens** — the tool for *seeing relationships* the angular
structure carries (law recovery; structure a non-equivariant fit would blur). For the **predictors
(Parts 2 & 3) equivariance is optional — included only where an ablation shows it helps, excluded
where it doesn't; where it helps, the model carries it.** The dragon that forces this: the **solution
shift IS the `0e`**; the `2e` is orthogonal to it and averages out under tumbling — so the tensor /
equivariance earns its keep where the **target is tensor-valued** (Part 1's DFT tensor, Part 2's
tensor), not where it is the isotropic solution shift (Part 3). **Both predictors are FIDO** —
everything-that-works goes in the pot, ablatable, chunky (now seasoned: the features are grounded and
cited). The four kernels are built regardless — responsibility-features, Part-1 analysis subjects,
and fido-bits for the models.

### Part 1 — The law / correlation study (DFT-anchored; standalone)

**Goal.** Characterise what *every sane metric we collect* carries against DFT shielding — the
signal-and-laws audit. This is **correlation and law**, not a predictor; keep it unconflated
with Parts 2–3.

**For each metric × stratum, produce:**
- **R² with the angular and without it** — T2-in vs T0-only. Does the tensor content add
  signal? (The thesis question — "the angular residual *is* the thesis" — made into a
  deliverable.)
- **PySR applied broadly** — AI-*designed* helpfully (good search space, sane operators), **not
  AI-*gated* and hand-fussed**. Let it search wide for closed-form laws.
- **Equivariant analysis — the relationship lens.** Equivariance here is for *seeing relationships*
  the angular structure carries (law recovery; does an equivariant fit expose a closed-form law),
  **not prediction**. **Deliberately consider what we feed it** — feature selection in the service of
  understanding, not throw-everything-in. Gated on a clean, defensible path. **Distinct deliverable
  from Part 2's predictor** (may reuse reader-side math).
- **Traditional stats alongside PySR** — partial / joint / confounder-controlling fits, the standalone-
  vs-joint analysis; the full statistical toolkit, not only the AI methods.

**Strata:** backbone-all · residue-backbone · sidechain-IUPAC.
**Targets:** the **751** (1P9J trajectory) **+ 720** (WL / WT static) DFTs.

### Part 2 — The shielding-*tensor* predictor (a real model)

**Goal.** Predict the shielding **tensor** (T2 preserved) as well as we can. The honest
question is *how far we get*: **30 % / 20 % / 10 %** of the tensor. Here **R² IS the metric** —
this is prediction, not explanation.

**Model.** **FIDO** — everything-that-works in the pot, ablatable, **as powerful as possible**
(probably a GNN). Inputs: the four kernels + electronic scalars + categoricals + anything that earns
its place. **Equivariance is not mandated — but the target here IS a tensor, so it will likely earn
its place; include it where the ablation says it helps.** (The `2e` is the target at the DFT-anchored
*static* level — it survives here, unlike Part 3's solution shift.) Earn each piece by **internal
ablation** — equivariance on/off, feature-by-feature (the fido cull), tensor-vs-scalar — so nothing
sits in the model by default. **NOT measured against the Stage-1 ridge:** that is a typed
(per-element, per-role), static, mutations-stage result on a different target — wrong context, wrong
quantity — and we do not optimize R² regardless.

**Test.** First on **held-out 1P9J frames + the 720 WL DFTs**; then on **a few hundred
purpose-run small-protein test DFTs**.

### Part 3 — The shift predictor (anything-goes; ablatable)

**Goal.** Predict **experimental chemical shifts against BMRB / RefDB** (isotropic, scalar).

**Data.** **~600 short (1–5 ns) ML-MD runs**, with **MOPAC and everything *but* DFT** (DFT too
costly for the fleet).

**Model.** **Everything is allowed in** — message passing or whatever makes a **nice, ablatable,
cool, buildable** system.

---

## Feature tiers (settled 2026-06-06)

Three tiers, by what we claim and how each is used:

- **The Four — our academic-responsibility kernels.** Four features we **take academic
  responsibility for** and **calculate ourselves** (perhaps from existing tools' outputs, but the
  kernel is ours). They go through the **full equivariant** T0/T1/T2 treatment (**required**), and
  may **also** be used in **scalar** form for the equation / PySR attempts. They are the
  **equation-extraction targets** of Part 1.
- **The Given set — taken as given, or scalar-led.** Established-tool outputs (MOPAC capture,
  AIMNet2, APBS), literature/scalar kernels we calculate but don't take novel responsibility for
  (Larsen H-bond, dispersion, π-quadrupole), and categoricals — participants in the Part-1
  statistics and **positive-control references** (the three-tool stool).
- **The Four + the Given set constitute Part 1 in whole.** Both also **participate** in Parts 2 & 3.
- **The models (Parts 2 & 3) are not limited to the above** — they take inputs at **any level, down
  to category and scalar, and *should*** (no-simplification-bias: use everything).

### The Three (refined 2026-06-07 from "the Four")

The academic-responsibility set we build and defend is **three** kernels. (Refines the
2026-06-06 "Four": Biot–Savart and Haigh–Mallion now count as **one** "ring" — the two-path
validation lives *inside* it, not as two headline members.)

- **Ring** = **Biot–Savart** (`BiotSavartResult`; `0e ⊕ 2e` Johnson–Bovey / Boyd–Skrynnikov 2002
  eq. 3 tensor lift) **+ Haigh–Mallion** (`HaighMallionResult`; decades-calibrated `0e`). Agreement
  *by different math* = the **two-path validation** (`feedback_two_path_validation`), kept, now
  internal to ring. (`RingSusceptibilityResult` is a candidate further cross-check, not a counted
  member.)
- **Charge / EFG** — `1o ⊕ 2e` (traceless-by-Laplace). One kernel, **source fork emergent**
  (FF14SB / MOPAC / AIMNet2 / APBS — map, don't decide). FF14SB stays (FullFat is the universal run
  shape → always-on); its cutoff is left as-is; the all-pairs/cutoff range difference is **not an EFG
  deliverable**. Code facts: `kernel_design/efg_reality_check.md`.
- **McConnell** — clean equivariant `D(r)·Q̂`, full even `0e ⊕ 1e ⊕ 2e`, two channels (fixed +
  Wiberg-BO). Form upgraded; scale **pinned** (held `Δχ`). **BUILT + adversarial-passed** —
  `kernel_design/mcconnell_spec.md` + `mcconnell_integration_addendum.md` are the **template for the
  other two**.

**Shared source-based pattern (ring + McConnell): emit the unit kernel; the physical scale is PINNED
from literature/physics, not learned.** Corrected 2026-06-07 ([[feedback_fittable_law_is_the_calibration]]):
our N — one protein within-instrument, 720-WT statics, thin strata — is too small to reliably *learn*
free coefficients, so the earlier "scale rides separately / fit downstream" framing overstates what the
data supports. **Very little is actually fit.** The scale is pinned where a law exists, and the pinned
value's **magnitude and sign must be right** — the fit will not rescue them.
- Ring emits the **unit-current** kernel `G` (`bs_shielding`, `bs_per_type_T0/T2` — all `I = 1 nA`); the
  **literature ring-current intensity** (`ring.Intensity()`, per ring type, cited, verified `I = −12` PHE
  → `+1.40 ppm`) is the **pinned** physical scale, applied via the per-type seam (`σ = I·G`). The gold
  case: a found fit space + de-circularised coefficient.
- McConnell emits the **unit-shape** `Q̂`; its scale `Δχ` is **pinned** from a held literature value, not
  fit. The 2–5× source scatter is handled by **picking a defensible value and owning it** (a sensitivity
  note if it bites), not by deferring to a fit. The rhombic C=O `Δχ` is now held (Hooper 1965 / Abraham
  1999 — `kernel_design/mcconnell_rhombic_spec.md`), so have-it-to-cite-it is **load-bearing physics**,
  not citation hygiene.
- The small amount that *is* learned is the residual ensemble combination (arc layer 3), **not** the
  per-kernel scales (arc layer 2 = the pinned laws). `feedback_emit_is_not_a_limiter` honoured as the
  *pattern* (emit unit + apply the physical scale); the scale is **pinned**, not "the direction-4
  hypothesis the fit refines." (The rediscover layer carries the literature-scaled `jb_T*` via
  `LiteratureIntensity()`, separate from the unit-current producer emit.)

The governing rule the three embody: **keep multiple paths only when they are independent routes to
the same quantity (ring: BS↔HM); replace when one strictly dominates (McConnell: lose the cruft);
leave open when undecided (charge/EFG: emergent source fork).**

**Caged — kept, not deleted, not updated, not featured as recovered laws**
(`project_three_kernels_and_cage`): **H-bond**, **π-quadrupole** (π-quad's EFG partition is an
ablation question, not ours to defend), **Larsen** (high complexity — a possible *late flourish if we
wish*, not now), **dispersion** (scalar-led).

**Standing roles, not darlings:** **water** (`WaterFieldResult` — physical EFG source, rides with
charge/EFG); **ORCA + Tripeptide** (DFT shielding tensors — the *targets/anchors*); **AIMNet2**
(feeder + EFG source).

**Caution carried into the ring build:** ring-current has **historical sign issues throughout**;
sign-convention verification (`σ = I·G`, `I < 0` diamagnetic, the header's worked example) is a
**first-class checkpoint**, not an afterthought. The three-kernel rework **will break** golden/smoke +
consumers — taken deliberately at the re-bless.

## The feature × value matrix

*Awaits Jessica's first grid.* **Rows = all features — the Four *and* the Given participants**
(the role columns are what distinguish a responsibility-row from a participant-row). Shape (to
confirm against the grid): **per feature** — the **T0 / T1 / T2 we want**, its **dependencies**,
the **calculation**, and its **role in each part separately**:

| feature | T0 want | T1 want | T2 want | dependencies | calculation | P1 role (metric to correlate) | P2 role (eq / scalar input) | P3 role (anything-goes input) |
|---|---|---|---|---|---|---|---|---|
| *(grid drops in here)* | | | | | | | | |

*(dependencies = what it consumes — upstream calculators, MOPAC quantities, AIMNet2, geometry;
calculation = how it's actually computed, down to drill-implementable.)*

A feature can serve the three differently (a metric to correlate in P1; an equivariant-or-scalar
input in P2; an anything-goes input in P3).

---

## Per-feature drill specs

*Filled per feature once the matrix lands.* Each becomes a "go build this calculator"
Result + NPYs drill. Start with the solid features (EFG, ring, McConnell — the wanted T0/T1/T2 +
method are nearly writable today); carry-alongs / partitions / which-MOPAC-quantity-becomes-
which-irrep wait on the clarifying examinations.
