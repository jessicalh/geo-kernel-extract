# The Controlling Spec — the forward build

*The spine of `doc/emerging/`. Per feature: the **T0/T1/T2 we want** + the **calculation
method, practically**, down to drill-implementable (hand to the standard "go build this
calculator" Result+NPYs drill). T2-preservation is the unit of the spec.*

The three Steps below **supersede the old three stages** (similar, but
re-framed as three complete deliverables; some run concurrently).

---

## Scope of the spec (pinned)

This document governs ONLY the **forward build** — the **new calculators, the revamped
(equivariant) model, and the revamped stats**. It is NOT laid back over existing settled work.

- **IN (the spec governs):** the MOPAC-extended capture; the revamp of the kept kernels
  (ring / McConnell / charge-EFG) into clean equivariant T0/T1/T2 emit; any genuinely new
  calculator (π-quadrupole if it earns a slot); the Step-1 statistical model and law work;
  the MatTen/e3nn tensor model; and the Step-3 shift model.
- **OUT (untouched, not re-governed):** the existing producer extractor as-is, the settled
  Stage-1 55-kernel ridge (weighted R² = 0.718), the current NPY/H5 emit surface, the blessed
  tests.
- **The mechanism keeping them apart is deprecate-and-add** — the old keeps running for the
  existing pipeline; the spec describes the *new* equivariant surface built **alongside** it.
  (MOPAC the one exception: run-once → consolidate in place.)

"Controlling document" = the spec for what we **newly build and emit**, never a retroactive
audit of what's done.

## Our novel geometry-encoding kernels are kept, not the defended contribution

**What these kernels are.** Our geometric kernels
— ring, McConnell, charge/EFG, and the others — are all the **same kind of object**: a low-order
(`l≤2`) **multipole moment of the weighted neighbourhood**, `K_i = Σ_j w_j (3 ρ̂ρ̂ − δ)/ρ³`. The
**geometry** is the substance; each source's chemistry enters **only** as the scalar weight `w_j` (a
charge for EFG, a ring-current intensity for ring, a `Δχ` for McConnell) — "thrown in." They are **not**
"physical-source tensors" separable from geometry. They **attempt to encode geometry**, and there are
many of them, not just McConnell.

**Disposition: our novel kernel architecture is not the defended contribution.** These geometry-encoding
kernels stay **built and emitted** (the producer surface is unchanged) — **not removed**. Defending the
maths of our novel encoding is too much work for this project's scope, so we do not feature it. They may
**find a use as we iterate** — reached for only if we are really reaching. This is **not** a claim that
Step 1 lacks tensor data or MOPAC (see below).

**What Step 1 actually relates.** Step 1's statistical model relates **our T0/T1/T2 tensor data to DFT
shielding**, with a **dihedral space**, **IUPAC backbone/sidechain identity**, and **pairwise residue
propinquity** as **modulators** (standard practice) — not as the primary inputs. **MOPAC is central — our
best predictor**: its electronic data and tensor contributions are core inputs. (Keep the lean predictive
MOPAC layer; the per-pose raw `O(N²)` matrices + the retained `.aux` are the disk bear, not the signal.)
Equivariant components are allowed where they help. Our novel `l≤2` geometric kernels are not featured —
feeding them into an `l≤2` e3nn is circular anyway (§7).

**PySR is unaffected.** The kernels are *lossy functions of the raw geometry* (§3), so the law work
(PySR, closed-form fits) on raw geometric inputs → the DFT shieldings has everything the kernels encode
**and more**, and is non-circular (handing PySR our pre-shaped kernels is the §7 circularity). Pulling
the kernels does not impair the law step.

## The bar — defensibility, not agent-confidence

A feature is **ADDED to the build only when the choice is defensible** — grounded in the literature,
the e3nn/equivariant field-standard, and our own code — *not* because an agent produced a confident
answer. The failure mode we will **not** repeat is **the back of the napkin**: choices that felt
right but couldn't be defended. A grounding agent's report is itself judged against this bar — a
thin, hand-wavy, or un-citable report is a *not-yet*, not a go. We must be able to **defend every one
of these choices**.

---

## The three Steps

Three **complete, independent** deliverables. The spec serves each independently and to a
**complete version** — no half-builds, no "fill in later." When this pass is done, writing
begins (modulo bug-fix runs). The three must not be conflated; they share machinery (the typed
model, the reader-side math, e3nn), not goals.

**Where equivariance lives.** Equivariance is **not mandated
as a predictor.** It is part of **Step 1's** statistical/law-analysis toolkit — a tool for *seeing
relationships* the angular structure carries (law recovery; structure a non-equivariant fit would blur).
For the **predictors (Steps 2 & 3)** equivariance is optional — included only where an ablation shows it helps, excluded
where it doesn't; where it helps, the model carries it. The dragon that forces this: the **solution
shift IS the `0e`**; the `2e` is orthogonal to it and averages out under tumbling — so the tensor /
equivariance earns its keep where the **target is tensor-valued** (Step 1's DFT tensor, Step 2's
tensor), not where it is the isotropic solution shift (Step 3). **The Step-2 and Step-3 predictors are
FIDO** — everything-that-works goes in the pot, ablatable. The geometry-encoding kernels are
kept and emitted, but not Step-1 inputs.

### Step 1 — The statistical model (primary deliverable; DFT-anchored)

**Goal.** A **statistical model** relating **our T0/T1/T2 tensor data to DFT shielding**, with a
**dihedral space, IUPAC backbone/sidechain identity, and pairwise residue propinquity as modulators**
(standard practice) — not as the primary inputs — and **equivariant components** allowed; **law work
folded in as appropriate**. **MOPAC is central — our best predictor** (its electronic data + tensor
contributions are core inputs; keep the lean predictive layer, not the per-pose raw matrices). It IS a
predictor (R² is the working metric); the law-finding within it is the explanation-flavoured part. Our
novel geometry-encoding kernel architecture is not featured here (see the note above).

**The law work + diagnostics it folds in, per stratum:**
- **R² with the angular and without it** — T2-in vs T0-only. Does the tensor content add
  signal? (The thesis question — "the angular residual *is* the thesis" — made into a
  deliverable.)
- **PySR applied broadly** — AI-*designed* helpfully (good search space, sane operators), **not
  AI-*gated* and hand-fussed**. Let it search wide for closed-form laws.
- **Equivariant analysis — the relationship lens.** Equivariance here is for *seeing relationships*
  the angular structure carries (law recovery; does an equivariant fit expose a closed-form law),
  **not prediction**. **Deliberately consider what we feed it** — feature selection in the service of
  understanding, not throw-everything-in. Gated on a clean, defensible path. **Distinct deliverable
  from Step 2's predictor** (may reuse reader-side math).
- **Traditional stats alongside PySR** — partial / joint / confounder-controlling fits, the standalone-
  vs-joint analysis; the full statistical toolkit, not only the AI methods.

**Strata:** backbone-all · residue-backbone · sidechain-IUPAC.
**Targets:** the **751** (1P9J trajectory) **+ 720** (WL / WT static) DFTs.

### Step 2 — The shielding-*tensor* predictor (a real model)

**Goal.** Predict the shielding **tensor** (T2 preserved) as well as we can. The honest
question is *how far we get*: **30 % / 20 % / 10 %** of the tensor. Here **R² IS the metric** —
this is prediction, not explanation.

**Model.** **FIDO** — everything-that-works in the pot, ablatable, **as powerful as possible**
(MatTen/e3nn is the live candidate). Inputs are earned by internal ablation; no kernel or scalar sits
in the model by default. **Equivariance is not mandated — but the target here IS a tensor, so it will likely earn
its place; include it where the ablation says it helps.** (The `2e` is the target at the DFT-anchored
*static* level — it survives here, unlike Step 3's solution shift.) Earn each piece by **internal
ablation** — equivariance on/off, feature-by-feature (the fido cull), tensor-vs-scalar — so nothing
sits in the model by default. **NOT measured against the Stage-1 ridge:** that is a typed
(per-element, per-role), static, mutations-stage result on a different target — wrong context, wrong
quantity — and we do not optimize R² regardless.

**Test.** First on **held-out 1P9J frames + the 720 WL DFTs**; then on **a few hundred
purpose-run small-protein test DFTs**.

### Step 3 — The shift predictor (anything-goes; ablatable)

**Goal.** Predict **experimental chemical shifts against BMRB / RefDB** (isotropic, scalar).

**Data.** **~600 short (1–5 ns) ML-MD runs**, with **MOPAC and everything *but* DFT** (DFT too
costly for the fleet).

**Model.** **Everything is allowed in** — message passing or whatever makes an **ablatable,
buildable** system.

---

## Feature roles

Three tiers, by what we claim and how each is used:

- **The Three — the featured geometry-encoding kernels.** Ring, charge/EFG, and McConnell are kept,
  emitted, and defended as calculations. They are **not** Step-1 inputs.
- **The Given set — taken as given, or scalar-led.** Established-tool outputs (MOPAC capture,
  AIMNet2, APBS), literature/scalar kernels we calculate but don't feature
  (Larsen H-bond, dispersion, π-quadrupole), and categoricals.
- **The cage** — H-bond, π-quadrupole, Larsen, dispersion — is kept and emitted, not updated, not
  featured, and not used in Step 1.
- **The models (Steps 2 & 3) are not limited to the above** — they take inputs at **any level, down
  to category and scalar, and *should*** (no-simplification-bias: use everything).

### The Three

The featured set is **three** kernels. Biot–Savart and HM-style count as **one** "ring"; the
two-path consistency check lives *inside* it, not as two headline members.

- **Ring** = **Biot–Savart** (`BiotSavartResult`; `0e ⊕ 2e` Johnson–Bovey / Boyd–Skrynnikov 2002
  eq. 3 tensor lift) **+ HM-style** (`HaighMallionResult`; `0e`). Agreement
  *by different math* is an internal consistency/regression guard, kept inside ring. (`RingSusceptibilityResult` is a candidate further cross-check, not a counted
  member.)
- **Charge / EFG** — `1o ⊕ 2e` (traceless-by-Laplace). One kernel, **source fork emergent**
  (FF14SB / MOPAC / AIMNet2 / APBS — map, don't decide). FF14SB stays (FullFat is the universal run
  shape → always-on); its cutoff is left as-is; the all-pairs/cutoff range difference is **not an EFG
  deliverable**. Code facts: `kernel_design/efg_reality_check.md`.
- **McConnell** — clean equivariant `D(r)·Q̂`, full even `0e ⊕ 1e ⊕ 2e`, two channels (fixed +
  Wiberg-BO). Form upgraded; scale **pinned** (held `Δχ`). `kernel_design/mcconnell_spec.md`
  + `mcconnell_integration_addendum.md` are the **template for the
  other two**.

**Shared source-based pattern (ring + McConnell): emit the unit kernel; the physical scale is PINNED
from literature/physics, not learned.** Our N — one protein within-instrument, 720-WT statics, thin strata — is too small to reliably *learn*
free coefficients. **Very little is actually fit.** The scale is pinned where a law exists, and the pinned
value's **magnitude and sign must be right** — the fit will not rescue them.
- Ring emits the **unit-current** kernel `G` (`bs_shielding`, `bs_per_type_T0/T2` — all `I = 1 nA`); the
  **literature ring-current intensity** (`ring.Intensity()`, per ring type, cited, verified `I = −12` PHE
  → `+1.40 ppm`) is the **pinned** physical scale, applied via the per-type seam (`σ = I·G`).
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

**Caged — kept, not deleted, not updated, not featured as recovered laws**:
**H-bond**, **π-quadrupole** (π-quad's EFG partition is an
ablation question, not ours to defend), **Larsen** (high complexity; possible later work),
**dispersion** (scalar-led).

**Standing roles:** **water** (`WaterFieldResult` — physical EFG source, rides with
charge/EFG); **ORCA + Tripeptide** (DFT shielding tensors — the *targets/anchors*); **AIMNet2**
(feeder + EFG source).

**Caution carried into the ring build:** ring-current has **historical sign issues throughout**;
sign-convention verification (`σ = I·G`, `I < 0` diamagnetic, the header's worked example) is a
**first-class checkpoint**, not an afterthought.

## The feature × value matrix

*Awaits Jessica's first grid.* **Rows = all candidate features**
(the role columns are what distinguish a responsibility-row from a participant-row). Shape (to
confirm against the grid): **per feature** — the **T0 / T1 / T2 we want**, its **dependencies**,
the **calculation**, and its **role in each Step separately**:

| feature | T0 want | T1 want | T2 want | dependencies | calculation | Step 1 role | Step 2 role | Step 3 role |
|---|---|---|---|---|---|---|---|---|
| *(grid drops in here)* | | | | | | | | |

*(dependencies = what it consumes — upstream calculators, MOPAC quantities, AIMNet2, geometry;
calculation = how it's actually computed, down to drill-implementable.)*

A feature can serve the three differently (statistical/law work in Step 1; tensor prediction in
Step 2; shift prediction in Step 3).

---

## Per-feature drill specs

*Filled per feature once the matrix lands.* Each becomes a "go build this calculator"
Result + NPYs drill. Start with the solid features (EFG, ring, McConnell — the wanted T0/T1/T2 +
method are nearly writable today); carry-alongs / partitions / which-MOPAC-quantity-becomes-
which-irrep wait on the clarifying examinations.
