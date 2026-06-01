# Rediscovery map — what we recover, how deep, and what's left

The **science** roadmap, kept separate from the engineering surface
(`SURFACE_DESIGN.md`) and the instrument framing (`PLAN.md`). Reads with `STATE.md`
(freshest numbers) + `analysis/FINDINGS.md` + `analysis/PATTERNS.md` (consumer discipline).

**The goal:** recover the classical NMR shielding kernels from DFT — show the
geometric kernels CARRY the shielding signal (correlate-not-match), per-element /
per-atom-type. Not prediction, not R²-optimisation: evidence the kernel set is
complete enough.

## Two depths of "rediscovered"

- **Depth A — a law falls out (symbolic).** The kernel's closed form is recovered
  FROM the data: a permutation-invariant / equivariant fit learns the per-source
  function, it is distilled to a closed form (PySR, or reading the learned radial
  functions), and validated NON-circularly (leave-atoms-out vs independent DFT;
  literature-coefficient-fixed). The strong claim — the equation came out, not imposed.
- **Depth B — the signal is captured (correlate).** The kernels demonstrably carry
  the per-atom-type signal: a credible per-stratum R² against DFT. The kernel set is
  complete enough for that atom-type; a closed form may not be distilled yet.

A cell can be at B without A. **Taking a B-cell to A is the symbolic frontier (item 1).**

## Status grid (2026-06-01)

| relationship | target | Depth B (signal captured) | Depth A (law fell out) |
|---|---|---|---|
| **ring current** | T0+T2 | ✅ k≈21, held-out R²=0.62; equiv-T2 R²=0.47 / \|T2\| 0.76 | ✅ PySR → `intensity·(cos²θ−0.224)/r³` ≈ Pople; angular 0.99, axial 0.98 |
| **mcconnell** | T0+T2 | ◐ scalar R²≈0.55 | ◐ form recovers R²≈0.85; producer-kernel gap (~0.55) |
| **broad backbone** (ring+bond+field) | T2, 8 strata | ✅ (this session) HN .76 / O .72 / HA .67 / N .65 / C .59 / CA .43; axes help all | ◐ PARTIAL — bond McConnell r⁻³ recovers (atom-split 0.48); ring null; charge partial. **De-circularised** (un-fitted literature kernel T2 vs DFT, component r): **N 0.69, O 0.53, C 0.51** — textbook physics predicts the heavy-atom + amide-N tensors un-fitted; CA/HN/HA weak (need the fit). Honest multi-mechanism mixture (`BACKBONE_LAW_EVIDENCE.md`) |
| charge_dipole → field | T0/T1 | ◐ field_z r=0.46, LOAO 0.21 (μ null) | ❌ |
| buckingham_efield (APBS) | l=1 | ❌ stub | ❌ |
| **efg** (APBS) | T2 | ❌ stub — *natural next; same equivariant machinery; APBS data present* | ❌ |
| charge_quadrupole | T2 | ❌ stub | ❌ |
| larsen_hbond | T0+T2 | ❌ stub | ❌ |
| AIMNet2 CRG | T0+T2 | ❌ stub | ❌ |
| AIMNet2 embedding | per-atom | ❌ stub | ❌ |

## MAP ITEM 1 — the symbolic-distillation frontier (take B → A) — COMPLETES THE TEMPLATE

**Do this BEFORE item 2** (Jessica, 2026-06-01). The backbone exemplar is at Depth B
(signal captured); the template only demonstrates HALF the rediscovery arc until the
law is distilled. Completing it (B→A on the backbone) makes the exemplar the FULL arc
— substrate → equivariant fit (B) → symbolic distillation (A) → non-circular
validation — so every item-2 calculator-cell copies a complete template, not a
capture-signal-only half.

Where signal is captured (B) but no closed form (A), distill the law:
- PySR / symbolic regression on the emitted substrate (as ring's T0 did).
- Read out the e3nn model's learned **per-source-type radial functions** — the
  equivariant model already *is* a learned `g(r, invariants)` per source; interpret it.
- Non-circular gates: leave-atoms-out vs independent DFT; the **literature-coefficient-
  FIXED** test (an un-fitted constant → DFT).
- **First target: the broad-backbone T2** (captured this session). Does a backbone law
  fall out the way ring's did? Now cheap — the substrate + a clean e3nn fitter exist
  (`equiv_t2_backbone_e3nn.py`, the exemplar).

## MAP ITEM 2 — the kernel set across the per-element dimensions (breadth at B)

The thesis claims per-element dimensionality (H≈20, C≈6, N≈3, O≈12); each kernel
carries a piece. Bring each to Depth B on 1P9J — a new relationship is a stratum +
source-kind + form through the engine (`PLAN.md` Axis 1; the substrate is
relationship-agnostic). Order of attack:
- **efg** (APBS EFG → T2) — natural next; reuses the backbone exemplar's equivariant
  machinery; APBS always-on so data is present.
- **buckingham_efield** (APBS E-field → l=1).
- **charge multipoles** (dipole / quadrupole; ff14sb / aimnet2 / mopac as a parameter sweep).
- **larsen_hbond** (exchangeable-H stratum, donor/acceptor geometry).
- **AIMNet2 CRG**, **AIMNet2 embedding**.
Each rides the unified engine (#29) and copies the backbone exemplar. Cost = a
relationship definition + a fitter run, not a rebuild *(plus, per cell, any quantity it
needs that the substrate lacks — a C++ emit-extension, spine-side, never a Python recompute).*

## The per-cell workflow — what every item-2 calculator copies (the FULL arc)

The goal is the PHYSICS, not just the signal: a cell is NOT done when the R² is good
(Depth B) — it's done when the law is examined (Depth A) and de-circularised. Every
relationship copies this arc:

1. **Use and support the functional interface.** Express the cell through the engine as
   a composed `Relationship` (curried closures — stratum / selectors / attachers /
   reducer), NOT a procedural walk; extend the interface where the cell needs it and keep
   it composable — cleanly and aesthetically where possible. Don't bypass it or
   proliferate sibling runners (the broad case forced one — unify them, #29). Using the
   functional API well *is* the deliverable, not just the number it produces.
2. **Emit what the fitter + distillation need — case by case, spine-side.** Most cells
   need a C++ emit-extension: the substrate must emit that cell's geometry, its source
   orientation vectors, AND the fixed-literature kernel T2 its de-circularising test
   needs. This is EXPECTED, recurring, proper work — the spine (C++) computes the
   physics and emits it (read the H5 / typed objects + canonical basis; additive;
   broad-specific so the ring/mc oracle stays intact), NEVER a Python recompute
   (`analysis/PATTERNS.md`). Exemplars this session: orientation vectors (#32),
   literature-kernel T2 (#34). Flag the gap → extend the emit; don't hand-roll in Python.
3. **Capture the signal (Depth B)** — the e3nn fitter; copy `equiv_t2_backbone_e3nn.py`
   (heterogeneous sources, per-source-type radial, scatter-pool, frozen-C, per-stratum,
   effective N).
4. **Distill the law (Depth A)** — read the learned per-source-type radials + PySR
   (`backbone_distill_evidence.py` / `backbone_pysr_distill.py`). Does the law fall out?
5. **De-circularise** — the literature-coefficient-FIXED test: plug the textbook
   coefficient (no fitting), correlate vs DFT (needs the emitted literature-kernel T2
   from step 2; see #34).
6. **Report honestly** — correlate-not-match; what fell out, what didn't, and WHY (the
   backbone: bond recovers, ring doesn't, charge partial — an honest mixture, not forced).

**The exemplar set to copy:** the composed `Relationship` (engine) + `equiv_t2_backbone_e3nn.py`
(step 3) + the distillation scripts (step 4) + the literature-fixed test (step 5) + the
emit-extension pattern (step 2). The signal is the first step; the physics is the destination.

## Dependencies / order
- **Item 1 COMPLETES the template and comes BEFORE item 2** (Jessica, 2026-06-01).
  Going broad on a Depth-B-only template copies a half-exemplar across every
  calculator. Finish the backbone B→A first.
- **#29 engine unification** gates item 2's BREADTH (each new relationship rides ONE
  engine, not a new sibling runner).
- The **completed backbone exemplar** (substrate → equivariant fit → symbolic
  distillation → non-circular validation; `equiv_t2_backbone_e3nn.py` +
  `analysis/PATTERNS.md`) is what every item-2 cell copies.

## Constraints carried (don't drift)
- Correlate-not-match; kernels carry the signal — not prediction, not R²-chasing.
- T2 is the thesis (the tensor, not the scalar); per-element/per-atom-type complexity
  is the story (no simplification bias).
- Instantaneous model (geometry → shielding; the trajectory is a geometry sampler).
- Breadth via MECHANISMS on 1P9J (depth B across kernels), NOT many proteins —
  transferability across proteins is `PLAN.md` Axis 2 / Stage 2, deferred.
- **Integrity = credibility:** the model lives in C++; Python fits via e3nn on the
  emitted substrate; no recompute (`analysis/PATTERNS.md`). A rediscovery that ran on
  a hand-rolled projection is not a rediscovery.
