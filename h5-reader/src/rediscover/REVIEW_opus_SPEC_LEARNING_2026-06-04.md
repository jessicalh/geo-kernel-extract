# Opus adversarial review — SPEC_LEARNING_MODEL_2026-06-04.md

> **Historical review record — not current truth (trued 2026-06-04).** Session
> provenance only; current truth is the relevant `SPEC_*`, `NOW.md`, and
> corrected `STATE.md`.

Reviewer: opus, adversarial-but-fair, against the codex brief and the repo.
Scope: docs-only. Verdict at top, evidence below. Branch git untouched.

## Verdict

**Spine: honored.** Equivariance / T2-sacred / no-imposed-local-frame /
e3nn-not-hand-rolled / no-Python-protein-model are all stated as hard
boundaries and matched against grounded code. No violations found.

**Ethos flip + three grounds: honest.** Stage-3-is-prediction / R²-is-the-metric
is clean (`:82-89`); the MLP-rejection-does-not-bind argument is correct and
load-bearing; the BMRB "in the dark, interpolate-to-validate-not-train" caveat
is stated three times and never softened (`:128`, `:134`, `:262`, `:279-280`).

**Grand-net drift: none — actively cut.** The spec explicitly refuses the single
all-calculator net (`:307-308`, `:312-315`) and stages v0→v5 (`:317-324`),
matching the physics-architecture review's "the .tex 'one network' is a picture"
(`PHYSICS_ARCHITECTURE_UNIFICATION_2026-06-04.md:68-71`).

This spec is precis-ready and accurate. The sharpenings below are refinements,
not corrections.

## 1. Spine check (line-by-line against code)

| Spine law | Spec statement | Grounded? |
|---|---|---|
| e3nn, never a 2nd hand-rolled model | `:45-49`, `:210-215` (`1o×1o→2e` `o3.TensorProduct`/`FullyConnected`, `Y2` from e3nn) | YES — matches `equiv_t2_e3nn.py:104-141` exactly (the live exemplar IS e3nn `o3.TensorProduct(1o,1o,2e)` + `o3.spherical_harmonics` + `index_add_` pool). |
| T2 sacred, never scalar-collapse | `:51-53`, `:246-249`, `:267-269` (T2 5-vec primary; T0 a comparison head "never as a replacement") | YES — consistent with `feedback_t2_sacred` and the exemplar's `(n,5)` target. |
| No imposed per-atom local frame | `:55-59`, `:251-254` (raw molecular-frame; physical source axes allowed *as source geometry*, not canonical target frames) | YES — matches `STATE.md:292-296` ("RAW lab-frame… NO imposed per-atom local frame") and `feedback_frames_from_physics_not_tests`. The distinction it draws (ring-normal/bond-axis are part of the source, not a frame attached to the target) is exactly right and is the subtle point a reviewer would attack — it survives. |
| No protein/spatial model in Python, not even a secondary aggregate | `:61-67` (Python never opens trajectory.h5, builds atoms/residues, KD-trees, discovers rings/bonds, decides self/bonded, assembles kernels; chewer binds C++ via pybind11) | YES — matches `feedback_model_is_spine`, `feedback_no_parallel_h5_in_python`, `GUIDANCE.md:24-33,54-56`. |
| The one legitimate Python transform | `:47-49` (`get_C()` is "a convention map on emitted 5-vectors, not a tensor projection") | YES — matches `change_of_basis.py` docstring verbatim in spirit; the orthogonality-below-`1e-12` gate it cites (`:286`) is real (`ALLATOM_FIT_SPEC_2026-06-03.md:476`). |
| No terabyte dump | `:69-73`, `:165-166` (lean row-aligned arrays + vectorized chewer query; no default source CSV) | YES — matches `PerAtomSubstrate.h:1-5` ("writes one row per (DFT frame slot, atom) and never emits a default source CSV"). |

**One latent spine risk to name (not a violation):** the radial-MLP boundary at
`:231-234` ("invariant quantities only… must not consume target values,
residuals, or held-out-derived statistics") is the place a future builder could
quietly break equivariance or leak — e.g. feeding a *vector* conditioner into
the radial net, or feeding a per-atom statistic computed on all frames. The spec
states the rule but does not name the failure mode. The exemplar enforces it
structurally: the radial MLP takes only `[r, ring_intensity]`
(`equiv_t2_e3nn.py:171`) and centering uses train rows only via `center_mask`
(`:132-152`). Recommend the spec point at that structural enforcement so "v2
conditioning" cannot drift into a vector-fed or leakage-fed radial gate. See
Sharpening 1.

## 2. Ethos flip + the MLP-rejection argument

Clean. `:82-89` states the flip without hedging and gives the correct reason the
Stage-1 MLP rejection does not bind: that MLP "was not an e3nn,
permutation-invariant, source-sum, T2-valued model over raw geometry… a
different model class with the wrong inductive bias." This is the right
argument and it matches `DESIGN.md:39-42` ("the Stage-1 MLP rejection says
nothing about it") and the project memory (`project_stage3_equivariant_gnn`:
"ETHOS FLIP… R² IS the metric").

Two honesty checks the spec passes that a skeptic would press:

- It does **not** quietly relax the project's "physics explanation, not
  prediction" cross-cutting rule into a blanket — it scopes the flip to Stage 3
  and keeps Stage 1-2 as explanation (`:77-82`). Good; that is the genuine
  position, not a contradiction of `CLAUDE.md`.
- It does **not** overclaim AIMNet2. `:135-136`, `:188` keep AIMNet2 as an
  allowed *conditioner* despite weak current T2 lift, explicitly because the
  goal changed from explanation to prediction — consistent with
  `project_aimnet2_win` without overselling it as a T2 rescue. Correct.

## 3. The three grounds

Honestly framed, and the dependency/uncertainty structure is correct.

- **DFT = train / last stable ground** (`:260`): correct anchor; the "do not
  train on dia/para split as the thesis target, total T2 is the stable target"
  caveat matches `GUIDANCE.md:139` and `STATE.md` audit verdict #1.
- **720-WT = transfer / the *only* between instrument** (`:261`): this is the
  single most important honesty point and the spec gets it exactly right. It
  cites `STATE.md:63-75`, which I verified says 1P9J true between-atom recovery
  is ~null (charge 0.036 / ring −1.0 / unified −105 overfit) and that
  "between / transferability / combine-DEPTH defers ENTIRELY to the 720-WT." The
  spec does not let 1P9J leave-atoms-out masquerade as transferability
  (`:276-277`: "do not call it transferability after the true-LOAO audit").
  This is precisely the trap the recent maths audit was about; the spec is
  audit-aware.
- **BMRB = shoot in the dark** (`:262`): the caveat is stated and never walked
  back — "Not clean validation: BMRB is an ensemble average, and the short MD
  ensemble is incomplete. Interpolate-to-validate against held-out DFT: yes.
  Train on BMRB: no." Matches the brief and `project_stage3_equivariant_gnn`
  ("shooting in the dark, informative not clean validation"). The split table
  reinforces it (`:279-280`: "no training split because it is not a training
  ground"). Honest.

## 4. MSc realism

Stays deliverable-sized; drift is cut, not merely disclaimed.

- `:299-310` names the deliverable as the per-mechanism dominance-clean
  equivariant fits + the conditioned predictor, with the 1P9J interpolation v0
  as the "runnable down-payment."
- `:312-315` explicitly rejects "a single all-calculator net over every emitted
  feature" as "too wide and too easy to over-read on current data."
- The v0→v5 tiering (`:317-324`) defers AIMNet2/Larsen/OF3/boosts to v3
  "only after explicit ablation and lead vet," and 720-WT to v4 — disciplined.

This is consistent with `PHYSICS_ARCHITECTURE_UNIFICATION_2026-06-04.md:68-71`
and the memory `project_stage3_equivariant_gnn` ("the deliverable-sized
version… NOT a single grand net"). No drift.

One realism gap worth flagging (Sharpening 2): the spec gives a build *order*
but no *stop condition* or expected-effort sense. For an MSc, "v2 conditioned
predictor with classical + MOPAC + selected geometry" is plausibly the thesis
ceiling; v3-v5 may not land. The spec should say which tier is the
minimum-viable thesis result so the build doesn't over-run chasing v5.

## 5. Honesty on what is BUILT vs DESIGNED vs OPEN

The BUILT/DESIGNED/OPEN/CAVEATED ledger (`:91-139`) is accurate and admirably
un-padded:

- It correctly marks the chewer and the GPU-scale source-set delivery as
  DESIGNED, not built (`:110-117`, `:328-341`).
- It correctly flags OF3-embedding and "boosts" as **OPEN / not grounded by any
  current column** (`:127`, `:192-193`, `:339`) and refuses to invent them —
  this is exactly right; those names came from the brief's feature wishlist and
  have no carrier in `PerAtomSubstrate.cpp`. Verified: grep finds no OF3/boost
  column.
- It honestly records that the two companion specs
  (`SPEC_720_STATIC_INGEST`, `SPEC_INTERP_1P9J_SCOPING`) are **not present
  locally** (`:36-41`, `:124`, `:335-338`) rather than pretending to
  cross-reference content it cannot see. Correct discipline.

Citation-accuracy spot check: the feature-pile column anchors (`:178-193`) and
the `ALLATOM_FIT_SPEC` leakage/basis refs (`:286`, `:291-292`) verify against
the actual files. Minor: some `PerAtomSubstrate.cpp` line numbers (e.g.
`:2870-2889`, `:3710-3733`) are approximate against the current 3979-line file,
but every *named* column/symbol they point at demonstrably exists (grep-checked
`dominant_fraction_*`, `gap_to_2nd_*`, `aimnet2_crg_scalar`, etc.). Treat the
column names as authoritative, the exact line numbers as soft.

## Sharpenings (recommended, none blocking)

1. **Name the equivariance/leakage failure mode at the radial boundary**
   (`:231-234`). State that the radial MLP takes *scalars only* (no vector
   conditioner) and is centered/normalized from train rows only, and point at
   the exemplar's structural enforcement (`equiv_t2_e3nn.py:115-118` radial on
   `[r,intensity]`; `:132-152` train-only `center_mask`). As written the rule is
   correct but a v2 builder could feed a vector conditioner or an all-frames
   statistic and silently break the spine. Make the enforcement structural, not
   a sentence.

2. **State the minimum-viable thesis tier.** The v0→v5 ladder (`:317-324`) is
   good but open-ended. Add one line: which tier is the floor the thesis stands
   on (likely v1 per-mechanism dominance-clean + v2 conditioned-on-DFT), and
   that v3-v5 are upside, not commitments. This keeps the build from chasing
   BMRB/720 before the DFT-ground result is the deliverable.

3. **Tie the angular-law claim to the data's verdict on *which* ring object.**
   `:217-220` ("shared angular `D_ab`, per-type radial") is well-grounded, but
   `PHYSICS_ARCHITECTURE_UNIFICATION_2026-06-04.md:62-65` notes BiotSavart/
   Haigh-Mallion are rank-≤1 outer products, NOT the symmetric-traceless `D_ab`
   shadow, and that the data already picked the clean ring path (jb/bs/hm agree
   +0.994, ringχ opposite). The spec's `:184` says "use JB as primary ring
   shadow, alternatives caveated by convention" — good — but the architecture
   section should explicitly note that the *shared angular `D_ab` basis* is the
   right shared object for the `D_ab`-family mechanisms (charge/McConnell/HBond),
   while ring-current enters through its own e3nn angular maps (Y2(r̂), Y2(n̂),
   cross term) as the exemplar already does. Otherwise "shared angular" can be
   misread as forcing the ring lobe object onto the `D_ab` basis. (Minor: the
   exemplar already does the right thing; the prose just under-specifies it.)

Optional nit: `:73` and `:165-166` both make the "no default source CSV" point;
the spec is slightly repetitive on the no-dump boundary (stated at `:69-73`,
`:165-166`, and implied in the architecture). Not wrong — these are load-bearing
— but it could consolidate.

## Precis-ready summary (4 sentences for the advisor's Stage-3 section)

> Stage 3 flips the project's stance from physics explanation to prediction: an
> equivariant, T2-valued sum-pooling network (e3nn) predicts the full shielding
> tensor of a held-out atom from its raw local source geometry plus
> "understood-or-not" conditioners (dihedrals, motion, the recovered classical
> kernel shadows, MOPAC/AIMNet2 fields), with held-out R² as the metric rather
> than a diagnostic. The architecture embodies the project's single-kernel
> insight — one shared angular tensor law, per-mechanism radial response — and
> keeps the spine intact: no imposed local frame (equivariance handles
> rotation), the spatial/protein model stays in C++ and is bound to GPU through
> the chewer, and the only Python transform is a frozen basis-convention
> constant. It is evaluated on three grounds of decreasing certainty: train on
> DFT (the last stable ground), transfer to the 720-WT static set (the only
> clean between-protein instrument, since 1P9J's between-atom signal is ~null),
> and finally a frozen-model shot at experimental BMRB shifts — an ensemble
> average our short MD never samples, so informative-but-not-clean validation
> (interpolate-to-validate, never train). The deliverable is the per-mechanism
> dominance-clean fits plus this conditioned predictor — a staged build, not a
> single grand network over every calculator — and it is currently DESIGNED, not
> built, blocked on the chewer (pybind11 over the live C++ model) and the 720
> static ingest.
