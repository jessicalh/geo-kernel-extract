# Codex brief — SPEC the Stage-3 learning model (the equivariant-conditioned GNN)

> **Historical session brief — not current truth (trued 2026-06-04).** Session
> provenance only; current truth is the relevant `SPEC_*`, `NOW.md`, and
> corrected `STATE.md`.

Status: **DRAFT — pending lead plan-vet, not yet fired.**

You own the grind; the lead vets + judges + owns ALL git. This is a **design-spec** pass. Deliverable: ONE
internal technical spec describing our Stage-3 learning model — the "ring-toss" equivariant-conditioned GNN.
It is the companion the advisor précis will summarize, so it must be precise and honest about built vs
designed vs open. Honest stakes: this is the forward arc of the thesis; a clean architecture spec is what
lets us build it deliberately once the chewer and the 720 ingest exist.

## CONTEXT (read first)
- `NEXT_SESSION_PROMPT.md` STEP 3 (the "ring-toss" predictor) + memory `project_stage3_equivariant_gnn`
  (the ethos flip + the three grounds) + memory `project_thesis_reporting_arc`.
- `DESIGN.md` + `GUIDANCE.md` (this dir) — why the substrate is method-agnostic: permutation-invariant over
  a variable source set, rotation-equivariant, tensor-valued (T2), a local sum σ_i = Σ_sources f(x_source);
  the equivariant sum-pooling shape; the un-summed source set + displacement vectors + target tensor.
- `PerAtomSubstrate.h` (the substrate contract) + `PER_ATOM_SUBSTRATE_SPEC_2026-06-02.md` (the column menu /
  the feature pile already emitted).
- The e3nn EXEMPLAR (the clean template to point at, NOT to re-author): `analysis/equiv_t2_efg_e3nn.py`,
  `analysis/equiv_t2_e3nn.py`, `analysis/change_of_basis.py`, `analysis/e3nn_protocol.py`; `analysis/FINDINGS.md`.
- `PHYSICS_ARCHITECTURE_UNIFICATION_2026-06-04.md` — the one-kernel `D_ab` picture (shared angular, per-type
  radial) that the architecture embodies.
- Memories (LAW): `feedback_frames_from_physics_not_tests` (e3nn equivariant → NO imposed per-atom local
  frame; emit raw geometry + tensor), `feedback_no_python_physics_except_labeled_integrity_test` (the model
  is e3nn, NEVER a hand-rolled second model), `feedback_t2_sacred` (never scalar-collapse), `project_aimnet2_win`,
  `feedback_model_is_spine` (per-source→GPU via the chewer / pybind11 over the C++ model).

## WHAT THE SPEC MUST DESIGN
- **The ethos flip, stated plainly:** Stages 1-2 are physics EXPLANATION (kernels carry signal; R² a
  diagnostic). Stage 3 is PREDICTION — **R² IS the metric** — while keeping the project's spine
  (equivariant, T2-valued, held-out, lean). Say why the MLP rejection (Stage 1) does NOT bind the
  equivariant path.
- **The feature pile (inputs / conditioners):** dihedrals + motion, the recovered kernel shadows (charge/
  ring/McConnell/field/H-bond), AIMNet2 embedding, MOPAC, OF3 embedding + OF3, geometry, all calculator
  outputs, Larsen DFTs (untried), boosts. "Understood-or-not" features are allowed here — prediction, not
  explanation. Map each to the substrate column that carries it (cite the spec).
- **The architecture:** equivariant sum-pooling (e3nn-class) over the variable source set; SHARED angular
  law, PER-TYPE radial (the `D_ab` one-object insight — independently re-derived by the physics-architecture
  review); conditioned on the feature pile; **T2-valued output** (the thesis; σ_iso/T0 alongside as the
  scalar comparison). **NO imposed per-atom local frame** — equivariance handles rotation. The angular map
  is e3nn (`1o×1o→2e` TensorProduct + scatter-add pooling), **never a hand-rolled second model**. Needs the
  **chewer** (per-source → GPU, pybind11 over the live C++ model) — name it as the dependency.
- **The three grounds (be honest about each):** TRAIN on **DFT** (the last stable ground; anchor here);
  TRANSFER on the **720-WT** (the statics transferability baseline — depends on the 720 ingest spec);
  SHOOT at **BMRB** experimental shift in the dark (an ensemble average our short ~15 ns MD never samples —
  informative, NOT clean validation; **interpolate-to-validate yes, to-train no**).
- **Protocol + scope:** held-out everywhere; R² the metric; lean. **MSc realism:** the deliverable-sized
  version is the per-mechanism dominance-clean equivariant fits + this conditioned predictor — NOT a single
  grand net over all calculators (the .tex "one network" is a picture, not the build).
- **Dependencies + status:** this is DESIGNED, not built; it needs the chewer + the 720 ingest first.
  Cross-reference `SPEC_720_STATIC_INGEST_2026-06-04.md` and `SPEC_INTERP_1P9J_SCOPING_2026-06-04.md`
  (the 1P9J interpolation model is the runnable down-payment / v0 of this).

## After you: an opus agent adversarially reviews this spec
An opus agent will check that the architecture honors equivariance / T2-sacred / no-imposed-local-frame /
e3nn-not-hand-rolled, that the ethos flip and the three grounds are honestly framed (especially the
BMRB-in-the-dark caveat), and will sharpen toward MSc realism (cut any drift to a grand all-calculator
net). Write it to survive that.

## HARD RULES
- **DOCS ONLY.** Write exactly ONE file:
  `/shared/2026Thesis/nmr-shielding/h5-reader/src/rediscover/SPEC_LEARNING_MODEL_2026-06-04.md`.
  Change NO code, run NO build / training, run NO `git`. Everything else READ-ONLY.
- **The model is e3nn; never propose a hand-rolled second model.** T2 is never scalar-collapsed. No imposed
  per-atom local frame. Python consumes the emitted substrate; the chewer binds the C++ model, not a Python
  re-model. **No protein/spatial model in Python — not even a secondary aggregate; the spatial work + the
  protein model live in C++, full stop, and never a terabyte Python dump.**
- Mark BUILT vs DESIGNED vs OPEN vs CAVEATED. Forks → lay out + recommend + flag for lead.
- Branch `h5-reader-pysr-spike` — never merge/switch/rebase/PR. **Lead owns ALL git.** Truthful, cited, no overclaim.
