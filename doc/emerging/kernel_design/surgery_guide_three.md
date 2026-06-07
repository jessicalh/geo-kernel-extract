# Surgery guide — the three (ring · charge/EFG · McConnell)

**What this is.** Not a physical-tensor audit — the **operating guide** for the three-kernel surgery:
*what we do differently now* that the 2026-06-07 session taught us. Each kernel's structured grounding
opens to this. Provenance: the law-study reframe ([[project_lawstudy_reframe]], directions 1 + 4), the
dragon ([[project_solution_shift_dragon_and_fido_reframe]]), the EFG reality-check
(`efg_reality_check.md`), the T2-physical-tensor audit, and the caging ([[project_three_kernels_and_cage]]).

The happy-tour claim: do these eight differently and the surgery is mostly **naming, switching on what's
already computed, and one careful sign pass** — not a rebuild. McConnell already proved the cadence
(reality-check + addendum → clean first-try); ring and EFG are *less* invasive than McConnell was.

## The eight cross-cutting rules (do differently)

1. **Emit the unit kernel; the physical scale rides separately — never fabricate a magnitude.** First
   classify each kernel **source-based** (ring, McConnell → unit geometric kernel + a scale that is
   literature-pinned or learned) vs **target-based** (EFG → *is* the physical field; no scale to ride).
   That classification decides whether there is a "scale surgery" at all.
2. **The kernel is a hypothesis to check the model against — not a feature to feed it.** Stop optimizing
   each kernel for "passes as a model input" (it doesn't — don't-double-feed) and optimize for "clean,
   honest, *checkable* closed-form claim." The model on raw geometry finds the law; the kernel is the
   closed-form hypothesis its ablation-discovery is compared to (direction 4).
3. **Report the relationship at honest confidence; drop the physical-pedigree theater** (direction 1).
   The bar is "this descriptor carries this much DFT signal, angular-in vs -out, with these confounds" —
   not "we recovered the textbook tensor."
4. **Per-part altitude (the dragon).** The three feed Part 1 + Part 2, both tensor-targeted → carry the
   full `0e ⊕ 1e ⊕ 2e`. T2 stays sacred as DFT-angular evidence, **never** sold as a solution-shift
   predictor.
5. **Confirm-and-refine ≠ rebuild.** Ring and EFG already compute the physics — the surgery is *naming,
   parity, dropped channels, sign, the two-path benchmark*, **not new physics**. **Front-load the code
   reality-check** (EFG's is done; ring's isn't) so the structured grounding writes against verified
   facts — that is what made McConnell go clean first-try.
6. **Verify, don't assert: sign, parity, frame are first-class checkpoints.** Ring sign especially
   (historical issues throughout). Parity per-quantity (EFG E-field is genuinely `1o`; shielding even,
   never `1o`). The mirror-reflection test is the runtime guard.
7. **Emit separably where cheap, lean on disk, decide downstream.** Turn the silently-dropped channels
   (sidechain EFG, partitioned E) into separable emits — you cannot conjure a channel you never wrote,
   and re-extract is <1 day — but stay <15 GB; drop downstream, not at the emit.
8. **Don't expand scope, don't prematurely resolve.** Source fork stays emergent (map, don't decide); the
   cage stays shut (H-bond / π-quad / Larsen / dispersion kept-not-updated); π-quad's partition is an
   *ablation*, never a code-grounds drop.

## Per-kernel: where each sits, and what the surgery actually is

- **McConnell — already built.** Surgery is light: the deliberate re-bless, the two tidy-ups (DRY
  manifest, layout-string), the X–H ablation (Part 1 decides keep/drop vs DFT), and — optional, later —
  the tabulated-`Δχ` physical hypothesis (rule 1's source-scale, scattered 2–5×). Stitches, not an incision.
- **EFG — confirm-and-refine; reality-check DONE** (`efg_reality_check.md`). Surgery: fix the
  `coulomb_shielding` misnaming; switch on the dropped channels (sidechain EFG, partitioned E) as
  separable; the catalog `1e→1o` parity-metadata fix; keep the source fork separate/emergent; **FF14SB
  untouched** (FullFat universal → always-on; cutoff as-is; range difference not an EFG deliverable). No
  new physics.
- **Ring — confirm-and-refine; reality-check NOT yet done, and the one that wants the most care.** Surgery:
  **sign-convention verification first** (against Boyd–Skrynnikov + the PHE worked example
  `I=−12 → +1.40 ppm`), then naming/parity; confirm the per-type seam pins-or-learns; the BS-`2e` ↔
  HM-`0e` two-path benchmark (HM `T0` vs the published bond-sum); unit-current emit + literature-intensity
  scale (rule 1). This is where "we will break things" and "sign issues throughout" both live → it gets
  the reality-check *and* a dedicated sign pass before anything moves.

## Resolved + corrected by codex confirm (2026-06-07)

A read-only codex pass (gpt-5.5/xhigh, ~402k tokens) confirmed the rules and resolved the to-verify items.

**Ring pin-vs-learn → likely LEARNS (not pins).** The top-level fitting code isn't present in `learn/`
(only `extract.py`, docs, logs); visible evidence (ridge coefficients as calibrated constants;
`calculator_params.toml` written by calibration) **leans learned**. The literature intensity is cited, is
the physical scale in `σ = I·G`, and the h5-reader viz uses it — but the calibration coefficient is *fit*,
not pinned. **So ring and McConnell are closer mirrors: both emit unit + learn the scale.** (Spec corrected.)
The rediscover layer *does* carry a literature-scaled `jb_T*` via `LiteratureIntensity()` — separate from
the unit-current producer emit; that's where `feedback_emit_is_not_a_limiter`'s pattern actually lives.

**Dropped-channel inventory (rule 7, full list).** FF14SB: sidechain EFG (computed, not written),
partitioned E (only total written), APBS solvent deltas, aromatic bond-projection + sidechain-source count.
MOPAC: same drops. AIMNet2: no E-field, no separate sidechain EFG. **Ring: per-type T1, per-ring
`B_field`/`B_cylindrical`, HM `hm_B_field`.** McConnell: near-field accepted/rejected counts. Re-extract is
cheap → turn the wanted ones into separable emits.

**Ring sign history → real, currently consistent in the producer.** Historical evidence exists
(`reviews/applied-math/biot-savart/plan.md` calls the minus "load-bearing"; HM reviews flagged a missing
minus). Producer code + tests are sign-consistent. **The live risk is in CONSUMERS/METADATA, not the
kernel:** the BS H5 time-series copies the unscaled unit kernel but labels units `ppm` (should be
`ppm_T_per_nA`); h5-reader viz multiplies `LiteratureIntensity()`; `_ring.py` describes HM as "intensity·H"
though C++ emits the unscaled `G`. Never mix these as the same physical quantity.

**Two-path (BS↔HM) → exists as a TEST, not a production emit.**
`tests/test_batch_biot_savart_haigh_mallion.cpp` computes per-ring + distance-binned BS/HM convergence and
the T2 cosine, explicitly "report, don't assert." There is **no emitted** two-path result/H5/group. If the
thesis wants it reported, that's a small to-build (surface the test's comparison as an emit).

**Build requirement (Jessica 2026-06-07): every two-path validation ships as an ASSERT test, in the same
build — never report-only.** An un-asserted two-path rots silently: a future sign/units/parity break prints
different numbers and still passes green → exactly the fast loop-back we must avoid. The tolerance is
**derived, not magic** — the test already prints `random=0.36, parallel=1.0`, so assert BS–HM `|cos|` sits
near *parallel* and well above the 0.36 random floor, and the negative controls (BS-vs-unrelated) stay at
~0.36. Flip `test_batch_biot_savart_haigh_mallion.cpp` + the other `test_batch_*` + `mopac_vs_ff14sb_recon`;
and #3 (HM↔published-bond-sum benchmark) and #4 (BS↔HM emit) each land their validation as an assert in the
same build, alongside the formal-tool positive controls ([[feedback_formal_tools_as_pipeline_references]]).
**These are CTest asserts on FIXTURES (the dev/CI gate) — never runtime asserts in the extractor.** A
regression blocks a *commit*, never a *fleet run*; the production `Compute`/`WriteFeatures` path stays
assert-free (emit only). The derived 0.36/1.0 band is deliberately *wide* → it fires on a real
sign/units/parity break, not on a protein's legitimately-odd geometry (so it won't flake). If a
production-side signal is ever wanted, it is a **logged `OperationLog` warning, never an abort.**

### Rule corrections
- **Rule 2 is FORWARD-looking, not current code.** The SDK/catalog still mark BS/HM/McConnell as *features*
  and `learn/` consumes them as ridge inputs. "Kernel = hypothesis, not feature" is the direction-4
  posture, not today's behavior — phrase it so.
- **Rule 4 has two legitimate exceptions.** Full `0e⊕1e⊕2e` holds for the shielding-style tensors
  (BS/HM/McConnell full arrays). **EFG is deliberately `2e`-only** (traceless-by-Laplace — no `0e/1e`; the
  E-field carries the `1o`), and **ring per-type arrays currently drop T1** (storage choice — candidate
  restore under rule 7). "Full irrep always" = for the even shielding tensors, not the EFG field.
- **Rule 6 — concrete trajectory-writer bugs:** BS/HM H5 time-series mark shielding parity `0e+1o+2e`, but
  the antisymmetric part of a shielding tensor is **`1e` (even), not `1o`** → should be `0e+1e+2e`; and the
  BS H5 units attr says `ppm` not `ppm_T_per_nA`. (The static-NPY catalog `parity="odd"`+`irreps="1e"` on
  E-fields is the separate, known item.)

### What the guide missed (add)
- **Ring scale/units consumer rule:** static BS catalog (per-nA) vs BS H5 (`ppm`) vs h5-reader viz
  (×literature intensity) vs learn-side (learned) — four places, one quantity; never conflate.
- **Ring-axis contract check before re-emit:** `TopologySidecar` says `ring_geometry`/`ring_contributions`
  are aromatic-ring-axis only, but `ConformationResult` writes `protein.RingCount()` rows and ring code has
  saturated-Pro at index 8 — confirm `RingCount()` is aromatic-only on a fixture before trusting the manifest.
- **SDK consumer-update item** for new EFG emits: `CoulombGroup`/`MopacCoulombGroup` lack sidechain-EFG /
  partitioned-E fields, and `EFGTensor` rejects old 9-wide EFGs — new emits need matching SDK entries.
- **Trajectory parity/units metadata pass** (not just the static catalog) — the BS/HM H5 attrs are the
  easiest place to get the right tensor with the wrong label.

## The evidence behind rule 1 — the T2-physical-tensor map (condensed)

- **Physical already (target-based, no scale surgery):** the EFG family (Coulomb / MOPAC / APBS / water /
  AIMNet2-EFG — genuine field-gradient tensors, traceless-projected); the DFT-resident tensors (ORCA
  shielding/dia/para, ProCS15 tripeptide, Larsen grids — the *targets/anchors*).
- **Geometric, source-based, cheaply upgradeable (the scale rides the frame):** **McConnell** (unit `Q̂` →
  literature `Δχ`, scattered), **ring** (unit-current `G` → literature ring intensity, solid/cited),
  π-quad (→ `Θ`, but caged). Tabulated first; Trp-cage QM as a reference check (verify its ORCA actually
  emits the moment before claiming it).
- **Feeders (no own `2e`):** MOPAC, AIMNet2 — charges / populations / Wiberg-BO / embedding feed the
  tensors above (charges → EFG, Wiberg-BO → McConnell strength); they do not carry their own physical T2.
