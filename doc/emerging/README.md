# doc/emerging/ — "Napkin Rescue"

The consolidated home for the **forward build**: the new calculators, the revamped equivariant model,
and the revamped stats — grounding the four kernels from back-of-the-napkin to defensible, *cited*
specs. **Keep this index fresh** — it's the map; update it when docs land. (Scope boundary + the three
NEW points + the bar: `CONTROLLING_SPEC.md`.)

## The spine

- **`CONTROLLING_SPEC.md`** — the controlling document. Scope; the three NEW points (law study ·
  tensor predictor · shift predictor); the Four (Biot–Savart · Haigh–Mallion · McConnell · charge/EFG);
  the feature tiers; the feature × value matrix; and **the bar** (defensibility, not agent-confidence —
  no going back to the napkin). *Now carries (2026-06-07) the **dragon** — the solution shift IS the
  `0e`, the `2e` is orthogonal to it — and the **fido reframe**: equivariance is Part 1's analysis lens,
  earned-only-where-it-helps in the predictors; both predictors are FIDO.*

## The narrative / account layer

- **`GEOMETRIC_KERNEL_ACCOUNT.md`** — what the geometric kernels literally are (low-order multipole moments).
- **`WHAT_ARE_WE_DOING.md`** — the open description of the pipeline through the model.
- **`STATS_PROGRAM.md`** — the stats guardrail (say-sensible-at-the-right-confidence; don't load the dice).
- **`DEFERRED_LEDGER.md`** — the living ledger of everything **deliberately deferred / "not ours" / parked / couponed** (what · why · owner · where). Every spent coupon lands here; nothing important lives only in chat.
- **`STAIRCASE_SOCIAL_HISTORY.md`** — the social history ("small committee in a basement").
- **`GEOMETRIC_KERNEL_MATH_LINEAGE.md`** — the kernel as an ℓ=2 steerable-kernel family member.
- **`equivariant_methodology_and_field_survey.md`** — the model's recognised-ground grounding: the cited
  e3nn methodology (linear algebra to modest depth) **+ the honest field survey** — every *piece* is
  recognised ground, the *whole combination* is genuinely unusual; a defensible MSc novelty if framed as
  integration with honest limits.

## `kernel_design/` — each kernel's design + its grounding, side by side

How many docs a kernel has shows its work-through state:

| kernel | design | first-stage grounding | structured grounding | shareable spec |
|---|---|---|---|---|
| **McConnell** | `mcconnell.md` | `mcconnell_grounding_agent1.md` | `mcconnell_structured_grounding.md` ✓ | **`mcconnell_spec.md` ✓** + **`mcconnell_integration_addendum.md` ✓** (all 5 decisions resolved: strict-backbone · ring-always-on→unconditional-zero · H-valency→ablate · metadata→JSON-manifest · 1o→1e contained). **BUILT + adversarial PASSED** (sound/contract-faithful; physics verified, PCS identity 1e-15, 6 checks real, golden untouched). Pending: deliberate re-bless + 2 tidy-ups (DRY manifest, layout-string) |
| **Biot–Savart + Haigh–Mallion** | `ring.md` | `bs_hm_grounding_agent1.md` (+ lit-sweep corroboration) | **`bs_hm_structured_grounding.md` ✓** | **`ring_spec.md` ✓ (draft — vet pending)** |
| **charge / EFG** | `charge_efg.md` | `efg_grounding_agent1.md` | reality-check ✓ (`efg_reality_check.md`, front-loaded) | **`efg_spec.md` ✓ (draft — vetting)** |

Supporting docs in `kernel_design/`: `pi_quadrupole.md`, `hbond_larsen.md`, `dispersion.md`,
`bond_anisotropy.md` (folds into McConnell), `mopac_extension.md` (the MOPAC capture spec),
`pipeline_adaptation.md`, `PIPELINE_SPINE.md` (the deep-vision architecture pass), `efg_reality_check.md` (the EFG code reality-check — **front-loaded** 2026-06-07, verified facts for the structured grounding), and **`CONTINUITY.md`**
(the live handoff / session-end state — **READ FIRST** when resuming).

## Work-through order

**McConnell (review + spec) → EFG → BS/HM.** The first-stage groundings are pre-staged so each
work-through starts warm. The bar holds throughout: a kernel is built only when its choices are
defensible and cited.
