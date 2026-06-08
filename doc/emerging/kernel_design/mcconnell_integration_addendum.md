# McConnell — Integration Addendum (the build contract)

*Drafted 2026-06-07 from the codex spec↔code reality-check. The physics is in `mcconnell_spec.md` +
`mcconnell_structured_grounding.md`; this is the **integration contract** that makes the spec buildable
**as written** against the current tree, so the coder doesn't discover these gaps mid-flight. Mechanical
parts are settled here; **DECISIONS marked for Jessica** are flagged and must be resolved before the
build fires.*

## 1. The output stems (concrete — settled)

Replace the old `mc_shielding` / `mc_category_T2` / `mc_scalars` emit. New emit: **one packed `(N,9)`
array per (category, channel)**, parity-even `[T0, T1[0..2], T2[0..4]]` per `SphericalTensor::PackFull9`:

```
mc_peptide_co_fixed.npy      mc_peptide_co_bo.npy
mc_peptide_cn_fixed.npy      mc_peptide_cn_bo.npy
mc_backbone_other_fixed.npy  mc_backbone_other_bo.npy
mc_sidechain_co_fixed.npy    mc_sidechain_co_bo.npy
mc_sidechain_other_fixed.npy mc_sidechain_other_bo.npy
mc_disulfide_fixed.npy       mc_disulfide_bo.npy
mc_aromatic_zeroed_fixed.npy mc_aromatic_zeroed_bo.npy        (zeros while ring active; kept for audit)
```
14 arrays, each `(N,9)`, units Å⁻³, irrep layout `0e,1e_x,1e_y,1e_z,2e_m-2..+2`.

## 2. The plumbing checklist (mechanical — settled)

- **C++ storage:** 7 categories × {fixed, bo} accumulators of `SphericalTensor` (or `Mat3` → decompose at
  emit). The current old aggregate buckets in `McConnellResult` are replaced.
- **Emit:** add `WriteFeatures` entries for the 14 stems via `ConformationResult::WriteAllFeatures` →
  `FrameNpyEmitter` → `NpyWriter` (which already writes arbitrary `rows×cols`).
- **SDK catalog:** one `ArraySpec` per new array in `python/nmr_extract/_catalog.py`; remove/deprecate the
  old McConnell array specs.
- **SDK loader/group:** `_protein.py` builds a new `McConnellGroup` exposing the 14 arrays as first-class
  outputs; `_tensors.py` wrapper exposes the irrep view.
- **Deprecate-and-add:** keep the old McConnell trajectory/time-series results compiling for the existing
  pipeline; the new emit is the forward surface.

## 3. The `1o → 1e` convention migration (correctness — settled fix, scope flagged)

The C++ `Types.cpp::Decompose` is correct (`1e`). Stale `1o` labels to fix:
- `python/nmr_extract/_tensors.py:39` (wrapper declares `1o`).
- H5 attrs: `McConnellShieldingTimeSeriesTrajectoryResult.cpp:101`,
  `MopacMcConnellShieldingTimeSeriesTrajectoryResult.cpp:143` (`0e+1o+2e`).

Fix = relabel to `0e+1e+2e` (even). The math doesn't change; only the declared parity. *(Same recurring
family as the BS/HM trajectory headers and `mcconnell.md` — already corrected there.)*

> **RESOLVED 2026-06-07 (refined — debt caught, risk contained; Jessica nervous about a 20-file producer
> sweep, rightly).** Even though it's a *label* (the math in `Decompose` is correct), changing the declared
> irrep strings could trip blessed/golden output and needs a per-file parity call. So:
> - **McConnell's OWN labels fixed inside the McConnell build** (in scope, reviewed there):
>   `McConnellShielding…`, `MopacMcConnellShielding…` H5 attrs, `_tensors.py:39` → `0e+1e+2e`.
> - **The rest = a dedicated, NAMED `1o`-parity-correctness pass** — caught now as a committed task (NOT
>   scattered to the work-throughs, NOT deferred-and-forgotten), done carefully: per-file judgment
>   (shielding-tensor → `1e`; genuine polar field — Coulomb/EFG E-field, AIMNet2 charge-response-gradient,
>   water field — KEEP `1o`), a golden-test impact check, and review. Approach gets Jessica's vet before
>   it fires.

## 4. Metadata home — **DISCUSS-FIRST (Jessica's declarative-format rule)**

The spec wants per-array metadata: `source_model = "unit susceptibility shape; scale pinned"`,
`bo_source = "MOPAC Wiberg bond order"`, `aromatic_zeroed_when_ring_active = true`, `irrep_layout`,
`units`. **`ArraySpec` has no fields for most of these today.** This is a new declarative surface, so per
your rule it is **discussed before the parser is written.** Options:
- **(a)** extend `ArraySpec` with an optional `metadata: dict` field (lightweight; lives in `_catalog.py`);
- **(b)** a small sidecar manifest beside the arrays — the **`TopologySidecar` / MOPAC
  `extraction_manifest.json` pattern** (NPY bulk + JSON manifest), which we already use elsewhere;
- **(c)** H5 attrs only (for the trajectory path) — but the static NPY path then has no metadata home.

> **RESOLVED 2026-06-07 — (b) the NPY-bulk + JSON-manifest pattern** (TopologySidecar / MOPAC
> `extraction_manifest.json`): bulk in the NPYs, a JSON manifest beside them carrying `source_model`,
> `bo_source`, `aromatic_zeroed_when_ring_active`, `irrep_layout`, `units`. One blessed seam, no new
> declarative format invented.

## 5. Ring-current-active policy — **DECISION**

The aromatic-zero needs a `ring_current_active` signal; McConnell has no such flag today. The runner
attaches Biot–Savart + Haigh–Mallion *before* McConnell (`OperationRunner.cpp:163`), so it's inferable.

> **RESOLVED 2026-06-07 — zero aromatic UNCONDITIONALLY, no flag.** `OperationRunner` attaches
> Biot–Savart + Haigh–Mallion **unconditionally** (no `--no-ring` toggle; core calculators run every
> production extraction, logging "no rings — nothing to compute" only when a protein has no aromatics).
> Ring-current is therefore *always* active → aromatic McConnell is *always* zeroed. No runtime flag /
> presence-check (matches no-gate-for-unconditional-sources). Document: "BS/HM always compute the aromatic
> ring-current; McConnell zeros aromatic to avoid the double-count." (No-aromatic-rings edge case → no
> aromatic bonds to zero anyway.)

## 6. Category semantics — **DECISION** (+ one settled)

Topology already classifies backbone / sidechain-co / sidechain-other / disulfide / aromatic
(`CovalentTopology.cpp`), and the sidecar writes `bond_category` + `is_aromatic`. Settled: **add
`Disulfide` to the McConnell category switch** (it's classified in topology but not handled in McConnell's
switch today).

> **RESOLVED 2026-06-07 — STRICT peptide backbone (N, Cα, C, O); CB → sidechain.** The settled
> protein-NMR / structural-biology convention (CB is the first sidechain atom). Apply it; where the
> topology helper (`CovalentTopology.cpp:101`) treats CB as backbone-ish, override to the strict
> convention for McConnell's categories.

## 7. The inherited open item — **JESSICA'S PHYSICS CALL** (separate)

`ALLBONDS` folds C–H/N–H/O–H into the bond-order/valency input (valency changed ~1116/1231 atoms).
**RESOLVED 2026-06-07 — don't fudge, MEASURE.** Build McConnell with the X–H sources **separable** (an
ablation scaffold, not a permanent hedge), then ablate **with-vs-without X–H against the DFT target**
(1P9J / 720) in Step 1; the DFT-match decides **A (keep)** or **B (drop)**, defended as a measurement, not
a physics fudge. **Build decision = keep X–H separable so the ablation is runnable;** the physics answer
is the Step-1 ablation result (collapse to A/B with data). (Lean on the likely outcome: small Δχ → probably
B — but the data rules.)

---

**Once §3–7 decisions land, the spec is buildable as written.** The coding brief then references
`mcconnell_spec.md` + `mcconnell_structured_grounding.md` + this addendum, and the core build is the
physics rewrite (old `M_ab` accumulation → the new `D(r)·Q̂_s` two-channel per-category source model) plus
this integration contract.
