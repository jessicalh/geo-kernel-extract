# Applied-mathematics audit — field/tensor handling + variance decomposition

> **Historical audit record — not current truth (trued 2026-06-04).** Superseded
> by the June 4 maths walk/checklist/fix postmortems; preserve as provenance.

Branch `h5-reader-pysr-spike` (NEVER MERGE). Authored by the Opus applied-maths
auditor, 2026-06-01. **Audit + corrected designs only; no C++ touched, no fit
re-run, nothing committed by me.** This sits beside (and does not supersede)
`analysis/VARIANCE_DECOMPOSITION_METHOD.md` — that doc is the
between/within design + the solvation footnote and is *correct*; this doc adds
the field/tensor MATH audit it deferred, and cross-references it for §B/§C.

Calibration to the lead's priority: I hunted for places the **applied maths is
wrong or incomplete**, not setup caveats. Low-confidence one-protein prototype
signal is fine (`correlate-not-match`, report effective N). What must survive a
methods-critical committee is the maths. Verdicts ("a law fell out") are the
lead's; mine is "is the maths right."

Sandbox note: Bash here is restricted to read-only `ls`; Python execution was
denied. The change-of-basis orthogonality below is verified **by hand** from the
pinned literal (it is exact rational/sign structure), not by running numpy.

---

## SUMMARY TABLE (read this first)

| # | location | severity | finding |
|---|---|---|---|
| F1 | `BroadBackbone.cpp:93-95` `bondKernelT2FromSources` | **MEDIUM-real, must-fix-or-justify** | bond-anisotropy tensor is an ad-hoc asymmetric expression, NOT a derived McConnell point-dipole shielding tensor; the `9cosθ d̂b̂ᵀ` term has no clean provenance. Its **symmetric** part is what survives `DecomposeLibrary`, so it is *a* valid l=2 object, but it is not the textbook neighbour-anisotropy tensor and the de-circularising "literature kernel" claim built on it (`bond→C component_r=+0.76`) is using a home-rolled, not literature, tensor. |
| F2 | `BroadBackbone.cpp:464-465` | **MEDIUM** | the "bond literature kernel" is the **H5 McConnell kernel** (`kernel_mc`) where present, else the F1 home-rolled tensor — two different objects under one column name, chosen per-row by a silent fallback. Per-stratum branch split unverified; a mixed-branch stratum's bond number is two mechanisms averaged. |
| F3 | `buckingham_efield_t0.py` whole-file | **LOW (incompleteness, disclose)** | the l=1 field is collapsed to two l=0 invariants (`E_z`, `|E|²`); the proper field→**T1** (antisymmetric l=1) equivariant target is emitted-but-never-fit. Correct for σ_iso; the dropped T1 channel is an *incomplete* treatment of "the electric field," not a wrong one. |
| F4 | `equiv_t2_efg_e3nn.py:154-160` + `EFG_ARC_EVIDENCE.md` | **LOW (correct, label tightening)** | the Schur scalar-gate `g(\|EFG\|)·EFG_T2` IS the correct 2e→2e equivariant map. But `g` is a function of `\|EFG\|`, itself a function of the multiplied feature; with a near-constant fitted `g` the "law" is a single scalar rescale per stratum — report it as a fitted linear coefficient (the evidence doc already does), not as a recovered field-gradient response law. |
| F5 | `backbone_literature_kernel_t2_corr.py:208` | **MEDIUM (axis-label, already flagged)** | `component_r` is **global-mean** Pearson on flattened 5-vectors over a 54-atom × 500-frame panel = **total/between-dominated** axis. The `+0.70` numbers are NOT the dynamic axis. (Confirmed; this is the core of `VARIANCE_DECOMPOSITION_METHOD.md` §0 — I concur, full stop.) |
| F6 | `chargeKernelT2FromSources` `BroadBackbone.cpp:107-127` | **OK — verified correct** | the charge→EFG tensor `Σ q(3d̂d̂ᵀ−I)/r³ · ke` is the correct traceless point-charge field-gradient; sign-invariant in `d`, origin = target atom (so origin-independent as an EFG *at the nucleus*), explicit re-trace-removal is redundant-but-harmless. No bug. |
| F7 | `change_of_basis.py` `_C_FROZEN` | **OK — verified correct by hand** | C is exactly orthogonal (det +1, norm-preserving); the (zz, xx−yy) 2×2 block is a clean −120° rotation `[[−½,−√3/2],[√3/2,−½]]`, the other three rows are unit selectors. `\|T2\|` is preserved. The pinned constant is right. |
| F8 | `DecomposeLibrary` / `reconstructLibraryTensor` | **OK — verified correct** | isometric real-SH basis, Frobenius-preserving (`Σ\|T2_m\|² = Σ S_ij²` for traceless S — checked); reconstruct is the exact inverse incl. the T1 antisymmetric re-add. |
| V1 | all fitters | **the incompleteness** | uniform per-atom de-mean = within-only; the between (static-environment) component is discarded for every mechanism. Design fix in §B (concurring with + extending `VARIANCE_DECOMPOSITION_METHOD.md`). |

Net: the **change-of-basis, the EFG charge tensor, the spherical decomposition,
the frames, the Schur gate, and the Buckingham σ_iso projection are
mathematically sound.** The two real maths problems are **F1 (the bond
anisotropy tensor is home-rolled, not derived)** and **V1 (variance
decomposition is one-sided)**. F2/F3/F4/F5 are labelling/completeness fixes.

---

## A. FIELD / TENSOR HANDLING (lead with these)

### F1 — the bond-anisotropy T2 tensor is ad-hoc, not a derived McConnell tensor

`BroadBackbone.cpp:79-105`, per bond source (d̂ = bond-midpoint→target unit
vector, b̂ = unit bond axis, cosθ = d̂·b̂, r = |displacement|):

```cpp
total += (9.0 * cosTheta * dHat * bHat.transpose()      // 9 cosθ (d̂ ⊗ b̂)
          - 3.0 * bHat * bHat.transpose()                //  −3 (b̂ ⊗ b̂)
          - (3.0 * dHat * dHat.transpose() - I)) / r3;   //  −(3 d̂⊗d̂ − I)
```

then `DecomposeLibrary(total)` → T2.

**Why this is a maths problem.** The neighbour-anisotropy (McConnell) shielding
*tensor* from an axially-symmetric group is a well-defined object: it is the
field-gradient-type response of the induced point dipole `m ∝ Δχ (b̂·B) b̂` at
the nucleus, giving a shielding tensor whose **trace** reproduces the McConnell
scalar `σ_iso ∝ (3cos²θ−1)/r³` (the very `dipolar` the cell already emits at
`McConnellNeighborhood.cpp:166`). A correctly derived tensor must therefore:

1. be a sum of the rank-2 dyads `{d̂⊗d̂, b̂⊗b̂, d̂⊗b̂_sym, I}` with coefficients
   that are **functions of cosθ only** (the field of a unit dipole has fixed
   angular structure), and
2. have an **isotropic part consistent with `(3cos²θ−1)`** so the tensor and the
   scalar kernel agree.

The code's expression fails the provenance test on both counts:

- The `9 cosθ (d̂⊗b̂)` term is **not symmetric** (d̂ ≠ b̂ in general), so only its
  symmetrised part `(9cosθ/2)(d̂⊗b̂ + b̂⊗d̂)` survives `DecomposeLibrary`. That a
  bug is *masked* by symmetrisation is itself a smell — whoever wrote it did not
  intend the antisymmetric half, and the coefficient `9` has no derivation in any
  comment or doc.
- **Trace check (the decisive one).** `tr(d̂⊗b̂) = cosθ`, `tr(b̂⊗b̂)=1`,
  `tr(3d̂⊗d̂−I) = 0`. So `tr(total)·r³ = 9cosθ·cosθ − 3·1 − 0 = 9cos²θ − 3 =
  3(3cos²θ − 1)`. So the trace IS `3(3cos²θ−1)/r³` — proportional to the
  McConnell scalar. **Good: the angular trace is right** (off by the overall
  factor 3, which a fit absorbs). So the tensor is not garbage — its isotropic
  part matches the scalar kernel. But `DecomposeLibrary` then *removes* the trace
  (T2 is the traceless-symmetric part), so the **T2 the fitter actually uses is
  the traceless part of an expression whose only verified property is its
  trace.** The traceless-symmetric structure — the part that is the entire point
  of a T2 kernel — is unverified against any literature form.

**What the correct derivation gives.** The shielding tensor from an
axially-symmetric neighbour group (the standard McConnell/Buckingham
point-dipole anisotropy, e.g. as written in the bond-anisotropy literature) is

```
σ_neighbour(r, d̂, b̂)  ∝  Δχ_ax / r³ · [ (3 d̂⊗d̂ − I) projected through the
                                          dipole induced along b̂ ] .
```

Concretely, the induced dipole is `m = α·b̂ (b̂·B)`; the dipole field at the
nucleus is `B_ind ∝ (3 d̂(d̂·m) − m)/r³`; collecting the linear map `B ↦ σ·B`:

```
σ_ij  ∝  (1/r³) [ 3 (d̂·b̂)(d̂_i b̂_j)  −  b̂_i b̂_j ]                       (★)
```

i.e. the coefficients are **3 (not 9) on the d̂⊗b̂ term and −1 (not −3) on b̂⊗b̂,
and there is NO standalone `(3d̂⊗d̂ − I)` term** — that term is the *isotropic
neighbour* (Δχ=0 / spherically-symmetric) contribution, which for an
anisotropic-susceptibility model should not be added separately. The code's (★)
is *not* the standard form; it is a different linear combination. Whether (★) or
the code's expression is "right" depends on which susceptibility model is
intended (full χ tensor vs axially-symmetric Δχ), and **that choice is
undocumented** — which is the audit finding: an l=2 kernel asserted as the
"bond literature kernel" must trace to a written susceptibility model, not to an
un-annotated dyad sum with a magic `9`.

**Correct treatment / required action (codex brief item 1):**
- Decide the susceptibility model explicitly (axially-symmetric Δχ is the
  standard McConnell choice) and **derive the tensor symbolically** to the form
  (★) (or the full-χ form if that is intended), with the coefficients written
  down. Put the derivation in a comment + `BACKBONE_LAW_EVIDENCE.md`.
- Replace lines 93-95 with the derived **symmetric** dyad sum (no asymmetric
  d̂⊗b̂; use the symmetrised `½(d̂⊗b̂ + b̂⊗d̂)`).
- Keep the trace-removal (T2 is traceless by construction) but now the traceless
  part is a *derived* object, so the de-circularising "un-fitted literature
  coefficient" claim is honest.
- Until then: **the `bond_literature_kernel` column is NOT a uniformly-literature
  kernel** for the source-built rows; relabel the F1 path
  `bond_anisotropy_T2_homerolled` and do not present `bond→C +0.76` as a
  literature-coefficient-fixed de-circularising result without first confirming
  (per F2) which branch those rows took. Where the rows are the producer
  `kernel_mc`, the +0.76 is a legitimate producer-kernel-vs-DFT correlation under
  that label; where they are F1, it inherits this "home-rolled, not literature"
  caveat.

### F2 — "bond literature kernel" is two different objects under one name

`BroadBackbone.cpp:460-466`:

```cpp
agg.bond_literature_kernel = h5KernelT2Local(body, ArrayId::KernelMc, atom, h5_row, fr); // producer MC
if (!agg.bond_literature_kernel.present)
    agg.bond_literature_kernel = bondKernelT2FromSources(fr, sources);                    // F1 home-rolled
```

So `broad_backbone_aggregated_bond_literature_kernel_T2.npy` is the **producer's
McConnell H5 kernel** (`kernel_mc`, a dense per-atom H5 tensor —
`Catalog.cpp:70-72`) rotated into the local frame *wherever that kernel is
present for the atom*, and the **F1 home-rolled tensor** wherever it is not.
These are different physical objects (the producer kernel is the library's full
McConnell result; F1 is the ad-hoc dyad sum). Mixing them per-row under one
column makes a per-mechanism correlation un-interpretable unless every row in a
stratum took the **same** branch.

- Which branch a stratum takes depends on whether the producer computed
  `kernel_mc` for that atom class. The McConnell kernel is typically a
  reporter-nucleus quantity (the oracle cell `McConnellNeighborhood` only strata
  on `IsBackboneAmideHydrogen` = HN), so the **C/O/N strata may be partly or
  wholly on the F1 home-rolled tensor while HN is on the producer kernel** — I
  did NOT verify the per-stratum `present()` split (would require reading the
  NPY values, out of scope for this audit). **Flag for the lead:** confirm, per
  stratum, which branch the `bond` rows took before quoting `bond→C +0.76` /
  `bond→HN −0.65` as a single mechanism. If a stratum is mixed-branch, its bond
  number is two mechanisms averaged.
- **Action:** split into two clearly-named columns — `bond_kernel_producer_T2`
  (present only where `kernel_mc` exists) and `bond_anisotropy_T2_homerolled`
  (the F1 form, after it is corrected per F1) — never fall one back to the other
  silently. The fitter/corr script then reports each under its true label, and
  the per-stratum branch split becomes visible (one column NaN, the other not).

### F3 — Buckingham collapses l=1 to l=0; the field→T1 channel is dropped

`buckingham_efield_t0.py` fits `σ_iso ≈ A·E_proj + B·|E|²` with `E_proj =
E_local.z` (`BuckinghamEfield.cpp:198`).

**The maths is correct for σ_iso.** A genuine rotation invariant cannot be a bare
vector component — but `E_z` here is the field projected onto a **molecule-fixed
reference axis** (z = the chemistry-defined bond direction: N→H for HN, CA→HA for
HA, O→C for carbonyl O, N→CA for N — all from `LocalFrameBasis.cpp`, all built
from typed atoms). Projecting a vector onto a molecule-fixed axis IS an invariant
(it co-rotates), so `A·E_z + B·|E|²` is the standard Buckingham scalar form and
is frame-consistent. The sign of `E_z` is fixed per stratum by the bond-axis
sign convention, so `A` has a consistent meaning within a stratum. **No bug.**

**The incompleteness:** the electric field is an l=1 object. Its proper
equivariant target is the antisymmetric/vector part of the shielding response —
the **T1** channel — and the linear field→σ_iso map is only the l=1⊗l=1→l=0
(quadratic, the `|E|²` and `E_z²`-type) plus the trivial scalar. The codebase
*emits* `dft_total_T1` but labels it `unverified_emit_only`
(`BuckinghamEfieldSink.cpp:151`) and **never fits it**. So "we treated the
electric field" is, today, "we treated the field's effect on the isotropic
scalar only." That is a defensible scope choice (T1 convention is genuinely
unverified — antisymmetric shielding is gauge/convention-fraught), but it must be
**disclosed as a scope limit**, not presented as a complete field treatment.

- **Action (codex brief item 3, optional/disclosable):** either (a) state
  explicitly in `BACKBONE_LAW_EVIDENCE.md` that the field arc is σ_iso-only and
  the T1 channel is out of scope pending a T1 convention verification, or (b) if
  T1 is wanted, fit it as a proper l=1→l=1 equivariant map (an `o3` `1o→1o`
  scalar-gate, the l=1 analogue of the EFG Schur gate) — **not** by reading raw
  T1 components. Do NOT hand-roll. (a) is the cheap, honest default.

### F4 — EFG Schur gate is correct; tighten the claim to "fitted linear coeff"

`equiv_t2_efg_e3nn.py:154-160`: `pred = g(mag_norm) · feature` then per-atom
de-mean. For a linear equivariant map between two copies of the `2e` irrep,
**Schur's lemma forces a single scalar multiplier** — so `g·EFG_T2` with `g` an
invariant scalar is exactly the right (and only) linear 2e→2e form. ✓ The
nonlinear `g(|EFG|)` generalises it to an invariant-gated map, still equivariant.
✓ De-meaning and `|EFG|`-standardisation are fine.

The only tightening: `EFG_ARC_EVIDENCE.md` already shows the nonlinear gate adds
≤0.011 R² over a fitted constant `g`, i.e. the supported form is `σ_T2 ≈ γ_stratum
· EFG_T2`. Since `|EFG|` is computed from the same feature that is multiplied, a
flat `g` makes the model a single per-stratum scalar rescale. **Report it as a
fitted linear coefficient (the evidence doc does), and do NOT claim a recovered
nonlinear field-gradient *response law* — `γ` is fitted, not literature-fixed**
(the doc's de-circularisation section already concedes no fixed literature γ
exists for these strata; concur). No code change; this is a wording guard.

### F5 — the `component_r` axis label (concur with the variance doc)

`backbone_literature_kernel_t2_corr.py:208`: `pearson(kernel[mask].reshape(-1),
target[mask].reshape(-1))` centres on the **global** mean over a 54-atom × 500-
frame flattened panel. That is the **total** axis, between-dominated. The
`charge→N +0.696`, `bond→C +0.764` numbers in
`broad_backbone_literature_kernel_t2_correlation.csv` are total-axis, NOT the
de-meaned dynamic axis. This is exactly `VARIANCE_DECOMPOSITION_METHOD.md` §0 —
**I concur fully**; it is the central reason the decomposition in §B is needed.
(Also: `mean_row_cosine` and `magnitude_r` are likewise pooled/total. The
effective-N column IS axis-honest, computed on the within target-|T2|
autocorrelation — good, but it is reporting a within N_eff next to a total-axis
r, which is a mismatch the §B table fixes by splitting the axes.)

---

## B. THE COMPLETE VARIANCE DECOMPOSITION (V1)

`VARIANCE_DECOMPOSITION_METHOD.md` already designs this correctly (within/between
panel transforms, NOT a random-effects model; per-component effective N;
per-component normalisation; the §1.4 reporting table). **I reviewed that design
as an applied mathematician and concur** — the within/between split is the
textbook panel decomposition, RE is correctly rejected (atom effect *is*
correlated with the static field a mechanism feels, violating RE's orthogonality
assumption), and the §2 normalisation logic (scale-free for single-mechanism R²/r;
per-axis z-scoring only under the future multi-kernel ridge) is right.

I add three **maths corrections/sharpenings** to that design:

**B1 — the between-axis tensor estimator must use the per-atom MEAN tensor in a
FIXED frame, and report |T2| r, not component_r, as primary.** §1.2 says correlate
the per-atom-mean kernel-T2 against `T̄2_i`. Correct — but note the per-atom mean
is taken over frames *in the local frame that rotates with the molecule each
frame*. Averaging a rotating-frame tensor over frames is well-defined (the local
frame co-rotates, so the components are already in a body-fixed basis — this is
why the emit rotates to local frame at all), so `T̄2_i` is a legitimate
body-frame mean. **But** the between-atom `component_r` then compares
body-frame-component patterns across atoms whose body frames are differently
*defined* (an N frame vs a C frame are different conventions) — within one
stratum the frame convention is identical, so it is fine **per stratum**, and
must NOT be pooled across strata. The §1.4 table is per-stratum, so this holds;
just never compute a cross-stratum between `component_r`. Make `|T2| r` (frame-
free) the headline between number and `component_r` the secondary, per stratum.

**B2 — between N is the unit count, and the between R²/r needs a leave-one-atom-
out twin to be defensible, not a permutation null.** §1.3/§1.5 say report the unit
count and flag THIN; §1.5 makes LOAO opt-in. I sharpen: for the between axis the
**LOAO R² is the honest generalisation number** (does the mechanism order
*unseen* atoms?), and on a 54-atom stratum it is cheap and should be the reported
between metric, not opt-in — the raw between R² on 54 fitted points overstates
(it is in-sample across environments). Keep permutation tests off (agreed — over-
thoroughness). So: **between column = LOAO-R² (or LOAO-|T2|-r) + unit count**, not
in-sample between R².

**B3 — the within-axis effective N is correct as coded but state the variance it
deflates.** `effective_n` uses `n·(1−ρ₁)/(1+ρ₁)` on the per-atom target-|T2|
series (`backbone_literature_kernel_t2_corr.py:148-180`). Mathematically this is
the standard AR(1) effective-sample-size for the *mean*; it is the right first-
order deflation. The every-other-frame DFT sampling gives ρ₁ ≈ 0.02–0.23 (CSV),
so within N_eff stays near the raw count — **and the frame-split random hold-out
is valid here precisely because ρ₁ is small** (lag-1 ≈ 0.05–0.23 ⇒ adjacent
held-out frames are near-independent; the lead's reasoning checks out). Two
honest caveats to state: (i) AR(1) is a *model* of the autocorrelation — HA/HN at
ρ₁≈0.23 have non-trivial deflation and a single-lag estimate can under-deflate if
the series has slow components; report the median ρ₁ (already done) so a reader
can judge. (ii) N_eff is computed on the **target** |T2| series, which is the
right denominator for a correlation's CI, but the *predictor's* autocorrelation
also matters for the regression's effective DoF — for single-mechanism r it is
fine; flag it for the multi-kernel ridge.

**The corrected reporting surface (per mechanism × target-channel × stratum):**

```
mechanism, target{σ_iso|T2}, stratum,
  n_atoms,                       # between unit count (FLAG THIN <10)
  within_N_eff, median_lag1_rho, # within axis honesty
  BETWEEN:  loao_R2 [σ_iso]  OR  loao_|T2|_r + loao_component_r [T2]
  WITHIN:   frame_split_R2 [σ_iso]  OR  frame_split_|T2|_r + component_r [T2]
  treatment_label,               # F-section: "magnetic" | "FF-Coulomb-mismatch" | ...
  thin_flag
```

No composite score, no ranking (SETI not editor). The point is *placement*: which
axis each mechanism's signal lives on. The counter-example the brief cites
(de-circularising charge→N T2 is on the within axis while Buckingham σ_iso is
near-null there, yet Stage-1 found the field strong on the between axis) is
exactly what this table surfaces — **do not assume "electrostatics are static";
decompose and read it off.**

---

## C. SOLVATION CAVEAT (demoted to one paragraph — disclosable, not a work item)

The DFT target is r²SCAN/def2-SVP **CPCM(Water)**, protein-only (846 atoms, no
explicit waters; confirmed in `VARIANCE_DECOMPOSITION_METHOD.md` §3 against the
ORCA `.out`). Our electrostatic kernels are all *different* treatments of the same
electrostatics: FF14SB vacuum-Coulomb (`field_z`), APBS Poisson–Boltzmann of FF
charges (EFG), explicit-MD-water field (hydration). None is the CPCM reaction
field the DFT actually responded to, and that reaction field is **not in the
retained artifacts** (only aggregate CPCM energies/surface-charge totals were
written; no per-segment charges or potential-at-nucleus). **This is a disclosable
caveat: "a possible reason the electrostatic kernels underperform is that all
three are treatment-mismatched to the DFT's implicit-continuum solvation."** Carry
a one-line treatment label on every electrostatic row (F-table / §B). The magnetic
mechanisms (ring current, McConnell anisotropy) are solvation-insensitive and
escape this entirely. **Do NOT extract the CPCM field; do NOT propose a re-run.**
(I spent zero effort here beyond confirming the prior agent's read, per the brief.)

---

## D. CODEX IMPLEMENTATION BRIEF

Front-load `analysis/PATTERNS.md`: ONE typed C++ model; Python READS the emitted
CSV/NPY substrate and decomposes/fits; never recompute a kernel/projection/field,
never open `trajectory.h5`, reuse the frozen `change_of_basis.get_C()`. e3nn for
equivariance (no hand-rolled SH). correlate-not-match, per-element/per-atom-type,
report effective N, model-not-physics. Never-merge spike; commit on the branch.

**Item 1 (C++, MEDIUM — the real maths fix). Fix the bond-anisotropy tensor.**
- In `BroadBackbone.cpp:79-105`, derive and document the neighbour-anisotropy
  shielding tensor for the chosen susceptibility model (axially-symmetric Δχ is
  the standard McConnell choice). The derived form is a **symmetric** dyad sum;
  write the coefficients and their provenance in a comment + `BACKBONE_LAW_EVIDENCE.md`.
  Verify the **trace** reproduces `c·(3cos²θ−1)/r³` (it must, to agree with the
  emitted `dipolar` scalar) and that the d̂⊗b̂ term is symmetrised. Replace the
  magic `9`/`−3`/standalone `(3d̂d̂−I)` with the derived coefficients.
- This is a `broad_backbone`-only path; the ring/mc oracle byte-parity must hold
  (no other C++ touched). Re-extract `/tmp/rdc-broad-backbone-axes` after.

**Item 2 (C++, MEDIUM — relabel, prevent the conflation). Split the bond column.**
- In `BroadBackbone.cpp:460-466` + `BroadBackboneSink` schema: emit two distinct
  columns/NPYs — `bond_kernel_producer_T2` (present only where `KernelMc` exists)
  and `bond_anisotropy_T2_homerolled` (the item-1-corrected source tensor). No
  silent fallback. Update `backbone_literature_kernel_t2_corr.py` KERNELS map.

**Item 3 (Python, the COMPLETE variance decomposition — §B). New script
`analysis/variance_decomposition.py`.** Reads the emitted substrate only
(`broad_backbone_aggregated.csv` + the target/kernel T2 NPYs + the buckingham
`*_aggregated.csv` + `efg_aggregated.csv`). For each (mechanism × target-channel
× stratum):
- WITHIN: per-atom de-mean predictor + target on train frames (reuse
  `centred_by_train_atom`), frame-split R² (σ_iso) or frame-split |T2| r +
  component_r (T2). This reuses the existing fitters' within path verbatim.
- BETWEEN: per-atom MEAN (over frames, body/local frame for tensors) → regress
  across atoms; report **LOAO-R² / LOAO-|T2|-r** (B2), NOT in-sample between R²,
  + the unit count, FLAG THIN <10.
- Effective N: within = AR(1)-deflated frame count (reuse `effective_n`,
  `:148-180`); between = unit count. Report median lag-1 ρ (B3).
- One row per (mechanism, target, stratum) in the §B table; carry the
  treatment_label (magnetic / FF-Coulomb-mismatch / PB-mismatch / explicit-water-
  mismatch). NO composite score, NO ranking.
- Tensor work goes through the frozen `get_C()` only; never reconstruct the
  5-vector from a 3×3 in Python.
- Print the median lag-1 ρ + the every-other-frame note so the frame-split
  validity is stated, not assumed.

**Item 4 (Python, LOW — label tightening, no new maths).**
- `equiv_t2_efg_e3nn.py` / `EFG_ARC_EVIDENCE.md`: ensure the reported EFG result
  is worded as "fitted linear coefficient γ_stratum" not "recovered field-gradient
  response law" (F4). One-line edit to the evidence doc + the script's header.
- `BuckinghamEfieldSink` / `BACKBONE_LAW_EVIDENCE.md`: state the field arc is
  σ_iso-only and the T1 channel is out of scope pending T1-convention
  verification (F3). Do NOT fit raw T1.

**Discipline gate for the codex run:** grep the new/edited Python for forbidden
recompute (`3cos`, `/r**3`, `q*`, `lib_T2`, `spherical_harmonics` outside e3nn,
`h5py`/`trajectory.h5`) — must be clean. Print `|CᵀC − I|max` from `get_C()` at
startup. Re-run the ring/mc oracle parity after any C++ change (items 1-2) and
confirm byte-parity. Report effective N per row. Nothing merged.
