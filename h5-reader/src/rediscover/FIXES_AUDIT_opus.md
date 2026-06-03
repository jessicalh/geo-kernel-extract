# Independent Opus audit of codex commit f299a01 (applied-maths fixes)

Branch `h5-reader-pysr-spike` (NEVER MERGE). Authored by the independent Opus
auditor, 2026-06-01. **Audit only. No C++ touched, no fit re-run, nothing
committed.** This is the *independent* review of codex's fixes + self-audit
(`APPLIED_MATHS_AUDIT_codex.md`), cross-read against the original Opus audit
(`APPLIED_MATHS_AUDIT.md`) and the substrate. I formed my own view and traced
F1 to first principles + the project's own canonical sources.

Sandbox: Bash limited to read-only `ls`; Python execution denied (same as the
prior auditor). The F1 tensor algebra below is verified **by hand** at an
explicit geometry; the variance-share identity is verified **algebraically**
plus by self-consistency of the emitted output (shares sum to exactly 1.0 and
match codex's independently-reported values).

---

## SUMMARY (read first)

| fix | verdict |
|---|---|
| **F1 — bond-anisotropy T2 (`BroadBackbone.cpp:93-95`)** | **codex is RIGHT.** The expression is the project's own canonical McConnell tensor, derived from first principles in `GEOMETRIC_KERNEL_CATALOGUE.md:914-966` and implemented byte-identically in the library (`McConnellResult.cpp:90`). It is `σ ∝ K·χ_aniso` with K the dipole-field operator `(3d̂d̂ᵀ−I)/r³` and χ_aniso the **traceless axial susceptibility** `(3b̂b̂ᵀ−I)`. The `9` = 3×3 (provenance documented). The asymmetry is real, physical, and `DecomposeLibrary` correctly routes the antisymmetric half to T1; the T2 the fitter uses is the symmetric-traceless part of a *derived* tensor. Codex's decision to decline the change was correct. Opus-1's "home-rolled, magic 9, the standalone (3d̂d̂−I) is an isotropic term that shouldn't be added" reading is **wrong** on the physics. **One honest caveat survives** (below): it is the *traceless-χ* convention, which differs from the bare-dyad McConnell form by a documented deviatoric choice — both are literature-defensible, this one is internally consistent and is the project canon. |
| **variance_decomposition.py** | **Correct.** Between=LOAO (whole-atom holdout), within=train-only-centred frame split, AR(1) N_eff, SS-identity shares — all implemented correctly and leakage-free. `|T2| r` is the between headline (not cross-stratum component_r). Two disclosable caveats: T2 variance share is pooled-Frobenius (not per-component), and the tensor between-LOAO **R²** column is catastrophically negative (single-scalar-γ predicting a 5-vector mean) — fine because the headline is `absT2_r`, but the R² column will mislead a casual reader. |
| **q/r² → q/r³ fix** | **Confirmed.** Charge T2/EFG kernel is `q(3d̂d̂ᵀ/r⁵ − I/r³)`, whose radial gate ∝ `q/r³`. The corrected comparator `gate_vs_q_r^-3` is the right power; charge-N R²=0.981 is in the emitted `radial_fit_summary.csv` and is believable (the EFG tensor genuinely IS a q/r³ object). The old `q/r²` comparator is gone. |
| **EFG blocked split** | **Confirmed correct AND honestly small — but it fixes only ONE of two leaks.** Blocked split has genuinely 0 adjacent train/test pairs (`cross_split_lag1_pairs=0` emitted) vs 170 for random; O moved 0.316→0.333, honestly small. BUT the EFG feature+target are **lab-frame** (`EfgFeature.cpp:83`), and subtracting a *constant* lab-frame mean from a *rotating* tensor leaves a rotation-co-movement confound in BOTH old and new within numbers. The blocked split fixes temporal autocorrelation; it does NOT fix the lab-frame-rotation confound (codex's own finding #4). The small movement is honest *about the temporal leak* but does not speak to the rotation confound. |

Net: codex's fixes are sound. The headline correction is that **F1 resolves in
codex's favour** — the bond→C kernel rests on a correct, project-canonical
tensor, not a home-rolled one. The remaining items are correct with disclosable
caveats; the one substantive open issue is the EFG lab-frame-rotation confound,
which neither the random nor the blocked split addresses (it needs a local-frame
EFG emit, flagged but not done).

---

## F1 — THE CENTREPIECE. Definitive verdict: **codex is right.**

### The disagreement
The emitted tensor (`BroadBackbone.cpp:93-95`, per bond) is
```
M = [ 9 cosθ (d̂⊗b̂) − 3 (b̂⊗b̂) − (3 d̂⊗d̂ − I) ] / r³
```
Opus-1 called it home-rolled/non-standard (magic 9, the standalone `(3d̂d̂−I)`
shouldn't be added, the "real" McConnell tensor is `∝ 3cosθ(d̂⊗b̂)_sym − b̂⊗b̂`).
Codex declined to change it, calling it `D(d̂)·ΔChi(b̂)`, the standard dipolar
shielding construction, with trace = 3(3cos²θ−1) matching McConnell.

### First-principles derivation (the right answer has a derivation)
McConnell neighbour anisotropy: a group with susceptibility χ in field **B₀**
acquires induced moment **m** = χ·**B₀**; its point-dipole field at the target
nucleus at displacement r=r·d̂ is **B_ind** = (1/r³)(3d̂d̂ᵀ − I)·**m**. So the
shielding tensor (the linear map **B₀** ↦ −**B_ind**) is

```
σ = (1/r³) K · χ ,   K ≡ 3 d̂d̂ᵀ − I   (the dipole-field operator; symmetric, traceless)
```

For an axially-symmetric group along bond axis b̂, the susceptibility is
χ = χ̄ I + (Δχ/3)(3 b̂b̂ᵀ − I). The isotropic part χ̄ I gives χ̄·K, which is
traceless (no McConnell isotropic shift) and orientation-of-d̂-only — the
*isotropic-susceptibility* mechanism, which the McConnell derivation drops. The
**anisotropic** part is the McConnell neighbour-anisotropy tensor:

```
σ_aniso = (Δχ/3)/r³ · K · (3 b̂b̂ᵀ − I)
        = (Δχ/3)/r³ · (3d̂d̂ᵀ−I)(3b̂b̂ᵀ−I)
```

Expand using d̂d̂ᵀ b̂b̂ᵀ = (d̂·b̂) d̂b̂ᵀ = cosθ d̂b̂ᵀ:

```
(3d̂d̂ᵀ−I)(3b̂b̂ᵀ−I) = 9 cosθ d̂b̂ᵀ − 3 d̂d̂ᵀ − 3 b̂b̂ᵀ + I
                     = 9 cosθ d̂b̂ᵀ − 3 b̂b̂ᵀ − (3 d̂d̂ᵀ − I)
```

**This is the code, term for term.** The `9` is 3 (from K's `3d̂d̂`) × 3 (from
χ's `3b̂b̂`). The standalone `−(3d̂d̂ᵀ−I)` is NOT a separately-added isotropic
neighbour term — it is `K·(−I)` arising **intrinsically** from the `−I` in the
*traceless anisotropic* χ. Opus-1's "that term is the isotropic contribution
which shouldn't be added" is the error: the isotropic-susceptibility term that
IS dropped is χ̄·K, and it is absent here.

### Provenance — it is NOT home-rolled
- `GEOMETRIC_KERNEL_CATALOGUE.md:914-966` ("RESOLVED: McConnell Full Shielding
  Tensor Derivation") gives exactly the derivation above, including the trace
  verification `Tr(σ) ∝ Δχ(3cos²θ−1)/r³` and the explicit note that the tensor
  is asymmetric (T1≠0 is physical).
- The library `src/McConnellResult.h:12-14` and `src/McConnellResult.cpp:84-95`
  implement the **identical** `9 cosθ d̂_a b̂_b − 3 b̂_a b̂_b − (3 d̂_a d̂_b − δ)`.
  The rediscover `bondKernelT2FromSources` is a faithful port of the project's
  own canonical McConnell tensor — same formula the producer's `kernel_mc` H5
  oracle computes. The `9` is documented, not magic.

### By-hand T2 check (geometry: b̂=ẑ, d̂=(sin60°,0,cos60°), cosθ=0.5)
Symmetric-traceless (T2) part of the **code/codex form** M_code:
`Sxx=−1.0, Syy=+1.25, Szz=−0.25, Sxz=+0.6495, Sxy=Syz=0`; trace(M)=−0.75 =
3·(3·0.25−1) ✓ (McConnell scalar × 3).
The **bare-dyad** form Opus-1 proposed (`3cosθ d̂b̂ᵀ − b̂b̂ᵀ`) gives
`Sxx=+0.083, Syy=+0.083, Szz=−0.167, Sxz=+0.6495`. These are **not** scalar
multiples — they differ by `K·(I/3)·Δχ`, i.e. the deviatoric-vs-bare-dyad
**convention** for χ. Algebraically `M_code = 3·[K·(b̂b̂ᵀ − I/3)] = 3·K·Δχ_traceless`
— precisely the paramagpy/PCS axial form (Δχ defined traceless). So codex's
literature citation (paramagpy, McConnell review) is apt; the project's choice
is the traceless-χ convention.

### (a)/(b)/(c) the brief's three sub-questions
- **(a) legitimate physical model or ad-hoc product?** Legitimate. `K·χ` is the
  derived dipolar shielding map; not an ad-hoc product. The matrix-product
  framing codex used is just the factored form of the derivation.
- **(b) is the `9cosθ d̂b̂ᵀ` asymmetry real, and does Decompose keep only the
  symmetric half?** Yes and yes — and that is **correct**, not a masked bug.
  `K·χ` (product of two symmetric matrices) is generally asymmetric; the
  catalogue states T1≠0 is physical. `DecomposeLibrary` (`SphericalBasis.cpp:15-25`)
  routes the antisymmetric half to T1 (`0.5(s_ij−s_ji)`) and keeps only the
  symmetric half in T2 (`0.5(s_ij+s_ji)`). The fitter's T2 is the symmetric-
  traceless part of the *derived* tensor — the correct observable shielding T2.
- **(c) does the emitted traceless-symmetric part equal the standard form?** It
  equals the standard form **under the traceless-χ (PCS) convention** and
  differs from the bare-dyad-χ McConnell form by the documented deviatoric term.
  Both are literature-defensible; the code's is the project canon and is
  internally consistent (matches the producer `kernel_mc`).

### Verdict
**Codex is right; the bond→C headline rests on a correct, project-canonical,
first-principles-derived kernel.** Do not change `BroadBackbone.cpp:93-95`.
F2 (Opus-1's "two objects under one name") is also substantially weakened: the
`kernel_mc` H5 branch and the `bondKernelT2FromSources` fallback compute the
**same tensor formula**; they can differ only in cutoff/source set, not in the
physical model. Worth a one-line note that the two branches share the formula.

**The one honest caveat to disclose** (not a fix): the kernel uses the
traceless-χ convention (implied χ_∥=2, χ_⊥=−1, Δχ=3 in the bare units before the
fitted Δχ/3 prefactor). This is a convention choice, documented in the
catalogue, and it differs from the bare-dyad McConnell tensor by an intrinsic
deviatoric term. A methods-critical reader may ask "why traceless χ and not the
bare anisotropy dyad?" — the answer is the standard separation of the isotropic
susceptibility, and the project should cite the catalogue derivation when this
kernel is presented. This is a labelling/citation point, not a maths error.

---

## variance_decomposition.py — correct, with two disclosable caveats

**Between = LOAO (whole-atom holdout): correct.** `scalar_between_loao`
(`:241-253`) and `tensor_between_loao` (`:285-296`) fit on per-atom MEAN
predictors/targets across all OTHER atoms and predict the held-out atom — the
atom is never in its own training fit. The unit is the atom, the right between
unit. ✓

**Within = train-only-centred frame split: correct and leakage-free.**
`centred_by_train_atom` (`:155-166`) subtracts the **train-frame** per-atom mean
(`mt = m & train`) from all rows; the fit uses train rows, scoring uses test
rows. Test-frame information never enters the centering — this is exactly the
leakage-free within fix codex's own audit finding #7 demanded, and it is a real
improvement over `equiv_t2_backbone_e3nn.py` (which de-meaned over all frames).
✓

**AR(1) N_eff: correct.** `effective_n_components` (`:201-234`) computes
`n·(1−ρ₁)/(1+ρ₁)` per atom per component on the de-meaned target series, clamped
to [1,n], summed over atoms. Standard AR(1) effective-sample-size. ✓

**Variance shares: correct (verified algebraically).** `variance_shares`
(`:183-198`) computes SS_between = Σ n_i(ȳ_i−μ)², SS_within = Σ(y_it−ȳ_i)²,
total = Σ(y_it−μ)², all over the **same** masked rows, so the ANOVA identity
SS_total=SS_between+SS_within holds exactly (the cross term vanishes since
Σ_t(y_it−ȳ_i)=0 per atom). Shares therefore sum to exactly 1.0.

**My spot-check** (Python denied, so algebraic + self-consistency):
- Emitted `variance_decomposition.md` shows every row's
  `between_share + within_share = 1.0000` (N T2 0.4647+0.5353; C T2 0.333+0.667;
  C sigma_iso 0.0811+0.9189; N sigma_iso 0.7707+0.2293). The identity is
  preserved row-by-row.
- These match the codex self-audit's independently-computed shares (N T2
  0.465/0.535, C T2 0.333/0.667, O T2 0.343/0.657; N sigma_iso 0.771/0.229,
  C sigma_iso 0.081/0.919). Two independent derivations agree → the
  implementation is computing the textbook decomposition correctly.
- Read-off that the design intended (and the numbers deliver): bond/charge T2
  on **C/O** are within-led (within share ~0.66), Buckingham σ_iso on **N** is
  between-led (0.771) while on **C** it is within-led (0.919) — i.e. "the field
  is static" is FALSE per-stratum, exactly the placement the decomposition was
  built to surface. Honest.

**`|T2| r` is the between headline, not cross-stratum component_r: confirmed.**
The markdown summary (`write_markdown` `:506-514`) shows `between_LOAO_absT2_r`
and omits `between_LOAO_component_r`; the print loop (`:575-579`) uses `absT2_r`
for the T2 between number. `component_r` is computed per-stratum into the CSV
but never pooled across strata. This honours Opus-1's B1/F5 (frame-free |T2| as
headline, no cross-stratum component_r). ✓

**Caveat 1 (disclosable): the tensor between-LOAO R² column is junk.**
`tensor_between_loao` predicts a held-out atom's mean 5-vector with a single
scalar `gamma·x̄_i` (`tensor_gamma`, `:269-273`). Per-atom mean tensors point in
different directions, not just different magnitudes, so a scalar rescale gives
R² = −39.5, −95, −108 (emitted). These are not a defect in the metric *as used*
(the headline is `absT2_r`, which is sane: N broad_total +0.771, N charge
+0.776), but the negative `between_LOAO_R2` values sit in the CSV and will
mislead a casual reader. Recommend: either drop the tensor between-R² column or
annotate it "scalar-γ rescale; use absT2_r for tensor between".

**Caveat 2 (disclosable): T2 variance share is pooled-Frobenius, not
per-component.** `variance_shares` sums squared deviations across all 5 T2
components into one scalar (`.sum()` at `:194-195`). That is the natural
Frobenius between/within share and is fine for *placement*, but it is a single
number, not the per-component split B1 also mentioned. Disclose; not wrong.

No leakage found. No recompute (the file reads emitted CSV/NPY, maps tensors
through frozen `get_C()`, prints `|CᵀC−I|max` at startup `:536-537`). Discipline
clean.

---

## q/r² → q/r³ fix (`backbone_distill_evidence.py:329`) — CONFIRMED

The emitted fixed charge T2 kernel is the EFG-like traceless tensor
`q(3 d̂d̂ᵀ/r⁵ − I/r³)·k_e` (`BroadBackbone.cpp:107-123`). Factoring r out:
`q(3d̂d̂ᵀ − I)/r³`, i.e. a unit-`Y2(d̂)` angular object times a radial gate
**∝ q/r³**. So the radial comparator for the learned charge gate must be `q·r⁻³`,
not the field's `q·r⁻²`. The fix at `:329` (`x = q·np.power(r,-3.0)`, fit
`gate_vs_q_r^-3`) is the **right power**; the old `q/r²` comparator is gone.

Believable, not an artifact: `radial_fit_summary.csv` shows
`N,charge,gate_vs_q_r^-3,−3.0,0.9814,30000` — charge-N R²=0.981. The charge
kernel genuinely IS a q/r³ EFG object, so a high R² against the matched radial
power is expected, not surprising; n=30000 sampled readouts. The companion
log-fit (`log_abs_gate_over_q_vs_log_r`) recovers power −4.03 for N — steeper
than −3, consistent with the learned gate carrying extra (angular/category)
structure beyond pure radial, which the linear `q·r⁻³` comparator still captures
at 0.981 because r⁻³ dominates the variance. The BACKBONE_LAW_EVIDENCE.md prose
(`:71`, "q·r⁻² R²=0.920") is now **stale** relative to the corrected run — flag
to refresh that doc to the q·r⁻³ / 0.981 number so the evidence doc and the
emitted CSV agree.

---

## EFG blocked split (`equiv_t2_efg_e3nn.py:43`) — CONFIRMED correct, but fixes
only one of two leaks

**The blocked split is correctly implemented and genuinely decorrelated.**
`make_split_masks` (`:204-211`) takes a contiguous test block, purges
`purge_frames` neighbours on each side from train, and `cross_split_lag1_pairs`
(`:220-225`) counts adjacent train/test frame pairs. Emitted output:
blocked rows show `cross_split_lag1_pairs = 0` (all 6 strata); random shows 170.
So with one-frame purge there are genuinely 0 adjacent train/test pairs. ✓
No other leakage in the EFG fitter: `mag_norm` standardisation uses train-only
mean/std (`:263-266`), `init_gate` uses train-only sums (`:268-272`), centering
is train-only (`:259-260`), forward de-mean uses `center_mask=train_mask`
(`:319`). Clean. ✓

**The small movement is honest about the temporal leak.** O within R² 0.316
(old random) → 0.333 (blocked); the other strata barely move. The lab-frame EFG
T2 lag-1 ρ ≈ 0.86 (emitted `median_lag1_rho` 0.80–0.86), so adjacent frames are
highly autocorrelated, yet the random-split R² was only marginally inflated —
because the signal is weak and the gate is a single scalar, the leak had little
to inflate. Honest.

**But a SECOND, unaddressed leak/confound remains** (codex's own finding #4,
NOT fixed by the blocked split). `EfgFeature.cpp:83` emits feature and target in
the **lab frame** (`LocalFrame labFrame; // invalid: BuildTarget emits lab-frame
T2 only`; manifest `t2_components: FRAME-ALIGNED` means feature/target are
mutually aligned per frame, still lab-frame). The within de-mean subtracts a
*constant* lab-frame per-atom mean from a tensor whose physical baseline
**rotates** with the molecule frame to frame. A constant subtraction does not
remove a rotating baseline; both feature and target retain rotational residual,
and they **co-rotate** (same neighbourhood, same lab frame), so `γ·feature` can
track the common rotation rather than dynamic physics. The blocked split removes
*temporal* coupling between train and test; it does NOT remove this lab-frame-
rotation confound. The correct fix — emit `efg_feature_local_T2` /
`efg_target_local_T2` in the shared backbone frame + a `local_frame_valid` flag
— is flagged in codex's finding #4 but was **not implemented**; only the split
strategy changed. So: the EFG within numbers (old and new) should be read as
lab-frame, rotation-confounded; the small old→new movement speaks only to the
(small) temporal leak, not to the confound. Disclose this; it is the one
substantive open issue among the fixes.

---

## Anything else the fixes got wrong / loose ends

1. **Stale evidence docs vs corrected runs.** `BACKBONE_LAW_EVIDENCE.md:71`
   still quotes the pre-fix `q·r⁻²` charge-N R²=0.920; the corrected substrate
   gives `q·r⁻³` R²=0.981. `EFG_ARC_EVIDENCE.md` reports the EFG numbers without
   stating the split is now blocked (the doc's "same deterministic frame split
   as the capture fitter" predates the blocked default). Both should be
   refreshed so docs match the emitted CSVs. (Documentation, not maths.)

2. **`bond_anisotropy_T2` in variance_decomposition is labelled
   "magnetic, solvation-insensitive"** (`variance_decomposition.py:419-421`).
   Correct given the F1 verdict — the McConnell tensor is a magnetic-
   susceptibility mechanism, solvation-clean — so this label is now *justified*
   (it would have been questionable under Opus-1's "home-rolled" reading).

3. **EFG `apbs_efg_T2_old_random` vs `apbs_efg_T2` both emitted.** Good practice
   — the script keeps the old random-split row beside the blocked row so the
   leak size is visible (170 vs 0 cross pairs side by side). Keep it.

4. **F1 not reflected in the doc surface.** The catalogue derivation
   (`GEOMETRIC_KERNEL_CATALOGUE.md:914`) is in the library tree, not the
   rediscover analysis docs. When the bond→C result is written up, cite that
   derivation explicitly so the "literature kernel" label is defensible on its
   own (the traceless-χ convention is the point a committee will probe).

---

## Discipline check on this audit
- No C++ edited, no fit re-run, nothing committed (this doc only; flagged).
- Read emitted substrate (`/tmp/rdc-broad-backbone-axes`, `/tmp/rdc-efg`,
  `/tmp/rediscover-variance-decomposition`, `/tmp/rdc-backbone-law-evidence-qfix`)
  and source only. No `trajectory.h5`. No recompute into a model path (the F1
  by-hand check is a one-off pinned-geometry algebra verification, the labelled-
  fixture exception, not reusable recompute code).
- F1 traced to first principles + `GEOMETRIC_KERNEL_CATALOGUE.md:914-966` +
  `src/McConnellResult.{h,cpp}` (read-only, for provenance).
