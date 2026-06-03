> PARTIALLY SUPERSEDED (2026-06-02): see `../STATE.md` top (OPEN: weak DYNAMIC field σ-response). §3's "solvation treatment mismatch" (CPCM) framing is a defensible BETWEEN/absolute footnote ONLY — it must NOT be used to excuse the weak WITHIN-axis field σ-response (R²≈0.10). The within-fluctuation is the protein's own charge motion (FF14SB/APBS capture it); CPCM is a geometry-set continuum reaction field, not a within-fluctuating source. Per STATE, the within gap is decided by the raw ΔField-vs-Δσ correlation (units/prefactor vs projection-axis vs genuinely-weak-driver), not by "treatment mismatch." The variance-axis method itself (between=LOAO, within=train-only-centred split, AR(1) N_eff) is current and correct.

# Between/within variance decomposition — method design

Branch `h5-reader-pysr-spike` (NEVER MERGE). Authored by the Opus methods agent,
2026-06-01. **Design only; not implemented, nothing committed.** The "did a law
fall out" verdict is reserved for the lead + Jessica; this doc fixes the
*defensible method* and the *re-interpretation framework*, calibrated to: survive
a stats-critical reader's due diligence, not lookit-me over-thoroughness. One
protein, ~50 coupled atoms/stratum — low-confidence prototype signal is fine and
expected (`correlate-not-match`); what must be airtight is the axis handling and
the effective-N honesty.

Discipline (front-loaded, also in the codex brief): ONE typed C++ model, no
second model in Python. The analysis **reads the emitted CSV/NPY substrate only**
— never recomputes a kernel/projection/field, never opens `trajectory.h5`,
reuses the frozen `change_of_basis.get_C()` for any tensor work. See
`analysis/PATTERNS.md`.

---

## 0. The confirmed methods/category error (verified against the code)

The data is panel / repeated-measures: rows are **(atom × frame)**. Total
variance of any quantity over the rows of one stratum splits exactly into

```
Var_total(x) = Var_between(x_atom_means) + Var_within(x - x_atom_mean)
```

The between part is the **static environment** axis (how atoms differ from each
other, averaged over the trajectory — the same axis as a Stage-1 WT-ALA
between-structure delta). The within part is the **dynamic** axis (frame-to-frame
motion of one atom about its own mean). A mechanism can carry signal on either,
both, or neither, and **the existing fitters do not agree on which axis they
report** — that is the error, confirmed in three files:

- `buckingham_efield_t0.py:119-130` (`centred_by_train_atom`) and
  `equiv_t2_backbone_e3nn.py:329-331` / `:280-284` (per-atom de-mean of target
  AND prediction) → **WITHIN-only**. The per-atom baseline tensor / scalar is
  subtracted before fitting. Buckingham σ_iso comes out near-null here (the
  doc's "best HN R²≈0.10"); that is a statement about the *dynamic* axis only.
- `backbone_literature_kernel_t2_corr.py:208` + `pearson()` `:131-132` →
  **POOLED (total)**. `component_r` is Pearson on `kernel[mask].reshape(-1)` vs
  `target[mask].reshape(-1)`, centered by the **global** mean, NOT per atom. So
  the de-circularising numbers in `broad_backbone_literature_kernel_t2_correlation.csv`
  (`charge→N component_r=+0.696`, `bond→C +0.764`, `charge→CA −0.073`) are
  **total-axis**, dominated by between-atom structure in a 54-atom × 500-frame
  panel.
- `efg_distill_evidence.py` / `EFG_ARC_EVIDENCE.md` fits the within axis (frame
  split, constant-g per stratum) and *already* reports an autocorrelation-based
  `N_eff_lag1` — the only fitter that is axis-honest about effective N today.

**Correction to the lead's read of the +0.70.** The brief states the `charge→N`
+0.70 "is on the DE-MEANED (dynamic) axis." It is not — it is the
`backbone_literature_kernel_t2_corr.py` **pooled/total** number (global-mean
Pearson, no per-atom de-mean; `:131-132,208`). The *fitted* equiv-T2 charge
channel is the de-meaned one, but the +0.696 figure specifically comes from the
total-axis literature-kernel CSV. The substantive point in the brief still holds
— charge clearly carries N-tensor signal and ring does not — but the axis label
on that specific number must be corrected before it is defended. The whole
purpose of this decomposition is to stop quoting numbers whose axis is implicit.

The fix is not to pick one axis. It is to **report every mechanism on BOTH
components** and let the placement be read off, per target, per stratum, with
effective N stated per component.

---

## 1. The decomposition design

### 1.1 Estimator choice — within/between transforms, NOT a mixed model

Recommendation: the classic econometric **within transform** and **between
transform**, run as two separate fits, NOT a random-effects / variance-components
mixed model.

Rationale (defensible, and the honest reason):

- The within and between transforms are the textbook panel decomposition
  (Wooldridge, *Econometric Analysis of Cross Section and Panel Data*, ch. 10;
  Mundlak 1978). A methods reader recognises them immediately and they make the
  axis of every reported number explicit — which is the entire point here.
- A random-effects model *re-pools* between and within through a single GLS
  weight (the quasi-demeaning parameter θ depends on σ²_between/σ²_within). That
  re-introduces exactly the implicit-axis ambiguity we are removing. RE also
  assumes the atom effect is uncorrelated with the regressors — false here:
  static environment (the atom effect) is *correlated with* the static field a
  mechanism feels. RE would be the wrong model and a critic would say so.
- The within transform is already what the current fitters do (per-atom
  de-mean). We are not changing the within fit; we are **adding its between
  twin** and labelling both. Minimal, clarifying, no new machinery.

So, per (mechanism, target, stratum):

```
WITHIN:   x_it - x̄_i.   regressed on   y_it - ȳ_i.     (frame-split gate)
BETWEEN:  x̄_i.          regressed on   ȳ_i.            (atoms ARE the units)
```

where `i` indexes atom, `t` indexes frame, `x` is the predictor (kernel /
field / literature-kernel T2), `y` is the DFT target.

### 1.2 What "predictor" and "target" are, per target channel

Two target channels, kept separate end-to-end (T2 is sacred — never collapse to
scalar):

- **σ_iso / T0** — scalar. Predictor is the scalar kernel/field
  (`field_z`, `|E|²`, `ring_sum_dipolar`, `bond_sum_dipolar`, EFG `|EFG|`,
  literature-kernel reductions). Report **R² per component** (within R², between
  R²). R² is scale-free within a component, so normalisation is moot here
  (see §2).
- **T2 tensor** — rank-2. Predictor is the emitted local-frame literature-kernel
  T2 (or the fitted equivariant readout). The tensor analogue of the within /
  between split:
  - the per-atom **mean tensor** `T̄2_i.` is the *between* unit; correlate the
    mechanism's per-atom-mean kernel-T2 against `T̄2_i.` across atoms.
  - the **deviation tensor** `T2_it - T̄2_i.` is the *within* unit.
  - report, on each component, the **mean-row cosine** and **|T2| Pearson r**
    (rotation-invariant magnitude) as the primary tensor correlations
    (`mean_row_cosine`/`magnitude_r` already exist in the corr script), plus the
    component-Pearson `component_r` — but now **computed per axis** instead of
    only pooled. |T2| r is the headline tensor number because it survives the
    frame-alignment caveat trivially; component_r is reported alongside since the
    frame alignment is verified (`manifest.json` `dft_frame_alignment` max angle
    2.4e-4°). Equivariant fits stay on the frozen `get_C()` path.

### 1.3 Effective N per component — the airtight part

Per-component effective N is where a critic will press, because between has
**~50 units** and within is **autocorrelated**. State both honestly:

- **Between effective N = number of coupled atoms in the stratum** (54 for
  N/CA/C/O, 52 HN, 50/58 HA, **4** for HA2/HA3). This is small. Report it as the
  unit count and FLAG thin strata (`THIN` < 10 atoms, as the fitters already do).
  A between R² on 4 atoms is not defensible as more than a flag — say so. Do not
  bootstrap-inflate it.
- **Within effective N = autocorrelation-deflated frame count**, summed over
  atoms, exactly as `backbone_literature_kernel_t2_corr.py:148-180` and
  `EFG_ARC_EVIDENCE.md` already do: `n_eff = n·(1−ρ₁)/(1+ρ₁)` per atom on the
  target-|T2| (or σ_iso) series, with the lag-1 ρ clamped to [−0.999, 0.999].
  Report the summed `effective_n` AND the median lag-1 ρ per stratum (the
  existing columns). The DFT sampling is **every-other-frame**, so consecutive
  rows are 2 frames / ~40 ps apart; ρ₁ is already mostly small (0.02–0.23 in the
  emitted CSV), which is *why* within N_eff stays near the raw row count — but it
  must be stated, not assumed, because HA/HN have the highest ρ (0.18–0.23).
- A single honest sentence the defense needs: *between asks whether ~50 distinct
  environments order the same way the mechanism predicts; within asks whether one
  atom's frame-to-frame wobble tracks the mechanism. They are different N, and we
  report each.*

### 1.4 Reporting surface (one table, per target channel)

Per (mechanism × target-channel × stratum), one row:

```
mechanism, target, stratum,
  n_atoms (between N_eff), within_N_eff, median_lag1_rho,
  between_R2 [or between |T2| r / cosine], within_R2 [or within |T2| r / cosine],
  thin_flag
```

The point of the table is **placement**, not a leaderboard: a reader sees at a
glance that (e.g.) field-σ_iso sits in the between column and charge-T2-on-N sits
in the within column. No composite score; no ranking; SETI, not editor.

### 1.5 Honest small-N handling (do not over-engineer)

- Between on a 54-atom stratum: report R²/|T2| r with the unit count; optionally
  a leave-one-atom-out R² (already an opt-in `--loao` in the equiv fitter) as the
  *honesty* number, not the gate. Do **not** add permutation tests or BCa
  bootstrap CIs as standard — that is the lookit-me failure mode for a one-protein
  prototype. A single permutation/atom-shuffle null is acceptable IF a number
  looks surprisingly strong and someone will quote it; offer it as opt-in, off by
  default.
- HA2/HA3 (4 atoms): between is uninterpretable — report the count and the flag,
  fit within only. The fitters already flag these THIN; keep that.
- Frame-split stays the within gate (test on held-out frames), as today.

---

## 2. Stage-1 normalisation — resolved (verify, not transplant)

**What I could and could not verify.** The producer copy of `learn/` is gitignored
and not on disk; the live copies are under `.claude/worktrees/**` which this
sandbox denies. I read the **authoritative plan** at
`learn/docs/stage1_plan.md` (worktree copy) and the calibration design docs
`/shared/2026Thesis/nmr-training/{LESSONS_FROM_CALIBRATION,CALIBRATION_VS_PREDICTOR}.md`.
I did **not** get to read the literal `Ridge(alpha=…)` / `StandardScaler` call in
the Stage-1 fitter. So the normalisation specifics below are stated as *what the
design requires* with that caveat flagged; the lead should confirm the one
`fit()` call before quoting it in the thesis. (I am flagging this rather than
reconstructing it — `feedback_huxley_data_discipline` / "I don't know" beats
handwaving.)

**What Stage-1 IS, confirmed from the plan.** Per-element (and per-AMBER-atom-type)
**ridge** on **723 WT-ALA mechanical-mutant pairs**, target = **WT−ALA delta**
tensors, ~55 kernels, R²=0.818 element-pooled (`stage1_plan.md`; `project_calibration_done`).
The WT−ALA delta is a **between-structure** comparison that **cancels the local
backbone baseline** — the static-environment axis. `CALIBRATION_VS_PREDICTOR.md:43`
makes the contrast explicit: the *delta* model cancels backbone conformation; only
the absolute RefDB predictor must carry it. Ridge regularisation (an L2 penalty on
coefficients) is the thing that requires features on a common scale — that is the
only reason a normalisation enters Stage-1 at all.

**The resolution (confirm or correct the lead's read — it is essentially right).**

1. **It does NOT transplant wholesale.** Stage-1's norm (whatever the exact
   per-element scaling/z-scoring is) was chosen so a *single regularised ridge*
   over many kernels on the *between-structure delta* axis had comparable
   coefficient scales. The trajectory decomposition runs the within and between
   axes at very different variance scales (within = frame wobble, often ×10–×100
   smaller than between-environment spread); a single Stage-1-style scaler fit on
   pooled rows would mis-scale one axis. So: **do not import the Stage-1 scaler.**
2. **For the R²-per-component reporting it is largely moot** (confirmed): R² and
   Pearson r and cosine are invariant to a positive rescaling of a single
   predictor *within a fixed regression*, and our headline reporting is
   single-mechanism R²/r per component. So §1's tables need NO normalisation to be
   defensible.
3. **It bites only under a multi-kernel ridge** (the eventual coupled fit). There
   it must be **component-aware**: standardise predictors (and ridge α) **separately
   within the within-transformed data and within the between-transformed data** —
   i.e. compute means/SDs on each axis's own rows, never on the pooled rows. State
   the scaler is fit on train frames / train atoms only (no leakage), exactly as
   `buckingham_efield_t0.py:centred_by_train_atom` already centres on the train
   mask.

So the defensible normalisation, stated per component:
- **σ_iso / T0, single-mechanism:** none needed (scale-free metrics). Center per
  axis (within = per-atom de-mean on train; between = atoms are already the
  units). 
- **T2, single-mechanism:** none beyond the frozen `get_C()` orthogonal map
  (|T2|-preserving). No per-component z-scoring of tensors.
- **Multi-kernel ridge (future):** z-score predictors per axis, on train rows
  only, with the ridge α tuned per axis. Never a pooled scaler. This is where
  Stage-1's lesson applies *in spirit* (regularised → standardise) but not its
  *parameters* (different axis, different scale).

---

## 3. Solvation treatment — the second category error, and what is extractable

**Confirmed against the actual ORCA outputs** (not reconstructed):

- Input line, both Stage-1 and 1P9J: `! r2SCAN def2-SVP def2/J NMR CPCM(Water)`
  (`…WT_nmr.out:223`; 1P9J job `…f000742…_nmr.out:223`).
- 1P9J QM region is **846 atoms, charge −1, protein-only** — no explicit waters
  (`:5389`, `:5517`). CPCM(Water) is an implicit continuum: ε=80.15, GAUSSIAN VDW
  surface, 62822 cavity points (`:5644`); the WT-ALA example shows the same CPCM
  block (ε=80.15, `:3809`).
- So the DFT shielding responded to a **CPCM reaction field** (the polarisation of
  a continuum dielectric by the solute), NOT to vacuum Coulomb, NOT to a
  PB-continuum of FF14SB point charges, NOT to explicit MD waters.

**Our electrostatic predictors are all the wrong treatment of the same physics:**

| emitted predictor | what it is | matches CPCM? |
|---|---|---|
| `field_z` / `charge_field_z` (FF14SB) | vacuum Coulomb of fixed FF point charges | no |
| APBS E-field/EFG | Poisson–Boltzmann continuum of FF point charges | no (PB ≠ CPCM; FF charges ≠ DFT density) |
| `water_efield`/`water_efg`/`hydration_shell` | explicit MD-water field | no (explicit ≠ continuum) |
| ring / McConnell (BS/HM) | through-space *magnetic* susceptibility | **solvation-insensitive — escapes the mismatch** |

The magnetic mechanisms (ring current, McConnell bond anisotropy) are
ring-current / susceptibility effects, essentially unchanged by the continuum, so
the ring/McConnell arc is **not** compromised by this. The mismatch is confined
to the electrostatic kernels (Buckingham field, APBS EFG, charge multipoles).

### 3.a Can the CPCM reaction field be extracted as the one matched predictor?

**No — not from the retained artifacts.** Verified by inspecting both an ORCA
`.out` and the on-disk job/consolidated directories:

- The text `.out` writes only **aggregate** CPCM properties: dielectric energy
  (`-1.119 Eh`, `…WT_nmr.out:3878`), total surface charge (`-2.064`, corrected
  `-2.0`, `:3893-3895`), cavity geometry. It does **NOT** write the per-segment
  surface charges, the segment coordinates, or the reaction-field electrostatic
  potential evaluated at the nuclei. (Grep for surface-charge / potential-at-nucleus
  blocks returns only these aggregates.)
- The job dirs retain **only** `*_nmr.out` + `.pdb` + `.xyz` (+ meta.json) — no
  `.gbw`, no `.cpcm`, no `.densities`/`.scfp`. (Both `1p9j-orcas/jobs/<frame>/`
  and `consolidated/<id>/` confirmed by `ls`.) The reaction-field operator and
  surface charges live in `.gbw`/`.cpcm` scratch that was discarded.

So the matched predictor is **not recoverable from the existing 1P9J campaign**.
Two honest options, neither buildable now:

1. **Re-run cost is prohibitive and unnecessary.** Re-running ORCA to print CPCM
   surface charges (`! KeepDens`, `%cpcm … end`, or `orca_vpot` on the saved
   density to get the reaction-field potential at nuclei) would be a new QM
   campaign — out of scope, and the every-other-frame 750-DFT set is already
   ~a-month of fleet time.
2. **The closest *available* matched-ish electrostatic signal** is the
   **Loewdin / Mulliken charges from the solvated SCF** (`…WT_nmr.out:9560` /
   `:4958`). These ARE polarised by the CPCM field (they are the solvated-density
   population charges), so a field built from *them* would be closer to the DFT's
   own electrostatics than FF14SB charges are. BUT: (a) building a field from
   them in Python is a **recompute** (forbidden — it would be the `q·d/r³` end-run
   the law bans); it would have to be a new C++ emit reading the producer's parse
   of those charges; (b) it is still a vacuum-Coulomb-of-DFT-charges, not the
   reaction field itself. So it is a *better-matched charge source*, not the
   matched reaction field. Flag as a possible future C++ emit, not a fix.

**Verdict for the defense:** the CPCM reaction field is the predictor the DFT
actually responded to, and it is **not in the data we kept**. Therefore the
electrostatic results must be **labelled as mismatched-treatment correlations**,
and hydration must **not** be reported as if it were the matched field.

### 3.b Honest labelling (the defensible framing)

For every electrostatic result, carry an explicit treatment label, e.g.:

- `field_z` → **"FF14SB vacuum-Coulomb proxy (treatment-mismatched to CPCM)"**
- APBS EFG → **"PB-continuum of FF charges (treatment-mismatched to CPCM)"**
- `water_efield`/hydration → **"explicit-MD-water field (treatment-mismatched;
  the DFT had no explicit waters)"** — and specifically: do **not** score
  hydration against the DFT target as if a positive/null result were evidence
  about the DFT's solvation; it is a different solvent model entirely.
- ring / McConnell → **"magnetic, solvation-insensitive"** (the one clean class).

The claim shape stays `correlate-not-match`: an electrostatic kernel correlating
with the DFT target is evidence the *geometry/charge pattern* carries signal, NOT
that we reproduced the reaction field. This keeps the electrostatic arc
defensible instead of silently overclaiming a matched mechanism. The
buckingham fitter already prints a PB-vs-Coulomb note
(`buckingham_efield_t0.py:223-226`) — generalise that into a one-line treatment
label on every electrostatic row of the §1 tables.

---

## 4. Re-interpretation framework (to be FILLED by the implementation, not a verdict)

Read off the **already-emitted** numbers (do not refit). The
`broad_backbone_literature_kernel_t2_correlation.csv` is **pooled/total** axis;
the equiv-T2 fits are **within**. Reading them together (with the axis caveat
that they are not yet the *same* decomposition — that is what §1 builds):

- **charge → N**: total component_r **+0.696**, |T2| r +0.51, cosine +0.63, lag-1
  ρ 0.047 (low → within N_eff ~23.8k of 27k rows). Strong on N's tensor. The
  *fitted* equiv charge channel also lights N. Framework expectation: charge
  carries N-T2 on **both** axes but the decomposition must say how the +0.70
  splits; the low ρ means within is well-sampled. **This is the de-circularising
  headline candidate and the axis split is exactly what must be shown.**
- **charge → CA**: total component_r **−0.073** (null), cosine +0.19. CA is
  local-bonding-dominated (`STATE.md`: CA σ_iso R²=0.055; EFG CA null). Framework
  expectation: CA target is near-**static and through-space-blind** — likely tiny
  on both axes. A clean negative control.
- **bond → C**: total component_r **+0.764**, cosine +0.66 — strongest single
  number; matches `BACKBONE_LAW_EVIDENCE.md` "strongest radial evidence is
  bond-like." Framework expectation: bond anisotropy on the carbonyl C, plausibly
  **between**-led (C=O orientation is a static-environment property) — the
  decomposition tests that.
- **field σ_iso (Buckingham)**: near-null on the **within** fit (HN R²≈0.10).
  Stage-1 found the field effect strong on the **between** (WT-ALA) axis.
  Framework expectation: **field lives on between, not within** — the headline
  case for why the uniform de-mean was wrong. The between fit (§1.1) is the test.
- **EFG → O / C**: within R² 0.33 / 0.12, |T2| r 0.31 / 0.13 (`EFG_ARC_EVIDENCE.md`).
  O carries within EFG signal; N/CA null within. Framework expectation: EFG has a
  *within* component on O/C (already shown) — check whether it ALSO has a between
  component (static EFG environment), which the current within-only fit discards.
- **ring → everything**: total component_r ≈ 0 across backbone strata (ring is
  the aromatic-H mechanism, not a backbone one). Clean null on backbone; the ring
  arc is its own (aromatic-H) stratum.

The framework's one sentence: **magnetic mechanisms (ring/McConnell) are
axis-agnostic and solvation-clean; electrostatic mechanisms split by axis —
field→between (Stage-1), charge-T2→within (N), EFG→within (O) with between
untested — and all electrostatic numbers carry a treatment-mismatch label.** The
implementation fills the between column that currently does not exist and reports
both, per stratum, with both effective Ns.

---

## 5. Codex implementation brief (the lead gates + fires this)

See the separate brief block below (also pasted into the lead summary). It is a
**Python analysis on the emitted substrate only**; no C++, no recompute, no H5.
