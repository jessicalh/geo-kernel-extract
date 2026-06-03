# Next-session kickoff — rediscover (current 2026-06-01, post-cleanup + calibration)

Continuing `rediscover` on `h5-reader-pysr-spike` — **NEVER MERGE** (experimental spike;
no merge / switch / PR / rebase, ever). Supersedes the older kickoff sections below.

## Read first
- **`STATE.md`** (top) — the applied-maths cleanup arc + the static calibration + the
  equation tally; freshest at top.
- Memory: **`project_thesis_reporting_arc`** (THE steering doc — the 4-layer arc),
  `feedback_applied_maths_over_methodology_caveats`, `feedback_token_economy_codex_codes`,
  `project_rediscover_state`, `feedback_no_python_physics_except_labeled_integrity_test`,
  `feedback_seti`, `feedback_model_is_spine`.
- The session record on disk: `APPLIED_MATHS_AUDIT.md`, `APPLIED_MATHS_AUDIT_codex.md`,
  `FIXES_AUDIT_opus.md`, `VARIANCE_DECOMPOSITION_METHOD.md`, `STATIC_ENVIRONMENT_CALIBRATION.md`,
  refreshed `BACKBONE_LAW_EVIDENCE.md` / `EFG_ARC_EVIDENCE.md`, the capstone charts
  (`analysis/rediscover_capstone_charts.py`). `analysis/PATTERNS.md`, `REDISCOVERY_MAP.md`.

## Where we are — the honest tally (equations, Depth-A / calibration)
- **Ring current (Pople)** σ ∝ k·Σintensity·(3cos²θ−1)/r³, k≈21 — **clean recovered law**,
  de-circularised, LOAO R²=0.62 (non-circular). The headline.
- **Charge field-gradient** σ_T2 ∝ Σq(3d̂d̂/r⁵−I/r³) — **form-recovered, scale-fitted**;
  charge→N T2 between |T2|r=0.78 (q/r³ fix → 0.981 comparator). Strongest between-axis.
- **McConnell** K·χ (canonical, library-identical — F1 was a false alarm) — form-recovered,
  scale-fitted; scalar R²≈0.85.
- **Buckingham** σ_iso≈A·E_proj+B·|E|² — form-recovered, scale-fitted; HN between R²≈0.40–0.48.
- **EFG → T2** — **deflated, not a law** (the lab-frame "O 0.34" was a tumbling rotation
  confound; clean local-frame ≈0).
- Bucket for the static three: **"form-recovered, scale-fitted"** (Jessica) — not yet
  de-circularised to literature coefficients (the γ-vs-units thread is open).
- **Confidence scope: within-protein only** (n=1 protein, ≈500 frames, ≈50 atoms/stratum) —
  jackknife/block-bootstrap, autocorr-aware; NO population (across-proteins) inference.

## The order (thesis-arc layers)
1. Signal — done (between/within variance decomposition).
2. Equation calibrations — done (the tally above; within-protein uncertainty).
3. **Ensemble + AIMNet2 — IN FLIGHT (codex cell, ~2026-06-01 late).** Wire AIMNet2 features
   (256-d embedding = the learnable CEILING; CRG = charge-response-gradient, NOT a
   polarizability — AIMNet2 has no true α; aimnet2-charge source-sweep), the ceiling
   diagnostic (esp. Cα), and the **ensemble honest-best-result** (good mechanisms ±
   AIMNet2, between/within, within-protein). **CHECK that cell's result first** (was it
   committed? what did the ensemble + the Cα ceiling show?).
4. **Equivariant transferability pilot — later, DATA-GATED.** e3nn pilot vs the 720 WT
   ("non-mutated mutant") backbones, statics as comparison. Needs per-backbone WT shielding
   targets through the pipeline (Stage-1 had mutation *deltas*) — scope feasibility first.

## Discipline
ONE typed model in C++; no second model in Python (read emitted substrate, no recompute,
no `trajectory.h5`, reuse frozen `get_C()`). **codex does the coding grind (C++ AND Python)**;
the lead reserves tokens for briefs / judgment / drift / cheap verification. **Applied-maths
errors > methodology caveats** (solvation/basis are disclosable footnotes). correlate-not-
match; within-protein confidence; report effective N + scatter, prototype-honest, no
parametric CI. **Captive self-audit** in every fix/build cell. Oracle byte-parity exit-0 if
C++ touched; ctest **scoped + headless** (`QT_QPA_PLATFORM=offscreen ctest -R
h5reader_rediscover` — NOT the full suite, it flashes GUI windows). codex builds via
`codex exec --dangerously-bypass-approvals-and-sandbox -C <h5-reader>`. NEVER MERGE.

## Data
The 1P9J **750-DFT full set lands ~3 days** (every-other-frame) → re-fire the capstone
charts (`analysis/rediscover_capstone_charts.py`, parameterized) + the variance
decomposition + the static calibration on it. Differencing (#35) stays PARKED (needs a
dense consecutive-frame DFT burst, not more every-other frames).

## Open threads
- The calibration **γ-vs-units** question (does γ≈1 in matched units = literature-predicts,
  or a units artifact) — bucket is "form-recovered, scale-fitted" until pinned.
- The AIMNet2 ensemble/ceiling result (the in-flight cell) — fold into the reporting arc.

---
(Older kickoff content below is historical — the four-layer arc + STATE.md supersede it.)
