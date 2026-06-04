# Stage-1 "learn" understanding — for the current rediscover / 720-ingest work

Comprehension pass over the archived `learn/` tree on the slow spinner
(`/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn`), 2026-06-04.
Read-only, nothing copied. Purpose: bring back the REAL Stage-1 numbers (the
advisor précis is stubbed), understand the 720-WT set ahead of the Stage-2
static-pose ingest, and trace the arc into Stage-2 (calculators-as-shadows) and
Stage-3 (equivariant model).

This is a reading of the archive; it does NOT change any conclusion in
`STATE.md`. Where the archive and project memory disagree on a headline number,
that disagreement is called out explicitly below (it matters for the précis).

---

## 0. The one reconciliation that matters for the précis

Memory / project CLAUDE.md says Stage 1 = **"720 proteins / 446K atoms, 55
kernels, R² = 0.818."** The archive shows that single sentence **fuses two
different runs**. Get this right in the précis:

- **R² = 0.818** is the **110-protein** headline (the "AzimuthalExtraction"
  set, 69,080 matched atoms). It is the **weighted-across-element** fair-set
  number from `twenty_eight_realities_2026-04-10.md` (per-element ridge,
  55 core kernels, per-protein normalised, + two fair categorical scalars:
  kernel-scale and mutation-type).
- **R² = 0.718** is the **720-protein / ~446K-atom** number — the same
  per-element fair-set pipeline re-run at full scale. The drop 0.818 → 0.718 is
  **cross-protein generalisation, not changed physics** (stated verbatim in
  `learn/CLAUDE.md` and `stage1_plan.md` step 6).
- **55 kernels** is a **core subset** label. The actual emitted kernel layout is
  larger: **126 T2 kernels** in the 720 run (`kernels.py` default `per_ring_k=6`
  → 60 + 6·11 = 126; `STAGE1_720_BASELINE` confirms "126 T2 tensor kernels, 316
  scalar columns"). "55 core" = the BS/HM/Disp/PQ per-type + bond-category +
  totals + EFG block without the per-ring-K expansion. Don't say "55 kernels
  gave 0.818 on 720 proteins" — that conflates the subset label, the 110-set
  R², and the 720 atom count.

**Correct précis sentence:** per-element, per-atom-type ridge on WT-ALA mutation
deltas; **R² = 0.818 weighted on 110 proteins (69K atoms)**, holding at
**R² = 0.718 on the full 720 proteins / 446K atoms** (the difference is
cross-protein generalisation). Use 0.718 whenever "720" is in the same sentence.

---

## 1. Stage-1 results — the real numbers and how they were obtained

### Method (one paragraph)
WT-ALA **mechanical mutants**: take a protein, replace each aromatic sidechain
(PHE/TYR/TRP/HIS) with alanine **with the backbone held fixed** (no relaxation).
The DFT shielding-tensor **delta** (WT minus ALA) is the regression target — i.e.
"the shielding effect of deleting the aromatic ring." The extractor computes
rank-2 (L=2, 5-component) geometric kernels per atom; ridge regression
(λ from `calibration.toml`, `ridge_lambda = 1e-2`) fits kernels → DFT-delta T2,
**stratified per element and per atom type**, **per-protein normalised** (strips
magnitude, isolates angular structure), with kernel-scale + mutation-type added
back as fair categorical scalars. Three continuous scalars (MOPAC valence, bond
order, molecular dipole) were tested and **each adds +0.000** — the kernels
already encode them.

### Headline per-element fair-set ridge (the 0.818 table, 110 proteins)
From `twenty_eight_realities_2026-04-10.md` ("progressive calibration"):

| Element | Base (55 norm) | +scales | +mutation-type | **Fair set** |
|---|---|---|---|---|
| H | 0.875 | +0.063 | +0.027 | **0.949** |
| C | 0.520 | +0.069 | +0.121 | **0.691** |
| N | 0.450 | +0.059 | +0.220 | **0.704** |
| O | 0.407 | +0.048 | +0.195 | **0.632** |
| **Weighted** | 0.684 | | | **0.818** |

Pooled (all atoms, one ridge) = **0.385**; weighted-per-element = 0.718 on this
same set's tensor-alignment section (the doc reports both 0.818 fair-set and a
0.718 "weighted per-element vs pooled 0.385, gap 0.333"). The **gap between
pooled 0.385 and per-element is the headline argument**: pooling hides the
physics.

### "Nitrogen is hard" = element-pooling artifact — the per-ATOM-TYPE table
This is the single most-cited Stage-1 finding. Source:
`stage1-mutations/notes/atom_type_stratification.md` (2026-04-15, 720 proteins,
446K atoms, per-atom-type ridge keyed on AMBER atom name from the prmtop).

**Normalised + progressive scalars (the "fair" column is the one to quote):**

| Atom type | n | Base | +Scales | +Mut | **Fair** |
|---|---|---|---|---|---|
| H (one group) | 230,135 | 0.848 | 0.924 | 0.861 | **0.928** |
| C all | 133,488 | 0.471 | 0.529 | 0.512 | 0.562 |
|  CA | 29,944 | 0.539 | 0.577 | 0.597 | 0.627 |
|  **C=O** | 29,944 | 0.372 | 0.411 | 0.430 | **0.463** |
|  CB | 27,429 | 0.544 | 0.597 | 0.604 | 0.647 |
|  **C side** | 46,171 | 0.646 | 0.690 | 0.694 | **0.729** |
| N all | 39,954 | 0.245 | 0.292 | 0.345 | 0.380 |
|  **N bb** | 29,944 | 0.254 | 0.301 | 0.351 | **0.387** |
|  **N side** | 10,010 | 0.720 | 0.762 | 0.870 | **0.887** |
| O all | 42,429 | 0.274 | 0.303 | 0.358 | 0.382 |
|  O bb | 29,944 | 0.303 | 0.338 | 0.395 | 0.422 |
|  **O side** | 12,485 | 0.334 | 0.373 | 0.543 | **0.566** |

The story: **backbone N (0.387) is hard; sidechain N (0.887) is the
second-best atom type after hydrogen.** The pooled N = 0.210 (raw) / 0.380
(fair) was 30K backbone atoms drowning 10K well-predicted sidechain atoms.
Symmetrically, **carbonyl C=O (0.463) is the hardest carbon** (sp2, C=O π*,
paramagnetic) and drags the carbon pooled number down from the 0.63-0.73 the
other carbons reach. Every element shows **sidechain > backbone** (sidechain
atoms sit closer to the deleted ring, so the kernels have more angular signal).

Mutation-type sensitivity (the `+Mut` delta) is itself atom-type-dependent and
is **Case-1995 ring-type intensity factors showing up in the T2 channel**:
N side +0.150, **O side +0.209** (largest), N bb +0.097; H only +0.013 (ring
current geometry alone suffices, ring identity irrelevant for protons).

### Effective-dimension inventory (H≈20 / C≈6 / N≈3 / O≈12)
Source: `physics_analysis_bridge.md` + `dimensionality_honest.md` (the INDEX
points at both). The number = independent predictive angular dimensions the
ridge extracts **after per-protein normalisation** (raw space is 3 for every
element — the shared "how far from a ring" magnitude axis). Normalisation
separates "how far" (magnitude, same for all kernels, protein-size-confounded)
from "which direction" (angular, element-specific).

| Element | dims (normalised) | R²_norm | What carries the diversity |
|---|---|---|---|
| H | **20** | 0.856 | ring current (17 kernels, 8-9 independent ring-type views); EFG adds 2-3; dispersion a few |
| C | **6** | 0.484 | EFG-dominated (MOPAC EFG 0.380 vs ff14SB 0.225); the +0.197 charge-polarisation gap is the carbon story |
| N | **3** | 0.267 | nothing dominates — 5 families each 0.015-0.08, 3 blurred mixtures; "sees N in the corner of its eye" |
| O | **12** | 0.304 | dispersion-driven (0.058 → 0.234 after normalisation — the single largest normalisation effect); r⁻⁶ ring-vertex angular structure |

(Raw R² for context: H 0.921, C 0.518, N 0.215, O 0.246, all at dim=3.)

### Other headline Stage-1 facts worth carrying
- **Distance slopes confirm multipole order to 3 sig figs** (un-fitted, encoded
  in geometry): Biot-Savart **-3.04**, Haigh-Mallion **-3.06** (both ring-current
  dipole, expect -3); Pi-Quadrupole **-5.05** (expect -5). This is a clean,
  citable instrument-validation result.
- **BS ≡ HM** — interchangeable. cos near-field 0.998, far-field 0.9999; R²(H)
  0.772 vs 0.784. Keep both for redundancy, not independent info.
- **Three effective T2 dimensions = three source geometries** (98% variance):
  ring-current magnetic dipole, EFG, bond magnetic anisotropy.
- **|cos(BS, EFG)| = 0.684** across 58K pairs — the two families share ring
  position but encode different angular patterns (proves angular independence).
- **BS exact additivity for fused rings**: T2(TRP5)+T2(TRP6)=T2(TRP9) at machine
  precision (wire-segment law is additive).
- **Mechanical-mutant negative controls pass**: HBond_total R² 0.002 (backbone
  H-bond unchanged → kernel sees nothing); DeltaAPBS_EFG 0.005 (solvation = 2% of
  direct aromatic EFG → the perturbation is local/through-space).
- **dia/para cancel**: mutation perturbs dia and para ~7 ppm each, they cancel to
  ~1-2 ppm total; kernels see the **net** (R² 0.70-0.73 H) not the channels
  (R²<0.05). Para/dia variance ratio 1.02→1.06→1.07→1.20 (H→C→N→O) confirms
  Saito 2010 in the T2 channel. **Echoed independently in current rediscover
  Build-3: "total-T2 is the target, dia/para gauge-dependent."**
- **Nonlinear (RF) headroom follows paramagnetic ordering**: N +0.169, C +0.128,
  O +0.013, H +0.002 — the "per-element gating for Stage 2" hook.
- **HIE (imidazole) is the worst ring type** (self-fit R² 0.062): the symmetric
  circular-loop BS model fails on the asymmetric 5-membered, two-N ring.
- **Per-protein ceiling**: 0.81 within-protein vs 0.35 across-protein. The gap =
  bulk electrostatic environment / protein shape that local kernels don't and
  shouldn't capture → the explicit motivation for Stage-2 ensembles.

### The 7-step argument chain (each step CITED or SHOWN) — `stage1_plan.md`
1. each mechanism has a known multipole order (Pople/JB/Buckingham/McConnell/
   Stone) — cited; 2. multipole order → angular symmetry (Stone Ch.3) — cited;
3. kernels confirm it (BS -3.04, PQ -5.05) — shown; 4. families angularly
independent (|cos|=0.684) — shown; 5. DFT target aligns to the right mechanism
per element (H=ring current Boyd 2002, C=EFG Sahakyan 2013) — shown+cited;
6. calibrated coefficients are physical constants (weighted R²=0.718, 720/446K)
— shown; 7. nonlinear follows paramagnetic ordering (Ramsey 1950, Saito 2010) —
shown+cited.

---

## 2. The 720-WT set itself (context for the Stage-2 static ingest)

**What it is.** A **BMRB-derived single-chain protein set**, each protein paired
as a **WT pose and an all-aromatics→ALA mechanical-mutant pose**. Each pose is a
mode-3-style AMBER-prepared structure (real `prmtop`, `xyz`, `nmr.out`) under
`calibration/{ID}/`. The DFT target is the **WT-ALA static shielding-tensor
delta** per matched atom. There are **723 mutant pairs** on disk; **720** produce
complete feature directories (3 are dropped for missing WT or ALA ORCA shielding
output: `A0A075FQU3`, `A0A7C4ZM98`, `A0A7J2L4W1`).

**The current clean run** (in the 3-week window — this is the live one):
- extractor features: `calibration/features/Stage1BMRB_20260513_topology`
- workup root: `learn/stage1-topology-workup-20260514/`
- Row contract: **720 OK dirs · 425,599 matched atom rows · 126 T2 kernels ·
  788 feature columns (this workup) / 316 scalar columns (baseline) · 23 target
  columns · 475,116 total topology atom rows.** Reproduction audit: 172 checks,
  0 failures. (Earlier same-shape run: `Stage1BMRB_20260509` / baseline
  `2026-05-10`, identical 720 / 425,599 / 126 counts.)

**How Stage 1 consumed it.** `extract.py` → `nmr_extract --mutant` per protein →
NPYs (≈53 arrays incl. `ring_contributions.npy`, `ring_geometry.npy`) → SDK
`nmr_extract.load()` returns a typed `Protein` (every tensor an e3nn-irrep SDK
type) → `mutation_set/{kernels,scalars,dataset}.py` assemble the (M, 126, 5)
kernel block + scalar blocks with **per-protein normalisation** → per-element /
per-atom-type ridge. **Ridge is the model; MLP was tested (0.61 pooled) and
rejected.**

**Matching policy — important for the ingest.** `MutationDeltaResult` treats a
residue-type change as a **hard boundary**: **all WT atoms in a mutation-site
residue are left unmatched** (C++ test `MutationSiteAtomsAreNotBound` pins this).
So the target surface is **"the shielding effect on the surrounding,
non-mutated protein,"** NOT the mutation site itself. Cost (reported, not hidden):
**17,346** chemically-shared skipped atoms (backbone + CB + termini) have an ALA
counterpart but are excluded by policy; **32,171** have no ALA counterpart
(aromatic ring atoms etc.). Matching is by **typed mechanical identity**, never
atom-name regex, never spatial. **If the Stage-2 720 static ingest re-uses these
poses, inherit this policy decision consciously** — the "static axis" the ingest
is built on is exactly this WT-ALA between-atom static surface.

**Why the 720 statics matter for the current work (from `STATE.md`):** the
rediscover Stage-2 1P9J results have a real **between-atom / static-axis
weakness** — 1P9J is one protein, ~24 clean cases per mechanism, 16 rings, and
the "LOAO" centering bug means **1P9J has NO clean between-atom axis**. The
720-WT set is the **only clean between/transferability instrument**: lots of
rings (fattens ring's thin between axis), cross-protein charge validation
(charge's strongest axis is the static one the 720 set is made of), and it is the
**same r²SCAN + same `.out` files** (confirmed: absolute σ present). This is
**arc layer 4 — the transferability/statics pilot** — and it is the next big
build named in `STATE.md`.

---

## 3. The arc — how Stage-1 grounds Stage-2 (shadows) and Stage-3 (equivariant)

The thesis reporting arc (memory `project_thesis_reporting_arc`):
**(1) signal → (2) equation calibrations (a law + coefficient, with confidence)
→ (3) ensemble model (the good mechanisms; AIMNet2 if magic) → (4) equivariant
transferability pilot (720 WT backbones, statics as baseline).** Stage-1 is the
calibration ground; the rediscover kit is doing (2)/(3); the 720 ingest is (4).

**Stage-1 → Stage-2 (calculators as shadows of one dipolar kernel).**
- **Per-element / per-atom-type is law, not preference.** Stage-1 proved pooling
  hides the physics (pooled 0.385 vs per-element). Current rediscover Build-2/3
  **re-learned this the hard way**: a single global ridge "sliced" by stratum
  mis-applies one coefficient set to H and produces wild anti-predictions
  (HN −54); **per-type FITTING rescues it** (HN → +0.578, O → 0.72). `STATE.md`
  cites `feedback_stage1_lessons` / `feedback_no_simplification` for exactly
  this. **Carry forward: fit per-type, do not just slice per-type.**
- **The mechanism rank order transfers.** Stage-1: H=ring current, C=EFG (the
  +0.197 charge-polarisation gap lives in MOPAC), N=blurred-multi-mechanism,
  O=dispersion+EFG. Rediscover's per-mechanism law fits land consistent with
  this: **charge q/r³ is the strongest recovered law** (the "cookie"), **ring
  (Pople) is form-recovered/scale-fitted and lives on aromatic-ring-facing H**,
  **MOPAC-field is the top contributor to the combine** (vindicates MOPAC>Amber,
  matching Stage-1's MOPAC-EFG > ff14SB-EFG by angular direction),
  **McConnell/H-bond can't stand alone → joint fit is their home.**
- **The "shadows" thesis is the Stage-1 redundancy structure made physical.**
  Stage-1: 126 kernels collapse to ~3 effective source geometries (ring dipole,
  EFG, bond anisotropy), with measured redundancy (BS≡HM cos 1.00, EFG families
  cos 0.99). Stage-2's claim — **the calculators are shadows of one classical
  dipolar object; the object recovers what the parts can't isolate** — is the
  same fact stated as physics. Current status: the **unified D_ab-sum combine**
  recovers within-axis (R² 0.43) but the between-axis depth is **PROVISIONAL,
  pending a joint-maths discussion** and the 720-WT between instrument.
- **Normalisation = a within-pose lens** transfers from per-protein (Stage 1) to
  per-frame (Stage 2): strip frame-to-frame magnitude to expose conformational
  angular signal. **total-T2 is the target** in both (dia/para cancel/gauge-vary).

**Stage-1 → Stage-3 (the equivariant learning model).**
- Stage-1's clean instrument (kernels see what theory says; coefficients are
  physical constants) is the **last stable ground** to anchor/train on (memory:
  DFT is the anchor; experimental BMRB shift = ensemble average our 15ns MD never
  samples). Stage-1 numbers (per-element/atom-type R², dims) are the
  **transferability yardstick** the Stage-3 GNN must beat or match.
- **No imposed per-atom frame**: Stage-3 is e3nn-equivariant; emit raw geometry +
  tensor and let equivariance handle rotation (memory `feedback_frames_from_
  physics_not_tests`). Stage-1's T2 spherical-tensor convention (isometric,
  √2 / √(3/2) factors, audited C++↔SDK↔analysis-clean) is the basis the
  equivariant path round-trips through (rediscover's frozen `get_C()` 5×5
  change-of-basis, orthogonal to 1.1e-16, is the C++-side embodiment).
- **AIMNet2 is a learnable ceiling, not a law.** Stage-1: AIMNet2 EFG occupies
  almost the same subspace as MOPAC/ff14SB EFG (canonical corr ~0.99) and "adds
  nothing beyond MOPAC for prediction." Current rediscover refines this for T2:
  **AIMNet2 does not lift on T2** (its magic was on σ_iso/T0; coherent
  hypothesis = AIMNet2 carries the isotropic/local part, T2 anisotropy is
  geometric and the through-space kernels already get it; a T0 companion fit
  would confirm the split). The producer keeps AIMNet2 **always-on, no flag**
  (`feedback_aimnet2_required_no_weasel`); "switchable" is analysis-side only.

---

## 4. Honest notes: stale / couldn't-find / divergences

- **The headline-number docs are "stale" by the 3-week filter but are the
  authoritative source.** `twenty_eight_realities_2026-04-10.md`,
  `element_physics_2026-04-10.md`, and the `stage1-mutations/notes/` working
  notes (atom-type, dimensionality, argument-chain) are dated **2026-04-10/12/15**
  — outside the 2026-05-14 window — but they ARE where the real Stage-1 numbers
  live, and the in-window `learn/CLAUDE.md` + `STAGE1_720_BASELINE_2026-05-10`
  point at them and confirm the 720/0.718 scale-up. `stage1-mutations/` is
  explicitly **"WRAPPED / frozen for the thesis"** (`stage1_plan.md`), so its
  not-recently-modified state is correct, not rot. The genuinely-recent
  (in-window) Stage-1 work is the **topology-workup-20260514** — a
  topology-aware *re-reading* of the same 720 outputs (named-atom signal
  visibility, blocked-stability resampling, dihedral sidecars), explicitly
  **"records which calculator outputs correlate, claims no new physics."** Treat
  it as the audit/figure layer over the frozen result, not a new result.
- **`stage1-mutations/` directory could not be `ls`-listed** (Bash was sandbox-
  denied on this thread after the first call — likely the spinner path tripped a
  block). I read its files directly by path via the INDEX/find listing; the
  `notes/` set named in `docs/INDEX.md` (atom_type_stratification,
  dimensionality_honest, dia_para_decomposition, master_chart, etc.) is the
  complete working-note set. I did **not** open every R script in the
  topology-workup (60+); the `STAGE1_COMPENDIUM.md` row-contract + provenance
  chain maps them and the instruction was not to trawl.
- **Nothing corrupt encountered.** Slow but clean reads. No truncated/zero-byte
  files in what I touched.
- **DIVERGENCE TO FLAG (don't resurrect):** the in-window `learn/stage2/*.py`
  scripts (`sigma_variance_vs_nh_s2.py`, `block_avg_convergence.py`,
  `sigma_essential_dynamics.py`, `markwick_overlay_1p9j.py`) **read
  `trajectory.h5` directly in Python** (the 1P9J σ-variance vs BMRB-5801 NH-S²
  cross-report, Wingens 2003). That is the **old learn workspace** and it
  **violates the current `feedback_no_parallel_h5_in_python` discipline** the
  rediscover kit was built to enforce (Python consumes the emitted substrate; if
  the reader reads the H5, extend the reader). These scripts are a useful record
  of the intended 1P9J dynamic-vs-experimental bridge, but **do not adopt their
  H5-in-Python pattern** — route any such cross-report through the rediscover
  substrate. (One science correction already baked into that script: per-atom
  σ_T0 **variance** in (ppm)² is NOT Lipari-Szabo S²; the earlier "σ-S²" label
  was retracted.)
- **"55 core kernels" has no single canonical list file** I could pin to an
  exact 55-name enumeration; it is the layout minus the per-ring-K expansion
  (`calibration.toml` [kernels] block: 32 ring-type + 10 bond-cat + 7 totals +
  9 EFG + 2 AIMNet2-dyadic = 60 "core-ish"; the round "55" predates the 2
  AIMNet2 dyadics and the 9th EFG). Report it as "≈55 core (per-type +
  bond-category + totals + EFG), full layout 126 with per-ring K=6."

---

## 5. Key artifacts and their spinner paths (point here, don't re-read)

Root: `/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/`

**Headline Stage-1 numbers (frozen, authoritative — outside 3wk window):**
- `docs/twenty_eight_realities_2026-04-10.md` — the 0.818 fair-set table; 28
  physics realities; distance slopes; BS≡HM; negative controls.
- `docs/element_physics_2026-04-10.md` — per-element R² by mechanism + forward
  selection + literature grounding.
- `stage1-mutations/notes/atom_type_stratification.md` — **the per-atom-type
  table** (N side 0.887 / N bb 0.387 / C=O 0.463 / C side 0.729 / O side 0.566).
- `docs/stage1_plan.md` — chapter plan, 7-step argument chain, H20/C6/N3/O12,
  the per-protein ceiling, lessons-for-Stage-2.
- `stage1-mutations/notes/physics_analysis_bridge.md` & `dimensionality_honest.md`
  — the effective-dimension derivation + bridging table.
- `docs/INDEX.md` — index to all 13 working notes (master_chart, dia_para, etc.).
- `EXPERIMENTS.md`, `docs/calibrated_weights_2026-04-10.md` — weight vectors /
  experiment log (not read in full; referenced).

**720-set definition + the live (in-window) Stage-1 re-reading:**
- `src/calibration.toml` — single config; `[paths].features = …/Stage1Results`;
  ridge λ=1e-2; per_ring_k=6; the 6 ring-presence strata.
- `src/mutation_set/kernels.py`, `scalars.py` — kernel/scalar assembly (126
  layout, the √2/√(3/2) T2 packing, per-protein normalisation).
- `stage1-atom-audit/STAGE1_720_BASELINE_2026-05-10.md` — the 720 baseline
  counts + AIMNet2/MOPAC/ff14SB EFG overlap (~0.99) + normalization audit.
- `stage1-topology-workup-20260514/docs/STAGE1_COMPENDIUM.md` — the in-window
  workup's full row contract + provenance + output entry points.
- `stage1-topology-workup-20260514/docs/MATCHING_POLICY_AUDIT.md` — the
  exclude-whole-mutated-residue policy (17,346 / 32,171 skipped).
- `stage1-topology-workup-20260514/derived/` — all the result packets (signal-
  visibility atlas, blocked-stability rollup, backbone multi-source, dihedral
  sidecars); `derived/reproduction_audit/…summary.csv` (172 checks, 0 fail).
- `stage1-atom-audit/README.md` (+ `EXTRACTION_CONTRACT.md`,
  `TOPOLOGY_SIDECAR_CONTRACT.md`) — the atom-level identity contract for the set.

**Stage-2 bridge scripts (in-window; record-only — do NOT adopt H5-in-Python):**
- `learn/stage2/sigma_variance_vs_nh_s2.py`, `block_avg_convergence.py`,
  `sigma_essential_dynamics.py`, `markwick_overlay_1p9j.py`.

**Calibration data (live, on /shared not the spinner):**
- `/shared/2026Thesis/nmr-shielding/calibration/{ID}/` — the 723 pose pairs
  (real prmtop/xyz/nmr.out). `calibration/features/Stage1BMRB_20260513_topology`
  is the current feature run.
