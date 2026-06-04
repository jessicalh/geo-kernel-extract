# McConnell Deep Audit — exhaustive self-error hunt

> **Historical audit record — not current truth (trued 2026-06-04).** Positive
> 1P9J leave-atoms-out/between claims below are superseded by the true-LOAO
> retraction; use `NOW.md` and corrected `STATE.md`.

Scope: read-only audit of the WHOLE McConnell path — calculator (`../src/`)
→ rediscover emit (`src/rediscover/`) → SDK contract
(`python/nmr_extract/_catalog.py`) → consumption
(`analysis/mcconnell_literature_decirc.py`, `analysis/mcconnell_dchi_calibration.py`).
Branch `h5-reader-pysr-spike`. No code change, no merge, no ORCA, no
`trajectory.h5` read. No git archaeology.

Governing principle (the reason this doc exists): when a result is null/weird
the DEFAULT attribution is OUR mistake (pipeline / wiring / frame / units /
exclusions / target choice), NOT "the mechanism is crap." McConnell was a
CONSISTENT Stage-1 contributor (720 proteins, R²=0.818) so the bar for an
"honest null" verdict is high. One of our mistakes was already found + fixed
(`1ceab65`: rediscover summed the target atom's OWN amide/near-field bonds; the
producer excludes them — the strong C "signal" was that own-bond artifact and
collapsed |T2| r 0.755→0.092 once removed). The FAIR (post-fix, `--mc-source-mode
valid`) test STILL does not converge. This audit KEEPS LOOKING for more of our
mistakes, calculator → emit → SDK → consumption, and reserves the final call for
the lead.

Status note on the committed reports: `analysis/MCCONNELL_LITERATURE_DECIRC.md`
and `analysis/MCCONNELL_DCHI_CALIBRATION.md` in this tree are the PRE-fix,
all-source runs (out-dir `/tmp/rediscover-mc-lit-fresh`; C within absT2_r=0.755,
q_CO sign-flipped). They predate the `1ceab65` C++ fix; the scripts now default
to `valid` mode but the committed reports were not re-emitted/re-run. Read the
brief's post-fix numbers (C |T2| r ≈ 0.092) as the current FAIR-test state, not
these committed tables.

---

## 0. What the consumption ACTUALLY eats (disambiguation — load-bearing)

There are TWO competing McConnell-tensor families emitted by the broad path.
Getting this wrong mis-attributes everything, so pin it first:

- **`mc_lit` family** (`mc_lit_T0`, `mc_lit_T2_local_0..4`, and the `_valid_`
  variants) — built by `McConnellSourceLiteratureKernelLocal`
  (`McConnellLiteratureKernel.cpp:50-88`): per-source, Δχ-scaled
  (`σ = -prefactor·q/3 · M_traceless`), summed by
  `sumMcLitKernel`/`sumValidMcLitKernel` (`BroadBackboneSink.cpp:52-87`) and by
  the calibration source-sum (`mcconnell_dchi_calibration.py:207-268`).
  **THIS is what both analysis scripts consume.** decirc reads aggregated
  `mc_lit_T2_local_valid_*` (`mcconnell_literature_decirc.py:28-29,609-614`);
  calibration source-sums per-source `mc_lit_T2_local_*` filtered by
  `mc_source_is_self_or_bonded==0` (`mcconnell_dchi_calibration.py:231-247`).

- **`bond_literature_kernel` family** (`bond_literature_kernel_T2_*`, NPY
  `..._bond_literature_kernel_T2.npy`, and the summed `literature_kernel`) —
  preference (a) read the PRODUCER's `mc_shielding`/`kernel_mc` H5 tensor and
  rotate to local (`h5KernelT2Local`, `BroadBackbone.cpp:66-78,485-489`), else
  (b) `bondKernelT2FromSources` raw-M reconstruction over all broad sources
  (`BroadBackbone.cpp:80-106`). **Neither analysis script consumes this family.**

Consequence: the de-circularisation / Δχ calibration ride the
**rediscover-reconstructed per-source Δχ kernels** (8 Å cutoff, rediscover
self/near-field flag), NOT the producer's canonical `mc_shielding`. The producer
kernel only ever reaches the consumed path indirectly, never directly. Several
candidates below turn on this.

---

## 1. CANDIDATE MISTAKES — ranked by likelihood-of-mattering

Each: location, what's suspect, severity, and whether it could explain the
post-fix non-convergence.

---

### C1 — TARGET ISOLATION: McConnell-only predictor vs the FULL DFT total T2 (THE DEEPEST CULPRIT)

**Location:** `mcconnell_literature_decirc.py:602-645` (target =
`broad_backbone_aggregated_target_local_T2.npy`), `mcconnell_dchi_calibration.py:911-943`
(same target); built by `BuildTarget` (`ExtractionSupport.cpp:45-67`,
`t.total_local = frame.TensorToLocal(s->total_raw)`); emitted at
`BroadBackboneSink.cpp:306-313`.

**What is suspect:** the target is `dft_total` — the FULL DFT shielding tensor
(dia + para, carrying EVERY mechanism: local dia/para, ring current, charge/EFG,
H-bond, AND McConnell bond anisotropy). The predictor is McConnell-ONLY
(`mc_lit`, with even aromatic zeroed — see C7). The scripts correlate a single
minority mechanism's tensor against the total, de-meaned per atom across frames
(`center_by_atom`). They never subtract the other mechanisms, never regress
against a residual, never use the para tensor alone, never restrict to a
McConnell-dominated atom subset. Everything McConnell does NOT model is variance
in the target — i.e. noise the fit cannot reduce.

**Why this is the deepest issue:** contrast the ring story, which DID
de-circularise (STATE.md:11-23). Ring-facing aromatic H are ring-DOMINATED: the
ring-current term is a large fraction of that atom's shielding modulation, so a
ring-only predictor competes against mostly-itself + thin noise. Backbone N/CA/C/O
T2 modulation is multi-mechanism — backbone shielding anisotropy is dominated by
local paramagnetic terms and is sensitive to φ/ψ/χ, H-bonding, and electrostatics,
NOT predominantly by neighbour-bond McConnell anisotropy. So a McConnell-only T2
buried under the total is exactly the "minority contributor drowned by
everything else" failure. The within-atom de-meaning removes each atom's static
baseline but NOT the per-frame modulation of the OTHER mechanisms, which co-varies
with geometry and is large.

**Direct evidence in our own numbers:** the post-fix |T2| r is small/scattered
across ALL strata (the brief's C 0.092; the committed table's small within-atom
component_r even in the all-source artifact: HN 0.177, N 0.077, CA 0.082, O 0.053
— `MCCONNELL_DCHI_CALIBRATION.md:36-46`). Within-atom R² ≈ 0–0.03 everywhere
except the (now-collapsed) C artifact. That is the fingerprint of a true minority
signal, not a wiring bug — a wiring bug usually produces a strong WRONG number
(as the own-bond artifact did), not uniformly-near-zero.

**Severity: HIGH. Could it explain non-convergence: YES — leading honest-null
candidate.** But it is only EARNED as a null after the fixable items below are
cleared, because some of them (C2 cutoff, C3 valid-set consistency) could ALSO
be suppressing what little McConnell signal exists.

**Concrete test that would confirm/refute (reserve for lead):** regress `mc_lit`
T2 against a McConnell-ISOLATED target, not the total. Three defensible routes,
all from already-emitted data: (i) target = DFT `para` T2 local
(`para_decomp`/`para_raw` are emitted to the record but NOT to the consumed
sidecar — needs an emit add of `..._target_para_local_T2.npy`); McConnell is a
paramagnetic-neighbour effect, so the para tensor is a less-diluted target.
(ii) Build a multi-mechanism design matrix from the ALREADY-emitted
`ring_literature_kernel_T2`, `charge_literature_kernel_T2`, and `mc_lit` valid
T2, joint-fit against total T2, and read McConnell's PARTIAL contribution /
unique variance — the present scripts fit McConnell ALONE, which conflates
"McConnell is weak" with "McConnell is collinear-with / dominated-by the others."
(iii) Pre-residualise: subtract the ring+charge literature-kernel prediction from
the total T2, then correlate McConnell against the residual. If McConnell's |T2|
r jumps under (ii)/(iii) but not (i), the non-convergence was target-isolation,
not a dead mechanism.

---

### C2 — CUTOFF MISMATCH: rediscover 8 Å vs producer 10 Å (genuine, fix left it at 8)

**Location:** rediscover broad bond cutoff defaults 8.0 Å
(`main_extract.cpp:183-185`, passed to `MakeBroadBackboneRelationship` →
`nearBackend(CloudKind::BondMidpoints, bond_cutoff_A)`,
`BroadBackbone.cpp:447-450`); standalone mc cutoff also 8.0 Å
(`McConnellNeighborhood.h:30`, `main_extract.cpp:168-170`). Producer
`mcconnell_bond_anisotropy_cutoff = 10.0` (`CalculatorConfig.cpp:50`), used at
`McConnellResult.cpp:140`.

**What is suspect:** the McConnell tensor falls as 1/r³, so the 8→10 Å shell
(r³ from 512 to 1000) carries small per-source magnitude — BUT the genuine
FAR-FIELD McConnell contribution is exactly the part the multipole expansion is
VALID for, and it is the part that is NOT contaminated by near-field self/bonded
artifacts. If McConnell's real, transferable signal lives in the far field
(which is the physically defensible regime), truncating at 8 Å cuts real signal
while keeping the near-field-dominated terms. The number of bonds in the 8–10 Å
shell scales with shell volume (~ (10³−8³)/8³ ≈ 0.95× the inner volume of bonds)
— this is NOT a tiny number of sources, even if each is individually small;
summed across many far-field carbonyls the T2 contribution is not negligible.

**Severity: MEDIUM.** Per-source magnitude is small, but it is the cleanest
(far-field) signal and the truncation is asymmetric (keeps the dirty near-field,
drops the clean far-field). This is itself OUR choice (the fix left it at 8).

**Could it explain non-convergence:** PARTIALLY — it cannot create signal that
isn't there, but it could be SUPPRESSING the only well-conditioned signal. Must
be cleared before C1's null is earned.

**Concrete test:** re-emit broad_backbone at `--bond-cutoff 10.0` (matching the
producer) and re-run decirc/calibration; compare |T2| r and γ_lit per stratum.
Cheap, one knob.

---

### C3 — VALID-SET / NEAR-FIELD FILTER FIDELITY (the "fix" — verify it is exactly the producer's two filters, no more, no less)

**Location:** rediscover flag `mc_source_is_self_or_bonded` set in
`makeBondAttacher` (`BroadBackbone.cpp:357-360`):
`endpointSelf = (st.atom == b.atomIndexA || st.atom == b.atomIndexB)`,
`nearField = (r <= axisNorm * mc_near_field_ratio)`, flag = `endpointSelf ||
nearField`. Producer: `SelfSourceFilter` (`KernelEvaluationFilter.h:216-241`,
`atom_index == source_atom_a/b`) + `DipolarNearFieldFilter`
(`KernelEvaluationFilter.h:174-197`, `distance > source_extent · 0.5` to ACCEPT,
i.e. reject when `distance ≤ 0.5·bond_length`), `near_field_exclusion_ratio=0.5`
(`CalculatorConfig.cpp:70`), `source_extent = bond_length`
(`McConnellResult.cpp:168`).

**What is suspect / verified:** the two filters DO match in form (`endpointSelf`
≙ SelfSourceFilter; `r ≤ 0.5·bond_length` ≙ DipolarNearFieldFilter with the
boundary on the same side, since rediscover excludes `r ≤ …` and producer
accepts `distance > …`). `mc_near_field_ratio` defaults 0.5
(`main_extract.cpp:186-189`), matching. The brief reports the near-field filter
fired on 0 non-endpoint sources, i.e. essentially only `endpointSelf` did the
work. **That "0" is itself suspect and IS a finding:** with `ratio=0.5` and a
~1.3 Å bond, the near-field radius is only ~0.65 Å from the MIDPOINT. The only
atoms within 0.65 Å of a bond midpoint that are NOT one of the two endpoints are
vanishingly rare in a real protein — so "near-field removed: 0" is the EXPECTED,
near-inert result, not evidence the filter is doing meaningful work. This means
the producer's near-field filter is ALSO nearly inert here, and the real
exclusion physics is entirely `SelfSourceFilter` (endpoints).

**The deeper issue the inert near-field exposes:** NEITHER producer nor
rediscover excludes a bond the target is ONE BOND away from but is not an
endpoint of (e.g. backbone N and the previous residue's C=O; or a backbone C and
the next residue's C–N where C is not the listed endpoint). Those sit at
r ≈ 2.0–2.6 Å from the midpoint — well outside 0.65 Å, so they survive both
filters — yet they are in the through-bond regime where the through-space
point-dipole McConnell model is least valid. This is consistent producer↔rediscover
(NOT a divergence), but it means the consumed `mc_lit` sum for N/C/O is still
dominated by near-neighbour amide-plane bonds whose tensor is physically the
weakest part of the model. The own-bond fix removed the endpoint case; the
adjacent-bond case is still there and is shared with the producer.

**Severity: MEDIUM (fidelity verified; the residual concern is shared physics,
not a rediscover-only bug).** One real sub-finding: the standalone
`McConnellNeighborhood` (`McConnellNeighborhood.cpp:213-215`) and composed
`mcReducer` (`ComposedRelationships.cpp:276-290`) do NOT set
`mc_source_is_self_or_bonded` at all ("no self/bonded for bonds") — so the
oracle/standalone `mcconnell` case is UNFILTERED. The consumed analyses use the
broad path (filtered), so this does not taint the decirc/calibration, but it
means the oracle parity check is NOT exercising the same exclusion as the
consumed path — a blind spot if anyone trusts the oracle to validate the broad
McConnell numbers.

**Could it explain non-convergence:** NO for the endpoint fix (it correctly
removed the artifact). The adjacent-bond inclusion could be DILUTING signal but
is shared with the producer, so it is a model-validity caveat, not a wiring bug.

**Concrete test:** verify decirc's aggregated `mc_lit_T2_local_valid_*`
(`sumValidMcLitKernel`) equals the calibration's per-source-summed valid set
component-by-component (the calibration's own self-audit RMS check at
`MCCONNELL_DCHI_CALIBRATION.md:80` was run in ALL mode — re-run it in VALID mode
to confirm the two valid paths agree to FP).

---

### C4 — LOCAL-FRAME CORRECTNESS (not just consistency): is `mc_lit`'s frame the one McConnell's signal lives in?

**Location:** backbone frames `backboneFrameFn` (`BroadBackbone.cpp:159-255`),
builders `BuildBackboneNFrame`/`Ca`/`CarbonylC`/`CarbonylO`/`Ha`
(`LocalFrameBasis.h:110-154`). Both predictor (`mc_lit`,
`McConnellLiteratureKernel.cpp` uses `disp_local`+`bond_axis_local`) and target
(`total_local = TensorToLocal(total_raw)`, `LocalFrameBasis.h:79-85`) live in the
SAME per-atom frame. `TensorToLocal` = `Rᵀ T R`, `R=[x y z]` — verified correct.

**What is suspect:** the vet checked frame CONSISTENCY (both sides use the same
frame) — that is necessary but NOT sufficient. A consistent-but-arbitrary frame
washes out real signal under per-atom de-meaning if the McConnell tensor's
principal axis is not stably oriented in that frame across frames. The backbone
N/CA/C/O frames are anchored on the atom's OWN backbone neighbours (N→CA, C→O,
bisectors). For the CARBONYL C, the frame z-axis is literally C→O — i.e. aligned
with the atom's own C=O bond (`LocalFrameBasis.h:139` even flags "this z-axis IS
the McConnell kernel's reference direction"). So for C, the dominant `mc_lit`
contributor (its own C=O, now excluded by C3; and the adjacent C–N) is highly
correlated with the frame axis — which is part of WHY C showed strong (but
artifactual) signal pre-fix. Post-fix, with the own C=O removed, the remaining
sources are off-axis and the de-meaned modulation in this frame may be small.

The concern is sharper for N/CA/HA: their frames are anchored on backbone
geometry that BARELY MOVES within a folded-protein trajectory (backbone bond
lengths/angles are stiff). If the local frame is nearly rigid across frames for a
given atom, then the per-frame VARIATION the within-atom fit relies on comes
almost entirely from the SOURCE positions moving in lab space (rings, distal
charges), not from anything McConnell-specific. The frame is "correct" (stable)
but may simply not be the basis in which McConnell's per-frame T2 modulation is
large — i.e. the modulation that survives de-meaning is dominated by far sources
moving, not by the near-bond anisotropy McConnell claims.

**Severity: MEDIUM.** Hard to rule out without the test below. Frame math is
correct; the question is whether it's the RIGHT correct frame for THIS mechanism.

**Could it explain non-convergence:** PARTIALLY — a stable-but-uninformative
frame interacts with C1 (target isolation): in a rigid backbone frame, the
McConnell T2's own per-frame variance may be tiny relative to the total T2's
per-frame variance from other mechanisms, so the correlation denominator is
dominated by non-McConnell variance.

**Concrete test:** for each stratum, report the per-atom across-frame STD of the
`mc_lit` valid T2 magnitude vs the across-frame STD of the target total T2
magnitude (both already emitted, both in local frame). If the McConnell T2 barely
moves across frames in this basis while the target moves a lot, the frame is
stable but the mechanism's modulation is genuinely a minority of the target's —
which is C1 confirmed, not a frame bug. If the McConnell T2 moves a lot but is
uncorrelated, look harder at frame ANCHOR drift (e.g. `prevResidueIndex` C anchor
flipping at chain breaks — see C6).

---

### C5 — DE-MEANING AXIS: within-atom centring may remove the McConnell signal itself

**Location:** `center_by_atom` (`mcconnell_literature_decirc.py:118-134`,
`mcconnell_dchi_calibration.py:285-299`); the within axis subtracts per-atom mean
of BOTH x (mc_lit) and y (DFT) across frames. The between axis uses per-atom
means (`atom_means`).

**What is suspect:** McConnell's contribution to a backbone atom is, to first
order, a near-static geometric term (the amide-plane neighbour bonds barely
reconfigure within a folded trajectory). The within-atom de-mean removes exactly
the per-atom STATIC component and keeps only the per-frame MODULATION. If
McConnell's signal is mostly the static across-atom offset (different atoms have
different fixed neighbour-bond geometry) and only weakly modulated frame-to-frame,
then the within-atom axis is fitting the WEAK part and discarding the strong part.
The between-atom axis (54 atom-means) is the one that would carry a static
McConnell signal — and note the committed table shows between-axis component_r is
consistently HIGHER than within (N between 0.963 vs within 0.077; C between 0.995
vs within 0.618 — `MCCONNELL_DCHI_CALIBRATION.md:36-46`). That pattern says the
signal IS more in the static/between structure than the within-frame modulation.

BUT the between axis with only ~54 points and a per-atom chemical baseline is
exactly where circularity/confounding lives (different atom TYPES have different
baselines that correlate with their fixed neighbour geometry). The scripts add an
intercept for between (`gamma_with_intercept`), which absorbs the global offset
but not type-structured confounds.

**Severity: MEDIUM.** This is a target-decomposition choice, partly entangled
with C1.

**Could it explain non-convergence:** PARTIALLY — within-axis may be the wrong
axis for a near-static mechanism; the high between-axis component_r (even
out-of-sample LOAO C between R²=0.30) is the most McConnell-looking number in the
whole table and is being under-weighted by the "within is the clean test" framing
inherited from ring (where within WAS clean because ring genuinely modulates
frame-to-frame via ring tumbling/stacking).

**Concrete test:** report both axes side by side per stratum with effective-N
(already done in the scripts) and explicitly ask whether the between-axis C/N
signal survives leave-one-atom-out with a TYPE-stratified null (shuffle atom
labels within type). If between survives a type-permutation null, McConnell has a
static across-residue signal the within axis was hiding.

---

### C6 — BOND CATEGORIZATION / FRAME-ANCHOR EDGE CASES (PRO, termini, cross-residue direction, sidechain amides)

**Location:** categorization `CovalentTopology.cpp:128-200` + `TagAromaticBonds`
247-270; frame anchors via `prevResidueIndex` (`BroadBackbone.cpp:177-188,232-238`);
`McConnellLiteratureCategory` (`McConnellLiteratureKernel.cpp:19-44`).

**What is suspect (enumerated):**
- **PRO:** no amide H, so HN stratum skips PRO (fine). PRO's ring is saturated
  (not aromatic) so PRO ring bonds are NOT tagged Aromatic and NOT indexed into
  the bond cloud — so PRO sidechain ring bonds contribute nothing to `mc_lit`
  (consistent producer↔rediscover; producer also only indexes the same
  categories). No divergence, but PRO N's frame and its peptide C–N (PRO N is the
  endpoint) are excluded by C3. OK.
- **Sidechain amides (Asn/Gln side-chain C(=O)–N, His ring N):** the Asn/Gln
  sidechain carbonyl is C–O with `dist < 1.35` → `SidechainCO` (q=2.41), correct;
  but the sidechain amide N–C(=O) bond is `SidechainOther` (not C–O) → NOT indexed,
  NOT in `mc_lit`. The producer ALSO does not give SidechainOther a category sum,
  but the producer's FULL M tensor (`mc_shielding`) DOES accumulate
  SidechainOther into `M_total` (`McConnellResult.cpp:235-237`). Divergence is
  only in the unconsumed `bond_literature_kernel` H5 path, so moot for the
  consumed analyses. Worth noting as a latent `bond_literature_kernel`↔producer
  mismatch (C9).
- **Cross-residue C–N direction (sign of b̂):** the McConnell M tensor's term 1
  `9 cosθ d̂⊗b̂` is LINEAR in b̂ (sign-dependent, asymmetric, T1-bearing); term 2
  `-3 b̂⊗b̂` is quadratic (sign-invariant). Producer b̂ = `Direction` = (B−A)/|B−A|
  with A=min(i,j), B=max(i,j) (`Bond.h:30-33`). Rediscover b̂ =
  `(posB − posA)/|posB − posA|` with `bond_atom_a/b` from `QtBond.atomIndexA/B`
  (`BroadBackbone.cpp:328`, `McConnellNeighborhood.cpp:170`). **If QtBond's A/B
  ordering does NOT match the producer's min/max convention, b̂ flips sign**, which
  flips term 1 and the T1 part of the reconstructed M. The traceless-symmetric T2
  that the analysis consumes is built AFTER trace removal of `scale·M`; T2 is
  even under b̂→−b̂ ONLY for the symmetric part of M. Term 1 (`9cosθ d̂⊗b̂`) is
  asymmetric, so its symmetric part `(d̂⊗b̂ + b̂⊗d̂)/2` is b̂-LINEAR and DOES enter
  T2. **So a b̂ sign flip DOES change the consumed T2.** This is a concrete,
  checkable correctness risk. (The `dipolar` scalar `(3cos²θ−1)/r³` is
  cosθ-squared so sign-invariant — only the tensor cares.)

**Severity for the b̂-sign sub-item: MEDIUM-HIGH if QtBond ordering differs from
producer min/max; LOW if it matches.** I could not confirm QtBond's A/B ordering
from the files in scope (it is loaded from the topology sidecar, owned by the
reader's loaders, outside this audit's file list). This is a NAMED gap.

**Could it explain non-convergence:** if b̂ is sign-inconsistent per-bond, the T2
would be partly scrambled and the correlation degraded — a plausible additional
suppressor, NOT a full explanation. The within-bond sign is at least consistent
frame-to-frame (same A/B each frame), so it would not destroy within-atom
correlation, only bias the absolute orientation — interacts with C1/C4.

**Concrete test:** confirm `QtBond.atomIndexA < atomIndexB` (min/max) matches the
producer's `Bond` ordering for every consumed bond; if not, the consumed `mc_lit`
T2 term-1 symmetric part is mis-signed relative to the producer's tensor.

---

### C7 — AROMATIC Δχ = 0 ZEROES MUCH OF THE BACKBONE-ADJACENT BOND SET

**Location:** `McConnellDeltaChiQ` returns 0.0 for Aromatic
(`McConnellLiteratureKernel.cpp:39-43`); calibration excludes category 4
(`mcconnell_dchi_calibration.py:44,237-241`).

**What is suspect:** the decision is physically motivated (RING carries the π
current; an aromatic McConnell term double-counts). But it means that for
backbone atoms NEAR an aromatic residue (His/Phe/Tyr/Trp), the aromatic-ring
bonds — which ARE indexed into the bond cloud and ARE close — contribute EXACTLY
ZERO to `mc_lit`. So the McConnell predictor for those atoms is built from only
the peptide CO/CN/sidechain-CO bonds, while the DFT target DOES feel the ring
current (which is in the total). This is correct de-circularisation IF the ring
term is handled separately — but the McConnell-ONLY fit (C1) does not add the
ring term back, so for ring-adjacent backbone atoms the predictor is missing a
real contributor to the target by construction. This compounds C1: not only is
McConnell a minority, but a chunk of the geometry it "sees" (aromatic bonds) is
deliberately silenced, widening the predictor/target mismatch.

**Severity: LOW-MEDIUM (correct decision, but it sharpens C1).**

**Could it explain non-convergence:** only as a contributor to the C1 dilution,
not on its own. The joint-fit test in C1 (ring + charge + McConnell kernels
together) is exactly what neutralises this.

---

### C8 — TRACE REMOVAL / SYMMETRIZATION ORDER & SIGN (verify across calculator ↔ kernel ↔ target)

**Location:** producer category-T2: symmetrize THEN remove trace
(`McConnellResult.cpp:271-281`); producer `mc_shielding` = `Decompose(M_total)`
with NO explicit symmetrization (Decompose handles sym/antisym internally,
`Types.cpp:25-60`). Rediscover `mc_lit`: scale M, remove trace
(`McConnellLiteratureKernel.cpp:83-84`), then `DecomposeLibrary` (which takes the
symmetric part internally). Target: `DecomposeLibrary(total_local).T2`
(`BroadBackboneSink.cpp:308`).

**What is verified:** `DecomposeLibrary` (`SphericalBasis.cpp:7-37`) is
byte-identical to producer `SphericalTensor::Decompose` (`Types.cpp:25-60`) —
same basis `[xy,yz,zz,xz,xx−yy]`, same isometric norms (√2, √(3/2), 1/√2). The
T2 extraction takes the symmetric part `(s_ij+s_ji)/2 − trace/3` regardless of
input symmetry, so the order of explicit trace-removal vs decomposition does not
change T2 (trace removal only shifts the diagonal uniformly, which the Sxx−T0
step re-does). `mc_lit` removes trace before decompose only to make `mc_lit_T0`
an explicit ≈0 audit channel (verified ≈1e-15, `MCCONNELL_LITERATURE_DECIRC.md:34-45`).
**Sign:** `scale = -prefactor·q/3` (`McConnellLiteratureKernel.cpp:82`); the
calibration de-weights by the SAME scale (`mcconnell_dchi_calibration.py:271-282`,
`SCALAR_PREF`) and the self-audit RMS = 1.8e-8 (`MCCONNELL_DCHI_CALIBRATION.md:80,82`).

**Severity: CLEAN.** No order/sign bug between calculator-style M, rediscover
`mc_lit`, and the target T2.

**Could it explain non-convergence:** NO.

---

### C9 — `bondKernelT2FromSources` / `bond_literature_kernel` RECONSTRUCTION FAITHFULNESS (NOT consumed — but a trap if anyone switches paths)

**Location:** `bondKernelT2FromSources` (`BroadBackbone.cpp:80-106`),
`h5KernelT2Local` (`BroadBackbone.cpp:66-78`), preference logic
(`BroadBackbone.cpp:485-489`).

**What is suspect:** `bondKernelT2FromSources` reconstructs the project-canonical
M (`9cosθ d̂⊗b̂ − 3b̂⊗b̂ − (3d̂⊗d̂ − I))/r³`, `BroadBackbone.cpp:94-96`) over ALL
broad bond sources within 8 Å with **NO self/near-field exclusion at all** (only
`r>1e-9`, axisNorm>1e-9) — this is the SAME all-source semantics the `1ceab65`
fix removed from the consumed path, but it is STILL LIVE here. `h5KernelT2Local`
reads the producer's `mc_shielding` (full M, 10 Å, producer-filtered, ALL
categories including BackboneOther/SidechainOther) and rotates to local. These
two reconstructions differ in: cutoff (8 vs 10), exclusion (none vs producer
filters), category set (4 indexed vs all), and `bondKernelT2FromSources` only
fires as a FALLBACK when `kernel_mc` H5 is absent. So `bond_literature_kernel`
is a self-inconsistent column (sometimes producer-faithful via H5, sometimes an
unfiltered 8 Å reconstruction). The oracle covers the standalone `mcconnell`
scalar path, NOT this broad raw-M reconstruction.

**Severity: LOW for the current question (NOT consumed by decirc/calibration),
HIGH as a latent trap.** If a future analysis switches to
`bond_literature_kernel` / `literature_kernel` / the
`..._bond_literature_kernel_T2.npy` sidecar thinking it is "the producer kernel,"
it will get a path-dependent quantity that is producer-faithful only when the H5
happens to carry `kernel_mc`, and an unfiltered all-source 8 Å sum otherwise.

**Could it explain the current non-convergence:** NO (not on the consumed path).
Flagged so it is not mistaken for a fix.

---

### C10 — SDK CONTRACT vs EMIT vs ASSUMPTION (units / irreps / sign / shape)

**Location:** `python/nmr_extract/_catalog.py:188-234` (mc_* and broad_backbone
specs); emit in `BroadBackboneSink.cpp`; consumption assumptions in the two
scripts.

**What is verified / suspect:**
- `mc_shielding` spec: 9 cols, `irreps="0e + 1e + 2e"`, units `Angstrom^-3`,
  sign `σ_ab = -dB_a^sec/dB_{0,b}` (`_catalog.py:189-190`) — matches
  `McConnellResult::WriteFeatures` `PackFull9` of `Decompose(M_total)`. OK.
- `mc_category_T2`: 25 cols (5 cat × 5 T2), `irreps="2e"`, `Angstrom^-3`
  (`_catalog.py:191-192`) — matches. OK.
- `broad_backbone_aggregated_target_local_T2`: shape (rows,5), `irreps="2e"`,
  units `ppm`, `mechanism=quantum_reference` (`_catalog.py:232-233`) — matches
  the sink write of `DecomposeLibrary(total_local).T2`. OK. **Unit mismatch in
  spirit:** the target is `ppm`; the consumed `mc_lit` columns are `ppm` (Δχ
  scaled, `McConnellNeighborhood.cpp:91-96` schema marks them `ppm`). So
  predictor and target are both ppm — γ_lit is dimensionless and the verdict
  "compatible with 1" is meaningful. OK.
- **The `mc_lit_*` and `..._valid_*` aggregated columns are NOT in `_catalog.py`
  at all** — they live only in the CSV (`BroadBackboneSink.cpp:148-187`
  `kAggregatedHeader`), not as ArraySpecs. The catalog's broad_backbone entries
  cover only the NPY sidecars (target_T2, target_local_T2, field_local,
  literature_kernel_T2 ×4 — `_catalog.py:230-245`). So the consumed `mc_lit`
  columns bypass the SDK contract entirely (they are read directly from the CSV
  by pandas). This is a CONTRACT GAP: the primary consumed feature has no
  ArraySpec, no declared units/irreps/sign in the single-source-of-truth. A
  future consumer cannot discover them via the catalog. Not a correctness bug for
  the current scripts (they hard-code the column names), but it means the
  load-bearing predictor is undocumented in the contract that CLAUDE.md says
  "every new output file needs an entry."
- **`literature_kernel`/`bond_literature_kernel` NPYs ARE in the catalog**
  (`_catalog.py:240-245`) with `units="ppm"`, `mechanism` mixed — but as C9
  notes these are path-dependent and NOT consumed. The catalog presents them as
  clean `ppm` tensors, masking the all-source-vs-producer ambiguity.

**Severity: LOW for correctness, MEDIUM for governance.** The consumed feature
(`mc_lit_*_valid_*`) being absent from `_catalog.py` is a real contract gap and a
discoverability hazard.

**Could it explain non-convergence:** NO.

---

### C11 — DFT TARGET CORRECTNESS (right tensor / sign / frame / basis) — secondary checks

**Location:** `BuildTarget` (`ExtractionSupport.cpp:45-67`); `total_raw` from
`run.dft.AtomShielding(atomIdx, originalIndex)`; `dftRows`/`originalIndex` map
(`BroadBackbone.cpp:522-523`); frame-alignment check noted in the vet at
`main_extract.cpp:262-270`.

**What is verified / suspect:**
- The target is the DFT TOTAL shielding tensor, decomposed in the SAME library
  basis as the kernel (`DecomposeLibrary(total_raw)`), and the local-frame
  variant uses the same `TensorToLocal`. Basis-consistent. OK.
- `present & local_frame_valid` gate (`mcconnell_literature_decirc.py:638-641`)
  correctly drops rows with no DFT or invalid frame. OK.
- **Frame/tumbling caveat (memory `project_rediscover_state` says "T2-frame
  caveat RESOLVED"):** the DFT tensor is in the lab frame of THAT frame's
  geometry, and the local frame is built from THAT frame's atoms, so the lab→local
  rotation is internally consistent per (atom, frame). The within-atom de-mean is
  over frames of the same atom, all expressed in their own per-frame local frame
  — so the de-meaning is comparing tensors that have each been rotated into a
  frame that itself reorients with the molecule. This is the SUBTLE point: if the
  local frame tracks the molecule's tumbling, the per-frame local T2 has the
  tumbling rotated OUT (good — that is the point of the local frame), and what
  remains is intramolecular reconfiguration. As long as the local frame is built
  from local atoms (it is), the residual is intramolecular. No lab-frame tumbling
  leak found. CLEAN, consistent with the vet's Suspect-2 = CLEAN and the memory
  "RESOLVED."
- **Sign of DFT shielding vs kernel:** the catalog sign convention
  `σ_ab = -dB^sec/dB_0` is the SHIELDING sign; the McConnell `mc_lit` PCS tensor
  carries the `-prefactor·q/3` shielding sign. Both are shieldings, so γ_lit ~ +1
  is the de-circularised expectation. No sign inversion found between target and
  predictor.

**Severity: CLEAN.**

**Could it explain non-convergence:** NO (target is the right tensor in the right
basis/frame; its PROBLEM is C1 — it is the TOTAL, not McConnell-isolated).

---

## 2. THE DEEPEST CULPRITS (called out)

1. **C1 — Target isolation (McConnell-only vs full DFT total T2).** The single
   most likely reason a real Stage-1 mechanism reads as non-convergent here.
   Backbone T2 is multi-mechanism; ring de-circularised because ring-facing H is
   ring-DOMINATED, McConnell does not get that luxury at backbone atoms. The fit
   has never been given a McConnell-isolated target or a joint multi-mechanism
   design matrix.
2. **C2 — 8 Å cutoff (vs producer 10 Å).** Asymmetrically truncates the CLEAN
   far-field while keeping the dirty near-field; could be suppressing the only
   well-conditioned McConnell signal. Cheap to test, must be cleared before C1's
   null is earned.
3. **C5 — Within-axis de-meaning may be the wrong axis for a near-static
   mechanism.** The between-axis component_r / LOAO numbers are the most
   McConnell-looking signal in the table and are being under-weighted by a
   "within is the clean test" heuristic imported from the ring case.
4. **C4/C6 — Frame informativeness + b̂-sign fidelity.** A stable-but-rigid
   backbone frame (C4) and an unconfirmed QtBond A/B ordering (C6) are
   second-order suppressors that interact with C1.

---

## 3. HONEST, EARNED VERDICT (final call reserved for the lead)

After this exhaustive calculator→emit→SDK→consumption hunt, the non-convergence
is MOST LIKELY a combination of one earned-null cause and at least two remaining
FIXABLE mistakes of ours, and the null is **NOT yet earned** because the fixable
items have not been cleared:

- **The leading remaining FIXABLE mistakes (ours):**
  (a) **C1 target isolation** — we are regressing a minority mechanism against
  the full total. The concrete confirming test: joint-fit the already-emitted
  ring + charge + McConnell literature-kernel T2 against the total T2 and read
  McConnell's PARTIAL/unique variance, AND/OR residualise the total by the
  ring+charge prediction before correlating McConnell. If McConnell's |T2| r
  rises materially under either, the "null" was target isolation, full stop.
  (b) **C2 cutoff** — re-emit at `--bond-cutoff 10.0` to match the producer; the
  fix left it at 8. One knob, cheap.
  (c) **C5 axis** — explicitly elevate the between-axis (static across-residue)
  read with a type-permutation null; the static McConnell signal may be real and
  hidden by the within-axis framing.

- **The honest-null candidate (C1, the minority-contributor / can't-isolate-the-
  mechanism reading):** is a LEGITIMATE possible outcome — backbone-heavy-atom
  T2 may simply be dominated by local paramagnetic / electrostatic terms with
  McConnell a true minority — BUT it can only be DECLARED a null AFTER (a)/(b)/(c)
  are run. Declaring it now would be the exact arrogance the governing principle
  forbids: attributing outward before exhausting our own errors. McConnell's
  Stage-1 record (720 proteins, R²=0.818 contributor) means the mechanism is
  real; the bar for "honest null on this one-protein per-frame local-T2 test" is
  high, and the test as currently wired (McConnell-only vs total, 8 Å,
  within-axis) is not yet a fair adjudication of it.

- **Confirmed CLEAN (not our bug):** T2 basis/normalisation (C8, C10 units),
  trace-removal order/sign (C8), DFT target tensor/frame/basis (C11), the
  endpoint self-source fix itself (C3 — verified faithful to the producer's
  SelfSourceFilter). The near-field filter is near-inert in BOTH producer and
  rediscover (C3) — a shared model-validity caveat, not a rediscover-only error.

- **Latent traps to NOT mistake for fixes:** the unfiltered all-source
  `bondKernelT2FromSources` / path-dependent `bond_literature_kernel` (C9, not
  consumed); the unfiltered standalone/composed `mcconnell` oracle path (C3); the
  consumed `mc_lit_*` columns being absent from the SDK catalog (C10).

Recommended order for the lead: run C1's joint/residualised test and C2's 10 Å
re-emit FIRST (both cheap, both could resurrect a fair signal), surface C5's
between-axis read, and confirm C6's b̂ ordering — THEN, only if McConnell's
isolated/partial T2 is still near-zero out-of-sample, is "honest null:
minority contributor at backbone atoms, not our bug" an EARNED conclusion.

---

## Appendix — file/line index used

Calculator: `../src/McConnellResult.cpp:54-98,105-302,360-395`,
`../src/MopacMcConnellResult.cpp:47-89,100-291`,
`../src/KernelEvaluationFilter.h:151-241`, `../src/CovalentTopology.cpp:128-270`,
`../src/Types.cpp:25-94`, `../src/CalculatorConfig.cpp:50,57,70`,
`../src/GeometryResult.cpp:19-27`, `../src/Bond.h:22-33`,
`../src/SpatialIndexResult.cpp:46`.
Emit: `BroadBackbone.cpp:66-106,159-414,418-498,500-559`,
`BroadBackboneSink.{h:45-139,cpp:34-330}`,
`McConnellLiteratureKernel.cpp:19-88`, `McConnellNeighborhood.cpp:64-221`,
`ComposedRelationships.cpp:225-354`, `LocalFrameBasis.h:79-154`,
`SphericalBasis.cpp:7-37`, `RediscoverTypes.h:33-157`,
`SpatialIndexSet.cpp:19-29,87-108`, `Catalog.cpp:67-72,219-223`,
`main_extract.cpp:165-189,431-526`.
SDK: `python/nmr_extract/_catalog.py:188-245,352-355`.
Consumption: `mcconnell_literature_decirc.py:28-29,118-134,602-698`,
`mcconnell_dchi_calibration.py:28-65,207-282,285-456,873-998`.
Context: `MCCONNELL_PIPELINE_VET.md`, `analysis/MCCONNELL_DCHI_CALIBRATION.md`,
`analysis/MCCONNELL_LITERATURE_DECIRC.md`, `analysis/FINDINGS.md:70-83`,
`analysis/PATTERNS.md`, `STATE.md:1-90`.
