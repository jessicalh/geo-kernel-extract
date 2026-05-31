# Rediscover — first findings (ring current + McConnell), 2026-05-31

Living doc. The question was never a correlation table; it was: can the emitted
substrate be stepped into PySR / an equivariant model and have the classical
kernel fall out as an equation? Short answer for the scalar (T0, σ_iso) case:
**yes for ring current, structurally yes for McConnell.** Honest caveats below.
No hyperbole; numbers are what the scripts printed (`analysis/*.py`,
substrate at `/tmp/rediscover-out/`).

## The modeling problem, stated honestly

Ring current / McConnell are **sums over a variable source set**, but DFT gives
one target per atom, not per source. Scalar PySR maps a fixed vector → scalar
and cannot sum over sources; no (atom,frame) group has a single source (~7
rings, ~30 bonds each). So the pipeline is:

1. a **permutation-invariant pooling model** (Deep Sets; the scalar precursor to
   the e3nn T2 model) assigns the per-source contribution and learns the
   per-source function `g` whose sum reconstructs the kernel;
2. **PySR distils `g`** into a closed-form equation.

Two targets, kept apart:
- **the producer's pure kernel** `bare_T0` (clean, analytic, ring-current-only)
  — well-posed; this is where the FORM is recovered;
- **DFT σ_iso** (the truth, but many contributors) — the physics validation,
  and a deliberate negative control.

## Ring current (aromatic ring-facing H; 41 atoms × 500 DFT frames)

**Form recovery (vs producer kernel).** Sum-pooling `g(r, cosθ, intensity)`, fed
RAW geometry only (never the precomputed dipolar), reconstructs the producer's
BS ring-current kernel to **test R² = 0.972** on held-out frames. Read-out of `g`:
- angular law `g(cosθ)` vs `(3cos²θ−1)`: **R² = 0.985–0.989**
- axial radial `g` vs `1/r³`: **R² = 0.978**
- in-plane radial: **R² = 0.67** — the point-dipole Pople form is known to break
  down in the ring plane (near-field / Johnson–Bovey); the model departs exactly
  where the textbook says it should. Figure: `ring_kernel_recovered.png`.
- far-field (r ≥ 4 Å) intensity-weighted recovery: R² = 0.71, coeff A = 0.139.

**Symbolic (PySR, distilling `g`, far-field).** Pareto front complexity-14:
`intensity·(cosθ² − 0.224)/(2.379·r³)` — the Pople form: intensity-weighted,
`(cos²θ − const)` over `r³`. The const 0.224 is near the ideal 1/3 of
`(3cos²θ−1)/3`; the textbook kernel scores R²=0.70 on the same samples, PySR's
0.68. `pysr_distill.py`.

**Physics validation (vs DFT).** Within-atom (the modulation the kernel claims):
- producer kernel vs DFT σ_iso: r = −0.65 (R² ≈ 0.42) — sign is convention.
- our recomputed through-space sum vs DFT: r = +0.58 → +0.60.
- concentrated where it should be: the atoms with a real moving neighbour ring
  (highest source-sum variance) are the strongly-coupled ones — res 24 HD2
  r=0.94, res 51 HD1 0.89, res 24 HE2 0.87, res 18 HD2 0.78, res 39 HE2 0.76.

**Self ring (Codex caveat, confirmed + characterised).** Nearest source is the
H's own ring: r ≈ 2.48 Å, in-plane, near-constant across frames (CV 0.03). It
carries ~0 within-atom signal vs DFT (r = −0.06) and ~12% of the sum variance;
removing it *sharpens* the through-space signal (0.575 → 0.604). NOTE: it is the
self ring that makes our self-inclusive sum match the producer kernel within-atom
at r=−0.97 — i.e. the producer's BS kernel appears to **include** the
self/bonded ring, contrary to the original Codex read. Resolve with `ring_index`
from the C++ re-run before asserting the producer's convention.

**Negative control (honest).** Sum-pooling rings directly against raw DFT σ_iso
gives **negative within-atom R²** and no kernel in the read-out. Correct and
expected: σ_iso has many contributors (local dia, neighbour bonds, H-bond,
E-field); rings alone cannot and should not explain all of it. This is why the
form is recovered against the producer kernel, then validated against DFT.

## McConnell (backbone amide HN; 52 atoms × 500 frames)

Same pooling recovery, contribution = `chi[bond_category] · h(r, cosθ_bondaxis)`:
- producer MC kernel reconstruction: **test R² = 0.547** (a real gap — see below).
- per-category `chi`: peptide C=O / C–N ≈ 1.0, sidechain / aromatic ≈ 0.5–0.6 —
  the carbonyl dominating the amide-HN anisotropy is chemically sensible.
- geometric read-out `h` vs `(3cos²θ−1)/r³`: **r = 0.92, R² = 0.85**; angular
  `h(cosθ)` vs `(3cos²θ−1)`: R² = 0.81. The bond-anisotropy form emerges.

**Gap to chase:** the producer's McConnell kernel is only ~55% reconstructable
from `(r, single bond-axis angle, category)`. Likely the producer uses a fuller
anisotropy (asymmetry η / a second angle) than the axial `(3cos²θ−1)` we feed.
Worth confirming against `src/McConnellResult.cpp` and, if so, emitting the
second angle in the substrate.

## Depth / effective-N (the honest limit)

- Scalar form recovery is **well-determined**: 110k far-field ring sources, 20.5k
  (atom,frame) ring groups, deterministic producer target. The kernel form is not
  depth-starved.
- DFT physics validation is **within-atom deep** (500 frames/atom) but
  **between-atom thin**: only ~10–13 aromatic H in 1P9J have a neighbour ring that
  moves enough to modulate, i.e. a handful of ring–ring (stacking) pairs.
  Within-1P9J modulation is solid; cross-environment transferability is not yet
  testable from this one protein. Stacked narrowing (single protein, DFT subset,
  C–H only, cutoffs) bites here, not on the form.

## Framing — instantaneous, not a process (credibility2_instantaneous.py)

Ring current has no dynamics: σ = f(geometry) is a STATIC map. The trajectory is
a geometry SAMPLER, not a process to fit; each ORCA frame is one independent
(geometry → shielding) sample. Autocorrelation is not evidence of a process
model — it only deflates the effective number of independent geometry draws
(~398 of 500), which touches the CI, not the model class. The trajectory earns
its keep (still instantaneously) by giving WITHIN-ATOM geometry spread, which
strips the per-atom chemical baseline; a single static structure gives each atom
one geometry and the baselines confound.

So the strong, autocorrelation-free generalization test is **leave-ATOMS-out**:
fit the universal coefficient on a set of atoms, predict the ring-current
modulation of atoms HELD OUT of the fit (independent units). Result: a single
**k ≈ 16.4 ppm·Å³** (stable across random splits, [13.3, 20.7]) predicts held-out
atoms at median within-atom R² = **0.28** (all) / **0.585** (the 7 coupled atoms).
One universal coefficient transfers to unseen atoms, against independent DFT, no
time axis. This is the headline credibility result; the frame-split numbers below
are consistent but the atom-split is the cleaner claim.

## Credibility audit — circular vs real (credibility_check.py)

Asked directly: is this indicative or a setup? Separated, with numbers.

**Circular half (sanity only, NOT physics).** Recovering the Pople form by
fitting the producer kernel `bare_T0` reverse-engineers the producer's own
formula — the `QtRing.LiteratureIntensity` values ARE Giessner-Prettre 1969
constants, so `bare_T0` is already a literature-Pople/Johnson-Bovey kernel.
The R²=0.97 reconstruction and the PySR "equation falling out" demonstrate the
NN→PySR pipeline can recover a known analytic function from input/output pairs.
That is a pipeline check, not evidence the physics is right. Do NOT headline it.

**Real half (non-circular): geometry → INDEPENDENT DFT, out-of-sample.** ORCA
computed the DFT tensors with no knowledge of ring-current models. Fit the
coefficient on TRAIN frames, predict HELD-OUT frames:
- a SINGLE shared ring-current coefficient k = +15.9 ppm·Å³ predicts unseen-frame
  within-atom shielding across all 41 atoms at **out-of-sample R² = 0.333**
  (≈ in-sample 0.330 — no overfitting).
- the 10 most-modulated atoms: **median out-of-sample R² = 0.51**, up to 0.89
  (res 24 HD2), 0.80 (res 51), 0.77 (res 24 HE2) — predicting frames never seen.
- **the per-atom coefficient is consistent across independent atoms**: k_atom
  median 20, range ~18–26 ppm·Å³ for the genuinely-coupled atoms (the flat atoms
  give k≈0 and oos R²≈0, as they should). A universal constant emerging from
  independent atoms is the strong, non-circular signal.
- **not an autocorrelation artifact**: MD lag-1 autocorr median 0.11 ⇒ effective
  N ≈ 398 of 500 frames/atom. The correlations are real samples, not inflated.
- **magnitude is sane, not absurd**: k≈16 ppm·Å³ ⇒ a proton 3 Å axially above a
  ring (kernel 2/27 Å³) shifts ≈1.2 ppm — the right order for a ring-current
  shift. Coefficient → literature-unit comparison still TODO.

**Verdict: indicative, and not absurd.** The geometry predicts independent DFT
out-of-sample with a consistent, physically-sane coefficient. Honest limits:
single protein; the signal concentrates in ~6–8 aromatic H (a few stacking
pairs); the coefficient consistency is across ~8 atoms in one protein;
cross-environment transferability is untested. The flashy PySR "equation" was
the circular half — the credible claim is the out-of-sample DFT prediction.

**To make it genuinely thesis-grade (not circular):** run the kernel against DFT
directly (not via `bare_T0`), ideally with the *literature* coefficient held
FIXED (un-fitted) — if a literature-constant kernel predicts independent DFT
out-of-sample, there is no circularity left at all.

**Confirmed on the canonical substrate (2026-05-31 eve).** Re-ran on the
identity-clean re-extraction (`/tmp/rediscover-out-v2/`, self/bonded rings
excluded by `ring_index`, not a distance proxy; typed CG/CD2 frame anchor). The
leave-atoms-out result on `sum_dipolar_producer_valid` holds and sharpens:
universal k ≈ 21 ppm·Å³ [15, 25], held-out-atom within-atom R² = 0.33 (all) /
**0.62 (the 7 coupled atoms)**. Identity-based self-exclusion didn't move the
verdict — good (it shouldn't, since the self ring de-means away). Still TODO: the
literature-coefficient-FIXED test (the last de-circularising step).

## Equivariant T2 (the tensor, not just the scalar) — equiv_t2.py

The T2 frame caveat resolved (rotation ~1e-4°), so the DFT T2 components are
comparable. Equivariant sum-pooling over through-space ring sources: per source,
three l=2 basis tensors in the LIBRARY basis — Y2(r̂), Y2(n̂_source), Y2sym(r̂,n̂)
— each scaled by a learned radial R(r, intensity), summed; target = the per-atom
de-meaned local-frame DFT T2 (same instantaneous logic: strip the local C–H
baseline tensor, fit the through-space modulation). Equivariant by construction
(l=2 features rotate correctly; radial weights invariant). Needed one additive
substrate field — `source_normal_local` (the dipole axis the l=2 tensor is
oriented by; the scalar cosθ can't reconstruct it).

- **basis check**: my lib_T2 projection == emitted `dft_total_T2` to 4.9e-8 (no
  convention bug between the model's l=2 and the DFT target).
- **T2 component fit (held-out frames)**: R² = 0.44 (5-vec), train 0.448 ≈ test
  0.438 (no overfit).
- **|T2| magnitude (rotation-invariant) modulation**: r = 0.75 (held-out frames).

The ring current is in the proton shielding TENSOR, not only the isotropic
scalar — supports [[feedback_t2_sacred]] (the angular residual is the argument).
Honest limits: this is the frame-split (out-of-sample in time); leave-atoms-out
is the sharper test still TODO. The model is the minimal 3-channel equivariant
ansatz, not a full e3nn gated network (which could extend it). Thin: ring-current
T2 lives in the same ~7 coupled aromatic H.

## What this answers

There IS something in the work: the classical ring-current law is recoverable
from the substrate as a function (R²=0.97 reconstruction; angular 0.99, radial
0.98) and as a symbolic equation (PySR → Pople form), and it explains a real,
spatially-concentrated fraction of the DFT shielding modulation. McConnell's form
recovers too, with a reconstruction gap to run down.

## Rate-of-change as an identifiability aid — tried on what we have (diag_differencing.py)

Hypothesis (instantaneous model, chain rule: Δσ ≈ Σ (∂fᵢ/∂gᵢ)Δgᵢ): differencing
across frame windows decorrelates competing sources that are collinear in levels,
aiding identifiability. Tested on the ring-current cell (the only competing
features we have yet: through-space rings on one atom).
- **Gate 1 (smoothness):** ring-geometry lag-1 autocorr 0.365 (a difference is a
  real derivative); but σ_iso autocorr 0.113 (near-white at the DFT stride →
  Δσ is noisy; the ~20-40 ps DFT spacing is the limit, and σ is only on the
  strided frames while geometry is on all MD frames).
- **Gate 2 (decorrelation):** per-ring-TYPE feature matrix max off-diagonal
  |corr| drops 0.33 → 0.19 (modest real gain). BUT the strongly-collinear pairs
  (|corr|>0.5) only go 0.93 → 0.89 — differencing can't separate rings co-moving
  in lockstep (same residue / fused), and 9% break below 0.3.

**Verdict (indicative, not definitive):** mild aggregate gain, fails on rigid
co-movers, real σ-noise cost. NOT refuted, though — what we can test is the
WORST case (two rings on one atom, often rigidly bonded). The idea's real target
is CROSS-MECHANISM separation (ring vs EFG vs bond), where independent sources
should decorrelate in velocity far more. That needs the heterogeneous substrate
(multiple source types per target atom), which doesn't exist yet — so the
decisive test is deferred to that build, not hand-waved away. Convergence: the
heterogeneous substrate that gives breadth-via-mechanisms ALSO settles the
differencing question. (Jessica's push: technique ≠ identity — keep it open.)

## Open / next (see STATE.md for the C++ side)

1. **Equivariant T2 (e3nn).** Predict the rank-2 tensor from source displacement
   vectors (the un-summed per-source rows carry `disp_local`). GATED on the T2
   Cartesian-frame caveat (ORCA-input vs H5 frame); T0 and |T2| invariants are
   safe meanwhile. This is the actual "equivariant model" half of the ask.
2. **Canonical C++ re-run.** The above used the geometry-proxy self exclusion
   (r ≥ 3 Å). Re-run with `ring_index` + identity-based self/bonded exclusion,
   McConnell bond-axis vector + 10 Å cutoff, typed CG/CD2 frame anchor — then
   re-confirm. The form recovery should be unchanged; per-type analysis gets clean.
3. **McConnell reconstruction gap** — check the producer's anisotropy model.
4. **Coefficient → literature.** Convert the recovered A and `chi` into the
   Pople susceptibility / Δχ units and compare to literature values (correlate,
   not match).

Scripts: `look01_ring_triangulate.py`, `look02_self_vs_throughspace.py`,
`look03_coefficient.py`, `sumpool_t0.py` (neg control), `sumpool_kernel.py`,
`refine_kernel.py`, `pysr_distill.py`, `sumpool_mcconnell.py`. venv with PySR:
`analysis/venv` (gitignored).
