# What Are We Doing — A Careful, Provisional Account of the Geometric-Kernel Work

*A working description, written in full, every term defined, ambiguity held open
on purpose. The companion to `doc/emerging/GEOMETRIC_KERNEL_ACCOUNT.md`: that document accounts
for **what the kernels are**; this one accounts for **what the operation does, what
its echo against DFT means and does not mean, and where the genuine ambiguities
live** as of the code and results read on 2026-06-06. It is a description, not a
verdict. Nothing here is a recommendation, a score, or a thumbs up/down. Where the
honest answer is "ambiguous," the document states the exact shape of the ambiguity
and stops there. Where the answer was settled by reading the source, that is
marked **[verified-from-code]**; where it was settled by reading our own results,
**[verified-from-results]**; where it remains genuinely open, **[open]**; where it
is borrowed standard theory we lean on but do not assert, **[standard-theory]**.*

*A note on register, kept from the companion. "Claim" = an assertion we would have
to defend as true of our data. "Account" = a true, situated description that stands
on its own correctness rather than on a contest with a skeptic. The requirement of
this part of the thesis is to **account**, not to **claim** and not to **explain**
(an MSc in methods for structural biology does not owe a mechanistic explanation of
why a result holds). The danger here is not overclaiming; it is failing to talk
about the work honestly. A prediction we cannot explain is fine; a result we cannot
honestly describe is not. This document is the honest description.*

---

## 0. The one-sentence frame, then the long version

The subject of the thesis is the **extraction methods** — what can be cheaply
pulled from a single structural frame and how much of the shielding it accounts
for — not a prediction system and not physics-discovery. So the question "what are
we doing?" is answered at the level of the operation and its evidential status, not
at the level of a final scientific claim. The long version is the rest of this
document, and it does not collapse to a tagline, because the honest content is
precisely in the distinctions that a tagline would erase.

---

## 1. The operation, stated literally

For a chosen *target atom* `i` (the atom whose nuclear magnetic shielding we are
about) and the *source atoms* `j` in a neighbourhood around it, each calculator
forms a small fixed set of numbers — a *kernel* — out of the geometry of the
neighbourhood. The geometry enters through the separation vectors
`ρ_ij = r_i − r_j` (length `ρ`, direction `ρ̂`); the chemistry of each source
enters only through a single scalar *weight* `w_j` (a charge, an induced magnetic
moment, a ring-current intensity). The output is the neighbourhood's low-order
*multipole moment*: a scalar (one number, called **T0**, the `l = 0` isotropic
part), a vector (three numbers, **T1**, `l = 1`), or a symmetric traceless rank-2
tensor (five numbers, **T2**, `l = 2`, the anisotropy). The companion document
derives why these are exactly the spatial derivatives of the Coulomb potential
`1/ρ`, why the `l = 2` member carries the `(3cos²θ − 1)/ρ³` "magic-angle" angular
shape, and what the compression keeps and discards. That derivation is the
mathematical spine and is not repeated here.

Read one way, the kernel is a **physical field** — the actual electric field
gradient, dipolar field, or ring-current field that perturbs the shielding at `i`.
Read the other way, it is a **lossy, fixed-basis geometric descriptor** of
structural context — an `l ≤ 2` projection of a weighted point cloud, in the manner
of a featurization. The companion's §5 establishes that these are not rival
hypotheses but **the same object**, because the low-order field of a source
distribution simply *is* its low-order geometric moment. That identity is the load-
bearing framing, and this document takes it as given and asks the next question:
**what does the comparison of this object against independently-computed DFT
shielding actually tell us?**

---

## 2. What the code actually computes, per calculator — verified, with corrections to the seed

The seed for this work proposed that each kernel is, cleanly, `Σ_j w_j ∂^(l)(1/ρ)`,
with ring/EFG/dipolar at `l = 2` and the Buckingham field at `l = 1`. Reading the
source confirms the *family resemblance* and the `l`-orders, but the literal
computation is richer than a single analytic derivative of `1/ρ` in three of the
four cases. The corrections matter for the account, because an honest description of
"what we are doing" has to describe what the code does, not an idealisation of it.

**Electric field gradient (APBS), `l = 2`. [verified-from-code]**
`src/ApbsFieldResult.cpp`. The kernel is the literal second spatial derivative of
the **Adaptive Poisson-Boltzmann Solver (APBS)** electrostatic potential `φ`:
`E = −∇φ` by central difference on the grid, then `EFG_ij = ∂E_i/∂r_j` by a second
central difference. The result is symmetrised (pinning T1 = 0) and trace-projected
(removing the atom's own self-potential delta-function), leaving a **pure T2**
five-vector emitted as `apbs_efg.npy`. So it is genuinely a field gradient — the
`l = 2` derivative-of-potential the seed names — but of the *solvated continuum*
potential computed on a grid, not an analytic `Σ_j q_j ∂²(1/ρ)` sum. The physical
field reading and the geometric-moment reading coincide here, with a solvation
treatment (CPCM/PB continuum) baked in that the companion notes as a disclosable
footnote.

**Ring current (Biot-Savart / Johnson-Bovey), emitted as full T0+T1+T2.
[verified-from-code]** `src/BiotSavartResult.cpp`. This is *not* an analytic
`(3cos²θ − 1)/ρ³` point dipole. It computes the actual magnetic field `B` of a
**Johnson-Bovey double current loop** (two loops offset from the ring plane, each
carrying half the ring current) via the **Biot-Savart law** summed over wire
segments, then forms the shielding kernel `G_ab = −n_b B_a · PPM_FACTOR` and
decomposes it into the full nine-component spherical tensor (T0, T1, T2). The
point-dipole `(3cos²θ − 1)/r³` form is the **far-field limit** of this object; in
the near field and in the ring plane the double-loop departs from it — which is
exactly what the early findings saw as "in-plane radial breaks down" (and which the
later units audit reframed; see §4). So the ring kernel is a genuine magnetic-field
computation whose far-field shadow is the textbook Pople form, not the Pople form
itself.

**McConnell bond anisotropy — richer than the bare D_ab. [verified-from-code]**
`src/McConnellResult.cpp` + `GEOMETRIC_KERNEL_CATALOGUE.md`. The full emitted
McConnell tensor per bond is

```
M_ab / r³ = [ 9 cosθ · d̂_a b̂_b   −  3 b̂_a b̂_b   −  (3 d̂_a d̂_b − δ_ab) ] / r³
```

three terms, where `d̂` is the (midpoint→atom) direction and `b̂` the bond axis.
Term 1 is generally **asymmetric** and can contribute **T1**; term 2 is symmetric
and contributes **T0 and T2**; term 3 is the pure symmetric-traceless dipolar
kernel `K_ab = (3 d̂d̂ − δ)/r³` — *that* term, and only that term, is the clean
`l = 2` `D_ab` the seed named. The familiar McConnell scalar
`f = (3cos²θ − 1)/r³` is a separate per-bond scalar the code also stores. So
McConnell is **not** reducible to a single `l = 2` member: it carries T0+T1+T2, and
the convention is the canonical pseudo-contact-shift form `σ_ab = (Δχ/3)·M_ab/r³`.
This is the direct, verified reason our own fitting found the McConnell kernel only
~55% reconstructable from a single bond-axis angle — the producer feeds more
angular structure (the `b̂` cross-term, the asymmetry) than a single `cosθ` can
carry. The seed's "McConnell = l = 2" is therefore *true of its dipolar sub-kernel
K and its scalar f, but not of its full emitted tensor M*.

**Buckingham electric-field effect, `l = 1`. [verified-from-code]**
`h5-reader/src/rediscover/BuckinghamEfield.cpp` + `GEOMETRIC_KERNEL_CATALOGUE.md`.
The rediscover feature reads the APBS **E-field** vector (`apbs_E`, V/Å), projects
it into the atom's local backbone frame, and the headline scalar is the axial
component `e_proj = E_local·ẑ`. This is the `l = 1` (vector) object — the first
derivative of the potential, the field itself, not its gradient. The seed's `l = 1`
claim for Buckingham holds exactly. (The Buckingham *shielding* response also uses
the `l = 2` EFG; the field-effect feature specifically is the `l = 1` field.)

**The per-calculator `l`-order picture therefore holds, with one real correction:**
EFG is `l = 2` (pure T2), Buckingham field is `l = 1`, the ring and McConnell
kernels *contain* `l = 2` content (and the ring's far-field shadow is the point-
dipole form) but are emitted as **fuller tensors** (T0+T1+T2) computed from real
field models, not as a bare analytic `l = 2` derivative. The companion's §4
"verification flag" — read `src/` before pinning any per-calculator line — is now
discharged for these four, with this correction recorded.

---

## 3. What the predictor actually ingests — the seed's central open question, resolved [verified-from-code]

The seed could not settle whether our equivariant predictor ingests **raw geometry**
or **our hand-computed kernel `K`**, and flagged that this "changes the whole
analysis." Reading the analysis code settles it, and the honest answer is **it
depends on the path, and the project runs both deliberately, with a sharp boundary
document for each.** The three paths under
`h5-reader/src/rediscover/analysis/`:

**Ring-current equivariant fitter (`equiv_t2_e3nn.py`): ingests RAW geometry.**
Its docstring and `load()` are explicit: inputs are the C++-emitted per-source
`disp_local` (displacement vector), `source_normal_local` (the ring-normal axis),
and invariant scalars `r`, `ring_intensity`, plus the `(atom, frame)` group index.
The `l = 2` features are built *inside* the model by **e3nn**
(`o3.spherical_harmonics` on `r̂` and `n̂`, an `o3` tensor product for the cross
term), summed over sources (Deep Sets permutation invariance). **No precomputed
dipolar kernel is read.** The target is the C++-emitted DFT T2 sidecar, required to
be present (it fails loud rather than recompute the projection in Python). So for
the ring path, the model is given coordinates and discovers the `l = 2` law; it is
not handed the law's value.

**Broad-backbone equivariant fitter (`equiv_t2_backbone_e3nn.py`): ingests RAW
geometry**, generalised to a heterogeneous neighbourhood — rings, anisotropic
bonds, and charge sites at once, each through the *same* shared e3nn angular
machinery `Y2(disp_hat)` but its own per-kind radial channel. Again: emitted
displacement vectors and invariants, e3nn-built `l = 2`, no precomputed kernel.

**EFG / Buckingham equivariant fitter (`equiv_t2_efg_e3nn.py`): ingests OUR
KERNEL.** This is the exception, and it is explicit: the inputs are the producer's
precomputed `efg_feature_T2.npy` (the APBS EFG rank-2 tensor) and its magnitude,
and the model is the "Schur-scalar" `l=2 → l=2` map `pred_T2 = g(|EFG|) · EFG_T2`
— a learned radial **rescaling of our own kernel** against the DFT T2 target. Here
the model *is* handed the kernel and asked only whether a scalar function of its
magnitude carries it onto the DFT anisotropy.

So the resolution is not "raw geometry" or "our kernel" but **both, on purpose, in
named separate paths** — and the distinction is exactly the distinction this whole
document turns on (§4). The ring/backbone paths measure whether a *generic* `l = 2`
basis, fed coordinates, lands on the DFT tensor. The EFG path measures whether *our
specific physical kernel*, rescaled, lands on it. Reading those two side by side is
the closest thing in the existing code to the form-ablation the seed asks for, and
§5 returns to that.

**What PySR fits, and against what. [verified-from-code/results]** The symbolic-
regression path (`pysr_distill.py`, `backbone_pysr_distill.py`) does **not** fit DFT
directly. It distils the per-source function `g` that a permutation-invariant
pooling model learned while reconstructing the **producer's own bare kernel**
(`bare_T0`, the clean analytic ring-current value). It recovered the Pople form
`intensity · (cos²θ − 0.224)/(2.379 r³)`. The project's own findings file is
unambiguous that this is a **pipeline check, not physics**: because
`QtRing.LiteratureIntensity` already holds Giessner-Prettre constants, `bare_T0` is
already a literature Pople/Johnson-Bovey kernel, so recovering the Pople form from
it "demonstrate[s] the NN→PySR pipeline can recover a known analytic function from
input/output pairs… Do NOT headline it." This is the self-echo the seed names, and
the project has already disciplined it.

---

## 4. The echo — what the comparison against DFT means, and does not

This is the live conceptual content. The deliverable's value is in stating it
exactly and leaving open what is open.

**The self-echo (K against our own K) is a tautology. [verified-from-results]** If
we feed `(3cos²θ − 1)/r³` in and a fit returns it, we have learned about our
pipeline, not about physics. The PySR ring result, the EFG-path form, the
"equation fell out" — these are self-echoes when scored against the producer's own
kernel. The project says so in its own words and explicitly de-headlines them. An
honest account *names* these as pipeline checks and moves on; they are not evidence
the physics is right, and we do not present them as such.

**The cross-echo (K against the DFT shielding σ) is not circular — but not free.**
The DFT shielding `σ` is computed by ORCA (r²SCAN/def2-SVP, CPCM water) with no
knowledge of our kernel, so a correlation between `K` and `σ` is a real, non-
circular statement. **But** `σ` is itself a function of the geometry, and `K` is a
summary of the geometry, so *any* geometric feature of the right order will
correlate with `σ` to some degree. There is a **cheap generic geometric baseline**.
The quantity that means something for mechanism is therefore not the raw
correlation but the **excess over what a generic same-order feature would already
give** — how surprised we would be by this echo if the specific physical kernel
(its `1/r³`, its physical weights, its specific axis) were *not* the thing
operating. The seed calls this the **severity** of the echo (a likelihood-ratio
intuition), and it is the right register; the project does not yet report a clean
numerical severity, and producing one is exactly what the form-ablation (§5) is
for. **[open — the severity is named and motivated; it is not yet a measured
number.]**

**The channel inversion (severity runs opposite to raw R²).** The seed's sharpest
observation, and the results bear out its *direction* even as the magnitudes have
since moved:

- The **scalar (T0) echo** against isotropic `σ` tends to be numerically the
  *stronger* correlation but the *more generic*: many scalar "how much stuff is
  nearby" summaries correlate with `σ_iso`, so a strong scalar R² mostly restates
  "shielding depends on environment," which we already knew. Low severity per unit
  correlation.
- The **T2 echo** against the shielding anisotropy tends to be *fainter* but *more
  specific*: a particular rank-2 angular pattern is hard to reproduce by accident,
  so a given amount of T2 correlation is harder to fake. Higher severity per unit
  correlation. The honest bound, which the seed states and we keep: `P(echo | no
  specific mechanism)` for T2 is **low but not zero** — `K`'s T2 and `σ`'s T2 still
  share the same underlying geometry. T2 specificity is a matter of degree, not a
  clean separation.

The consequence is that **the faintness of T2 is not the weakness it looks like,
and the strength of the scalar is not the solidity it looks like** — for mechanism
purposes. This is why the project's reporting standard fixes `|T2| r` (the rotation-
invariant tensor-magnitude correlation) as the headline metric for a rank-2 feature
and treats the scalar-gamma R² as "a diagnostic, not the headline."

**Where the echo actually came back, by stratum. [verified-from-results]** As of
the current, corrected results (`NOW.md`, the trued top of `STATE.md`, the
2026-06-04 audit), on the single protein 1P9J across its 15 ns trajectory:

- On the **within-atom (dynamic) axis** — how one atom's shielding moves as its own
  neighbourhood moves, with the per-atom chemical baseline stripped — the echoes
  that stand are: **charge `q/r³`** (within R² ≈ 0.28), **ring current** (within R²
  ≈ 0.28, concentrated in ~5-7 aromatic ring-facing protons), and a **unified
  `D_ab`-sum combine** (within R² ≈ 0.43, carried by the MOPAC-field and McConnell
  channels, not by charge alone). The equivariant T2 backbone fits report frame-
  split T2 R² in the 0.43-0.76 range across the eight backbone strata, with `|T2| r`
  up to ~0.86 for HN — these are within/frame-split numbers.
- The **standalone nulls** on this cut: the **standalone field** (≈ 0.03-class —
  nonzero coefficient but null-class recovery; its value is as the top *contributor
  to the combine*, not as a standalone law), **standalone McConnell**, and
  **standalone H-bond** (confidence intervals span zero).
- The **EFG → DFT-T2 echo collapsed under a frame correction.**
  `EFG_ARC_EVIDENCE.md`: the apparently strong O and C EFG signals were a **lab-
  frame orientation artifact** (target lag-1 autocorrelation ~0.75-0.86); after
  rotating into the backbone local frame before decomposition, the autocorrelation
  dropped to ~0.05-0.23 and the held-out predictive signal went weak/null. This is a
  concrete case of an echo that *looked* like mechanism and turned out to be a
  rotation confound — exactly the kind of own-bug-first vetting the project holds as
  discipline, caught and recorded.

**The between-atom axis was retracted, and this reframes the strongest-sounding
claims. [verified-from-results]** The single most consequential current fact: the
2026-06-04 "true-LOAO" audit found that the path previously reporting "leave-atoms-
out" / between-atom transfer was **centering each held-out atom by its own mean**,
which measures *within*-modulation, not between-atom recovery. The genuine between-
atom recovery on 1P9J is **approximately null** for the current kernels (charge
0.036, ring −1.0, unified −105/overfit). Every prior positive between/LOAO number
(charge 0.38, ring 0.17, unified 0.26) is **retracted as between evidence** and
reread as a mislabeled within number. The within results stand; the between axis —
the transferability question, "does a coefficient fit on some atoms predict atoms
held out of the fit" — is now declared thin **by construction** on one protein and
is deferred entirely to the planned 720-WT static-pose corpus. This is a place where
"what are we doing" genuinely shifted as technique developed: the work is, right
now, honestly a **within-instrument** demonstration on one protein, with the
between/transferability claim explicitly *not yet made*.

---

## 5. The e3nn caveat, and what would measure the excess [verified-from-code on the mechanism; open on the result]

The seed's e3nn caveat is correct and the code makes it sharp. e3nn builds `l = 2`
features as a **generic, complete basis coordinate** (Weiler et al.'s completeness:
equivariant convolutions are the most general equivariant linear maps between fields
over R³ — `doc/emerging/GEOMETRIC_KERNEL_MATH_LINEAGE.md`). And `σ`-as-a-function-
of-geometry has `l = 2` content generically. So **"the equivariant model predicts
the DFT T2 well and uses `l = 2`" is a prediction win, not mechanism evidence** —
it measures the raw correlation (the within-axis R² and `|T2| r` of §4), never the
excess over a generic feature. The ring and broad-backbone paths, which ingest raw
geometry and let e3nn build the basis (§3), are precisely these generic-basis
prediction measurements. They are real and reportable as **how much of the DFT
anisotropy a cheap equivariant geometric feature accounts for** — a methods result —
but they are not, on their own, evidence that the *specific physical mechanism* (the
ring's Biot-Savart field, the bond's `Δχ`) is the operative one.

The thing that *would* measure the excess is the **form-ablation**: does the
**specific physical kernel** — its `1/r³`, its physical per-source weights, its
specific orienting axis — beat a **generic same-order feature** at predicting the
independent DFT `σ`? The EFG path's `g(|EFG|)·EFG_T2` (§3) is the nearest existing
instance of feeding the specific kernel, and its post-frame-correction result is
weak/null (§4) — i.e., on the EFG channel the specific kernel does *not* currently
beat its own baseline once the rotation confound is removed. The clean head-to-head
— specific kernel vs generic `l = 2` basis, both against DFT, the difference being
the severity — is **named and partially instanced but not yet run as a single
controlled ablation**. **[open — this is the test that would convert "predicts
well" into "the mechanism is visible," and it is the natural next move, not a
present result.]**

---

## 6. The company — how others in the field describe the same ambiguity [verified-from-local-corpus]

Jessica asked for honest neighbours: people whose framing is humble and methods-
shaped, and who leave the same descriptor-vs-mechanism ambiguity open in plain
vocabulary. The local reference corpus (`references-meta/`) supplies several, with
on-disk, page-anchored quotes:

- On the **descriptor side**, the equivariant-ML and GP-kernel literatures use our
  exact `l = 2` channel as a *basis coordinate* and say so flatly: Thomas et al.
  2018 ("rotation orders `l = 0, 1, 2` correspond to scalars, vectors… and
  symmetric traceless matrices"), Geiger & Smidt 2022 (`l = 2` *is* the symmetric
  rank-2 decomposition), Ben Mahmoud et al. 2024 (irreducible tensor components
  "transform independently, hence they can be targeted by independent ML models").
  These license saying our `l = 2` is a mathematically distinct, non-redundant
  channel — a featurization claim — **without** asserting it is the physical
  mechanism. The math-lineage doc is careful that none of these "prove that our
  specific `(3cos²θ − 1)/r³` traversal improves a shielding predictor… that remains
  an empirical claim."
- On the **mechanism side**, the molecular-NMR tradition treats the same forms as
  physics: Case 1995 fits ring-current and electrostatic *parameters* to DFT shifts
  ("slopes near unity, intercepts near zero"), Sahakyan & Vendruscolo 2013 fit the
  RC/EF *hierarchy* per nucleus class, Larsen et al. 2015 (ProCS15) add geometry
  terms as named physical contributions. These are the same `(3cos²θ − 1)/r³` forms
  read as mechanism.
- On the **law-recovery side** — the company for the "an equation falls out" move —
  the staircase doc collects Cranmer et al. 2020 (force laws recovered from learned
  representations), Udrescu & Tegmark 2020 (AI-Feynman), Iten et al. 2020, all with
  the honest note that several are "toy examples" or "controlled physics cases,"
  i.e. faint or in-vitro recoveries, not field-strength claims. Our self-echo (§4)
  sits in exactly this company and is held to the same modesty.

The field, in short, **leaves the descriptor-vs-mechanism question unresolved by
using the same object in both registers and labelling which register it is in.**
The staircase doc's own placement of this project is the honest neighbour Jessica
wanted: *"this project is one more flight of the same stairs. The part that is more
ours is keeping the full tensor, especially T2 anisotropy, on the way up instead of
collapsing immediately to the isotropic shift, and anchoring the fit on DFT
shielding tensors rather than fitting measured shifts directly."* That sentence is a
true, methods-shaped, non-discovery-shaped account of what we are doing, and it is
ours.

---

## 7. So — what are we doing? (the description, ambiguity held open)

We are **building and characterising a method that extracts, from a single
structural frame, the low-order (`l ≤ 2`) multipole moments of an atom's weighted
neighbourhood — cheaply, equivariantly, per frame — and asking how much of the
independently-computed DFT shielding (scalar and, centrally, the rank-2 anisotropy
T2) those features account for.** That much is settled and is a clean methods
statement. The genuine ambiguities, held open and given their exact shape:

1. **Descriptor or mechanism — same object, two readings, and the reading is not
   decided by the correlation alone.** The kernel is at once a physical field and a
   geometric descriptor (companion §5). A strong echo against DFT is consistent with
   both, because `σ` is a function of geometry and so is `K`. The discriminating
   quantity is the *severity* — the excess over a generic same-order feature — which
   is **named and motivated but not yet measured as a number** (§4, §5). Until the
   form-ablation is run, "the mechanism is visible" is *not* something the
   correlations alone establish; "a cheap geometric feature accounts for this much
   of the DFT anisotropy" *is*. We are, right now, entitled to the second sentence
   and not yet to the first.

2. **The severity inversion is real and reframes the headline numbers.** The strong
   scalar echo is the generic, low-severity one; the faint T2 echo is the specific,
   higher-severity one (`P(echo | no-mechanism)` low but not zero). So we do not
   read a large scalar R² as solidity nor a small T2 R² as weakness — we report
   `|T2| r` as the tensor headline and treat scalar R² as a diagnostic. This is
   settled as a *reporting stance*; the exact severity it implies per channel is the
   open number above.

3. **The code-reading resolved what the predictor ingests — and it is split on
   purpose.** Ring and broad-backbone e3nn ingest **raw geometry** and build `l = 2`
   themselves (generic-basis prediction); the EFG path ingests **our kernel** and
   rescales it (specific-kernel test). PySR distils against the **producer's bare
   kernel**, not DFT — a pipeline check the project already de-headlines. This is no
   longer open. **[resolved]**

4. **The per-calculator `l`-orders held up, with one correction.** EFG is pure
   `l = 2`; Buckingham field is `l = 1`; ring and McConnell *carry* `l = 2` but are
   emitted as fuller T0+T1+T2 tensors from real field models (ring = Johnson-Bovey
   Biot-Savart, far-field shadow = Pople; McConnell = a three-term tensor richer
   than the bare `D_ab`). The seed's clean "`l = 2` derivative of `1/ρ`" is the
   *shadow* these compute, not the literal computation. **[resolved-with-
   correction]**

5. **The strongest-sounding result has been narrowed, honestly.** The between-atom /
   transferability echo on 1P9J is **retracted** (it was mislabeled within-
   modulation); 1P9J is now a **within-instrument** and the transferability claim is
   explicitly deferred to the 720-WT corpus. The EFG-T2 echo **collapsed** under a
   lab-frame→local-frame correction. These are not failures to hide; they are the
   own-bug-first vetting working, and they sharpen what we can presently say to: a
   within-protein demonstration, on one protein, that charge-`q/r³`, ring current,
   and a unified through-space combine account for a real, spatially-concentrated
   fraction of the DFT shielding modulation — with the between-protein
   generalisation an open, instrument-gated question, not a present claim.

6. **How the understanding will shift as technique develops.** The form-ablation
   (specific kernel vs generic `l = 2` basis, both against DFT) converts §4's
   "predicts well" into a measured severity — and could land either way per channel;
   the EFG channel's current weak/null specific-kernel result is the honest preview
   that it will not all come back strong. The 720-WT static corpus turns the thin,
   case-study between-axis into a population the between-atom question can be asked
   of with probability rather than anecdote. The full e3nn run expected later will
   add the generic-basis prediction numbers across more strata. Each of these is
   expected to *change the description in specific, named ways* — which is why this
   account is provisional by construction: it is true as of the code and results on
   2026-06-06, and it is meant to be reread, and partly revised, the moment those
   land.

This is what we are doing. The settled parts are stated as settled; the open parts
are left open with their exact shape; nothing is reduced to a verdict, because the
honesty of the methods thesis lives precisely in not reducing it.

---

*Provenance for this document: the per-calculator computations were read from
`src/ApbsFieldResult.cpp`, `src/BiotSavartResult.cpp`, `src/McConnellResult.cpp`,
`h5-reader/src/rediscover/BuckinghamEfield.cpp`, and
`GEOMETRIC_KERNEL_CATALOGUE.md`. The predictor-ingest question was read from
`h5-reader/src/rediscover/analysis/{equiv_t2_e3nn.py, equiv_t2_backbone_e3nn.py,
equiv_t2_efg_e3nn.py, pysr_distill.py}`. Results and their corrections are from
`h5-reader/src/rediscover/{NOW.md, STATE.md (trued 2026-06-04 top),
analysis/FINDINGS.md, analysis/EFG_ARC_EVIDENCE.md, analysis/BACKBONE_LAW_EVIDENCE.md,
analysis/VARIANCE_DECOMPOSITION_METHOD.md}`. The literature company and borrowed
theorems are from `doc/emerging/GEOMETRIC_KERNEL_MATH_LINEAGE.md` and
`doc/emerging/STAIRCASE_SOCIAL_HISTORY.md`. The framing rests on
`doc/emerging/GEOMETRIC_KERNEL_ACCOUNT.md §§2-7.*
