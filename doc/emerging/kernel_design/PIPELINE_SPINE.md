# Pipeline spine — extraction → reader → model, held whole

*The foundational knowledge-spine for the numeric pipeline of the NMR
chemical-shielding work. Written to be read, then to **become** the knowledge that
informs how we build the pipeline. It holds the whole flow in one frame — the
feature set, the equivariant architecture that puts it together, the two ends
(DFT-anchored shielding emulator; experiment-facing shift predictor) — and its
central job is to **hold every ambiguity at once and triage each: cruft (accidental
complexity we can clear) or realism (a genuine open question we must honour)**, with
the reason. The worth of the pass is the **delta**, not the restatement: the
surprises, the additional feature categories and scalars, the MOPAC quantities we
should start using, and whether the architecture itself changes once a rich
local-electronic scalar stream sits beside the geometric tensors.*

*Companions held whole while writing: the `kernel_design/*.md` set (ring, mcconnell,
charge_efg, bond_anisotropy, dispersion, hbond_larsen, pi_quadrupole,
mopac_extension, pipeline_adaptation), `CONTINUITY.md`, `doc/emerging/GEOMETRIC_KERNEL_ACCOUNT.md`,
`doc/emerging/WHAT_ARE_WE_DOING.md`, `doc/emerging/GEOMETRIC_KERNEL_MATH_LINEAGE.md`,
`doc/emerging/STAIRCASE_SOCIAL_HISTORY.md`, `doc/emerging/STATS_PROGRAM.md`, the Stage-3 vision
(`project_stage3_equivariant_gnn`), the four-layer reporting arc, and the actual
e3nn fitter code (`equiv_t2_backbone_e3nn.py`, `equiv_t2_efg_e3nn.py`,
`aimnet2_ceiling_ensemble.py`, `e3nn_protocol.py`). Method note: WebFetch was denied
in this pass; literature deltas are grounded on search summaries + the held corpus,
and that limit is stated where it bears. No outcome is forecast anywhere in this
document — the sweep determines outcomes; "would help" always means **structurally
carries signal the others cannot**, never a performance prediction.*

---

## 0. The premise, in one breath

A MOPAC-extended set of calculators provides a coherent **e3nn input on a defensible
foundation**: ~five kept geometric kernels emitted as clean equivariant
spherical-tensor features (`0e`/`1o`/`2e`); the **full MOPAC capture** as the
local-electronic complement the through-space kernels structurally cannot reach; and
AIMNet2. MOPAC is the default source; ours overrides only where ours is demonstrably
better (geometry and the through-space kernels are ours). The two ends are a
DFT-anchored, T2-preserving **shielding emulator** (Model 1) and a wild-west
**shift predictor** against experiment (Model 2, "fido & biscuits"). This document
is the connected account of that flow and the triage of its open seams.

The spine has a shape worth naming before the parts, because the shape is the
intelligent thing and the parts hang off it:

> **The through-space geometric tensors and the local-electronic scalars are
> complements by construction, not competitors.** The geometric kernels are
> `l ≤ 2` multipole moments of the weighted neighbourhood — they carry exactly the
> part of shielding that is *through-space and low-angular-order*, and they are
> structurally blind to the local-quantum term (hybridisation, the paramagnetic
> excited-state sum) that dominates heavy-atom shielding
> (`doc/emerging/GEOMETRIC_KERNEL_ACCOUNT.md` §6). MOPAC's electronic descriptors are a cheap
> proxy for exactly that local structure. So the architecture is not "more features";
> it is **two orthogonal coverage regions stitched at the atom** — tensors where the
> field reaches, scalars where the electrons sit. Everything below is the literal
> form of that stitch.

---

## 1. The complete feature set into the model

Three parallel streams enter at every target atom. They are kept distinct on purpose
(parallel channels, fit-weighed — `doc/emerging/STATS_PROGRAM.md`; `pipeline_adaptation.md` §0),
and the architecture (§2) is what combines them.

### 1.1 Stream A — the geometric through-space tensors (the kernels), as irreps

Each kept kernel is emitted as a clean equivariant spherical tensor, raw in the
molecular frame, geometry unscaled, with the uncertain physical scale carried as a
calibratable coefficient or a parallel channel. The irrep content per kernel, read
off the design docs:

| Kernel | Primary irrep | Also | Physical object | Scale knob |
|---|---|---|---|---|
| **Ring current** | `1x2e` (CSA anisotropy) | `1x0e` (iso shift, the Haigh–Mallion number, a built-in sanity anchor) | distributed-source J-B/H-M field + Boyd–Skrynnikov full tensor | per-ring-type intensity (Giessner-Prettre / Case 1995) |
| **McConnell** (bond anisotropy folded in) | `1x2e` per source category | `0e`/`1o` only where a real trace/asymmetry exists | symmetric-traceless propagator `D(r)=(3r̂⊗r̂−I)/r³` contracted with axial (or rhombic for C=O/C=C) Δχ tensor | per-category Δχ (Williamson–Asakura / Abraham; ~2–5× source scatter → fitted) |
| **Charge / EFG** | `1o` (E-field) **⊕** `1x2e` (EFG) | `0e` invariants \|E\|², EFG asymmetry η | Coulomb field + field-gradient; EFG is *pure* `2e` by Laplace (the single cleanest claim in the family) | charge-source channels (ff14SB / AIMNet2 / **MOPAC**) × screening (vacuum / APBS-PB / explicit-water) |
| **H-bond / Larsen** *(carry-along)* | `0e` scalar (Barfield/Larsen ppm) | `1o` H-bond direction + `0e` distance powers; `2e` CSA change where a tensor readout exists | the odd sibling: an empirical scalar law in ppm, **already a calibration** | per-acceptor-type table; 2.07 ppm water offset |
| **π-quadrupole** *(carry-along)* | `1x2e` (quadrupole-sourced EFG) | `1o`/`0e` field-order (the "A-term") | ring's *permanent* electrostatic quadrupole, Stone T-tensor | per-ring-type Θ; **load-bearing partition vs charge/EFG** |
| **Dispersion** *(carry-along)* | `1x0e` (⟨E²⟩, R⁻⁶ contact scalar) | weak `2e` tail | London ⟨E²⟩ near-contact sum; **no DFT anchor** (campaign wrote total σ only) | per-element C₆·B scale |

The solid core is **ring, McConnell, charge/EFG**. The three carry-alongs ride per
the two-context inclusion test (`CONTINUITY.md`): a feature is cut only when it earns
nothing in *either* the shielding-tensor work *or* the shift predictor.

Two cross-cutting facts the substrate already records and the model must honour
(`pipeline_adaptation.md`): every tensor is **raw molecular-frame** (no imposed
per-atom local frame; equivariance handles rotation — `feedback_frames_from_physics_not_tests`),
and parity is declared correctly to e3nn (a shielding/susceptibility tensor is
even, `2e`; an outer-product-of-displacements tensor is not — a correctness point,
not a physics ambiguity).

### 1.2 Stream B — the MOPAC local-electronic descriptors (the complement)

This is the stream the geometric kernels structurally cannot reach, and the
`mopac_extension.md` capture is the contract. The descriptors, by form:

**Scalars (`0e`), per atom — the workhorses of the complement:**
- **Mulliken net atomic charge** — conformation-responsive; feeds the charge/EFG
  field as a parallel charge channel *and* rides as a bare `0e`.
- **s and p shell populations** (and d/f/shell-resolved where printed) — the direct
  read on local hybridisation, the quantity the paramagnetic term turns on and the
  geometric kernels are blind to.
- **per-atom electron density** and **per-atom dipole-contribution** fields.
- **MOPAC diagonal valency** (separate from the project's derived valency).
- **AO-resolved populations and overlap populations** (per AO, and per AO-pair) —
  the finest-grained local-electronic descriptors in the capture.

**Bond/pair-axis quantities:**
- **Wiberg bond orders** (full directed object, valency diagonals, hydrogen bonds
  under `ALLBONDS`) — the McConnell source weight, conformation-responsive.

**MO/LMO and matrices (high-dimensional, optional):**
- LMO energies/occupations/bonding contributions, MO coefficient matrices, density
  and overlap matrices when emitted.

**Graph-level scalars (conditioning metadata, not per-atom):**
- heat of formation, electronic / core-core energies, **ionization potential /
  frontier-orbital summary**, the full molecular dipole table (POINT-CHG / HYBRID /
  SUM).

The honest shape: most of Stream B is `0e` (local-electronic descriptors are
rotation-invariant scalars by nature). It has **no native `2e` content of its own** —
which is exactly why it is a *complement* to Stream A and not a competitor: it adds
the scalar coverage region, the tensors stay Stream A's job. The one place MOPAC
becomes tensorial is *indirectly*, when its charges/bond-orders source a geometric
kernel (charge → EFG `2e`; bond order → McConnell `2e`).

### 1.3 Stream C — AIMNet2

The pretrained AIM vector: a learned, charge-equilibrated per-atom electronic
embedding (`256x0e`-class), plus AIMNet2 charges and the charge-response-gradient
(CRG, labelled as CRG, **not** a polarizability —
`aimnet2_ceiling_ensemble.py`). It is reported as a **learnable ceiling / transfer
diagnostic**, not a recovered law (`project_thesis_reporting_arc` layer 3): a
representation trained for energies/forces that demonstrably carries shielding
signal — including the local-bonding residual the geometric kernels cannot reach
(Cα σ_iso within 0.27→0.60, the atom no through-space kernel touches). AIMNet2 is a
`0e` stream like MOPAC, but a *learned* one.

---

## 2. The equivariant architecture that puts them together

The architecture is **already half-built and proven in the rediscover fitters**;
the spine names it concretely from that code, then states where it grows. The
governing idea, with literature warrant: **scalar electronic streams condition the
combination of geometric tensors; they do not get summed in as fake tensors.** This
is the SEGNN result — "node and edge attributes are not restricted to invariant
scalars but can contain covariant information… the more geometric and physical
quantities are injected the better [it] performs" — with the mechanism being scalar
*gates* and Clebsch–Gordan products interleaved
([Brandstetter et al., SEGNN, ICLR 2022](https://arxiv.org/abs/2110.02905)). e3nn's
own `Gate` is the canonical scalar-gated nonlinearity for mixing irreps while
preserving equivariance ([e3nn](https://docs.e3nn.org/)).

### 2.1 The proven core (what the code already does)

From `equiv_t2_backbone_e3nn.py`, the working heterogeneous equivariant pool:

1. **Shared angular machinery.** `o3.spherical_harmonics(2, disp_hat)` projects each
   source displacement to a `2e` (and, parity-safe, `Y2(axis_hat)` for ring normals
   / bond axes). This `Y2` projection is **shared across all source kinds** — it is
   the geometry, and the geometry is one law.
2. **Per-kind radial channels.** Each source kind (ring / bond / charge) has its own
   invariant MLP `w_kind(invariants_kind)` producing the gate weight on its `2e`
   contribution, fed *only* its own emitted invariants (`ring: r, cosθ, intensity`;
   `bond: r, cosθ, category`; `charge: r, q`). "A ring at 4 Å and a charge at 4 Å
   must be allowed to scale their `l=2` contribution differently" — a single shared
   radial gate would force a law they do not physically obey.
3. **Scatter-pool (Deep Sets).** `pooled.index_add_(group, contrib)` sums per-source
   `2e` contributions into the `(atom, frame)` group accumulator — permutation
   invariance over a variable, heterogeneous neighbourhood.
4. **Per-atom de-meaning.** The per-atom chemical baseline tensor is stripped using
   **train rows only** (`centred_by_train_atom`) — the fit is the through-space
   *modulation*, the within-atom signal the trajectory exists to expose.
5. **One frozen change-of-basis.** `change_of_basis.get_C()` maps the library `2e`
   5-vector into e3nn's `(y,z,x)` 2e convention, pinned, orthogonality-checked at
   startup. `|T2|` is invariant under it.
6. **Honest split.** Blocked temporal frame-split with purged neighbour frames
   (`e3nn_protocol.make_split_masks`); target and feature normalisation centred/scaled
   from train rows only; leave-atoms-out opt-in, effective N printed per stratum,
   thin strata flagged not force-fit.

This is a real equivariant model: `2e` features rotate correctly (Wigner-D tested),
radial weights invariant, prediction de-meaned, target = emitted local-frame DFT T2.
It already fits eight backbone strata.

### 2.2 Where the architecture grows for the full feature set — concrete irreps

The full Stage-3 model is "the exemplar scaled" (`project_stage3_equivariant_gnn`):
an e3nn `radius_graph` message-passing GNN with invariant conditioning and
scatter-pool, **not a giant model**. Concretely:

- **Node irreps.** Per atom, a steerable feature `Nx0e + Mx1o + Kx2e`. The `0e`
  multiplicity is large and carries Streams B and C: MOPAC scalars (charge, s/p
  populations, valency, electron density, AO populations as optional high-dim),
  AIMNet2 embedding (`256x0e`), graph-level conditioning broadcast to nodes (energy,
  IP, dipole magnitude), and the kernels' own `0e` projections (ring iso shift,
  dispersion ⟨E²⟩, H-bond ppm). The `1o` carries the E-field (charge/EFG) and the
  H-bond direction. The `2e` carries the kernel anisotropies (ring CSA, McConnell,
  EFG, π-quadrupole), summed by mechanism but kept as distinguishable channels so the
  fit can weigh them.
- **The combination — gates, then Clebsch–Gordan.** The MOPAC/AIMNet2 `0e` stream
  enters as **steerable node attributes that gate the tensor message-passing**
  (the SEGNN move): the scalar electronic state of an atom modulates *how its
  geometric-tensor neighbourhood is read*. Mechanistically: `0e` features pass
  through `Gate`-style nonlinearities and multiply the `1o`/`2e` channels; the
  message function forms tensor products (`o3.FullyConnectedTensorProduct` /
  Clebsch–Gordan) between the displacement `Y(r̂)` of an edge and the neighbour's
  steerable feature, with weights produced by an MLP of the **invariant** edge and
  node scalars (the radial+electronic conditioning). This is exactly the
  per-kind-radial pattern of §2.1 generalised: the radial MLP's input set grows from
  `(r, cosθ, intensity)` to `(r, cosθ, intensity, charge, s_pop, p_pop, bond_order,
  AIMNet2_scalars…)`. The electronic scalars enter where they belong — **conditioning
  the weights**, not masquerading as tensors.
- **The prediction head, T2 preserved.** The output is `1x0e + 1x2e` per nucleus:
  the isotropic shielding (`0e`) and the symmetric-traceless anisotropy (`2e`), even
  parity. **T2 is preserved end-to-end** (`feedback_t2_sacred`): the rank-2 target is
  the thesis argument, never collapsed to a scalar. The head is the standard
  equivariant tensor head the field uses — decompose to spherical irreps, predict each
  by irrep ([Venetos et al. 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10026072/),
  MAE 1.05 ppm on the full ²⁹Si tensor; [Ben Mahmoud et al. 2024](https://arxiv.org/abs/2412.15063)).

### 2.3 The chewer — the one genuinely new engineering need

Per-source geometry → e3nn tensors → GPU is **the chewer** (parked,
`project_unified_stats_engine`; `project_stage3_equivariant_gnn`). The proven fitters
load emitted CSV/NPY per relationship; the full GNN needs the typed model's per-source
neighbourhood exposed to the GPU at scale (pybind11 binding the spine, not a second
Python model — `feedback_model_is_spine`, `h5-reader/CLAUDE.md`). The chewer and the
720-WT static-pose B-path are the same extension family. **This is the load-bearing
new build; everything else is wiring what we have.**

---

## 3. The two ends

### 3.1 Model 1 — the shielding emulator, anchored on DFT (the stable ground)

Target: the per-frame ORCA r²SCAN/def2-SVP DFT shielding tensor (T0 + T2),
the **last stable ground** (`project_stage3_equivariant_gnn`). The metric against DFT
is honest because ORCA computed σ with no knowledge of our kernels — a non-circular
cross-echo (`doc/emerging/WHAT_ARE_WE_DOING.md` §4). This is the explanation/attribution end: the
partition between kernels matters, the rank-2 content is load-bearing, ablation is
the instrument (funky ones off the top, then everything —
`doc/emerging/STATS_PROGRAM.md`). Reporting is per-element **and** per-IUPAC, always, even at
small n (Stage 1's lesson: backbone N R²≈0.387 vs sidechain N R²≈0.887 — pooling
hides the finding). 1P9J is **within-instrument** (transferability deferred to
720-WT); the between-atom axis was retracted as a frame/centering artifact and is
honoured as deferred, not as loss.

The discipline that makes this end honest is already codified and must carry into the
full model: **form-recovery against our own kernel is tautological** (we fed the
gradient field in); **only the match to independent DFT is real**
(`doc/emerging/GEOMETRIC_KERNEL_ACCOUNT.md` §7). The form-ablation — specific physical kernel vs a
generic same-order `l=2` basis, both against DFT — is the test that converts "predicts
well" into "the mechanism is visible," and it is named-but-not-yet-run
(`doc/emerging/WHAT_ARE_WE_DOING.md` §5).

### 3.2 Model 2 — the shift predictor against experiment ("fido & biscuits")

Target: experimental BMRB shift. The **ethos flips** (`project_stage3_equivariant_gnn`):
prediction not explanation, **R² is the metric**, take the boost understood-or-not,
overlap does not matter, a useful `0e` scalar earns its seat regardless of
attribution. This is where the carry-along scalars (dispersion ⟨E²⟩, H-bond ppm) and
the whole MOPAC/AIMNet2 `0e` pile ride freely. What stays: e3nn (never a hand-rolled
second model), T2 preserved (predict the tensor), held-out honesty, correlate-not-match
on validation (no fabrication), lean.

The honesty that governs this end: the BMRB shift is an **ensemble average over states
our 15 ns MD never samples** — comparing short-MD-averaged predicted shifts to it is
"shooting in the dark," informative (movement/charge/predicted-shift trends vs
one-or-more BMRB shifts) but **not a clean validation**. Interpolate-to-validate
(predict held-out frames) yes; interpolate-to-train (fabricate labels) never.

**The bridge between the two ends is the inclusion test, and it is the reason the
feature set is a superset, not a committee vote:** a feature confounded or redundant
for *attribution* (Model 1) can still ride as a *scalar* in Model 2 if it helps
there. The model reads the superset and chooses; the old fields stay for cross-check
(`pipeline_adaptation.md`, deprecate-and-add).

---

## 4. The ambiguity triage — held, then called cruft or realism, with the reason

This is the heart. Each ambiguity is stated, then triaged. **CRUFT** = accidental
complexity a decision can clear; the spine gets simpler. **REALISM** = a genuine open
question the pipeline is built to honour; it becomes a disclosed open seam, not a
papered-over one. The discipline cuts both ways: do not collapse a realism ambiguity
into a tidy label (the tool's named failure mode — `CONTINUITY.md`), and do not
dignify a bookkeeping double-count as deep physics.

### 4.1 Descriptor-vs-mechanism (the kernel is a field AND a geometric descriptor)
**REALISM.** The companion proves these are the *same object* — the low-order field of
a source distribution *is* its low-order geometric moment
(`doc/emerging/GEOMETRIC_KERNEL_ACCOUNT.md` §5), so it is not a choice to be resolved. What is open
is the **severity**: the excess of the specific kernel's echo over a generic
same-order feature, which is named and motivated but not yet a measured number
(`doc/emerging/WHAT_ARE_WE_DOING.md` §4). The pipeline honours it by running the form-ablation and
by reporting `|T2| r` as the tensor headline with scalar R² as a diagnostic — never
asserting "mechanism visible" from correlation alone. **Why realism, not cruft:** the
field uses the same object in both registers and labels which register it is in
(`doc/emerging/STAIRCASE_SOCIAL_HISTORY.md`); collapsing it would be exactly the over-tidy move the
guardrail forbids.

### 4.2 The round-trip circularity (we feed the gradient field in)
**REALISM — but with a crisp cruft-removal inside it.** The realism: form-recovery is
partly tautological; only the DFT-target match is non-circular. That is permanent and
disclosed. The cruft *inside* it: any pipeline path that scores a model against the
**producer's own bare kernel** and headlines it (PySR recovering the Pople form from
`bare_T0`, which already holds Giessner-Prettre constants) is an accidental
self-echo — already caught and de-headlined ("Do NOT headline it" — `FINDINGS.md`,
`doc/emerging/WHAT_ARE_WE_DOING.md` §3). **Action:** keep the circularity disclosure (realism);
keep every distillation-against-our-own-kernel explicitly labelled a pipeline check
and out of the headline (cruft, swept). The architecture already enforces the clean
side: the broad/ring fitters ingest **raw geometry** and let e3nn build `l=2`, so the
model is given coordinates and discovers the law, not handed its value.

### 4.3 The partition web (which feature overlaps are real)
**MIXED — triage each edge.** Overlap is fine for the predictor (Model 2), a confound
for attribution (Model 1). The edges, called individually:

- **bond-anisotropy ↔ McConnell:** **CRUFT.** They are the *same machine* over the
  same bonds (`bond_anisotropy.md` §2/§4). Resolution is decided: **extend the one
  McConnell calculator** (richer per-bond Δχ, rhombic C=O/C=C, the S–S/C–C categories),
  do not stand up a second kernel. Swept.
- **aromatic Δχ ↔ ring current:** **CRUFT.** Aromatic susceptibility anisotropy *is*
  the π ring current re-expressed; zero aromatic bonds while the ring kernel is active
  (`mcconnell.md` §1.6). Decided, swept.
- **π-quadrupole ↔ charge/EFG:** **REALISM until a partition decision is made, then
  CRUFT.** If charge/EFG already sums the ring's partial atomic charges, the
  quadrupole is the same charge re-expanded — a double-count. The term carries
  independent signal *only* under an explicit charge-partition (delegate ring atoms to
  the quadrupole, OR source Θ from a better QM moment). This is "the single most
  important design decision for the kernel" (`pi_quadrupole.md` §3.1) and it is **not
  yet settled in the code**. **Action:** make the decision explicit at the composition
  (the new MOPAC-extended calculator's shared substrate is exactly where "each
  partition is decided once" — `CONTINUITY.md`). Realism now because undecided; becomes
  cruft the moment it is decided.
- **dispersion ↔ ring / charge / H-bond:** **REALISM, leaning drop-candidate.** The
  R⁻⁶ contact scalar is explicitly hard to separate from anisotropy, electrostatics,
  and H-bonding, and has **no DFT anchor** (the campaign wrote total σ only —
  `dispersion.md` §1.5). It rides as a labelled `0e` channel the fit may reject; "may
  not earn its place" is an acceptable, designed-for outcome. Honour it as an open
  question; the sweep decides, not us.
- **H-bond ↔ charge/EFG:** **REALISM (genuinely partial).** The long-range
  electrostatic part of the H-bond effect *is* inside the charge/EFG field already; the
  H-bond channel's distinct content is the short-range **covalent/Pauli** depletion the
  point-charge field cannot represent (the proton sits *inside* the polarised density —
  `hbond_larsen.md` §1c). No clean numeric partition exists at this scale. **Action:**
  emit both, disclose the overlap at the feature boundary, let the fit weigh the shared
  vs unique content. Do not claim a partition.
- **MOPAC-charge ↔ AIMNet2-charge ↔ ff14SB-charge:** **REALISM as parallel channels;
  no partition needed.** These are three estimates of the same physical field; the
  defensible move is to emit all and let the model weigh them (`charge_efg.md` §3) —
  not to pick one. Realism that the architecture *absorbs* rather than resolves.
- **MOPAC-bond-order McConnell ↔ literature-Δχ McConnell:** **REALISM as channels.**
  Conformation-responsive Wiberg weighting beside the fixed-literature Δχ weighting —
  channels the fit weighs, overlapping for attribution. Disclose, do not collapse.

### 4.4 MOPAC ↔ geometry / AIMNet2 redundancy
**REALISM — and the spine's most important open seam to name precisely.** Three
distinct questions, often conflated:
- **MOPAC charges vs ff14SB/AIMNet2 charges** → §4.3, parallel channels, absorbed.
- **MOPAC electronic descriptors (s/p populations, AO/overlap populations, valency,
  IP) vs AIMNet2's learned embedding.** This is the genuinely open one. AIMNet2's AIM
  vector is *itself* an iterated, charge-equilibrated electronic embedding that
  "incorporate[s] rich electronic structure information" and refines the local
  chemical environment each message pass
  ([AIMNet2, Chem. Sci. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12057637/)). So
  MOPAC's explicit PM7 populations and AIMNet2's learned scalars may cover overlapping
  ground. **But they differ in kind:** MOPAC's are *named, interpretable, physically
  grounded* quantities (an s-population is an s-population); AIMNet2's are *learned,
  opaque, higher-capacity*. For Model 1 (attribution) the named MOPAC scalars are the
  ones you can *say something about*; for Model 2 (prediction) AIMNet2's capacity may
  dominate. The inventory the project already queues — *"what does MOPAC tell us that
  our geometry and work tell us differently?"* run **after** the full MOPAC capture
  lands (`CONTINUITY.md`) — is exactly the instrument for this seam. **Action:** honour
  as open; the inventory is upstream-of-build (which source feeds which calculation),
  distinct from the downstream partition. Do not assume AIMNet2 subsumes MOPAC or
  vice-versa; the sweep ablates both.
- **MOPAC geometry vs ours** → **CRUFT, decided.** Geometry is ours (MOPAC doesn't
  frame positions/neighbourhoods); through-space kernels are ours (MOPAC doesn't
  compute a ring current); local-electronic is MOPAC's. The default-source rule (MOPAC
  default, override only where ours is *shown* better) resolves the clean cases by
  construction (`CONTINUITY.md`).

### 4.5 The equivariant-architecture choices
- **Scalars as gates/conditioning vs scalars as concatenated channels.** **REALISM
  leaning decided.** The SEGNN literature says inject physical/geometric quantities as
  steerable conditioning and they help; e3nn's `Gate` is the standard mechanism. The
  spine adopts conditioning (§2.2) as the principled default, but whether the gain over
  plain `0e` concatenation is structural for *our* streams is for the sweep. Honoured,
  not forced.
- **Per-kind radial channels vs one shared radial law.** **CRUFT, decided** by the
  working code: heterogeneous sources must scale their `2e` differently
  (`equiv_t2_backbone_e3nn.py`); a shared gate forces a law they do not obey.
- **`2e`-only vs `1o ⊕ 2e` for the electrostatic family.** **REALISM, emit both.** The
  field (`1o`) and gradient (`2e`) are different multipole orders, not derivable from
  each other (different `l`, different distance law — `charge_efg.md` §2). Both earn a
  place; the fit weighs.
- **Frame/parity correctness.** **CRUFT (a correctness obligation, not an ambiguity).**
  Parity must be declared right or equivariance is silently wrong; raw molecular frame,
  no imposed local frame. Get it right and state it — already the discipline.

### 4.6 Which scalars belong / which MOPAC quantities to use
Triaged as the delta — see §5.2 and §5.3. Short form: the `0e` MOPAC pile **belongs**
because it is the complement region (REALISM the architecture is *for*); *which*
specific MOPAC quantities carry independent signal is the open inventory (REALISM,
sweep-decided), with the spine's recommended starting set in §5.3.

### 4.7 The two contexts (shielding-tensor vs shift)
**REALISM, and it is the organising principle, not a problem.** The same feature can be
cruft for attribution and a darling for prediction. The pipeline honours this by being
*two models off one superset*, with the inclusion test as the bridge (§3). Not a thing
to resolve — a thing to build around. Already decided.

---

## 5. The delta — the worth of the pass

A polished restatement is worth nothing here. These are the surprises, additions, and
changes that research or first principles surfaced.

### 5.1 Surprises

1. **The architecture itself changes once the rich `0e` stream is real — and the
   change has a name and a citation.** Before this pass, the implicit picture was
   "concatenate MOPAC/AIMNet2 scalars as extra `0e` channels next to the tensors." The
   SEGNN result reframes it: the scalar electronic stream should **condition the tensor
   message-passing as steerable node attributes** ("the more geometric and physical
   quantities are injected the better" — [Brandstetter et al. 2022](https://arxiv.org/abs/2110.02905)),
   via scalar gates and interleaved Clebsch–Gordan products. The mechanistic upshot is
   concrete (§2.2): the per-kind radial MLP that already gates the `2e` contributions
   **grows its input set to include the electronic scalars**, so an atom's MOPAC/AIMNet2
   state modulates *how its geometric-tensor neighbourhood is read*. This is a small,
   principled change to a model we already have running — not a new architecture — and
   it is the structurally correct place for Stream B to enter. **(Stated structurally,
   no outcome forecast: this is where the electronic scalars *belong* in an equivariant
   model, not a claim they will raise R².)**

2. **MOPAC has essentially no native tensor content — and that is the point, not a
   gap.** Holding the whole feature set at once makes it obvious: Stream B is almost
   entirely `0e`. The instinct to find a "MOPAC `2e`" is a trap. MOPAC's tensorial
   contribution is *indirect and already covered* — its charges source the EFG `2e`,
   its bond orders source the McConnell `2e`. The complement is clean precisely because
   the streams divide by irrep content as well as by physics: tensors from geometry,
   scalars from electrons. This sharpens the architecture (no wasted machinery hunting
   a MOPAC tensor) and the triage (the MOPAC↔geometry overlap is a *scalar-conditions-
   tensor* relationship, not a competing-tensor double-count).

3. **The dispersion kernel's "no DFT anchor" is a structural fact, not a tuning
   problem.** The fixed-DFT campaign wrote total σ only; there is no separable DFT
   dispersion shielding to validate against, ever (`dispersion.md` §1.5). This is a
   harder constraint than the other carry-alongs face and it changes how dispersion is
   honoured: it can *only* be a geometric surrogate offered to the fit, never a
   pointwise-validated readout. It belongs in Model 2 (prediction) far more
   comfortably than in Model 1 (attribution).

4. **AIMNet2 is not "just another charge source" — it is a learned electronic
   embedding that overlaps MOPAC's descriptor role.** The pass surfaced that AIMNet2's
   AIM vector is an iterated charge-equilibrated electronic representation
   ([AIMNet2 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12057637/)), so the real
   redundancy seam is **MOPAC-populations ↔ AIMNet2-embedding**, not MOPAC-charge ↔
   ff14SB-charge. The two differ in *kind* (named/interpretable vs learned/opaque),
   which maps cleanly onto the two ends (MOPAC for Model-1 attribution, AIMNet2 capacity
   for Model-2 prediction). This is the redundancy to instrument, and it has been
   mis-slotted as a charge question.

### 5.2 Additional feature categories we should carry

1. **Graph-level MOPAC scalars as broadcast conditioning.** Heat of formation,
   electronic/core-core energies, **ionization potential / frontier-orbital summary**,
   and the full molecular dipole are in the capture (`mopac_extension.md`) but not
   slotted as model inputs. Broadcast to nodes as `0e` conditioning, they are cheap
   global context (the IP is a molecule-level frontier-orbital read; the dipole
   magnitude a global polarity). They are conditioning metadata, not per-atom physics —
   carry them as such.

2. **AO-resolved and overlap populations as optional high-dimensional local
   descriptors.** The capture preserves per-AO populations and AO/AO overlap populations
   (`mopac_extension.md`). These are the finest-grained local-electronic descriptors
   available and are a natural high-dim `0e` block for the Model-2 (prediction) end,
   where capacity is welcome. Carry as an optional channel the sweep can switch on.

3. **The McConnell rhombic component for C=O / C=C.** The current McConnell path is
   axial-by-construction; the literature tabulates a genuine *two-component* (axial +
   rhombic) Δχ for the carbonyl and the C=C double bond (`bond_anisotropy.md` §1.3,
   ApSimon/Abraham). Carrying the rhombic part where it is tabulated is more `2e`
   structure for the bonds where it physically exists — an addition the fold-into-
   McConnell already scopes.

### 5.3 Additional scalars — and the MOPAC quantities to start using, with why each
carries signal

The capture print (`mopac_extension.md`) is vastly richer than the four legacy arrays
(`mopac_charges`, `mopac_scalars`, `mopac_bond_orders`, `mopac_global`). The
quantities to **start using**, each with the structural reason it carries signal the
others cannot:

| Quantity | Form | Why it carries signal the geometric kernels structurally cannot |
|---|---|---|
| **s / p shell populations** | `0e` per atom | Direct read on local hybridisation — the quantity the paramagnetic σ term turns on. The through-space `l≤2` kernels are *blind* to it (`doc/emerging/GEOMETRIC_KERNEL_ACCOUNT.md` §6); this is the heavy-atom complement, the Cα/backbone-N gap the tensors leave on the table. |
| **AO populations + overlap populations** | `0e` per atom / per AO-pair | Finest-grained local electronic structure; an overlap population is a direct bonding-character measure between two AOs — local-quantum information no neighbourhood multipole encodes. |
| **MOPAC diagonal valency** | `0e` per atom | MOPAC's own valency (distinct from derived) — a conformation-responsive local-bonding scalar. |
| **Wiberg bond orders (full directed object)** | bond/pair axis | The conformation-responsive McConnell source weight; the directed object retains valency diagonals and H-bonds (`ALLBONDS`) the legacy symmetric projection drops. |
| **Ionization potential / frontier-orbital summary** | `0e` graph-level | Frontier-orbital energy is the molecule-level handle on the excited-state manifold the paramagnetic term sums over; a global conditioning scalar for the local-quantum coverage region. **Honest limit:** the search this pass did not surface a *direct* IP↔shielding correlation in the literature — it is a structurally-motivated conditioning scalar, sweep-decided, not a literature-validated predictor. Have-it-to-cite-it: carried as a hypothesis channel, not a claimed law. |
| **Per-atom electron density / dipole-contribution** | `0e` per atom | Local charge environment beyond the single Mulliken number. |
| **Full molecular dipole (POINT-CHG / HYBRID / SUM)** | `1o`/`0e` graph-level | The HYBRID/SUM split carries more than the legacy SUM vector; global polarity conditioning. |

The honest framing on *why* this pile is worth the capture cost: MOPAC "tested
strong" (`CONTINUITY.md`), and the plausible-but-held account is that the geometric
kernels are structurally blind where shielding is local-electronic-dominated, and
MOPAC's descriptors are a cheap proxy for exactly that. **The test is whether MOPAC's
lift concentrates on the heavy, local-dominated atoms — verify it, do not claim it**
(`CONTINUITY.md`). The literature corroborates the *direction*: PM7-geometry-based
descriptors are used for near-DFT NMR prediction at scale, and the standard descriptor
sets include hybridisation, bond counts, and aromaticity — the local-bonding flavour
([arXiv:2510.05623](https://arxiv.org/abs/2510.05623); search-summary grounded,
WebFetch denied this pass — stated per have-it-to-cite-it).

### 5.4 Does the equivariant architecture change once the rich scalar stream sits
beside the geometric tensors?

**Yes, structurally — and the change is the §5.1.1 surprise made into the build
instruction.** Without Stream B, the model is a tensor pool with invariant radial
gates (the proven fitter). With a rich `0e` electronic stream, the principled
architecture is **SEGNN-style**: the electronic scalars become steerable node
attributes that condition the tensor message-passing (gates + Clebsch–Gordan), and
the per-kind radial MLP's input set grows to include them. This is not a different
model family — it is the same e3nn scatter-pool/message-passing GNN with a richer
conditioning input — but it *is* a real architectural commitment: the electronic
stream enters as **conditioning on the geometric combination**, not as side-channel
scalars bolted on at the head. The prediction head and the T2 preservation are
unchanged. **(Structural statement only: this is where the streams belong in an
equivariant model; whether it lifts R² is the sweep's to say.)**

---

## 6. The connected flow, end to end (the spine, in one pass)

Extraction (the C++ library, frozen during reader work; producer changes
user-scheduled): the typed protein model computes, per frame, the through-space
geometric kernels as clean `0e`/`1o`/`2e` spherical tensors in the molecular frame
(ring, McConnell-with-bond-anisotropy, charge/EFG, plus the carry-alongs), and runs
**MOPAC per-frame on everything** — held whole in the extended `MopacResult` buffer —
emitting the full local-electronic capture (charges, s/p and AO/overlap populations,
Wiberg bond orders, valency, MO/LMO, energies, IP, dipole) plus AIMNet2. The
default-source rule keeps the split clean: geometry and through-space are ours,
local-electronic is MOPAC's, charges are parallel channels. The emit is
deprecate-and-add: legacy arrays stay, the canonical irrep features and the full MOPAC
sidecar are added (`pipeline_adaptation.md`).

Reader (the trust/vetting boundary; rides the same typed model): MOPAC's per-bond and
per-atom facts propagate to `QtAtom`/`QtBond` so "what MOPAC knows" shows on selection,
as time-series over the trajectory. The reader builds the equivariant feature
substrate from the resident body and emits it as per-relationship CSV + NPYs, read
through the resident `Catalog` (`valueVec3`/`valueT2`/`valueTensor` already carry
`1o`/`2e`/full). Python consumes the substrate only — never `trajectory.h5`
(`feedback_no_parallel_h5_in_python`).

Model (the far end): an e3nn message-passing GNN, scatter-pooled, T2-preserving. Node
irreps `Nx0e + Mx1o + Kx2e` carry the three streams: the kernel tensors (`1o`/`2e`),
the MOPAC/AIMNet2 electronic scalars (`0e`), the graph-level conditioning. The
electronic scalars **condition** the tensor combination (SEGNN gates + CG products);
the per-kind radial MLPs grow their input set to include them. Two ends off one
superset: **Model 1** anchored on DFT σ (T0+T2), within-instrument on 1P9J,
transferability to 720-WT, attribution-disciplined (ablation, form-recovery vs DFT-match
kept apart, per-element/per-IUPAC always); **Model 2** against BMRB shift, prediction
ethos, R²-the-metric, the carry-along scalars riding free, validation honest about the
MD-ensemble gap. The chewer (per-source-geometry → e3nn → GPU, pybind11 on the spine) is
the one new build that makes the full model run at scale.

That is the spine. The settled parts are settled; the realism seams — descriptor-vs-
mechanism severity, the π-quadrupole/charge partition decision, the H-bond covalent
overlap, the MOPAC-populations/AIMNet2-embedding redundancy, dispersion's no-anchor
status — are the honest open questions the pipeline is **built to honour**, not hide;
the cruft seams — bond-anisotropy↔McConnell, aromatic↔ring, MOPAC↔geometry,
per-kind-radial, the self-echo headline — are swept by named decisions, and the spine
is simpler for it.

---

### Sources (this pass; WebFetch denied — grounded on search summaries + held corpus,
stated per have-it-to-cite-it)
- [Brandstetter et al., SEGNN — Geometric and Physical Quantities improve E(3) Equivariant Message Passing, ICLR 2022](https://arxiv.org/abs/2110.02905) — steerable node attributes; scalar gates + Clebsch–Gordan; "the more geometric and physical quantities are injected the better."
- [Venetos et al., ML full ²⁹Si NMR shift tensors with equivariant GNNs, JPCA 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10026072/) — symmetric-spherical-tensor target, MAE 1.05 ppm; tensorial handling is the key.
- [Ben Mahmoud et al., GNN solid-state NMR parameters via spherical-tensor decomposition, arXiv:2412.15063](https://arxiv.org/abs/2412.15063) — decompose to irreps, predict each by irrep.
- [AIMNet2, Chem. Sci. 2025 / PMC12057637](https://pmc.ncbi.nlm.nih.gov/articles/PMC12057637/) — learned AIM vector, neural charge equilibration, iterated electronic embedding.
- [Neighborhood-Informed Representations for AIM NMR shielding ML, arXiv:2510.05623](https://arxiv.org/abs/2510.05623) — neighbourhood electronic/structural enrichment of shielding ML (search-summary grounded).
- [e3nn docs — irreps, Gate, tensor products](https://docs.e3nn.org/) — the scalar-gated equivariant combination mechanism.
- Held corpus / project docs: the `kernel_design/*.md` set, `doc/emerging/GEOMETRIC_KERNEL_ACCOUNT.md`, `doc/emerging/WHAT_ARE_WE_DOING.md`, `doc/emerging/GEOMETRIC_KERNEL_MATH_LINEAGE.md` + `doc/emerging/STAIRCASE_SOCIAL_HISTORY.md`, `doc/emerging/STATS_PROGRAM.md`, and the rediscover fitter code (`equiv_t2_backbone_e3nn.py`, `equiv_t2_efg_e3nn.py`, `aimnet2_ceiling_ensemble.py`, `e3nn_protocol.py`).
