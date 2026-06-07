# MatTen equivariant shielding predictor — adaptation & implementation plan

Status: reviewed (2026-06-07). Scope decided in session; codex code-verified
review folded in — see `MATTEN_CODEX_REVIEW_2026-06-07.md`. Net of the review:
the two reader emits and the equivariance framing stand; the **reader-path
element check is reinstated as a required guard** (it does not exist today, only
the producer-side path has it); several precision fixes landed (row≠molecule
regrouping, exact reflection rule, hidden-lmax as a sensitivity choice).
**Correction to codex (Jessica, 2026-06-07): the 720 DFT targets DO exist — they
are the mutant set** (the WT/ALA poses that calibrated the Stage-1 ridge; codex
hit a location/naming miss, not absence). So the transferability pilot is
**data-backed**; the remaining task is to **locate and characterize** that DFT,
not to establish it exists (§7.1).

## 0. Decision and scope

We adopt **MatTen** (e3nn-based equivariant GNN, `github.com/wengroup/matten`)
as the architecture for an equivariant per-atom NMR shielding-tensor
predictor. It is the safest defensible choice: a published, peer-reviewed,
maintained e3nn model whose irreducible-tensor target machinery already
predicts full NMR shielding tensors (Venetos et al., silicon oxides,
*J. Phys. Chem. A* 2023 / PMC10026072), parity-correct, and capable of the
full `0e ⊕ 1e ⊕ 2e` target including the antisymmetric part.

**There is exactly one version: v1.** v1 predicts the DFT shielding tensor
from geometry alone — per-atom Cartesian coordinates + atomic number → the
per-atom shielding tensor. We do **not** feed the calculator kernels as input
features. A kernel-conditioned variant ("v2") is explicitly out of scope and
considered indefensible: the kernels *are* the low-order (ℓ≤2) multipole
moments of the weighted neighbourhood, and an ℓ≤2 steerable network computes
the same projections internally — feeding the kernels in would let the model
relay them and would destroy the one thing the predictor exists to show
(whether geometry carries the signal). Geometry-only is both the cleaner
equivariance choice and the only honest experiment.

### What v1 is

v1 is the **equivariant transferability pilot** of the thesis arc: an equivariant
network, given only geometry, reproduces DFT shielding tensors (including the
anisotropy T2 and the antisymmetric part T1) on **held-out proteins** — trained
**leave-proteins-out on the 720-WT static poses**, with the Stage-1 ridge as the
invariant baseline. This is the thesis arc's "720 WT backbones, statics
baseline."

This is **data-backed**: the full-tensor per-atom DFT targets exist as the
**mutant set** — the WT/ALA poses that calibrated the Stage-1 720 ridge (§7.1).
The remaining work is to **locate and characterize** that DFT (coverage, confirm
absolute-σ), not to establish that it exists. 1P9J's 751 frames remain available
as a within-protein / conformational companion, but the headline is the
transferable pilot.

### What v1 is not

- It does **not** predict experimental chemical shift. It predicts and is
  validated against **DFT shielding**. Experimental solution shift is the
  ensemble average ⟨σ⟩ that a 15 ns trajectory never samples — the dragon —
  and is out of v1's claim entirely.
- It is **not** a calibration. It is a model (signal-capture / R²), not a
  fittable law.

## 1. Why MatTen over the alternatives

- **e3nn / MatTen** — output irreps are declared explicitly, so the full
  `0e ⊕ 1e ⊕ 2e` parity-correct target is a config line; the tensor-readout
  head is published and NMR-validated. Heavier CG/spherical machinery, but
  that cost is the library's, not ours.
- **TensorNet** (Cartesian rank-2) — best object-model match (its hidden
  state *is* a general 3×3 = our T0/T1/T2) and cheaper, but the released
  TorchMD-Net ships no rank-2 *output head* and no arbitrary-feature *input
  port*; we would build both ends. A "TensorNet-MatTen" hybrid is a category
  error — they are alternative backbones, not composable.
- **PaiNN + tensor channels** — symmetric-rank-2 only (T0+T2, **no T1**),
  parity implicit/polar, and the tensor-channel extension ships no code. It
  silently drops the antisymmetric part we refuse to drop.

MatTen is the one route that is (a) proven on real DFT NMR tensors, (b)
parity-correct for a magnetic (axial) response tensor, (c) T1-capable by
configuration, and (d) reusable rather than green-field.

### Parity note (why "never 1o")

Shielding relates an induced magnetic field to an applied one — a response of
an axial vector to an axial vector — so σ is parity-**even**. Its irreducible
pieces are `0e ⊕ 1e ⊕ 2e`: isotropic 0e, antisymmetric **1e** (axial
pseudovector, *not* 1o), symmetric-traceless 2e. Input geometry is polar
(positions are `1o`, `Y₁` is `1o`, `Y₂` is `2e`); the even-parity target is
reachable from polar geometry via even products (`1o ⊗ 1o → 0e ⊕ 1e ⊕ 2e`).
e3nn enforces this **iff** the irreps are typed correctly — one wrong label is
a silent correctness bug.

Two concrete guards (both from the code review):
- **The adapter must declare the output irreps explicitly** and not infer parity
  from existing metadata. There is a stale-`1o` hazard: older docs/comments still
  imply `1o` for shielding T1 (`python/doc/index.rst`,
  `QtTimeSeriesBuffers.h`, an older codex brief), and some Vec3 catalog metadata
  pairs `irreps="1e"` with `parity="odd"`. The live contract is right
  (`_SHIELD_IRREPS = "0e + 1e + 2e"`, `SphericalBasis.cpp:13-17` names T1 the
  antisymmetric pseudovector), but the adapter must not pick up the stale labels.
- **Reflection test, stated exactly:** an axial T1 pseudovector does **not** flip
  under inversion; under a general improper rotation it transforms as
  `v' = det(R) · R · v` (e.g. reflecting x→−x acts by `diag(1,−1,−1)`). The test
  must encode that exact rule, plus a rotation-equivariance check, not a vague
  "pseudo-components flip."

## 1b. Assumptions from Venetos/MatTen we do NOT inherit

Venetos et al. (the precedent, and MatTen's home paper) is silicon-oxide
*crystals*. Read in full from the corpus (2026-06-07,
`references/venetos-2023-ml-nmr-shift-tensors-silicon-oxides-egnn.pdf`). Four
assumptions are baked in that are invalid here. The **model (equivariance)
transfers; these do not.**

1. **Symmetric-only target.** Their best config drops the antisymmetric part:
   *"A symmetric spherical tensor target was found to yield the best loss"* (p.2393);
   *"only the symmetric part of the tensor is used"* (p.2389). That is our T1.
   MatTen off-the-shelf = symmetric-only → silently drops T1. We **force output
   irreps `0e+1e+2e` and keep T1** (overriding their default, not adopting it).

2. **Point-group site symmetry drives their accuracy — proteins have none.**
   *"the shift tensor is very closely linked to the structural point group… grouped
   by Qⁿ into Td, C3, C2"* (p.2393); accuracy tracks symmetry (Q4/Td 1.52 ppm vs
   Q1 8.58 ppm, Table 4). Every protein atom is ~C1 — our regime is their *worst*.
   Distinction to hold: **model equivariance transfers (their 53% win over
   invariant models); data site-symmetry does not (their 1.05 ppm absolute).** Do
   NOT quote their accuracy as a target — the honest bar is "does equivariance beat
   the ridge," not "1 ppm."

3. **Crystal graph + periodic boundary conditions** (p.2391) — but PBC was doing
   *real maths work*, not plumbing: it guaranteed **environment completeness** (no
   surfaces; σ is a complete function of the periodic environment the DFT also
   saw). That is the actual "can the model see the signal" argument — equivariance
   only fixes the hypothesis class; completeness is what makes the target a
   function of the input. **We get completeness by a different route —
   finite-cluster consistency** — and *more cleanly*: a crystal's infinite lattice
   sum is always truncated at `rcut`, while a finite molecule fed in full has no
   tail to truncate.
   Verified DFT convention (mutant set, e.g. `A0A822ILJ2_WT_nmr.out`):
   `! r2SCAN def2-SVP def2/J NMR CPCM(Water)`, **543 atoms = whole domain, all
   atoms shielded, implicit CPCM solvent, no explicit waters, no fragment caps.**
   So completeness holds: we feed exactly the atom set the DFT saw. What remains of
   the PBC worry is two concrete items, not a blocker:
   - **Receptive field ≥ molecule extent.** `n_layers × rcut` must span the domain
     (~30–40 Å) to capture long-range ring-current/field contributions. `rcut` is
     therefore a *signal-visibility* parameter, not just a physics-range knob.
   - **CPCM is the one non-local piece of σ.** The implicit-solvent reaction field
     depends on the whole molecular cavity/surface — deterministic from geometry
     (so completeness and equivariance hold; isotropic dielectric) but hardest for
     a *local* GNN, most so on solvent-exposed atoms. Prime suspect if surface σ
     predicts poorly; disclose as a known limit.
   Build the radial graph **without periodic images** (finite molecule).

4. **Single target nucleus (²⁹Si), Si-only readout** (p.2391), predicting absolute
   *shift* δ via a reference. We have **H/C/N/O all as targets** → per-element
   readout/standardisation, and we stay in absolute *shielding* σ (no referencing).

Net: MatTen is the right model; we consciously reject these four and frame the
thesis around the harder, symmetry-free protein regime.

## 2. Architecture: the bridge

The data path is **not** "Python reads NPYs." It is:

```
trajectory.h5 + 5-NPY topology sidecar + per-frame ORCA .out
      │   (all named by the .LGS manifest — no file discovery)
      ▼
h5-reader CLI prep path  (rediscover, headless: h5reader_extract)
      │   NPY/H5 → fully-resident typed C++ model (QtProtein spine)
      │   + indexes (per-frame KD-trees, frame/DFT-row maps, typed atom index)
      │   + functional surface (Catalog.value / verbs / RunTraversal)
      ▼
lean substrate  (row-keyed CSV + NPY sidecars, lab frame, <15 GB)
      ▼
MatTen adapter  (new Python consumer) → pymatgen Molecule + per-atom target
      ▼
MatTen / e3nn  (GPU)
```

This honours the standing discipline: the model is the spine (C++ holds the
gigabytes and the indexes), Python gets a lean pre-digested input, no parallel
H5 in Python, lean disk / generous RAM.

Key code anchors (read-only review, 2026-06-07):

- CLI entry: `h5-reader/src/rediscover/main_extract.cpp:145`
  (binary `h5reader_extract`), the `all_atom_equivariant` case at
  `:438`/`:454`.
- Immutable carrier `RunData` (`RunData.h:125`), functional `Body{run, idx,
  catalog}` (`AnalysisBody.h:11`), `Catalog` (`Catalog.h:88`),
  `ResidentIndexes` (per-frame nanoflann KD-trees, frame/DFT-row maps,
  chemistry-typed atom index) (`ResidentIndexes.h:14`).
- ORCA target: the reader **re-parses the raw ORCA `.out`** itself
  (`io/OrcaShieldingParser.cpp:46`), deliberately keeping more than the
  producer NPY (dia/para split, verbatim asymmetric 3×3, ORCA-input coords).
  Target built in `ExtractionSupport.cpp:45`; T1 computed in
  `SphericalBasis.cpp:15-17`.
- The local-frame target (`*_target_local_T2`) lives in a **different** sink
  (`BroadBackboneSink.cpp:198`, the scalar/ridge path), correctly quarantined
  away from the equivariant sink.

## 3. Frame and atom-order consistency — settled

The per-frame **DFT PDB was emitted by the extractor during trajectory
traversal**, and ORCA ran on exactly that geometry. So:

- the H5 node positions and the ORCA tensor are the **same coordinates in the
  same lab frame** by construction, and
- ORCA nucleus order == topology atom order by construction.

Consistency has been tested for 1P9J, so for the **existing 1P9J campaign** the
invariant holds and neither guard is on the critical path.

But the review corrected an over-generous downgrade: the **reader DFT path does
not actually check element identity**. `DftShieldingLoader` validates atom count,
parser holes, and `total == dia + para`
(`io/DftShieldingLoader.cpp:78-105`) but never compares
`DftAtomShielding.element` to the topology atom. The hard element check exists
only on the **producer/static** side (`nmr-720/src/OrcaShieldingResult.cpp:181-204`)
— which is probably why the by-construction story feels safe. So:

- **Reinstate the element check in the reader path as a required fail-loud guard**
  before any 720-WT or re-prepped trajectory data is trusted. It is cheap and
  closes the only way a tensor lands on the wrong atom undetected (a
  same-cardinality permutation passes the count check silently).
- **Kabsch alignment:** keep `CheckDftFrameAlignment` (`ExtractionSupport.cpp:101-162`)
  as diagnostic-only for 1P9J if the audited angle/RMSD stats are near zero; for
  future campaigns, record a threshold in the substrate build log / coverage
  manifest and reject above it.

## 4. Reader-side changes (the only code we add to the producer/prep path)

Two additive emits in the equivariant sink (`AllAtomEquivariantSink`). Both
read data already resident in RAM; neither is a redesign, neither touches the
extractor library, neither touches `MopacResult`/analysis (the lead's
concurrent-edit surface — confirmed independent of the ORCA-target path).

1. **Raw per-atom point cloud.** Emit `positions.npy` (lab-frame Cartesian,
   `(rows, 3)`) keyed to the same `row_id` as the targets, sourced from the
   resident `Conformation::atomPosition` / `QtPositionsTimeSeries`. This is the
   central gap: the current substrate is *relational* (per-source displacement
   vectors under mechanism-specific cutoffs, rings/bonds as aggregated
   entities), which suits the classical sum-pool consumer but is the wrong
   shape for MatTen. MatTen wants a homogeneous atom graph and builds its own
   radial edges from node positions + Z. `Z` is already on every target row
   (`element`).

   *Optional:* atom–atom edges could be dumped from the existing
   `CloudKind::Atoms` KD-tree, but the simpler and more standard choice is to
   emit positions only and let MatTen construct the radial graph.

2. **`target_T1.npy`** (`(rows, 3)`), from the already-computed
   `total_decomp.T1`, plus `dft_total_T1_{0,1,2}` CSV columns. Today the sink
   writes `target_raw` (3×3), `target_T2` (5), `target_sigma_iso` (T0) but
   **drops T1** — so a consumer training on the named irrep sidecars silently
   gets a symmetric-only target. T1 is recoverable from `target_raw`, but it
   must be surfaced as a first-class irrep array so the full-irrep target is
   the obvious one to use, not the easy-to-miss one. (Full-irrep discipline:
   never collapse to symmetric-only.)

**Presence gate:** today the gate is split — `dftPresent` is
`target.present && finiteT2(...)` (`AllAtomEquivariantSink.cpp:183-187`) but
`sigma_iso`/`raw` are appended under only `target.present`
(`:210-215`, `:246-250`). Unify all of `target_raw`, `target_T2`,
`target_sigma_iso`, `target_T1` on one predicate —
`target.present && finite(total_raw) && finite(T0/T1/T2)` — so a row is uniformly
usable or uniformly NaN.

## 5. The MatTen adapter (new Python consumer)

A new module reading the substrate (NPYs + the row-keyed CSV — never the H5):

- **Row ≠ molecule.** The traversal is DFT-rows-outer, atoms-inner
  (`RelationshipEngine.h:79-86`; all-atom stratum is every atom,
  `AllAtomEquivariant.cpp:220-223`), so a `row_id` is one **(atom, frame)**
  target row, not one molecule. The adapter must **group rows by `h5_row` /
  `original_index`, sort by `atom_index`, and assemble one whole molecule per
  frame**. Assert this grouping invariant in the adapter tests — it is the most
  likely place to silently scramble a structure.
- Per frame: assemble `positions[N,3]`, `Z[N]`, and the per-atom target
  tensor. Reconstruct the full 3×3 from `target_raw` (or from
  T0 ⊕ T1 ⊕ T2) — **keep T1**.
- Build a non-periodic pymatgen `Molecule` (or a large-vacuum-box `Structure`)
  per frame/pose — MatTen's data layer is crystal/PBC-oriented, so we feed it
  finite systems.
- Configure MatTen output irreps to `0e + 1e + 2e`. Predict the per-atom
  tensor directly (no readout-only outer-product shortcut).

Precedent in-tree: `analysis/equiv_t2_e3nn.py` is a bespoke e3nn sum-pool over
the *sources* CSV; the MatTen adapter is a different consumer that needs the
point-cloud emit (§4.1) instead of the mechanism source sets.

## 6. Model configuration (defaults, all justifiable as physics)

- **Output irreps:** fixed at `0e + 1e + 2e` (full shielding tensor,
  parity-correct). This is genuinely pinned by physics — a 3×3 tensor decomposes
  into exactly 0/1/2.
- **Hidden lmax = 2** is a **default / sensitivity-and-compute choice, not a
  physics fact.** The *output* is ℓ≤2, but the angular complexity of a many-atom
  protein environment is not guaranteed exhausted by hidden ℓ=2. lmax=2 is
  defensible and cheap; treat higher intermediate irreps as a sensitivity knob to
  check, not as physically irrelevant.
- **Cutoff radius** is both a physics statement and a *completeness /
  signal-visibility* parameter (§1b.3): `n_layers × rcut` sets how far a nucleus
  "sees," and σ is only a function of the input within that reach. **Set it from
  the LONGEST relevant range, not the shortest — this is a trap.** The DFT target
  has *no* spatial cutoff (whole-molecule QM + CPCM), and the extractor's own TOML
  already acknowledges long electrostatic reach: Coulomb/e-field **20 Å**, EEQ-CN
  **25 Å**, H-bond dipolar **50 Å** — far beyond ring-current 15 / McConnell 10 /
  rediscover charge 6 (the repo even warns "charge cutoff 6 truncates the
  long-range 1/r² field"). MatTen predicts the *total* σ, so rcut must be the max
  over mechanisms. **Default: `n_layers × rcut` ≈ domain extent (≥20–25 Å)** for
  the ~543-atom (~30–40 Å) domains — this covers the QM target's reach and the
  longest classical cutoff at once, and puts every atom in reach for the non-local
  CPCM term. Do NOT inherit a short calculator cutoff as rcut (would be incomplete
  vs both the QM target and the extractor's own electrostatic ranges). The
  extractor cutoffs thus give a transparent, derived *floor* for rcut, not a
  guess. The Stage-1 ridge baseline (§7.3) was built with these cutoffs, so
  MatTen-vs-ridge also tests "does seeing past the kernel cutoffs help" — name it,
  don't misread it as pure architecture.
- **Layers / body order:** modest depth for v1; each added layer is a
  longer-range claim, stated as such.

## 7. What the *predictor* honestly needs (beyond the two emits)

The C++ is small; the defensibility is mostly here.

### 7.1 DFT targets — locate and characterize (they exist: the mutant set)

**The 720 DFT targets exist — located authoritatively (Jessica + `MUTANT720_RECONCILE_2026-06-04.md`).**
Clean canonical source: **`/shared/2026Thesis/shielding-calcsets/data/orca-alphafold-and-mutants/`**
(exclude `disqualified/`) = **720/720 complete WT+ALA pairs**, each with `_nmr.out`
+ geometry for both halves. Mirrored DFT-only in
**`/shared/2026Thesis/nmr-shielding/calibration/`** (exclude `features/`, ~2.4 GiB,
in-repo). A sample `_nmr.out` carries 463 "Total shielding tensor" blocks →
**full per-atom tensors confirmed**. (codex reported "no corpus" because it
searched `nmr-720/calibration/` and the now-deleted `/shared/2026Thesis/consolidated`,
and globbed the wrong names — a broad `find -name '*_nmr.out'` finds ~14.5k of them.
The lesson: do the broad search, don't trust a narrow one.)

**Two structural facts that shape the training set and the split:**
- These are **AlphaFold structures**, keyed by UniProt IDs (e.g. `A0A062V9G2`),
  static single pose per half — *predicted* geometry, not experimental. A real
  data-quality caveat to disclose.
- The "ALA mutant" is **not a point mutation**: it knocks *all aromatic residues*
  (PHE/TYR/TRP/HIS) to ALA — often 6–8 residues per protein (built to isolate
  ring-current via the WT−ALA delta). So WT and ALA share a backbone scaffold but
  differ at many sidechains; the ALA halves are aromatic-stripped, artificial
  structures.

**Training-set decision:** MatTen trains on the **~720 WT AlphaFold poses** (real
proteins, full-tensor σ), leave-proteins-out. The **ALA halves are not run by
default** — they exist for the Stage-1 ridge (the WT−ALA delta that isolates
ring-current); they are not the predictor's data.

**ALA as an *optional* augmentation/probe (allowed — still pure geometry→σ, not a
v2).** The aromatic-knockout structures are free, valid `(geometry, DFT σ)`
supervision already on disk, so two uses are legitimate if we want them:
1. *Extra training data* — ~720 more structures; more is what a data-limited
   transferability pilot wants. Report whether the aromatic-depleted distribution
   shift helps or hurts held-out WT (empirical).
2. *Built-in physics probe* — a WT-trained model that still predicts the ALA
   shieldings (rings removed) is showing it learned the **geometric** ring-current
   dependence rather than memorizing "aromatic ⇒ shift." Evidence for the thesis,
   for free.
Either use keeps it pure geometry→σ (no kernel conditioning), so it does not
become the forbidden v2. Only care: an ALA shares its WT's backbone → **WT+ALA in
the same fold**. Fold by sequence/structure non-redundancy.

1P9J's 751 frames remain a within-protein conformational companion;
`fleet_calibration_dft/` holds per-frame DFT for ~10 further proteins (per-ns
sampling), another within-protein source.

**Remaining audit (characterize, not discover):** confirm absolute full-tensor σ
per WT pose across the 720 (effective ~676 after disulfide drops, CLAUDE.md),
count usable proteins, produce the coverage manifest (§7.7.1). The targets are
already spent and on disk — we characterize existing data, buy no new DFT.

### 7.2 Split — leakage-free, matched to the claim

A random split leaks (adjacent frames near-identical; same-protein frames
correlated) and yields a flattering, meaningless R².
- Single-protein → temporal/conformational holdout with a real gap.
- Multi-protein → **leave-proteins-out** (the Stage-1 "fair set" instinct).

### 7.3 Baselines — or the number means nothing

- Predict-the-per-element-mean (the floor).
- **Invariant baseline = the Stage-1 ridge-on-kernels.** MatTen-equivariant vs
  ridge is the clean "does carrying the tensor through beat scalarizing it"
  experiment — the silicate paper's 2× argument, on our own data. That
  comparison is the result; the raw MatTen number is not.

### 7.4 Loss and per-element scaling

The isotropic 0e is hundreds of ppm; T2 is smaller, T1 smaller still. Naive MSE
lets 0e swamp the anisotropy — and the anisotropy is the thesis. Use per-irrep
weighting (or standardization) so the model is forced to learn T2/T1, plus
per-element target normalization (H/C/N/O/**S** ranges and effective
dimensionalities differ — H≈20, C≈6, N≈3, O≈12; no simplification bias). Note S
is present in the DFT (3 atoms in the sample domain) and **H is the majority**
(277/543) — and H is exactly where the receptive-field and CPCM-non-locality
risks (§1b.3) bite hardest, being ring-current-dominated and long-range-sensitive.

### 7.5 Variance, not a point estimate

Small effective N → report across seeds and folds with intervals. Here R²/MAE
genuinely *is* the metric (the Stage-3 ethos flip), which is exactly why it
must be earned on a leak-free split with baselines and variance.

### 7.6 Framing held throughout

Predicts DFT shielding, validated against DFT. Not experimental shift. Stated
once, plainly, up front.

### 7.7 Required artifacts before the first number is defensible

From the review — concrete deliverables, not just principles:

1. **Target-coverage manifest** (machine-readable), produced *before* training:
   parseable frames, per-protein/per-element/per-atom counts, NaN/gap counts, and
   the Kabsch alignment stats. This is also how §7.1 gets pinned.
2. **A committed split artifact** beside the substrate. For 1P9J, temporal blocks
   with gaps between train/test; random per-row splits are invalid (adjacent
   frames are near-identical).
3. **An autocorrelation / nearby-frame baseline** in addition to the
   per-element-mean and ridge baselines — e.g. a block-aware or nearest-frame
   predictor — so the MatTen number is shown to beat "adjacent frames look alike,"
   not just exploit it.
4. **Metrics by irrep and element:** T0, T1 (vector + norm), T2 (component +
   norm), raw-3×3 Cartesian MAE/RMSE/R², and full-tensor reconstruction error
   after the inverse transform — reported per element (H/C/N/O), never pooled.
5. **Round-trip + equivariance tests:** raw-3×3 → `0e+1e+2e` → raw-3×3 identity;
   rotation equivariance; reflection parity (`v'=det(R)Rv`, §1); and row-grouped
   molecule reconstruction from the emitted substrate (§5 invariant).
6. **A compact cutoff/depth grid, not a large sweep.** DFT is the limiting
   resource; the first result should be robust, not exhaustively tuned. A spent
   Trp-cage run is an external sanity check, never hyperparameter data.

## 8. Non-goals

- No kernel-conditioned variant (v2). Geometry-only, full stop.
- No experimental-shift comparison dressed as validation.
- No changes to the extractor library or to `MopacResult`/analysis.
- No re-encoding of the "not nice" orientation features (sign-ambiguous ring
  normals / bond axes) — they are simply not inputs in v1, which dissolves the
  parity footgun rather than engineering around it.

## 9. Open questions

1. **Where does the mutant-set DFT live, and what is the coverage?** Locate the
   consolidated tree / mutant-pair extraction outputs; confirm absolute
   full-tensor σ per WT (and ALA) pose (`orca_total`, incl. T1); count usable
   proteins (~676 effective after the disulfide drops). This is
   characterization, not existence (§7.1).
2. The point-cloud + `target_T1` emit must land where the **720-WT mutant-set
   substrate** is built — i.e. the `nmr-720` worktree's rediscover path — since
   that is what MatTen trains on. Confirm that tree's equivariant sink matches
   the one reviewed here.
3. Cutoff radius: derive from ring-current range, or sensitivity-sweep?

## 10. Build order

1. **Locate + characterize the mutant-set DFT (§7.1) → coverage manifest
   (§7.7.1)** — read-only. Find the consolidated / mutant-pair tree, confirm
   absolute full-tensor σ per WT pose, count usable proteins. This sets the
   leave-proteins-out split (WT+ALA in the same fold; non-redundant folds); it
   does **not** gate whether the transferable claim exists — it does.
2. **Reinstate the reader-path element check** (§3) as a fail-loud guard — cheap,
   and required before any non-1P9J data is trusted.
3. Add the two emits (§4) — `positions` + `target_T1` — to
   `AllAtomEquivariantSink`, in whichever tree feeds training; unify the presence
   gate.
4. Write the MatTen adapter (§5): row→molecule regrouping (assert the invariant),
   explicit `0e+1e+2e` output irreps, non-periodic structure handling.
5. Commit the split artifact (§7.7.2); wire baselines (§7.3): ridge +
   per-element-mean + nearby-frame (§7.7.3) on that split.
6. Train with the §7.4 loss; report §7.5 variance and §7.7.4 per-irrep/element
   metrics.
7. Round-trip + rotation + reflection tests (§7.7.5 / §1).
