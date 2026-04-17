# Investigation Layers

Each layer applies one operation to the data and observes what
emerges. No layer depends on a model. Each produces figures and
findings independently. The GNN comes after, informed by what
the layers reveal.

---

## Layer 0: Census and join

**Question:** How many atoms can we match between extraction and RefDB?

**Method:**
- For each fleet protein, match extraction atom indices to RefDB
  shifts via (protein_id, residue_seq, atom_name)
- The chain_reconciliation table in protein_catalog handles
  PDB/BMRB numbering differences
- Count matched atoms per protein, per nucleus type (H, HA, HN,
  CA, CB, CO, N)
- Flag proteins with <50% match rate for investigation

**Output:** Match table. Per-protein coverage. Quality filter.

**Assumption tested:** The atom-matching join works at scale.

---

## Layer 1: Ensemble kernel statistics

**Question:** What do the kernel distributions look like across
5005 frames?

**Method:**
- For each protein, stream through frames (or use the 6
  initial-samples as pilot)
- Per atom: accumulate mean, variance, skewness of each kernel
  T0 (and optionally T2 magnitude) across frames
- The EnsembleConformation accumulator (spec/ENSEMBLE_MODEL.md)
  does exactly this

**Output:** Per-atom kernel moment tables. 685 proteins x ~1500
atoms x ~20 kernel T0 features x 3 moments = ~60M values, ~500 MB.

**What to look at:**
- Which kernels have the highest variance? (conformationally
  sensitive physics)
- Is variance correlated with ring distance? (expected: near-ring
  atoms have high BS variance)
- Do different elements show different variance patterns?
- Are there atoms with near-zero variance in ALL kernels?
  (structurally locked -- these are the easy cases)

**Assumption tested:** Kernel values actually vary meaningfully
across the ensemble. If everything is flat, ensemble analysis adds
nothing over single-frame.

---

## Layer 2: Residue-type centering of experimental shifts

**Question:** After removing the baseline, what's left?

**Method:**
- From RefDB, compute global mean shift per (residue_type,
  atom_name) across all 685 proteins
- Subtract. The residual delta_exp(i) = shift(i) - mean(residue,
  atom_name) is the geometry-dependent part
- Characterise the residual: what's its magnitude by nucleus?
  What's the per-protein variance?

**Output:** Centred shift residuals for every matched atom.

**Expected magnitudes (from literature):**
- HA residuals: ~0.5-1.5 ppm range (ring currents move HA by ~1 ppm)
- HN residuals: ~1-2 ppm range (H-bond + ring current)
- CA residuals: ~2-5 ppm range (backbone conformation dominates)
- N residuals: ~5-15 ppm range (multi-mechanism)

**Assumption tested:** The centering removes enough baseline that
geometric perturbations are above the noise. If the residuals are
pure noise, no kernel analysis will help.

---

## Layer 3: T0 kernel vs shift residual (the 4% finding)

**Question:** Do ensemble-averaged kernel T0 values explain any of
the centred shift residual?

**Method:**
- Per element, per nucleus type: correlate each kernel's mean T0
  against the centred shift residual
- This is a univariate correlation, not a fit
- Report per-kernel correlation coefficient and p-value

**The finding we're looking for:**
"Ensemble-averaged Biot-Savart T0 explains X% of residual HA
variance after residue-type centering, concentrated within 8A of
aromatic rings, following the 1/r^3 falloff that the calibration
independently verified."

If X = 4%, that's a three-sentence finding with a figure that has
never been shown at this scale. If X = 0%, that's also a finding
-- it means the ensemble average washes out the single-frame signal.

**Assumption tested:** T0 kernels computed from geometry carry
information about experimental shifts. This bridges the calibration
(engine works on mutants) to nature (engine sees real physics).

---

## Layer 4: Distance stratification

**Question:** Does the kernel-shift correlation follow the expected
distance dependence?

**Method:**
- Stratify atoms by distance to nearest aromatic ring (from
  ring_geometry data)
- Recompute Layer 3 correlation within each distance bin
- Bins: 0-4A, 4-6A, 6-8A, 8-12A, >12A

**Expected:** Ring current correlation strongest at <6A, following
r^-3 (BS/HM) or r^-5 (PQ). EFG correlation may not depend on ring
distance (it's about backbone electrostatics).

**The finding:** Distance-dependent kernel-shift correlation that
matches the theoretical falloff. Or doesn't -- and that's
interesting too.

**Assumption tested:** The r^-3 falloff validated on mutants
replicates on experimental ensemble data.

---

## Layer 5: Three-view convergence and kernel covariance

**Question:** Do BS, HM, and EFG_aro co-vary across conformations
(confirming geometric equivalence), and where do they diverge?

**Method:**
- Per protein, per atom: compute covariance matrix of kernel T0
  values across frames (or the 6 poses as pilot)
- Specifically: track BS-HM correlation, BS-EFG correlation,
  HM-EFG correlation per atom across frames
- Stratify by ring distance and ring type
- Average covariance matrices across atoms in the same stratum
- Look for block structure: ring current + EFG should form one
  correlated block (same geometry). McConnell, H-bond, backbone
  EFG should form separate blocks (different physics).

**The calibration predicts:**
BS, HM, and EFG_aro carry the same angular signal (cos ~0.93)
in single conformations. Across conformations they SHOULD co-vary
perfectly in the far-field. Near-field divergence is the signal:
- BS-HM divergence: circular-loop vs surface-integral model
- BS-EFG divergence: magnetic vs electric mechanism
- HM-EFG divergence: isolates the charge-vs-geometry question

**What this reveals:**
- Correlated block = same structural motion drives these kernels
  (expected for BS/HM/EFG_aro from the calibration finding)
- Divergence within the block = the three models see the geometry
  differently in that conformation (near-field, asymmetric rings)
- Independent blocks = different physics, different motions
- Off-diagonal surprise = unexpected coupling between mechanisms

**Assumption tested:** The three-view geometric equivalence holds
under conformational variation, and the divergence pattern is
physically interpretable (not noise).

---

## Layer 6: Per-ring tracking (conformational gating)

**Question:** Do individual ring-atom interactions gate on/off
across the ensemble?

**Method:**
- Use ring_contributions sparse array (59 columns per atom-ring pair)
- Track each ring's BS/HM T0 contribution to each atom across frames
- Classify rings as: always-on (contribution stable and significant),
  always-off (always negligible), gating (large variance, transitions
  between on and off)

**What this reveals:**
- Structurally locked ring currents vs dynamically modulated ones
- Which ring types gate? (HIE expected: small, mobile, variable
  protonation state. PHE expected: less mobile, more stable.)
- Gating frequency relative to NMR timescale

**Assumption tested:** Conformational gating of ring currents is
observable at the MD timescale. If all rings are either always-on or
always-off, gating is not a significant effect.

---

## Layer 7: Conditional structure (emergent)

**Question:** What patterns exist that we didn't predict?

**Method:**
- PCA on per-atom kernel moment profiles (from Layer 1)
- Cluster atoms by their kernel-distribution signatures
- For each cluster: what is it? Element? Residue type? Secondary
  structure? Ring distance? Or something else?
- Cross-reference clusters with shift residuals from Layer 2
- Look for clusters that have coherent shift residuals but don't
  correspond to obvious categories

**What this reveals:**
- If clusters are just element + residue_type: the kernels see
  what we already knew
- If clusters cut across residue types: the kernels see geometric
  structure that chemistry doesn't predict
- If clusters correlate with shift residuals: the geometric structure
  matters for NMR

**Assumption tested:** There exists interpretable structure in the
kernel-ensemble space beyond what element/residue categorisation
provides.

---

## Sequencing

Layers 0-2 are prerequisites (census, centering). Can be done on
the 6 initial-samples poses before committing to full trajectory
extraction.

Layer 3 is the key experiment. If T0-vs-residual shows signal, every
subsequent layer has something to decompose. If it shows nothing,
we need to understand why before proceeding.

Layers 4-5 deepen the physics. Layers 6-7 are discovery.

The GNN work (nmr-training/) runs in parallel. It doesn't wait for
the analytical layers. But the analytical findings INFORM the GNN
architecture: which kernels to weight, which couplings to allow,
what the learned residual should look like.

---

## What connects to the GNN

The analytical findings become the GNN's design spec:

- Layer 3 tells the GNN which kernel features carry signal per
  nucleus type (feature importance prior)
- Layer 4 tells the GNN the effective range of each kernel
  (radius graph cutoff, distance weighting)
- Layer 5 tells the GNN which tensor product channels to support
  (kernel couplings that co-vary)
- Layer 6 tells the GNN about temporal averaging (whether to
  predict per-frame or per-ensemble)
- Layer 7 tells the GNN what structure exists beyond the classical
  kernels (the learned residual should capture this)

The thesis statement becomes: "The classical kernels explain X% of
the geometry-dependent shift variation, decomposed into [ring
current / EFG / bond anisotropy] contributions that follow the
theoretically predicted [distance / angular / element] dependence.
The equivariant GNN recovers an additional Y% by learning [specific
residual patterns], corroborating the analytical decomposition."
