# Analysis Trajectory Mode — Tentative Design

**Date:** 2026-04-14
**Status:** TENTATIVE. Working design, will evolve with the data.

---

## What this is

A 10-protein workspace that carries us from exploratory analysis
through to runnable models. The directory at
`/shared/2026Thesis/fleet_calibration/` holds 10 proteins, each
with a 25 ns plain MD trajectory (1,250 frames, ~335-5000 protein
atoms, full-system XTC with explicit water and ions).

This is an integration step. The analysis H5 file is an event
collector whose schema will evolve as we learn what signals the
trajectory contains. The 10 proteins are the proving ground for
every downstream consumer: R analysis, DFT point selection, GNN
message design, time series model training, and UI visualization.

**The worst mistake we can make is conflating all downstream
consumers into one design.** Each phase has its own questions and
its own outputs. We advance one phase at a time, using the same
10-protein workspace, refining the schema as we go.

---

## Phases (sequential, not parallel)

### Phase 1: Analysis for DFT selection (NOW)

**Question:** Where in the trajectory should we run DFTs?

**Method:** Write exhaustive per-frame data to analysis H5.
Explore in R. Find frames where the physics dimensions diverge
most from their means, where structural events cluster, where
new features (water, charges) show unexpected correlations.
Pick ~10 DFT frames per protein.

**Output:** `{protein_id}_analysis.h5` with full per-atom and
per-residue time series. This is the only H5 we design now.

### Phase 2: DFT validation and signal refinement

**Question:** Do the trajectory kernel variations predict the
DFT shielding variations?

**Method:** Run ORCA on the selected frames. Compare DFT T2
against kernel T2 per frame. Identify which dimensions carry
conformational signal and which are noise.

**Output:** Refined understanding of which channels matter.
May trigger re-extraction with modified parameters.

### Phase 3: GNN message design

**Question:** What discrete events should the GNN receive as
messages?

**Method:** Using the analysis H5 + DFT results, define event
types (rotamer flips, H-bond changes, water shell exchange,
charge jumps). Write event detectors in Python. Populate the
`events/` section of the H5.

**Output:** Event schema, Python detectors, populated events.

### Phase 4: Production time series

**Question:** What per-frame data does the time series model
need, and what does the UI need?

**Method:** Using everything learned in phases 1-3, define
the production `--trajectory --timeseries` output. This is
leaner than the analysis H5 — only channels that proved useful.

**Output:** Production H5 schema, C++ writer, SDK reader.

---

## Phase 1 design: Analysis H5

### Run parameters

- **Frame selection:** uniform stride to ~600 frames from 1,250
  (stride 2, skip odd frames via `Skip()`)
- **Calculators:** ALL (including APBS, AIMNet2). No cheap subset.
  Estimated ~100 min per protein on 5090.
- **MOPAC:** skipped (10 min/frame, not affordable per-frame)
- **Coulomb:** skipped (APBS replaces it)
- **Buffering:** entire protein buffered in memory, H5 written
  once at end. ~560 MB per protein at float64. Fine for 64 GB.
- **Overwrite:** re-extract from scratch when anything changes.
  No incremental writes. Forward always forward.

### CLI

    nmr_extract --trajectory --analysis \
                --tpr topology.tpr --xtc trajectory.xtc \
                --output /path/to/output

Single pass over XTC with stride. All calculators run on each
sampled frame. Analysis H5 written at end.

### H5 schema

The file is self-contained. R reads it with `rhdf5` without
needing the Python SDK or understanding the C++ object model.
Groups are for R memory management — load one group at a time.
Names match the physics, not the tool.

```
{protein_id}_analysis.h5

meta/
  frame_times          (T,)     float64   # ps timestamps for sampled frames
  frame_indices        (T,)     int32     # original XTC frame indices
  attrs:
    protein_id         string             # e.g. "1B1V_4292" (from directory name)
    n_atoms            int                # protein atom count
    n_frames           int                # harvested frame count
    n_residues         int                # residue count
    stride             int                # frame stride used (e.g. 2)

atoms/
  element              (N,)     int32     # atomic number: 1=H, 6=C, 7=N, 8=O, 16=S
  residue_index        (N,)     int32     # 0-based, which residue
  atom_name            (N,)     string    # "CA", "CB", "N", "HB2" (display only)
  atom_role            (N,)     int32     # AtomRole enum
  hybridisation        (N,)     int32     # sp2/sp3/...
  n_bonded             (N,)     int32     # covalent neighbors
  graph_dist_ring      (N,)     int32     # bond hops to nearest ring
  is_backbone          (N,)     int8      # boolean
  is_conjugated        (N,)     int8

residues/
  residue_name         (R,)     string    # "ALA", "GLY", "PRO"
  residue_number       (R,)     int32     # PDB numbering
  chain_id             (R,)     string    # chain identifier

topology/
  # Bond graph, ring membership, hydrogen→heavy atom mapping.
  # This is what forces the Python consumer to build a protein model.
  bond_atoms           (B,2)    int32     # atom index pairs
  bond_category        (B,)     int32     # BondCategory enum
  bond_order           (B,)     int32     # BondOrder enum
  parent_atom_index    (N,)     int32     # heavy atom parent (-1 = not hydrogen)
  # Rings: CSR ragged array (ring_offsets[i]:ring_offsets[i+1] into ring_atom_indices)
  ring_type            (n_rings,)    int32     # RingTypeIndex enum
  ring_residue         (n_rings,)    int32     # parent residue index
  ring_fused_partner   (n_rings,)    int32     # fused ring index (-1 = none)
  ring_offsets         (n_rings+1,)  int32     # CSR offsets into ring_atom_indices
  ring_atom_indices    (K,)          int32     # flat atom indices for all rings
  attrs:
    n_bonds            int
    n_rings            int

positions/
  xyz                  (T,N,3)  float64   # Angstroms, PBC-fixed

ring_current/
  # Biot-Savart per ring type (8 types)
  bs_T0_per_type       (T,N,8)     float64   # isotropic per ring type
  bs_T2_per_type       (T,N,8,5)   float64   # T2 components per ring type
  bs_shielding         (T,N,9)     float64   # total SphericalTensor (T0,T1[3],T2[5])
  # Haigh-Mallion per ring type
  hm_T0_per_type       (T,N,8)     float64
  hm_T2_per_type       (T,N,8,5)   float64
  hm_shielding         (T,N,9)     float64
  # Ring susceptibility
  rs_shielding         (T,N,9)     float64   # ringchi_shielding_contribution
  # Ring proximity
  n_rings_3A           (T,N)       int16
  n_rings_5A           (T,N)       int16
  n_rings_8A           (T,N)       int16
  n_rings_12A          (T,N)       int16
  mean_ring_dist       (T,N)       float64
  nearest_ring_atom    (T,N)       float64
  # Exponential-weighted sums
  G_iso_exp_sum        (T,N)       float64
  G_T2_exp_sum         (T,N,5)     float64
  G_iso_var_8A         (T,N)       float64
  # B-field
  total_B_field        (T,N,3)     float64

efg/
  # Coulomb EFG (ff14SB charges) — decomposed by source
  coulomb_total        (T,N,9)     float64   # SphericalTensor
  coulomb_backbone     (T,N,9)     float64
  coulomb_aromatic     (T,N,9)     float64
  # Coulomb E-field vectors
  E_total              (T,N,3)     float64
  E_backbone           (T,N,3)     float64
  E_sidechain          (T,N,3)     float64
  E_aromatic           (T,N,3)     float64
  E_solvent            (T,N,3)     float64
  # Coulomb scalars
  E_magnitude          (T,N)       float64
  E_bond_proj          (T,N)       float64
  E_backbone_frac      (T,N)       float64
  # APBS solvated
  apbs_efg             (T,N,9)     float64   # SphericalTensor
  apbs_efield          (T,N,3)     float64
  # AIMNet2 charge-derived EFG — decomposed
  aimnet2_total        (T,N,9)     float64
  aimnet2_backbone     (T,N,9)     float64
  aimnet2_aromatic     (T,N,9)     float64
  aimnet2_shielding    (T,N,9)     float64

bond_aniso/
  # McConnell per bond category
  mc_shielding         (T,N,9)     float64   # total SphericalTensor
  T2_backbone          (T,N,9)     float64   # backbone bonds
  T2_sidechain         (T,N,9)     float64
  T2_aromatic          (T,N,9)     float64
  T2_CO_nearest        (T,N,9)     float64
  T2_CN_nearest        (T,N,9)     float64
  # McConnell scalars
  co_sum               (T,N)       float64
  cn_sum               (T,N)       float64
  sidechain_sum        (T,N)       float64
  aromatic_sum         (T,N)       float64
  co_nearest           (T,N)       float64
  nearest_CO_dist      (T,N)       float64
  nearest_CN_dist      (T,N)       float64
  nearest_CO_midpoint  (T,N,3)     float64
  dir_nearest_CO       (T,N,3)     float64

quadrupole/
  pq_shielding         (T,N,9)     float64   # total SphericalTensor
  pq_T0_per_type       (T,N,8)     float64
  pq_T2_per_type       (T,N,8,5)   float64

dispersion/
  disp_shielding       (T,N,9)     float64   # total SphericalTensor
  disp_T0_per_type     (T,N,8)     float64
  disp_T2_per_type     (T,N,8,5)   float64

hbond/
  hbond_shielding      (T,N,9)     float64   # SphericalTensor
  nearest_spherical    (T,N,9)     float64
  nearest_dist         (T,N)       float64
  nearest_dir          (T,N,3)     float64
  inv_d3               (T,N)       float64
  count_3_5A           (T,N)       int16
  is_donor             (T,N)       int8
  is_acceptor          (T,N)       int8
  is_backbone          (T,N)       int8

sasa/
  sasa                 (T,N)       float64   # Shrake-Rupley (A^2)
  normal               (T,N,3)     float64   # outward surface normal

water/
  # E-field from all waters
  efield               (T,N,3)     float64
  efg                  (T,N,9)     float64   # SphericalTensor
  # First shell only
  efield_first         (T,N,3)     float64
  efg_first            (T,N,9)     float64
  # Shell counts
  n_first              (T,N)       int16
  n_second             (T,N)       int16
  # Hydration shell geometry (COM-based, HydrationShellResult)
  half_shell_asymmetry (T,N)       float64
  dipole_cos           (T,N)       float64
  nearest_ion_dist     (T,N)       float64
  nearest_ion_charge   (T,N)       float64
  # Hydration geometry (SASA-normal, HydrationGeometryResult)
  dipole_vector        (T,N,3)     float64
  surface_normal       (T,N,3)     float64
  sasa_asymmetry       (T,N)       float64
  sasa_dipole_align    (T,N)       float64
  sasa_dipole_cohere   (T,N)       float64
  sasa_first_shell_n   (T,N)       int16

charges/
  aimnet2_charge       (T,N)       float64   # Hirshfeld (e)
  eeq_charge           (T,N)       float64   # Caldeweyher (e)
  eeq_cn               (T,N)       float64   # coordination number

aimnet2_embedding/
  aim                  (T,N,256)   float64   # learned electronic embedding

dihedrals/
  # Per-residue, in radians. NaN for undefined (terminal, GLY/ALA chi).
  phi                  (T,R)       float64
  psi                  (T,R)       float64
  omega                (T,R)       float64   # peptide planarity — NEW
  chi1                 (T,R)       float64
  chi2                 (T,R)       float64
  chi3                 (T,R)       float64
  chi4                 (T,R)       float64
  # Also store cos/sin pairs for circular analysis
  chi1_cos             (T,R)       float64
  chi1_sin             (T,R)       float64
  chi2_cos             (T,R)       float64
  chi2_sin             (T,R)       float64
  chi3_cos             (T,R)       float64
  chi3_sin             (T,R)       float64
  chi4_cos             (T,R)       float64
  chi4_sin             (T,R)       float64

dssp/
  ss8                  (T,R)       int8      # 8-class secondary structure
  hbond_energy         (T,R)       float64   # strongest acceptor energy

predictions/
  # Ridge applied to this frame's kernels. Weights from Stage 1
  # calibration (per-element, 55 kernels). Additive: total = sum
  # of group contributions. We are SETI — collect the signal,
  # do not editorialize. Python checks our work from the raw.
  #
  # Raw: ridge applied to unnormalized kernel features.
  # Norm: ridge applied after per-frame per-protein normalization.
  #
  # GROUNDING (Stage 1, 720 proteins, 446K atoms):
  # The groups below are the physics groups from the dimension
  # tree (stage1-mutations/notes/dimension_tree.md). Each has
  # measured element-specific R² and measured cross-group
  # independence (cosines). These are not ad hoc buckets —
  # they are the calibrated decomposition of the shielding
  # perturbation into independent physical dimensions.
  #
  #   Group         | H R²_n | C R²_n | N R²_n | O R²_n | Independence
  #   ring_current  | .782   | .121   | .124   | .195   | cos(rc,efg)=0.44
  #   efg           | .701   | .380   | .064   | .150   | cos(efg,ba)=0.39
  #   bond_aniso    | .107   | .025   | .029   | .041   | cos(ba,rc)=0.40
  #   quadrupole    | .030   | .013   | .063   | .038   | cos(pq,rc)=0.45
  #   dispersion    | .387   | .115   | .068   | .234   | cos(d,efg)=0.53
  #   hbond         | .021   | .017   | .006   | .018   | zero for mutants
  #
  # The additive sum of group contributions IS the prediction.
  # "These two together means more" is grounded in the measured
  # independence — ring_current + efg add because cos=0.44
  # (nearly random in 5D), not because we decided they should.
  # New feature groups (water, charges, SASA) were zero for
  # mutants. They ride along with zero weights. If they light
  # up against DFT in Phase 2, they earn their weights.
  #
  # NONLINEAR SIGNAL (Stage 1 RF delta):
  #   H: +0.002 — linear. Ridge is complete.
  #   C: +0.128 — nonlinear. Ridge misses 26% relative.
  #   N: +0.169 — nonlinear. Ridge misses 63% relative.
  #   O: +0.013 — linear. Ridge is complete.
  #
  # For N and C, the gated MLP (KernelMixingHead: MLP → per-kernel
  # weights, gated by kernel self-reported magnitude) captures
  # structure the ridge cannot. This is per-element gating for
  # Stage 2. Both predictions are written — ridge (linear) and
  # gated MLP (nonlinear, N and C). The difference between them
  # is the nonlinear signal.
  #
  # DO NOT add editorial metrics (confidence, agreement,
  # stability). The contributions are the signal. The stats
  # say what "more" means. Ad hoc reasoning does not.
  raw_T0               (T,N)       float64   # sum of weighted kernel T0
  raw_T2               (T,N,5)     float64   # sum of weighted kernel T2
  norm_T0              (T,N)       float64
  norm_T2              (T,N,5)     float64
  # Per-group additive contributions (what adds up to the total).
  # Each is the weighted sum of kernels in that group only.
  # New features (water, charges) get their own groups with
  # initially zero weights — they ride along for free.
  raw_ring_current     (T,N,5)     float64   # BS+HM+RS contribution
  raw_efg              (T,N,5)     float64   # Coulomb+APBS+AIMNet2
  raw_bond_aniso       (T,N,5)     float64   # McConnell
  raw_quadrupole       (T,N,5)     float64
  raw_dispersion       (T,N,5)     float64
  raw_hbond            (T,N,5)     float64
  norm_ring_current    (T,N,5)     float64
  norm_efg             (T,N,5)     float64
  norm_bond_aniso      (T,N,5)     float64
  norm_quadrupole      (T,N,5)     float64
  norm_dispersion      (T,N,5)     float64
  norm_hbond           (T,N,5)     float64
  # New feature groups — zero weights initially, signal rides along
  raw_water            (T,N,5)     float64   # water EFG contribution
  raw_charges          (T,N,5)     float64   # charge-derived contribution
  raw_sasa             (T,N,5)     float64   # surface contribution
  norm_water           (T,N,5)     float64
  norm_charges         (T,N,5)     float64
  norm_sasa            (T,N,5)     float64
  # Gated MLP predictions (nonlinear, for N and C)
  # KernelMixingHead: scalar context → per-kernel weights × T2.
  # Written for all elements but the N/C signal is where it matters.
  mlp_T0               (T,N)       float64
  mlp_T2               (T,N,5)     float64
  # Per-kernel gating weights from the MLP (what the model chose)
  mlp_kernel_weights   (T,N,55)    float64   # 55 kernel weights per atom

projections/
  # Mathematical projections from the calibrated model. These are
  # dot products and distances — no editorial judgment.
  #
  # The model carries the Stage 1 eigenvectors, residual directions,
  # and training distribution statistics. It applies them per atom
  # per frame via the exported TorchScript model.
  #
  # 1. Element-specific subspace coordinates.
  # PCA on normalized kernels identified per-element predictive
  # subspaces (H=20, C=6, N=3, O=12 dims). Project each atom's
  # kernel vector into its element's subspace. R sees the protein
  # moving through shielding-relevant space over time.
  subspace_coords      (T,N,20)    float64   # padded to max_dims=20, zeros beyond element's count
  subspace_n_dims      (N,)        int8      # actual dims per atom (from element)
  #
  # 2. Residual projection of new features.
  # The ridge residual (what 55 kernels can't explain) has a
  # direction in kernel space. Project new trajectory features
  # onto it: "how much does this new feature align with the
  # unexplained variance?"
  residual_water_efg   (T,N)       float64   # water EFG projected onto residual
  residual_water_first (T,N)       float64   # first-shell water EFG onto residual
  residual_aimnet2_efg (T,N)       float64   # AIMNet2 EFG onto residual
  residual_hydration   (T,N)       float64   # hydration geometry onto residual
  residual_charges     (T,N)       float64   # EEQ/AIMNet2 charge signal onto residual
  #
  # 3. Training distribution distance.
  # Mahalanobis distance of each atom's kernel vector from the
  # 720-protein calibration distribution. Not confidence — distance.
  # High = trajectory exploring uncalibrated geometry.
  # Those frames are where DFTs add the most new information.
  mahalanobis          (T,N)       float64
  #
  # 4. MLP gating decomposition (N and C nonlinear signal).
  # Per-kernel weights chosen by the gated MLP. Shows WHERE
  # the nonlinear signal lives at each atom in each frame.
  mlp_gate_weights     (T,N,55)    float64   # core kernel gating weights

events/
  # PLACEHOLDER. Populated by Python event detectors in Phase 3.
  # Schema TBD after exploratory analysis reveals what events matter.
```

### Dimensions

| Symbol | Meaning | Example (1B1V) |
|--------|---------|-----------------|
| T | sampled frames | 600 |
| N | protein atoms | 335 |
| R | residues | ~22 |

### Precision: float64 throughout

All datasets are float64 (double). The C++ calculators compute in
double. The signal is faint — kernel T2 values live around 0.01-0.04
ppm, and we need to see whether group contributions add or cancel
at the fourth decimal place. float64's 7 significant digits would
wash out exactly the subtlety we need. float64's 15 digits preserve
the full computation. The size cost (2x vs float64) is irrelevant
for 10 proteins.

### Estimated sizes (1B1V, 335 atoms, 600 frames, float64)

| Group | Channels per atom | Size |
|-------|------------------|------|
| positions | 3 | 4.8 MB |
| ring_current | ~130 | 210 MB |
| efg | ~70 | 112 MB |
| bond_aniso | ~76 | 122 MB |
| quadrupole | ~57 | 92 MB |
| dispersion | ~57 | 92 MB |
| hbond | ~30 | 48 MB |
| sasa | 4 | 6.4 MB |
| water | ~30 | 48 MB |
| charges | 3 | 4.8 MB |
| aimnet2_embedding | 256 | 412 MB |
| predictions | ~40 | 64 MB |
| **Total per-atom** | **~756** | **~1.2 GB** |
| dihedrals (per-residue) | 15 × R | ~1.6 MB |
| dssp (per-residue) | 2 × R | ~0.2 MB |

~1.2 GB per protein. Buffers fine on 64 GB (one protein at a time).
10 proteins = ~12 GB on disk. R loads one group at a time. Largest
single load: aimnet2_embedding at 412 MB. Fine on 64 GB.

---

## What is NOT in this spec

### Production time series (Phase 4)

The `--trajectory --timeseries` flag on nmr_extract exists in the
architecture but its schema is undefined. It will be a strict
subset of the analysis H5, determined by what phases 1-3 prove
useful. Do not design it now.

### GNN message format (Phase 3)

The `events/` group is a placeholder. Event types, message schema,
and detection logic depend on what we find in the analysis. The
Python session that fills this section will be informed by the R
analysis from Phase 1 and the DFT results from Phase 2.

### UI trajectory visualization (Phase 4+)

The viewer needs a digestible time series for interactive display.
That is a further reduction of the production time series, driven
by what the UI can render meaningfully. Not designed here.

---

## Implementation scope (Phase 1 only)

### C++ changes

1. **AnalysisWriter** class (new). Buffers per-frame data from
   ConformationAtom. Writes the analysis H5 via HighFive at end.
   One method per H5 group: `HarvestRingCurrent(conf)`,
   `HarvestEFG(conf)`, etc. Each reads fixed fields from
   ConformationAtom / FindResult<T>. No runtime configuration.

2. **RidgePredictor** class (new). Loads calibrated ridge weights
   from a file (exported from Python calibration — per-element
   weight vectors, normalization statistics). Per frame: assembles
   the 55 kernel T2 vectors from ConformationAtom, applies
   normalization, computes weighted sums per group and total.
   Pure linear algebra — dot products on the data we already have.
   Called by AnalysisWriter after harvest. New feature groups
   (water, charges, SASA) have zero weights initially and ride
   along — their raw kernel values are projected through the same
   machinery so they appear in the decomposition.

3. **DsspResult**: add omega dihedral computation (CA-C-N-CA).
   Store phi, psi, omega, chi1-4 in radians on the per-residue
   struct. Currently only phi/psi stored; chi computed at write
   time. Move chi computation to Compute() so it is available
   for the analysis writer.

4. **GromacsFrameHandler**: new `--analysis` code path. Single
   pass with stride. Calls AnalysisWriter::Harvest() after each
   frame's calculators. Calls AnalysisWriter::WriteH5() at end.
   Does NOT run the two-pass scan/extract pattern.

5. **JobSpec**: add `--analysis` flag (mutually exclusive with
   `--timeseries` when that eventually has content).

6. **Model export** (Python session builds this):

   a. **Ridge weights** — per-element weight vectors for the
      trajectory kernel layout (~50 kernels). TOML or JSON.
      New feature slots have zero weights. ~4 KB.

   b. **Normalization stats** — std floor, per-kernel reference
      statistics from the 720-protein calibration set.

   c. **Gated MLP** — TorchScript `.pt` for KernelMixingHead.
      Input: kernel magnitudes + kernel scales + new scalar
      features (SASA, water counts, charges, hydration).
      Output: per-kernel gating weights.
      Includes the N/C nonlinear signal. Small model (~50 KB).

   d. **PCA eigenvectors** — per-element predictive subspace
      from Stage 1 (H=20, C=6, N=3, O=12 directions in
      normalized kernel space). Serialized as arrays.

   e. **Residual directions** — per-element residual vector
      (what the ridge can't explain), for projecting new
      features onto. From Stage 1 analysis.

   f. **Training distribution** — mean + inverse covariance
      of kernel vectors across 720 proteins, per element.
      For Mahalanobis distance computation.

   g. **Kernel manifest** — ordered list of kernel names mapping
      positions in the weight vector to ConformationAtom fields.
      This is the contract between Python and C++. Both sides
      assemble/interpret kernels in this order.

   The Python calibration pipeline is the source of truth. C++
   loads these artifacts and applies them. Python checks C++'s
   work from the raw kernels in the same H5.

### Python changes

5. **SDK**: `load_analysis(path)` reader for the analysis H5.
   Returns typed groups with named accessors. Not urgent for
   Phase 1 — R reads the H5 directly.

### Not changed

- GromacsProtein, GromacsProteinAtom, Welford accumulators
- OperationRunner, calculator dispatch
- Existing two-pass trajectory mode
- Production H5 master file
- Any existing use case (A-E)

---

## Design constraints

1. **The 10-protein workspace is the integration vehicle.**
   Every downstream consumer — R analysis, DFT selection, GNN,
   time series model, UI — will be developed against these 10
   proteins. The analysis H5 is the common data layer.

2. **One phase at a time.** Do not design GNN messages while
   building the analysis writer. Do not design production time
   series while exploring in R. Each phase informs the next.

3. **Forward always forward.** Re-extract from scratch when
   anything changes. No incremental updates, no layered outputs.

4. **The H5 is the interface.** It must be self-contained and
   self-evident. A cold Python or R session with no context
   should be able to open it and understand what each dataset is.

5. **Write everything the conformation knows.** The analysis
   is not reductionist. We are looking for signals we don't
   expect. Editorial decisions about what matters come after
   we have looked, not before.

6. **Buffer, don't stream.** ~575 MB per protein fits in memory
   on all 4 machines. Write the H5 once at end. Either you get
   a complete file or nothing.

7. **Per-atom and per-residue are separate axes.** Do not
   broadcast per-residue quantities to atoms. Dihedrals and
   SS are (T, R). Everything else is (T, N).

8. **Dihedrals in radians, not cosines.** The full angle is
   needed for event detection. Store cos/sin pairs alongside
   for circular statistics.

---

## Build and test plan

### Build sequence

Three workstreams, partially parallel:

**Workstream A: Python model session (learn/ side)**

Builds and exports the model artifacts. Works entirely against the
existing 720-protein calibration data. No C++ changes needed.

1. **Define trajectory kernel layout.** Modify KernelLayout to
   produce a trajectory variant:
   - Drop 10 MOPAC/DeltaAPBS slots (always zero in trajectory)
   - Add 5 new T2 slots: AIMNet2 EFG (total, backbone, aromatic),
     water EFG (total, first shell)
   - Keep per-ring kernels at K=6 (11 calcs × 6 nearest = 66)
   - Total: ~111 kernels (45 surviving core + 5 new + 66 per-ring)
     minus whatever MOPAC per-ring slots drop
   - Write the kernel manifest (ordered name list → JSON)

2. **Fit ridge per element** on the trajectory layout. The
   calibration data fills core + per-ring slots; new T2 slots are
   zero (calibration proteins have no trajectory water/AIMNet2 EFG
   data). Ridge weights for new slots will be zero — they ride along.

3. **Train gated MLP** (KernelMixingHead) on the trajectory layout.
   Scalar inputs: kernel magnitudes, kernel scales, plus new scalar
   features (SASA, water counts, charges, hydration) — all zero in
   calibration data, but the architecture accommodates them so they
   activate when trajectory data arrives.

4. **Extract Stage 1 artifacts:**
   - PCA eigenvectors per element (from full_space_analysis.py
     output, recomputed on the trajectory kernel layout)
   - Residual directions per element (ridge residual vector in
     normalized kernel space)
   - Training distribution per element (mean + covariance of
     kernel vectors across 720 proteins, for Mahalanobis)

5. **Export everything:**
   - `trajectory_ridge_weights.json` — per-element weight vectors
   - `trajectory_norm_stats.json` — normalization reference stats
   - `trajectory_model.pt` — TorchScript gated MLP
   - `trajectory_pca.npz` — eigenvectors per element
   - `trajectory_residual.npz` — residual directions per element
   - `trajectory_distribution.npz` — mean + inv_cov per element
   - `trajectory_kernel_manifest.json` — ordered kernel names

6. **Self-test:** apply all exports to the 720 calibration proteins
   in Python. Verify ridge predictions match the Stage 1 results
   (core kernels identical, new kernel contribution = zero). Verify
   PCA projections, Mahalanobis distances are sane. Save reference
   predictions for 10 test proteins.

**Workstream B: C++ implementation (src/ side)**

Builds the AnalysisWriter and model loading. Can start on the
harvester before the model exports are ready.

1. **AnalysisWriter** — buffer allocation, per-group harvest
   methods, H5 writing via HighFive. Test with just raw data
   groups (no predictions) initially.

2. **DsspResult** — add omega dihedral, move chi to Compute(),
   store radians + cos/sin.

3. **Kernel assembler** — reads ConformationAtom fields into a
   flat kernel T2 vector in manifest order. This is the C++ side
   of the kernel manifest contract. Must exactly match the Python
   assembly order.

4. **RidgePredictor** — loads ridge weights + normalization stats.
   Applies per-frame: normalize, weight, decompose by group.

5. **TorchScript model loader** — loads `trajectory_model.pt`,
   runs per-frame inference for MLP predictions and gating weights.
   Same pattern as AIMNet2 model loading.

6. **Projections** — loads PCA eigenvectors, residual directions,
   distribution stats. Computes subspace coords, residual
   projections, Mahalanobis per atom per frame.

7. **GromacsFrameHandler** — `--analysis` code path. Single pass,
   stride, calls AnalysisWriter per frame, writes H5 at end.

8. **JobSpec** — `--analysis` flag.

**Workstream C: Verification (Python, after A+B)**

Three independent checks that together prove correctness.

1. **Python self-consistency (720 proteins).**
   Already done in Workstream A step 6. The model exports applied
   in Python reproduce the Stage 1 calibration results. This tests
   model correctness independent of C++.

2. **C++ kernel assembly vs existing NPY output.**
   Run `nmr_extract --trajectory` on 1B1V (one frame) using the
   existing two-pass mode, which writes per-calculator NPY files.
   Run `nmr_extract --trajectory --analysis` on the same frame.
   Load both in Python. Verify that the raw kernel values in the
   analysis H5 match what the SDK loads from the NPY files.
   This tests that the C++ kernel assembler agrees with the Python
   SDK's kernel assembly.

   ```
   # Pseudocode
   from nmr_extract import load, load_analysis
   p = load("1B1V_output/frame_0000/")
   a = load_analysis("1B1V_analysis.h5")

   # Assemble kernels from NPY via SDK (Python ground truth)
   K_python = assemble_kernels(p, all_atoms, trajectory_layout)

   # Read kernels from analysis H5 (C++ assembled)
   K_cpp = reconstruct_kernels_from_h5_groups(a, frame=0)

   assert np.allclose(K_python, K_cpp, atol=1e-12)
   ```

3. **End-to-end: Python recomputes from H5 raw data.**
   Load the analysis H5. From the raw physics groups (ring_current/,
   efg/, bond_aniso/, etc.), reassemble the kernel T2 matrix.
   Apply the Python ridge weights, normalization, MLP, projections.
   Compare against the predictions/ and projections/ datasets in
   the same H5. This is the ongoing verification — runs on every
   extraction.

   ```
   verify_analysis.py {protein}_analysis.h5

   For each frame t:
     K = assemble_kernels_from_raw_groups(h5, t)  # (N, n_kernels, 5)
     # Ridge
     raw_pred = apply_ridge(K, weights)
     assert close(raw_pred, h5["predictions/raw_T2"][t])
     # Normalized ridge
     K_norm, scales = normalize(K)
     norm_pred = apply_ridge(K_norm, weights)
     assert close(norm_pred, h5["predictions/norm_T2"][t])
     # MLP
     mlp_pred = apply_mlp(K, scales, scalars)
     assert close(mlp_pred, h5["predictions/mlp_T2"][t])
     # Projections
     coords = project_subspace(K_norm, eigenvectors, elements)
     assert close(coords, h5["projections/subspace_coords"][t])
     # Mahalanobis
     maha = mahalanobis(K_norm, mean, inv_cov, elements)
     assert close(maha, h5["projections/mahalanobis"][t])

   PASS: all frames, all datasets, max |diff| < 1e-10
   ```

### Bootstrap sequence

The workstreams have dependencies:

```
A1 (kernel manifest) ──→ B3 (C++ kernel assembler)
A2-A5 (model export) ──→ B4-B6 (load + apply)
B1 (AnalysisWriter)  ──→ B7 (frame handler integration)
B1-B7 complete       ──→ C2 (kernel assembly check)
A6 + B7 complete     ──→ C3 (end-to-end verification)
```

**What can run in parallel:**
- A1-A6 (Python model session) and B1-B2 (AnalysisWriter + DSSP)
  are independent. Start both immediately.
- B3 needs A1 (the manifest). B4-B6 need A2-A5 (the exports).

**First milestone:** AnalysisWriter writes raw data groups only
(no predictions, no projections). Verify the H5 is loadable in R,
groups are correctly shaped, values are non-zero. This tests the
harvest + buffer + write pipeline without the model.

**Second milestone:** Model exports loaded, predictions and
projections populated. Run verify_analysis.py. This tests the
full pipeline.

**Third milestone:** Run on all 10 proteins. Load in R. Look at
the signal. Phase 1 is complete.
