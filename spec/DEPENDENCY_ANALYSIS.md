# Dependency Analysis — Complete Pipeline Map

## Pipeline Stages

```
Stage 1: ORCA DFT (external, pre-existing)
  → data/orca_runs/{protein}/

Stage 2: Delta Computation (batch_compute_deltas, v2 library)
  → data/deltas/{protein}_delta.json  (554 proteins)

Stage 3: Protonation + APBS (scripts/prep_protein.py + tleap + APBS)
  → data/prepared/{protein}/protonated.pdb, efield_solvated.bin

Stage 4: Feature Extraction (bs-equivariant, v1 binary)
  Input: data/deltas/ + data/ff14sb_params.dat + data/orca_runs/ + data/prepared/
  → data/equivariant/  (7 files: L0_features.bin, L1_features.bin, L2_features.bin, targets.bin, manifest.json, protein_ids.txt, elements.txt)

Stage 5: Training (python/train.py)
  Input: data/equivariant/
  → models/best_bilinear.pt, models/best_e3nn.pt

Stage 6: System Evaluation (python/evaluate_system.py)
  Input: data/equivariant/ + models/*.pt
  → results/system_eval_*.json

Stage 7: Inference (C++ eval module + trained models — TO BUILD)
  Input: PDB file + trained models
  → per-atom {tensor, confidence, tier}
```

## Two Extractors (MUST RECONCILE)

### Working extractor (v1 binary): 136 L0 + 17 L1 + 35 L2
- Binary: /mnt/extfast/ring/archive/biot-savart/build-linux/bs-equivariant
- Source: old worktrees (agent-afadf593/src/features/EquivariantExtractor.cpp, 1740 lines)
- Dynamically linked against glibc/libstdc++ only, no APBS/Eigen/etc .so deps
- Contains: AtomSite precomputation, nearest-ring geometry (14 features), molecular
  graph (BFS, hybridization, bond perception), bond anisotropy decomposition,
  McConnell projections, structural context, APBS E-field integration

### v2 extractor: 94 L0 + 6 L1 + 22 L2
- Binary: build-fresh/extract_equivariant
- Source: src/Main/ExtractEquivariant.cpp (734 lines)
- Produces 94 L0 features per the v10 manifest
- Better architecture (v2 library APIs) but missing 42 L0, 11 L1, 13 L2 features
- Missing: nearest-ring geometry, molecular graph features, bond anisotropy
  decomposition, APBS fields, structural context, McConnell projections

### Missing features (in v1 but not v2):
- near1_*: 14 nearest-ring geometry features (r, rho, z, theta, G_iso, McConnell, B-field components)
- near2_*: 4 second-ring features
- ring_*: 12 multi-ring context features (counts within shells, distances, normal products)
- graph_*: 9 molecular graph features (BFS distance, hybridization, electronegativity)
- struct_*: 11 structural context (H-bond, C=O distance, packing density, aromatic detection)
- bond_aniso_*: 3 decomposed bond anisotropy (backbone/sidechain/aromatic)
- mcconnell_*: 5 McConnell carbonyl features
- disp_*: 2 dispersion contact features (count, sum 1/r^6)
- hbond_*: 6 H-bond features (distance, angle, count, donor/acceptor flags)
- coulomb_*: 4 Coulomb E-field projections (magnitude, ring/bond/backbone projections)
- apbs_*: 3 APBS solvated E-field features
- Additional L1 vectors: 11 extra (ring normals x2, B-field components, E-field, directions)
- Additional L2 tensors: 13 extra (per-type sums, nearest-ring, bond aniso decomposition)

## Training Configuration

### Bilinear Model (primary)
- 84,284 params, 7s/epoch on DGX Spark CUDA
- Best: T0 R²=0.862, T2 R²=0.689, Sign=84.2% (test set, epoch 43)
- System: REPORT sign=91.1%, false confidence=1.87%, coverage=95.2%

### e3nn Model (secondary)
- 6,788,808 params, 31s/epoch on DGX Spark (train on 5090 instead)
- Not yet converged with current data

## System Evaluation Rubric

Three tiers based on classical ring sum magnitude:
- REPORT: |classical| >= threshold_high AND element in report_elements → ML prediction
- PASS: |classical| >= threshold_low → classical prediction only
- SILENT: below threshold_low → no prediction

Tuned thresholds: threshold_high=0.01, threshold_low=0.001
Report elements: {H, C, N, O, S} (all), no pass overrides
Note: pass_override_elements mechanism allows forcing specific elements into PASS tier
even when they qualify for REPORT, useful for elements where ML is unreliable

## Configuration for pyproject.toml

```toml
[paths]
data_dir = "data/features"
delta_dir = "data/deltas"
charge_file = "data/ff14sb_params.dat"
archive_dir = "data/orca_runs"
prep_dir = "data/prepared"
model_dir = "models"

[extraction]
min_ring_current_ppm = 0.01
max_target_T0_ppm = 30.0

[training.bilinear]
hidden_dims = [256, 128, 64]
lr = 1e-3
weight_decay = 1e-4
batch_size = 512
epochs = 200
patience = 30
sign_weight = 0.1
seed = 42

[training.e3nn]
hidden_irreps = "32x0e+8x1e+4x1o+8x2e"
compact_irreps = "16x0e+4x1e+4x1o+8x2e"
lr = 1e-4
weight_decay = 1e-5
batch_size = 512
epochs = 300
patience = 50
warmup_epochs = 20
seed = 42

[evaluation]
threshold_high_sweep = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
threshold_low_sweep = [0.001, 0.005, 0.01, 0.05]
report_element_options = [["H","C"], ["H","C","N"], ["H","C","N","O","S"]]

[tools]
# Empty = rely on PATH; set explicitly if not on PATH
# Note: ToolPaths.cpp reads ~/.nmr_tools.toml; this [tools] section mirrors that format
propka3 = ""
pdb2pqr = ""
pdb4amber = ""
tleap = ""
amberhome = ""
apbs = ""
```

## Known Inconsistency

Ring current intensities differ between modules:
- irreps.py physics filter: [-12.00, -11.28, -12.48, -6.72, -5.16, -5.16, -5.16] (Giessner-Prettre 1969)
- ring_relevance.py: [-25.26, -20.62, -4.34, -4.34, -30.47, -30.47, -30.47] (fitted)
- C++ RingType table: 8 per-type effective intensities including TRP9 (DFT-calibrated, different from both Python sets which have only 7 values)

These affect the physics filter threshold. The training results were produced with
the Giessner-Prettre values in irreps.py. Changing them changes which atoms enter training.
