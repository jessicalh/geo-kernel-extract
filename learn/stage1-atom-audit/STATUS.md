# Stage 1 Atom Audit Status

Date: 2026-04-25 (status block from 2026-05-13 below)

## 2026-05-13 Update — Topology sidecar landed

Commits `f2781da` + `dc50917`. Closes most of the "C++ Side In
Progress" items in this doc; remaining items below are
re-annotated.

C++ side now emits:
- `atoms_category_info.npy` — extended with 6 new topology-projection
  fields (chain_id, residue_number, insertion_code, parent_atom_index,
  ff_atom_type_string, equivalence_class).
- `residues.npy` (new — codex Residue Table).
- `bonds.npy`, `rings.npy`, `ring_membership.npy` (new).
- `extraction_manifest.json` (new — axis_sizes + axis_alignment +
  topology summary).

Python SDK:
- Strict load: required files + post-load invariant validation
  (axis sizes match, bond endpoints valid, ring refs valid,
  residue atom_count sum == n_atoms). Malformed exports fail loud.
- `ArraySpec` extended with `native_axis` / `irreps` / `units` /
  `sign_convention` / `tensor_rank` / `parity` / `mechanism`,
  populated for all ~108 catalog entries (resolves OI-016).

Still deferred:
- OF3 retention contract (685-fleet stopped per CLAUDE.md, OF3
  not flowing).
- R-side regex-mechanism refactor (downstream consumer of the new
  catalog metadata; OI-120).

## Current Decision State

- NPY remains the stable Stage 1 analysis contract.
- Small HDF5 output may change during the extractor-side rebuild.
- Chain id is not required for this Stage 1 pass.
- WT/ALA atom matching should be symbolic, not nearest-neighbour geometry.
- NMR atom identity is being added extractor-side now.

## C++ Side In Progress

- Add NMR atom identity assignment.
- Export `atom_identity.npy`.
- Export sidecar vocab/metadata.
- Replace mutation compare matching with `residue_index + atom_position`.
- Export `mutation_match.npy`.
- Export `removed_rings.npy`.
- Add `manifest.json`.
- Optionally export ORCA dia/para deltas.

## Analysis Side After Next Run

- Verify atom identity coverage for every atom.
- Verify Stage 1 coarse strata reproduction.
- Verify exact BMRB/RefDB atom-position slicing.
- Audit WT/ALA symbolic match coverage.
- Inspect match-distance distribution as a diagnostic.
- Rebuild Stage 1 around atom identity, raw/normalised kernels, removed-ring source geometry, and ORCA tensor deltas.

## Python Side Added 2026-05-08

- `learn/extract.py --stage1-audit` now requires WT/ALA ORCA
  `_nmr.out` inputs and verifies the new Stage 1 arrays after each
  extraction:
  - `atoms_category_info.npy`
  - WT/mut/delta diamagnetic and paramagnetic shielding tensors.
- `learn/src/actual_physics/bmrb_identity_analysis.py` consumes
  `atoms_category_info.npy` directly.  It emits per-matched-atom rows
  and BMRB-level summaries using BMRB/IUPAC names, residue context,
  mechanical identity, backbone role, pseudoatom/ring/prochirality
  flags, naming provenance, match distance, removed-ring distance, and
  dia/para DFT component cancellation.

Smoke command:

```bash
cd learn/src
python3 -m actual_physics.bmrb_identity_analysis \
    --features /shared/2026Thesis/nmr-shielding/calibration/features/Stage1Smoke_20260508 \
    --output-dir output/actual_physics/bmrb_identity_smoke \
    --max-proteins 1
```

## Stage1BMRB Run 2026-05-09

- Background extraction completed: 720 proteins OK.
- Expected missing ORCA inputs:
  - `A0A075FQU3`: missing WT NMR output.
  - `A0A7C4ZM98`: missing ALA NMR output.
  - `A0A7J2L4W1`: missing WT NMR output.
- One extra failed log entry, `features`, came from the pre-fix job-list
  scan including `calibration/features` as if it were a protein.

## BMRB Atom R Export And R Analysis 2026-05-09

- `bmrb_export_for_r.py` now exports isotropic `T0` targets alongside
  anisotropic `T2` targets for total/dia/para shielding.
- 100-protein iso/T2 export:
  `learn/src/output/actual_physics/stage1_bmrb_r_export_100_iso`
- Per-BMRB atom R analysis:
  `learn/src/output/actual_physics/r_stage1_bmrb_atom_dimensions_100`
- Main tables:
  - `bmrb_atom_ridge.csv`
  - `bmrb_atom_dimensions.csv`
  - `bmrb_atom_dimension_sources.csv`
  - `bmrb_atom_dimension_top_loadings.csv`
