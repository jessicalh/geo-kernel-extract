# Stage 1 Atom Audit Status

Date: 2026-04-25

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

