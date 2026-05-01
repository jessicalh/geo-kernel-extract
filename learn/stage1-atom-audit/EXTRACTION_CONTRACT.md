# Stage 1 Atom-Audit Extraction Contract

## NMR Atom Identity

Export `atom_identity.npy` with stable integer columns:

- `atom_index`
- `residue_index`
- `residue_type`
- `prev_residue_type`
- `next_residue_type`
- `element`
- `atom_position`
- `nmr_class`
- `locant`
- `ring_atom_role`
- `ring_membership`
- `symmetry_class_id`
- `parent_atom_index`
- `parent_atom_position`
- `partner_atom_index`
- `prochirality`
- `terminal_flag`

Export sidecar vocab/metadata:

- canonical atom names
- original loaded atom names
- pseudo atom names
- enum names
- column names

## Mutation Compare

Match WT/ALA atoms by:

```text
residue_index + atom_position
```

Keep geometry distance as audit output only.

Export `mutation_match.npy`:

- `matched`
- `mutant_atom_index`
- `match_distance`

## Removed Rings

Export `removed_rings.npy`:

- `removed_ring_id`
- `wt_ring_index`
- `ring_type`
- `source_residue_index`

Existing `ring_geometry.npy` supplies ring center, normal, and radius by
`wt_ring_index`.

## Stage 1 Identity Checks

Verify the new identity fields reproduce the old late Stage 1 strata:

- `CA`
- `C=O`
- `CB`
- backbone `N`
- sidechain `N`
- backbone `O`
- sidechain `O`
- `H`

Also verify exact BMRB/RefDB-facing atom-position slicing:

- backbone: `N`, `HN`, `CA`, `HA`, `C`, `O`
- beta atoms: `CB`, `HB*`
- sidechain heavy atoms by exact canonical atom position
- sidechain hydrogens by exact canonical atom position or pseudo group
- aromatic ring atoms by exact position plus `ring_atom_role`
- methyl and methylene groups by `symmetry_class_id` and `pseudo_atom_name`

## Manifest

Export `manifest.json`:

- schema version
- protein id
- WT root
- ALA root
- extractor version/config
- column names for packed arrays
- enum vocab references

## Optional NPY Additions

- `delta_diamagnetic.npy`
- `delta_paramagnetic.npy`

## NPY-Side Calculator Comparison

If ALA calculator arrays are exported, Stage 1 analysis can compute
calculator deltas from NPY:

```text
calculator_delta[wt_i] = wt_array[wt_i] - ala_array[mutant_atom_index[wt_i]]
```

The extractor does not need per-calculator mutation-delta classes for this pass.

