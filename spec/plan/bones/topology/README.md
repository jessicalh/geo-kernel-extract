# Topology Semantic Substrate Tests

These tests are the SDET guardrail for the LegacyAmberTopology semantic
substrate. They are intentionally isolated from the broader unit and structure
test buckets because they validate a generated taxonomy, not molecule loading.

## Targets

- `topology_semantic_api_tests`
  - Links `libnmr_shielding`.
  - Exercises only the public runtime surface in
    `src/generated/LegacyAmberSemanticTables.h`.
  - Guards `LookupBy`, `LookupCap`, invalid-residue and invalid-variant
    fail-fast behavior, variant chemistry dispatch, and `ApplyCapDelta`.

- `topology_semantic_table_tests`
  - Does not link `libnmr_shielding`.
  - Includes `src/generated/LegacyAmberSemanticTables.cpp` into the test
    translation unit so the internal constexpr arrays can be audited without
    making them production API.
  - Guards row-level invariants: element population, mechanical-identity
    uniqueness except explicit equivalent-H sets, cap chemistry, backbone-role
    conventions, variant inventory, variant charge override ordering, and CCD
    standard-residue element scope.

## What This Catches

- A standard residue or variant routed to the wrong table.
- Invalid variant indices silently falling back to the base residue.
- Unknown residues returning a default row instead of `nullptr`.
- Cap composition clobbering identity fields or RDKit-derived chain fields.
- Terminal caps drifting from the residue-reference Section 4 chemistry.
- Variant additions/removals drifting from the AMBER protonation table.
- Parent synthesis charge overrides not being patched back in variants.
- `Element::Unknown` appearing in emitted standard/variant/cap rows.
- A future CCD standard-residue element outside runtime `Element` scope.

## Corpus Hook

These tests are deterministic invariants over the generated tables. Corpus
grounding should build on top of the same public API: generate or ingest AMBER
fixtures, compute each atom's `AtomMechanicalIdentity`, call `LookupBy` and
`LookupCap`, then persist an audit artifact with residue, variant index,
terminal state, atom name, identity tuple, table hit, charge, polar-H kind, and
cap-composition outcome.
