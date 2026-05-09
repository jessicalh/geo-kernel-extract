# tests/bones/topology — archaeology

These three files (`README.md`, `test_legacy_amber_semantic_api.cpp`,
`test_legacy_amber_semantic_tables.cpp`) were authored 2026-05-07 as
SDET guardrails on the substrate tables. They captured a useful audit
shape — row-level invariants on `LegacyAmberSemanticTables`, public
LookupBy / LookupCap surface, fail-fast on invalid residue/variant —
but predate **Bundle C Slice B** (commit `6ec9bff`, 2026-05-07
afternoon, "rings move to LegacyAmberTopology via RingTopology") which
moved ring data into the substrate proper.

When the slice landed, several of the tests' assertions started
failing against the new substrate shape: 1/6 in the API tests
(`ApplyCapDeltaChangesOnlyDeclaredDeltaFields` — the Bundle C ring
migration touches `ring_position` during cap composition) and 10/14
in the table tests (similar substrate-shape drift across multiple
audit invariants).

The 2026-05-09 post-topology doc-cleanup decision: **retire to
archaeology, do not re-write the assertions in this pass**. The live
substrate coverage is `tests/topology/test_legacy_amber_semantic_integration.cpp`
(the tracked target `topology_semantic_integration_tests`), which loads
real AMBER-prepared protein PDBs through `BuildFromProtonatedPdb`,
runs `Protein::FinalizeConstruction`, and audits the per-atom
`AtomSemanticTable` substrate populated on `LegacyAmberTopology`. That
target is green on current `master` and is the authoritative substrate
guardrail.

The bones'd files preserve the audit shape for future re-write if
anyone wants the fine-grained constexpr-array invariants surfaced
again. Updating the assertions would be a separate slice.

Provenance: see `spec/plan/post-topology-doc-cleanup-2026-05-09.md`
Tier 2.
