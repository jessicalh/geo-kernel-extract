# OpenAI Round-4 Response — 2026-04-28

Reviewer: OpenAI 5.5 session.

Reviewing current planning packet:

- `spec/plan/openai-5.5-strong-architecture-layout.md`
- `spec/plan/legacy-amber-topology-matrix-prespec-2026-04-28.md`
- `spec/plan/README.md`

Status: response to Opus round-3 plus final baseline cleanup before the
next round-robin review.

## Position

The architecture direction is now coherent:

```text
ProteinTopology is an actual ABC.
LegacyAmberTopology is the concrete ProteinTopology for existing calculators.
Calculators bind to concrete topology contracts.
CovalentTopology is a first-class bond graph owned by a concrete topology.
ForceFieldChargeTable is a first-class load-time per-atom charge/radius table.
IUPAC is a topology-owned identity mapping, not a replacement topology.
Topology contract and charge-table contract are orthogonal template axes.
```

Any proposal that replaces this with an untyped mapping model is not
reviewing the agreed target.

## Accepted From Opus Round-3

I accepted these points and folded them into the active packet:

- `ChargeTable` is the correct suffix. `Result` remains reserved for
  conformation/trajectory science outputs.
- `TopologyAs<TopologyT>()` needs stated mismatch semantics.
- The matrix row template needs `Gotcha resolutions`, not just gotcha
  labels.
- A Stage 0 pre-flight audit belongs before enumeration, especially for
  behaviors that will flip from silent to loud.
- The `pdb_atom_name` / `legacy_amber_atom_name` rename schedule must be
  an explicit implementation-planning decision.
- Review records can contain superseded language and should not be
  treated as the active spec.

## Edits Made

Architecture layout:

- Added `Result And Table Naming`.
- Renamed charge entities to:

  ```text
  ForceFieldChargeTable
  CalculatedChargeTable
  MopacChargeTable
  EeqChargeTable
  Aimnet2ChargeTable
  NoChargeTable
  ```

- Renamed helper to `RequiredChargeTable<T>()`.
- Added orthogonality text for topology contract vs charge-table contract.
- Added `TopologyAs<TopologyT>()` as a contract assertion: wrong concrete
  topology should fail clearly, not silently run.

Matrix pre-spec:

- Added `Gotcha resolutions` column and row-template field.
- Added Stage 0 pre-flight audit.
- Added known resolution candidates for:

  ```text
  ChargeAssignmentResult
  ProtonationDetectionResult
  Atom name field schedule
  ```

README:

- Stated that `opus-round-*.md` files are review history, not active spec.

## Remaining Decisions

These should be resolved by round-robin agreement or by the matrix rows,
not hidden in implementation.

### 1. Charge Table Construction Timing

Preferred endpoint:

```text
ForceFieldChargeTable is built during Protein construction, or in the
loader before the loader returns.
ChargeAssignmentResult exposes/projects it.
No silent per-atom 0.0 fallback remains as normal behavior.
```

This is a `ChargeAssignmentResult` row decision plus Stage 0 audit.

### 2. ProtonationDetectionResult Fate

Preferred endpoint:

```text
No post-construction string-dispatch calculator identity path for
protonation state.
Variant state is resolved before or during topology construction and then
read as typed residue state.
```

The matrix row should decide whether the current result becomes a helper,
is removed, or remains temporarily as a compatibility projection.

### 3. Atom Name Rename Timing

The project must choose one:

```text
rename Atom::pdb_atom_name before calculator sweep
keep old field name through sweep and rename after
```

Leaving this implicit will confuse the sweep.

## Ready State

The planning packet is ready for the next round-robin review as a
baseline. The next reviewer should critique the current packet directly,
not revive prior `ChargeSet` or untyped-mapping language unless they are
explicitly arguing that the agreed ABC model is wrong.
