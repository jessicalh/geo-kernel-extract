# nmr-extract I/O and Schema — 2026-05-09 consolidation

**Status:** PENDING-WORK consolidation doc. Captures I/O-surface and
schema items migrated from `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md`
§C as part of the 2026-05-09 post-topology doc-cleanup pass. Each
section is the canonical pending design for that I/O surface —
DESIDERATA retires after this consolidation completes; this doc
becomes the home for I/O / schema work.

**Scope.** This doc covers contracts the extractor consumes (input
loaders) and emits (output schemas + metadata conventions). It is
*not* a calculator catalogue — calculators sit in
`spec/PLANNED_CALCULATORS_2026-04-22.md` and
`spec/plan/comprehensive-calculator-inventory-2026-04-30.md`. Three
buckets:

- **Input loaders** (§1 ORCA bond orders, §2 ExperimentalReference,
  §3 ORCA magnetizability).
- **Output schema conventions** (§4 irrep tags, §5 unit / sign
  metadata) — cross-cutting attributes per field.
- **Output schema groups** (§6 validation-bench slices, §7
  TrajectoryResult formalisation) — pre-organised H5 slices for
  downstream consumption.

**Companion docs:** `spec/PLANNED_CALCULATORS_2026-04-22.md`
(calculators that produce or consume the surfaces here — A.5 PCS,
A.7 Green-Kubo, D.1 per-SS CSA there; A.11 BenchmarkBackCalculation
arriving as Amendment 2026-05-09(i) in this same cleanup);
`spec/plan/comprehensive-calculator-inventory-2026-04-30.md`;
`spec/EXTRACTION_SDK.md` + `python/nmr_extract/_catalog.py` (every
output-schema item lands as a corresponding catalog registration);
`OBJECT_MODEL.md`, `spec/CONSTITUTION.md` (architectural
invariants — in particular the "Calculator Shielding Contribution
Contract").

**Architectural-state caveats (2026-05-09):** items below were
authored 2026-04-22 against a pre-substrate pipeline. References to
"per-atom typed identity / categorical record / atom-name
canonicalisation" now resolve through the landed substrate
(`atoms_category_info.npy`, `LegacyAmberTopology::SemanticAt(ai)`,
`NamingApplicator`). Where prose below assumed pre-substrate
framing, an architectural-state note flags the update; the pending
design intent is preserved verbatim from §C.

---

## 1. ORCA Wiberg / Mayer bond orders (input loader)

*Migrated from DESIDERATA §C.1.*

### What
Read Wiberg / Mayer bond orders from ORCA `.out` files (with
`%output Print[P_Mayer] 1` in the input deck). Wherever ORCA is
already run for shielding (the 260-pose calibration set, future
μs-harvester DFT poses), populate the same per-bond order field
that the MOPAC pipeline currently produces, skipping the redundant
MOPAC run on those frames.

### Why
On any frame where DFT shielding is already being computed, ORCA
can also report Wiberg / Mayer per-bond orders cheaply (negligible
cost on top of the shielding calculation). Currently the pipeline
runs MOPAC separately for bond orders even when ORCA output is
available — that's a ~45 s-per-protein MOPAC run avoidable on
DFT-already-computed frames. Reduces calibration-path cost; MOPAC
remains the bond-order source on the non-DFT pipeline (per-frame
trajectory output, where ORCA is not run).

### Implementation sketch
- Add `OrcaBondOrderLoader` parsing `Mayer bond orders larger than
  X` blocks from ORCA `.out` files; emit per-bond-pair (atom_a,
  atom_b, order) records.
- Wire into the existing ORCA driver alongside `OrcaShieldingResult`
  ingestion. On frames where the loader produces results, populate
  `mopac_bond_orders` (or a parallel `orca_bond_orders` if the
  unit / convention differs enough to warrant separate storage)
  from those records instead of triggering a MOPAC run.
- Same NPY contract as `mopac_bond_orders` (3 cols: atom_a,
  atom_b, order). Wrapper class `BondOrders` already exists.

### Dependencies
- ORCA input decks regenerated with `%output Print[P_Mayer] 1`
  (one-line addition; future calibration / harvester runs only).
- Existing `mopac_bond_orders` NPY contract + `BondOrders` wrapper.
- Existing ORCA driver (`OrcaRunLoader`, `OrcaShieldingResult`).

### Cost
Trivial in the DFT calculation itself; loader is parsing, ~0
runtime. Eliminates one MOPAC run (~45 s / 4000-atom protein) on
DFT-paired frames.

### Origin
DESIDERATA §C.1.

### Status
**PENDING.** No `OrcaBondOrderLoader` in `src/`; ORCA input decks
do not currently set `Print[P_Mayer]`. Loader work + input-deck
update + driver wiring all pending. Has no upstream blockers; can
land any time the calibration-path cost becomes a friction point.

---

## 2. ExperimentalReferenceLoader (typed BMRB / RefDB / CSA / PCS / CCR loader)

*Migrated from DESIDERATA §C.2.*

### What
Typed loader that reads experimental reference datasets — BMRB and
RefDB chemical shifts, L6–L12 CSA tensor tables, N5–N7 PCS tables,
M10–M12 CCR datasets — into a typed `ExperimentalReference` object.
The object attaches per ProteinBuildContext or as per-atom lookup
via the categorical-record join (see architectural-state note
below).

### Why
Currently per-atom experimental comparison happens in Python
post-passes (the `learn/` calibration pipeline, the NMR forensics
work for the 10-protein calibration set). Promoting it into a
typed C++ load step means:

- Per-atom residual (predicted – experimental) is computable inside
  the C++ pipeline, not only after Python loads the H5.
- The benchmark-back-calculation calculator (PLANNED_CALCULATORS
  Amendment 2026-05-09(i), formerly DESIDERATA A.11) can attach
  predicted-vs-experimental directly per bench, no Python join.
- The 0.9 validation-benches table in `spec/PHYSICS_FOUNDATIONS.md`
  becomes consumable pipeline output, not an external spreadsheet.

Pairs with §6 below (validation-bench H5 slices) — the loader
provides the experimental side; the bench slice carries the
predicted + experimental + residual triplet in H5.

### Implementation sketch
- New `ExperimentalReference` class holding per-atom typed columns:
  `experimental_shift_ppm`, `experimental_csa_principal_components`,
  `experimental_pcs_ppm`, `experimental_ccr_rate`.
- Per-source loader hierarchy: `BmrbStarLoader`, `RefDbCsvLoader`,
  `CsaTensorTableLoader` (Yao-Bax 2010 / Wylie 2011 / Loth 2005),
  `PcsTableLoader` (Tharayil 2021), `CcrRateLoader` (Loth 2005
  ubiquitin 64-bond rates).
- Each loader emits records keyed by `(residue_index, atom_name)`
  on the canonical AMBER-name vocabulary. The categorical record
  (`atoms_category_info.npy`) is the join key; unmatched rows are
  logged-and-dropped (Huxley discipline per memory
  `feedback_huxley_data_discipline` — never silently translated).
- Attached per ProteinBuildContext rather than per-conformation;
  experimental references are conformation-invariant.

### Architectural-state note (2026-05-09)
DESIDERATA §C.2 was authored before the 2026-05-08 substrate slice
landed. The "per-atom lookup" path it described now resolves through
the categorical record + `LegacyAmberTopology::SemanticAt(ai)`:
join key is the typed semantic, not a string. The loader prose
above reflects that — match against typed semantic identity, log
unmatched rows, never translate silently.

### Dependencies
- Existing per-atom typed identity (`atoms_category_info.npy`,
  `LegacyAmberTopology` substrate). Already landed.
- Per-source file formats; BMRB STAR is the trickiest (multi-line
  records, optional fields). The existing 10-protein forensics
  work has a Python BMRB / RefDB loader that can serve as the
  shape reference (`project_nmr_forensics` memory entry).
- Datasets themselves: BMRB / RefDB downloads, plus the literature
  CSA / PCS / CCR tables digitised once into a versioned data
  directory (similar to `data/calculator_params.toml`).

### Cost
Loader code is moderate (one parser per source + one merge step).
Bigger cost is digitising the literature tables into machine-
readable form once. The 10-protein forensics work
(`project_nmr_forensics`) already did this for BMRB / RefDB shifts;
extending to CSA tensors / PCS / CCR is incremental.

### Origin
DESIDERATA §C.2.

### Status
**PENDING.** No `ExperimentalReference` class in `src/`; the BMRB /
RefDB / literature reads currently happen Python-side. Pairs with
PLANNED_CALCULATORS Amendment 2026-05-09(i)
(BenchmarkBackCalculationResult); the loader is the upstream
half. Can land calculator-by-calculator (start with BMRB / RefDB
shifts which already have a Python implementation to mirror).

---

## 3. ORCA-computed magnetizability extraction (input loader, validation)

*Migrated from DESIDERATA §C.4.*

### What
Read the ORCA-computed molecular magnetizability tensor from `.out`
files into a typed per-conformation object. ORCA computes
magnetizability natively (it's a separate top-level property like
shielding); the extractor currently does not consume it.

### Why
A3 `BulkSusceptibilityAccumulator` (DESIDERATA A.3, pending in
PLANNED_CALCULATORS amendment) accumulates protein-level χ
anisotropy from per-atom kernel sums (McConnell + RingSusceptibility
+ HBond — Babaei 2017 N3 hierarchical approach). Where ORCA
already computed shielding, ORCA can also compute magnetizability
on the same DFT — providing a same-DFT cross-check on the kernel-
accumulator output at the protein scale. Same DFT means the
comparison isolates kernel-vs-DFT residual, not DFT-vs-DFT
methodological drift.

### Implementation sketch
- Add `OrcaMagnetizabilityLoader` parsing ORCA's magnetizability
  output blocks from `.out` files (same driver, same files as
  `OrcaShieldingResult`).
- Emit a per-conformation `OrcaMagnetizabilityResult` with a Mat3
  + SphericalTensor field for the magnetizability tensor.
- One value per protein per conformation, not per-atom (it's a
  global molecular property).

### Dependencies
- ORCA input decks set to compute magnetizability (one-line
  addition to the ORCA `%nmr` block or analogous).
- Existing ORCA driver.
- A3 BulkSusceptibilityAccumulator landed (the consumer; without
  it the magnetizability sits on the H5 unused).

### Cost
Trivial in DFT cost (already running ORCA for shielding). Loader
is parsing.

### Origin
DESIDERATA §C.4.

### Status
**PENDING.** No magnetizability loader; A3 BulkSusceptibilityAccumulator
also pending. The magnetizability on its own without A3 is
reduced-utility — has a DFT value but no kernel-derived comparand.
Land both together when the validation pair becomes a thesis
priority.

---

## 4. Irrep-tagged H5 metadata (output schema convention)

*Migrated from DESIDERATA §C.5.*

### What
Per-tensor-field H5 attribute carrying the e3nn `Irreps`
representation: `0e` for scalar, `1o` for vector, `2e` for rank-2
symmetric traceless, `1x0e+1x1o+1x2e` for the full SphericalTensor.
Stored as an H5 string attribute on each tensor dataset.

### Why
Downstream learners (the equivariant model in `learn/c_equivariant/`,
future GNN consumers, the Python SDK's `ShieldingTensor` /
`EFGTensor` / `VectorField` wrappers) need to know how each
tensor field decomposes into irreps. Currently this is documented
in `spec/EXTRACTION_SDK.md` prose and encoded in the wrapper class
selection in `_catalog.py`. Tagging the H5 itself means consumers
that don't go through the SDK (e.g. an upstream model loading H5
directly via h5py + e3nn) can still type-check.

The Python SDK already encodes this — `ShieldingTensor` carries
`Irreps("1x0e+1x1o+1x2e")`. The H5 attribute is the byte-level
record of the same fact, in machine-readable form, no Python
import required.

### Implementation sketch
- Extend the H5 writer in `fileformat/analysis_file.cpp` to emit
  an `irreps` H5 attribute on each tensor dataset.
- Convention: `0e` (scalar), `1o` (polar vector), `1e` (axial
  vector), `2e` (rank-2 symmetric), and direct-sum forms like
  `1x0e+1x1o+1x2e`.
- Per-dataset emission driven by `_catalog.py` wrapper class:
  `ShieldingTensor` and `EFGTensor` → `1x0e+1x1o+1x2e`,
  `VectorField` → `1x1o`. Add an `irreps: str` field to
  `ArraySpec` so the catalog is the single source of truth.
- Round-trip test in `fileformat/roundtrip_test.cpp`.

### Dependencies
- `fileformat/analysis_file.cpp` — central H5 writer.
- `python/nmr_extract/_catalog.py` — adding the irrep field to
  `ArraySpec`.
- e3nn convention reference (E2 Geiger-Smidt — already covered in
  bibliography).

### Cost
Trivial schema change. Adds ~one string attribute per tensor
dataset; storage cost negligible. Implementation is one writer
extension + one catalog field + a round-trip test.

### Origin
DESIDERATA §C.5.

### Status
**PENDING.** No `irreps` attribute on tensor datasets currently;
SDK encodes the convention via wrapper class selection. Zero-cost
schema change; can land any time. Should land before the
equivariant model (`learn/c_equivariant/`) goes into thesis-figure
production — the tagging is most useful as a self-describing format
for that consumer.

---

## 5. Machine-readable units + sign-convention metadata (output schema)

*Migrated from DESIDERATA §C.6.*

### What
Per-H5-field unit and sign-convention attributes. Currently the
units are documented in `_catalog.py` `description` strings (e.g.
"Coulomb total E-field") and in source comments; the actual unit
("V/Å") and sign convention (e.g. "shielding σ, positive = more
shielded; chemical shift δ = σ_ref − σ, positive = deshielded")
are scattered.

### Why
Eliminates the "is this ppm or ppm × 10⁶" / "is the sign convention
σ or δ" class of bug across downstream Python and viewer consumers.
Anyone reading the H5 directly (h5py, MATLAB, R, Python without
the SDK) gets the units off the dataset, not by tribal knowledge.
Particularly load-bearing for the cross-tool work
(`h5-reader/`, the calibration pipeline, future thesis-figure
scripts) where tribal knowledge does not transfer.

### Implementation sketch
- H5 attributes per dataset:
  - `units: str` — SI unit string (e.g. `"V/Å"`, `"e"`, `"ppm"`,
    `"kJ/mol"`, `"Å²"`).
  - `sign_convention: str` — when sign is conventional, a short
    label naming which convention applies (e.g.
    `"shielding_sigma_positive_more_shielded"`,
    `"chemical_shift_delta_positive_deshielded"`).
  - `description: str` — copy of the `_catalog.py` description
    string (mostly already-present guidance on what the field is).
- Emission driven by extending `ArraySpec` in `_catalog.py` with
  `units: str` and `sign_convention: Optional[str]` fields. The
  `_catalog.py` becomes the single source of truth; the C++ writer
  reads from a parallel C++ table populated from the same canonical
  list.
- Round-trip test in `fileformat/`.

### Dependencies
- `fileformat/analysis_file.cpp`.
- `python/nmr_extract/_catalog.py` (add units / sign columns to
  `ArraySpec`).
- C++-side table mirroring the catalog.
- One-time prose pass — every existing field gets a units +
  sign-convention review.

### Cost
Schema change is trivial; the prose pass is the work. Each of the
~75 catalog entries needs ~one minute of "what's the unit, what's
the sign convention" review. ~2 hours total. Worth it as a
once-and-done.

### Origin
DESIDERATA §C.6.

### Status
**PENDING.** No `units` / `sign_convention` attributes currently.
Pairs with §4 (irrep tagging) — both are zero-cost schema
extensions deserving one combined slice. Lands ahead of the
external consumer phase (other groups' agents reading our H5
without our Python SDK).

---

## 6. Validation-bench H5 slices (output schema, paired with bench-back-calc)

*Migrated from DESIDERATA §C.7.*

### What
Pre-computed H5 groups `/benches/{bench_name}/`, each carrying
per-atom predicted, experimental, and residual values for one
named validation bench. Initial roster:

- `/benches/loth_2005/` — ubiquitin DD/DD CCR rates per N–H bond.
- `/benches/yao_bax_2010/` — α-helix vs β-sheet ¹⁵N CSA principal
  components.
- `/benches/babaei_2017/` — bulk protein χ anisotropy.
- `/benches/tharayil_2021/` — ubiquitin / GB1 PCS per atom.
- `/benches/wylie_2006/` — δ₂₂(¹³C') vs CO···HN distance.

Each group carries: `atom_index`, `prediction`, `experimental`,
`residual`, plus per-bench provenance attributes (paper DOI, dataset
hash). Bench-scoped slicing means thesis-plot scripts open only
what they need.

### Why
Turns the `spec/PHYSICS_FOUNDATIONS.md` §0.9 validation-benches
table into concrete, consumed, tested pipeline output. Without
this, every thesis figure that compares pipeline-vs-experiment
re-derives the join in Python — same pipeline output, same
experimental data, every figure script does the join again.
Pre-computed once at extraction time, the H5 slice is the
contract. Consumers (thesis scripts, paper supplementary tables)
read the residual directly.

Pairs tightly with PLANNED_CALCULATORS Amendment 2026-05-09(i)
(`BenchmarkBackCalculationResult`, formerly DESIDERATA A.11). The
calculator computes the predicted side per bench (per-bench filter
on which atoms apply, per-bench prediction formula); the
ExperimentalReferenceLoader (§2 above) provides the experimental
side; this schema is the H5 layout for the joined result.

### Implementation sketch
- Schema: per-bench H5 group at `/benches/{bench_name}/`. Datasets
  per group: `atom_index: int32`, `residue_id: bytes`,
  `prediction: float64`, `experimental: float64`,
  `residual: float64`, `experimental_error: float64`.
- Per-bench attributes: `paper_doi`, `dataset_version`,
  `prediction_formula`, `experimental_source`.
- Driven by `BenchmarkBackCalculationResult.WriteFeatures`, one
  call per registered bench. Per-bench prediction logic is
  bench-specific (CCR rate formula for Loth 2005, CSA
  principal-component computation for Yao-Bax 2010); the schema
  is uniform.
- Pairs with the per-bench pipeline mode (`nmr_extract --bench
  loth-2005`, DESIDERATA E.6 routed to the diagnostics doc) —
  per-bench mode runs only the calculators that bench needs.

### Dependencies
- §2 ExperimentalReferenceLoader (provides experimental side).
- PLANNED_CALCULATORS Amendment 2026-05-09(i)
  `BenchmarkBackCalculationResult` (the producer).
- Specific calculators per bench: A.1 CSAPrincipalAxis (for
  Yao-Bax / Wylie), A.10 CCRRate (for Loth), A.5 PCS (for
  Tharayil), A.3 BulkSusceptibility (for Babaei). Most are
  PENDING; benches land bench-by-bench as the calculators land.
- §4 / §5 schema metadata (units / sign / irrep) — bench datasets
  carry the same metadata convention as other H5 datasets.

### Cost
Calculator cost dominated by the per-bench prediction logic
(separate slice per calculator). The schema layer is
straightforward — one writer per bench, parallel structure.
Total for the slicing layer once calculators exist: ~hours per
bench.

### Origin
DESIDERATA §C.7.

### Status
**PENDING.** No `/benches/` group emission currently. Lands
incrementally as each upstream calculator lands. Loth 2005 +
Yao-Bax 2010 are the highest-priority benches per the thesis
narrative (W3 handle in `spec/IDENTITY_AND_DYNAMICS_ROLLUP_2026-04-22.md`
section 13.2; α/β CSA split per the per-SS stratification
diagnostic).

---

## 7. TrajectoryResult serialisation formalisation (output schema)

*Migrated from DESIDERATA §C.8.*

### What
Once-and-done schema decisions for trajectory-level NPY / H5
output, decided before any of the trajectory calculators (A.7
GreenKubo, A.8 MemoryKernel, A.10 CCRRate per PLANNED_CALCULATORS
amendments) are implemented. Decisions cover:

- **Autocorrelation lag layout** — `(N_atoms, N_lags, N_components)`
  vs alternates. The `BsT0AutocorrelationTrajectoryResult` exemplar
  (landed) chose one layout; formalise that as the convention.
- **Covariance tensor ordering** — for kernel-kernel covariance
  per rollup 13.5, axis ordering choice across atoms / kernels.
- **Per-field block-SD fields** — A.9 ErgodicityMetric carries a
  block-averaged SD per trajectory field; sibling NPY
  (`<field>_block_sd.npy`) vs same-NPY column vs H5 attribute.
- **Irrep tagging** (§4 above) per trajectory tensor field.
- **Lag-axis units** — picoseconds vs nanoseconds vs frames.
- **Mean / normalisation convention** — already locked for the ACF
  family in `spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md`;
  this slice reaffirms / cross-references.

### Why
Every subsequent trajectory calculator depends on these decisions.
Implementing one calculator and then changing the layout for the
next is a) a multi-calculator rewrite, b) a multi-consumer rewrite
(the SDK, the viewer, any thesis-figure script). Once-and-done
schema is the existing project posture (per memory
`feedback_no_export_pipelines` + the
`spec/TRAJECTORY_WRITE_SURFACE_2026-04-24.md` prior pattern).

### Implementation sketch
- A small `spec/TRAJECTORY_NPY_LAYOUT_2026-XX-YY.md` doc locks each
  decision with rationale. Single-page, decisions enumerated. (Not
  this doc — this doc records the *need* for that decision pass;
  the locking pass is its own slice.)
- Update `python/nmr_extract/_catalog.py` with new layout-encoding
  fields on `ArraySpec` (e.g. `axes: tuple[str, ...]` —
  `("atoms", "lags", "components")`), so consumers can introspect.
- Round-trip tests in `fileformat/` for each layout.
- Cross-reference §4 (irrep) and §5 (units / signs).

### Architectural-state note (2026-05-09)
DESIDERATA §C.8 was authored before the 2026-04-24 trajectory-result
exemplar landed. Some decisions are now de facto locked by the
exemplar and by the ACF-family discipline in
`spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md`. The
formalisation slice should **read the exemplar layouts and lock
them** rather than re-debate. Open questions remain only for
not-yet-implemented categories: covariance ordering, kernel-kernel
cross-correlation ordering, block-SD storage location.

### Dependencies
- Existing `BsT0AutocorrelationTrajectoryResult` exemplar (read
  for current convention).
- `spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md` (ACF-family
  discipline already locked there).
- `spec/TRAJECTORY_RESULT_PLAN_2026-04-24.md` and
  `spec/TRAJECTORY_WRITE_SURFACE_2026-04-24.md` (companion
  trajectory-architecture docs).

### Cost
Decision pass: ~2 hours (read exemplar, enumerate categories,
write the layout doc). Code: writer-side enforcement is part of
the next trajectory calculator slice (whichever lands first); the
layout doc itself is the immediate output.

### Origin
DESIDERATA §C.8.

### Status
**PENDING.** Blocks PLANNED_CALCULATORS amendments for A.7
GreenKubo, A.8 MemoryKernel, A.10 CCRRate (formal-amendments
pending in this same 2026-05-09 cleanup slice). Most likely
landing path: the next trajectory-calculator slice gates itself
on this decision pass, runs it as a 2-hour preamble, then proceeds.

---

## Cross-references

### Items in this doc paired with calculator-side designs

| I/O / schema item | Paired calculator (origin) | Calculator landing home |
|---|---|---|
| §1 ORCA bond orders | None directly; replaces MOPAC bond-order runs | Pre-existing `mopac_bond_orders` consumer; no new calculator |
| §2 ExperimentalReferenceLoader | A.11 BenchmarkBackCalculation | PLANNED_CALCULATORS Amendment 2026-05-09(i) |
| §3 ORCA magnetizability | A.3 BulkSusceptibilityAccumulator | PLANNED_CALCULATORS amendment (pending in this cleanup) |
| §6 Validation-bench H5 slices | A.11 BenchmarkBackCalculation | Same as §2 |
| §7 TrajectoryResult schema | A.7 / A.8 / A.10 trajectory calculators | PLANNED_CALCULATORS amendments (pending in this cleanup) |

§4 (irrep tagging) and §5 (units / signs) are cross-cutting schema
extensions — they do not pair with a specific calculator; every
existing tensor field acquires the new attributes.

### DESIDERATA §C items not in this doc

DESIDERATA §C.3 (lanthanide-tag position + χ-tensor input) is
**paired with PLANNED_CALCULATORS §2** (`PseudocontactShiftResult`,
A.5). It lives with that calculator's spec — a calculator-specific
input is captured alongside the calculator, not in this
cross-cutting I/O doc. See PLANNED_CALCULATORS §2 "Dependencies"
for the lanthanide-tag input.

### Section E items handled elsewhere

`spec/post-topology-doc-cleanup-2026-05-09.md` Tier 2.5 routes
DESIDERATA §E items to: E.2 event menu and E.3 probe-point API →
the sibling `DIAGNOSTICS_AND_WORKFLOWS_*` doc; E.5 calibration-vs-
raw separation → `PATTERNS.md` as a discipline rule; E.6 per-bench
pipeline mode → the diagnostics doc, paired with §6 here.

---

## How this doc is used

Forward I/O-or-schema work touches one of these sections, picks up
its pending design intent, and either implements OR refines the
design here. Updates land as amendments at the bottom of the
relevant section, same convention as PLANNED_CALCULATORS:

- New constraint or decision → append-only amendment under the
  section.
- Implementation lands → flip Status to `LANDED` with commit hash;
  preserve the design prose above.
- Design refinement / scope change → amendment with rationale; the
  original prose above stays as the historical record.

A future session forward-working on (say) BenchmarkBackCalculation
should open PLANNED_CALCULATORS Amendment 2026-05-09(i) for the
calculator-side design plus this doc §2 + §6 for the I/O / schema
halves, and build all three together — calculator + loader +
schema. A benchmark back-calculation calculator without the loader
is incomplete; without the schema slice the H5 has nowhere honest
to write.

## Amendments

*(None yet. Append below, never rewrite above.)*
