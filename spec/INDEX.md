# Document Index — NMR Shielding Tensor Prediction System

**READ THIS FIRST** when starting a new context.

## Reading Order

### Tier 1 — Start here (what the system is)

1. **doc/ARCHITECTURE.md** — structure, provenance paths, pipeline,
   kernel equations, diagrams. Points at code with line numbers.
2. **PATTERNS.md** — what holds the system together, what breaks it.
   Anti-patterns, numerical stability, C++ rules, lessons learned.

### Tier 2 — Open while working (field-level detail)

3. **OBJECT_MODEL.md** — every class, property, type, unit. The
   canonical concrete specification. Code against this.
4. **EXTRACTION_ORDER.md** — dependency graph, feature manifest
   (189 features), what each ConformationResult computes.
5. **CALCULATOR_PARAMETER_API.md** — 93 tuneable parameters with
   equations, defaults, literature references.

### Tier 3 — Design decisions (ask before diving in)

These documents are deep. Reading them burns significant context
for significant gain — they ground you in the physics and sign
conventions that every calculator depends on. Ask the user whether
the task requires this depth. The user is aware of these tiers and
may direct you to dive deep on one or more — if so, read thoroughly.
Skimming these is worse than not reading them.

6. **spec/CONSTITUTION.md** — supreme constraint. Sign conventions,
   inviolable rules, anti-simplification rules, minimum output
   contracts per calculator.
7. **spec/MATHS_GOALS.md** — the three mathematical pillars: T0
   parameter corrections, T2 angular residual, classical calculation
   correctness. Validation checklist.
8. **GEOMETRIC_KERNEL_CATALOGUE.md** — full derivations for the 8
   classical + 2 MOPAC-derived calculators. McConnell tensor from
   first principles.

### When implementing a calculator

Read an existing calculator first (e.g., McConnellResult.cpp).
You must:
- Produce full rank-2 tensor output (Mat3 + SphericalTensor)
- Hold a KernelFilterSet and call AcceptAll() before every evaluation
- Record GeometryChoices via GeometryChoiceBuilder during Compute(),
  showing every inclusion, exclusion, and triggered event with
  entities, roles, outcomes, and named numbers with units
- Document the unit chain in comments
- Verify sign with a known physical scenario

See PATTERNS.md "Numerical Stability" section.

### When adding or changing WriteFeatures output

Validate with `./smoke_tests` (~80s), not individual calculator tests.
See **spec/TEST_FRAMEWORK.md** for the full test tier decision tree.

Any new NPY file must be registered in the Python reader SDK so
consumers can load it.  See **spec/EXTRACTION_SDK.md** for the
design and **python/API.md** for the full column-level reference.

The reader lives at `python/nmr_extract/`.  The format contract is
`python/nmr_extract/_catalog.py` — one `ArraySpec` per NPY file.
If you add a column to `ring_contributions.npy` or create a new
array, update `_catalog.py` and the corresponding wrapper class,
then run `python -m pytest python/tests/` to verify.

## All Documents

### Top-level
- **OBJECT_MODEL.md** — canonical concrete specification
- **EXTRACTION_ORDER.md** — dependency graph + feature manifest
- **CALCULATOR_PARAMETER_API.md** — model-to-calculator parameter interface
- **PATTERNS.md** — what holds the system together (living doc)
- **GEOMETRIC_KERNEL_CATALOGUE.md** — calculator physics, derivations

### doc/
- **ARCHITECTURE.md** — system map with diagrams
- **diagrams/** — mermaid sources + rendered SVG/PNG
- **generated/doxygen/** — full class UML, collaboration graphs

### spec/
- **INDEX.md** — this file
- **CONSTITUTION.md** — supreme constraint
- **MATHS_GOALS.md** — mathematical validation plan
- **USE_CASES.md** — the 5 use cases (--pdb, --orca, --mutant, --trajectory, --trajectory --analysis)
- **GEOMETRY_CHOICE_BRIEF.md** — GeometryChoice recording spec
- **DEPENDENCIES.md** — external library list
- **DIRECTORY_SET.md** — directory structure (historical)
- **TEST_FRAMEWORK.md** — test tiers, what to run when, smoke test design
- **ENSEMBLE_MODEL.md** — GromacsProtein trajectory-based model (revised 2026-04-12, replaces old EnsembleConformation design)
- **TRAJECTORY_EXTRACTION.md** — full-system trajectory path: explicit solvent E-field, water packing, ion field calculators (design, 2026-04-12)
- **TIMING_4876_ATOMS.md** — measured per-calculator times on 4876-atom protein (2026-04-12)
- **OUTSTANDING_GROMACS_PATH.md** — open items for trajectory extraction path (GeometryChoice, KernelFilterSet, TOML, SDK tests)
- **POLARISABILITY_ROADMAP_2026-04-13.md** — charge polarisation proxy approaches: SASA normal, HydrationGeometryResult, water-embedded AIMNet2, EEQ, E-field variance. What was rejected and why.
- **ANALYSIS_TRAJECTORY_2026-04-14.md** — analysis trajectory mode: exhaustive per-frame H5, per-ring K=6, ridge+MLP predictions, projections. 10-protein workspace design (tentative, evolving).
- **IDENTITY_AND_DYNAMICS_ROLLUP_2026-04-22.md** — tentative rollup spec: NmrAtomIdentity per-atom invariants, TrajectoryResult continuous summaries, threshold-parameterised event menu with hookable C++/Python forms. Living document; re-read before acting.
- **PHYSICS_FOUNDATIONS.md** — IN PROGRESS. Properly-referenced physics theory for every mechanism the extractor touches. Session 0 landscape underway 2026-04-22; three drafting sessions follow. Must exist before rollup implementation per IDENTITY_AND_DYNAMICS_ROLLUP section 13.11.
- **MICROSECOND_MD_HARVESTER_2026-04-22.md** — DECIDED design, NOT YET BUILT. Pose-harvester scanner (`nmr_extract --harvest-poses`) + AI-curator DFT submission for μs-MD extension of the 10 calibration proteins during a 41-day 5090 rental window. Au-courant ~2026-04-24 when current 25 ns fleet run completes; queued after Session 1 drafting.
- **PLANNED_CALCULATORS_2026-04-22.md** — NOTE, not spec. Five calculator / diagnostic / model-architecture ideas surfaced during Session 0 literature pass: GreenKuboSpectralDensity (J(ω) from σ-autocorrelation), PseudocontactShift (same K_ab as RingSusceptibility with external χ tensor), per-SS CSA stratification diagnostic, Lie-group GP on SE(3), FNO for volumetric shielding field. None committed; captured so ideas survive between sessions.
- **NMR_EXTRACT_DESIDERATA_2026-04-22.md** — NOTE, not spec. Consolidated library-wide pass of every calculator / variation / input-output surface / diagnostic / architectural idea surfaced during Session 0 — organised, not triaged. Sections A-F: new ConformationResults (A1-A11 + 2 cross-ref to PLANNED), variations on existing kernels (B1-B7), I/O surfaces (C1-C8), diagnostics (D1-D3), architectural (E1-E6), scope boundaries (F incl. explicit xTB exclusion). Starting place before the user organises the work against use cases and thesis narrative.
- **pending_include_trajectory_scope_2026-04-22.md** — was WIP_OBJECT_MODEL.md during the design pass; renamed after the trajectory-scope shape landed in src/ and folded (tight form) into OBJECT_MODEL.md + PATTERNS.md §§13-18. Holds working-note material that doesn't belong in the control docs: rejected alternatives (Appendix B — TrajectoryBond first-class store), design-option discussions, pending appendices awaiting user review (NmrAtomIdentity §2 + Appendix A, full TrajectoryResult catalog Appendix F, H5 metadata schema §7), and the anti-pattern warnings §3 that the `feedback_trajectory_scope_gotchas` memory entry references. **Not authoritative for anything landed** — the code + OBJECT_MODEL.md + PATTERNS.md win. Read only when activating a pending-review section.
- **TRAJECTORY_LANDING_STATE_2026-04-23.md** + **TRAJECTORY_REFACTOR_GAPS_2026-04-23.md** — session-dated records of the first landing push. Historical; useful for tracing how the trajectory-scope refactor got from WIP to landed.
- **meta-docs-review/** — 2026-04-03 documentation audit artifacts

### ui/
- **CLAUDE.md** — UI-specific rules, library API quick-reference,
  current state, what you own vs read-only
- **UI_ROADMAP.md** — forward-looking per-calculator visualization plan
  (stable — nothing here is needed for correctness)
- **GEOMETRY_CHOICE_UI_NOTES_20260407.md** — library-to-UI handoff
  for GeometryChoice display (what each calculator records, UI ideas)
- **src/REST_INTERFACE_SPEC.md** — all 20+ REST commands with parameters
- **doc/generated/doxygen/** — viewer class diagrams, collaboration graphs

### python/ — Extraction Reader SDK
- **spec/EXTRACTION_SDK.md** — design intent, trust boundary, what it covers
- **python/API.md** — consumer reference: every type, field, column
- **python/nmr_extract/_catalog.py** — format contract (one entry per NPY)
- **python/doc/** — Sphinx autodoc source (build with `sphinx-build`)
- **python/tests/test_load.py** — 79 tests against real extractions

### Historical (learn/bones/)
Session notes, early specs, analysis passes, resolved feedback.
Reference only.

### Source
- **src/** — 74 headers, 62 .cpp files
- **ui/src/** — 14 viewer source files (Qt6/VTK)
- **tests/** — GTest suite (287 tests, 40 test files, ~2.5 hours with MOPAC)
- **extern/** — header-only libraries (nanoflann, sphericart)

### External
- **/shared/2026Thesis/consolidated/** — 723 WT+ALA protein pairs
- **tests/data/fleet** — GROMACS CHARMM36m ensemble

## Current Status (2026-04-15)

20 calculators (8 classical, 2 MOPAC-derived, AIMNet2, SASA,
WaterField, HydrationShell, GromacsEnergy, HydrationGeometry, EEQ,
BondedEnergy, plus foundation: Geometry, SpatialIndex, Enrichment).
DSSP extended with 8-class SS, H-bond energies, chi1-4. Calibration
settled at R²=0.818 (per-element ridge, 55 kernels).

**Single TPR parse (2026-04-14):** FullSystemReader::ReadTopology
parses the TPR once — atom ranges, bonded parameters, and protein
construction all from one `read_tpx_state` call. Previously three
independent parses. `--trajectory DIR` takes a protein directory
(must contain md.tpr, md.xtc, md.edr — all required). EDR is no
longer optional. BuildFromTrajectory(dir_path) derives all file
paths from the directory convention.

Two-pass trajectory streaming: GromacsProtein (adapter +
accumulators) + GromacsFrameHandler (streaming XTC reader, PBC fix,
frame lifecycle) + GromacsRunContext (bonded params, EDR energy
frames, cursor state). Per-calculator timing via OperationLog::Scope
on every Compute call. 75 NPY arrays registered in SDK catalog.

**Analysis trajectory mode (2026-04-14):** AnalysisWriter harvests
all per-atom and per-residue data from each sampled frame, writes
`{protein_id}_analysis.h5` at end. 21 H5 groups: raw physics
(ring_current, efg, bond_aniso, quadrupole, dispersion, hbond,
sasa, water, charges, aimnet2_embedding, per_ring K=6), per-residue
dihedrals/ (phi, psi, omega, chi1-4 with cos/sin), dssp/ (ss8,
hbond_energy), bonded_energy/ (per-atom CHARMM36m bond/angle/UB/
proper/improper/CMAP decomposition), energy/ (42 EDR terms).
PDB snapshots at ~1ns intervals — these ARE the ORCA DFT inputs.
Predictions/ and projections/ groups pending (need Python model export).

**Adversarial validation (2026-04-15):** 1CBH_192 analysis H5
validated. NPY-vs-H5 bit-identical at 0ns and 12ns. No NaN/Inf.
Ring current signs correct. Bonded energy sums exact. Charge
conservation to machine precision. Dihedrals cos²+sin²=1 to
machine epsilon. MolProbity: 0 twisted peptides across all PDBs,
Ramachandran outliers typical for MD snapshots (not minimised).

**Fleet extraction (2026-04-15):** 10 proteins at
`/shared/2026Thesis/fleet_calibration-working/`. ~25 min per protein at
~4000 atoms (APBS dominates at ~2.7s/frame). Two proteins run in
parallel. ORCA launched on PDB snapshots.

**AIMNet2 aim embedding:** float32 (native torch precision) through
the full path: ConformationAtom, NPY, H5. Fixed 2026-04-15 — prior
H5 files have float64 (upshifted).
