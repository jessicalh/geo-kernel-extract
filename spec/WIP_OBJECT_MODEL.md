# WIP Object Model — Trajectory Scope (2026-04-22)

**Status:** Work in progress. This is the source of truth for the
trajectory-scope object model + attending patterns during this design pass.
When it settles, sections here fold into `OBJECT_MODEL.md` and `PATTERNS.md`;
this file is then retired. Until then, header notes at the superseded
sections of those two documents point readers here.

**Purpose.** The library has three model scopes: per-protein (Protein),
per-conformation (ProteinConformation + ConformationAtom + ConformationResult),
and per-trajectory. The third is currently partially realised — `GromacsProtein`
is the adapter that ended up holding per-atom trajectory state, accumulation
is baked into a hardcoded `GromacsProteinAtom` struct as a flat bag of Welford
instances, and serialisation flows through `AnalysisWriter` as a buffered
holder outside the model. This document designs the third scope properly as
a peer of the first two, extending the same discipline (typed entities,
modular calculators, self-serialising results, singleton-per-type) into the
multi-frame case.

**Scope of this pass:** (1) define `TrajectoryProtein` and `TrajectoryAtom`
as canonical physical-world entities, (2) define `TrajectoryResult` as the
modular calculator unit at trajectory scope, (3) define `Trajectory` as the
process entity that runs a traversal, (4) introduce `NmrAtomIdentity` on
Protein (topology-like typed identity), (5) name the explicit renames,
dissolutions, and bones-migrations. Code changes follow the design; this
document comes first.

---

**NON-REGRESSION CONTRACT (load-bearing, preservation commitment).**
This rollup is trajectory-scope. **The static, non-trajectory extraction
paths must continue to work unchanged through the activation.** Specifically:

1. **Static single-PDB load** (protonated or unprotonated) with classical
   kernels: unchanged. `Protein` built from a single structure;
   `ProteinConformation` is the single conformation; `OperationRunner::Run`
   sequences the per-conformation calculators; NPY / H5 output emitted
   per existing conventions.
2. **Static single-PDB + ORCA shielding load** (calibration path):
   unchanged. `OrcaShieldingResult::Compute` loads DFT tensors for
   the single conformation; stored on `ConformationAtom.orca_shielding_*`;
   consumed by calibration scripts in `learn/`.
3. **Static WT + ALA mutant pair with ORCA + MutationDeltaResult**:
   unchanged. Two Protein instances, one conformation each, ORCA loaded
   on both, `OperationRunner::RunMutantComparison` produces
   `MutationDeltaResult` on the WT conformation. This is the 720-pair
   calibration pattern; it is the heart of the thesis and must not
   regress.

These paths use `Protein + ProteinConformation + ConformationResult +
OperationRunner::Run` (and `::RunMutantComparison`) as they do today.
The only whole-library change that touches them is `NmrAtomIdentity`
on Protein (§2) — which is additive and does not change what any
existing calculator reads or writes.

No trajectory-scope additions — TrajectoryProtein, TrajectoryAtom,
TrajectoryResult, Trajectory, RunConfiguration, RunContext — change
the static paths. Static paths don't use any of these. `Protein` and
`ProteinConformation` keep their existing shape; `OperationRunner`
keeps its existing two entry points.

**Test enforcement: concrete regression sentinels.** The "existing
static-path test suite must pass" commitment resolves to specific test
categories that exercise the static paths. These must remain green
throughout every stage of the rollup (see Appendix H for the stage
plan); any failure is a signal that a trajectory-scope change has
leaked into non-trajectory scope and must be contained.

1. **Static single-PDB load tests.** Build Protein from PDB input —
   one per residue type plus variants (HIE/HID/HIP, CYS/CYX), at
   minimum. Verify atom count, bond detection output, ring detection,
   residue assignment. `tests/test_protein_*` or equivalent current
   test names.
2. **Conformation construction tests.** ProteinConformation built
   from positions, ConformationAtom field initialisation, spatial-
   index build, neighbourhood queries. `tests/test_protein_conformation*`.
3. **OperationRunner::Run tests.** Per-calculator attach order,
   skip-flag behaviour, RunOptions dispatch, dependency validation.
   `tests/test_operation_runner*`.
4. **Per-calculator smoke tests.** Each of the ~20 classical
   calculators has an individual smoke test (e.g.,
   `tests/test_biot_savart_result*`, `tests/test_mcconnell_result*`,
   `tests/test_coulomb_result*`, ...). Each computes against a
   canonical input and asserts output fields match expected values.
5. **OrcaShieldingResult tests.** Load ORCA shielding from a
   reference file; verify per-atom tensor parsing, diamagnetic /
   paramagnetic decomposition, atom-matching during load.
6. **OperationRunner::RunMutantComparison tests.** Per the 720-pair
   calibration pattern: WT conformation + mutant conformation + ORCA
   on both + MutationDeltaResult production. These specifically
   exercise the static-mutant-pair path in the non-regression
   contract above.
7. **NPY output tests.** Per-calculator WriteFeatures + Python SDK
   round-trip (`python -m pytest python/tests/test_load.py`).
   Downstream consumers' primary contract is the NPY/H5 schema.

**Enforcement mechanism.** Run the full test suite at every stage
boundary (Appendix H). A regression in any of the above categories
blocks progression to the next stage — no stage progresses until
the sentinel tests are green again.

**Companion docs:**
- `OBJECT_MODEL.md` — the canonical object model. Superseded sections
  carry header notes pointing here until this document folds in.
- `PATTERNS.md` — the pattern language holding the system together.
  A header note in `PATTERNS.md` points here for trajectory-scope patterns.
- `spec/CONSTITUTION.md` — supreme constraint on the code. Not
  superseded; trajectory-scope entities inherit its rules (one sign
  convention, typed identity, no string dispatch, full tensors).
- `src/ConformationResult.h`, `src/OperationRunner.h` — the
  conformation-scope pattern being extended here. Read these as the
  reference shape.
- `spec/ChangesRequiredBeforeProductionH5Run.md` — the coordinated
  pre-fleet activation discipline. This rollup lands under the same
  discipline: single activation, fleet-wide.

---

## 1 — Entity categories: physical-world model and process model

The library has two categories of model entity. Both are real; both carry
typed fields and methods; both participate in the object model. The
distinction is what each one models.

**Physical-world model entities** describe the protein and the things we
compute about it. They are the library's representation of the molecule and
its properties. Rules: permanence and identity — once an entity exists and
is attached, it stays; once a field is written by its owning computation,
it stays; identity is typed and stable; one writer per field.

**Process model entities** describe the workflow we do. A trajectory
traversal is a process with a beginning, middle, and end. The thing that
drives it exists because we're doing computational work, not because there's
a protein with that structure. Rules: lifecycle — create, run, complete;
may hold transient state during the active phase and typed records after;
may persist as a queryable run-record.

**Both categories coexist and can refer to each other.** A Trajectory
(process) references a TrajectoryProtein (physical-world) it traversed.
A TrajectoryResult (physical-world, modular calculator) is attached by a
Trajectory during its run. Keeping the rules distinct while allowing
mutual reference is what makes the model readable.

### Inventory

**Physical-world model:**
- `Protein` — topology and identity. What the molecule IS. Unchanged in
  role; extended with `NmrAtomIdentity` per §8. Has `Atom` objects
  (identity-only).
- `ProteinConformation` — one geometric instance of a Protein. Unchanged.
  Has `ConformationAtom` objects (per-conformation computed data store).
- `ConformationResult` — modular calculator at conformation scope.
  Unchanged.
- `TrajectoryProtein` — a protein in trajectory context. Wraps a Protein,
  adds the per-atom running buffer and the attached TrajectoryResults.
  **New in this pass; replaces the existing `GromacsProtein` role.** Has
  `TrajectoryAtom` objects (per-trajectory computed data store).
- `TrajectoryResult` — modular calculator at trajectory scope. Parallel
  to ConformationResult in role but with multi-phase compute (per frame
  + finalize). **New in this pass.**
- `NmrAtomIdentity` — typed identity lenses attached to `Protein`.
  Topology-level, geometry-independent. **New in this pass.**

**Process model:**
- `Trajectory` — the process entity representing one traversal run.
  Owns source paths, holds frame metadata and selection records during
  and after the run, drives the traversal via `Run()`. **New in this
  pass.**
- `GromacsFrameHandler` — per-frame mechanics. Unchanged in role; fed
  by Trajectory, feeds TrajectoryResults.
- `RunConfiguration` — typed description of a specific run shape.
  Three concrete configurations: `ScanForDftPointSet`,
  `PerFrameExtractionSet`, `FullFatFrameExtraction`. **New in this
  pass as a first-class typed object.**
- `RunContext` — caller-constructed object carrying configuration and
  extras for a run. Not a struct bag; a class with constructor and
  methods. **New in this pass; absorbs some of the role the current
  `GromacsRunContext` played, with ownership redistributed.**

---

## 2 — Protein + Atom (extended): `NmrAtomIdentity`

**⚠ STATUS: PROPOSAL PENDING USER REVIEW (flagged 2026-04-22 pm).**

The *need* for NMR-aware typed atom identity is not in question —
user confirmed 2026-04-22: *"I do and did need NMR IUPAC and
OpenBabel."* Both the IUPAC-typed atom identity and the OpenBabel
hybridisation source are wanted in the library, independent of this
rollup.

What is pending is the *specific design* of NmrAtomIdentity as
described in this section (IupacAtomPosition enum inventory, NmrClass
taxonomy, the semantic-lens field set, the orthogonal-index API), plus
the generator-script approach in Appendix A and the OpenBabel migration
plan in Appendix E. User flagged: *"I did not design this with you you
parachuted it in with some complexity. I probably agree but haven't
yet."*

Read §2 + Appendix A + Appendix E as a **proposal to discuss**, not as
settled design. The shape of the NmrClass taxonomy, the exact enum
inventory, the decision to generate from TOML vs hand-write, the
EnrichmentResult-becomes-shim migration — each is open for user
revision. The §2 content stays in place as the proposal's current
form; it will evolve during user review before implementation begins.

Do not implement §2 / Appendix A / Appendix E until user sign-off.

---

`Protein` is unchanged in its existing role — topology, residue list,
atom list, ring list, bond list, ProteinBuildContext. It gets one typed
extension in this pass: every atom gains `NmrAtomIdentity`.

### Why NmrAtomIdentity lives on Protein: topology-invariant rule

Per `spec/CONSTITUTION.md` (section "The Object Model"): *"Protein is
identity and topology only. Per-atom geometry lives on
ProteinConformation."* IUPAC atom position, NmrClass, MethylGroup
pairings, RingAtomRole, ChiParticipation — all determined by
(residue type, bond connectivity, IUPAC convention), all
geometry-independent. They are protein topology. They belong on
Protein.

`ConformationAtom` today carries categorical booleans (`is_amide_H`,
`is_methyl`, `is_aromatic_H`, `parent_is_sp2`, etc.) written by
`EnrichmentResult` per-frame. These are topology-invariants attached
to the per-conformation data store — a longstanding violation of the
rule above. Per this revision's OpenBabel decision (§"OpenBabel and
EnrichmentResult" below + Appendix E), those booleans stay on
ConformationAtom for migration-cost reasons. `NmrAtomIdentity` on
Protein provides the correct architectural home for new code; typed
queries through it (`protein.AtomIdentity(i).nmr_class ==
NmrClass::BackboneHN`) supersede the flat boolean reads over time as
callers migrate.

### Why the library itself needs typed identity (not just downstream)

The typed identity is load-bearing for physics work the library
performs directly, independently of any SDK-side consumption. Specific
literature-anchored uses, from the bibliography pass:

- **Per-atom-type / per-stratum calibration pooling (Stage 1, settled
  at R² = 0.818 per `CLAUDE.md`).** Today stratifies by AMBER atom
  name within each element via string-filtering. `NmrClass` slicing
  via `Protein::AtomsByNmrClass(cls)` is the typed version of what
  calibration already does.
- **Havlin-Oldfield (B12, B13): Cα CSA principal values individually
  sensitive to (φ, ψ, χ₁) across six amino acids, extended to all 20
  in the 2001 paper with carbonyl context.** `ChiParticipation`
  bitfield identifies the χₖ-defining atoms for dihedral-dependent
  analysis; `NmrClass::BackboneCA` + `NmrClass::BackboneCO_C` give
  the stratification keys. The carbonyl-context extension means
  residue-category stratification also matters — supplied by
  `ResidueCategory` broadcast per atom.
- **Yao-Bax 2010 (L7): α-helix 15N CSA ≈ −173 ± 7 ppm vs β-sheet
  ≈ −162 ± 6 ppm — 11 ppm secondary-structure split.** Per-stratum
  diagnostic needs `AtomsByNmrClass(BackboneN)` typed slice crossed
  with DSSP secondary structure; without the typed selection, this
  becomes string-matching on atom_name columns.
- **Loth-Pelupessy-Bodenhausen 2005 (L12): ubiquitin amide CSA from
  CCR, 64 amide bonds, all three nuclei (C', 15N, HN).** Per-amide
  triplet identification via NmrAtomIdentity primary tuples
  (`(AminoAcid, IupacAtomPosition)`); backbone-triplet construction
  via `NmrClass::BackboneCO_C` + `BackboneN` + `BackboneHN` typed
  selection.
- **Kasinath-Wand (M7): methyl entropy meter with pro-R/pro-S
  resolution.** VAL γ₁↔γ₂, LEU δ₁↔δ₂, ILE γ₂ / δ₁, THR γ₂, MET ε.
  `MethylGroup` enum + `methyl_partner_atom_index` preserve pro-R /
  pro-S at the typed level; flat per-atom booleans cannot.
- **Brender-Taylor-Ramamoorthy 2001 (L11): amide 15N CSA tensor
  orientation relative to N-H bond vector.** Per-atom NmrClass
  (`BackboneN` + `BackboneHN` pairing) + typed N-H bond geometry
  gives the frame-construction input without string matching.
- **Pauling 1979 (N2), Babaei et al. 2017 (N3): per-peptide and bulk
  χ anisotropy.** `ResidueCategory` (Aromatic contributes differently
  by a significant margin) broadcast per-atom gives O(1) stratification
  during kernel summation in the bulk-susceptibility accumulator.

These are library-internal consumers — kernels, diagnostics, per-stratum
calibration, validation-bench back-calculation. The SDK-side
consumption (BMRB shift binding via `IupacAtomPosition` as primary
key) is a parallel consumer downstream; both matter, both read the
same typed identity, neither subsumes the other. NmrAtomIdentity is
infrastructure the library needs to do its own work; the fact that
the SDK uses the same identity for external-data integration is a
consequence of typing invariants correctly, not the primary case.

### Scope framing

`NmrAtomIdentity` is a whole-library change. The fleet activation
window is where it lands because (a) coordinated library changes
happen there, (b) trajectory-scope TrajectoryResults depend on typed
per-atom slicing for their attachment and query patterns, (c) the
existing calibration work and the bibliography-anchored validation
benches both depend on the typed identity being available in the
library, and (d) the architectural rule (topology invariants on
Protein) has been violated longer than necessary.

If the bundle decomposes for any reason, `NmrAtomIdentity` ships as
a pre-bundle (additive, testable independently) rather than staying
coupled to the trajectory refactor.

### What `NmrAtomIdentity` is

Topology-like typed identity for an atom: the things about an atom that
are determined by (residue type, position in residue, bond connectivity)
and do not change when the protein moves. It's called "NmrAtomIdentity"
rather than "topology" because it extends beyond pure bond connectivity
to include IUPAC position, NMR-specific semantic categorisation, and
per-residue invariants broadcast per-atom for fast slicing.

Lives on `Protein` (not on `ProteinConformation` or `TrajectoryProtein`).
Populated once during Protein construction via an EnrichmentResult-adjacent
process. Const thereafter. Accessible via `protein.AtomIdentity(atom_idx)`.

### Primary identity tuple

The key that external NMR data binds to:

- `AminoAcid residue_type` (existing enum)
- `IupacAtomPosition` — new per-residue enum. Every valid (residue, IUPAC
  atom name) pair enumerated once. E.g., `IupacAtomPosition::ILE_CD1`,
  `IupacAtomPosition::ILE_HD11`, `IupacAtomPosition::ALA_HB3`. Zero
  runtime strings past load.

### Semantic lenses

Derived at load from topology + IUPAC convention. All typed enums; no
string dispatch:

- `NmrClass` — unified NMR classification absorbing the scattered
  booleans on ConformationAtom (`is_amide_H`, `is_alpha_H`, `is_methyl`,
  `is_aromatic_H`, ...). Values: `BackboneHN`, `BackboneHA`, `BackboneN`,
  `BackboneCA`, `BackboneCO_C`, `BackboneCO_O`, `SidechainAmideNH`,
  `GuanidiniumNH`, `AmmoniumNH`, `RingNH`, `HydroxylOH`, `CarboxylOH`,
  `ThiolSH`, `ThioetherS`, `DisulfideS`, `SidechainC_sp3`,
  `SidechainC_sp2_carbonyl`, `SidechainC_aromatic`,
  `SidechainC_sp2_guanidinium`, `MethylCH3`, `MethyleneCH2`, `MethineCH`,
  `RingCH`, and their non-H peers. ~30 values.
- `Locant` — sidechain depth letter: α, β, γ, δ, ε, ζ, η, `Backbone`,
  `NotApplicable`.
- `SidechainDepthBonds` — integer bond-count from nearest backbone atom.
  Redundant with `Locant`; numerically useful for distance-decay features.
- `MethylGroup` — enumerated methyl identity: `ALA_β`, `VAL_γ1`, `VAL_γ2`,
  `LEU_δ1`, `LEU_δ2`, `ILE_γ2`, `ILE_δ1`, `THR_γ2`, `MET_ε`, or `None`.
  pro-R / pro-S encoded by which enum slot the atom lands in.
- `methyl_partner_atom_index` — for twin methyls (VAL γ1↔γ2, LEU
  δ1↔δ2), partner atom's index. `SIZE_MAX` otherwise.
- `MethyleneGroup` — enum for methylenes, pro-R/pro-S encoded by slot.
- `methylene_partner_atom_index` — prochiral partner's atom index.
- `ChiParticipation` — bitfield, bit k set if atom is one of the four
  dihedral-defining atoms for `chi_k`.
- `RingAtomRole` — ring-member role. PHE/TYR: `Ipso`, `Ortho`, `Meta`,
  `Para`. TRP benzene: `C4/C5/C6/C7`. Indole perimeter + TrpPyrrole
  positions including `Nε1`. HIS variants: `Nδ1`, `Cε1`, `Nε2`, `Cδ2`,
  `Cγ`. `NotInRing` otherwise.
- `ring_membership` — bitfield of which `protein.Rings()` indices this
  atom belongs to (fused rings have atoms in multiple).

### Residue-level invariants, broadcast per atom

Stored on `NmrAtomIdentity` so per-atom slicing avoids a residue
indirection:

- `ResidueCategory`: `Hydrophobic`, `Polar`, `AcidicCharged`,
  `BasicCharged`, `Aromatic`, `GLY`, `PRO`, `CYS_disulfide_bonded`,
  `CYS_free`.
- `ResidueHydropathy` (Kyte-Doolittle), `ResidueVolume` (Chothia),
  `ResidueIsoelectricPoint` — scalars.
- `ResidueHbondDonorCount`, `ResidueHbondAcceptorCount` — ints.

### Symmetry / pseudo-atom

- `SymmetryClassId` — atoms that are chemically equivalent under
  topology and local symmetry (methyl Hs, PHE Hε1/Hε2, Hδ1/Hδ2) share
  an ID. For downstream pooling where topologically equivalent atoms
  should be averaged.
- `PseudoAtomName` — BMRB-compliant collapsed name (HB#, HG#). String
  at the edge only; not used in C++ logic.

### Orthogonal indices on Protein (pre-built, O(1) slice)

- `Protein::AtomsByNmrClass(NmrClass) -> const vector<size_t>&`
- `Protein::AtomsByLocant(Locant)`
- `Protein::AtomsByMethylGroup(MethylGroup)`
- `Protein::AtomsByRingAtomRole(RingAtomRole)`
- `Protein::AtomsByResidueCategory(ResidueCategory)`
- `Protein::AtomsByResidueAndPosition(AminoAcid, IupacAtomPosition)`
- `Protein::MethylGroups() -> vector<MethylGroupInfo>` (three atoms +
  twin-pair linkage)
- `Protein::MethyleneGroups()`

### Enum completeness is load-bearing; see Appendix A

The sketch above lists ~30 `NmrClass` values and suggests an
`IupacAtomPosition` enum with "example values" (e.g., `ILE_CD1`,
`ILE_HD11`, `ALA_HB3`). The actual enum inventory is substantially
larger and covers variants the sketch glosses over:

- HIS tautomers: `HID`, `HIE`, `HIP` (δ-protonated, ε-protonated,
  doubly-protonated — all three appear in the fleet depending on pKa /
  protonation assignment)
- CYS variants: `CYS` (free thiol), `CYX` (disulfide-bonded — SG atom
  role changes, no SG-H atom present)
- N-terminal caps: `H1`, `H2`, `H3` on backbone N when ammonium;
  ACE cap atoms (acetyl group on N-terminus)
- C-terminal caps: `OXT` (second carbonyl oxygen, when -COOH or -COO⁻);
  NME cap atoms (N-methyl amide on C-terminus)
- Protonation variants: `ASH` (protonated ASP), `GLH` (protonated GLU),
  `LYN` (neutral LYS). In practice the fleet uses ff14SB defaults with
  charged states; variants must still be representable.

**"Zero runtime strings past load" is a rule** — a single uncovered
atom drops through to a string fallback and the rule breaks silently.
Appendix A describes the generator-driven approach: a canonical
dictionary (TOML, one entry per `(residue, IUPAC position)` pair
including variants, with CHARMM / PDB / BMRB aliases) and a build-step
script that emits the C++ enum + mapping tables. Hand-writing the enum
would be error-prone and would drift from reality; generating it from
a reviewable source of truth gives us a place to fix bugs when they
surface.

### Relation to existing EnrichmentResult

`EnrichmentResult` currently writes categorical booleans
(`is_amide_H`, `is_methyl`, `is_aromatic_H`, ...) onto ConformationAtom.
That pattern was correct at the conformation level but awkward because
these booleans are topology-determined, not geometry-determined — they
don't belong on a per-conformation data store; they belong on Protein
via `NmrAtomIdentity`.

Clean path: EnrichmentResult continues to run during the per-frame
ConformationResult pipeline but writes no new fields itself (it becomes
a validator that the topology has been enriched correctly). `NmrAtomIdentity`
is populated once during Protein construction, with the same derivation
logic that EnrichmentResult currently runs. The ConformationAtom booleans
are kept (as redundant convenience) or deprecated in favour of
`NmrAtomIdentity` predicates — implementation call.

Consumers query `protein.AtomIdentity(i).nmr_class == NmrClass::BackboneHN`
instead of `conf.atoms[i].is_amide_H`. Same information, typed at the
right scope.

### Derivation algorithm: per-residue full-template matching

User preference, recorded 2026-04-22 pm: `NmrAtomIdentityBuilder`
derives identity via **full-residue pattern matching** against typed
templates, one residue at a time, with full match required. This is
the preferred algorithm, in explicit preference to any more efficient
or opportunistic scheme.

**Algorithm.**

1. For each residue in the Protein (ordered by sequence), look up the
   canonical atom template for this `(AminoAcid, variant)` pair in the
   generated residue-template database (see Appendix A — typed-enum
   generator sourced from PDB CCD or AMBER residue library, committed
   as repo data).
2. Match every expected atom in the template against the residue's
   actual atoms. The match criteria are joint:
   - Atom name in the residue matches the template's IUPAC name
     (via the NamingRegistry translation for CHARMM↔IUPAC gaps).
   - Bond connectivity to expected neighbours matches the template's
     edge list.
   - Stereochemistry slot (for twin methyls, methylenes, ring roles)
     matches the template's pro-R / pro-S / ring-position assignment.
3. **All expected atoms must match.** Missing atoms, extra atoms,
   connectivity mismatches — any of these flag the whole residue as
   non-canonical.
4. On full match: populate `NmrAtomIdentity` for every atom in this
   residue from the template slots. On mismatch: log loudly with the
   residue ID, the expected-vs-found difference, and mark the residue
   as `NonCanonical` with downstream code respecting that marker
   (typically: falls back to string-name access, does not produce
   typed-identity slicing results for that residue, but does not
   crash the extraction).
5. Move to next residue.

**Rationale.**

- **Robustness beats efficiency.** Atom-level opportunistic
  derivations can miss edge cases that only surface later as
  per-stratum statistical anomalies — you see a weird spike in
  calibration residuals at some residue subset and trace it back to
  a mis-labelled methyl partner. Full-residue matching fails at
  residue granularity immediately.
- **Typed top-down, not typed bottom-up.** The template is canonical;
  atoms inherit their typing from the template slot they fill. No
  atom is ever assigned a type "because it looks like one"; atoms
  get types because the whole residue matched a canonical template
  that says "this slot is pro-R γ₁ of VAL."
- **Debuggable by scope.** Mismatch logs are residue-local with
  specific expected-vs-found differences. Reproducing a label error
  means running the builder on one residue in isolation.
- **The slower path is fast enough.** 20–400 residues per protein,
  ~10–30 atoms per residue, template-match runs in milliseconds per
  Protein at construction time. Negligible vs OpenBabel's per-frame
  cost (retained per §"OpenBabel and EnrichmentResult") and vs the
  per-frame calculator pipeline.

**Anti-patterns this preference rules out:**

- No SMARTS pattern matching on individual atoms.
- No heuristic name-similarity fallbacks ("probably a methyl H
  because it's bonded to a C with three other Hs").
- No "assign what you can figure out; leave the rest as Other."
- No silent partial success where some atoms of a residue get
  typed identities and others get fallbacks, within the same residue.

The only acceptable outcomes per residue are: full match → fully
typed, or mismatch → `NonCanonical` marker + loud log.

### Validation discipline for NmrAtomIdentity

The derivation above is still "C++ pattern-matching code applied to
bond graphs" and can have bugs. External tests anchor the labels to
reference sources, caught via the test suite before fleet activation.

**External ground-truth sources, priority-ordered.**

1. **PDB Chemical Component Dictionary (CCD) — PRIMARY reference.**
   Covers every standard residue + every residue the PDB has ever
   seen, including the mutants (ALA, all 20 standard residues).
   Authoritative for atom names, bond connectivity, stereochemistry.
   The residue-template database (Appendix A) is generated from CCD
   in the first place; testing against CCD is testing the generation
   output against its source.
2. **AMBER ff14SB residue library — CROSS-CHECK.** Redundant with
   CCD for standard residues; useful as an independent check that the
   library's typed identity is consistent with the force field's
   residue definitions. Covers all muties.
3. **IUPAC rules (Markley 1998 and subsequent recommendations) —
   BASELINE.** Canonical source for pro-R / pro-S stereochemistry
   conventions that CCD / AMBER / BMRB all derive from. Used to
   verify stereochemistry cases (methyl / methylene twin pairing,
   ring role assignment) where CCD / AMBER encode the convention
   but didn't originate it.
4. **BMRB — SECONDARY, protein-specific.** Useful for the 10
   calibration proteins that have BMRB entries — the `nmr_forensics`
   work already demonstrated BMRB binding as a way to surface
   NamingRegistry gaps (`ChangesRequiredBeforeProductionH5Run`
   2026-04-20). **Not applicable to the mutant rerun** (no BMRB for
   artificial ALA mutants). Valuable for the downstream
   H5-binding-correctness check (does our typed identity resolve
   published shifts?) but not as a ground-truth source for the typed
   labels themselves.

The reordering reflects user clarification 2026-04-22 pm: *"This WILL
apply to the rerun of the muties, because we need these distinctions
for our stats, so bmrb is not probably our first best choice."* The
mutant pair rerun goes through this library; calibration stats depend
on typed per-atom stratification; PDB CCD (plus AMBER cross-check)
covers muties uniformly; BMRB doesn't. CCD is primary.

**Test categories.**

1. **Coverage (no atom is "Unknown").** Every atom in every
   calibration protein must produce a typed `IupacAtomPosition` from
   `NmrAtomIdentityBuilder`; every residue must match a canonical
   template (full match). No fallbacks to `Other`. No residues marked
   `NonCanonical` in the calibration set. Fails loudly on any atom
   that falls through or any residue that doesn't match its template.
   Catches enum-inventory gaps + template-database gaps.
2. **Geometry invariance.** Load the same Protein twice from different
   conformations (frame 0 vs frame 100 of an MD trajectory). Assert
   `protein_a.AtomIdentity(i) == protein_b.AtomIdentity(i)` for every
   atom. Catches accidental use of geometry in the derivation.
3. **Variant correctness.** Load a protein containing HIS residues
   three times, once each as HIE / HID / HIP. Assert the
   tautomer-specific atoms (`HE2`, `HD1`, both) resolve to the right
   `IupacAtomPosition`; `NmrClass::RingNH` attaches to the correct
   nitrogen in each variant. Same pattern for CYS / CYX, ASP / ASH,
   GLU / GLH, LYS / LYN. Catches variant-dispatch bugs.
4. **PDB CCD agreement (primary cross-check).** For each residue
   type in its canonical variant, load the PDB CCD entry via cifpp
   (or gemmi if we go that route). For each atom named in CCD, assert
   `NmrAtomIdentityBuilder`'s label agrees with CCD's atom attributes:
   name, element, bond connectivity, stereo descriptor, aromatic
   membership. This is the broadest cross-check; catches template-
   generation bugs that pass the coverage test but disagree with
   the authoritative source.
5. **Stereochemistry correctness via published benchmarks.** For
   ubiquitin (1UBQ — BMRB-backed + Loth-Pelupessy-Bodenhausen-L12-
   backed + a well-known reference protein), assert MethylGroup slot
   assignments for VAL / LEU / ILE / THR / MET match the IUPAC
   convention. Reference assertion list committed to
   `tests/nmr_identity/ubiquitin_methyl_stereo.txt`. Same for
   methylene stereochemistry (HB2/HB3 etc.) at β-methylene residues.
6. **BMRB cross-check (secondary, 10-protein only).** For each of
   the 10 calibration proteins, load its BMRB entry. For every shift
   row, resolve to a library atom via NmrAtomIdentity primary tuple.
   Assert 100% resolution rate. Every remaining ambiguity is either
   expected (BMRB's own ambiguity flag set) or a bug. Cheaper-to-run
   than CCD cross-check; catches CHARMM↔IUPAC translation gaps.

Tests 1, 2, 6 run on every CI pass. Test 3 requires the variant
handling to be complete, so runs on the pre-activation validation
pass. Test 4 requires cifpp ChemComp access (or gemmi); runs on
pre-activation. Test 5 requires the curated stereochemistry
assertion list, committed once; runs on every CI pass.

**Activation gate.** The rollup does not activate for the fleet
until all six test categories pass on all 10 calibration proteins,
with zero fallbacks to `NonCanonical` and zero unbound BMRB rows
beyond the expected ambiguity flags. A failure in any test category
is a label bug; label bugs do not ship to the fleet.

### OpenBabel and EnrichmentResult: no migration (decided 2026-04-22)

**Decision:** leave OpenBabel and `EnrichmentResult` unchanged. User
preference (verbatim): *"I am happy to pay for openbabel per frame in
exchange for not mucking with everything on that level."*

`EnrichmentResult` continues to run per-frame in the OperationRunner
sequence, invoking OpenBabel for hybridisation and populating
`ConformationAtom.hybridisation` plus the categorical booleans
(`is_amide_H`, `is_methyl`, `is_aromatic_H`, etc.) as it does today.
No migration to Protein construction; no shim phase; no cross-cutting
calculator changes.

**What this means for NmrAtomIdentity.** It becomes a purely additive
layer on `Protein` alongside the existing per-conformation enrichment
machinery. New code paths (particularly trajectory-scope TrajectoryResults
that want typed atom slicing — `Protein::AtomsByNmrClass(BackboneHN)`)
consume NmrAtomIdentity. Existing calculators that read
`conf.atoms[i].is_amide_H` keep reading the boolean as they do today.
No deprecation of the ConformationAtom booleans; they stay as per-
conformation convenience fields even though they're redundant with
NmrAtomIdentity queries — paying for redundancy is the trade for not
mucking with existing pipelines.

**NmrAtomIdentity's derivation is independent of EnrichmentResult.**
`NmrAtomIdentityBuilder` runs once at `Protein::FinalizeConstruction()`,
reads the same topology inputs (bond graph, residue list, ring
detection output), produces `NmrAtomIdentity` per atom. It doesn't
call OpenBabel; hybridisation is not part of NmrAtomIdentity
(hybridisation stays on ConformationAtom via EnrichmentResult as
today). The builder uses only: (residue type, IUPAC atom name, bond
connectivity from OpenBabel-already-run-by-EnrichmentResult on frame
0, ring type assignments from `Protein::DetectAromaticRings`). Since
these are all topology-fixed, the builder runs once per Protein,
never again.

Appendix E reflects this decision (it had proposed the migration path;
now explicitly declined).

---

## 3 — `TrajectoryProtein`: the running buffer at trajectory scope

`TrajectoryProtein` is the model of a protein in a trajectory. It wraps a
`Protein` (identity + topology) and adds the trajectory-scope running
buffer: a vector of `TrajectoryAtom` objects (per-atom, parallel to
`protein.atoms` in indexing) plus attached `TrajectoryResult` objects.

It replaces the existing `GromacsProtein` class — same role in the system,
corrected in content.

---

## 🛑 READ THIS FIRST — ANTI-PATTERNS THE CURRENT CODE EMBODIES THAT MUST NOT SURVIVE THE REFACTOR

This refactor is not a rename. It is the replacement of an
architecturally wrong pattern with a different pattern. A future
session reading the existing code under time pressure, or trained on
conventional C++ idioms, will gravitate toward preserving the existing
structure under new names. **That is not the refactor.** The points
below exist specifically to prevent that outcome.

### Anti-pattern 1: accumulator state baked into the per-atom struct

**What the current code does** (see `src/GromacsProteinAtom.h` for the
live version):

```cpp
// WRONG — THIS PATTERN DOES NOT SURVIVE THE REFACTOR
struct GromacsProteinAtom {
    Welford bs_T0;           // Welford = mean, M2, min, max, frame indices
    Welford bs_T2mag;
    Welford hm_T0;
    DeltaTracker bs_T0_delta;
    TransitionCounter chi1_transitions;
    Welford water_emag;
    // ... 50+ accumulator-state instances as fields
};
```

**Why this is wrong:**

1. The `Welford` / `DeltaTracker` / `TransitionCounter` objects carry
   *implementation state of the accumulation process* — a running
   mean, a sum-of-squared-deviations, a ring buffer, a last-frame
   value, a threshold-crossing count. That is the guts of how
   accumulation computes, not the outputs of accumulation. Baking it
   into the per-atom data struct makes the struct a god-struct
   holding every physics subsystem's internal accumulator state
   alongside every other's.

2. Adding a new accumulator requires modifying the struct
   *and* updating a parallel enumeration method (`AllWelfords`)
   *and* updating the per-frame dispatcher (`AccumulateFrame`).
   Three edits in three files for one new accumulator kind. The
   struct doesn't own its fields; the world mutates them.

3. No singleton-per-type discipline. Any code with access to the
   struct can read or write any field. There is nothing that says
   "the `bs_T0` Welford is owned by the BS accumulator alone."

4. Heterogeneous accumulators can't coexist without the struct
   growing additional typed fields. A windowed autocorrelation
   accumulator has a rolling buffer of 120 per-atom tensors; an
   FFT accumulator has frequency-domain state; a transition
   counter has a last-bin integer. The struct can't express this
   variability without exploding in size and losing type coherence.

5. Output layout depends on hand-coded enumeration order in
   `AllWelfords`. Reorder the enumeration by accident and downstream
   schema shifts silently.

**What must replace it:**

`TrajectoryAtom` is a **declared-upfront bag of finalized output
fields** — typed summary values only. `double bs_t0_mean`, `double
bs_t0_std`, `int chi1_transitions_count`. These are OUTPUTS, not
accumulator state. Accumulator state lives *inside `TrajectoryResult`
subclasses* and never leaks out onto the per-atom struct.

Each field on `TrajectoryAtom` has **exactly one writer** (the
`TrajectoryResult` subclass whose `Compute` or `Finalize` writes it).
The attach discipline enforces singleton-per-type.

There is no `AllWelfords()` or equivalent. Serialization is by each
`TrajectoryResult`'s `WriteH5Group` — no central enumeration, no
parallel list to maintain.

**The refactor temptation to avoid:**

> *"I'll just make `Welford` / `DeltaTracker` / `TransitionCounter`
> inherit from an abstract `Accumulator` base class with a virtual
> `Update()` method. Then the struct can hold
> `std::vector<std::unique_ptr<Accumulator>>`. Much cleaner C++."*

**This is the wrong direction.** It generalises the anti-pattern
instead of replacing it. The accumulator state would still live on
the per-atom struct; the polymorphism would just hide the god-struct
problem behind virtual dispatch. The correct direction is:
accumulator state moves **off the per-atom struct entirely** and onto
typed `TrajectoryResult` objects that own their own state.

### Anti-pattern 2: monolithic per-frame dispatcher (`AccumulateFrame`)

**What the current code does** (see
`src/GromacsProtein.cpp::AccumulateFrame`):

```cpp
// WRONG — THIS PATTERN DOES NOT SURVIVE THE REFACTOR
void GromacsProtein::AccumulateFrame(
    const ProteinConformation& conf, size_t frame_idx, double time)
{
    for (size_t i = 0; i < atom_count; ++i) {
        auto& ga = atoms_[i];  // GromacsProteinAtom
        const auto& ca = conf.AtomAt(i);
        // Monolithic update: every field, every frame, hand-written
        ga.bs_T0.Update(ca.bs_shielding_contribution.T0, frame_idx);
        ga.bs_T2mag.Update(T2Magnitude(ca.bs_shielding_contribution), frame_idx);
        ga.hm_T0.Update(ca.hm_shielding_contribution.T0, frame_idx);
        // ... 50+ more lines per atom, hand-edited when fields added
    }
}
```

**Why this is wrong:**

- Monolithic per-frame work. A single method knows every accumulator
  that exists in the library.
- `GromacsProtein` (the trajectory-context adapter) has to import
  every accumulator's knowledge about which ConformationAtom field it
  reads. Cross-subsystem coupling at the dispatcher layer.
- Adding a new accumulator requires editing this method, adding a
  line per accumulator. The new accumulator doesn't own its own
  per-frame compute; the central dispatcher owns it.
- Per-atom inner loop hand-written. No ordering discipline. No
  dependency declaration. No typed dispatch.

**What must replace it:**

The frame handler iterates attached `TrajectoryResult`s and calls
each one's virtual `Compute(conf, tp, frame_idx, time_ps)`. Each
`TrajectoryResult` reads the `ConformationAtom` fields it cares about
and writes to its own `TrajectoryAtom` fields. Polymorphic dispatch
through the base class; no central dispatcher-of-everything.

**The refactor temptation to avoid:**

> *"I'll just keep `AccumulateFrame` but have it iterate a list of
> accumulator pointers stored on GromacsProtein. Each accumulator has
> a method I call in the loop. Same shape, a bit cleaner."*

**This is half the refactor.** It adds polymorphism over the
accumulators but keeps the central dispatcher holding the list, keeps
the per-atom struct as the output location, keeps the god-struct
pattern. The correct direction: the frame handler (`GromacsFrameHandler`)
iterates `tp.ResultsInAttachOrder()` — which is a method on
`TrajectoryProtein`, which *is the adapter formerly known as
GromacsProtein*, which holds attached `TrajectoryResult`s through the
same typed-result mechanism `ProteinConformation` uses for
`ConformationResult`. No `AccumulateFrame` method survives.

### Anti-pattern 3: `AllWelfords()` / hand-enumerated field lists for I/O

**What the current code does** (see
`src/GromacsProteinAtom.h::AllWelfords`):

```cpp
// WRONG — THIS PATTERN DOES NOT SURVIVE THE REFACTOR
std::vector<std::pair<std::string, const Welford*>> AllWelfords() const {
    return {
        {"bs_T0",          &bs_T0},
        {"bs_T2mag",       &bs_T2mag},
        {"hm_T0",          &hm_T0},
        // ... hand-maintained parallel list
    };
}
```

**Why this is wrong:**

- Central enumeration couples output order and schema to hand-
  maintained list contents.
- Can only enumerate `Welford` instances (note the return type).
  `DeltaTracker` and `TransitionCounter` need their own parallel
  enumeration methods. Every new accumulator type = new enumeration
  method + new output logic elsewhere that consumes it.
- Adding a field and forgetting to add to `AllWelfords()` silently
  drops it from H5 output. No compile-time check.

**What must replace it:**

Each `TrajectoryResult` has `WriteH5Group(tp, file)` that knows its
own fields and emits them. The top-level writer iterates
`tp.ResultsInAttachOrder()` and calls `WriteH5Group` on each. No
central enumeration. No hand-maintained parallel list. Adding a new
`TrajectoryResult` = new subclass with `WriteH5Group` override; that
class's output lands in the H5 without touching any central code.

### Anti-pattern 4: parallel `Gromacs*` typed hierarchy as a substitute for canonical types

**What the current code does:** `GromacsProtein` holds `atoms_`
(vector of `GromacsProteinAtom`) as a parallel per-atom store
alongside the canonical `Protein`'s topology. `GromacsRunContext`
holds trajectory-level state parallel to `ProteinBuildContext`.
`GromacsProteinAtom` is a parallel data store to `ConformationAtom`.

**Why this is wrong:**

- Duplicates data stores without unifying them. Every atom has a
  canonical `Atom` (identity, on Protein), a `ConformationAtom`
  (per-conformation computed, on ProteinConformation), *and* a
  `GromacsProteinAtom` (trajectory-scope accumulated, on
  GromacsProtein). Three stores per atom. The third's existence is
  a workaround, not a modelling commitment.
- Turns "trajectory vs static" into a type-dispatch concern. Code
  that works with a static conformation uses `Protein +
  ProteinConformation`. Code that works with a trajectory uses
  `GromacsProtein + parallel machinery`. Same physics, two code
  paths.
- "Gromacs" in the name is lying about the scope. The concepts
  (trajectory-scope running buffer, per-atom accumulation) are not
  Gromacs-specific; they're trajectory-scope generic. The prefix
  pretends the machinery is format-specific when it isn't.

**What must replace it:**

`TrajectoryProtein` replaces `GromacsProtein` *by name and by
content*. It wraps the canonical `Protein` (as `GromacsProtein` does
today) but its per-atom store is `TrajectoryAtom`, which is the
trajectory-scope peer of `ConformationAtom` at the same logical level
of the object model — not a third store invented to work around a
naming mismatch.

Per-atom identity still flows through the canonical `Protein` via
back-pointer (no triplicated identity fields). Trajectory-scope
computed data lives on `TrajectoryAtom`, same typing discipline as
`ConformationAtom`.

**The refactor temptation to avoid:**

> *"I'll rename `GromacsProtein` to `TrajectoryProtein`,
> `GromacsProteinAtom` to `TrajectoryAtom`, and I'm done. The
> existing code works; the names just needed updating."*

**This is the rename-but-not-replace failure mode**. It preserves the
anti-patterns (1, 2, 3) under new names. A file named
`TrajectoryAtom.h` that contains 50 `Welford` instances as fields is
the thing the refactor specifically exists to prevent. If the
implementing session produces that file, the refactor hasn't
happened — only a rename has.

---

### What the successful refactor looks like (test for "am I doing it right?")

After the refactor:

- `grep -r "Welford" src/` returns hits only inside specific
  `TrajectoryResult` subclass sources (e.g.,
  `BsWelfordTrajectoryResult.cpp`) where the accumulation
  implementation lives. **No `Welford` fields on any per-atom
  struct.**
- `grep -r "AllWelfords\|AccumulateFrame" src/` returns zero hits.
  Those methods do not survive.
- `src/TrajectoryAtom.h` declares typed summary fields only
  (doubles, ints, SphericalTensors, small arrays). No
  accumulator-state types as fields.
- `src/TrajectoryProtein.h` holds `std::vector<TrajectoryAtom>
  atoms_`, a result map
  `std::unordered_map<std::type_index,
  std::unique_ptr<TrajectoryResult>>`, wrapped `Protein`, and
  nothing else structural.
- `GromacsFrameHandler::Next` iterates
  `tp.ResultsInAttachOrder()` and calls `result->Compute(...)` on
  each. No other per-frame accumulator-level calls.
- `GromacsProtein` / `GromacsProteinAtom` / `GromacsRunContext` do
  not exist in `src/`. They moved to `learn/bones/` with a short
  note explaining why.

If any of those doesn't hold, the refactor is incomplete — regardless
of what compiles and what tests pass. Rename-without-replace compiles
and passes existing tests. It's still wrong.

### Recommended file-header comment for the new classes

To help the implementing session keep the pattern straight while
writing, the following comment blocks go at the top of the new files.
They embed the architectural commitment in the code itself.

**`src/TrajectoryAtom.h`** (top-of-file, before any code):

```cpp
// TrajectoryAtom: per-atom trajectory-scope DATA STORE.
//
// ⚠ ARCHITECTURAL NOTE — READ BEFORE MODIFYING ⚠
//
// This class replaces GromacsProteinAtom. The old struct held
// Welford / DeltaTracker / TransitionCounter instances as fields
// (accumulator IMPLEMENTATION STATE baked into a per-atom data
// struct). That pattern was wrong and must not be reintroduced under
// new names.
//
// TrajectoryAtom holds finalized OUTPUT fields only — typed summary
// values (double, int, SphericalTensor, small arrays).
//
// DO NOT add as fields on this class:
//   - Welford, DeltaTracker, TransitionCounter, or any other
//     accumulator-state object.
//   - std::vector<std::unique_ptr<Accumulator>> or similar
//     "polymorphic accumulator list" — that is the anti-pattern
//     under virtual dispatch.
//   - Any field whose writer is not clearly identifiable as a
//     specific TrajectoryResult subclass (singleton-per-type).
//
// DO add as fields on this class:
//   - Finalized means / stds / deltas / counts (doubles, ints).
//   - Finalized tensor summaries (SphericalTensor, Mat3, Vec3).
//   - Rich per-source structured vectors (paralleling
//     ConformationAtom::ring_neighbours) populated at Finalize.
//
// Each field has EXACTLY ONE WRITER: the TrajectoryResult subclass
// that owns the physics producing it. Enforced by attach discipline
// at TrajectoryProtein::AttachResult.
//
// For the full rationale, see
// spec/WIP_OBJECT_MODEL.md §3 anti-patterns subsection.
```

**`src/TrajectoryResult.h`** (top-of-file):

```cpp
// TrajectoryResult: base class for per-trajectory modular calculators.
// Parallel to ConformationResult at conformation scope.
//
// ⚠ ARCHITECTURAL NOTE — READ BEFORE MODIFYING ⚠
//
// TrajectoryResult subclasses OWN their per-frame accumulator state.
// Welford, DeltaTracker, TransitionCounter, rolling windows, FFT
// buffers — all internal state of the subclass, populated during
// Compute, dissolved or transferred at Finalize.
//
// They do NOT write accumulator state onto TrajectoryAtom. They
// write finalized OUTPUT fields onto TrajectoryAtom, one writer per
// field, enforced by singleton-per-type discipline at attach time.
//
// If you find yourself writing
//   struct MyTrajectoryResult { std::vector<Welford> per_atom_; };
// that is CORRECT (accumulator state owned by the result).
//
// If you find yourself writing
//   struct TrajectoryAtom { Welford bs_T0; /* ... */ };
// that is WRONG (accumulator state on the per-atom struct). Stop;
// re-read §3 of the WIP object model doc.
//
// For the full pattern, see spec/WIP_OBJECT_MODEL.md §4 and the two
// worked examples (BsWelfordTrajectoryResult,
// BsT0AutocorrelationTrajectoryResult).
```

These comments are redundant with the spec doc on purpose. The spec
is not always in the implementing session's context; the comments
always are.

---

### The role

During streaming: the running buffer that TrajectoryResults write to as
they compute per-frame updates. TrajectoryResult::Compute(conf, tp, idx,
time) reads from the frame's ConformationAtoms, reads and writes this
TrajectoryProtein's TrajectoryAtoms, reads and writes its own internal
state.

After streaming: the finalised record of what the trajectory produced.
Queryable. Holds the finalised per-atom fields, the attached
TrajectoryResults with their query methods, and the wrapped Protein for
identity.

### Structure

```cpp
class TrajectoryProtein {
public:
    explicit TrajectoryProtein(std::unique_ptr<Protein> protein);

    // Identity delegation — TrajectoryProtein does not duplicate
    // topology. Ask the wrapped Protein.
    const Protein& Protein() const;
    size_t AtomCount() const;  // == protein().AtomCount()

    // Per-atom trajectory data store.
    const TrajectoryAtom& AtomAt(size_t atom_idx) const;
    TrajectoryAtom& MutableAtomAt(size_t atom_idx);
    const std::vector<TrajectoryAtom>& Atoms() const;

    // Attached TrajectoryResults. Typed singleton-per-class.
    bool AttachResult(std::unique_ptr<TrajectoryResult> result);

    template <typename T> T& Result() const;        // throws if missing
    template <typename T> bool HasResult() const;

    const std::unordered_map<std::type_index,
          std::unique_ptr<TrajectoryResult>>& AllResults() const;

    // Iteration order matches attach order (so per-frame dispatch
    // respects the order the caller attached, which respects
    // declared dependencies by the singleton + deps check at attach).
    const std::vector<TrajectoryResult*>& ResultsInAttachOrder() const;

    // Dense buffers transferred from TrajectoryResults at their Finalize
    // (e.g., per-atom-per-frame time series stored contiguously for
    // cache-friendly access). Typed by the donating TrajectoryResult.
    // See §4 "Dense buffers" for the pattern.
    template <typename T>
    void AdoptDenseBuffer(std::unique_ptr<DenseBuffer<T>> buffer,
                          std::type_index owner);

    // State for gated reads.
    bool IsFinalized() const { return finalized_; }

private:
    std::unique_ptr<class Protein> protein_;
    std::vector<TrajectoryAtom> atoms_;
    std::unordered_map<std::type_index,
          std::unique_ptr<TrajectoryResult>> results_;
    std::vector<TrajectoryResult*> results_attach_order_;
    std::unordered_map<std::type_index,
          std::unique_ptr<DenseBufferBase>> dense_buffers_;
    bool finalized_ = false;
    // ... sources-linked metadata as §5 Trajectory determines.
};
```

### What goes on `TrajectoryProtein` vs `TrajectoryAtom` vs TrajectoryResult internal state

The partition is the same as at conformation scope, with one trajectory-
specific addition:

- **On `TrajectoryProtein` (protein-level, not per-atom):** identity
  references (the wrapped Protein), attached TrajectoryResults, dense
  buffers whose indexing is (atom, frame) or similar and who are too
  large to fit as per-atom `std::vector` fields without wasteful
  per-atom allocation.
- **On `TrajectoryAtom` (per-atom, declared upfront, one writer per
  field):** running scalar fields (always-valid mid-stream where the
  accumulator shape allows it, finalised-only where it doesn't),
  small-array fields (e.g., per-atom 8-type Welford arrays), rich
  structured data as `std::vector<StructuredPerSource>` parallel to
  `ConformationAtom::ring_neighbours`.
- **Internal state on `TrajectoryResult`:** rolling windows, FFT
  buffers, intermediate matrices for block-averaging diagnostics — anything
  whose shape is not per-atom or whose lifetime is only during streaming.

The partition question per field is: does this value have a single atom
as its primary index, and is it small enough to not require central
buffer pooling? If yes → `TrajectoryAtom` field. If it's inherently
multi-indexed (atom × frame, atom × lag, atom × frequency, ring × frame,
pair × pair × frame), consider a dense buffer on TrajectoryProtein with
typed access methods. If it's not per-atom at all (e.g., FFT working
state), it's internal to the TrajectoryResult.

### Memory math for the declared-upfront discipline

The "all fields present even if unused" discipline on ConformationAtom
allocates one field slot per atom per field; for a 4000-atom protein
with ~216 fields (current ConformationAtom), that's ~864K field slots
per conformation, roughly 150–300 MB per materialised frame (position
Vec3s + tensor Mat3s + SphericalTensor components + ring_neighbours
vectors + bond_neighbours vectors). At trajectory scope the discipline
applies differently because the scope of the store is different:

- **One TrajectoryProtein per protein per run**, not one per frame.
  A fully-populated TrajectoryAtom with the proposed field set —
  summary Welfords + DeltaTrackers + TransitionCounters + rich
  per-source vectors like `ring_neighbour_stats` — at 4000 atoms
  totals ~80–200 MB, allocated once per run.
- **Dense buffers absorb the per-atom × frame / × lag multiplicands
  separately** (see "Dense buffers on TrajectoryProtein" in §4).
  A `PositionsTimeSeriesTrajectoryResult` holding (4000 atoms × 1200
  frames × Vec3) is ~115 MB; a `BsShieldingTimeSeriesTrajectoryResult`
  with SphericalTensor per atom per frame is ~345 MB. Handful of such
  buffers fits comfortably.
- **Running total during a run:** one materialised ProteinConformation
  (~200 MB) + one TrajectoryProtein with all TrajectoryAtoms (~100 MB)
  + sum of attached dense buffers (~2–4 GB under `PerFrameExtractionSet`
  for a 25 ns protein). Against 128 GB available per trajectory worker,
  this is comfortable headroom.
- **μs harvester phase** (20–40× frame count): dense buffers scale
  linearly → 40–160 GB per run. Still fits, but at the upper end; worth
  a perf-test run at μs length before fleet activation of that phase.
  Selective TrajectoryResult attachment (runs that don't need every
  time series) is the lever.

No memory concern at current fleet scale. Flag for re-check at μs
harvester activation.

### Per-bond trajectory state: option (b) decided 2026-04-22

**STATUS: DECIDED — option (b).** User picked (b) 2026-04-22: per-bond
state lives internal to bond-scope TrajectoryResults (e.g.,
`BondLengthStatsTrajectoryResult` holds
`std::vector<PerBondWelford> per_bond_` indexed by
`protein.BondAt(i).index()`; queries route through the Result's
methods; H5 emission by the Result itself). No first-class
`TrajectoryBond` store on TrajectoryProtein.

Rationale (short): matches §3's anti-pattern principle that
accumulator state lives on the Result, not on the physical-world
data store. Bond counts are small and bond-scope calculators are few
(2–3 total: length stats, length-delta stats, Wiberg stats in FullFat),
so a parallel typed entity for symmetry-with-TrajectoryAtom doesn't
earn its weight. Trade-off: loss of symmetry with TrajectoryAtom
accepted — the honest asymmetry is that per-atom running state is
shared across multiple Results, per-bond state is not.

Appendix B below preserves the full option-(a) design as the
rejected alternative. It should not be implemented; it is context for
why option (b) was chosen.

---

`GromacsProtein::bonds_accum_` today carries per-bond Welfords (bond
length mean / variance / min / max; per-bond length-fluctuation rate).
Those observables remain wanted; the per-field-per-bond Welfords now
live inside `BondLengthStatsTrajectoryResult` internal state, emitted
at Finalize into the Result's own H5 group (e.g., `/bonds/length/*`).

### Identity-through-back-pointer

`TrajectoryAtom` holds no identity data. Element, bonds, residue,
NmrAtomIdentity all flow through the wrapped Protein:

```cpp
const TrajectoryAtom& ta = tp.AtomAt(42);
const Atom& identity_atom = tp.Protein().AtomAt(42);
Element elem = identity_atom.element;
NmrClass cls = tp.Protein().AtomIdentity(42).nmr_class;

// Trajectory-scope data
double mean_bs_t0 = ta.bs_t0_mean;  // always-valid
SphericalTensor autocorr_result =
    tp.Result<BsAutocorrelationTrajectoryResult>().LagAtAtom(42, 10);
```

One index, three tiers of information (identity, trajectory-scope
summary, trajectory-scope rich query). No duplication. Parallel to how
ConformationAtom's identity flows through `conf.Protein().AtomAt(i)`.

### `TrajectoryAtom`: per-atom running buffer

`TrajectoryAtom` is the trajectory-scope analog of `ConformationAtom`.
Private constructor, owned exclusively by TrajectoryProtein's
`atoms_` vector, constructed once, never resized. Public typed fields
declared upfront. Singleton-per-class discipline on writers enforced by
which TrajectoryResult type owns which field.

No `Position()` member — per-atom position is conformation-scoped, not
trajectory-scoped. Frame-wise positions during streaming are on the
current `ProteinConformation`; if the trajectory wants to preserve a
time series of positions, that lives on a `FramePositionsTrajectoryResult`
(via a dense buffer attached to TrajectoryProtein — see §4).

The field set is declared upfront per the "ALL fields present even if
unused" discipline on ConformationAtom. Current inventory is defined by
the accumulator-kin fields in `GromacsProteinAtom` today, redistributed
as the TrajectoryResult subclasses land. Running means / standard
deviations / deltas / transition counts are the natural shape; per-lag
functions are not (those go to dense buffers).

Concrete example of what `TrajectoryAtom` looks like (narrative, not
literal code — actual field list emerges from TrajectoryResult declarations):

```cpp
class TrajectoryAtom {
    friend class TrajectoryProtein;
public:
    // Written by BsWelfordTrajectoryResult. Always-valid mid-stream:
    // Compute() updates these fields directly each frame.
    double bs_t0_mean = 0.0;
    double bs_t0_m2 = 0.0;        // sum of squared deviations; std derived
    double bs_t0_min = std::numeric_limits<double>::infinity();
    double bs_t0_max = -std::numeric_limits<double>::infinity();
    size_t bs_n_frames = 0;

    // Written by BsDeltaTrajectoryResult. Also updated each frame;
    // Finalize divides by (n-1) for variance of deltas.
    double bs_t0_delta_mean = 0.0;
    double bs_t0_delta_m2 = 0.0;

    // Written by McWelfordTrajectoryResult. Same pattern.
    double mc_t0_mean = 0.0;
    double mc_t0_m2 = 0.0;
    // ...

    // Written by RotamerTransitionTrajectoryResult. Updated when a
    // threshold crossing is detected.
    int chi1_transitions = 0;
    int chi2_transitions = 0;

    // Written by HBondLifetimeTrajectoryResult. Finalize-only fields.
    // Set only at Finalize after the full stream has been observed.
    double hbond_mean_lifetime_ps = 0.0;
    int hbond_break_form_events = 0;

    // Rich structured per-source data (parallel to
    // ConformationAtom::ring_neighbours) written by
    // RingNeighbourhoodTrajectoryResult. Per atom, per near-by ring,
    // the trajectory-averaged neighbourhood descriptor.
    std::vector<RingNeighbourhoodTrajectoryStats> ring_neighbour_stats;

private:
    explicit TrajectoryAtom() = default;
};
```

The declared-upfront pattern means adding a new TrajectoryResult that
produces per-atom fields requires adding fields here. That's the same
pattern as adding fields on ConformationAtom when adding a new
ConformationResult — one place to declare, one writer per field,
typed readers.

---

## 4 — `TrajectoryResult`: the modular calculator at trajectory scope

`TrajectoryResult` is the trajectory-scope analog of `ConformationResult`.
The modularity properties of ConformationResult (see `PATTERNS.md` and
`src/ConformationResult.h`) transfer directly; the structural difference
is multi-phase compute across frames rather than one-shot compute on a
conformation.

### Base class

```cpp
class TrajectoryResult {
public:
    virtual ~TrajectoryResult() = default;

    // Human-readable name, for logging and AttachResult diagnostics.
    virtual std::string Name() const = 0;

    // Declared dependencies. Type indices of:
    //   - other TrajectoryResult types that must be attached first, AND/OR
    //   - ConformationResult types that must be run per frame (i.e., the
    //     RunConfiguration's per-frame calculator set must include them).
    // Attach-time check fires if any declared dependency is missing.
    virtual std::vector<std::type_index> Dependencies() const = 0;

    // THE work method. Called once per frame during streaming.
    //
    // Arguments:
    //   conf       — this frame's ProteinConformation, with its
    //                ConformationResults already attached and its
    //                ConformationAtoms populated.
    //   tp         — the TrajectoryProtein running buffer. Reads from
    //                and writes to tp.AtomAt(i) fields this result
    //                owns. Reads from tp.AtomAt(i) fields owned by
    //                earlier-attached TrajectoryResults.
    //   frame_idx  — zero-based frame index in the current traversal.
    //   time_ps    — simulation time of this frame.
    //
    // No return value. Errors are logged via OperationLog, not
    // thrown; a failing Compute either skips its update for this
    // frame or marks its own output invalid via a state flag.
    virtual void Compute(const ProteinConformation& conf,
                         TrajectoryProtein& tp,
                         size_t frame_idx,
                         double time_ps) = 0;

    // End-of-stream synthesis. Called once after the last frame.
    // Default no-op for results whose TrajectoryAtom fields are
    // always-valid mid-stream (pure Welfords, transition counters).
    // Override for results that need end-of-stream work:
    //   - FFT over a buffered window
    //   - division by (n-1) for unbiased variance
    //   - transferring dense-buffer ownership to TrajectoryProtein
    //   - computing derived fields from accumulated state
    virtual void Finalize(TrajectoryProtein& tp) {}

    // Self-serialisation. Same discipline as ConformationResult —
    // each result knows its own fields and writes them. The top-level
    // writer iterates attached results and calls this on each.
    // Default returns 0 (no arrays written).
    virtual int WriteFeatures(const TrajectoryProtein& tp,
                              const std::string& output_dir) const { return 0; }

    // H5 group emission — parallel to WriteFeatures but targeting the
    // analysis H5 schema. Each result writes its own group.
    virtual void WriteH5Group(const TrajectoryProtein& tp,
                              HighFive::File& file) const {}
};
```

Four virtuals on the base instead of ConformationResult's three. The
extra one (`Compute`) is virtual-on-base because per-frame iteration
dispatches through the base class — the frame loop doesn't know
concrete types. ConformationResult's Compute is a static factory on
the subclass because conformation-scope compute is one-shot and attaches
a fully-populated result. For multi-phase streaming compute, the result
is constructed empty at attach and called polymorphically per frame.

### Concrete example: `BsWelfordTrajectoryResult`

A full worked example, parallel to `BiotSavartResult` at conformation
scope.

```cpp
// BsWelfordTrajectoryResult: running mean, variance, min, max of the
// BiotSavart T0 shielding contribution per atom, accumulated across
// all frames of a trajectory.
//
// Source: each frame, after BiotSavartResult::Compute has run on the
// frame's ProteinConformation, this result reads
// conf.atoms[i].bs_shielding_contribution.T0 and ticks a Welford
// update per atom. Running mean / M2 / min / max live as fields on
// TrajectoryAtom; they are always valid mid-stream.
//
// Finalize divides M2 by (n-1) and writes the final std to
// TrajectoryAtom.bs_t0_std.

class BsWelfordTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override { return "BsWelfordTrajectoryResult"; }

    std::vector<std::type_index> Dependencies() const override {
        // Per-frame dependency: BiotSavartResult must run each frame.
        return { std::type_index(typeid(BiotSavartResult)) };
    }

    // Factory: create empty with running Welford state zeroed.
    static std::unique_ptr<BsWelfordTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp) override;

    int WriteFeatures(const TrajectoryProtein& tp,
                      const std::string& output_dir) const override;

    // Post-finalize query methods. Callers can also read the fields
    // directly from TrajectoryAtom; these exist for convenience.
    double MeanAtAtom(size_t atom_idx) const;
    double StdAtAtom(size_t atom_idx) const;  // valid only after Finalize
    size_t NumFrames() const { return n_frames_; }

private:
    size_t n_frames_ = 0;
    bool finalized_ = false;
};
```

`Compute` implementation, concrete:

```cpp
void BsWelfordTrajectoryResult::Compute(
    const ProteinConformation& conf,
    TrajectoryProtein& tp,
    size_t frame_idx,
    double time_ps)
{
    // Read from frame's ConformationAtoms; write to tp's TrajectoryAtoms.
    // Welford update in place, field-by-field.
    const size_t N = tp.AtomCount();
    for (size_t i = 0; i < N; ++i) {
        const double x = conf.AtomAt(i).bs_shielding_contribution.T0;
        TrajectoryAtom& ta = tp.MutableAtomAt(i);

        // Running min/max.
        if (x < ta.bs_t0_min) ta.bs_t0_min = x;
        if (x > ta.bs_t0_max) ta.bs_t0_max = x;

        // Welford update. Fields are read, updated, written back.
        // After this frame, ta.bs_t0_mean is the mean over frames [0..idx].
        const size_t n = ta.bs_n_frames + 1;
        const double delta = x - ta.bs_t0_mean;
        const double new_mean = ta.bs_t0_mean + delta / n;
        const double new_m2 = ta.bs_t0_m2 + delta * (x - new_mean);
        ta.bs_t0_mean = new_mean;
        ta.bs_t0_m2 = new_m2;
        ta.bs_n_frames = n;
    }
    ++n_frames_;
}
```

`Finalize` implementation:

```cpp
void BsWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp) {
    const size_t N = tp.AtomCount();
    for (size_t i = 0; i < N; ++i) {
        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        // Unbiased std from M2.
        if (ta.bs_n_frames > 1) {
            ta.bs_t0_std = std::sqrt(ta.bs_t0_m2 / (ta.bs_n_frames - 1));
        } else {
            ta.bs_t0_std = 0.0;
        }
    }
    finalized_ = true;
}
```

Points from this example:

- The Compute method writes directly to TrajectoryAtom fields.
  `bs_t0_mean` is updated in place each frame and is always the mean
  over frames seen so far. Mid-stream reads are valid.
- `bs_t0_std` is a Finalize-only field. Before Finalize, its value
  is undefined (default 0.0, but meaningless as std). Readers consult
  `tp.IsFinalized()` or the Result's own `finalized_` flag to know.
  In practice, results only care about Finalized-state data, so the
  flag is rarely checked — but it's there when needed.
- The field name `bs_t0_m2` (Welford's running sum of squared deviations)
  is an intermediate variable exposed as a field. Could be internal to
  the Result, but since Welford's M2 is small (one double per atom) and
  reading it is cheap, keeping it on TrajectoryAtom is simpler than
  indexing into result-internal storage. Calls this decision out: fields
  vs internal state is a per-field judgment.
- The Result class itself has minimal internal state (`n_frames_`,
  `finalized_`). All per-atom data lives on TrajectoryAtom. This matches
  ConformationResult's pattern — the concrete Result is a thin orchestration
  class, not a data container.
- Dependencies declare `BiotSavartResult` — a ConformationResult type.
  The attach-time check reads this declaration; the RunConfiguration
  must include BiotSavartResult in its per-frame calculator set or
  attach of BsWelfordTrajectoryResult fails.

### Concrete example: `BsT0AutocorrelationTrajectoryResult`

A contrasting example: windowed, rich internal state, Finalize-only
output, needs a dense buffer.

```cpp
// BsT0AutocorrelationTrajectoryResult: per-atom autocorrelation of the
// BiotSavart T0 shielding contribution over the trajectory, at a set
// of named lags.
//
// Source: each frame reads conf.atoms[i].bs_shielding_contribution.T0
// and appends to a rolling window of the last N_LAGS frames per atom.
// After the stream ends, compute the autocorrelation function per atom
// per lag and write it to a dense buffer owned by TrajectoryProtein.
//
// No always-valid mid-stream fields; all output is Finalize-only.

class BsT0AutocorrelationTrajectoryResult : public TrajectoryResult {
public:
    static constexpr size_t N_LAGS = 120;  // ~12 ns of lag at 100 ps cadence

    std::string Name() const override { return "BsT0AutocorrelationTrajectoryResult"; }

    std::vector<std::type_index> Dependencies() const override {
        return { std::type_index(typeid(BiotSavartResult)) };
    }

    static std::unique_ptr<BsT0AutocorrelationTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp) override;

    int WriteFeatures(const TrajectoryProtein& tp,
                      const std::string& output_dir) const override;

    // Post-finalize query: per atom, per lag, autocorrelation value.
    double AutocorrelationAt(size_t atom_idx, size_t lag_idx) const;
    // ... convenience methods for spectral density etc.

private:
    // Per-atom rolling window of recent T0 values. std::array of
    // fixed size (N_LAGS), one per atom. Rolling index per atom.
    // This is the biggest internal state.
    std::vector<std::array<double, N_LAGS>> windows_;
    std::vector<size_t> window_head_indices_;  // per atom; (frame_idx % N_LAGS)
    std::vector<double> sum_per_atom_;          // running mean
    std::vector<double> sum_sq_per_atom_;        // running M2

    // Accumulated per-atom per-lag correlation sums.
    // After Finalize, transferred to tp as a dense buffer and
    // normalised to autocorrelation function values.
    // Size: atoms × N_LAGS doubles. Kept internal during streaming,
    // transferred at Finalize.
    std::vector<std::array<double, N_LAGS>> lag_sums_;
    size_t n_frames_ = 0;
    bool finalized_ = false;
};
```

Points from this example:

- The internal state (`windows_`, `lag_sums_`) is large — atoms × N_LAGS
  doubles — but it lives on the Result, not on TrajectoryAtom. Reason: it's
  per-atom per-lag, not per-atom. Per-atom per-lag is a dense 2D structure
  better suited to contiguous storage than to per-atom `std::vector<double>`
  fields.
- Compute updates the internal state each frame, reads only from
  ConformationAtom. It writes nothing to TrajectoryAtom during streaming —
  all its output is Finalize-only.
- Finalize does the real synthesis: compute autocorrelation from lag sums,
  normalise, transfer to a typed dense buffer owned by TrajectoryProtein.
  The Result's `lag_sums_` internal state is moved out (ownership transfer)
  into TrajectoryProtein's dense-buffer pool, keyed by the Result's type.
  Post-Finalize, the Result's `lag_sums_` is empty; the data lives on tp.
- Post-Finalize queries go through the Result's methods (`AutocorrelationAt`),
  which dereference into the dense buffer on tp. Alternatively the caller
  can grab the buffer directly from tp if it wants raw contiguous access.
- The `WriteFeatures` or `WriteH5Group` method knows its own group layout:
  `/trajectory/bs_t0_autocorrelation/lag_values` as a (N_atoms, N_LAGS)
  array, lag time units in a metadata attribute, etc. Schema lives on
  the Result that knows what it's emitting.

### Dense buffers on TrajectoryProtein

When a TrajectoryResult's output is per-atom × (frame | lag | frequency |
...) and the total data is large enough that per-atom separate allocation
is wasteful, the Result transfers ownership of a contiguous buffer to
TrajectoryProtein at Finalize. TrajectoryProtein stores it keyed by the
owning Result's `type_index`; query methods on the Result dereference
into it.

```cpp
template <typename T>
class DenseBuffer {
public:
    DenseBuffer(size_t atom_count, size_t stride_per_atom);
    T& At(size_t atom_idx, size_t offset);
    const T& At(size_t atom_idx, size_t offset) const;
    std::span<T> AtomSlice(size_t atom_idx);
private:
    std::vector<T> storage_;  // atom-major: [atom_i * stride .. (atom_i+1)*stride)
    size_t stride_per_atom_;
};
```

Atom-major layout — per-atom time series is contiguous, good for
per-atom autocorrelation access patterns. Frame-major layouts (for
per-frame-across-atoms analyses) can coexist if needed; same class
template with a different stride shape.

TrajectoryProtein's `AdoptDenseBuffer<T>(buffer, owner_type)` takes
ownership. Queries go through the owning TrajectoryResult's methods, not
direct access to TrajectoryProtein's buffer map (that map is an
implementation detail).

**Non-scalar payloads: see Appendix C.** The `DenseBuffer<T>` shown
above carries `T = double`. The real fleet uses `Vec3` (E-field,
positions), `Mat3` (EFG tensors), `SphericalTensor` (9-component irrep
decomposition). Appendix C commits to the layout convention:
`DenseBuffer<T>` instantiated with the native C++ type
(`DenseBuffer<SphericalTensor>`, `DenseBuffer<Vec3>`, `DenseBuffer<Mat3>`)
since these are memory-contiguous structs with no indirection; at H5
emission, flatten to trailing dimensions with metadata attribute
describing the irrep layout (e.g., `(atom, frame, 9)` for SphericalTensor
with attribute `irrep_layout = "T0,T1_m-1,T1_m0,T1_m1,T2_m-2..T2_m2"`).
This matches how ConformationAtom stores SphericalTensor fields as the
typed struct and emits to NPY as 9 contiguous doubles.

### Attach discipline (mirrors ConformationResult)

```cpp
bool TrajectoryProtein::AttachResult(std::unique_ptr<TrajectoryResult> result) {
    if (!result) return false;

    std::string name = result->Name();
    std::type_index tid(typeid(*result));

    // Singleton check.
    if (results_.find(tid) != results_.end()) {
        OperationLog::Log(OperationLog::Level::Warning, LogResultAttach,
                          "AttachResult",
                          "rejected " + name + ": already attached");
        return false;
    }

    // Dependency check. Dependencies may refer to other TrajectoryResult
    // types (attached here) OR ConformationResult types (which must be in
    // the per-frame RunConfiguration). The RunConfiguration's declared
    // per-frame set is cross-checked at AttachResult when it's available;
    // otherwise, ConformationResult-typed dependencies are validated at
    // the start of Run by Trajectory itself.
    for (const std::type_index& dep : result->Dependencies()) {
        // ... check attached traj results + cross-check run config per-frame set
        // (see Trajectory::Run pre-loop validation in §5)
    }

    results_[tid] = std::move(result);
    results_attach_order_.push_back(results_[tid].get());

    OperationLog::Log(OperationLog::Level::Info, LogResultAttach,
                      "AttachResult",
                      "attached " + name + " to trajectory protein");
    return true;
}
```

Same pattern as `ProteinConformation::AttachResult`. Identical signature,
identical failure modes, identical log levels. A session that knows the
ConformationResult pattern knows this one.

### Selection emission from scan-mode: see Appendix D

In `ScanForDftPointSet` mode, specific TrajectoryResults detect events
(rotamer transitions, RMSD spikes, χ₁ bin crossings, ring-flip hints)
and need to record frame-selection events that the Trajectory entity
collects post-run for downstream DFT submission.

Appendix D specifies the pattern: `SelectionEmittingTrajectoryResult` as
a separate interface (mixin) with a `std::vector<FrameSelectionRecord>`
member + accessor. TrajectoryResults that emit selections inherit both
`TrajectoryResult` and `SelectionEmittingTrajectoryResult`; those that
don't inherit only `TrajectoryResult`. At end of `Trajectory::Run`'s
Finalize phase, Trajectory iterates attached results checking for the
interface via `dynamic_cast` and merges their records into `selections_`.

Key property: the base `TrajectoryResult::Compute` signature stays
clean — no callback parameter, no coupling to Trajectory. Selection is
an opt-in interface, not a cross-cutting concern.

### Per-frame dispatch overhead

Every frame, `Trajectory::Run` iterates `tp.ResultsInAttachOrder()` and
calls each attached TrajectoryResult's virtual `Compute`. Under
`PerFrameExtractionSet` there are ~20–30 attached TrajectoryResults,
one virtual call each, loop body per call mostly doing atom-loop work
inside the concrete Compute. Dispatch overhead is O(10 µs) per frame —
utterly dominated by the 1–10 s that the per-frame ConformationResult
pipeline spends on calculator work (APBS at ~2.7 s, AIMNet2 at ~0.5 s,
ring current + McConnell + Coulomb at ~seconds combined).

Dispatch is not a performance concern. Hot-path considerations stay on
the individual Compute implementations, most of which already
parallelise across atoms (OpenMP / `#pragma omp parallel for`). Adding
a 31st TrajectoryResult to the attach list is free at the dispatch level.

### Catalog of concrete TrajectoryResult subclasses: see Appendix F

§6's `RunConfiguration` factories name the TrajectoryResults attached
in each configuration (13 in `PerFrameExtractionSet`, a handful in
`ScanForDftPointSet`, ~3 additions in `FullFatFrameExtraction`). This
section's two worked examples cover the two patterns (always-valid
Welford and Finalize-only windowed autocorrelation); the other ~25
subclasses follow one of those two patterns.

Appendix F gives the full table: for each TrajectoryResult subclass,
its name, source ConformationResult dependency, emission shape
(TrajectoryAtom fields + dense buffers), always-valid-mid-stream vs
Finalize-only, brief description. The catalog is the implementation
checklist; each row is roughly one class to write.

### Adding a new TrajectoryResult — the modularity checklist

What it takes to add a new trajectory-scope calculator:

1. New `.h/.cpp` pair with `FooTrajectoryResult : public TrajectoryResult`.
2. Override `Name()`, `Dependencies()`, `Compute()`, optionally `Finalize()`,
   optionally `WriteFeatures()` / `WriteH5Group()`.
3. Declare any new fields on TrajectoryAtom this result writes to (the
   upfront-declaration discipline). Or, if the output is a dense buffer,
   register the buffer type.
4. Add the factory method to one or more `RunConfiguration` static
   factories (`ScanForDftPointSet`, `PerFrameExtractionSet`,
   `FullFatFrameExtraction`) — by design choice, which config sets see
   this accumulator.
5. If the new result's WriteFeatures produces a new NPY / H5 group,
   update the Python SDK's `_catalog.py` to register it.

Five points. Analogous to adding a new ConformationResult (one new
`.h/.cpp`, one line in OperationRunner, new fields on ConformationAtom,
SDK catalog update). The modularity property holds.

### Removing a TrajectoryResult — the modularity check in reverse

What it takes to drop a calculator from the library:

1. Remove from the relevant `RunConfiguration` factories.
2. Delete the source files.
3. Remove its fields from TrajectoryAtom (optional — they can live as
   dead fields if we don't want to break serialisation compatibility
   with old H5s; the declared-upfront pattern makes dead fields
   harmless).
4. Remove its SDK catalog entry.

Calculators that depended on the removed one will fail at AttachResult
with a logged warning, be skipped, and their downstream outputs go
unwritten. Callers who read `tp.Result<RemovedType>()` will throw; that's
the contract.

### Taking a TrajectoryResult in and out per run

Within a single binary, a caller controls which TrajectoryResults run
by choosing a RunConfiguration and optionally adding extras via
RunContext. See §6.

---

## 5 — `Trajectory`: the process entity

`Trajectory` is the first-class representation of a traversal run. It
holds both process-role state (what's being traversed, how far through,
which frames have been selected) and record-role state (what was
traversed, what was produced). Process entities can carry both; they
just need to keep the two kinds of state distinct via flags / lifecycle.

`Trajectory` owns no per-atom state — that's TrajectoryProtein's job.
`Trajectory` owns what's about the run itself.

### Structure

```cpp
class Trajectory {
public:
    // Constructor takes source file paths. Records them; does not
    // open them until Run() is called.
    Trajectory(std::filesystem::path xtc_path,
               std::filesystem::path tpr_path,
               std::filesystem::path edr_path);

    // The main method. Drives the traversal:
    //   - Open source files via internal GromacsFrameHandler.
    //   - Read frame 0, invoke gp.FinalizeProtein() to add conf0 to
    //     the wrapped Protein's conformations_ vector.
    //   - Attach TrajectoryResults declared by ctx.Configuration().
    //   - Validate TrajectoryResult dependencies against the per-frame
    //     ConformationResult set named by ctx.Configuration().
    //   - Enter the per-frame loop.
    //   - After loop, Finalize each attached TrajectoryResult.
    //   - Fill in this Trajectory's post-run record fields (frame_count_,
    //     frame_times_, selections_, etc.).
    //
    // Updates tp in place. Post-Run, tp holds trajectory-scope data and
    // this Trajectory holds run metadata.
    void Run(TrajectoryProtein& tp, RunContext& ctx);

    // === Record fields (valid post-Run) ===
    size_t FrameCount() const;
    double TotalTimePs() const;
    const std::vector<double>& FrameTimes() const;  // (frame_count,)
    const std::vector<size_t>& FrameIndices() const; // original XTC indices if strided
    const std::vector<FrameSelectionRecord>& Selections() const;
    const std::filesystem::path& XtcPath() const;
    const std::filesystem::path& TprPath() const;
    const std::filesystem::path& EdrPath() const;

    // State flags.
    bool IsComplete() const { return state_ == State::Complete; }

    // Serialisation. Emits trajectory-scope metadata to H5:
    //   /trajectory/source/xtc_path, tpr_path, edr_path (attributes)
    //   /trajectory/frames/time_ps (T,)
    //   /trajectory/frames/original_index (T,)
    //   /trajectory/selections/... (records of which frames were selected)
    // Does NOT emit per-atom data — that's TrajectoryProtein's
    // WriteH5(), called separately or composed in a higher-level write.
    void WriteH5(HighFive::File& file) const;

private:
    enum class State { Constructed, Running, Complete };
    State state_ = State::Constructed;

    // Sources.
    std::filesystem::path xtc_path_;
    std::filesystem::path tpr_path_;
    std::filesystem::path edr_path_;

    // Record fields.
    size_t frame_count_ = 0;
    std::vector<double> frame_times_;
    std::vector<size_t> frame_indices_;
    std::vector<FrameSelectionRecord> selections_;

    // Process/active state (during Run only).
    std::unique_ptr<GromacsFrameHandler> handler_;
    std::unique_ptr<GromacsRunContextReplacement> active_context_;  // cursor etc.
};

struct FrameSelectionRecord {
    size_t frame_idx;
    double time_ps;
    std::string reason;             // e.g., "rotamer_transition_chi1_LYS_42"
    std::string selector;           // which Accumulator flagged it
};
```

### The Run() loop, concrete

```cpp
void Trajectory::Run(TrajectoryProtein& tp, RunContext& ctx) {
    if (state_ == State::Complete) {
        throw std::logic_error("Trajectory::Run called on already-completed run");
    }
    state_ = State::Running;

    // === Phase 1: attach TrajectoryResults per the RunConfiguration ===
    for (auto& factory : ctx.Configuration().TrajectoryResultFactories()) {
        auto result = factory(tp);
        if (!tp.AttachResult(std::move(result))) {
            // AttachResult logged the reason. Caller's config is wrong;
            // fail loudly.
            throw std::runtime_error("Failed to attach TrajectoryResult");
        }
    }
    // Any ad-hoc extras from RunContext
    for (auto& extra : ctx.MoveOutExtraResults()) {
        if (!tp.AttachResult(std::move(extra))) {
            throw std::runtime_error("Failed to attach extra TrajectoryResult");
        }
    }

    // === Phase 2: validate ConformationResult deps in the run config ===
    for (auto* result : tp.ResultsInAttachOrder()) {
        for (auto& dep : result->Dependencies()) {
            if (ctx.Configuration().RequiresConformationResult(dep)) {
                continue;  // ConformationResult dep declared by config
            }
            if (tp.AllResults().count(dep) > 0) {
                continue;  // TrajectoryResult dep attached
            }
            // Unknown dep — fail loud.
            throw std::runtime_error(
                std::string("TrajectoryResult ") + result->Name() +
                " declares unmet dependency");
        }
    }

    // === Phase 3: open source, read frame 0, finalize protein ===
    handler_ = std::make_unique<GromacsFrameHandler>(tp, xtc_path_, tpr_path_);
    handler_->Open(ctx.PerFrameRunOptions());
    // handler_->Open() reads frame 0 internally, builds conf0, adds it
    // to tp.Protein().conformations_, runs this frame's ConformationResults,
    // and exposes it via handler_->LastConformation().
    const ProteinConformation& conf0 = handler_->LastConformation();

    // Tick TrajectoryResults on frame 0, recording metadata.
    for (auto* result : tp.ResultsInAttachOrder()) {
        result->Compute(conf0, tp, 0, handler_->LastTimePs());
    }
    frame_times_.push_back(handler_->LastTimePs());
    frame_indices_.push_back(0);
    frame_count_ = 1;

    // === Phase 4: per-frame loop ===
    while (handler_->Next(ctx.PerFrameRunOptions())) {
        const ProteinConformation& conf = handler_->LastConformation();
        const size_t idx = frame_count_;
        const double time_ps = handler_->LastTimePs();

        // Each TrajectoryResult acts on this frame. Iteration is
        // attach-order, which respects declared dependencies.
        for (auto* result : tp.ResultsInAttachOrder()) {
            result->Compute(conf, tp, idx, time_ps);
        }

        // Record per-frame metadata on Trajectory.
        frame_times_.push_back(time_ps);
        frame_indices_.push_back(handler_->LastXtcIndex());
        ++frame_count_;

        // Any selections triggered by accumulators during this frame
        // are pulled from tp or from accumulator callbacks — design
        // decision; see RunContext's SelectionCollector for how
        // accumulators report DFT-pose candidates.
    }

    // === Phase 5: Finalize ===
    for (auto* result : tp.ResultsInAttachOrder()) {
        result->Finalize(tp);
    }
    // tp.finalized_ = true (tp's own field)

    state_ = State::Complete;
    handler_.reset();  // Active state is cleared; record fields remain.
}
```

Five clean phases. The loop body of Phase 4 is the "I am on a trajectory,
I am at a frame" phrase made executable: each TrajectoryResult acts on the
frame using the frame's ConformationAtom data + the trajectory-scope
TrajectoryAtom running buffer + its own internal state.

### Frame 0 is permanent; frames 1..N are ephemeral

The existing discipline stays: frame 0 is added to `tp.Protein().
conformations_` via `gp.FinalizeProtein()` during `handler_->Open()`.
Frames 1..N are free-standing `ProteinConformation` objects that live in
`handler_->last_conf_` for one iteration and die at the next `handler_->
Next()` call or at end of loop.

This preserves:
- conf0 as a usable single-conformation reference after the run (e.g.,
  for topology queries, for "the first frame's state").
- Memory bounded at one materialised frame at a time.
- Every ConformationResult that the RunConfiguration declares runs on
  every frame via `OperationRunner::Run` internally to `handler_->Next()`.

### Ordering guarantees

The Run() phases impose a specific order that Compute dispatch and
Finalize dispatch both rely on:

1. **Attach order equals dispatch order.** `ResultsInAttachOrder()`
   returns TrajectoryResults in the order they were attached in Phase 1.
   Phase 4's per-frame loop iterates in that order; Phase 5's Finalize
   loop iterates in that order.
2. **RunConfiguration factories attach first, in source-code order.**
   The factories vector in each `RunConfiguration` static factory
   method is deterministic at compile time, and Phase 1 iterates it
   in that order.
3. **Extras attach after config factories, in call-site order.**
   `RunContext::MoveOutExtraResults()` returns extras in the order
   the caller attached them via `AttachExtraResult`; Phase 1 attaches
   those after the config factories finish.
4. **Dependency validation at attach time.** If a TrajectoryResult
   declares a dependency on another TrajectoryResult type, the
   dependency must already be attached (Phase 1's earlier iterations
   or an earlier extra). If a TrajectoryResult declares a dependency
   on a ConformationResult type, that's validated against
   `RunConfiguration.RequiresConformationResult()` in Phase 2.

This means a config-default TrajectoryResult can depend on another
config-default (the factories are ordered to satisfy deps). An extra
can depend on a config-default. A config-default depending on an extra
would fail attach — but that's an authorial error (extras should be
extras, not dependencies); the design surfaces the problem immediately.

Callers inspecting `tp.ResultsInAttachOrder()` post-Run see exactly
what ran, in exactly the order it ran. Serialisation orders follow
suit.

### EDR preload timing (GromacsRunContext dissolution)

With `GromacsRunContext` dissolved (§9), the preloaded EDR energy
frames and the bonded parameters from TPR need a new home. Concrete
assignment:

- **EDR energy frames** (~100 MB for a 25 ns trajectory, preloaded as
  `std::vector<GromacsEnergy>`): owned by `Trajectory` as a private
  member. Loaded in `Trajectory`'s constructor, not during `Run()` —
  opening the file and preloading happens before any traversal starts,
  so `Run()` has the data available from phase 3 onward.
- **Bonded parameters from TPR**: these describe protein topology
  (bond / angle / dihedral / CMAP terms with force-field constants).
  They belong naturally on TrajectoryProtein (topology-scope, set once
  per protein). Loaded by `GromacsEnsembleLoader::BuildFromTpr` during
  TrajectoryProtein construction, not during `Run()`.
- **Cursor state** (current frame index, current energy pointer): lives
  inside `GromacsFrameHandler` as private state during Run() only. No
  persistent owner.

During Run()'s per-frame loop, `GromacsFrameHandler` asks Trajectory
for the current frame's energy via a typed accessor
(`traj.EnergyAtTime(time_ps)` — O(log T) binary search into the
preloaded vector). The handler then populates
`RunOptions.frame_energy` before calling `OperationRunner::Run` on
the frame, so `GromacsEnergyResult::Compute` receives the energy
through the usual RunOptions channel. No changes to the
ConformationResult pipeline; only the plumbing behind
`frame_energy`'s source changes.

---

## 6 — `RunConfiguration` and `RunContext`

### `RunConfiguration` as a first-class typed object

Three named run shapes. Each is a concrete instance of `RunConfiguration`
produced by a static factory method. Each names its per-frame
ConformationResult set (via `RunOptions` skip flags), its
TrajectoryResult factory list, and any frame-selection policy.

```cpp
class RunConfiguration {
public:
    using TrajectoryResultFactory =
        std::function<std::unique_ptr<TrajectoryResult>(const TrajectoryProtein&)>;

    // The three named shapes.
    static RunConfiguration ScanForDftPointSet();
    static RunConfiguration PerFrameExtractionSet();
    static RunConfiguration FullFatFrameExtraction();

    // Accessors used by Trajectory::Run and validators.
    const std::vector<TrajectoryResultFactory>& TrajectoryResultFactories() const;
    const RunOptions& PerFrameRunOptions() const;
    bool RequiresConformationResult(std::type_index tid) const;
    const std::string& Name() const;  // "ScanForDftPointSet" etc., for logs

private:
    std::string name_;
    RunOptions per_frame_opts_;
    std::vector<TrajectoryResultFactory> traj_factories_;
    std::unordered_set<std::type_index> required_conf_result_types_;
};
```

The three named configurations materialise as follows. Concrete, not
abstract — each one is what the library actually does for that mode
today (mapping the current scan-mode vs analysis-mode code paths onto
the new typed object).

```cpp
RunConfiguration RunConfiguration::ScanForDftPointSet() {
    RunConfiguration c;
    c.name_ = "ScanForDftPointSet";

    // Per-frame: very cheap. No APBS, no MOPAC, no Coulomb, no AIMNet2.
    c.per_frame_opts_.skip_dssp = false;   // need DSSP for Ramachandran bins
    c.per_frame_opts_.skip_mopac = true;
    c.per_frame_opts_.skip_apbs = true;
    c.per_frame_opts_.skip_coulomb = true;
    c.per_frame_opts_.aimnet2_model = nullptr;
    // Calculators that DO run each frame (not controlled by skip flags
    // because they're always-on unless their gate fails): Geometry,
    // SpatialIndex, Enrichment, DSSP, BiotSavart (for ring proximity),
    // HaighMallion, McConnell, SASA.

    c.required_conf_result_types_ = {
        typeid(GeometryResult),
        typeid(SpatialIndexResult),
        typeid(DsspResult),
        typeid(SasaResult),
    };

    c.traj_factories_ = {
        // Cheap accumulators for scan-mode: dihedral bins, rotamer
        // transitions, Ramachandran-region transitions, RMSD tracking,
        // ring-flip detection.
        [](const TrajectoryProtein& tp) {
            return DihedralBinTransitionTrajectoryResult::Create(tp);
        },
        [](const TrajectoryProtein& tp) {
            return RmsdTrackingTrajectoryResult::Create(tp);
        },
        [](const TrajectoryProtein& tp) {
            return DftPoseSelectionTrajectoryResult::Create(tp);
        },
        // etc.
    };

    return c;
}

RunConfiguration RunConfiguration::PerFrameExtractionSet() {
    RunConfiguration c;
    c.name_ = "PerFrameExtractionSet";

    // Per-frame: the full classical stack, no MOPAC.
    c.per_frame_opts_.skip_dssp = false;
    c.per_frame_opts_.skip_mopac = true;
    c.per_frame_opts_.skip_apbs = false;
    c.per_frame_opts_.skip_coulomb = false;
    // AIMNet2 set by caller if available.

    c.required_conf_result_types_ = {
        typeid(GeometryResult), typeid(SpatialIndexResult),
        typeid(EnrichmentResult), typeid(DsspResult),
        typeid(ChargeAssignmentResult), typeid(ApbsFieldResult),
        typeid(BiotSavartResult), typeid(HaighMallionResult),
        typeid(McConnellResult), typeid(CoulombResult),
        typeid(HBondResult), typeid(RingSusceptibilityResult),
        typeid(PiQuadrupoleResult), typeid(DispersionResult),
        typeid(SasaResult), typeid(WaterFieldResult),
        typeid(HydrationShellResult), typeid(HydrationGeometryResult),
        typeid(EeqResult), typeid(BondedEnergyResult),
    };

    c.traj_factories_ = {
        // The analysis-H5-producing family. Time series of every
        // per-atom tensor field + Welford accumulators + per-residue
        // dihedral time series + DSSP time series.
        [](const TrajectoryProtein& tp) { return PositionsTimeSeriesTrajectoryResult::Create(tp); },
        [](const TrajectoryProtein& tp) { return BsWelfordTrajectoryResult::Create(tp); },
        [](const TrajectoryProtein& tp) { return BsShieldingTimeSeriesTrajectoryResult::Create(tp); },
        [](const TrajectoryProtein& tp) { return McWelfordTrajectoryResult::Create(tp); },
        [](const TrajectoryProtein& tp) { return McShieldingTimeSeriesTrajectoryResult::Create(tp); },
        [](const TrajectoryProtein& tp) { return CoulombFieldTimeSeriesTrajectoryResult::Create(tp); },
        [](const TrajectoryProtein& tp) { return ApbsFieldTimeSeriesTrajectoryResult::Create(tp); },
        [](const TrajectoryProtein& tp) { return WaterEnvironmentTimeSeriesTrajectoryResult::Create(tp); },
        [](const TrajectoryProtein& tp) { return DihedralTimeSeriesTrajectoryResult::Create(tp); },
        [](const TrajectoryProtein& tp) { return DsspTimeSeriesTrajectoryResult::Create(tp); },
        [](const TrajectoryProtein& tp) { return BondedEnergyTimeSeriesTrajectoryResult::Create(tp); },
        [](const TrajectoryProtein& tp) { return GromacsEnergyTimeSeriesTrajectoryResult::Create(tp); },
        // Rollup statistics
        [](const TrajectoryProtein& tp) { return SasaWelfordTrajectoryResult::Create(tp); },
        [](const TrajectoryProtein& tp) { return HBondCountTrajectoryResult::Create(tp); },
        // etc. — one factory per TrajectoryResult class that this
        // configuration wants attached.
    };

    return c;
}

RunConfiguration RunConfiguration::FullFatFrameExtraction() {
    RunConfiguration c = PerFrameExtractionSet();  // inherit
    c.name_ = "FullFatFrameExtraction";

    // Add MOPAC.
    c.per_frame_opts_.skip_mopac = false;
    c.required_conf_result_types_.insert(typeid(MopacResult));
    c.required_conf_result_types_.insert(typeid(MopacCoulombResult));
    c.required_conf_result_types_.insert(typeid(MopacMcConnellResult));

    // Add MOPAC-family TrajectoryResults.
    c.traj_factories_.push_back(
        [](const TrajectoryProtein& tp) { return MopacCoulombTimeSeriesTrajectoryResult::Create(tp); });
    c.traj_factories_.push_back(
        [](const TrajectoryProtein& tp) { return MopacMcConnellTimeSeriesTrajectoryResult::Create(tp); });
    c.traj_factories_.push_back(
        [](const TrajectoryProtein& tp) { return MopacVsFf14SbReconciliationTrajectoryResult::Create(tp); });

    return c;
}
```

**MOPAC is sparse-frame, not every-frame.** PM7 + MOZYME takes ~45 s
per ~889-atom protein per conformation (per `CLAUDE.md`). Every-frame
MOPAC on a 25 ns trajectory (~1200 frames) is 15 hours per protein; on
the 685-protein fleet it's hundreds of core-years. **`FullFatFrameExtraction`
is the configuration for selected frames only** — the 260 DFT pose set
for calibration today, checkpoint frames during the μs harvester phase
later. The selected-frame discipline is enforced at the caller layer
via one of two mechanisms:

1. *Extract a frame-filtered XTC.* A preceding `ScanForDftPointSet` run
   produces `FrameSelectionRecord`s; the calling code extracts a new XTC
   containing only selected frames and runs `FullFatFrameExtraction` on
   that. Clean separation; current fleet uses this pattern.
2. *Selected-frame skip in the frame handler.* `RunContext::SetSelectedFrameIndices(
   const std::set<size_t>&)` tells `GromacsFrameHandler` to skip frames
   not in the set (read, discard without running the calculator
   pipeline). Less I/O overhead, no intermediate file. Needs a small
   addition to GromacsFrameHandler; flagged for implementation as an
   open decision (§12 item 6, added during this revision).

The per-frame cost note applies to MOPAC and to
`MopacVsFf14SbReconciliationTrajectoryResult` specifically; the other
two `PerFrameExtractionSet` calculators this config inherits stay
per-frame (they're cheap — milliseconds per frame even with rich tensor
output).

Each configuration is a small, readable static factory. "What does this
run do?" is answered by reading one function. Changing a configuration
means changing one function. Adding a new configuration means writing a
new static factory.

See **Appendix F** for the full catalog of TrajectoryResult subclasses
referenced by these factories: one row per class with its dependencies,
emission shape, and lifecycle (always-valid-mid-stream vs Finalize-only).

### `RunContext` as a class with a constructor

Not a struct. Takes its required arguments via constructor; lets caller
attach ad-hoc extras via methods.

```cpp
class RunContext {
public:
    // Caller supplies a configuration and an output destination.
    // Optional extras can be added after construction.
    RunContext(RunConfiguration config,
               std::filesystem::path output_dir);

    // Add an ad-hoc TrajectoryResult beyond the config's defaults.
    // Common for tests, for ad-hoc analyses, for the μs harvester's
    // custom selectors.
    void AttachExtraResult(std::unique_ptr<TrajectoryResult> extra);

    // Set AIMNet2 model (optional).
    void SetAimnet2Model(AIMNet2Model* model);

    // Accessors.
    const RunConfiguration& Configuration() const;
    const std::filesystem::path& OutputDir() const;
    const RunOptions& PerFrameRunOptions() const;  // delegates to config, with aimnet2 patched in

    // Moved out by Trajectory::Run when it attaches extras.
    std::vector<std::unique_ptr<TrajectoryResult>> MoveOutExtraResults();

    // Room to grow: validation, derived field lookups, AttachSelector,
    // SetFrameSelectionPolicy, etc.

private:
    RunConfiguration config_;
    std::filesystem::path output_dir_;
    std::vector<std::unique_ptr<TrajectoryResult>> extras_;
    AIMNet2Model* aimnet2_model_ = nullptr;
    // No mutable cache. PerFrameRunOptions() returns by value,
    // composed fresh at call time from config_ + aimnet2_model_. Called
    // once per Run() setup (Phase 1 / 2 use), plus once per frame as
    // Trajectory passes it into the handler — the copy cost is trivial
    // (a struct of pointers and booleans). Avoids the const-correctness
    // smell of a mutable cache.
};
```

Construction reads as intent, not assembly:

```cpp
RunContext ctx(RunConfiguration::PerFrameExtractionSet(),
               "output/1ubq");
ctx.SetAimnet2Model(aimnet2_model);
// Optional extras:
ctx.AttachExtraResult(
    std::make_unique<CustomSelectionTrajectoryResult>(custom_criterion));

Trajectory traj(xtc_path, tpr_path, edr_path);
TrajectoryProtein tp = GromacsEnsembleLoader::BuildFromTpr(tpr_path);
traj.Run(tp, ctx);
```

---

## 7 — Writer semantics and the model-is-source-of-truth discipline

There is no separate writer class holding buffers. The model IS the
source of truth during and after streaming. A writer is a function
that traverses the model and emits.

### The pattern

1. During streaming, TrajectoryResults write their data to TrajectoryAtom
   fields (always-valid or Finalize-only), to their own internal state,
   and via Finalize to dense buffers owned by TrajectoryProtein.
2. After streaming, the writer asks the model to serialise itself:
   - `tp.WriteH5(file)` writes the TrajectoryProtein — iterates
     `tp.ResultsInAttachOrder()` and calls each result's
     `WriteH5Group(tp, file)`. Plus writes TrajectoryAtom-level fields
     across all atoms as a single "atoms" group.
   - `traj.WriteH5(file)` writes the Trajectory record — source paths,
     frame times, frame indices, selections.
   - Top-level wrapper function calls both:
     ```cpp
     void WriteTrajectoryH5(const Trajectory& traj,
                            const TrajectoryProtein& tp,
                            const std::filesystem::path& output_h5_path);
     ```

### Layers leave things for other layers — preserved

The discipline holds: a layer 2 TrajectoryResult that depends on layer 1
TrajectoryResult's output reads from the model (TrajectoryAtom fields
or tp's dense buffers), not from some writer's buffer. Example:

```cpp
// Layer 1: BsWelfordTrajectoryResult writes bs_t0_mean per atom.
// Layer 2: BsAnomalousShiftTrajectoryResult reads bs_t0_mean per atom
// and flags atoms whose value is > 2 std from the protein's median.
// Layer 2 doesn't need a writer — it reads layer 1's field directly
// in its Compute (per frame) or in its Finalize (once).
```

Both layers write to the model; the writer at end traverses the model.

### AnalysisWriter dissolves

The current `AnalysisWriter` buffers per-frame `(T, N, ...)` arrays in
its own vectors. Those buffers migrate into TrajectoryResult objects:
`PositionsTimeSeriesTrajectoryResult` owns a dense buffer of size
(atoms × frames × 3); it's populated each frame by Compute (reading
conf.AtomAt(i).position()); at Finalize it transfers ownership to tp;
at WriteH5Group it writes `/positions/xyz` as a 3D dataset.

AnalysisWriter the class is not needed. Its logic migrates into:
- Multiple `*TimeSeriesTrajectoryResult` classes, one per current H5 group.
- Each owns its dense buffer during streaming.
- Each emits its own H5 group.

The top-level `WriteTrajectoryH5` orchestrates by iterating attached
results.

### Structure-ID emission (per-protein provenance)

Beyond per-atom `NmrAtomIdentity`, the H5 output commits to emitting
per-protein structure provenance in fixed group paths that downstream
consumers (SDK, h5-reader, learn/ calibration, external analysts)
can rely on. The atom-names-and-structure-ids pair is the foundation
of every downstream typed binding.

**H5 group `/metadata/source/` — per-Protein provenance.** Populated
from `ProteinBuildContext` (see `OBJECT_MODEL.md`). Attributes emitted
at group level:

- `pdb_source` — string, PDB deposition ID if known (e.g., `"1UBQ"`)
  or file path
- `deposition_date` — string, ISO format (`YYYY-MM-DD`), if known
- `organism` — string, if known
- `crystal_resolution` — double, Angstroms (NaN if not crystallographic)
- `r_factor` — double (NaN if not crystallographic)
- `temperature_kelvin` — double (NaN if not known)
- `protonation_tool` — string, e.g., `"PROPKA 3.5.1"`, `"KaML-CBTrees"`,
  `"Manual"`
- `protonation_pH` — double (NaN if not pH-based)
- `force_field` — string, e.g., `"Amber_ff14SB"`, `"CHARMM36m"`
- `stripped` — string array, e.g., `["waters", "heteroatoms"]`
- `assumptions` — string array, e.g., `["missing_loops_rebuilt"]`
- `extractor_version` — string, commit SHA of the nmr-extract
  library at extraction time
- `extraction_date` — string, ISO timestamp of the extraction run

Emission is handled by a small writer function called once per
extraction (not a ConformationResult — this is protein-level and
pre-calculator). The contract: every H5 file emitted by the library
carries this group; absence of a field is represented by NaN /
empty-string, never by dataset absence.

**H5 group `/atoms/identity/` — per-atom typed identity.** Populated
from `NmrAtomIdentity` at Protein construction (unchanged across
conformations; emitted once). Datasets indexed by atom index; one
dataset per identity field:

- `residue_index` (int) — into `protein.Residues()`
- `residue_type` (int) — `AminoAcid` enum value
- `residue_variant` (int) — HIE / HID / HIP / CYS / CYX / ALA / ... enum
- `residue_sequence_number` (int) — PDB sequence number
- `chain_id` (string) — PDB chain
- `insertion_code` (string) — PDB column 27
- `pdb_atom_name` (string) — for display / provenance only; not the
  binding key
- `iupac_atom_position` (int) — `IupacAtomPosition` enum value; **the
  primary binding key for external data**
- `nmr_class` (int) — `NmrClass` enum value
- `locant` (int) — α / β / γ / δ / ε / ζ / η / Backbone / NotApplicable
- `sidechain_depth_bonds` (int)
- `methyl_group` (int) — `MethylGroup` enum value or `None`
- `methyl_partner_atom_index` (int64) — SIZE_MAX if not applicable
- `methylene_group` (int)
- `methylene_partner_atom_index` (int64)
- `chi_participation` (uint32) — bitfield for χ₁..χ₄
- `ring_atom_role` (int) — `RingAtomRole` enum value
- `ring_membership` (uint64) — bitfield over `protein.Rings()`
- `symmetry_class_id` (int)
- `pseudo_atom_name` (string) — BMRB-collapsed form (e.g., `HB#`), edge
  only; not used in C++ logic
- `residue_category` (int) — `ResidueCategory` enum
- `residue_hydropathy` (double)
- `residue_volume` (double)
- `residue_isoelectric_point` (double)
- `residue_hbond_donor_count` (int)
- `residue_hbond_acceptor_count` (int)
- `is_canonical` (bool) — false if the residue was flagged
  `NonCanonical` during NmrAtomIdentity derivation (per §2's
  per-residue full-template-matching discipline); downstream code
  respects this flag

Each enum-valued dataset carries a sidecar attribute `enum_values`
listing value names in order (e.g., on `nmr_class`:
`"BackboneHN=0,BackboneHA=1,BackboneN=2,..."`) so consumers can
decode without hardcoded enum mappings. This makes the SDK and any
external consumer robust to future enum additions (enum additions
append, never renumber).

**Schema stability commitment.** These group paths, attributes, and
dataset names are committed as stable across the rollup activation.
Schema changes post-activation require an amendment entry
documenting the change and are propagated to the SDK's `_catalog.py`
in the same activation commit. The byte-parity test strategy in §10
guards value-level stability; this commitment guards path-level
stability.

### The library ↔ SDK boundary

With the rollup activated, the library commits to a specific output
contract that downstream work reads against. Naming the boundary
explicitly, rather than leaving it implicit in scattered
WriteFeatures calls, prevents "what does the SDK see?" drift.

**The library (nmr-extract) produces, per extraction:**

1. Typed per-atom identity (`/atoms/identity/` — see subsection above).
2. Typed per-Protein provenance (`/metadata/source/`).
3. Per-conformation ConformationResult outputs (for static paths:
   the existing calculators' NPY / H5 output, unchanged — preserved
   by §0 non-regression contract).
4. Per-trajectory TrajectoryResult outputs (for trajectory paths:
   TrajectoryAtom fields + dense buffers + TrajectoryResult-group
   H5 emissions, per Appendix F catalog).
5. `Trajectory` record metadata (`/trajectory/` — frame times, frame
   indices, `FrameSelectionRecord`s from scan-mode runs).
6. Topology (`/topology/` — bonds, rings, parent-atom mappings,
   existing, unchanged).

**The library does not:**

- Ingest external experimental data (BMRB `.str` files, RefDB shift
  tables, published CSA tables, lanthanide-tag χ tensors). Library
  emits identity + provenance; SDK performs joins.
- Run ML training or inference beyond `AIMNet2Result`'s per-atom
  charges + embedding (which is a typed ConformationResult, not an
  outside pipeline).
- Render anything (`ui/`, `h5-reader/` are downstream consumers).
- Host Python extensions as part of the extraction path (hard rule
  from project discipline; no subprocess, no cpython binding into
  the library).
- Attempt to infer typed identity in the presence of disagreement
  with reference sources. If the validation discipline (§2, Appendix A)
  fails, the extraction fails — it does not silently fall through.

**What SDK work this enables (downstream of the rollup):**

- **BMRB / RefDB shift binding** via per-atom
  `iupac_atom_position` tuple + `residue_sequence_number` +
  `chain_id`. The identity is already typed; the SDK just reads the
  H5's `/atoms/identity/` group, opens a `.str` file, joins on the
  tuple, returns typed shift records per atom.
- **Validation-bench back-calculation.** L7 Yao-Bax α/β 15N split,
  L12 Loth ubiquitin 64-amide CCR, M7 Kasinath-Wand methyl entropy,
  N3 Babaei bulk χ, etc. Each bench needs typed per-atom slicing
  (by `nmr_class`, `methyl_group`, `residue_category`) + per-atom
  data (from `/atoms/identity/` + calculator output groups). SDK
  writes these back-calculations as analysis modules.
- **Cross-protein pooling for calibration.** `learn/stage1-mutations/`
  and successors read multiple proteins' H5s, slice each by typed
  identity, accumulate per-stratum statistics. Stratum keys are
  `iupac_atom_position` or `nmr_class` or `methyl_group`, all typed
  integers the SDK reads once per protein.
- **Per-atom residual computation** (predicted shift - experimental).
  SDK holds both the library's predicted shielding tensor + the
  BMRB/RefDB experimental shift after binding; residual computation
  is SDK-native arithmetic.

**What the SDK doesn't do:**

- Doesn't compute physics. Everything computable is in the library's
  output; SDK reads + joins + correlates.
- Doesn't modify H5 files it reads (read-only consumer).
- Doesn't second-guess `is_canonical`. If a residue is flagged
  NonCanonical by the library's NmrAtomIdentity derivation, SDK sees
  the flag and respects it (usually: exclude from typed-slicing; still
  include in raw-numeric access if the consumer is doing something
  name-based).
- Doesn't run against H5 files from pre-activation library versions
  without an explicit compatibility mode (post-activation, the
  schema in §0–§7 is the contract; older files may be readable with
  a compat mode but the main SDK path targets the activated schema).

**Post-activation roadmap (SDK-side, out of rollup scope but
enabled by it):**

- `python/nmr_extract/experimental/bmrb.py` — BMRB `.str` ingestion
  + join against `/atoms/identity/`.
- `python/nmr_extract/benches/` — one module per validation bench
  (Yao-Bax, Loth, Kasinath-Wand, Babaei, ...), each doing its
  back-calculation against library output + BMRB / published data.
- `python/nmr_extract/stratify/` — typed-slicing utilities for
  per-stratum calibration and analysis.

None of these SDK modules exist yet. The rollup's library output is
what enables them; the SDK work happens post-fleet-activation in a
separate sequence.

---

## 8 — Per-frame vs per-trajectory, concrete

Not a section of new design; a consolidation of how the pattern handles
the existing scan-mode vs analysis-mode distinction.

### Scan mode (ScanForDftPointSet)

Per-frame: cheap ConformationResults (Geometry, SpatialIndex,
Enrichment, DSSP, BiotSavart, HaighMallion, McConnell, SASA). No APBS,
no MOPAC, no Coulomb, no AIMNet2. Fast.

Attached TrajectoryResults: DihedralBinTransitions (per-residue χ/φ/ψ
transitions), RmsdTracking (running RMSD from frame 0 or reference),
DftPoseSelection (aggregator that fires when a transition / novelty
threshold is crossed, records the frame_idx and reason into the
Trajectory's `selections_` vector via a callback).

Output: TrajectoryProtein has Welford summaries on TrajectoryAtom for
the running stats. Trajectory has `selections_` full of records. Writer
emits a "scan summary" H5 + a list of selected frames to re-extract
under FullFatFrameExtraction.

### Analysis mode (PerFrameExtractionSet)

Per-frame: full classical stack. Takes ~seconds per frame at 4000 atoms.

Attached TrajectoryResults: the big family of `*TimeSeriesTrajectoryResult`
plus the corresponding Welford rollups.

Output: exhaustive per-frame per-atom H5, ~1 GB per 25 ns trajectory at
~4000 atoms.

### Full-fat mode (FullFatFrameExtraction)

Per-frame: analysis + MOPAC (sparse selected frames — or across all if
time allows; typically 260 DFT pose set).

Attached TrajectoryResults: analysis set + MopacTimeSeries family +
MopacVsFf14SbReconciliation diagnostic (per spec desiderata D5).

Output: analysis H5 with MOPAC fields populated for the subset of frames
where MOPAC ran.

---

## 9 — Migration: renames, dissolutions, bones

What moves where.

### Renames (bold; git is the safety net)

| Current | Renamed | Reason |
|---|---|---|
| `src/GromacsProtein.{h,cpp}` | `src/TrajectoryProtein.{h,cpp}` | Generic role, not format-specific |
| `class GromacsProtein` | `class TrajectoryProtein` | |
| `src/GromacsProteinAtom.h` | `src/TrajectoryAtom.h` | |
| `class GromacsProteinAtom` | `class TrajectoryAtom` | |

Call sites (`GromacsEnsembleLoader`, `GromacsFrameHandler`, `nmr_extract.cpp`,
tests) update en masse.

### Stays Gromacs-named (format-specific, role is XTC/TPR/EDR)

- `src/GromacsFrameHandler.{h,cpp}` — per-frame XTC reader.
- `src/FullSystemReader.{h,cpp}` — TPR parser + frame splitter.
- `src/GromacsEnsembleLoader.{h,cpp}` — builds Protein + TrajectoryProtein
  from Gromacs sources.
- `src/GromacsEnergyResult.{h,cpp}` — parses EDR frames.
- `src/SolventEnvironment.h` — per-frame water/ion snapshot, used by
  water calculators.

### Dissolutions (these classes/methods/patterns are removed)

- `GromacsProteinAtom`'s flat bag of Welford fields — replaced by
  TrajectoryAtom with fields written by attached TrajectoryResults.
- `GromacsProteinAtom::AllWelfords()` — gone. The per-field enumeration
  for the rollup H5 becomes a per-TrajectoryResult `WriteH5Group`.
- `GromacsProtein::AccumulateFrame(conf, idx, time)` — gone. The frame
  loop calls each attached TrajectoryResult's `Compute` directly.
- `GromacsProtein::atoms_`, `GromacsProtein::bonds_accum_` — gone.
  Per-atom state lives on `TrajectoryProtein::atoms_` (of type
  `TrajectoryAtom`). Per-bond accumulators become a Bond-indexed
  TrajectoryResult if retained.
- `AnalysisWriter` the class — gone. Its per-frame buffers migrate into
  individual `*TimeSeriesTrajectoryResult` classes.
- `GromacsRunContext` — gone. Its contents redistribute:
  - Bonded parameters (from TPR) → owned by TrajectoryProtein or
    Trajectory as a preload, accessed by ConformationResults that need
    them (BondedEnergyResult) via RunOptions per frame.
  - Preloaded EDR frames → owned by Trajectory as a preload, accessed
    by GromacsEnergyResult via RunOptions.frame_energy per frame.
  - Cursor state (frame_index, frame_time_ps) → lives in
    GromacsFrameHandler's internal state during Run only; not a
    persistent member of anyone.

### Bones migrations

The dissolved classes' source files move to `learn/bones/` (the project's
historical-code graveyard) with brief migration notes attached. Bonus:
reading a future session's bones is how they see "here is the shape we
evolved away from."

### Non-trivial downstream updates

- Python SDK `python/nmr_extract/_catalog.py` — every current
  array catalog entry updates to match the new H5 group layout. Much
  of the layout doesn't change (`/positions/xyz` stays at the same
  path); what changes is the provenance ("written by PositionsTimeSeriesTrajectoryResult"
  instead of "written by AnalysisWriter").
- Python SDK wrapper classes — unchanged where fields are unchanged;
  updated import paths.
- h5-reader — reads H5; the schema changes outlined here will update
  what the reader sees. Its parallel QtProtein hierarchy stays but as
  a Qt wrapper, not a parallel type hierarchy (per earlier design
  discussions).
- Test fixtures that reference `GromacsProtein` or `GromacsProteinAtom`
  by name need updating.

---

## 10 — Coordinated activation

This rollup activates as a single coordinated library update before the
685-protein fleet extraction, per
`spec/ChangesRequiredBeforeProductionH5Run.md`. Bundled with:
- NamingRegistry activation (the 8 CHARMM↔IUPAC rules and the ALA
  wildcard narrowing, already deferred and documented).
- NmrAtomIdentity introduction on Protein (§2 above).
- Accumulator / TrajectoryResult refactor (§3, §4, §5, §6).
- AnalysisWriter dissolution (§7).
- Renames (§9).
- Python SDK `_catalog.py` update (§9).

No partial activation. The 10-protein calibration H5s stay as-is (not
regenerated in this pass — per the existing discipline, regeneration is
physically possible later because the ORCA output is preserved, but not
on the table for this activation).

Validation:
1. Build library with all changes.
2. Run test suite.
3. Run extraction on 1 calibration protein under PerFrameExtractionSet;
   verify output H5 has expected groups + fields + values match old
   extraction. See **Test strategy** below for the acceptance criteria.
4. Spot-check on 2-3 more calibration proteins.
5. Activate for fleet.

### Test strategy for byte-parity validation

The activation criterion is: **bit-identical output for Welford-family
fields, byte-identical output for time-series fields, explicit
documented changes for schema evolution.** Any ε-tolerance exemption
needs justification and a specific numerical reason.

The reasoning: Welford's M2 is associative enough that a frame-order-
identical input produces numerically-identical output — mean and M2
commute under the Welford update. But any reordering of frame tick
calls (if the new dispatch iterates a different order than the old
`GromacsProtein::AccumulateFrame` sweep) surfaces floating-point
differences at the ULP level. That's a reordering bug, not a
correctness improvement, and must be caught before fleet activation.

Test approach:

1. **Baseline capture.** Seed-stable calibration-protein extraction
   under the old pipeline → capture full H5 output + any NPY outputs +
   any CSV rollups. Archive as read-only fixtures.
2. **New-pipeline extraction.** Same seed, same inputs, new pipeline
   → capture equivalent outputs.
3. **Strict diff.** `h5diff` with 0 tolerance; `numpy.array_equal` on
   NPY arrays; byte-diff on CSV. Any mismatch reports as either:
   a. A reordering bug (must fix to preserve parity — most common case).
   b. An intentional schema change (documented in Appendix G diff,
      justified in the activation artifact, accepted).
   c. A numerical-precision change from a recomputation whose new form
      is more numerically stable — acceptable iff justified with a
      specific reason (e.g., "new pipeline uses Welford on the whole
      series; old pipeline averaged in batches and accumulated roundoff
      differently"). Requires review.
4. **Multi-protein spot-check.** Repeat on 2–3 additional calibration
   proteins. Same bit-identical expectation.
5. **Test-suite additions.** Two new test categories:
   a. **Welford reordering test.** Feed a TrajectoryResult synthetic
      frame data in two different attach orders; assert the final
      per-atom field values are bit-identical to each other. Catches
      accidental order-sensitivity.
   b. **Dense-buffer roundtrip test.** TrajectoryResult populates a
      dense buffer → Finalize → TrajectoryProtein adopts → WriteH5 →
      reload into memory → byte-compare. Catches serialisation layout
      bugs early.

### Schema changes — documented where they happen, not as a separate artifact

Decision 2026-04-22: no dedicated schema-diff artifact for fleet
activation. User's framing verbatim: *"yes it can come out on the
diff. Good idea, not for us. If I wanted file ver I would have spec'd
it. Thought about it but honestly I know what version I launch for a
3 week run."* A formal schema-diff process (originally Appendix G,
now removed) is corporate overhead that two-collaborator + git doesn't
need.

Where schema changes actually live during activation:
- Commit messages for the rollup commits describe schema changes
  inline.
- WIP Amendments section records each schema change with a one-line
  entry (path changed, reason).
- The existing test suite's dataset-existence + shape-assertion tests
  catch accidental drops / shape mutations, same as for any library
  change.
- The byte-parity test on a calibration protein (above) catches value
  changes where old pipeline → new pipeline should produce identical
  output.

No separate `spec/fleet_activation/schema_diff_*.md` artifact. No gate
on "reviewer signs off on the diff document." The activation commit
and its message are the record.

### Scope and timeline

This rollup is **12–20 implementation sessions at ~60 hour/week pace**
per user estimate 2026-04-22. That cost is cheap relative to the
downstream work it enables: every library-level incorrectness
discovered post-fleet costs weeks of re-extraction across 685 proteins
(at ~25 min/protein per `CLAUDE.md`, a full re-extraction is ~285
hours of compute); every library-level correctness captured pre-fleet
compounds across all remaining thesis work (calibration, Stage 2/3
model training, μs harvester). The design-level investment is
deliberate and up-front.

**No scope-tightening**; activation bundle stays whole. The fleet
timing adapts to the bundle, not the other way around. User
explicitly framed this as *"our working on this hard a week is nothing
compared to all the stuff that has to run."* That framing is the
scoping principle for this rollup.

---

## 11 — Patterns this document defines (for PATTERNS.md update)

Once this design settles and folds, the following patterns move into
`PATTERNS.md` as peer sections to the existing "ConformationAtom:
private construction, typed fields", "Typed results with dependency
declarations", etc.:

- **TrajectoryAtom: private construction, typed fields, running buffer.**
  Same pattern as ConformationAtom at conformation scope, scope-shifted
  to trajectory. Private ctor; only TrajectoryProtein constructs.
  Public typed fields declared upfront. One writer per field, enforced
  by the singleton-per-type attachment of TrajectoryResults. No
  identity duplication (flow through wrapped Protein).
- **TrajectoryResult: multi-phase modular calculator with declared
  dependencies.** Four virtuals (Name, Dependencies, Compute, Finalize;
  plus optional WriteFeatures/WriteH5Group). Compute called per frame
  with (conf, tp, frame_idx, time_ps). Finalize called once at end.
  Dependencies mix TrajectoryResult and ConformationResult types.
  Add = write subclass + register in RunConfiguration. Remove = drop
  from RunConfiguration.
- **RunConfiguration as named typed object, not switch.** Three concrete
  configs produced by static factories. Each encodes (per-frame RunOptions,
  TrajectoryResult factory list, required ConformationResult set). New
  configurations are new factories; existing ones are read as ordinary
  functions.
- **Trajectory::Run as the orchestrator.** Five phases: attach
  TrajectoryResults, validate deps, open and read frame 0 + tick, per-frame
  loop, Finalize. No intermediate sinks; model-is-source-of-truth.
- **Writer as model traverser, not owner of buffers.** tp.WriteH5(file)
  iterates attached TrajectoryResults and calls each one's WriteH5Group.
  Trajectory.WriteH5(file) writes source metadata and frame records.
  Top-level WriteTrajectoryH5 combines.
- **Dense buffers owned by TrajectoryProtein, queried via TrajectoryResult
  methods.** For per-atom × N output, TrajectoryResult transfers
  ownership of a contiguous buffer at Finalize; query methods on the
  Result dereference into it. Avoids per-atom std::vector allocation for
  dense data.
- **NmrAtomIdentity as topology-like on Protein.** Typed identity
  lenses populated at Protein construction from topology + IUPAC
  convention. Zero-string after load. Orthogonal indices on Protein
  for fast per-lens slicing.

---

## 12 — Open questions and decisions to make during implementation

Preserved here so implementation doesn't silently close them.

1. **The frame-0 ConformationResults on tp.Protein().conformations_[0].**
   After Run completes, conf0 has its ConformationResults attached and
   queryable. Does the analysis H5 emit them as "conf0's data" alongside
   per-frame time series? Today's AnalysisWriter emits `/positions/xyz`
   at (T, N, 3) which implicitly includes conf0 as frame 0. That's fine.
   But conf0 also has e.g. an OrcaShieldingResult if the run includes
   DFT comparison — does that go in the trajectory H5 or stay on a
   sibling conformation H5? Design call during implementation.

2. **Per-frame NPY snapshots.** The current analysis mode writes per-frame
   NPY snapshots at 1 ns intervals via `ConformationResult::WriteAllFeatures`.
   In the new model, those snapshots are the ConformationResults on
   selected ProteinConformations — but those conformations die after
   their frame. To preserve the snapshot, we'd need a mechanism to
   either (a) write the NPYs during Compute of a SnapshotTrajectoryResult,
   or (b) keep a small ring of "interesting" conformations alive for
   later snapshot emission. Design call.

3. **Selection callback shape.** ~~Options: a selection-collector
   pointer on RunContext; a Selection facility on TrajectoryProtein;
   a callback on Trajectory. Implementation call; pattern will stabilise.~~
   **RESOLVED (revision 2026-04-22 after fresh-agent review): see
   Appendix D.** `SelectionEmittingTrajectoryResult` mixin, collected by
   `Trajectory::Run` at end of Finalize via `dynamic_cast`. Keeps
   `TrajectoryResult::Compute` signature clean; opt-in interface.

4. **Dependency validation for ConformationResult types.** ~~The design
   in §5 Phase 2 does this in Trajectory::Run after attach. Alternatively
   it could happen at AttachResult itself by passing the RunConfiguration
   to AttachResult. Design call; small.~~
   **RESOLVED (revision 2026-04-22): stays in §5 Phase 2.** Reason:
   attach doesn't need to know about run configuration; Phase 2 is the
   one place with both attached results and run config in scope; keeping
   attach signature clean matches ConformationResult's pattern exactly.

5. **Scope of MopacVsFf14SbReconciliationTrajectoryResult.** Per
   desiderata D5, this is a diagnostic that would fit under
   FullFatFrameExtraction. Design the diagnostic before or after the
   accumulator refactor lands. Deferred.

6. **Frame-selection mechanism for FullFatFrameExtraction (NEW,
   2026-04-22).** Two options noted in §6: (a) extract a frame-filtered
   XTC as an intermediate and run FullFatFrameExtraction on that, or
   (b) `RunContext::SetSelectedFrameIndices()` makes `GromacsFrameHandler`
   skip frames not in the set. Fleet currently uses pattern (a).
   Pattern (b) is cleaner but needs a small addition to
   `GromacsFrameHandler::Next()`. Decide during implementation;
   non-blocker.

---

## Notes for future sessions reading this document

- **This file supersedes sections of OBJECT_MODEL.md and PATTERNS.md
  pertaining to trajectory scope until it folds in.** Header notes in
  those files point here. Trust this file over those two for
  trajectory-scope guidance during the design window.

- **Do not treat this file as committed design.** It is the design
  document for a coordinated library update. Specific details
  (virtual method signatures, exact field names, exact dependency
  declarations) are subject to refinement during the implementation
  pass. The architectural pattern and entity model are the load-bearing
  parts.

- **Concrete > general.** Every claim here is grounded in an existing
  code pattern (`ConformationResult.h`, `ProteinConformation.cpp`,
  `OperationRunner.cpp`, `GromacsProtein.h/cpp`, `GromacsFrameHandler.cpp`,
  `BiotSavartResult.h/cpp`) or in a specific desiderata item from
  `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md` or
  `spec/PLANNED_CALCULATORS_2026-04-22.md`. Deviations from those
  patterns are called out and justified.

- **Appendices A–H close specific gaps** surfaced by cold-read
  reviews on 2026-04-22. Main body sections point at the relevant
  appendix where they previously hand-waved a detail. Reading order:
  main body for architecture, appendices for the specific
  inventories and decisions. Appendix H (implementation staging
  roadmap) sequences the 12–20 session work into seven
  checkpointable stages.

- **Amendments below, never rewrite above.** Same discipline as other
  WIP specs in the project.

---

## Appendix A — `IupacAtomPosition` enum inventory and generator

**⚠ STATUS: PROPOSAL PENDING USER REVIEW (see §2).** The *need* for
NMR-aware IUPAC-typed atom identity is confirmed by user. The
specific approach below — TOML-driven canonical dictionary, build-step
generator script, the particular enum / mapping shapes — is one
concrete design. It is not yet user-approved. Alternatives that
deserve discussion: hand-written enum with a smaller core set
(covering just the 20 residues without tautomer / cap expansion, with
an "other" escape hatch); a hybrid where the core set is typed and
edge cases fall through to string lookup with logging; a generator
driven by a source other than TOML (e.g., a Python module with the
dictionary in native Python, or a directly generated from a
reference like the AMBER prepi files). User sign-off on this
appendix is a prerequisite to implementation. Do not begin building
`tools/generate_iupac_enum.py` or the TOML dictionary until the
approach is confirmed.

---

### Scope

`IupacAtomPosition` is a typed enum with one value per valid
(AminoAcid residue type, IUPAC atom position) pair, including
protonation-state variants and terminal caps. It is the primary
per-atom identity key for external-data binding (BMRB, RefDB, PCS
datasets, published CSA tensors) — everything that names an atom in
the outside world ends up passing through this enum at the library
boundary.

### Source: PDB CCD (primary) or AMBER residue library (cross-check)

The TOML dictionary is not authored — it is generated at library
build time from a canonical reference source. Priority:

1. **PDB Chemical Component Dictionary (CCD)** via `cifpp` (already a
   project dependency) if `cifpp` exposes per-residue ChemComp
   templates with stereochemistry attributes; else via `gemmi` (new
   build-time dep only, not runtime). CCD covers every standard
   residue plus every variant the PDB has seen — including all 20
   standard amino acids used in mutant calibration.
2. **AMBER ff14SB residue library** (`.prepi` / `.lib` / `.off`
   files) as a cross-check source. For standard residues this
   redundantly covers CCD; useful when the library's force field
   representation needs to match the residue template exactly.
3. **IUPAC rules** (Markley 1998 + subsequent recommendations)
   consulted for stereochemistry assignments (pro-R / pro-S slots
   in methyl / methylene groups). Encoded in the template via
   stereo-descriptor fields rather than derived at generation time.

No hand-authored residue templates. Every template entry traces to
CCD or AMBER; any deviation from those sources requires an explicit
justification in the TOML (a per-entry comment pointing at the
reason — typically a protonation variant or a cap that one source
doesn't express).

### Why generate, not hand-write

The enum is large (several hundred values), every value is
significant, and the "zero runtime strings past load" rule means a
single uncovered atom silently breaks the rule in production. Three
concrete classes of drift if hand-written:

1. New variant added to the protonation-state dictionary (e.g., a
   special ASN tautomer for a reviewer's paper) — easy to forget the
   enum addition, easy for a caller's atom to fall through.
2. BMRB / RefDB aliases for the same position (e.g., HA2 vs 1HA for
   glycine alpha Hs) — if the mapping table isn't kept in sync with
   the enum, external-data binding silently uses the wrong slot.
3. Typos in hand-written enum values carrying CHARMM naming
   conventions where we meant IUPAC.

Generation eliminates all three.

### The canonical dictionary

A TOML file (working name `data/iupac_atom_dictionary.toml`) is the
reviewable source of truth. One top-level table per residue type,
including variants:

```toml
[residues.ALA]
iupac_name = "ALA"
variants = ["ALA"]  # Only one variant (no protonation state ambiguity).

[residues.ALA.atoms]
N = { role = "BackboneN", locant = "backbone" }
H = { role = "BackboneHN", locant = "backbone",
      charmm_alias = "HN", bmrb_alias = "H" }
CA = { role = "BackboneCA", locant = "backbone" }
HA = { role = "BackboneHA", locant = "backbone" }
C = { role = "BackboneC", locant = "backbone" }
O = { role = "BackboneO", locant = "backbone" }
CB = { role = "SidechainC_sp3", locant = "beta" }
HB1 = { role = "MethylCH3", locant = "beta",
        methyl_group = "ALA_beta", methyl_slot = 1 }
HB2 = { role = "MethylCH3", locant = "beta",
        methyl_group = "ALA_beta", methyl_slot = 2 }
HB3 = { role = "MethylCH3", locant = "beta",
        methyl_group = "ALA_beta", methyl_slot = 3 }

[residues.HIS]
# HIS expands to three variants by default.
iupac_name = "HIS"
variants = ["HID", "HIE", "HIP"]

[residues.HIS.variants.HID]
# delta-protonated: ND1 carries H, NE2 does not
extra_atoms = { HD1 = { role = "RingNH", locant = "delta" } }
absent_atoms = ["HE2"]
ring_atom_role = { ND1 = "Ndelta1", CE1 = "Cepsilon1",
                   NE2 = "Nepsilon2", CD2 = "Cdelta2",
                   CG = "Cgamma" }

[residues.HIS.variants.HIE]
extra_atoms = { HE2 = { role = "RingNH", locant = "epsilon" } }
absent_atoms = ["HD1"]
ring_atom_role = { ND1 = "Ndelta1", CE1 = "Cepsilon1",
                   NE2 = "Nepsilon2", CD2 = "Cdelta2",
                   CG = "Cgamma" }

[residues.HIS.variants.HIP]
# Doubly-protonated: both ND1 and NE2 carry H
extra_atoms = { HD1 = { role = "RingNH", locant = "delta" },
                HE2 = { role = "RingNH", locant = "epsilon" } }
absent_atoms = []
ring_atom_role = { ND1 = "Ndelta1", CE1 = "Cepsilon1",
                   NE2 = "Nepsilon2", CD2 = "Cdelta2",
                   CG = "Cgamma" }

[residues.CYS]
iupac_name = "CYS"
variants = ["CYS", "CYX"]
# CYS = free thiol (has SH); CYX = disulfide-bonded (no SH)

[residues.CYS.variants.CYS]
extra_atoms = { HG = { role = "ThiolSH", locant = "gamma" } }

[residues.CYS.variants.CYX]
absent_atoms = ["HG"]

# ... (17 more standard residues)

[caps.ACE]
# Acetyl N-terminal cap.
atoms = {
    CH3 = { role = "MethylCH3", locant = "cap" },
    HH31 = { ... }, HH32 = { ... }, HH33 = { ... },
    C = { role = "BackboneC", locant = "cap" },
    O = { role = "BackboneO", locant = "cap" },
}

[caps.NME]
# N-methyl amide C-terminal cap.
atoms = { ... }

[caps.NTerminus]
# Add to any residue at N-terminus: extra H1/H2/H3 (ammonium) or
# cap residue ACE instead.
extra_atoms_if_ammonium = {
    H1 = { role = "BackboneHN", locant = "backbone", terminus = "N" },
    H2 = { role = "BackboneHN", locant = "backbone", terminus = "N" },
    H3 = { role = "BackboneHN", locant = "backbone", terminus = "N" },
}

[caps.CTerminus]
extra_atoms_if_carboxyl = {
    OXT = { role = "CarboxylO", locant = "backbone", terminus = "C" },
}
```

Variant handling: each `(residue, variant)` pair generates its own set
of `IupacAtomPosition` enum values. For ALA, which has one variant,
the values are `ALA_N`, `ALA_H`, `ALA_CA`, `ALA_HA`, ..., `ALA_HB3`.
For HIS, values are generated per tautomer: `HID_N`, `HID_H`, ...,
`HID_HD1` (but no `HID_HE2`); `HIE_N`, `HIE_H`, ..., `HIE_HE2` (but
no `HIE_HD1`); `HIP_N`, ..., `HIP_HD1`, `HIP_HE2`.

### Generator script

A small Python script (working name `tools/generate_iupac_enum.py`,
not yet written) reads the TOML at build time and emits:

- `src/generated/IupacAtomPosition.h` — the `enum class
  IupacAtomPosition` with all ~400 values, alphabetised for
  diff-stability.
- `src/generated/IupacAtomPositionMappings.cpp` — lookup tables from
  `(AminoAcid variant, pdb name)` to `IupacAtomPosition`, from CHARMM
  alias names, from BMRB alias names. Also `IupacAtomPosition` →
  `NmrClass`, `Locant`, `MethylGroup` / slot, etc.
- `src/generated/IupacAtomPositionCount.h` — integer constant with
  the count, for array sizing if needed.

Build integration: the generator runs at CMake configure time if the
TOML changes; emits files into `build/generated/`; those are added to
the build target. Editor tooling (clangd) sees them via
`compile_commands.json`.

### Coverage checklist

Before first generation, the dictionary must cover:

- [ ] All 20 standard amino acids with default protonation (pH 7)
- [ ] HIS variants: HID, HIE, HIP
- [ ] CYS variants: CYS (free), CYX (disulfide-bonded)
- [ ] Caps: ACE, NME, N-terminus ammonium (H1/H2/H3), C-terminus OXT
- [ ] Protonation variants: ASH (protonated ASP), GLH (protonated
      GLU), LYN (neutral LYS) — optional, but include for completeness
- [ ] Deuterated variants: **deferred**. If deuterium appears, treat
      as H for now; re-visit if isotope labelling becomes first-class.
- [ ] Non-standard residues (seleno-Met, phosphoserine, etc.):
      **deferred**. Fleet is standard-20 only; re-visit when expanded.

A coverage test in the test suite walks every `(residue, variant)`
combination, queries the generated mapping, and asserts every atom in
a canonical reference PDB of that residue maps to a defined
`IupacAtomPosition`. Fails loudly on uncovered atoms.

### Interaction with `NmrAtomIdentityBuilder`

`NmrAtomIdentityBuilder` consumes the generated mappings. Given a
loaded Protein's atom list (element, residue type, variant, pdb name),
the builder:

1. Looks up `IupacAtomPosition` from `(variant, pdb_name)`.
2. Looks up `NmrClass`, `Locant`, `MethylGroup`, `ChiParticipation`,
   `RingAtomRole`, `ring_membership` from the generated tables.
3. Populates `protein.atom_identities_[i]` with the full
   `NmrAtomIdentity` struct.

No string matching, no hardcoded tables in code — all data lives in
the TOML.

---

## Appendix B — `TrajectoryBond`: per-bond trajectory store (rejected alternative)

**STATUS: DECIDED AGAINST 2026-04-22 — option (b) chosen.** See §3's
"Per-bond trajectory state" subsection for the decision and rationale.
Per-bond accumulator state lives internal to bond-scope
TrajectoryResults, not in a first-class `TrajectoryBond` store on
TrajectoryProtein. This appendix is preserved as the rejected option
(a) design — the shape that was considered and declined — so a future
session considering this design space sees the path that was
considered and why it was passed over. Do not implement this
appendix.

---

### The role

Parallel to `TrajectoryAtom`, scoped to bonds. `TrajectoryProtein`
holds `std::vector<TrajectoryBond>` of length `protein.BondCount()`,
indexed the same as `protein.Bonds()`. Per-bond fields are declared
upfront (same discipline), written by specific TrajectoryResults,
singleton-per-writer.

### Structure

```cpp
class TrajectoryBond {
    friend class TrajectoryProtein;
public:
    // Identity through back-pointer, same as TrajectoryAtom.
    // Element, atom_a / atom_b indices, order, category → access via
    // tp.Protein().BondAt(bond_idx).

    // === Written by BondLengthWelfordTrajectoryResult ===
    double length_mean = 0.0;    // Å
    double length_m2 = 0.0;      // Welford sum-of-sq-dev; std at Finalize
    double length_min = std::numeric_limits<double>::infinity();
    double length_max = -std::numeric_limits<double>::infinity();
    size_t length_n_frames = 0;
    double length_std = 0.0;     // Finalize-only

    // === Written by BondLengthDeltaTrajectoryResult ===
    // Frame-to-frame bond-length fluctuation rate.
    double length_delta_mean = 0.0;
    double length_delta_m2 = 0.0;
    double length_delta_std = 0.0;  // Finalize-only

    // === Written by BondOrderWelfordTrajectoryResult (FullFat only) ===
    // MOPAC Wiberg bond orders, accumulated only when MOPAC runs.
    // Populated sparse — only frames where MOPAC computed a bond order.
    double wiberg_mean = 0.0;
    double wiberg_m2 = 0.0;
    double wiberg_std = 0.0;
    size_t wiberg_n_frames = 0;

    // (additional per-bond fields as TrajectoryResults are added)

private:
    explicit TrajectoryBond() = default;
};
```

### Accessors on TrajectoryProtein

```cpp
size_t BondCount() const;  // == protein().BondCount()
const TrajectoryBond& BondAt(size_t bond_idx) const;
TrajectoryBond& MutableBondAt(size_t bond_idx);
const std::vector<TrajectoryBond>& Bonds() const;
```

### Writers

`BondLengthWelfordTrajectoryResult` is the canonical example. Its
`Compute` iterates `tp.Protein().Bonds()`, computes current-frame
length from `conf.AtomAt(atom_a).Position()` and `conf.AtomAt(atom_b).Position()`,
and ticks the Welford on `tp.MutableBondAt(bond_idx)`. Same pattern as
`BsWelfordTrajectoryResult` (§4), just indexed by bond instead of
atom.

### H5 emission

Natural H5 group: `/bonds/length_mean`, `/bonds/length_std`, etc.
Each emitting TrajectoryResult owns its own group within `/bonds/`,
written via `WriteH5Group`. Metadata (bond atom indices, bond order,
category) written once in `/bonds/topology/` by a small
`BondTopologyTrajectoryResult` that attaches unconditionally (always
in every RunConfiguration's factory list).

### Note

No `TrajectoryRing` at this pass. Per-ring trajectory fields currently
attach to per-atom storage (e.g., `ring_neighbours` on ConformationAtom
becomes `ring_neighbour_stats` on TrajectoryAtom). If per-ring
trajectory quantities become first-class later (e.g., ring-centre
time series, ring-normal autocorrelation), a sibling `TrajectoryRing`
store follows the same pattern. Deferred.

---

## Appendix C — Dense buffer layouts for non-scalar payloads

### Template convention

`DenseBuffer<T>` is instantiated with the native C++ type (`Vec3`,
`Mat3`, `SphericalTensor`) directly — these are memory-contiguous
structs whose `sizeof(T)` corresponds to their natural serialisation
size. No intermediate flattening to double arrays in memory.

```cpp
template <typename T>
class DenseBuffer {
public:
    DenseBuffer(size_t atom_count, size_t stride_per_atom);

    // Typed element access: atom_idx, offset within the stride.
    T& At(size_t atom_idx, size_t offset);
    const T& At(size_t atom_idx, size_t offset) const;

    // Typed atom slice: span of length stride_per_atom for one atom.
    std::span<T> AtomSlice(size_t atom_idx);
    std::span<const T> AtomSlice(size_t atom_idx) const;

    // Raw storage access for H5 emission. Layout is atom-major,
    // contiguous per-atom.
    const T* RawData() const { return storage_.data(); }
    size_t TotalElementCount() const { return storage_.size(); }

    size_t AtomCount() const;
    size_t StridePerAtom() const;

private:
    std::vector<T> storage_;  // [atom_0_stride, atom_1_stride, ...]
    size_t atom_count_;
    size_t stride_per_atom_;
};
```

### H5 emission conventions by payload type

**`DenseBuffer<double>`** (e.g., per-atom per-frame T0 shielding scalar
time series). H5 dataset shape `(atom_count, stride)`. Row-major.
Metadata: `"units"` attribute (e.g., `"ppm"`), `"stride_meaning"`
attribute (e.g., `"frame"`, `"lag"`, `"frequency"`).

**`DenseBuffer<Vec3>`** (e.g., per-atom per-frame position). H5 dataset
shape `(atom_count, stride, 3)`. Emission flattens the `Vec3`'s three
doubles to the trailing dimension. Metadata: `"units"` (e.g., `"Å"`),
`"components"` attribute (e.g., `"[x,y,z]"`).

**`DenseBuffer<Mat3>`** (e.g., per-atom per-frame EFG tensor). H5
dataset shape `(atom_count, stride, 3, 3)`. Trailing dimensions are
the row/column indices. Metadata: `"units"`, `"layout"` (e.g.,
`"row_major"`).

**`DenseBuffer<SphericalTensor>`** (e.g., per-atom per-frame shielding
tensor irrep decomposition). H5 dataset shape `(atom_count, stride, 9)`.
Trailing dimension is `[T0, T1[-1], T1[0], T1[1], T2[-2], T2[-1], T2[0],
T2[1], T2[2]]` — matches the existing SphericalTensor memory layout.
Metadata: `"irrep_layout"` attribute (the string above verbatim),
`"units"`, `"sign_convention"` (e.g., `"NMR shielding positive"`).

### Why the native-type template works

The memory layout of each of these types is fixed and contiguous in
C++:

- `Vec3` (Eigen `Matrix<double, 3, 1>`): 3 doubles, 24 bytes, no
  padding.
- `Mat3` (Eigen `Matrix<double, 3, 3>`): 9 doubles, 72 bytes, no
  padding, column-major by default (Eigen's convention — note: H5
  emission needs to respect or transpose).
- `SphericalTensor`: `{double T0; array<double, 3> T1; array<double, 5>
  T2;}` → 9 contiguous doubles. `sizeof(SphericalTensor) == 72`.

For H5 emission, cast `RawData()` (a `T*`) to `const double*` and emit
as `(atom_count × stride × elements_per_T)` shape. Elements-per-T = 3
for Vec3, 9 for Mat3 and SphericalTensor. Emission has no memcpy,
no reinterpretation — just declared shape in the H5 dataset
definition.

### Single caveat: Eigen Mat3 column-major

Eigen defaults to column-major storage. If the caller expects
row-major in H5, transpose at emission time (or declare the H5
dataset as column-major via the HighFive C-vs-F convention). The
choice is library-wide; default to column-major in H5 to avoid
unnecessary transposes, and document the convention in the schema.

---

## Appendix D — Selection callback pattern

### The problem

In `ScanForDftPointSet` mode, specific TrajectoryResults detect
events during per-frame Compute — rotamer transitions, RMSD spikes,
χ₁ bin crossings, ring-flip candidates — and need to record them for
downstream DFT submission. The records live on `Trajectory`'s
`selections_` vector as `FrameSelectionRecord`s. The question: how
does data flow from a TrajectoryResult's Compute back up to
Trajectory.

### The chosen pattern: opt-in interface

A separate interface (`SelectionEmittingTrajectoryResult`) implemented
by TrajectoryResults that emit selections. Trajectory collects from
them at end of Run via `dynamic_cast`. TrajectoryResults that don't
emit selections implement only the base `TrajectoryResult` interface
and have no awareness of selection machinery.

```cpp
struct FrameSelectionRecord {
    size_t frame_idx;
    double time_ps;
    std::string reason;            // e.g., "rotamer_transition_chi1_LYS_42"
    std::string selector;          // which TrajectoryResult flagged it
    // Free-form metadata for downstream consumers:
    std::map<std::string, std::string> metadata;
};

// Mixin interface.
class SelectionEmittingTrajectoryResult {
public:
    virtual ~SelectionEmittingTrajectoryResult() = default;
    virtual std::span<const FrameSelectionRecord> SelectionRecords() const = 0;
};

// Example: multiple-inheritance at the subclass.
class ChiRotamerSelectionTrajectoryResult
    : public TrajectoryResult,
      public SelectionEmittingTrajectoryResult {
public:
    std::string Name() const override { return "ChiRotamerSelectionTrajectoryResult"; }
    std::vector<std::type_index> Dependencies() const override { ... }
    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 size_t frame_idx, double time_ps) override {
        // Detect chi1 bin crossings by reading conf.AtomAt(...) chi
        // values and comparing to this result's internal "last frame's
        // bins" state. When a transition is detected, push a record.
        // ... detection logic ...
        if (transition_detected) {
            selections_.push_back({frame_idx, time_ps,
                "chi1_rotamer_transition_" + residue_name,
                Name(), metadata});
        }
    }
    void Finalize(TrajectoryProtein& tp) override { /* no-op or cleanup */ }

    // From SelectionEmittingTrajectoryResult:
    std::span<const FrameSelectionRecord> SelectionRecords() const override {
        return selections_;
    }

private:
    std::vector<FrameSelectionRecord> selections_;
    // + internal state for detection
};
```

### Collection in Trajectory::Run

At end of Finalize phase:

```cpp
// Phase 6: collect selections.
for (auto* result : tp.ResultsInAttachOrder()) {
    if (auto* emitter =
            dynamic_cast<SelectionEmittingTrajectoryResult*>(result)) {
        auto records = emitter->SelectionRecords();
        selections_.insert(selections_.end(),
                           records.begin(), records.end());
    }
}
```

One pass over attached results; `dynamic_cast` on each to check for
the interface; merge records into `Trajectory::selections_`. Cheap,
explicit, no virtual dispatch added to the hot-path `Compute` call.

### Why not a callback parameter on Compute

Three reasons:

1. Most TrajectoryResults don't emit selections — adding a parameter
   they must accept but will never use is signature pollution.
2. The base `TrajectoryResult` class should not depend on
   `FrameSelectionRecord`, `Trajectory`, or any process-scope type.
   The type hierarchy for physical-world results should not reach
   across into process-scope entities at the base-class level.
3. `dynamic_cast` in the collection loop is the explicit "who knows
   about selections" check. It surfaces at a single point in the
   code — easy to audit, easy to find with grep, not scattered across
   every TrajectoryResult subclass.

### Why not a facility on TrajectoryProtein

Would couple the physical-world container (TrajectoryProtein) to a
process concept (selection). TrajectoryProtein is the protein as
observed; selections are bookkeeping about the run. They belong on
Trajectory.

---

## Appendix E — EnrichmentResult + OpenBabel disposition

**STATUS: DECIDED 2026-04-22 — NO MIGRATION. OpenBabel and
EnrichmentResult stay unchanged.** User preference verbatim: *"I am
happy to pay for openbabel per frame in exchange for not mucking with
everything on that level."* The migration plan that previously lived
here (moving hybridisation inference from per-frame EnrichmentResult
to one-shot Protein construction, EnrichmentResult becomes a
validator + shim, eventual removal of ConformationAtom hybridisation
field, etc.) is explicitly declined.

### What this means

`src/EnrichmentResult.{h,cpp}` continues to run per-frame in the
OperationRunner sequence, exactly as it does today:
- Invokes OpenBabel for hybridisation inference.
- Populates `ConformationAtom::hybridisation`.
- Sets all existing categorical booleans on ConformationAtom
  (`is_backbone`, `is_amide_H`, `is_alpha_H`, `is_methyl`,
  `is_aromatic_H`, `is_on_aromatic_residue`, `is_hbond_donor`,
  `is_hbond_acceptor`, `parent_is_sp2`).

No shim, no validator mode, no retirement path. The per-frame cost
(milliseconds per frame, absorbed in the existing calculator pipeline
timing) is accepted in exchange for zero disruption to any existing
calculator that reads from ConformationAtom.

### Relationship to NmrAtomIdentity (Appendix A + §2)

NmrAtomIdentity is an **additive layer on Protein**, independent of
EnrichmentResult:

- `NmrAtomIdentityBuilder` runs once at `Protein::FinalizeConstruction`,
  reading residue type, IUPAC atom name, bond graph, and ring detection
  output (all topology-fixed).
- Does NOT invoke OpenBabel. Hybridisation is not part of NmrAtomIdentity.
- Populates `protein.atom_identities_[i]` with the typed identity
  fields (IupacAtomPosition, NmrClass, MethylGroup, RingAtomRole, etc.).

New consumers (trajectory-scope TrajectoryResults, Python SDK
typed-slicing queries) use `protein.AtomIdentity(i).nmr_class`,
`protein.AtomsByNmrClass(...)`, etc. Existing consumers that read
`conf.atoms[i].is_amide_H` continue to do so; nothing changes for
them. Redundancy accepted; consistency by construction (both sources
derive from the same topology, so they don't disagree).

### Why this shape

- **No cross-cutting changes during the activation window.** Every
  calculator that reads from ConformationAtom's booleans / hybridisation
  keeps reading unchanged. Static-path calculators (non-regression
  contract in §0) are therefore trivially preserved — there's no
  migration for them to miss.
- **The "redundancy" cost is tiny.** Per-frame OpenBabel is
  milliseconds; the booleans are microseconds. Absorbed by the
  per-frame calculator work that dominates wall time (APBS, AIMNet2,
  ring current, McConnell all individually cost more).
- **The rollup scope narrows.** `EnrichmentResult`, ConformationAtom's
  booleans, hybridisation-storage disposition — all unchanged. The
  rollup touches fewer files, surfaces fewer migration risks, adds
  less validation burden.

### What Appendix E had and no longer does

The migration plan previously documented here is preserved in git
history if ever wanted. It is not the current plan. Future sessions
reading this appendix should take this decision (not the previous
plan) as authoritative.

---

## Appendix F — TrajectoryResult catalog

Concrete classes referenced by `RunConfiguration` factories in §6.
Each row is approximately one class to write in implementation.

Column conventions:
- **Source**: the ConformationResult type(s) that must run per frame
  for this TrajectoryResult to read valid data. Enforced by the
  TrajectoryResult's `Dependencies()` return value.
- **Lifecycle**: `AV` = always-valid mid-stream; reads during
  streaming return current partial values. `FO` = finalize-only;
  reads before Finalize are undefined.
- **Emission**: where output lives. `TA` = TrajectoryAtom fields.
  `DB` = dense buffer owned by TrajectoryProtein post-Finalize.
  `RI` = result-internal (bond-scope Results; see §3 option-(b)
  decision + Appendix B rejected alternative). May be multiple.

### `PerFrameExtractionSet` (the default fleet extraction)

| Class | Source | Lifecycle | Emission | Brief |
|---|---|---|---|---|
| `PositionsTimeSeriesTrajectoryResult` | (none — reads positions) | FO | DB | Per-atom per-frame `Vec3` positions. Large. |
| `BsWelfordTrajectoryResult` | BiotSavartResult | AV | TA | Per-atom Welford of `bs_shielding_contribution.T0` + per-type T0 sums. |
| `BsShieldingTimeSeriesTrajectoryResult` | BiotSavartResult | FO | DB | Per-atom per-frame `SphericalTensor` for BS shielding. |
| `HmWelfordTrajectoryResult` | HaighMallionResult | AV | TA | Per-atom Welford of HM T0 + per-type T0 sums. |
| `HmShieldingTimeSeriesTrajectoryResult` | HaighMallionResult | FO | DB | Per-atom per-frame SphericalTensor for HM. |
| `McConnellWelfordTrajectoryResult` | McConnellResult | AV | TA | Per-atom Welford of mc_shielding T0 + per-category (CO/CN/sidechain/aromatic) sums. |
| `McConnellShieldingTimeSeriesTrajectoryResult` | McConnellResult | FO | DB | Per-atom per-frame SphericalTensor. |
| `CoulombFieldTimeSeriesTrajectoryResult` | CoulombResult | FO | DB | Per-atom per-frame `coulomb_E_total` (Vec3) + `coulomb_EFG_total` (Mat3) + decomposed (backbone/sidechain/aromatic). Large. |
| `CoulombFieldWelfordTrajectoryResult` | CoulombResult | AV | TA | Per-atom Welford of E-field magnitude + EFG Frobenius norm. |
| `ApbsFieldTimeSeriesTrajectoryResult` | ApbsFieldResult | FO | DB | Per-atom per-frame APBS E-field (Vec3) + EFG (Mat3). |
| `WaterEnvironmentTimeSeriesTrajectoryResult` | WaterFieldResult | FO | DB | Per-atom per-frame water E-field, EFG, shell counts. |
| `WaterEnvironmentWelfordTrajectoryResult` | WaterFieldResult | AV | TA | Per-atom Welford of water E-field magnitude + first-shell count. |
| `HydrationShellTimeSeriesTrajectoryResult` | HydrationShellResult | FO | DB | Per-atom half_shell_asymmetry, mean_water_dipole_cos, nearest_ion_distance time series. |
| `HydrationGeometryTimeSeriesTrajectoryResult` | HydrationGeometryResult | FO | DB | Per-atom SASA-normal water dipole + alignment + coherence per frame. |
| `AIMNet2ChargeTimeSeriesTrajectoryResult` | AIMNet2Result | FO | DB | Per-atom per-frame AIMNet2 Hirshfeld charge. Small. |
| `AIMNet2EmbeddingTimeSeriesTrajectoryResult` | AIMNet2Result | FO | DB | Per-atom per-frame 256-dim AIMNet2 embedding. Very large. Optional. |
| `EeqChargeWelfordTrajectoryResult` | EeqResult | AV | TA | Per-atom Welford of EEQ charge + coordination number. |
| `SasaTimeSeriesTrajectoryResult` | SasaResult | FO | DB | Per-atom SASA time series + SASA normal time series. |
| `SasaWelfordTrajectoryResult` | SasaResult | AV | TA | Per-atom Welford of SASA, delta-SASA. |
| `HBondTimeSeriesTrajectoryResult` | HBondResult | FO | DB | Per-atom HBond nearest-distance + count-within-3.5Å time series. |
| `HBondCountWelfordTrajectoryResult` | HBondResult | AV | TA | Per-atom Welford of HBond count. |
| `DihedralTimeSeriesTrajectoryResult` | DsspResult | FO | DB | Per-residue per-frame φ, ψ, χ1-4 (cos, sin, exists). |
| `DihedralBinTransitionTrajectoryResult` | DsspResult | AV | TA | Per-residue running transition counts (Ramachandran + χ rotamer). |
| `Dssp8TimeSeriesTrajectoryResult` | DsspResult | FO | DB | Per-residue per-frame 8-class SS code + H-bond energies. |
| `Dssp8TransitionTrajectoryResult` | DsspResult | AV | TA | Per-residue SS transition counts. |
| `BondedEnergyTimeSeriesTrajectoryResult` | BondedEnergyResult | FO | DB | Per-atom per-frame energy decomposition (bond/angle/UB/proper/improper/CMAP). |
| `GromacsEnergyTimeSeriesTrajectoryResult` | GromacsEnergyResult | FO | DB | Per-frame aggregate energy terms (scalar per frame, not per atom). |
| `PiQuadrupoleShieldingTimeSeriesTrajectoryResult` | PiQuadrupoleResult | FO | DB | Per-atom per-frame PQ shielding contribution. |
| `RingSusceptibilityShieldingTimeSeriesTrajectoryResult` | RingSusceptibilityResult | FO | DB | Per-atom per-frame RingChi shielding contribution. |
| `DispersionShieldingTimeSeriesTrajectoryResult` | DispersionResult | FO | DB | Per-atom per-frame dispersion shielding contribution. |
| `RingNeighbourhoodTrajectoryStats` | BiotSavartResult + HaighMallionResult + PiQuadrupoleResult + RingSusceptibilityResult + DispersionResult | FO | TA (`vector<...>` per atom) | Rich per-atom-per-ring stats accumulated across trajectory. |
| `BondTopologyTrajectoryResult` | (none — topology) | AV | TA (writes to `/bonds/topology/`) | Bond atom indices, order, category. Stable across frames; just needs one write. |
| `BondLengthWelfordTrajectoryResult` | (none — reads positions) | AV | RI | Per-bond Welford of bond length. |
| `BondLengthDeltaTrajectoryResult` | (none — reads positions) | AV | RI | Per-bond Welford of frame-to-frame length change. |

~30 classes in the full set. Most follow the `BsWelfordTrajectoryResult`
(AV) or `BsT0AutocorrelationTrajectoryResult` (FO + DB) templates.

### `ScanForDftPointSet` (scan mode)

Strict subset of PerFrameExtractionSet with additional selection-
emitting classes:

| Class | Source | Lifecycle | Emission | Brief |
|---|---|---|---|---|
| `DihedralBinTransitionTrajectoryResult` | DsspResult | AV | TA | (same as in PerFrameExtractionSet; shared implementation) |
| `Dssp8TransitionTrajectoryResult` | DsspResult | AV | TA | (shared) |
| `RmsdTrackingTrajectoryResult` | (none — reads positions) | AV | TA (per-frame vector on tp, not per-atom) | Running RMSD of protein vs reference (e.g., frame 0). |
| `SasaWelfordTrajectoryResult` | SasaResult | AV | TA | (shared) |
| **`ChiRotamerSelectionTrajectoryResult`** | DsspResult | AV | TA + emits selections | Implements `SelectionEmittingTrajectoryResult`. Detects χ₁/χ₂ transitions; records FrameSelectionRecord. |
| **`RmsdSpikeSelectionTrajectoryResult`** | (RmsdTrackingTrajectoryResult) | AV | emits selections | Dependency on RmsdTracking's output. Detects RMSD > threshold; records. |
| **`DftPoseCoordinatorTrajectoryResult`** | (all selection emitters) | FO | emits selections | At Finalize: deduplicates, budgets, selects final DFT pose set. Reads other selection emitters' SelectionRecords at Finalize. |

The bold ones implement `SelectionEmittingTrajectoryResult`; the rest
are shared with PerFrameExtractionSet where the behaviour matches.

### `FullFatFrameExtraction` (sparse full-fat mode)

PerFrameExtractionSet plus MOPAC-family additions:

| Class | Source | Lifecycle | Emission | Brief |
|---|---|---|---|---|
| (inherits PerFrameExtractionSet set) | | | | |
| `MopacChargeWelfordTrajectoryResult` | MopacResult | AV (sparse) | TA | Per-atom Welford of MOPAC Mulliken charge across frames where MOPAC ran. |
| `MopacBondOrderWelfordTrajectoryResult` | MopacResult | AV (sparse) | RI | Per-bond Welford of MOPAC Wiberg bond order. |
| `MopacCoulombShieldingTimeSeriesTrajectoryResult` | MopacCoulombResult | FO (sparse) | DB | Per-atom per-(sparse)-frame MOPAC Coulomb shielding contribution. |
| `MopacMcConnellShieldingTimeSeriesTrajectoryResult` | MopacMcConnellResult | FO (sparse) | DB | Per-atom per-(sparse)-frame MOPAC McConnell shielding. |
| `MopacVsFf14SbReconciliationTrajectoryResult` | MopacCoulombResult + CoulombResult (+ OrcaShieldingResult for fleet-calibration frames) | FO | TA | Per-atom |cos| between (MopacCoulomb T2 - Coulomb T2) and (DFT delta T2 - Coulomb T2). Computed at Finalize from the sparse-frame data of the two MOPAC TrajectoryResults. |

"Sparse" means populated only on frames where MOPAC ran; off-frame
values are not present.

### Using the catalog

This table is the implementation checklist. Each row is approximately
one C++ class to write. The two worked examples in §4
(`BsWelfordTrajectoryResult` — AV scalar Welford; `BsT0AutocorrelationTrajectoryResult`
— FO dense-buffer) cover the two primary templates; all other rows
clone one of those two templates with source / emission / lifecycle
per the table.

---

## Appendix G — Scope and timeline framing (preserved verbatim)

User framing, 2026-04-22 afternoon:

> "This is maybe 12–20 sessions, a sixty hour week, and it is well
> worth it. The constraints on the rollup functionality are real in
> terms of dependencies — our working on this hard a week is nothing
> compared to all the stuff that has to run."

This sets the design-level investment as deliberately up-front. The
fleet extraction on 685 proteins (~25 min/protein = ~285 hours of
compute), followed by calibration, followed by Stage 2 / Stage 3 model
work, followed by the μs harvester phase on rented 5090s (41-day
window), followed by thesis writing — all of that downstream work
depends on the library being correct at activation time. A library-
level incorrectness discovered post-fleet costs weeks of re-extraction
plus schedule disruption to every downstream phase. A library-level
correctness captured pre-fleet compounds through all of it.

**Sequencing consequence:** this rollup does not decompose to hit a
tighter activation window. The fleet timing adapts to the rollup's
readiness, not the other way around. The session-count estimate
(12–20) drives the scheduling, not the other way.

**Quality consequence:** the appendices above (A–H) are not optional.
Test strategy, enum inventory, TrajectoryBond design, catalog,
implementation staging — these are what the 12–20 sessions produce.
If implementation happens without them, the rollup ships without the
correctness guarantees the investment was meant to buy.

---

## Appendix H — Implementation staging roadmap

**Revised 2026-04-22 pm after user direction:** original 7-stage plan
consolidated into **two primary stages** (trajectory refactor as one
push + activation), with the `NmrAtomIdentity` pre-bundle (original
Stage A) **deferred to a post-rollup effort**. Rationale in two parts:

1. **Shim drift is a known failure mode for AI-assisted refactors.**
   A multi-stage plan with an intermediate shim state (old and new
   pattern coexisting) risks the situation where a future session —
   working from conventional C++ patterns and partial context —
   settles into the shim state as a local optimum, rationalising it
   as "it works, don't touch it." The critical replacement work (kill
   GromacsProtein, kill GromacsProteinAtom, kill AccumulateFrame, kill
   AllWelfords) never happens. The rollup nominally completed with
   the worst of both patterns coexisting indefinitely.

2. **Adding typed names to atoms before cleaning up the existing
   atom-level mess is out of order.** `NmrAtomIdentity` is genuinely
   needed (§2 motivation holds) but adding typed identity fields on
   top of the current `GromacsProteinAtom` mess means the atom-level
   data store grows new responsibilities before the old responsibilities
   are cleaned up. Cleanup first; additive typed identity second.

The consolidated plan below keeps the seven original stages' content
as sub-checkpoints within a single push — the grouping is tighter but
no content is lost.

---

### 🔥 Stage 1 — Trajectory-scope refactor (single focused push)

**Scope:** the whole trajectory-scope replacement in one push, from
base infrastructure through rename through AnalysisWriter dissolution.
No intermediate shim state that exits the push commit. Original
Stages B + C + D + E + F as sub-checkpoints.

**Design refs:** §3 (anti-patterns + TrajectoryProtein), §4
(TrajectoryResult), §5 (Trajectory), §6 (RunConfiguration + RunContext),
§7 (Writer semantics + §7 structure-ID emission + library↔SDK
boundary), §9 (migration + renames), Appendix F (catalog).

**Proposal-pending items gating the push:** §3 TrajectoryBond design
(first-class store vs `BondStatsTrajectoryResult` internal). Must be
resolved before the push begins. `NmrAtomIdentity`-typed slicing is
NOT used by TrajectoryResults in this push (Stage A deferred);
TrajectoryResults iterate atoms by index, as ConformationResults do
today.

**Session budget:** 2–4 focused claude-code sessions, depending on
the implementing session's batching discipline in the per-field
migration sub-checkpoint.

#### Sub-checkpoint 1a — TrajectoryResult base + TrajectoryAtom + TrajectoryProtein (new infrastructure, zero removal yet)

Files:
- New: `src/TrajectoryResult.h/.cpp` — base class with virtuals per §4.
- New: `src/TrajectoryAtom.h` — declared-upfront field bundle with
  the §3-anti-patterns file-header comment block. Starts empty
  (field additions accumulate in 1b); no Welford/DeltaTracker
  instances ever.
- New: `src/TrajectoryProtein.h/.cpp` — wraps Protein; holds `atoms_`
  vector of TrajectoryAtom; attach discipline mirroring
  `ProteinConformation::AttachResult`; `ResultsInAttachOrder()`;
  `AdoptDenseBuffer<T>` template.
- New: `src/DenseBuffer.h` — template per Appendix C; supports T =
  double, Vec3, Mat3, SphericalTensor.
- Tests: attach discipline (singleton + dependency check), dense-
  buffer adoption + ownership transfer, basic Compute dispatch.

Checkpoint: library compiles. Existing tests unchanged (no existing
code touched yet). New infrastructure tested in isolation.

#### Sub-checkpoint 1b — Migrate every GromacsProteinAtom field into a TrajectoryResult, in batches

The grind. Per Appendix F catalog. Batch by source ConformationResult
(all BS accumulators, all MC accumulators, etc.) for natural commit
granularity.

Per batch, a fixed procedure:
1. Write the TrajectoryResult subclass(es) for this batch, following
   `BsWelfordTrajectoryResult` (AV) or `BsT0AutocorrelationTrajectoryResult`
   (FO + dense buffer) template from §4.
2. Add their output fields to `TrajectoryAtom` upfront. **Only finalized
   output fields — doubles, ints, SphericalTensors. No Welford
   instances, no accumulator-state objects. Per the §3 anti-pattern
   warning.**
3. Register their factories in a temporary "attach-all" helper
   (RunConfiguration infrastructure comes in 1d).
4. **Remove the corresponding old code path:** the specific Welford
   fields on `GromacsProteinAtom`, the specific lines in
   `AccumulateFrame` that updated them, the corresponding entries in
   `AllWelfords`. The migration is **replacement, not addition**.
   Parallel paths do not coexist post-commit.
5. Byte-parity test on 1 calibration protein: old output (from
   pre-push commit) vs new output (from this batch's commit) for
   the specific fields just migrated. Bit-identical for Welford-shaped;
   byte-identical for time series.

Batch ordering (smallest blast-radius first):
1. Scalar Welfords: BS T0, MC T0, Coulomb T0
2. Per-type Welford arrays: BS per-type, HM per-type
3. DeltaTracker family: BS T0 delta, AIMNet2 charge delta, SASA
   delta, water_n_first delta
4. TransitionCounter family: chi transitions, SS transitions
5. Water-environment Welfords
6. Hydration-geometry Welfords
7. Bond-length Welfords (via TrajectoryBond if that proposal lands,
   or via `BondLengthStatsTrajectoryResult` internal state if the
   alternative wins)
8. Dense-buffer time-series classes — positions, per-frame tensors
   for each calculator. Longest sub-batch; these are the current
   `AnalysisWriter::*_` buffers moving into typed Results.
9. Rich structured per-atom vectors (`RingNeighbourhoodTrajectoryStats`)
10. Scan-mode-specific classes (`DihedralBinTransitionTrajectoryResult`,
    `RmsdTrackingTrajectoryResult`, `ChiRotamerSelectionTrajectoryResult`
    + `SelectionEmittingTrajectoryResult` mixin per Appendix D)
11. FullFat-mode MOPAC family + `MopacVsFf14SbReconciliationTrajectoryResult`

**Checkpoint per batch:** old and new outputs agree byte-for-byte on
the migrated fields. Earlier-migrated batches remain green. Existing
test suite remains green. `GromacsProteinAtom` shrinks; `AccumulateFrame`
shrinks; `AllWelfords` shrinks. At end of 1b, `GromacsProteinAtom`
has zero Welford fields left and `AccumulateFrame` is an empty loop.

#### Sub-checkpoint 1c — Dissolve GromacsRunContext + rename

Now that the old accumulator content is gone, the structural rename
is mechanical.

Files:
- Delete: `src/GromacsRunContext.h/.cpp`. Bonded params → loaded by
  `GromacsEnsembleLoader::BuildFromTpr` into TrajectoryProtein.
  EDR preload → Trajectory (built in 1d). Cursor → internal to
  `GromacsFrameHandler` during Run only.
- Rename files: `src/GromacsProtein.h/.cpp` →
  `src/TrajectoryProtein.h/.cpp`; `src/GromacsProteinAtom.h` →
  `src/TrajectoryAtom.h` (the latter is actually an already-written
  new file; the rename here is: the old `GromacsProteinAtom.h`
  content is gone to bones, `TrajectoryAtom.h` exists with new
  content).
- Rename classes (grep-and-replace): `GromacsProtein` →
  `TrajectoryProtein`; `GromacsProteinAtom` → `TrajectoryAtom`.
- Update call sites in `src/`, `tests/`, `ui/`, `h5-reader/`.
- Move old `GromacsProtein.*` and `GromacsProteinAtom.h` content to
  `learn/bones/` with brief migration notes per §9 ("this was wrong,
  replaced by TrajectoryProtein / TrajectoryAtom, see
  spec/WIP_OBJECT_MODEL.md §3 anti-patterns subsection for the
  rationale").

Checkpoint: `grep -r "GromacsProtein\|GromacsProteinAtom\|GromacsRunContext\|AllWelfords\|AccumulateFrame" src/` returns zero hits. Tests pass.

#### Sub-checkpoint 1d — Trajectory + RunConfiguration + RunContext typed process entities

Files:
- New: `src/Trajectory.h/.cpp` — 5-phase `Run()` per §5. Owns
  preloaded EDR + cursor during Run. Holds `selections_` and frame
  metadata post-Run. `WriteH5()` emits `/trajectory/` group +
  `/metadata/source/` per §7.
- New: `src/RunConfiguration.h/.cpp` — class with three static
  factories per §6.
- New: `src/RunContext.h/.cpp` — class (not struct), constructor
  per §6 (no `mutable` cache).
- New: `src/SelectionEmittingTrajectoryResult.h` — mixin per
  Appendix D.
- Modified: `src/nmr_extract.cpp` (driver) — uses `Trajectory::Run`
  instead of direct GromacsFrameHandler dance.
- Tests: each RunConfiguration factory, extras-attach,
  selection-record collection.

Checkpoint: library drives extraction through Trajectory::Run. Full
test suite passes.

#### Sub-checkpoint 1e — AnalysisWriter dissolution + top-level H5 writer

Files:
- New: `src/WriteTrajectoryH5.h/.cpp` — orchestrator that calls
  `traj.WriteH5(file)` + iterates `tp.ResultsInAttachOrder()`
  calling `WriteH5Group` on each.
- Delete: `src/AnalysisWriter.h/.cpp` (to `learn/bones/`). Per-frame
  buffer logic is already in typed TrajectoryResults from 1b.
- Modified: `python/nmr_extract/_catalog.py` — update provenance /
  schema strings to match new group ownership.
- Tests: H5 output byte-parity vs archived old-pipeline output
  across all 10 calibration proteins.

Checkpoint at end of Stage 1: AnalysisWriter is gone. TrajectoryProtein
is the serialization source of truth. `grep` for every anti-pattern
token returns zero hits. Full test suite passes. 10-protein byte-
parity test passes.

---

### 🎯 Stage 2 — Activation

**Scope:** validation + fleet commit.

1. Run full test suite — all existing tests + new trajectory tests +
   validation-gate tests pass.
2. 10-protein byte-parity check against archived old-pipeline output.
3. Run 720-pair mutant rerun under the static single-PDB + ORCA +
   MutationDelta path (non-regression contract from §0). Verify
   output matches previous calibration baseline.
4. NamingRegistry fix activated in the same commit per
   `spec/ChangesRequiredBeforeProductionH5Run.md`.
5. Commit activation tag. Fleet extraction (685 proteins) runs from
   this commit.
6. Update WIP Amendments with activation date + commit SHA. Retire
   proposal-pending flags (TrajectoryBond resolution recorded).
   Fold WIP back into `OBJECT_MODEL.md` and `PATTERNS.md` per §0.

**Session budget:** 1 session.

---

### Deferred — `NmrAtomIdentity` (was Stage A)

**Moved to post-rollup** per user direction 2026-04-22 pm. Rationale:
adding new typed identity fields on atoms before cleaning up the
atom-level mess (GromacsProteinAtom's Welford fields) is out of
order. The Stage 1 refactor removes that mess; `NmrAtomIdentity` can
then land cleanly on top.

The design in §2 + Appendices A + E stays valid; it just happens
after activation (Stage 2), as a separate pre-fleet-ish additive
landing. Same activation-discipline commitments apply when that work
happens — byte-parity test, proposal-pending sign-off on specific
design, etc.

**This deferral does not block fleet extraction.** The 685-protein
fleet can run against the Stage 2 activation commit without
NmrAtomIdentity. When NmrAtomIdentity lands later, the fleet's H5
output does not need regeneration — the new `/atoms/identity/` group
can be computed retroactively from the (Protein, conformation) data
already in the H5 via a one-shot upgrade pass.

---

### Total estimated effort (revised)

| Stage | Work | Sessions |
|---|---|---|
| 1 (sub-checkpoints 1a–1e) | Trajectory-scope refactor, full replacement | 2–4 |
| 2 | Activation + validation | 1 |
| Deferred A | NmrAtomIdentity, separate later | 2–3 (whenever scheduled) |

Stage 1's 2–4 session range depends on batching discipline and
whether the implementing session preserves full context across
sittings (1M-context Opus runs suggest 2–3 sessions is typical for
this scope if design is locked).

### Staging discipline (unchanged from prior revision)

- **No commit preserves an anti-pattern-shaped state.** A commit where
  `TrajectoryAtom` has Welford fields, or where `GromacsProtein` and
  `TrajectoryProtein` coexist, or where `AccumulateFrame` is still
  called — is not a valid checkpoint. Either finish 1b / 1c / 1d in
  one commit group or defer the whole group.
- **No stage breaks earlier checkpoints.** Git as safety net.
- **Proposal-pending items block relevant stages.** TrajectoryBond
  gates the per-bond batch in 1b; NmrAtomIdentity is deferred so no
  sub-checkpoint in Stage 1 depends on it.
- **Test-suite health is the gate, not the goal.** Passes required
  at every checkpoint.

---

### Original stage labels (preserved as reference)

The original seven-stage plan (A through G) is preserved below for
reference; the consolidated Stage 1 / Stage 2 above is the authoritative
plan. Original stages map to sub-checkpoints:

- Original Stage A → Deferred (post-rollup)
- Original Stage B → Sub-checkpoint 1a
- Original Stage C → Sub-checkpoint 1b
- Original Stage D → Sub-checkpoint 1c
- Original Stage E → Sub-checkpoint 1d
- Original Stage F → Sub-checkpoint 1e
- Original Stage G → Stage 2

### Stage A — `NmrAtomIdentity` pre-bundle (purely additive)

**Design refs:** §2, Appendix A, Appendix E (which declines a
migration but leaves the scope of Stage A unchanged).

**Files:**
- New: `src/NmrAtomIdentity.h` — struct declaration + typed enums
  (`IupacAtomPosition`, `NmrClass`, `Locant`, `MethylGroup`,
  `MethyleneGroup`, `RingAtomRole`, `ResidueCategory`,
  `ChiParticipation` bitfield semantics, `SymmetryClassId`).
- New: `src/NmrAtomIdentityBuilder.h/.cpp` — derivation via
  per-residue full-template matching (§2 algorithm).
- New: `data/iupac_atom_dictionary.toml` — canonical residue
  template database, generated at library build time from PDB CCD
  (via cifpp) or AMBER residue library (Appendix A).
- New: `tools/generate_iupac_enum.py` — build-step generator
  producing `src/generated/IupacAtomPosition.h/.cpp` from the TOML.
- Modified: `src/Protein.h/.cpp` — add `atom_identities_` member,
  `AtomIdentity(i)` accessor, orthogonal indices
  (`AtomsByNmrClass`, `AtomsByMethylGroup`, `AtomsByRingAtomRole`,
  etc.).
- Modified: `src/Protein.cpp::FinalizeConstruction` — call
  `NmrAtomIdentityBuilder` after bond detection + ring detection.
- New: `tests/nmr_identity/` — six test categories per §2 validation
  discipline (coverage, geometry invariance, variant correctness,
  PDB CCD agreement, stereochemistry benchmarks, BMRB cross-check).
  Reference assertion lists committed as repo data.

**Checkpoint:** library compiles; all existing 293+ tests pass
unchanged (NmrAtomIdentity is additive, EnrichmentResult unchanged
per Appendix E); all six identity test categories pass on 10
calibration proteins. No behaviour changes for existing calculators.

**Estimated effort:** 2–3 sessions. Enum inventory + generator +
builder is the work; tests run the generator output against the
external references.

### Stage B — `TrajectoryResult` base + `TrajectoryAtom` + shim on `GromacsProtein`

**Design refs:** §3, §4.

**Files:**
- New: `src/TrajectoryResult.h/.cpp` — base class with virtuals
  `Name()`, `Dependencies()`, `Compute(conf, tp, idx, time_ps)`,
  `Finalize(tp)`, `WriteFeatures(tp, dir)`, `WriteH5Group(tp, file)`.
- New: `src/TrajectoryAtom.h` — declared-upfront typed field bundle
  with private constructor, friend `TrajectoryProtein`. Starts empty
  (field additions accumulate as TrajectoryResults land in Stage C).
- New: `src/DenseBuffer.h` — template per Appendix C with native T
  support for `Vec3` / `Mat3` / `SphericalTensor`.
- Modified: `src/GromacsProtein.h/.cpp` — **still named
  GromacsProtein at this stage** (rename in Stage D). Add:
  - `std::unordered_map<std::type_index, std::unique_ptr<TrajectoryResult>> traj_results_`
  - `AttachResult()` method mirroring `ProteinConformation::AttachResult`
  - template `Result<T>()`, `HasResult<T>()`, `ResultsInAttachOrder()`
  - dense-buffer adoption (`AdoptDenseBuffer<T>`)
  - Retain all existing `atoms_`, `bonds_accum_`, `AccumulateFrame()`
    code paths unchanged (shim mode — dual paths coexist).
- Tests: attach discipline (singleton + dependency check), dense-
  buffer ownership transfer, TrajectoryResult base dispatch.

**Checkpoint:** library compiles with new types, all existing tests
pass, new TrajectoryResult infrastructure is exercised by the new
tests but no TrajectoryResult subclasses exist yet. Old accumulation
path still drives all current output.

**Estimated effort:** 1–2 sessions.

### Stage C — Per-field migration of `GromacsProteinAtom` to `TrajectoryResult` subclasses (the grind)

**Design refs:** §4 (worked examples for the two patterns), Appendix F
(catalog of ~30 TrajectoryResult classes).

**Per-field subtask** (repeats ~30 times, one per Appendix F row):

1. Write the concrete TrajectoryResult subclass following
   `BsWelfordTrajectoryResult` template (for always-valid-mid-stream
   cases) or `BsT0AutocorrelationTrajectoryResult` template (for
   Finalize-only dense-buffer cases). Per §4.
2. Add that class's output fields to `TrajectoryAtom` upfront. Per §3.
3. Register its factory in the relevant `RunConfiguration`
   static factory (Stage E builds these; during Stage C use a
   temporary "attach everything" helper that enumerates factories
   for testing).
4. In `GromacsProtein::AccumulateFrame`, remove the old code path
   that wrote to `GromacsProteinAtom`'s Welford for that field. The
   new TrajectoryResult's `Compute` writes to TrajectoryAtom directly
   during the frame loop.
5. Byte-parity test on 1 calibration protein for this field group:
   old-pipeline output vs new-pipeline output, for the specific
   fields migrated in this subtask. Bit-identical required for
   Welford-shaped; byte-identical for time series.

**Sub-ordering** (smallest blast-radius first):
1. Scalar Welfords: BS T0, MC T0, Coulomb T0 (~3 subclasses)
2. Per-type Welford arrays: BS per-type, HM per-type (~2 subclasses)
3. DeltaTracker family: BS T0 delta, AIMNet2 charge delta, SASA
   delta, water_n_first delta (~4 subclasses)
4. TransitionCounter family: chi transitions, SS transitions,
   ring flip hints (~3 subclasses)
5. Water-environment Welfords (~2 subclasses)
6. Hydration geometry Welfords (~2 subclasses)
7. Bond-length Welfords (via `TrajectoryBond` per Appendix B — if
   that proposal is accepted; alternative design per Appendix B
   reviews if not)
8. Dense-buffer time-series classes (~15 subclasses — positions,
   per-frame tensors for every calculator). Longest subtask; these
   are the current `AnalysisWriter::*_` buffers moving into typed
   Results.
9. Rich-structured per-atom vectors (`RingNeighbourhoodTrajectoryStats`)
10. Scan-mode-specific classes (`DihedralBinTransitionTrajectoryResult`,
    `RmsdTrackingTrajectoryResult`, `DftPoseSelectionTrajectoryResult`
    + `SelectionEmittingTrajectoryResult` mixin per Appendix D)
11. FullFat-mode-specific classes (MOPAC time series,
    `MopacVsFf14SbReconciliationTrajectoryResult`)

**Checkpoint per migrated field group:** old calibration-protein
output and new output agree byte-for-byte on the migrated fields.
Earlier-migrated groups remain green. Existing test suite remains
green.

**Estimated effort:** 6–10 sessions. This is the majority of the
rollup's session budget. Granularity helps — checkpointable every
1–3 subtasks.

### Stage D — Dissolve `GromacsRunContext` and rename

**Design refs:** §5 (EDR preload timing), §6 (RunContext dissolution),
§9 (rename list).

**Files:**
- Delete: `src/GromacsRunContext.h/.cpp`.
  - Bonded parameters → owned by TrajectoryProtein (loaded by
    `GromacsEnsembleLoader::BuildFromTpr`).
  - Preloaded EDR frames → owned by Trajectory (built by Stage E).
    Placeholder location during Stage D: a `std::vector<GromacsEnergy>`
    member on GromacsProtein is fine as an intermediate step.
  - Cursor state → private to `GromacsFrameHandler` during Run only.
- Rename files: `src/GromacsProtein.h/.cpp` → `src/TrajectoryProtein.h/.cpp`;
  `src/GromacsProteinAtom.h` → `src/TrajectoryAtom.h`.
- Rename classes: `GromacsProtein` → `TrajectoryProtein`;
  `GromacsProteinAtom` → `TrajectoryAtom`.
- Update all call sites in `src/`, `tests/`, `ui/`, `h5-reader/`.
- Move `src/GromacsProteinAtom.h`'s historic content (dissolved
  Welfords) to `learn/bones/` with a brief migration note per §9.

**Checkpoint:** no `GromacsProtein` / `GromacsProteinAtom` /
`GromacsRunContext` remain in `src/`. All tests pass. Renames
are mechanical; grep for remaining references should return zero.

**Estimated effort:** 1–2 sessions.

### Stage E — `Trajectory` + `RunConfiguration` + `RunContext` (new typed process objects)

**Design refs:** §5, §6.

**Files:**
- New: `src/Trajectory.h/.cpp` — the process entity; 5-phase `Run()`
  method per §5; owns preloaded EDR + cursor during Run; holds
  `selections_` and frame metadata post-Run; `WriteH5()` emits
  `/trajectory/` H5 group.
- New: `src/RunConfiguration.h/.cpp` — the class with three static
  factories: `ScanForDftPointSet()`, `PerFrameExtractionSet()`,
  `FullFatFrameExtraction()`. Each factory's body enumerates the
  per-frame RunOptions + the TrajectoryResult factory list + the
  required ConformationResult type set.
- New: `src/RunContext.h/.cpp` — class (not struct), constructor
  takes `RunConfiguration` + output_dir; methods `AttachExtraResult`,
  `SetAimnet2Model`, accessors. Replaces the "attach-everything"
  helper from Stage C.
- New: `src/SelectionEmittingTrajectoryResult.h` — mixin interface
  per Appendix D.
- Modified: `src/nmr_extract.cpp` (driver) — main replaces direct
  `GromacsProtein` + `GromacsFrameHandler` dance with
  `Trajectory traj(...); tp = ...; traj.Run(tp, ctx);`.
- Tests: `Trajectory::Run` with each `RunConfiguration`, extras-attach
  path, dependency validation (ConformationResult types validated at
  Phase 2 per §5), EDR-preload timing, selection-record collection
  via `SelectionEmittingTrajectoryResult` dynamic_cast pass.

**Checkpoint:** library drives extraction through `Trajectory::Run`;
output H5 matches Stage D's expected shape. Full test suite passes.

**Estimated effort:** 2–3 sessions.

### Stage F — `AnalysisWriter` dissolution + top-level H5 writer

**Design refs:** §7.

**Files:**
- New: `src/WriteTrajectoryH5.h/.cpp` — top-level function
  `WriteTrajectoryH5(traj, tp, path)` that opens the H5 file, calls
  `traj.WriteH5(file)` for `/trajectory/*` + `/metadata/source/*`
  content, and iterates `tp.ResultsInAttachOrder()` calling each
  attached TrajectoryResult's `WriteH5Group(tp, file)`.
- Delete: `src/AnalysisWriter.h/.cpp` (move to `learn/bones/`).
  Per-field buffer logic migrated to individual
  `*TimeSeriesTrajectoryResult` classes in Stage C.
- Modified: `python/nmr_extract/_catalog.py` — update any
  `provenance` / schema-metadata strings to match the new group
  ownership. Most array paths don't change per §9.
- Tests: H5 output byte-parity vs archived old-pipeline output
  across all 10 calibration proteins.

**Checkpoint:** `AnalysisWriter` class no longer in the tree;
library's H5 output is entirely model-traversal-based (typed
entities self-serialising). Full test suite passes.

**Estimated effort:** 1–2 sessions.

### Stage G — Activation

**Design refs:** §10, `spec/ChangesRequiredBeforeProductionH5Run.md`.

1. Build library with all rollup changes.
2. Run full test suite — all 293+ existing tests + all new
   identity / trajectory tests + all validation-gate tests pass.
3. Run extraction on 1 calibration protein under
   `PerFrameExtractionSet`; byte-parity check against archived
   old-pipeline output.
4. Run on 2–3 more calibration proteins; same discipline.
5. Run the 720-pair mutant rerun under the static single-PDB +
   ORCA + MutationDelta path (non-regression contract from §0).
   Verify output matches previous calibration baseline.
6. Activate the NamingRegistry fix per
   `ChangesRequiredBeforeProductionH5Run.md` — same activation
   commit, bundled per that doc's discipline.
7. Commit activation tag. Fleet extraction (685 proteins) runs
   from this commit.
8. Update WIP Amendments with activation date + commit SHA. Retire
   proposal-pending flags on §2 / §3 / Appendix A / Appendix B
   (sign-off recorded in the activation commit). Fold WIP back
   into `OBJECT_MODEL.md` and `PATTERNS.md` per the §0 merge
   commitment.

**Checkpoint:** rollup is live; fleet extraction proceeds; WIP
retires.

**Estimated effort:** 1–2 sessions.

### Total estimated effort

| Stage | Sessions |
|---|---|
| A — NmrAtomIdentity pre-bundle | 2–3 |
| B — TrajectoryResult base + shim | 1–2 |
| C — Per-field migration (grind) | 6–10 |
| D — Rename + RunContext dissolve | 1–2 |
| E — Trajectory + RunConfiguration | 2–3 |
| F — AnalysisWriter dissolution | 1–2 |
| G — Activation | 1–2 |
| **Total** | **14–24** |

Fits the user's 12–20 session estimate (upper bound). Tighter if
Stage C's field list turns out shorter than Appendix F's catalog
suggests (some classes may merge); longer if proposal-pending items
(§2 NmrAtomIdentity design, §3 TrajectoryBond design, Appendix A
generator approach) need significant redesign after user review.

### Staging discipline

- **No stage breaks earlier stages' checkpoints.** If regression
  is detected between stages, the offending change is fixed in-stage
  before moving on. Git as safety net (rollback discipline); commit
  granularity aligned with checkpointable milestones.
- **Proposal-pending items block the stages that depend on them.**
  §2 NmrAtomIdentity (proposal-pending) gates Stage A; §3
  TrajectoryBond (proposal-pending) gates the per-bond subtask in
  Stage C. User sign-off required before those stages begin.
- **Test-suite health is a gate, not a goal.** The test suite is
  what enforces the non-regression contract (§0). Any stage boundary
  that fails its tests is not a stage boundary until the tests are
  green.

---

## Amendments

**2026-04-22 — revision pass after fresh-agent cold-read review.**
Cold-read review by a general-purpose agent (see conversation
transcript) surfaced specific gaps at the edges of the main-body
design. Main architectural direction accepted; edge gaps closed via:

- §2 (NmrAtomIdentity): added scope-framing paragraph; added
  enum-completeness discussion; added OpenBabel migration paragraph;
  pointers to Appendices A and E.
- §3 (TrajectoryProtein): added memory-math paragraph; added
  `TrajectoryBond` sub-section pointing at Appendix B.
- §4 (TrajectoryResult): added dense-buffer non-scalar-payload
  paragraph pointing at Appendix C; added selection-emission
  sub-section pointing at Appendix D; added per-frame-dispatch
  overhead paragraph; added catalog pointer to Appendix F.
- §5 (Trajectory::Run): added ordering-guarantees sub-section; added
  EDR preload timing sub-section.
- §6 (RunConfiguration / RunContext): replaced `mutable` cache in
  RunContext with comment explaining compute-fresh-each-call;
  added MOPAC sparse-frame note on `FullFatFrameExtraction`; pointer
  to Appendix F catalog.
- §10 (Coordinated activation): added test-strategy sub-section with
  byte-parity discipline; added schema-diff-as-prerequisite
  sub-section pointing at Appendix G; added scope-and-timeline
  sub-section per user framing, full quote preserved in Appendix H.
- §12 (Open questions): items 3 (selection callback shape) and 4
  (ConformationResult dependency validation) marked RESOLVED via
  Appendix D and §5 Phase 2 respectively. Item 6 added for
  frame-selection mechanism in `FullFatFrameExtraction` (filter-XTC
  vs skip-in-handler).
- **Appendices A–G added** at the end of the main body: A
  (IupacAtomPosition enum generator), B (TrajectoryBond), C
  (dense buffer layouts for non-scalar payloads), D (selection
  emission pattern), E (EnrichmentResult + OpenBabel disposition),
  F (TrajectoryResult catalog, ~30 classes), G (scope and timeline
  framing). (Original Appendix G on schema-diff process was removed
  during a subsequent revision pass — see below — and former
  Appendix H was renumbered to G.)

No architectural changes to §1, §7, §8, §9, §11. The design direction
and entity model are preserved; only specific gaps closed.

**2026-04-22 (continued) — further edits from same-day user review.**
Three prominent-note additions and two scope decisions:

- §0 (Status): added **non-regression contract** — static PDB paths
  (protonated/unprotonated single-structure, single-structure + ORCA,
  WT + ALA + ORCA + MutationDelta) must continue to work unchanged
  through the activation. Preservation commitment with test-enforcement
  discipline. This is load-bearing, not background.
- §2 (NmrAtomIdentity) + Appendix A: flagged **PROPOSAL PENDING USER
  REVIEW**. User confirmed the need for IUPAC-typed atom identity and
  OpenBabel; flagged that the specific design ("parachuted in with
  some complexity") is not yet user-approved. §2 content stays as
  proposal; pending sign-off before implementation.
- §3 (TrajectoryBond subsection) + Appendix B: flagged **PROPOSAL
  PENDING USER REVIEW**. User expressed concern; alternative design
  (per-bond state internal to a `BondStatsTrajectoryResult`, no
  parallel `TrajectoryBond` store) has equal weight and deserves
  discussion before implementation.
- §2 + Appendix E: **OpenBabel migration DECLINED.** User preference
  verbatim: *"I am happy to pay for openbabel per frame in exchange for
  not mucking with everything on that level."* EnrichmentResult stays
  unchanged; per-frame OpenBabel stays; ConformationAtom booleans
  stay. NmrAtomIdentity becomes a purely additive layer on Protein
  with zero disruption to existing calculators. Rollup scope narrows
  accordingly.
- Appendix G (schema-diff process) **removed**. User direction: formal
  schema-diff artifact + activation gate is corporate overhead that
  doesn't apply to our two-collaborator + git workflow. Schema
  changes documented in commit messages + Amendments. §10 test
  strategy updated. Former Appendix H (scope and timeline) renumbered
  to G.

**2026-04-22 (continued) — §2 derivation discipline + validation
discipline added after user review.**

- **Per-residue full-template matching** established as the preferred
  derivation algorithm for `NmrAtomIdentityBuilder` in §2. User
  preference: full-residue-at-a-time, full-match-required, in
  explicit preference to more efficient or opportunistic schemes.
  Anti-patterns enumerated (no SMARTS, no name-similarity
  heuristics, no partial assignment).
- **Validation discipline** for NmrAtomIdentity added in §2: four
  external reference sources priority-ordered (PDB CCD primary,
  AMBER cross-check, IUPAC rules baseline, BMRB secondary), six
  test categories including CCD agreement + stereochemistry
  correctness + BMRB cross-check, activation gate ties rollup to
  test-category pass before fleet runs.
- **Reference priority reordered from BMRB-first to CCD-first.** User
  clarification that the mutant pair rerun (720 pairs, ALA mutants
  do not have BMRB entries) goes through this library; calibration
  stats depend on typed per-atom stratification for muties; CCD
  covers muties uniformly, BMRB does not. CCD becomes the primary
  reference; BMRB becomes useful-but-secondary (10 calibration
  proteins only).
- **Runtime dependency discipline.** Brief C++ lib scan confirms
  cifpp (already a dep) and gemmi (new) both cover CCD access.
  Recommendation: pre-generate residue-template database at library
  build time from CCD (via cifpp) or AMBER residue library; ship as
  committed TOML/JSON repo data; no new runtime dep. Appendix A's
  generator approach extends cleanly.

**2026-04-22 (continued) — pre-second-cold-read revision pass.**
User-requested additions before a second fresh-agent review:

- **§0 non-regression contract — concrete test sentinels enumerated.**
  Replaced "existing test suite must pass" with a specific enumeration
  of test categories (static PDB load, conformation construction,
  OperationRunner::Run, per-calculator smoke tests, ORCA loading,
  RunMutantComparison, NPY output + SDK round-trip) and named the
  test-suite-at-stage-boundary enforcement mechanism per Appendix H.
- **§7 structure-ID emission subsection added.** Makes explicit the
  H5 group layout for per-Protein provenance (`/metadata/source/`)
  and per-atom identity (`/atoms/identity/`). Schema stability
  committed; enum-valued datasets carry `enum_values` sidecar
  attributes for decoder-free consumer access.
- **§7 library ↔ SDK boundary subsection added.** Names the
  contract explicitly — what the library produces, what it does not
  do, what SDK work it enables downstream (BMRB, validation benches,
  cross-protein pooling, per-atom residuals), what the SDK does not
  do. Post-activation SDK roadmap flagged (BMRB module, benches
  directory, stratify utilities) as out-of-rollup-scope but enabled
  by the rollup's output contract.
- **Appendix H added — implementation staging roadmap.** Seven
  checkpointable stages (A: NmrAtomIdentity pre-bundle, B:
  TrajectoryResult base + shim, C: per-field migration, D: rename +
  RunContext dissolve, E: Trajectory + RunConfiguration, F:
  AnalysisWriter dissolution, G: activation). Per-stage files +
  checkpoint + estimated effort. Stage ordering respects
  dependencies; proposal-pending items gate their relevant stages.
- **Cleanup: duplicate Appendix G labels resolved.** Former
  "Appendix G — REMOVED" stub removed (its history is captured in
  this Amendments section). Scope/timeline framing stays as
  Appendix G; implementation staging is the new Appendix H.

**2026-04-22 (continued) — anti-pattern warnings added; staging
consolidated; NmrAtomIdentity deferred.** User concern:
*"It is VERY possible the next new session will talk to me for an hour,
and then just go implement exactly what is there now with different
names and C++ off stack exchange because they know better. I am quite
scared actually — you actually get it but it is not in tune enough
with your training the next you will."*

Three substantive additions to head off that failure mode:

1. **New prominent subsection at the top of §3: "🛑 READ THIS FIRST —
   ANTI-PATTERNS THE CURRENT CODE EMBODIES THAT MUST NOT SURVIVE
   THE REFACTOR."** Four anti-patterns enumerated in detail with
   their current-code form, explicit reasoning for why each is
   wrong, what must replace each, and the specific "refactor
   temptation to avoid" that a future session is likely to
   rationalise. Plus a "test for am-I-doing-it-right" checklist
   (`grep` commands that must return zero hits post-refactor) plus
   recommended file-header comment blocks for `TrajectoryAtom.h`
   and `TrajectoryResult.h` that embed the architectural commitment
   in the code itself — comments that survive into the implementing
   session's context when the spec is not loaded.
2. **Appendix H staging consolidated into two stages.** Original
   seven-stage plan compressed into: Stage 1 = trajectory-scope
   refactor as one focused push (sub-checkpoints 1a–1e covering
   original B through F), Stage 2 = activation (original G). Driver:
   avoid an intermediate shim state that becomes a local optimum
   and never completes the replacement. No commit within Stage 1
   preserves an anti-pattern-shaped state (coexisting old + new
   paths, Welford fields on TrajectoryAtom, etc.). Original stage
   labels A–G preserved as a reference map.
3. **NmrAtomIdentity (original Stage A) deferred to post-rollup.**
   User's reasoning: adding typed identity fields on top of the
   current atom-level mess (GromacsProteinAtom's Welford fields)
   is out of order. Clean up first; additive typed identity
   second. The §2 + Appendix A + Appendix E design stays valid
   for when the deferred work happens; the deferral does not block
   fleet extraction (NmrAtomIdentity can be computed retroactively
   post-activation from data already in the H5).

Session budget revised: Stage 1 is 2–4 focused claude-code sessions
for the trajectory refactor push; Stage 2 is 1 session for
activation. Deferred A is 2–3 sessions whenever scheduled.

**Explicit guidance to the implementing session.** This revision
makes the WIP a tool for preventing re-emergence of the old pattern,
not just a design for the new pattern. The implementing session is
expected to re-read §3 anti-patterns before each sub-checkpoint and
to verify the test-for-am-I-doing-it-right `grep` commands after
each commit within Stage 1.

**2026-04-22 (continued) — §2 motivation reframed after user review.**
User correction: I had framed NmrAtomIdentity as "infrastructure for
SDK-side BMRB binding" in a proposal revision. User pushed back:
*"This is protein topology. Open babel aside, we have explicit rules
about introducing invariants to conformations. It is not only
important to downstream, it is also an important part of some of the
tensors from literature we talked about."* §2 rewritten to lead with
the topology-invariants-on-Protein architectural rule (per
CONSTITUTION) and with concrete library-internal physics uses —
per-atom-type calibration stratification (Stage 1, settled),
Havlin-Oldfield CSA sensitivities to (φ, ψ, χ₁), Yao-Bax α/β 15N
split, Loth ubiquitin CCR bench, Kasinath-Wand methyl entropy,
Brender-Taylor-Ramamoorthy tensor orientations, Pauling / Babaei
bulk χ. SDK consumption is named as a parallel downstream consumer,
not as the primary case. The library needs the typed identity for
its own physics work independent of any integration-with-external-data
concern. Proposal-pending flag stays; specific enum design / generator
approach still open.
