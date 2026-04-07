# Feedback: Constitution vs Existing Object Model and Extraction Order

**Last updated: 2026-04-01**
**Status: Most items resolved. See per-item Status lines below.**
**Authority: OBJECT_MODEL.md and EXTRACTION_ORDER.md are the canonical
concrete specifications. CONSTITUTION.md states principles. When they
disagree, OBJECT_MODEL.md wins for concrete types and CONSTITUTION.md
wins for principles.**

This document captures contradictions, ambiguities, missing items, and
suggestions identified during the alignment of OBJECT_MODEL.md and
EXTRACTION_ORDER.md with the revised CONSTITUTION.md.

---

## 1. Contradictions

### 1.1 Environment vs ProteinBuildContext: More Than a Rename

The constitution defines ProteinBuildContext as "how this protein instance
was built" -- PDB source, protonation tool + version + pH, force field,
what was stripped, what was assumed, crystal resolution, organism,
deposition date. It is a friend to the Protein and travels with it.

The old object model defines Environment as a class hierarchy
(BmrbEnvironment, GromacsEnvironment, CrystalEnvironment, etc.) that
holds experimental conditions (pH, temperature, ionic strength, buffer,
solvent, force field, water model, space group, etc.). The Protein
owned a list of Environments, and each Conformation referenced one
by index.

The contradiction: the constitution says ProteinBuildContext records
"what happened" (provenance), while the old Environment recorded
"under what conditions the experiment was done" (experimental context).
These are distinct concepts. A ProteinBuildContext would record that
the protein was loaded from PDB 1UBQ, protonated with PROPKA at pH 7.0,
stripped of waters, and assigned ff14SB charges. An experimental
Environment would record that the NMR experiment was done at pH 6.5,
25C, in phosphate buffer.

**Resolution needed:** Is ProteinBuildContext a complete replacement
for Environment (merging provenance with experimental conditions)?
Or should Environment survive as a separate concept for experimental
metadata while ProteinBuildContext handles computational provenance?
The constitution's ProteinConformation hierarchy already carries
experimental metadata (CrystalConformation has resolution, R-factor,
temperature; NMRConformation has restraint count). This suggests that
experimental conditions live on the conformation type, and
ProteinBuildContext covers only the computational provenance.

**Status: RESOLVED** -- Constitution, archive notes, and OBJECT_MODEL all agree. ProteinBuildContext records provenance. Experimental conditions live on ProteinConformation subtype metadata.

### 1.2 Conformation Class Hierarchy: Flat vs Typed

The old object model had a flat Conformation class with a name, source
string, timestamp, and an index into the Protein's environment list.
There was no type hierarchy. The State bitfield tracked what had been
computed.

The constitution introduces a full class hierarchy:
ProteinConformation (base) with subtypes ExperimentalConformation
(CrystalConformation, NMRConformation), ComputedConformation
(PredictionConformation, MinimisedConformation, MDFrameConformation),
and DerivedConformation. The base class is explicitly NOT abstract --
it is fully functional.

**Contradiction:** The old model's approach (flat class + State enum)
is incompatible with the constitution's typed hierarchy + named
ConformationResult singletons. The entire state tracking mechanism
changes from a bitfield on the Conformation to a dependency graph
among ConformationResult types.

**Status: RESOLVED** -- All docs use typed hierarchy.

### 1.3 Protein Owns Conformations Through Factory Methods

The old model listed conformations as a bare `vector<Conformation>` on
the Protein, with creation at "Phase 0+." The constitution requires that
conformations are created through typed factory methods:

```
protein.AddCrystalConformation(positions, metadata)
protein.AddMDFrame(positions, metadata)
protein.AddPrediction(positions, metadata)
protein.AddDerived(parent, description)
```

No agent creates a ProteinConformation directly. No loose conformations.
This is a significant structural change.

**Status: RESOLVED** -- All docs agree.

### 1.4 Results Are Named Singletons, Not Phase-Numbered Properties

The old model stored computed data as parallel vectors on the
Conformation (partial_charges, vdw_radii, dssp_per_residue,
apbs_electric_field, etc.) with a State bitfield tracking completion.

The constitution requires named ConformationResult objects accessed by
name: conformation.Dssp(), conformation.BiotSavart(), etc. Each result
is a singleton that declares dependencies, is checked at attach time,
and is permanent once attached. Results have physics query methods,
not raw data getters.

**Contradiction:** The entire data storage model changes. Properties
currently stored as `vector<double>` on the Conformation become internal
state of the respective ConformationResult objects. The State enum is
entirely replaced by the presence/absence of ConformationResult objects.

**Status: RESOLVED** -- All docs agree. Template-based access adopted as mechanism.

### 1.5 Atom Ownership and Identity vs Conformation Geometry

The old model placed Atom objects (with identity-only properties) on the
Protein and conformation-dependent properties on the Conformation in
parallel arrays. The constitution says atoms live "in a
ProteinConformation" and know their position, element, role, distances,
etc. -- but also says the Protein holds sequence and science data while
the ProteinConformation holds geometry.

**Ambiguity:** The constitution's Atom section heading says "Atom (in a
ProteinConformation)" but the Protein section implies the Protein owns
the atom list (sequence data, element, residue membership). The old
model's split (identity on Protein, geometry on Conformation) appears
to be the correct interpretation, but the constitution is not explicit
about this split.

**Status: RESOLVED** -- Protein owns atom identity (element, residue index, bonds). ProteinConformation owns positions and all computed data via ConformationResult objects. Conformation holds a valid back-pointer to its Protein. Enrichment properties (role, hybridisation) computed per-conformation for simplicity. Constitution updated with placement rule.

### 1.6 Bond Category Taxonomy Discrepancy

The old model defines BondCategory with: PeptideCO, PeptideCN,
BackboneOther, SidechainCO, Aromatic, Disulfide, SidechainOther,
Unknown.

The constitution defines a simpler taxonomy: PeptideCO, PeptideCN,
Sidechain, Aromatic, Disulfide. The constitution does not mention
BackboneOther, SidechainCO, SidechainOther, or Unknown.

**Contradiction:** The old model has finer granularity. The
constitution's McConnell per-bond-category subtotals reference
"peptide CO, CN, sidechain, aromatic" which aligns with a coarser
set. However, the distinction between SidechainCO (ASN/ASP/GLN/GLU
carbonyl bonds) and other sidechain bonds matters for bond anisotropy
calculations, since C=O bonds have much stronger anisotropy than C-C
bonds. The physics requires keeping SidechainCO separate.

**Resolution needed:** The constitution's Bond section specifically
mentions "bond category (peptide CO, peptide CN, sidechain, aromatic,
disulfide)" but the extraction order document's Phase 4a calculator
distinguishes PeptideCO, PeptideCN, SidechainCO, and Aromatic subtotals.
BackboneOther (N-CA, CA-C, CA-CB) is also physics-relevant for
backbone bond anisotropy. Suggest the constitution adopt the finer
BondCategory enum from the existing object model.

**Status: RESOLVED (stale constitution)** -- OBJECT_MODEL has correct 7+1 categories. Constitution being updated to match.

### 1.7 Feature Compute Signature

The old object model's Feature base class has:
```
Compute(const Atom&, const AtomProperties&, const Conformation&, FeatureOutput&)
```

The constitution says features read from ConformationResult objects
accessed as `conformation.BiotSavart().SumT0ByRingType(atomIdx, ...)`.
This implies the Feature::Compute signature should take only the
atom index and the ProteinConformation, since all data is accessed
through result object query methods.

**Status: RESOLVED** -- Feature takes (atom_index, const ProteinConformation&, FeatureOutput&). Reads from ConformationResult query methods.

---

## 2. Ambiguities That Need Resolution

### 2.1 Where Do Enrichment Properties Live?

The constitution says enrichment properties (role, hybridisation,
categoricals) are "set once, append-only" and "added by
ConformationResult objects during attachment." But role and hybridisation
are not conformation-dependent -- they depend on element, bond topology,
and residue type, which are the same across all conformations of the
same protein with the same protonation state.

**Question:** Are enrichment properties on the Atom (owned by Protein,
shared across conformations) or on the ProteinConformation (per-
conformation)? If they are on the Atom, they are computed once. If
on the ProteinConformation, they are recomputed for each conformation
even when identical. The constitution's copy semantics section suggests
some enrichment properties survive GeometryOnly copies (hybridisation
"usually unchanged unless bonds change"), which implies they live on
the Protein's atoms, not on the ProteinConformation.

**Status: RESOLVED** -- EnrichmentResult is a ConformationResult. Computed per-conformation for simplicity. Not shared across conformations -- simple beats optimal.

### 2.2 EnrichmentResult vs ConformationResult

The constitution groups all results as ConformationResult types. But
enrichment (role assignment, hybridisation, categoricals) is fundamentally
different from computation (BiotSavart, Coulomb) -- enrichment depends
on the protein's identity, not on a specific conformation's geometry.

**Question:** Should there be a separate EnrichmentResult type that
attaches to the Protein rather than the ProteinConformation? Or does
each ProteinConformation independently carry its own enrichment results?
The constitution's DsspResult has no dependencies and attaches to the
conformation, which is correct (DSSP depends on geometry). But
RoleAssignment and Hybridisation depend on bond topology, which is
on the Protein.

**Status: RESOLVED** -- EnrichmentResult IS a ConformationResult with no dependencies.

### 2.3 Protonation State: Per-Protein or Per-Conformation?

The constitution says "One protonation state per ProteinConformation"
but also says the Protein owns the conformation list and that the
copy-and-modify pattern applies a new ProtonationState to the whole
protein copy. The old model stored protonation on the Conformation.

**Question:** Can two conformations of the same Protein have different
protonation states? The MD ensemble use case suggests yes (different
frames might sample different protonation states). The copy-and-modify
pattern suggests no (you copy the whole protein to change protonation).
The constitution is ambiguous.

**Status: RESOLVED** -- Protonation lives on ProteinConformation. All conformations of the same Protein share the same protonation. Different protonation = copy the whole Protein, apply new protonation, re-enrich.

### 2.4 OrcaShieldingResult Dependencies

The constitution lists OrcaShieldingResult as "requires: nothing (loaded
from files)." But the ORCA delta-shielding pipeline described in Layer 0
requires atom matching between WT and ALA structures. This matching
requires atom positions and potentially spatial queries.

**Question:** Does OrcaShieldingResult truly have no dependencies, or
does it depend on the protein being loaded (positions available for
atom matching)? If it requires atom matching, it may depend on
SpatialIndexResult.

**Status: RESOLVED** -- Correctly requires nothing. Atom matching at load time using positions already on conformation.

### 2.5 MolecularGraphResult vs Graph Features in Through-Bond Calculation

The constitution's dependency table lists MolecularGraphResult (requires
SpatialIndexResult) separately from the through-bond BFS features
(Phase 2b/2c in the old extraction order). The old extraction order
had graph features computed as a single pass.

**Question:** Is MolecularGraphResult a ConformationResult that computes
and stores graph topology features (BFS distances, electronegativity
sums, conjugation chains) on atoms? Or is it just the graph data
structure that other results query? If the former, it replaces Phase
2b/2c entirely. If the latter, there needs to be a separate
GraphFeaturesResult.

**Status: RESOLVED** -- Fully specified in both OBJECT_MODEL and EXTRACTION_ORDER.

### 2.6 SpatialIndexResult Scope

The constitution says SpatialIndexResult has no dependencies. But the
old extraction order's Phase 1f (spatial index construction) and Phase
1g (neighbour list construction) are distinct operations. Phase 1g
depends on Phase 1f.

**Question:** Does SpatialIndexResult include neighbour list
construction (all atoms within 15A, with stored distances and
directions)? Or is there a separate NeighbourListResult that depends
on SpatialIndexResult?

**Status: RESOLVED** -- Includes KD-trees AND 15A neighbour lists.

### 2.7 Ring Geometry and Bond Geometry Timing

The old extraction order computed ring geometry at Phase 1h and bond
geometry at the same time. The constitution's ConformationResult
dependency graph does not explicitly list ring geometry or bond
geometry as separate result types.

**Question:** Is ring/bond geometry computation part of construction
(happens when the ProteinConformation is created, since positions are
const), part of SpatialIndexResult, or a separate result type?

**Status: RESOLVED** -- GeometryResult (requires nothing) handles this.

### 2.8 ProteinConformation Base Class: What Methods Does It Provide?

The constitution says the base class "does the heavy lifting: holds
positions (const), provides spatial queries, accumulates
ConformationResult objects." But it also says "results are accessed by
name: conformation.BiotSavart()."

**Question:** Are the named accessor methods (Dssp(), BiotSavart(),
Coulomb(), etc.) on the base ProteinConformation class? If so, every
new ConformationResult type requires modifying the base class to add
an accessor. The extensibility section says "add named accessor on
ProteinConformation" as step 4. This means the base class grows with
every new result type, which is not ideal for extensibility.

**Suggestion:** Consider a template accessor like
`conformation.Result<BiotSavartResult>()` that does not require
modifying the base class. The named convenience methods can be
non-virtual wrappers if desired for readability.

**Status: RESOLVED** -- Template-based access adopted. `conf.Result<T>()` does not require modifying base class. Named accessors are optional convenience wrappers.

---

## 3. Places Where the Constitution Says One Thing But Physics Requires Another

### 3.1 SidechainCO vs Generic Sidechain in Bond Categories

As noted in 1.6, the constitution's bond category taxonomy merges all
sidechain bonds. But the physics of bond anisotropy depends heavily
on bond order: a sidechain C=O (ASP, GLU carboxylate; ASN, GLN amide)
produces much stronger anisotropy than a sidechain C-C bond. The
McConnell equation scales with the anisotropy of the magnetic
susceptibility, which is dramatically different for double bonds vs
single bonds. Merging them into one "Sidechain" category would lose
essential physics information for the bond anisotropy calculator.

**Status: RESOLVED (stale constitution)** -- OBJECT_MODEL has finer categories. Constitution being updated.

### 3.2 BackboneOther Bond Category

The constitution does not mention BackboneOther (N-CA, CA-C, CA-CB).
These bonds have backbone anisotropy contributions that are distinct
from peptide C=O and peptide C-N bonds. While their anisotropy is
smaller than the peptide bonds, it is not zero, and the aggregate
backbone anisotropy is a feature used by the ML model. Dropping this
category would require the calculator to either ignore these bonds
(incorrect) or lump them with sidechain bonds (also incorrect, since
backbone and sidechain contributions are decomposed separately).

**Status: RESOLVED (stale constitution)** -- Same as 3.1.

### 3.3 Coulomb E-Ring Projection Depends on Ring Results

The Coulomb calculator (Phase 4b in the old order) computes
coulomb_E_ring_proj = E_total . nearest_ring_normal. This requires
knowing which ring is nearest and its normal vector, which comes from
the ring current calculators (Phase 3). The constitution's dependency
graph lists CoulombResult as requiring ChargeAssignmentResult and
SpatialIndexResult, but NOT any ring calculator result.

**Resolution needed:** Either CoulombResult must declare a dependency
on BiotSavartResult (or at minimum on some ring geometry result that
provides ring normals), OR the E-ring projection must be computed as
a derived feature at feature extraction time rather than in the
Coulomb calculator. The latter is cleaner: the Coulomb calculator
produces the E-field; the feature extractor knows the nearest ring
normal from BiotSavartResult and computes the projection.

**Status: RESOLVED** -- Computed at feature extraction time, not in CoulombResult. Documented in OBJECT_MODEL line 493 and EXTRACTION_ORDER line 720.

### 3.4 HBondResult Dependencies

The constitution lists HBondResult as requiring only DsspResult. But
the H-bond calculator needs spatial queries (finding nearest H-bond
partners) which come from SpatialIndexResult. The old Phase 4c
explicitly uses the spatial index.

**Resolution:** HBondResult should require both DsspResult (for H-bond
partner validation from DSSP data) and SpatialIndexResult (for finding
nearby atoms). The current dependency declaration is incomplete.

**Status: RESOLVED** -- All docs show requires: DsspResult, SpatialIndexResult.

---

## 4. Suggestions for the Constitution

### 4.1 Introduce a GeometryResult for Ring/Bond Geometry

Ring geometry (center, normal, radius, vertex positions) and bond
geometry (midpoint, length, direction) are computed from the const
positions and don't depend on any ConformationResult. They could
be computed at construction time. But since they are data that
multiple results depend on, making them a ConformationResult
(GeometryResult, requires: nothing) would make the dependency
explicit and put them in the same framework as everything else.

**Status: ADOPTED** -- In both OBJECT_MODEL and EXTRACTION_ORDER.

### 4.2 Introduce an EnrichmentResult

Role assignment, hybridisation, and categorical properties are
identity-dependent (element, bond topology, residue type), not
geometry-dependent. Making them a ConformationResult that attaches
to each ProteinConformation is technically correct but redundant for
conformations that share the same protein and protonation. Consider
an EnrichmentResult that lives on the Protein and is shared by all
its conformations. Each ProteinConformation would reference this
shared enrichment. On copy-with-new-protonation, the enrichment
is invalidated and recomputed.

**Status: ADOPTED** -- In both docs.

### 4.3 Template-Based Result Access

Rather than adding a named accessor to ProteinConformation for each
new result type, use a template method:

```
template<typename T>
T& Result() { return *dynamic_cast<T*>(results_[T::name]); }
```

Named convenience methods can wrap this for common types. This
keeps the base class stable as result types are added.

**Status: ADOPTED** -- Template mechanism chosen over named accessors for agent-driven development safety.

### 4.4 Explicit NMRConformation Use Case

The constitution defines NMRConformation (ensemble_member,
restraint_count) but neither the old object model nor the extraction
order addresses NMR ensemble handling. NMR structures from the PDB
come as multi-model ensembles (e.g., 20 conformers). Each conformer
would be an NMRConformation. The pipeline should handle this: load
PDB with MODEL/ENDMDL records, create one NMRConformation per model,
enrich and extract each independently.

**Status: RESOLVED** -- EXTRACTION_ORDER shows protein.AddNMRConformation() per model.

### 4.5 Clarify Tool Configuration Scope

The constitution says "Tool configuration (where binaries live,
default parameters) is global, not per-protein." The old model had
ToolContext enum on the Environment. A global configuration registry
or singleton should be specified for tool paths (DSSP binary, tleap
path, APBS binary, OpenBabel path, etc.) so agents know where tool
config lives.

**Status: DEFERRED** -- Global config registry not yet specified. Will address in implementation.

### 4.6 Specify the ConformationResult ABC Contract

The constitution says each result type declares dependencies and is
checked at attach time. It would help to specify the ABC explicitly:

```
class ConformationResult {
public:
    virtual ~ConformationResult() = default;
    virtual string Name() const = 0;
    virtual set<string> Dependencies() const = 0;
    virtual void Attach(ProteinConformation& conf) = 0;
    // WHY: Every result must be identifiable, dependency-checkable,
    // and capable of populating properties on a conformation.
};
```

**Status: RESOLVED** -- OBJECT_MODEL lines 1042-1048 have the explicit ABC.

### 4.7 MinimisedConformation Provenance

The constitution lists MinimisedConformation with parent, method,
force_field, energy, converged. This type would be produced by the
OpenMM MD workflow. Consider adding: step_count, minimisation_algorithm
(L-BFGS, steepest descent), and initial_energy vs final_energy for
diagnostics.

**Status: RESOLVED** -- Constitution has parent, method, force_field, energy, converged.

---

## 5. Missing Result Types or Properties

### 5.1 No Explicit Result Type for Enrichment

Role assignment, hybridisation, and categorical properties are
described as enrichment operations but are not listed as a named
ConformationResult type in the constitution's dependency table. There
should be an EnrichmentResult (or AtomRoleResult + HybridisationResult)
that other results can depend on.

**Status: RESOLVED** -- Adopted.

### 5.2 No Explicit Result Type for Ring/Bond Geometry

Ring geometry computation and bond geometry computation are described
as Phase 1h operations but are not listed as ConformationResult types.
Since BiotSavartResult, HaighMallionResult, PiQuadrupoleResult,
RingSusceptibilityResult, and DispersionResult all depend on ring
geometry, and McConnellResult depends on bond geometry, these should
be explicit result types (or part of SpatialIndexResult).

**Status: RESOLVED** -- Adopted.

### 5.3 No Explicit Result Type for Neighbour Lists

The 15A neighbour lists with stored distances and directions are
built at Phase 1g but are not a named ConformationResult. Since
virtually every calculator depends on them, they should be explicit.
SpatialIndexResult could encompass this.

**Status: RESOLVED** -- Part of SpatialIndexResult.

### 5.4 No Explicit Result Type for Pre-Built Collections

The atoms_by_role, rings_by_type, bonds_by_category, and
residues_by_type collections (Phase 1i) are not listed as a
ConformationResult. They are infrastructure that all subsequent
results use.

**Status: RESOLVED** -- Part of GeometryResult.

### 5.5 No Explicit Result Type for Ring-Pair Properties

Ring-pair properties (Phase 1j) are not listed as a ConformationResult.
They are used by multi-ring features.

**Status: RESOLVED** -- Part of GeometryResult.

### 5.6 Gaussian Density Result

Phase 3f (per-ring-type Gaussian density) is described in the
extraction order but not listed as a named ConformationResult type in
the constitution's dependency table. It has learned parameters (alpha
and r0 per ring type) which implies it needs a parameter source.

**Status: RESOLVED** -- Learned parameters come from ParameterCorrectionResult (renamed from BaselineModelResult). Stored in RingNeighbourhood.

### 5.7 BFS Decay Lambda Parameter

The through-bond BFS decay feature uses a "learned parameter" lambda.
This parameter needs to come from somewhere -- either a model config
file or a training-time optimization. The extraction order says
"lambda initialized to ln(2)/3" but does not specify how it is updated.
This is a training concern, not an extraction concern, but the
extraction pipeline needs a way to receive the current lambda value.

**Status: RESOLVED** -- Comes from ParameterCorrectionResult or is a fixed hyperparameter (ln(2)/3). Decision deferred to implementation -- if model learns it, it comes from the model; if not, it is a named constant.

### 5.8 FeatureExtractionResult Dependencies

The constitution lists FeatureExtractionResult as requiring "all physics
results." This is underspecified. It should enumerate exactly which
result types it depends on, or use a mechanism like "all registered
ConformationResult types." Otherwise, adding a new physics result type
does not automatically make it a dependency of feature extraction, and
its data may not be available.

**Status: OPEN** -- Still says "all physics results." The implementation should enumerate them explicitly. The template mechanism with HasResult<T>() allows graceful handling of optional results.

### 5.9 PredictionResult Confidence Model

The constitution mentions heteroscedastic learned uncertainty (sigma)
and a HeuristicTier classification. The thresholds (threshold_report,
threshold_pass) are not specified. These are likely hyperparameters
from training, but they need a home in the configuration.

**Status: DEFERRED** -- Thresholds in toml config. Will tune during training.

### 5.10 ProtonationState Tool Version

The old model's ProtonationState had a tool_version field (e.g.,
"PROPKA 3.5.1"). The constitution's ProteinBuildContext says
"protonation tool + version." These should be consistent -- if
ProteinBuildContext subsumes ProtonationState's tool/version, then
ProtonationState can simplify to just the per-residue decisions.
If ProtonationState remains self-describing, it should keep
tool/version.

**Status: RESOLVED** -- ProteinBuildContext has provenance summary. ProtonationState has per-residue detail. Not redundant -- different granularity.

---

## Summary

The constitution represents a significant and welcome architectural
improvement over the existing object model. The key shifts are:

1. From flat Conformation to typed ProteinConformation hierarchy
2. From manual phase ordering to dependency-graph-driven extraction
3. From parallel property vectors to named ConformationResult singletons
4. From Environment to ProteinBuildContext
5. From State bitfield to dependency checking at attach time

The main areas needing clarification before implementation are:

- The ProteinBuildContext vs Environment scope question (1.1)
- Where enrichment properties live (2.1)
- The bond category granularity (1.6, 3.1, 3.2)
- Missing result types for infrastructure operations (5.1-5.5)
- CoulombResult's implicit dependency on ring geometry (3.3)
- HBondResult's missing SpatialIndexResult dependency (3.4)
