# Constitution: Proposed Edits

**Date**: 2026-04-03
**Status**: Draft for review. Each section is replacement text for the
corresponding section in CONSTITUTION.md. Notes explain the reasoning.

---

## 1. Protonation State (replaces lines 49-52)

**Note**: The original says "one protonation state per
ProteinConformation." The code and protonation design history
established that protonation determines the atom list, which is
per-Protein. The ProteinConformation hierarchy is about geometry
and provenance, not chemistry.

### Proposed text:

> ### Protonation State
>
> Protonation determines the atom list. Which titratable groups are
> charged, which histidine tautomer exists, which cysteines form
> disulfide bonds — these decisions change which atoms are present
> (adding/removing hydrogens changes the flat index array). Therefore:
>
> **Different protonation = different Protein.** Not a different
> conformation of the same protein. The atom indices, ring vertex
> indices, bond endpoints, and every ConformationAtom field are
> indexed into the flat atom array. Changing that array invalidates
> everything.
>
> A ProteinConformation holds geometry for a *fixed* atom list. Many
> conformations can share the same Protein (same protonation, different
> poses — this is MD ensembles). But two protonation states require
> two Protein instances, each with their own conformations, each
> independently enriched and extracted.
>
> Protonation state is recorded in the ProteinBuildContext (which tool,
> which pH, which decisions) and in the Residue variant indices (which
> variant was applied to each titratable residue).

---

## 2. Protein-to-Conformation Relationship (new subsection, after
   ProteinConformation hierarchy)

**Note**: This makes the back-pointer pattern idiomatic rather than
suspicious. A future agent seeing `conf.ProteinRef().AtomAt(i).element`
should understand this is the designed access path, not a violation of
encapsulation.

### Proposed text:

> ### Accessing the Protein from a Conformation
>
> A ProteinConformation always holds a valid pointer to its owning
> Protein. This is the designed access path for identity information:
>
> ```cpp
> // Element, bonds, residue type — identity, not geometry
> const auto& protein = conf.ProteinRef();
> Element elem = protein.AtomAt(i).element;
> const Residue& res = protein.ResidueAt(protein.AtomAt(i).residue_index);
> const Ring& ring = protein.RingAt(ri);
>
> // Position, computed fields — geometry, per-conformation
> Vec3 pos = conf.PositionAt(i);
> const auto& ca = conf.AtomAt(i);  // ConformationAtom
> ```
>
> The Protein is non-copyable and non-movable. It lives on the heap
> via unique_ptr. Conformations are owned by the Protein and cannot
> outlive it. The raw pointer is safe because the lifetime guarantee
> is structural, not reference-counted.
>
> **Do not** duplicate identity information onto ConformationAtom to
> avoid the back-pointer traversal. Element, residue type, bond
> connectivity, ring membership — these live on Protein and are
> accessed through ProteinRef(). This is not a code smell. This is
> the separation between what the molecule IS (Protein) and where
> its atoms ARE (ProteinConformation).

---

## 3. Rebuild Pattern (replaces Copy Semantics section, lines 551-635)

**Note**: The original described a copy-and-modify pattern with
GeometryOnly/Full copy policies. Protein is non-copyable and the
protonation design history explains why: different protonation =
different atom list = different indices = rebuild, not copy. The
*principle* (same geometry, different chemistry, delta is the
experiment) is preserved. The *mechanism* is corrected.

### Proposed text:

> ## Protein Rebuild Pattern
>
> ### The principle
>
> It must be possible, tested, and routine to analyse the same
> protein geometry under different chemical conditions:
>
> - Different protonation (pH scanning)
> - Different mutation (WT vs ALA)
> - Different charge model (ff14SB vs xTB vs CHARMM36m)
> - Different force field parameters
>
> The delta between two conditions IS the experiment. This is how
> the training data works: WT protein → full extraction. ALA mutant
> → full extraction. Delta in DFT shielding → training target. Delta
> in classical calculations → what the model should learn to predict.
>
> ### The mechanism: rebuild, not copy
>
> Different protonation changes the atom list (adding/removing H
> atoms). Different atom list means different flat indices. Every
> ring vertex, every bond endpoint, every ConformationAtom, every
> spatial neighbour list is indexed into this array. Patching indices
> is fragile and error-prone. The correct mechanism is a clean
> rebuild:
>
> 1. Take the source protein's heavy-atom data (symbolic: residue
>    types, atom names within residues, sequence)
> 2. Apply new protonation decisions (from PROPKA, KaML, or explicit)
> 3. Run tleap (or equivalent) to produce the new atom list with
>    hydrogens placed according to the new protonation
> 4. Build a new Protein with fresh indices
> 5. Create conformations on the new Protein using heavy-atom
>    position transfer (symbolic matching by residue + atom name)
> 6. Run the full extraction pipeline on the new conformations
>
> The heavy-atom skeleton survives because symbolic identity (residue
> 17, atom CA) maps to a position in both old and new index spaces.
> The name-within-residue lookup bridges them. This is a dictionary
> lookup, not geometric matching.
>
> ### What exists now
>
> - **MutationDeltaResult**: two separate Proteins (WT and mutant),
>   atom matching by position + element, delta tensors computed.
>   Attached to WT conformation as a ConformationResult.
> - **GromacsEnsembleLoader**: builds a fresh Protein from TPR
>   topology with CHARMM naming → canonical translation.
> - **OrcaRunLoader**: builds from prmtop + XYZ with AMBER naming.
> - **ProtonationState, Protonator interface, NamingRegistry**: the
>   infrastructure for making protonation decisions and translating
>   between naming universes exists and is tested.
>
> ### What the builder will do (TBD, 2-3 sessions)
>
> A ProteinBuilder that unifies the loading paths: takes source data
> (from any loader), takes protonation decisions (from any protonator),
> runs tleap to produce the protonated atom list, and constructs the
> Protein through the shared FinalizeConstruction pipeline. The three
> existing loaders (PdbFileReader, OrcaRunLoader, GromacsEnsembleLoader)
> become Layer 1 sources feeding into the builder.
>
> ### Why copy was wrong
>
> The original design called for copying a Protein and applying a new
> ProtonationState. This required a copy constructor that could
> add/remove atoms, renumber every index, and selectively invalidate
> computed results. The two indexing systems (symbolic identity and
> geometric flat array) are fused at FinalizeConstruction time. There
> is no safe way to patch one without rebuilding the other. A clean
> rebuild from symbolic data is simpler, correct, and already proven
> by the three existing loaders.
>
> See spec/PROTONATION_DESIGN_HISTORY.md for the full analysis.

---

## 4. Enforcement: Logging and Tracing (replaces "Framework Stores"
   section lines 846-877 and "UDP Logging" section lines 971-1015)

**Note**: The original described a StoreRingContribution() API with
per-property UDP logging that doesn't exist. What DOES exist: the
KernelFilterSet logs every rejection with specific context (atom,
source, distance, filter name, reason). OperationLog provides
operation-level tracing. The ConformationResult singleton guarantee
provides the structural tracing. The proposed text describes what
works and where to extend it.

### Proposed text:

> ## Enforcement: Structural Guarantees and Tracing
>
> ### The singleton guarantee
>
> Each ConformationResult type is a singleton on a ProteinConformation.
> BiotSavartResult is computed once, attached once, and never replaced.
> Each ConformationAtom field has exactly one writer, determined by
> which ConformationResult type owns it. The dependency graph ensures
> results attach in the correct order. AttachResult checks:
>
> - Not null
> - Not already attached (singleton violation → logged error, rejected)
> - All dependencies present (missing dependency → logged error, rejected)
>
> This is the structural guarantee that prevents overwrites and
> ordering errors. It replaces per-property overwrite detection with
> a stronger invariant: the scenario cannot arise.
>
> ### Filter rejection logging (per-evaluation)
>
> Every KernelFilterSet logs per-evaluation rejections when enabled:
>
> ```
> RingBondedExclusionFilter: atom=47 is bonded to ring 3
> DipolarNearFieldFilter: atom=123 distance=0.82A < threshold=1.15A
>   (source_extent=2.30A)
> SelfSourceFilter: atom=15 IS source_atom=15
> ```
>
> Each rejection record includes: filter name, atom index, source
> indices, distance, source extent, and a human-readable reason.
> Per-filter rejection counts are reported at computation end.
> This is the "what happened to atom 42" diagnostic.
>
> Enable per-rejection logging for single-protein diagnostics:
> `filter_set.SetLogRejections(true)`. Disable for batch performance.
>
> ### Operation logging (per-result)
>
> OperationLog (OperationLog.h) provides structured JSON logging
> over UDP or stderr with channel bitmask control. Every result
> attachment, external tool invocation, and computation boundary
> emits a log entry. Scoped logging tracks elapsed time:
>
> ```cpp
> OperationLog::Scope scope("BiotSavartResult::Compute",
>     "protein=" + protein_name + " atoms=" + std::to_string(n));
> // ... computation ...
> // destructor emits END with elapsed time
> ```
>
> 22 log channels cover bond classification, ring detection,
> calculator dispatch, file I/O, DSSP, APBS, xTB, charges,
> atom mapping, protonation, conformation management, and result
> attachment.
>
> ### Where to extend
>
> The filter rejection logging demonstrates the right pattern:
> specific context (which atom, which source, what distance, why
> rejected) at the evaluation level. The same pattern should be
> applied to calculator algorithm steps — not every field write,
> but the key decision points:
>
> - Which rings were evaluated for this atom (and which excluded)
> - Which bonds contributed to the McConnell sum
> - What the accumulated T0 and T2 are after this ring/bond
>
> This would be opt-in per-calculator logging (like filter rejection
> logging), not a framework-level store intercept.

---

## 5. Planned Result Types (replaces the entries for
   ParameterCorrectionResult, FeatureExtractionResult,
   PredictionResult in the dependency graph, lines 369-371)

**Note**: These are the design for the next phase. They follow the
same ConformationResult pattern as all existing results. They should
be clearly marked as planned rather than implemented.

### Proposed text (for the dependency graph):

> ```
> --- Planned (not yet implemented) ---
> FeatureExtractionResult    requires: all physics results
> ParameterCorrectionResult  requires: FeatureExtractionResult
> PredictionResult           requires: ParameterCorrectionResult
> ```
>
> These follow the same ConformationResult singleton pattern as all
> existing results. Each reads from what prior results stored on the
> conformation, computes its own contribution, writes to its own
> fields on ConformationAtom, and attaches as a singleton.
>
> - **FeatureExtractionResult**: reads the shielding_contribution
>   fields from all 8 calculators plus enrichment, DSSP, charges,
>   and graph features. Organises them into the 189-feature vector
>   (137x0e + 6x1e + 10x1o + 36x2e) as queryable runtime state.
>   WriteFeatures already serializes this data to NPY; this result
>   makes it available for in-process queries by downstream results.
>
> - **ParameterCorrectionResult**: trained e3nn model that predicts
>   corrected parameters for classical calculators from environment
>   features. See CALCULATOR_PARAMETER_API.md for the 93-parameter
>   interface. Provides residual tracking that updates as classical
>   results attach.
>
> - **PredictionResult**: final ML inference producing per-atom T0,
>   T2 predictions with confidence and heuristic tier (REPORT/PASS/
>   SILENT).
>
> The ConformationAtom fields for these results already exist as
> default-initialised placeholders (ConformationAtom.h).

---

## 6. ProteinConformation Hierarchy (replaces lines 94-112)

**Note**: The original listed intermediate grouping classes
(ExperimentalConformation, ComputedConformation) and types that
don't exist (NMRConformation, MinimisedConformation). The flat
hierarchy was the right call — the intermediate classes added
no functionality.

### Proposed text:

> **Class hierarchy:**
>
> ```
> ProteinConformation (base -- NOT abstract, fully functional)
> +-- CrystalConformation
> |     resolution_angstroms, r_factor, temperature_kelvin, pdb_id
> +-- PredictionConformation
> |     method, confidence
> +-- MDFrameConformation
> |     walker, time_picoseconds, boltzmann_weight, rmsd_nm, rg_nm
> +-- DerivedConformation
>       derivation_description
> ```
>
> The hierarchy is deliberately flat. The base class does the heavy
> lifting: holds positions (const), provides ConformationResult
> attachment and template access, holds the ConformationAtom wide
> table. Subclasses add only provenance metadata. No intermediate
> grouping classes (Experimental, Computed) — they were considered
> and added no value.
>
> Future types (NMRConformation, MinimisedConformation) can be added
> as direct subclasses of ProteinConformation when needed, following
> the same pattern: metadata only, inheriting all base class machinery.
>
> **Primary access**: `protein.Conformation()` returns the first
> conformation as ProteinConformation& regardless of subtype. This
> is what calculators and tests use. Typed accessors (CrystalConf(),
> MDFrameAt(), PredictionAt()) exist for metadata consumers.

---

These six edits cover the Category 2 judgement calls. Each preserves
the original principle while correcting the mechanism to match reality.
