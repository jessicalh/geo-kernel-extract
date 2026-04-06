# Design Brief — DRAFT 1 of N

This is the first of several drafts. It captures constraints and process,
not solutions. Solutions come from the categorical analysis and adversarial
review that follow.

## Principles (non-negotiable)

1. **One library.** No v1/v2. No shadow types. No adapters.
2. **Categories are typed.** If you categorise an object, that category
   is embedded in the object's representation, part of a recognised class
   hierarchy, and DRY. Not strings. Not inline conditionals.
3. **The protein model is the single source of truth.** Calculators and
   feature extractors query the protein. They do not build parallel
   representations.
4. **Third-party tools are authorities.** DSSP for secondary structure.
   APBS for solvation. Force fields for charges. We represent their
   answers; we do not reimplement their logic.
5. **Sign conventions, normalization, and ring types are chosen once,
   documented with worked examples, and used everywhere.**
6. **Features are subclassed, not configured.** Adding a feature is
   writing a class, not editing a manifest. The base class enforces
   name, irrep, and compute method.
7. **The extraction uses the same protein model as the calculators.**
   Same Protein. Same Conformation. Same AromaticRing. One truth.

## Process (three stages before code)

### Stage 1: Categorical analysis

Enumerate every categorical decision in the v1 EquivariantExtractor
(1740 lines, ~16 methods). For each decision:

- What is the decision? (e.g., "is this atom an amide hydrogen?")
- What object does it belong to? (Atom, Bond, Residue, Ring, Interaction)
- What structural fact does it derive from? (bond connectivity, element, etc.)
- Is it already on a v2 object?
- Is it duplicated?

This analysis requires multiple adversarial agents because a single pass
will miss implicit categorisations (e.g., "if nearest ring has nitrogen"
is really a ring type query, not a feature decision).

Also enumerate categorical decisions we WANT but don't have:
- Through-bond distance (BFS from ring atoms)
- Protonation-variant-specific ring properties
- Insertion code handling
- Multi-chain sequence separation

### Stage 2: Object model breakout (shitty first draft)

Express the FULL set of categories as a class hierarchy. Not simplified.
Not implemented. Just the structure:

- What objects exist?
- What does each object know about itself?
- What queries can you ask?
- Where do categorical properties live?

This will be overdesigned. That's intentional. You need the full picture
before you can simplify. Subclassing a ring type vs adding a property:
that's a design judgement that can only be made when you see ALL the
ring-type-dependent decisions in one place.

Include: everything the current extractor does.
Include: everything we WANT it to do (mutant isolation features,
through-bond proxy, ensemble statistics, confidence field).

Then simplify. Remove what's redundant. Flatten what's over-hierarchical.
The simplified version is the design.

### Stage 3: Operations flow

Given the object model, define:

- What happens when a PDB arrives?
- What happens per conformation?
- What gets computed, in what order, with what dependencies?
- What serializes, what stays live?
- What does the viewer read? What does the upstream model read?
- What does the mutant comparison query?

The operations flow is NOT a pipeline configuration. It's a loop
with well-typed interfaces between stages.

## What we have now (the ingredients)

### From v2 (good object model, incomplete physics):
- Protein / Atom / Residue / Conformation / Environment
- AromaticRing / CovalentBond / RingType table
- 9 Calculator classes with clean interfaces
- ProtonationState (not fully wired to ring types)
- DSSP, APBS, OpenBabel integration
- Qt/VTK viewer with REST interface and working isosurfaces

### From v1 (good feature engineering, bad object model):
- AtomSite: distance-sorted rings, parent atom, structural context
- MolecularGraph: BFS, hybridization, electronegativity
- EquivariantExtractor: 136 features with known physics meaning
- Category-decomposed McConnell (backbone CO/CN/sidechain)
- Aromatic sidechain Coulomb field (Buckingham isolation)
- Nearest-ring cylindrical coordinate features

### From this session (architecture insights):
- Gated physics correction (ML modulates classical baseline)
- Heteroscedastic confidence (learned per-atom uncertainty)
- Three-tier heuristic (REPORT/PASS/SILENT)
- Through-bond distance as a feature (BFS from MolecularGraph)
- Curriculum training (sign → magnitude → refinement)
- Asymmetric nodal loss
- Per-ring-type Gaussian density (learned spatial envelope)

### Wanted but not yet built:
- Mutant isolation features (charge-flip, ring modification, salt-bridge)
- Ensemble statistics across conformations
- Per-calculator decomposed output for upstream model
- Insertion code handling
- Multi-chain correctness throughout
- Automated physics verification tests

## What this brief does NOT contain

- Class diagrams (premature before categorical analysis)
- Code (premature before object model)
- Feature list (premature before "what we want" is complete)
- Architecture choices (premature before adversarial review)

These come from the stages above, not from this document.
