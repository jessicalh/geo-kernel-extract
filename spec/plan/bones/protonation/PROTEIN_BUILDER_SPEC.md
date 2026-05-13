# ProteinBuilder: Spec Change

**Status**: Design. No code changes yet. Evaluation required first.

**Why this matters**: The thesis requires what-if NMR protonation — run
the calculators, change protonation, run again, compare. Without a
reprotonation path, the project cannot do this. The pieces exist
(ProtonationState, PROPKA, KaML, NamingRegistry, two working loaders)
but the unifying factory that makes reprotonation a single operation
does not.

---

## The concept

ProteinBuilder is the single point of entry for getting a working
Protein. Every Protein in the system comes from the builder. There
is no other way to get one.

### Inputs

- **Sequence + heavy atom positions.** From a PDB file, from an
  existing Protein's heavy atoms, from an XYZ file — the builder
  doesn't care where they come from. It needs the residue sequence
  and the heavy atom coordinates.

- **Protonation instructions.** One of:
  - "leave as-is" — the input already has H atoms, keep them.
    Detect protonation from what's present. Fail if unprotonated
    and calculators need H.
  - "protonate with [tool] at pH [X]" — strip any existing H,
    run PROPKA/KaML to get decisions, run the force field tool
    (tleap, pdb2gmx) to place H atoms with correct geometry.
  - "apply this ProtonationState" — strip H, add H per these
    specific decisions. For when you already have decisions and
    don't need to re-predict.

- **Force field context.** ff14SB, CHARMM36m, etc. Determines how
  H atoms are placed and what the atom names are. NamingRegistry
  translates at the tool boundary.

### What it does

1. If protonation requested: strip any existing H from the atom list.
2. If a protonation tool is requested: run it (PROPKA/KaML) to get
   decisions. If a ProtonationState was provided directly, use that.
3. Run the force field tool (tleap, pdb2gmx) to place H atoms with
   correct geometry. The tool gets residue names translated by
   NamingRegistry into its naming universe.
4. Build Protein: atoms (heavy + H), residues (with variant indices
   set from protonation decisions), backbone cache, ring detection
   (HIS variant determines ring type), bond detection.
5. Create first ProteinConformation from positions (heavy atom
   positions unchanged, H positions from the force field tool).
6. Set ProteinBuildContext recording: source, protonation tool,
   pH, force field, timestamp.

### What it produces

- `unique_ptr<Protein>` — complete, one conformation, ready for
  `RunClassicalCalculators()`.
- `ProtonationState` — the decisions that were made. Value type.
  Travels separately from the Protein. Available for comparison,
  for passing to a different tool, for recording.

### The three use cases

**First load from PDB (most common):**
```
builder.SetSource(pdb_path);
builder.SetProtonation(Protonate::WithTool, propka, 7.0);
builder.SetForceField(ForceField::Amber_ff14SB);
auto result = builder.Build();
// result.protein: complete, protonated, one conformation
// result.protonation_state: PROPKA decisions at pH 7
```

**Reprotonate an existing protein (the thesis use case):**
```
// protein_a exists, was protonated with PROPKA at pH 7
// I want to try KaML at pH 5

builder.SetSource(protein_a);  // takes sequence + heavy atom positions
builder.SetProtonation(Protonate::WithTool, kaml, 5.0);
builder.SetForceField(ForceField::Amber_ff14SB);
auto result = builder.Build();
// result.protein: NEW Protein, different H atoms, different rings
// protein_a still exists with all its conformations and results
```

**Load pre-protonated structure (ORCA/fleet data):**
```
builder.SetSource(prmtop_path, xyz_path);  // already has H
builder.SetProtonation(Protonate::AsIs);   // keep what's there
auto result = builder.Build();
// result.protein: loaded from prmtop+XYZ, protonation detected
```

---

## What exists today

Two loaders that each do PART of what the builder should do:

### LoadProtein (PdbFileReader.cpp)
- Parses PDB via cif++
- Builds atoms, residues from whatever's in the file
- Calls FinalizeConstruction (backbone, rings, bonds)
- Creates CrystalConformation
- Sets ProteinBuildContext (minimal)
- Does NOT protonate. Does NOT strip H. Does NOT run any tool.

### LoadOrcaRun (OrcaRunLoader.cpp)
- Reads prmtop for atom list (protonated, from tleap)
- Reads XYZ for positions
- Sets protonation_variant_index from AMBER residue labels (HID/HIE/etc)
- Calls FinalizeConstruction
- Creates CrystalConformation
- Sets ProteinBuildContext (tleap, ff14SB, prmtop path)
- Does NOT protonate. Loads an already-protonated structure.

### ProtonationDetectionResult
- Runs after loading. Reads H atoms present, sets variant_index
  on Residue via const_cast. Retroactive annotation.
- NOT a builder step. Confirms what the loader already built.

### PROPKA/KaML Protonators
- Take loaded Protein + Conformation + pH
- Return ProtonationState (decisions)
- Do NOT modify the protein. Do NOT add/remove atoms.

### NamingRegistry
- Translates residue names between tool contexts
  (HIS/HID/HIE/HIP ↔ HSD/HSE/HSP, ASP/ASH ↔ ASPP, etc.)
- Translates atom names (H ↔ HN, HB2 ↔ HB1)
- The translation layer between our canonical model and each tool

### ProtonationVariant (on AminoAcidType)
- Name, description, formal charge, registry_key for each variant
- The registry_key is how NamingRegistry looks up tool-specific names
- No atom list delta — doesn't say which H atoms the variant has.
  That knowledge lives in the force field tool (tleap knows what
  HID has), not in our data tables. This is probably correct —
  we don't need to replicate tleap's knowledge.

---

## What changes in the design

### Constitution changes

1. **Remove "one protonation state per ProteinConformation".**
   Protonation determines the atom list. It's per Protein.
   Different protonation = new Protein.

2. **Replace copy-and-modify with builder pattern.** The old text
   says "copy protein, apply new ProtonationState." The real
   pattern is: builder takes existing protein's heavy atoms as
   source, builds a new Protein with different protonation.
   No copy constructor needed. The builder IS the factory.

3. **ProteinBuilder is the single entry point.** Every Protein
   comes from the builder. LoadProtein and LoadOrcaRun become
   convenience wrappers that configure the builder and call Build().

4. **ProtonationState is a builder input/output, not a Protein
   field.** It's the decisions. The Protein records what happened
   in ProteinBuildContext. The ProtonationState is returned
   alongside the Protein for anyone who needs the structured
   decisions.

### What does NOT change

- Protein owns atoms including H. Correct. Approved.
- Protein is non-movable (conformations hold raw back-pointers).
  Correct. Keep.
- ProteinConformation holds positions parallel to Protein's atoms.
  Correct.
- ConformationResult singleton pattern. All of it. Untouched.
- All 8 calculators. All 274 tests. All batch validation.
- NamingRegistry. ProtonationVariant. AminoAcidType tables.
- ProtonationState as a value type with decisions.
- PROPKA and KaML protonators.
- The Pipeline (RunClassicalCalculators).

### What might be removable (evaluate first)

- `ProtonationDetectionResult` — if the builder always sets
  variant_index during construction, retroactive detection may
  be redundant. BUT: it's useful for validating what a loader
  built, and for analysing pre-protonated structures where the
  builder wasn't involved. Keep unless it causes confusion.
- The const_cast in ProtonationDetectionResult — if variant_index
  is set by the builder at construction, this shouldn't be needed.
  But if detection is kept as a validation tool, it still needs
  write access to Residue.
- `Protein(const Protein&) = delete` — the builder doesn't need
  a copy constructor. It takes heavy atom positions and sequence,
  not a Protein reference. The Protein stays non-copyable.

---

## Evaluation needed before implementation

The agent's job is to evaluate what building ProteinBuilder breaks
or changes. **Change nothing.** Report only.

### Questions to answer

1. **What calls LoadProtein() today?** Every call site needs to
   go through the builder instead, or LoadProtein becomes a
   builder convenience wrapper. List them all.

2. **What calls LoadOrcaRun() today?** Same analysis.

3. **What reads Residue::protonation_variant_index?** This field
   gets set by detection or by the loader. Under the builder, it
   gets set at build time. Does anything depend on it being set
   AFTER construction?

4. **What depends on Protein being default-constructible?** The
   loaders create empty Proteins and fill them. The builder would
   do the same internally. But if test code creates empty Proteins,
   does that break?

5. **What is the test impact?** Run all 274 tests. Identify which
   ones create Proteins directly (not through loaders) and would
   need updating if the builder becomes the only entry point.

6. **tleap integration.** We call tleap for ORCA prep (tleap scripts
   exist in the data). How hard is it to call tleap programmatically
   with a specific ProtonationState? The NamingRegistry already
   translates to AMBER names. What's missing?

7. **Can the two existing loaders become builder configurations?**
   LoadProtein = builder with PDB source + as-is protonation.
   LoadOrcaRun = builder with prmtop+XYZ source + as-is protonation.
   What breaks if we restructure them this way?

### What NOT to do

- Do not add complexity. The builder is ONE class with ONE Build()
  method. Not a strategy pattern. Not a template. Not a plugin
  architecture.
- Do not create intermediate types (UnprotonatedProtein,
  HeavyAtomSkeleton, etc.). The builder takes inputs and produces
  a Protein. Internally it can have steps. Externally it's one call.
- Do not touch calculator code. The calculators work. They will
  continue to receive complete Proteins from the builder exactly
  as they receive them from the loaders today.
- Do not move NamingRegistry, ProtonationVariant, or AminoAcidType.
  They're fine where they are.

---

## For the UI

The viewer's surface for reprotonation:

```cpp
// In Pipeline.h or a new ProtonationPipeline.h:

struct ReprotonationResult {
    std::unique_ptr<Protein> protein;
    ProtonationState protonation_state;
    std::string error;
    bool Ok() const { return protein != nullptr && error.empty(); }
};

// Reprotonate an existing protein at a new pH.
// Takes heavy atom positions from the specified conformation.
// Returns a NEW Protein ready for RunClassicalCalculators.
ReprotonationResult Reprotonate(
    const Protein& source,
    const ProteinConformation& source_conf,
    Protonator& tool,
    double pH,
    ForceField ff = ForceField::Amber_ff14SB);
```

The viewer calls `Reprotonate()`, gets back a new Protein, runs
Pipeline on it, displays results alongside the original. The viewer
never sees the builder internals. It doesn't know about tleap or
NamingRegistry or variant indices.
