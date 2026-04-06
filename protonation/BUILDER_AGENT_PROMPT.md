# Agent Task: Understand the Protein Builder Gap

**Output**: A written analysis document. No code. No lists of files
to change. No implementation plan.

**You succeed if**: your analysis shows you understood the object
model, the original design intent, where reality diverged, and
what constraints and freedoms apply to fixing it.  You succeed by genuine engagement as a thoughtful C++ developer who understands why Stroustrop built C++ as a modeling language, and how a simple clean design can reflect physics and model

**You fail if**: you propose a generic factory pattern, add
unnecessary abstractions, list code changes, or demonstrate that
you skimmed the docs instead of reading them.

---

## Your reading list (in order, all of it, not the first 200 lines)

1. `spec/CONSTITUTION.md` — the supreme constraint. Read the object
   model section, the copy-and-modify section, the protonation tools
   section, the static vs dynamic properties section. All of it.
2. `OBJECT_MODEL.md` — the concrete types. Every class, property, type.
3. `PATTERNS.md` — how things are built, what gets rejected.
4. `spec/PROTEIN_BUILDER_SPEC.md` — the design we arrived at through
   a long conversation. Read it AFTER the constitution and object model
   so you understand what it's correcting.
5. Source code (read fully, not grep):
   - `src/Protein.h` and `src/Protein.cpp`
   - `src/Atom.h`
   - `src/Residue.h`
   - `src/Ring.h` and `src/Ring.cpp`
   - `src/ProteinBuildContext.h`
   - `src/ProtonationState.h` and `src/ProtonationState.cpp`
   - `src/Protonator.h`
   - `src/PropkaProtonator.h` and `src/PropkaProtonator.cpp`
   - `src/ProtonationDetectionResult.h` and `src/ProtonationDetectionResult.cpp`
   - `src/AminoAcidType.h` and `src/AminoAcidType.cpp`
   - `src/NamingRegistry.h` and `src/NamingRegistry.cpp`
   - `src/ChargeSource.h`
   - `src/PdbFileReader.h` and `src/PdbFileReader.cpp`
   - `src/OrcaRunLoader.h` and `src/OrcaRunLoader.cpp`
   - `src/ProteinConformation.h`
   - `src/ConformationAtom.h`
   - `src/Pipeline.h` and `src/Pipeline.cpp`
   - `tests/test_protonation_pipeline.cpp`

---

## The problem

In structural biology, protonation is always an issue when modeling and working with proteins.  Most PDBs do not have it.  Model results do.
Systems like this one, which do a  one-way analysis of a protein's physics properties, must be able to analyise the protein based on protonation
via various other modeling systems -- propka or KaML for example -- and then recaclculate the physics properties based on that change.  But since
modeling happens in an instantited set of structures with implicit geometry and physics, you either create a generality we do not want to support
or become locked to a single path.  And while a single path is fine, it also fails to collect or make orgothonal, as methods are added to proponate,
the properties like naming which must be tracked.

The core vision of this sytem is a protein object, which contains geometry invariant properties, and conformations, which contain the 
real geometry, whether of a pdb or MD frame poses.  It was considered unlikely we would do CpHMD in this project, so it was OK if protonation
stuck to the protein object -- but there must be a pattern to support deprotonation and reprotonation, with good C++ semantics.  And ideally -- even if
things turned from 1 protein instance with 1 protonation enforced and 100 poses to 100 protein instances with different protonations enforced, well
that was OK too if we ended up handling CpMHD results.

The idea was that if you wanted to make a new protonation state, you would use a builder. The builder would take the old protein, apply a ctor
that gave you back a hydrogen free copy, create a protonation state using kaml or proka or whatever, even a poses' old protonation state (classed -- 
the fact this was not clearly defined may have contrinuted to the mess) in the state according to the how, with the state holding receipts (not
hydrogens).  You would also supply the old confomration (or conformations) and get them back with your geometry but absolutely no calculation results, and any
geometry changes required for the hydrogens -- you have to run all the calculatins again but yu have your start state.

Instead, as happens, we have what we have. We are asking you to take a very very deep dive. If we try to fix this, what will it break?  Do not ask what is the
first case I find and how to fix it and done.  Ask how would it be fixed and what would _that_ break.  We are specifying a fix pass.  It must give idential
results at the highest level of test (see all the tests, there's no finessing this!).
Unfortunately it was simply necessary to get the basics working first,
and along the way the vision got locked into a single load-and-go path
where the protein arrives already protonated and that's that.

This is a physics thesis project that computes geometric kernels —
the spatial shape of classical electromagnetic effects at every atom
in a protein. 8 calculators, 274 tests, 723 proteins batch-validated.
The calculators work. The thesis requires what-if NMR protonation:
run the calculators, change protonation (different pH, different tool),
run again, compare. Without a reprotonation path, the project cannot
do this.

## The original design vision

Read the constitution's copy-and-modify section and protonation tools
section. The builder concept described above was the intent. The
constitution also specified ProtonationState as "one per
ProteinConformation" — that turned out to be wrong (protonation
determines the atom list, which is per Protein, not per conformation).
The constitution has other stale text around copy-and-modify. Don't
try to make the stale text work. Understand what it was TRYING to say
in light of the problem description above.

## What actually got built (and why)

Two loaders got built instead of a builder:

- `LoadProtein` (PdbFileReader): parses PDB, builds atoms from
  whatever's in the file, finalises (bonds, rings), creates one
  conformation. No protonation step.

- `LoadOrcaRun` (OrcaRunLoader): reads an AMBER prmtop for the
  protonated atom list, reads XYZ for positions, sets protonation
  variant indices from AMBER residue labels. Protonation happened
  externally (tleap) before our code ever sees it.

Both loaders do part of what the builder was supposed to do. Neither
offers protonation as a step you can request.

The protonation infrastructure exists: ProtonationState (value type
with decisions), Protonator interface, PROPKA and KaML implementations,
NamingRegistry (translates between tool naming universes), per-variant
metadata on AminoAcidType (name, charge, registry_key). All tested.
But nothing connects "I have protonation decisions" to "build me a
Protein with the right H atoms."

## Where the original design didn't survive

1. **"One ProtonationState per ProteinConformation"** — wrong.
   Protonation determines which atoms exist. Different protonation =
   different atom list = different Protein. A conformation's atom
   vector is parallel to its Protein's atom list. Two conformations
   on the same Protein must have the same atoms. Protonation is per
   Protein, not per conformation.

2. **Copy-and-modify** — the concept was right (keep the old protein,
   produce a new one with different protonation) but the mechanism
   was wrong. A copy constructor that "applies" a ProtonationState
   would need to add/remove H atoms, which changes the atom count,
   which invalidates every index into the atom list. The real
   mechanism is: the builder takes heavy-atom data from the existing
   protein and builds a new one from scratch. Not a copy — a rebuild.

3. **Protein is non-copyable** — `= delete` on copy and move. This
   was done for pointer safety (conformations hold raw Protein*
   back-pointers). The design intended a Copy() factory method but
   it was never built because the early agents deleted copy/move
   entirely instead of providing a proper factory.

These aren't bugs. They're the result of real constraints changing
the design during implementation. The atom list IS the protein's
identity — that's correct and approved. It just means protonation
can't be a state change on an existing Protein.

## What we DON'T care about

- That the original copy-and-modify language is wrong. It's being
  corrected. Don't try to make copy-and-modify work.
- That ProtonationVariant doesn't have per-variant atom lists. tleap
  knows what atoms each variant has. We don't need to replicate that.
- That Protein is non-copyable. The builder doesn't need to copy a
  Protein. It takes heavy-atom positions and sequence data.

## What your analysis must cover

### 1. Show you understand the object model

Explain in your own words: what is a Protein, what does it own, why
does it own atoms including hydrogens, what is a ProteinConformation's
relationship to its Protein, why are positions const, why is Protein
non-movable.

### 2. Show you understand the protonation lifecycle

Trace the path from "I have a PDB file" to "I have a protonated
Protein with a conformation ready for calculators." What happens at
each step. Where does protonation enter. What does NamingRegistry do.
What does ProtonationDetectionResult do vs what do the Protonators do.

### 3. Show you understand the builder gap

What does the builder need to do that neither loader does? Be
specific about the protonation step — what tool gets called, what
data flows in, what data flows out, where NamingRegistry translates.

### 4. Identify the real constraints

What CANNOT change without breaking things? Be specific. Think about:
atom indices, conformation atom vectors, ring vertex indices, bond
indices, the spatial index, the filter framework, all 8 calculators.

### 5. Identify the real freedoms

What CAN change safely? What is genuinely not coupled to anything
else? Think about: the loader functions, ProteinBuildContext, how
Proteins are created, what happens before FinalizeConstruction.

### 6. What would the builder need from tleap (or equivalent)?

The force field tool places H atoms. What does it need as input?
What does it produce? How does NamingRegistry fit? Be specific
about the data flow across the tool boundary.

### 7. What happens to existing tests?

274 tests pass. Which ones create Proteins? Through which path?
Would routing them through a builder change their behavior? Don't
list every test — identify the PATTERNS of how tests create proteins
and what the impact categories are.

---

## Rules

- Read ALL the files listed above. Not summaries. Not first 200
  lines. Read them.
- Do not propose code changes.
- Do not propose new classes, interfaces, or abstractions.
- Do not list files that would need to change.
- Write your analysis as prose, not bullet lists. Show your
  reasoning.
- If you find contradictions between the constitution and the code,
  say so. The code wins — we know the constitution has stale text.
- If something seems wrong in the code, check if a test exercises
  it before calling it wrong.
- Your document should be 2000-4000 words. Shorter means you
  skimmed. Longer means you're padding.
