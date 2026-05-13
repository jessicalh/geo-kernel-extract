# Protonation Design: Analysis and Decision Record

This directory contains the artifacts of a deep design analysis
conducted on 2026-04-02 across a single long session. Read the
design history first, then the analysis chain if you need the
reasoning.

## Start here

- **DESIGN_HISTORY.md** — the full story: original vision, what got
  built, where it diverged, the two indexing systems (symbolic vs
  geometric), findings, and design conclusions. This is the
  authoritative document.

## The decision

**Protonate on load.** The system is a oneshot analysis engine: load,
protonate, calculate, export. If you want different protonation, load
again from source with different protonation decisions. This is not a
limitation — it's what the physics requires. Changing protonation
changes the atom list, which invalidates every geometric index.

**Clean up the divergence.** Two loaders (LoadProtein, LoadOrcaRun)
descended from an unbuilt builder concept. They share a pattern but
don't share code. The constitution describes copy-and-modify
semantics that were never implemented. ProtonationState exists as a
decision value type but nothing consumes it for construction.

**The three real architecture issues to fix:**

1. ProtonationDetectionResult writes to Protein via const_cast
   through a conformation-level result. Identity data should be
   set at construction, not by a retroactive annotation.

2. Everything loads as CrystalConformation, even AlphaFold
   predictions. Loses metadata. The type system lies. Fix required
   before MD frame loading (which needs MDFrameConformation).

3. CoulombResult is the only calculator that bypasses
   KernelFilterSet. Inline self-exclusion is invisible to filter
   rejection logging.

## Analysis chain (for the reasoning)

Read in order if you need to understand HOW we arrived at the
decision:

1. `BUILDER_AGENT_PROMPT.md` — the task given to agents (includes
   the human's description of the original vision)
2. `BUILDER_ANALYSIS.md` — round 1 analysis (3400 words)
3. `BUILDER_ANALYSIS_CRITIQUE.md` — what round 1 missed
4. `BUILDER_ANALYSIS_ROUND2.md` — supplementary analysis (2500 words)
5. `PROTEIN_BUILDER_SPEC.md` — the builder design that emerged

## Key insight: two indexing systems

A Protein has symbolic identity (residue sequence, atom names, ring
types) and geometric identity (flat atom index 0..N-1). Every ring
vertex, bond endpoint, conformation atom, spatial neighbour, and
calculator result is indexed into the geometric array. Changing
protonation changes which atoms exist, which changes the geometric
array, which invalidates everything indexed into it. The symbolic
identity survives. The geometric identity is rebuilt.

This is why reprotonation is a reload, not a patch. The loader
(or future builder) fuses symbolic identity with a new geometric
index space each time.

## What's NOT in this directory

A fragility audit was attempted and reverted. It analysed enterprise
data-flow risks in a one-way pipeline that doesn't have them. The
system's ConformationResult pattern is clean: each result owns its
data, writes to its conformation once, no hidden globals, no
cross-contamination. The reverted analysis would have been misleading
to future agents.

## Related code

- `src/CalculationRunner.h` — the three allowed calculation sequences
  (PDB only, single DFT, WT+ALA comparison). This is the "run stuff"
  entry point that ensures correct ordering.
- `src/Pipeline.h` — lower-level convenience that CalculationRunner
  may eventually replace (Pipeline was written earlier in the session
  before the CalculationRunner design crystallised).
- `src/ProtonationState.h` — the decision value type
- `src/PropkaProtonator.h`, `src/KamlProtonator.h` — protonation tools
- `src/NamingRegistry.h` — translation between tool naming universes
- `src/PdbFileReader.cpp`, `src/OrcaRunLoader.cpp` — the two loaders
