---
name: Rewrite decisions from 2026-03-31 session
description: Clean-room rewrite decisions including mutant strategy, feature design, and architecture principles. This is the plan.
type: project
---

## Core Decision
One final rewrite. No v1/v2. No dead code. No convenient fictions. Last rewrite before thesis submission unless something goes fundamentally wrong.

## Two Pillars

### Pillar 1: What IS a protein
Atoms, residues, chains, conformations, protonation states, aromatic rings, covalent bonds, secondary structure. Geometry and identity. Third party tools (dssp, OpenBabel, APBS) are authorities where they apply. Our object model represents their answers.

### Pillar 2: What do the calculators COMPUTE
Each calculator takes a point and a protein and returns a tensor. Equations from textbooks. Sign conventions chosen once. Outputs in irrep form.

## Feature Design (two tiers)

### Tier 1: Raw physics (for e3nn)
Raw calculator outputs in natural irrep form, per source. No projections, no ratios, no pre-computed combinations. Let the equivariant TP layers discover which products matter. This is where unexpected physics emerges.

### Tier 2: Engineered features (for bilinear)
Projections, ratios, nearest-ring geometry, McConnell decomposition. The bilinear model needs these because it can't construct tensor products freely.

The upstream model gets both tiers plus predictions.

## Mutant Strategy (NEW - for isolation training)

### Current: 600 more WT→ALA proteins running, capacity for 800 total

### Three new mutant categories (~200 each):

1. **Charge-flip (ASP→ASN, GLU→GLN)**: Changes electrostatic without changing aromatic ring. Isolates Buckingham from ring current. The delta between WT-ASP and WT-ASN is pure electrostatic correction.

2. **Ring modification**: Changes ring electronic structure without changing backbone geometry. Options from Gemini: remove ring oxygens (TYR→PHE equivalent), remove paired ring from TRP (indole→pyrrole or indole→benzene). Tests whether the model correctly separates intensity I from geometry G in I×G.

3. **Salt-bridge breakers (LYS→MET, ARG→LEU)**: Changes solvation environment by removing charged residues near rings. Isolates the APBS/solvation contribution. Tests whether solvated vs vacuum E-field matters.

### Why mutants matter
Each mutant type gives the model a controlled experiment. If charge-flip mutants show the model's learned electrostatic correction matches the DFT delta, we've VALIDATED the Buckingham decomposition experimentally. That's a thesis result.

### Tool recommendations needed
- DFT tool for running mutant calculations (ORCA r2SCAN, same as current pipeline)
- Structure preparation: need to validate mutant PDB construction (which tool? PyMOL mutagenesis? Modeller? pdb-tools?)
- Sanity checks: are the mutant structures physically reasonable before running DFT?

## Architecture Principles
- No string-driven identity in physics code
- Protonation state is a first-class object, detected from structure, authorities consulted
- Ring type determined by the type system + protonation state, not by grepping atom names
- Every calculator's equation documented with literature reference
- One sign convention, stated once, used everywhere
- Features documented with equations, not just names

## What went wrong with v1/v2
- v1 had the feature engineering but old physics conventions
- v2 had corrected physics but couldn't extract features
- I (Opus) kept using the v1 binary while developing v2, hiding the divergence
- Neither codebase was the truth. Both were partly wrong.
- The byte comparison exposed real v2 library bugs (histidine, graph, charges)
- The target values don't match any decomposition from the JSON -- the provenance of the old training data is unclear
- This is why a clean-room rewrite from first principles is necessary

## Gemini's ring mutant suggestions
- Pulling oxygens off rings (TYR hydroxyl removal → effectively TYR→PHE)
- Removing second paired ring from TRP (indole → single ring)
- Worth considering alongside the three categories above
- These test ring-type-specific intensity parameters directly
