# Assessment: Layer 0 Foundation (2026-04-01)

Honest assessment of what was built, what works, what's suspect.
Written at the end of the session that built it.

## What's solid

- **SphericalTensor**: Mathematically verified, 10 tests, isometric
  normalization. Adapted from hard-won old code. Correct.
- **Ring type hierarchy**: 8 classes, virtual properties, factory.
  Calculator code is ring-type-agnostic. Clean.
- **ConformationResult framework**: Template mechanism, dependency
  checking, singleton guarantee, diagnostic messages. The architectural
  spine holds.
- **EnrichmentResult**: Zero string comparisons. Every role from typed
  properties. The PATTERNS.md vision realised in code.
- **MolecularGraphResult**: Clean BFS on typed bond graph. Fixes old
  Bug 2 (different bond thresholds).
- **Bond classification**: Fixed after the string machine debacle.
  Uses backbone index cache. Correct pattern.

## What works but needs validation

- **ApbsFieldResult**: Real PB solve, real grid interpolation, traceless
  EFG. Differs from vacuum Coulomb (good). But E-field magnitudes not
  cross-checked against published values. Being different is not the
  same as being right.
- **XtbChargeResult**: Calls real binary, parses JSON. But only runs
  on atom subsets. Whole-protein fragment decomposition not implemented.
  HOMO-LUMO gap plausible but not validated.
- **SpatialIndexResult**: nanoflann wired, 15A neighbours. Mean 292
  neighbours for 1UBQ seems high but geometrically plausible. Spot-check
  distances.
- **DsspResult**: Correct SS for 1UBQ. Clunky temp-file approach but
  works. Helix/sheet assignments match known structure.
- **PDB reader**: Works for crystal structures. Thin — doesn't fully
  use NamingRegistry. Needs attention for non-crystal structures.

## What's broken

- **ChargeAssignmentResult -75.46 for 1UBQ**: RESOLVED. Not a lookup
  bug — 1UBQ crystal structure has no hydrogens. ff14SB charges sum to
  zero per residue only when all atoms (including H) are present. The
  fix is upstream: protonate first, then charge. The protonation
  pipeline (ProtonationState, PropkaProtonator, KamlProtonator, typed
  ChargeSource) was built to address this. For protonated structures
  (ORCA runs from prmtop+XYZ, GROMACS fleet from .tpr), total charge
  is now exact integer. ARG was also corrected: was `is_titratable=false,
  charged_formal_charge=0`, now `true, +1` with ARN deprotonated variant.

## What the agents taught us

- Agents default to strings unless given strong reasons not to. The
  first agent built a string machine for bond classification. The
  correction agent fixed it.
- Agents bail on OS-level issues (library linking, tool paths). We did
  the plumbing ourselves (OpenBabel, APBS bridge, ToolPaths).
- Agents produce structurally competent code. They do NOT validate
  physics. Every number needs checking.
- Feeding PATTERNS.md inline in the prompt works better than "go read
  this file."
- The correction pass pattern works: implement, review, fix. Each pass
  improves the code and updates PATTERNS.md.
- Protein-by-value, gratuitous move constructors, and exception
  hierarchies are agent defaults that we had to correct.

## What's next

1. ~~Fix the charge bug~~ DONE — protonation pipeline built, charges
   correct on protonated structures (ORCA: -4.0 exact, CHARMM fleet: -2.0 exact)
2. tleap topology builder (open PDB path, prmtop regeneration for ~95 early runs)
3. MutantProteinConformationComparison (two complete Proteins, atom matching)
4. MD trajectory loading (GROMACS walker XTC files)
5. 4-way maths critique using MATHS_GOALS.md
6. Classical calculators (Layer 1): BiotSavart first
7. ParameterCorrectionResult training

## Numbers to verify before trusting

- ff14SB charges: must sum to integer net charge per protein
  **VERIFIED**: PrmtopChargeSource on ORCA prmtop → -4.0000 exact
  **VERIFIED**: GmxTprChargeSource on CHARMM .tpr → -2.0000 exact
  **UNDERSTOOD**: ff14SB on heavy-atom-only crystal PDB → -75.46 is
  CORRECT for that input (missing H charges). Not a bug — incomplete input.
- APBS E-field: cross-check against published values for standard proteins
- xTB charges: must sum to net charge, compare to ff14SB
- ORCA shielding: total = dia + para verified at every matrix element
  **VERIFIED**: atom 0 (N) iso = 227.975 ppm matches raw ORCA output
- Biot-Savart (when built): proton 3A above PHE → sigma > 0, known value
- McConnell (when built): known geometry → known tensor value
- All T2 components: non-zero near relevant sources
