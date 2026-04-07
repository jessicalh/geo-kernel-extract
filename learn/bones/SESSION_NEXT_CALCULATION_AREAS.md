# Session prompt: CalculationArea — first three concrete types

## What this project is

NMR shielding tensor prediction for proteins. 10 classical calculators
compute geometric kernels (ring currents, bond anisotropy, Coulomb EFG,
H-bonds, dispersion, etc.) at every atom in a protein conformation.
Each calculator makes geometric decisions — "is this source close enough?",
"is this atom inside the source distribution?", "is this atom bonded to
the ring?" — using hardwired constants. There are 41 of these decision
points. We are reifying them into CalculationArea objects so we can:
document them, calibrate them against mutation-delta R², ask if MOPAC
electronic structure can inform them, dump what's included/excluded,
and draw them as isosurfaces in the UI.

## What you must read before writing any code

Read these files IN FULL. Do not skim. Do not summarise. You need the
physics associations between the geometry and the tensors, not a checklist.

1. `nmr-shielding/GEOMETRIC_KERNEL_CATALOGUE.md` — the mathematical
   foundation. Every calculator's physics, tensor structure, and why
   the geometric parameters matter. This is the document that explains
   WHY each CalculationArea exists.

2. `nmr-shielding/OBJECT_MODEL.md` (first 200 lines minimum) — the
   ConformationResult pattern, KernelEvaluationFilter ABC, FieldValue,
   SphericalTensor. This is how results attach to conformations.

3. `nmr-shielding/src/KernelEvaluationFilter.h` — the existing filter
   framework. CalculationAreas will eventually replace or wrap this.
   Understand what it does before designing the replacement.

4. `nmr-shielding/src/PhysicalConstants.h` — every named constant.

5. `nmr-shielding/src/Ring.h` — ring type hierarchy with Intensity()
   and JBLobeOffset() per type. These are CalculationAreas too.

6. `nmr-shielding/src/ConformationAtom.h` — THE ANCHOR POINT.
   Every CalculationArea evaluation stands on a ConformationAtom.
   It has position, topology, enrichment, partial_charge, mopac_charge,
   s_pop, p_pop, valency, Wiberg bond orders, ring_neighbours,
   bond_neighbours. The area sees all of this, not just xyz.

7. `nmr-shielding/src/BiotSavartResult.cpp` — a complete calculator
   that uses RingHorizon + MultipoleInnerBoundary + RingBondedExclusion.
   Read the Compute() method to see how the three areas interact.

8. `nmr-shielding/spec/CALCULATION_AREA_CATALOGUE.csv` — the
   complete catalogue of all 41 CalculationAreas with names, values,
   code locations, abc_group assignments, and physics descriptions.

## What was agreed in the previous session

- ABC is minimal: Name(), Description(), DomainKind(). Nothing else.
- Subtypes own their query interfaces. No lowest-common-denominator API.
- Evaluation is ConformationAtom-anchored.
- Results attach to conformations (same pattern as ConformationResult).
- UI queries subtype-specific geometry (spheres, atom sets, fields).
- Templates where the answer type varies, but we follow the complexity
  minimum — don't force a template if a simple virtual works.
- The code is CamelCase throughout.
- MOPAC electronic properties on ConformationAtom mean areas can be
  charge/orbital-aware, not purely geometric.

## What to do this session

1. Come back and TALK about the physics and the object model. Tell me
   what you understand about how the calculators use geometry to gate
   kernel evaluation, why the tensor structure (T0/T1/T2) matters for
   calibrating these areas, and how ConformationAtom's electronic
   properties could inform future area boundaries. I will correct you
   where needed. Do not start writing code until we have talked.

2. Write the ABC (`CalculationArea`) and the first three concrete
   subtypes as header files with full physics documentation in comments.
   NO implementation bodies — just class definitions, constructors,
   method signatures, and English descriptions of what each method does
   and why. The three types are:
   
   - `RadialThreshold` — spatial sphere inclusion/exclusion (covers
     RingHorizon, BondAnisotropyHorizon, SingularityGuard, shells, etc.)
   - `SourceRelativeExclusion` — inner boundary scaled by source extent
     (covers MultipoleInnerBoundary)
   - `TopologicalExclusion` — bond-graph-based exclusion (covers
     RingBondedExclusion, SelfSourceExclusion, DispersionBondedExclusion)

3. After writing, do NOT implement. Instead, write a critique: what
   doesn't fit? Which of the other 38 areas would be awkward under
   this hierarchy? Where does the "here is your ConformationAtom"
   pattern break down?

4. We will then plan a binary comparison run: extract features with
   the current code, refactor to use CalculationArea objects, extract
   again, diff. Zero change in output = correct refactor.

## What NOT to do

- Do not write implementation .cpp files
- Do not modify any existing calculator code
- Do not create "utility" or "helper" abstractions
- Do not add features beyond what is described above
- Do not write a summary document of the physics — the kernel catalogue
  IS the physics document, point to it
