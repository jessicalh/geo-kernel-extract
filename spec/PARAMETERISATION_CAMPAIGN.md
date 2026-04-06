# Parameterisation Campaign: CalculationArea + Ring Refactor

Session date: 2026-04-06 (afternoon session, follows MOPAC calculator and
Round 2 extraction sessions earlier today).

This document is the complete record of the CalculationArea design session.
It preserves the reasoning, decisions, and actual conversation text so that
the next session can pick up without re-deriving anything.

---

## 1. What this project is and why we are here

NMR shielding tensor prediction for proteins. 10 classical calculators
compute geometric kernels at every atom. Each calculator makes geometric
decisions using hardcoded constants scattered across PhysicalConstants.h,
Ring.h, DispersionResult.cpp, and individual calculator implementations.
There are 39 of these decision points (catalogued in
spec/CALCULATION_AREA_CATALOGUE.csv, previously thought to be 41 — some
were consolidated).

We are reifying them into CalculationArea objects so we can: document
them, give them provenance, calibrate them against mutation-delta T2 R^2,
ask if MOPAC electronic structure can inform them, dump what's included/
excluded, and draw them as isosurfaces in the UI.

The current system achieves R^2 > 0.60 on T2 deltas with these
undocumented cutoffs. The hypothesis: if sloppy constants get 0.60,
properly documented and eventually calibrated constants might do better.
But calibration is a LATER goal. The NOW goal is:

> "We go from some stuff in the code with no provenance and an AI'd
> number run to a specific object of a specific type which fully
> documents all choices and on its first pass does EXACTLY what the
> current code does."

---

## 2. Session briefing

The session was driven by `spec/SESSION_NEXT_CALCULATION_AREAS.md`.
Required reading (all read in full):

1. GEOMETRIC_KERNEL_CATALOGUE.md — mathematical foundation
2. OBJECT_MODEL.md (first 250 lines) — ConformationResult pattern
3. KernelEvaluationFilter.h — existing filter framework
4. PhysicalConstants.h — every named constant
5. Ring.h — ring type hierarchy
6. ConformationAtom.h — the anchor point
7. BiotSavartResult.cpp — a calculator using all three patterns
8. CALCULATION_AREA_CATALOGUE.csv — all 39 named areas

---

## 3. Design constraints agreed (from previous session + this session)

1. ABC is minimal: Name(), Description(), DomainKind(). Nothing else.
2. Subtypes own their query interfaces — no lowest-common-denominator API.
3. At MOST two abstract layers (ABC + subtype). May be pushing it.
4. No templates to start.
5. All areas in one CalculationArea.h and CalculationArea.cpp (CamelCase).
6. This is NOT an evaluation interface. CalculationAreas do not Accept()
   or Reject(). KernelEvaluationFilter does that. CalculationAreas are the
   CONSTANTS that filters and calculators consume.
7. The CalculationArea type is the PATTERN, not the calculator.
8. Ring parameters (currents, lobe heights) MUST be in the hierarchy — they
   are choice points too. "Rings need to count."
9. The Ring.h refactor (class hierarchy -> constexpr table) must produce
   binary identical output. It is "a class that takes a conformation
   pointer in ctor, gets back the stuff through wrapper calls, and if we
   want to fiddle the magic numbers, we can."

---

## 4. Physics discussion: how calculators use geometry

### How calculators gate kernel evaluation

Every calculator evaluates a kernel at atom-source pairs. The 39
CalculationAreas are the decision boundaries:

- **Outer radial horizons** (RingHorizon 15A, BondAnisotropyHorizon 10A):
  where kernel contribution becomes negligible. Physics: 1/r^3 at 15A is
  0.03% of the 3A value.

- **Source-relative inner boundary** (MultipoleInnerBoundary, 0.5 x extent):
  multipole expansion diverges inside the source distribution. Without it,
  HBond T2 max goes from 0.78 to 1908 A^-3.

- **Topological exclusions** (RingBondedExclusion, SelfSource): bond graph,
  not distances. Through-bond regime. The TRP additivity proof (T0(TRP5) +
  T0(TRP6) / T0(TRP9) = 1.000) only holds AFTER RingBondedExclusion removes
  the shared CD2-CE2 edge atoms.

- **Singularity guard** (0.1A): numerical, not physical.

### Why T0/T1/T2 matters for future calibration

- T0 constrains magnitudes (I, Delta_chi, Buckingham coefficients).
- T2 constrains geometry (lobe offset, midpoint shift, boundary placement).
- T1 is the magnetic dipole signature (exists in McConnell, ring
  susceptibility, H-bond, BS, HM; absent in Coulomb EFG).

For calibration against mutation-delta R^2: T2 of the delta isolates the
angular structure of the ring's contribution. If an area boundary is
misplaced, the T2 residual at boundary atoms will be systematically wrong.

### How MOPAC electronic properties could inform future area boundaries

ConformationAtom carries mopac_charge, s_pop, p_pop, valency, Wiberg bond
orders. Three future directions:

1. Soft topological exclusion: Wiberg bond order is continuous. Hard
   exclusion may throw away legitimate signal at weakly-bonded neighbours.
2. Charge-aware Buckingham coupling: orbital populations affect
   polarizability.
3. Valency-dependent horizon: saturated atoms may be less responsive to
   distant through-space effects.

All speculative. Geometric framework first, electronic enrichment later.

---

## 5. Full hierarchy design process

### First proposal (full tree for all 39 areas)

One abstract layer, 12 concrete types, 39 named instances:

```
CalculationArea (ABC: Name, Description, DomainKind)

RadialThreshold          Radius(), Sense()     14 instances
SourceRelativeExclusion  Factor()               1 instance
TopologicalExclusion     (mechanism varies)      3 instances
SequenceGate             MinSeparation()         1 instance
SwitchingFunction        Onset(), Cutoff()       1 instance
DecayFunction            DecayLength()           2 instances
RingParameter            Value(), RingType()    10 instances
NumericalAccuracy        Threshold()             2 instances
ValueClamp               Limit()                 1 instance
ValueGate                Floor()                 2 instances
Sentinel                 Value()                 1 instance
WholeDomain              (no cutoff)             1 instance
```

### Critique agent findings (7 points)

1. **RingParameter breaks ConformationAtom-anchored evaluation.** Ring
   currents and lobe heights are properties of Ring types, not atom-source
   decisions. "No ConformationAtom evaluation can answer 'given this atom,
   what is the PHE ring current?'"

2. **Dispersion switching pair is one decision.** DispersionVertexHorizon
   (5.0A hard cutoff) + DispersionTaperOnset (4.3A smooth onset) should be
   one SwitchingFunction, not two separate areas in different types.

3. **TopologicalExclusion needs per-instance classes.** RingBonded needs
   Protein at construction (graph walk), SelfSource needs nothing (index
   check), DispersionBonded rebuilds the same sets. "These share nothing
   besides 'topological check'."

4. **Missing: NEAR_ZERO_NORM (1e-10)** — normalisation guard in 7+ files.

5. **Ring shells are bin edges, not cutoffs.** Sense enum doesn't capture
   "cumulative counter boundary."

6. **DomainKind redundancy.** SpatialSphere/SpatialExclusion/SpatialShell
   too fine; RadialThreshold instances all share the same interface.

7. **DispersionVertexHorizon listed as RadialThreshold but catalogue says
   SwitchingFunction.** Calling a switching-function endpoint a standalone
   threshold is misleading.

### Revised hierarchy after critique

Ring parameters stay IN (user decision: "rings need to count"). Dispersion
pair merged. Topological types split into separate classes.
DispersionVertexHorizon moved out of RadialThreshold.

```
CalculationArea (ABC: Name, Description, DomainKind)

RadialThreshold             9 instances (outer horizons + guards)
ShellBoundary               4 instances (ring count bin edges)
SourceRelativeExclusion     1 instance (MultipoleInnerBoundary)
RingBondedExclusion         own class, needs Protein
SelfSourceExclusion         own class, pure identity
DispersionBondedExclusion   own class, needs Protein
SequenceGate                1 instance (HBondSequenceGate)
SwitchingFunction           1 instance (DispersionTaper 4.3/5.0)
DecayFunction               2 instances
RingCurrent                 6 instances (per ring type)
LobeOffset                  4 instances (per ring size group)
NumericalAccuracy           2 instances
ValueClamp                  1 instance
ValueGate                   3 instances (+ NEAR_ZERO_NORM)
Sentinel                    1 instance
WholeDomain                 1 instance
Infrastructure              1 instance (SpatialIndexHorizon)
```

### The ring orthogonality tension

Key observation from the session:

> "My sense is we are getting more orthogonal with the calculation area
> and less orthogonal with rings."

The CalculationArea hierarchy slices by PATTERN — how a decision works.
Rings slice by SOURCE — what a ring IS. The more orthogonal the pattern
axis gets, the more a ring's identity scatters:

- RingHorizon -> RadialThreshold
- MultipoleInnerBoundary -> SourceRelativeExclusion
- RingBondedExclusion -> its own topological class
- PheRingCurrent -> RingCurrent
- SixMemberedLobeHeight -> LobeOffset
- RingShellStrong -> ShellBoundary

Six different types for one ring doing its job in BiotSavart. This
fragmentation IS orthogonal — it just takes us further into Ring.h.
The ring-centric view is a query over the index, not a separate hierarchy.

Ring.h stays the home for "what IS a ring." CalculationAreas index the
CHOICES, organized by choice pattern.

---

## 6. Ring object model proposal

An agent analysed all 12 ring-related source files and produced
`spec/RING_OBJECT_MODEL_PROPOSAL.md`. Key findings:

### Current Ring.h is a lookup table disguised as a class hierarchy

8 leaf classes, 3 intermediate classes, 1 base class, 1 factory function.
Every leaf class returns 5-6 hardcoded constants via virtual methods. No
calculator ever downcasts. The virtual dispatch adds nothing.

### Proposed: constexpr table + single concrete Ring class

```cpp
struct RingTypeDescriptor {
    RingTypeIndex   index;
    const char*     name;
    double          intensity;       // nA
    double          jb_lobe_offset;  // Angstroms
    int             ring_size;
    int             nitrogen_count;
    RingAromaticity aromaticity;
};

constexpr std::array<RingTypeDescriptor, 8> RING_TYPE_TABLE = {{
    { PheBenzene,   "PHE",  -12.00, 0.64, 6, 0, Full    },
    { TyrPhenol,    "TYR",  -11.28, 0.64, 6, 0, Full    },
    { TrpBenzene,   "TRP6", -12.48, 0.64, 6, 0, Full    },
    { TrpPyrrole,   "TRP5",  -6.72, 0.52, 5, 1, Reduced },
    { TrpPerimeter, "TRP9", -19.20, 0.60, 9, 1, Full    },
    { HisImidazole, "HIS",   -5.16, 0.50, 5, 2, Weak    },
    { HidImidazole, "HID",   -5.16, 0.50, 5, 2, Weak    },
    { HieImidazole, "HIE",   -5.16, 0.50, 5, 2, Weak    },
}};
```

Ring becomes concrete (no ABC, no virtual methods). Same public API:
Intensity(), JBLobeOffset(), TypeName() — inline delegation to table.

### What this gives CalculationAreas

Concrete addresses. "PheRingCurrent" is RING_TYPE_TABLE[0].intensity.
All 8 values for any property visible in one place. No indirection
through virtual methods on subclasses.

### Additional improvements in the proposal

- `fused_perimeter_index` on Ring: TRP sub-rings navigate to TRP9
- Bonded exclusion sets computed once on Protein, not rebuilt per calculator
- Intensity() noted as not consumed at kernel evaluation time — the
  geometric kernel G is intensity-independent (evaluated at unit current).
  Whether I stays as a literature constant, becomes a learnable weight, or
  gets replaced by a calculated estimate from surrounding CalculationArea
  values is an open question for the calibration campaign. Finding that
  PHE intensity in a specific environment is better predicted by a function
  of 6 nearby area values than by Case (1995) would be a genuine finding.

### Constraint from the user

The Ring.h refactor must produce binary identical output. The ring class
takes a conformation pointer in the constructor, reaches back to protein
definitions through wrapper calls. If we want to fiddle the magic numbers,
we can — that is the whole point.

The coupling between intensity and lobe offset (change d and you need
different I to match experiment) is "a ball of worms but one we need to
open just not today."

---

## 7. What was written this session

### CalculationArea.h (src/CalculationArea.h)

ABC + 4 concrete types. Header only, no implementation bodies.

**CalculationArea (ABC)**
- Name() -> const char*
- Description() -> const char*
- Domain() -> DomainKind

**DomainKind enum**: 13 values (Spatial, Shell, SourceRelative, Topological,
Sequence, Switching, Decay, RingMagnitude, RingGeometry, Numerical,
ValueThreshold, Sentinel, WholeDomain).

**RadialThreshold**: Radius(), ThresholdSense(). Parameterised by name,
description, domain, radius, sense. 14 instances documented in comments
(9 spatial + 4 shell + 1 infrastructure). DispersionVertexHorizon
deliberately excluded (it's a switching function endpoint, not a standalone
threshold).

**SourceRelativeExclusion**: Factor(). One instance: MultipoleInnerBoundary
(factor=0.5). Full physics documentation on multipole convergence. Cites
demonstrated effect numbers (HBond T2 max with/without).

**RingBondedExclusion**: Takes Protein& in constructor. Builds per-ring
exclusion sets. ExcludedAtoms(ring_index) returns the pre-computed set.
Documents TRP additivity proof. Notes DispersionResult.cpp reimplements
the same walk independently.

**SelfSourceExclusion**: Pure identity — Name, Description, Domain only.
No query method. The comparison logic lives on SelfSourceFilter. This was
corrected during critique: the original version had IsExcluded() which
broke the "areas are constants" rule.

### Critique and fixes applied

1. Removed IsExcluded() from SelfSourceExclusion (broke "areas are
   constants" rule — the filter evaluates, the area holds identity).
2. Removed DispersionVertexCutoff from RadialThreshold instance list
   (it's a switching function endpoint, not a standalone threshold).
3. Added <climits> include.

### Ring object model proposal (spec/RING_OBJECT_MODEL_PROPOSAL.md)

Full proposal document produced by agent. Covers: what a ring IS (3
categories: type constants, topology, per-conformation geometry), current
model strengths/weaknesses, proposed RingTypeDescriptor table + concrete
Ring class, calculator needs table, complexity budget, migration path.

---

## 8. What was completed after initial writeup

All 15 concrete types written as class definitions in CalculationArea.h
(later reduced to 14 — InfrastructureArea killed, SpatialIndexHorizon
becomes a RadialThreshold instance). Full implementations in
CalculationArea.cpp for all 14 types. RingBondedExclusion topology walk
mirrors existing RingBondedExclusionFilter constructor.

## 9. What was NOT done this session
- No modifications to existing calculator code.
- No Ring.h refactor (proposal written, not implemented).
- No remaining CalculationArea types written (SequenceGate, SwitchingFunction,
  DecayFunction, RingCurrent, LobeOffset, NumericalAccuracy, ValueClamp,
  ValueGate, Sentinel, WholeDomain, ShellBoundary, Infrastructure).
- No factory function creating named instances.
- No binary comparison run (extract features, refactor, extract again, diff).
- No calibration work (T2 R^2 optimisation over area parameters).

---

## 9. Next session plan

### Phase 1: Complete the header (CalculationArea.h)

Write the remaining concrete types. Priority order:

1. **RingCurrent** and **LobeOffset** — depends on Ring.h refactor decision.
   If we do the refactor first, these reference RING_TYPE_TABLE entries.
   If not, they reference Ring virtual methods. Either way, they're
   parameterised by RingTypeIndex + value.

2. **SequenceGate** (MinSeparation=2, HBond)
3. **SwitchingFunction** (DispersionTaper onset=4.3, cutoff=5.0)
4. **DecayFunction** (RingProximityDecay 4.0A, GraphDistanceDecay 4.0 hops)
5. **ShellBoundary** (4 ring count shells)
6. **NumericalAccuracy** (HmRefineNear 2.0, HmRefineClose 1.0)
7. **ValueClamp** (EFieldSanityClamp 100.0 V/A)
8. **ValueGate** (MopacBondOrderFloor 1e-6, ChargeNoiseFloor 1e-15,
   NearZeroNorm 1e-10)
9. **Sentinel** (NoDataMarker 99.0)
10. **WholeDomain** (CoulombWholeProtein)
11. **Infrastructure** (SpatialIndexHorizon 15.0)

### Phase 2: Ring.h refactor (separate commit)

1. Add RingTypeDescriptor and RING_TYPE_TABLE to Ring.h (or new file).
2. Replace class hierarchy with concrete Ring + table lookup.
3. Add fused_perimeter_index.
4. Move bonded exclusion sets to Protein.
5. Verify binary identical extraction output (full 725-protein run, diff
   all .npy files).

### Phase 3: Integration (TBD — user decides how areas enter the pipeline)

How CalculationArea objects get instantiated and consumed by calculators
is an open design question. No factory function assumed.

### Phase 4: Calibration campaign (the reason all of this exists)

1. Systematic sweep over area parameter values.
2. Measure T2 R^2 for each parameter setting.
3. Look for the Gaussian in results — optimal values should show a peak.
4. The current "what the heck cutoffs" are R^2 > 0.60. How much is
   recoverable by getting the geometry right?

---

## 10. Files created/modified this session

| File | Status | What |
|------|--------|------|
| src/CalculationArea.h | NEW | ABC + 4 concrete types, no implementations |
| spec/RING_OBJECT_MODEL_PROPOSAL.md | NEW | Ring refactor proposal |
| spec/PARAMETERISATION_CAMPAIGN.md | NEW | This document |
| spec/SESSION_NEXT_CALCULATION_AREAS.md | EXISTS | Session briefing (not modified) |
| spec/CALCULATION_AREA_CATALOGUE.csv | EXISTS | All 39 areas (not modified) |

---

## 11. Key quotes from this session (for context recovery)

On the goal:
> "We go from some stuff in the code with no provenance and an AI'd
> number run to a specific object of a specific type which fully
> documents all choices and on its first pass does EXACTLY what the
> current code does. The only reason the learning stuff is in there is
> that we are at just over 0.6 on the last T2 runs — our simple
> calculators with what the heck cutoffs manage that. So maybe doing
> this right has some possibilities."

On ring orthogonality:
> "My sense is we are getting more orthogonal with the calculation area
> and less orthogonal with rings."
>
> "That kind of fragmentation IS orthogonal it just takes us further in
> to rings.h. Give an agent the job of an object model for rings."

On ring parameters and coupling:
> "I think that coupled is a ball of worms but one we need to open just
> not today."

On hierarchy depth:
> "At MOST two abstract layers and that may end up pushing it."

On the Ring.h refactor:
> "Provided that conversion ends up in binary identical output, sure:
> but it is a class with that in it that takes a conformation pointer
> in ctor, and gets back the stuff through wrapper calls, and if we
> want to fiddle the magic numbers, we can."

On what CalculationArea types represent:
> "Type is pattern not calculator."

On ring inclusion:
> "Rings need to count and that means this needs to include our approach
> to ring geometry or it won't do what we want."

On the concept:
> "If there are six types of ring and six things or 20 that might exist
> per ring, that's fair dinkum — we have described what we do math on
> in words and the code is 10x clearer. But only if the concept is
> sensible and maps well."
