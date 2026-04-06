# GeometryChoice: Library-to-UI Handoff Notes

**Date:** 2026-04-07
**From:** Library/calculation agent sessions
**To:** UI agent sessions

The library now records every geometric decision made by every calculator
as a `GeometryChoice` object on the `ProteinConformation`. This document
describes what the UI receives and how to use it.

---

## What lives on the conformation

```cpp
// ProteinConformation.h
std::vector<GeometryChoice> geometry_choices;  // flat, append-only
```

After `OperationRunner::Run()` completes, this vector contains every
geometric decision from all 10 calculators. On 1ubq (1231 atoms, 4 rings,
1237 bonds, 161 H-bonds) this is on the order of a few thousand entries.

## GeometryChoice structure (sealed class)

Read-only accessors:
- `Label()` -> `const string&` -- 1-3 words: "ring current", "filter exclusion", "E-field clamp", etc.
- `Calculator()` -> `CalculatorId` -- which calculator created this
- `GroupKey()` -> `size_t` -- groups choices by source (ring index, bond index, atom index)
- `Entities()` -> `const vector<GeometryEntity>&` -- the bag of model objects
- `Numbers()` -> `const vector<NamedNumber>&` -- named numeric values
- `HasSampler()` -> `bool` -- can this choice evaluate its field at any 3D point?
- `SampleAt(Vec3)` -> `SphericalTensor` -- live field evaluation (only if HasSampler)

## GeometryEntity -- one entry in the bag

```cpp
struct GeometryEntity {
    const ConformationAtom* atom;     // non-null if this is an atom ref
    const Ring*             ring;     // non-null if this is a ring ref
    const Bond*             bond;     // non-null if this is a bond ref
    size_t                  atom_index;  // into Protein's atom list (SIZE_MAX if not atom)

    EntityRole    role;     // Source, Target, Context
    EntityOutcome outcome;  // Included, Excluded, Triggered, NotTriggered

    std::string   filter_name;  // which filter rejected (empty if not rejected)
};
```

Exactly one of {atom, ring, bond} is populated per entry. The pointers
are into the live model -- follow them directly for positions, ring
geometry, bond endpoints, etc.

## NamedNumber -- a value with semantics

```cpp
struct NamedNumber {
    std::string name;   // "horizon", "intensity", "distance", "bond_order"
    double      value;
    std::string unit;   // "A", "nA", "V/A", "count", ""
};
```

## Samplers

BiotSavart "ring current" choices carry a sampler lambda that captures
the ring geometry and evaluates the Johnson-Bovey B-field + shielding
kernel at any 3D point. This is the same physics the calculator runs.

```cpp
for (const auto& gc : conf.geometry_choices) {
    if (gc.HasSampler()) {
        SphericalTensor st = gc.SampleAt(cursor_world_pos);
        // st.T0 is the isotropic shielding at that point
        // Use for isosurfaces, field probes, tooltips
    }
}
```

HaighMallion samplers are TODO (need AccumulateTensor refactor to be
callable from a lambda).

## What each calculator records

### BiotSavartResult (CalculatorId::BiotSavart)

| Label | When | Entities | Numbers | Notes |
|-------|------|----------|---------|-------|
| "ring current" | once per ring | Ring* (Source) | intensity (nA), lobe_offset (A) | Has sampler. Deduplicated -- one per ring, not per atom. |
| "filter exclusion" | atom-ring pair rejected | Ring* (Source), Atom (Target/Excluded) | distance (A), source_extent (A) | filter_name tells you MinDistance vs DipolarNearField vs RingBonded |
| "ring shells" | per atom, after ring loop | Atom (Target) | n_within_3A, n_within_5A, n_within_8A, n_within_12A (count) | Proximity histogram |

### HaighMallionResult (CalculatorId::HaighMallion)

| Label | When | Entities | Numbers | Notes |
|-------|------|----------|---------|-------|
| "surface integral" | once per ring | Ring* (Source) | -- | Sampler TODO |
| "filter exclusion" | atom-ring pair rejected | Ring* (Source), Atom (Target/Excluded) | distance (A), source_extent (A) | Same filter_name pattern as BS |
| "adaptive refinement" | atom within 2.0 A of ring | Ring* (Source), Atom (Target/Triggered) | distance (A), L1_threshold (A), L2_threshold (A) | Subdivision was needed |

### McConnellResult (CalculatorId::McConnell)

| Label | When | Entities | Numbers | Notes |
|-------|------|----------|---------|-------|
| "filter exclusion" | atom-bond pair rejected | Bond* (Source), Atom (Target/Excluded) | distance (A), source_extent (A) | filter_name: MinDistance, SelfSource, or DipolarNearField |
| "bond anisotropy" | per atom, after bond loop | Atom (Target) | nearest_CO_dist (A), nearest_CN_dist (A) | 99.0 = sentinel (no bond found) |

### CoulombResult (CalculatorId::Coulomb)

| Label | When | Entities | Numbers | Notes |
|-------|------|----------|---------|-------|
| "E-field clamp" | E_mag > 100 V/A | Atom (Target/Triggered) | actual_E_magnitude (V/A), scale_factor | Fires rarely -- clash geometries |

### HBondResult (CalculatorId::HBond)

| Label | When | Entities | Numbers | Notes |
|-------|------|----------|---------|-------|
| "hbond resolution" | H-bond rejected at resolution | -- | donor_residue, acceptor_residue (index), seq_sep or distance, rejection reason | Not an atom-level decision |
| "filter exclusion" | atom-hbond pair rejected | donor Atom (Context), acceptor Atom (Context), field Atom (Target/Excluded) | distance (A), source_extent (A) | Three atoms in the bag |
| "hbond neighbourhood" | per atom | Atom (Target) | hbond_count_within_3_5A (count) | |

### PiQuadrupoleResult (CalculatorId::PiQuadrupole)

| Label | When | Entities | Numbers | Notes |
|-------|------|----------|---------|-------|
| "filter exclusion" | atom-ring pair rejected | Ring* (Source), Atom (Target/Excluded) | distance (A), source_extent (A) | Same pattern as BS |

### RingSusceptibilityResult (CalculatorId::RingSusceptibility)

| Label | When | Entities | Numbers | Notes |
|-------|------|----------|---------|-------|
| "filter exclusion" | atom-ring pair rejected | Ring* (Source), Atom (Target/Excluded) | distance (A), source_extent (A) | Same pattern as BS |

### DispersionResult (CalculatorId::Dispersion)

| Label | When | Entities | Numbers | Notes |
|-------|------|----------|---------|-------|
| "near-field exclusion" | ring-level DipolarNearField | Ring* (Source), Atom (Target/Excluded) | distance (A), source_extent (A) | |
| "through-bond exclusion" | atom bonded to ring vertex | Ring* (Source), Atom (Target/Excluded) | distance (A) | Dispersion's own bonded check, not KernelFilterSet |
| "dispersion taper" | once per ring | Ring* (Source) | switch_onset (A), cutoff (A) | CHARMM switching function parameters |
| "switching noise floor" | S < 1e-15 in taper zone | Ring* (Source), Atom (Target/Excluded) | vertex_distance (A) | Vertex contribution dropped |

### MopacCoulombResult (CalculatorId::MopacCoulomb)

| Label | When | Entities | Numbers | Notes |
|-------|------|----------|---------|-------|
| "mopac E-field clamp" | E_mag > 100 V/A | Atom (Target/Triggered) | actual_E_magnitude (V/A), scale_factor | Same as Coulomb |
| "mopac charge floor" | any zero-charge sources skipped | Atom (Target) | zero_charge_skipped (count) | Per-atom count, not per-pair |

### MopacMcConnellResult (CalculatorId::MopacMcConnell)

| Label | When | Entities | Numbers | Notes |
|-------|------|----------|---------|-------|
| "filter exclusion" | atom-bond pair rejected | Bond* (Source), Atom (Target/Excluded) | distance (A), source_extent (A) | |
| "mopac bond anisotropy" | per atom | Atom (Target) | nearest_CO_dist (A), nearest_CN_dist (A), zero_bo_skipped (count) | Bond order floor count included |

## UI visualization ideas

### Atom inspector panel
When user picks an atom, walk `geometry_choices` and collect every entry
where that atom appears as Target. Group by Calculator. Show:
- Which rings/bonds affected this atom
- Which were excluded and why (filter_name)
- Shell counts, nearest distances
- Whether any clamps or guards fired

### Ring/bond inspector
When user picks a ring or bond, walk `geometry_choices` and collect every
entry where that ring/bond appears as Source. Show:
- How many atoms it reached
- How many were excluded by each filter
- Its parameters (intensity, lobe offset, taper range)

### Field sampler overlay
For choices with `HasSampler()`, evaluate on a grid around the source
and render as:
- Isosurface at a chosen T0 level (shielded/deshielded regions)
- Color-mapped slice plane through the ring
- Field lines from the B-field

### Exclusion zone overlay
For DipolarNearFieldFilter exclusions, draw a sphere at the ring/bond
center with radius = source_extent/2 (the exclusion boundary). Color
excluded atoms differently from included ones.

### Shell overlay
For "ring shells" choices, draw concentric transparent shells at 3/5/8/12 A
from each ring center. Color-code by how many rings fall in each shell.

## Filter set composition per calculator

For reference, the filters in each calculator's KernelFilterSet:

| Calculator | Filters (in order) |
|---|---|
| BiotSavart | MinDistance, DipolarNearField, RingBonded |
| HaighMallion | MinDistance, DipolarNearField, RingBonded |
| McConnell | MinDistance, SelfSource, DipolarNearField |
| Coulomb | MinDistance, SelfSource |
| HBond | MinDistance, SelfSource, DipolarNearField |
| PiQuadrupole | MinDistance, DipolarNearField, RingBonded |
| RingSusceptibility | MinDistance, DipolarNearField, RingBonded |
| Dispersion | DipolarNearField only (ring-level); vertex-level has own checks |
| MopacCoulomb | MinDistance, SelfSource |
| MopacMcConnell | MinDistance, SelfSource, DipolarNearField |

## Key files

- `src/GeometryChoice.h` -- the sealed class, builder, free functions
- `src/GeometryChoice.cpp` -- builder implementation
- `src/KernelEvaluationFilter.h` -- MinDistanceFilter, LastRejectorName()
- `src/ProteinConformation.h` -- the `geometry_choices` vector
