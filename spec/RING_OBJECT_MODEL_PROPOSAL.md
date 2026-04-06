# Ring Object Model Proposal

What a ring IS, what it carries, and how to restructure Ring.h so that
type constants, topology, and per-conformation geometry each live in the
right place -- without adding abstraction that does not serve the
geometry or mathematics.

---

## 1. What a ring IS in this system

Every ring property falls into exactly one of three categories:

### Type constants (compile-time, 8 values total per property)

| Property         | PHE    | TYR    | TRP6   | TRP5   | TRP9   | HIS   | HID   | HIE   |
|------------------|--------|--------|--------|--------|--------|-------|-------|-------|
| Intensity (nA)   | -12.00 | -11.28 | -12.48 | -6.72  | -19.20 | -5.16 | -5.16 | -5.16 |
| JB lobe d (A)    | 0.64   | 0.64   | 0.64   | 0.52   | 0.60   | 0.50  | 0.50  | 0.50  |
| Ring size         | 6      | 6      | 6      | 5      | 9      | 5     | 5     | 5     |
| N count           | 0      | 0      | 0      | 1      | 1      | 2     | 2     | 2     |
| Aromaticity       | Full   | Full   | Full   | Reduced| Full   | Weak  | Weak  | Weak  |

Source: literature values (Case 1995, Johnson & Bovey 1958, Osapay & Case 1991).

### Topology (set once at protein construction, const thereafter)

- `atom_indices`: vertex atom indices in the Protein (5, 6, or 9 atoms)
- `type_index`: which of the 8 RingTypeIndex values
- `parent_residue_index`: which Residue owns this ring
- `parent_residue_number`: PDB residue number (for diagnostics)
- `fused_partner_index`: for TRP5/TRP6, the index of the other sub-ring (SIZE_MAX otherwise)

### Per-conformation geometry (recomputed per MD frame by GeometryResult)

- `center`: centroid of vertex positions
- `normal`: SVD-derived normal to best-fit plane
- `radius`: mean vertex distance from center
- `vertices`: ordered vertex positions (Vec3 array)

Stored on `ProteinConformation::ring_geometries[ri]`, NOT on Ring.

---

## 2. Current model strengths

The current Ring.h gets several things right:

1. **Calculator code is ring-type-agnostic.** Every calculator calls
   `ring.Intensity()`, `ring.JBLobeOffset()`, `ring.type_index` -- never
   downcasts or switches on type. This is the correct abstraction
   boundary.

2. **Geometry is separated from identity.** RingGeometry lives on
   ProteinConformation, not on Ring. A Ring's identity does not change
   between conformations; only vertex positions do.

3. **RingAccumulated tracks diagnostic state** without polluting the
   per-type constants or the per-conformation geometry.

4. **The factory pattern** (CreateRing from RingTypeIndex) keeps
   construction centralised.

---

## 3. Current model weaknesses

1. **8 classes returning hardcoded scalars.** PheBenzeneRing, TyrPhenolRing,
   etc. are entire classes whose only job is returning 5-6 constants. A
   class hierarchy is the wrong tool for a lookup table. The virtual
   dispatch adds nothing: no calculator ever holds a PheBenzeneRing*; they
   all hold Ring* and call the same virtual methods.

2. **SixMemberedRing/FiveMemberedRing/FusedRing intermediate classes are
   load-free.** SixMemberedRing overrides only RingSizeValue. That is a
   type constant (derivable from the ring size in the table). The
   intermediate classes add an inheritance level that carries no data and
   enables no polymorphism that the leaf classes do not already provide.

3. **Intensity() and LiteratureIntensity() are always equal.** Every
   subclass returns the same value for both. The intent was to allow a
   fitted intensity to differ from literature, but this never materialised.
   Two virtual methods returning the same number is dead design.

4. **No explicit fused-ring relationship.** TRP's three rings (TRP5 +
   TRP6 + TRP9) and the additivity property I(TRP9) = I(TRP5) + I(TRP6)
   are implicit knowledge. The model stores `fused_partner_index` but
   has no way to navigate from a sub-ring to the perimeter ring, or to
   express that the perimeter ring's atom_indices is the union of the
   sub-rings' vertex sets (minus shared edge).

5. **CalculationAreas have no clean reference point.** The CATALOGUE.csv
   lists "PheRingCurrent = -12.00 nA, defined at Ring.h:93" -- pointing
   at a virtual method on a specific subclass. There is no single place
   where all 8 type constants for a given property are visible together.

6. **Dispersion rebuilds bonded-exclusion sets** from scratch in its
   Compute() via a local BondedToVertices() function, duplicating what
   RingBondedExclusionFilter builds from the same topology. The ring
   model does not own its exclusion set.

---

## 4. Proposed object model

### 4a. RingTypeDescriptor -- the type constant table

Replace 8 leaf classes with a struct-of-constants, one instance per
RingTypeIndex, stored in a static constexpr array.

```cpp
struct RingTypeDescriptor {
    RingTypeIndex   index;
    const char*     name;            // "PHE", "TRP5", etc.
    double          intensity;       // nA (literature, calibration starting point)
    double          jb_lobe_offset;  // Angstroms
    int             ring_size;       // 5, 6, or 9
    int             nitrogen_count;
    RingAromaticity aromaticity;
};

// Compile-time table, indexed by RingTypeIndex.
constexpr std::array<RingTypeDescriptor, 8> RING_TYPE_TABLE = {{
    { RingTypeIndex::PheBenzene,   "PHE",  -12.00, 0.64, 6, 0, RingAromaticity::Full    },
    { RingTypeIndex::TyrPhenol,    "TYR",  -11.28, 0.64, 6, 0, RingAromaticity::Full    },
    { RingTypeIndex::TrpBenzene,   "TRP6", -12.48, 0.64, 6, 0, RingAromaticity::Full    },
    { RingTypeIndex::TrpPyrrole,   "TRP5",  -6.72, 0.52, 5, 1, RingAromaticity::Reduced },
    { RingTypeIndex::TrpPerimeter, "TRP9", -19.20, 0.60, 9, 1, RingAromaticity::Full    },
    { RingTypeIndex::HisImidazole, "HIS",   -5.16, 0.50, 5, 2, RingAromaticity::Weak    },
    { RingTypeIndex::HidImidazole, "HID",   -5.16, 0.50, 5, 2, RingAromaticity::Weak    },
    { RingTypeIndex::HieImidazole, "HIE",   -5.16, 0.50, 5, 2, RingAromaticity::Weak    },
}};

inline const RingTypeDescriptor& RingTypeInfo(RingTypeIndex t) {
    return RING_TYPE_TABLE[static_cast<int>(t)];
}
```

This is the single reference point for CalculationAreas. "PheRingCurrent"
is `RING_TYPE_TABLE[0].intensity`. All 8 values for JB lobe offset are
visible in one table. No virtual dispatch, no class hierarchy, no factory.

### 4b. Ring -- the topology instance

Ring becomes a concrete class (no ABC, no virtual methods). One per ring
in the protein. Owns topology, delegates type queries to the descriptor.

```cpp
class Ring {
public:
    // Topology (set at construction, const thereafter)
    std::vector<size_t> atom_indices;
    RingTypeIndex       type_index = RingTypeIndex::PheBenzene;
    size_t              parent_residue_index = 0;
    int                 parent_residue_number = 0;
    size_t              fused_partner_index = SIZE_MAX;
    size_t              fused_perimeter_index = SIZE_MAX;  // NEW: TRP5/TRP6 -> TRP9

    // Type queries (inline, no virtual dispatch)
    const RingTypeDescriptor& Type() const { return RingTypeInfo(type_index); }
    double      Intensity()      const { return Type().intensity; }
    double      JBLobeOffset()   const { return Type().jb_lobe_offset; }
    int         RingSizeValue()  const { return Type().ring_size; }
    int         NitrogenCount()  const { return Type().nitrogen_count; }
    RingAromaticity Aromaticity()const { return Type().aromaticity; }
    const char* TypeName()       const { return Type().name; }
    int         TypeIndexAsInt() const { return static_cast<int>(type_index); }

    // Topology queries
    bool IsFused()    const { return fused_partner_index != SIZE_MAX; }
    bool IsSubRing()  const { return fused_perimeter_index != SIZE_MAX; }
    bool IsPerimeter()const { return type_index == RingTypeIndex::TrpPerimeter; }

    // Geometry (delegates to conformation -- Ring does not own positions)
    RingGeometry ComputeGeometry(const std::vector<Vec3>& positions) const;

    // Diagnostic accumulation (mutable, written by post-pass)
    RingAccumulated accumulated;
};
```

### 4c. RingGeometry -- unchanged

```cpp
struct RingGeometry {
    Vec3              center = Vec3::Zero();
    Vec3              normal = Vec3::Zero();
    double            radius = 0.0;
    std::vector<Vec3> vertices;
};
```

Lives on `ProteinConformation::ring_geometries[ri]`. No change needed.

### 4d. TRP fused ring navigation

A TRP residue produces three Ring objects. The relationships:

```
TRP5 (pyrrole, 5 atoms) --fused_partner_index--> TRP6
TRP6 (benzene, 6 atoms) --fused_partner_index--> TRP5
TRP5 --fused_perimeter_index--> TRP9
TRP6 --fused_perimeter_index--> TRP9
TRP9 (perimeter, 9 atoms) --fused_partner_index--> SIZE_MAX (no partner)
TRP9 has no fused_perimeter_index (it IS the perimeter)
```

The additivity property is a statement about the type table:
`RING_TYPE_TABLE[TrpPerimeter].intensity == RING_TYPE_TABLE[TrpPyrrole].intensity + RING_TYPE_TABLE[TrpBenzene].intensity`

This can be a static_assert if the table is constexpr.

### 4e. Bonded exclusion set -- owned by Protein, not rebuilt per calculator

The per-ring exclusion set (ring vertices + bonded neighbours) is pure
topology. It should be computed once when rings are added to the Protein
and stored alongside them.

```cpp
// On Protein:
const std::set<size_t>& RingBondedExclusion(size_t ring_index) const;
```

RingBondedExclusionFilter and DispersionResult::BondedToVertices() both
query this instead of rebuilding it. The filter constructor becomes trivial.

---

## 5. What calculators need from a ring

| Property             | BS | HM | RingSusc | PiQuad | Disp |
|----------------------|----|----|----------|--------|------|
| type_index           | Y  | Y  | --       | Y      | Y    |
| Intensity()          | *  | *  | --       | --     | --   |
| JBLobeOffset()       | Y  | --  | --       | --     | --   |
| atom_indices         | -- | --  | --       | --     | Y    |
| TypeIndexAsInt()     | Y  | Y  | --       | Y      | Y    |
| fused_partner_index  | -- | --  | --       | --     | --   |
| **RingGeometry:**    |    |    |          |        |      |
| center               | Y  | Y  | Y        | Y      | Y    |
| normal               | Y  | Y  | Y        | Y      | --   |
| radius               | Y  | Y  | Y        | Y      | Y    |
| vertices             | Y  | Y  | --       | --     | Y    |
| **Filter inputs:**   |    |    |          |        |      |
| ring diameter (2*r)  | Y  | Y  | Y        | Y      | Y    |
| bonded exclusion set | Y  | Y  | Y        | Y      | Y    |

`*` Intensity is not consumed at kernel evaluation time. The geometric
kernel G is intensity-independent. Intensity is a learnable weight
applied at the feature level: sigma = I * G. Calculators evaluate G at
unit current (I=1 nA for BS) or unit susceptibility.

Key observation: no calculator calls Aromaticity(), NitrogenCount(), or
RingSizeValue() during computation. These are extracted as features for
the learning model but do not enter any geometric kernel. They belong on
the type descriptor for feature extraction, not in the kernel interface.

---

## 6. How CalculationAreas reference ring parameters

With the type table, each CalculationArea entry has a concrete address:

| CalculationArea name   | abc_group     | Reference                                     |
|------------------------|---------------|-----------------------------------------------|
| PheRingCurrent         | RingParameter | `RING_TYPE_TABLE[PheBenzene].intensity`        |
| Trp9LobeHeight         | RingParameter | `RING_TYPE_TABLE[TrpPerimeter].jb_lobe_offset` |
| RingHorizon            | RadialThreshold | `RING_CALC_CUTOFF` (PhysicalConstants.h)     |
| MultipoleInnerBoundary | SourceRelativeExclusion | `DipolarNearFieldFilter` (0.5 * extent)|
| RingBondedExclusion    | TopologicalExclusion | `Protein::RingBondedExclusion(ri)`       |

The type table IS the parameter catalogue for ring-specific constants.
No indirection through virtual methods. A CalculationArea that says
"the lobe offset for HIS" can point at a constexpr table entry.

---

## 7. Complexity budget

### Essential complexity (physics demands it, must be in the model)

- **8 distinct ring types** with different intensities and lobe offsets.
  The physics is different: PHE and TRP5 produce measurably different
  shielding patterns. Collapsing types loses information.

- **TRP fused ring triple** (TRP5 + TRP6 + TRP9). The sub-ring
  decomposition is load-bearing: the model sees both component and
  whole-system features, and the additivity property is verified on
  720 proteins. This is not redundancy; it is structure.

- **Per-conformation geometry** separate from type identity. A ring's
  type never changes; its vertex positions change every MD frame.

- **Bonded exclusion topology**. Through-bond vs through-space is a
  real physics boundary. Removing it breaks TRP additivity by 13%.

- **5 different calculators** consuming overlapping but distinct subsets
  of ring properties. The RingNeighbourhood struct on ConformationAtom
  is the right place to accumulate per-calculator results per ring.

### Accidental complexity (current implementation, not physics)

- **8 leaf classes** for 8 rows of a lookup table. Replace with
  RingTypeDescriptor table.

- **3 intermediate classes** (SixMembered, FiveMembered, Fused) that
  carry no data and enable no polymorphism. Eliminate.

- **Virtual dispatch** for accessing compile-time constants. Replace
  with inline delegation to constexpr table.

- **Intensity() and LiteratureIntensity()** returning the same value.
  Keep one. If calibration later produces fitted intensities, that is
  a per-calculator parameter (a CalculationArea), not a ring property.

- **Per-calculator rebuilding of bonded exclusion sets.** Compute once
  on Protein, query by ring index.

---

## 8. Migration path

The proposed Ring keeps the same public API (Intensity(), JBLobeOffset(),
TypeName(), etc.) with identical return values. The only change is the
mechanism: inline table lookup instead of virtual dispatch. This means:

1. No calculator code changes. Every call site compiles and returns the
   same value.
2. RingTypeDescriptor and RING_TYPE_TABLE can be added first, Ring
   refactored second.
3. The factory function CreateRing() can be replaced by direct Ring
   construction with a RingTypeIndex.
4. Tests that construct specific ring types (PheBenzeneRing, etc.) change
   to `Ring r; r.type_index = RingTypeIndex::PheBenzene;`.

The bonded exclusion set migration is independent and lower priority.

---

## 9. What this proposal does NOT do

- Does not add ring-type-specific behaviour that varies at runtime.
  If future physics requires it (e.g., different HM vs BS intensities
  per type), that belongs on the calculator's parameter struct, not on
  Ring. Ring is identity; parameters are calibration.

- Does not merge RingGeometry into Ring. Geometry is per-conformation
  and lives on ProteinConformation. This separation is correct.

- Does not add an ABC or template parameter to Ring. There is exactly
  one kind of ring in this system. The variation is in the type
  constants, which are data, not behaviour.

- Does not change RingNeighbourhood on ConformationAtom. That struct
  is the per-atom, per-ring accumulator for all 5 ring calculators.
  Its layout is driven by what calculators write, not by ring identity.
