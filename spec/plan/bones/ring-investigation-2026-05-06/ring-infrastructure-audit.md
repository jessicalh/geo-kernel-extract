# Ring infrastructure audit — 2026-05-07

Read-only investigation of the Ring class hierarchy, the
`Protein::DetectAromaticRings` migration target, the
`Protein::FinalizeConstruction` ordering, and the consumer surface
that Bundle C's `ConstructRingsFromSubstrate` must satisfy. Scope is
substrate-side only: no calculator updates are recommended, the
`if (ti < 8)` guards stay as-is per the locked Bundle C decision.

This document deepens `inventory.md` for the
`Infrastructure-Ring-Detection-and-Class-Hierarchy.md` slot of the
per-calculator dispatch list. The other six per-calculator reports
remain pending; this one focuses on the Ring object itself, the
construction site, and the contact-surface consumers impose.


## 1. Executive summary

- The Ring class hierarchy at `src/Ring.h:50-202` is a virtual
  3-tier dispatch (`Ring` → `SixMemberedRing` / `FiveMemberedRing` /
  `FusedRing` → leaf class per RingTypeIndex) with seven required
  virtuals: `Intensity()`, `LiteratureIntensity()`, `JBLobeOffset()`,
  `NitrogenCount()`, `Aromaticity()`, `RingSizeValue()`, `TypeName()`.
  `ProPyrrolidineRing` plugs into this hierarchy as a new
  `FiveMemberedRing` leaf with all seven virtuals overridden.
- `Ring::atom_indices` cyclic ordering is load-bearing: it determines
  the ring normal sign via the cross-product orientation fix at
  `src/Ring.cpp:32-37`. Bundle C must preserve current ordering for
  the 8 aromatic ring types to keep shielding-NPY bit-identity.
- `DetectAromaticRings` (`src/Protein.cpp:548-646`) reads atom-name
  strings at exactly four sites (568, 575, 604, 605); all four are
  Bundle C clearance targets. The substrate's typed
  `LegacyAmber().SemanticAt(ai).ring_position.primary.ring` already
  covers the chemistry these strings discriminate.
- `FinalizeConstruction` orders ring detection BEFORE bond
  resolution and BEFORE substrate composition. Bundle C requires
  reordering: the substrate-driven ring construction must run
  AFTER the substrate is composed.
- Consumer surface: 5 ring-current calculators (BS, HM, RingChi,
  PiQuad, Disp) read `protein.RingAt(ri)`; `Ring::JBLobeOffset()`,
  `Ring::Intensity()`, `Ring::TypeIndexAsInt()`, `Ring::type_index`,
  `Ring::atom_indices`. 2 indirect Coulomb consumers
  (`CoulombResult.cpp:77-81`, `MopacCoulombResult.cpp:73-75`) read
  `RingAt(ri).atom_indices` for membership-only `is_aromatic_atom`
  masks. 1 infrastructure consumer (`GeometryResult.cpp:17-19, 63-64`)
  builds `conf.ring_geometries[]` and `conf.rings_by_type[]` from
  `Ring::ComputeGeometry` and `Ring::type_index`. 1 filter constructor
  (`KernelEvaluationFilter.cpp:22-41`) walks `Ring::atom_indices`
  to build the per-ring exclusion set.
- TRP perimeter is currently constructed at runtime in
  `DetectAromaticRings` from a static atom-name list at
  `src/AminoAcidType.cpp:188`. The substrate encodes only
  `Indole_Trp_5` and `Indole_Trp_6`; perimeter is synthesised.
- DispersionResult locally reimplements `RingBondedExclusionFilter`'s
  topology walk at `DispersionResult.cpp:136-149,202-204`. The
  reimplementation is functionally equivalent (same protein bond
  graph walk, same per-ring set semantics); it is accidental
  duplication, not a science divergence.


## 2. Ring class hierarchy map

Source: `src/Ring.h:1-205`, `src/Ring.cpp:1-65`.

### 2.1 Polymorphism mechanism

Pure virtual functions on the abstract `Ring` base, with three
intermediate non-leaf classes (`SixMemberedRing`, `FiveMemberedRing`,
`FusedRing`) that supply only `RingSizeValue()` defaults (6, 5, and
none respectively). Leaf classes provide all seven virtuals (the
intermediates' `RingSizeValue` is re-overridden in `FusedRing`'s
sole leaf `IndolePerimeterRing` to return 9).

A non-virtual integer-tag dispatch is also present:
`Ring::TypeIndexAsInt()` (line 72) returns
`static_cast<int>(type_index)`; calculators use it to index into
per-type-array fields on `ConformationAtom`. The two dispatch paths
co-exist: virtuals for ring-physics queries
(`Intensity`/`JBLobeOffset`), integer tag for per-type accumulation.

`CreateRing(RingTypeIndex)` factory at `Ring.cpp:50-62` returns
`std::unique_ptr<Ring>`; switch on RingTypeIndex.

### 2.2 Class tree

```text
Ring (abstract base, src/Ring.h:50-79)
|
+-- SixMemberedRing                             (src/Ring.h:86-89)
|   |   RingSizeValue() = 6
|   +-- PheBenzeneRing       I=-12.0, JB=0.64A, Aromaticity::Full,  N_count=0
|   +-- TyrPhenolRing        I=-11.28, JB=0.64A, Aromaticity::Full, N_count=0
|   +-- TrpBenzeneRing       I=-12.48, JB=0.64A, Aromaticity::Full, N_count=0
|
+-- FiveMemberedRing                            (src/Ring.h:129-132)
|   |   RingSizeValue() = 5
|   +-- TrpPyrroleRing       I=-6.72, JB=0.52A, Aromaticity::Reduced, N_count=1
|   +-- HisImidazoleRing     I=-5.16, JB=0.50A, Aromaticity::Weak,    N_count=2
|   +-- HidImidazoleRing     I=-5.16, JB=0.50A, Aromaticity::Weak,    N_count=2
|   +-- HieImidazoleRing     I=-5.16, JB=0.50A, Aromaticity::Weak,    N_count=2
|   +-- [Bundle C ADD] ProPyrrolidineRing  I=0, JB=0,   Aromaticity::None, N_count=1
|
+-- FusedRing                                   (src/Ring.h:183)
    |   (no RingSizeValue default; leaf overrides)
    +-- IndolePerimeterRing  I=-19.2, JB=0.60A, Aromaticity::Full, N_count=1, RingSize=9
```

`Aromaticity::None` does NOT currently exist in `RingAromaticity`
enum (`src/Types.h:168` defines only `Full, Reduced, Weak`). Bundle
C adds `None` as a new enum value.

### 2.3 Per-class fields

All Rings inherit from base:

| Field                   | Type                  | Set by                 |
| ----------------------- | --------------------- | ---------------------- |
| `atom_indices`          | `std::vector<size_t>` | DetectAromaticRings (Bundle C: ConstructRingsFromSubstrate) |
| `type_index`            | `RingTypeIndex`       | Constructor (per leaf class) |
| `parent_residue_index`  | `size_t`              | DetectAromaticRings    |
| `parent_residue_number` | `int`                 | DetectAromaticRings    |
| `fused_partner_index`   | `size_t` (SIZE_MAX default) | DetectAromaticRings post-pass for TRP |
| `accumulated`           | `RingAccumulated`     | calculator post-passes |

### 2.4 Per-class virtual surface

Required virtuals on every leaf class:

| Virtual                       | Signature                          | Purpose                                         | Sites read                                                                 |
| ----------------------------- | ---------------------------------- | ----------------------------------------------- | -------------------------------------------------------------------------- |
| `Intensity()`                 | `double`                           | Calibration-time current (nA/T) via CalculatorConfig | BS:185 (record only); not used in kernel       |
| `LiteratureIntensity()`       | `double`                           | Literature reference value                      | (no calculator reads; reference for diagnostics) |
| `JBLobeOffset()`              | `double`                           | Johnson-Bovey lobe offset (Å) via CalculatorConfig | BS:181, 186, 223, 380, 406                  |
| `NitrogenCount()`             | `int`                              | Heteroatom count for diagnostic emission        | (no calculator reads at present)                                          |
| `Aromaticity()`               | `RingAromaticity`                  | Categorical aromaticity                         | (no calculator reads at present)                                          |
| `RingSizeValue()`             | `int`                              | 5/6/9                                           | (no calculator reads at present; available for substrate-extension calcs) |
| `TypeName()`                  | `const char*`                      | Short name for diagnostics                      | (logging only)                                                            |

Plus non-virtual:

| Method                     | Signature                          | Purpose                                                             |
| -------------------------- | ---------------------------------- | ------------------------------------------------------------------- |
| `IsFused()`                | `bool`                             | True iff `fused_partner_index != SIZE_MAX`                          |
| `TypeIndexAsInt()`         | `int`                              | `static_cast<int>(type_index)` for indexing                         |
| `ComputeGeometry(positions)` | `RingGeometry`                   | SVD normal + centroid + radius + vertices; cyclic-order sign rule   |

`ComputeGeometry` (`Ring.cpp:6-47`) reads `atom_indices` to gather
vertex positions; the orientation fix at lines 34-37 enforces
right-hand-rule on the first three vertices (`edge01 × edge02`).
This is where `atom_indices` cyclic ordering converts to ring-normal
sign, which calculators use as the multipole axis. Reversing
`atom_indices` reverses the normal.

### 2.5 ProPyrrolidineRing API requirement (Bundle C addition)

The class must derive from `FiveMemberedRing` and override every
virtual:

```cpp
class ProPyrrolidineRing : public FiveMemberedRing {
public:
    ProPyrrolidineRing() { type_index = RingTypeIndex::ProPyrrolidine; }
    double Intensity() const override { return 0.0; }                     // saturated; no ring current
    double LiteratureIntensity() const override { return 0.0; }
    double JBLobeOffset() const override { return 0.0; }
    int NitrogenCount() const override { return 1; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::None; }
    // RingSizeValue() inherited from FiveMemberedRing (= 5); locked
    const char* TypeName() const override { return "PRO"; }
};
```

Three substrate-prerequisite changes accompany the class:

1. `RingTypeIndex::ProPyrrolidine = 8`, `Count = 9`
   (`src/Types.h:184`).
2. `RingAromaticity::None` enum value (`src/Types.h:168`).
3. `CreateRing` factory case for `ProPyrrolidine`
   (`src/Ring.cpp:50-62`).

`Intensity()` returning `0.0` should be the literal value, not
`CalculatorConfig::Get(...)` — Pro Ring intensity is physically zero
by saturated heterocycle chemistry (Joule & Mills 2010 ch. 7), not
a calibration parameter. Other Rings reach `CalculatorConfig` because
their literature values are calibration starting points.


## 3. DetectAromaticRings full audit

Source: `src/Protein.cpp:548-646`. 99-line function.

### 3.1 Function structure

| Lines    | Block                                                                          |
| -------- | ------------------------------------------------------------------------------ |
| 548-549  | Function signature; `rings_.clear()`                                           |
| 551-557  | `TrpRingIndices` local struct + `trp_rings` map for fused-pair post-pass       |
| 559-636  | Per-residue loop                                                               |
| 561-563  | Read `aatype = res.AminoAcidInfo()`; skip if `!aatype.is_aromatic`             |
| 565-569  | **STRING-READ SITE 1**: build `name_to_idx[atoms_[ai]->pdb_atom_name] = ai`    |
| 571-583  | **STRING-READ SITE 2**: per `ring_def`, match `ring_def.atom_names[k]` strings |
| 585-615  | RingTypeIndex resolution (HIS tautomer logic)                                  |
| 600-614  | **STRING-READ SITES 3 & 4**: `name_to_idx.find("HD1")`, `find("HE2")` fallback |
| 617-623  | Construct Ring via `CreateRing(effective_type)`; populate identity fields     |
| 625-635  | Track TRP ring indices for fused-partner post-pass                            |
| 638-645  | Set `fused_partner_index` on the TRP 5-ring + 6-ring pair                     |

### 3.2 String-read sites — verbatim with substrate replacement

**Site 1: line 568.** Build per-residue name→index map.

```cpp
for (size_t ai : res.atom_indices) {
    name_to_idx[atoms_[ai]->pdb_atom_name] = ai;
}
```

Substrate replacement: walk `res.atom_indices`; for each `ai`, query
`LegacyAmber().SemanticAt(ai).ring_position.primary.ring`; collect
indices grouped by `(residue_index, RingSystemKind)`. The string map
is replaced by a typed-key map:

```text
std::map<RingSystemKind, std::vector<size_t>> ring_atoms;
// populated from substrate; each atom that has primary.ring != NotInRing
// joins the appropriate group; bridge atoms also join their secondary.ring
```

**Site 2: lines 571-582.** Per-ring iteration plus atom-name lookup.

```cpp
for (const auto& ring_def : aatype.rings) {
    std::vector<size_t> atom_indices;
    bool all_present = true;
    for (const char* aname : ring_def.atom_names) {
        auto it = name_to_idx.find(aname);
        if (it == name_to_idx.end()) {
            all_present = false;
            break;
        }
        atom_indices.push_back(it->second);
    }
```

Substrate replacement: drop the `aatype.rings[]` outer loop entirely.
Iterate `RingSystemKind` values in canonical order; for each that
the residue's substrate populates, walk the atoms in canonical
`RingPositionLabel` order to produce ordered `atom_indices`. For
PHE: `Ipso → Ortho1 → Meta1 → Para → Meta2 → Ortho2`. For TRP-5:
`Ipso → PyrroleBeta → Heteroatom_NH → BridgeFusion → BridgeFusion`
(two bridge atoms; substrate's secondary RingMembership disambiguates
the ordering). Per the locked decision Pro adopts
`N → Cα → Cβ → Cγ → Cδ` residue-walk order.

The `aatype.rings[]` const-char-pointer arrays at
`src/AminoAcidType.cpp` (PHE line 147; TYR line 201; HIS line 92;
TRP lines 186-188) become deletable, as does the entire
`AminoAcidRing` struct definition at `src/AminoAcidType.h:58-61`.
Verified single-site consumer: `grep -rn "rings\[\]\|aatype.rings"
src/` returns only `src/Protein.cpp:571` (per inventory §3.1).

**Site 3 & 4: lines 604-605.** HIS-tautomer hydrogen detection
fallback.

```cpp
bool has_HD1 = name_to_idx.find("HD1") != name_to_idx.end();
bool has_HE2 = name_to_idx.find("HE2") != name_to_idx.end();
```

Substrate replacement: this branch runs only when
`res.protonation_variant_index < 0` (substrate has no resolved
variant). After Bundle C lands at the post-FinalizeConstruction
ordering, substrate composition has populated the per-N
`RingPositionLabel` (`Heteroatom_NH` vs `Heteroatom_NoH`) before
`ConstructRingsFromSubstrate` runs; the variant index is already
resolved. The fallback path becomes structurally unreachable. Bundle
C deletes the fallback branch outright.

Per the science README §11 (fail-loud), if a residue arrives with
HIS chemistry but no resolved variant index by the time
`ConstructRingsFromSubstrate` runs, the function must FATAL-abort
rather than silently default to `HisImidazole` (the current line
598 default). This is a stronger invariant than today.

### 3.3 Behavioural items beyond the strings

- Line 583: `if (!all_present) continue` — silent-skip of residues
  with incomplete ring atoms (e.g. backbone-only residues, mutation
  caps). Bundle C should fail-loud here per README §11; substrate-
  driven construction sees the substrate's empty `ring_position` and
  must distinguish "not a ring residue" from "ring residue with
  missing atoms."
- Lines 638-645: post-pass to set `fused_partner_index` on TRP 5-ring
  and 6-ring pairs. Substrate-driven construction at the same
  function would reuse this post-pass identically; the post-pass
  reads only typed `Ring*` references, no strings.
- The TRP perimeter ring is materialised in the same outer loop
  (line 588: `effective_type = ring_def.type_index`) because
  `aatype.rings[]` for TRP includes a perimeter entry at
  `src/AminoAcidType.cpp:188`. After substrate-driven construction,
  the perimeter is either (a) synthesised at runtime by
  `ConstructRingsFromSubstrate` from union of `Indole_Trp_5` and
  `Indole_Trp_6` atoms, or (b) synthesised by an extended substrate
  with `RingSystemKind::Indole_Trp_9` per the README's strongly
  recommended default.


## 4. FinalizeConstruction ordering

Source: `src/Protein.cpp:215-374`.

### 4.1 Current ordering

```text
Step                                         Line     What it does
---------------------------------------------------------------------------------------
ResolveResidueTerminalStates                 224      Per-chain N-terminus / C-terminus tagging
CacheResidueBackboneIndices [STRING]         225      First pass: PDB-name match for N/CA/C/O/H/HA/CB
ResolveProtonationStates(bonds=nullptr)      226      First pass: HIS / LYS / ASP / GLU / TYR variants
                                                       from explicit-H presence (no bond graph yet)
DetectAromaticRings  ◄── current ring-build  227      Reads aatype.rings[].atom_names + name_to_idx
                                                       lookups; HIS variant from variant_index OR
                                                       fallback string lookup at lines 604-605
CovalentTopology::Resolve                    232      Bond graph (covalent radius check + ring info)
OverrideDisulfides (if invariants)           246-255  TPR authority overlay on disulfides
ResolveProtonationStates(bonds=non-null)     261      Second pass: CYS → CYX from disulfides
Copy bonds → Atom.bond_indices etc.          267-270  Per-atom connectivity surface
Pass-2 NamingApplicator                      292-356  Re-apply naming rules with resolved variants
ComposeAtomSemantic                          360-361  Build per-atom AtomSemanticTable substrate
LegacyAmberTopology construction             364-366  Move substrate + invariants into topology
CacheResidueBackboneIndices_Typed            373      Overwrite backbone cache from substrate
```

### 4.2 What the current ordering implies for ring construction

The current `DetectAromaticRings` runs at step 4 (line 227), which
is BEFORE `CovalentTopology::Resolve` (line 232) and BEFORE
`ComposeAtomSemantic` (line 360-361). It does NOT depend on the
substrate; it consumes the `aatype.rings[].atom_names` static table
and the per-atom `pdb_atom_name` strings.

`CovalentTopology::Resolve` at line 232-233 takes `rings_` as input:

```cpp
auto bonds = CovalentTopology::Resolve(atoms_, rings_, residues_,
                                       positions, bond_tolerance);
```

Bond resolution uses ring identity to refine bond categorisation
(aromatic bonds, ring-closing bonds). So rings must exist before
bonds for the current pipeline.

`ComposeAtomSemantic` at line 361 takes `(atoms_, residues_, *bonds)`
and produces the substrate. It does NOT take `rings_` — the substrate's
`ring_position` per atom is computed from the residue tables in
`src/generated/LegacyAmberSemanticTables.cpp` (the residue-level
encoding), not from runtime ring detection.

### 4.3 Reordering required for Bundle C

For `ConstructRingsFromSubstrate` to read the substrate, the
substrate must exist when the function runs. The substrate is
created at line 360-366 (after both bond resolution and second
protonation pass). So `ConstructRingsFromSubstrate` must move from
line 227 to AFTER line 366.

But step 232's `CovalentTopology::Resolve(atoms_, rings_, ...)`
takes `rings_` as input. If we move ring construction to post-line-366,
bond resolution sees an empty `rings_`. Two options:

**Option A: split bond-ring dependency.** Examine whether
`CovalentTopology::Resolve` actually needs `rings_` at construction
time, or only later for refinement. If only refinement, move the
refinement pass to post-rings.

**Option B: temporal dependency reversal.** The substrate's per-atom
`ring_position` is residue-table-encoded — it does not require
runtime bond information. The rings can be built from the substrate
without a complete bond graph; bond categorisation can then refine
post-rings. This requires a substrate-only ring construction pass
followed by bond resolution.

The shape question for Bundle C: does
`CovalentTopology::Resolve` USE `rings_` for chemistry decisions, or
only for diagnostic emission?

Verification needed (substrate-side investigation, not in this
pass): grep `CovalentTopology::Resolve` body for reads of
`rings_` / Ring fields. If reads are diagnostic-only, the reorder is
"move ring construction below substrate composition; bond
resolution unchanged." If reads are chemistry-load-bearing, the
reorder is harder and may require splitting bond resolution into
two passes.

### 4.4 Sequencing constraint summary

Per the locked Bundle C decisions:

```text
Step                                         Required ordering
---------------------------------------------------------------------------------------
ResolveResidueTerminalStates                 unchanged
CacheResidueBackboneIndices                  unchanged (string pass; overwritten later)
ResolveProtonationStates(nullptr)            unchanged
   [DetectAromaticRings deleted]
CovalentTopology::Resolve(atoms_, residues_, positions, bond_tolerance)
                                             rings argument REMOVED if bond resolution
                                             does not require rings; else split.
OverrideDisulfides                           unchanged
ResolveProtonationStates(bonds)              unchanged
Copy bonds → Atom.bond_indices               unchanged
Pass-2 NamingApplicator                      unchanged
ComposeAtomSemantic                          unchanged
LegacyAmberTopology construction             unchanged
   [ConstructRingsFromSubstrate ← NEW]       reads LegacyAmber().SemanticAt;
                                             builds rings_ + sets fused_partner_index
CacheResidueBackboneIndices_Typed            unchanged
```

The reorder is the only sequencing change Bundle C requires.


## 5. Consumer-surface table

This table tabulates what each Ring consumer reads from the Ring
class. Source coverage:

- `BiotSavartResult.cpp:124-503`
- `HaighMallionResult.cpp:190-432`
- `RingSusceptibilityResult.cpp:103-270`
- `PiQuadrupoleResult.cpp:109-292`
- `DispersionResult.cpp:163-423`
- `CoulombResult.cpp:43-405` (indirect: aromatic mask)
- `MopacCoulombResult.cpp:39-405` (indirect: aromatic mask)
- `GeometryResult.cpp:7-93` (infrastructure: ring_geometries + rings_by_type)
- `KernelEvaluationFilter.cpp:22-41` (filter: per-ring exclusion sets)

### 5.1 Ring fields consumed (read-side)

| Consumer            | atom_indices | type_index | parent_residue_* | fused_partner_index | accumulated |
| ------------------- | ------------ | ---------- | ---------------- | ------------------- | ----------- |
| BiotSavartResult    | indirect (via filter) | yes (302, 248) | no | no | no |
| HaighMallionResult  | indirect (via filter) | yes (336, 307) | no | no | no |
| RingSusceptibility  | indirect (via filter) | yes (189) | no | no | no |
| PiQuadrupoleResult  | indirect (via filter) | yes (213, 192) | no | no | no |
| DispersionResult    | **direct** (139, 275) | yes (336, 313) | no | no | no |
| CoulombResult       | **direct** (78) | no | no | no | no |
| MopacCoulombResult  | **direct** (74) | no | no | no | no |
| GeometryResult      | indirect (via ComputeGeometry) | yes (64, 86) | no | yes (86) | no |
| RingBondedExclusion | **direct** (30) | no | no | no | no |

Direct readers of `atom_indices`:

- `DispersionResult` reads `ring.atom_indices.size()` and iterates
  for vertex-by-vertex kernel summation (lines 275-298).
- `CoulombResult` and `MopacCoulombResult` walk
  `RingAt(ri).atom_indices` once at compute-prep to mark atoms as
  aromatic in a bool vector (lines 77-81 / 73-75).
- `RingBondedExclusionFilter` walks `Ring::atom_indices` plus
  per-vertex `bond_indices` to build per-ring exclusion sets.
- `Ring::ComputeGeometry` (in GeometryResult's flow at line 19) reads
  `atom_indices` to gather vertex positions for SVD.
- DispersionResult's local `BondedToVertices` (lines 136-149) is the
  same walk as `RingBondedExclusionFilter` — see §7.

Indirect readers (BS, HM, RingChi, PiQuad) access `atom_indices`
ONLY through filters or geometry; they never iterate vertices.

### 5.2 atom_indices ordering dependence

| Consumer            | Ordering matters?       | Why                                                                 |
| ------------------- | ----------------------- | ------------------------------------------------------------------- |
| BiotSavartResult    | YES (transitively)      | Reads `geom.normal` (line 222); SVD orientation set by atom_indices first-three cyclic walk |
| HaighMallionResult  | YES (transitively)      | Reads `geom.normal` (line 285); same dependency                      |
| RingSusceptibility  | YES (transitively)      | Reads `geom.normal` (line 155); same dependency                      |
| PiQuadrupoleResult  | YES (transitively)      | Reads `geom.normal` (line 156-157); same dependency                  |
| DispersionResult    | YES (transitively)      | Reads `geom.normal` (line 220-319); per-vertex iteration order does NOT matter (kernel is vertex-symmetric) |
| CoulombResult       | NO                      | Membership-only check; order of `atom_indices` walk is irrelevant   |
| MopacCoulombResult  | NO                      | Same                                                                |
| GeometryResult      | YES                     | `Ring::ComputeGeometry` orientation-fix at `Ring.cpp:34-37` reads first three atoms |
| RingBondedExclusion | NO                      | Insert-into-set is unordered                                         |

Net implication: Bundle C MUST preserve the cyclic walk direction
of `atom_indices` for the 8 aromatic ring types. Reversing
atom_indices for any aromatic ring flips the normal sign and the
shielding sign for every nearby probe atom. This is the
bit-identity gate of bless-compare on the standard fixtures.

The cyclic walk for each ring type today (read from
`src/AminoAcidType.cpp`):

| Ring type        | atom_names array (from AminoAcidType.cpp) |
| ---------------- | ------------------------------------------ |
| PheBenzene       | `{CG, CD1, CE1, CZ, CE2, CD2}`            |
| TyrPhenol        | `{CG, CD1, CE1, CZ, CE2, CD2}`            |
| HisImidazole/HID/HIE/HIP | `{CG, ND1, CE1, NE2, CD2}`        |
| TrpBenzene       | `{CD2, CE2, CZ2, CH2, CZ3, CE3}`          |
| TrpPyrrole       | `{CG, CD1, NE1, CE2, CD2}`                |
| TrpPerimeter     | `{CG, CD1, NE1, CE2, CZ2, CH2, CZ3, CE3, CD2}` |

`ConstructRingsFromSubstrate`'s typed walk must produce identical
sequences for these ring types. The substrate's
`RingPositionLabel` walk has an opportunity for handedness-mismatch:
if the substrate emits `Ortho1=CD1, Ortho2=CD2, Meta1=CE1,
Meta2=CE2, Para=CZ`, the canonical walk
`Ipso→Ortho1→Meta1→Para→Meta2→Ortho2` reproduces the existing PHE
order (`CG, CD1, CE1, CZ, CE2, CD2`). If the substrate's
labeling is opposite-handed, the order becomes `CG, CD2, CE2, CZ,
CE1, CD1` — same cyclic walk, opposite direction, opposite normal
sign. This must be verified empirically per ring type at Bundle C
draft time using fixture-level inspection.

For Pro: no current consumer reads Pro Ring atom_indices for
chemistry (Coulomb consumers exclude Pro because `is_aromatic =
false` on AminoAcidType). The chosen walk N→Cα→Cβ→Cγ→Cδ is
locally definable and need not match any current convention. The
constraint comes from FUTURE puckering descriptors needing the
dihedral N-Cα-Cβ-Cγ-Cδ.

### 5.3 Iteration-mode and aromaticity filtering

| Consumer            | Iteration mode                                                                    | Filtered by aromaticity? | Pro Ring iterated? |
| ------------------- | ---------------------------------------------------------------------------------- | ------------------------ | ------------------ |
| BiotSavartResult    | per-atom × per-nearby-ring (spatial cutoff 15Å), `RingsWithinRadius`               | NO; `if (ti < 8)` guards `per_type_*` accumulation | YES (would fall through `< 8` guard, total still includes Pro raw kernel) |
| HaighMallionResult  | per-atom × per-nearby-ring                                                          | NO; same `if (ti < 8)` pattern | YES (same as BS) |
| RingSusceptibility  | per-atom × per-nearby-ring                                                          | NO per-type emission     | YES (no per-type filter; total includes Pro) |
| PiQuadrupoleResult  | per-atom × per-nearby-ring                                                          | NO; same `if (ti < 8)`   | YES (same as BS) |
| DispersionResult    | per-atom × per-nearby-ring × **per-vertex** (`for vi in 0..ring.atom_indices.size()`) | NO; same `if (ti < 8)`   | YES (Pro contributes non-trivially per-vertex; through-bond exclusion catches Pro residue's own atoms) |
| CoulombResult       | per-atom × residue-walk × per-ring-walk (mark-aromatic only)                        | NO                       | YES (Pro atoms get marked, but Pro is NOT aromatic per AminoAcidType so... wait — Coulomb iterates ALL ring atoms; Pro Ring atoms WOULD get marked is_aromatic_atom=true) |
| MopacCoulombResult  | same as CoulombResult                                                               | NO                       | YES (same — Pro atoms would be marked as aromatic source) |
| GeometryResult      | per-ring (ring_geometries + rings_by_type)                                          | NO                       | YES — `conf.rings_by_type[ProPyrrolidine]` populated naturally |
| RingBondedExclusion | per-ring (build per-ring exclusion set)                                              | NO                       | YES |

**Critical Pro-specific finding for indirect Coulomb consumers:**

`CoulombResult.cpp:77-81` and `MopacCoulombResult.cpp:73-75` walk
EVERY Ring's atom_indices (the inventory's analysis at
inventory.md:974-984 stated "Pro has `is_aromatic = false` so the
source class is unchanged"). That analysis was incorrect on the
narrow point: `is_aromatic_atom[]` is set per-atom via
`RingAt(ri).atom_indices` walk, NOT via residue's `is_aromatic`
flag. Adding a Pro Ring would mark Pro's N, Cα, Cβ, Cγ, Cδ atoms as
`is_aromatic_atom[ai] = true`, which would change the source
classification for those 5 atoms per Pro residue from sidechain to
aromatic in the Coulomb decomposition.

This is a Bundle C compatibility constraint, not a calculator
update. Two options surface:

- **Option 1: filter at consumer.** Coulomb consumers read
  `Ring::Aromaticity()` and skip non-aromatic rings: `if
  (RingAt(ri).Aromaticity() != RingAromaticity::None)` before the
  inner walk. This is one line, not a calculator-physics change. It
  arguably IS a calculator-side update though — out of Bundle C
  scope per the locked decision.
- **Option 2: filter at consumer via type_index.** `if
  (RingAt(ri).type_index != RingTypeIndex::ProPyrrolidine)` is also
  one line. Same out-of-scope concern.
- **Option 3: substrate gates the iteration via `is_aromatic`
  RingMembership field.** Substrate already records
  `RingMembership::aromatic` (`src/SemanticEnums.h:604`); a
  Bundle-C-substrate-side filter could be Coulomb consumers walking
  `LegacyAmber().SemanticAt(ai).ring_position.primary.aromatic`
  instead of `RingAt(ri).atom_indices`. This is also a calculator
  change.

The locked Bundle C decision says "calculator outputs unchanged" and
"Pro Ring's `Intensity = 0` produces zero contribution to ring-current
shielding regardless." For Coulomb the Pro contribution is NOT zero
— it's a categorical mask change that would re-route 5 atoms per
Pro residue from sidechain to aromatic source bucket in the per-atom
Coulomb decomposition.

This is an **open question for Bundle C**. The substrate-side
ConstructRingsFromSubstrate produces a Pro Ring object;
CoulombResult would naturally pick it up via `RingCount()` /
`RingAt(ri).atom_indices` walk and reclassify Pro atoms. There are
three paths out:

(a) Accept the reclassification (Bundle C lands; Coulomb's per-atom
    aromatic-source decomposition shifts for Pro residues; bless
    NPY drift documented but not bit-identity).
(b) Add a `Ring::Aromaticity()` skip at Coulomb's mark-aromatic walk
    (1-line calculator change, out of Bundle C scope per locked
    decision).
(c) Defer Pro Ring construction in Bundle C (build it in a later
    slice with the Coulomb update bundled).

This is documented for surfacing — the audit deliberately does not
recommend.

### 5.4 RingGeometry fields consumed

This is conformation-dependent geometry, populated by
`GeometryResult::Compute` (`GeometryResult.cpp:17-19`) at calculator
prep time, before any ring-current calculator runs. Calculator
reads:

| Field                       | BS | HM | RingChi | PiQuad | Disp |
| --------------------------- | --- | --- | ------- | ------ | ---- |
| `RingGeometry::vertices`    | yes (179, 222) | yes (282) | no | yes (157) | yes (276) |
| `RingGeometry::center`      | yes (198, 251, 286, 369) | yes (236, 282) | yes (152-155) | yes (157) | yes (224, 315) |
| `RingGeometry::normal`      | yes (180, 222, 235, 255, 287) | yes (285, 313) | yes (155, 195) | yes (157) | yes (220, 319) |
| `RingGeometry::radius`      | yes (203) | yes (232) | yes (160) | yes (160) | yes (229) |

All five calculators read all four geometry fields. Bundle C does
not change `RingGeometry`. `GeometryResult` rebuilds it per frame.

### 5.5 Per-RingTypeIndex array contracts

Five fields on `ConformationAtom` (`src/ConformationAtom.h:131-134, 244-249`) sized as
`std::array<double, 8>` and `std::array<std::array<double, 5>, 8>`:

```cpp
std::array<double, 8> per_type_G_T0_sum = {};
std::array<std::array<double, 5>, 8> per_type_G_T2_sum = {};
std::array<double, 8> per_type_hm_T0_sum = {};
std::array<std::array<double, 5>, 8> per_type_hm_T2_sum = {};
std::array<double, 8> per_type_pq_scalar_sum = {};
std::array<std::array<double, 5>, 8> per_type_pq_T2_sum = {};
std::array<double, 8> per_type_disp_scalar_sum = {};
std::array<std::array<double, 5>, 8> per_type_disp_T2_sum = {};
```

These are sized as 8 (the current `RingTypeIndex::Count`). The
locked Bundle C decision is to NOT extend these to 9; Pro's
`type_index = ProPyrrolidine` (= 8) falls outside the existing
`if (ti < 8)` guards (`BiotSavart:303`, `HaighMallion:337`,
`PiQuadrupole:214`, `Dispersion:337`), so Pro contributions land in
total-G accumulators (line 298, line 343 etc.) but NOT in
`per_type_*` arrays. This is the explicit Bundle C scope decision.

Calculator-side migration to 9-arrays is a future per-calculator
slice, not Bundle C.


## 6. TRP perimeter representation status

### 6.1 Current state

The TRP perimeter (the 9-atom indole π-current loop) is materialised
at runtime by `DetectAromaticRings`. The trigger is the
`AminoAcidType::rings[]` table at `src/AminoAcidType.cpp:188`:

```cpp
{RingTypeIndex::TrpPerimeter, {"CG","CD1","NE1","CE2","CZ2","CH2","CZ3","CE3","CD2"}}
```

(Plus the TrpBenzene and TrpPyrrole entries on the lines above.)

`DetectAromaticRings` walks `aatype.rings[]` for each TRP residue
and produces THREE Ring instances per TRP: `TrpBenzeneRing` (6 atoms),
`TrpPyrroleRing` (5 atoms), `IndolePerimeterRing` (9 atoms — same
atom set as the union of the prior two minus the shared edge). The
post-pass at lines 638-645 sets `fused_partner_index` ONLY between
benzene and pyrrole; the perimeter is NOT linked to either via
`fused_partner_index`.

Net: every TRP residue produces 3 Ring objects in `protein.rings_`.
Calculators iterate all 3 and the per-type accumulators record
per-type contributions separately. The 5-ring + 6-ring + 9-ring
sums are physically meaningful (PATTERNS.md §24: ratios verified
unity within `RingBondedExclusionFilter`).

### 6.2 Substrate state

The substrate at `src/SemanticEnums.h:496-504` defines:

```cpp
RingSystemKind::Indole_Trp_5    ///< Trp pyrrole 5-ring
RingSystemKind::Indole_Trp_6    ///< Trp benzene 6-ring fused with pyrrole
```

There is NO `RingSystemKind::Indole_Trp_9`.

Per `src/SemanticEnums.h:516-518`, fused atoms (TRP CD2, CE2)
populate `primary` (smaller ring per the convention) and `secondary`
(larger ring) — the perimeter atoms are individually addressable via
`primary.ring + secondary.ring` membership checks.

### 6.3 Bundle C target (per locked decision)

The locked decision is "substrate-extension as strongly recommended
default." This means:

1. Add `RingSystemKind::Indole_Trp_9` to the enum.
2. Generate substrate to populate `Indole_Trp_9` membership on
   the 9 perimeter atoms (CG, CD1, NE1, CE2, CZ2, CH2, CZ3, CE3, CD2).
   For atoms in two rings + perimeter (CD2, CE2): primary remains
   the 5-ring per the smaller-ring convention; secondary becomes one
   of {6-ring, 9-ring}; or extend `RingPosition` to support
   tertiary/list of memberships. This is a substrate-side schema
   question.
3. `ConstructRingsFromSubstrate` reads the substrate and builds the
   `IndolePerimeterRing` instance with atom_indices in the correct
   cyclic walk order.

The simplest substrate extension preserves the current 9-atom
ordering. The cyclic walk
`CG → CD1 → NE1 → CE2 → CZ2 → CH2 → CZ3 → CE3 → CD2` traces the
indole perimeter; this is the order required to match today's
shielding NPY values. The substrate's `RingPositionLabel` walk
order for the perimeter must match. Substrate generator extension
plus a generator re-run is needed.

**Diff:** today's runtime synthesis becomes substrate-encoded.
Calculator behaviour is identical (calculator iterates 3 TRP rings
either way; Pro pyrrolidine adds a 4th ring object at TRP-adjacent
residues that's a Pro ring of a different residue, unrelated to
TRP).

The alternative — runtime synthesis at
`ConstructRingsFromSubstrate` from the union of `Indole_Trp_5` +
`Indole_Trp_6` atoms — keeps the substrate identical and pushes the
9-atom logic into Bundle C runtime code. The locked decision favours
substrate-extension; this audit documents the alternative for
completeness only.


## 7. DispersionResult local reimplementation finding

`DispersionResult.cpp:136-149` and `DispersionResult.cpp:202-204`
implement a function `BondedToVertices` and a per-ring set-vector
`std::vector<std::set<size_t>> ring_bonded(n_rings)` that mirrors
exactly what `RingBondedExclusionFilter` does at
`src/KernelEvaluationFilter.cpp:22-41`.

### 7.1 Side-by-side inspection

DispersionResult local `BondedToVertices` (lines 136-149):

```cpp
static std::set<size_t> BondedToVertices(
        const Ring& ring, const Protein& protein) {
    std::set<size_t> bonded;
    for (size_t vi : ring.atom_indices) {
        bonded.insert(vi);
        const auto& atom = protein.AtomAt(vi);
        for (size_t bi : atom.bond_indices) {
            const auto& bond = protein.BondAt(bi);
            bonded.insert(bond.atom_index_a);
            bonded.insert(bond.atom_index_b);
        }
    }
    return bonded;
}
```

`RingBondedExclusionFilter` constructor body (lines 22-41):

```cpp
RingBondedExclusionFilter::RingBondedExclusionFilter(const Protein& protein) {
    size_t n_rings = protein.RingCount();
    ring_bonded_.resize(n_rings);

    for (size_t ri = 0; ri < n_rings; ++ri) {
        const Ring& ring = protein.RingAt(ri);
        auto& bonded = ring_bonded_[ri];

        for (size_t vi : ring.atom_indices) {
            bonded.insert(vi);
            const auto& atom = protein.AtomAt(vi);
            for (size_t bi : atom.bond_indices) {
                const auto& bond = protein.BondAt(bi);
                bonded.insert(bond.atom_index_a);
                bonded.insert(bond.atom_index_b);
            }
        }
    }
}
```

The walks are byte-equivalent except for outer-loop framing (the
filter loops over rings at construction; the helper is called
per-ring by the calculator).

### 7.2 Use site in DispersionResult

DispersionResult lines 247-258:

```cpp
if (ring_bonded[ri].count(ai)) {
    // ---- GeometryChoice: through-bond exclusion ----
    choices.Record(...);
    bonded_exclusions++;
    continue;
}
```

Compare to other ring-current calculators using the filter:
`BiotSavart:151` adds `RingBondedExclusionFilter` to the filter set
and then `filters.AcceptAll(ctx)` at line 206 returns false when
ai is in any ring's exclusion set. Same physics, different surface.

### 7.3 Verdict — accidental duplication, not science divergence

Both implementations:
- Walk `Ring::atom_indices` to seed each ring's exclusion set.
- Walk per-vertex `bond_indices` to add bonded neighbours.
- Read the same `Protein::BondAt(bi)` for `atom_index_a`/`b`.

The walks return identical sets (modulo the outer-loop framing).

DispersionResult's GeometryChoice at line 248-255 records
"through-bond exclusion" with role explicitly named `"ring_bonded"`,
whereas `KernelFilterSet::AcceptAll` records the rejection with the
filter's `LastRejectorName()` (a `const char*` that is
`"RingBondedExclusion"`). The strings differ; the topology decision
is identical.

DispersionResult does NOT add `RingBondedExclusionFilter` to its
`KernelFilterSet`; it adds only `DipolarNearFieldFilter`
(line 189). The local helper exists because the per-vertex
iteration at line 275 needs the exclusion set referenced by index
into `ring_bonded[ri]` rather than as a per-call kernel context.

This is a refactoring opportunity for a future calculator-side
slice (out of Bundle C scope per the locked decision). The
substrate-side investigation question is resolved: the divergence is
implementation, not science.


## 8. Implications for Bundle C (substrate-side only)

### 8.1 Ring class API requirements for ProPyrrolidineRing

ProPyrrolidineRing must implement all 7 base virtuals plus inherit
intermediate-class virtual via `FiveMemberedRing`. Required
implementations:

| Virtual                       | Pro return value          | Rationale                                               |
| ----------------------------- | ------------------------- | ------------------------------------------------------- |
| `Intensity()`                 | `0.0` literal             | Saturated heterocycle; no π current. Joule & Mills 2010 ch. 7 |
| `LiteratureIntensity()`       | `0.0` literal             | Same                                                    |
| `JBLobeOffset()`              | `0.0` literal             | No JB loop physics; zero is the unambiguous null value |
| `NitrogenCount()`             | `1`                       | Pro N is in the ring                                    |
| `Aromaticity()`               | `RingAromaticity::None`   | New enum value (Bundle C extends `RingAromaticity` enum) |
| `RingSizeValue()`             | `5` (inherited)           | FiveMemberedRing default                                |
| `TypeName()`                  | `"PRO"`                   | Short diagnostic name                                   |

Required type-system extensions:

1. `RingTypeIndex::ProPyrrolidine = 8`, `Count = 9` (`Types.h:175-185`).
2. `RingAromaticity::None` enum value (`Types.h:168`).
3. `RingTypeName(ProPyrrolidine)` switch case returning `"PRO"`.
4. `CreateRing(ProPyrrolidine)` factory case in `Ring.cpp:50-62`
   returning `std::make_unique<ProPyrrolidineRing>()`.

### 8.2 ConstructRingsFromSubstrate input/output spec

**Inputs:**

- `LegacyAmber().SemanticAt(ai).ring_position` per `ai` in
  `protein.atoms_`. Specifically:
  - `primary.ring` (RingSystemKind) — per-atom primary ring
    membership.
  - `primary.position` (RingPositionLabel) — per-atom position
    within primary ring.
  - `secondary.ring`, `secondary.position` — for fused-bridge
    atoms.
- `protein.residues_` for parent-residue lookup.
- (The protonation_variant_index for HID/HIE/HIP discrimination is
  resolved before this function runs.)

**Output:**

- Populates `protein.rings_` with one `Ring*` per detected ring.
- Sets `Ring::atom_indices` in canonical cyclic walk order.
- Sets `Ring::type_index`, `Ring::parent_residue_index`,
  `Ring::parent_residue_number`.
- Post-pass sets `Ring::fused_partner_index` for TRP fused
  benzene+pyrrole pair (and perimeter if substrate-extended).

**Algorithm sketch (substrate-driven):**

```text
For each residue ri in residues_:
    Group atoms by RingSystemKind (primary + secondary memberships).
    For each (ring_system_kind, atom_set):
        Determine RingTypeIndex (from RingSystemKind + protonation_variant_index for HIS).
        Sort atoms by canonical RingPositionLabel walk-order for that ring_system_kind.
        Construct Ring via CreateRing(ring_type).
        Populate atom_indices, parent_residue_index, parent_residue_number.
        Append to rings_.
        If TRP, track for fused-partner post-pass.
        Record GeometryChoice (per "always record" decision).
For each TRP fused ring pair (benzene + pyrrole):
    Set fused_partner_index on both.
For each Pro residue:
    Construct ProPyrrolidineRing with atom_indices = [N, Cα, Cβ, Cγ, Cδ].
    Append to rings_.
    Record GeometryChoice.
```

The canonical walk-order per RingSystemKind is encoded in the
RingPositionLabel sequencing logic of the function. Substrate
emits per-atom labels; the function consumes them.

### 8.3 Compatibility constraints

**Constraint A — atom_indices cyclic walk preservation.** For
PHE/TYR/HIS/HID/HIE/HIP/TrpBenzene/TrpPyrrole/TrpPerimeter, the
substrate-driven walk must produce `atom_indices` in an order whose
first-three-atoms cross product matches today's order's first-three-
atoms cross product (sign agreement). Bless-compare on shielding
NPYs is the test. Verify per-ring-type at fixture level.

**Constraint B — RingTypeIndex::ProPyrrolidine = 8 (NOT lower).**
Bundle C scope decision: per-type `[8]` arrays remain size 8; Pro's
type_index 8 falls outside `if (ti < 8)` guards on purpose. Pro's
contribution does NOT enter `per_type_*` accumulators. (If
ProPyrrolidine took index 7 or earlier, it would shift other
ring types' assignments and break bit-identity.)

**Constraint C — Pro Ring atom_indices ordering N→Cα→Cβ→Cγ→Cδ.**
Per locked decision. Future puckering descriptors (planned
calculators consuming Pro Ring) require the dihedral
N-Cα-Cβ-Cγ-Cδ in this fixed walk order.

**Constraint D — `is_aromatic_atom` mask in Coulomb consumers.**
Current behaviour walks `RingAt(ri).atom_indices` for ALL ring
types. Adding Pro Ring would mark Pro N, Cα, Cβ, Cγ, Cδ as
aromatic-source. This is the open finding from §5.3 — needs
explicit Bundle C decision (accept reclassification, defer Pro
construction, or out-of-scope calculator filter).

**Constraint E — `accumulated.fused_partner_index` for TRP perimeter.**
If substrate-extension extends `RingSystemKind::Indole_Trp_9`,
ConstructRingsFromSubstrate must establish three-way fused
relationships (benzene ↔ pyrrole ↔ perimeter). Today's
`fused_partner_index` is a single SIZE_MAX-default index. Either
extend the field to a small vector, or maintain only the binary
benzene↔pyrrole link as today and let perimeter stand
separately (current behaviour).

**Constraint F — fail-loud on substrate-miss.** Per README §11. If
a residue is `AminoAcid::HIS` with no resolved `protonation_variant_index`
when ConstructRingsFromSubstrate runs, FATAL+abort, do NOT default
to `HisImidazole`. If a residue's substrate has no `ring_position`
on atoms expected to be ring members, FATAL+abort, do NOT silently
skip. Stronger invariant than today's `if (!all_present) continue`.

### 8.4 GeometryChoice records expected during ConstructRingsFromSubstrate

Per the locked "always record" decision. Each ring constructed
emits a GeometryChoice:

```cpp
choices.Record(CalculatorId::ProteinFinalize, ri, "ring constructed",
    [&ring, /* parent_residue ref */](GeometryChoice& gc) {
        AddRing(gc, ring.get(), EntityRole::Source, EntityOutcome::Included);
        AddNumber(gc, "ring_size", static_cast<double>(ring->RingSizeValue()), "atoms");
        AddNumber(gc, "intensity_literature", ring->LiteratureIntensity(), "nA");
        AddNumber(gc, "n_heteroatoms", static_cast<double>(ring->NitrogenCount()), "atoms");
        // Aromaticity recorded as a string-tag via SetTextField if available; else
        // map to integer-via-enum-cast.
    });
```

For Pro rings, the record is identical with intensity_literature =
0.0, NitrogenCount = 1, RingSize = 5. Per the locked decision
"always record"; the Pro entries are documentary value, not noise.

A second record category for skipped rings (residues whose substrate
has no ring chemistry but expected one — e.g. malformed input):

```cpp
choices.Record(CalculatorId::ProteinFinalize, ri, "no ring substrate",
    [&residue](GeometryChoice& gc) {
        // Document the residue, the failure mode, and the FATAL trigger
    });
```

(The fail-loud rule means this branch typically doesn't return; the
record is emitted before abort.)


## 9. Open questions

These remain for Bundle C drafting:

### Q1. `CovalentTopology::Resolve(rings_)` actual dependency on rings

The current `FinalizeConstruction` at line 232 passes `rings_` to
bond resolution. To move ring construction post-substrate-composition
(which Bundle C requires), bond resolution must either run without
rings (preferred) or be split into a pre-rings pass (covalent edges)
plus a post-rings pass (categorisation refinement). This audit did
not deep-dive `CovalentTopology::Resolve`; that's a quick follow-up
investigation.

### Q2. Coulomb consumer reclassification on Pro Ring (§5.3)

Three paths: accept the per-atom aromatic-mask shift on Pro residues,
defer Pro Ring construction, or update Coulomb consumers (out of
Bundle C scope). User decides.

### Q3. TRP perimeter substrate schema for atoms in 3 rings

Currently `RingPosition` is `{primary, secondary}`. With Indole_Trp_9
substrate extension, CD2 and CE2 are in 3 rings. Schema options:
(a) keep CD2/CE2 with primary=Indole_Trp_5, secondary=Indole_Trp_6,
omit perimeter membership from substrate but synthesise the
perimeter Ring at runtime from union of 5+6 atoms; (b) extend
`RingPosition` to `{primary, secondary, tertiary}`; (c) change
`RingPosition` to `std::vector<RingMembership>`. (a) keeps
substrate-extension minimal; (c) is most flexible. Bundle C scope
is to choose.

### Q4. ProPyrrolidineRing values: literal 0.0 or CalculatorConfig?

The audit recommends literal `0.0` for Intensity / LiteratureIntensity
/ JBLobeOffset on Pro because saturated heterocycles don't carry π
current as a matter of physics (Joule & Mills 2010 ch. 7), not as a
calibration parameter. Other rings use CalculatorConfig::Get because
their values are fitted. Bundle C should commit explicitly.

### Q5. RingPositionLabel Pro extensions

Per the README addition (Pro per-atom RingPositionLabel
enrichment): five Pro-specific labels (`ProRingNitrogen`,
`ProRingAlphaCarbon`, `ProRingBeta`, `ProRingPuckerPivot`,
`ProRingDelta`). Verified against current code: these enum values
are NOT yet present in `src/SemanticEnums.h:528-589`. Pro atoms
currently emit `RingPositionLabel::Saturated` for all 5 ring atoms
(`src/generated/LegacyAmberSemanticTables.cpp:323-329` — all show
either `Saturated` or `Heteroatom_NoH` for the N). Bundle C must
extend the enum and re-run the substrate generator.

### Q6. fused_partner_index extension for perimeter

Today `Ring::fused_partner_index` is a single SIZE_MAX index pointing
to the partner ring. With perimeter as a substrate-encoded ring, the
benzene↔pyrrole↔perimeter triangle does not fit. Options:
keep binary linkage, extend to small-vector, or leave perimeter
unlinked (today's behaviour — perimeter has fused_partner_index =
SIZE_MAX). Bundle C should commit.
