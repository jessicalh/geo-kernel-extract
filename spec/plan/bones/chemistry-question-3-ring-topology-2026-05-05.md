# Chemistry Question 3 — Ring topology (2026-05-05)

## Summary

Ring detection in this codebase is currently a string-keyed pass:
`Protein::DetectAromaticRings` (`src/Protein.cpp:431-529`) builds a
per-residue `name → atom_index` map by reading `Atom.pdb_atom_name`,
then matches that map against `const char*` literals stored in
`AminoAcidType::rings[].atom_names` (`src/AminoAcidType.cpp:92, 147,
186-188, 201`). The `Ring` object that calculators consume
(`src/Ring.h`) carries `atom_indices`, `RingTypeIndex`,
`parent_residue_index`, `parent_residue_number`, and
`fused_partner_index` — none of those fields demand a string lookup;
the string lookup is purely the *populator* of the index list.

Meanwhile, the chemistry-substrate work that landed at commit
`ee1f1b4` (per topology-encoding-dependencies §H.9) emits typed
`RingSystemKind` and `RingPositionLabel` per-atom membership in
the generated `AtomSemanticTable` rows. No production calculator
reads those fields yet. The substrate has the typed answer for
ring membership; `DetectAromaticRings` is reaching past it to
strings that pre-date its existence.

The recommended Phase-1 architecture (§4-5 below) is **substrate-
driven ring construction**: after `LookupBy` populates per-atom
substrate rows at `FinalizeConstruction`, a new
`ConstructRingsFromSubstrate()` pass reads `ring_position.primary`
and `ring_position.secondary`, groups atoms by `(residue index,
RingSystemKind)`, and builds the same `Ring` objects calculators
already consume. `Protein::DetectAromaticRings` dissolves into
substrate population. `AminoAcidType::rings[].atom_names` becomes
unreachable from the runtime and is removed; ring atom membership
is encoded once, in the generated substrate table, populated by the
generator from CCD + RDKit.

No calculator behaviour changes in Phase 1. The same `Ring`
objects appear, with the same `atom_indices` ordering preserved
(or canonicalised in a deliberate way — see §7). Phase 2 (§6)
opens up new direct-from-substrate calculator paths
(per-position ring-current stratification, ring-position-aware
diagnostics).

A non-trivial dependency: substrate population must run at
`FinalizeConstruction` for any of this to work. Per the audit's
closing note 2, the LookupBy-based populator is not yet wired
into FinalizeConstruction; this work and that work are linked.

---

## 1. The chemistry question

**How does the model identify which atoms belong to which ring within a residue?**

### 1.1 The rings, per-residue

The standard 20 amino acids contribute six ring topologies, all
sidechain-bound:

| Residue | Ring | Size | Atoms (canonical AMBER names) | Aromaticity |
|---|---|---|---|---|
| PHE | benzene | 6 | CG, CD1, CE1, CZ, CE2, CD2 | full |
| TYR | phenol (benzene + para-OH) | 6 | CG, CD1, CE1, CZ, CE2, CD2 | full |
| HIS | imidazole | 5 | CG, ND1, CE1, NE2, CD2 | weak (variant-dependent) |
| TRP | pyrrole (5-ring of indole) | 5 | CG, CD1, NE1, CE2, CD2 | reduced |
| TRP | benzene (6-ring of indole) | 6 | CD2, CE2, CZ2, CH2, CZ3, CE3 | full |
| TRP | indole perimeter (fused) | 9 | CG, CD1, NE1, CE2, CZ2, CH2, CZ3, CE3, CD2 | full |
| PRO | pyrrolidine | 5 | N, CA, CB, CG, CD | none (saturated) |

TRP carries three rings simultaneously: the 5-ring, the 6-ring,
and a synthetic 9-atom perimeter ring used by some ring-current
formulations (Case 1995). The perimeter is not a chemical ring —
it's the atom set traversed by the conjugated π current that
encircles the indole — but it lives in the same data structure
because calculators want it the same way.

PRO is included for completeness. It is **not** in the current
`AminoAcidType::rings[]` table. Its pyrrolidine ring is invisible
to `DetectAromaticRings` (the function early-returns on
`!aatype.is_aromatic`). The substrate's `RingSystemKind::
Pyrrolidine_Pro` and `RingPositionLabel::Saturated` exist for it,
but no `Ring` object is constructed. Whether it should be is one
of the §7 user-decision points.

### 1.2 What calculators use ring information

Five production calculators consume `Protein::RingAt(ri)` /
`RingCount()`:

- **BiotSavartResult** (`src/BiotSavartResult.cpp`) — ring-current
  shielding via Biot-Savart-Johnston-Bovey two-loop construction.
  Uses `ring.atom_indices` to compute centre, normal, and the
  loops; `ring.Intensity()`, `ring.JBLobeOffset()`,
  `ring.TypeIndexAsInt()`.
- **HaighMallionResult** (`src/HaighMallionResult.cpp`) —
  Haigh-Mallion ring-current shielding. Uses
  `ring.atom_indices`, `ring.type_index`, `ring.TypeIndexAsInt()`.
- **RingSusceptibilityResult** (`src/RingSusceptibilityResult.cpp`)
  — Pople point-dipole ring-current shielding.
  Uses `ring.atom_indices`, `ring.type_index`.
- **DispersionResult** (`src/DispersionResult.cpp`) — uses ring
  atoms as dispersion-source points (via `ring.atom_indices`).
- **EnrichmentResult** / **MolecularGraphResult** /
  **CoulombResult** / **MopacCoulombResult** /
  **AIMNet2Result** / **CovalentTopology** — incidental ring
  iterations for atom-role tagging, neighbour graphs, etc.

Plus **KernelEvaluationFilter** which traverses
`Protein → RingAt(ri) → atom_indices → AtomAt(vi) →
bond_indices` to scope geometry choices.

What none of them currently consume: the typed RingSystemKind /
RingPositionLabel substrate. Calculators ask `ring.Intensity()`,
which uses `RingTypeIndex` indirection to a CalculatorConfig
parameter. Per-position discrimination (ipso vs ortho vs meta vs
para — different chemistry, very different observed shifts) is
not possible from the current `Ring` interface. Phase 2 §6 picks
this up.

### 1.3 Why the question matters

Aromatic rings are one of the dominant signals in protein NMR
shielding. Ring-current shifts for atoms near phenylalanine /
tyrosine / tryptophan rings can run 1-3 ppm for ¹H, larger for
¹³C — comparable to or larger than the secondary-structure
contribution. The thesis's calibration depends on ring detection
working correctly across every protein in the fleet.

The current implementation works, on standard fixtures. It works
because the atom-name strings produced by the loaders agree with
the atom-name literals in the AminoAcidType table. That agreement
is brittle: it depends on every loader producing PDB-canonical
names (NamingRegistry's job), the strings round-tripping through
H5 emit/read paths, and the `AminoAcidType` table tracking any
naming-convention drift (e.g. a hypothetical ff15 update). When
the agreement breaks, ring detection produces zero rings for the
affected residue and `BiotSavartResult` etc. silently produce
zero ring-current contribution at every nearby probe — no error,
no warning, just absent kernel signal.

The risk is not theoretical. The audit calls this out at Hotspot
3 (§2.3 of `mechanical-identity-model-and-audit-2026-05-05.md`):
"A miss means zero ring current, silently."

---

## 2. What the current code does

### 2.1 `Protein::DetectAromaticRings`

Implementation, reading from `src/Protein.cpp:431-529` (the lines
the audit calls Hotspot 3):

```
for ri in 0..residues_.size():
    res = residues_[ri]
    aatype = res.AminoAcidInfo()
    if (!aatype.is_aromatic) continue          # <-- skips PRO entirely

    name_to_idx = {}
    for ai in res.atom_indices:
        name_to_idx[atoms_[ai]->pdb_atom_name] = ai

    for ring_def in aatype.rings:               # AminoAcidRing entries
        atom_indices = []
        for aname in ring_def.atom_names:        # const char* literals
            if name_to_idx has aname:
                atom_indices.push_back(name_to_idx[aname])
            else:
                all_present = false; break
        if !all_present: continue

        effective_type = ring_def.type_index
        if res.type == HIS:
            # HIS variant resolution
            if res.protonation_variant_index >= 0:
                # typed path: 0 -> HID, 1 -> HIE, 2 -> HIP
            else:
                # PDB LOADING BOUNDARY string fallback:
                has_HD1 = name_to_idx.find("HD1") != end
                has_HE2 = name_to_idx.find("HE2") != end
                if has_HD1 && has_HE2: HIP
                elif has_HD1: HID
                elif has_HE2: HIE

        ring = CreateRing(effective_type)
        ring->atom_indices = atom_indices
        ring->parent_residue_index = ri
        ring->parent_residue_number = res.sequence_number
        rings_.push_back(ring)
```

Plus a TRP fused-pair fixup pass at the end that links the
benzene ring to its pyrrole partner via `fused_partner_index`.

**What this code is doing chemistry-wise:** for each residue
flagged `is_aromatic`, look up the ring atom set defined in
`AminoAcidType`, find the protein's atoms by string matching
against `pdb_atom_name`, build the `Ring` object. The ring
*chemistry* — which atoms participate, in what order, with what
ring system, position labels — is encoded entirely in:

1. `AminoAcidType::is_aromatic` (a bool, residue-level filter).
2. `AminoAcidType::rings[].type_index` (typed `RingTypeIndex`).
3. `AminoAcidType::rings[].atom_names` (`const char*` array).

Items 1 and 2 are typed. Item 3 is the string surface this audit
pivot is removing from the codebase.

The function is run from `Protein::FinalizeConstruction` at line
220, immediately after `CacheResidueBackboneIndices` and
`ResolveProtonationStates`. Same ordering on every load path
(PDB, ORCA, AMBER trajectory).

### 2.2 The `Ring` object

`src/Ring.h` declares a small inheritance hierarchy, eight types
in three size categories (`SixMemberedRing`, `FiveMemberedRing`,
`FusedRing`). The base class fields are:

```cpp
class Ring {
public:
    std::vector<size_t> atom_indices;        // populated by DetectAromaticRings
    RingTypeIndex       type_index;          // typed, set in subclass ctor
    size_t              parent_residue_index;
    int                 parent_residue_number;
    size_t              fused_partner_index;  // SIZE_MAX if not fused

    virtual double Intensity() const = 0;
    virtual double LiteratureIntensity() const = 0;
    virtual double JBLobeOffset() const = 0;
    virtual int NitrogenCount() const = 0;
    virtual RingAromaticity Aromaticity() const = 0;
    virtual int RingSizeValue() const = 0;
    virtual const char* TypeName() const = 0;

    bool IsFused() const;
    int TypeIndexAsInt() const;
    RingGeometry ComputeGeometry(const std::vector<Vec3>& positions) const;

    RingAccumulated accumulated;             // mutable: filled by passes
};
```

Subclasses override the virtual properties to return per-type
constants (PHE intensity = -12.0, TYR = -11.28, etc., parameterised
through `CalculatorConfig::Get(...)` for thesis calibration).

**Observation:** the `Ring` object is *already* the right shape.
`atom_indices` is just `vector<size_t>` — populating it from
typed substrate is a different code path producing the same
data structure. Calculators do not need to change.

What the `Ring` object does NOT carry today: per-atom
RingPositionLabel within the ring. There is no
`ring.PositionLabelOfAtom(ai)` accessor. That information is
present in the substrate but not surfaced through the `Ring`.
Phase 2 §6 picks this up.

### 2.3 `AminoAcidType::rings[]`

The string-shape data lives in `src/AminoAcidType.h:57-61`:

```cpp
struct AminoAcidRing {
    RingTypeIndex            type_index;     // TYPED
    std::vector<const char*> atom_names;     // STRINGS
};
```

and is populated for the four aromatic residues:

```
PHE: {RingTypeIndex::PheBenzene,   {"CG","CD1","CE1","CZ","CE2","CD2"}}   (1 ring)
TYR: {RingTypeIndex::TyrPhenol,    {"CG","CD1","CE1","CZ","CE2","CD2"}}   (1 ring)
HIS: {RingTypeIndex::HisImidazole, {"CG","ND1","CE1","NE2","CD2"}}        (1 ring; variant resolved at detection)
TRP: {RingTypeIndex::TrpBenzene,   {"CD2","CE2","CZ2","CH2","CZ3","CE3"}}
     {RingTypeIndex::TrpPyrrole,   {"CG","CD1","NE1","CE2","CD2"}}
     {RingTypeIndex::TrpPerimeter, {"CG","CD1","NE1","CE2","CZ2","CH2","CZ3","CE3","CD2"}}   (3 rings)
```

This is a hybrid: typed `type_index`, string `atom_names`. The
`type_index` is what the calculator dispatches on; the
`atom_names` only purpose is to be string-matched against
`pdb_atom_name` to populate the index list. The ordering of
`atom_names` *is* load-bearing — the SVD-based ring-normal
calculation in `Ring::ComputeGeometry` (`src/Ring.cpp:6-47`)
uses the cross product of "first two edges" to fix orientation,
which means consistent vertex ordering matters for sign
conventions on the normal vector. Whatever replaces this
string list has to preserve a defined ordering.

PRO is absent: `is_aromatic = false` in PRO's row, no `rings[]`
entry. The pyrrolidine ring exists chemically but not in this
table.

---

## 3. Where typed information already lives

### 3.1 `AminoAcidType::rings[].ring_type_index` (typed)

The `RingTypeIndex` enum (`PheBenzene / TyrPhenol / TrpBenzene /
TrpPyrrole / TrpPerimeter / HisImidazole / HidImidazole /
HieImidazole`) is the dispatch axis the calculators use. It is
typed already and survives any pivot.

### 3.2 Per-atom RingSystemKind / RingPositionLabel substrate

In `src/SemanticEnums.h`:

```cpp
enum class RingSystemKind : uint8_t {
    NotInRing       = 0,
    Benzene_Phe     = 1,    // Phe sidechain six-ring aromatic
    Benzene_Tyr     = 2,    // Tyr sidechain six-ring aromatic
    Imidazole_His   = 3,    // His imidazole; variant-dependent aromaticity
    Indole_Trp_5    = 4,    // Trp pyrrole 5-ring
    Indole_Trp_6    = 5,    // Trp benzene 6-ring fused
    Pyrrolidine_Pro = 6,    // Pro saturated 5-ring
};

struct RingMembership {
    RingSystemKind    ring;
    RingPositionLabel position;
    uint8_t           ring_size;
    bool              aromatic;
    bool              planar;
    uint8_t           n_heteroatoms;
};

struct RingPosition {
    RingMembership primary;     // primary ring (smaller-ring rule for fused atoms)
    RingMembership secondary;   // fused partner; NotInRing when single-ring atom
};
```

Each `AtomSemanticTable` row carries `RingPosition ring_position`
populated by the table generator from CCD topology + RDKit
aromaticity perception + the synthesised Trp-perimeter convention
(see `topology-encoding-dependencies-2026-05-05.md` §A.3 and §H).
The substrate already encodes:

- Which ring system the atom belongs to (`ring_position.primary.ring`)
- Per-position chemistry label (Ipso / Ortho1 / Ortho2 / Meta1 /
  Meta2 / Para / PyrroleAlpha / PyrroleBeta / BridgeFusion /
  Heteroatom_NH / Heteroatom_NoH / Heteroatom_OH / Saturated)
- Ring size, aromatic flag, planar flag, heteroatom count
- Fused-ring secondary membership (Trp CD2 / CE2 are bridgeheads
  with primary = 5-ring, secondary = 6-ring)

The granularity is finer than what `AminoAcidType::rings[]`
encodes. The substrate distinguishes Tyr's para-Cζ (Para) from
Phe's para-Cζ (Para in `Benzene_Phe`); both are "CZ" string-wise.
Ring-position-stratified shielding statistics (which the thesis
methodology calls for, given the ~30 ppm ¹³C spread between Phe
ortho/meta/para) are blocked behind that distinction today.

PRO Cα through Cδ + N are encoded with `Pyrrolidine_Pro/Saturated/5`
in the substrate. There is no corresponding `Ring` object for
them.

### 3.3 `LookupBy` runtime lookup

The generator emits `LookupBy(residue, variant_idx, identity)`
returning `const AtomSemanticTable*` (per `src/generated/Legacy
AmberSemanticTables.h:39-72` and topology-encoding §H.2). One
call per atom at protein construction populates the per-atom
substrate. After this runs, every atom carries its typed
ring_position field.

**Status check (per audit closing note 2):** the generator and
LookupBy exist; the wiring at `Protein::FinalizeConstruction`
that calls LookupBy per atom and writes the substrate row onto
the runtime atom record is not yet in place. Substrate population
is a prerequisite for substrate-driven ring construction. They
land together or not at all (§5).

### 3.4 CovalentTopology bond graph (option C ingredient)

`src/CovalentTopology.h` exposes the bond list resolved from
geometry at `FinalizeConstruction`. Algorithmic ring detection
(SSSR — Smallest Set of Smallest Rings) over the bond graph
produces a residue-agnostic ring inventory. Some external
chemistry libraries do this. We do not currently have an SSSR
implementation, and adopting one introduces:

- An external ring-perception library (e.g. RDKit at runtime —
  blocked by the string-barrier rule for libnmr_shielding;
  topology-encoding §H.7) or a hand-rolled algorithm.
- A separate problem: mapping detected rings to `RingTypeIndex`
  (the calculators' dispatch axis) requires recognising which
  topology corresponds to which residue's ring. The substrate
  already encodes this relationship; SSSR would re-derive it.
- Performance: O(atoms) per protein, plus a deterministic
  ordering choice for SSSR results so vertex orderings match
  expectations downstream.

Option C is mentioned in §4 below for completeness; it is not
the recommended path.

---

## 4. Architectural choice

Three paths were considered. The architecture rule we are
evaluating against: **objects answer questions about themselves**;
**no string dispatch**; **the substrate already encodes the
chemistry that pre-substrate code re-derives**; **calculator
behaviour preserved in Phase 1**.

### Path A — Substrate-driven ring construction

After per-atom substrate population at `FinalizeConstruction`,
walk atoms; for each atom with `ring_position.primary.ring !=
NotInRing`, group by `(residue index, RingSystemKind)`; build
one `Ring` object per group. Map `RingSystemKind` → `RingTypeIndex`
deterministically, using `Residue.protonation_variant_index` to
disambiguate the three His variants (HID/HIE/HIP). Populate
`atom_indices` from the group; resolve order from the typed
substrate (the `RingPositionLabel` enum carries the canonical
ordering: Ipso, Ortho1, Meta1, Para, Meta2, Ortho2 for benzene
rings).

**Pros:**
- Reads chemistry from the typed authority (the substrate). The
  pre-substrate hybrid string/typed table is bypassed and can
  be removed.
- Per-position chemistry labels are present on every atom for
  free; Phase 2 calculators can stratify by ring position
  without further work.
- Removes the duplicate HIS-tautomer string fallback at
  `Protein.cpp:485-497`. Variant index is already on `Residue`;
  substrate variant-table selection during `LookupBy` does the
  work.
- TRP fused-partner pairing is straightforward: groups for
  `Indole_Trp_5` and `Indole_Trp_6` within one residue are
  partners by definition.
- PRO pyrrolidine ring, if we want it, falls out automatically
  (substrate carries `Pyrrolidine_Pro`); no special-case code.
- Same `Ring` object on the calculator side; no calculator
  changes.

**Cons:**
- Depends on substrate population being live at
  `FinalizeConstruction`. That work is not yet sequenced (per
  audit closing note 2). This refactor blocks on that one.
- TRP perimeter ring is *synthetic* — it's the 9-atom union of
  the 5-ring and 6-ring atom lists in a defined sequential
  order. The substrate doesn't have a `RingSystemKind::Indole_
  Trp_Perimeter` (and shouldn't — perimeter is not a chemical
  ring, it's a calculation construct). The perimeter `Ring` has
  to be built post-hoc from the union of the 5-ring and 6-ring
  groups, with explicit ordering.
- Atom ordering within a group: `RingPositionLabel` enum
  provides position labels (Ipso, Ortho1, Meta1, Para, Meta2,
  Ortho2 for benzene; PyrroleAlpha, PyrroleBeta, BridgeFusion,
  Heteroatom_NH/NoH for imidazole/pyrrole) but does not by
  itself impose a traversal order on the ring. We have to define
  one (canonical: walk in increasing position label, starting
  from Ipso) and verify it matches the current order produced
  by the string-list path so the SVD normal-vector orientation
  doesn't flip silently.

### Path B — Typed atom-list in AminoAcidType

Replace `std::vector<const char*> atom_names` with
`std::vector<AtomMechanicalIdentity> atom_identities` in
`AminoAcidRing`. `DetectAromaticRings` does typed-equality
matching using `Atom.identity` (added in audit Phase 1) instead
of `pdb_atom_name`-keyed maps.

**Pros:**
- Smallest-impact change. The `DetectAromaticRings` shape is
  preserved; only the *keying* changes from string to typed
  identity.
- Doesn't depend on substrate population; AminoAcidType is the
  authority directly.
- Consistent with the audit's Hotspot 3 recommendation when
  read literally: "encode ring atom membership as
  `AtomMechanicalIdentity` in `AminoAcidType` and resolve via
  `res.AtomWithIdentity(...)`."

**Cons:**
- Two authorities for the same chemistry: `AminoAcidType.rings`
  duplicates information already in the generated substrate
  table. The substrate has *more* (RingPositionLabel per atom);
  the new typed AminoAcidType list has *less* (just the atoms).
  This is the opposite of what the project documents say
  (PATTERNS.md "single authority for amino acid chemistry").
- Calculators still cannot ask for per-position chemistry
  labels. Phase 2 stratification is unblocked by Path A but
  not by Path B; Path B forces a second refactor to surface
  the substrate fields.
- HIS variant disambiguation still needs a special branch (the
  `protonation_variant_index → RingTypeIndex` switch) but Path
  A handles it for free at the substrate level.
- The TRP perimeter ring is still hand-listed as 9
  AtomMechanicalIdentity entries, just like it's hand-listed
  as 9 strings now. Net duplication.

### Path C — Graph-algorithm SSSR

Detect rings from `CovalentTopology` bond graph using SSSR;
classify each ring by (residue type, ring size, heteroatom
count) into `RingTypeIndex`. Build `Ring` objects from the
detected vertex sets.

**Pros:**
- Residue-agnostic. Would work for non-standard residues out
  of the box if their bond graphs are populated (modified amino
  acids, ligands).
- No data table needed for ring atom lists at all.

**Cons:**
- Requires either an external library or a hand-rolled SSSR
  algorithm. Both add complexity — SSSR is non-trivial to
  implement correctly; choosing the canonical "smallest" ring
  set has multiple definitions (Horton's, Vismara's) with
  different downstream behaviours.
- TRP perimeter ring is not minimal; SSSR will find the 5-ring
  and 6-ring, not the 9-atom perimeter. We'd have to compose
  the perimeter on top, exactly as in Path A.
- Mapping detected ring → `RingTypeIndex` requires either
  (a) re-deriving Tyr-vs-Phe distinction from para-OH presence
  (the residue type already says it directly), or (b) consulting
  the residue type, which makes SSSR redundant — we already
  knew Tyr has a phenol ring, the question was always which
  atoms.
- His variant assignment still needs the variant_index switch.
- Calculator semantics depend on consistent atom ordering;
  SSSR doesn't give that for free.
- The substrate already encodes the answer SSSR would compute,
  with strictly more information (ring position labels). Using
  SSSR is computing what we already have.

### Recommendation: Path A

Path A is the only one that:

1. Reads chemistry from the existing typed authority (the
   substrate);
2. Removes a string surface rather than relocating it;
3. Surfaces ring-position chemistry on every atom for Phase 2
   to consume directly;
4. Leaves Phase 1 calculator behaviour bit-identical (modulo
   the ordering question, which is testable).

Path B is acceptable as an *interim* if substrate population
slips its sequencing (e.g. if we want ring detection migrated
before substrate population lands). It is not the destination
shape; the destination has the substrate as authority.

Path C is rejected.

---

## 5. Phase 1 design (existing calculators must work)

The smallest refactor that lands Path A.

### 5.1 Sequencing dependency

Substrate population at `Protein::FinalizeConstruction` must
exist before substrate-driven ring construction can run. The
audit's Phase 1 (per `mechanical-identity-model-and-audit-2026
-05-05.md` §2.4) adds `Atom.identity` and the per-atom semantic
record. Once that lands, every atom carries:

- `Atom.identity` (typed AtomMechanicalIdentity)
- `Atom.semantic` (or pointer to const AtomSemanticTable row;
  contains `RingPosition ring_position`)

Substrate-driven ring construction depends on `Atom.semantic`
being populated; the audit's Phase 1 lands that. This work is
slot-in after audit Phase 1, before audit Phase 2.

### 5.2 The new construction pass

Replace `Protein::DetectAromaticRings` with `Protein::Construct
RingsFromSubstrate` (or absorb it; the function name is
secondary). Pseudocode:

```
ConstructRingsFromSubstrate():
    rings_.clear()

    # Group atoms by (residue_index, RingSystemKind, primary-or-secondary)
    map<(ri, RingSystemKind, primary_or_secondary), vector<atom_indices>> groups
    for atom_idx in 0..atoms_.size():
        rp = atoms_[atom_idx].semantic.ring_position
        if rp.primary.ring != NotInRing:
            ri = atoms_[atom_idx].residue_index
            groups[(ri, rp.primary.ring, primary)].push(atom_idx)
        if rp.secondary.ring != NotInRing:
            ri = atoms_[atom_idx].residue_index
            groups[(ri, rp.secondary.ring, secondary)].push(atom_idx)

    # Build a Ring per (residue, RingSystemKind)
    # When a group is the secondary ring of a fused atom, we still record
    # the atom in that ring's index list — bridgehead atoms appear in both.
    for (ri, ring_system, _), atom_idx_list in groups:
        ring_type = MapRingSystemToRingTypeIndex(ring_system, residues_[ri])
        ring = CreateRing(ring_type)
        ring->atom_indices = OrderedAtomIndices(atom_idx_list, ring_system, residues_[ri])
        ring->parent_residue_index = ri
        ring->parent_residue_number = residues_[ri].sequence_number
        rings_.push_back(ring)

    # Trp perimeter: synthesise from 5-ring + 6-ring atom unions.
    for ri where residues_[ri].type == TRP:
        five_ring = find ring with parent ri and type_index TrpPyrrole
        six_ring  = find ring with parent ri and type_index TrpBenzene
        if both present:
            perimeter_atoms = OrderedPerimeterFromFusedPair(five_ring, six_ring)
            perim = CreateRing(TrpPerimeter)
            perim->atom_indices = perimeter_atoms
            perim->parent_residue_index = ri
            perim->parent_residue_number = residues_[ri].sequence_number
            # set fused_partner_index for 5/6 rings to each other (NOT the perimeter)
            five_ring->fused_partner_index = six_ring_idx
            six_ring->fused_partner_index = five_ring_idx
            rings_.push_back(perim)
```

`MapRingSystemToRingTypeIndex(RingSystemKind, Residue)`:

```
Benzene_Phe                     -> PheBenzene
Benzene_Tyr                     -> TyrPhenol
Indole_Trp_5                    -> TrpPyrrole
Indole_Trp_6                    -> TrpBenzene
Imidazole_His + variant 0 (HID) -> HidImidazole
Imidazole_His + variant 1 (HIE) -> HieImidazole
Imidazole_His + variant 2 (HIP) -> HisImidazole
Imidazole_His + variant -1      -> HisImidazole  (ambiguous default)
Pyrrolidine_Pro                 -> [no current RingTypeIndex; see §7.4]
```

The HIS branch reads `Residue::protonation_variant_index`,
which by the time `FinalizeConstruction` runs has been resolved
either by the GROMACS readback block (trajectory loads) or by
`Protein::ResolveProtonationStates` (PDB / ORCA loads). The
`-1` (ambiguous) case maps to `HisImidazole` exactly as the
current code does.

`OrderedAtomIndices(atom_list, RingSystemKind, Residue)`:

The atoms in the group are unordered (insertion order = atom
construction order, which for AMBER is roughly canonical but not
guaranteed). Define a canonical traversal using
`RingPositionLabel`:

- 6-ring (Benzene_Phe / Benzene_Tyr / Indole_Trp_6):
  Ipso → Ortho1 → Meta1 → Para → Meta2 → Ortho2
- 5-ring imidazole (Imidazole_His):
  Ipso → PyrroleBeta(CD2) → Heteroatom_{NH or NoH}(NE2)
  → PyrroleAlpha(CE1) → Heteroatom_{NH or NoH}(ND1)
  (or: walk Ipso → Heteroatom_X → PyrroleAlpha → Heteroatom_Y
   → PyrroleBeta — must match the current order in
   AminoAcidType.cpp:92)
- 5-ring pyrrole (Indole_Trp_5):
  Ipso → PyrroleBeta(CD1) → Heteroatom_NH(NE1) → BridgeFusion(CE2)
  → BridgeFusion(CD2)
- Saturated 5-ring (Pyrrolidine_Pro):
  Cα → Cβ → Cγ → Cδ → N (or as defined; Pro is the only saturated
  case so the convention is local)

The ordering must match the existing `AminoAcidType.cpp` ring
atom-list ordering bit-identically for Phase 1 to be a no-op
on calculator output. A test that compares the new output to a
golden snapshot of the old `Ring::atom_indices` order on
representative fixtures (1UBQ, 1Z9B, 1P9J) catches any
inadvertent reordering. If a difference is unavoidable due to
substrate-side ordering choices, it goes through `golden/
blessed/bless_policy.toml`-style review, with the SVD-normal
sign convention re-checked against a known case.

### 5.3 `AminoAcidType::rings[]`

Phase 1 deletes `AminoAcidType::AminoAcidRing::atom_names`. It
*may* keep the struct around with just `type_index` if any code
benefits from "this residue type has a ring of this kind"
metadata at residue scope (the `is_aromatic` bool answers this
already, as does iterating `RingAt(ri)` after construction).
Recommendation: drop `AminoAcidRing` entirely and remove the
`rings` field from `AminoAcidType`. The `is_aromatic` bool
stays.

### 5.4 What the calculators see

Same `Ring` interface, same `RingTypeIndex` dispatch, same
ordering. `BiotSavartResult`, `HaighMallionResult`,
`RingSusceptibilityResult`, `DispersionResult` continue to read
`ring.Intensity()`, `ring.JBLobeOffset()`, `ring.atom_indices`.
No calculator-side change in Phase 1.

### 5.5 What goes away

- The `pdb_atom_name`-keyed `name_to_idx` map at
  `Protein.cpp:449-452`.
- The string match against `AminoAcidType::rings[].atom_names`
  at `Protein.cpp:454-465`.
- The duplicate HIS-tautomer string fallback at
  `Protein.cpp:485-497`.
- `AminoAcidType::AminoAcidRing` struct + `rings` field.
- The `atom_names` literals in every aromatic residue's
  `AMINO_ACID_TYPES` row.

### 5.6 Validation

The same fixtures the audit calls out for Phase 1 validation:
1UBQ, 1Z9B, 1P9J, plus any PHE/TYR/TRP/HIS-rich case. The
verifying assertion: for every ring in the old detection output
and the new substrate-driven output, the `(parent_residue_index,
type_index, sorted(atom_indices))` triple matches; ordered
`atom_indices` matches per the canonical convention; ring count
matches.

The smoke-test bit-identity contract on classical calculator
output (BiotSavart, HaighMallion, RingSusceptibility) is the
end-to-end check. Drift acceptable on the same terms as the
2026-05-04 bless policy.

### 5.7 Estimated scope

- `src/Protein.cpp`: `DetectAromaticRings` removed (~100 lines);
  `ConstructRingsFromSubstrate` added (~80-120 lines depending
  on how clean the ordering helpers come out).
- `src/AminoAcidType.h`: remove `AminoAcidRing` struct + `rings`
  field (~10 lines).
- `src/AminoAcidType.cpp`: remove `rings: {...}` initialisers
  from PHE, TYR, HIS, TRP rows (~10 lines).
- New helper functions (in Protein.cpp or a new RingConstruction
  unit): `MapRingSystemToRingTypeIndex`, `OrderedAtomIndices`,
  `OrderedPerimeterFromFusedPair` (~80 lines).
- Tests: `tests/test_object_model.cpp` and any ring-construction
  test fixtures need their AminoAcidType-fixture writes
  removed. Identity-driven test fixtures replace them.

Total: ~250-350 lines. Comparable to one of the audit's phases.

---

## 6. Phase 2 considerations

Once Phase 1 lands, the substrate is the authority, the strings
are gone, and `Atom.semantic.ring_position` is populated for
every atom. Phase 2 can use this directly.

### 6.1 Per-position ring-current kernel stratification

Currently `BiotSavartResult` and friends produce a single
shielding contribution per (probe atom, ring) pair, parameterised
by `RingTypeIndex`. The literature observation: ring-current
shifts on the *para* carbon of a Tyr ring are systematically
different from those on the *ortho* carbons (~30 ppm spread on
Tyr ¹³C; the para-OH electron-donation perturbs the ring
current on top of geometric factors). Today this is hidden in
the calibration (Tyr para C is just one of the six ring carbons
in the kernel sum).

With per-atom `ring_position.primary.position`, calculators can
emit separate kernel columns for ipso/ortho/meta/para
contributions. The kernel set grows; the ridge calibration sees
finer-grained features; the per-element-and-position
stratification matches the per-element stratification already
in place (Stage 1 `learn/stage1-mutations`).

This is upstream of any data — it changes the calculator's
output schema. PerFrameExtractionSet would need new NPY columns.
The Python SDK would need new ArraySpecs.

### 6.2 Direct ring-position queries from calculators

```cpp
// Today (chain through Ring.atom_indices, no position info):
for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
    const Ring& ring = protein.RingAt(ri);
    for (size_t ai : ring.atom_indices) {
        // know which atom; do not know ipso/ortho/etc.
    }
}

// Phase 2:
for (size_t ai = 0; ai < protein.AtomCount(); ++ai) {
    const RingPosition& rp = protein.AtomAt(ai).semantic.ring_position;
    if (rp.primary.ring == RingSystemKind::Benzene_Tyr &&
        rp.primary.position == RingPositionLabel::Para) {
        // know this is Tyr-Cζ; can apply different kernel.
    }
}
```

The atom-loop form often replaces the ring-loop form when the
question is per-atom (per-position kernel contribution, per-
position statistic, per-position calibration). The ring-loop
form remains right when the question is per-ring (centre,
normal, total dipole moment).

### 6.3 Pyrrolidine ring (Pro)

Phase 2 is when "should there be a Pro Ring object?" gets
answered by the data. Pro's ring puckering (endo vs exo at Cγ)
is a known degree of freedom that affects backbone angles ψ
and chemical shielding (Cα and Cβ shifts vary by 1-2 ppm
between puckers). A `Ring` object for Pro would surface ring-
puckering measurement as a `RingGeometry` operation — `Ring::
ComputeGeometry` already returns a vertex list and centre; an
additional ring-puckering pseudo-rotation diagnostic could read
from there.

If Phase 2 builds a Pro ring: extend `RingTypeIndex` with
`ProPyrrolidine`, add a `PyrrolidineRing` subclass (or a new
`SaturatedRing` base). Substrate already returns
`Pyrrolidine_Pro` so the construction code from §5.2 picks it
up automatically — only the `MapRingSystemToRingTypeIndex`
table grows.

If Phase 2 doesn't: the substrate carries the chemistry without
needing a `Ring` object to surface it; per-atom queries
(§6.2) work for ring-puckering diagnostics with the substrate
ring_position directly.

### 6.4 Ring-flip diagnostics for symmetric rings

PHE / TYR rings flip on the NMR timescale; ortho1/ortho2 pairs
average to a single observable. The substrate distinguishes
them (Ortho1 / Ortho2); aware analysis can compute symmetric vs
antisymmetric kernel contributions and then take the symmetric
combination as the observable kernel. This is downstream
calibration work; the substrate enables it.

### 6.5 Per-residue ring statistics in TrajectoryResult

A trajectory-scope ring-puckering statistic (Pro endo/exo
fraction over the trajectory; Trp Indole-perimeter dihedral
distribution) attaches naturally as a `*TrajectoryResult` if a
suitable per-frame `Ring` accessor exists. `Ring::Compute
Geometry(positions)` already provides centre/normal; rolling
that into a per-frame angular signal is a few lines on top.
This is queue-worthy for the planned-calculators backlog
(spec/PLANNED_CALCULATORS_2026-04-22.md and successors), not
Phase 2 directly.

---

## 7. Choices that need a user decision

These are decisions the user is the right authority on; they
shape Phase 1 specifics.

### 7.1 Does `Protein::DetectAromaticRings` stay as a function?

**Question:** post-Phase-1, does the function still exist with
the substrate-driven body, or does ring construction dissolve
into the `LookupBy` substrate-population pass that runs over
every atom?

**Position favouring "function stays":** discoverability. A
function named `ConstructRingsFromSubstrate` (or kept as
`DetectAromaticRings` with an updated docstring) called
explicitly from `FinalizeConstruction` makes the construction
order obvious. The line that today reads
`DetectAromaticRings()` becomes `ConstructRingsFromSubstrate()`
and the developer reading FinalizeConstruction knows what's
happening.

**Position favouring "dissolves":** the substrate population
loop already touches every atom. Constructing rings inside the
same loop (group as we go, materialise `Ring` objects at end)
saves a second pass at the cost of coupling.

Recommendation: **stays as a function**. The two passes are
conceptually different — substrate population fills typed
fields per atom; ring construction groups across atoms — and
the two-pass shape matches `FinalizeConstruction`'s existing
`CacheResidueBackboneIndices → ResolveProtonationStates →
DetectAromaticRings → CovalentTopology::Resolve` sequence.

But the user should bless this. The dissolved form is also
defensible.

### 7.2 Removal of `AminoAcidType::AminoAcidRing` entirely?

**Question:** does `AminoAcidType` keep an empty `rings` slot
with just `RingTypeIndex` entries (saying "this residue has
rings of these kinds"), or is the field removed entirely?

The substrate is the authority on which atoms are in which ring
of which kind. The residue-level "what kinds of rings does this
residue have?" question answers from iterating
`Protein::RingAt(ri)` and filtering by `parent_residue_index`,
or from a derived helper. We don't need a static residue-level
table for it.

**Recommendation: remove entirely.** Single authority is the
substrate. PATTERNS.md "single authority for amino acid
chemistry" reads the table as an authority for chemistry the
substrate doesn't have (chi-angle definitions, atom rosters);
the substrate now has ring chemistry, so the table loses that
role.

The user should bless. The conservative alternative (keep the
typed half, drop the strings) is structurally fine but adds an
authority site we'd have to keep in sync.

### 7.3 Atom-ordering convention within `Ring::atom_indices`

**Question:** does Phase 1 preserve the current ordering exactly
(as listed in `AminoAcidType.cpp` for each ring), or adopt the
substrate's typed ordering (Ipso → Ortho1 → Meta1 → Para → Meta2
→ Ortho2 for benzene; analogous for others)?

Calculator behaviour: the SVD ring-normal calculation in
`Ring::ComputeGeometry` orients the normal via a cross product
of the first two edges. Reordering changes which two atoms
those edges go between. For a planar ring the normal direction
flips only if the *cyclic order* reverses (clockwise → counter-
clockwise viewing the ring from one side). The current
PHE ordering `CG → CD1 → CE1 → CZ → CE2 → CD2` is one cyclic
walk; the substrate's Ipso → Ortho1 → Meta1 → Para → Meta2 →
Ortho2 with the right Ortho1/Ortho2 assignment is the same walk.

But the actual mapping AMBER-name → RingPositionLabel for PHE
is: `CG=Ipso, CD1=Ortho1, CE1=Meta1, CZ=Para, CE2=Meta2,
CD2=Ortho2` — assuming the substrate's Ortho1/Ortho2 follows
the chi2-priority convention from Markley 1998 (per
`SemanticEnums.h:540`). If that mapping holds, current order
and substrate order are the same. If Ortho1 vs Ortho2 are
assigned the other way around, current order is `Ipso → Ortho1
→ Meta1 → Para → Meta2 → Ortho2` but the cyclic ring is
traversed the other direction, which flips the SVD normal sign.

**Recommendation:** verify the mapping at fixture level by
running the new construction code, comparing ring-normal
direction against the old code, and adjusting the
`OrderedAtomIndices` traversal definition (Ortho1-first vs
Ortho2-first) so signs match. Lock it in a test; document the
convention in `Ring.h`.

### 7.4 Should there be a Pro `Ring` object?

**Question:** Phase 1 — does the new construction code produce
a `Ring` for Pro's pyrrolidine, or skip it?

The substrate emits `Pyrrolidine_Pro` for Pro N + Cα + Cβ + Cγ
+ Cδ. There is no `RingTypeIndex` for it today; no `Ring`
subclass; no calculator currently consumes a Pro ring. Adding
one in Phase 1 introduces:
- A new `RingTypeIndex::ProPyrrolidine`.
- A new `Ring` subclass `ProPyrrolidineRing` with
  `Aromaticity() = None`, `RingSizeValue() = 5`,
  `Intensity() = 0` (no ring current; saturated ring).
- New `CreateRing(ProPyrrolidine)` factory case.
- Empty-but-present effect: `BiotSavartResult` etc. iterate
  every ring; encountering a ring with `Intensity() = 0`
  produces zero contribution and no behavioural change.

Skipping Pro in Phase 1 is the alternative — keep it in the
substrate (no work to remove it), don't materialise a `Ring`
object. Phase 2 reconsiders if a calculator wants it.

**Recommendation: skip in Phase 1.** Phase-1 goal is calculator-
behaviour-preserving; adding a ring with non-trivial geometry
that no calculator yet uses is yak-shaving. If Phase 2 needs
Pro ring puckering, the substrate-driven construction path
trivially extends to it.

The user should bless. Adding it in Phase 1 is also defensible
(no behavioural change; surfaces the chemistry symmetrically;
Phase 2 starting point is shorter). Decision is style.

### 7.5 Ordering and naming: does the function get renamed?

**Question:** post-Phase-1, is the function still
`DetectAromaticRings` (now misnamed because it includes
non-aromatic Pro if §7.4 lands), `ConstructRings`,
`ConstructRingsFromSubstrate`, or another name?

Recommendation: **`ConstructRingsFromSubstrate`**. Names that
describe the shape of the operation rather than the chemistry
result are clearer in this codebase (audit's hotspot 3
discussion uses similar shape-vocabulary). The "Aromatic"
qualifier is wrong post-Pro; the "Detect" verb is the wrong
verb for "look up in substrate."

User decision is style.

### 7.6 Docstring convention

The current function has a "PDB LOADING BOUNDARY" comment
(lines 420-424) marking it as the place where strings cross
into typed objects. Post-Phase-1 there are no strings here at
all. The docstring should be replaced. Recommended convention:

```cpp
// ConstructRingsFromSubstrate -- typed authority chain.
//
// Per-atom RingPosition substrate (populated by LookupBy at
// FinalizeConstruction from the generated AtomSemanticTable
// rows) is the authority for which atoms belong to which ring
// in which residue. This function groups atoms by
// (residue, RingSystemKind), maps to RingTypeIndex via
// residue.protonation_variant_index for HIS, and constructs
// the typed Ring objects calculators consume.
//
// No string surface. AminoAcidType::rings[] (formerly the
// const char* list of atom names) was removed in <commit>.
```

User decision: tone, length.

### 7.7 What about ChiAngleDef?

The audit's Hotspot 2 mentions chi angles use the same
`const char*[4]` shape (`AminoAcidType::ChiAngleDef::atoms`);
the recommendation there mirrors this one (encode as
identity). Phase 1 of *this* refactor should leave chi angles
alone — they are a separate string surface, called out in a
separate hotspot, with separate consumers. The user has already
seen many in-flight scopes balloon out from "while we're
here..." additions; this question intentionally stops at rings.

Question to confirm: do chi-angle changes go in their own
slice, or are they coupled to ring changes? Recommendation:
their own slice (separate audit hotspot, separate calculator
consumers).

---

## 8. Risks and unknowns

### 8.1 Substrate population not yet wired

The single largest risk: substrate population at
`Protein::FinalizeConstruction` is the prerequisite this work
depends on, and it is not yet sequenced (audit closing note 2).
If substrate population slips, ring detection refactor blocks.

Mitigation paths:
- Path B (typed atom-list in AminoAcidType) is implementable
  *before* substrate population, against `Atom.identity` only.
  It's strictly worse than Path A but unblocks the ring
  detection refactor on its own schedule. If we accept Path B
  as interim, the destination shape (Path A) becomes a Phase
  1.5 follow-up after substrate population lands.
- Recommendation: do Path A. Sequence after substrate
  population. The rings refactor is small enough (~250 lines)
  that it doesn't need to ride alongside substrate population's
  own work.

### 8.2 Atom-ordering verification

Per §7.3, Ortho1-vs-Ortho2 assignment in the substrate must
match the current AminoAcidType ordering, *or* the
OrderedAtomIndices function inverts the direction. If neither,
SVD ring-normal sign flips silently and downstream sign
conventions in BiotSavart / HaighMallion break.

Mitigation: test the ring-normal direction on a fixture
explicitly, with a known reference vector (e.g., for 1UBQ's
PHE45 ring, we know which side of the plane the protein
backbone is on). Lock in the result.

### 8.3 TRP perimeter synthesis

The synthetic 9-atom perimeter ring is currently constructed by
listing atoms in a defined order in `AminoAcidType.cpp:188`.
The substrate doesn't have a `RingSystemKind::Indole_Trp_
Perimeter`. The new code synthesises the perimeter from the
union of the 5-ring and 6-ring atom lists; the canonical order
of that union (CG → CD1 → NE1 → CE2 → CZ2 → CH2 → CZ3 → CE3 →
CD2 in current code) must be reproducible from the typed atom
identities.

Mitigation: define `OrderedPerimeterFromFusedPair(five_ring,
six_ring)` to produce exactly this sequence. Verifiable by
comparing the resulting `atom_indices` to the old code's
output on a 1UBQ-class fixture with at least one Trp.

### 8.4 Non-standard residues

The substrate generator runs against the standard 20. Modified
amino acids, ligands, terminal caps don't have substrate rows
and don't carry typed `ring_position`. If a future load path
delivers a phosphorylated tyrosine or a modified Trp with a
substrate-undefined identity, ring construction silently
produces nothing.

Mitigation:
- For now, our load paths handle only standard 20 plus AMBER
  variants (HID/HIE/HIP/CYX/etc.); this is not a gap in
  practice.
- Phase 2 substrate generator extension to non-standard
  residues (per topology-residue-reference plans) addresses
  this naturally.
- A diagnostic: if `Protein::RingCount()` for a residue with
  `is_aromatic == true` is zero post-construction, log it.
  Catches generator gaps; does not paper over them.

### 8.5 The "Aromaticity()" virtual on Ring

`Ring::Aromaticity()` returns `RingAromaticity::Full / Reduced
/ Weak / None`. This is a per-RingTypeIndex constant in the
subclass, not derived from substrate. Phase 1 leaves it alone;
the values map cleanly. But the substrate carries
`RingMembership::aromatic` (bool), `RingMembership::planar`
(bool), `RingMembership::n_heteroatoms` (uint8_t). Phase 2
might surface these per-atom and the per-ring `Aromaticity()`
becomes derivable, removing one more piece of constant data
from the subclass override.

Not a Phase 1 concern. Note for Phase 2.

### 8.6 The `accumulated` field on Ring

`Ring::accumulated` is a mutable struct (`RingAccumulated` —
total B at center, intensity used, etc.) populated during the
extraction passes. Substrate-driven construction does not
affect this — the field is set later, by calculators, on the
`Ring` object regardless of how it was constructed. Phase 1
preserves this contract.

### 8.7 Tests as a pivot signal

`tests/test_object_model.cpp` currently builds tiny proteins
with `n->pdb_atom_name = "N";` etc. The audit has these as
Category F (test-code reflection of the surface). Phase 1 of
this refactor doesn't strictly need them migrated — the rings
in those tests are constructed by the production code path,
which (post-pivot) reads substrate, which (in tests) is
populated by the same `LookupBy` running over the test atoms.
The test atoms need to have `Atom.identity` set (so LookupBy
can find them) and that's it; the strings can stay until audit
Phase 7.

But: any test that previously asserted ring detection by
checking `ring.atom_indices` against a hand-written list of
atom indices will continue to work iff the new ordering
matches the old. Per §7.3.

### 8.8 Coupling to audit phasing

This refactor sits between audit Phase 1 (add Atom.identity +
substrate population) and audit Phase 2 (migrate other internal
Protein.cpp builders: CacheResidueBackboneIndices,
ResolveProtonationStates). It IS the rings part of audit Phase
2, broken out for separate analysis.

If audit Phase 2 lands as a single commit doing all three
(backbone cache, protonation, rings) — the audit's recommended
shape — then this analysis applies but the slicing is one
commit not three. If audit Phase 2 lands per-builder, this is
its own commit. User decides.
