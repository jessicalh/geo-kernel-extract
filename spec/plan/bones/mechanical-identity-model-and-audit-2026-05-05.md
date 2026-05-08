# AtomMechanicalIdentity — model definition + audit of `pdb_atom_name` (2026-05-05)

## Reviewer

Audit-and-define agent commissioned by the user to (1) state precisely
what `AtomMechanicalIdentity` IS in the model, and (2) walk every
existing `pdb_atom_name` call site against that definition. The user
has been burned by prior sessions agreeing in principle and then
wiring strings into calculators anyway. This audit is written with
the explicit instruction to call uses as they are, and to flag
anything borderline rather than absolve it.

The audit was performed by reading: `OBJECT_MODEL.md`, `PATTERNS.md`,
`spec/CONSTITUTION.md`, `spec/plan/topology-encoding-dependencies-2026-05-05.md`
(esp. §H), `src/SemanticEnums.h`, `src/generated/LegacyAmberSemanticTables.h`,
`src/Atom.h`, and the 24 `.cpp/.h` files containing `pdb_atom_name`.

---

## Part 1 — Definition

### 1.1 What is `AtomMechanicalIdentity` in the model?

**Shape (`src/SemanticEnums.h:848-863`).** A 5-tuple of typed enums /
typed structs:

```cpp
struct AtomMechanicalIdentity {
    Element             element;       // H, C, N, O, S
    Locant              locant;        // None, Alpha, Beta, ..., Eta
    BranchAddress       branch;        // {outer ∈ 0..2, inner ∈ 0..2}
    DiastereotopicIndex di_index;      // None, Position2, Position3
    BackboneRole        backbone_role; // None, N, CA, C, O, H, HA
};
```

All five fields are typed, none is a string, and the tuple is
constructible from typed parser output (`ParseAtomName` /
`ComputeAtomMechanicalIdentity` lifted in
`LegacyAmberSemanticTables`). Equality is defined as field-wise
equality.

**Role.** `AtomMechanicalIdentity` is a **typed lookup key into the
chemistry-substrate tables** built at code-generation time from
CCD + RDKit + Markley conventions. Its sole runtime use is
`LookupBy(residue, variant_idx, identity)` /
`LookupCap(state, identity)` returning a pointer to the
`AtomSemanticTable` row for that atom. The substrate row is the
answer surface; the identity is the key.

Per §H.1 / §H.7 of `topology-encoding-dependencies-2026-05-05.md`
(emphasis added): "NOT atom_local_idx (index spaces don't align
across CCD, AmberAminoAcidVariantTable, and other producers) and
**NOT atom_id strings (the string wall is sacred)**." The five
fields are typed precisely because the requirement was a typed
substitute for the indexing role atom names had previously played.

Concretely: identity is **(c) a graph-navigation primitive +
(d) a substrate-side identity** in the sense of §3.3.5 of the
spec — it tells you "this atom IS the prochiral-2 beta-hydrogen
on this residue" using only typed enum values. It is NOT a typed
substitute for `pdb_atom_name`; the name carries a per-atom string
the substrate row does not. Identity is the locator, not the
chemistry; the chemistry comes from the row. (§1.4 below settles
this.)

This is consistent with the project conventions:
- **`PATTERNS.md` "Objects answer questions about themselves"**:
  identity asks the substrate for the answer; it does not embed
  the answer.
- **`PATTERNS.md` "Translations between naming systems happen in
  one place"**: identity is the post-translation typed object.
  The string lived at the load boundary; identity replaces it
  going forward.
- **`OBJECT_MODEL.md` substrate-vs-conformation discipline**:
  identity is an INVARIANT atom-level property. It does not change
  between conformations (Element, Locant, etc. are all geometry-
  free). It belongs on the chemistry-substrate side.

### 1.2 What questions does identity legitimately answer?

Identity answers **graph-navigation and substrate-lookup
questions**:

- **"Look up this atom's substrate row."**
  `LookupBy(residue.type, residue.protonation_variant_index,
  identity) → AtomSemanticTable*`. This is the *only* primary
  consumer per §H.5. Justified.

- **"Find the atom in residue R that is at locant gamma, branch 1,
  pro-R."** A typed iteration over `res.atom_indices` testing
  identity equality. Useful for: "does this residue have an HG2 vs.
  HG?"-shaped questions, where the answer is structural. A typed
  query — `Residue::AtomWithIdentity(identity) → atom_index` —
  is the right surface to expose; calculators do NOT iterate
  atom names directly. Justified.

- **"Are these two atoms chemically equivalent (do they belong to
  the same equivalence class within the residue)?"** Methyl Hs
  share an identity by design (the docstring on
  `AtomMechanicalIdentity` makes this explicit); arginine
  guanidinium Hs do not. Identity equality is the right test for
  NMR-equivalence-by-fast-exchange. Justified — but note this is
  identity equality, NOT atom equality (§1.6).

- **"Find the next residue's `BackboneRole::AlphaCarbon`."** When
  a calculator needs cross-residue navigation, the right form is
  `next_res.AtomWithBackboneRole(BackboneRole::AlphaCarbon)`,
  which is `next_res.CA` already if the cache exists, or a typed
  lookup using a partial identity if it doesn't. (§1.5 expands.)

### 1.3 What questions does identity NOT answer?

Identity does NOT answer chemistry questions. Those go through
the substrate row (`AtomSemanticTable`) populated via `LookupBy`
at construction time:

- **"Is this atom polar?"** → `atom.semantic().polar_h !=
  PolarHKind::NotPolar`, **not** "is the locant Alpha and the
  element H" or any name-based proxy.

- **"Is this atom in a peptide-amide group?"** →
  `atom.semantic().planar_group == PlanarGroupKind::PeptideAmide`,
  not "is the residue not Pro and the role
  `BackboneRole::AmideHydrogen`."

- **"Is this atom aromatic?"** → `atom.semantic().aromatic`, not
  "is the residue PHE and the locant Delta-or-Epsilon-or-Zeta."

- **"What's this atom's formal charge?"** →
  `atom.semantic().formal_charge` (Lewis-localised per the
  encoding-dependencies §C.1), not any name-based pattern matching.

- **"What pseudoatom group does this atom belong to?"** →
  `atom.semantic().pseudoatom`, not a methyl-H-name regex.

The principle: **identity locates; the row answers.** If a
calculator branches on identity fields directly to make a
chemistry decision, it has cut out the substrate and is doing
synthesis it shouldn't be. The exception is when the decision IS
graph-navigational ("which atom is HB2 vs HB3?"), and even then
the right surface is a typed accessor that hides the lookup.

This matches `PATTERNS.md`'s anti-string-dispatch discipline.
The pre-substrate version of this rule was "no `atom.pdb_atom_name
== "HB2"` in calculators." The post-substrate version is "no
`atom.identity().di_index == DiastereotopicIndex::Position2 &&
atom.identity().locant == Locant::Beta` in calculators either,
unless the calculator is doing graph navigation specifically." A
typed disguise of the same string-dispatch hack is still
string-dispatch.

### 1.4 Is identity a calculator surface, or only substrate-side?

**Position: identity is substrate-side machinery; calculators
consume the typed semantic fields, not identity itself, except for
narrow graph-navigation use.**

The argument FOR identity-as-calculator-surface:
- Calculators need to find specific atoms (e.g., "the HB2 of this
  residue, for a methylene-asymmetry calculation"). If they can't
  use identity, they need name-keyed maps, which is the very
  string-dispatch the pivot is removing.
- Identity is typed; using it doesn't violate the string wall.

The argument AGAINST:
- Once identity is a calculator surface, calculators will branch
  on its fields (`if (id.locant == Locant::Beta && id.di_index ==
  Position2)`) and the string-dispatch antipattern returns under
  a typed disguise. The user has stated this is the failure mode
  to prevent.
- The five identity fields are *deliberately not* the right
  vocabulary for chemistry decisions. `Locant::Beta` does not
  mean "polar" or "aromatic" or "hydrogen-bond donor"; the
  substrate fields do. A calculator branching on `Locant` is
  re-deriving chemistry the substrate already encoded, badly.
- Equivalent chemically-distinct atoms collide in identity by
  design (methyl Hs); calculators using identity directly to
  pick "atom HB2" will sometimes get "atom HB3" instead, depending
  on which the row holds.

**Resolution.** Identity is exposed to calculator code through ONE
surface: typed accessors on `Residue` (or on `Atom` itself for
identity equality testing). The accessors hide `LookupBy` and
return `atom_index` or `nullopt`. The calculator never assembles
an `AtomMechanicalIdentity` literal in its body. This keeps the
mechanism for graph navigation while preventing field-by-field
chemistry dispatch.

The right shape is something like:

```cpp
// On Residue (or accessible from it):
std::optional<size_t> AtomWithRole(BackboneRole role) const;
std::optional<size_t> AtomWithLocant(Locant l, BranchAddress b) const;
std::optional<size_t> AtomWithIdentity(const AtomMechanicalIdentity& id) const;

// On Atom or via Protein:
const AtomSemanticTable& Semantic(size_t ai) const;
const AtomMechanicalIdentity& Identity(size_t ai) const;
bool SameRole(size_t a, size_t b) const;  // identity-equality test
```

Calculators consume `Semantic(ai).polar_h` etc. They use
`AtomWithRole` for graph navigation; they do not construct
identities themselves. Cross-residue navigation goes through
typed accessors on the parent objects, not through raw identity
construction in calculator code.

### 1.5 Cross-residue graph navigation

A calculator may legitimately need "the Cα of the next residue."
Options:

- **Typed cache field** (current shape for backbone): `res.CA`
  is a `size_t`. For backbone atoms this is the right answer
  because the cache exists.
- **`BackboneRole::AlphaCarbon` lookup**: if the cache field is
  empty (e.g., a residue type doesn't have a CA — which doesn't
  arise in this project, but would for non-standard residues),
  this reduces to the same answer.
- **Partial identity**: `res.AtomWithBackboneRole(AlphaCarbon)` —
  same answer, slightly different shape.

For **side-chain cross-residue navigation** (e.g., "Asn HD21 of
the i+1 residue"), the cache is not populated; a typed accessor
is needed. The right surface is:

```cpp
std::optional<size_t>
Residue::AtomWith(Locant l, BranchAddress b,
                  DiastereotopicIndex d = None,
                  Element e = Unknown) const;
```

Or, for the rare case where every field matters (NMR-equivalence
testing), `AtomWithIdentity(const AtomMechanicalIdentity&)`. The
calculator does not build a string and does not iterate
`atom_indices` checking `pdb_atom_name`. The accessor hides the
linear scan or any future indexing.

For backbone: `res.CA`, `res.N`, etc. continue to work; they ARE
the cache. The post-pivot question is whether those indices are
populated by string equality (`if (name == "CA")`, current code)
or by identity equality (`if (Identity(ai).backbone_role ==
AlphaCarbon)`). Both produce the same `size_t` field; only the
former requires the string.

### 1.6 Equality semantics

Two distinct issues:

**Identity equality vs. atom equality.**
- `AtomMechanicalIdentity` `==` is field-wise typed equality.
  Two atoms with `id_a == id_b` may or may not be the same atom.
  Methyl Hs collide deliberately (per the docstring). Across
  residues, identity equality is meaningless (identities are
  local within a residue table).
- "Same atom" is `atom_index` equality (a `size_t`) or `Atom*`
  pointer equality (since `Protein` is non-movable per the
  back-pointer-safety lesson, pointers are stable).
- "Same role within different residues" is identity equality
  PLUS residue-equivalence (same `AminoAcid` type and
  `protonation_variant_index`).

**When to use which.**

| Question | Right test |
|---|---|
| "Is this the same physical atom?" | `ai == aj` (size_t) |
| "Are these atoms NMR-equivalent within a residue?" | identity equality |
| "Is this atom the backbone Cα?" | `Identity(ai).backbone_role == AlphaCarbon` (typed accessor) OR `ai == res.CA` (cache) |
| "Are these two residues' Cα atoms 'the same role'?" | identity equality on backbone_role field (or just check `BackboneRole::AlphaCarbon` on each separately) |

Calculators need both. The atom-index test is the dominant one
(every spatial calculation, every neighbour list); the identity
test is rare and goes through typed accessors that wrap it.

### 1.7 Post-pivot architecture

**Where strings live (the only legitimate sites):**

1. **PDB / AMBER / TPR loaders** (PdbFileReader, FullSystemReader,
   OrcaRunLoader). Strings enter from external file bytes,
   immediately produce typed identity + populate substrate rows
   via `LookupBy`. The string can be preserved temporarily for
   diagnostics within the loader, but does not flow out onto
   the typed `Atom`.

2. **File-format emitters** (AmberLeapInput, DSSP-PDB writer in
   DsspResult, KamlProtonator-PDB writer, PropkaProtonator-PDB
   writer). PDB ATOM records require column-13-16 atom names —
   this is wire format. Per the **Crystal Projection Rule**
   (`amber-implementation-plan-2026-04-29.md`): names are PURE
   FUNCTIONS on the typed substrate, computed at the emit site
   from typed identity, never cached. The emitter's local helper
   `AmberAtomNameFor(identity, residue.type, terminal_state)` is
   the right shape; it lives next to the writer.

3. **External-table key projection** (ChargeSource
   ParamFileChargeSource, AmberChargeResolver verdict, AmberPrepared
   PRMTOP atom-name match). The flat ff14SB dat file is keyed by
   atom-name strings; the PRMTOP format carries atom-name strings.
   These are wire-boundary lookups; the typed-to-string projection
   happens at the lookup site. Internal storage stays typed.

**What lives on `Atom` / `Protein` instead:**

- `Atom.identity` of type `AtomMechanicalIdentity` (typed,
  populated at construction by `ComputeAtomMechanicalIdentity`).
- `Atom.semantic` of type `AtomSemanticTable` (or a pointer to
  the canonical row for this atom's residue × variant × identity,
  populated at `Protein::FinalizeConstruction` via `LookupBy` /
  `LookupCap` / `ApplyCapDelta`).
- Backbone caches on `Residue` (already exist; populated from
  `Identity(ai).backbone_role` instead of `pdb_atom_name == "CA"`).
- Typed accessors on `Residue` for cross-residue navigation when
  the backbone cache isn't enough.

**What comes off `Atom`:**

- `pdb_atom_name`. Quarantined to load boundary, then deleted
  from the type. Calculator code never sees it. Tests asserting
  the surface go away with the surface.

**Phasing for the pivot** (sketched, more in §2.4):
1. Add `Atom.identity` populated from existing `pdb_atom_name`
   via `ComputeAtomMechanicalIdentity` at construction time
   (parallel surface).
2. Add `Atom.semantic` populated via `LookupBy` at
   `FinalizeConstruction` (parallel surface).
3. Migrate internal `Protein.cpp` index-builders
   (`CacheResidueBackboneIndices`, `DetectAromaticRings`,
   `ResolveProtonationStates`) from string match to identity /
   semantic match. These are the single hardest sites because
   the chemistry decisions are fused with the string keying.
4. Migrate calculators: replace identity-proxy hacks with
   semantic-row reads (Category C below). This is straightforward
   once the substrate is populated.
5. Migrate file-format emitters to compute names from typed
   identity at emit time (Crystal Projection Rule).
6. Migrate external-table lookups (charges) to project at the
   lookup site.
7. **Delete `Atom.pdb_atom_name`**. Loaders use it only on the
   stack inside their own scope; once ingested, the typed
   identity is the only surface.

The user has stated firmly: deletion is the goal, not quarantine.
The phasing makes the deletion deferred but inevitable; the
field comes off the type at the end of step 7.

---

## Part 2 — Audit

### 2.1 Inventory

`grep -rn pdb_atom_name --include='*.cpp' --include='*.h' src/
tests/ tools/` returns **49 hits across 22 files**, including the
declaration on `Atom.h` and the documentation hit in
`PdbFileReader.h`. The original brief estimated ~67 hits / ~24
files; the actual numbers are slightly lower because some of the
hits projected from EVIL_STRING_AUDIT (2026-04-28) have been
collapsed by intermediate refactors (e.g., `MutableLegacyAmber`
removal, `AmberFFData` removal). The four `tests/bones/` hits
are quarantined CHARMM-path code per
`project_charmm_retired_amber_only_2026-05-02`.

**Files (excluding `tests/bones/` quarantine):**

Source (15 files):
- `src/Atom.h` — type declaration
- `src/PdbFileReader.h` — header doc only
- `src/PdbFileReader.cpp` — load-write site
- `src/FullSystemReader.cpp` — load-write site (TPR / AMBER trajectory)
- `src/OrcaRunLoader.cpp` — load-write site (XYZ)
- `src/TrajectoryProtein.cpp` — H5 emit-read site
- `src/Protein.cpp` — internal index/ring/protonation builders (4 sites)
- `src/AmberLeapInput.cpp` — file-format emit + disulfide detect (3 sites)
- `src/AmberChargeResolver.cpp` — verdict (2 sites)
- `src/ChargeSource.cpp` — ff14SB lookup (2 sites)
- `src/AmberPreparedChargeSource.cpp` — PRMTOP match (4 sites)
- `src/PropkaProtonator.cpp` — PDB writer (3 sites)
- `src/KamlProtonator.cpp` — PDB writer (3 sites)
- `src/DsspResult.cpp` — PDB writer (3 sites)
- `src/AIMNet2Result.cpp` — error-message use (1 site)
- `src/ApbsFieldResult.cpp` — error-message use (1 site)

Tests (5 files):
- `tests/test_object_model.cpp` (5 sites)
- `tests/test_amber_prepared_charge_source.cpp` (2 sites)
- `tests/test_amber_charge_resolver.cpp` (2 sites)
- `tests/test_amber_leap_input.cpp` (3 sites)
- `tests/test_foundation_results.cpp` (1 site)
- `tests/test_traversal_dump.cpp` (2 sites: comment + display printf)

Quarantined (do not pivot, but recorded):
- `tests/bones/src/GromacsEnsembleLoader.cpp` (4 sites; CHARMM
  path, retired 2026-05-02 per memory entry).

**Category overview** (anticipating §2.2):

| Category | Count |
|---|---|
| A. Load-boundary write | 4 sites |
| B. Boundary-out string-key projection | 9 sites |
| C. Identity-proxy hack | **8 sites** (the rot) |
| D. Diagnostic / log / error message | 6 sites |
| E. Internal model use that should not have been allowed | **5 sites** (Protein.cpp) |
| F. Test-code reflection of the surface | 13 sites |

The 8 Category-C sites + 5 Category-E sites = **13 sites that
should not exist post-pivot**. The 9 Category-B sites need
projection-at-the-call-site refactoring. Categories A, D, F
mostly survive in form (their inputs change).

### 2.2 Per-file classification

#### `src/Atom.h:23` — type declaration

| | |
|---|---|
| Role | Header declaring `pdb_atom_name` field on `Atom` |
| Use | `std::string pdb_atom_name;` |
| Classification | **The surface itself.** Not a use site, but the source of every other site. |
| Risk | HIGH (the existence of this field IS the risk) |
| Recommendation | **Delete after the pivot.** Add `AtomMechanicalIdentity identity;` and (eventually) `AtomSemanticTable* semantic;` (or stored-by-value `AtomSemanticTable semantic;`). Until deletion is sequenced, marking the field deprecated and `[[deprecated]]`-ing direct reads is a useful in-flight signal. |

#### `src/PdbFileReader.h:12` — header documentation

| | |
|---|---|
| Use | Comment explaining "Atom-name strings from cif++ are assigned raw to `Atom::pdb_atom_name` — no NamingRegistry translation on this path." |
| Classification | D (documentation; not load-bearing) |
| Risk | LOW |
| Recommendation | Update post-pivot: "atom names from cif++ are translated to `AtomMechanicalIdentity` via `ComputeAtomMechanicalIdentity` and stored as the typed identity field; the load-time string is local-scoped." |

#### `src/PdbFileReader.cpp:133` — load-boundary write

| | |
|---|---|
| Use | `new_atom->pdb_atom_name = atom_name;` where `atom_name = atom.get_label_atom_id();` from cifpp |
| Classification | **A (load-boundary write).** Strings legitimately enter at this exact site. |
| Risk | LOW (legitimate boundary) |
| Recommendation | Post-pivot, this becomes: keep `atom_name` local; produce `new_atom->identity = ComputeAtomMechanicalIdentity(atom_name, residue_type, ...)` and discard the string. Document that the string never escapes the function. |

#### `src/FullSystemReader.cpp:876,884` — TPR/AMBER load + backbone cache

| | |
|---|---|
| Use line 876 | Writes `atom->pdb_atom_name = registry.TranslateAtomName(charmm_atom_name, canonical_residue, ToolContext::Charmm, ToolContext::Standard);` |
| Use line 884-896 | Reads the just-written name to populate `res_ref.N`, `res_ref.CA`, etc. via `name == "N"` etc. — string-keyed backbone cache. |
| Classification | Line 876: **A** (load-boundary write, with explicit NamingRegistry translation — this is the pattern PATTERNS.md endorses). Line 884-896: **E** (internal model use; backbone-cache string keying that should be identity keying). |
| Risk | MEDIUM (the cache shape is right; the keying isn't typed) |
| Recommendation | Line 876: Move the typed-identity computation here too (`new_atom->identity = ComputeAtomMechanicalIdentity(...)`); name string is local. Line 884-896: Replace `if (name == "N")` with `if (new_atom->identity.backbone_role == BackboneRole::Nitrogen)` etc. The backbone cache stays; its keying becomes typed. |

#### `src/OrcaRunLoader.cpp:206` — load-boundary write

| | |
|---|---|
| Use | `atom->pdb_atom_name = atom_names[ai];` where `atom_names` came from ORCA output parser |
| Classification | A (load-boundary write) |
| Risk | LOW |
| Recommendation | Same as PdbFileReader: produce `identity` from the string; discard the string. |

#### `src/TrajectoryProtein.cpp:257,269,274` — H5 emit-read

| | |
|---|---|
| Use 257 (comment) | Documentation about emitting `pdb_atom_name` |
| Use 269 | `atom_names[i] = a.pdb_atom_name;` builds an emit vector |
| Use 274 | `atoms.createDataSet("pdb_atom_name", atom_names);` writes an H5 dataset |
| Classification | **B (boundary-out string-key projection).** The H5 file is a wire format; downstream Python SDK consumers parse `pdb_atom_name` from it. |
| Risk | MEDIUM (emit lock-in: SDK contract baked into a string column) |
| Recommendation | Post-pivot: the OBJECT_MODEL §H5 schema notes "/atoms/legacy_amber/ with the typed `LegacyAmberTopology` semantic fields … plus naming projections (AMBER native, IUPAC, BMRB) emitted as separate H5 columns derived from the typed enums via `LegacyAmberTopology`'s projection surface; no label strings are stored authoritatively." This is the right answer. The emit path computes `AmberNameFor(identity, residue.type, terminal_state)` at the H5-emit site (Crystal Projection Rule); the names are derived columns, not the source. The Python SDK reads them as before. The substrate columns are the new authoritative surface. Sequencing: emit BOTH during transition, deprecate name column when SDK readers migrate. |

#### `src/Protein.cpp:329` — `ResolveProtonationStates`

| | |
|---|---|
| Use | Builds `name_to_idx` map from `pdb_atom_name`, then looks for "HD1", "HE2", "HD2", "HZ1", "HZ2", "HZ3", "HG", "HH", "SG" to detect HID/HIE/HIP/ASH/GLH/CYX/LYN/TYM. |
| Classification | **E (internal model use that should not have been allowed).** This is chemistry inference from atom names — exactly the antipattern PATTERNS.md prohibits. The signal "this residue's protonation variant" is being read off PDB atom names. |
| Risk | HIGH. The list of names is hardcoded; if a load path produces a name that the encoding doesn't enumerate, the variant is silently misclassified. CHARMM's HSD/HSE → HID/HIE translation in NamingRegistry happens earlier; if the translation misses, this falls through to "ambiguous." |
| Recommendation | Post-pivot, this function should branch on **substrate semantics** populated from the source: GROMACS-readback (already authoritative for trajectory loads, see `GromacsToAmberReadbackBlock`) sets `protonation_variant_index` directly; PDB-load detection becomes "atom-presence" via typed identity (`res.AtomWithIdentity({Element::H, Locant::Delta, {1, 0}, ...})` returns the HD1 atom or nullopt). The branch is then on identity-presence, not on name match. Even better: this whole function moves into the loader (PdbFileReader specifically) where the load context is fresh, and the resulting protonation state is set on `Residue` once during construction. The `Protein.cpp` call site becomes "trust the resolved state." |

#### `src/Protein.cpp:451` — `DetectAromaticRings` ring-atom matching

| | |
|---|---|
| Use | Builds `name_to_idx` from `pdb_atom_name`, then matches ring definitions from `aatype.rings[].atom_names` (which are `const char*` literals like `"CG"`, `"CD1"`, ...). |
| Classification | **E (internal model use).** Matching `AminoAcidType` ring-definition strings against `pdb_atom_name`. The PATTERNS.md "PDB LOADING BOUNDARY" comment is present (line 422-424); the function is *labelled* as the boundary, but the string boundary is being repeated every call to `FinalizeConstruction`. |
| Risk | MEDIUM. The string keys are stable; the names in `AminoAcidType` are project-controlled (canonical AMBER names). But the entire RING structure of every aromatic residue depends on string equality with these literals — a minor regression in any loader's name normalisation breaks ring detection silently. |
| Recommendation | The right shape is to encode ring atom membership as `AtomMechanicalIdentity` (or a partial identity, specifically `(Element, Locant, BranchAddress)`) in `AminoAcidType`, and resolve via `res.AtomWithIdentity(...)` instead of a name table. This is a non-trivial refactor of `AminoAcidType` itself; sequenced after the per-atom identity surface lands. **Do not call this acceptable** because of the boundary comment; the comment marks intent, not status. The string match runs every load. |

#### `src/Protein.cpp:545,563` — `CacheResidueBackboneIndices`

| | |
|---|---|
| Use 545 | `if (name == "N") res.N = ai; else if (name == "CA") res.CA = ai;` ... seven backbone names string-matched. |
| Use 563 | Chi-angle atom matching — `if (atoms_[ai]->pdb_atom_name == def.atoms[j])` against the `ChiAngleDef::atoms` `const char*[4]` from `AminoAcidType`. |
| Classification | **E (internal model use).** Backbone-cache populator + chi-angle index resolver. Both run at every protein construction. |
| Risk | HIGH for line 545 — every calculator that uses `res.CA` depends on this matching working. Silent failure means `res.CA = NONE` and downstream calculators see `Residue::NONE`. The "`HN`" alias for "`H`" and "`HA2`" alias for "`HA`" at line 550-551 is exactly the string-juggling PATTERNS.md says should be one-shot at NamingRegistry. |
| Recommendation | Post-pivot: this function becomes typed: `if (Identity(ai).backbone_role == BackboneRole::Nitrogen) res.N = ai;` etc. The aliases (`H` vs `HN`, `HA` vs `HA2`) collapse — both map to the same `BackboneRole`. For chi angles: encode `ChiAngleDef` with `AtomMechanicalIdentity` instead of `const char*`. |

#### `src/AmberLeapInput.cpp:285` — PDB ATOM record emit

| | |
|---|---|
| Use | `WriteAtomRecord(pdb_out, serial, atom.pdb_atom_name, ambr_name, ...);` writing AMBER-conforming PDB for tleap input. |
| Classification | **B (boundary-out string-key projection).** PDB column 13-16 is wire format; tleap will reject input it doesn't recognise. |
| Risk | MEDIUM (the string is the wire format, but going through `pdb_atom_name` instead of computing it at emit time means the emitter trusts the loader's string preservation — a fragile contract). |
| Recommendation | Post-pivot, emit via `AmberAtomNameFor(identity, residue.type, terminal_state)`. The function lives next to `WriteAtomRecord`; it knows AMBER's atom-name conventions (HG vs HG1 for SER, etc.) and computes from typed identity. The Crystal Projection Rule from `amber-implementation-plan-2026-04-29.md` mandates this exact shape. |

#### `src/AmberLeapInput.cpp:320,328` — `DetectDisulfides` SG matching

| | |
|---|---|
| Use 320 (comment) | `// and atom.pdb_atom_name == "SG"` — describes the function |
| Use 328 | `if (atom.pdb_atom_name == "SG" && atom.element == Element::S)` — finds the Sγ atom of each CYS residue |
| Classification | **C (identity-proxy hack)** in non-trivial form. The function is gated on `res.type == AminoAcid::CYS`, so the residue context is typed; but inside the residue the SG atom is identified by string. The right typed primitive is `Identity(ai).element == Element::S && Identity(ai).locant == Locant::Gamma` — there's only one S in CYS (`SG`), so element alone narrows it. |
| Risk | MEDIUM. CYS only has one S; element + residue-type sufficiency means the bug surface is small in practice. But the call shape is "find atom-named-X in residue-typed-Y," which is exactly the antipattern. |
| Recommendation | Replace with: `if (atom.element == Element::S)` (CYS only has one) OR `if (Identity(ai).locant == Locant::Gamma && atom.element == Element::S)` (more specific). Even cleaner: `cys_residue.AtomWithLocant(Locant::Gamma, BranchAddress{})` returns the SG. Note that the AMBER-LEAP-INPUT calling-site already lives behind `BondCategory::Disulfide` in production (CovalentTopology resolves disulfides via the GROMACS readback, see `feedback_readback_block_is_a_compiler_trace`). Audit the whether `DetectDisulfides` is even needed post-readback-block; it may be deletable rather than refactorable. |

#### `src/AmberChargeResolver.cpp:306,321` — flat-table coverage verdict

| | |
|---|---|
| Use 306 | `const auto candidates = AtomNameCandidates(atom.pdb_atom_name, res.terminal_state);` — produces alternate spellings (`H` vs `H1`) for terminal-N hydrogens. |
| Use 321 | Records the failed atom name into the verdict for human-readable error message. |
| Classification | **B (boundary-out string-key projection).** The flat ff14SB dat file is keyed by atom-name string triples `(terminal, resname, atom)`. The verdict tests whether the protein's atoms are covered by the flat table's keys. This IS a wire-format lookup. |
| Risk | LOW for line 306 (legitimate projection at the lookup site). MEDIUM for line 321 (the verdict captures the string for diagnostics; if that string is also used by downstream code as identity, propagation grows). |
| Recommendation | Post-pivot, projection happens at the verdict computation: `AmberNameFor(atom.identity, res.type, res.terminal_state)` → string → `AtomNameCandidates` → flat-table lookup. The atom never carries the AMBER-specific name. The verdict's diagnostic string is computed at error time, not stored on the atom. |

#### `src/ChargeSource.cpp:222,243` — `ParamFileChargeSource::LoadCharges`

| | |
|---|---|
| Use 222 | `const auto atom_candidates = AtomNameCandidates(identity.pdb_atom_name, res.terminal_state);` — same pattern as 306 above; per-atom flat-table key projection at load time. |
| Use 243 | Error-message printf when the verdict-loader contract is violated (FATAL abort path). |
| Classification | 222: **B**. 243: **D** (diagnostic only). |
| Risk | 222: LOW (legitimate). 243: LOW. |
| Recommendation | Post-pivot, line 222 calls `AmberNameFor(atom.identity, res.type, res.terminal_state)` → projects → looks up. The function is at the AMBER-table boundary (ff14SB flat file is the wire format) so projection is the right shape. The diagnostic at 243 stays — it's a FATAL abort message; carry the typed identity into the error output too if desired. |

#### `src/AmberPreparedChargeSource.cpp:376,378,382,427` — PRMTOP atom-name matching

| | |
|---|---|
| Use 376-382 | Builds a per-residue `(name → atom_index)` map from `pdb_atom_name`, deduplicating via `m.count(...)`. |
| Use 427 | Error-message use (PRMTOP atom name not found in extractor protein). |
| Classification | **B**. PRMTOP is wire format (Fortran 5E16.8 binary-ish text); its atom names are the keys. The matching is "find this PRMTOP atom in the extractor's residue." |
| Risk | MEDIUM. The match is per-residue scoped (no cross-residue collision), but the dedup check (`m.count(...) > 0`) treats two atoms with the same `pdb_atom_name` in the same residue as a fatal error — which is correct only if `pdb_atom_name` is unique per residue, which the substrate's identity-uniqueness contract gives (within a residue table). |
| Recommendation | Post-pivot: build the map keyed by `AtomMechanicalIdentity` (typed) on the extractor side; on the PRMTOP side, project `prmtop.atom_names[prmtop_ai]` → `AtomMechanicalIdentity` via `ComputeAtomMechanicalIdentity(name, residue_type)` and look up. Wire-format strings stay local to the parser; the cross-side comparison is typed. |

#### `src/PropkaProtonator.cpp:45,47,50` — PDB writer for propka input

| | |
|---|---|
| Use | Formats `pdb_atom_name` into PDB columns 13-16 |
| Classification | **B**. PDB wire format. |
| Risk | LOW (legitimate emit). |
| Recommendation | Post-pivot, compute the AMBER name from typed identity at emit time (Crystal Projection Rule). The function-local `atomField` is built from a function-local string. This was originally `D` (display) but on reflection it's wire — propka reads it. |

#### `src/KamlProtonator.cpp:42,44,47` — PDB writer for kaml input

Same as `PropkaProtonator.cpp`. **B**. Same recommendation.

#### `src/DsspResult.cpp:57,59,62` — PDB writer for DSSP input

| | |
|---|---|
| Use | DSSP requires PDB input; this function writes a temp PDB with the conformation's atoms. |
| Classification | **B**. DSSP is the consumer; PDB is wire. |
| Risk | LOW (DSSP itself is at the boundary). The DsspResult field surface (ss8, hbond_energy, etc.) is typed; the boundary's input is the only string-touching site. |
| Recommendation | Same Crystal-Projection-Rule recommendation. Compute atom name from typed identity at emit. |

#### `src/AIMNet2Result.cpp:169` — error message

| | |
|---|---|
| Use | Error message when an atom has `Element::Unknown` (AIMNet2 has no Z=0 embedding). Identifies the offending atom via `pdb_atom_name + " in residue " + residue_index`. |
| Classification | **D (diagnostic / error message).** |
| Risk | LOW. |
| Recommendation | Survives. Post-pivot, render the typed identity to a string at the error site (e.g., `RenderIdentity(atom.identity, residue_type)` returning a human-readable like "HG (Cys42)"). The render function is for human consumption; calculator code does not consume strings from it. |

#### `src/ApbsFieldResult.cpp:140` — error message

| | |
|---|---|
| Use | Error message when a PB radius is missing for an atom; identifies via `pdb_atom_name`. |
| Classification | **D (diagnostic).** |
| Risk | LOW. |
| Recommendation | Same as AIMNet2: render typed identity to string at error time. **Note that this calculator (APBS) does NOT consume `pdb_atom_name` for any chemistry decision — only for the diagnostic** — confirming that the substrate-side typed surface (`partial_charge`, `pb_radius`) is what calculator logic reads. The error path is the only string touch, and it's display. |

#### `tests/test_object_model.cpp` (5 sites: 28, 29, 30, 31, 44, 69)

| | |
|---|---|
| Use | Test fixture builds a tiny protein with `n->pdb_atom_name = "N";` etc. |
| Classification | **F (test-code reflection of the surface).** |
| Risk | LOW. |
| Recommendation | Post-pivot, construct via typed identity: `n->identity = AtomMechanicalIdentity{Element::N, Locant::None, {}, DiastereotopicIndex::None, BackboneRole::Nitrogen};`. Or, more ergonomic: a helper `MakeBackboneAtom(BackboneRole::Nitrogen)` that fills the fields. The tests will need to migrate when the field comes off `Atom`. |

#### `tests/test_amber_prepared_charge_source.cpp` (45, 104)

| | |
|---|---|
| Use 45 | Test fixture: `atom->pdb_atom_name = templ.name;` from AminoAcidType templates. |
| Use 104 | Test injects `stray->pdb_atom_name = "ZZZZ";` to verify error path. |
| Classification | **F**. |
| Risk | LOW. |
| Recommendation | Post-pivot: AminoAcidType templates carry `AtomMechanicalIdentity`, not strings. Stray-atom test becomes "atom whose identity has no entry in the residue's table" — which the verdict already handles (via the substrate's missing-row signal). The test surface migrates with the production surface. |

#### `tests/test_amber_charge_resolver.cpp` (45, 385)

Same as above. **F**. Same recommendation.

#### `tests/test_amber_leap_input.cpp` (35, 93, 294)

Same. **F**. Same recommendation.

#### `tests/test_foundation_results.cpp` (68)

Same. **F**. Same recommendation.

#### `tests/test_traversal_dump.cpp` (6, 74)

| | |
|---|---|
| Use 6 (comment) | `// pdb_atom_name is used for printing only, never for comparison.` — documentation establishing the test's discipline. |
| Use 74 | `<< " (" << atom.pdb_atom_name << ") has no bonds";` — display-only printf in EXPECT message. |
| Classification | **D**. The whole point of this test is to demonstrate the typed surface; the only `pdb_atom_name` appearance is for human error messages. |
| Risk | LOW. |
| Recommendation | Post-pivot, the printf renders typed identity. The test's discipline-comment becomes "identity is used for navigation; substrate fields for chemistry; strings only for human display." |

#### `tests/bones/src/GromacsEnsembleLoader.cpp` (338, 345, 529, 540)

| | |
|---|---|
| Status | **Quarantined.** CHARMM/GROMACS-via-XTC path retired 2026-05-02 per `project_charmm_retired_amber_only_2026-05-02`. Code lives under `tests/bones/` and is not built into production paths. |
| Classification | A (write) + E (backbone-cache string-keying), same shape as `FullSystemReader.cpp:876,884` |
| Risk | None for production. |
| Recommendation | No action. If the path is ever resurrected (unlikely given the AMBER-only stance), the same recommendations as for `FullSystemReader.cpp` apply. |

### 2.3 Risk hotspots

The 3 highest-priority sites, where the rot is most concentrated:

#### Hotspot 1: `Protein::ResolveProtonationStates` (`Protein.cpp:325-413`)

**The bad thing.** This function infers protonation chemistry
(HID/HIE/HIP, ASH, GLH, CYX/CYS, LYN, TYM) from PDB atom-name
presence. Every titratable residue is tested via:

```cpp
std::map<std::string, size_t> name_to_idx;
for (size_t ai : res.atom_indices) {
    name_to_idx[atom.pdb_atom_name] = ai;
}
// ...
if (res.type == AminoAcid::HIS) {
    bool has_HD1 = name_to_idx.find("HD1") != name_to_idx.end();
    bool has_HE2 = name_to_idx.find("HE2") != name_to_idx.end();
    if (has_HD1 && has_HE2) variant_idx = 2;  // HIP
    else if (has_HD1) variant_idx = 0;        // HID
    else if (has_HE2) variant_idx = 1;        // HIE
    // ...
}
```

There are nine such per-residue branches. Each is a tiny taxonomic
encoding of "what name does this protonation state's H carry?"
spread across the function. The substrate's `LookupBy` with a
partial identity (`Element::H, Locant::Delta, BranchAddress{1, 0}`)
already returns "is HD1 in this residue's table" — typed.

**Why high priority.** Protonation determines variant index, which
drives charge assignment, ring type, and AMBER residue label. A
silent misclassification here propagates through every calculator
downstream. The string match is brittle: any loader producing a
slightly-off name (CHARMM's `HD1` vs `HD1+` if such a thing
existed; future PDB-with-explicit-charge-state) breaks this.

**Status of mitigation.** The trajectory path now sets
`protonation_variant_index` directly via the GROMACS readback
block (committed 2026-05-02); this function's output is overridden
on that path. The PDB and ORCA paths still rely on it.

**Recommendation.** Move the detection logic into the loader (where
the load context is fresh). The loader sets variant_index on
`Residue` directly. `ResolveProtonationStates` becomes
"already-resolved by load; consistency-check only." Even the
consistency-check should branch on identity, not on names.

#### Hotspot 2: `Protein::CacheResidueBackboneIndices` (`Protein.cpp:542-571`)

**The bad thing.** Every backbone cache field (`res.N`, `res.CA`,
`res.C`, `res.O`, `res.H`, `res.HA`, `res.CB`) is populated by
string equality on `pdb_atom_name`. The `H` cache treats
`"H" || "HN"` as equivalents (NamingRegistry alias seeping into
the cache builder). The `HA` cache treats `"HA" || "HA2"` as
equivalents (Gly's HA is HA2 in the ff14SB convention; everyone
else's is HA). Chi angles are resolved by string match against
`ChiAngleDef::atoms[]`.

**Why high priority.** Every spatial-geometry calculator that
asks for the alpha carbon, the backbone N, the carbonyl O — all
of them — depends on these caches being populated. Silent failure
means the calculator sees `Residue::NONE` and either skips the
atom or, worse, produces zeroed kernels.

The aliases are PATTERNS.md's "naming boundary cross once" rule
broken: NamingRegistry.cpp does the canonical translation, then
the alias re-emerges in this function as a string disjunction.
The single source of truth was supposed to be NamingRegistry;
having a second site gives two places where the alias set must
stay synchronised.

**Recommendation.** Replace string match with `BackboneRole`
match: `if (Identity(ai).backbone_role == BackboneRole::Nitrogen)`,
etc. Aliases collapse: `H || HN` → `BackboneRole::AmideHydrogen`
because NamingRegistry already mapped them through. `HA || HA2`
→ `BackboneRole::AlphaHydrogen`. For chi angles, encode
`ChiAngleDef` with identity tuples instead of name strings.

#### Hotspot 3: `Protein::DetectAromaticRings` (`Protein.cpp:431-529`)

**The bad thing.** Every aromatic residue's ring atoms are matched
by string. The `AminoAcidType::rings[].atom_names` is a
`const char*[6]` (or [5] / [9]) of literal AMBER names. The
function builds `name_to_idx` from `pdb_atom_name` and intersects.

Compounding issue: when a HIS residue has no resolved
`protonation_variant_index`, the function falls through to a PDB
LOADING BOUNDARY name match on "HD1" / "HE2" (lines 487-496),
duplicating the protonation-detection logic from
`ResolveProtonationStates`. Two sites, same fragile string match.

**Why high priority.** Ring detection is the foundation of
ring-current calculations (BiotSavartResult, HaighMallionResult,
RingSusceptibilityResult). A miss means zero ring current,
silently. The hardcoded ring atom-names list in `AminoAcidType`
is a separate string surface that has to stay synchronised with
every loader's name normalisation.

**Recommendation.** Encode ring atom membership in `AminoAcidType`
as `AtomMechanicalIdentity` tuples (or partial identities —
`(Element, Locant)` is often sufficient: PHE ring is six C atoms
at locants Gamma / Delta(1,2) / Epsilon(1,2) / Zeta). Resolve via
`res.AtomWithIdentity(...)` instead of `name_to_idx`. The HIS
fallback string match goes away because `LookupBy` already
distinguishes HID/HIE/HIP rows by variant_idx.

### 2.4 Recommended pivot phases

Given the audit, the pivot is sequenceable in seven phases. Each
phase is independently testable; each leaves the codebase in a
working state with a parallel typed surface alongside the existing
string surface, until the final phase deletes the string.

**Phase 1 — Add the typed identity field on `Atom`.**
- New: `Atom.identity` of type `AtomMechanicalIdentity`.
- Populate at construction in every loader: PdbFileReader,
  FullSystemReader, OrcaRunLoader, all three Protonator PDB
  writers (which currently only read the field — they don't
  construct atoms — so no change needed there).
- Helper: `ComputeAtomMechanicalIdentity(string, residue_type,
  is_terminal)` already lifted in
  `LegacyAmberSemanticTables`; use it.
- `pdb_atom_name` stays for the moment; this phase is purely
  additive.
- Estimated lines: ~80 src + ~30 tests.
- Validation: `Identity(ai).backbone_role` agrees with the
  current `name == "N"` / etc. across the standard 20 across all
  test fixtures.

**Phase 2 — Migrate internal `Protein.cpp` builders to typed.**
- `CacheResidueBackboneIndices`: typed.
- `DetectAromaticRings`: typed (ring atom-name → identity tuple
  in `AminoAcidType`; this is its own internal change).
- `ResolveProtonationStates`: typed (LookupBy partial identity
  for atom-presence detection).
- `pdb_atom_name` continues to exist but stops being the keying
  surface for any internal builder.
- Estimated lines: ~150 src + the AminoAcidType ring-data change.
- Validation: every existing test passes; the 1UBQ traversal-dump
  test continues to count rings correctly.

**Phase 3 — Migrate calculator-side identity-proxy hacks
(Category C).**
- `AmberLeapInput::DetectDisulfides:328` → typed S-by-locant.
- AIMNet2 + ApbsField error messages → keep `pdb_atom_name`
  use as render-only OR migrate to `RenderIdentity` helper.
- Anything else discovered during audit that this document missed
  — cycle through Category C grep.
- Estimated lines: ~30 src.
- Validation: smoke tests bit-identical (the calculator
  computations don't change — just how they identify atoms).

**Phase 4 — Migrate file-format emitters (Category B
write-side).**
- AmberLeapInput PDB emitter: compute name from typed identity
  via `AmberAtomNameFor(identity, residue.type, terminal_state)`
  at emit time.
- Three Protonator PDB writers: same.
- DsspResult PDB writer: same.
- TrajectoryProtein H5 emit: emit BOTH the typed substrate
  columns AND the legacy `pdb_atom_name` column (transition
  period); flag the latter deprecated. Update SDK catalog.
- Estimated lines: ~80 src.
- Validation: PDB outputs byte-identical for standard fixtures
  (1Z9B, 1P9J, 1UBQ); tleap accepts the regenerated PDB.

**Phase 5 — Migrate external-table lookups (Category B
read-side).**
- AmberChargeResolver verdict: project at lookup site.
- ChargeSource ParamFileChargeSource: project at lookup site.
- AmberPreparedChargeSource PRMTOP match: project at the
  cross-side comparison.
- The flat ff14SB dat file and PRMTOP are unchanged; only the
  extractor's usage of them changes.
- Estimated lines: ~60 src.
- Validation: charge assignment bit-identical on all standard
  fixtures (this is the existing `ChargeFF14SBVariantTest` and
  similar suites).

**Phase 6 — Migrate loader writes to be local-scoped only.**
- PdbFileReader, FullSystemReader, OrcaRunLoader: the
  `atom_name` string is used to compute identity, then discarded
  at function scope (does not write to `Atom.pdb_atom_name`).
- This is when the internal load-side strings stop existing on
  Atom.
- Estimated lines: ~40 src.
- Validation: all existing tests pass; H5 emit still works
  because it now reads the typed identity.

**Phase 7 — Delete `Atom.pdb_atom_name`.**
- Remove field declaration.
- Remove all remaining test-fixture writes (Category F sites).
- Remove the `tests/bones/` GromacsEnsembleLoader writes (or
  leave; it's quarantined and out of build).
- Update PdbFileReader.h header documentation.
- Update OBJECT_MODEL.md `Atom` section — replace
  `pdb_atom_name` with `identity` of type
  `AtomMechanicalIdentity`.
- Estimated lines: ~50 deletions + ~30 doc updates.
- Validation: project builds; full ctest passes.

**Total estimated diff:** ~600 lines across the seven phases, of
which 250-300 lines are core typed substitution and the rest is
test-fixture migration + doc updates. The substrate already
exists; this audit found no missing typed primitives — just
missing pivot.

**Sequencing notes:**

- Phases 1-2 must precede 3-5; the typed surface has to exist
  before consumers migrate to it.
- Phases 3-5 can run in parallel if the substrate population is
  stable.
- Phase 6 is gated on Phase 4 (emitters need the surface, not
  the field).
- Phase 7 is the final consequence; the field comes off only
  when no consumer remains.
- Each phase is committable independently. The repo's git
  discipline (commits per topic) supports this.

---

## Closing notes

A few tensions that surfaced during this audit and that the user
should resolve before pivot starts:

1. **Identity vs partial identity for graph lookup.** Constructing
   a full `AtomMechanicalIdentity` literal in calculator code is
   verbose: 5 fields, most of which are `None` for any specific
   query. The right surface is probably partial-identity accessors
   (`AtomWithRole`, `AtomWithLocant`, `AtomWithDiastereotopicIndex`)
   on `Residue`, with `AtomWithIdentity` reserved for the
   identity-equality comparison case. This document's
   recommendations assume that surface; the user should bless or
   reshape it before Phase 2 lands.

2. **Substrate population timing.** §H.5 of the
   topology-encoding-dependencies doc specifies that `LookupBy`
   runs at `Protein::FinalizeConstruction`. The current
   implementation does not yet do this — the substrate tables
   exist (`LegacyAmberSemanticTables.cpp` is regenerated and
   landed at `ee1f1b4`), but `FinalizeConstruction` does not call
   `LookupBy` to populate per-atom semantic rows. This document
   has assumed the population is in scope of the pivot; the user
   should confirm. If population is its own piece of work, Phase
   1 above grows by ~50 lines (the population loop in
   `FinalizeConstruction`) and ~20 lines (the `Atom.semantic`
   field declaration, plus a `Residue.AtomWithIdentity` accessor).

3. **Ring `atom_names` in `AminoAcidType`.** The change to encode
   ring membership as identities rather than strings is a
   structural change to `AminoAcidType` (a frequently-read
   table). Phase 2 assumes this lands; the alternative is to
   keep `atom_names` as strings and project at use site, which
   leaves the strings inside the typed-table data. The user
   should decide whether `AminoAcidType` is part of the pivot or
   stays string-keyed (with the rule that the strings never
   escape it).

4. **Chi-angle definitions in `AminoAcidType`.** Same shape as
   above — `ChiAngleDef::atoms[4]` is `const char*`. Phase 2
   assumes typed identity; alternative is per-site projection.

5. **`tests/bones` GromacsEnsembleLoader.** Quarantined and out
   of build; this audit notes its existence but recommends no
   action. The user should confirm the quarantine status is
   stable (it is, per memory entry
   `project_charmm_retired_amber_only_2026-05-02`); if the
   CHARMM path is ever revived, its load-side string match
   needs the same treatment as `FullSystemReader.cpp`.

6. **SDK contract for `pdb_atom_name` H5 column.** The Python
   SDK reader has `pdb_atom_name` as a known column. Removing it
   in Phase 4 (or 7) requires updating the SDK catalog. The
   substrate columns (`/atoms/legacy_amber/`) are the new
   answer per OBJECT_MODEL §H5; sequencing the SDK migration
   alongside is a non-trivial cross-project commit. The user
   should plan this — it is downstream consumer-impact, not
   library-internal.

7. **The "AmberAtomNameFor" emitter helper.** Phase 4 introduces
   `AmberAtomNameFor(identity, residue_type, terminal_state)`
   per the Crystal Projection Rule. AMBER has its own atom-name
   conventions (HG vs HG1, OH vs HO, etc. for SER; CYS HG vs
   omitted for CYX) that diverge from CCD/IUPAC. The function
   needs careful spec; it is essentially the inverse of
   `ComputeAtomMechanicalIdentity`. A single `AtomNameProjections`
   header housing IUPAC, AMBER, BMRB rendering functions is the
   natural shape. The user should decide whether one helper
   per format or a unified `Render(identity, format_enum)`
   surface — both are typed, the difference is ergonomics.

These are not blockers; they are decisions the user should make
before the pivot begins, so the next agent has a stable target.
The audit's job is done: it has named every site, classified
each, and proposed a phasing. The next session implements; this
document is what they read first.
