# Chemistry Question 2 — Backbone topology indexing (2026-05-05)

## Summary

The backbone is a typed concept: `BackboneRole` enumerates six slots
(Nitrogen, AlphaCarbon, CarbonylCarbon, CarbonylOxygen, AmideHydrogen,
AlphaHydrogen) plus the `Locant::Alpha` shape for Gly's HA2/HA3 and
the no-H Pro case. Yet across the runtime library, indexing INTO that
typed concept is done by string equality on `Atom::pdb_atom_name`,
with hand-rolled `H || HN` and `HA || HA2` (and in OrcaRunLoader,
`H || HN || H1`) alias disjunctions repeated across **four** sites:
three loaders plus `Protein::CacheResidueBackboneIndices`. The typed
substitute (`BackboneRole`) is fully designed and lives in
`SemanticEnums.h`, but no `Atom::backbone_role` field exists yet at
runtime — only the build-time generator in `tools/topology/` computes
it. Phase 1 must (a) lift `BackboneRole` derivation into the runtime
library so every loader can compute it from a name on entry, (b)
store it on `Atom`, (c) re-key the four backbone-cache populators on
the typed field. Phase 2 may unify backbone access through typed
accessors (`res.alpha_carbon()`, `Residue::AtomWithRole(...)`).

## 1. The chemistry question

**How does the model find a residue's backbone N / CA / C / O / H /
HA atoms?**

The peptide-amide unit (Pauling, Corey, Branson, PNAS 37 (1951)
205-211) defines six canonical backbone slots. Every standard residue
has all six EXCEPT:

- Proline lacks the backbone amide H (secondary amine — chain N
  carries no H in internal Pro). `BackboneRole::AmideHydrogen` is
  absent.
- Glycine has two prochiral alpha hydrogens HA2/HA3 instead of a
  single HA. By Markley 1998 convention these carry `Locant::Alpha`
  and `DiastereotopicIndex::Position2` / `Position3`, with
  `BackboneRole::None` — the (Locant, DiIndex) pair already
  disambiguates them. (This is documented inline in `SemanticEnums.h`
  lines 86-88 and reproduced in `tools/topology/build_semantic_tables.cpp:519-521`.)

The chemistry question — "which atom in this residue is the alpha
carbon?" — has a typed answer: the atom whose `BackboneRole ==
AlphaCarbon`. The substrate row produced by the build-time generator
already carries that answer (`AtomSemanticTable::backbone_role` field
at `SemanticEnums.h:790`). What does NOT yet exist is any runtime-
library code that puts that answer onto `Atom` instances when a
loader runs. So when `Protein::CacheResidueBackboneIndices` and the
three load-site mini-caches need to find the alpha carbon, they fall
back to the only field they have: the raw PDB-style string name.

The gap, sharply: typed backbone identity is fully defined; typed
backbone identity is not yet on `Atom` at construction time; the
backbone-cache populator runs at construction time. Calculators
downstream then read `res.CA` (a `size_t` index) as their typed
surface — and that surface IS typed; it's the populator that hides
strings under the hood.

## 2. What the current code does

### 2.1 `Protein::CacheResidueBackboneIndices` (src/Protein.cpp:542-571)

The populator. Called from `FinalizeConstruction` (line 218) after
`ResolveResidueTerminalStates` and before `ResolveProtonationStates`.
Iterates each residue's atom list; for each atom, walks a chain of
`name == "<literal>"` comparisons:

```cpp
const std::string& name = atoms_[ai]->pdb_atom_name;
if      (name == "N")   res.N  = ai;
else if (name == "CA")  res.CA = ai;
else if (name == "C")   res.C  = ai;
else if (name == "O")   res.O  = ai;
else if (name == "H" || name == "HN")  res.H  = ai;
else if (name == "HA" || name == "HA2") res.HA = ai;
else if (name == "CB")  res.CB = ai;
```

Two alias disjunctions: `H || HN` (CHARMM uses `HN` for the backbone
amide H; AMBER and IUPAC use `H`) and `HA || HA2` (Glycine's
diastereotopic HA2/HA3, where HA2 is conventionally the IUPAC
"forward" pro-R H). These are alias sets that NamingRegistry
(`src/NamingRegistry.cpp:138-139`) is supposed to collapse at the
load boundary:

```cpp
AddAtomNameRule("H",  "HN", "*", ToolContext::Standard, ToolContext::Charmm);
AddAtomNameRule("HN", "H",  "*", ToolContext::Charmm, ToolContext::Standard);
```

That NamingRegistry rule should mean post-load every backbone amide
H carries the canonical `H` name, and the `|| name == "HN"` clause
in `CacheResidueBackboneIndices` would never fire. **That the clause
is still there indicates the alias collapse does not reach all load
paths** (or did not, at the time the line was written, and the
clause is defensive scar tissue). This is a finding; it is not
benign. PATTERNS.md §"The naming boundary: cross once, then typed
objects forever" forbids exactly this duplication.

The chi-angle resolver in the same function (lines 559-569) does
similar string matching against `ChiAngleDef::atoms[]` (a
`const char*[4]` literal table). Out of scope for backbone topology
specifically, but adjacent — flagged in audit Hotspot 2.

### 2.2 Backbone-cache code at the three loader sites

**`src/FullSystemReader.cpp:884-896`** (TPR / GROMACS-via-AMBER trajectory):
duplicates the same backbone-cache populator INLINE in the load loop.
Reads the just-written `pdb_atom_name` (which has been NamingRegistry-
translated CHARMM→Standard at line 876) and matches:

```cpp
if      (name == "N"  && res_ref.N  == Residue::NONE) res_ref.N  = idx;
else if (name == "CA" && res_ref.CA == Residue::NONE) res_ref.CA = idx;
else if (name == "C"  && res_ref.C  == Residue::NONE &&
         protein->AtomAt(idx).element == Element::C &&
         res_ref.CA != Residue::NONE) res_ref.C = idx;
else if (name == "O"  && res_ref.O  == Residue::NONE) res_ref.O  = idx;
else if ((name == "H" || name == "HN") &&
         res_ref.H == Residue::NONE) res_ref.H = idx;
else if ((name == "HA" || name == "HA2") &&
         res_ref.HA == Residue::NONE) res_ref.HA = idx;
else if (name == "CB" && res_ref.CB == Residue::NONE) res_ref.CB = idx;
```

**`src/OrcaRunLoader.cpp:242-252`** (ORCA XYZ):

```cpp
if (name == "N" && res.N == Residue::NONE) res.N = idx;
else if (name == "CA" && res.CA == Residue::NONE) res.CA = idx;
else if (name == "C" && res.C == Residue::NONE &&
         protein->AtomAt(idx).element == Element::C &&
         res.CA != Residue::NONE) res.C = idx;
else if (name == "O" && res.O == Residue::NONE) res.O = idx;
else if ((name == "H" || name == "HN" || name == "H1") &&
         res.H == Residue::NONE) res.H = idx;
else if ((name == "HA" || name == "HA2") &&
         res.HA == Residue::NONE) res.HA = idx;
else if (name == "CB" && res.CB == Residue::NONE) res.CB = idx;
```

OrcaRunLoader carries an EXTRA alias the others don't: `H || HN ||
H1`. The `H1` is the N-terminal first ammonium H — at an internal
residue this should never match; at the N-terminus the ammonium H1
is a cap atom that `BackboneRole::AmideHydrogen` does NOT cover (it's
distinct chemistry — see `kCapNtermCharged` table). Treating H1 as
the backbone amide H is a category error that probably works on most
fixtures but encodes confusion.

**`tests/bones/src/GromacsEnsembleLoader.cpp:347-348, 542-543`**:
quarantined CHARMM/XTC path (retired 2026-05-02). Same pattern;
out-of-scope but completeness-noted.

**`src/PdbFileReader.cpp:130-138`**: writes `pdb_atom_name` from
cifpp (line 133); does NOT populate the backbone cache directly.
Relies on `FinalizeConstruction` → `CacheResidueBackboneIndices`
being called at line 148 after all atoms are added. So the PDB path
goes through `Protein::CacheResidueBackboneIndices` (§2.1); the TPR
and ORCA paths populate the cache *also* during their load loops AND
also again at FinalizeConstruction. The redundancy is a separate
finding (audit didn't flag it; worth flagging here): the load-site
caches and the post-load `CacheResidueBackboneIndices` overlap.
Either the load-site caches are dead code, or `Residue::NONE`-guarded
double-write makes the second invocation a no-op. Either way: the
single source of truth principle is broken three ways over.

### 2.3 `src/DsspResult.cpp` — pdb_atom_name reads

`DsspResult.cpp` has three `pdb_atom_name` reads, all in
`WriteTempPdb` (lines 30-78), and **all of them are file-format
emit, not backbone navigation**:

- Line 57: `if (atom.pdb_atom_name.size() <= 3)` — width branch for
  PDB column 13-16 formatting.
- Line 59: `snprintf(... " %-3s", atom.pdb_atom_name.c_str())` —
  emits the name into the temp PDB.
- Line 62: `snprintf(... "%-4s", atom.pdb_atom_name.c_str())` —
  same, four-character form.

These are Category B in the audit's taxonomy (boundary-out string-
key projection). DSSP itself, the cifpp-and-libdssp wrapper that
reads the temp PDB, lives at the PDB-file boundary — once the temp
PDB is read by cifpp, DSSP's residue-level results are mapped back
to typed `residue_index` via `(chain, sequence_number)` (line 144)
WITHOUT any further atom-name reading. The hydrogen-bond donor/
acceptor results come back as residue-pair tuples from libdssp's
internal navigation (which DOES use atom names internally — but
that's libdssp's business, not ours). DSSP at the calculator
surface here does NOT use `pdb_atom_name` for backbone navigation;
all its "find the H-bond donor H of this residue" logic happens
inside libdssp on cifpp's structure.

This is important to call out separately. The audit's mention of
"DSSP processor that uses backbone atom name lookups" was at the
PDB write boundary, not at the C++ navigation layer. Phase 1 work
on backbone topology indexing does NOT need to change `DsspResult.cpp`'s
internal logic; it only needs to change the AMBER/IUPAC name being
written into the temp PDB (which becomes "compute name from typed
identity" per the Crystal Projection Rule, but that's a Phase 4
emitter concern, not Phase 1 backbone indexing).

The take: DSSP reads `pdb_atom_name` only as a wire-format value;
post-pivot it computes the wire format from typed identity at emit
time. NOT a Phase 1 blocker for the backbone-topology question.

## 3. Where typed information already lives

- **`BackboneRole` enum** (`src/SemanticEnums.h:89-99`). Six values
  plus `None`. Cited to Pauling 1951 + Ramachandran 1968. This IS
  the typed answer to "which backbone slot."

- **`AtomSemanticTable::backbone_role`** field (`src/SemanticEnums.h:790`).
  Every atom in every per-residue / per-variant / per-cap table
  carries its `BackboneRole` value. The build-time generator
  (`tools/topology/build_semantic_tables.cpp:519-535`) computes it
  from the canonical AMBER atom name + parent-atom map.

- **`AtomMechanicalIdentity`** (`src/SemanticEnums.h:848-863`). The
  5-tuple lookup key including `backbone_role`. Used by
  `LookupBy`/`LookupCap` at construction time to find a substrate
  row given typed identity.

- **`Residue` backbone-cache fields** (`src/Residue.h:48-55`):
  `static constexpr size_t NONE`, `size_t N, CA, C, O, H, HA, CB`.
  These are typed `size_t` indices (a thin typed surface).
  Calculators consume them by name: `EnrichmentResult.cpp:36-41`,
  `CoulombResult.cpp:68-74`, `CovalentTopology.cpp:94-96`,
  `AmberLeapInput.cpp:266,296`. The downstream surface is already
  typed; only the population side leaks strings.

- **`AminoAcidType::atoms`** chain inventory (`src/AminoAcidType.cpp`).
  The `BB` macro defines six backbone atoms by AMBER-canonical name:
  `BB = N, CA, C, O, H, HA`. `BB_PRO` drops `H`. `BB_GLY` swaps
  HA for HA2 + HA3. Today these are `const char*` literals; they
  are project-internal canonical names, not "atom-name strings the
  user provided." The macro is the single inventory authority for
  what backbone atoms exist per residue type.

- **`CovalentTopology` bond graph** (`src/CovalentTopology.{h,cpp}`).
  Resolved in `FinalizeConstruction` AFTER `CacheResidueBackboneIndices`.
  In principle the backbone N/CA/C/O could be discovered topologically
  (find the N bonded to a heavy atom across the prior residue's
  carbonyl C; the CA bonded to that N and to a carbonyl C; the
  carbonyl C bonded to that CA and to a carbonyl O and to the next
  N; etc.). This is a fallback path, not the primary one.

- **`NamingRegistry`** (`src/NamingRegistry.{h,cpp}`). Has rules
  for `H ↔ HN` collapsing (lines 138-139). Has known coverage gaps
  documented in a multi-finding comment block (lines 153-200) for
  ILE δ-methyl, GLY α-methylene, ARG/LYS/PRO δ-methylene, LYS ε-
  methylene, plus an ALA β-methyl wildcard bug. The comment
  documents these are deferred pending fleet-wide vetting. The
  H/HN canonicalization specifically IS active and IS supposed to
  reach FullSystemReader.cpp (which calls it explicitly at line
  876). The fact that `Protein::CacheResidueBackboneIndices` still
  has `name == "H" || name == "HN"` after NamingRegistry should
  have collapsed them is the smoking-gun of incomplete naming-
  boundary discipline. (Possibility: PdbFileReader does not call
  NamingRegistry at all — just writes the cifpp-provided name raw,
  per `src/PdbFileReader.h:12`. So "HN" can still arrive on a PDB
  load if the source PDB used the CHARMM convention.)

What's missing:

- **No runtime `ComputeAtomMechanicalIdentity` function.** The
  parser logic that derives `BackboneRole` from a string name lives
  ONLY in `tools/topology/build_semantic_tables.cpp` — the
  build-time generator binary, which doesn't link into the runtime
  library. The audit at line 768 says "use it" referring to a
  helper in `LegacyAmberSemanticTables`; that helper does not
  exist. `LegacyAmberSemanticTables.h` contains only `LookupBy`,
  `LookupCap`, `ApplyCapDelta`. The parser is not lifted.

- **No `Atom::backbone_role` (or `Atom::identity`) field** on
  the runtime `Atom` (`src/Atom.h:20-45`). The class has `element`,
  `pdb_atom_name` (the string we want gone), `residue_index`,
  `bond_indices`, `parent_atom_index`. Nothing typed about identity.

- **No `Protein::FinalizeConstruction` call to `LookupBy`/`LookupCap`.**
  §H.5 of the topology dependencies doc (the LookupBy invocation
  example at line 467-490) describes this composition but
  `FinalizeConstruction` (Protein.cpp:213-262) does not currently
  do it. The substrate population timing is the §H.5 spec target
  but is not implemented.

These three gaps mean: **the typed answer exists in tables, but no
code on the runtime side reads those tables to put a `BackboneRole`
on an `Atom`.** Phase 1 must close these.

## 4. The order-of-operations question

Here is the chicken-and-egg, named.

Today's order at `FinalizeConstruction`:

1. Loader writes `Atom::pdb_atom_name` from external bytes (PDB,
   TPR, ORCA XYZ).
2. (TPR/ORCA paths only) Load-site backbone cache populates
   `res.N`, `res.CA`, etc. from string match.
3. Loader calls `protein->FinalizeConstruction(positions, ...)`.
4. `ResolveResidueTerminalStates` — uses chain_id/sequence; no
   strings.
5. `CacheResidueBackboneIndices` — string match again, redundant
   with step 2 for TPR/ORCA paths, primary for PDB path.
6. `ResolveProtonationStates(false)` — string match on `HD1`/`HE2`/
   `SG`/`HG`/etc. (out of scope for THIS chemistry question; called
   out as Hotspot 1 in the audit).
7. `DetectAromaticRings` — string match against `AminoAcidType::rings[].atom_names`
   (audit Hotspot 3).
8. `CovalentTopology::Resolve` — geometric/connectivity inference,
   uses element + position + bond perception. NO string use.
9. (Conditional) `OverrideDisulfides` from Amber readback.
10. Construct `LegacyAmberTopology`.
11. `ResolveProtonationStates(true)` — same as step 6 with
    disulfide info.

**The chicken-and-egg.** To populate `Atom::backbone_role` from a
typed parser at step 1, the parser has to exist in the runtime
library. Today it does not. Once it does, step 5 becomes:

```cpp
for (auto& res : residues_) {
    for (size_t ai : res.atom_indices) {
        switch (atoms_[ai]->backbone_role) {  // typed
            case BackboneRole::Nitrogen:        res.N  = ai; break;
            case BackboneRole::AlphaCarbon:     res.CA = ai; break;
            case BackboneRole::CarbonylCarbon:  res.C  = ai; break;
            case BackboneRole::CarbonylOxygen:  res.O  = ai; break;
            case BackboneRole::AmideHydrogen:   res.H  = ai; break;
            case BackboneRole::AlphaHydrogen:   res.HA = ai; break;
            default: break;
        }
        // CB falls out of the AlphaCarbon-bonded-Cβ rule via
        // CovalentTopology, OR via Locant::Beta + Element::C +
        // BranchAddress {0,0} on a sidechain non-Gly residue.
        // (CB is not in BackboneRole; it sits between regimes.)
    }
}
```

For this to work, `Atom::backbone_role` has to be populated BEFORE
`CacheResidueBackboneIndices` runs.

There are two places that population can happen:

**Option A — at the loader (per-atom, name-keyed parse).** Each
loader, after reading the external atom name and before calling
`protein->AddAtom`, calls a runtime helper `BackboneRoleFromName(name,
residue_type, terminal_state) → BackboneRole`. The helper lives in
the runtime library (lifted from `tools/topology/build_semantic_tables.cpp:511-535`)
and is just a switch on canonical names. The string-to-role
translation happens once per loader; the result is typed onto
`Atom`. **`CacheResidueBackboneIndices` then runs typed.**

**Option B — at FinalizeConstruction via the substrate.**
`FinalizeConstruction` walks atoms, computes a partial
`AtomMechanicalIdentity` from the name (using the lifted parser),
calls `LookupBy(res.type, res.protonation_variant_index, identity)`,
copies the row's `backbone_role` field onto the atom. Then
`CacheResidueBackboneIndices` runs typed.

**Tension.** Option A puts the parse at every load site (three
sites now); Option B puts it once in `FinalizeConstruction`. Option
B is more aligned with the audit's "loaders compute identity once
at boundary" framing AND with §H.5 of the topology dependencies
doc. But Option B requires `LookupBy` to actually run — and §H.5
notes that runtime composition is "the next-session deliverable;
the substrate data is ready for it" (§H.11). So Option B has a
prerequisite (substrate population) that Option A does not.

**Resolution proposal.** Phase 1 backbone-topology refactor uses
Option A: a runtime `BackboneRoleFromName(name, residue_type)` helper
that loaders call when they read each atom. This gives backbone-
indexing typing WITHOUT requiring the full substrate-population
machinery (`LookupBy` at FinalizeConstruction) to land first. When
substrate population DOES land (the §H.5 work, separately
sequenced), `Atom::backbone_role` becomes redundant with
`Atom::semantic.backbone_role` and can either be retained as a
fast-access copy OR resolved through the substrate. That is a
Phase 1.5 concern; it doesn't block Phase 1.

Concretely: lift the `BackboneRole` derivation from
`tools/topology/build_semantic_tables.cpp:519-535` into a small
runtime header (e.g. `src/AtomNameParser.h` or extend
`src/NamingRegistry.h`). Three lines at the top of each loader's
atom-add loop produce the typed value; `Atom::backbone_role`
stores it; `CacheResidueBackboneIndices` switches on it.

The chi-angle resolver (lines 559-569) is a separate question —
its inputs are `ChiAngleDef::atoms` strings, not Atom names, so
re-keying it requires changing `AminoAcidType` table data. Audit
Hotspot 2 flags this; it can be deferred from this Phase-1 surface
without blocking the backbone-cache typing.

## 5. Phase 1 design

The smallest refactor that keeps every existing calculator working
and removes string-keying from backbone indexing:

**Step 1.** Add `BackboneRole backbone_role = BackboneRole::None;`
to `Atom` (`src/Atom.h`). Field defaults to None for safety.

**Step 2.** Lift `BackboneRoleFromName` into the runtime library. A
new file (or addition to `NamingRegistry`) declares:

```cpp
BackboneRole BackboneRoleFromName(const std::string& name,
                                   AminoAcid residue_type);
```

The implementation is a small switch matching the canonical-name
patterns used in `tools/topology/build_semantic_tables.cpp:519-535`
(for the standard Pauling six) plus the Gly HA2/HA3 special case
(which return `BackboneRole::None` because they carry the locant).
NTERM-cap H1/H2/H3 return `BackboneRole::None` (they are NOT the
backbone amide H — distinct chemistry per `kCapNtermCharged`).

**Step 3.** Each loader, when adding an atom, calls
`atom->backbone_role = BackboneRoleFromName(atom_name, residue_type);`
right after writing `atom->pdb_atom_name`. This adds one call per
loader at:

- `src/PdbFileReader.cpp:133` (after `pdb_atom_name = atom_name`).
- `src/FullSystemReader.cpp:876` (after the NamingRegistry
  translation).
- `src/OrcaRunLoader.cpp:206` (XYZ atom add).
- (Quarantined: `tests/bones/src/GromacsEnsembleLoader.cpp` —
  out-of-build; no action.)

The callable shape is uniform across the three; the loaders no
longer differ in their alias-handling because `BackboneRoleFromName`
is the single source of truth. The `H1` extra-alias bug in
OrcaRunLoader (line 248) goes away — `BackboneRoleFromName("H1",
...)` returns `None`, which is the correct chemistry call (H1 is a
cap atom, not a backbone amide H).

**Step 4.** Rewrite `Protein::CacheResidueBackboneIndices`
(`src/Protein.cpp:542-554`) to switch on
`atoms_[ai]->backbone_role` instead of the name string. The
function body becomes:

```cpp
void Protein::CacheResidueBackboneIndices() {
    for (auto& res : residues_) {
        for (size_t ai : res.atom_indices) {
            switch (atoms_[ai]->backbone_role) {
                case BackboneRole::Nitrogen:        res.N  = ai; break;
                case BackboneRole::AlphaCarbon:     res.CA = ai; break;
                case BackboneRole::CarbonylCarbon:  res.C  = ai; break;
                case BackboneRole::CarbonylOxygen:  res.O  = ai; break;
                case BackboneRole::AmideHydrogen:   res.H  = ai; break;
                case BackboneRole::AlphaHydrogen:   res.HA = ai; break;
                default: break;
            }
        }
        // CB: not in BackboneRole. Either fall back to identity-
        // based partial-match (Locant::Beta + Element::C + branch
        // {0,0}) or keep one defensive name-match line. See §7
        // user-decision below.

        // Chi angles: still strings against ChiAngleDef. Out of
        // scope for Phase 1; audit Hotspot 2 covers separately.
        const AminoAcidType& aatype = res.AminoAcidInfo();
        for (int ci = 0; ci < aatype.chi_angle_count && ci < 4; ++ci) {
            const ChiAngleDef& def = aatype.chi_angles[ci];
            for (int j = 0; j < 4; ++j) {
                for (size_t ai : res.atom_indices) {
                    if (atoms_[ai]->pdb_atom_name == def.atoms[j]) {
                        res.chi[ci].a[j] = ai;
                        break;
                    }
                }
            }
        }
    }
}
```

**Step 5.** Delete the redundant load-site backbone caches in
`FullSystemReader.cpp:884-896` and `OrcaRunLoader.cpp:242-252`.
After step 4 the post-FinalizeConstruction populator handles every
load path uniformly; the per-loader inline caches were defensive/
redundant.

**Step 6.** Verify `CovalentTopology::Resolve` (which runs AFTER
`CacheResidueBackboneIndices` in FinalizeConstruction) gets the
backbone indices it expects (`CovalentTopology.cpp:94-96` reads
`res.N`, `res.CA`, `res.C`). No code change here; it already reads
the typed `size_t` indices.

**What does NOT change in Phase 1.**

- `Atom::pdb_atom_name` stays — it's still used by the chi-angle
  resolver, by the protonation-state inference (Hotspot 1), by the
  ring-detection (Hotspot 3), by AmberLeapInput's PDB emitter,
  by DsspResult's PDB emitter, and so on. Those are separate
  pivots; this Phase 1 only kills the backbone-indexing string use.

- All downstream calculators that read `res.N`, `res.CA`, `res.C`,
  etc. (CoulombResult, EnrichmentResult, CovalentTopology, AmberLeapInput,
  HBondResult, AIMNet2Result, MopacCoulombResult, BsT0Autocorrelation
  trajectory) continue to work unchanged. Their typed surface
  (`size_t`-indexed cache fields) is the same.

- DsspResult's internal logic is untouched. Its `pdb_atom_name`
  reads are wire-format PDB emit, not backbone navigation; they
  don't move in Phase 1.

**Estimated diff.** ~50 lines added (the runtime helper +
`Atom::backbone_role` field + three loader call sites), ~25 lines
deleted (two inline caches + the seven name-match branches in
`CacheResidueBackboneIndices`). Net ~+25 lines. One new symbol on
`Atom`; one new symbol in the runtime library (`BackboneRoleFromName`
or `ParseBackboneRole`).

**Validation.** Existing tests pass: the typed indexing produces
the same `size_t` values as the string indexing on every standard
fixture (1Z9B, 1P9J, 1UBQ, 1A6J, etc.). Add a test asserting that
on a CHARMM-named PDB load (atom name `HN`), `res.H` is populated;
on an AMBER-named load (atom name `H`), `res.H` is populated. Both
go through the typed channel; neither requires the alias disjunction
in `CacheResidueBackboneIndices`. Add a test on Glycine residues
asserting `res.HA == residue's HA2 atom index` (which the typed
parser must enforce — Gly has BackboneRole::None on HA2/HA3, so
the case statement must NOT match them; instead Gly's HA cache
lookup needs the Locant::Alpha + DiastereotopicIndex::Position2
identity, OR a hand-coded "if Gly, use HA2" check). This is a
sub-decision flagged in §7.

## 6. Phase 2 considerations (flagged, not solved)

Once `Atom::backbone_role` is typed and the cache populator is
typed, calculators continue to read `res.CA` (a `size_t`). A future
walk might:

- **Add typed accessors on Residue.** `res.alpha_carbon() →
  size_t` (or `optional<size_t>`), `res.carbonyl_oxygen()`, etc.
  These hide the cache-field naming convention. Calculators read
  through accessors, not raw fields. Migration is mechanical;
  benefit is consistency with `Residue::AtomWithRole(BackboneRole)`
  for non-cache cases.

- **Introduce a `ResidueBackbone` struct.** A small typed bundle
  `struct ResidueBackbone { size_t N, CA, C, O, H, HA; }` with
  constexpr `NONE` slots. `Residue::backbone() → const
  ResidueBackbone&` returns it. Calculators do `r.backbone().CA`
  instead of `r.CA`. Pro and Gly variants are visible at type
  level (separate typed structs, or const-references).

- **Generalise to `AtomWithRole(BackboneRole)`.** When the
  substrate population (§H.5) lands, every atom has a typed
  `backbone_role` (or `semantic.backbone_role`); the accessor
  scans `res.atom_indices` for the matching role. This makes the
  cache redundant — the cache becomes an optimisation, not the
  primary surface. Probably not worth pursuing because the cache
  IS fast and there are only ~30K atoms even on large fixtures;
  the linear scan is fine.

- **Cross-residue backbone navigation.** Calculators that need
  the i-1 carbonyl C or the i+1 N currently have no typed primitive;
  they'd read `protein.ResidueAt(ri+1).N`. Phase 2 might add
  `Protein::PreviousResidue(size_t ri) → const Residue*` and
  similar, which is independent of backbone-role typing.

- **Gly HA disambiguation cleanup.** Today `res.HA` for Gly holds
  the index of HA2 (per the `name == "HA2"` alias in the current
  code). Phase 2 might split this: `res.HA` for non-Gly, `res.HA2
  / res.HA3` for Gly, surfacing the chemistry difference. Or
  `res.HA` always returns the pro-R alpha hydrogen (HA in non-Gly,
  HA3 in Gly per Markley Fig 1). This is a chemistry decision
  about what calculators should see, not a typing question. Flag.

- **Pro special-case visibility.** `res.H` is `Residue::NONE` for
  Pro, which is correct chemistry (no backbone amide H). Calculators
  that read `res.H` as if it's always present have latent bugs.
  Phase 2 might tighten the type: `optional<size_t>` instead of
  `size_t` with NONE sentinel. This is a wider refactor (touches
  every consumer); flag.

These are notes; do not solve.

## 7. Choices that need a user decision

**Decision 1: Field on Atom — `backbone_role` only, full
`identity`, or `semantic` row?**

Three shapes of typed surface on `Atom`:

(a) `BackboneRole backbone_role` only. Smallest, most targeted.
    Solves the backbone-indexing question and nothing else. Other
    chemistry questions (protonation, rings, polar-H taxonomy)
    keep using `pdb_atom_name` until their respective Phase 1s
    land.

(b) `AtomMechanicalIdentity identity` (the full 5-tuple). Solves
    backbone indexing AND prepares for the §H.5 substrate
    population (when that lands, `LookupBy` keys on this struct
    directly). More machinery up front; demands lifting the full
    `ParseAtomName` parser into the runtime library, not just the
    backbone-role-name switch.

(c) `const AtomSemanticTable* semantic` (or stored-by-value).
    Solves everything, but requires `LookupBy` at
    FinalizeConstruction, which §H.5 says is "next-session
    deliverable" but not yet built.

This document recommends (a) for Phase 1 (smallest, unblocked,
satisfies the backbone-topology question). User should confirm or
direct otherwise.

**Decision 2: Does `pdb_atom_name` stay on `Atom` post-Phase-1?**

Phase 1 of THIS chemistry question doesn't delete `pdb_atom_name`
— other consumers (chi-angle resolver, ring detection, protonation
state, charge-table lookup, PDB emitters) still use it. The
audit's overall plan deletes the field at Phase 7 after all those
consumers migrate. User should confirm: Phase 1 of backbone-topology
is purely additive (add `Atom::backbone_role`, switch `CacheResidue
BackboneIndices` to typed, leave `pdb_atom_name` alone), with the
field deletion deferred until every chemistry-question Phase 1
lands.

**Decision 3: Where does `BackboneRoleFromName` live?**

Three reasonable homes:
- New file `src/AtomNameParser.{h,cpp}` housing only the parser.
  Cleanest separation; small new module.
- Extend `NamingRegistry` with a `BackboneRoleFor(canonical_name)`
  method. Treats it as part of the naming-boundary discipline.
  Conceptually aligned with PATTERNS.md "naming boundary cross
  once": NamingRegistry already collapses `H ↔ HN`, so giving it
  the typed-output API "what role is this name?" is consistent.
- Inline in `src/Atom.{h,cpp}`. Simplest. Probably wrong because
  Atom is an identity-only class (no chemistry methods).

User should decide. The recommendation here is the second:
NamingRegistry is already the typed-translation gatekeeper. Adding
`BackboneRole BackboneRoleFor(const std::string& canonical_name)`
keeps the parsing-and-typing in one place. Also opens the path for
future `LocantFor`, `DiastereotopicIndexFor` if the full identity
surface is needed.

**Decision 4: Glycine HA cache.**

Today `res.HA` for Gly is set to the HA2 atom (the `name == "HA2"`
alias at `Protein.cpp:551`). Per the typed semantics, Gly's HA2
and HA3 carry `BackboneRole::None` (their identity comes from
`Locant::Alpha` + `DiastereotopicIndex`). So a strict `switch
(backbone_role)` DOES NOT MATCH Gly's HA2/HA3 — meaning post-
Phase 1, `res.HA` for Gly would be `Residue::NONE`. That breaks
every calculator that reads `res.HA` for Gly residues today.

Three options:
- (i) Hand-add a `if (res.type == AminoAcid::GLY) { ... }` branch
  in `CacheResidueBackboneIndices` to find HA2 by typed identity
  (Locant::Alpha, DiastereotopicIndex::Position2).
- (ii) Change `BackboneRoleFromName("HA2", AminoAcid::GLY)` to
  return `BackboneRole::AlphaHydrogen` (deviating from the
  substrate generator's convention but matching the existing
  cache shape).
- (iii) Phase-2 the question: `res.HA` for Gly becomes `NONE` in
  Phase 1; downstream consumers learn to check `res.type ==
  AminoAcid::GLY` and read both HA2/HA3 via a separate accessor.
  This breaks downstream code in Phase 1 — not acceptable.

The recommendation is (i): special-case Gly in the cache populator,
with a typed lookup (LocantSetWithDiastereotopic) hidden behind
the special case. Option (ii) muddies the substrate convention;
option (iii) breaks Phase 1's "calculators continue working" rule.

User should bless or redirect.

**Decision 5: CB cache.**

CB is in `Residue.h:55` (with NONE for Gly) but is NOT in the
`BackboneRole` enum (it's a sidechain Cβ, Locant::Beta, Element::C,
no branch). Phase 1's typed `switch (backbone_role)` will not
populate `res.CB`. Three options:
- (i) Keep one defensive `if (atoms_[ai]->pdb_atom_name == "CB")
  res.CB = ai;` line for CB only; deal with it later.
- (ii) Use a typed identity check: `Locant::Beta + Element::C +
  branch{0,0} + di_index == None`. Requires `Atom::locant` field
  too — bigger surface than backbone alone.
- (iii) Move CB out of `Residue` entirely; it's not a backbone
  atom and lives more naturally as a Locant::Beta sidechain
  lookup. Bigger refactor; affects every consumer of `res.CB`.

Recommendation: (i) for Phase 1 (single defensive name match,
flagged for later, isolated). Audit's Phase 2 unifies through the
substrate.

User should bless.

## 8. Risks and unknowns

**Risk 1.** PDB loads that bring atoms with non-canonical CHARMM
names. PdbFileReader does NOT call NamingRegistry (per
`src/PdbFileReader.h:12`), so a PDB authored under CHARMM
conventions could land with `HN` in `Atom::pdb_atom_name`. If the
runtime `BackboneRoleFromName` accepts ONLY canonical AMBER/IUPAC
names (`H` not `HN`), the H cache for that load is empty.
Mitigation: `BackboneRoleFromName` accepts both `H` and `HN`,
collapsing the alias internally — this puts the alias rule in ONE
place (the parser), not in four (the cache populator + three
loaders). Same for `HA` and `HA2` (Gly). The H1 alias from
OrcaRunLoader is correctly NOT accepted (H1 is a cap atom, not the
backbone amide H — fixing the latent bug in current OrcaRunLoader).

**Risk 2.** The redundant cache populators in load sites are
"defensive": maybe they exist because someone debugging once saw
`Protein::CacheResidueBackboneIndices` not get called on a path,
or saw it overwritten downstream. Removing them in step 5 of Phase
1 might surface a load path that doesn't call FinalizeConstruction.
Mitigation: audit-grep for `protein->FinalizeConstruction` /
`protein->AddAtom` pairs to confirm every loader calls Finalize.
If a loader doesn't, fix the loader, don't keep the redundant
cache.

**Risk 3.** `BackboneRole` doesn't cover CB. Today's cache code
treats CB as a special seventh slot; Phase 1's typed switch
silently drops CB. The Decision 5 question above is real and must
not be skipped — silently leaving `res.CB == NONE` post-pivot is a
calculator regression.

**Risk 4.** Substrate population (the §H.5 work) might land before
Phase 1 of THIS question. If so, Phase 1 should be re-evaluated:
the typed `Atom::semantic.backbone_role` would already exist and
`Atom::backbone_role` becomes a redundant field. The cleaner path
is to switch the cache populator to read
`atoms_[ai]->semantic.backbone_role` directly (or whichever field
shape lands). Mitigation: this Phase 1 plan is small enough (~25
net lines) that it can be done EITHER as standalone (now) OR
folded into the substrate-population landing. User should confirm
sequencing.

**Risk 5.** Tests asserting `pdb_atom_name == "..."` on backbone
atoms (some of the Category F sites in the audit) will continue to
work in Phase 1 (the field still exists). When `pdb_atom_name`
deletes (audit Phase 7), test fixtures need to migrate — but
that's not a Phase 1 concern.

**Unknown 1.** The `OrcaRunLoader.cpp:248` `H1` alias — is it
load-bearing? Does some ORCA-output flow actually emit `H1` for
the backbone amide H? If yes, the runtime `BackboneRoleFromName`
must accept it; if no (it's just paranoia), it should not. Quick
check: scan tests/data/ for ORCA outputs and verify atom-name
conventions. Probable answer: H1 is the N-terminal cap H, and the
current OrcaRunLoader is silently mis-classifying it. The Phase 1
typed lookup correctly drops it. But this is a behaviour change
worth confirming on the existing fleet.

**Unknown 2.** Are there other load paths the audit doesn't list?
Check `src/protonation/` and `src/PropkaProtonator.cpp` /
`src/KamlProtonator.cpp` — they read PDB they wrote, so cycle
back into PdbFileReader; no new load surface. The only loaders
are the four named.

**Unknown 3.** Does the chi-angle resolver in
`CacheResidueBackboneIndices` (lines 559-569) need to move with
the backbone cache? It uses `pdb_atom_name` against
`ChiAngleDef::atoms[]`. Out of scope for THIS chemistry question
but flagged: anyone touching `CacheResidueBackboneIndices` will
see both. Audit Hotspot 2 covers chi separately. Recommendation:
leave chi alone in this Phase 1; it's a separate atomic refactor.
