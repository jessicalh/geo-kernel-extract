# Chemistry Question 4 — External-model atom typing (2026-05-05)

## Summary

`AIMNet2Result.cpp` and `ApbsFieldResult.cpp` are the two
calculators flagged in the audit as touching `pdb_atom_name`
(Part 2.2). Reading them shows the diagnosis is narrower and
cleaner than the audit's worst-case framing:

- **AIMNet2Result has exactly one `pdb_atom_name` read** — at
  `src/AIMNet2Result.cpp:169`, inside an error message that fires
  when the protein contains an `Element::Unknown` atom. The
  calculator's *physics* surface — the model input tensor — uses
  only `Element` (lines 224-231) and Cartesian positions; the
  neural net does not see atom names. **AIMNet2 needs no
  string-keyed lookup at all** (Category D in the audit:
  diagnostic-only).
- **ApbsFieldResult has exactly one `pdb_atom_name` read** — at
  `src/ApbsFieldResult.cpp:140`, inside an error message that
  fires when an atom's `pb_radius` is missing. The calculator's
  *physics* surface — the C bridge call to APBS — consumes only
  `partial_charge` and `pb_radius` from `ConformationAtom` (read
  at lines 134, 136). The string-keyed projection happens
  *upstream* of this calculator, in `ChargeAssignmentResult`'s
  ChargeSource implementations (`ParamFileChargeSource::LoadCharges`,
  `AmberPreparedChargeSource::LoadCharges`); APBS itself reads
  pre-projected typed numbers (Category D: diagnostic-only).

So the audit's Category-C "identity-proxy hack" framing is wrong
for both of these calculators. The architecture is already
substantially correct: **physics through typed numerical fields,
strings only at error-rendering time**. The string-keyed boundary
exists, but it is enforced one architectural layer up
(ChargeSource), not at the APBS calculator surface. This
document explains where the projection happens, why that's the
right shape, and what should change for Phase 1 + Phase 2.

The two cases are different in *shape* but identical in
*conclusion for the calculator*: drop both `pdb_atom_name` reads,
substitute typed-identity rendering for the error-message use,
let the upstream ChargeSource layer carry the wire-format
projection responsibility for APBS's input data.

---

## 1. The chemistry question

External models and external parameter tables have their own
atom-typing schemes. AIMNet2 takes elements + topology + positions
and outputs predicted charges. APBS takes positions + charges +
radii and outputs a potential grid. PARSE-format radii tables are
keyed by `(residue_name, atom_name)` strings. AMBER ff14SB charge
tables are keyed by `(terminal_state, residue_name, atom_name)`
strings. PRMTOP files store atom-name strings as the chain-residue
inventory.

The audit (`mechanical-identity-model-and-audit-2026-05-05.md`
Part 2.2) classified the two calculators here as Category C —
identity-proxy hacks where calculator code branches on
`pdb_atom_name` to make a chemistry decision. Reading the code
shows that classification is wrong for both: there is no chemistry
decision flowing through `pdb_atom_name` in either calculator's
physics path. Both reads are render-only diagnostic uses inside
error messages.

The *real* chemistry question for these two calculators is
slightly different than the audit framed:

> Where does the wire-format string projection (typed Atom →
> AMBER atom-name string for ff14SB lookup; typed Atom →
> AMBER atom-name string for PRMTOP match) live, and is that
> location architecturally sound?

The post-audit answer is yes: the projection happens in the
ChargeSource implementations (`ParamFileChargeSource`,
`AmberPreparedChargeSource`), and the calculator surface
(`ApbsFieldResult::Compute`) sees only typed numerical results
populated on `ConformationAtom` by `ChargeAssignmentResult` (the
adapter that mounts a `ForceFieldChargeTable` onto per-atom
`partial_charge` + `pb_radius`). This is the **Crystal Projection
Rule** in action: typed identity throughout the runtime, the
projection-to-string only at the wire boundary itself, with the
calculator consuming the typed result the boundary produced.

The two `pdb_atom_name` reads at the calculator surface are
Category D (diagnostic) leftovers, not chemistry-bearing
substrate violations. Phase 1 swaps them for typed-identity
renderers; Phase 2 may or may not move the wire-boundary
projection point itself.

---

## 2. AIMNet2Result analysis

### 2.1 What the neural net needs

AIMNet2 is a graph neural network: per-atom embeddings keyed by
atomic number, message passing on a neighbour graph (short-range
+ long-range cutoffs), and a charge-readout head. The input
contract is set by the `.jpt` TorchScript model:

| Input key | Type | Shape | Source |
|---|---|---|---|
| `coord` | float32 | (N+1, 3) | atom Cartesian positions |
| `numbers` | int64 | (N+1,) | atomic numbers Z |
| `charge` | float32 | (1,) | net molecular charge (zero) |
| `mol_idx` | int64 | (N+1,) | molecule grouping (single mol = zeros) |
| `nbmat` | int32 | (N+1, max_nb) | short-range half-neighbour list |
| `nbmat_lr` | int32 | (N+1, max_nb_lr) | long-range half-neighbour list |
| `cutoff_lr` | float32 | (1,) | LR cutoff radius |

Atom names do not appear in this contract. AIMNet2's per-atom
representation is **(Z, position, neighbours)** — element +
geometry only. Methylene-vs-methyl distinction, prochirality,
diastereotopic indices, IUPAC locant — none of these ground-truth
atom-typing concepts is visible to the network. They are
re-discovered (or not) from the geometry by the network's own
features.

This is consistent with the AIMNet2 paper's design: it is meant
to generalise across chemistry without a pre-curated atom typing
scheme. The model carries its own atom-typing implicitly in the
embedding tables and message-passing parameters.

### 2.2 What the current code does

`grep "pdb_atom_name" src/AIMNet2Result.cpp` returns one site:

- **`src/AIMNet2Result.cpp:169`** — inside the
  `Element::Unknown` guard that fires before any neural-net
  input is built. The error message reads
  `"Atom " + std::to_string(i) + " (" + protein.AtomAt(i).pdb_atom_name + " in residue " + std::to_string(protein.AtomAt(i).residue_index) + ") has Element::Unknown..."`.
  The string is consumed by `OperationLog::Error` and the
  function returns `nullptr`. The string never enters the model
  input tensor; it never participates in any chemistry decision.

The rest of the calculator's chemistry path:

- Lines 224-231: build the `numbers` tensor by **switching on
  `protein.AtomAt(i).element`** (the typed `Element` enum) into
  atomic-number integers. No atom name involved.
- Line 286: store output charges by `conf.MutableAtomAt(i).aimnet2_charge = ...`.
  The `aimnet2_charge` field on `ConformationAtom` is the typed
  result.
- Lines 314-318: store the AIM embedding into `aimnet2_aim`.
- Lines 348-441 (`ComputeCoulombEFG`): downstream EFG
  computation reads `conf.AtomAt(j).aimnet2_charge` (typed) and
  classifies sources via the **typed backbone cache**
  (`res.N`, `res.CA`, `res.C`, `res.O`, `res.H`, `res.HA`,
  `res.CB`, lines 377-378) and the typed ring inventory
  (`protein.RingAt(ri).atom_indices`, lines 381-383). These
  reads are typed; no `pdb_atom_name` involvement. (The backbone
  cache itself is populated by string match in
  `Protein::CacheResidueBackboneIndices` per the audit's
  Hotspot 2 — that's an upstream issue, not an AIMNet2-side
  issue.)

So AIMNet2's Category-C audit classification is empirically
unsupported. There is no identity-proxy hack here; there is one
diagnostic string in an error message.

### 2.3 Phase 1 design

Replace the single error-message use with a typed-identity
renderer:

```cpp
// Today (line 167-172):
"Atom " + std::to_string(i) + " (" +
 protein.AtomAt(i).pdb_atom_name + " in residue " +
 std::to_string(protein.AtomAt(i).residue_index) +
 ") has Element::Unknown..."

// Phase 1 target:
"Atom " + std::to_string(i) + " (" +
 RenderAtomIdentity(protein.AtomAt(i).identity, ...) +
 " in residue " + ...
 ") has Element::Unknown..."
```

The shape of `RenderAtomIdentity` is a project-wide decision (see
the audit's Closing Note 7) — it could be a free function, a
method on `Atom`, or a method on `Protein`. The point for
AIMNet2's Phase 1 is the same regardless: nothing about the
calculator's *physics* surface needs to change. The neural net
already consumes `Element` + position + neighbours.

The `Element::Unknown` guard itself stays — AIMNet2 has no
embedding for Z=0, and silently producing zero charges for those
atoms would corrupt downstream calculations. The guard is correct
physics. Only the message rendering needs to switch from string
field to typed identity.

**Flag for Phase 2.** AIMNet2 is a candidate for replacement.
Future neural-net charge models (Schnet variants, EquiformerV2,
DimeNet++, etc.) follow the same input contract: element +
geometry + neighbours, no atom names. So the Phase 1 design — no
string at the calculator surface — is also the Phase 2 design.
The calculator surface is stable across model swaps; only the
TorchScript loading path and the per-model attribute reads
(`module.attr("cutoff")`) would change.

---

## 3. ApbsFieldResult analysis

### 3.1 What APBS needs

APBS is a Poisson-Boltzmann solver. Its input is purely numerical:

| Input | Type | Source |
|---|---|---|
| atom positions (x, y, z) | doubles | `conf.PositionAt(i)` |
| partial charge per atom | double | `conf.AtomAt(i).partial_charge` |
| Poisson-Boltzmann radius per atom | double | `conf.AtomAt(i).pb_radius` |
| dielectric constants (pdie, sdie) | scalars | hardcoded (4.0, 78.54) |
| temperature, ionic strength | scalars | hardcoded (298.15 K, 0.15 M) |
| grid sizing (fine/coarse extents, dim) | scalars | computed from atom bounding box |

This is a pure-data contract. APBS does not consume atom names,
residue names, or any string-keyed identity. The C bridge
`apbs_solve` (declared in `apbs_bridge.h`) takes raw double
arrays. The output is a potential grid; the calculator extracts
E-field and EFG via central differences on the grid.

### 3.2 What the current code does

`grep "pdb_atom_name" src/ApbsFieldResult.cpp` returns one site:

- **`src/ApbsFieldResult.cpp:140`** — inside the `pb_radius`
  validation loop, error message when `r <= 0.0`:
  `"missing PB radius for atom " + std::to_string(i) + " (" + protein.AtomAt(i).pdb_atom_name + ")"`.
  String is consumed by `OperationLog::Error` and the
  function returns `false`. No model decision flows through it.

The rest of the calculator's path (lines 126-144, 222-260):

- Lines 130-133: read positions from `conf.PositionAt(i)` — typed
  `Vec3`.
- Line 134: read `conf.AtomAt(i).partial_charge` — typed `double`.
- Line 136: read `conf.AtomAt(i).pb_radius` — typed `double`.
- Lines 184-194: pass raw `double*` arrays to the C bridge.
- Lines 222-260: extract E-field and EFG by grid interpolation;
  store into `conf.MutableAtomAt(i).apbs_efield` (Vec3),
  `apbs_efg` (Mat3), `apbs_efg_spherical` (SphericalTensor). All
  typed. No strings.

So `ApbsFieldResult` itself is already typed end-to-end on the
physics path. Its single string read is a diagnostic. Same shape
as AIMNet2.

### 3.3 Where the radius lookup goes

This is the architecturally interesting part. APBS needs
`pb_radius` per atom; that radius comes from somewhere; the
"somewhere" is the wire-boundary projection.

Tracing upstream from `ApbsFieldResult::Compute`:

1. **`ApbsFieldResult` declares `Dependencies()` →
   `ChargeAssignmentResult`** (line 39 of `ApbsFieldResult.h`).
2. **`ChargeAssignmentResult::Compute`** (`src/ChargeAssignmentResult.cpp`
   lines 51-52): mounts `ForceFieldChargeTable` onto
   `ConformationAtom` by reading `table.PartialChargeAt(ai)` and
   `table.PbRadiusAt(ai)` and writing them onto
   `ca.partial_charge` and `ca.pb_radius`. **Pure numerical
   transfer; no string involvement.**
3. **`ForceFieldChargeTable`** is built from a `ChargeSource`
   implementation. The string-to-typed projection happens here,
   inside the source.

The `ChargeSource` hierarchy (`src/ChargeSource.h`):

- **`ParamFileChargeSource`** (ff14SB flat-table path).
  `LoadCharges` (`src/ChargeSource.cpp` lines 174-259) is the
  string-keyed projection site. For each protein atom, it reads
  `identity.pdb_atom_name`, builds candidate AMBER names via
  `AtomNameCandidates(...)` (handles N-terminal `H` ↔ `H1`
  alias), calls `Ff14sbVariantResidueName(...)` to project the
  residue's typed `(AminoAcid, protonation_variant_index)` into
  an AMBER residue-name string, and looks up
  `(terminal_token, ff_resname, atom_name)` in the parsed flat
  file. Returns numerical `AtomChargeRadius{partial_charge,
  pb_radius}`. **THIS is where the wire-boundary projection lives
  for the ff14SB case.** It is at the lookup site, not at the
  calculator surface.

- **`PrmtopChargeSource`** (upstream PRMTOP from `--orca` /
  `--mutant`). `LoadCharges` (`src/ChargeSource.cpp` lines 328-375)
  reads the PRMTOP `%FLAG CHARGE` and `%FLAG RADII` sections by
  byte position — purely numerical. It assumes
  PRMTOP-atom-order matches Protein-atom-order (justified
  upstream by `OrcaRunLoader` which constructs the Protein from
  the PRMTOP atom inventory). **No atom-name string projection
  in this source** — the upstream loader's atom-order invariant
  carries the identity coupling, no wire-boundary string lookup
  per atom.

- **`AmberPreparedChargeSource`** (PRMTOP regenerated by
  in-process tleap when the flat ff14SB table can't represent
  the protein). `LoadCharges`
  (`src/AmberPreparedChargeSource.cpp` lines 259-441) DOES use
  `pdb_atom_name` (lines 376-382, 427) to build a per-residue
  `(name → atom_index)` map for cross-walking PRMTOP atom names
  against extractor atom names. **This IS a wire-boundary string
  projection — but it's at the source, not at APBS.**

So for the radius lookup question:

- The string-keyed lookup happens inside the ChargeSource
  implementation (whichever source is active for this Protein's
  build path).
- The result of that lookup is a typed pair `(partial_charge,
  pb_radius)` per atom.
- That pair is mounted onto `ConformationAtom` by
  `ChargeAssignmentResult`.
- `ApbsFieldResult` reads the typed pair off `ConformationAtom`
  and feeds it to the C bridge as `double*`.

The architecture is **already substantially correct** by the
Crystal Projection Rule's standard. The string-keyed projection
exists — it has to, the flat ff14SB file and the PRMTOP atom-name
column are wire formats — but it lives at the wire boundary
(inside the ChargeSource), not at the APBS calculator surface.

The audit's Category-C framing of `ApbsFieldResult` was based on
the existence of one `pdb_atom_name` read in error message text.
That read is Category D (diagnostic-only). The architectural
question — "where does the wire-boundary projection happen?" —
has the right answer already: in the source.

### 3.4 Phase 1 design

For `ApbsFieldResult` itself: replace the one error-message read
with `RenderAtomIdentity(...)`, same shape as AIMNet2. The
calculator's physics path is already typed end-to-end; nothing
else needs to change.

For the upstream string projection (the actual chemistry-bearing
boundary): this is inherited from the audit's Phase 5 ("Migrate
external-table lookups (Category B read-side)"). The three
projection sites are:

- `AmberChargeResolver.cpp:306,321` — verdict computation.
- `ChargeSource.cpp:222,243` — `ParamFileChargeSource::LoadCharges`.
- `AmberPreparedChargeSource.cpp:376,378,382,427` — PRMTOP atom-name match.

Phase 5's Phase-1 design (per audit): keep the projection at the
lookup site, but project from typed identity to AMBER atom-name
string at that site, instead of reading `Atom.pdb_atom_name`
directly. The function signature is something like:

```cpp
std::string AmberAtomNameFor(const AtomMechanicalIdentity& identity,
                              AminoAcid residue_type,
                              ResidueTerminalState terminal_state);
```

This function lives next to the wire-format readers/writers
(`ChargeSource.cpp`, `AmberLeapInput.cpp`, `AmberPreparedChargeSource.cpp`,
`AmberChargeResolver.cpp`, the Protonator PDB writers). It is
the Crystal Projection Rule's function form: pure projection from
typed substrate fields to AMBER's wire convention, never cached.

For `ApbsFieldResult` specifically, the Phase 1 deliverable is
trivial:

1. Replace `protein.AtomAt(i).pdb_atom_name` (line 140) with the
   project-wide error-rendering helper that takes typed identity.
2. No other changes; the C bridge call and the per-atom radius
   read stay.

---

## 4. Crystal Projection Rule, applied

The Crystal Projection Rule (`spec/plan/amber-implementation-plan-2026-04-29.md`
§ "Crystal projection rule") states:

> A projection in this design is a PURE FUNCTION on the typed
> substrate fields. It is never a stored string. Every call
> recomputes the projected name from the typed fields.
>
> Output writers (H5 emitters, methods text generators, BMRB
> joiners) call the projection at write time and consume the
> returned string locally. Nothing stores it on the topology side.

For the external-model-typing chemistry question, the rule maps as
follows:

**Wire boundary definition.** A wire boundary is any interface
where bytes leave or enter the process under an external schema
the project does not control. For the two calculators here:

| Calculator | Wire boundary | Schema | Strings live here? |
|---|---|---|---|
| AIMNet2Result | `module.forward({input_dict})` to TorchScript | numerical-tensor input dict (Element/positions/neighbours) | NO — pure numerical |
| ApbsFieldResult | `apbs_solve(...)` C bridge | `double*` arrays (positions, charges, radii) | NO — pure numerical |

Neither calculator's *direct* wire boundary uses strings. The
indirect wire boundary for APBS is the ff14SB flat file, the
PRMTOP file, or the prepared-PRMTOP atom-name column. Those wire
boundaries DO use strings — but they live one architectural layer
upstream, in the ChargeSource implementations, where projection-to-string
is enforced at the lookup site:

```
[Atom typed identity]
    ↓ ParamFileChargeSource::LoadCharges (projects typed → AMBER name string)
[ff14SB flat file lookup; numerical AtomChargeRadius per atom]
    ↓ ForceFieldChargeTable::Build (numerical accumulation)
[Per-atom partial_charge, pb_radius numerical]
    ↓ ChargeAssignmentResult::Compute (mounts onto ConformationAtom)
[ConformationAtom.partial_charge, .pb_radius (typed numerical)]
    ↓ ApbsFieldResult::Compute (reads typed numerical, calls C bridge)
[apbs_solve(double*); typed E/EFG output]
```

The string-keyed lookup is bracketed between two typed surfaces.
The calculator never sees the string. The wire-format reader
projects at its own boundary — which is correct.

**Calculator-side typed; utility-side projection; no strings on
the calculator's input/output surface.** This is what the rule
prescribes. For these two calculators it is already substantially
in place; the only deviation is the two error-message uses of
`pdb_atom_name` at the calculator surface. Those are diagnostic
strings, not chemistry surface, and they get rendered from typed
identity by a project-wide helper in Phase 1.

**Why this matters for the substrate refactor.** When
`Atom.pdb_atom_name` is deleted (audit's Phase 7), the
calculator-side error-message renderer must already exist and
take typed identity. The renderer does not know about ff14SB or
AMBER's quirks (HG vs HG1 for SER, etc.); it produces a
human-readable form like `"HB2 (Cys42)"` or
`"Hβ (CYS-42 prochiral-2)"`. The wire-boundary functions
(`AmberAtomNameFor`, `IupacAtomNameFor`, `BmrbAtomNameFor`,
etc.) are separate — each format-specific, each pure projection
from typed identity, each used at its own wire-format
read/write site.

The clean separation is: error rendering is a *display*
function; wire-format projection is a *protocol* function.
Both project from typed identity; they project to different
strings.

---

## 5. Phase 2 considerations

The Phase 1 design above keeps the calculator surface stable
across future model swaps. Phase 2 considerations:

**AIMNet2 → AIMNet3 / Schnet / Equiformer.** The input contract
(element + position + neighbours) is shared across this entire
class of GNN-based atom-property predictors. A Phase 2 swap would
change the TorchScript loading path and the output extraction
(reading `aim` embedding, etc.), but the calculator's input
construction stays. The `Element` switch at lines 224-231 is the
load-bearing element-typing mapping and is stable. No
`pdb_atom_name` involvement at any stage. Phase 2-safe.

**APBS → DelPhi / Zap / numerical PB solver in libgromacs.**
Different solvers consume different input formats:

- DelPhi: PQR file (atom-name-keyed; atom-name from typed identity).
- Zap: ZAP-format inputs.
- libgromacs PB: would need TPR + parameters.

The C bridge interface (`apbs_solve(positions, charges, radii,
...)`) is solver-specific. A swap would replace the bridge call
but keep the input data shape (positions / charges / radii are
numerical in every PB solver). The calculator side still reads
typed `partial_charge` and `pb_radius` from `ConformationAtom`.
Phase 2-safe.

**Charge-source diversification.** The current ChargeSource
hierarchy already absorbs:

- ff14SB flat table (string-keyed wire format).
- AMBER PRMTOP (positional, atom-order-coupled).
- Pre-loaded charges (caller-supplied; e.g. from libgromacs).
- AMBER prepared PRMTOP (run-time tleap; atom-name-keyed
  cross-walk).

A future ff19SB / OPLS / FEP / RESP source would add another
ChargeSource subclass. Each subclass owns its own wire-boundary
projection. The calculator (APBS) sees only the projected typed
result. Phase 2-safe.

**The radii table specifically.** ff14SB ships PB radii inline
with charges in the same flat file (mbondi2 set, per the file
header). PARSE-format radii (used by some PB workflows) are a
separate file; if the project ever supports a PARSE-radii path
distinct from the ff14SB-radii path, that would be a new
ChargeSource subclass (or a new RadiiSource sibling abstraction).
The Crystal Projection Rule still applies at that subclass's
boundary. Phase 2-safe; the architectural shape extends cleanly.

The general claim: the typed calculator surface is the stable
contract. Wire-boundary subclasses come and go beneath it. The
substrate refactor (typed `AtomMechanicalIdentity` on `Atom`,
typed `AtomSemanticTable` on each atom) makes that contract
explicit. Future model swaps reshape the wire-boundary subclasses
but never propagate up to the calculator.

---

## 6. Choices that need a user decision

The substantive questions for the user before Phase 1 lands:

**Q1. Where does the typed-identity error-rendering helper live?**

Two `pdb_atom_name` reads at the calculator surface (AIMNet2:169,
ApbsField:140) need a typed-identity replacement. Three plausible
homes:

- **(a) Free function in a new header `AtomIdentityRender.h`.**
  Signature `std::string RenderAtomIdentity(const AtomMechanicalIdentity&, AminoAcid residue_type, size_t residue_index, int sequence_number);`.
  Lives next to `SemanticEnums.h`. Pro: clear separation between
  the rendering helper and the typed substrate. Con: adds a
  new header.

- **(b) Method on `Atom`.** Signature
  `std::string Atom::RenderForLog(const Residue&) const;`.
  Pro: object-answers-questions-about-itself convention.
  Con: requires `Residue` reference at call site, slightly more
  ceremony.

- **(c) Method on `Protein`.** Signature
  `std::string Protein::RenderAtomForLog(size_t atom_index) const;`.
  Pro: single argument from any calculator; Protein has all the
  context. Con: `Protein` accumulates methods that aren't really
  about the protein itself.

This choice has project-wide implications because every
diagnostic string that currently contains `atom.pdb_atom_name`
will use this helper. Five Category-D sites in the audit
(`AIMNet2Result.cpp:169`, `ApbsFieldResult.cpp:140`,
`PdbFileReader.h:12` doc, `tests/test_traversal_dump.cpp:6,74`,
plus implicit uses in error returns from charge sources) all
adopt the same helper.

**Q2. Does ChargeSource own a typed-to-AMBER-name projection
function, or does it call a separate utility?**

Today's `ParamFileChargeSource::LoadCharges` reads
`Atom.pdb_atom_name` directly. Phase 5 of the audit says project
at the lookup site. The shape question is:

- **(a) Inline projection.** `ParamFileChargeSource::LoadCharges`
  contains a private helper `AmberAtomNameFromIdentity(...)`
  that turns typed `AtomMechanicalIdentity + Residue` into the
  AMBER atom-name string for flat-table lookup. Pros: locality.
  Cons: every wire-format reader (`AmberPreparedChargeSource`,
  `AmberChargeResolver`, `AmberLeapInput`, the three Protonator
  writers) reimplements the same projection.

- **(b) Shared utility.** A single
  `AmberAtomNameFor(identity, residue_type, terminal_state)`
  function in a new file `AmberAtomNames.{h,cpp}` (or as static
  methods on `LegacyAmberTopology` per the existing architecture
  plan). All wire-format reader/writer sites call it. Pros:
  single source of truth for AMBER's quirks (HG vs HG1 for SER,
  the H1 alias for N-terminal hydrogens, etc.). Cons: a new
  file.

The architecture plan (`openai-5.5-strong-architecture-layout.md`)
already specifies (b) — naming projections as pure functions on
`LegacyAmberTopology` (or in a sibling header). The audit's
Phase 4 + Phase 5 also point at (b). User confirms the existing
plan, or proposes (a) as alternative.

**Q3. Does `AmberPreparedChargeSource`'s PRMTOP cross-walk stay
atom-name-string-keyed or move to typed-identity-keyed?**

Today `AmberPreparedChargeSource::LoadCharges` (lines 370-382)
builds a per-residue `map<string, size_t>` from
`Atom.pdb_atom_name`, then projects PRMTOP atom names against
that map. Both sides are AMBER's wire format. Two options:

- **(a) Stay atom-name-string-keyed.** When `Atom.pdb_atom_name`
  is removed in Phase 7, project from typed identity at the
  extractor side: `AmberAtomNameFor(atom.identity, ...)` produces
  the AMBER name; the resulting string is the map key. Symmetric
  to the PRMTOP side, which already has AMBER-format strings.

- **(b) Move to typed-identity-keyed.** Project the PRMTOP atom
  name into typed identity via `ComputeAtomMechanicalIdentity(...)`
  at parse time; build a `map<AtomMechanicalIdentity, size_t>`
  on the extractor side. The cross-walk happens between two
  typed identities, regardless of what AMBER's name conventions
  do. Pros: typed throughout; the AMBER-name projection function
  is needed only for write paths, not read. Cons: more churn in
  the cross-walk code.

The audit (Phase 5 recommendation) leans toward (b): "build the
map keyed by `AtomMechanicalIdentity` (typed) on the extractor
side." User confirms (b) is the target, or chooses (a) as
expedient transition.

**Q4. ChargeAssignmentResult's role.**

Today `ChargeAssignmentResult` is the *adapter* between
`ForceFieldChargeTable` (per-Protein typed table) and
`ConformationAtom` (per-conformation per-atom fields). It does
not touch strings. After the substrate refactor, this is
unchanged: it stays a pure numerical adapter. Worth confirming
with user that this layer is stable and not part of the audit's
phasing.

---

## 7. Risks and unknowns

**R1. The audit miscategorised these calculators.**

The audit's Part 2.2 listed `AIMNet2Result.cpp:169` and
`ApbsFieldResult.cpp:140` under Category C (identity-proxy
hack) in the count overview, and then in the per-file
classification correctly identified both as Category D
(diagnostic). The count and the per-file analysis disagree.
Reading the code confirms Category D for both: there is no
chemistry decision flowing through `pdb_atom_name` in either
calculator. Phase 1 work is correspondingly small (2 sites,
both error-message-only), not large.

This isn't a defect in the audit — the per-file analysis is
correct and detailed — but the count overview misled. User
should not budget time as if these calculators had identity-proxy
hacks; they don't.

**R2. The wire-boundary projection in ChargeSource is in scope
but understated.**

The real chemistry-bearing string sites for the
external-model-typing question are NOT in the calculators —
they're in the ChargeSource implementations:

- `ParamFileChargeSource::LoadCharges`: 1 site
  (`identity.pdb_atom_name` line 222).
- `AmberPreparedChargeSource::LoadCharges`: 4 sites
  (lines 376, 378, 382, 427).
- `AmberChargeResolver::AnalyzeFlatTableCoverage`: 2 sites
  (lines 306, 321).

These are wire-boundary projections — Category B in the audit's
classification — and the audit's Phase 5 covers them. They are
the substantive work for the external-model-typing chemistry
question. The two calculator-surface sites (Category D) are
trivial.

**R3. PB radii source provenance.**

The current ff14SB flat-table embeds PB radii (mbondi2 set) inline
with charges. The PRMTOP embeds PB radii in `%FLAG RADII`. So
neither path requires a separate radii source. If the project
ever adopts a PARSE-format radii path (where charges come from
one source and radii from another), the architecture will need a
sibling `RadiusSource` abstraction parallel to `ChargeSource`.
This is not on the current roadmap; flag for future awareness.

The PB radius unit chain is undocumented at the calculator
surface — `pb_radius` on `ConformationAtom` carries no unit
comment. Cross-checking against `ChargeSource.h:88-89` shows
"Angstroms, electrostatic/PB radius" but the comment lives on
`AtomChargeRadius`, not on `ConformationAtom`. Worth tightening
the unit-chain documentation when the substrate refactor lands.

**R4. The "non-authoritative PB radii" warning.**

`ApbsFieldResult.cpp:117-123` checks
`charge_result.ChargeTable().NonAuthoritativePbRadiusCount()` and
warns when it's > 0. This relates to the GROMACS/CHARMM TPR path
(now quarantined per `project_charmm_retired_amber_only_2026-05-02`)
that supplied charges but used a placeholder PB radius
(`kCompatibilityPlaceholderPbRadiusAngstrom = 1.5`,
`ChargeSource.h:97`). This is dead code on the AMBER-only
project plan (CHARMM path retired); flag for confirmation that
it can be removed alongside the GmxTprChargeSource path. Not
strictly part of the chemistry-question-4 pivot but adjacent.

**R5. Element::Unknown handling.**

Both calculators guard against `Element::Unknown` — AIMNet2 by
returning nullptr (line 167-174), APBS not directly (it'll
propagate a degenerate result through the C bridge). The
guard-and-render-error pattern is the right shape. The error
message rendering using `pdb_atom_name` is a quality-of-life
issue (humans read the error), not a chemistry issue. The
substrate refactor's renderer must produce something at least
as informative as `pdb_atom_name + " in residue N"`. A typed
identity string like `"H_beta_2 (LYS-42)"` is more informative
than `"HB2 (LYS-42)"` for a human reader; either is acceptable.
Don't degrade the diagnostic.

**R6. Test fixtures.**

`test_object_model.cpp`, `test_amber_prepared_charge_source.cpp`,
`test_amber_charge_resolver.cpp` build test fixtures by writing
`atom->pdb_atom_name = "..."` directly. These are Category F in
the audit (test-code reflection of the surface) and migrate when
the surface migrates. No risk; flag for sequencing — the
calculator-side helper renderer needs to exist before the field
comes off `Atom`, and the test-fixture migration happens
alongside.

---

## Closing

For chemistry question 4, the headline is: **the calculator
surface is already typed; the wire-boundary projection is one
architectural layer up; the audit's two calculator-surface
flagged sites are diagnostic-only; Phase 1 is two trivial
edits**. The substantive substrate work for these calculators
(the wire-boundary projection in ChargeSource implementations)
is captured by Phase 5 of the audit, not Phase 3.

The two questions for the user that gate Phase 1 progress:
- Q1 — error-renderer location (free function vs Atom method vs
  Protein method).
- Q2 — `AmberAtomNameFor` location (inline in ChargeSource vs
  shared utility vs static method on `LegacyAmberTopology`).

Both are project-wide decisions; once made they enable the audit's
Phase 1 (typed identity field on Atom) and Phase 5 (wire-boundary
projection migration) cleanly. Phase 2 is conceptually a no-op
for these calculators — same typed surface, possibly different
external model behind it.
