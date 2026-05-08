# Topology Substrate + Identity Pivot — synthesis and forward plan (2026-05-05)

## Executive summary

In one long working session today (six commits on master, plus a stashed
runtime-integration attempt that we pivoted away from), the project
landed Phase 5 of the chemistry-substrate model end-to-end and then
ran a deliberate audit that surfaced both a structural gap (atom-name
strings on the runtime `Protein` model) and — by cross-cutting
investigation — a serendipity payoff: **the calculators are already
typed-surface-clean**. The string contamination is concentrated at the
load and protein-construction layers, not throughout the calculator
floor. The pivot scope shrinks from the audit's first-cut "20 files,
~600 lines" to roughly 5-7 files of substantive refactor, with several
of the chemistry questions essentially dissolving on contact.

This document is the synthesis. It captures **how we got here**, **what
we discovered**, **what we are thinking about**, and **how to pick up**.

It supersedes neither the evening session handoff
(`session-handoff-20260505-evening.md`) nor the audit
(`mechanical-identity-model-and-audit-2026-05-05.md`) — those remain
authoritative on their respective scopes. This doc is the integrating
narrative + the forward plan.

---

## Part I — How we got here

### 1. The starting point this morning

Phase 5 of the substrate work was the morning's deliverable: define
14 typed semantic fields per atom on `LegacyAmberTopology`; emit them
as a generated constexpr table; expose runtime API for lookup. PHE
was proven end-to-end (commit `721e681`); the architecture was sound
on one residue.

Reading order at start-of-day: `session-handoff-20260505.md` (morning
handoff). Status: substrate architecture LANDED on PHE; 19 standard
residues + 10 protonation variants pending; runtime integration
pending; coverage tests pending.

### 2. The substrate model was built (commits `bb4584d` → `e7e47e6`)

The aim throughout: a **model** in the Stroustrup-replacing-Simula
sense. Typed structure that lets programs say what they mean. Not a
grand edifice, not a religious enterprise — a workable model
where types carry symbolic truth in code and algorithms can do
science rather than power through to a result. The substrate is one
component of that larger model; the runtime objects (`Atom`,
`Residue`, `Protein`, `LegacyAmberTopology`) are the others. Today's
work was building the substrate component cleanly so the rest of
the model can lean on it.

Six commits this evening, in sequence:

```
bb4584d  Phase 5 substrate: residue reference + critical review +
         dependencies + enum extensions
56fb6b9  Phase 5 substrate: encode standard 20 residues; variants pending
ab3bd5e  Phase 5 substrate: emit AMBER variant tables via parent-CCD +
         delta architecture
60f9a11  Session handoff: update for evening Phase 5 substrate completion
3f1d03a  Phase 5 substrate: critique-driven fixes + structural-matching
         architecture spec
ee1f1b4  Phase 5 substrate: structural-matching lookup + cap-table separation
57db9f6  Phase 5 substrate: OpenAI-eval fix bundle (6 findings + §H.10 landed)
e7e47e6  Phase 5 substrate: cap-override field-level rule + runtime API
         surface + stale-comment fixes
```

The work spanned several agent runs under supervision plus two rounds
of OpenAI external evaluation. Key milestones:

- **Reference + review + dependencies**: a deep research pass produced
  the per-(residue, variant, atom) reference doc, an empowered LLM
  review caught two real bugs (Asn/Gln amide H E/Z inversion; LYN HZ
  atom names), the dependencies file captured application-time
  decisions in a way that won't drop.
- **Standard 20 + variants**: encoded chemistry-substrate fields for
  20 residues + 10 AMBER protonation variants. A first agent run
  failed because of a CCD-3-letter-code coincidence trap (CCD's
  `data_HID`, `data_ASH`, etc. are unrelated small molecules); we
  stripped the bogus variant tables and re-encoded via a parent-CCD +
  delta-patch architecture.
- **Structural-matching lookup architecture**: the audit-discovered
  "atom_local_idx mismatch" between CCD atom order and runtime
  `AmberAminoAcidVariantTable` was dissolved by typed `AtomMechanicalIdentity`
  matching + cap-table separation (`§H` of the dependencies file).
- **OpenAI eval rounds**: surfaced six findings (cap deltas, formal
  charges, fail-fast variant_idx, PRO backbone H, identity uniqueness,
  stale comments) — all fixed in one bundle. Then a second round
  found the cap-override-as-whole-row issue, the runtime-API-not-
  exposed issue, and stale comments — all fixed. The substrate is now
  signed off by the SDET (codex) for runtime integration.

**Final substrate state** at HEAD `e7e47e6`:
- 30 typed-enum constexpr `AtomSemanticTable` arrays (20 standard + 10
  variants) in `src/generated/LegacyAmberSemanticTables.cpp`.
- 4 cap tables (`kCapNtermCharged`, `kCapNtermNeutral`, `kCapCtermDeprotonated`,
  `kCapCtermProtonated`) with both *added cap atoms* (H1/H2/H3, OXT,
  HXT) AND *backbone overrides* (N for NTERM, C+O for CTERM) so that
  `ApplyCapDelta` field-level composition produces correct chemistry
  at residue termini.
- Runtime API in `src/generated/LegacyAmberSemanticTables.h`:
  `LookupBy(residue, variant_idx, identity)`, `LookupCap(state, identity)`,
  inline `ApplyCapDelta(chain, cap)`, sentinel `kBaseVariantIdx`.
- The lookup uses **typed structural matching** — `AtomMechanicalIdentity`
  is `(Element, Locant, BranchAddress, DiastereotopicIndex, BackboneRole)`
  — never a string.
- String barrier intact at the linker: `nm libnmr_shielding.a | grep
  -c _ZN5RDKit` returns 0; cifpp confined to `PdbFileReader.cpp` +
  `DsspResult.cpp`; gemmi/RDKit absent.
- ctest 352/352 pass; StringBarrier 5/5; runtime build clean.

The substrate is the model the rest of the project will lean on.

### 3. The runtime integration attempt — and why we pivoted

After the substrate signoff, an agent was empowered to wire
`Protein::FinalizeConstruction` to call `LookupBy` + `ApplyCapDelta`
for every atom. The agent did the wiring; the integration compiled
and tested green; the result *worked* in the narrow sense.

But it worked by reaching for `Atom.pdb_atom_name` (a pre-existing
PDB-name string field on `Atom`) at FinalizeConstruction time, parsing
that string into an `AtomMechanicalIdentity`, and then querying the
substrate. That parse-from-string path is structurally available to
every other piece of code that holds an `Atom*`.

The user's response, the heart of the pivot:

> We absolutely cannot have atom_name on the protein side of the
> barrier. If it exists, calculators will use it: no documentation
> can stop them. We need a solution which is well considered,
> preserves our excellent work today, and which leaves us with a
> *model*, made of types which can become objects or at the least
> have symbolic truth in code. With that string on there every
> calculator is at risk of becoming an epistemic hack.

And the deeper point:

> The issue is not some nicety, it is that you code based on name,
> name is meaning, and some choices are virtually automatic. It is
> our job, together, in planning mode, to make sure we get a model
> where our algorithms can do science, if we are careful and lucky,
> not just power to a result.

The string carries chemistry meaning implicitly — name patterns
encode chemistry (HD1 → HID protonation; CG → ipso ring carbon;
HB2 → prochiral methylene). Code that reads the string and branches
on its shape is doing chemistry inference via name pattern. The fix
is *structural*: make the string structurally inaccessible from the
runtime model, so that even coding-mode reflexes have to ask the
typed substrate.

Action: `git stash` the agent's integration (preserved as `stash@{0}`
if any pieces are reusable); master returns to `e7e47e6` (substrate
complete, no name-pattern-reading integration).

### 4. The mechanical-identity audit

The first analytical step was a deep audit of what
`AtomMechanicalIdentity` IS in the model, and where strings (specifically
`pdb_atom_name`) actually live in the codebase today.

Output: `spec/plan/mechanical-identity-model-and-audit-2026-05-05.md`
(944 lines). Key results:

**Definition (Part 1)**: `AtomMechanicalIdentity` is a typed *lookup
key* into the chemistry-substrate tables and an *equality predicate*
for atom-as-role. It is NOT a typed substitute for `pdb_atom_name`.
Calculators consume the *substrate data* (`AtomSemanticTable` fields
like `polar_h`, `planar_group`, `formal_charge`, `ring_position`) —
not identity itself, except via narrow typed accessors on Residue/Protein
for graph navigation.

The user's correction (which sharpens the audit's framing):

> Mechanical identity comparisons are fair and a clean translation
> layer.

Identity equality, identity-as-lookup-key, identity-as-substrate-table-
locator are all legitimate model concepts. The antipattern is
*chemistry inference encoded in name shape* — whether the shape is a
`std::string` like `"HB2"` or a typed-tuple like `(Beta, {0,0},
Position2)`. In both cases the algorithm reads the spelling of an
atom's label and concludes chemistry from it. The model should keep
chemistry in the substrate's typed semantic fields; calculators ask
the substrate, not the name.

**Audit (Part 2)**: 67 `pdb_atom_name` hits across 24 files. Categorised:
- A. Load-boundary write: 4
- B. Boundary-out string-key projection: 9
- C. Identity-proxy hack (calculator-side chemistry from name): 8
- D. Diagnostic / log / error: 6
- E. Internal model use that should not have been allowed: 5
- F. Test-code reflection: 13

Three top-3 hotspots: `Protein::ResolveProtonationStates`,
`Protein::CacheResidueBackboneIndices`, `Protein::DetectAromaticRings`.

### 5. The four chemistry questions — agents in parallel

Rather than one comprehensive analysis, we split into four independent
chemistry questions and spawned an agent per question. The user's
reasoning:

> We will have a better time working through this by subcontracting
> out the apprehension of the problem context and the underlying
> physics and chemistry question if we tell one agent for each of
> our four uses where to look for the actual chemistry information
> we have in topology, and what the story is.

The four:

1. **Protonation state determination** (`Protein::ResolveProtonationStates`)
2. **Backbone topology indexing** (`Protein::CacheResidueBackboneIndices` + `DsspResult.cpp`)
3. **Ring topology** (`Protein::DetectAromaticRings`)
4. **External-model atom typing** (`AIMNet2Result.cpp` + `ApbsFieldResult.cpp`)

Each agent: read project foundation + substrate work + audit + the
specific call sites; articulate the chemistry question rigorously;
identify what typed information our model already carries for that
question; propose Phase 1 design (existing calculators must work)
and Phase 2 considerations (the future walk).

Outputs at `spec/plan/chemistry-question-{1,2,3,4}-...-2026-05-05.md`.

---

## Part II — What we discovered

### 6. The substrate model works

Substrate work at `e7e47e6` is structurally sound:
- Standard 20 + 10 variants emit correctly.
- 4 cap tables carry both added atoms AND backbone overrides for the
  cap-on-chain composition rule.
- `ApplyCapDelta` is the canonical composition primitive; whole-row
  override is structurally prevented by the helper's enumeration of
  cap-controlled vs chain-controlled fields.
- The runtime API is exposed via the hand-written
  `LegacyAmberSemanticTables.h`; the cpp compiles into
  `libnmr_shielding.a` (verified via `nm`); the string barrier holds.
- Two rounds of external review found and fixed all the issues that
  were findable in the substrate model at this scope.

The substrate is *not yet wired into runtime*. That's the integration
work this synthesis is the prelude to.

### 7. Calculators are already clean (the serendipity payoff)

The four chemistry questions revealed that **the calculator floor of
the model is already typed-surface-clean**. Specifically:

**AIMNet2** (the neural-net charge calculator): the audit's first-pass
classified it Category C (identity-proxy hack). On deep read by Q4's
agent, the classification was wrong. AIMNet2's input contract is
`(Element, Cartesian position, bond topology)` — atom names are not
in the contract and never could be. There is exactly one `pdb_atom_name`
read in `AIMNet2Result.cpp` (line 169) — in a `Element::Unknown` error
message (Category D, diagnostic only).

**APBS** (the Poisson-Boltzmann electrostatic-field calculator):
similarly miscategorised. APBS's input is `(x, y, z, radius, charge)` —
typed numerical, not labels. There is exactly one `pdb_atom_name` read
in `ApbsFieldResult.cpp` (line 140) — in a "missing PB radius" error
message (Category D, diagnostic only).

**Ring-current calculators** (BiotSavart, HaighMallion, RingSusceptibility,
Dispersion, plus one more): consume `Ring::atom_indices` (typed) and
`Ring::Intensity()` (typed). They never see atom names. The Ring object
itself has a typed API; only the populator (`DetectAromaticRings`) uses
strings.

**DSSP**: the agent for Q2 found that all three `pdb_atom_name` reads
in `DsspResult.cpp` are PDB-emit wire-format formatting (Category B,
boundary-out projection). DSSP's actual hydrogen-bond donor/acceptor
chemistry happens inside libdssp on cifpp's structure; results map
back via `(chain, sequence_number)`, not atom names. The DSSP
calculator surface is typed.

The pattern: **chemistry-bearing strings are concentrated at LOAD and
PROTEIN-CONSTRUCTION sites, not in calculators**. The model's
load-time rooms are messy; its calculator floors are well-floored.

### 8. Where the chemistry-bearing strings actually concentrate

The Category C + Category E sites — the ones that should not exist
post-pivot — number 13. They cluster into four chemistry areas:

1. **Protonation-state inference** in `Protein::ResolveProtonationStates`
   — three sites encode the same chemistry-from-name-presence pattern
   (the audit named one; Q1's agent found that
   `Protein::DetectAromaticRings` and
   `FullSystemReader::DetectHisVariantFromAtoms` repeat the same
   HD1/HE2 inference logic).

2. **Backbone-atom indexing** in `Protein::CacheResidueBackboneIndices`
   — the H||HN, HA||HA2 alias-disjunction pattern duplicated FOUR
   times (Protein.cpp, FullSystemReader.cpp, OrcaRunLoader.cpp,
   quarantined GromacsEnsembleLoader.cpp). Q2's agent found this
   is the audit's second hotspot underplayed: NamingRegistry exists
   and was supposed to canonicalise H/HN, but PdbFileReader bypasses
   it entirely — so each consumer defensively re-aliases.

3. **Ring-atom enumeration** in `Protein::DetectAromaticRings` —
   string-match against `AminoAcidType::rings[].atom_names` (`const
   char*` literals like `"CG"`, `"CD1"`, plus a duplicate fallback
   for HIS-tautomer). The substrate already carries `RingPosition`
   per atom but the populator doesn't consume it.

4. **Diagnostic error messages** in the two Category-D calculator
   reads (AIMNet2, APBS). Trivial.

### 9. The substrate already has typed answers for most chemistry questions

Each of the four chemistry questions, the substrate or the typed
project model already carries the answer. The gap is *not* missing
substrate data; the gap is **calculators not consuming it because
the population path isn't wired and the runtime parser doesn't
exist in the runtime library**.

Concretely:
- `BackboneRole` enum exists; `AtomSemanticTable.backbone_role`
  exists; the 6-slot backbone (Nitrogen / AlphaCarbon /
  CarbonylCarbon / CarbonylOxygen / AmideHydrogen / AlphaHydrogen)
  is fully encoded. But: nobody computes `backbone_role` on a
  runtime `Atom` because `ParseAtomName` lives only in the
  build-time generator binary.
- `RingPosition` (RingSystemKind + RingPositionLabel + structural
  flags) is populated in every emitted `AtomSemanticTable`. But:
  no production calculator currently consumes it; ring construction
  uses string match.
- `protonation_variant_index` (typed) lives on `Residue`. Three of
  four load paths populate it correctly upstream. PROPKA and KaML
  protonators produce typed `ResidueProtonation::variant_index`.
  But: PdbFileReader doesn't wire them in — there's a `TODO` at
  `PdbFileReader.cpp:245`.

The model's chemistry rooms are populated; the corridors that
let runtime code walk them are not yet built.

### 10. The four chemistry questions, summarised

#### Q1 — Protonation state determination

Three sites encode the same chemistry-from-H-presence pattern:
`Protein::ResolveProtonationStates`, `Protein::DetectAromaticRings`
HIS-tautomer fallback, `FullSystemReader::DetectHisVariantFromAtoms`.

Three of four load paths already produce typed `protonation_variant_index`:
- `FullSystemReader` via `GromacsToAmberReadbackBlock` (landed
  2026-05-02).
- `OrcaRunLoader` via residue-label-string-match against the typed
  variant table.
- `--protonated-pdb` similarly.

The fourth — `PdbFileReader::BuildFromPdb` (the `reduce`-based path) —
leaves `variant_index = -1`, forcing `ResolveProtonationStates` to do
the inference work. This is the single point the inference logic
exists for.

`PropkaProtonator` and `KamlProtonator` already produce typed
`ResidueProtonation::variant_index` from their pKa predictions
(`PropkaProtonator.cpp:174-179` for HIS). They are **built but
unused** — `TODO` at `PdbFileReader.cpp:245`. Wiring them eliminates
the inference logic at its source.

The sweet path: wire PROPKA/KaML at PdbFileReader, then
`ResolveProtonationStates` shrinks to ~10 lines (CYX-from-disulfide
branch + per-titratable-residue diagnostic), the
`DetectAromaticRings` HIS fallback dissolves (substrate-driven ring
construction makes it irrelevant), and `DetectHisVariantFromAtoms`
becomes redundant or moves entirely behind the GromacsToAmberReadbackBlock
authority.

#### Q2 — Backbone topology indexing

The H||HN, HA||HA2 alias-disjunction is duplicated four times across
the project. NamingRegistry's H↔HN canonicalisation
(`NamingRegistry.cpp:138-139`) is *either not reaching all load paths
or being defensively worked around*. PdbFileReader (line 133) writes
the cifpp-provided name raw — never calls NamingRegistry — so PDB
loads with CHARMM-style HN bypass the registry entirely. Each
consumer that needs backbone indexing has independently re-aliased.

OrcaRunLoader has a latent bug: `H||HN||H1` conflates N-terminal cap
H1 (which is in `kCapNtermCharged`, chemically distinct) with backbone
amide H. Different chemistry; the typed lookup silently fixes this.

Critical finding: **`ComputeAtomMechanicalIdentity` does NOT exist in
the runtime library.** The parser (`ParseAtomName`) lives only in
`tools/topology/build_semantic_tables.cpp:511-535` — the build-time
generator binary, which doesn't link to the runtime. This is the
chicken-and-egg: typed identity is fully designed but no runtime
code computes it. Phase 1 must lift the parser (or at minimum the
`BackboneRoleFromName` switch) into the runtime library — a single
source of truth between generator and runtime.

DSSP's three `pdb_atom_name` reads are PDB-emit wire-format formatting
(Category B), NOT calculator-side navigation. DSSP migrates with
Phase 4 PDB-emitter Crystal-Projection-Rule work, decoupled from
this chemistry question.

#### Q3 — Ring topology

The `Ring` object's calculator-facing API (`atom_indices`,
`Intensity()`, `parent_residue_index`, `fused_partner_index`) is
already typed. Five calculators consume it correctly. **The string
surface is purely the populator** (`DetectAromaticRings`), not part
of the calculator API.

Recommended path: substrate-driven ring construction. Replace
`DetectAromaticRings` with `ConstructRingsFromSubstrate` — walk
atoms with `ring_position.primary.ring != NotInRing`; group by
`(residue_index, RingSystemKind)`; build `Ring` objects from the
typed groups. No string matching anywhere. Phase 1 calculator
behaviour preserved bit-identical (atom-ordering convention is
the verification gate).

PRO is currently invisible to ring detection (early return on
`is_aromatic = false`); the substrate carries `Pyrrolidine_Pro`
for it. Phase 1: skip Pro `Ring` materialisation, defer to Phase
2 if a calculator wants saturated-ring Ring objects.

Critical dependency: substrate population at `Protein::FinalizeConstruction`
(currently not wired) must land first; ring refactor follows it.

#### Q4 — External-model atom typing (the dissolution)

The audit's count-overview overcounted. Both AIMNet2 and APBS are
Category D (diagnostic-only). The audit's per-file classification
states this correctly; the count-overview was the artifact.

AIMNet2's neural net contract: Element + topology + position. Atom
names are not in the contract. The single `pdb_atom_name` read is in
an `Element::Unknown` error message.

APBS's contract: per-atom (x, y, z, radius, charge) — typed numerical.
The single `pdb_atom_name` read is in a "missing PB radius" error
message.

The actual chemistry-bearing string-keyed projection lives one
architectural layer up in the `ChargeSource` implementations
(`ParamFileChargeSource::LoadCharges`, `AmberPreparedChargeSource::LoadCharges`,
`AmberChargeResolver::AnalyzeFlatTableCoverage`). These are
wire-boundary-out projections (Crystal Projection Rule already
substantially obeyed): typed throughout the runtime, projection-to-string
at the wire boundary inside the ChargeSource subclass, calculator
surface consumes typed numerical fields.

Phase 1 cost for Q4: two trivial error-message edits.

The substantive work for Q4 is Phase 5 of the audit's pivot scope
(wire-boundary projection migration in ChargeSource subclasses) —
not the Q4 calculator surface. That's a separate cleanup with its
own scope.

---

## Part III — What we're thinking about

### 11. The model definition (carry-forward from the audit)

`AtomMechanicalIdentity` is, in the model:

- A **typed lookup key** for the chemistry-substrate tables. The
  primary use is `LookupBy(residue, variant_idx, identity)` →
  `AtomSemanticTable*`. This is the legitimate purpose; identity is
  the locator, the substrate is the answer.
- A **typed equality predicate** for atom-as-role. Two atoms with
  the same identity are the same role within the residue (modulo
  chemically-equivalent collisions for methyl Hs and cap Hs, which
  are deliberate per the docstring).
- A **clean translation layer**: strings live at the load boundary
  (PDB / AMBER topology / ORCA output), become typed identity
  through the model, and project back to strings only at output
  boundaries (file-format emitters, external parameter tables we
  don't own).

It is NOT:

- A typed substitute for `pdb_atom_name`. Calculators that branch
  on identity field-by-field (`if (identity.locant == Beta &&
  identity.di_index == Position2)`) are doing the same chemistry-
  from-name-shape inference under a typed disguise. The right
  question is "is this atom prochiral?" via `ProchiralStereo`, not
  "is this called HB2?" via locant + di_index inspection.

- A direct calculator surface. Calculators consume the substrate's
  typed semantic fields (`PolarHKind`, `PlanarGroupKind`, `BackboneRole`,
  `RingPosition`, etc.) — not identity. Identity is the lookup
  machinery.

### 12. The Phase 1 / Phase 2 framing

The pivot has two phases, explicitly separated by the user:

**Phase 1 (immediate refactor)**: existing calculators must continue
to work. The string `pdb_atom_name` comes off the runtime model.
Inference logic that reads name-shape moves upstream (to load
boundaries, to typed-input authorities like PROPKA). The runtime
parser is lifted from the generator into the library so `Atom.identity`
or `Atom.backbone_role` can be populated at load.

**Phase 2 (future walk)**: every existing calculator is revisited.
Calculator-by-calculator, we confirm typed surfaces work under
refactor and pursue deeper structural improvements: time-dependent
protonation under MD; ring-pucker via substrate; substrate-aware
shielding stratification; etc.

The Phase 1 budget is "minimum viable to remove the string from the
runtime model AND keep calculators behaving identically." The Phase 2
budget is "everything else, deliberately."

### 13. The Phase 1 plan, sequenced (TENTATIVE)

Six bullets, in the order they need to land. Total realistic footprint
~300-400 lines across 5-7 files. Each bullet below names the file(s),
the function(s) involved, the expected change shape, the verification
gate, and the dependencies between steps. This is the tentative plan;
the order and shapes are subject to refinement when execution begins.

**Cross-cutting verification gates** that apply at every step (every
bullet must pass these before moving to the next):
- `cmake --build /shared/2026Thesis/nmr-shielding/build -j 8` clean.
- `ctest --test-dir build -j 8` at full count green (currently 352).
- `ctest --test-dir build -R StringBarrier` 5/5.
- `nm build/libnmr_shielding.a | grep -c '_ZN5RDKit'` returns 0.
- Each step's specific spot-check (named per bullet).

**Cross-cutting principle** that applies at every step: never
introduce a new code path that reads `pdb_atom_name` outside the
load-boundary site that wrote it. If a step can't avoid it,
escalate to user — that's a model gap, not a tactical detail.

#### 13.1 Lift `ParseAtomName` into the runtime library

**Files**: source location `tools/topology/build_semantic_tables.cpp:511-535`
(the existing parser, build-time only). Target location: either
`src/NamingRegistry.{h,cpp}` (extends the existing string-canonicalisation
utility — Q2's recommendation) or `src/generated/LegacyAmberSemanticTables.h`
(single source between generator and runtime).

**Decision pending**: which home. Q2 recommended NamingRegistry on the
basis of "alias collapse already lives there." LegacyAmberSemanticTables.h
has the merit that the generator and runtime would literally share one
function. User to decide; both are tractable.

**Function signature** (regardless of home):
```cpp
ParsedAtomName ParseAtomName(std::string_view name,
                             std::string_view parent_name);
struct ParsedAtomName {
    Element             element;
    Locant              locant;
    BranchAddress       branch;
    DiastereotopicIndex di_index;
    BackboneRole        backbone_role;
    bool                is_n_terminus;
    bool                is_c_terminus;
};
```

**Verification gate**: byte-identical generator output (the parser's
generator-side behaviour must not change). Run the generator after the
lift; diff the regenerated `LegacyAmberSemanticTables.cpp` against the
checked-in version; expect zero diff.

**Dependencies**: none upstream. Blocks 13.6 (FinalizeConstruction
needs to compute `AtomMechanicalIdentity` at runtime) and 13.4
(substrate-driven ring construction queries identity).

**Estimated size**: 50-80 lines moved + ~20 lines of supporting struct
definitions if the home is NamingRegistry.

#### 13.1a String boundary discipline check

After the lift, verify with grep that `ParseAtomName` is the **only**
runtime entry point that reads atom-name strings outside loaders. If
other paths exist, surface them — they may be additional Category-C
hits the audit missed.

#### 13.2 Wire PROPKA/KaML at `PdbFileReader.cpp:245`

**Files**: `src/PdbFileReader.cpp` (the TODO site), `src/PropkaProtonator.{h,cpp}`,
`src/KamlProtonator.{h,cpp}`, `src/Residue.h` (the
`protonation_variant_index` field), `src/Protein.cpp:325-413`
(`ResolveProtonationStates` — shrinks).

**The TODO at `PdbFileReader.cpp:245`**: this was always going to be
wired; this is when. PROPKA is the install-target environment default;
KaML is the alternative (selected per project config). Both already
produce typed `ResidueProtonation::variant_index` (e.g.
`PropkaProtonator.cpp:174-179`).

**Wiring shape**:
```cpp
// In PdbFileReader::BuildFromPdb, after structure load:
auto protonation = config.protonator->Predict(structure);
for (auto& residue : residues) {
    auto pred = protonation.find(residue.canonical_id());
    if (pred != protonation.end()) {
        residue.protonation_variant_index = pred->variant_index;
    }
    // Else: leaves variant_index = -1 for ResolveProtonationStates'
    // shrunken consistency-check pass to flag (no longer infer).
}
```

**`Protein::ResolveProtonationStates` shrinks to**:
- CYX-from-disulfide branch (genuinely needs bond-graph detection;
  no name reads).
- Per-titratable-residue consistency diagnostic: if
  `variant_index == -1` AND the residue is titratable, log a typed
  diagnostic naming the residue and recommended action (run PROPKA
  / set explicit annotation). Do NOT infer from name patterns.

**Knock-on effect 1**: HIS-tautomer fallback in
`Protein::DetectAromaticRings` (`src/Protein.cpp:484-497`) dissolves
because step 13.4 replaces the function entirely with substrate-driven
construction.

**Knock-on effect 2**: `FullSystemReader::DetectHisVariantFromAtoms`
(`src/FullSystemReader.cpp:74-86`) dissolves: the
GromacsToAmberReadbackBlock (landed 2026-05-02) already produces
typed variants; this is dead inference.

**Verification gate**: corpus regression test. Load `1ubq_protonated.pdb`
+ `fleet_amber/1Z9B` + `fleet_amber/1P9J`; verify each residue's
`protonation_variant_index` matches what `ResolveProtonationStates`
would have inferred (bit-identical) AND that the new path doesn't
read `pdb_atom_name` for chemistry inference (verify with grep on
the new code paths only).

**Dependencies**: none upstream (independent of 13.1; can land in
parallel). Doesn't block any subsequent step but unblocks the
shrinking of three Category-C sites.

**Estimated size**: 30-50 lines added to PdbFileReader; ~80 lines
removed/shrunk from Protein.cpp + FullSystemReader.cpp.

#### 13.3 Centralise H/HN, HA/HA2 alias collapse in `NamingRegistry`

**Files**: `src/NamingRegistry.{h,cpp}` (extend), `src/PdbFileReader.cpp:133`
(start calling the registry), `src/Protein.cpp:542-571`
(`CacheResidueBackboneIndices` — reduce alias-disjunction),
`src/FullSystemReader.cpp:892-895` (TPR loader inline cache —
reduce), `src/OrcaRunLoader.cpp:248-251` (XYZ loader inline cache —
reduce + fix the H1 latent bug).

**Why centralise**: the H||HN, HA||HA2 alias-disjunction is duplicated
four times across loaders (per Q2's audit). Each consumer defensively
re-aliases because PdbFileReader bypasses NamingRegistry entirely.
Centralisation eliminates the duplication and one latent bug.

**API extension** on NamingRegistry (suggested):
```cpp
// Returns canonical AMBER ff14SB name for the given input + residue
// context. Handles H↔HN, HA↔HA2 (Gly), HZ1↔HZ (some legacy), etc.
// Returns nullopt if the input doesn't match any known alias and no
// canonicalisation is needed (i.e. input is already canonical).
std::optional<std::string> CanonicalAtomName(
    std::string_view input,
    AminoAcid residue_type,
    uint8_t variant_idx,
    TerminalState terminal_state);
```

**Wiring change** in `PdbFileReader.cpp:133`:
```cpp
// BEFORE: writes the cifpp-provided name raw.
atom.pdb_atom_name = cifpp_name;
// AFTER: canonicalise at the load boundary.
atom.pdb_atom_name = NamingRegistry::CanonicalAtomName(cifpp_name,
    residue.type, residue.variant_idx, residue.terminal_state)
    .value_or(std::string(cifpp_name));
```

**Reduction at consumers**: each `if (atom.pdb_atom_name == "H" ||
atom.pdb_atom_name == "HN")` becomes `if (atom.pdb_atom_name == "H")`
(canonical) — or, after 13.6 lands, `if (atom.backbone_role ==
BackboneRole::AmideHydrogen)`.

**OrcaRunLoader latent bug fix**: today's code has
`H || HN || H1` which conflates the N-terminal cap H1 with the
backbone amide H. After canonicalisation, H is canonical for backbone
amide; H1 is a cap atom (`kCapNtermCharged`). The bug fix is
emergent — the centralised registry simply doesn't alias H1 to H.

**Verification gate**: a discipline test that asserts NamingRegistry
is called at every loader's atom-name write site. New gtest:
`tests/test_naming_registry_discipline.cpp` (or extend an existing
test). The test greps loader source files for `pdb_atom_name = `
assignments and asserts each is preceded by a NamingRegistry call.

**Dependencies**: independent of 13.1 and 13.2; can land in parallel.
Helps 13.4 + 13.6 by reducing the alias-noise surface.

**Estimated size**: 30-50 lines added to NamingRegistry; ~5-10 lines
modified at each loader (4 sites = ~20-40 lines net change); ~40
lines reduced from the four duplicated alias-disjunction code paths.

#### 13.4 Replace `DetectAromaticRings` with substrate-driven `ConstructRingsFromSubstrate`

**Files**: `src/Protein.cpp:431-529` (`DetectAromaticRings` — replaces),
`src/Ring.{h,cpp}` (the calculator-facing API stays as-is),
`src/AminoAcidType.h:60`/`AminoAcidType.cpp:92,147,186-188,201`
(`rings[].atom_names` — likely removed).

**Algorithm shape**:
```cpp
std::vector<Ring> ConstructRingsFromSubstrate(const Protein& p) {
    std::vector<Ring> rings;
    std::map<std::pair<size_t, RingSystemKind>, std::vector<size_t>> groups;
    for (size_t a = 0; a < p.atoms_.size(); ++a) {
        const auto& sem = p.legacy_amber_.SemanticAt(a);
        if (sem.ring_position.primary.ring != RingSystemKind::NotInRing) {
            groups[{p.atoms_[a].residue_index,
                    sem.ring_position.primary.ring}].push_back(a);
        }
        // Trp bridgeheads also live in the secondary ring:
        if (sem.ring_position.secondary.ring != RingSystemKind::NotInRing) {
            groups[{p.atoms_[a].residue_index,
                    sem.ring_position.secondary.ring}].push_back(a);
        }
    }
    for (auto& [key, atom_indices] : groups) {
        auto type_idx = RingSystemToRingTypeIndex(key.second,
            p.residues_[key.first].variant_idx);
        rings.push_back(Ring{type_idx, atom_indices, key.first, ...});
    }
    return rings;
}
```

**RingSystemKind → RingTypeIndex map** (small typed switch, no strings):
| `RingSystemKind` | `RingTypeIndex` |
|---|---|
| `Benzene_Phe` | `PheRing` |
| `Benzene_Tyr` | `TyrRing` |
| `Imidazole_His` (variant_idx 0/HIE) | `HieImidazoleRing` |
| `Imidazole_His` (variant_idx HID) | `HidImidazoleRing` |
| `Imidazole_His` (variant_idx HIP) | `HipImidazoleRing` |
| `Indole_Trp_5` | `TrpPyrroleRing` |
| `Indole_Trp_6` | `TrpBenzeneRing` |
| `Pyrrolidine_Pro` | (skip Phase 1 — no calculator consumes saturated
    ring as Ring object today; Phase 2 if needed) |

**Verification gate**: bit-identical `Ring.atom_indices` ordering
versus the string-matched populator's output. Run against
`1ubq_protonated.pdb` + `fleet_amber/1Z9B` + `fleet_amber/1P9J` — for
every Ring object produced by the old populator, the new populator
must produce a Ring with identical `atom_indices` (potentially up to
permutation, but ideally byte-identical — Q3's agent named the
atom-ordering convention as the gate).

The five ring-current calculators (BiotSavart, HaighMallion,
RingSusceptibility, Dispersion, +1) consume `Ring.atom_indices`
typed; they should produce bit-identical output if the new populator
preserves ordering. Smoke tests verify.

**`AminoAcidType::rings[].atom_names` disposition**: removed entirely
(user decision per §15). The substrate carries the chemistry; the
const-char-literal list becomes redundant. `RingTypeIndex` stays on
`AminoAcidType::rings[]` (already typed).

**Dependencies**: blocks on 13.6 (substrate must be populated before
`ConstructRingsFromSubstrate` can read `SemanticAt`). 13.4 lands
after 13.6.

**Estimated size**: ~80 lines new (`ConstructRingsFromSubstrate` +
the type map); ~100 lines removed (`DetectAromaticRings` body +
`AminoAcidType::rings[].atom_names` literal arrays).

#### 13.5 Two error-message edits in `AIMNet2Result.cpp` and `ApbsFieldResult.cpp`

**Files**: `src/AIMNet2Result.cpp:169` (one line),
`src/ApbsFieldResult.cpp:140` (one line). Plus a small new helper.

**Pattern**: each calculator has a single diagnostic error message
that today renders the offending atom by `atom.pdb_atom_name`. After
13.6 lands and atoms carry typed identity, the message renders by
typed-identity-to-string projection at the diagnostic boundary.

**Helper signature** (one suggested location: free function in
`src/SemanticEnums.h` adjacent to the `AtomMechanicalIdentity` struct,
or static method on `LegacyAmberTopology`):

```cpp
// Render an AtomMechanicalIdentity to a human-readable string for
// diagnostic output. NEVER read at calculator surface for chemistry.
std::string RenderIdentity(const AtomMechanicalIdentity& id,
                           AminoAcid residue_type,
                           uint8_t variant_idx);
```

The rendering is the inverse of `ParseAtomName`: given the typed
identity, produce the canonical AMBER name. This is a
typed-to-string projection at a diagnostic-output boundary
(Crystal Projection Rule); calculators consume the typed identity,
the diagnostic-render utility produces the string only when a
human needs to read it.

**User decision pending**: where does the helper live?
- Free function in `SemanticEnums.h` neighbourhood (closest to the
  identity types).
- Static method on `LegacyAmberTopology` (substrate-side projection).
- Standalone utility in `src/AtomNameProjection.{h,cpp}` (its own
  module).

Q4's agent flagged this as user-decision Q1.

**Verification gate**: AIMNet2 + APBS unit tests pass with diagnostic
messages now using rendered identity; no `pdb_atom_name` read
in either file.

**Dependencies**: blocks on 13.1 (the helper inverts the parser; both
must exist) and 13.6 (typed identity must be populated on Atom).

**Estimated size**: 2 lines changed in calculators + ~30 lines for
the helper.

#### 13.6 Wire `Protein::FinalizeConstruction` to call `LookupBy` + `ApplyCapDelta`

**Files**: `src/Protein.cpp` (`FinalizeConstruction`),
`src/LegacyAmberTopology.{h,cpp}` (per-atom semantic store +
accessor), `src/Atom.h` (typed identity field added; `pdb_atom_name`
NOT removed yet — Phase 1 keeps it until all consumers migrate, then
removes in a follow-up step at end of Phase 1).

**The structural prevention point**: per the user's framing —
strings at the load boundary, typed at the runtime model.
`Atom.pdb_atom_name` stays for now (still populated by loaders,
read by load-boundary code) but is structurally inaccessible from
calculator-side code at the END of Phase 1 — either by being moved
to a load-only namespace or by being deleted with all consumers
already migrated. The decision: delete entirely (per user's
direction); the migration must therefore complete before the field
is removed.

**Composition algorithm** in `FinalizeConstruction` (per dependencies
§H.5):

```cpp
void Protein::FinalizeConstruction(/* existing args */) {
    // ... existing finalisation (LegacyAmberInvariants, etc.) ...

    // Populate per-atom semantic store via structural matching.
    std::vector<AtomSemanticTable> semantic;
    semantic.resize(atoms_.size());

    for (size_t i = 0; i < residues_.size(); ++i) {
        const Residue& res = residues_[i];
        const TerminalState n_state = ComputeNTermState(res);
        const TerminalState c_state = ComputeCTermState(res);

        for (size_t a : res.atom_indices()) {
            // Atom carries typed identity from the loader (set via
            // ParseAtomName at load boundary, step 13.1 + loader
            // changes).
            const auto& ident = atoms_[a].identity();

            const AtomSemanticTable* base =
                topology_generated::LookupBy(res.type, res.variant_idx, ident);

            if (base != nullptr) {
                semantic[a] = *base;
            } else if (IsCapOnlyAtom(ident)) {
                const auto* cap = topology_generated::LookupCap(
                    AtomTerminalState(ident, n_state, c_state), ident);
                if (cap == nullptr) {
                    Fail("Unresolved cap-only atom in residue index " +
                         std::to_string(i));
                }
                semantic[a] = *cap;
            } else {
                Fail("Unresolved standard atom in residue index " +
                     std::to_string(i));
            }

            // Backbone-cap composition: terminal-residue backbone N (NTERM)
            // or backbone C/O (CTERM) get cap-delta overlay on top of chain.
            if (IsBackboneCapOverlayAtom(ident)) {
                const auto* cap = topology_generated::LookupCap(
                    AtomTerminalState(ident, n_state, c_state), ident);
                if (cap != nullptr) {
                    topology_generated::ApplyCapDelta(semantic[a], *cap);
                }
            }
        }
    }

    legacy_amber_.SetAtomSemantic(std::move(semantic));
}
```

**Fail-loudly** on unresolved standard atoms — codex SDET signoff
condition. Use the project's existing `Fail(...)` invariant-violation
pattern (per `feedback_no_attach_lifecycle_for_invariant_data`); do
not throw `std::runtime_error` bare.

**Atom typed-identity field**: `Atom` gains
`AtomMechanicalIdentity identity_;` populated at load. The string
`pdb_atom_name` stays for the duration of Phase 1; deletion happens
at the end-of-Phase-1 cleanup once all consumers have migrated.

**Loader changes** (one line per loader): each loader that today
writes `atom.pdb_atom_name = X;` also writes
`atom.identity_ = ParseAtomName(X, parent_name);`. Five loaders to
update (`PdbFileReader`, `FullSystemReader`, `TrajectoryProtein`,
`OrcaRunLoader`, `KamlProtonator`/`PropkaProtonator`).

**Reusable from stash@{0}**: the corpus/audit test scaffolding
(load representative AMBER fixtures and assert per-atom semantic
record) is reusable; the parser-lifted-into-header logic is
reusable. Do NOT reuse the FinalizeConstruction wiring as-is — it
parsed `pdb_atom_name` *at FinalizeConstruction* (calculator-adjacent
code path), not at the loader (boundary). The principle was right;
the implementation crossed the barrier.

**Verification gates**:
- Corpus/audit test passes against 1ubq + fleet_amber 1Z9B + 1P9J.
- Spot-check: NTERM N has `formal_charge=+1` after cap composition;
  CTERM C/O have `planar_group=Carboxylate`; HID Nδ1 has
  `Heteroatom_NH`; LYN HZ2/HZ3 have `AmineNH`; TYM Oη has
  `AromaticOxide`; GLY HA2 has `ProchiralStereo::ProR`.
- ctest 352/352 + StringBarrier 5/5.
- `nm libnmr_shielding.a | grep _ZN5RDKit | wc -l` returns 0.

**Dependencies**: blocks on 13.1 (parser lifted) and 13.3 (alias
collapse centralised so loaders write canonical names). Unblocks
13.4 (substrate-driven ring construction).

**Estimated size**: ~120-180 lines (FinalizeConstruction extension +
helpers + semantic store on LegacyAmberTopology + identity field on
Atom + 5 loader one-liner additions + corpus/audit test).

#### 13.7 End-of-Phase-1 cleanup: delete `pdb_atom_name` from `Atom`

**Files**: `src/Atom.h` (delete the field), every consumer that the
audit identified, the SDK/H5 boundary if it reads from Atom directly.

**Order**: this step lands LAST. Until every Category-C and
Category-E consumer has migrated to typed identity (steps 13.1-13.6),
the field stays. Only when no calculator-side code reads it does the
field come off `Atom`.

**Migration verification**: `grep -rn "atom.pdb_atom_name\|atom_p->pdb_atom_name"
--include='*.cpp' --include='*.h' src/` returns zero hits *outside*
loader files (and the loaders are the only legitimate writers — they
will be reading their own load-boundary string into typed identity
via 13.1's parser).

**Once green**: delete the field. ctest stays green; StringBarrier
stays 5/5.

**Estimated size**: deletion is mechanical; verification is the work.

#### 13.x Cross-step risk: substrate-population vs backbone-cache ordering

`Protein::CacheResidueBackboneIndices` runs in `FinalizeConstruction`
ahead of where 13.6 would populate the substrate. Two options:

1. **Reorder**: cache backbone indices AFTER substrate is populated.
   The cache code reads `atom.semantic.backbone_role` typed; no
   dependency on the string.
2. **Loader-side BackboneRole**: each loader computes
   `BackboneRole` at load time (parser output, not substrate query)
   and stores on Atom directly. Cache reads `atom.backbone_role`
   typed before substrate population.

Q2's agent recommended option 2 (loader-side). It's smaller and
removes the chicken-and-egg without reordering existing init code.
Adopt unless 13.6's wiring naturally orders before backbone caching.

### 14. The Phase 2 calculator walk

Phase 2 is the deeper revisit. Q1, Q2, Q3 each flagged Phase 2
considerations that don't dissolve under Phase 1 alone:

- **Time-dependent protonation under constant-pH MD**: the substrate
  is invariant per assumption; per-frame protonation state would
  need a per-frame ConformationResult companion. Today's design
  doesn't support this.
- **Tautomer ambiguity in crystal-only structures**: the AMBER LEAP
  default of HIS → HIE quietly resolves, but the "confidence" of
  that resolution isn't surfaced.
- **Ring puckering and ring-flip dynamics**: substrate-side
  RingPosition + per-frame `PlanarGeometryResult` ConformationResult
  was the design (per dependencies §H.5); the per-frame companion is
  pending.
- **Atom-typing audit across 13 calculators**: confirm typed surfaces
  hold under the substrate-and-no-atom_name regime. The expected
  outcome is "they already hold" (per the serendipity finding) but
  Phase 2 is the audit pass that confirms.

Phase 2 is *less work than feared* because of the calculator-floor
serendipity. But it's not zero.

### 15. User decisions still open

Aggregated from the four chemistry questions and the audit:

**Architectural**:
- Where does the runtime `ParseAtomName` live: in
  `LegacyAmberSemanticTables.h` (single source of truth with the
  generator), or in `NamingRegistry` (centralises alias collapse)?
  Q2 recommended NamingRegistry; both have merit.
- Does `Atom` carry typed `backbone_role` / `identity` directly, or
  do these live as parallel arrays on `LegacyAmberTopology`? Q2's
  recommendation is `BackboneRole` directly on `Atom` for fastest
  access in backbone-dense calculators.
- Is `pdb_atom_name` deleted in Phase 1, or quarantined behind a
  load-only friend-restricted accessor with deletion in Phase 2?
  The user has stated firmly: **deletion**. Q2's recommendation was
  "yes, deferred to Phase 7"; the user's answer narrows that.

**Tactical**:
- Where does the typed-identity error-rendering helper live (Q4)?
- Where does `AmberAtomNameFor(identity, residue_type, terminal_state)`
  live (Q4)?
- Does `AmberPreparedChargeSource`'s PRMTOP cross-walk move to
  typed-identity-keyed maps, or stay atom-name-string-keyed at the
  ChargeSource boundary? Q4's recommendation: stay string-keyed at
  the boundary; the wire format is its authority.
- Glycine HA cache disposition (Q2): special-case Gly with typed
  Locant lookup, since Gly's HA2/HA3 carry `Locant::Alpha` not
  `BackboneRole::AlphaHydrogen`.
- CB cache disposition (Q2): defensive name-match as Phase 1
  workaround, fixed in Phase 2 via Locant::Beta substrate.
- Does `DetectAromaticRings` survive as a function post-Phase-1, or
  dissolve into substrate population (Q3)?
- Should `AminoAcidType::rings[].atom_names` be removed entirely
  (Q3)?
- Should PROPKA wiring be the Phase 1 baseline (Q1) — or a typed
  fallback inferring from atom-presence — or both?

**Cross-cutting**:
- Substrate population timing: `LookupBy` at `FinalizeConstruction`
  is specified but not yet called. Phase 1 step 13.6 wires this.
- SDK catalog migration: when `pdb_atom_name` comes off the H5 column
  (or is renamed to a clear load-time-only artefact), downstream
  Python consumers update.
- AminoAcidType::rings[].atom_names + ChiAngleDef::atoms[4] (currently
  `const char*`): become typed identities or stay string-keyed with
  projection-at-use?

### 16. Risks and unknowns

- **Substrate population vs cache population ordering**: backbone
  caching today happens in `FinalizeConstruction` ahead of where
  substrate population would land (because the cached indices are
  used by other invariant-data setup). Either reorder, or have the
  loaders compute `BackboneRole` at load time (so the cache can
  use typed role at populate time, before substrate is populated).
  Q2 recommends the latter.
- **The stashed agent-integration code**: useful pieces include the
  parser-lifted-into-header and the corpus/audit test scaffolding.
  Do not reuse the FinalizeConstruction wiring as-is — it consumed
  `pdb_atom_name` to derive identity. The principle was right; the
  implementation crossed the barrier.
- **Phase 1 must not break the SDET's signoff conditions**: use
  `ApplyCapDelta` not whole-row; fail loudly on unresolved standard
  atoms; corpus/audit test against AMBER fixtures.
- **OrcaRunLoader latent bug**: the H||HN||H1 conflation will
  silently change behaviour after centralisation (H1 properly
  resolves to cap atom, not backbone H). Verify with corpus tests
  before declaring Phase 1 complete.
- **The SDK reads `pdb_atom_name` from H5**: the H5 emission boundary
  is a wire-format projection. Phase 1 keeps the string at that
  emission point (it's a wire boundary); the runtime model surface
  is what changes. SDK readers continue to work.

---

## Part IV — How to pick up

### 17. Reading order for a fresh session

To rejoin this work, read in this order:

1. **`spec/INDEX.md`** — project foundation pointer.
2. **`OBJECT_MODEL.md`** + **`PATTERNS.md`** + **`spec/CONSTITUTION.md`**
   — sacred docs.
3. **This document** — the synthesis (you're reading it).
4. **`spec/plan/mechanical-identity-model-and-audit-2026-05-05.md`**
   — the audit's full definition + 67-call-site classification.
5. **`spec/plan/chemistry-question-1-protonation-state-2026-05-05.md`**
   — Q1 detailed plan.
6. **`spec/plan/chemistry-question-2-backbone-topology-2026-05-05.md`**
   — Q2 detailed plan.
7. **`spec/plan/chemistry-question-3-ring-topology-2026-05-05.md`**
   — Q3 detailed plan.
8. **`spec/plan/chemistry-question-4-external-model-typing-2026-05-05.md`**
   — Q4 detailed plan (note: Q4 mostly dissolves).
9. **`spec/plan/topology-encoding-dependencies-2026-05-05.md`**
   §H — substrate architecture.
10. **`session-handoff-20260505-evening.md`** — the substrate's end-of-day
    state.

### 18. Verification commands

```bash
# 1. HEAD invariant.
git -C /shared/2026Thesis/nmr-shielding log --oneline -10
# Expect: e7e47e6 most recent; 57db9f6 / cc480c8 / ee1f1b4 / ab3bd5e /
# 56fb6b9 / bb4584d above (substrate landings); b079d4a (morning handoff).

# 2. Working tree clean — substrate state.
git -C /shared/2026Thesis/nmr-shielding status --short
# Expect: empty (the agent's runtime integration is in stash@{0}).

# 3. Stash preserved if any pieces are reusable.
git -C /shared/2026Thesis/nmr-shielding stash list
# Expect: stash@{0} present, message "agent-integration-WIP: substrate
# wired into Protein::FinalizeConstruction via pdb_atom_name parsing
# — pivoted away pending mechanical-identity audit"

# 4. Substrate tests pass.
ctest --test-dir /shared/2026Thesis/nmr-shielding/build -j 8
# Expect: 100% tests passed, 0 tests failed out of 352

# 5. StringBarrier holds.
ctest --test-dir /shared/2026Thesis/nmr-shielding/build -R StringBarrier
# Expect: 5/5 pass.

# 6. Linker-level barrier intact.
nm /shared/2026Thesis/nmr-shielding/build/libnmr_shielding.a | grep -c '_ZN5RDKit'
# Expect: 0.

# 7. Generated table integrity.
grep -c "^// === " /shared/2026Thesis/nmr-shielding/src/generated/LegacyAmberSemanticTables.cpp
# Expect: 30 (20 standard + 10 variants); plus 4 cap tables visible
# under "// === Terminal-state cap tables ===".

# 8. Plan documents.
ls spec/plan/topology-*-2026-05-05.md spec/plan/chemistry-question-*-2026-05-05.md
spec/plan/mechanical-identity-model-and-audit-2026-05-05.md
spec/plan/topology-and-identity-pivot-synthesis-2026-05-05.md
spec/plan/session-handoff-20260505*.md
# Expect: this file plus 4 chemistry questions plus audit plus 5 topology
# docs plus 2 session handoffs.
```

### 19. The HEAD invariant

Master is at **`e7e47e6`** — Phase 5 substrate model complete,
no string-on-Protein integration. The agent's runtime integration
(which crossed the barrier via `pdb_atom_name`) is preserved in
`stash@{0}` for reference but should NOT be applied as-is.

### 20. The Phase 1 work order, when ready

In the order they need to land, with rough size:
- 13.1 Lift parser into runtime: ~50-80 lines (in NamingRegistry or
  generated header).
- 13.3 Centralise alias collapse in NamingRegistry: ~30-50 lines.
- 13.2 Wire PROPKA/KaML at PdbFileReader: ~30 lines (the TODO
  resolution).
- 13.6 Wire `LookupBy` + `ApplyCapDelta` in `FinalizeConstruction`:
  ~80-120 lines (with corpus/audit test).
- 13.4 Replace `DetectAromaticRings` with `ConstructRingsFromSubstrate`:
  ~60-80 lines.
- 13.5 AIMNet2/APBS error-message edits: ~10 lines combined.
- Then: `ResolveProtonationStates` shrinks; `DetectHisVariantFromAtoms`
  shrinks/dissolves; backbone-cache code reduces; `pdb_atom_name`
  comes off the runtime model.

Total: ~300-400 lines net change across 5-7 files. Plus tests.

### 21. The shape of "done"

Phase 1 done when:
- `nm libnmr_shielding.a | grep pdb_atom_name` returns 0 calls from
  any non-loader file. (The string lives at load boundary, dies at
  load.)
- `Atom` does not have a public `pdb_atom_name` field.
- All 8 Category-C and 5 Category-E sites from the audit are
  refactored to typed-identity or substrate-data consumption.
- ctest stays green; StringBarrier stays 5/5; corpus/audit test
  passes against AMBER fixtures.
- The four chemistry questions' Phase 1 designs are landed.

Phase 2 begins after Phase 1 lands and the SDET signs off again.

---

## Closing note

The journey today: from "scale the substrate to 19 residues" through
two OpenAI eval rounds, through a runtime-integration attempt that
revealed the string-on-Protein structural risk, through an
identity-audit that defined the model concept, through four parallel
chemistry-question agents that revealed the calculator floor is
already typed-clean.

The user's framing throughout: structure-not-discipline; types-not-strings;
algorithms doing science not powering through; the model carrying
meaning explicitly via typed structure rather than implicitly via
name patterns.

The substrate model is built. The pivot is well-considered, the
phases are sequenced, the user decisions are explicit, and the
calculator floor is more cooperative than feared. Phase 1 is a small
focused refactor; Phase 2 is the deeper walk that builds on it.

When picking this up, the work is well-scoped enough that one or two
agent runs (under supervision) should cover Phase 1 substantively.
The calculators that worried us most (AIMNet2, APBS, DSSP, the ring
currents) are already typed-surface-clean. The string concentration
is at four named hotspots; the fix per hotspot is well-articulated.

This is good ground to stop on.
