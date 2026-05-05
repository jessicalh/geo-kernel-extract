# Topology encoding — application-time dependencies (2026-05-05)

Items in this file are decision points that the encoder (the next session
generating `src/generated/LegacyAmberSemanticTables.cpp` from the residue
reference) MUST consider, AND/OR that downstream calculator consumers MUST
know about when reading the corresponding `LegacyAmberTopology` fields.

Companions:
- `spec/plan/topology-residue-reference-2026-05-05.md` — the per-(residue,
  variant, atom) reference. The "what to encode."
- `spec/plan/topology-reference-review-2026-05-05.md` — the critical-review
  signoff. The "we have agreed to encode it."
- This file — the "watch-outs at the time of doing the encoding, and at
  the time downstream code reads the values."

The point of a separate file: the reference is a static document of values;
this file captures the decisions that *depend* on something at the moment
of application. Without it, the dependencies hide inside the choices and
get rediscovered painfully later.

---

## A. Vocabulary preconditions — apply to `src/SemanticEnums.h` BEFORE encoding

The reference document uses two enum values that don't exist in
`SemanticEnums.h` yet. The C++ generator cannot reference them until the
header is extended. Apply these first; commit; THEN regenerate.

### A.1 `PolarHKind::AmineNH`

Add as a new enum value. Suggested ordinal: 11 (insert before `OtherPolarH`,
renumbering it to 12). Docstring:

> Neutral primary-amine N-H. Lys neutral form (LYN: HZ2, HZ3) and the
> NTERM_NEUTRAL state (H1, H2) on the backbone N. Distinguished from
> `AmmoniumNH` because the chemistry differs (pKa ~9-10 for primary amine
> vs always-deprotonated water-exchange for ammonium); SHIFTX2 (Han 2011)
> trains separate models for these environments. Wüthrich, "NMR of
> Proteins and Nucleic Acids" (1986) ch. 2; Englander, Annu. Rev.
> Biophys. Biomol. Struct. 39 (2008) 289-307 for exchange framework.

Consumers of `PolarHKind::AmmoniumNH` must check whether they meant
"protonated ammonium specifically" or "any -NH-on-N polar H" — if the
latter, they need to add `AmineNH` to their match-set.

### A.2 `PlanarGroupKind::AromaticOxide`

Add as a new enum value. Suggested ordinal: 9 (after `AromaticHydroxyl`).
Docstring:

> Aryloxide. The deprotonated phenolate -O- on TYM (Tyr deprotonated
> form). The phenolate conjugates more strongly with the aromatic ring
> π system than neutral phenol-OH does, redistributing ring current and
> shifting the para-Cζ environment by 5-10 ppm relative to neutral Tyr.
> Distinguished from `AromaticHydroxyl` because the H is absent and the
> conjugation is qualitatively different.

### A.3 Documentation comments (cosmetic, no behaviour change)

Per review item 11d and cross-cutting observations:

- `RingPosition` struct comment: "An atom can be in a ring AND have
  Locant::None — the two fields are orthogonal: Locant tells you the
  backbone-vs-side-chain Greek-letter position; RingPosition tells you
  ring membership. Pro Cα is in the pyrrolidine ring with Locant::None."
- `RingPositionLabel::Ortho1` / `Ortho2`: extend comment to cover the
  Trp 6-ring perimeter convention (Ortho1=Cε3, Ortho2=Cζ2 — see Section D).
- `PseudoatomKind` comment block: "N-terminal H1/H2/H3 are encoded as
  Q-with-Locant=None, treating them as a 3-H equivalent group. Markley
  Table 1 does not list a pseudoatom for these terminus-specific Hs;
  the Q taxonomy is the closest fit and is used here without a new
  letter."

---

## B. External-authority dependencies — encoder MUST verify against the live source

Values in the reference document are not always self-sufficient: some
require cross-checking against an authority that may have evolved since
the reference was written. The encoder verifies each before emitting.

### B.1 ARN atom inventory

**Authority**: `src/AminoAcidType.cpp::AMINO_ACID_TYPES` (the
`AmberAminoAcidVariantTable`) at code-gen time.

**Reference's commitment**: removes HE for ARN (chemistry-inferred —
the lone pair after deprotonation goes to Nε, the formerly-protonated
nitrogen).

**Encoder rule**: read the variant table. If ARN has an explicit atom
list, the table wins and the reference is overridden. If ARN is absent
or has the chemistry-inferred default, the reference's HE-removal stands.
If the project later adopts a propka-style on-the-fly inference or a
custom ARN template, regenerate the ARN row from that authority.

**Caveat**: AMBER ff14SB has no canonical ARN template at present (per
review item 2 verdict). The reference's HE-removal is the encoder's
commit-rule unless and until that changes.

### B.2 LYN HZ atom names

**Authority**: AMBER ff14SB at `data/ff14sb_params.dat` lines 449-451.

**Authoritative atom inventory**: HZ2, HZ3 (no HZ1).

**Reference's commitment**: needs correction — currently says "HZ3
removed, leaving HZ1 + HZ2"; should say "HZ1 removed, leaving HZ2 + HZ3."

**Encoder rule**: trust ff14sb_params.dat. The reference will be patched
in the same edit pass that adds AmineNH (see Section A.1).

### B.3 PHE C++ acceptance-gate dead code

**Authority**: `tools/topology/build_semantic_tables.cpp::SynthesisedForPhe`
lines 779-842. PHE has no atom named "HG" — CG is the ipso carbon and
has no H because it's substituted.

**Encoder action**: remove `|| atom_id == "HG"` from line 794 (PHE has
no HG atom; dead code with future generalisation risk). After removal,
re-run the PHE acceptance check on the residue reference — should still
match.

---

## C. Reading-convention dependencies — consumer calculators MUST know

These are conventions baked into the substrate that downstream code reading
the `LegacyAmberTopology` fields must be aware of. Mis-aligned consumers
will silently produce wrong answers.

### C.1 Charge placement is Lewis-localised, not delocalised

The substrate's `formal_charge` is a per-atom int8_t. For groups where
real-world charge distribution is delocalised, the substrate exposes the
**Lewis-structure** value:

- Asp / Glu / CTERM_DEPROTONATED carboxylate: -1 on Oδ2 / Oε2 / O'';
  Oδ1 / Oε1 / O' = 0.
- Arg guanidinium: +1 on Nε; Cζ, Nη1, Nη2 = 0.
- HIP imidazolium: +1 on Nε2; Nδ1, Cε1, Cδ2, Cγ = 0.
- Lys ammonium / NTERM_CHARGED: +1 on N; Hζ1/2/3 (or H1/2/3) = 0.
- CYM thiolate: -1 on Sγ.
- TYM aryloxide: -1 on Oη.

**What this means for consumers**: a calculator that wants the time-
averaged or symmetrised charge across a delocalised group MUST apply its
own transformation. The substrate is the topology label, not the live
electron distribution.

**Documented in**: §2 of `topology-residue-reference-2026-05-05.md`
(Charge-distribution conventions).

### C.2 E/Z mapping is residue-specific (Asn/Gln vs Arg)

After the corrections in Section E, the substrate's `PlanarStereo`
follows BMRB `atom_nom.tbl` conventions:

- Asn: HD21 = E, HD22 = Z.
- Gln: HE21 = E, HE22 = Z.
- Arg: HH11 = Z, HH12 = E, HH21 = Z, HH22 = E.

The asymmetry (Asn/Gln have "1 = E"; Arg has "1 = Z") is a real
chemistry consequence: for Asn/Gln, the heavy atom 3 bonds in (Cβ for
Asn, Cγ for Gln) is the *lower-priority* substituent on the planar
carbon (the carbonyl O is higher); for Arg, the heavy atom 3 bonds in
(Nε) IS the higher-priority substituent on Cζ. Markley's caption gives
the *labelling rule* (1 = cis to heavy atom 3 bonds in); CIP determines
the Z/E from priority.

**What this means for consumers**: do NOT assume a single rule like
"odd-numbered = Z" or "1 = always cis = Z." Read the substrate's
`PlanarStereo` field directly per atom — it is residue-aware.

### C.3 Trp 6-ring perimeter labels are synthesised, not Markley-canonical

Markley does not publish ipso/ortho/meta/para for the indole 6-ring
perimeter. The substrate uses Ortho1=Cε3, Ortho2=Cζ2, Meta1=Cζ3,
Meta2=Cη2 (synthesised; per the dossier's convention). No literature
contradiction.

**What this means for consumers**: if a calculator stratifies Trp ring
atoms by ring position, the labels above are stable and consistent
within the substrate but may not align with publications using
explicitly-labelled "perimeter" or "C2/C3/C4/C5" conventions. Document
the Markley-extension convention in the calculator's notes.

### C.4 HIP +1 placement is on Nε2 specifically

HIP imidazolium charge +1 is localised on Nε2. Alternative conventions
(localise on Nδ1, distribute) exist; the substrate has chosen Nε2.

**What this means for consumers**: a calculator reading
`formal_charge` across HIP atoms will see +1 on Nε2 and 0 elsewhere.
For symmetrised treatment, transform; or reach into the substrate-side
PlanarGroupKind=Imidazole atoms and apply a delocalisation rule.

---

## D. Locked architectural decisions

These are decisions that affect encoding but should not be re-litigated
without explicit user direction. Listed here to immunise against
context-to-context drift (the project lead has flagged this risk).

### D.1 AMBER ff14SB default for "HIS" = HIE

A residue parsed as the bare 3-letter code "HIS" without a tautomer
override has the substrate row identical to the HIE row (Nε2-
protonated, no Hδ1, neutral). HID and HIP variant blocks list deltas
vs that HIE/HIS-default reference.

**What changes if this is revisited**: the entire HIS default block in
the reference document, plus the HID/HIE/HIP variant deltas, would
need to be regenerated against the new default tautomer.

### D.2 Pro Cα is in the pyrrolidine ring with Locant::None

Locant is for sidechain atoms (Greek-letter rule); RingPosition is
ring-membership. The two fields are orthogonal. Pro Cα has Locant::None
AND ring_primary = Pyrrolidine_Pro/Saturated/5/f/1.

This is documented in SemanticEnums.h per Section A.3 above.

### D.3 Markley pseudoatom taxonomy is M / Q / R only

XPLOR / CARA / DIANA introduce additional letters (X, Y, S, T) but
those are not IUPAC and are not adopted. For the terminal H1/H2/H3
case (where Markley's Table 1 doesn't define a pseudoatom), use
`Q/Locant::None/0/false` as the closest-fit encoding.

### D.4 String barrier at linker level

The runtime library (`libnmr_shielding.a`) does not link RDKit. The
generator (`tools/topology/build_semantic_tables`) is the only build
artefact that does. The chemistry-strings boundary is enforced by
`tests/test_string_barrier.cpp`; consumers must not bypass.

---

## E. Reference-document corrections to apply

The reviewer flagged corrections to apply to the residue reference
document itself. Listed here for completeness; these are *fixes*, not
*ongoing dependencies*. Apply once.

1. **Asn HD21**: planar_stereo Z → E.
2. **Asn HD22**: planar_stereo E → Z.
3. **Gln HE21**: planar_stereo Z → E.
4. **Gln HE22**: planar_stereo E → Z.
5. **LYN variant block**: change "HZ3 REMOVED, leaving HZ1+HZ2" →
   "HZ1 REMOVED, leaving HZ2+HZ3" (per `data/ff14sb_params.dat`).
6. **LYN HZ2 / HZ3 rows**: polarH AmmoniumNH → AmineNH (after enum added).
7. **NTERM_NEUTRAL H1 / H2 rows**: polarH AmmoniumNH → AmineNH.
8. **TYM Oη row**: planar_group AromaticHydroxyl → AromaticOxide.
9. **§2 charge-distribution conventions block**: rewrite the 1=Z/2=E
   sub-clause for Asn/Gln; explain the residue-specific Z/E mapping
   per Section C.2 above. Cite BMRB atom_nom.tbl as the canonical
   E/Z source.
10. **Disagreement log entries 1, 2**: invert resolution. Caption
    gives the labelling rule (1 = cis to heavy atom 3 bonds in); CIP
    determines the Z/E from priority on each end of the planar bond.
    BMRB atom_nom.tbl is the authoritative cross-check.
11. **Disagreement log entries 8, 9**: tighten — the AmineNH and
    AromaticOxide additions are now ACTION items in this dependencies
    file, not just open flags.

---

## F. Open question: ARN remains a flag, not a verdict

Per review item 2: AMBER has no canonical ARN template, and the
reference's chemistry-inference (HE removed) is defensible but not
literature-authoritative. The encoder commits this rule UNLESS a future
project decision lands a propka-style inference or custom ARN template.
If that happens, regenerate the ARN row from the new authority.

This is not an unresolved verdict — it's an explicitly-flagged degree
of freedom that the encoder commits to with its own rule when it
generates `LegacyAmberSemanticTables.cpp`.

---

## G. Variant-table synthesis architecture (PENDING — discovered 2026-05-05 evening)

**Status**: the ten AMBER protonation variants (HID, HIE, HIP, ASH,
GLH, CYX, CYM, LYN, ARN, TYM) are NOT yet emitted into
`src/generated/LegacyAmberSemanticTables.cpp`. The standard 20 are.
Per-variant `SynthesisedFor<Variant>` patch functions exist in the
generator (compiled but unused, marked `[[maybe_unused]]`) and carry
the correct delta logic; they need a different upstream feeder.

**The discovered constraint**: AMBER's protonation-variant 3-letter
codes (HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN/ARN/TYM) are an AMBER-internal
naming convention. They are NOT canonical wwPDB Chemical Component
Dictionary entries for those chemistries. The CCD does happen to
contain `data_X` blocks for several of these codes — but those blocks
refer to UNRELATED small molecules that share the 3-letter code by
historical accident:

| AMBER variant | CCD's `data_X` actually contains |
|---|---|
| HID | (5-hydroxy-1H-indol-3-yl)acetic acid |
| HIE | tetrahydroindazol benzamide (60-atom small molecule) |
| HIP | ND1-phosphonohistidine (a phosphorylated histidine, not HIP) |
| ASH | imidazopyrazine compound |
| GLH | cyclohexyl-glutamine |
| CYX | (3-formyl-but-3-enyl)-phosphonic acid |
| CYM | S-methylcysteine |
| LYN | lysine amide |
| ARN | 1-imino-5-pentanone |
| TYM | tryptophanyl-AMP |

A naive `grep "^data_HID" components.cif` will report a hit; the entry
exists. But the chemistry inside it is not histidine.

**Implication**: variants must be synthesised, not looked up. The
generator needs a parent-CCD-plus-delta architecture:

1. For each variant, identify the parent (HID/HIE/HIP → HIS;
   ASH → ASP; GLH → GLU; CYX/CYM → CYS; LYN → LYS; ARN → ARG;
   TYM → TYR).
2. Read the parent's CCD entry and run the existing pipeline
   (cifpp + RDKit perception).
3. Apply the variant's per-atom delta from the reference doc Section 3
   variant block:
   - **Atoms ABSENT** in variant but PRESENT in parent: drop the
     entry from the table (e.g. LYN drops HZ1; TYM drops HH; CYX/CYM
     drop HG).
   - **Atoms PRESENT** in variant but ABSENT in parent: synthesise an
     atom record. The parent CCD does NOT carry these atoms (HID's
     Hδ1, HIP's Hδ1, ASH's Hδ2, GLH's Hε2). The synthesis needs:
     - Element (always H for the variant adds in scope).
     - Bond graph: the new H bonds to the variant's parent heavy atom
       (Hδ1 to Nδ1 in HID; Hδ2 to Oδ2 in ASH; Hε2 to Oε2 in GLH).
     - Hybridisation, formal_charge, equivalence_class: derive from
       chemistry (RDKit on the variant's heavy-atom topology).
   - **Field-delta atoms**: existing entries change values
     (Nδ1's RingPositionLabel flips Heteroatom_NoH → Heteroatom_NH
     in HID; Oδ2's formal_charge flips -1 → 0 in ASH).

4. Pipe the synthesised atom inventory through `BuildAtomSemanticEntry`
   with the existing `SynthesisedFor<Variant>` function as the
   chemistry-source plug-in.

**Authority for per-variant atom additions**: AMBER ff14SB at
`data/ff14sb_params.dat` is canonical. Cross-check against
`src/AminoAcidType.cpp::AMINO_ACID_TYPES`.

**Estimated effort**: one focused agent run with a corrected brief
(see `spec/plan/topology-variants-agent-brief-2026-05-05.md` if it
exists; otherwise the brief is in the session-handoff doc).

**Why this wasn't caught earlier**: the prior research dossier
(`topology-fields-research-2026-05-05.md`) and prior session handoffs
contained the implicit assumption "use the CCD entry where it exists"
without verifying chemistry. The bare presence-test (grep) was treated
as confirmation. Lesson recorded in the next-agent brief.

---

## H. Structural-matching lookup architecture (PENDING — discovered 2026-05-05 evening, OpenAI critique)

**Status**: the comment at `tools/topology/build_semantic_tables.cpp:2144`
promises "(AminoAcid, atom_local_idx)" lookup at runtime, but the
generator iterates raw CCD atom order at line 2085. Runtime
`AminoAcidType` (in `src/AminoAcidType.cpp`) carries chain-residue
inventory; CCD's free-amino-acid form is a superset (extra cap atoms
OXT/HXT/H2 plus, for ASP/GLU/HIS, the protonated forms HD2/HE2/HD1).
Direct atom_local_idx integration would silently mis-assign semantics.

**Architectural fix — locked decision (do NOT relitigate)**: the
runtime lookup MUST use **typed structural matching**, not raw
atom_local_idx and not atom_id strings. The string wall is sacred.

### H.1 Mechanical-identity tuple

Each `AtomSemanticTable` entry exposes its mechanical identity. The
parser already computes `Locant + BranchAddress + DiastereotopicIndex`;
add `Element` (already typed in `src/Types.h`) and a `BackboneRole`
enum to disambiguate the six backbone slots (which all share
`Locant::None`):

```cpp
enum class BackboneRole : uint8_t {
    None              = 0,   // Sidechain atom or cap; identity comes from
                              // (Locant, Branch, DiIndex, Element).
    Nitrogen          = 1,   // Backbone N.
    AlphaCarbon       = 2,   // Backbone CA.
    CarbonylCarbon    = 3,   // Backbone C.
    CarbonylOxygen    = 4,   // Backbone O.
    AmideHydrogen     = 5,   // Backbone H / HN.
    AlphaHydrogen     = 6,   // Backbone HA (Gly's HA2/HA3 keep BackboneRole::None
                              // because they carry Locant::Alpha).
};

struct AtomMechanicalIdentity {
    Element             element;
    Locant              locant;
    BranchAddress       branch;
    DiastereotopicIndex di_index;
    BackboneRole        backbone_role;
};
```

The `AtomMechanicalIdentity` is unique within a residue's table.

### H.2 Lookup function (runtime)

```cpp
// Linear scan over the table; N is small (≤ ~25 entries per residue).
const AtomSemanticTable*
LookupBy(AminoAcid residue, uint8_t variant_idx,
         const AtomMechanicalIdentity& identity);
```

Lookup time ≈ 25 ns per call; one call per atom at protein
construction; entirely inert at calculator runtime (lookup happens
once at `Protein::FinalizeConstruction`).

### H.3 Benign-extras tolerance

The substrate's per-residue tables MAY contain extra atoms not
present in the runtime's chain-residue inventory (CCD-form caps,
protonated atoms in CCD ASP/GLU/HIS, etc.). These are benign extras:
runtime queries the lookup function for the atoms it actually has;
extras sit unqueried in the table.

This is robustness by design: the substrate is a superset of any
runtime use-case; the runtime composes its semantic record by
querying the typed identity for each atom it owns. If the substrate
table is later filtered to chain-form (per H.4 below), benign extras
disappear; the lookup mechanism is unaffected.

### H.4 Option C: cap separation (paired with structural matching)

The standard 20 tables today carry CCD's cap atoms (OXT, HXT, H2)
with default-no-chemistry entries — the existing `SynthesisedFor<Residue>`
functions only handle `is_n_terminus` / `is_c_terminus` polar_h via
per-atom branching; other cap-specific fields (formal_charge,
pseudoatom membership, planar_group context) get default values that
are wrong for an actual terminus.

**Decision**: emit four cap tables alongside the 30 residue tables:

```cpp
constexpr std::array<AtomSemanticTable, N> kCapNtermCharged    = { ... };
constexpr std::array<AtomSemanticTable, N> kCapNtermNeutral    = { ... };
constexpr std::array<AtomSemanticTable, N> kCapCtermDeprotonated = { ... };
constexpr std::array<AtomSemanticTable, N> kCapCtermProtonated   = { ... };
```

The reference doc Section 4 (`spec/plan/topology-residue-reference-2026-05-05.md`)
has the chemistry. Generator partitions CCD atoms: chain atoms go to
per-residue tables; cap atoms go to terminal-state tables; variant
atoms stay in variant tables (already done).

### H.5 Composition at protein construction

```cpp
void Protein::FinalizeConstruction() {
    for (Residue& res : residues_) {
        for (Atom& atom : res.atoms()) {
            // Compute mechanical identity from atom name (parser).
            auto ident = ComputeAtomMechanicalIdentity(atom);

            // Standard semantic record from per-residue / per-variant table.
            const AtomSemanticTable* base = LookupBy(res.type, res.variant, ident);
            if (base) ApplySemantic(atom, *base);

            // If this residue is at a terminus, also look up cap atoms.
            if (atom.is_n_terminal_atom) {
                const auto* cap = LookupCap(res.n_term_state, ident);
                if (cap) ApplySemantic(atom, *cap);
            }
            if (atom.is_c_terminal_atom) {
                const auto* cap = LookupCap(res.c_term_state, ident);
                if (cap) ApplySemantic(atom, *cap);
            }
        }
    }
}
```

`ApplySemantic` writes the typed semantic fields onto the runtime
atom record (or a parallel `LegacyAmberTopology` substrate slot).

**Override-composition rule (per Fix 6 of the 2026-05-05 OpenAI
evaluation; live in the cap-table comment block of the generated
.cpp)**: when a residue is at a terminus, `Protein::FinalizeConstruction`
MUST query both `LookupBy` (chain) and `LookupCap` (terminal-state).
If `LookupCap` returns non-null, the cap entry's fields override the
chain entry's fields for that atom. The cap table MAY contain entries
for backbone atoms (N for NTERM, C/O for CTERM) whose chemistry
differs at the terminus; those entries have priority.

Backbone overrides have non-None `BackboneRole` matching the chain's
same-role entry; terminus-added Hs / Os (H1, H2, H3, OXT, HXT) carry
`BackboneRole::None` and only match via the cap-table lookup. The
`ApplySemantic` ordering above (cap last, after chain) implements the
override automatically: the cap-table fields land on top of any chain
values that were written first.

After Fix 1 (2026-05-05), the cap tables carry override entries:
- `kCapNtermCharged` — N (BackboneRole::Nitrogen, formal_charge=+1) +
  H1, H2, H3 (4 entries).
- `kCapNtermNeutral` — N (BackboneRole::Nitrogen, formal_charge=0) +
  H1, H2 (3 entries).
- `kCapCtermDeprotonated` — C, O (Carboxylate planar group, both
  formal_charge=0) + OXT (formal_charge=-1) (3 entries).
- `kCapCtermProtonated` — C, O (Carboxylate planar group) + OXT, HXT
  (4 entries).

### H.6 What this means for the generator

The generator must:

1. **Partition CCD atoms per residue** into:
   - Chain atoms (those in `AmberAminoAcidVariantTable`'s standard set
     for that residue) → per-residue table.
   - Cap atoms (OXT, HXT, H1, H2, H3) → terminal-state table.
   - Variant-specific atoms — already handled by variant tables.
2. **Add `Element` and `BackboneRole` to `AtomSemanticTable`**
   (and the parser computes them).
3. **Emit four cap tables** carrying the reference doc Section 4 chemistry.
4. **Re-anchor the `Lookup` function comment** at line 2144 to describe
   structural matching, not atom_local_idx.

### H.7 Rationale (why this architecture, not the alternatives)

- **NOT atom_local_idx** because index spaces don't align across CCD,
  AmberAminoAcidVariantTable, and any other producer of atom inventories.
- **NOT atom_id strings** because the string wall is sacred; the
  substrate is typed-only.
- **YES typed mechanical identity** because the parser already computes
  these fields; reusing them as the lookup key adds no string surface
  and gives unambiguous matching within a residue.
- **YES cap-table separation** because cap chemistry is residue-
  independent (NTERM_CHARGED chemistry is the same across all 20
  residues); centralising it eliminates per-residue duplication and
  makes the chain table semantically clean.

### H.8 Estimated effort

One focused agent run:
- ~50 lines new types in `SemanticEnums.h` (`BackboneRole`, `AtomMechanicalIdentity`,
  Element field on `AtomSemanticTable`).
- ~150 lines generator changes (partitioning, cap-table emission,
  parser BackboneRole computation, lookup-function emission).
- Reference doc Section 4 already has the cap chemistry.
- Tests updated.

The `Lookup` function lives in the generated `.cpp` (typed enums only;
no strings); the runtime composition lives in `Protein.cpp` (next-
session integration, not in this work).

This entry takes precedence over earlier "atom_local_idx" comments
in the codebase. Update those comments at the same time the
architectural change lands.

### H.9 LANDED 2026-05-05 evening at commit `ee1f1b4`

The structural-matching architecture is implemented:
- `BackboneRole` + `TerminalState` enums and `AtomMechanicalIdentity`
  struct in `src/SemanticEnums.h`.
- `Element` + `BackboneRole` fields on `AtomSemanticTable`.
- Chain-form filtering in `BuildResidueEntries` drops cap atoms
  (OXT/HXT/H1/H2/H3) from per-residue tables.
- Four cap tables emitted (`kCapNtermCharged`, `kCapNtermNeutral`,
  `kCapCtermDeprotonated`, `kCapCtermProtonated`) per Section 4 of
  the residue-reference doc.
- `LookupBy(AminoAcid, variant_idx, AtomMechanicalIdentity)` and
  `LookupCap(TerminalState, AtomMechanicalIdentity)` functions in
  the generated `.cpp`.
- Stale "atom_local_idx" comments refreshed.
- Tests 352/352 pass; StringBarrier 5/5.

### H.10 LANDED 2026-05-05 (Fix 4 of OpenAI evaluation): chain-atom whitelist filter

The §H landing partitioned out cap atoms but did NOT partition out
variant-specific atoms (HD2 in ASP, HE2 in GLU, HD1 in HIS) from the
standard-residue tables, and missed PRO's backbone H. CCD's free-
amino-acid form has these atoms (CCD ASP carries the protonated HD2;
CCD HIS carries both HD1 and HE2; CCD PRO carries the backbone amide
H) but AMBER's chain-residue inventory for the standard charged /
HIE-default / Pro-as-secondary-amine forms does not.

The benign-extras tolerance from §H.3 made this functionally inert
at runtime, but for "clean" encoding the standard-residue tables
should match the runtime data shape exactly.

**Implementation**: a `ChainAtomWhitelist(residue_3letter)` function
in `tools/topology/build_semantic_tables.cpp` returns the per-residue
chain-atom set, mirroring `AmberAminoAcidVariantTable` from
`src/AminoAcidType.cpp::AMINO_ACID_TYPES`. (Duplicated rather than
linked — the build-time generator can't pull AminoAcidType.cpp
without significant runtime deps.) `BuildResidueEntries` filters CCD
atoms against this whitelist on the standard-20 path; the variant
path passes an empty whitelist (no filtering) so variant tables
retain their full atom set.

**After Fix 4 (counts in regenerated `LegacyAmberSemanticTables.cpp`)**:

  ALA: 10 (unchanged)
  ARG: 24 (unchanged)
  ASN: 14 (unchanged)
  ASP: 12 (was 13; HD2 partitioned to ASH)
  CYS: 11 (unchanged)
  GLN: 17 (unchanged)
  GLU: 15 (was 16; HE2 partitioned to GLH)
  GLY: 7  (unchanged)
  HIS: 17 (was 18; HD1 partitioned to HID/HIP)
  ILE: 19 (unchanged)
  LEU: 19 (unchanged)
  LYS: 22 (unchanged)
  MET: 17 (unchanged)
  PHE: 20 (unchanged)
  PRO: 14 (was 15; backbone H dropped — Pro is a secondary amine)
  SER: 11 (unchanged)
  THR: 14 (unchanged)
  TRP: 24 (unchanged)
  TYR: 21 (unchanged)
  VAL: 16 (unchanged)

Variant tables (kHisAtoms_HID, kAspAtoms_ASH, kGluAtoms_GLH, etc.)
retain CCD's full atom set; the dropped atoms appear there.

### H.11 No outstanding follow-ups from the 2026-05-05 OpenAI evaluation

Six fixes from the OpenAI evaluation have all landed:

1. Cap-table N/C/O semantic deltas (Fix 1) — `BuildCapTables` now
   emits backbone overrides; cap counts 4 / 3 / 3 / 4.
2. Default residue formal charges (Fix 2) — ARG NE=+1; ASP OD2=-1;
   GLU OE2=-1; HIS ND1=0 (HIS=HIE per AMBER convention).
3. LookupBy fail-fast on invalid indices (Fix 3) — explicit
   `case kBaseVariantIdx:` for chain table; `default: return nullptr`.
4. Chain-atom whitelist filter (Fix 4) — this section.
5. Identity uniqueness docstring relaxation (Fix 5) — docstring on
   `AtomMechanicalIdentity` in `src/SemanticEnums.h` acknowledges
   chemically-equivalent H collisions.
6. Override-composition rule documentation (Fix 6) — see §H.5 above
   and the `LookupCap` comment block in the generated `.cpp`.

Tests: 352/352 pass; StringBarrier 5/5; generator builds zero
warnings. Runtime composition (`Protein::FinalizeConstruction`
applying both `LookupBy` and `LookupCap`) is the next-session
deliverable; the substrate data is ready for it.

---

## How this file should be used

- **Before encoding** (running `tools/topology/build_semantic_tables`):
  1. Confirm Section A vocabulary patches have landed in `src/SemanticEnums.h`.
  2. Apply Section E reference-document corrections.
  3. Walk Section B authorities and verify each value the encoder reads
     against the live source.
  4. Carry Section C reading-conventions into any consumer-side notes.
  5. Encoder builds; runs; emits `src/generated/LegacyAmberSemanticTables.cpp`.

- **Before consumer code** is written or modified:
  1. Read Section C carefully — these are the substrate's documented
     conventions; consumer code must respect them.
  2. Section D protects locked decisions.

- **When this file is updated**: append-only. New dependencies move to
  the bottom of their relevant section, dated. Old entries get
  resolved-by line if the dependency was concluded (e.g. an ARN custom
  template lands).
