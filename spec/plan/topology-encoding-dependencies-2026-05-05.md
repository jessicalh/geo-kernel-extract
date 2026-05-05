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
