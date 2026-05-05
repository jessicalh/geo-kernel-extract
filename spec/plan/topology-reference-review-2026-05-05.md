# Critical review — topology residue reference (2026-05-05)

## Reviewer

LLM agent (Claude Opus 4.7) executing chemistry signoff under empowerment from the project lead. Verdicts below are decisive; the reference document and its downstream C++ encoding should be updated per the actions enumerated. Authority basis: Markley 1998 IUPAC nomenclature recommendations (primary, local references-text/), BMRB atom_nom.tbl (secondary cross-check, fetched live), AMBER ff14SB params (data/ff14sb_params.dat, for variant atom inventory), Cahn-Ingold-Prelog rules for E/Z / R/S derivation.

## Verdict summary table

| # | Item | Verdict | Specific action |
|---|------|---------|-----------------|
| 1 | Charge-placement convention (carboxylate, guanidinium, ammonium) | ACCEPT | Lewis-localised int8 is correct for typed substrate |
| 2 | ARN encoding (HE removed) | ACCEPT-WITH-FLAG | Defensible; no canonical ff14SB ARN; keep regenerate-from-rtp clause |
| 3 | LYN as AmmoniumNH (no AmineNH bucket) | MODIFY | Add `PolarHKind::AmineNH`; LYN HZ2/HZ3 use it (and BMRB neutral-N H1/H2) |
| 4 | TYM Oη as AromaticHydroxyl | MODIFY | Add `PlanarGroupKind::AromaticOxide`; TYM Oη uses it |
| 5 | Trp 6-ring perimeter labels | ACCEPT | Ortho1/Ortho2/Meta1/Meta2 synthesis is internally consistent |
| 6 | HIP +1 on Nε2 specifically | ACCEPT | Defensible Lewis localisation; document convention as in §2 |
| 7 | Asn/Gln HD21/HE21 = Z assignment | REJECT | HD21=E, HD22=Z; HE21=E, HE22=Z (BMRB-consistent CIP derivation) |
| 8 | Arg Hη two-level branch + Z/E assignments | ACCEPT | HH11=Z, HH12=E, HH21=Z, HH22=E correct per Markley + BMRB + CIP |
| 9 | Glycine HA2=ProR, HA3=ProS | ACCEPT | Markley Fig 1 + BMRB atom_nom.tbl agree |
| 10 | Branched heavy atoms prochirality (Val/Leu/Ile) | ACCEPT | Val Cγ1=R, Leu Cδ1=R; Ile Cγ1 methylene Hγ13=R, Hγ12=S |
| 11a | Add `PolarHKind::AmineNH` | ADD | Required by LYN, NTERM_NEUTRAL |
| 11b | Add `PlanarGroupKind::AromaticOxide` | ADD | Required by TYM |
| 11c | Pseudoatom letter for terminal H1/H2/H3 | MODIFY | Use `PseudoatomKind::Q` with `Locant::None`; document convention |
| 11d | `Locant::None` for Pro Cα-in-ring | ACCEPT | Cosmetic; current behavior consistent |
| 12 | C++ `SynthesisedForPhe` correctness | ACCEPT-WITH-FLAG | Logic correct for CCD free-PHE entry; spurious "HG" check is dead code |

## Per-item details

### 1. Charge-placement convention for delocalised groups

**Verdict: ACCEPT**

Lewis-localised int8 placement is the correct choice for a typed-substrate where `formal_charge` is per-atom int8_t. The substrate has to commit to integer values; fractional placement (-0.5 / +0.333) cannot be expressed in int8_t and would require a wider field with no corresponding gain in calculator-side usefulness. Calculators that need delocalised/partial-charge information already get that from `LegacyAmberInvariants` (the AMBER per-atom partial charges from ff14SB) — that is the right home for the physically-realistic delocalised charge distribution. The substrate's `formal_charge` is the bookkeeping label, not the physical charge; calculators never compute electrostatics from it.

Conventions chosen (Asp/Glu OD2/OE2, CTERM O''; Arg Nε; HIP Nε2; Lys Nζ; CYM Sγ; TYM Oη) are all standard Lewis assignments and are documented in §2 of the reference. Calculators that stratify by formal_charge==-1 will pick up the carboxylate oxygen carrying the label; this is the intended semantic. Symmetrisation in MD/DFT is a property of the *physical* charge distribution from `LegacyAmberInvariants`, not the substrate label.

### 2. ARN encoding (HE removed)

**Verdict: ACCEPT-WITH-FLAG**

ARN is not in `data/ff14sb_params.dat` (line 7: "Terminal ASH/CYM/GLH/LYN are absent in this ff14SB source and must fail explicitly"; ARN absent entirely from the file). AMBER ff14SB does not ship a standard ARN INTERNAL template; ARN is a non-standard variant requiring user-supplied parameters. There is therefore no ff14SB-canonical answer to "which H is removed" — the reference's choice of removing HE (giving the neutral guanidine with the lone pair on Nε) is chemistry-defensible: it preserves the four guanidinium NH's in symmetric Z/E pattern and places the lone pair on the most-CIP-prominent N (the one bonded to Cδ). The chemistry-formal alternative (remove one of HH11/HH12/HH21/HH22) would break the guanidinium's Z/E symmetry asymmetrically.

The reference correctly flags this with "encoder MUST regenerate from `AmberAminoAcidVariantTable["ARN"]` at code-gen time" — but `AmberAminoAcidVariantTable` (in `src/AminoAcidType.cpp` line 35) only contains a name/description/formal_charge tuple, not an atom-removal list. The encoder will need to either (a) keep the HE-removal default and document, (b) refuse ARN in production until ff14SB ships a canonical template, or (c) use a propka-style on-the-fly inference. Status quo encoding is fine; flag stays.

### 3. LYN as AmmoniumNH (no AmineNH bucket)

**Verdict: MODIFY**

Add `PolarHKind::AmineNH` to SemanticEnums.h. The chemistry distinction is real and load-bearing for shielding statistics: a neutral primary amine N-H (LYN -NH2; NTERM_NEUTRAL -NH2) has a pKa ~9-10, exchanges slower than ammonium, and exhibits a chemical shift environment ~1-2 ppm distinct from the charged ammonium. Wuethrich 1986 ch. 2 distinguishes them; SHIFTX2 trains separate models on neutral-amine vs ammonium environments. Calculators stratifying shielding by polar-H-environment will incorrectly bucket LYN with LYS HZ if these collapse onto the same enum value.

The change is small (one new enum value, one docstring). It also makes the convention internally consistent — `BackboneAmide` and `SidechainPrimaryAmide` already distinguish chemically-distinct N-H environments, so distinguishing AmineNH from AmmoniumNH follows the same pattern. After adding: LYN HZ2/HZ3 → AmineNH; LYS HZ1/HZ2/HZ3 → AmmoniumNH; NTERM_CHARGED H1/H2/H3 → AmmoniumNH; NTERM_NEUTRAL H1/H2 → AmineNH.

**Caveat to fix in the same pass**: the reference's LYN row says "HZ3: REMOVED (Lys neutral has -NH2 with two ζ-Hs; AMBER convention removes Hζ3 to leave HZ1 and HZ2)." This is **wrong per `data/ff14sb_params.dat` lines 449-451**: AMBER ff14SB LYN has **HZ2 and HZ3** (no HZ1). The reference must correct this to "HZ1: REMOVED". Same logic applies to NTERM_NEUTRAL — the H atom names follow ff14SB's NTRESIDUE convention.

### 4. TYM Oη as AromaticHydroxyl (no H present)

**Verdict: MODIFY**

Add `PlanarGroupKind::AromaticOxide` to SemanticEnums.h. The post-deprotonation Tyr -O- is a fundamentally distinct chemical group from -OH: the lone pair on phenolate O conjugates strongly with the ring π system, redistributes the ring current, and shifts the para-Cζ chemical environment by 5-10 ppm relative to neutral Tyr. Calling Oη `AromaticHydroxyl` when there is no hydroxyl proton is semantically wrong — the enum docstring on `AromaticHydroxyl` literally states "C-zeta-O-H rotation is a real conformation degree of freedom captured by PlanarGeometryResult", which doesn't apply when the H doesn't exist.

After adding: TYM Oη → AromaticOxide; Tyr (default) Oη → AromaticHydroxyl. The two-value distinction is small but clean and is consistent with the project's "first-class typed semantic, not a string-equality test" philosophy.

### 5. Trp 6-ring perimeter labels (Ortho1=Cε3, Ortho2=Cζ2, Meta1=Cζ3, Meta2=Cη2)

**Verdict: ACCEPT**

The labeling is internally consistent (BMRB atom_nom.tbl shows Trp atoms without positional descriptors, so this synthesis goes beyond BMRB; Markley does not publish ipso/ortho/meta/para for indole 6-ring). The reference's choice — atoms adjacent to bridgeheads as "ortho", atoms two-bonds-from-bridgehead as "meta" — is the standard organic-chemistry intuition for fused-ring perimeters. Symmetry is preserved (Ortho1/Cε3 paired with bridgehead Cδ2 closer to Cα; Ortho2/Cζ2 paired with bridgehead Cε2 closer to Nε1).

The alternative encodings considered (PyrroleAlpha/PyrroleBeta-style; raw Markley atom names alone) are weaker: PyrroleAlpha/Beta are 5-ring concepts that don't transfer cleanly to a 6-ring perimeter; raw atom names defeat the point of having typed RingPositionLabel. Document the synthesis explicitly in the SemanticEnums.h docstring on `RingPositionLabel::Ortho1` and add a Trp-specific note that "Ortho1/Ortho2 in Trp 6-ring refer to perimeter positions adjacent to bridgehead atoms Cδ2 and Cε2 respectively; this extends Vollhardt benzene convention to fused-ring nomenclature".

### 6. HIP +1 on Nε2 specifically

**Verdict: ACCEPT**

The chemistry rationale: HIP arises from HIE (AMBER's default) by protonating Nδ1 — i.e. the "extra proton" that introduces the +1 charge is the new Hδ1, which sits on Nδ1. Lewis-localising the +1 on the second-protonated N is therefore *the* canonical interpretation; placing it on Nδ1 is also defensible; placing it on Nε2 (the one that was already protonated in HIE) is somewhat arbitrary. The reference picks Nε2 and documents the choice; either Nδ1 or Nε2 is fine, and the documentation is what matters.

For NMR-shielding stratification purposes, the choice does not matter: HIP imidazolium delocalises across the ring, and the per-atom AMBER partial charges in `LegacyAmberInvariants` (HIP HD1=0.3866, HE2=0.3911 per `data/ff14sb_params.dat` lines 442/447) are already symmetric across the two NHs. The substrate label is bookkeeping. Minor docstring suggestion: in §2 of the residue reference, change the rationale from "second-protonated nitrogen" to "the canonical AMBER ff14SB HIP template structure protonates Nδ1 to convert HIE → HIP, so we localise the +1 on the same Nε2 position used for the HIE neutral state to keep ring atom labels consistent across the HIE/HIP variant transition." That is the principled reason.

### 7. Asn / Gln amide H E/Z

**Verdict: REJECT**

The reference encodes Hδ21=Z, Hδ22=E (Asn) and Hε21=Z, Hε22=E (Gln) on the basis of "1 = cis = Z, 2 = trans = E". This derivation is wrong. The Markley caption says "1 or 2 by relationship (cis or trans, respectively) to the heavy atom three bonds closer to the main chain (Cβ for Asn, Cγ for Gln, Nε for Arg)" — i.e. "1 = cis to Cβ" for Asn, "1 = cis to Cγ" for Gln. But "cis to the low-priority substituent" does not equal Z. Z (zusammen) requires the **higher-priority** substituents on each end of the planar bond to be cis.

Working through Asn: planar amide CG=ND2 (resonance double bond). On CG, substituents are OD1 (priority O) and CB (priority C-bonded). OD1 > CB. On ND2, two H's labeled 1/2 by Markley convention. Hδ21 cis to Cβ means Hδ21 on the SAME side as Cβ on the C=N bond (i.e. opposite the OD1 side, since OD1 and Cβ are on opposite sides of CG in the planar amide). Therefore Hδ21 is *trans* to OD1 (the higher-priority CG substituent) → Hδ21 = E. Hδ22 = Z.

BMRB atom_nom.tbl confirms this: ASN HD21 = E, HD22 = Z; GLN HE21 = E, HE22 = Z; ARG HH11 = Z, HH12 = E, HH21 = Z, HH22 = E. The Arg encoding in the reference is correct (because for Arg, Nε is BOTH the "heavy atom 3 bonds in" AND the higher-priority substituent on Cζ — so 1=cis-to-Nε does coincide with 1=cis-to-high-priority=Z). Only Asn and Gln are wrong.

**Action**: in the residue reference, Asn HD21 row → planar_stereo=E (was Z); HD22 row → planar_stereo=Z (was E). Gln HE21 row → planar_stereo=E (was Z); HE22 row → planar_stereo=Z (was E). Update §2 convention note: "For Asn/Gln, the '1' index labels the H cis to the heavy atom 3 bonds closer to main chain (Cβ for Asn, Cγ for Gln); BECAUSE that heavy atom is the LOWER-priority substituent on the planar carbon (the carbonyl O is higher), 'cis to Cβ/Cγ' means 'trans to the high-priority O' which is E by CIP. For Arg, '1' = cis to Nε which IS the higher-priority substituent on Cζ (because Nε is bonded to Cδ, while the other Nη is bonded only to H atoms), so '1=Z' there. The Z/E mapping is residue-specific in the way that matters for the substrate label." Update disagreement-log entries 1, 2: "Resolution: BMRB atom_nom.tbl is canonical and the reference's earlier claim 'caption wins → 1=Z' was a conflation of 'cis to Cβ' with 'cis to high priority'."

### 8. Arg Hη two-level branch and PlanarStereo

**Verdict: ACCEPT**

Per BMRB atom_nom.tbl confirmed live: HH11 = Z, HH12 = E, HH21 = Z, HH22 = E. Per Markley caption (text-2:357-372): outer index = which Nη (Nη1 cis to Cδ → 1; Nη2 trans → 2); inner index = cis (1) or trans (2) to Nε (heavy atom 3 bonds in). Per CIP on Cζ=Nη bond: Nε > the-other-Nη (because Nε is bonded to Cδ, which leads to the rest of side chain; the other Nη is bonded to H's only). So "cis to Nε" = "cis to high priority" = Z. Therefore inner=1 → Z, inner=2 → E.

The reference correctly encodes:
- Hη11: outer=1 (Nη1), inner=1 → Z. ✓
- Hη12: outer=1, inner=2 → E. ✓
- Hη21: outer=2 (Nη2), inner=1 → Z. ✓
- Hη22: outer=2, inner=2 → E. ✓

The two-level BranchAddress encoding (outer=Nη, inner=H-within-Nη) is exactly the structural shape required by Markley caption + the existing `BranchAddress` struct in SemanticEnums.h lines 97-103 (which already documents the Arg case explicitly: "For 99% of atoms only `outer` is non-zero; for Arg side-chain Hs both indices apply"). Encoding lines up cleanly with substrate vocabulary.

### 9. Glycine HA inversion: HA2 = ProR, HA3 = ProS

**Verdict: ACCEPT**

BMRB atom_nom.tbl confirms: HA2 = pro-R, HA3 = pro-S. Markley Fig 1 (text-2 OCR'd region) marks HA2 as (R). Reference encodes HA2=ProR, HA3=ProS, correct. Note the BMRB SC column does say "pro-R/pro-S" explicitly here, distinct from the E/Z convention used for amide H's — Markley's Glycine treatment IS R/S because Hα atoms are tetrahedral CHHC (with N + C as the other substituents), not on a planar centre.

The "INVERTED vs typical alternation" comment in the reference correctly captures that Gly is the only standard residue where DiastereotopicIndex 2 maps to ProR (everywhere else, 3=R, 2=S). The inversion happens because Cα in Gly has substituents N, H, H, C; reversed CIP priority order vs the standard L-amino acid Cα with N, H, side-chain-C, C(carbonyl).

### 10. Branched heavy atoms prochirality (Val / Leu / Ile)

**Verdict: ACCEPT**

Per Markley Fig 1 + research dossier curated table (lines 298-340 of topology-fields-research-2026-05-05.md):
- Val Cγ1=ProR, Cγ2=ProS. ✓ (reference matches)
- Leu Cδ1=ProR, Cδ2=ProS. ✓
- Ile Cγ1 methylene Hγ13=ProR, Hγ12=ProS. ✓
- Ile Cγ2 (γ-methyl): per Markley convention this is the "ProS" branch — Ile has a γ-methylene (Cγ1) AND a γ-methyl (Cγ2). The diastereotopy is between the methylene (HG12/HG13) at Cγ1 and the methyl-bearing Cγ2. The "ProS"-on-Cγ2 label is the chain-branch prochirality, not Hγ-prochirality. The reference encodes Cγ2=ProS, which is correct under the Markley Cγ1>Cγ2 priority numbering (Cγ1 has H+H+CD1 going further; Cγ2 is a terminal methyl with H+H+H — Cγ1 wins by CIP, so Cγ1=ProR, Cγ2=ProS).
- Ile Cδ1=NotProchiral. ✓ (only one δ-carbon; nothing to be prochiral with)

The methyl Hs inheriting their parent C's prochiral label (e.g., Val HG11/HG12/HG13 all on Cγ1=ProR) is the standard convention and the reference encodes this correctly. The methyl Hs themselves are NOT prochiral within the methyl group (they are equivalent under rotation), so their `prochiral` field is `NotProchiral`; the parent C carries the ProR label.

### 11a. PolarHKind::AmineNH — ADD

**Verdict: ADD**

Required by item 3. Add to SemanticEnums.h line ~358 with docstring naming LYN, NTERM_NEUTRAL as occupants. Cite Wuethrich 1986 ch. 2 + SHIFTX2 separate-model evidence.

### 11b. PlanarGroupKind::AromaticOxide — ADD

**Verdict: ADD**

Required by item 4. Add to SemanticEnums.h line ~232 with docstring naming TYM as occupant; describe the conjugation distinction vs AromaticHydroxyl.

### 11c. Pseudoatom letter for terminal H1/H2/H3 — MODIFY

**Verdict: MODIFY (use Q/None/0/false)**

Markley Table 1 doesn't define a pseudoatom letter for backbone-amide H or for N-terminal H1/H2/H3. The reference's current shorthand "Q/None/0/false" is the right answer: H1/H2/H3 in NTERM_CHARGED form an equivalent-H group (3 H's on a single N, rapidly exchanging in solution NMR), which is exactly what Q is for ("group of equivalent Hs that aren't methyl"). Locant=None correctly captures the absence of a Greek-letter side-chain position. Don't introduce a new pseudoatom letter ("T") — that would diverge from the IUPAC M/Q/R taxonomy without chemistry justification.

Document the convention in SemanticEnums.h::PseudoatomKind comment block: "N-terminal H1/H2/H3 are encoded as Q-with-Locant=None, treating them as a 3-H equivalent group. Markley Table 1 does not list a pseudoatom for these terminus-specific Hs; the Q taxonomy is the closest fit and is used here without a new letter."

### 11d. Locant::None clarification for Pro Cα-in-ring — ACCEPT

**Verdict: ACCEPT (cosmetic only; document but don't change)**

Pro Cα is in the pyrrolidine ring (per the Pro entry: "ring atoms are N, Cβ, Cγ, Cδ (and Cα for completeness)"). Markley convention reserves Locant for side-chain atoms only; Cα and other backbone atoms have Locant::None. This is consistent with the docstring on `Locant::None` ("Backbone or terminal atom; no Greek-letter locant"). For Pro specifically, Cα happens to also be a ring member, but that's a RingPosition fact, not a Locant fact. The reference's encoding (Cα with Locant=None and ring_primary=Pyrrolidine_Pro/Saturated/5/f/1) correctly separates the two orthogonal facts.

A clarifying comment in SemanticEnums.h on `RingPosition` documenting "An atom can be in a ring AND have Locant::None — the two fields are orthogonal: Locant tells you the backbone-vs-side-chain Greek-letter position; RingPosition tells you ring membership. Pro Cα is in the pyrrolidine ring with Locant=None." would help. Not a behavioural change.

### 12. Did the existing C++ `SynthesisedForPhe` itself encode PHE correctly?

**Verdict: ACCEPT-WITH-FLAG**

Logic is correct for the CCD free-PHE entry. Ring tagging (CG=Ipso, CD1/CD2=Ortho1/Ortho2, CE1/CE2=Meta1/Meta2, CZ=Para) matches both BMRB convention and the residue reference. PolarH classification (H/HN → BackboneAmide; n_terminus → AmmoniumNH; c_terminus → CarboxylOH) handles the free-PHE CCD entry correctly. PeptideAmide assignment to N/C/O/H + Carboxylate fall-through for C/O/OXT (only when planar_group not already set) preserves the right ordering: peptide-context PHE has C/O as PeptideAmide (set first), free-PHE adds Carboxylate to OXT only.

Pseudoatom assignments (HB2/HB3 → Q/Beta/0/false; HD1/HD2 → Q/Delta/0/super=true; HE1/HE2 → Q/Epsilon/0/super=true; HZ → R/Zeta/0/super=true) match the reference exactly.

**Dead-code flag**: line 794 reads `if (atom_id == "CG" || atom_id == "HG") set_ring(...)`. PHE has no atom named "HG" (CG is the ipso carbon, attached to CB; CG has no H because it's substituted). The "HG" check is harmless dead code that suggests the function was templated against a generic ring-residue pattern (like Tyr where HG might exist as a typo for HH the para-OH proton, though that's also not standard). Either remove the "|| atom_id == \"HG\"" clauses (lines 794-799) or document that they're forward-compatibility scaffolding. Recommendation: remove them, since silent acceptance of a non-existent atom is the kind of thing that breaks future TYR/TRP generalisation.

**Function correctness conclusion**: PHE encoding from `SynthesisedForPhe` is chemically correct and matches the residue reference's PHE rows. The C++ acceptance gate works.

## Cross-cutting observations

**The Z/E bug for Asn/Gln is the most important single finding.** The reference document conflates two different "cis" relationships: (a) cis to the heavy atom 3 bonds closer to the main chain (Markley's labeling rule) and (b) cis to the higher-CIP-priority substituent on the planar partner (CIP's E/Z rule). For Arg these coincide; for Asn/Gln they don't. The reference's blanket "1 = cis = Z" generalisation is wrong for Asn/Gln. BMRB atom_nom.tbl is the canonical cross-check that catches this. Disagreement-log entries 1 and 2 in the reference need to be inverted: caption gives the *labeling rule* (1=cis-to-Cβ), but the *Z/E mapping* requires CIP priority on each end of the planar bond.

**Vocabulary patches (AmineNH, AromaticOxide) are small and high-value.** Both are chemistry-real distinctions (neutral amine vs ammonium; phenol vs phenolate) that downstream calculators can stratify on. Cost: two new enum values, two docstring entries, regenerating `LegacyAmberSemanticTables.cpp` for LYN/NTERM_NEUTRAL/TYM rows. Benefit: substrate stops conflating chemically-distinct environments under the same enum value, eliminating a sliding-bug class.

**The LYN HZ atom-name typo should be fixed in the same pass as the AmineNH addition.** AMBER ff14SB LYN has HZ2 + HZ3 (no HZ1). The reference says "HZ3 removed, leaving HZ1 + HZ2" — wrong. Fix to "HZ1 removed, leaving HZ2 + HZ3" per `data/ff14sb_params.dat` lines 449-451.

**The reference is otherwise high-quality.** Out of 12 substantive items, 7 are pure ACCEPTs, 1 is ACCEPT-WITH-FLAG (ARN, defensibly so), 1 is REJECT (the Z/E bug), 4 are MODIFY (LYN HZ name, AmineNH addition, AromaticOxide addition, terminal-H pseudoatom doc). The decision tree, citation discipline, and self-audit checklist are well-constructed. The disagreement log (Section 5) is honest and complete except for the Z/E inversion error.

## Final encoding action list

1. **Residue reference document** (`spec/plan/topology-residue-reference-2026-05-05.md`):
   - Asn HD21 row: planar_stereo Z → E.
   - Asn HD22 row: planar_stereo E → Z.
   - Gln HE21 row: planar_stereo Z → E.
   - Gln HE22 row: planar_stereo E → Z.
   - LYN variant block: change "HZ3 REMOVED, leaving HZ1+HZ2" → "HZ1 REMOVED, leaving HZ2+HZ3" (per ff14SB).
   - LYN HZ2/HZ3 rows: polarH AmmoniumNH → AmineNH (after enum value added).
   - NTERM_NEUTRAL H1/H2 rows: polarH AmmoniumNH → AmineNH.
   - TYM Oη row: planar_group AromaticHydroxyl → AromaticOxide.
   - §2 charge-distribution conventions block: rewrite the 1=Z/2=E sub-clause for Asn/Gln; clarify that BMRB E/Z labels and CIP-derived Z/E from Markley's labeling rule require knowing the priority of the heavy-atom-3-bonds-in vs the OTHER substituent on the planar carbon. Arg is the special case where 1=Z; Asn/Gln have 1=E.
   - Disagreement log entries 1, 2: invert resolution. Caption gives the labeling rule; CIP gives the Z/E. BMRB atom_nom.tbl is the canonical Z/E source. Acknowledge the residue-specific Z/E mapping.
   - Disagreement log entry 8: tighten — the AmineNH addition is now ACTION, not just a flag.
   - Disagreement log entry 9: tighten — the AromaticOxide addition is now ACTION.

2. **SemanticEnums.h** (`src/SemanticEnums.h`):
   - Add `PolarHKind::AmineNH = 12` (or insert before OtherPolarH and renumber). Docstring: "Neutral primary-amine N-H (LYN HZ2/HZ3; NTERM_NEUTRAL H1/H2). pKa ~9-10; slower exchange than ammonium and ~1-2 ppm chemical-shift offset. Wuethrich 1986 ch. 2; SHIFTX2 (Han 2011) trains separate models for neutral-amine vs ammonium environments."
   - Add `PlanarGroupKind::AromaticOxide = 9`. Docstring: "Aryloxide (TYM Oη after deprotonation). The phenolate -O- conjugates strongly with the ring π system, redistributing ring current and shifting the para-Cζ environment by 5-10 ppm relative to neutral Tyr. Distinguished from `AromaticHydroxyl` because the H is absent and the conjugation is qualitatively different."
   - Add comment on `RingPosition` struct documenting Pro Cα as in-ring with Locant::None.
   - Add comment on `RingPositionLabel::Ortho1`/`Ortho2` documenting the Trp 6-ring extension (perimeter atoms adjacent to bridgehead = Ortho).
   - Add comment on `PseudoatomKind` documenting the Q-with-Locant::None convention for terminal H1/H2/H3.

3. **C++ generator function** (`tools/topology/build_semantic_tables.cpp::SynthesisedForPhe`):
   - Remove `|| atom_id == "HG"` from line 794 (PHE has no HG atom; dead code with risk of silently mis-handling future TYR generalisations).
   - Otherwise unchanged. PHE encoding is correct.

4. **PHE acceptance gate**: re-run after removing the dead "HG" check; the PHE rows in the residue reference still match.

## Anything you couldn't verdict

**ARN atom-removal canonicalisation**: the choice of which Hη/HE to remove is documented as a flag rather than verdicted to a specific final answer because there is no AMBER ff14SB canonical template (ARN absent from `data/ff14sb_params.dat`). The reference's choice (remove HE) is chemistry-defensible and is accepted with the existing flag in place. If the project later adopts a propka-style on-the-fly inference or a custom ARN template, the reference row should regenerate from that authority. This is not an unresolved verdict — it's an explicitly-flagged degree of freedom in the substrate that the encoder must commit to with its own rule when it generates the LegacyAmberSemanticTables.cpp.
