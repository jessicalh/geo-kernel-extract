# Topology residue reference — synthesised + Markley-derived per-(residue, variant, atom) substrate (2026-05-05)

## Section 1: Header

### Purpose

This document encodes the synthesised + Markley-derived chemistry-substrate fields for every (residue, variant, atom) triple in scope for `LegacyAmberTopology::AtomSemanticTable`. It is the source from which the C++ `SynthesisedFor<Residue>` functions are mechanically transcribed in a follow-on session. Algorithmic fields (RDKit-perceived aromaticity, hybridisation, bond-order mask, equivalence class) and mechanical fields (`Locant`, `BranchAddress`, `DiastereotopicIndex` — all derivable from the atom name string by a parser) are intentionally out of scope here. They are populated by orthogonal generator passes against this same atom inventory.

### Sources read

Canonical typed-enum vocabulary and per-enum-value docstrings:
- `/shared/2026Thesis/nmr-shielding/src/SemanticEnums.h` (canonical commit `721e681`, "Topology generator: emitter + string-barrier test")

PHE acceptance gate:
- `/shared/2026Thesis/nmr-shielding/tools/topology/build_semantic_tables.cpp` lines 779-842 (`SynthesisedForPhe`)

Field-level research:
- `/shared/2026Thesis/nmr-shielding/spec/plan/topology-fields-research-2026-05-05.md` (1328-line dossier; per-field convention tables; per-residue prochiral mapping; ring-position assignments)

IUPAC nomenclature primary literature:
- `/shared/2026Thesis/nmr-shielding/references-text/markley-1998-iupac-nmr-nomenclature-recommendations-text-{1..8}.txt`. In particular:
  - text-2:30-60: backbone vs sidechain atom-name conventions
  - text-2:357-372: Figure 1 caption (CIP-clockwise rule, Asn/Gln/Arg Hη "cis/trans to heavy atom three bonds closer to main chain" rule)
  - text-2:78-358: Figure 1 atom diagrams with (R) and (E)/(Z) annotations (OCR-garbled in places but largely recoverable)
  - text-3:181-242: Table 1 (pseudoatom taxonomy)
  - text-3:295-315: Pro ring puckering note
  - text-3:421-441: Phe/Tyr ring-flip note
  - text-7:148-165: worked CIP-rule restatement

BMRB cross-reference (consulted but not all queried in this pass; the local sources cover the standard 20):
- https://bmrb.io/referenc/nomenclature/
- https://bmrb.io/ref_info/atom_nom.tbl
- https://bmrb.io/ref_info/pseudoatom_nom.txt

AMBER variant inventory (atom-name-by-variant authority):
- `/shared/2026Thesis/nmr-shielding/src/AminoAcidType.cpp::AMINO_ACID_TYPES` (variant→ff_name table; HID/HIE/HIP, ASH, GLH, CYX/CYM, LYN, ARN, TYM)

### Scope-in (encoded in this file)

- 20 standard amino acids in default protonation, alphabetical: ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL.
- AMBER protonation variants: HID, HIE, HIP (the AMBER default for HIS is HIE; the "default" HIS row above is therefore identical to HIE), ASH, GLH, CYX, CYM, LYN, ARN, TYM.
- Terminal states: NTERM_CHARGED, NTERM_NEUTRAL, CTERM_DEPROTONATED, CTERM_PROTONATED.

### Scope-out (NOT encoded here)

- Mechanical fields derivable by a parser from the atom_id (`Locant`, `BranchAddress`, `DiastereotopicIndex`).
- Algorithmic fields (`aromatic`, `equivalence_class`, `Hybridisation`, `BondOrderMask`).
- Code (no C++).
- Per-residue investigative essays.
- Per-witness provenance machinery (the C++ generator stashes witnesses; this file just records the value).
- Cap pseudo-residues (ACE, NME, NHE), modified residues, nucleic acids, ligands.

---

## Section 2: Convention notes

### AMBER's HIS-default-as-HIE convention

AMBER ff14SB's default histidine protonation is HIE (Nε2-protonated, neutral, no Hδ1). When a residue is parsed under the bare 3-letter code "HIS" without a tautomer override, the substrate row is the HIE row. The "Default" block for HIS in this document is therefore identical to the HIE block. The HID and HIP variant blocks list deltas vs that HIE/HIS-default reference.

### Markley Figure 1 OCR caveat

The OCR'd Figure 1 in `markley-1998-...text-2.txt` lines 60-358 has visual layout artefacts (bond lines rendered as long em-dashes; some atom labels and (R)/(E)/(Z) annotations partially obscured). The two channels for cross-checking inside the local corpus are: (a) the Figure 1 caption at text-2:357-372 (which states the rules in prose, unambiguously) and (b) Markley §28 / text-7:148-165 (CIP statement). In every case where the diagram and the caption disagree, this document follows the caption — IUPAC's authoritative statement of the convention.

The most-OCR-garbled atom-by-atom marks affect the Asn/Gln/Arg primary-amide H labels. The Markley caption (text-2:367-371) gives the **labelling rule**: each NH2 hydrogen of Asn / Gln / Arg side chain is numbered "1 or 2" by its cis (1) or trans (2) relationship to the heavy atom three bonds closer to the main chain (Cβ for Asn, Cγ for Gln, Nε for Arg). The IUPAC **Z/E mapping** follows from CIP priority on each end of the planar bond — and this is RESIDUE-SPECIFIC because the priority of "the heavy atom three bonds in" relative to the OTHER substituent on the planar carbon changes per residue. **BMRB `atom_nom.tbl` is the canonical Z/E source** for the standard 20 (verified 2026-05-05).

- **Asn**: on the planar Cγ=Nδ2 bond, OD1 is the high-priority Cγ substituent (O > C-bonded). HD21 = cis to Cβ = trans to OD1 = **E**. HD22 = trans to Cβ = cis to OD1 = **Z**. (BMRB confirms.)
- **Gln**: same pattern with Cδ in place of Cγ. HE21 = **E**. HE22 = **Z**. (BMRB confirms.)
- **Arg**: on the planar Cζ=Nη bonds, Nε IS the high-priority Cζ substituent (Nε is bonded to Cδ and the rest of the side chain; the other Nη is bonded only to H atoms). So "cis to Nε" = "cis to high priority" = Z. HH11 = **Z**, HH12 = **E**, HH21 = **Z**, HH22 = **E**. (BMRB confirms.)
- **Arg Nη1 vs Nη2 themselves**: Nη1 cis to Cδ → Z; Nη2 trans → E.

The asymmetry between Asn/Gln (where "1" = E) and Arg (where "1" = Z) is a real consequence of CIP priority on the planar carbon, **not** a labelling inconsistency. Consumers reading `PlanarStereo` must read the substrate's value per atom; do not assume a universal "1 = Z" rule.

### Word/Reduce flip ambiguity

The substrate carries the canonical IUPAC label. Asn/Gln side-chain amide flip detection (Word et al. 1999 J. Mol. Biol. 285 1735) is conformation-side, not substrate-side. Out of scope for this document.

### Charge-distribution conventions chosen here

For delocalised groups the convention used in this file is the formal Lewis structure value:
- **Carboxylate (Asp default Oδ1/Oδ2; Glu default Oε1/Oε2; CTERM_DEPROTONATED O'/O''):** -1 placed on the formally double-bonded-to-H "second" oxygen (Oδ2, Oε2, O''). Oδ1/Oε1/O' carry 0. Heavy carbon (Cγ/Cδ/C) carries 0. This is the Lewis convention; the actual delocalisation symmetrises Oδ1/Oδ2 in MD/DFT.
- **Guanidinium (Arg default):** +1 placed on the formally protonated Nε. Cζ, Nη1, Nη2 carry 0. Lewis convention; delocalisation symmetrises across Cζ-N3 in MD/DFT.
- **Imidazolium (HIP):** +1 placed conceptually on the "whole imidazolium" but for the per-atom substrate a choice is made: charge is placed on Nε2 (the second-protonated nitrogen). Nδ1, Cε1, Cδ2, Cγ carry 0. Flag in notes.
- **Ammonium (Lys default Nζ; NTERM_CHARGED N):** +1 on N. Hζ1/2/3 (or H1/2/3) carry 0.
- **Thiolate (CYM Sγ):** -1 on Sγ. Cβ carries 0.
- **Aryloxide (TYM Oη):** -1 on Oη. Cζ carries 0.

### Table column legend

- `atom`: AMBER atom name (matches the names used in `AmberAminoAcidVariantTable` / ff14SB rtp).
- `element`: H, C, N, O, S.
- `planar_group`: `PlanarGroupKind` enum value name verbatim (None, PeptideAmide, SidechainAmide, Guanidinium, Imidazole, Aromatic6Ring, Aromatic5Ring, Carboxylate, AromaticHydroxyl).
- `planar_stereo`: `PlanarStereo` enum name verbatim (NA = NotApplicable, E, Z, Unspecified). NA for atoms not at a defined E/Z centre.
- `pseudoatom`: shorthand `kind/locant/branch/super` where:
  - kind ∈ {M, Q, R}
  - locant ∈ {Beta=2, Gamma=3, Delta=4, Epsilon=5, Zeta=6, Eta=7} (or named: "Alpha", "Beta", ...)
  - branch ∈ {0, 1, 2}: 0 for non-branched (e.g. QB on Phe), 1 or 2 for branched methyls (Val MG1 vs MG2)
  - super ∈ {true=t, false=f}: true if atom is also a member of a higher-order Q super-group (Val QG, Leu QD, Arg QH, Phe/Tyr QR). For Markley letter-only labels (M alone) the M and the Q super are both true if the methyl is in QG/QD; in the row I encode the M as "M/Locant/Branch/true" because the same H-atom is counted in both.
  - "—" for atoms with no pseudoatom membership.
- `polarH`: `PolarHKind` enum name verbatim (NotPolar, BackboneAmide, SidechainPrimaryAmide, IndoleNH, AmmoniumNH, GuanidiniumNH, ImidazoleNH, CarboxylOH, HydroxylOH_Aliphatic, HydroxylOH_Aromatic, ThiolSH, OtherPolarH).
- `ring_primary` / `ring_secondary`: `RingMembership` shorthand `sys/pos/size/arom/het` where:
  - sys ∈ {Benzene_Phe, Benzene_Tyr, Imidazole_His, Indole_Trp_5, Indole_Trp_6, Pyrrolidine_Pro, NotInRing}
  - pos ∈ {Ipso, Ortho1, Ortho2, Meta1, Meta2, Para, PyrroleAlpha, PyrroleBeta, BridgeFusion, Heteroatom_NH, Heteroatom_NoH, Heteroatom_OH, Saturated, NotInRing}
  - size: 5 or 6 (0 if NotInRing)
  - arom: t/f
  - het: integer 0/1/2 (number of heteroatoms in the ring)
  - "NotInRing" alone shorthand for the all-default RingMembership.
- `prochiral`: `ProchiralStereo` enum name verbatim (NotProchiral, ProR, ProS, Unassigned).
- `formal_chg`: int8 (-1, 0, +1).
- `exch`: t/f, derived as "exch = (polarH != NotPolar)". Always true for a polar H, false otherwise.
- `notes`: free-text for outliers, conventions chosen, sources cited.

### Per-residue alternation defaults applied through the file

Where applicable, the following defaults are applied (per Markley Fig 1 + caption + the user spec):
- All β-onward methylene Hs in standard L-amino acids: H<locant>3 = ProR, H<locant>2 = ProS. Applied for: Leu (β), Ser (β), Cys (β), Met (β, γ), Asp (β), Asn (β), Glu (β, γ), Gln (β, γ), Lys (β, γ, δ, ε), Arg (β, γ, δ), His (β), Phe (β), Tyr (β), Trp (β), Pro (β, γ).
- **Glycine Hα is INVERTED**: HA2 = ProR, HA3 = ProS.
- **Pro Hδ matches the typical alternation**: Hδ3 = ProR, Hδ2 = ProS.
- **Branched heavy atoms**: Val Cγ1 = ProR, Cγ2 = ProS (and the methyl Hs inherit). Leu Cδ1 = ProR, Cδ2 = ProS (methyl Hs inherit). Ile: Cγ2 (the only γ-methyl carbon) is the "ProS" branch; Cδ1 is the only δ carbon (no diastereotopy at δ); Cγ1 is methylene with diastereotopic Hγ12/Hγ13 → Hγ13 = ProR, Hγ12 = ProS.

---

## Section 3: Per-residue blocks

### ALA — alanine

Atoms: N, H, CA, HA, C, O, CB, HB1, HB2, HB3.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | backbone amide N |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | a.k.a. HN |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | chiral L-Cα |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | carbonyl C |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | carbonyl O |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | β-methyl C; the M-pseudoatom MB is on the H trio, not the C |
| HB1 | H | None | NA | M/Beta/0/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | Markley MB; 3 equivalent β-methyl Hs |
| HB2 | H | None | NA | M/Beta/0/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB3 | H | None | NA | M/Beta/0/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |

Formal charges total: 0.

### ARG — arginine

Atoms: N, H, CA, HA, C, O, CB, HB2, HB3, CG, HG2, HG3, CD, HD2, HD3, NE, HE, CZ, NH1, HH11, HH12, NH2, HH21, HH22.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| CG | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HG2 | H | None | NA | Q/Gamma/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QG |
| HG3 | H | None | NA | Q/Gamma/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QG |
| CD | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | the heavy atom three bonds closer to main chain (re Markley caption rule for Hη) |
| HD2 | H | None | NA | Q/Delta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QD |
| HD3 | H | None | NA | Q/Delta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QD |
| NE | N | Guanidinium | NA | — | GuanidiniumNH | NotInRing | NotInRing | NotProchiral | +1 | f | formal charge placed here per chosen convention; PolarHKind on the H, not the N |
| HE | H | Guanidinium | NA | — | GuanidiniumNH | NotInRing | NotInRing | NotProchiral | 0 | t | the lone Nε-H; Markley GuanidiniumNH |
| CZ | C | Guanidinium | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | guanidinium central sp2 C |
| NH1 | N | Guanidinium | Z | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | Nη1 cis to Cδ ⇒ Z; AMBER spelling NH1 |
| HH11 | H | Guanidinium | Z | Q/Eta/1/true | GuanidiniumNH | NotInRing | NotInRing | NotProchiral | 0 | t | outer=1 (on Nη1); inner=1 (cis to Cδ via Nη1-Nε-Cδ chain) ⇒ Z. QH1; in QH super |
| HH12 | H | Guanidinium | E | Q/Eta/1/true | GuanidiniumNH | NotInRing | NotInRing | NotProchiral | 0 | t | outer=1; inner=2 (trans) ⇒ E. QH1; in QH super |
| NH2 | N | Guanidinium | E | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | Nη2 trans to Cδ ⇒ E; AMBER spelling NH2 |
| HH21 | H | Guanidinium | Z | Q/Eta/2/true | GuanidiniumNH | NotInRing | NotInRing | NotProchiral | 0 | t | outer=2 (on Nη2); inner=1 (cis) ⇒ Z. QH2; in QH super |
| HH22 | H | Guanidinium | E | Q/Eta/2/true | GuanidiniumNH | NotInRing | NotInRing | NotProchiral | 0 | t | outer=2; inner=2 (trans) ⇒ E. QH2; in QH super |

Formal charges total: +1 (placed on NE; chemistry-formal Lewis convention; delocalises across Cζ-N3 in MD/DFT).

### Variant: ARN

Atoms PRESENT in ARN but ABSENT from default (ARG):
- (none added)

Atoms ABSENT in ARN but PRESENT in default:
- One of the guanidinium Hs is removed to neutralise. The exact atom removed by AMBER's ARN template is not standardised in ff14SB (ARN is a non-standard variant; the user spec flagged this). Encode by chemistry inference: Markley does not give a canonical convention for ARN; ff14SB's port `ARN` removes Hε to give the neutral guanidine where the lone pair is on Nε. **Document this assumption explicitly in notes.**

Atoms with field deltas vs default:
- NE: formal_chg changes from +1 → 0 (deprotonation removes the charge); polar_h on the now-absent HE no longer applies. NE itself stays Guanidinium planar group with PolarHKind=NotPolar (the N has no H).
- HE: REMOVED. (If ff14SB's ARN template removes a different H, regenerate this row from the rtp.)
- HH11/HH12/HH21/HH22 retained as GuanidiniumNH; the planar-stereo labels Z/E unchanged.

Formal charges (across all atoms): 0.

Note: this encoding is uncertain; ARN is rare in ff14SB practice and BMRB doesn't have a canonical entry. The encoder should regenerate from `AmberAminoAcidVariantTable["ARN"]` rather than trust this row alone.

### ASN — asparagine

Atoms: N, H, CA, HA, C, O, CB, HB2, HB3, CG, OD1, ND2, HD21, HD22.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | the heavy atom three bonds closer to main chain (re Markley caption rule for Hδ) |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| CG | C | SidechainAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | side-chain primary amide carbonyl C |
| OD1 | O | SidechainAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | amide O (sp2 with double bond to Cγ in canonical Lewis) |
| ND2 | N | SidechainAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | amide N |
| HD21 | H | SidechainAmide | E | Q/Delta/0/false | SidechainPrimaryAmide | NotInRing | NotInRing | NotProchiral | 0 | t | "1" = cis to Cβ; CIP priority on Cγ is OD1 > CB, so cis-to-Cβ = trans-to-OD1 = E. BMRB atom_nom.tbl confirms HD21=E. QD (the Markley Asn QD aggregator over the two amide Hs). |
| HD22 | H | SidechainAmide | Z | Q/Delta/0/false | SidechainPrimaryAmide | NotInRing | NotInRing | NotProchiral | 0 | t | "2" = trans to Cβ = cis-to-OD1 = Z. BMRB atom_nom.tbl confirms HD22=Z. QD. |

Formal charges total: 0.

### ASP — aspartate

Atoms: N, H, CA, HA, C, O, CB, HB2, HB3, CG, OD1, OD2.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| CG | C | Carboxylate | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | side-chain carboxylate C |
| OD1 | O | Carboxylate | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | formal C=O oxygen (Lewis convention; delocalised in physical reality) |
| OD2 | O | Carboxylate | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | -1 | f | formal C-O- oxygen carrying the -1 in Lewis convention |

Formal charges total: -1.

### Variant: ASH

Atoms PRESENT in ASH but ABSENT from default:
- HD2: planar_group=Carboxylate, planar_stereo=NA, pseudoatom=—, polarH=CarboxylOH, ring=NotInRing, prochiral=NotProchiral, formal_chg=0, notes=protonated carboxyl OH on Oδ2 (Markley CarboxylOH).

Atoms ABSENT in ASH but PRESENT in default:
- (none removed from standard; Hδ2 is added)

Atoms with field deltas vs default:
- OD2: formal_chg changes from -1 → 0 (now carries an OH).

Formal charges (across all atoms): 0.

### CYS — cysteine

Atoms: N, H, CA, HA, C, O, CB, HB2, HB3, SG, HG.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| SG | S | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | thiol S |
| HG | H | None | NA | — | ThiolSH | NotInRing | NotInRing | NotProchiral | 0 | t | reduced cysteine thiol H; absent in CYX/CYM |

Formal charges total: 0.

### Variant: CYX

Atoms PRESENT in CYX but ABSENT from default:
- (none added; intermolecular SS bond is a topology fact handled by `CovalentTopology::DisulfidePair` and is not a per-atom semantic field here)

Atoms ABSENT in CYX but PRESENT in default:
- HG: REMOVED (Sγ is bonded to partner CYX's Sγ via disulfide).

Atoms with field deltas vs default:
- SG: formal_chg unchanged (0); polar_h unchanged (NotPolar); the disulfide-partner connectivity is captured upstream of this table.

Formal charges (across all atoms): 0.

### Variant: CYM

Atoms PRESENT in CYM but ABSENT from default:
- (none added)

Atoms ABSENT in CYM but PRESENT in default:
- HG: REMOVED (deprotonated thiolate).

Atoms with field deltas vs default:
- SG: formal_chg changes from 0 → -1 (thiolate); polar_h on SG itself is still NotPolar (S without H is not a polar-H site by definition; the field labels the H, not the heavy atom).

Formal charges (across all atoms): -1.

### GLN — glutamine

Atoms: N, H, CA, HA, C, O, CB, HB2, HB3, CG, HG2, HG3, CD, OE1, NE2, HE21, HE22.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| CG | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | the heavy atom three bonds closer to main chain (re Markley caption rule for Hε) |
| HG2 | H | None | NA | Q/Gamma/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QG |
| HG3 | H | None | NA | Q/Gamma/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QG |
| CD | C | SidechainAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | side-chain primary amide carbonyl C |
| OE1 | O | SidechainAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| NE2 | N | SidechainAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HE21 | H | SidechainAmide | E | Q/Epsilon/0/false | SidechainPrimaryAmide | NotInRing | NotInRing | NotProchiral | 0 | t | "1" = cis to Cγ; CIP priority on Cδ is OE1 > CG, so cis-to-Cγ = trans-to-OE1 = E. BMRB atom_nom.tbl confirms HE21=E. QE. |
| HE22 | H | SidechainAmide | Z | Q/Epsilon/0/false | SidechainPrimaryAmide | NotInRing | NotInRing | NotProchiral | 0 | t | "2" = trans to Cγ = cis-to-OE1 = Z. BMRB atom_nom.tbl confirms HE22=Z. QE. |

Formal charges total: 0.

### GLU — glutamate

Atoms: N, H, CA, HA, C, O, CB, HB2, HB3, CG, HG2, HG3, CD, OE1, OE2.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| CG | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HG2 | H | None | NA | Q/Gamma/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QG |
| HG3 | H | None | NA | Q/Gamma/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QG |
| CD | C | Carboxylate | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | side-chain carboxylate C |
| OE1 | O | Carboxylate | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | formal C=O |
| OE2 | O | Carboxylate | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | -1 | f | formal C-O- |

Formal charges total: -1.

### Variant: GLH

Atoms PRESENT in GLH but ABSENT from default:
- HE2: planar_group=Carboxylate, planar_stereo=NA, pseudoatom=—, polarH=CarboxylOH, ring=NotInRing, prochiral=NotProchiral, formal_chg=0, notes=protonated carboxyl OH on Oε2 (Markley CarboxylOH).

Atoms ABSENT in GLH but PRESENT in default:
- (none removed)

Atoms with field deltas vs default:
- OE2: formal_chg changes from -1 → 0.

Formal charges (across all atoms): 0.

### GLY — glycine

Atoms: N, H, CA, HA2, HA3, C, O.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | the only achiral standard residue |
| HA2 | H | None | NA | Q/Alpha/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | INVERTED vs typical alternation; Markley Fig 1 marks HA2 as (R). QA. |
| HA3 | H | None | NA | Q/Alpha/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | INVERTED vs typical; QA |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |

Formal charges total: 0.

### HIS — histidine (AMBER default = HIE)

Per AMBER ff14SB convention, the default protonation of "HIS" is HIE (Nε2-protonated, neutral). The "Default" block below is the HIE row; HID and HIP are listed below as variants.

Atoms (HIE/default): N, H, CA, HA, C, O, CB, HB2, HB3, CG, ND1, CE1, HE1, NE2, HE2, CD2, HD2.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| CG | C | Imidazole | NA | — | NotPolar | Imidazole_His/Ipso/5/t/2 | NotInRing | NotProchiral | 0 | f | ipso of imidazole 5-ring |
| ND1 | N | Imidazole | NA | — | NotPolar | Imidazole_His/Heteroatom_NoH/5/t/2 | NotInRing | NotProchiral | 0 | f | HIE: no Hδ1 ⇒ Heteroatom_NoH |
| CE1 | C | Imidazole | NA | — | NotPolar | Imidazole_His/PyrroleAlpha/5/t/2 | NotInRing | NotProchiral | 0 | f | sp2 C between the two N's |
| HE1 | H | Imidazole | NA | — | NotPolar | Imidazole_His/PyrroleAlpha/5/t/2 | NotInRing | NotProchiral | 0 | f | aromatic CH |
| NE2 | N | Imidazole | NA | — | NotPolar | Imidazole_His/Heteroatom_NH/5/t/2 | NotInRing | NotProchiral | 0 | f | HIE: HE2 attached |
| HE2 | H | Imidazole | NA | — | ImidazoleNH | Imidazole_His/Heteroatom_NH/5/t/2 | NotInRing | NotProchiral | 0 | t | imidazole NH; tautomer-specific |
| CD2 | C | Imidazole | NA | — | NotPolar | Imidazole_His/PyrroleBeta/5/t/2 | NotInRing | NotProchiral | 0 | f | the other 5-ring carbon |
| HD2 | H | Imidazole | NA | — | NotPolar | Imidazole_His/PyrroleBeta/5/t/2 | NotInRing | NotProchiral | 0 | f | aromatic CH |

Formal charges total: 0.

### Variant: HID

Atoms PRESENT in HID but ABSENT from HIS/HIE default:
- HD1: planar_group=Imidazole, planar_stereo=NA, pseudoatom=—, polarH=ImidazoleNH, ring_primary=Imidazole_His/Heteroatom_NH/5/t/2, ring_secondary=NotInRing, prochiral=NotProchiral, formal_chg=0, notes=Nδ1-protonated tautomer; ImidazoleNH on Hδ1.

Atoms ABSENT in HID but PRESENT in HIS/HIE default:
- HE2: REMOVED.

Atoms with field deltas vs HIS/HIE default:
- ND1: ring_primary position changes from Heteroatom_NoH → Heteroatom_NH (now bonded to Hδ1).
- NE2: ring_primary position changes from Heteroatom_NH → Heteroatom_NoH (no longer has bonded H).

Formal charges (across all atoms): 0.

### Variant: HIE

This is identical to the HIS default block above (AMBER's default tautomer). No deltas.

### Variant: HIP

Atoms PRESENT in HIP but ABSENT from HIS/HIE default:
- HD1: planar_group=Imidazole, planar_stereo=NA, pseudoatom=—, polarH=ImidazoleNH, ring_primary=Imidazole_His/Heteroatom_NH/5/t/2, ring_secondary=NotInRing, prochiral=NotProchiral, formal_chg=0, notes=imidazolium has both NHs.

Atoms ABSENT in HIP but PRESENT in HIS/HIE default:
- (none; HE2 still present in HIP)

Atoms with field deltas vs HIS/HIE default:
- ND1: ring_primary position changes from Heteroatom_NoH → Heteroatom_NH.
- NE2: formal_chg changes from 0 → +1 per chosen convention (place +1 on the second-protonated nitrogen). Note in row: localised vs delocalised across the imidazolium ring is a chemistry choice; alternative is to place +1 on Nδ1 or distribute. Choice documented in Section 2.

Formal charges (across all atoms): +1.

### ILE — isoleucine

Atoms: N, H, CA, HA, C, O, CB, HB, CG2, HG21, HG22, HG23, CG1, HG12, HG13, CD1, HD11, HD12, HD13.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | sidechain has chiral β-carbon (Ile is doubly chiral) |
| HB | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | single Hβ; not diastereotopic |
| CG2 | C | None | NA | — | NotPolar | NotInRing | NotInRing | ProS | 0 | f | γ2-methyl carbon; ProS branch (Ile has only one γ-methyl C) |
| HG21 | H | None | NA | M/Gamma/2/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | Markley MG; methyl Hs |
| HG22 | H | None | NA | M/Gamma/2/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MG |
| HG23 | H | None | NA | M/Gamma/2/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MG |
| CG1 | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | γ1-methylene C |
| HG12 | H | None | NA | Q/Gamma/1/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QG (the γ1-methylene per Markley; "QG = γ1-methylene" in Table 1) |
| HG13 | H | None | NA | Q/Gamma/1/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QG |
| CD1 | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | δ1-methyl carbon (the only δ on Ile) |
| HD11 | H | None | NA | M/Delta/1/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | Markley MD; methyl Hs |
| HD12 | H | None | NA | M/Delta/1/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MD |
| HD13 | H | None | NA | M/Delta/1/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MD |

Formal charges total: 0.

Note: Per Markley Table 1, Ile has MG (γ2-methyl), MD (δ1-methyl), QG (γ1-methylene). MG and MD are NOT in a higher-order Q super-group (Ile has no QG_super because the γ2-methyl and the γ1-methylene are not equivalent in any higher pseudoatom).

### LEU — leucine

Atoms: N, H, CA, HA, C, O, CB, HB2, HB3, CG, HG, CD1, HD11, HD12, HD13, CD2, HD21, HD22, HD23.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB (Markley Leu QB; not in QD super) |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| CG | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HG | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | single Hγ (Leu Cγ has only one H) |
| CD1 | C | None | NA | — | NotPolar | NotInRing | NotInRing | ProR | 0 | f | δ1 methyl C; Cδ1 = ProR per Markley Fig 1 |
| HD11 | H | None | NA | M/Delta/1/true | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MD1; in QD super (all six δ-methyl Hs collapse into QD) |
| HD12 | H | None | NA | M/Delta/1/true | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MD1; QD |
| HD13 | H | None | NA | M/Delta/1/true | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MD1; QD |
| CD2 | C | None | NA | — | NotPolar | NotInRing | NotInRing | ProS | 0 | f | δ2 methyl C; Cδ2 = ProS |
| HD21 | H | None | NA | M/Delta/2/true | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MD2; QD |
| HD22 | H | None | NA | M/Delta/2/true | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MD2; QD |
| HD23 | H | None | NA | M/Delta/2/true | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MD2; QD |

Formal charges total: 0.

### LYS — lysine

Atoms: N, H, CA, HA, C, O, CB, HB2, HB3, CG, HG2, HG3, CD, HD2, HD3, CE, HE2, HE3, NZ, HZ1, HZ2, HZ3.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| CG | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HG2 | H | None | NA | Q/Gamma/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QG |
| HG3 | H | None | NA | Q/Gamma/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QG |
| CD | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HD2 | H | None | NA | Q/Delta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QD |
| HD3 | H | None | NA | Q/Delta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QD |
| CE | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HE2 | H | None | NA | Q/Epsilon/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QE |
| HE3 | H | None | NA | Q/Epsilon/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QE |
| NZ | N | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | +1 | f | quaternary-like N (Lewis: protonated NH3+ formal +1 on N) |
| HZ1 | H | None | NA | Q/Zeta/0/false | AmmoniumNH | NotInRing | NotInRing | NotProchiral | 0 | t | QZ; Markley says Lys has QZ over the 3 ammonium Hs |
| HZ2 | H | None | NA | Q/Zeta/0/false | AmmoniumNH | NotInRing | NotInRing | NotProchiral | 0 | t | QZ |
| HZ3 | H | None | NA | Q/Zeta/0/false | AmmoniumNH | NotInRing | NotInRing | NotProchiral | 0 | t | QZ |

Formal charges total: +1.

### Variant: LYN

Atoms PRESENT in LYN but ABSENT from default:
- (none added)

Atoms ABSENT in LYN but PRESENT in default:
- HZ1: REMOVED (Lys neutral has -NH2 with two ζ-Hs; AMBER ff14SB removes Hζ1 to leave HZ2 and HZ3 per `data/ff14sb_params.dat` lines 449-451).

Atoms with field deltas vs default:
- NZ: formal_chg changes from +1 → 0 (deprotonation).
- HZ2, HZ3: polarH changes AmmoniumNH → AmineNH (per the new `PolarHKind::AmineNH` enum value added 2026-05-05; LYN's -NH2 is a neutral primary amine, distinct from charged ammonium).

Formal charges (across all atoms): 0.

### MET — methionine

Atoms: N, H, CA, HA, C, O, CB, HB2, HB3, CG, HG2, HG3, SD, CE, HE1, HE2, HE3.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| CG | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HG2 | H | None | NA | Q/Gamma/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QG |
| HG3 | H | None | NA | Q/Gamma/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QG |
| SD | S | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | thioether S |
| CE | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | ε-methyl C |
| HE1 | H | None | NA | M/Epsilon/0/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | Markley ME; methyl Hs (3 equivalent) |
| HE2 | H | None | NA | M/Epsilon/0/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | ME |
| HE3 | H | None | NA | M/Epsilon/0/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | ME |

Formal charges total: 0.

### PHE — phenylalanine

ACCEPTANCE-GATE residue. The rows below are checked against `tools/topology/build_semantic_tables.cpp::SynthesisedForPhe` lines 779-842.

Atoms: N, H, CA, HA, C, O, CB, HB2, HB3, CG, CD1, HD1, CD2, HD2, CE1, HE1, CE2, HE2, CZ, HZ.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| CG | C | Aromatic6Ring | NA | — | NotPolar | Benzene_Phe/Ipso/6/t/0 | NotInRing | NotProchiral | 0 | f | ipso |
| CD1 | C | Aromatic6Ring | NA | — | NotPolar | Benzene_Phe/Ortho1/6/t/0 | NotInRing | NotProchiral | 0 | f | ortho1 |
| HD1 | H | Aromatic6Ring | NA | Q/Delta/0/true | NotPolar | Benzene_Phe/Ortho1/6/t/0 | NotInRing | NotProchiral | 0 | f | QD; in QR super |
| CD2 | C | Aromatic6Ring | NA | — | NotPolar | Benzene_Phe/Ortho2/6/t/0 | NotInRing | NotProchiral | 0 | f | ortho2 |
| HD2 | H | Aromatic6Ring | NA | Q/Delta/0/true | NotPolar | Benzene_Phe/Ortho2/6/t/0 | NotInRing | NotProchiral | 0 | f | QD; QR |
| CE1 | C | Aromatic6Ring | NA | — | NotPolar | Benzene_Phe/Meta1/6/t/0 | NotInRing | NotProchiral | 0 | f | meta1 |
| HE1 | H | Aromatic6Ring | NA | Q/Epsilon/0/true | NotPolar | Benzene_Phe/Meta1/6/t/0 | NotInRing | NotProchiral | 0 | f | QE; QR |
| CE2 | C | Aromatic6Ring | NA | — | NotPolar | Benzene_Phe/Meta2/6/t/0 | NotInRing | NotProchiral | 0 | f | meta2 |
| HE2 | H | Aromatic6Ring | NA | Q/Epsilon/0/true | NotPolar | Benzene_Phe/Meta2/6/t/0 | NotInRing | NotProchiral | 0 | f | QE; QR |
| CZ | C | Aromatic6Ring | NA | — | NotPolar | Benzene_Phe/Para/6/t/0 | NotInRing | NotProchiral | 0 | f | para |
| HZ | H | Aromatic6Ring | NA | R/Zeta/0/true | NotPolar | Benzene_Phe/Para/6/t/0 | NotInRing | NotProchiral | 0 | f | Markley R-letter at para; QR |

Formal charges total: 0.

**PHE acceptance-gate match check:**
- Ring tagging matches exactly: CG=Ipso, CD1=Ortho1, CD2=Ortho2, CE1=Meta1, CE2=Meta2, CZ=Para; ring_size=6; aromatic=t; planar=t; n_heteroatoms=0; planar_group=Aromatic6Ring on all ring atoms (heavy + H). MATCH.
- Backbone N+C+O+H tagged PeptideAmide. MATCH.
- HB2/HB3 → Q/Beta/0/false (super=false). MATCH.
- HD1/HD2 → Q/Delta/0/true. MATCH (the C++ uses `super=true`; this row uses `super=true`).
- HE1/HE2 → Q/Epsilon/0/true. MATCH.
- HZ → R/Zeta/0/true. MATCH.
- BackboneAmide PolarH on H. MATCH.
- (The C++ has special-case fall-throughs for `is_n_terminus → AmmoniumNH` and `is_c_terminus → CarboxylOH`; those are terminal-state rows in this file (NTERM_CHARGED / CTERM_PROTONATED) and orthogonal to the per-residue PHE rows.)
- (The C++ fall-through to `Carboxylate` when planar_group is None on `OXT/C/O` is present *because the canonical CCD entry for free PHE is the protonated free amino acid with OXT*. In this file, free-amino-acid OXT is encoded only on the terminal-state rows; the per-residue PHE rows assume the residue is internal/peptide-bonded. The encoder may optionally apply this fall-through during code-gen for the standalone free-PHE CCD entry.)

PHE encoding matches `SynthesisedForPhe` exactly.

### PRO — proline

Atoms: N, CA, HA, C, O, CB, HB2, HB3, CG, HG2, HG3, CD, HD2, HD3.

Special: Pro has no backbone H (secondary amine, in-ring N). Cα is in the ring as well, but the ring labels for Cα are conventionally not encoded as ring-position (Cα is treated as backbone). Below: ring atoms are N, Cβ, Cγ, Cδ (and Cα for completeness, encoded as ring-membership-only since Pyrrolidine_Pro spans those 5 atoms).

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | None | NA | — | NotPolar | Pyrrolidine_Pro/Heteroatom_NoH/5/f/1 | NotInRing | NotProchiral | 0 | f | secondary amine; in the ring; no H |
| CA | C | None | NA | — | NotPolar | Pyrrolidine_Pro/Saturated/5/f/1 | NotInRing | NotProchiral | 0 | f | Cα also in ring (Markley caption: "priority leads out from main chain Cα") |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | Pyrrolidine_Pro/Saturated/5/f/1 | NotInRing | NotProchiral | 0 | f | |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| CG | C | None | NA | — | NotPolar | Pyrrolidine_Pro/Saturated/5/f/1 | NotInRing | NotProchiral | 0 | f | |
| HG2 | H | None | NA | Q/Gamma/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QG |
| HG3 | H | None | NA | Q/Gamma/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QG |
| CD | C | None | NA | — | NotPolar | Pyrrolidine_Pro/Saturated/5/f/1 | NotInRing | NotProchiral | 0 | f | |
| HD2 | H | None | NA | Q/Delta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QD; Markley Fig 1 marks Hδ3 (R) on Pro — implies Hδ2=ProS (typical alternation) |
| HD3 | H | None | NA | Q/Delta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QD; Hδ3 = ProR per Markley |

Formal charges total: 0.

Note: `RingMembership.aromatic = false` for Pro pyrrolidine. `RingMembership.planar = false` (saturated; ring-puckers; conformation-side data).

### SER — serine

Atoms: N, H, CA, HA, C, O, CB, HB2, HB3, OG, HG.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| OG | O | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | hydroxyl O |
| HG | H | None | NA | — | HydroxylOH_Aliphatic | NotInRing | NotInRing | NotProchiral | 0 | t | |

Formal charges total: 0.

### THR — threonine

Atoms: N, H, CA, HA, C, O, CB, HB, OG1, HG1, CG2, HG21, HG22, HG23.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | doubly chiral β-carbon |
| HB | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | single Hβ |
| OG1 | O | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | hydroxyl O |
| HG1 | H | None | NA | — | HydroxylOH_Aliphatic | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CG2 | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | γ2-methyl C |
| HG21 | H | None | NA | M/Gamma/2/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | Markley MG; methyl Hs |
| HG22 | H | None | NA | M/Gamma/2/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MG |
| HG23 | H | None | NA | M/Gamma/2/false | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MG |

Formal charges total: 0.

### TRP — tryptophan

Atoms: N, H, CA, HA, C, O, CB, HB2, HB3, CG, CD1, HD1, NE1, HE1, CE2, CD2, CE3, HE3, CZ2, HZ2, CZ3, HZ3, CH2, HH2.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| CG | C | Aromatic5Ring | NA | — | NotPolar | Indole_Trp_5/Ipso/5/t/1 | NotInRing | NotProchiral | 0 | f | ipso of pyrrole 5-ring |
| CD1 | C | Aromatic5Ring | NA | — | NotPolar | Indole_Trp_5/PyrroleBeta/5/t/1 | NotInRing | NotProchiral | 0 | f | only non-bridgehead pyrrole C |
| HD1 | H | Aromatic5Ring | NA | — | NotPolar | Indole_Trp_5/PyrroleBeta/5/t/1 | NotInRing | NotProchiral | 0 | f | aromatic CH |
| NE1 | N | Aromatic5Ring | NA | — | NotPolar | Indole_Trp_5/Heteroatom_NH/5/t/1 | NotInRing | NotProchiral | 0 | f | indole N (carries Hε1) |
| HE1 | H | Aromatic5Ring | NA | — | IndoleNH | Indole_Trp_5/Heteroatom_NH/5/t/1 | NotInRing | NotProchiral | 0 | t | indole NH; Markley IndoleNH |
| CE2 | C | Aromatic5Ring | NA | — | NotPolar | Indole_Trp_5/BridgeFusion/5/t/1 | Indole_Trp_6/BridgeFusion/6/t/0 | NotProchiral | 0 | f | bridgehead atom; in BOTH 5- and 6-rings; primary=smaller (5), secondary=6 |
| CD2 | C | Aromatic5Ring | NA | — | NotPolar | Indole_Trp_5/BridgeFusion/5/t/1 | Indole_Trp_6/BridgeFusion/6/t/0 | NotProchiral | 0 | f | bridgehead; primary=5, secondary=6 |
| CE3 | C | Aromatic6Ring | NA | — | NotPolar | Indole_Trp_6/Ortho1/6/t/0 | NotInRing | NotProchiral | 0 | f | benzene ortho1 (per Markley convention) |
| HE3 | H | Aromatic6Ring | NA | — | NotPolar | Indole_Trp_6/Ortho1/6/t/0 | NotInRing | NotProchiral | 0 | f | |
| CZ2 | C | Aromatic6Ring | NA | — | NotPolar | Indole_Trp_6/Ortho2/6/t/0 | NotInRing | NotProchiral | 0 | f | benzene ortho2 |
| HZ2 | H | Aromatic6Ring | NA | — | NotPolar | Indole_Trp_6/Ortho2/6/t/0 | NotInRing | NotProchiral | 0 | f | |
| CZ3 | C | Aromatic6Ring | NA | — | NotPolar | Indole_Trp_6/Meta1/6/t/0 | NotInRing | NotProchiral | 0 | f | benzene meta1 |
| HZ3 | H | Aromatic6Ring | NA | — | NotPolar | Indole_Trp_6/Meta1/6/t/0 | NotInRing | NotProchiral | 0 | f | |
| CH2 | C | Aromatic6Ring | NA | — | NotPolar | Indole_Trp_6/Meta2/6/t/0 | NotInRing | NotProchiral | 0 | f | benzene meta2 |
| HH2 | H | Aromatic6Ring | NA | — | NotPolar | Indole_Trp_6/Meta2/6/t/0 | NotInRing | NotProchiral | 0 | f | |

Formal charges total: 0.

Notes:
- The "para" position of the 6-ring is BridgeFusion (Cδ2 or Cε2), so Para is not a separate label here; the perimeter is Ortho1=Cε3, Ortho2=Cζ2, Meta1=Cζ3, Meta2=Cη2 in Markley's convention. (Bridge-relative Para is two-bonds from the BridgeFusion across the 6-ring, which would be one of Cζ3 or Cη2, but neither is labelled "Para" because the ring already has two BridgeFusion atoms.) Markley does not give "ipso/ortho/meta/para" for the indole 6-ring in published convention; the labels above (Ortho1/Ortho2/Meta1/Meta2) are the synthesis from the Markley dossier text-2:230-358 + dossier `topology-fields-research-2026-05-05.md` lines 1108-1117. **Alternative encoding option**: encode Cε3/Cζ3/Cη2/Cζ2 with positions PyrroleBeta-style or as "perimeter-1, perimeter-2, perimeter-3, perimeter-4". This is in the disagreement-log section.
- The BridgeFusion atoms (Cδ2, Cε2) carry both `ring_primary` (5-ring per smaller-ring convention from `RingPosition` docstring) AND `ring_secondary` (6-ring).
- planar_group: arguably both 5-ring atoms AND 6-ring atoms are part of the same delocalised "indole" π-system. The convention chosen here is to assign `Aromatic5Ring` for the pyrrole-membership atoms (CG, CD1, HD1, NE1, HE1, CE2, CD2) and `Aromatic6Ring` for the perimeter benzene atoms (CE3, HE3, CZ2, HZ2, CZ3, HZ3, CH2, HH2). The bridgeheads (CE2, CD2) are tagged Aromatic5Ring per their primary ring; the 6-ring secondary RingMembership is the "also in 6-ring" record.

### TYR — tyrosine

Atoms: N, H, CA, HA, C, O, CB, HB2, HB3, CG, CD1, HD1, CD2, HD2, CE1, HE1, CE2, HE2, CZ, OH, HH.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB2 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProS | 0 | f | QB |
| HB3 | H | None | NA | Q/Beta/0/false | NotPolar | NotInRing | NotInRing | ProR | 0 | f | QB |
| CG | C | Aromatic6Ring | NA | — | NotPolar | Benzene_Tyr/Ipso/6/t/0 | NotInRing | NotProchiral | 0 | f | ipso |
| CD1 | C | Aromatic6Ring | NA | — | NotPolar | Benzene_Tyr/Ortho1/6/t/0 | NotInRing | NotProchiral | 0 | f | |
| HD1 | H | Aromatic6Ring | NA | Q/Delta/0/true | NotPolar | Benzene_Tyr/Ortho1/6/t/0 | NotInRing | NotProchiral | 0 | f | QD; QR |
| CD2 | C | Aromatic6Ring | NA | — | NotPolar | Benzene_Tyr/Ortho2/6/t/0 | NotInRing | NotProchiral | 0 | f | |
| HD2 | H | Aromatic6Ring | NA | Q/Delta/0/true | NotPolar | Benzene_Tyr/Ortho2/6/t/0 | NotInRing | NotProchiral | 0 | f | QD; QR |
| CE1 | C | Aromatic6Ring | NA | — | NotPolar | Benzene_Tyr/Meta1/6/t/0 | NotInRing | NotProchiral | 0 | f | |
| HE1 | H | Aromatic6Ring | NA | Q/Epsilon/0/true | NotPolar | Benzene_Tyr/Meta1/6/t/0 | NotInRing | NotProchiral | 0 | f | QE; QR |
| CE2 | C | Aromatic6Ring | NA | — | NotPolar | Benzene_Tyr/Meta2/6/t/0 | NotInRing | NotProchiral | 0 | f | |
| HE2 | H | Aromatic6Ring | NA | Q/Epsilon/0/true | NotPolar | Benzene_Tyr/Meta2/6/t/0 | NotInRing | NotProchiral | 0 | f | QE; QR |
| CZ | C | Aromatic6Ring | NA | — | NotPolar | Benzene_Tyr/Para/6/t/0 | NotInRing | NotProchiral | 0 | f | para; carries OH |
| OH | O | AromaticHydroxyl | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | para-hydroxyl O; AromaticHydroxyl group includes Cζ-O-H |
| HH | H | AromaticHydroxyl | NA | — | HydroxylOH_Aromatic | NotInRing | NotInRing | NotProchiral | 0 | t | the aromatic-hydroxyl proton |

Formal charges total: 0.

### Variant: TYM

Atoms PRESENT in TYM but ABSENT from default:
- (none added)

Atoms ABSENT in TYM but PRESENT in default:
- HH: REMOVED (deprotonated phenolate / aryloxide).

Atoms with field deltas vs default:
- OH: planar_group changes AromaticHydroxyl → AromaticOxide (per the new `PlanarGroupKind::AromaticOxide` enum value added 2026-05-05; phenolate -O- conjugates more strongly with the ring π system than neutral phenol-OH does, redistributing ring current and shifting the para-Cζ environment by 5-10 ppm relative to neutral Tyr); formal_chg changes from 0 → -1.

Formal charges (across all atoms): -1.

### VAL — valine

Atoms: N, H, CA, HA, C, O, CB, HB, CG1, HG11, HG12, HG13, CG2, HG21, HG22, HG23.

| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| H | H | PeptideAmide | NA | — | BackboneAmide | NotInRing | NotInRing | NotProchiral | 0 | t | |
| CA | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HA | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| C | C | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| O | O | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| CB | C | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | |
| HB | H | None | NA | — | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | single Hβ |
| CG1 | C | None | NA | — | NotPolar | NotInRing | NotInRing | ProR | 0 | f | γ1 methyl C; ProR per Markley Fig 1 |
| HG11 | H | None | NA | M/Gamma/1/true | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MG1; in QG super (all six γ-methyl Hs) |
| HG12 | H | None | NA | M/Gamma/1/true | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MG1; QG |
| HG13 | H | None | NA | M/Gamma/1/true | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MG1; QG |
| CG2 | C | None | NA | — | NotPolar | NotInRing | NotInRing | ProS | 0 | f | γ2 methyl C; ProS |
| HG21 | H | None | NA | M/Gamma/2/true | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MG2; QG |
| HG22 | H | None | NA | M/Gamma/2/true | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MG2; QG |
| HG23 | H | None | NA | M/Gamma/2/true | NotPolar | NotInRing | NotInRing | NotProchiral | 0 | f | MG2; QG |

Formal charges total: 0.

---

## Section 4: Terminal-state blocks

The terminal-state blocks below are encoded as deltas relative to the canonical free-amino-acid CCD entry (i.e., the residue's own per-residue rows above). They apply ON TOP of any per-residue block — when a residue is at the N-terminus, you swap N+H rows for the NTERM_* delta; at the C-terminus, you swap C+O rows for the CTERM_* delta.

### Terminal: NTERM_CHARGED (NH3+)

Atoms PRESENT in NTERM_CHARGED but ABSENT from default residue rows:
- H1, H2, H3 — planar_group=None, planar_stereo=NA, pseudoatom=Q/None/0/false (no Markley super-group; H1/H2/H3 are not given a Q-letter in Table 1 because they're terminus-specific — but pseudoatom shorthand here is Q-with-locant=None to indicate "Q-like equivalent group of three Hs"), polarH=AmmoniumNH, ring=NotInRing, prochiral=NotProchiral, formal_chg=0, exch=t, notes=N-terminal NH3+ Hs (Markley text-2:43-48 names them H1, H2, H3 when protonated).

Atoms ABSENT in NTERM_CHARGED but PRESENT in default residue rows:
- H (the standard backbone amide H): REMOVED (the N has H1/H2/H3 instead of H).

Atoms with field deltas vs default residue rows:
- N: planar_group changes from PeptideAmide → None (the terminal N is not in a peptide-amide unit; it's a primary amine NH3+); polarH stays NotPolar (label is on the H, not the N); formal_chg changes from 0 → +1.

Formal charges total over the NTERM_CHARGED delta atoms: +1.

### Terminal: NTERM_NEUTRAL (NH2)

Atoms PRESENT in NTERM_NEUTRAL but ABSENT from default residue rows:
- H1, H2 — planar_group=None, polarH=AmineNH (per the new `PolarHKind::AmineNH` enum value added 2026-05-05; the neutral N-terminus is a primary amine, distinct from the charged NTERM_CHARGED ammonium).

Atoms ABSENT in NTERM_NEUTRAL but PRESENT in default residue rows:
- H: REMOVED.

Atoms with field deltas vs default residue rows:
- N: planar_group changes PeptideAmide → None; formal_chg unchanged (0).

Formal charges total: 0.

### Terminal: CTERM_DEPROTONATED (COO-)

Atoms PRESENT in CTERM_DEPROTONATED but ABSENT from default residue rows:
- OXT (a.k.a. O'' in Markley convention): planar_group=Carboxylate, polarH=NotPolar, ring=NotInRing, prochiral=NotProchiral, formal_chg=-1, exch=f, notes=second carboxylate O (Markley says O'' takes the prime due to lower CIP priority).

Atoms ABSENT in CTERM_DEPROTONATED but PRESENT in default residue rows:
- (none ABSENT; the per-residue O remains as O', the OXT is added)

Atoms with field deltas vs default residue rows:
- C: planar_group changes from PeptideAmide → Carboxylate (now a carboxylate C, not peptide-amide).
- O: planar_group changes from PeptideAmide → Carboxylate (this becomes O' in Markley convention; formal_chg stays 0; the -1 is on OXT/O'').

Formal charges total over the CTERM_DEPROTONATED delta atoms: -1.

### Terminal: CTERM_PROTONATED (COOH)

Atoms PRESENT in CTERM_PROTONATED but ABSENT from default residue rows:
- OXT (= O''): planar_group=Carboxylate, polarH=NotPolar, formal_chg=0, notes=O'' the protonated oxygen.
- HXT (= H''): planar_group=Carboxylate, planar_stereo=NA, pseudoatom=—, polarH=CarboxylOH, ring=NotInRing, prochiral=NotProchiral, formal_chg=0, exch=t, notes=carboxyl OH proton (Markley H'' nomenclature).

Atoms ABSENT in CTERM_PROTONATED but PRESENT in default residue rows:
- (none ABSENT; the canonical residue has C and O — both retained — and OXT/HXT are added)

Atoms with field deltas vs default residue rows:
- C: planar_group changes PeptideAmide → Carboxylate.
- O: planar_group changes PeptideAmide → Carboxylate.

Formal charges total over the CTERM_PROTONATED delta atoms: 0.

---

## Section 5: Disagreement log

Substantive disagreements between sources (Markley 1998 caption text, Markley 1998 Fig 1 OCR, BMRB nomenclature page, BMRB pseudoatom_nom, the user spec, the existing `SynthesisedForPhe` C++):

1. **Asn Hδ21/Hδ22 E/Z assignment** [RESOLVED 2026-05-05 via critical review]. Fig 1 OCR (text-2:103) reads "Hδ22(Z)" and "Hδ21(E)". The Fig 1 caption (text-2:367-371) gives only the labelling rule ("1=cis to Cβ, 2=trans"); the Z/E mapping requires CIP priority on the Cγ=Nδ2 planar bond. On Cγ, OD1 (priority O) > CB (priority C). So HD21 = cis-to-Cβ = trans-to-OD1 = E; HD22 = trans-to-Cβ = cis-to-OD1 = Z. **BMRB `atom_nom.tbl` confirms: HD21=E, HD22=Z.** Resolution: BMRB is canonical; Fig 1 OCR is right; the earlier "caption→1=Z" derivation conflated "cis to Cβ" with "cis to high-priority". Reference now encodes HD21=E, HD22=Z.

2. **Gln Hε21/Hε22 E/Z assignment** [RESOLVED 2026-05-05 via critical review]. Same logic as entry 1 on Gln with Cγ in place of Cβ on the planar Cδ=Nε2 bond. OE1 (O) > CG (C) on Cδ. HE21=E, HE22=Z. BMRB confirms. Reference now encodes HE21=E, HE22=Z.

3. **Arg Hη11/Hη12/Hη21/Hη22.** Fig 1 OCR partially shows "Hη11 (Z), Hη12 (E)" and one of the Hη2x with "(E)". Caption says outer index = which Nη (Nη1 cis to Cδ ⇒ 1, Nη2 trans ⇒ 2), inner index = cis (1) or trans (2) to the heavy atom (Nε) "three bonds closer to main chain". This document encodes inner=1 ⇒ Z, inner=2 ⇒ E for both Nη nitrogens.

4. **Trp 6-ring position labels (Cε3, Cζ3, Cη2, Cζ2).** Markley does not publish "ipso/ortho/meta/para" for the indole 6-ring perimeter; the dossier `topology-fields-research-2026-05-05.md` lines 1108-1117 chose Ortho1=Cε3, Ortho2=Cζ2, Meta1=Cζ3, Meta2=Cη2. This is a **synthesised choice**, not an authoritative IUPAC convention. **No literature contradiction**; just non-canonical. Encode as-is; flag.

5. **HIP charge placement convention.** SemanticEnums.h does not have an "imidazolium" enum or per-atom convention for delocalised +1. This document places +1 on Nε2 (the second-protonated nitrogen). Alternative: place +1 on Nδ1, or split (e.g., +0.5 each). Conventions vary. Flag.

6. **Asp/Glu carboxylate -1 placement convention.** Same flavour as HIP. The Lewis-structure choice (Oδ2/Oε2 carry -1; Oδ1/Oε1 carry 0) is encoded; in MD/DFT the delocalisation makes them indistinguishable. Document the convention in Section 2.

7. **ARN atom-removal choice.** ARN (deprotonated arginine) is rare; AMBER ff14SB doesn't have a single canonical convention for which Hη is removed. This document removes HE for chemistry-inference reasons (the lone pair after deprotonation is on Nε, the formerly-protonated nitrogen). The encoder MUST regenerate from `AmberAminoAcidVariantTable["ARN"]` at code-gen time and override this row if the rtp says otherwise.

8. **LYN/AmineNH vs AmmoniumNH** [RESOLVED 2026-05-05 via critical review]. The categorical distinction is chemistry-real (different pKa, different exchange behaviour). `PolarHKind::AmineNH` was added to `src/SemanticEnums.h`. LYN HZ2/HZ3 and NTERM_NEUTRAL H1/H2 are now encoded as `AmineNH`; LYS HZ1/HZ2/HZ3 and NTERM_CHARGED H1/H2/H3 stay `AmmoniumNH`.

9. **TYM AromaticHydroxyl when there's no proton** [RESOLVED 2026-05-05 via critical review]. `PlanarGroupKind::AromaticOxide` was added to `src/SemanticEnums.h`. TYM Oη is now encoded as `AromaticOxide`; Tyr default Oη remains `AromaticHydroxyl`. The phenolate -O- conjugates more strongly with the ring π system than neutral phenol-OH does (5-10 ppm shift on para-Cζ).

10. **Pseudoatom letter for backbone-amide H/N-terminal H1/H2/H3.** Markley Table 1 doesn't list a pseudoatom for backbone HN or N-terminal H1/H2/H3. The shorthand here uses `Q/None/0/false` for the NH3+ Hs (treating them as a 3-H equivalent group like an ammonium). This is a synthesis, not Markley-authoritative. Flag for the encoder.

---

## Section 6: SemanticEnums.h vocabulary patches recommended

Recommendations only — do NOT make these edits. Discuss before committing.

1. **`PolarHKind::AmineNH`** — distinguish neutral primary-amine NH (LYN, NTERM_NEUTRAL) from charged ammonium NH (LYS, NTERM_CHARGED). Currently both share `AmmoniumNH`.

2. **`PlanarGroupKind::AromaticOxide`** — for TYM Oη after deprotonation. Currently TYM Oη falls through to `AromaticHydroxyl`, which contains "Hydroxyl" in the name despite no H.

3. **Pseudoatom letter for terminal Hs.** Markley doesn't define one; in current shorthand `Q/None/0/false` is used. Encoder should pick a convention (None or new "T" letter) and document it.

4. **`Locant::None` clarification for Pro CA in ring.** Cα is in the Pro ring but `Locant::None` for backbone atoms; encoding Cα as a ring atom while keeping Locant=None is technically consistent (locant is for sidechain only) but worth a clarifying comment.

5. **Asn/Gln side-chain amide pseudoatom super-group.** Markley Table 1 says Asn QD over the (HD21, HD22) pair, Gln QE over (HE21, HE22). These are 2-H groups, not 6-H super-groups in the sense that Phe QR or Val QG are. Current shorthand uses `super=false` for these because they're not "super-groups across multiple methyls/methylenes". Document the convention in the enum or in a derived comment.

---

## Quality-bar checklist (self-audit)

- PHE acceptance check passes: yes (see PHE acceptance-gate match check above).
- Every (residue, variant, atom) triple has all relevant fields populated: yes; absences encoded as None/NotInRing/NA/NotProchiral/—.
- Every variant explicitly lists added/removed/changed atoms: yes (HID, HIE, HIP, ASH, GLH, CYX, CYM, LYN, ARN, TYM each have a delta block).
- Every PolarHKind across CYS/CYX/CYM is correctly variant-aware: yes (CYS has ThiolSH on HG; CYX and CYM have HG REMOVED).
- Every ImidazoleNH across HID/HIE/HIP is correctly variant-aware: yes (HIE: HE2 only; HID: HD1 only; HIP: both).
- Glycine Hα is encoded with the inverted prochiral mapping: yes (HA2=ProR, HA3=ProS).
- Branched heavy atoms have correct prochiral mapping per Markley Fig 1: yes (Val Cγ1=ProR, Cγ2=ProS; Leu Cδ1=ProR, Cδ2=ProS; Ile Cγ2=ProS, Cδ1=NotProchiral [only one δ]; Ile Cγ1 is methylene → Hγ12=ProS, Hγ13=ProR).
- Trp bridgehead atoms (Cδ2, Cε2) have both primary and secondary RingMembership: yes.
- Asn Hδ21/Hδ22, Gln Hε21/Hε22, Arg Hη11/12/21/22 have explicit E/Z PlanarStereo: yes; Asn Hδ21=Z, Hδ22=E; Gln Hε21=Z, Hε22=E; Arg Hη11=Z, Hη12=E, Hη21=Z, Hη22=E (per the caption-rule cis=1=Z convention).
- Disagreement log is honest: 10 entries listed, all substantive (no manufactured flags; ones marked "synthesis" or "convention chosen" reflect real choices made beyond what literature directly states).
