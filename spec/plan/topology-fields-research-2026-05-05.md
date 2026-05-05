# Topology fields — research synthesis (2026-05-05)

Research dossier for the eight pending atom-level semantic fields on
`LegacyAmberTopology`. Built from the local NMR-shielding reference corpus
(text under `references-text/`, summaries under `references-meta/`) plus
external authoritative sources (CCD / mmCIF dictionary, BMRB nomenclature
tables, RDKit, Open Babel). Citations of local material use absolute
file:line where the line carries the claim. External sources are linked.

This is a research synthesis, not a design spec. Recommendations are
framed as enums and shapes; the user makes the integration call.

The eight fields:
1. **Locant** — Greek-letter position label
2. **BranchIndex** — disambiguator when atoms share a locant
3. **DiastereotopicIndex** — HB2 vs HB3 for prochiral methylene
4. **ProchiralStereo** — pro-R / pro-S per Markley 1998
5. **PlanarStereo** — sp2 / planar group classification + the conformation question
6. **PseudoatomClass** — Markley pseudoatom taxonomy
7. **PolarHKind** — functional-group classification of polar hydrogens
8. **RingPosition** — ring-atom positional information

---

## Field 1: Locant (Greek-letter position label)

### What it represents

The Greek-letter position label that places an atom along the side chain,
counting outward from Cα. It is a **graph-distance role on the sidechain
backbone**, not a free index. Settled by IUPAC-IUB tentative rules of 1969
(referenced in Markley 1998 §2.1.1 as "ref. 20"), and still the thing the
1998 IUPAC-IUBMB-IUPAB recommendations for protein NMR adopt verbatim.

### Markley / IUPAC convention — full enum

From `/shared/2026Thesis/nmr-shielding/references-text/markley-1998-iupac-nmr-nomenclature-recommendations-text-2.txt:25-56`:

> "These rules designate all atoms by Greek letters (or Roman counterparts)
> and employ the main-chain precedence rule for numbering prochiral sites.
> [...] In the peptide backbone, the nitrogen is denoted by N, its attached
> hydrogen by H, and the carbonyl carbon and oxygen by C and O,
> respectively." (text-2:32-37)
>
> "In situations where Greek letters are unavailable [...] they may be
> replaced by upper case Roman letters (α = A, β = B, γ = G, δ = D,
> ε = E, ζ = Z, η = H)." (text-2:50-56)

The full enum, in canonical order:

| Greek | Roman | What it labels |
|-------|-------|---------------|
| (none) | N | Backbone amide nitrogen |
| (none) | H or HN | Backbone amide hydrogen — Markley recommends `HN` over `H` (text-2:33-37) |
| (none) | Cα, Hα, C, O | Backbone Cα + its hydrogen, carbonyl C, carbonyl O |
| α | A | Cα; in glycine and as ring-system reference |
| β | B | Cβ — first sidechain heavy atom |
| γ | G | Cγ |
| δ | D | Cδ |
| ε | E | Cε |
| ζ | Z | Cζ |
| η | H | Cη (Arg, Tyr-OH proton has special treatment) |

Backbone atoms have **no locant** — they are addressed by name (N, H, CA,
HA, C, O), not Greek letter. Only sidechain atoms carry locants.

### Special considerations baked into the convention

- **Proline**: `markley...text-2:24-28` — "For cyclic amino acids such as
  proline, the priority leads out from the main chain Cα atom rather than
  from the N atom (i.e., the priority is given by Cα > Cβ > Cγ > ... > N)."
  So Pro's Cδ connects to N but its locant is δ, not anything implying
  proximity to N. The cycle does not break the locant scheme; it just
  forces the walk to terminate at N.
- **N-terminal amine**: `text-2:43-48` — H1, H2, H3 (protonated) or H1, H2
  (unprotonated). These are NOT pseudo-locants — they are positional
  numbers without Greek letters. They label atoms on the backbone N, not
  on Cα.
- **C-terminal carboxyl**: `text-2:48-50` — O' and O'' for the two oxygens,
  H'' for the protonated form (because lower CIP priority gets the prime —
  see BMRB nomenclature page).

### What shift predictors actually use

- **SHIFTX2 (Han 2011)**: predicts 40 atom types (`han-2011...text-3:165-168`),
  named per IUPAC. Side-chain table in Table 2 (`han-2011...text-4:1-39`)
  uses CD, CD1, CD2, CE, CE1, CE2, CG, CG1, CG2, CZ, HB, HB2, HB3, HD1,
  HD2, HD3, HE, HE1, HE2, HE3, HG, HG1, HG12, HG13, HG2, HG3, HZ. These
  are exactly Greek-letter locants (B/G/D/E/Z) with branch indices.
- **ProCS15 (Larsen 2015)**: backbone + Cβ only, but parameterizes
  side-chain effects via χ1...χ4 dihedrals which are themselves locant-
  based (Cα-Cβ-Cγ-Cδ etc.). Larsen text-3:64-68 mentions per-side-chain
  scans on χ1 through χ4.
- **CH3Shift (Sahakyan 2011)**: predicts methyl groups labeled by
  Cβ, Cγ1, Cγ2, Cδ1 etc. (see `sahakyan-2011...text-2:621-672` Table 1).
  The carbon-side label IS the locant.

### Most-granular reasonable enum / shape

```
enum class Locant : uint8_t {
    Alpha   = 0,   // Cα, Hα — only on Gly (αHs) and as ring reference
    Beta    = 1,
    Gamma   = 2,
    Delta   = 3,
    Epsilon = 4,
    Zeta    = 5,
    Eta     = 6,
    None    = 255  // backbone atoms, terminal hydrogens, caps
};
```

The Greek vs Roman issue is purely presentation. Internal storage is
the integer; rendering picks Greek (UTF-8) or Roman (ASCII) per output.

### Open questions / known controversies

None substantive. The 1969 rules are stable and adopted by every
reference predictor in the corpus. The only friction is **typographic**:
PDB stores ASCII (CG1, CG2), BMRB and journal text use Greek (Cγ1, Cγ2),
some old code uses 'B' as both backbone (legacy Wüthrich) and β (Roman
counterpart). Internal representation should be locant-as-int and never
attempt to round-trip through strings.

---

## Field 2: BranchIndex (disambiguator when atoms share a locant)

### What it represents

When two or more atoms share a Greek-letter locant, BranchIndex
distinguishes them. Three distinct chemistry contexts use this:

1. **Branched aliphatic** — Val Cγ1/Cγ2, Leu Cδ1/Cδ2, Ile Cγ1/Cγ2 + Cδ1
   (Ile has only one δ, but Cγ1 vs Cγ2 split is real). These are
   diastereotopic at the **carbon** level, distinct heavy atoms.
2. **Aromatic ring positions** — Phe/Tyr Cδ1/Cδ2 and Cε1/Cε2; His Nδ1
   vs Cδ2 / Cε1 vs Nε2; Trp Cδ1 vs Cδ2, Cε2 vs Cε3, Cζ2 vs Cζ3.
3. **Diastereotopic methylene hydrogens** — HB2/HB3, HG2/HG3 etc.
   This case is what `DiastereotopicIndex` (Field 3) is for, NOT
   `BranchIndex`. Separation matters because the rules for assigning
   `1` vs `2` differ between the heavy-atom and prochiral-H contexts.

### Markley / IUPAC convention

The numbering rule for branched heavy atoms is the **CIP-priority
clockwise rule**, given in the Figure 1 caption
(`markley-1998...text-2:357-372`):

> "The Cα or the substituent closer to Cα (in the order Cα, Cβ, Cγ, ...)
> takes precedence over atoms in branches in defining stereochemical
> relationships. For example, if tetrahedral carbon C has four substituents
> X, Y, Z, and Z' (with priority X > Y > Z = Z'; i.e., Z and Z' are
> diastereotopic substituents designated provisionally as unprimed and
> primed), their numbering is derived as follows: if one sights down the
> X—C axis (with the X atom toward the viewer), the equivalent atoms,
> Z and Z', are designated Z2 and Z3, such that Y, Z2, and Z3 follow a
> clockwise orientation."

So the rule is: **sight from the highest-priority substituent toward the
chiral center; number the remaining diastereotopic atoms 2 and 3 such
that Y, Z2, Z3 are clockwise**. This determines BranchIndex per IUPAC.

For aromatic rings (`text-2:368-371`):

> "Numbering of Phe and Tyr rings gives higher priority to the atom with
> the smaller absolute value of the χ2 torsion angle (ref. 20). For
> example, the ring carbons of Phe and Tyr lying in the plane with the
> smaller χ2 torsion angle are designated as Cδ1 and Cε1."

For arginine sidechain NH2 nitrogens (`text-2:362-368`):

> "The side-chain-NH2 nitrogens of Arg are designated as Nη1 and Nη2 by
> their relationship (cis or trans, respectively) to Cδ. The hydrogen
> atoms of the side-chain -NH2 groups of Asn, Gln, and Arg are
> distinguished by numbers (1 or 2) on the basis of their relationship
> (cis or trans, respectively) to the heavy atom three bonds closer to
> the main chain (Cβ for Asn, Cγ for Gln, Nε for Arg)."

So BranchIndex is **derived from a physical configuration**: CIP-clockwise
for sp3, χ2-magnitude for aromatic, cis/trans for E/Z planar groups.

### What shift predictors actually use

SHIFTX2 has separate models for HB2 vs HB3 (`han-2011...text-4:23-25`):
correlations 0.9817 vs 0.9785, RMSDs 0.142 vs 0.156 — they ARE measurably
different chemically. For Val Cγ1 vs Cγ2 (`text-4:17-18`): correlations
0.8851 vs 0.8131 — even more different. For Ile Cγ1/Cγ2 + Cδ1: each
needs its own slot.

Sahakyan 2011 explicitly fits separate models for Leu Hδ1 vs Hδ2, Val
Hγ1 vs Hγ2, Ile Hγ2 (only one γ-methyl in Ile because Cγ1 is methylene)
and Ile Hδ1 (`sahakyan-2011...text-2:621-672`). The diastereotopic
chemical-shift differences encode rotamer information — Mulder 2009's
ΔδC(δ1−δ2) = −5 + 10·p_trans rotamer relation cited in
`/shared/2026Thesis/nmr-shielding/references-meta/sahakyan-2011-methyl-chemical-shifts-proteins-summary-qwen-test.txt:7`
puts a number on it.

### Most-granular reasonable enum / shape

```
enum class BranchIndex : uint8_t {
    NotApplicable = 0,   // atom does not need a branch disambiguator
    Branch1       = 1,
    Branch2       = 2,
    Branch3       = 3,   // for cases like ARG Hη11/η12/η21/η22 — see below
};
```

**Caveat for Arg side-chain hydrogens**: `markley-1998...text-2:367-371`
says η-hydrogens carry **two** numbers. So Hη11 has nitrogen-branch=1
and H-branch=1 within that nitrogen. This is a **two-level branch
hierarchy**, not a single index. The cleanest representation is:

```
struct BranchAddress {
    uint8_t  outer;   // which Nη (Arg), or which methyl carbon (Val/Leu)
    uint8_t  inner;   // which H within that group
};
```

Outer=1 for Nη1, inner=1 for the proton with E geometry to the heavy
atom three bonds in (per Markley convention). For 99% of atoms only
one of the two indices is non-zero; for Arg side-chain Hη atoms both
are.

### Open questions / known controversies

**Ring-flip ambiguity in Phe/Tyr**: from `markley-1998...text-3:421-441`,
when ring flip is slow on the chemical-shift time scale, signals from δ-
and ε-atoms can be non-equivalent, and the |χ2|-priority rule that
assigns 1 vs 2 requires knowing χ2 — which often is not known. The
literature handles this with the slash rule (`Hδ1/ε1/δ2/ε2`). For the
substrate, **BranchIndex is a static topology field** — it gets the
canonical Markley value computed from the standard residue topology,
NOT from the live conformation. Ring-flip is a conformation phenomenon
that should not flip the substrate label.

For **Arg cis/trans η** and **Asn/Gln/Arg amide hydrogen E/Z**: these
also depend on a definition of cis/trans (see Markley) but the topology
field stores the canonical IUPAC convention which is itself a
configuration choice. ProtonationDetectionResult or AmideOrientationResult
on the conformation side handles real-time E/Z (Asn/Gln amide flips are
a thing — see Word 1999 cited in the SHIFTX2 paper).

---

## Field 3: DiastereotopicIndex (HB2 vs HB3 for prochiral methylene)

### What it represents

The 2/3 designation on prochiral methylene hydrogens (HB2/HB3, HG2/HG3,
HD2/HD3, HE2/HE3 plus glycine HA2/HA3) and on certain ribose hydrogens
(H5'/H5'' etc.). This is **distinct** from `BranchIndex` because:

1. The **rule for assigning 2 vs 3** is the explicit "X-toward-the-viewer,
   Y/Z2/Z3 clockwise" convention from Markley Figure 1 caption.
2. The atoms involved are **hydrogens**, not heavy atoms.
3. The **mapping to pro-R/pro-S** is residue-specific (next field).

### Markley / IUPAC convention — full rule

From `markley-1998...text-2:357-372` (Figure 1 caption, quoted in
Field 2 above): if sighting down X→C with X toward the viewer, Z' (the
duplicate of Z) takes precedence as if it were the heavier isotope.
Then Y, Z', Z (projecting away) clockwise → Z' is pro-R, Z is pro-S
(see Markley §28 / `markley-1998...text-7:155-165`).

Specifically for amino acid Hβ atoms:
- **Sight from Cα toward Cβ**.
- The two diastereotopic Hβ atoms are Hβ2 and Hβ3.
- **Hβ3 = pro-R, Hβ2 = pro-S** is the alternation that Markley's rule
  produces for an L-amino acid with a heavy sidechain branch (because
  the branched sidechain heavy atom takes Y priority).

But — and this is the per-residue subtlety — the rule produces
**different alternations** when the sidechain branches differently or
when there's no β-heavy substituent (Gly). For **glycine HA**:
- Sight from N toward Cα.
- HA2 and HA3 are diastereotopic.
- `markley-1998...text-2:78-82`: Figure 1 marks HA2 explicitly as (R)
  and HA3 as (S). So for **Gly: HA2 = pro-R, HA3 = pro-S** (opposite
  of the typical β alternation).

For **Pro Hδ**: `markley-1998...text-2:79-82`: Hδ3 marked (R), Hδ2
unmarked. So for **Pro: Hδ3 = pro-R, Hδ2 = pro-S** (matches the
typical β alternation).

The general rule: **always derive from CIP, never assume an alternation**.
A topology-field initializer should run CIP at construction or read from
a curated per-residue table that has been CIP-checked.

### Curated table — diastereotopic protons in standard amino acids

This list is built from `markley-1998...text-2:Figure 1` annotations and
the canonical IUPAC convention. The (R) marks in Figure 1 encode pro-R.

| Residue | Atom positions | pro-R | pro-S |
|---------|---------------|-------|-------|
| Gly     | HA2, HA3      | HA2   | HA3   |
| Ala     | (no diastereotopic Hβ — methyl is symmetric) | – | – |
| Val     | Cγ1, Cγ2      | Cγ1 (R per Fig 1) | Cγ2 |
| Val     | Hβ — single H, no diastereotopy | – | – |
| Leu     | Cδ1, Cδ2      | Cδ1 (R per Fig 1) | Cδ2 |
| Leu     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |
| Leu     | Hγ — single H | – | – |
| Ile     | Hβ — single H | – | – |
| Ile     | Hγ12, Hγ13    | Hγ13 (R) | Hγ12 (S) |
| Ile     | Cγ2 (methyl), Cδ1 (methyl) — chain branch, see Field 2 | – | – |
| Pro     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |
| Pro     | Hγ2, Hγ3      | Hγ3 (R) | Hγ2 (S) |
| Pro     | Hδ2, Hδ3      | Hδ3 (R per Fig 1) | Hδ2 (S) |
| Ser     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |
| Thr     | Hβ — single H | – | – |
| Cys     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |
| Met     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |
| Met     | Hγ2, Hγ3      | Hγ3 (R) | Hγ2 (S) |
| Asp     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |
| Asn     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |
| Asn     | Hδ21, Hδ22 (E/Z, NOT pro-R/S — planar group) | – | – |
| Glu     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |
| Glu     | Hγ2, Hγ3      | Hγ3 (R) | Hγ2 (S) |
| Gln     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |
| Gln     | Hγ2, Hγ3      | Hγ3 (R) | Hγ2 (S) |
| Gln     | Hε21, Hε22 (E/Z) | – | – |
| Lys     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |
| Lys     | Hγ2, Hγ3      | Hγ3 (R) | Hγ2 (S) |
| Lys     | Hδ2, Hδ3      | Hδ3 (R) | Hδ2 (S) |
| Lys     | Hε2, Hε3      | Hε3 (R) | Hε2 (S) |
| Arg     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |
| Arg     | Hγ2, Hγ3      | Hγ3 (R) | Hγ2 (S) |
| Arg     | Hδ2, Hδ3      | Hδ3 (R) | Hδ2 (S) |
| Arg     | Hη11, Hη12, Hη21, Hη22 (E/Z, not pro-R/S — guanidinium plane) | – | – |
| His     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |
| His     | Hδ1 (NH proton on imidazole) — not prochiral | – | – |
| Phe     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |
| Phe     | ring atoms 1/2 — see RingPosition | – | – |
| Tyr     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |
| Tyr     | ring atoms 1/2 + OH proton — see RingPosition | – | – |
| Trp     | Hβ2, Hβ3      | Hβ3 (R) | Hβ2 (S) |

**Authoritative source for this table**: Markley 1998 Figure 1 + caption
(`text-2:357-372`). Cross-check via BMRB nomenclature page (which maps
HB1/HB2 in some conventions to HB2/HB3 in IUPAC); see <https://bmrb.io/referenc/nomenclature/>:
the BMRB explicitly notes "If one or two protons are missing in a methyl
or amino group...the protons will be numbered starting from 1" — so
PDB files with HB1 only have it because someone ran a missing-H prep.
In the IUPAC convention, complete methylenes are HB2/HB3, never HB1/HB2.

The user-flagged warning ("agents have to verify pro-R/pro-S in the next
pass" per memory `project_iupac_topology_landed_20260426.md`) applies
here. **Don't trust this table without independent CIP verification per
residue**. The table is consistent with Markley Figure 1 marks but those
marks are 1998-vintage and small typos in the original could propagate.

### What shift predictors actually use

- **SHIFTX2**: separate models for HB2 vs HB3 that achieve r=0.9817 vs
  r=0.9785 — they are CHEMICALLY distinct at the shift level, not just
  the label level (`han-2011...text-4:23-25`).
- **Sahakyan 2011 CH3Shift**: filters out non-stereospecifically-assigned
  Val/Leu records during database curation
  (`sahakyan-2011...text-1:127-141`): "we included only chemical shift
  entries with stereospecific assignment for Val and Leu residues. Cases
  for which chemical shifts were flagged as stereospecifically assigned
  but no difference between the two methyl chemical shifts were
  discarded." — i.e. the diastereotopic distinction IS the signal.

### Most-granular reasonable enum / shape

```
enum class DiastereotopicIndex : uint8_t {
    NotProchiral = 0,
    Position2    = 2,
    Position3    = 3,
};
```

This is the IUPAC numeric label (2 vs 3). Mapping 2→ProS / 3→ProR is
done by `ProchiralStereo` (Field 4); they are separate fields because:
- DiastereotopicIndex is a topology label (which atom IS this).
- ProchiralStereo is a chirality classification (R or S).
- Most residues alternate (3=R, 2=S) but Gly inverts (2=R, 3=S).

### Open questions / known controversies

**Stereospecific-assignment rate in BMRB is incomplete** —
`markley-1998...text-3:402-419` (the slash-rule examples) shows that
when assignments aren't stereospecific, the convention is to write
`Hβ2/β3 = 3.12/2.44 ppm` joined by slash. The substrate field is
deterministic per-residue (always 2 or 3 by IUPAC) but the experimental
shift table may carry ambiguity that flows through calibration.

**HA2/HA3 in Gly** — the inversion described above means generic code
that assumes "3=R" will mislabel Gly. Per-residue tables are unavoidable.

---

## Field 4: ProchiralStereo (pro-R / pro-S per Markley 1998)

### What it represents

The Cahn-Ingold-Prelog R/S designation of a prochiral atom. See
`markley-1998...text-7:148-165`:

> "The 1983 recommendations (ref. 21) specified the use of numbers to
> designate atoms and the use of the R/S rules of Cahn, Ingold and Prelog
> (refs. 25 and 26, see also refs. 27–29) to define stereochemistry.
> According to the R/S rules, a tetrahedral carbon X with prochiral
> substituents A > B > C = C' (where the priority is according to
> decreasing atomic mass and C and C' have equivalent mass and are
> labeled arbitrarily) is labeled as follows: one sights down the A–X
> axis and assumes that C' takes precedence over C (as it would if it
> were replaced, for example, with a heavier isotope). Then, if B, C',
> and C (projecting away from the viewer) are clockwise, C' is labeled
> as pro-R and C is labeled as pro-S; if they are counterclockwise,
> C' is labeled as pro-S and C is labeled as pro-R."

This is the formal definition. The application of it to the 20 standard
amino acids gives the per-residue table in Field 3.

### Markley / IUPAC convention — full enum

```
enum class ProchiralStereo : uint8_t {
    NotProchiral = 0,
    ProR         = 1,   // CIP rectus: Cahn-Ingold-Prelog right-handed
    ProS         = 2,   // CIP sinister: Cahn-Ingold-Prelog left-handed
    Unassigned   = 3,   // residue is prochiral but no CIP determination
                        // available — should not happen for std amino
                        // acids where the convention is in Markley
};
```

This corresponds 1:1 with the CCD's `_chem_comp_atom.pdbx_stereo_config`
field which has allowed values `R / S / N` (none) — see
<https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_atom.pdbx_stereo_config.html>.

### Mapping table — DiastereotopicIndex ↔ ProchiralStereo

For the standard amino acids:

- **β-onward methylenes (most residues)**: index=3 → ProR; index=2 → ProS.
- **Gly Hα methylene**: index=2 → ProR; index=3 → ProS. **Inverted.**
- **Pro Hδ methylene**: index=3 → ProR; index=2 → ProS (per Markley
  Figure 1 marker on Hδ3 — same alternation as the typical β-onward).
- **Branched heavy atoms** (Val Cγ1/Cγ2, Leu Cδ1/Cδ2, Ile Cγ2/Cδ1):
  Cγ1/Cδ1 = ProR, Cγ2/Cδ2 = ProS per Markley Figure 1 marks.

### What shift predictors actually use

The chemical-shift correlate of pro-R vs pro-S is real and substantial.
Mulder 2009 (cited in `sahakyan-2011...text-1:71-72`) gave:
`ΔδC(δ1-δ2) = -5 + 10·p_trans` for leucine.

This means knowing which methyl is pro-R and which is pro-S is
**necessary to interpret the chemical-shift difference as a rotamer
population**. Without ProchiralStereo, the per-conformation shift
information cannot be tied to the residue rotamer state.

### Most-granular reasonable enum / shape

The 4-value enum above. For the topology substrate, the Unassigned
state should never appear in canonical residues but is needed for
non-standard residues (DFT model compounds, modified residues).

The CCD / RDKit / Open Babel agreement makes this field portable —
all three use R/S/N.

### Open questions / known controversies

- **Naming-convention divergence in legacy code**: the BMRB nomenclature
  page (<https://bmrb.io/referenc/nomenclature/>) explicitly notes:
  "C-terminal oxygens are labelled O' and O''. If protonated the proton
  will be labelled O'' because of it's lower priority" — i.e. CIP
  drives the prime-vs-double-prime convention. PDB files don't always
  follow this.
- **Markley 1998 Figure 1 has the (R) marks but not (S) marks** —
  unmarked diastereotopic atom is implicitly pro-S. A few annotations
  in Figure 1 are illegible from OCR (text-2:79-228 was garbled chunk).
  An independent verification pass against e.g. Hanson 1966 (Markley
  ref. 23) or recent IUPAC compilations is wise.
- **Variant residues** (HID, HIE, HIP for histidine; CYX for disulfide
  cysteine) — the CIP designations of HD methylenes etc. can in
  principle change with the variant if the priority ranking changes.
  In practice they don't (the side-chain heavy atoms outrank the
  protonation difference), but a careful check per-variant is wise.

---

## Field 5: PlanarStereo (sp2 / planar group classification + the conformation question)

### What it represents

This is the field with the **substrate-vs-conformation question** the
user flagged ("we need conformation, conformation is not just binary but
angle"). The right answer is: **planar groups have a topology side
(which atoms are in the planar group, what the canonical/idealized
configuration is) and a conformation side (the actual angle in this
particular conformation, including deviations and flips).**

The substrate field carries the **topology side**. The conformation
side belongs on `PlanarGroupResult` or similar per-frame calculator
output.

### Markley / IUPAC convention — full enum

The places in Markley 1998 where E/Z or planar-stereo annotations appear:

1. **Peptide bond (ω torsion)**: `markley-1998...text-3:200-260` —
   the ω torsion is defined as Cαi-1 → Ci-1' → Ni → Cαi. Idealized
   trans is 180°, cis is 0°. Markley Figure 7 (`text-4:140-256`)
   gives the Klyne-Prelog conventions (±sp = synperiplanar, ±sc =
   synclinal, ±ap = antiperiplanar, ±ac = anticlinal) and the
   spectroscopist convention (cis = 0°±30°, trans = 180°±30°,
   gauche± = ±60°±30°). Real proteins are 99.97% trans for non-Pro
   peptides and ~95% trans for Xaa-Pro (Craveur 2013, cited via
   <https://pubmed.ncbi.nlm.nih.gov/23728840/>).
2. **Asn/Gln amide H E/Z**: `markley-1998...text-2:367-371` (quoted
   above) — H21/H22 distinguished by cis/trans relationship to a heavy
   atom three bonds closer to the main chain. So (E) and (Z) are written
   in Figure 1 — see `text-2:80-228` for Asn Hδ22(Z), Hδ21(E) and
   Gln Hε21(E), Hε22(Z) etc. The configuration is **implied** by IUPAC
   atom labels.
3. **Arg side-chain Nη / Hη E/Z**: `text-2:362-371` — Nη1 cis to Cδ,
   Nη2 trans. Each H on each N is further marked E or Z relative to
   the heavy atom. So Hη11 vs Hη12 differ by E/Z assignment.
4. **Phe/Tyr ring 1/2**: `text-2:368-371` — Cδ1/Cε1 = smaller |χ2|,
   so the ring orientation in the conformation IS the assignment
   criterion. Slow ring flip means the assignment is real chemistry;
   fast ring flip means the assignment is degenerate at NMR
   timescales (the two atoms have identical observed shift).
5. **CCD `_chem_comp_bond.pdbx_stereo_config`**: allowed values
   E / Z / N (see <https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_bond.pdbx_stereo_config.html>).
   This is for **double bonds** (C=N in amide, C=N in guanidinium).

### The substrate-vs-conformation split for this field

The user's question — "we need conformation, conformation is not just
binary but angle" — gets a careful answer. Here's the split that
matches both Markley and the way reference predictors treat the data:

**Substrate (topology) side — what `LegacyAmberTopology` should hold:**
- Membership: which atoms are in a planar/sp2 group. Asn Cβ-Cγ-Oδ1-Nδ2-Hδ21-Hδ22
  is the Asn amide planar group. Phe Cγ-Cδ1-Cδ2-Cε1-Cε2-Cζ ring atoms.
  His five-membered ring atoms. Arg guanidinium planar group. **This is
  determined by AMBER force-field types or by the residue chemistry —
  unchanged across conformations.**
- Canonical / idealized configuration: trans peptide bond is the IUPAC
  default; Asn/Gln/Arg E/Z labels come from Markley Figure 1; ring 1/2
  numbering is a convention choice. **Topology fixes which atom is
  E and which is Z.**
- Group type: peptide ω, amide (Asn/Gln side-chain), guanidinium (Arg),
  imidazole (His), aromatic 6-ring (Phe/Tyr/Trp benzene part), aromatic
  5-ring (Trp pyrrole part, His imidazole), carboxyl (Asp/Glu), Tyr OH
  on aromatic.

**Conformation side — what `PlanarGroupResult` / `BondTorsionResult`
should hold per-frame:**
- Actual ω angle in radians (not just trans/cis classification).
- Planarity-deviation angle of sp2 N in Asn/Gln side chain (small
  pyramidalization).
- Phe/Tyr χ2 angle (which determines whether the Cδ1/Cε1 labels match
  the live geometry — see Markley caption).
- Ring puckering for non-aromatic ring (Pro: UP / DOWN per
  `markley-1998...text-3:295-315`; ribose: pseudorotation phase).

**This split is consistent with the `feedback_object_model_scope_discipline`
memory** — invariant facts on the substrate, conformation-dependent on
per-frame results.

### What shift predictors actually use

- **ProCS15 (Larsen 2015)**: the backbone term σ_BB^i is parameterized
  on a tripeptide grid scanning φ, ψ, χ1, χ2, χ3, χ4. **The ω angle
  is FIXED at 180°** (`larsen-2015...text-3:65-67`):

  > "The ω dihedral angle was fixed at 180°. The φ/ψ backbone angles on
  > the N and C-termini alanine residues were fixed at -140° and 120°
  > corresponding to typical β-sheet residue backbone angles."

  This is a strong simplifying assumption. For the 99.97% of non-Pro
  residues with trans ω, it's fine; for cis-Pro it would be wrong, and
  ProCS15 doesn't distinguish.
- **SHIFTX2 (Han 2011)**: includes ω as a feature with weights 13.9%
  for C', 10.4% for Cα, 5.7% for Cβ, 5.3% for HN, 3.8% for Hα, 7.1%
  for N (`han-2011...text-4:147-157` Table 5). **ω is the third- or
  fourth-most influential feature for backbone shifts.** It is
  treated explicitly as a continuous feature, not a binary cis/trans.
  And SHIFTX2 includes ω_{i-1}, ω_i, ω_{i+1}, ω_{i-2} all separately
  in the ML model — so neighborhood ω values matter.
- **ω varies in real proteins**: a "trans" ω is rarely exactly 180°
  — typical fluctuations span ±10-15°. The ML model captures these
  fluctuations. Cis ω (around 0°) is genuinely rare except at
  Xaa-Pro junctions.

So PlanarStereo on the topology layer should NOT lock ω to its
idealized value; it should record only the **canonical assignment**
(trans is the default for substrate analysis), while leaving the
**actual angle** to the per-frame conformation result.

### Most-granular reasonable enum / shape

Two pieces. First, the classification of which planar-group an atom
belongs to — this is topology:

```
enum class PlanarGroupKind : uint8_t {
    None              = 0,
    PeptideAmide      = 1,   // Cα(i-1)-C(i-1)-N(i)-Cα(i) backbone
    SidechainAmide    = 2,   // Asn/Gln side chain Cγ/Cδ-O/N
    Guanidinium       = 3,   // Arg Cζ + Nε/Nη1/Nη2
    Imidazole         = 4,   // His side chain
    Aromatic6Ring     = 5,   // Phe + Tyr + Trp benzene component
    Aromatic5Ring     = 6,   // His imidazole + Trp pyrrole component
    Carboxylate       = 7,   // Asp/Glu side-chain CO-O
    AromaticHydroxyl  = 8,   // Tyr Cζ-OH (the OH is in the ring plane)
};
```

Second, the canonical E/Z designation if it's an E/Z-distinguishable
atom — this is the configuration label, independent of live angle:

```
enum class PlanarStereo : uint8_t {
    NotApplicable = 0,
    E             = 1,   // entgegen / opposite
    Z             = 2,   // zusammen / together
    Unspecified   = 3,
};
```

This 4-value enum aligns with CCD `_chem_comp_bond.pdbx_stereo_config`'s
{E, Z, N} so it's portable.

The **angle** is conformation, not topology. Per-frame conformation
result owns it. To answer the user's "we need angle" — yes, but the
angle goes on the conformation, not on `LegacyAmberTopology`. The
substrate carries:
- Which atoms are in the planar group (topology).
- The canonical E/Z label per IUPAC (topology — defined by Markley
  Figure 1, doesn't change per frame).
- The PlanarGroupKind taxonomy.

### Open questions / known controversies

- **Asn/Gln amide flip ambiguity**: the Word 1999 paper (cited by
  SHIFTX2 — `han-2011...text-1:88` ref. and `han-2011...text-3:155-157`)
  Word et al. 1999 J. Mol. Biol. 285:1735 — "Asparagine and glutamine:
  using hydrogen atom contacts in the choice of sidechain amide
  orientation" is the canonical paper showing that ~30% of Asn/Gln
  side-chain amide groups in PDB files have **the flipped (wrong)
  assignment** of Oδ/Nδ vs Hδ21/Hδ22 because crystallography can't
  distinguish O from N at typical resolution. **Reduce/Word's flip is
  a chemistry change**, not a conformation change — it changes which
  atom IS Hδ21 vs Hδ22. So the substrate label can be wrong if the
  input PDB is mis-assigned.
  This is a known data-quality issue and is the kind of thing AIMNet2-
  in-the-loop should catch (chemistry + geometry together).
- **Cis-trans peptide bonds at Pro junctions are stable on the NMR
  time scale** (`Wikipedia: "Peptide bond" + Craveur 2013`). NMR can
  resolve cis-Pro from trans-Pro as separate populations. The substrate
  layer should support both (e.g. a peptide-bond-i field with E or Z),
  not assume trans.
- **Tyr OH proton in the ring plane vs out of plane**: the Tyr OH is
  typically described as in-plane with the ring, but Cζ-O-H can rotate.
  Most ring-current models (Sahakyan 2013, Larsen 2015) treat the ring
  as the planar group and ignore the OH-rotation as a separate planar
  question.

---

## Field 6: PseudoatomClass (Markley pseudoatom taxonomy)

### What it represents

The user's framing was "Markley pseudoatom taxonomy: M, Q, X, Y, S — we
need the lot." After reading the corpus carefully, the actual Markley
taxonomy is **smaller** than that — only **Q, M, R** (and one residue-
specific letter, A, used only for the Gly α-methylene). X, Y, S, T are
**not** in Markley 1998. They appear in some other contexts (XPLOR
nomenclature, CARA), but are NOT IUPAC pseudoatoms.

### Markley / IUPAC convention — full enum from Table 1

`markley-1998...text-3:181-242` — Table 1 of the paper. Verbatim:

| Residue | Pseudoatom | What it represents |
|---------|-----------|--------------------|
| Gly | QA | α-methylene |
| Ala | MB | β-methyl |
| Val | MG1, MG2 | γ1-, γ2-methyl |
| Val | QG | All six γ-methyl Hs |
| Ile | MG, MD | γ2-, δ1-methyl |
| Ile | QG | γ1-methylene |
| Leu | MD1, MD2 | δ1-, δ2-methyl |
| Leu | QB | β-methylene |
| Leu | QD | All six δ-methyl Hs |
| Pro | QB, QG, QD | β-, γ-, δ-methylene |
| Ser, Asp, Cys, His, Trp | QB | β-methylene |
| Thr | MG | γ2-methyl |
| Asn | QB | β-methylene |
| Asn | QD | δ2-amino (the side-chain amide H pair) |
| Glu | QB, QG | β-, γ-methylene |
| Gln | QB, QG | β-, γ-methylene |
| Gln | QE | ε2-amido |
| Lys | QB, QG, QD, QE | β-, γ-, δ-, ε-methylene |
| Lys | QZ | ζ-amino (NH3+ trio) |
| Arg | QB, QG, QD | β-, γ-, δ-methylene |
| Arg | QH1, QH2 | η11+η12, η21+η22 (each Nη's H pair) |
| Arg | QH | All four η-guanidino Hs |
| Met | QB, QG | β-, γ-methylene |
| Met | ME | ε-methyl |
| Phe, Tyr | QB | β-methylene |
| Phe, Tyr | QD, QE | δ1+δ2 ring, ε1+ε2 ring |
| Phe, Tyr | QR | All ring Hs |

Plus from `markley-1998...text-3:184-189`:

> "M describes the location of methyl groups, and Q is used in all
> other situations [...] R is used to denote a ring."

So the taxonomy is:

- **M = methyl** (3 equivalent Hs averaged)
- **Q = unassigned/group of equivalent Hs** (general fallback)
- **R = ring** (used as the locant index for ring-atom pseudoatoms,
  e.g. QR = Q across all ring atoms)

Plus the specific letter (B/G/D/E/Z/H/A) that identifies the locant of
the pseudoatom. So `MD1` = methyl-at-δ-position-1, `QB` = group-at-β,
`QR` = group-at-ring.

### What shift predictors actually use

- The Markley pseudoatom labels are used in **structure constraints**
  (NOE distance constraints), not in chemical-shift prediction directly.
  Predictors operate on individual atoms.
- However, when a stereo-specific assignment is missing for, say, the
  Hβ pair, the shift assignment table reports them with the slash rule
  `Hβ2/β3 = 3.12/2.44`. If both have the same observed shift,
  the QB pseudoatom can be associated with that single shift.

### Most-granular reasonable enum / shape

For a **per-atom topology field** on `LegacyAmberTopology`, the
question is: "for the standard residue that this atom belongs to, what
pseudoatom group(s) does this atom belong to?" An atom can belong to
multiple pseudoatoms (e.g., an Hβ in Lys belongs to QB; a δ1 methyl H
in Leu belongs to MD1 AND QD).

The cleanest representation:

```
enum class PseudoatomKind : uint8_t {
    None         = 0,
    M            = 1,   // methyl (3 Hs)
    Q            = 2,   // group/methylene (typically 2 Hs)
    R            = 3,   // ring (used as a locant suffix, not a pseudoatom by itself)
};

// Per atom: which pseudoatom group(s) include this atom?
struct PseudoatomMembership {
    uint8_t  primary_kind;        // M, Q, or 0
    uint8_t  primary_locant;      // β=1, γ=2, δ=3, ε=4, ζ=5, η=6, A=7, R=8
    uint8_t  primary_branch;      // 1 or 2 for MD1/MD2 etc., 0 if none
    bool     in_super_group;      // true if also part of QG/QD/QH/QR aggregate
};
```

### Open questions / known controversies

- **Wüthrich 1983 (ref. 30 in Markley)** introduced the underlying
  pseudoatom formalism. Markley 1998 Table 1 uses Wüthrich's set
  almost unchanged. The phrase "all pseudoatoms needed to describe
  proteins and nucleic acids" (`text-3:181-185`) suggests this is
  considered closed for standard residues.
- **The "X, Y, S, T" pseudoatom letters mentioned in the task brief
  are NOT in Markley 1998.** Looking at the BMRB nomenclature pages
  and the CARA wiki (which couldn't load), they appear in:
  - **XPLOR pseudoatom convention**: HB*, HB#, HB%, HB+ — these are
    shorthand for `2/3 split` or `complete group` rather than separate
    letter names. Per the BMRB conversion table summary, XPLOR uses
    suffix characters (*, #, %, +) attached to atomic labels.
  - **CARA (Computer-Aided Resonance Assignment)**: defines its own
    pseudoatom set adapted from Wüthrich; some letters (X, Y) appear
    as "any atom in this position" placeholders, not as canonical
    pseudoatoms.
  - **DIANA / older tooling**: legacy notation for ambiguous cases.

  None of these are IUPAC. **The Markley pseudoatom set (M, Q, R) is
  the authoritative IUPAC taxonomy. X/Y/S/T should not be added to
  the substrate field.**
- **For non-standard residues** (modified amino acids, ligands), the
  CCD has its own atom names but NOT a parallel pseudoatom convention
  beyond gross M/Q. New pseudoatoms should be defined per the Markley
  rule (M for methyl, Q for everything else, plus a locant) when
  needed.

---

## Field 7: PolarHKind (functional-group classification of polar hydrogens)

### What it represents

The classification of a hydrogen atom by what functional group it sits
on, restricted to polar hydrogens (those that exchange or hydrogen-bond).
This is **not** the same as just "is this a polar H" — the chemistry
matters. An amide NH and a hydroxyl OH and an indole NH all behave
differently in NMR.

### IUPAC / chemistry convention

There is no single IUPAC table. The relevant physical chemistry distinguishes:

1. **Backbone amide N-H** (HN/H) — the canonical hydrogen-bond donor of
   the polypeptide. Chemical shift typically 6-10 ppm, very sensitive
   to H-bond geometry and secondary structure.
2. **Side-chain amide N-H**:
   - Asn Hδ21, Hδ22 — primary amide on Cγ.
   - Gln Hε21, Hε22 — primary amide on Cδ.
   - Trp Hε1 — indole secondary amide.
3. **Side-chain ammonium N-H** (Lys Hζ1/2/3, N-terminal H1/2/3) — fully
   exchanged with water typically; rarely observable.
4. **Guanidinium N-H** (Arg Hε, Hη11, Hη12, Hη21, Hη22) — extensively
   delocalized; intermediate exchange.
5. **Imidazolium / imidazole N-H** (His Hδ1 in HID, Hε2 in HIE; both in
   HIP/HIH) — sensitive to protonation state.
6. **Carboxyl O-H** (Asp Hδ2 in ASH, Glu Hε2 in GLH, C-terminal OH).
   Very acidic — typically not observed in solution NMR at neutral pH.
7. **Hydroxyl O-H** (Ser Hγ, Thr Hγ1, Tyr Hη). Slow exchange, observable
   under some conditions.
8. **Thiol S-H** (Cys Hγ in reduced state, NOT in CYX disulfide). Slow
   exchange.

### What shift predictors actually use

- **SHIFTX2 (Han 2011)** — predicts amide HN explicitly with high
  accuracy (`han-2011...text-2:148`: r=0.9714 RMSD=0.171 ppm). It also
  has separate slots for HD21/HD22/HE21/HE22 (the side-chain amide Hs)
  and HH11/HH12/HH21/HH22 (Arg guanidinium Hs)
  (`han-2011...text-2:21-23`):

  > "Too few shifts were available for 9 side chain 15N atoms, eight 1H
  > atoms (HD21, HD22, HE21, HE22, HH11, HH12, HH21 and HH22) and four
  > carbon atoms (CE3, CH2, CZ2, CZ3)."

  i.e. SHIFTX2 wanted to predict each polar-H type separately but the
  training database was too small for these 8. This tells us the
  granularity matters chemically.
- **ProCS15 (Larsen 2015)** — predicts only backbone HN, with a
  hydrogen-bond term Δσ_HB^i parameterized on N-methylacetamide model
  systems (`larsen-2015...text-3:1-10`). The H-bond geometry is
  parameterized as (rOH, θ, ρ) — distance plus two angles.
- **Christensen 2011** (the predecessor to ProCS15) — focused
  exclusively on backbone amide HN and ring-current effects on it.
- **Yi-McDermott 2024 SSNMR study** —
  `/shared/2026Thesis/nmr-shielding/references-meta/yi-mcdermott-2024-temperature-shifts-conformational-dynamics-summary-qwen-test.txt:7`:
  finds that "Hydrogen-bond angle θ correlated strongly with ψ, whereas
  distance showed weak correlation with shifts" — this is for backbone
  HN. The angle dependence is not yet exploited in classical predictors.

### Most-granular reasonable enum / shape

```
enum class PolarHKind : uint8_t {
    NotPolar             = 0,    // C-H or aliphatic H
    BackboneAmide        = 1,    // HN/H — backbone N-H
    SidechainPrimaryAmide = 2,   // Asn Hδ2, Gln Hε2 — C(=O)-NH2
    IndoleNH             = 3,    // Trp Hε1 — secondary aromatic NH
    AmmoniumNH           = 4,    // Lys Hζ, N-term H1/2/3 — quaternary nitrogen
    GuanidiniumNH        = 5,    // Arg Hε, Hη11/12/21/22 — delocalized
    ImidazoleNH          = 6,    // His Hδ1 (HID) or Hε2 (HIE) or both (HIP)
    CarboxylOH           = 7,    // Asp Hδ2 (ASH), Glu Hε2 (GLH), C-term OH
    HydroxylOH_Aliphatic = 8,    // Ser Hγ, Thr Hγ1
    HydroxylOH_Aromatic  = 9,    // Tyr Hη
    ThiolSH              = 10,   // Cys Hγ (reduced only — not in CYX)
    OtherPolarH          = 11,   // catch-all for non-standard residues
};
```

This 12-class enum is **substantially finer** than just "polar / nonpolar"
or "donor / acceptor" — it carries the chemistry that drives shielding
behavior. SHIFTX2's Table 2 (`han-2011...text-4:1-39`) implicitly uses
roughly this granularity by training separate models per residue+atom.

### Open questions / known controversies

- **Protonation-state-dependent assignments**: a Cys H-γ exists only
  in reduced Cys (not in CYX disulfide); a His H-δ1 exists only in
  HID/HIP; a His H-ε2 exists only in HIE/HIP. The PolarHKind for an
  atom IS coupled to the protonation variant. The substrate field
  has to either (a) carry the variant-dependent value or (b) carry a
  variant-independent value with the variant indexed separately.
  Option (b) is cleaner: PolarHKind is determined by atom name +
  variant, computed from both at construction time.
- **ASH / GLH / CYM**: alternative protonation states with the
  proton present (ASH = Asp-protonated) need their own PolarHKind.
  See the AMBER variant index recently consolidated in
  `src/AminoAcidType.cpp::VariantIndexFromForceFieldName`.
- **N-terminus / C-terminus**: the ammonium N-H1/H2/H3 set on the
  N-terminus, the carboxyl O-H'' on the C-terminus, and the amine
  versions when neutralized — these are all PolarH variants of the
  base residue. Cap pseudo-residues (NHE, NME, ACE) introduce more.
  Currently scoped per `feedback_object_model_scope_discipline` —
  cap chemistry on the substrate.

---

## Field 8: RingPosition (ring-atom positional information)

### What it represents

The user said "we want all the specific ring-position information we
can get — hugely informative for NMR and shielding statistics." This
is about giving each ring atom a **position label that captures**:
- Which ring(s) it's in (Phe has 1 ring, Trp has 2 rings + perimeter).
- Where it sits in that ring (ipso/ortho/meta/para; α/β positions in
  five-rings; bridgehead in fused systems).
- Distance from heteroatom (in heterocycles).
- Heteroatom type if the position itself is a heteroatom.

### Survey of conventions

Different conventions exist. None is universal in NMR-shielding work,
but several are used.

#### Phe / Tyr (benzene-like)

The standard NMR convention places the ring along Cβ-Cγ-Cδ1/2-Cε1/2-Cζ.
Markley defines Cδ1/Cε1 (vs Cδ2/Cε2) by the smaller |χ2| convention
(`text-2:368-371`).

In the **chemistry / IUPAC organic-numbering** view, the ring positions
relative to Cγ (the **ipso** carbon — the one bonded to the rest of the
sidechain) are:
- ipso = Cγ
- ortho = Cδ1, Cδ2
- meta = Cε1, Cε2
- para = Cζ

(For Tyr, Cζ also carries the OH, making para the donor position.)

#### Trp (indole)

Indole has a five-membered pyrrole ring fused to a six-membered benzene
ring, sharing two atoms. Markley `text-2:230-358` (Figure 1, Trp
section) labels the atoms:

- Pyrrole part: Cβ(connector)-Cγ-Cδ1-Nε1-Cδ2(fusion)-Cε2(fusion)
  - Cδ1 is the only non-bridgehead pyrrole carbon
  - Nε1 carries the indole NH
  - Cδ2 and Cε2 are the bridgehead atoms
- Benzene part: Cδ2(fusion)-Cε2(fusion)-Cζ2-Cη2-Cζ3-Cε3
  - Cε3, Cζ3, Cη2, Cζ2 are the perimeter carbons

#### His (imidazole)

Five-membered ring with two nitrogens. Connectivity Cβ-Cγ-Nδ1-Cε1-Nε2-Cδ2.
Markley `text-2:163-225` shows: Cγ-Nδ1-Cε1-Nε2-Cδ2-Cγ (5-ring).
- Cγ = ipso (connector to backbone)
- Nδ1 = NH in HID, deprotonated in HIE (variant-dependent)
- Cε1 = sp2 C between two nitrogens
- Nε2 = NH in HIE, deprotonated in HID
- Cδ2 = the other sp2 ring carbon

#### Pro (pyrrolidine — non-aromatic)

Pro is special because it's saturated. Ring puckering DOWN/UP
(`markley-1998...text-3:295-315` and Figure 4 caption at text-4:1-30).

### What shift predictors actually use

- **ProCS15 (Larsen 2015)** ring-current term
  `larsen-2015...text-3:38-47`:

  > "ΔσRC^i = i·B·(1 − 3cos²θ) / r³. The model depends on the parameters
  > i, which is the side-chain-specific ring-current intensity relative
  > to benzene, B, which is a constant in the model, and the vector r,
  > which is the vector from the proton to the center of the aromatic
  > ring."

  So ProCS15 uses **per-ring intensity factors i** (one per ring type:
  Phe, Tyr, His, Trp 5-ring, Trp 6-ring) calibrated against
  Christensen-Sauer-Jensen 2011. The ring is treated as a single point
  dipole at the ring center; **individual ring-atom positions are not
  used as predictor features** in this term, only the ring itself.

- **Sahakyan 2011 CH3Shift** uses Haigh-Mallion ring-current geometric
  factors (`sahakyan-2011...text-2:1-50`) which DO use individual ring
  atom positions:

  > "where Sij is the algebraic (signed) triangle area formed by the O'
  > projection of the query point O onto the ring plane and the ring
  > atoms i and j. [...] The summation goes over all the adjacent ij
  > atom pairs forming the ring, that is over the number of bonds in
  > the conjugated ring."

  So each ring atom is a vertex in the Haigh-Mallion summation. The
  per-atom position is used, but **as a coordinate**, not as a
  classification label.

  Sahakyan 2011 also uses "rings" plurally for Trp:
  `sahakyan-2011...text-2:32-37`:

  > "For tryptophan residues, if one of the two rings satisfy the above
  > mentioned criterion, the second ring is included as well."

  So Trp's two rings are treated as separate ring objects, both
  contributing.

- **SHIFTX2 (Han 2011)** has ring-current feature weights of 11.5%
  for HN, 11.2% for Hα, and ~0% for backbone carbons
  (`han-2011...text-4:147-156` Table 5). Ring-current is the second-
  largest feature for HN and Hα. The implementation uses Case 1995
  parameterization (Haigh-Mallion). Individual ring atom positions
  are summed via the Haigh-Mallion formula.

- **Sahakyan-Vendruscolo 2013** RNA bases study extends ring-current
  modeling to nucleic acids, classifying inter-ring arrangements as
  ADJ (adjacent, stacked), SPT (spatial, diagonal) or HBD (hydrogen-
  bonded, coplanar) (`sahakyan-vendruscolo-2013...text-2:80-96`).
  This is **inter-ring relative geometry** rather than ring-position
  labels per se, but it captures the "what kind of ring environment
  is this atom in" question.

- **None of the current predictors use ipso/ortho/meta/para labels
  explicitly.** The ring-current physics is geometric (Haigh-Mallion
  triangle areas, Pople point-dipole). Position-label features are
  only implicit via per-atom-name training data.

But the **statistics within each ring position** are real:

- Phe Cδ1/Cδ2 (ortho) at ~131 ppm vs Cε1/Cε2 (meta) at ~129 ppm vs
  Cζ (para) at ~128 ppm in random coil. The differences are small
  but consistent.
- Tyr Cδ at ~133 ppm, Cε at ~118 ppm, Cζ at ~157 ppm — much bigger
  differences because of the OH.
- Trp Hε1 (indole NH) at ~10 ppm — completely separate from amide
  HN region (8 ppm).
- His Hδ2/Hε1 distinguishable; protonation state changes them by
  several ppm.

These are the kinds of "shielding statistics" the user wants stratified
by RingPosition.

### Most-granular reasonable enum / shape

The right thing is a **structured ring address**, not a flat enum.
A flat enum would have to enumerate every (ring-system, position)
combination, which is brittle.

```
enum class RingSystemKind : uint8_t {
    NotInRing       = 0,
    Benzene_Phe     = 1,    // Phe sidechain 6-ring
    Benzene_Tyr     = 2,    // Tyr sidechain 6-ring
    Imidazole_His   = 3,    // His sidechain 5-ring
    Indole_Trp_5    = 4,    // Trp pyrrole 5-ring
    Indole_Trp_6    = 5,    // Trp benzene 6-ring (fused with the pyrrole)
    Pyrrolidine_Pro = 6,    // Pro saturated 5-ring (non-aromatic)
};

enum class RingPositionLabel : uint8_t {
    NotInRing       = 0,
    Ipso            = 1,    // bonded to the rest of the sidechain (Cγ)
    Ortho1          = 2,    // ortho, branch-1
    Ortho2          = 3,    // ortho, branch-2
    Meta1           = 4,    // meta, branch-1
    Meta2           = 5,    // meta, branch-2
    Para            = 6,    // para
    PyrroleAlpha    = 7,    // five-ring α-position (next to N)
    PyrroleBeta     = 8,    // five-ring β-position (Trp Cδ1, His Cδ2)
    BridgeFusion    = 9,    // shared atom in fused-ring system (Trp Cδ2/Cε2)
    Heteroatom_NH   = 10,   // ring N with attached H (Trp Nε1, His Nδ1/HID, His Nε2/HIE)
    Heteroatom_NoH  = 11,   // ring N without H (His Nε2/HID, His Nδ1/HIE)
    Heteroatom_OH   = 12,   // ring O with attached H (none in standard 20)
    Saturated       = 13,   // Pro Cβ/Cγ/Cδ — sp3 ring carbon
};

struct RingMembership {
    RingSystemKind   ring;
    RingPositionLabel position;
    uint8_t          ring_size;     // 5 or 6
    bool             aromatic;      // false for Pro pyrrolidine
    bool             planar;        // tied to aromatic + His charge state
    uint8_t          n_heteroatoms; // 0 (Phe, Tyr-mostly), 1 (Trp pyrrole), 2 (His)
};

// For atoms in fused-ring systems (Trp), an atom may belong to BOTH rings:
struct RingPosition {
    RingMembership primary;     // the smaller ring if shared
    RingMembership secondary;   // the second ring; .ring == NotInRing if not shared
};
```

### Specific assignments per residue

#### Phe (Benzene_Phe, 6-ring)
- Cγ = Ipso, Cδ1 = Ortho1, Cδ2 = Ortho2, Cε1 = Meta1, Cε2 = Meta2, Cζ = Para
- Hδ1, Hδ2 = on Ortho1, Ortho2; Hε1, Hε2 = on Meta1, Meta2; Hζ = on Para

#### Tyr (Benzene_Tyr, 6-ring)
- Cγ = Ipso, Cδ1 = Ortho1, Cδ2 = Ortho2, Cε1 = Meta1, Cε2 = Meta2, Cζ = Para
- Tyr Cζ has the hydroxyl. The OH oxygen is **not in the ring** but
  on the para position.

#### Trp (fused: Indole_Trp_5 + Indole_Trp_6)
Pyrrole ring (5):
- Cγ = Ipso, Cδ1 = PyrroleBeta, Nε1 = Heteroatom_NH, Cε2 = BridgeFusion,
  Cδ2 = BridgeFusion
- Hδ1 on Cδ1, Hε1 on Nε1
Benzene ring (6):
- Cδ2 = BridgeFusion, Cε2 = BridgeFusion, Cζ2 = Ortho2, Cη2 = Meta2,
  Cζ3 = Meta1, Cε3 = Ortho1
- Hε3, Hζ3, Hη2, Hζ2 on Cε3, Cζ3, Cη2, Cζ2

#### His (Imidazole_His, 5-ring)
- Cγ = Ipso, Nδ1 = Heteroatom_NH (HID/HIP) or Heteroatom_NoH (HIE),
  Cε1 = PyrroleAlpha (between the two N's), Nε2 = Heteroatom_NH (HIE/HIP)
  or Heteroatom_NoH (HID), Cδ2 = PyrroleBeta
- Hδ1 on Nδ1 (only in HID/HIP), Hε1 on Cε1, Hε2 on Nε2 (only in HIE/HIP),
  Hδ2 on Cδ2

#### Pro (Pyrrolidine_Pro, 5-ring, non-aromatic, sp3)
- Cα, Cβ, Cγ, Cδ, N — all in the ring (Cα is the only one shared with
  the backbone)
- Position labels are Saturated for Cβ, Cγ, Cδ
- N is a Heteroatom (no H in the normal Pro, since Pro's nitrogen IS
  the secondary amine — already part of the ring)

### What this gives you for shielding statistics

If RingPosition is recorded per atom, you can stratify per-element
shielding distributions by:
- ring system (Phe/Tyr/His/Trp-5/Trp-6/Pro)
- ring position (ipso/ortho/meta/para/bridgehead/heteroatom)
- heteroatom presence and type
- protonation variant (HID vs HIE vs HIP for His)

This goes substantially beyond what current shift predictors do as
explicit categorical features, but matches what the underlying physics
expresses (ring currents are geometric; ring-position label gives
structural context).

### Open questions / known controversies

- **Trp 5-ring + 6-ring overlap**: Cδ2 and Cε2 are in both rings
  (the bridgehead). The structure-of-ring-membership has to handle
  multi-ring atoms. The `primary` / `secondary` shape above does this.
- **Aromatic vs non-aromatic discrimination**: His protonation state
  affects aromaticity. HID/HIE imidazole is aromatic (6 π electrons,
  Hückel 4n+2). HIP imidazolium has 6 π electrons distributed over
  more atoms, still aromatic. The aromatic-flag stays True for all His
  variants. But His being "less aromatic" than benzene is a real
  thing — ring-current intensity is ~30% of benzene per ProCS15
  (Christensen-Sauer-Jensen 2011 ratio).
- **Tyr ring is ring-current intensity** is similar to Phe (because
  the para-OH only weakly perturbs the π system); Trp 6-ring is
  similar to benzene; Trp 5-ring (pyrrole) is around 60% of benzene.
  Per Christensen-Sauer-Jensen 2011 calibrated factors at 1.074 vs
  benzene baseline.
- **Ring-flip dynamics**: Phe/Tyr rings flip on the NMR timescale.
  In Markley `text-3:421-441`, slow-flip allows separate δ1 vs δ2
  shifts; fast-flip averages them. The substrate label is fixed
  (RingPosition is topology); the flip dynamics belong on a per-
  conformation result.

---

## Cross-cutting findings

### The conformation question — where does angle data live?

The user's framing was right: "conformation is not just binary but
angle." The honest split is:

- **Topology (substrate)**: which atoms are in a planar/sp2/ring/methylene
  group, what the canonical IUPAC label is for E/Z and pro-R/S, the
  Markley pseudoatom membership. **All static, all derivable from
  the 3-letter residue code + protonation variant alone.**
- **Conformation (per-frame)**: actual ω angle, χ angles, ring-pucker
  phase, sp2 pyramidalization, ring-flip state classification, planar-
  group dihedral deviation. **Continuous angles, fluctuating each frame.**

The literature evidence:

- ProCS15 fixes ω at 180° in the substrate-side scan (`larsen-2015...text-3:65-67`).
- SHIFTX2 uses ω as a continuous feature of the conformation
  (`han-2011...text-4:147-157` Table 5: ω accounts for 5-14% of the
  shielding signal across atoms).
- Yi-McDermott 2024 finds ψ correlates with H-bond θ, both are
  per-frame conformation features that drive 25 ppm of 15N shift
  variation (yi-mcdermott summary line 7).

So the substrate label of "this atom is in a peptide-amide planar
group with canonical trans (E) configuration" is correct and stable.
The actual ω angle goes on the per-frame data. **Don't put angles
on the substrate.**

### Sources and citations

Local corpus:
- Markley 1998 (the IUPAC NMR nomenclature paper):
  `/shared/2026Thesis/nmr-shielding/references-text/markley-1998-iupac-nmr-nomenclature-recommendations-text-{1..8}.txt`
- ProCS15 (Larsen 2015):
  `/shared/2026Thesis/nmr-shielding/references-text/larsen-2015-procs15-dft-chemical-shift-predictor-text-{1..7}.txt`
- SHIFTX2 (Han 2011):
  `/shared/2026Thesis/nmr-shielding/references-text/han-2011-shiftx2-protein-chemical-shift-prediction-text-{1..5}.txt`
- CH3Shift (Sahakyan 2011):
  `/shared/2026Thesis/nmr-shielding/references-text/sahakyan-2011-methyl-chemical-shifts-proteins-text-{1..6}.txt`
- Sahakyan-Vendruscolo 2013 RNA ring-currents:
  `/shared/2026Thesis/nmr-shielding/references-text/sahakyan-vendruscolo-2013-ring-current-electric-field-contributions-text-{1..4}.txt`
- Yi-McDermott 2024:
  `/shared/2026Thesis/nmr-shielding/references-meta/yi-mcdermott-2024-temperature-shifts-conformational-dynamics-summary-qwen-test.txt`

External authoritative:
- mmCIF dictionary (CCD per-atom and per-bond fields):
  - <https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_atom.pdbx_aromatic_flag.html> (Y/N)
  - <https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_atom.pdbx_stereo_config.html> (R/S/N)
  - <https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_atom.pdbx_leaving_atom_flag.html> (Y/N)
  - <https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_atom.pdbx_ordinal.html>
  - <https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_atom.pdbx_align.html>
  - <https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_bond.value_order.html> (sing/doub/trip/quad/arom/delo/pi/poly)
  - <https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_bond.pdbx_stereo_config.html> (E/Z/N)
- Westbrook 2015 CCD overview: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4393513/>
- BMRB nomenclature pages:
  - <https://bmrb.io/referenc/nomenclature/>
  - <https://bmrb.io/ref_info/atom_nom.tbl> (BMRB/SC/PDB/UCSF/MSI/XPLOR/SYBYL/MIDAS/DIANA cross-table)
  - <https://bmrb.io/ref_info/pseudoatom_nom.txt> (BMRB pseudoatom table)
- RDKit:
  - <https://www.rdkit.org/docs/cppapi/classRDKit_1_1Atom.html> (ChiralType enum)
  - <https://www.rdkit.org/docs/cppapi/classRDKit_1_1RingInfo.html> (ring queries)
  - <https://www.rdkit.org/docs/RDKit_Book.html> (CIP labeler, SSSR vs SymmSSSR)
- Open Babel: <https://open-babel.readthedocs.io/en/latest/Stereochemistry/stereo.html>
  (Tetrahedral / CisTrans / SquarePlanar)
- Craveur 2013 cis-trans isomerization survey: <https://pubmed.ncbi.nlm.nih.gov/23728840/>

### Recommendations beyond what was asked

A few items that the literature suggests would be worth carrying alongside
the eight fields:

1. **Aromatic flag (binary)** — `pdbx_aromatic_flag` from CCD. Binary
   per-atom; cheap; CCD already has it. Useful as a fast filter
   ("is this atom in an aromatic ring") without needing to traverse
   the RingPosition field. SHIFTX2 ring-current calculation uses
   exactly this.

2. **Hybridization (sp/sp2/sp3)** — RDKit and Open Babel both expose
   this; the Markley conventions implicitly use it (sp2 for amides,
   guanidinium, carboxylates, aromatic; sp3 for everything else). Two
   bits per atom. Distinguishes "Cα is sp3" from "C' is sp2" without
   needing to walk the bond graph. Sahakyan 2011's distance-based
   shielding terms explicitly merge atoms by sp2/sp3 hybridization
   (`sahakyan-2011...text-2:101-115`).

3. **Bond type to neighbors (single/double/aromatic)** — `_chem_comp_bond.value_order`
   from CCD: sing/doub/trip/quad/arom/delo/pi/poly. The "delo" value
   (delocalized) is used in CCD for guanidinium and carboxylate where
   the formal Lewis structure has a double bond but the actual electron
   distribution is symmetric. Carrying this as a per-bond-relative-to-
   atom field on the topology lets calculators classify atoms by their
   bonded-electron environment without re-running the perception.

4. **Charge formal state** — for distinguishing Asp vs ASH, etc. Currently
   handled via the protonation_variant_index integer. A per-atom formal
   charge (-1, 0, +1) gives the same information at the atom level
   without needing to look up the variant.

5. **An explicit "is this atom an exchangeable hydrogen" flag** — beyond
   PolarHKind. Many datasets (BMRB, DSS-corrected solution NMR) only
   have observable shifts for non-exchangeable Hs. A binary flag here
   is cheaper than re-deriving from PolarHKind != NotPolar each time.

6. **Atom symmetry equivalence class within the residue** — for a
   methyl, all three Hs are NMR-equivalent at room temperature (Markley
   M-pseudoatom). For a ring-flip-fast Phe, Cδ1 and Cδ2 are equivalent.
   A per-atom integer "equivalence-class index within this residue's
   topology" would let calculators compute averages once per class.

These are all consistent with the substrate-vs-conformation discipline:
all of (1) through (6) are static topology, no angles, no per-frame
state.

### What's NOT settled in the literature

Three areas where the corpus and the external sources disagree or are
silent:

- **Pseudoatom letters beyond Markley's M/Q/R**: the user asked about
  X/Y/S/T. These appear sporadically in legacy NMR tooling (CARA,
  XPLOR, DIANA) but are not IUPAC. They are NOT in Markley 1998 Table 1.
  Recommendation: do not adopt them; follow Markley.
- **Per-residue pro-R/pro-S table**: Markley 1998 Figure 1 marks (R)
  on specific atoms but the marks were OCR-recovered with some
  garbling in our text-2 chunk. Independent verification against
  Hanson 1966 (CIP rule statement) or the BMRB nomenclature page
  (<https://bmrb.io/referenc/nomenclature/>) is wise per residue.
  This is the named-debt-1 from `project_iupac_topology_landed_20260426`
  memory.
- **Ring-position labels**: ipso/ortho/meta/para is well-defined for
  benzene-like Phe/Tyr but the literature does NOT have a universally
  adopted vocabulary for the 5-ring (His imidazole) or fused rings
  (Trp indole). The PyrroleAlpha/PyrroleBeta/BridgeFusion proposal
  here is synthesized from chemistry conventions, not quoted from
  Markley or any predictor paper. **Use with caution**; it's a
  defensible choice but not authoritative.

---

## Summary table

| Field | Authority | Convention | Granularity choice | Notes |
|-------|-----------|------------|---------------------|-------|
| Locant | Markley 1998 / IUPAC 1969 | Greek letter outward from Cα | 7-value enum + None | Stable; both Greek and Roman are presentation only |
| BranchIndex | Markley 1998 Fig 1 caption | CIP-clockwise rule | (outer, inner) two-level | Arg η has both; everything else just outer |
| DiastereotopicIndex | Markley 1998 Fig 1 (R) marks | IUPAC numeric 2/3 | 3-value enum | Use per-residue table |
| ProchiralStereo | CIP 1966 + Markley 1998 | R/S | 4-value enum {N,R,S,Unassigned} | Maps to CCD pdbx_stereo_config |
| PlanarStereo | Markley 1998 + CCD | E/Z + group kind | Split: PlanarGroupKind enum + PlanarStereo enum (E/Z/None) | Angle goes on conformation |
| PseudoatomClass | Markley 1998 Table 1 (Wüthrich 1983) | M / Q / R + locant | Membership struct | Only IUPAC letters; X/Y/S/T are not IUPAC |
| PolarHKind | Synthesis from chem groups | Functional-group | 12-value enum | Coupled to protonation variant |
| RingPosition | Synthesis from chem topology | Ring system + position | RingSystemKind + RingPositionLabel + ring_size + aromatic + n_het | Trp atoms can have two rings |

The fields above use roughly 15-20 bytes of substrate per atom. For
the standard 20 amino acids that's about 200-400 bytes per residue
table. Insignificant against the kernel-output and conformation data
sizes. Worth carrying.
